//! Memoization of the pure combinatorial kernels.
//!
//! Every table here caches the result of a **pure function** of its key, so the
//! library stays referentially transparent and thread-safe: two threads racing
//! on the same key compute the same value, and the worst case is duplicated
//! work, never a wrong answer. That is categorically different from
//! Symmetrica's mutable global scratch state — nothing here is observable in
//! the results.
//!
//! Every table also **accounts for its bytes**: each insert adds the heap
//! behind the key and value ([`HeapSize`]) and the bucket array to a counter
//! the table owns, and [`cache_stats`] reads the counters without taking a
//! lock. The estimate is calibrated against the allocator in
//! `tests/memory.rs`, so a caller reading it is reading bytes and not a guess.
//!
//! Caches are unbounded by default, which suits interactive research (degrees
//! stay modest and reuse is high). [`clear_caches`] releases them if a
//! long-running session wants the memory back, and [`set_cache_budget`] holds
//! them under a byte count from then on, clearing whole tables in the tier
//! order `docs/plans/cache-budget.md` fixes.

use std::any::{Any, TypeId};
use std::collections::{BTreeMap, HashMap};
use std::hash::{BuildHasher, Hash};
use std::sync::atomic::{AtomicUsize, Ordering::Relaxed};
use std::sync::{Arc, OnceLock, RwLock, RwLockReadGuard, RwLockWriteGuard, TryLockError};

use crate::bh::Rat;
use crate::coeff::Ring;
use crate::fasthash;
use crate::frac::Frac;
use crate::guard::GuardedRat;
use crate::partition::{partitions_of, Partition};
use crate::qt::QtPoly;
use crate::sym::{PowerSum, Schur, SymFn};

/// Bytes a value owns on the heap, beyond its own `size_of`.
///
/// What a table charges for an entry is `size_of` the key and value plus this
/// for each, so an inline type answers zero and a `Vec` answers its capacity
/// times its element size plus what the elements own. Capacity, not length:
/// the allocator holds the capacity.
pub(crate) trait HeapSize {
    fn heap_bytes(&self) -> usize;
}

macro_rules! inline_heap {
    ($($t:ty),*) => {$(
        impl HeapSize for $t {
            #[inline]
            fn heap_bytes(&self) -> usize {
                0
            }
        }
    )*};
}

inline_heap!(u32, u64, u128, i64, i128, TypeId, GuardedRat);

impl<A: HeapSize, B: HeapSize> HeapSize for (A, B) {
    #[inline]
    fn heap_bytes(&self) -> usize {
        self.0.heap_bytes() + self.1.heap_bytes()
    }
}

impl<A: HeapSize, B: HeapSize, C: HeapSize> HeapSize for (A, B, C) {
    #[inline]
    fn heap_bytes(&self) -> usize {
        self.0.heap_bytes() + self.1.heap_bytes() + self.2.heap_bytes()
    }
}

impl<T: HeapSize> HeapSize for Vec<T> {
    fn heap_bytes(&self) -> usize {
        self.capacity() * std::mem::size_of::<T>()
            + self.iter().map(HeapSize::heap_bytes).sum::<usize>()
    }
}

impl<T: HeapSize> HeapSize for Arc<T> {
    fn heap_bytes(&self) -> usize {
        // Two reference counts precede the value in an `ArcInner`.
        2 * std::mem::size_of::<usize>() + std::mem::size_of::<T>() + (**self).heap_bytes()
    }
}

impl<K: HeapSize, V: HeapSize> HeapSize for BTreeMap<K, V> {
    fn heap_bytes(&self) -> usize {
        btree_node_bytes::<K, V>(self.len())
            + self
                .iter()
                .map(|(k, v)| k.heap_bytes() + v.heap_bytes())
                .sum::<usize>()
    }
}

/// The node storage of a B-tree holding `len` entries.
///
/// A leaf holds up to 11 entries beside a parent pointer and two `u16`s, and
/// an internal node adds 12 child pointers; a tree built by insertion runs
/// about three-quarters full. `len / 8` leaves of the full size is that
/// estimate, and `tests/memory.rs` is what holds it to the allocator.
fn btree_node_bytes<K, V>(len: usize) -> usize {
    let leaf = 11 * (std::mem::size_of::<K>() + std::mem::size_of::<V>()) + 16;
    len.div_ceil(8) * leaf
}

/// The heap behind one coefficient, which is zero for every fixed-width ring
/// and the digit vectors for a bignum one.
///
/// Answered by type inspection rather than by a bound on [`Ring`], because a
/// bound would leak this crate-private trait into the public signatures of
/// every generic entry point above a cache. The inspection runs once per
/// coefficient at insert time and never on the read path.
pub(crate) fn coeff_heap_bytes<C: 'static>(c: &C) -> usize {
    let any: &dyn Any = c;
    #[cfg(feature = "bignum")]
    {
        if let Some(b) = any.downcast_ref::<num_bigint::BigInt>() {
            return bigint_heap_bytes(b);
        }
        if let Some(r) = any.downcast_ref::<num_rational::BigRational>() {
            return bigint_heap_bytes(r.numer()) + bigint_heap_bytes(r.denom());
        }
    }
    let _ = any;
    0
}

/// A `BigInt` holds its magnitude as 64-bit digits in a `Vec`; the capacity is
/// not observable, so the length is what is charged.
#[cfg(feature = "bignum")]
fn bigint_heap_bytes(b: &num_bigint::BigInt) -> usize {
    usize::try_from(b.bits().div_ceil(64) * 8).unwrap_or(usize::MAX)
}

impl HeapSize for Schur<QtPoly<i128>> {
    fn heap_bytes(&self) -> usize {
        self.terms().heap_bytes()
    }
}

impl HeapSize for PowerSum<GuardedRat> {
    fn heap_bytes(&self) -> usize {
        self.terms().heap_bytes()
    }
}

impl HeapSize for Arc<dyn Any + Send + Sync> {
    /// Unknowable through the erased type; [`transition_cached`] measures the
    /// concrete table before erasing it and charges that instead.
    fn heap_bytes(&self) -> usize {
        0
    }
}

/// One memo table: the map, and the bytes it holds.
///
/// `heap` is the sum of [`HeapSize::heap_bytes`] over every key and value
/// present; `slots` is the bucket array at the map's current capacity. Both
/// are read without the lock, which is what lets [`cache_stats`] and the
/// budget check run on a path that never waits on a reader.
pub(crate) struct Table<K, V, S = std::hash::RandomState> {
    map: OnceLock<RwLock<HashMap<K, V, S>>>,
    heap: AtomicUsize,
    slots: AtomicUsize,
}

impl<K, V, S: Default + BuildHasher> Table<K, V, S> {
    const fn new() -> Self {
        Table {
            map: OnceLock::new(),
            heap: AtomicUsize::new(0),
            slots: AtomicUsize::new(0),
        }
    }

    fn map(&self) -> &RwLock<HashMap<K, V, S>> {
        self.map
            .get_or_init(|| RwLock::new(HashMap::with_hasher(S::default())))
    }

    fn bytes(&self) -> usize {
        self.heap.load(Relaxed) + self.slots.load(Relaxed)
    }

    fn entries(&self) -> usize {
        rd(self).len()
    }

    /// The bucket array behind `capacity` entries: hashbrown keeps its buckets
    /// a power of two and fills them to 7/8, at one control byte per bucket
    /// beside the entry itself.
    fn slot_bytes(capacity: usize) -> usize {
        if capacity == 0 {
            return 0;
        }
        let buckets = (capacity * 8).div_ceil(7).next_power_of_two();
        buckets * (std::mem::size_of::<(K, V)>() + 1)
    }

    /// Record what an insert that grew the map to `capacity` and added
    /// `heap` bytes changed.
    fn note_insert(&self, heap: usize, capacity: usize) {
        self.heap.fetch_add(heap, Relaxed);
        self.slots.store(Self::slot_bytes(capacity), Relaxed);
    }

    /// Insert `key → value`, charging `heap` for the pair.
    ///
    /// A key already present is replaced and its bytes are not recovered;
    /// that happens only when two threads miss on the same key, and the
    /// overcharge is one entry.
    fn insert(&self, key: K, value: V, heap: usize)
    where
        K: Eq + Hash,
    {
        let mut m = wr(self);
        m.insert(key, value);
        let capacity = m.capacity();
        drop(m);
        self.note_insert(heap, capacity);
        enforce_budget();
    }

    /// Drop every entry and the bucket array with it.
    ///
    /// A fresh map rather than `HashMap::clear`, which keeps the buckets: a
    /// cleared character memo would otherwise still hold its array at the
    /// size the largest run reached.
    fn clear_locked(&self, m: &mut HashMap<K, V, S>) {
        *m = HashMap::with_hasher(S::default());
        self.heap.store(0, Relaxed);
        self.slots.store(0, Relaxed);
    }

    /// Clear, waiting for the lock.
    fn clear_blocking(&self) {
        let mut m = wr(self);
        self.clear_locked(&mut m);
    }

    /// Clear if the lock is free right now; `false` if a reader or writer
    /// holds it. A character recursion holds its mask table's read guard for
    /// its whole run, and a caller freeing memory must not wait behind it.
    fn clear(&self) -> bool {
        match self.map().try_write() {
            Ok(mut m) => {
                self.clear_locked(&mut m);
                true
            }
            Err(TryLockError::Poisoned(p)) => {
                self.clear_locked(&mut p.into_inner());
                true
            }
            Err(TryLockError::WouldBlock) => false,
        }
    }
}

/// Read guard, ignoring poisoning; [`wr`] is the same for writes.
///
/// A poisoned lock means some thread panicked while holding a guard — not that
/// the table is unsound. Every value here is a pure function of its key
/// (module doc), and `compute` runs *outside* the guard, so a panicking
/// computation cannot leave a half-built entry behind; what a panic can leave
/// is the lock flag. Propagating that flag would convert one thread's failure
/// into a permanent crash of every cached path in the process, which is a
/// strictly worse outcome and not a correctness one — a caller who overflowed
/// an `i128` and caught it would find the library dead afterwards. So the flag
/// is cleared and the table used (`docs/policies/failure.md`, R2: a panic is
/// for a violated contract, and a neighbor's panic is not this call's
/// contract).
fn rd<K, V, S: Default + BuildHasher>(t: &Table<K, V, S>) -> RwLockReadGuard<'_, HashMap<K, V, S>> {
    t.map().read().unwrap_or_else(|e| e.into_inner())
}

fn wr<K, V, S: Default + BuildHasher>(
    t: &Table<K, V, S>,
) -> RwLockWriteGuard<'_, HashMap<K, V, S>> {
    t.map().write().unwrap_or_else(|e| e.into_inner())
}

/// Look up `key`, computing and inserting it on a miss.
///
/// The read guard is released before `compute` runs, so a `compute` that
/// recurses back into the same table cannot deadlock.
fn lookup<K, V>(table: &Table<K, V>, key: &K, compute: impl FnOnce() -> V) -> V
where
    K: Eq + Hash + Clone + HeapSize,
    V: Clone + HeapSize,
{
    if let Some(v) = rd(table).get(key) {
        return v.clone();
    }
    let v = compute();
    let key = key.clone();
    let heap = key.heap_bytes() + v.heap_bytes();
    table.insert(key, v.clone(), heap);
    v
}

/// What [`cache_stats`] reports about one table.
#[derive(Clone, Copy)]
struct TableInfo {
    name: &'static str,
    tier: u8,
    bytes: fn() -> usize,
    entries: fn() -> usize,
    clear_blocking: fn(),
    clear: fn() -> bool,
}

macro_rules! hasher {
    () => {
        std::hash::RandomState
    };
    ($h:ty) => {
        $h
    };
}

/// Declares every table and the registry that lists them, in one place so a
/// table cannot exist without a row in [`cache_stats`] and a line in
/// [`clear_caches`].
///
/// The tier is the eviction order `docs/plans/cache-budget.md` fixes: 0 is
/// never cleared by a budget, 2 goes first, then 3, then 1. Declaration order
/// is tier, then name, and is the order [`cache_stats`] reports.
macro_rules! tables {
    ($($tier:literal $label:literal $name:ident : $key:ty => $val:ty $(, $hasher:ty)?;)*) => {
        $(
            fn $name() -> &'static Table<$key, $val, hasher!($($hasher)?)> {
                static T: Table<$key, $val, hasher!($($hasher)?)> = Table::new();
                &T
            }
        )*
        static REGISTRY: &[TableInfo] = &[$(TableInfo {
            name: $label,
            tier: $tier,
            bytes: || $name().bytes(),
            entries: || $name().entries(),
            clear_blocking: || $name().clear_blocking(),
            clear: || $name().clear(),
        }),*];
    };
}

tables! {
    0 "flip_rows" flip_row_table: u32 => Arc<Vec<(Partition, i64)>>;
    0 "lex_partitions" lex_parts_table: u32 => Arc<Vec<Partition>>;
    0 "partitions" partitions_table: u32 => Arc<Vec<Partition>>;
    0 "power_in_h_rows" power_in_h_row_table: u32 => Arc<Vec<(Partition, i128)>>;
    1 "bh_ell" bh_ell_table: (Partition, Partition) => Rat<i128>;
    1 "bh_pieri" bh_pieri_table: (Partition, Partition) => Rat<i128>;
    1 "character" character_table: (Partition, Partition) => i128;
    1 "character_masks" character_mask_table: (u64, u64) => i128,
        std::hash::BuildHasherDefault<fasthash::MixHasher>;
    1 "character_wide_masks" character_wide_mask_table: (u128, u128) => i128,
        std::hash::BuildHasherDefault<fasthash::MixHasher>;
    1 "kostka" kostka_table: (Partition, Partition) => u128;
    1 "lr" lr_table: (Partition, Partition, Partition) => u128;
    2 "bold_p" bold_p_table: Partition => Arc<PowerSum<GuardedRat>>;
    2 "ht_to_st" ht_to_st_table: Partition => Arc<Vec<(Partition, i128)>>;
    2 "htilde" htilde_table: u32 => Arc<Vec<(Partition, Schur<QtPoly<i128>>)>>;
    2 "inverse_kostka_rows" inverse_kostka_row_table: Partition => Arc<Vec<i128>>;
    2 "jt_rows" jt_row_table: Partition => Arc<Vec<(Partition, i64)>>;
    2 "products" product_table: (Partition, Partition) => Arc<Vec<(Partition, u128)>>;
    2 "reduced_kronecker" reduced_kronecker_table: (Partition, Partition)
        => Arc<Vec<(Partition, i128)>>;
    2 "schur_to_st" schur_to_st_table: Partition => Arc<Vec<(Partition, i128)>>;
    2 "skews" skew_table: (Partition, Partition) => Arc<Vec<(Partition, u128)>>;
    2 "st_to_ht" st_to_ht_table: Partition => Arc<Vec<(Partition, i128)>>;
    2 "st_to_schur" st_to_schur_table: Partition => Arc<Vec<(Partition, i128)>>;
    3 "transitions" transition_store: (TypeId, u32) => Arc<dyn Any + Send + Sync>;
}

/// One table's share of the memory the caches hold, as [`cache_stats`]
/// reports it.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct CacheStat {
    /// The table's name, fixed for the life of the crate: `"products"`,
    /// `"skews"`, `"character_masks"`, and so on.
    pub name: &'static str,
    /// The eviction tier `docs/plans/cache-budget.md` assigns the table; 0 is
    /// structural and is never evicted by a budget.
    pub tier: u8,
    /// Entries held right now.
    pub entries: usize,
    /// Bytes held right now, counting the bucket array and everything the
    /// keys and values own on the heap.
    pub bytes: usize,
}

/// What every memo table holds right now, one row per table.
///
/// Rows come in a fixed order — by tier and then by name — and every table
/// has a row whether or not it holds anything, so a caller can diff two
/// readings. Bytes are the crate's own accounting, calibrated against the
/// allocator in `tests/memory.rs`; reading them takes no lock and stalls no
/// computation. The sum over the rows is what a long session has kept, and
/// [`clear_caches`] is what returns it.
///
/// ```
/// use symfn::{cache_stats, clear_caches, partitions_of};
///
/// clear_caches();
/// let before = cache_stats();
/// assert!(before.iter().all(|row| row.entries == 0 && row.bytes == 0));
///
/// // `partitions_of` is uncached; the cached path is what the kernels use.
/// let _ = symfn::kostka(&partitions_of(4)[1], &partitions_of(4)[2]);
/// let after = cache_stats();
/// let kostka = after.iter().find(|row| row.name == "kostka").unwrap();
/// assert_eq!(kostka.entries, 1);
/// assert!(kostka.bytes > 0);
/// assert_eq!(after.len(), before.len());
/// ```
pub fn cache_stats() -> Vec<CacheStat> {
    REGISTRY
        .iter()
        .map(|t| CacheStat {
            name: t.name,
            tier: t.tier,
            entries: (t.entries)(),
            bytes: (t.bytes)(),
        })
        .collect()
}

/// The byte budget, with `usize::MAX` for unbounded — the default, so no
/// Rust caller sees eviction unless it asked for it.
static BUDGET: AtomicUsize = AtomicUsize::new(usize::MAX);

/// The order a budget clears tiers in: the large, cheaply rebuilt expansions
/// first, then the whole-degree matrices, then the scalar tables that make
/// characters and Kostka numbers usable at all. Tier 0 is never here.
const EVICTION_ORDER: [u8; 3] = [2, 3, 1];

/// The byte budget the caches are held to, or `None` for unbounded.
pub fn cache_budget() -> Option<usize> {
    match BUDGET.load(Relaxed) {
        usize::MAX => None,
        b => Some(b),
    }
}

/// Hold the caches to `budget` bytes, or lift the bound with `None`.
///
/// The bound is enforced at every insert: when the sum over [`cache_stats`]
/// passes it, whole tables are cleared — the largest first within tier 2,
/// then tier 3, then tier 1, and never tier 0 — until the sum is under. A
/// table another thread is reading at that moment is skipped and tried again
/// at the next insert, so enforcement never waits on a computation. Results
/// are unaffected, because every table holds a pure function of its key;
/// what a tighter budget costs is recomputation, measured in
/// `docs/record/memory.md`.
///
/// `Some(0)` is a legal budget and means every insert clears what it can.
/// Setting a budget below the current holdings takes effect at the next
/// insert, not immediately.
///
/// ```
/// use symfn::{cache_budget, cache_stats, clear_caches, set_cache_budget, SymFn};
///
/// assert_eq!(cache_budget(), None);
/// set_cache_budget(Some(1 << 20));
/// assert_eq!(cache_budget(), Some(1 << 20));
///
/// // A product that fits in a megabyte stays cached; a sweep that does not
/// // is trimmed back under the line at each insert.
/// let s = symfn::Schur::<i64>::monomial(symfn::Partition::new([2, 1]), 1);
/// let _ = &s * &s;
/// let held: usize = cache_stats().iter().map(|r| r.bytes).sum();
/// assert!(held <= 1 << 20);
///
/// set_cache_budget(None);
/// clear_caches();
/// ```
pub fn set_cache_budget(budget: Option<usize>) {
    BUDGET.store(budget.unwrap_or(usize::MAX), Relaxed);
}

/// Bytes held across every table, from the counters alone.
fn total_bytes() -> usize {
    REGISTRY.iter().map(|t| (t.bytes)()).sum()
}

/// Bring the caches back under the budget, if there is one and they are over.
///
/// Called by every insert after its own lock is released, so the loop below
/// can take each table's write lock in turn with nothing else held. Within a
/// tier it clears the table holding the most bytes and re-reads the total,
/// so a single large table goes before many small ones and the loop stops
/// as soon as the sum is under. A table it cannot lock right now is left for
/// the next insert.
fn enforce_budget() {
    let budget = BUDGET.load(Relaxed);
    if budget == usize::MAX || total_bytes() <= budget {
        return;
    }
    for tier in EVICTION_ORDER {
        loop {
            let largest = REGISTRY
                .iter()
                .filter(|t| t.tier == tier)
                .map(|t| (t, (t.bytes)()))
                .filter(|(_, b)| *b > 0)
                .max_by_key(|(_, b)| *b);
            let Some((table, _)) = largest else {
                break;
            };
            if !(table.clear)() {
                // Skipped: another thread holds it. Try the rest of the tier.
                break;
            }
            if total_bytes() <= budget {
                return;
            }
        }
    }
}

/// h_n in the e-basis (equivalently e_n in the h-basis), as integer terms.
///
/// Tier 0 with [`partitions_cached`]: a few hundred bytes per degree, and the
/// row for `n` is built from every row below it, so evicting it would cost
/// the whole ladder.
pub fn flip_row_cached(
    n: u32,
    compute: impl FnOnce() -> Vec<(Partition, i64)>,
) -> Arc<Vec<(Partition, i64)>> {
    lookup(flip_row_table(), &n, || Arc::new(compute()))
}

/// p_n in the h-basis, as integer terms; see [`flip_row_cached`].
pub fn power_in_h_row_cached(
    n: u32,
    compute: impl FnOnce() -> Vec<(Partition, i128)>,
) -> Arc<Vec<(Partition, i128)>> {
    lookup(power_in_h_row_table(), &n, || Arc::new(compute()))
}

/// `H̃_μ` in the Schur basis for a whole degree — the modified (q,t)-Kostka
/// coefficients, cached at `i128` and converted by the caller.
///
/// Cached because the recursion's unit of work is the **degree**, while
/// `qt_kostka` and `qt_kostka_column` are asked for one value or one column.
/// Without this, each of those rebuilds the whole table — the mistake
/// [`kostka_table`](mod@crate::kostka) documents, arrived at from the other
/// direction (`docs/record/qt-kostka.md`).
///
/// `i128` rather than a type parameter, following [`character_cached`]: a
/// `static` cannot be generic, and the values are integers. That is safe here
/// by a bound and not just a measurement. `K̃_{λμ}` has non-negative
/// coefficients (Haiman) summing to `K̃_{λμ}(1,1) = f^λ`, and `Σ_λ (f^λ)² =
/// n!`, so no coefficient exceeds `√(n!)` — past `i128` only around degree 57,
/// which the enumeration never reaches. Measured, they are 9 bits at degree 12
/// and growing about 1 per degree (`docs/record/qt-kostka.md`).
pub fn htilde_cached(
    n: u32,
    compute: impl FnOnce() -> Vec<(Partition, Schur<QtPoly<i128>>)>,
) -> Arc<Vec<(Partition, Schur<QtPoly<i128>>)>> {
    lookup(htilde_table(), &n, || Arc::new(compute()))
}

/// The `s → J` transition matrix of a whole degree, cached **per coefficient
/// ring**.
///
/// Cached because the projection's unit of work is the degree, while
/// [`schur_to_macdonald_j`](crate::schur_to_macdonald_j) is asked for one
/// element: without this, expanding p(n) shapes one at a time rebuilds the
/// table p(n) times, which is measured at 17.4× the whole-degree call at
/// degree 8 (`docs/record/qt-kostka.md`).
///
/// **Keyed by the ring, not stored at one concrete ring** — the opposite of
/// [`htilde_cached`], and the reason is escalation. `schur_in_macdonald_j` and
/// `schur_to_macdonald_j` run a guarded fixed-width pass and re-run over
/// `BigRational` when it reports; a single `i128`-shaped cache would hand the
/// wide pass the narrow pass's values and make the escalation a lie. The
/// values here are also not integers — they are `Frac`s over ℚ(q,t) — so
/// [`htilde_cached`]'s bound argument would not apply anyway.
///
/// **The store is conditional**, which no other table here is except
/// [`bold_p_peek`]'s pair, and for that table's reason. Over a reporting ring a
/// value computed by an overflowing call is garbage, and caching it is worse
/// than recomputing it: a later reader inside a clean [`guarded`] window sees
/// an untouched counter and accepts it. So the counter is read on both sides of
/// `compute` and the entry is stored only if it did not move. Over a ring that
/// cannot report — `Rational`, `BigRational` — the counter never moves and
/// every entry is stored.
///
/// [`guarded`]: crate::guard::guarded
///
/// # Panics
///
/// Panics if an entry stored under `C`'s [`TypeId`] does not hold a table over
/// `C`, which is a bug in this function rather than a reachable state.
pub fn schur_in_j_cached<C: Ring + Send + Sync + 'static>(
    n: u32,
    compute: impl FnOnce() -> Vec<Vec<Frac<C>>>,
) -> Arc<Vec<Vec<Frac<C>>>> {
    transition_cached(n, compute)
}

/// The `m → P` transition of a whole degree, cached on the same terms as
/// [`schur_in_j_cached`] and for the same reason: the unit of work is the
/// degree, while [`monomial_to_macdonald_p`](crate::monomial_to_macdonald_p)
/// is asked for one element (`docs/record/macdonald.md`).
///
/// # Panics
///
/// Panics if an entry stored under this table's [`TypeId`] does not hold a
/// table of that type, which is a bug in [`transition_cached`] rather than a
/// reachable state.
pub fn mac_p_inverse_cached<C: Ring + Send + Sync + 'static>(
    n: u32,
    compute: impl FnOnce() -> Vec<std::collections::BTreeMap<Partition, Frac<C>>>,
) -> Arc<Vec<std::collections::BTreeMap<Partition, Frac<C>>>> {
    transition_cached(n, compute)
}

/// The Jack `m → P` transition of a whole degree, on the same terms as
/// [`mac_p_inverse_cached`]. A separate entry point only because the table is
/// over [`AFrac`](crate::AFrac) rather than [`Frac`]; the store is shared, and
/// the two cannot collide because the key carries the table's type.
///
/// # Panics
///
/// Panics if an entry stored under this table's [`TypeId`] does not hold a
/// table of that type, which is a bug in [`transition_cached`] rather than a
/// reachable state.
pub fn jack_p_inverse_cached<C: crate::coeff::Integral + Send + Sync + 'static>(
    n: u32,
    compute: impl FnOnce() -> Vec<std::collections::BTreeMap<Partition, crate::AFrac<C>>>,
) -> Arc<Vec<std::collections::BTreeMap<Partition, crate::AFrac<C>>>> {
    transition_cached(n, compute)
}

/// One degree's transition table, keyed by its own Rust type and the degree.
///
/// The type parameter carries both the shape of the table and the coefficient
/// ring it is over, so two tables of different shape — and one table over two
/// rings — never share a key. That is what makes the key enough on its own.
fn transition_cached<T: HeapSize + Send + Sync + 'static>(
    n: u32,
    compute: impl FnOnce() -> T,
) -> Arc<T> {
    let key = (TypeId::of::<T>(), n);
    if let Some(v) = rd(transition_store()).get(&key) {
        return Arc::clone(v)
            .downcast::<T>()
            .expect("a transition cache is keyed by the type it stores");
    }
    let before = crate::guard::overflow_count();
    let value = Arc::new(compute());
    if crate::guard::overflow_count() == before {
        let heap = value.heap_bytes();
        transition_store().insert(key, Arc::clone(&value) as Arc<dyn Any + Send + Sync>, heap);
    }
    value
}

/// The Bergeron–Haiman Pieri coefficient `c⁽ʳ⁾_{μν}`, and `L_{μν} = ⟨H̃_μ,
/// h_ν⟩`.
///
/// Cached **across degrees**, which is the point: computing degree `n` needs
/// both at every size below `n`, so a degree-12 run rebuilds most of what a
/// degree-11 run already knew. Sharing them speeds up a walk up the degrees for
/// **no extra memory at all** (`docs/record/llt.md`): the degree-12 call was
/// building that cache inside itself either way, so all the sharing does is
/// stop the smaller degrees rebuilding it.
///
/// It buys nothing for a single cold degree, which is what
/// `bench_qt_kostka.py` measures, so the headline comparison against Sage is
/// unaffected either way.
///
/// `i128` for the same reason as [`htilde_cached`]: a `static` cannot be
/// generic. Safe by measurement rather than by a bound here, since these are
/// intermediate rational functions and not the coefficients Haiman's theorem
/// constrains — 25 bits at degree 12, growing about 3 per degree, so `i128`
/// holds past degree 45 (`docs/record/qt-kostka.md`).
pub fn bh_pieri_cached(
    key: &(Partition, Partition),
    compute: impl FnOnce() -> Rat<i128>,
) -> Rat<i128> {
    lookup(bh_pieri_table(), key, compute)
}

/// `L_{μν}`; see [`bh_pieri_cached`].
pub fn bh_ell_cached(
    key: &(Partition, Partition),
    compute: impl FnOnce() -> Rat<i128>,
) -> Rat<i128> {
    lookup(bh_ell_table(), key, compute)
}

/// The partitions of `n`, shared rather than regenerated.
///
/// Conversions and the coproduct enumerate `partitions_of` inside nested
/// loops, where the same list would otherwise be rebuilt thousands of times.
pub fn partitions_cached(n: u32) -> Arc<Vec<Partition>> {
    lookup(partitions_table(), &n, || Arc::new(partitions_of(n)))
}

/// Memoized χ^λ(μ), keyed on the partitions themselves. The Murnaghan–Nakayama
/// recursion has heavily overlapping subproblems, so this is the difference
/// between exponential and near-linear on repeated use.
///
/// This is the table for the pairs whose β-sets fit neither mask width —
/// past |λ| = 127. Everything below runs on [`character_masks_read`], keyed
/// on β-masks, and never reaches here.
///
/// `compute` returning `None` means the character overflows `i128`. That is
/// **not** cached: it is the absence of a representable value rather than a
/// value, and caching it would make the table unable to distinguish "not yet
/// computed" from "cannot be represented". Overflow is confined to |λ| ≳ 58, so
/// recomputing it is not a hot path.
pub fn character_cached(
    lambda: &Partition,
    mu: &Partition,
    compute: impl FnOnce() -> Option<i128>,
) -> Option<i128> {
    let key = (lambda.clone(), mu.clone());
    if let Some(&v) = rd(character_table()).get(&key) {
        return Some(v);
    }
    let v = compute()?;
    let heap = key.heap_bytes();
    character_table().insert(key, v, heap);
    Some(v)
}

/// A β-mask width that has a character memo of its own: `u64` through degree
/// 63, `u128` through 127 (`character::fits_mask`). Two tables rather than one
/// wide one so the narrow recursion, which is the hot one, hashes two words
/// and not four.
pub(crate) trait MaskKey: crate::convert::Beta {
    fn table(
    ) -> &'static Table<(Self, Self), i128, std::hash::BuildHasherDefault<fasthash::MixHasher>>;
}

impl MaskKey for u64 {
    fn table(
    ) -> &'static Table<(u64, u64), i128, std::hash::BuildHasherDefault<fasthash::MixHasher>> {
        character_mask_table()
    }
}

impl MaskKey for u128 {
    fn table(
    ) -> &'static Table<(u128, u128), i128, std::hash::BuildHasherDefault<fasthash::MixHasher>>
    {
        character_wide_mask_table()
    }
}

/// The β-mask-keyed character memo, read-locked for the duration of one
/// Murnaghan–Nakayama recursion.
///
/// Keys are `(β-mask of λ, β-mask of μ)` in the canonical form
/// `convert::Beta::canonical` fixes, so a value stored by one caller is found
/// by every other. This is the table [`character_cached`] would be if a
/// partition were a word: the recursion looks up once per node, and a
/// `Partition`-keyed table charged two heap clones and a SipHash per lookup
/// (`docs/record/transitions.md`).
///
/// **The guard is held across the whole recursion, and the recursion never
/// writes.** New values go to a local map and land here through
/// [`character_masks_store`] after the guard is dropped — one lock round trip
/// per top-level character rather than one per node. Holding a read guard
/// while asking for another can deadlock against a waiting writer, so a
/// recursion must not call back into anything that takes this lock.
pub(crate) fn character_masks_read<M: MaskKey>(
) -> RwLockReadGuard<'static, fasthash::Map<(M, M), i128>> {
    rd(M::table())
}

/// See [`character_masks_read`]: merge one recursion's new values.
pub(crate) fn character_masks_store<M: MaskKey>(entries: fasthash::Map<(M, M), i128>) {
    let table = M::table();
    let mut t = wr(table);
    t.reserve(entries.len());
    t.extend(entries);
    let capacity = t.capacity();
    drop(t);
    table.note_insert(0, capacity);
    enforce_budget();
}

/// Memoized Kostka number K_{λμ}. A `compute` that declines stores nothing,
/// so a later caller asks again rather than reading a refusal as a value
/// (`docs/policies/failure.md`, R7).
pub fn kostka_cached(
    lambda: &Partition,
    mu: &Partition,
    compute: impl FnOnce() -> Option<u128>,
) -> Option<u128> {
    let table = kostka_table();
    let key = (lambda.clone(), mu.clone());
    if let Some(v) = rd(table).get(&key) {
        return Some(*v);
    }
    let v = compute()?;
    let heap = key.heap_bytes();
    table.insert(key, v, heap);
    Some(v)
}

/// Memoized Littlewood–Richardson coefficient c^λ_{μν}.
pub fn lr_cached(
    lambda: &Partition,
    mu: &Partition,
    nu: &Partition,
    compute: impl FnOnce() -> u128,
) -> u128 {
    lookup(
        lr_table(),
        &(lambda.clone(), mu.clone(), nu.clone()),
        compute,
    )
}

/// The partitions of `n` in **decreasing-lex** order, the order in which the
/// Kostka matrix is upper-unitriangular.
pub fn lex_parts_cached(n: u32) -> Arc<Vec<Partition>> {
    lookup(lex_parts_table(), &n, || {
        let mut v = partitions_of(n);
        v.sort_by(|a, b| b.parts().cmp(a.parts()));
        Arc::new(v)
    })
}

/// Memoized **one row** of the inverse Kostka matrix, indexed against
/// [`lex_parts_cached`].
///
/// Keyed by μ rather than by degree because m → s needs a single row per term.
/// Caching the whole matrix instead meant one conversion paid for p(n)² tableau
/// counts to read p(n) of them.
pub fn inverse_kostka_row_cached(
    mu: &Partition,
    compute: impl FnOnce() -> Vec<i128>,
) -> Arc<Vec<i128>> {
    lookup(inverse_kostka_row_table(), mu, || Arc::new(compute()))
}

/// One Jacobi–Trudi row: `s_λ` expanded in `h` (or, on the conjugate index, in
/// `e`), which is `det(h_{λ_i − i + j})` read as terms.
///
/// Keyed by the index alone, and ring-free, because the determinant is: the
/// coefficient ring enters only when a caller scales the row into its own
/// output. That is what makes one table serve every caller.
///
/// The gap this closes was visible from Sage rather than from here. `s → m`
/// was already memoized through [`kostka_cached`]; `s → h` and `s → e` were
/// not. A workload making thousands of repeated small conversions paid the
/// determinant every time — `sage.combinat.sf.character`'s peel is one, at one
/// conversion per term removed. On a path where a single conversion wins, that
/// lost to Symmetrica by 20x.
pub fn jt_row_cached(
    index: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, i64)>,
) -> Arc<Vec<(Partition, i64)>> {
    lookup(jt_row_table(), index, || Arc::new(compute()))
}

/// `𝐩_γ` (OZ Eq 24) in the power-sum basis — the image of `p_γ` under the map Γ
/// that carries `s_λ` to `s̃_λ`.
///
/// Keyed by γ alone and shared across every λ and every degree, because that is
/// what it is: a function of γ. One `s̃_λ` of degree n needs all p(n) of them,
/// and every other λ of that degree needs the same ones again.
///
/// **Peek and store are separate on purpose**, unlike every other table here.
/// The values are fixed-width rationals that *report* overflow rather than
/// wrapping, and a value computed by a call that overflowed is garbage. Caching
/// it would be worse than recomputing it: the overflow counter is checked
/// around the call that produced it, so a later reader of the poisoned entry
/// would see a clean counter and accept a wrong answer. The caller stores only
/// what it has confirmed clean; see `character_basis::bold_guarded`.
pub fn bold_p_peek(gamma: &Partition) -> Option<Arc<PowerSum<GuardedRat>>> {
    rd(bold_p_table()).get(gamma).cloned()
}

/// See [`bold_p_peek`]. Storing an entry asserts it was computed without
/// overflow.
pub fn bold_p_store(gamma: &Partition, value: PowerSum<GuardedRat>) {
    let value = Arc::new(value);
    let heap = gamma.heap_bytes() + value.heap_bytes();
    bold_p_table().insert(gamma.clone(), value, heap);
}

/// Memoized `s̃_λ` in the Schur basis, and its inverse `s_λ` in the `s̃` basis.
///
/// Cached at `i128` with the generic conversion at the edges — the same shape
/// as [`htilde_cached`], for the same reason: a `static` cannot be generic, and
/// both transitions are integral (OZ Thm 1(2) makes the `s → s̃` direction a
/// matrix of *non-negative* integers, being multiplicities in a restriction).
pub fn st_to_schur_cached(
    lambda: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, i128)>,
) -> Arc<Vec<(Partition, i128)>> {
    lookup(st_to_schur_table(), lambda, || Arc::new(compute()))
}

/// See [`st_to_schur_cached`]; this is the other direction.
pub fn schur_to_st_cached(
    nu: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, i128)>,
) -> Arc<Vec<(Partition, i128)>> {
    lookup(schur_to_st_table(), nu, || Arc::new(compute()))
}

/// Memoized `h̃_μ` in the `s̃` basis, and its inverse `s̃_λ` in the `h̃` basis.
///
/// The same shape as [`st_to_schur_cached`]. `ht_to_st_row` sweeps every
/// partition of every size up to |μ|, and `st_to_ht_row` runs a back
/// substitution that calls it once per pivot. A miss is therefore expensive,
/// and a row is small.
pub fn ht_to_st_cached(
    mu: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, i128)>,
) -> Arc<Vec<(Partition, i128)>> {
    lookup(ht_to_st_table(), mu, || Arc::new(compute()))
}

/// See [`ht_to_st_cached`]; this is the other direction.
pub fn st_to_ht_cached(
    lambda: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, i128)>,
) -> Arc<Vec<(Partition, i128)>> {
    lookup(st_to_ht_table(), lambda, || Arc::new(compute()))
}

/// Memoized `s̃_λ · s̃_μ` — one whole column of reduced Kronecker coefficients.
///
/// Whole expansion rather than single coefficients, following
/// [`product_cached`]: the engine produces every ν at once, so caching per
/// coefficient would recompute the column once per ν asked for.
pub fn reduced_kronecker_cached(
    lambda: &Partition,
    mu: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, i128)>,
) -> Arc<Vec<(Partition, i128)>> {
    lookup(
        reduced_kronecker_table(),
        &(lambda.clone(), mu.clone()),
        || Arc::new(compute()),
    )
}

/// Memoized full expansion of s_μ · s_ν. Caching the whole product (rather than
/// individual coefficients) is what makes repeated multiplication cheap.
pub fn product_cached(
    mu: &Partition,
    nu: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, u128)>,
) -> Arc<Vec<(Partition, u128)>> {
    lookup(product_table(), &(mu.clone(), nu.clone()), || {
        Arc::new(compute())
    })
}

/// Memoized full expansion of the skew Schur function s_{outer/inner}.
///
/// Keyed on the *shape* rather than on (λ, μ, ν), because one traversal of the
/// shape produces every ν at once. A caller sweeping over ν — `skew_schur`, the
/// coproduct, a single `lr_coeff` in a loop — therefore pays for exactly one
/// enumeration instead of p(n) of them.
pub fn skew_cached(
    outer: &Partition,
    inner: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, u128)>,
) -> Arc<Vec<(Partition, u128)>> {
    lookup(skew_table(), &(outer.clone(), inner.clone()), || {
        Arc::new(compute())
    })
}

/// The cached expansion of s_{outer/inner}, if some earlier call computed it.
///
/// A read-only peek — never computes. This is what lets a coefficient query
/// answer from an expansion a previous caller already paid for (e.g. the
/// whole-product expansion behind c^λ_{μν}) instead of starting a fresh
/// traversal of its own.
pub fn skew_cache_peek(
    outer: &Partition,
    inner: &Partition,
) -> Option<Arc<Vec<(Partition, u128)>>> {
    rd(skew_table())
        .get(&(outer.clone(), inner.clone()))
        .cloned()
}

/// Drop every cached table, releasing the memory.
///
/// Every table, tier 0 included, and waiting for each lock in turn: a timing
/// run wants a cold start, and a caller returning memory wants all of it.
/// Every row of [`cache_stats`] reads zero afterwards.
pub fn clear_caches() {
    for t in REGISTRY {
        (t.clear_blocking)();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A table computed while the overflow counter moved is not stored.
    ///
    /// The positive control is in the same test: without it, a `schur_in_j_cached`
    /// that never stored anything would pass the negative half.
    #[test]
    fn a_reported_overflow_is_not_cached() {
        use crate::guard::Guarded;

        let _lock = crate::guard::serial();
        let calls = std::cell::Cell::new(0);
        let clean = || {
            calls.set(calls.get() + 1);
            vec![vec![Frac::<Guarded>::from_poly(QtPoly::zero())]]
        };
        let dirty = || {
            calls.set(calls.get() + 1);
            // Reported, not wrapped — so the table this stands for is garbage.
            let _ = Guarded(i128::MAX).mul(&Guarded(2));
            vec![vec![Frac::<Guarded>::from_poly(QtPoly::zero())]]
        };

        // Degrees no computation reaches, so the entries are this test's alone.
        schur_in_j_cached(u32::MAX, clean);
        schur_in_j_cached(u32::MAX, clean);
        assert_eq!(calls.get(), 1, "a clean table was not cached");

        calls.set(0);
        schur_in_j_cached(u32::MAX - 1, dirty);
        schur_in_j_cached(u32::MAX - 1, dirty);
        assert_eq!(calls.get(), 2, "a table computed past the width was cached");

        wr(transition_store()).clear();
    }

    #[test]
    fn stats_rows_come_by_tier_then_name_and_cover_every_table() {
        let rows = cache_stats();
        assert_eq!(rows.len(), REGISTRY.len());
        let order: Vec<(u8, &str)> = rows.iter().map(|r| (r.tier, r.name)).collect();
        let mut sorted = order.clone();
        sorted.sort();
        assert_eq!(order, sorted, "registry order is tier, then name");
        assert!(rows.iter().any(|r| r.name == "flip_rows" && r.tier == 0));
        assert!(rows.iter().any(|r| r.name == "transitions" && r.tier == 3));
    }

    // A table of the test's own: the shared ones are process-wide and the
    // test threads race on them, so an absolute reading there proves nothing.
    #[test]
    fn an_insert_charges_the_key_and_value_and_the_buckets() {
        static OWN: Table<(Partition, Partition), u128> = Table::new();
        let lam = Partition::new([13, 11, 9, 7, 5, 3, 1]);
        let mu = Partition::new([7, 7, 7, 7, 7, 7, 7]);
        assert_eq!((OWN.entries(), OWN.bytes()), (0, 0));
        lookup(&OWN, &(lam.clone(), mu.clone()), || 12345);
        let keys = lam.heap_bytes() + mu.heap_bytes();
        let slot = std::mem::size_of::<((Partition, Partition), u128)>() + 1;
        assert_eq!(OWN.entries(), 1);
        assert!(
            OWN.bytes() >= keys + slot,
            "{} bytes charged for {keys} of keys and at least one {slot}-byte slot",
            OWN.bytes()
        );
        lookup(&OWN, &(lam.clone(), mu.clone()), || 0);
        assert_eq!(OWN.entries(), 1, "a hit inserts nothing");
        OWN.clear_blocking();
        assert_eq!(
            (OWN.entries(), OWN.bytes()),
            (0, 0),
            "clear releases the buckets too"
        );
    }

    #[test]
    fn partitions_are_shared_and_correct() {
        let a = partitions_cached(6);
        let b = partitions_cached(6);
        assert_eq!(a.len(), 11); // p(6) = 11
        assert!(
            Arc::ptr_eq(&a, &b),
            "second call should reuse the cached Arc"
        );
    }

    #[test]
    fn lookup_computes_once_then_reuses() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        // A key no other test touches: `clear_caches` is global, and the test
        // harness runs tests in parallel, so isolation comes from the key.
        let calls = AtomicUsize::new(0);
        let lam = Partition::new([9, 7, 5, 3, 1]);
        let mu = Partition::new([5, 5, 5, 5, 5]);
        for _ in 0..5 {
            kostka_cached(&lam, &mu, || {
                calls.fetch_add(1, Ordering::SeqCst);
                Some(12345)
            });
        }
        assert_eq!(
            calls.load(Ordering::SeqCst),
            1,
            "should compute exactly once"
        );
    }
}
