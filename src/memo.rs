//! Memoization of the pure combinatorial kernels.
//!
//! Every table here caches the result of a **pure function** of its key, so the
//! library stays referentially transparent and thread-safe: two threads racing
//! on the same key compute the same value, and the worst case is duplicated
//! work, never a wrong answer. That is categorically different from
//! Symmetrica's mutable global scratch state — nothing here is observable in
//! the results.
//!
//! Caches are unbounded, which suits interactive research (degrees stay modest
//! and reuse is high). [`clear_caches`] releases them if a long-running session
//! wants the memory back.

use std::any::{Any, TypeId};
use std::collections::HashMap;
use std::hash::Hash;
use std::sync::{Arc, OnceLock, RwLock, RwLockReadGuard};

use crate::bh::Rat;
use crate::coeff::Ring;
use crate::fasthash;
use crate::frac::Frac;
use crate::guard::GuardedRat;
use crate::partition::{partitions_of, Partition};
use crate::qt::QtPoly;
use crate::sym::{PowerSum, Schur};

type Table<K, V, S = std::hash::RandomState> = RwLock<HashMap<K, V, S>>;

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
fn rd<K, V, S>(t: &Table<K, V, S>) -> RwLockReadGuard<'_, HashMap<K, V, S>> {
    t.read().unwrap_or_else(|e| e.into_inner())
}

fn wr<K, V, S>(t: &Table<K, V, S>) -> std::sync::RwLockWriteGuard<'_, HashMap<K, V, S>> {
    t.write().unwrap_or_else(|e| e.into_inner())
}

/// Look up `key`, computing and inserting it on a miss.
///
/// The read guard is released before `compute` runs, so a `compute` that
/// recurses back into the same table cannot deadlock.
fn lookup<K, V>(table: &Table<K, V>, key: &K, compute: impl FnOnce() -> V) -> V
where
    K: Eq + Hash + Clone,
    V: Clone,
{
    if let Some(v) = rd(table).get(key) {
        return v.clone();
    }
    let v = compute();
    wr(table).insert(key.clone(), v.clone());
    v
}

macro_rules! table {
    ($name:ident, $key:ty, $val:ty) => {
        table!($name, $key, $val, std::hash::RandomState);
    };
    ($name:ident, $key:ty, $val:ty, $hasher:ty) => {
        fn $name() -> &'static Table<$key, $val, $hasher> {
            static T: OnceLock<Table<$key, $val, $hasher>> = OnceLock::new();
            T.get_or_init(|| RwLock::new(HashMap::with_hasher(<$hasher>::default())))
        }
    };
}

table!(partitions_table, u32, Arc<Vec<Partition>>);
table!(transition_store, (TypeId, u32), Arc<dyn Any + Send + Sync>);
table!(character_table, (Partition, Partition), i128);
table!(
    character_mask_table,
    (u64, u64),
    i128,
    std::hash::BuildHasherDefault<fasthash::MixHasher>
);
table!(
    character_wide_mask_table,
    (u128, u128),
    i128,
    std::hash::BuildHasherDefault<fasthash::MixHasher>
);
table!(kostka_table, (Partition, Partition), u128);
table!(lr_table, (Partition, Partition, Partition), u128);
table!(lex_parts_table, u32, Arc<Vec<Partition>>);
table!(inverse_kostka_row_table, Partition, Arc<Vec<i128>>);
table!(jt_row_table, Partition, Arc<Vec<(Partition, i64)>>);
table!(
    product_table,
    (Partition, Partition),
    Arc<Vec<(Partition, u128)>>
);
table!(
    skew_table,
    (Partition, Partition),
    Arc<Vec<(Partition, u128)>>
);
table!(
    htilde_table,
    u32,
    Arc<Vec<(Partition, Schur<QtPoly<i128>>)>>
);
table!(bh_pieri_table, (Partition, Partition), Rat<i128>);
table!(bh_ell_table, (Partition, Partition), Rat<i128>);
table!(bold_p_table, Partition, Arc<PowerSum<GuardedRat>>);
table!(st_to_schur_table, Partition, Arc<Vec<(Partition, i128)>>);
table!(schur_to_st_table, Partition, Arc<Vec<(Partition, i128)>>);
table!(ht_to_st_table, Partition, Arc<Vec<(Partition, i128)>>);
table!(st_to_ht_table, Partition, Arc<Vec<(Partition, i128)>>);
table!(
    reduced_kronecker_table,
    (Partition, Partition),
    Arc<Vec<(Partition, i128)>>
);

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
fn transition_cached<T: Send + Sync + 'static>(n: u32, compute: impl FnOnce() -> T) -> Arc<T> {
    let key = (TypeId::of::<T>(), n);
    if let Some(v) = rd(transition_store()).get(&key) {
        return Arc::clone(v)
            .downcast::<T>()
            .expect("a transition cache is keyed by the type it stores");
    }
    let before = crate::guard::overflow_count();
    let value = Arc::new(compute());
    if crate::guard::overflow_count() == before {
        wr(transition_store()).insert(key, Arc::clone(&value) as Arc<dyn Any + Send + Sync>);
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
    wr(character_table()).insert(key, v);
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
    let mut t = wr(M::table());
    t.reserve(entries.len());
    t.extend(entries);
}

/// Memoized Kostka number K_{λμ}.
pub fn kostka_cached(lambda: &Partition, mu: &Partition, compute: impl FnOnce() -> u128) -> u128 {
    lookup(kostka_table(), &(lambda.clone(), mu.clone()), compute)
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
    wr(bold_p_table()).insert(gamma.clone(), Arc::new(value));
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
pub fn clear_caches() {
    wr(htilde_table()).clear();
    wr(transition_store()).clear();
    wr(bh_pieri_table()).clear();
    wr(bh_ell_table()).clear();
    wr(partitions_table()).clear();
    wr(character_table()).clear();
    wr(character_mask_table()).clear();
    wr(character_wide_mask_table()).clear();
    wr(kostka_table()).clear();
    wr(lr_table()).clear();
    wr(lex_parts_table()).clear();
    wr(inverse_kostka_row_table()).clear();
    wr(jt_row_table()).clear();
    wr(product_table()).clear();
    wr(skew_table()).clear();
    wr(bold_p_table()).clear();
    wr(st_to_schur_table()).clear();
    wr(ht_to_st_table()).clear();
    wr(st_to_ht_table()).clear();
    wr(schur_to_st_table()).clear();
    wr(reduced_kronecker_table()).clear();
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
                12345
            });
        }
        assert_eq!(
            calls.load(Ordering::SeqCst),
            1,
            "should compute exactly once"
        );
    }
}
