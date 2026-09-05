//! Whole-shape Littlewood–Richardson expansion.
//!
//! The other two backends answer one question at a time: "how many LR tableaux
//! of shape λ/μ have content ν?". Getting a whole expansion out of them costs
//! one independent search per candidate answer. This backend answers the
//! *shape-level* question instead —
//!
//! ```text
//!   s_{λ/μ} = Σ_ν c^λ_{μν} · s_ν
//! ```
//!
//! — by walking the skew diagram once and letting each completed filling report
//! its own content. Since the content of an LR tableau is automatically a
//! partition, every ν with a nonzero coefficient is discovered by the walk; no
//! list of candidates is ever consulted, so the cost of an expansion is
//! independent of `p(n)`.
//!
//! One enumeration is necessary but nowhere near sufficient: a shape can have
//! astronomically more LR tableaux than it has distinct contents — the shape
//! behind `s_{8,7,6,5,4,3}²` has 2.1 × 10⁸ of them across 164 037 terms — so
//! visiting tableaux one at a time would be hopeless. The traversal therefore
//! advances a *layer* of merged partial fillings rather than a stack of
//! individual ones: two partial fillings that agree on the row above and on the
//! content so far are interchangeable, so they collapse into one weighted
//! state. That is where the asymptotic win lives; the single enumeration only
//! removes the `p(n)` factor on top of it.
//!
//! Sage computes the same expansion as `s[outer].skew_by(s[inner])` in the
//! Schur basis (`scripts/compare_sage.py`).
//!
//! Products reduce to the same primitive. See [`SkewLr::schur_product`].

// Key packing. Every element serialized into a `Key` is at most the cell count
// of the shape (see `elem_width`, which chooses the byte width from exactly
// that bound), so each narrowing is inside the width the same function picked;
// a `PackedKey` bit position is a part index or a column, both bounded by the
// shape's `rows + width`, which `LayerKey::fits` has already checked.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
// The crate denies `unsafe_code`; this module carries one block, the
// `&[u8]` → `&KeyBytes` cast in `KeyBytes::new`, with its reason beside it.
#![allow(unsafe_code)]

use std::borrow::Borrow;
use std::collections::hash_map::Entry;
use std::hash::{Hash, Hasher};
use std::sync::{Arc, OnceLock};

use crate::fasthash::Map;
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};

use crate::lr::LrBackend;
use crate::memo::skew_cached;
use crate::partition::Partition;

/// Backend that expands an entire skew shape in one traversal.
///
/// Implements [`LrBackend`] on top of [`expand_skew`]; both `lr_coeff` and
/// `schur_product` are a memoized expansion plus a lookup, never a sweep over
/// candidate partitions.
#[derive(Clone, Copy, Debug, Default)]
pub struct SkewLr;

/// The full Schur expansion of the skew Schur function s_{outer/inner}.
///
/// Returns every ν with c^outer_{inner,ν} ≠ 0, paired with that coefficient,
/// sorted by ν and free of zero terms. Yields the empty vector when
/// `inner ⊄ outer`, and `[(∅, 1)]` when the two shapes are equal.
///
/// Memoized on the shape (see `memo::skew_cached`), so a
/// caller sweeping many ν against one (outer, inner) pays for one traversal.
///
/// # Panics
///
/// Panics if a layer multiplicity exceeds `u128`. The accumulator retries the
/// whole traversal in `u128` when `u64` overflows, and refuses loudly above
/// that. Multiplicities are tableau counts and dwarf the final coefficients, so
/// the refusal is a wall on the traversal rather than on the answer. Memory
/// arrives first on every shape measured: `[24,20,16,12]²` has 5 313 471 terms
/// (`docs/record/littlewood-richardson.md`), and its expansion is already past
/// what fits.
pub fn expand_skew(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
    (*expand_skew_shared(outer, inner)).clone()
}

/// [`expand_skew`] without the copy: the memoized expansion itself.
///
/// The memo holds an `Arc<Vec<_>>` and `expand_skew` hands back a deep clone of
/// it, which for a caller that only reads is the whole expansion stored twice —
/// **one allocation per term**, since every `Partition` owns a `Vec`. On
/// `[8,7,6,5,4,3]²` that is 164 037 allocations and 14.1 MB per call, on top of
/// an identical 14.1 MB already sitting in the cache. On `[24,20,16,12]²`, at
/// 5.3M terms, it is a growing share of peak RSS and the reason a shape that
/// fits can still fail to run.
///
/// Every caller inside the crate iterates and drops, so they take this. The
/// owned version stays for callers that want to mutate or keep the vector past
/// a [`clear_caches`](crate::clear_caches).
///
/// # Panics
///
/// Panics where [`expand_skew`] does: on a layer multiplicity past `u128`.
pub fn expand_skew_shared(outer: &Partition, inner: &Partition) -> Arc<Vec<(Partition, u128)>> {
    if !outer.contains(inner) {
        return Arc::new(Vec::new());
    }
    skew_cached(outer, inner, || expand_skew_uncached(outer, inner))
}

/// A partially-filled diagram, reduced to what the rest of the fill can see,
/// as the layer map stores it.
///
/// Two partial fillings that agree on this state are interchangeable:
/// everything still to be decided depends on the row above (column strictness)
/// and on the content accumulated so far (the ballot condition), and on nothing
/// else. So they are merged and carry a multiplicity, which is what keeps the
/// traversal from degenerating into an enumeration of individual tableaux.
///
/// The state is the pair `(content, above)`: `content[v - 1]` is the number of
/// `v`s placed so far (no trailing zeros), and `above` holds the previous
/// row's entries, clipped to the columns the next row overlaps. Peak memory
/// is (number of live states) × (bytes per state), and the traversal is bound
/// by how much of the layer fits in cache, so the serialization is the
/// second factor and it is chosen per shape:
///
/// * [`PackedKey`] holds each half as a lattice-path bitmap in one machine
///   word, and carries every shape whose `rows + width` fits the word — two
///   `u64`s up to 64, two 128-bit words up to 128. It has no length header
///   and never spills to the heap, and its fixed size lets the hash table
///   probe it in one or two word compares.
/// * [`Key`] serializes the elements at a byte or more each and is the
///   fallback past 128, where a bitmap would be longer than the bytes.
///
/// [`expand_oriented`] decides once per expansion, from the shape alone
/// ([`packed_bits`]); every worker of that expansion then uses the same
/// representation, which is what makes equal states equal keys.
trait LayerKey: Sized + Eq + Hash + Send + Sync + 'static {
    /// Per-worker scratch the serialization needs; `()` when it needs none.
    type Scratch: Default + Send;

    /// Whether this representation can hold every state of the shape.
    fn fits(outer: &Partition, inner: &Partition) -> bool;

    /// The state before any row is filled: empty content, no row above.
    fn root(elem_width: usize) -> Self;

    /// Split the state back into content and the row above.
    ///
    /// `elem_width` is the byte width [`Key`] serializes at (see
    /// [`elem_width`]); `above_len` is the number of cells in the row above,
    /// which the geometry knows and the packed form checks itself against.
    fn decode(
        &self,
        elem_width: usize,
        above_len: usize,
        content: &mut Vec<u32>,
        above: &mut Vec<u32>,
    );

    /// Build the successor of the state in `ctx` from the row just completed
    /// — content folded in, row clipped to the next overlap — and merge it
    /// into its shard of `ctx.out`.
    fn commit<C: Acc>(ctx: &mut RowCtx<'_, C, Self>);
}

/// A machine word wide enough for one lattice-path bitmap.
///
/// `u64` for shapes with `rows + width ≤ 64`; [`W128`] up to 128. The trait
/// is the handful of bit operations the encoding needs, so the encode and
/// decode loops are written once.
trait Word: Copy + Eq + Hash + Send + Sync + 'static {
    /// The number of bit positions the word holds.
    const BITS: usize;
    const ZERO: Self;
    fn set(self, pos: u32) -> Self;
    fn is_zero(self) -> bool;
    fn count_ones(self) -> u32;
    fn trailing_zeros(self) -> u32;
    /// The word with its lowest set bit cleared.
    fn clear_lowest(self) -> Self;
    /// A 64-bit digest for shard routing; any fixed mixing will do.
    fn fold(self) -> u64;
}

impl Word for u64 {
    const BITS: usize = 64;
    const ZERO: Self = 0;
    #[inline]
    fn set(self, pos: u32) -> Self {
        self | (1u64 << pos)
    }
    #[inline]
    fn is_zero(self) -> bool {
        self == 0
    }
    #[inline]
    fn count_ones(self) -> u32 {
        u64::count_ones(self)
    }
    #[inline]
    fn trailing_zeros(self) -> u32 {
        u64::trailing_zeros(self)
    }
    #[inline]
    fn clear_lowest(self) -> Self {
        self & self.wrapping_sub(1)
    }
    #[inline]
    fn fold(self) -> u64 {
        self
    }
}

/// 128 bits as two `u64`s, low word first.
///
/// Not `u128`: that type is 16-byte aligned, which would pad a
/// `(PackedKey<u128>, u64)` layer entry from 40 to 48 bytes, and entry size
/// is what the layer's speed follows. Two `u64`s keep it at 40.
#[derive(Clone, Copy, PartialEq, Eq)]
struct W128([u64; 2]);

impl Hash for W128 {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u64(self.0[0]);
        state.write_u64(self.0[1]);
    }
}

impl Word for W128 {
    const BITS: usize = 128;
    const ZERO: Self = W128([0, 0]);
    #[inline]
    fn set(self, pos: u32) -> Self {
        let W128([lo, hi]) = self;
        if pos < 64 {
            W128([lo | (1u64 << pos), hi])
        } else {
            W128([lo, hi | (1u64 << (pos - 64))])
        }
    }
    #[inline]
    fn is_zero(self) -> bool {
        self.0 == [0, 0]
    }
    #[inline]
    fn count_ones(self) -> u32 {
        self.0[0].count_ones() + self.0[1].count_ones()
    }
    #[inline]
    fn trailing_zeros(self) -> u32 {
        if self.0[0] != 0 {
            self.0[0].trailing_zeros()
        } else {
            64 + self.0[1].trailing_zeros()
        }
    }
    #[inline]
    fn clear_lowest(self) -> Self {
        let W128([lo, hi]) = self;
        if lo != 0 {
            W128([lo & lo.wrapping_sub(1), hi])
        } else {
            W128([lo, hi & hi.wrapping_sub(1)])
        }
    }
    #[inline]
    fn fold(self) -> u64 {
        self.0[0] ^ self.0[1].rotate_left(32)
    }
}

/// A layer state as two lattice-path bitmaps.
///
/// Both halves of the state are sequences of bounded, monotone integers, and
/// such a sequence is a set of distinct bit positions:
///
/// * `content` is a partition — weakly decreasing, and with `k` nonzero parts
///   part `i` (0-indexed) is stored as bit `content[i] + (k − 1 − i)`. The
///   positions strictly decrease with `i`, so the word has exactly `k` bits
///   set, and its highest one sits at `content[0] + k − 1`.
/// * `above` is weakly increasing, and cell `j` of `d` is stored as bit
///   `above[j] + j`, strictly increasing in `j`.
///
/// Both maps are injective (the count of set bits recovers `k`; `d` is the
/// row geometry), so equal states are equal words. The bound that makes a
/// word suffice: `content[0]` counts 1's, which form a horizontal strip and so
/// number at most `width`, while `k` grows by at most one per row and
/// `above`'s values are at most one more than `k`. So every position is
/// below `rows + width` — [`packed_bits`], the quantity
/// [`fits`](LayerKey::fits) tests. Shifts are checked in every profile, so a
/// state outside the bound is a panic, not a wrong coefficient.
///
/// The bitmap is shorter than a byte-per-element serialization wherever the
/// layer is large: the conjugate walk of `[20,16,12,8]²` carries content of
/// up to 40 parts in a 48-bit word, and its whole state in 16 bytes.
#[derive(Clone, Copy, PartialEq, Eq, Hash)]
struct PackedKey<W: Word> {
    content: W,
    above: W,
}

/// `rows + width` of the skew shape: one more than the largest bit position a
/// [`PackedKey`] of the shape can need.
fn packed_bits(outer: &Partition, _inner: &Partition) -> usize {
    outer.len() + outer.part(0) as usize
}

/// Decode a partition bitmap (see [`PackedKey`]) into its parts, largest first.
#[inline]
fn unpack_partition<W: Word>(mut w: W, out: &mut Vec<u32>) {
    let k = w.count_ones() as usize;
    out.clear();
    out.resize(k, 0);
    // The j-th lowest set bit is part k − 1 − j, stored at part + j.
    let mut j = 0u32;
    while !w.is_zero() {
        out[k - 1 - j as usize] = w.trailing_zeros() - j;
        w = w.clear_lowest();
        j += 1;
    }
}

/// Decode a row bitmap (see [`PackedKey`]) into its cells, left to right.
#[inline]
fn unpack_row<W: Word>(mut w: W, out: &mut Vec<u32>) {
    out.clear();
    let mut j = 0u32;
    while !w.is_zero() {
        out.push(w.trailing_zeros() - j);
        w = w.clear_lowest();
        j += 1;
    }
}

impl<W: Word> LayerKey for PackedKey<W> {
    type Scratch = ();

    fn fits(outer: &Partition, inner: &Partition) -> bool {
        packed_bits(outer, inner) <= W::BITS
    }

    fn root(_elem_width: usize) -> Self {
        PackedKey {
            content: W::ZERO,
            above: W::ZERO,
        }
    }

    #[inline]
    fn decode(
        &self,
        _elem_width: usize,
        above_len: usize,
        content: &mut Vec<u32>,
        above: &mut Vec<u32>,
    ) {
        unpack_partition(self.content, content);
        unpack_row(self.above, above);
        debug_assert_eq!(above.len(), above_len);
    }

    #[inline]
    fn commit<C: Acc>(ctx: &mut RowCtx<'_, C, Self>) {
        let clen = ctx.content.len();
        // `added` has exactly one slot past the old content (the row-wide
        // value cap); a new value appears there or nowhere.
        let new = ctx.added[clen];
        let k = clen + usize::from(new > 0);
        let mut content = W::ZERO;
        for (i, (&c, &a)) in ctx.content.iter().zip(ctx.added.iter()).enumerate() {
            content = content.set(c + a + (k - 1 - i) as u32);
        }
        if new > 0 {
            // Part `clen` of `k = clen + 1`: shift zero.
            content = content.set(new);
        }
        let mut above = W::ZERO;
        for j in 0..ctx.dn_hi - ctx.dn_lo {
            above = above.set(ctx.row[ctx.dn_lo - ctx.lo + j] + j as u32);
        }
        let key = PackedKey { content, above };

        // Any fixed mixing of the two words routes every copy of a state to the
        // same shard whichever worker produced it, which is all the parallel
        // merge needs. R4's modular-by-definition case, as in `shard_of`.
        let mixed =
            (content.fold() ^ above.fold().rotate_left(29)).wrapping_mul(0x9e37_79b9_7f4a_7c15);
        let shard = &mut ctx.out[((mixed >> 32) as usize) % ctx.out.len()];
        // The key is a word or two, so it is built before the probe and the
        // entry API hashes it once for merge and insert alike.
        match shard.entry(key) {
            Entry::Occupied(mut e) => match e.get().checked_add(ctx.mult) {
                Some(sum) => *e.get_mut() = sum,
                None => *ctx.overflow = true,
            },
            Entry::Vacant(e) => {
                e.insert(ctx.mult);
            }
        }
    }
}

/// The byte serialization of a layer state: the fallback [`LayerKey`] for
/// shapes too wide for a [`PackedKey`].
///
/// The state is the sequence `[len, content.., above..]`: `content[v - 1]` is
/// the number of `v`s placed so far (no trailing zeros, so the leading `len`
/// makes the split unambiguous), and `above` holds the previous row's entries,
/// clipped to the columns the next row overlaps.
///
/// * Every element — the length header, each content count, each cell value —
///   is at most the cell count of the shape (see [`elem_width`]), so the
///   sequence is serialized at the narrowest sufficient byte width, fixed per
///   expansion.
/// * Keys of ≤ [`INLINE`] bytes are stored inline in the enum, so such a *new*
///   state costs no heap allocation. Longer keys spill to a box.
///
/// [`INLINE`] is where two costs meet, and the meeting point is sharp. Making
/// keys wide enough that nothing spills removes every one of those
/// allocations and runs slower, because the layer is the working set and more
/// bytes per entry outweigh a malloc and a free per state; making them
/// narrower spills everything and is worse still
/// (`docs/record/littlewood-richardson.md`). That measurement is what the
/// bitmap form answers: it needs neither the header nor the box.
///
/// Merges (the common case) are probed with a borrowed scratch buffer via
/// [`KeyBytes`], so the hot path allocates nothing either way.
enum Key {
    /// `(length, bytes)`; only `bytes[..length]` is meaningful.
    Inline(u8, [u8; INLINE]),
    Heap(Box<[u8]>),
}

/// Inline capacity, chosen so `size_of::<Key>()` is 32: tag + 1 + 30 on one
/// side, a 16-byte box on the other. That makes a layer entry 40 bytes, the
/// tuned point between spilling and widening
/// (`docs/record/littlewood-richardson.md`).
const INLINE: usize = 30;

const _: () = assert!(std::mem::size_of::<Key>() == 32);
const _: () = assert!(std::mem::size_of::<PackedKey<u64>>() == 16);
const _: () = assert!(std::mem::size_of::<PackedKey<W128>>() == 32);
const _: () = assert!(std::mem::align_of::<PackedKey<W128>>() == 8);

/// Per-worker scratch for assembling a byte key: the `u32` sequence and its
/// serialization.
#[derive(Default)]
struct ByteScratch {
    key: Vec<u32>,
    kbuf: Vec<u8>,
}

impl LayerKey for Key {
    type Scratch = ByteScratch;

    fn fits(_outer: &Partition, _inner: &Partition) -> bool {
        true
    }

    fn root(elem_width: usize) -> Self {
        Key::from_bytes(&vec![0u8; elem_width])
    }

    #[inline]
    fn decode(
        &self,
        elem_width: usize,
        above_len: usize,
        content: &mut Vec<u32>,
        above: &mut Vec<u32>,
    ) {
        let bytes = self.bytes();
        // The header is one element; decode it alone to find the split.
        decode_into(content, &bytes[..elem_width], elem_width);
        let split = elem_width * (1 + content[0] as usize);
        decode_into(content, &bytes[elem_width..split], elem_width);
        decode_into(above, &bytes[split..], elem_width);
        debug_assert_eq!(above.len(), above_len);
    }

    #[inline]
    fn commit<C: Acc>(ctx: &mut RowCtx<'_, C, Self>) {
        // Assemble the successor key in the scratch buffer: [len, content, above].
        let ByteScratch { key, kbuf } = &mut *ctx.scratch;
        key.clear();
        key.push(0); // length, patched below
        key.extend(
            ctx.content
                .iter()
                .zip(ctx.added.iter())
                .map(|(&c, &a)| c + a),
        );
        // `added` has exactly one slot past the old content (the row-wide value
        // cap); a new value appears there or nowhere, so no zero-trimming is
        // ever needed.
        if let Some(&n) = ctx.added.get(ctx.content.len()) {
            if n > 0 {
                key.push(n);
            }
        }
        key[0] = (key.len() - 1) as u32;
        if ctx.dn_lo < ctx.dn_hi {
            key.extend_from_slice(&ctx.row[ctx.dn_lo - ctx.lo..ctx.dn_hi - ctx.lo]);
        }
        encode_into(kbuf, key, ctx.width);

        // Probe with the borrowed scratch; a merge (the common case) allocates
        // nothing, and a genuinely new state is copied inline unless oversized.
        let shard = &mut ctx.out[shard_of(kbuf, ctx.out.len())];
        match shard.get_mut(KeyBytes::new(kbuf)) {
            Some(w) => match w.checked_add(ctx.mult) {
                Some(sum) => *w = sum,
                None => *ctx.overflow = true,
            },
            None => {
                shard.insert(Key::from_bytes(kbuf), ctx.mult);
            }
        }
    }
}

impl Key {
    /// Canonical constructor: inline iff it fits, so equal byte strings always
    /// get the same representation (though `Eq`/`Hash` don't rely on that).
    #[inline]
    fn from_bytes(b: &[u8]) -> Key {
        if b.len() <= INLINE {
            let mut buf = [0u8; INLINE];
            buf[..b.len()].copy_from_slice(b);
            Key::Inline(b.len() as u8, buf)
        } else {
            Key::Heap(b.into())
        }
    }

    #[inline]
    fn bytes(&self) -> &[u8] {
        match self {
            Key::Inline(len, buf) => &buf[..*len as usize],
            Key::Heap(b) => b,
        }
    }
}

impl PartialEq for Key {
    #[inline]
    fn eq(&self, other: &Key) -> bool {
        self.bytes() == other.bytes()
    }
}

impl Eq for Key {}

impl Hash for Key {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write(self.bytes());
    }
}

/// Borrowed view of a [`Key`], so a layer probe can use a scratch slice.
///
/// `#[repr(transparent)]` makes the `&[u8]` → `&KeyBytes` cast sound; `Hash`
/// and `Eq` agree with [`Key`]'s exactly, which is what `Borrow` requires.
#[derive(PartialEq, Eq)]
#[repr(transparent)]
struct KeyBytes([u8]);

impl KeyBytes {
    #[inline]
    fn new(s: &[u8]) -> &KeyBytes {
        // SAFETY: `KeyBytes` is `repr(transparent)` over `[u8]`, so the two
        // have the same layout, alignment and pointer metadata, and the
        // returned reference borrows `s` for the same lifetime.
        unsafe { &*(s as *const [u8] as *const KeyBytes) }
    }
}

impl Borrow<KeyBytes> for Key {
    #[inline]
    fn borrow(&self) -> &KeyBytes {
        KeyBytes::new(self.bytes())
    }
}

impl Hash for KeyBytes {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write(&self.0);
    }
}

/// High-water mark of live layer states, for measurement harnesses.
///
/// Peak memory is (states) × (bytes per state); this records the first factor,
/// which — unlike RSS — is not perturbed by the allocator. Sampled once per
/// row at the point both the old and new layer are fully populated, so it
/// is the true in-traversal maximum of live map entries. Monotone across
/// expansions until read: [`take_peak_layer_states`] returns and resets it.
static PEAK_LIVE_STATES: AtomicUsize = AtomicUsize::new(0);

/// Returns the high-water mark of live layer states and resets it to zero.
///
/// The count is the number of live map entries, monotone across expansions
/// until read. A measurement hook, not part of the semantic API.
pub fn take_peak_layer_states() -> usize {
    PEAK_LIVE_STATES.swap(0, Ordering::Relaxed)
}

/// Running total of row fillings committed to a layer, for measurement
/// harnesses.
///
/// Every completed row of every state is one production, merged or not, so
/// this counts the work the traversal does where [`PEAK_LIVE_STATES`] counts
/// what it holds. Their ratio to the number of LR tableaux is the layer's
/// compression, the quantity that decides between orientations
/// (`docs/record/littlewood-richardson.md`). Accumulated once per chunk from
/// a worker-local count, so the hot path pays one register increment.
static PRODUCTIONS: AtomicU64 = AtomicU64::new(0);

/// Returns the number of row fillings committed since the last call and
/// resets it to zero. A measurement hook, not part of the semantic API.
pub fn take_productions() -> u64 {
    PRODUCTIONS.swap(0, Ordering::Relaxed)
}

/// The byte width every element of a layer key fits in, for this shape.
///
/// The bound is the cell count `n = |outer| − |inner|`: content counts total
/// exactly the cells filled so far; any placed value `v` has all of `1..v`
/// present in the content (ballot), so cell values — hence also the content
/// length and its header — never exceed `n`. Row count does not enter.
///
/// Narrowing on a proven bound is exactly where an off-by-one becomes a wrong
/// coefficient, so [`encode_into`] debug-asserts every element against the
/// chosen width, and the width-boundary tests below cross n = 255 and 65535.
fn elem_width(outer: &Partition, inner: &Partition) -> usize {
    let n = outer.size() - inner.size();
    if n <= 0xFF {
        1
    } else if n <= 0xFFFF {
        2
    } else {
        4
    }
}

/// Serialize `src` little-endian at `w` bytes per element.
#[inline]
fn encode_into(dst: &mut Vec<u8>, src: &[u32], w: usize) {
    dst.clear();
    match w {
        1 => {
            debug_assert!(src.iter().all(|&x| x <= 0xFF));
            dst.extend(src.iter().map(|&x| x as u8));
        }
        2 => {
            debug_assert!(src.iter().all(|&x| x <= 0xFFFF));
            for &x in src {
                dst.extend_from_slice(&(x as u16).to_le_bytes());
            }
        }
        _ => {
            for &x in src {
                dst.extend_from_slice(&x.to_le_bytes());
            }
        }
    }
}

/// Deserialize a key back to `u32` elements; inverse of [`encode_into`].
#[inline]
fn decode_into(dst: &mut Vec<u32>, src: &[u8], w: usize) {
    dst.clear();
    match w {
        1 => dst.extend(src.iter().map(|&b| b as u32)),
        2 => dst.extend(
            src.chunks_exact(2)
                .map(|c| u16::from_le_bytes([c[0], c[1]]) as u32),
        ),
        _ => dst.extend(
            src.chunks_exact(4)
                .map(|c| u32::from_le_bytes([c[0], c[1], c[2], c[3]])),
        ),
    }
}

/// A layer multiplicity: `u64` for the fast pass, `u128` for the fallback.
///
/// Multiplicities are *tableau counts*, which dwarf the final coefficients
/// (2.1 × 10⁸ tableaux behind 10⁵-ish coefficients on `s[8,7,6,5,4,3]²`), so
/// no fixed narrow type is provably safe. Instead every add is checked: the
/// `u64` pass detects saturation and [`expand_skew_uncached`] transparently
/// reruns in `u128` — correctness never rests on an unproven bound, and the
/// 8-bytes-per-state saving is kept on every shape that stays under 2⁶⁴.
trait Acc: Copy + Send + Sync {
    const ONE: Self;
    fn checked_add(self, other: Self) -> Option<Self>;
    fn widen(self) -> u128;
}

impl Acc for u64 {
    const ONE: Self = 1;
    #[inline]
    fn checked_add(self, other: Self) -> Option<Self> {
        u64::checked_add(self, other)
    }
    #[inline]
    fn widen(self) -> u128 {
        self as u128
    }
}

impl Acc for u128 {
    const ONE: Self = 1;
    #[inline]
    fn checked_add(self, other: Self) -> Option<Self> {
        u128::checked_add(self, other)
    }
    #[inline]
    fn widen(self) -> u128 {
        self
    }
}

fn expand_skew_uncached(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
    expand_walk(outer, inner, prefer_conjugate(outer, inner))
}

/// The expansion of `outer/inner`, walked on the diagram itself or, with
/// `conjugate`, on its transpose with the terms conjugated back; sorted.
fn expand_walk(outer: &Partition, inner: &Partition, conjugate: bool) -> Vec<(Partition, u128)> {
    let (co, ci);
    let (o, i) = if conjugate {
        co = outer.conjugate();
        ci = inner.conjugate();
        (&co, &ci)
    } else {
        (outer, inner)
    };
    let mut out = expand_oriented::<u64>(o, i, conjugate)
        .or_else(|| expand_oriented::<u128>(o, i, conjugate))
        .expect("LR tableau multiplicity exceeded u128");
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

/// Whether to walk the transposed diagram instead, for a skew expansion.
///
/// Products do not come here: [`product_walk`] chooses their orientation
/// from the two factors, which say more than the juxtaposed shape does.
///
/// c^λ_{μν} = c^{λ'}_{μ'ν'}, so the expansion may run on either orientation
/// and conjugate its terms back. Peak layer size is close to
/// orientation-independent — within ±25% on every case measured, since the
/// states carry the same information either way. *Time* is not: the per-row
/// run fill enumerates fillings whose count grows combinatorially with row
/// width and with the number of distinct values. So wide-and-deep diagrams
/// walk far faster on their side, and the gap widens with the shape
/// (`docs/record/littlewood-richardson.md`).
///
/// The thresholds are empirical. Diagrams with fewer rows than this stay
/// direct, as do diagrams no wider than tall (a staircase's conjugate is
/// itself; tall shapes already walk well) and small shapes (the conjugate's
/// extra rows cost more than they save — the known losses to this rule are
/// rectangles like `[12⁶]²`). Everything the rule fires on was measured at a
/// clear win or a tie.
///
/// The rectangle losses do not reach here through the default backend:
/// [`AutoLr`](crate::strip_lr::AutoLr) routes a rectangle-times-rectangle
/// product to [`crate::rect`], which has a closed form. They still matter for
/// callers that name `SkewLr` directly, and for rectangular *skew* shapes.
///
/// `SKEW_ORIENT=direct` or `SKEW_ORIENT=conj` in the environment overrides
/// the rule, so a calibration harness can time both walks of one shape in
/// one binary. A measurement hook, read once per uncached expansion.
fn prefer_conjugate(outer: &Partition, inner: &Partition) -> bool {
    match std::env::var_os("SKEW_ORIENT") {
        Some(v) if v == "conj" => return true,
        Some(v) if v == "direct" => return false,
        _ => {}
    }
    let rows = outer.len();
    let width = outer.part(0) as usize;
    let cells = outer.size() - inner.size();
    // Only where the transpose compresses. The rows of `outer` below the
    // last row of `inner` are full rows; the direct walk meets them last and
    // the transposed walk meets them first, and with more than four of them
    // the transposed layer commits as many fillings as the direct one or
    // more, at 1.3–2x the cost each for its longer content. With four or
    // fewer — a band, the inner shape reaching nearly to the bottom — the
    // transposed walk merges where the direct one cannot, up to 6x. An outer
    // shape descending by at most one cell per row transposes at any depth;
    // that is the family the rule was first calibrated on, and it holds
    // (`docs/record/littlewood-richardson.md`, "mixed-step outer shapes").
    let gentle = outer.parts().windows(2).all(|w| w[0] - w[1] <= 1);
    let band = outer.len().saturating_sub(inner.len()) <= 4;
    rows >= 8 && width > rows && cells >= 60 && (gentle || band)
}

/// One layered traversal with multiplicities in `C`; `None` means some merge
/// overflowed `C` and the caller should retry wider.
///
/// Picks the layer key from the shape: the narrowest [`PackedKey`] whose word
/// holds [`packed_bits`], else the byte [`Key`].
///
/// `conjugate_terms` reports each term as the conjugate of the content the
/// walk found — the form [`prefer_conjugate`] needs — one term at a time. So
/// no intermediate vector of unconjugated partitions is ever materialized:
/// each would own a heap block of up to `rows` parts beside the output's own,
/// and on a multi-million-term expansion that transient is the size of the
/// answer again (`docs/record/memory.md`, "The unconjugated-term
/// transient").
fn expand_oriented<C: Acc>(
    outer: &Partition,
    inner: &Partition,
    conjugate_terms: bool,
) -> Option<Vec<(Partition, u128)>> {
    let width = elem_width(outer, inner);
    if PackedKey::<u64>::fits(outer, inner) {
        expand_layer::<C, PackedKey<u64>>(outer, inner, width, conjugate_terms)
    } else if PackedKey::<W128>::fits(outer, inner) {
        expand_layer::<C, PackedKey<W128>>(outer, inner, width, conjugate_terms)
    } else {
        expand_layer::<C, Key>(outer, inner, width, conjugate_terms)
    }
}

/// Whether `SKEW_TRACE` is set, read from the environment once per process.
///
/// Every Schur product passes through [`expand_layer`], and
/// `std::env::var_os` takes the process environment lock and scans it on
/// each call: 60 ns, against about 2 µs for the smallest cold product
/// (`docs/record/littlewood-richardson.md`, "Reading `SKEW_TRACE` once"). The
/// variable is set before a traced run starts, never during one, so a
/// per-process read is the same facility.
fn skew_trace() -> bool {
    static TRACE: OnceLock<bool> = OnceLock::new();
    *TRACE.get_or_init(|| std::env::var_os("SKEW_TRACE").is_some())
}

/// [`expand_oriented`] with the layer key type and the byte element width
/// chosen by the caller. Split out so tests can force each representation
/// onto the same shape and pin their agreement; `K` must fit the shape and
/// `width` must be sufficient for it (any width is, when ≥ the chosen one).
fn expand_layer<C: Acc, K: LayerKey>(
    outer: &Partition,
    inner: &Partition,
    width: usize,
    conjugate_terms: bool,
) -> Option<Vec<(Partition, u128)>> {
    debug_assert!(
        K::fits(outer, inner),
        "layer key too narrow for {outer}/{inner}"
    );
    let rows = outer.len();
    let trace = skew_trace();

    // The layer of the traversal: reduced state -> number of ways to reach
    // it. Between rows it is held as a plain `Vec`: the hash table is only
    // needed on the side being merged *into*, and a table's bucket array
    // (power-of-two, reserved ahead) can run 2–4× the entry payload. Draining
    // each finished table into an exactly-sized vector caps the steady-state
    // layer at real entries only, and reading it back is a linear scan
    // instead of a table walk.
    let mut cur: Vec<(K, C)> = vec![(K::root(width), C::ONE)];

    // Scratch now lives in `fill_chunk`, which is per worker; these two are
    // for decoding the finished layer below.
    let mut content: Vec<u32> = Vec::new();
    let mut above: Vec<u32> = Vec::new();
    let mut overflow = false;

    // On the row loop and not inside `fill_chunk`: the fill runs on scoped
    // workers, whose join treats any panic as a bug, so a cancellation raised
    // there would arrive as that bug rather than as itself
    // (`crate::interrupt`, "Which thread may poll").
    for r in 0..rows {
        crate::interrupt::poll();
        let lo = inner.part(r) as usize;
        let hi = outer.part(r) as usize;
        // Columns of row r-1 that row r sits under, and the ones of row r that
        // row r+1 will sit under. Clipping to these is what makes distinct
        // histories collapse.
        // `r == 0` wraps to `usize::MAX`, which is the sentinel for "there is no
        // row above row 0": `overlap` rejects any index past `outer.len()` and
        // returns the empty span. Wrapping is the encoding, not an accident —
        // an unsigned row index has no −1 to hold (R4).
        let (up_lo, up_hi) = overlap(inner, outer, r.wrapping_sub(1), r);
        let (dn_lo, dn_hi) = overlap(inner, outer, r, r + 1);

        let geom = RowGeom {
            lo,
            hi,
            up_lo,
            up_hi,
            dn_lo,
            dn_hi,
            width,
        };
        let next: Vec<(K, C)> = fill_row(&cur, &geom, &mut overflow);
        // Both layers are momentarily live here; record the sum (once per
        // row, so the cost is nil).
        PEAK_LIVE_STATES.fetch_max(cur.len() + next.len(), Ordering::Relaxed);
        if trace {
            eprintln!("skew_lr row {r}: {} -> {} states", cur.len(), next.len());
        }
        if overflow {
            return None;
        }
        cur = next;
        if cur.is_empty() {
            break;
        }
    }

    // Every row's down-overlap with a nonexistent next row is empty, so the
    // final keys are pure content — already merged, one term each.
    Some(
        cur.into_iter()
            .map(|(key, c)| {
                key.decode(width, 0, &mut content, &mut above);
                debug_assert!(above.is_empty());
                // The ballot condition forces the content to be weakly
                // decreasing at every prefix, so the final one is already a
                // partition.
                let parts = if conjugate_terms {
                    crate::partition::conjugate_parts(&content)
                } else {
                    content.clone()
                };
                (Partition::from_sorted(parts), c.widen())
            })
            .collect(),
    )
}

/// The columns shared by the skew parts of rows `a` and `b` (`a` may underflow
/// to mean "no row", giving an empty range).
fn overlap(inner: &Partition, outer: &Partition, a: usize, b: usize) -> (usize, usize) {
    if a >= outer.len() || b >= outer.len() {
        return (0, 0);
    }
    let lo = (inner.part(a) as usize).max(inner.part(b) as usize);
    let hi = (outer.part(a) as usize).min(outer.part(b) as usize);
    if lo >= hi {
        (0, 0)
    } else {
        (lo, hi)
    }
}

/// The column geometry of one row: fixed for the whole row, so it is computed
/// once and shared by every state (and every thread) processing that row.
#[derive(Clone, Copy)]
struct RowGeom {
    lo: usize,
    hi: usize,
    up_lo: usize,
    up_hi: usize,
    dn_lo: usize,
    dn_hi: usize,
    width: usize,
}

/// Below this many states a row is filled on one thread. Spawning costs tens of
/// microseconds and the merge is not free, so on small layers the parallel
/// path is pure loss — and most rows of most shapes are small even when the
/// peak is not.
const PARALLEL_MIN_STATES: usize = 24_576;

/// Which shard a byte [`Key`] belongs to; a [`PackedKey`] mixes its two words
/// the same way when it commits.
///
/// Cheap on purpose: hashbrown will hash the key again on insert, so anything
/// thorough here is paid twice. Keys are `[len, content.., row..]`, and the row
/// suffix is what actually varies between states, so mixing the tail bytes
/// spreads them well. Balance only affects merge parallelism, never
/// correctness — every copy of a key lands in the same shard whichever worker
/// produced it, which is the property the parallel merge needs.
#[inline]
fn shard_of(bytes: &[u8], shards: usize) -> usize {
    let n = bytes.len();
    let mut w = n as u64;
    for &b in bytes.iter().rev().take(8) {
        w = (w << 8) | b as u64;
    }
    // R4's modular-by-definition case: this is a hash mixer, where the product
    // mod 2^64 *is* the intended operation rather than a truncated one. The
    // constant is the 64-bit golden-ratio odd multiplier.
    w = w.wrapping_mul(0x9e37_79b9_7f4a_7c15);
    ((w >> 32) as usize) % shards
}

/// Fill one row: every state in `cur` expanded into a fresh layer.
///
/// The row is a barrier — row r+1 cannot start until row r is complete — so
/// this is bulk-synchronous, and the only question is how to split the states
/// within a row. Each worker owns private layers and they are combined at
/// the end, rather than sharing one behind a lock: the layer is written on
/// *every* emitted filling, so a shared table would serialize the hot path
/// exactly where the work is.
///
/// **The layer is sharded, so the combine is parallel too.** Merging every
/// worker's table into one was a third to a half of wall time on the large
/// shapes, which is an Amdahl ceiling no core count can lift
/// (`docs/record/littlewood-richardson.md`). Routing each key to a shard by a
/// cheap hash puts every copy of a key in the same shard whoever produced it,
/// so shard `j` can be combined from all workers independently of shard `k`.
fn fill_row<C: Acc, K: LayerKey>(
    cur: &[(K, C)],
    geom: &RowGeom,
    overflow: &mut bool,
) -> Vec<(K, C)> {
    // Many small chunks claimed from a shared counter, rather than one slice per
    // worker. This machine — like most now — is heterogeneous: 4 performance
    // cores and 6 efficiency cores, the latter roughly a third the throughput.
    // With an even split the row barrier waits on whichever chunk landed on the
    // slowest core, so a static partition gives up much of the parallelism
    // before any of it is used. Claiming work on demand lets a fast core take
    // three chunks while a slow one takes one.
    const CHUNK: usize = 2_048;
    let threads = worker_count(cur.len());
    if threads <= 1 {
        let mut out = vec![Map::with_capacity_and_hasher(cur.len(), Default::default())];
        // Chunked only so there is somewhere to poll: a row of a large shape is
        // the longest stretch this path runs without returning, and states are
        // independent, so splitting the slice changes nothing but that.
        for lo in (0..cur.len()).step_by(CHUNK) {
            crate::interrupt::poll();
            fill_chunk(
                &cur[lo..(lo + CHUNK).min(cur.len())],
                geom,
                &mut out,
                overflow,
            );
        }
        // One shard, so flattening is the shard — and unlike `pop` it needs no
        // claim about how many there are.
        return out.into_iter().flatten().collect();
    }
    let shards = threads;
    let nchunks = cur.len().div_ceil(CHUNK);
    let cursor = AtomicUsize::new(0);
    let parts: Vec<(Vec<Map<K, C>>, bool)> = std::thread::scope(|s| {
        let handles: Vec<_> = (0..threads)
            .map(|_| {
                let cursor = &cursor;
                s.spawn(move || {
                    let cap = cur.len() / (threads * shards) + 1;
                    let mut out: Vec<Map<K, C>> = (0..shards)
                        .map(|_| Map::with_capacity_and_hasher(cap, Default::default()))
                        .collect();
                    let mut of = false;
                    loop {
                        let i = cursor.fetch_add(1, Ordering::Relaxed);
                        if i >= nchunks {
                            break;
                        }
                        let lo = i * CHUNK;
                        let hi = (lo + CHUNK).min(cur.len());
                        fill_chunk(&cur[lo..hi], geom, &mut out, &mut of);
                    }
                    (out, of)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .expect("an LR fill worker panicked; its panic is the bug")
            })
            .collect()
    });

    // True peak: before the merge the same key can exist once per worker, so
    // the live entry count here exceeds the merged layer. Sampling only
    // after the merge (as the row loop does) cannot see that, and would report
    // the parallel path as free when it is not.
    let pre: usize = parts
        .iter()
        .map(|(ms, _)| ms.iter().map(|m| m.len()).sum::<usize>())
        .sum();
    PEAK_LIVE_STATES.fetch_max(cur.len() + pre, Ordering::Relaxed);

    // Transpose: shard j gathers its table from every worker.
    let mut columns: Vec<Vec<Map<K, C>>> =
        (0..shards).map(|_| Vec::with_capacity(threads)).collect();
    for (maps, of) in parts {
        if of {
            *overflow = true;
        }
        for (j, m) in maps.into_iter().enumerate() {
            columns[j].push(m);
        }
    }
    // …and each shard is combined independently, in parallel.
    let merged: Vec<(Vec<(K, C)>, bool)> = std::thread::scope(|s| {
        let handles: Vec<_> = columns
            .into_iter()
            .map(|mut group| {
                s.spawn(move || {
                    let mut of = false;
                    // Into the largest, so the biggest table is never reinserted.
                    let best = group
                        .iter()
                        .enumerate()
                        .max_by_key(|(_, m)| m.len())
                        .map(|(i, _)| i)
                        .unwrap_or(0);
                    let mut acc = group.swap_remove(best);
                    for m in group {
                        for (k, v) in m {
                            // The key is owned here, so one hash serves both
                            // outcomes; a probe followed by an insert would
                            // hash every key new to `acc` twice.
                            match acc.entry(k) {
                                Entry::Occupied(mut e) => match e.get().checked_add(v) {
                                    Some(sum) => *e.get_mut() = sum,
                                    None => of = true,
                                },
                                Entry::Vacant(e) => {
                                    e.insert(v);
                                }
                            }
                        }
                    }
                    (acc.into_iter().collect::<Vec<_>>(), of)
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| {
                h.join()
                    .expect("an LR fill worker panicked; its panic is the bug")
            })
            .collect()
    });

    let mut out = Vec::with_capacity(merged.iter().map(|(v, _)| v.len()).sum());
    for (v, of) in merged {
        if of {
            *overflow = true;
        }
        out.extend(v);
    }
    out
}

/// How many workers to use for a layer of `states`.
///
/// Returns 1 whenever the row is too small to pay for the split.
fn worker_count(states: usize) -> usize {
    if states < PARALLEL_MIN_STATES {
        return 1;
    }
    let avail = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    // Never more workers than there is work to give them.
    avail.min(states / (PARALLEL_MIN_STATES / 2)).max(1)
}

/// Expand every state in `chunk` into `out`. All scratch is local, so this is
/// what a worker runs.
fn fill_chunk<C: Acc, K: LayerKey>(
    chunk: &[(K, C)],
    geom: &RowGeom,
    out: &mut [Map<K, C>],
    overflow: &mut bool,
) {
    let RowGeom {
        lo,
        hi,
        up_lo,
        up_hi,
        dn_lo,
        dn_hi,
        width,
    } = *geom;
    let mut row: Vec<u32> = vec![0; hi.saturating_sub(lo)];
    let mut added: Vec<u32> = Vec::new();
    let mut cut: Vec<usize> = Vec::new();
    let mut gap: Vec<u32> = Vec::new();
    let mut scratch = K::Scratch::default();
    let mut content: Vec<u32> = Vec::new();
    let mut above: Vec<u32> = Vec::new();
    let mut produced: u64 = 0;

    for (state, mult) in chunk.iter() {
        state.decode(width, up_hi - up_lo, &mut content, &mut above);
        let clen = content.len();
        let top = clen as u32 + 1;

        // cut[v]: the first column whose cell above holds a value ≥ v, or
        // `hi` when no column does. `above` is weakly increasing, so the
        // columns closed to a value v are exactly [cut[v], up_hi) — one
        // merge scan answers every column-strictness question for this
        // state's row.
        cut.clear();
        cut.push(0); // v = 0, unused
        let mut j = 0;
        for v in 1..=top {
            while j < above.len() && above[j] < v {
                j += 1;
            }
            cut.push(if j < above.len() { up_lo + j } else { hi });
        }

        // gap[v]: how many v's this row may add before the ballot condition
        // #(v-1) ≥ #v breaks. Only content *before* the row counts, because
        // the reading word takes the row right to left (see `fill_runs`).
        gap.clear();
        gap.extend_from_slice(&[0, 0]); // v = 0, 1: never constrained
        for v in 2..=top {
            let vi = v as usize - 1;
            gap.push(content[vi - 1] - content.get(vi).copied().unwrap_or(0));
        }

        added.clear();
        added.resize(clen + 1, 0);
        let mut ctx = RowCtx {
            lo,
            hi,
            up_hi,
            dn_lo,
            dn_hi,
            top,
            content: &content,
            cut: &cut,
            gap: &gap,
            added: &mut added,
            row: &mut row,
            scratch: &mut scratch,
            width,
            mult: *mult,
            out,
            overflow,
            produced: &mut produced,
        };
        fill_runs(lo, 1, &mut ctx);
    }
    PRODUCTIONS.fetch_add(produced, Ordering::Relaxed);
}

/// Scratch for filling one row of one layer state.
struct RowCtx<'a, C, K: LayerKey> {
    lo: usize,
    hi: usize,
    up_hi: usize,
    dn_lo: usize,
    dn_hi: usize,
    /// Largest value this row may use: one past the content length (a value v
    /// needs a v-1 banked by an *earlier* row, so the cap is fixed row-wide).
    top: u32,
    /// Content *before* this row; fixed while the row is filled.
    content: &'a [u32],
    /// `cut[v]`: first column whose cell above is ≥ v (`hi` if none).
    cut: &'a [usize],
    /// `gap[v]`: how many v's the ballot condition lets this row add.
    gap: &'a [u32],
    /// Occurrences of each value contributed by this row so far; exactly one
    /// slot past `content`, for the one new value a row may introduce.
    added: &'a mut [u32],
    row: &'a mut [u32],
    /// Whatever the key representation needs to assemble a successor.
    scratch: &'a mut K::Scratch,
    /// Bytes per serialized key element (see [`elem_width`]).
    width: usize,
    mult: C,
    out: &'a mut [Map<K, C>],
    /// Set when a merge overflows `C`; the expansion is then abandoned and
    /// rerun with a wider accumulator.
    overflow: &'a mut bool,
    /// Row fillings committed by this worker's chunk (see [`PRODUCTIONS`]).
    produced: &'a mut u64,
}

/// Fill columns `a ..` of the current row with runs of equal values, the runs
/// strictly increasing in value left to right.
///
/// A weakly increasing row *is* a sequence of such runs, so this enumerates
/// exactly the fillings a cell-at-a-time recursion would, but decides a whole
/// run per stack frame. Each constraint costs O(1) per run:
///
/// * **Column strictness.** The row above is weakly increasing, so the columns
///   whose cell above blocks a value v form the suffix [cut[v], up_hi). A run
///   of v starting at `a < up_hi` may extend to `cut[v]` and no further, and
///   `cut[v] = hi` when nothing blocks — which also lets the run spill into
///   the overhang [up_hi, hi) where there is no cell above. Columns at or past
///   `up_hi` are never blocked.
///
/// * **Ballot.** The reading word takes a row right to left, so every v in the
///   row is read before every v-1 in it: the v-1s available to justify this
///   row's v's are only those banked *before* the row, and the whole-row
///   condition collapses to `#v added ≤ gap[v] = before[v-1] − before[v]`
///   (values 1-indexed). That per-value cap is exact — see the module docs for
///   why nothing else is needed.
fn fill_runs<C: Acc, K: LayerKey>(a: usize, vmin: u32, ctx: &mut RowCtx<C, K>) {
    if a == ctx.hi {
        *ctx.produced += 1;
        K::commit(ctx);
        return;
    }
    for v in vmin..=ctx.top {
        // Furthest column a run of v could reach from here.
        let emax = if a >= ctx.up_hi {
            ctx.hi
        } else {
            ctx.cut[v as usize]
        };
        if emax <= a {
            // Blocked immediately — but a larger value may still fit: both
            // `cut` and `gap` are non-monotone in v.
            continue;
        }
        let mut cap = emax - a;
        if v >= 2 {
            let g = ctx.gap[v as usize] as usize;
            if g == 0 {
                continue;
            }
            cap = cap.min(g);
        }
        let vi = (v - 1) as usize;
        // `n` is a run length, not a degree — this is not a V6 sweep floor. An
        // empty run is the `continue` above, and `row[a - lo + n - 1]` has no
        // n = 0 form.
        for n in 1..=cap {
            ctx.row[a - ctx.lo + n - 1] = v;
            ctx.added[vi] = n as u32;
            fill_runs(a + n, v + 1, ctx);
        }
        ctx.added[vi] = 0;
    }
}

impl LrBackend for SkewLr {
    fn lr_coeff(&self, lambda: &Partition, mu: &Partition, nu: &Partition) -> u128 {
        if lambda.size() != mu.size() + nu.size() || !lambda.contains(mu) {
            return 0;
        }
        // A caller sweeping many λ against one (μ, ν) — the natural way to
        // read off a product — may already have paid for the whole product
        // expansion. That lives in the skew cache under the juxtaposed shape
        // (see `schur_product`); answering from it costs a lookup, where the
        // λ/μ route below would run one fresh traversal *per λ*. Peek only:
        // for a one-shot query the λ/μ shape (|ν| cells) is the cheaper
        // expansion, so nothing is computed speculatively here.
        let (outer, inner, _) = product_walk(mu, nu);
        if let Some(product) = crate::memo::skew_cache_peek(&outer, &inner) {
            return product
                .binary_search_by(|(p, _)| p.cmp(lambda))
                .map(|i| product[i].1)
                .unwrap_or(0);
        }
        // Expand whichever side is cheaper. |λ/μ| = |ν| and |λ/ν| = |μ|, and
        // c^λ_{μν} = c^λ_{νμ}, so peeling off the *larger* factor leaves the
        // smaller diagram to walk. Always expanding λ/μ is badly wrong when
        // |ν| ≫ |μ|: on c^λ_{μν} with λ = [13,12..2], μ = [5,4,3,2,1] it built
        // all 42 335 terms of a 75-cell expansion to read one coefficient
        // (91ms, against lrcalc's 4.4ms) where the 15-cell side answers at once.
        let (inner, want) = if mu.size() >= nu.size() {
            (mu, nu)
        } else {
            (nu, mu)
        };
        // Shared, not cloned: this reads one coefficient out of an expansion that
        // may hold millions of terms.
        let expansion = expand_skew_shared(lambda, inner);
        // `expand_skew` is sorted by content, which is what the sort is for.
        expansion
            .binary_search_by(|(p, _)| p.cmp(want))
            .map(|i| expansion[i].1)
            .unwrap_or(0)
    }

    /// The full Schur expansion of the product s_μ · s_ν.
    ///
    /// A product is a skew expansion of one disconnected shape. Place one
    /// factor up and to the right of the other so the two diagrams share no
    /// row and no column; the result is a shape whose fillings are exactly a
    /// filling of one alongside a filling of the other, hence
    /// s_{shape} = s_μ · s_ν. Expanding that one shape yields every λ in the
    /// product at once — no candidate sweep, and no per-λ call to `lr_coeff`.
    /// Which factor goes where, and whether the walk runs on the transposed
    /// diagram, is `product_walk`'s choice, calibrated in
    /// `docs/record/littlewood-richardson.md`. Returns only the nonzero terms,
    /// sorted by λ.
    ///
    /// Memoized under that shape, so a repeat costs the copy; the shared
    /// form is [`schur_product_shared`](LrBackend::schur_product_shared).
    fn schur_product(&self, mu: &Partition, nu: &Partition) -> Vec<(Partition, u128)> {
        (*self.schur_product_shared(mu, nu)).clone()
    }

    fn schur_product_shared(&self, mu: &Partition, nu: &Partition) -> Arc<Vec<(Partition, u128)>> {
        memoized_product(mu, nu, || None)
    }
}

/// The expansion of s_μ · s_ν, memoized under the shape [`product_walk`]
/// chooses; on a miss, `shortcut`'s answer, or the layer's walk of that shape
/// when it declines.
///
/// One table entry per product, whatever route computed it: the memo is
/// keyed by the shape, and a closed form or a counting route yields the same
/// expansion the walk would. That is what lets [`SkewLr::lr_coeff`]'s peek
/// answer a sweep of coefficients from a product any backend built, and what
/// makes a repeated product cost a lookup — `AutoLr` routes rectangle,
/// two-row and three-row products through here for that reason.
///
/// A `shortcut` that answers must return what
/// [`LrBackend::schur_product`] promises: nonzero terms, sorted by λ.
pub(crate) fn memoized_product(
    mu: &Partition,
    nu: &Partition,
    shortcut: impl FnOnce() -> Option<Vec<(Partition, u128)>>,
) -> Arc<Vec<(Partition, u128)>> {
    let (outer, inner, conjugate) = product_walk(mu, nu);
    skew_cached(&outer, &inner, || {
        shortcut().unwrap_or_else(|| expand_walk(&outer, &inner, conjugate))
    })
}

/// The skew shape a product `s_μ · s_ν` is expanded as, and whether the walk
/// runs on its transpose.
///
/// The shape is [`juxtapose`]d with the larger factor (by cells, then
/// lexicographically — a total order, so `s_μ·s_ν` and `s_ν·s_μ` share one
/// cache entry) on top: the upper block has a single canonical filling, so
/// the layer enumerates the smaller factor, with the larger as the ballot
/// offset. Four walks compute the same product — either factor enumerated,
/// on the diagram or its transpose, and transposing swaps the roles — and
/// they differ in how many row fillings the layer commits and what each one
/// costs (`examples/calibrate_orientation.rs`).
///
/// The transpose pays only where the direct fill is loose enough that there
/// is compression left to gain — a square or near-square, whose products carry
/// many tableaux per term — and costs elsewhere, because its content is longer
/// and each production dearer. Measured, the boundary is: the factors have
/// the same number of rows, the smaller is at least 0.7 of the larger by
/// cells, the product has at least 72 cells, and the shape is at least eight
/// rows and wider than tall (the same two conditions as
/// [`prefer_conjugate`]). Every asymmetric pair measured wants the direct
/// walk, by up to 2.9x; the fitted thresholds and their cases are in
/// `docs/record/littlewood-richardson.md`.
fn product_walk(mu: &Partition, nu: &Partition) -> (Partition, Partition, bool) {
    let (top, bottom) = if (mu.size(), mu.parts()) >= (nu.size(), nu.parts()) {
        (mu, nu)
    } else {
        (nu, mu)
    };
    let (outer, inner) = juxtapose(top, bottom);
    let rows = outer.len();
    let width = outer.part(0) as usize;
    let cells = top.size() + bottom.size();
    let conjugate = top.len() == bottom.len()
        && 10 * bottom.size() >= 7 * top.size()
        && cells >= 72
        && rows >= 8
        && width > rows;
    (outer, inner, conjugate)
}

/// The disconnected skew shape whose Schur function is s_μ · s_ν.
///
/// μ is translated right by ν₁ and stacked above ν; the inner shape is the
/// ν₁-wide block that translation leaves empty.
fn juxtapose(mu: &Partition, nu: &Partition) -> (Partition, Partition) {
    let shift = nu.part(0);
    let mut outer = Vec::with_capacity(mu.len() + nu.len());
    outer.extend(mu.parts().iter().map(|&p| p + shift));
    outer.extend_from_slice(nu.parts());
    let inner = vec![shift; if shift == 0 { 0 } else { mu.len() }];
    (Partition::from_sorted(outer), Partition::from_sorted(inner))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lr::NaiveLr;
    use crate::partition::partitions_of;

    fn p(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    /// What `expand_skew` must equal, computed the slow honest way.
    fn reference_skew(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
        if !outer.contains(inner) {
            return Vec::new();
        }
        let mut v: Vec<(Partition, u128)> = partitions_of(outer.size() - inner.size())
            .into_iter()
            .filter_map(|nu| {
                let c = NaiveLr.lr_coeff(outer, inner, &nu);
                (c != 0).then_some((nu, c))
            })
            .collect();
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }

    #[test]
    fn degenerate_shapes() {
        // inner ⊄ outer
        assert!(expand_skew(&p(&[2, 1]), &p(&[3])).is_empty());
        assert!(expand_skew(&p(&[]), &p(&[1])).is_empty());
        // outer == inner: s_{λ/λ} = 1
        for sh in [&[][..], &[1], &[3, 2, 1], &[4, 4]] {
            let l = p(sh);
            assert_eq!(expand_skew(&l, &l), vec![(Partition::default(), 1)]);
        }
        // empty inner: s_{λ/∅} = s_λ
        for sh in [&[1][..], &[3, 2, 1], &[2, 2, 2], &[5]] {
            let l = p(sh);
            assert_eq!(expand_skew(&l, &Partition::default()), vec![(l.clone(), 1)]);
        }
        // one row and one column
        assert_eq!(expand_skew(&p(&[5]), &p(&[2])), vec![(p(&[3]), 1)]);
        assert_eq!(
            expand_skew(&p(&[1, 1, 1, 1]), &p(&[1])),
            vec![(p(&[1, 1, 1]), 1)]
        );
    }

    /// A single column deeper than any small type would hold: the entries are
    /// forced to 1, 2, 3, … down the column, so this pins the cell-value range.
    #[test]
    fn very_deep_column() {
        let col = Partition::from_sorted(vec![1; 300]);
        assert_eq!(
            expand_skew(&col, &Partition::default()),
            vec![(col.clone(), 1)]
        );
        // …and with an inner shape lopping off the top of it.
        let inner = Partition::from_sorted(vec![1; 40]);
        assert_eq!(
            expand_skew(&col, &inner),
            vec![(Partition::from_sorted(vec![1; 260]), 1)]
        );
    }

    #[test]
    fn disconnected_shape() {
        // λ/μ = [4,4,2,2]/[2,2] splits into two 2x2 blocks that share no row
        // and no column, so s_{λ/μ} = s_{22}·s_{22}.
        let got = expand_skew(&p(&[4, 4, 2, 2]), &p(&[2, 2]));
        let want = NaiveLr.schur_product(&p(&[2, 2]), &p(&[2, 2]));
        assert_eq!(got, want);
        assert_eq!(got, reference_skew(&p(&[4, 4, 2, 2]), &p(&[2, 2])));
    }

    #[test]
    fn agrees_with_naive_on_every_small_skew_shape() {
        let mut checked = 0;
        for n in 0..=8u32 {
            for outer in partitions_of(n) {
                for k in 0..=n {
                    for inner in partitions_of(k) {
                        if !outer.contains(&inner) {
                            assert!(expand_skew(&outer, &inner).is_empty());
                            continue;
                        }
                        assert_eq!(
                            expand_skew(&outer, &inner),
                            reference_skew(&outer, &inner),
                            "s_{{{outer}/{inner}}}"
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 300, "expected a real sweep, got {checked}");
    }

    #[test]
    fn agrees_with_naive_on_every_small_product() {
        let mut checked = 0;
        for a in 0..=7u32 {
            for b in 0..=(7 - a) {
                for mu in partitions_of(a) {
                    for nu in partitions_of(b) {
                        assert_eq!(
                            SkewLr.schur_product(&mu, &nu),
                            NaiveLr.schur_product(&mu, &nu),
                            "s{mu}*s{nu}"
                        );
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 200, "expected a real sweep, got {checked}");
    }

    /// Shapes chosen to drive the run-based fill into its corners, each checked
    /// against the independent naive backend.
    ///
    /// The run fill decides a maximal block of equal values at once, so what
    /// can go wrong is a run that stops one column early or late. These shapes
    /// put a run against each boundary in turn:
    ///
    /// * an *overhang* — columns past the end of the row above, where nothing
    ///   blocks and a run may spill to the row end;
    /// * a row above that blocks in the middle (`cut` interior);
    /// * a row above that blocks at its very first column;
    /// * rows that share no column at all — empty overlap, so the layer's
    ///   `above` half is empty.
    #[test]
    fn run_fill_boundaries() {
        for (o, i) in [
            // Overhang: row 1 is far wider than row 0, so most of it has no
            // cell above and a single run may cover the tail.
            (&[9, 2][..], &[][..]),
            (&[9, 9, 2], &[]),
            // Inner shape shifts the overlap right, so `cut` starts mid-row.
            (&[9, 7], &[4]),
            (&[9, 7, 5], &[4, 2]),
            // Disjoint rows: no shared column, `above` is empty every step.
            (&[8, 3], &[5]),
            (&[8, 4, 2], &[6, 3]),
            // A single wide row, and a wide row over a wide row.
            (&[12], &[]),
            (&[12, 12], &[]),
            (&[12, 11], &[3]),
            // Narrow-over-wide and wide-over-narrow, both orientations of the
            // overlap clipping.
            (&[7, 3], &[]),
            (&[3, 3, 7], &[]),
            // Ragged, so every row has a different overlap on both sides.
            (&[10, 8, 5, 1], &[6, 3, 1]),
        ] {
            let (outer, inner) = (p(o), p(i));
            assert_eq!(
                expand_skew(&outer, &inner),
                reference_skew(&outer, &inner),
                "s_{{{outer}/{inner}}}"
            );
        }
    }

    /// Rows of width zero, which reach `finish_row` without entering a run.
    ///
    /// `outer_r == inner_r` makes row `r` empty; the traversal must still carry
    /// the layer through it rather than dropping or duplicating states.
    #[test]
    fn empty_rows_are_traversed() {
        for (o, i) in [
            (&[4, 3, 3, 1][..], &[0, 3, 0, 0][..]), // row 1 empty, mid-shape
            (&[4, 4], &[4, 0]),                     // row 0 empty
            (&[5, 3, 3], &[2, 3, 3]),               // two trailing empty rows
        ] {
            let (outer, inner) = (Partition::new(o.iter().copied()), p(i));
            assert_eq!(
                expand_skew(&outer, &inner),
                reference_skew(&outer, &inner),
                "s_{{{outer}/{inner}}}"
            );
        }
    }

    /// The packed layer key must distinguish states the old two-`Vec` key
    /// did, in particular a content/`above` split that could be read two ways.
    ///
    /// `[len, content.., above..]` is only unambiguous because of the leading
    /// length, so this pins that a shape whose `above` values collide with
    /// plausible content values still separates its states. A shape wide enough
    /// to carry several columns forward, over enough rows to grow the content
    /// past the `above` width, exercises exactly that.
    #[test]
    fn packed_key_separates_states() {
        for (o, i) in [
            (&[6, 5, 4, 3][..], &[2, 1][..]),
            (&[5, 5, 5, 5], &[]),
            (&[7, 6, 4, 2], &[3, 2, 1]),
        ] {
            let (outer, inner) = (p(o), p(i));
            let got = expand_skew(&outer, &inner);
            assert_eq!(got, reference_skew(&outer, &inner), "s_{{{outer}/{inner}}}");
            // Total tableau count is the sum of coefficients; a key collision
            // would merge states and lose some, a spurious split would keep
            // duplicates. Cross-check against the naive count.
            let total: u128 = got.iter().map(|(_, c)| c).sum();
            let want: u128 = reference_skew(&outer, &inner).iter().map(|(_, c)| c).sum();
            assert_eq!(total, want, "tableau total for {outer}/{inner}");
        }
    }

    #[test]
    fn lr_coeff_matches_naive_including_zeros() {
        for lambda in partitions_of(7) {
            for mu in partitions_of(3) {
                for nu in partitions_of(4) {
                    assert_eq!(
                        SkewLr.lr_coeff(&lambda, &mu, &nu),
                        NaiveLr.lr_coeff(&lambda, &mu, &nu),
                        "c^{lambda}_{{{mu},{nu}}}"
                    );
                }
            }
        }
        // Mismatched degree and non-containment are both 0.
        assert_eq!(SkewLr.lr_coeff(&p(&[3, 1]), &p(&[2]), &p(&[1])), 0);
        assert_eq!(SkewLr.lr_coeff(&p(&[2, 2]), &p(&[3]), &p(&[1])), 0);
    }

    /// `lr_coeff` may answer from a previously computed product expansion
    /// (the skew cache under the juxtaposed shape) instead of expanding λ/μ.
    /// Both routes must give identical answers, including the zeros — sweep
    /// every λ of the right degree cold, then again with the product warm.
    #[test]
    fn lr_coeff_agrees_before_and_after_product_is_cached() {
        let mu = p(&[4, 2, 1]);
        let nu = p(&[3, 2]);
        let n = mu.size() + nu.size();
        let cold: Vec<u128> = partitions_of(n)
            .iter()
            .map(|l| SkewLr.lr_coeff(l, &mu, &nu))
            .collect();
        let product = SkewLr.schur_product(&mu, &nu); // warms the cache
        assert_eq!(product, NaiveLr.schur_product(&mu, &nu));
        for (lambda, was) in partitions_of(n).iter().zip(&cold) {
            assert_eq!(
                SkewLr.lr_coeff(lambda, &mu, &nu),
                *was,
                "cached-product route diverged on c^{lambda}_{{{mu},{nu}}}"
            );
            // The peek finds the product under either argument order, since
            // the walk is a function of the unordered pair.
            assert_eq!(
                SkewLr.lr_coeff(lambda, &nu, &mu),
                *was,
                "cached-product route diverged on c^{lambda}_{{{nu},{mu}}}"
            );
            let want = product
                .iter()
                .find(|(l, _)| l == lambda)
                .map_or(0, |(_, c)| *c);
            assert_eq!(*was, want, "vs product term for {lambda}");
        }
    }

    /// Shapes that fire the orientation dispatch must give the same expansion
    /// through the conjugated walk as through the direct one.
    ///
    /// `[14,13,…,7]/[6,5,…,1]` (8 rows, width 14, 63 cells, every step one)
    /// fires the rule, is cheap even in debug builds, and its
    /// direct-orientation expansion is computed here explicitly as the
    /// reference.
    #[test]
    fn orientation_dispatch_preserves_expansions() {
        let outer = p(&[14, 13, 12, 11, 10, 9, 8, 7]);
        let inner = p(&[6, 5, 4, 3, 2, 1]);
        assert!(prefer_conjugate(&outer, &inner), "test shape must dispatch");
        let dispatched = expand_skew(&outer, &inner);
        let mut direct = expand_oriented::<u64>(&outer, &inner, false).expect("no overflow");
        direct.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(dispatched, direct);

        // The rule's stated boundaries, pinned so a future edit is deliberate.
        // Transposed: the step-1 staircases it was first measured on, at any
        // depth of the inner shape; and bands — the inner shape reaching to
        // within four rows of the bottom — whatever the steps, including the
        // 6.06x and 3.95x cases.
        assert!(prefer_conjugate(
            &p(&[16, 15, 14, 13, 12, 11, 10, 9]),
            &p(&[8, 7, 6, 5, 4, 3, 2, 1])
        ));
        assert!(prefer_conjugate(
            &p(&[12, 11, 10, 9, 8, 7, 6, 5, 4, 3]),
            &p(&[5, 4, 3, 2, 1])
        ));
        assert!(prefer_conjugate(
            &p(&[26, 23, 20, 17, 14, 11, 8, 5, 2]),
            &p(&[14, 12, 10, 8, 6, 4, 2])
        ));
        assert!(prefer_conjugate(
            &p(&[26, 23, 20, 17, 14, 11, 8, 5, 2]),
            &p(&[12, 10, 8, 6, 4])
        ));
        assert!(prefer_conjugate(
            &p(&[16, 15, 14, 12, 11, 10, 9, 8]),
            &p(&[7, 6, 5, 4, 3, 2, 1])
        ));
        // Direct: five or more full rows under a steep outer shape, and the
        // size floor.
        assert!(!prefer_conjugate(
            &p(&[20, 17, 14, 11, 8, 5, 2, 2]),
            &p(&[8, 5, 2])
        ));
        assert!(!prefer_conjugate(
            &p(&[16, 14, 12, 10, 8, 6, 4, 2]),
            &p(&[6, 4, 2])
        ));
        assert!(!prefer_conjugate(
            &p(&[30, 26, 22, 17, 13, 9, 5, 3, 1]),
            &p(&[20, 16, 12, 8])
        ));
        assert!(!prefer_conjugate(
            &p(&[24, 20, 16, 12, 8, 4, 2]),
            &p(&[12, 10, 8, 6, 4, 2])
        ));
        assert!(!prefer_conjugate(&p(&[40, 36, 32]), &Partition::default()));
        assert!(!prefer_conjugate(
            &p(&[12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12]),
            &Partition::default()
        ));
        assert!(!prefer_conjugate(
            &p(&[9, 7, 5, 3, 2, 2, 1, 1]),
            &Partition::default()
        ));
    }

    /// A product's four walks — either factor enumerated, on the diagram or
    /// its transpose — must agree, and [`product_walk`] must pick the one it
    /// says it picks. One pair from each regime: a near-square of rectangles
    /// (transposed, and cheap because the product is multiplicity-free) and
    /// an asymmetric pair (direct). The dispatched product is compared with
    /// every walk it did not take.
    #[test]
    fn product_walk_agrees_with_every_other_walk() {
        for (a, b, conj) in [
            (&[10u32, 10, 10, 10][..], &[9u32, 9, 9, 9][..], true),
            (&[10, 8, 6, 4], &[8, 6, 4, 2], false),
        ] {
            let (a, b) = (p(a), p(b));
            let (outer, inner, chosen) = product_walk(&a, &b);
            assert_eq!(chosen, conj, "dispatch for s{a}·s{b}");
            let got = SkewLr.schur_product(&a, &b);
            let (swapped_outer, swapped_inner) = juxtapose(&b, &a);
            for (o, i, c) in [
                (&outer, &inner, !conj),
                (&swapped_outer, &swapped_inner, false),
                (&swapped_outer, &swapped_inner, true),
            ] {
                assert_eq!(
                    expand_walk(o, i, c),
                    got,
                    "walk {o}/{i} conj={c} on s{a}·s{b}"
                );
            }
        }
    }

    /// The walk is a function of the unordered pair — the same shape, so the
    /// same cache entry, whichever way the product is written — and its
    /// transpose decision must match the calibration it was fitted to
    /// (`docs/record/littlewood-richardson.md`, "product orientation").
    /// Squares from the earlier calibration of `prefer_conjugate` are pinned
    /// with the outcome they measured today.
    #[test]
    fn product_walk_is_symmetric_and_matches_the_calibration() {
        let cases: &[(&[u32], &[u32], bool)] = &[
            // Asymmetric: direct, by 1.5–2.9x over the transpose.
            (&[16, 13, 10, 7], &[8, 6, 4, 2], false),
            (&[20, 16, 12, 8], &[8, 6, 4, 2], false),
            (&[18, 15, 12, 9], &[9, 7, 5, 3], false),
            (&[24, 20, 16, 12], &[10, 8, 6, 4], false),
            (&[20, 16, 12, 8], &[10, 8, 6, 4], false), // equal rows, ratio 0.5
            (&[16, 13, 10, 7], &[10, 8, 6, 4], false), // equal rows, ratio 0.61
            (&[16, 13, 10, 7], &[5, 4, 3, 2, 1], false),
            (&[12, 11, 10, 9, 8, 7], &[6, 5, 4, 3], false),
            (&[16, 13, 10, 7], &[8, 7, 6, 5, 4, 3], false),
            (&[12, 10, 8, 6], &[9, 8, 7, 6, 5], false), // ratio 0.97, rows differ
            (&[16, 13, 10, 7], &[10, 9, 8, 7, 6, 5], false),
            // Near-squares: transpose, by 1.06–1.55x.
            (&[16, 13, 10, 7], &[12, 10, 8, 6], true),
            (&[10, 9, 8, 7, 6, 5], &[9, 8, 7, 6, 5, 4], true),
            (&[10, 9, 8, 7, 6, 5], &[8, 7, 6, 5, 4, 3], true),
            // Near-squares under 72 cells: direct (a' loses 1.31–1.49x).
            (&[12, 10, 8, 6], &[10, 8, 6, 4], false),
            (&[11, 9, 7, 5], &[10, 8, 6, 4], false),
            (&[9, 8, 7, 6, 5], &[9, 7, 5, 3, 1], false),
            // Squares: transpose from 72 cells (1.13–2.25x); direct below
            // (`[12,9,6,3]²` at 60 cells loses 1.39x transposed).
            (&[16, 13, 10, 7], &[16, 13, 10, 7], true),
            (&[12, 10, 8, 6], &[12, 10, 8, 6], true),
            (&[10, 9, 8, 7, 6, 5], &[10, 9, 8, 7, 6, 5], true),
            (&[12, 9, 6, 3], &[12, 9, 6, 3], false),
            (&[10, 8, 6, 4], &[10, 8, 6, 4], false),
            // Too few rows, or taller than wide: direct regardless of size.
            (&[40, 32, 24], &[40, 32, 24], false),
            (
                &[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                &[3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
                false,
            ),
        ];
        for (a, b, want) in cases {
            let (a, b) = (p(a), p(b));
            let (outer, inner, conj) = product_walk(&a, &b);
            assert_eq!(conj, *want, "s{a}·s{b}");
            let (outer_r, inner_r, conj_r) = product_walk(&b, &a);
            assert_eq!(
                (outer, inner, conj),
                (outer_r, inner_r, conj_r),
                "s{b}·s{a}"
            );
        }
        // The larger factor sits on top, so the smaller is the enumerated
        // block: the outer shape's first rows are the larger factor shifted.
        let (outer, inner, _) = product_walk(&p(&[3, 2]), &p(&[6, 4, 2]));
        assert_eq!(outer, p(&[9, 7, 5, 3, 2]));
        assert_eq!(inner, p(&[3, 3, 3]));
    }

    /// Keys are serialized at the narrowest element width the cell count
    /// allows, so the dangerous inputs are the ones that sit exactly on a
    /// width boundary. Pieri gives the exact answer independently of any LR
    /// machinery: s_a · s_1 = s_{a+1} + s_{a,1}.
    ///
    /// n = 255 is the last one-byte shape (a content count hits 0xFF exactly),
    /// n = 256 the first two-byte one; 65535/65536 likewise for two → four.
    #[test]
    fn width_boundaries_match_pieri() {
        for a in [254u32, 255, 256, 65534, 65535] {
            let got = SkewLr.schur_product(&p(&[a]), &p(&[1]));
            let mut want = vec![(p(&[a + 1]), 1), (p(&[a, 1]), 1)];
            want.sort_by(|x, y| x.0.cmp(&y.0));
            assert_eq!(got, want, "s[{a}]·s[1]");
        }
    }

    /// The traversal reports terms in map order, which legitimately varies
    /// with the key representation; sort before comparing.
    fn sorted(v: Option<Vec<(Partition, u128)>>) -> Vec<(Partition, u128)> {
        let mut v = v.expect("no overflow at u64");
        v.sort_by(|a, b| a.0.cmp(&b.0));
        v
    }

    /// The three byte serializations must be interchangeable: any width wide
    /// enough for the shape yields the same expansion. Forcing 2 and 4 bytes
    /// onto one-byte shapes exercises every encode/decode pair on layers with
    /// real merging, where a mis-split key would corrupt coefficients.
    #[test]
    fn widths_agree_on_merging_shapes() {
        for (o, i) in [
            (&[10, 8, 5, 1][..], &[6, 3, 1][..]),
            (&[6, 5, 4, 3], &[2, 1]),
            (&[7, 6, 4, 2], &[3, 2, 1]),
        ] {
            let (outer, inner) = (p(o), p(i));
            let narrow = sorted(expand_layer::<u64, Key>(&outer, &inner, 1, false));
            let mid = sorted(expand_layer::<u64, Key>(&outer, &inner, 2, false));
            let wide = sorted(expand_layer::<u64, Key>(&outer, &inner, 4, false));
            assert_eq!(narrow, mid, "widths 1 vs 2 on {outer}/{inner}");
            assert_eq!(narrow, wide, "widths 1 vs 4 on {outer}/{inner}");
            assert_eq!(narrow, reference_skew(&outer, &inner), "vs naive");
        }
    }

    /// The bitmap keys and the byte key must be interchangeable: every
    /// representation that fits a shape yields the same expansion as every
    /// other, and as the reference where the reference is affordable.
    ///
    /// The shapes carry real merging (a mis-encoded state would merge with
    /// the wrong neighbors and move a coefficient rather than fail), and the
    /// last two fire the conjugate dispatch and a large layer respectively.
    #[test]
    fn key_representations_agree_on_merging_shapes() {
        let small: &[(&[u32], &[u32])] = &[
            (&[10, 8, 5, 1], &[6, 3, 1]),
            (&[6, 5, 4, 3], &[2, 1]),
            (&[7, 6, 4, 2], &[3, 2, 1]),
            (&[9, 7, 5], &[4, 2]),
        ];
        for &(o, i) in small {
            let (outer, inner) = (p(o), p(i));
            let w = elem_width(&outer, &inner);
            let bytes = sorted(expand_layer::<u64, Key>(&outer, &inner, w, false));
            let narrow = sorted(expand_layer::<u64, PackedKey<u64>>(
                &outer, &inner, w, false,
            ));
            let wide = sorted(expand_layer::<u64, PackedKey<W128>>(
                &outer, &inner, w, false,
            ));
            assert_eq!(narrow, bytes, "u64 bitmap vs bytes on {outer}/{inner}");
            assert_eq!(wide, bytes, "128-bit bitmap vs bytes on {outer}/{inner}");
            assert_eq!(bytes, reference_skew(&outer, &inner), "vs naive");
        }
        for m in [
            &[6u32, 5, 4, 3, 2][..],
            &[7, 6, 5, 4, 3],
            &[10, 10, 10, 10, 10],
        ] {
            let (outer, inner) = juxtapose(&p(m), &p(m));
            let w = elem_width(&outer, &inner);
            let bytes = sorted(expand_layer::<u64, Key>(&outer, &inner, w, false));
            let narrow = sorted(expand_layer::<u64, PackedKey<u64>>(
                &outer, &inner, w, false,
            ));
            let wide = sorted(expand_layer::<u64, PackedKey<W128>>(
                &outer, &inner, w, false,
            ));
            assert_eq!(narrow, bytes, "u64 bitmap vs bytes on s{}²", p(m));
            assert_eq!(wide, bytes, "128-bit bitmap vs bytes on s{}²", p(m));
            // …and through the dispatch, which may conjugate first.
            assert_eq!(
                SkewLr.schur_product(&p(m), &p(m)),
                bytes,
                "dispatched s{}²",
                p(m)
            );
        }
    }

    /// Both bitmaps are bijections onto their words: encoding then decoding
    /// gives the sequence back, and distinct sequences never share a word.
    /// The extremes place a bit at the top position of each word width, which
    /// is where an off-by-one in the bound would wrap or panic.
    #[test]
    fn packed_bitmaps_round_trip_and_separate() {
        fn pack_partition<W: Word>(parts: &[u32]) -> W {
            let k = parts.len();
            let mut w = W::ZERO;
            for (i, &c) in parts.iter().enumerate() {
                w = w.set(c + (k - 1 - i) as u32);
            }
            w
        }
        fn pack_row<W: Word>(cells: &[u32]) -> W {
            let mut w = W::ZERO;
            for (j, &v) in cells.iter().enumerate() {
                w = w.set(v + j as u32);
            }
            w
        }
        let mut buf = Vec::new();

        // Every partition in a 6 × 6 box, on both words.
        let mut seen64 = std::collections::HashSet::new();
        let mut seen128 = std::collections::HashSet::new();
        let mut count = 0;
        for n in 0..=36u32 {
            for q in partitions_of(n) {
                if q.len() > 6 || q.part(0) > 6 {
                    continue;
                }
                let w64: u64 = pack_partition(q.parts());
                unpack_partition(w64, &mut buf);
                assert_eq!(buf, q.parts(), "u64 round trip of {q}");
                assert!(seen64.insert(w64), "u64 word collision at {q}");
                let w128: W128 = pack_partition(q.parts());
                unpack_partition(w128, &mut buf);
                assert_eq!(buf, q.parts(), "128-bit round trip of {q}");
                assert!(seen128.insert(w128.0), "128-bit word collision at {q}");
                count += 1;
            }
        }
        assert_eq!(count, 924, "partitions in a 6 × 6 box");

        // Every weakly increasing row of up to 4 cells with values in 1..=5.
        fn rows_of(d: usize, min: u32, cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
            if cur.len() == d {
                out.push(cur.clone());
                return;
            }
            for v in min..=5 {
                cur.push(v);
                rows_of(d, v, cur, out);
                cur.pop();
            }
        }
        let mut rows = Vec::new();
        for d in 0..=4 {
            rows_of(d, 1, &mut Vec::new(), &mut rows);
        }
        assert_eq!(
            rows.len(),
            1 + 5 + 15 + 35 + 70,
            "multisets of size ≤ 4 from 5 values"
        );
        let mut seen = std::collections::HashSet::new();
        for cells in &rows {
            let w: u64 = pack_row(cells);
            unpack_row(w, &mut buf);
            assert_eq!(&buf, cells, "row round trip of {cells:?}");
            assert!(
                seen.insert((cells.len(), w)),
                "row word collision at {cells:?}"
            );
        }

        // The top bit of each word: 40 parts of 24 put bit 63 in play, 64
        // parts of 64 put bit 127, and eight cells of 56 end at bit 63.
        let tall = vec![24u32; 40];
        unpack_partition::<u64>(pack_partition(&tall), &mut buf);
        assert_eq!(buf, tall);
        let taller = vec![64u32; 64];
        unpack_partition::<W128>(pack_partition(&taller), &mut buf);
        assert_eq!(buf, taller);
        let row = vec![56u32; 8];
        unpack_row::<u64>(pack_row(&row), &mut buf);
        assert_eq!(buf, row, "row bits up to position 63");
    }

    /// The representation is chosen from `rows + width` alone, so the shapes
    /// on either side of each boundary must dispatch as stated and still
    /// expand correctly — a straight shape has one filling, `s_{λ/∅} = s_λ`,
    /// and its single term's content puts a bit at the word's top position.
    #[test]
    fn packed_key_dispatch_boundaries() {
        let e = Partition::default();
        // rows + width = 64: the last u64 shape; content [62, 3] is bit 63.
        let l = p(&[62, 3]);
        assert!(PackedKey::<u64>::fits(&l, &e));
        assert_eq!(expand_skew(&l, &e), vec![(l.clone(), 1)]);
        // 65: the first 128-bit shape.
        let l = p(&[63, 3]);
        assert!(!PackedKey::<u64>::fits(&l, &e));
        assert!(PackedKey::<W128>::fits(&l, &e));
        assert_eq!(expand_skew(&l, &e), vec![(l.clone(), 1)]);
        // 128: the last 128-bit shape; content [126, 3] is bit 127.
        let l = p(&[126, 3]);
        assert!(PackedKey::<W128>::fits(&l, &e));
        assert_eq!(expand_skew(&l, &e), vec![(l.clone(), 1)]);
        // 129: bytes.
        let l = p(&[127, 3]);
        assert!(!PackedKey::<W128>::fits(&l, &e));
        assert_eq!(expand_skew(&l, &e), vec![(l.clone(), 1)]);

        // Pieri astride the u64 boundary, with real branching in the last row:
        // s_a · s_1 = s_{a+1} + s_{a,1}, and the juxtaposed shape has
        // rows + width = a + 3.
        for a in [60u32, 61, 62, 63] {
            let got = SkewLr.schur_product(&p(&[a]), &p(&[1]));
            let mut want = vec![(p(&[a + 1]), 1), (p(&[a, 1]), 1)];
            want.sort_by(|x, y| x.0.cmp(&y.0));
            assert_eq!(got, want, "s[{a}]·s[1]");
        }
    }

    /// The byte fallback keeps working, both sides of its inline capacity.
    ///
    /// `[140, 136]/∅` carries a 136-wide clipped row (all heap) and has
    /// exactly one filling, so the answer is pinned: s_{λ/∅} = s_λ. `[33,31]²`
    /// mixes inline and heap byte keys in one layer *with* merging: forced
    /// onto the byte key, it is checked against the dispatched expansion
    /// (a 128-bit bitmap, its `rows + width` being 70) and against the
    /// conjugate orientation, an independent traversal of the same
    /// coefficients via c^λ_{μν} = c^{λ'}_{μ'ν'}.
    #[test]
    fn byte_keys_merge_and_agree_with_conjugate_orientation() {
        let wide = p(&[140, 136]);
        assert!(!PackedKey::<W128>::fits(&wide, &Partition::default()));
        assert_eq!(
            expand_skew(&wide, &Partition::default()),
            vec![(wide.clone(), 1)]
        );

        let m = p(&[33, 31]);
        let (outer, inner) = juxtapose(&m, &m);
        let w = elem_width(&outer, &inner);
        let bytes = sorted(expand_layer::<u64, Key>(&outer, &inner, w, false));
        let direct = SkewLr.schur_product(&m, &m);
        assert_eq!(direct, bytes, "dispatched vs forced byte keys");
        let mc = m.conjugate();
        let mut via_conjugate: Vec<(Partition, u128)> = SkewLr
            .schur_product(&mc, &mc)
            .into_iter()
            .map(|(l, c)| (l.conjugate(), c))
            .collect();
        via_conjugate.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(direct, via_conjugate);
    }

    // A `u8` accumulator, so the saturation boundary is cheap to reach. Written
    // here rather than inside the test body: an `impl` is never scoped, even
    // nested in a function, so the in-body form implemented `Acc for u8` across
    // the whole test module while looking local.
    impl Acc for u8 {
        const ONE: Self = 1;
        fn checked_add(self, other: Self) -> Option<Self> {
            u8::checked_add(self, other)
        }
        fn widen(self) -> u128 {
            self as u128
        }
    }

    /// Layer multiplicities are tableau counts, so no narrow accumulator is
    /// provably safe — the engine must *detect* saturation and retry wider.
    /// The `u8` accumulator above makes the boundary cheap to reach: the pass
    /// must report overflow (not wrap), and the two production widths must
    /// agree.
    #[test]
    fn accumulator_overflow_is_detected_not_wrapped() {
        let sorted = |v: Option<Vec<(Partition, u128)>>| {
            v.map(|mut v| {
                v.sort_by(|a, b| a.0.cmp(&b.0));
                v
            })
        };

        // s[6,5,4,3,2]² has a coefficient of 644, and a final state's
        // multiplicity *is* its coefficient, so the u8 pass must saturate…
        let (outer, inner) = juxtapose(&p(&[6, 5, 4, 3, 2]), &p(&[6, 5, 4, 3, 2]));
        assert_eq!(expand_oriented::<u8>(&outer, &inner, false), None);
        let full = sorted(expand_oriented::<u64>(&outer, &inner, false));
        assert!(full.is_some(), "u64 must not saturate here");
        assert_eq!(full, sorted(expand_oriented::<u128>(&outer, &inner, false)));

        // …while a tiny shape stays under 255, so u8 must also *succeed* when
        // nothing overflows (the detection is not a blanket refusal).
        let (outer, inner) = juxtapose(&p(&[2, 1]), &p(&[2, 1]));
        let small = sorted(expand_oriented::<u8>(&outer, &inner, false));
        assert!(small.is_some());
        assert_eq!(
            small,
            sorted(expand_oriented::<u128>(&outer, &inner, false))
        );
    }
}
