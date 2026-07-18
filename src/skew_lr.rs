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
//! advances a *frontier* of merged partial fillings rather than a stack of
//! individual ones: two partial fillings that agree on the row above and on the
//! content so far are interchangeable, so they collapse into one weighted
//! state. That is where the asymptotic win lives; the single enumeration only
//! removes the `p(n)` factor on top of it.
//!
//! Products reduce to the same primitive. See [`SkewLr::schur_product`].

use std::borrow::Borrow;
use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hash, Hasher};
use std::sync::atomic::{AtomicUsize, Ordering};

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
/// Memoized on the shape (see [`crate::memo::skew_cached`]), so a
/// caller sweeping many ν against one (outer, inner) pays for one traversal.
pub fn expand_skew(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
    if !outer.contains(inner) {
        return Vec::new();
    }
    (*skew_cached(outer, inner, || expand_skew_uncached(outer, inner))).clone()
}

/// A partially-filled diagram, reduced to what the rest of the fill can see.
///
/// Two partial fillings that agree on this state are interchangeable:
/// everything still to be decided depends on the row above (column strictness)
/// and on the content accumulated so far (the ballot condition), and on nothing
/// else. So they are merged and carry a multiplicity, which is what keeps the
/// traversal from degenerating into an enumeration of individual tableaux.
///
/// The state is the sequence `[len, content.., above..]`: `content[v - 1]` is
/// the number of `v`s placed so far (no trailing zeros, so the leading `len`
/// makes the split unambiguous), and `above` holds the previous row's entries,
/// clipped to the columns the next row overlaps.
///
/// Peak memory is (number of live states) × (bytes per state), and this type
/// is the second factor, so it is packed hard:
///
/// * Every element — the length header, each content count, each cell value —
///   is at most the cell count of the shape (see [`elem_width`]), so the
///   sequence is serialized at the narrowest sufficient byte width, fixed per
///   expansion. For every practically computable shape that is one byte per
///   element, a 4× cut over the old `u32` words.
/// * Keys of ≤ [`INLINE`] bytes — all of them, in practice — are stored inline
///   in the enum, so a *new* state costs no heap allocation at all. The old
///   representation paid a malloc per state, and on shapes with multi-million
///   frontiers the allocator (and the kernel behind it) was a measured 38% of
///   wall time. Longer keys spill to a box and everything still works.
///
/// Merges (the common case) are probed with a borrowed scratch buffer via
/// [`KeyBytes`], so the hot path allocates nothing either way.
enum Key {
    /// `(length, bytes)`; only `bytes[..length]` is meaningful.
    Inline(u8, [u8; INLINE]),
    Heap(Box<[u8]>),
}

/// Inline capacity, chosen so `size_of::<Key>()` is 32: tag + 1 + 30 on one
/// side, a 16-byte box on the other. At one byte per element this holds a
/// header plus content plus a 20-wide clipped row with room to spare.
const INLINE: usize = 30;

const _: () = assert!(std::mem::size_of::<Key>() == 32);

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

/// Borrowed view of a [`Key`], so a frontier probe can use a scratch slice.
///
/// `#[repr(transparent)]` makes the `&[u8]` → `&KeyBytes` cast sound; `Hash`
/// and `Eq` agree with [`Key`]'s exactly, which is what `Borrow` requires.
#[derive(PartialEq, Eq)]
#[repr(transparent)]
struct KeyBytes([u8]);

impl KeyBytes {
    #[inline]
    fn new(s: &[u8]) -> &KeyBytes {
        // SAFETY: KeyBytes is repr(transparent) over [u8].
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

/// A multiply-xor-rotate hasher.
///
/// The frontier map is the hot data structure — several million short integer
/// keys per expansion — and SipHash's per-key setup dominates there. Keys are
/// small and structured, so a cheap mixing step is a large net win; measured at
/// roughly 1.3× on `s[8,7,6,5,4,3]²`. Not cryptographic, and nothing here is
/// exposed to adversarial input.
#[derive(Default)]
struct MixHasher(u64);

const MIX: u64 = 0x9e37_79b9_7f4a_7c15;

impl MixHasher {
    #[inline]
    fn add(&mut self, w: u64) {
        self.0 = (self.0.rotate_left(5) ^ w).wrapping_mul(MIX);
    }
}

impl Hasher for MixHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        let mut chunks = bytes.chunks_exact(8);
        for c in &mut chunks {
            self.add(u64::from_le_bytes(c.try_into().unwrap()));
        }
        let rest = chunks.remainder();
        if !rest.is_empty() {
            let mut buf = [0u8; 8];
            buf[..rest.len()].copy_from_slice(rest);
            self.add(u64::from_le_bytes(buf));
        }
    }
    #[inline]
    fn write_u8(&mut self, i: u8) {
        self.add(i as u64);
    }
    #[inline]
    fn write_u32(&mut self, i: u32) {
        self.add(i as u64);
    }
    #[inline]
    fn write_u64(&mut self, i: u64) {
        self.add(i);
    }
    #[inline]
    fn write_usize(&mut self, i: usize) {
        self.add(i as u64);
    }
    #[inline]
    fn finish(&self) -> u64 {
        // Avalanche (the murmur3/splitmix-style finalizer). The per-chunk mix
        // above only moves entropy *upward* (multiply and a small left
        // rotation), which was fine when a key spanned many chunks and the
        // rotation wrapped, but byte-packed keys fit in 2–4 chunks and left
        // the low bits — exactly the ones hashbrown takes the bucket index
        // from — barely mixed. That showed up as probe-chain clustering worth
        // ~2.7× on `[20,16,12]²`; these two multiply-xorshift rounds fix the
        // dispersion for a few cycles per key.
        let mut x = self.0;
        x ^= x >> 33;
        x = x.wrapping_mul(0xff51_afd7_ed55_8ccd);
        x ^= x >> 33;
        x = x.wrapping_mul(0xc4ce_b9fe_1a85_ec53);
        x ^ (x >> 33)
    }
}

type Map<K, V> = HashMap<K, V, BuildHasherDefault<MixHasher>>;

/// High-water mark of live frontier states, for measurement harnesses.
///
/// Peak memory is (states) × (bytes per state); this records the first factor,
/// which — unlike RSS — is not perturbed by the allocator. Sampled once per
/// row at the point both the old and new frontier are fully populated, so it
/// is the true in-traversal maximum of live map entries. Monotone across
/// expansions until read: [`take_peak_frontier_states`] returns and resets it.
static PEAK_LIVE_STATES: AtomicUsize = AtomicUsize::new(0);

/// Read and reset the peak live frontier-state count (see
/// [`PEAK_LIVE_STATES`]). A measurement hook, not part of the semantic API.
pub fn take_peak_frontier_states() -> usize {
    PEAK_LIVE_STATES.swap(0, Ordering::Relaxed)
}

/// The byte width every element of a frontier key fits in, for this shape.
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

/// A frontier multiplicity: `u64` for the fast pass, `u128` for the fallback.
///
/// Multiplicities are *tableau counts*, which dwarf the final coefficients
/// (2.1 × 10⁸ tableaux behind 10⁵-ish coefficients on `s[8,7,6,5,4,3]²`), so
/// no fixed narrow type is provably safe. Instead every add is checked: the
/// `u64` pass detects saturation and [`expand_skew_uncached`] transparently
/// reruns in `u128` — correctness never rests on an unproven bound, and the
/// 8-bytes-per-state saving is kept on every shape that stays under 2⁶⁴.
trait Acc: Copy {
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
    let mut out = expand_oriented::<u64>(outer, inner)
        .or_else(|| expand_oriented::<u128>(outer, inner))
        .expect("LR tableau multiplicity exceeded u128");
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

/// One frontier traversal with multiplicities in `C`; `None` means some merge
/// overflowed `C` and the caller should retry wider.
fn expand_oriented<C: Acc>(outer: &Partition, inner: &Partition) -> Option<Vec<(Partition, u128)>> {
    expand_with_width::<C>(outer, inner, elem_width(outer, inner))
}

/// [`expand_oriented`] at an explicit element width. Split out so tests can
/// force a wider serialization than [`elem_width`] picks and pin agreement
/// across the width-specific encode/decode paths; `width` must be sufficient
/// for the shape (any width is, when ≥ the chosen one).
fn expand_with_width<C: Acc>(
    outer: &Partition,
    inner: &Partition,
    width: usize,
) -> Option<Vec<(Partition, u128)>> {
    let rows = outer.len();
    let trace = std::env::var_os("SKEW_TRACE").is_some();

    // The frontier of the traversal: reduced state -> number of ways to reach
    // it. Between rows it is held as a plain `Vec`: the hash table is only
    // needed on the side being merged *into*, and a table's bucket array
    // (power-of-two, reserved ahead) can run 2–4× the entry payload. Draining
    // each finished table into an exactly-sized vector caps the steady-state
    // frontier at real entries only, and reading it back is a linear scan
    // instead of a table walk.
    let mut cur: Vec<(Key, C)> = vec![(Key::from_bytes(&vec![0u8; width]), C::ONE)];

    // Per-row scratch, reused across states so the inner loop never allocates.
    let mut row: Vec<u32> = Vec::new();
    let mut added: Vec<u32> = Vec::new();
    let mut cut: Vec<usize> = Vec::new();
    let mut gap: Vec<u32> = Vec::new();
    let mut key: Vec<u32> = Vec::new();
    let mut kbuf: Vec<u8> = Vec::new();
    let mut st: Vec<u32> = Vec::new();
    let mut overflow = false;

    for r in 0..rows {
        let lo = inner.part(r) as usize;
        let hi = outer.part(r) as usize;
        // Columns of row r-1 that row r sits under, and the ones of row r that
        // row r+1 will sit under. Clipping to these is what makes distinct
        // histories collapse.
        let (up_lo, up_hi) = overlap(inner, outer, r.wrapping_sub(1), r);
        let (dn_lo, dn_hi) = overlap(inner, outer, r, r + 1);

        let mut next: Map<Key, C> =
            Map::with_capacity_and_hasher(cur.len(), Default::default());
        row.clear();
        row.resize(hi.saturating_sub(lo), 0);
        for (state, mult) in cur.iter() {
            decode_into(&mut st, state.bytes(), width);
            let clen = st[0] as usize;
            let content = &st[1..1 + clen];
            let above = &st[1 + clen..];
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
                content,
                cut: &cut,
                gap: &gap,
                added: &mut added,
                row: &mut row,
                key: &mut key,
                kbuf: &mut kbuf,
                width,
                mult: *mult,
                out: &mut next,
                overflow: &mut overflow,
            };
            fill_runs(lo, 1, &mut ctx);
        }
        // Both frontiers are momentarily live here; record the sum (once per
        // row, so the cost is nil).
        PEAK_LIVE_STATES.fetch_max(cur.len() + next.len(), Ordering::Relaxed);
        if trace {
            eprintln!("skew_lr row {r}: {} -> {} states", cur.len(), next.len());
        }
        if overflow {
            return None;
        }
        cur = next.into_iter().collect();
        if cur.is_empty() {
            break;
        }
    }

    // Every row's down-overlap with a nonexistent next row is empty, so the
    // final keys are pure content — already merged, one term each.
    Some(
        cur.into_iter()
            .map(|(key, c)| {
                decode_into(&mut st, key.bytes(), width);
                debug_assert_eq!(st[0] as usize + 1, st.len());
                // The ballot condition forces the content to be weakly
                // decreasing at every prefix, so the final one is already a
                // partition.
                (Partition::from_sorted(st[1..].to_vec()), c.widen())
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

/// Scratch for filling one row of one frontier state.
struct RowCtx<'a, C> {
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
    /// Occurrences of each value contributed by this row so far.
    added: &'a mut [u32],
    row: &'a mut [u32],
    /// Scratch buffer for the candidate output key.
    key: &'a mut Vec<u32>,
    /// Scratch for the serialized form of `key`.
    kbuf: &'a mut Vec<u8>,
    /// Bytes per serialized key element (see [`elem_width`]).
    width: usize,
    mult: C,
    out: &'a mut Map<Key, C>,
    /// Set when a merge overflows `C`; the expansion is then abandoned and
    /// rerun with a wider accumulator.
    overflow: &'a mut bool,
}

/// Fill columns `a ..` of the current row with runs of equal values, the runs
/// strictly increasing in value left to right.
///
/// A weakly increasing row *is* a sequence of such runs, so this enumerates
/// exactly the fillings the old cell-at-a-time recursion did, but decides a
/// whole run per stack frame. Each constraint costs O(1) per run:
///
/// * **Column strictness.** The row above is weakly increasing, so the columns
///   whose cell above blocks a value v form the suffix [cut[v], up_hi) — a run
///   of v starting at `a < up_hi` may extend to `cut[v]` and no further
///   (`cut[v] = hi` when nothing blocks, which also lets the run spill into
///   the overhang [up_hi, hi) where there is no cell above). Columns at or past
///   `up_hi` are never blocked.
///
/// * **Ballot.** The reading word takes a row right to left, so every v in the
///   row is read before every v-1 in it: the v-1s available to justify this
///   row's v's are only those banked *before* the row, and the whole-row
///   condition collapses to `#v added ≤ gap[v] = before[v-1] − before[v]`
///   (values 1-indexed). That per-value cap is exact — see the module docs for
///   why nothing else is needed.
fn fill_runs<C: Acc>(a: usize, vmin: u32, ctx: &mut RowCtx<C>) {
    if a == ctx.hi {
        finish_row(ctx);
        return;
    }
    for v in vmin..=ctx.top {
        // Furthest column a run of v could reach from here.
        let emax = if a >= ctx.up_hi { ctx.hi } else { ctx.cut[v as usize] };
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
        for n in 1..=cap {
            ctx.row[a - ctx.lo + n - 1] = v;
            ctx.added[vi] = n as u32;
            fill_runs(a + n, v + 1, ctx);
        }
        ctx.added[vi] = 0;
    }
}

/// Commit a completed row: fold its content in, clip the frontier, merge.
fn finish_row<C: Acc>(ctx: &mut RowCtx<C>) {
    // Assemble the successor key in the scratch buffer: [len, content, above].
    let key = &mut *ctx.key;
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
    encode_into(ctx.kbuf, key, ctx.width);

    // Probe with the borrowed scratch; a merge (the common case) allocates
    // nothing, and a genuinely new state is copied inline unless oversized.
    match ctx.out.get_mut(KeyBytes::new(ctx.kbuf)) {
        Some(w) => match w.checked_add(ctx.mult) {
            Some(sum) => *w = sum,
            None => *ctx.overflow = true,
        },
        None => {
            ctx.out.insert(Key::from_bytes(ctx.kbuf), ctx.mult);
        }
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
        let (a, b) = if mu >= nu { (mu, nu) } else { (nu, mu) };
        let (outer, inner) = juxtapose(a, b);
        if let Some(product) = crate::memo::skew_cache_peek(&outer, &inner) {
            return product
                .binary_search_by(|(p, _)| p.cmp(lambda))
                .map(|i| product[i].1)
                .unwrap_or(0);
        }
        let expansion = expand_skew(lambda, mu);
        // `expand_skew` is sorted by ν, which is the whole point of the sort.
        expansion
            .binary_search_by(|(p, _)| p.cmp(nu))
            .map(|i| expansion[i].1)
            .unwrap_or(0)
    }

    /// A product is a skew expansion in disguise.
    ///
    /// Place μ up and to the right of ν so the two diagrams share no row and no
    /// column; the result is a skew shape whose fillings are exactly a filling
    /// of μ alongside one of ν, hence s_{shape} = s_μ · s_ν. Expanding that one
    /// shape yields every λ in the product at once — no candidate sweep, and no
    /// per-λ call to `lr_coeff`.
    fn schur_product(&self, mu: &Partition, nu: &Partition) -> Vec<(Partition, u128)> {
        // c^λ_{μν} is symmetric, so fix an orientation: s_μ·s_ν and s_ν·s_μ
        // then land on the same shape and share one cache entry.
        let (a, b) = if mu >= nu { (mu, nu) } else { (nu, mu) };
        let (outer, inner) = juxtapose(a, b);
        expand_skew(&outer, &inner)
    }
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
            assert_eq!(
                expand_skew(&l, &Partition::default()),
                vec![(l.clone(), 1)]
            );
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
    /// The run fill decides a maximal block of equal values at once, so what can
    /// go wrong is a run that stops one column early or late. These shapes put a
    /// run against each boundary in turn: an *overhang* (columns past the end of
    /// the row above, where nothing blocks and a run may spill to the row end),
    /// a row above that blocks in the middle (`cut` interior), a row above that
    /// blocks at its very first column, and rows that share no column at all
    /// (empty overlap, so the frontier's `above` half is empty).
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
    /// the frontier through it rather than dropping or duplicating states.
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

    /// The packed frontier key must distinguish states the old two-`Vec` key
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
            let want = product
                .iter()
                .find(|(l, _)| l == lambda)
                .map_or(0, |(_, c)| *c);
            assert_eq!(*was, want, "vs product term for {lambda}");
        }
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

    /// The three serializations must be interchangeable: any width wide enough
    /// for the shape yields the same expansion. Forcing 2 and 4 bytes onto
    /// one-byte shapes exercises every encode/decode pair on frontiers with
    /// real merging, where a mis-split key would corrupt coefficients.
    #[test]
    fn widths_agree_on_merging_shapes() {
        // The traversal reports terms in map order, which legitimately varies
        // with the serialization; sort before comparing.
        let sorted = |v: Option<Vec<(Partition, u128)>>| {
            let mut v = v.expect("no overflow at u64");
            v.sort_by(|a, b| a.0.cmp(&b.0));
            v
        };
        for (o, i) in [
            (&[10, 8, 5, 1][..], &[6, 3, 1][..]),
            (&[6, 5, 4, 3], &[2, 1]),
            (&[7, 6, 4, 2], &[3, 2, 1]),
        ] {
            let (outer, inner) = (p(o), p(i));
            let narrow = sorted(expand_with_width::<u64>(&outer, &inner, 1));
            let mid = sorted(expand_with_width::<u64>(&outer, &inner, 2));
            let wide = sorted(expand_with_width::<u64>(&outer, &inner, 4));
            assert_eq!(narrow, mid, "widths 1 vs 2 on {outer}/{inner}");
            assert_eq!(narrow, wide, "widths 1 vs 4 on {outer}/{inner}");
            assert_eq!(narrow, reference_skew(&outer, &inner), "vs naive");
        }
    }

    /// Keys longer than the inline capacity spill to the heap; both paths must
    /// coexist and merge correctly.
    ///
    /// `[40, 36]/∅` carries a 36-wide clipped row (39-byte keys, all heap) and
    /// has exactly one filling, so the answer is pinned: s_{λ/∅} = s_λ.
    /// `[33, 31]²` mixes inline and heap keys in one frontier *with* merging;
    /// its expansion is checked against the conjugate orientation, which is an
    /// independent traversal (2-wide keys, all inline) of the same
    /// coefficients via c^λ_{μν} = c^{λ'}_{μ'ν'}.
    #[test]
    fn heap_keys_merge_and_agree_with_conjugate_orientation() {
        let wide = p(&[40, 36]);
        assert_eq!(
            expand_skew(&wide, &Partition::default()),
            vec![(wide.clone(), 1)]
        );

        let m = p(&[33, 31]);
        let direct = SkewLr.schur_product(&m, &m);
        let mc = m.conjugate();
        let mut via_conjugate: Vec<(Partition, u128)> = SkewLr
            .schur_product(&mc, &mc)
            .into_iter()
            .map(|(l, c)| (l.conjugate(), c))
            .collect();
        via_conjugate.sort_by(|a, b| a.0.cmp(&b.0));
        assert_eq!(direct, via_conjugate);
    }

    /// Frontier multiplicities are tableau counts, so no narrow accumulator is
    /// provably safe — the engine must *detect* saturation and retry wider.
    /// A `u8` accumulator makes the boundary cheap to reach: the pass must
    /// report overflow (not wrap), and the two production widths must agree.
    #[test]
    fn accumulator_overflow_is_detected_not_wrapped() {
        impl Acc for u8 {
            const ONE: Self = 1;
            fn checked_add(self, other: Self) -> Option<Self> {
                u8::checked_add(self, other)
            }
            fn widen(self) -> u128 {
                self as u128
            }
        }

        let sorted = |v: Option<Vec<(Partition, u128)>>| {
            v.map(|mut v| {
                v.sort_by(|a, b| a.0.cmp(&b.0));
                v
            })
        };

        // s[6,5,4,3,2]² has a coefficient of 644, and a final state's
        // multiplicity *is* its coefficient, so the u8 pass must saturate…
        let (outer, inner) = juxtapose(&p(&[6, 5, 4, 3, 2]), &p(&[6, 5, 4, 3, 2]));
        assert_eq!(expand_oriented::<u8>(&outer, &inner), None);
        let full = sorted(expand_oriented::<u64>(&outer, &inner));
        assert!(full.is_some(), "u64 must not saturate here");
        assert_eq!(full, sorted(expand_oriented::<u128>(&outer, &inner)));

        // …while a tiny shape stays under 255, so u8 must also *succeed* when
        // nothing overflows (the detection is not a blanket refusal).
        let (outer, inner) = juxtapose(&p(&[2, 1]), &p(&[2, 1]));
        let small = sorted(expand_oriented::<u8>(&outer, &inner));
        assert!(small.is_some());
        assert_eq!(small, sorted(expand_oriented::<u128>(&outer, &inner)));
    }
}
