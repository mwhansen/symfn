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
/// The state is packed into a single `u32` buffer, `[len, content.., above..]`:
/// `content[v - 1]` is the number of `v`s placed so far (no trailing zeros, so
/// the leading `len` makes the split unambiguous), and `above` holds the
/// previous row's entries, clipped to the columns the next row overlaps. One
/// buffer means one allocation per *new* state and a plain slice compare on a
/// hit — and hits can be probed with a borrowed scratch buffer, so the hot path
/// allocates nothing at all (see [`commit`]).
#[derive(PartialEq, Eq)]
struct Key(Box<[u32]>);

/// Borrowed view of a [`Key`], so a frontier probe can use a scratch slice.
///
/// `#[repr(transparent)]` makes the `&[u32]` → `&KeySlice` cast sound; `Hash`
/// and `Eq` agree with [`Key`]'s exactly, which is what `Borrow` requires.
#[derive(PartialEq, Eq)]
#[repr(transparent)]
struct KeySlice([u32]);

impl KeySlice {
    #[inline]
    fn new(s: &[u32]) -> &KeySlice {
        // SAFETY: KeySlice is repr(transparent) over [u32].
        unsafe { &*(s as *const [u32] as *const KeySlice) }
    }
}

impl Borrow<KeySlice> for Key {
    #[inline]
    fn borrow(&self) -> &KeySlice {
        KeySlice::new(&self.0)
    }
}

impl Hash for KeySlice {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        write_words(state, &self.0);
    }
}

impl Hash for Key {
    #[inline]
    fn hash<H: Hasher>(&self, state: &mut H) {
        write_words(state, &self.0);
    }
}

/// Feed a `u32` slice to a hasher two elements at a time.
///
/// Frontier keys are hashed by the million and are only a handful of words
/// long, so the per-element trait round trip is a real fraction of the cost.
#[inline]
fn write_words<H: Hasher>(state: &mut H, words: &[u32]) {
    let mut it = words.chunks_exact(2);
    for c in &mut it {
        state.write_u64((c[0] as u64) << 32 | c[1] as u64);
    }
    if let [last] = it.remainder() {
        state.write_u64(*last as u64);
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
        self.0
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

fn expand_skew_uncached(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
    let mut out = expand_oriented(outer, inner);
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

fn expand_oriented(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
    let rows = outer.len();

    // The frontier of the traversal: reduced state -> number of ways to reach it.
    let mut states: Map<Key, u128> = Map::default();
    states.insert(Key(vec![0u32].into_boxed_slice()), 1);

    // Per-row scratch, reused across states so the inner loop never allocates.
    let mut row: Vec<u32> = Vec::new();
    let mut added: Vec<u32> = Vec::new();
    let mut cut: Vec<usize> = Vec::new();
    let mut gap: Vec<u32> = Vec::new();
    let mut key: Vec<u32> = Vec::new();

    for r in 0..rows {
        let lo = inner.part(r) as usize;
        let hi = outer.part(r) as usize;
        // Columns of row r-1 that row r sits under, and the ones of row r that
        // row r+1 will sit under. Clipping to these is what makes distinct
        // histories collapse.
        let (up_lo, up_hi) = overlap(inner, outer, r.wrapping_sub(1), r);
        let (dn_lo, dn_hi) = overlap(inner, outer, r, r + 1);

        let mut next: Map<Key, u128> =
            Map::with_capacity_and_hasher(states.len() * 2, Default::default());
        row.clear();
        row.resize(hi.saturating_sub(lo), 0);
        for (state, mult) in states.iter() {
            let clen = state.0[0] as usize;
            let content = &state.0[1..1 + clen];
            let above = &state.0[1 + clen..];
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
                mult: *mult,
                out: &mut next,
            };
            fill_runs(lo, 1, &mut ctx);
        }
        // Both frontiers are momentarily live here; record the sum (once per
        // row, so the cost is nil).
        PEAK_LIVE_STATES.fetch_max(states.len() + next.len(), Ordering::Relaxed);
        states = next;
        if states.is_empty() {
            break;
        }
    }

    // Every row's down-overlap with a nonexistent next row is empty, so the
    // final keys are pure content — already merged, one term each.
    states
        .into_iter()
        .map(|(key, c)| {
            debug_assert_eq!(key.0[0] as usize + 1, key.0.len());
            // The ballot condition forces the content to be weakly decreasing
            // at every prefix, so the final one is already a partition.
            (Partition::from_sorted(key.0[1..].to_vec()), c)
        })
        .collect()
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
struct RowCtx<'a> {
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
    mult: u128,
    out: &'a mut Map<Key, u128>,
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
fn fill_runs(a: usize, vmin: u32, ctx: &mut RowCtx) {
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
fn finish_row(ctx: &mut RowCtx) {
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

    // Probe with the borrowed scratch; only a genuinely new state allocates.
    match ctx.out.get_mut(KeySlice::new(key)) {
        Some(w) => *w += ctx.mult,
        None => {
            ctx.out.insert(Key(key.as_slice().into()), ctx.mult);
        }
    }
}

impl LrBackend for SkewLr {
    fn lr_coeff(&self, lambda: &Partition, mu: &Partition, nu: &Partition) -> u128 {
        if lambda.size() != mu.size() + nu.size() || !lambda.contains(mu) {
            return 0;
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
}
