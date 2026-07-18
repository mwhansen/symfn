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

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hash, Hasher};

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
/// Two partial fillings that agree on this pair are interchangeable: everything
/// still to be decided depends on the row above (column strictness) and on the
/// content accumulated so far (the ballot condition), and on nothing else. So
/// they are merged and carry a multiplicity, which is what keeps the traversal
/// from degenerating into an enumeration of individual tableaux.
///
/// `content[v - 1]` is the number of `v`s placed so far, trailing zeros
/// trimmed; `above` holds the previous row's entries, clipped to the columns
/// the next row actually overlaps.
#[derive(Clone, PartialEq, Eq)]
struct Frontier {
    content: Vec<u32>,
    above: Vec<u32>,
}

impl Hash for Frontier {
    fn hash<H: Hasher>(&self, state: &mut H) {
        state.write_u32(self.content.len() as u32);
        write_words(state, &self.content);
        write_words(state, &self.above);
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

fn expand_skew_uncached(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
    let rows = outer.len();

    // The frontier of the traversal: reduced state -> number of ways to reach it.
    let mut states: Map<Frontier, u128> = Map::default();
    states.insert(
        Frontier {
            content: Vec::new(),
            above: Vec::new(),
        },
        1,
    );

    for r in 0..rows {
        let lo = inner.part(r) as usize;
        let hi = outer.part(r) as usize;
        // Columns of row r-1 that row r sits under, and the ones of row r that
        // row r+1 will sit under. Clipping to these is what makes distinct
        // histories collapse.
        let (up_lo, up_hi) = overlap(inner, outer, r.wrapping_sub(1), r);
        let (dn_lo, dn_hi) = overlap(inner, outer, r, r + 1);

        let mut next: Map<Frontier, u128> =
            Map::with_capacity_and_hasher(states.len() * 2, Default::default());
        let mut row = vec![0u32; hi.saturating_sub(lo)];
        let mut added: Vec<u32> = Vec::new();
        for (state, mult) in states.iter() {
            added.clear();
            added.resize(state.content.len() + 1, 0);
            let mut ctx = RowCtx {
                lo,
                hi,
                up_lo,
                up_hi,
                dn_lo,
                dn_hi,
                content: &state.content,
                above: &state.above,
                added: &mut added,
                row: &mut row,
                mult: *mult,
                out: &mut next,
            };
            fill_row(0, 1, &mut ctx);
        }
        states = next;
        if states.is_empty() {
            break;
        }
    }

    let mut totals: HashMap<Vec<u32>, u128> = HashMap::new();
    for (state, mult) in states {
        *totals.entry(state.content).or_insert(0) += mult;
    }
    let mut out: Vec<(Partition, u128)> = totals
        .into_iter()
        .filter(|(_, c)| *c != 0)
        // The ballot condition forces the content to be weakly decreasing at
        // every prefix, so in particular the final one is already a partition.
        .map(|(content, c)| (Partition::from_sorted(content), c))
        .collect();
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
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
    up_lo: usize,
    up_hi: usize,
    dn_lo: usize,
    dn_hi: usize,
    /// Content *before* this row; fixed while the row is filled.
    content: &'a [u32],
    /// The previous row's entries, indexed from column `up_lo`.
    above: &'a [u32],
    /// Occurrences of each value contributed by this row so far.
    added: &'a mut [u32],
    row: &'a mut [u32],
    mult: u128,
    out: &'a mut Map<Frontier, u128>,
}

/// Fill columns `lo + i ..` of the current row, left to right.
///
/// Left-to-right is the order in which all three constraints are already
/// decided: the row entry to the left bounds this one below (weak increase),
/// and the entry above bounds it below strictly.
///
/// The ballot condition also collapses to a per-row test, which is what lets a
/// row be treated as one step. The reading word takes a row right to left, so
/// within a row it is weakly *decreasing*; a prefix of the word that stops
/// partway through the row has therefore taken all entries greater than some
/// value `t` and some of the `t`s. For a pair `(i, i+1)` with `i+1 > t` both
/// counts are already final for the row, and for `i+1 < t` neither has moved.
/// Only `i+1 = t` is a new constraint, and it is worst when every `t` in the
/// row has been read while none of the `t-1`s have. So over the whole row the
/// condition is exactly
///
/// ```text
///   for every value v ≥ 2 :  before[v-1]  ≥  before[v] + (number of v in row)
/// ```
///
/// where `before` is the content accumulated by earlier rows. This implies the
/// `i+1 > t` cases, so it is the only check needed — and it is a statement about
/// counts alone, hence checkable as each cell is committed.
fn fill_row(i: usize, prev: u32, ctx: &mut RowCtx) {
    let width = ctx.hi.saturating_sub(ctx.lo);
    if i == width {
        finish_row(ctx);
        return;
    }
    let col = ctx.lo + i;
    // Column strictness: beat the cell above, if this column has one.
    let floor = if col >= ctx.up_lo && col < ctx.up_hi {
        ctx.row_above(col) + 1
    } else {
        1
    };
    let start = prev.max(floor);
    // A value v may only be used while some v-1 is already banked, so the
    // reachable alphabet is exactly one wider than the content so far. (This
    // also recovers the familiar bound "entries in row r are at most r+1".)
    let top = ctx.content.len() as u32 + 1;

    for v in start..=top {
        let vi = (v - 1) as usize;
        if v >= 2 {
            // Ballot: after this row, #(v-1) must still dominate #v. The v-1s
            // available are only those banked *before* this row, because the
            // reading word takes this row's cells right to left — every v in
            // the row is read before every v-1 in it.
            let have = ctx.content[vi - 1];
            let want = ctx.content.get(vi).copied().unwrap_or(0) + ctx.added[vi] + 1;
            if have < want {
                // Only this value is blocked: the condition is about the gap
                // content[v-1] − content[v], which is not monotone in v, so a
                // larger value may still be admissible.
                continue;
            }
        }
        ctx.row[i] = v;
        ctx.added[vi] += 1;
        fill_row(i + 1, v, ctx);
        ctx.added[vi] -= 1;
    }
}

/// Commit a completed row: fold its content in and clip the frontier.
fn finish_row(ctx: &mut RowCtx) {
    let mut content = ctx.content.to_vec();
    for (vi, &n) in ctx.added.iter().enumerate() {
        if n == 0 {
            continue;
        }
        if vi >= content.len() {
            content.resize(vi + 1, 0);
        }
        content[vi] += n;
    }
    while content.last() == Some(&0) {
        content.pop();
    }

    let above = if ctx.dn_lo < ctx.dn_hi {
        ctx.row[ctx.dn_lo - ctx.lo..ctx.dn_hi - ctx.lo].to_vec()
    } else {
        Vec::new()
    };

    *ctx.out.entry(Frontier { content, above }).or_insert(0) += ctx.mult;
}

impl RowCtx<'_> {
    /// The entry of the previous row in column `col`.
    #[inline]
    fn row_above(&self, col: usize) -> u32 {
        // `above` was clipped to [up_lo, up_hi) when the previous row finished.
        debug_assert!(col >= self.up_lo && col < self.up_hi);
        self.above[col - self.up_lo]
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
