//! Kostka numbers K_{λμ} = number of semistandard Young tableaux of shape λ and
//! content μ (μ_i entries equal to i, rows weakly increasing, columns strictly
//! increasing).
//!
//! These are the s → m transition: s_λ = Σ_μ K_{λμ} m_μ. The reverse, m → s,
//! inverts the (unitriangular) Kostka matrix; see `convert`.
//!
//! ## How it is computed
//!
//! Not by enumerating tableaux. Grouping an SSYT's cells by value gives the
//! standard bijection with **chains of horizontal strips**
//!
//! ```text
//!   ∅ = λ⁰ ⊆ λ¹ ⊆ … ⊆ λ^ℓ(μ) = λ,   λⁱ/λⁱ⁻¹ a horizontal strip of size μᵢ
//! ```
//!
//! (the cells holding i form a horizontal strip, since two i's in one column
//! would break column-strictness). So K_{λμ} counts chains, and chains through
//! the same intermediate shape **merge** — which is the whole difference
//! between counting and enumerating. This is `strip_lr`'s dynamic program
//! without the lattice condition, which is exactly what makes it count SSYT
//! rather than LR tableaux.
//!
//! Every intermediate shape is pruned to λ, so the state space is bounded by the
//! partitions inside λ rather than by the tableaux of shape λ.
//!
//! Enumerating SSYT one cell at a time instead is exponential, and at degree 20
//! a single K_{λμ} that way costs more than this DP costs for a whole
//! conversion (`docs/record/transitions.md`).

// Kostka numbers are counts of tableaux, so the `i128 → u128` on the way out is
// non-negative by definition; the rest are shape indices.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::HashMap;

use crate::coeff::Ring;
use crate::convert::{beta_mask, pieri_trie, unit_mask_layer, MASK_LIMIT};
use crate::fasthash::Map;
use crate::memo::kostka_cached;
use crate::partition::Partition;

/// The Kostka number K_{λμ}, the number of semistandard Young tableaux of
/// shape λ and content μ.
///
/// Requires nothing of the arguments beyond being partitions. Returns 0 unless
/// |λ| = |μ| and λ ⊵ μ in dominance order, which together are exactly the
/// nonzero cases. Returns 1 when both are empty, counting the empty tableau.
///
/// Sage computes it as `SemistandardTableaux(λ, μ).cardinality()`, which
/// `scripts/compare_sage.py` drives. Symmetrica's entry point is
/// `kostka_number` (`scripts/compare_symmetrica.py`).
///
/// **Range.** The largest value at degree n is K_{λ,1ⁿ} = f^λ ≈ √(n!), which
/// passes `u128` near n ≈ 58. A whole degree walls earlier on memory; see
/// [`kostka_table`].
pub fn kostka(lambda: &Partition, mu: &Partition) -> u128 {
    if lambda.size() != mu.size() {
        return 0;
    }
    if lambda.is_empty() {
        return 1; // both empty: the empty tableau
    }
    // K_{λμ} ≠ 0 iff λ dominates μ, so this is an exact O(rows) early exit.
    // Worth only ~1.1x in practice, not the large win the sparsity of dominance
    // suggests: the chain DP's capacity pruning already rejected these quickly,
    // so the test mostly replaces a fast zero with a faster one. Kept because it
    // is free and states the fact outright, not because it is a real speedup.
    if !dominates(lambda, mu) {
        return 0;
    }
    kostka_cached(lambda, mu, || kostka_uncached(lambda, mu))
}

/// Whether λ ⊵ μ in dominance order: every prefix sum of λ is at least μ's.
pub(crate) fn dominates(lambda: &Partition, mu: &Partition) -> bool {
    let (mut a, mut b) = (0u32, 0u32);
    for i in 0..lambda.len().max(mu.len()) {
        a += lambda.part(i);
        b += mu.part(i);
        if a < b {
            return false;
        }
    }
    true
}

fn kostka_uncached(lambda: &Partition, mu: &Partition) -> u128 {
    if mu.is_empty() {
        return 0;
    }
    let bound = lambda.parts();
    // Layer of the chain: intermediate shape -> number of ways to reach it.
    let mut cur: HashMap<Vec<u32>, u128> = HashMap::new();
    cur.insert(Vec::new(), 1);
    let mut next: HashMap<Vec<u32>, u128> = HashMap::new();
    let mut buf: Vec<u32> = Vec::new();

    for &r in mu.parts() {
        next.clear();
        for (shape, &ways) in cur.iter() {
            buf.clear();
            buf.extend_from_slice(shape);
            buf.resize(bound.len(), 0);
            grow(0, r, u32::MAX, &mut buf, bound, &mut |grown: &[u32]| {
                let end = grown.iter().rposition(|&x| x > 0).map_or(0, |i| i + 1);
                *next.entry(grown[..end].to_vec()).or_insert(0) += ways;
            });
        }
        std::mem::swap(&mut cur, &mut next);
        if cur.is_empty() {
            return 0;
        }
    }
    cur.get(bound).copied().unwrap_or(0)
}

/// Add a horizontal strip of `left` cells to `shape` in place, staying inside
/// `bound`, calling `emit` on each result.
///
/// A horizontal strip means the interlacing `shape_{i-1} ≥ new_i ≥ shape_i`, so
/// row `i` may grow only up to the *previous row's old* value — and never past
/// `bound_i`, which is what keeps the layer proportional to the partitions
/// inside λ rather than to all partitions.
fn grow(
    i: usize,
    left: u32,
    prev_old: u32,
    shape: &mut Vec<u32>,
    bound: &[u32],
    emit: &mut impl FnMut(&[u32]),
) {
    if i == shape.len() {
        if left == 0 {
            emit(shape);
        }
        return;
    }
    // What all remaining rows can still absorb; bail when `left` cannot fit.
    let capacity: u32 = bound[i..]
        .iter()
        .zip(shape[i..].iter())
        .map(|(b, s)| b.saturating_sub(*s))
        .sum();
    if capacity < left {
        return;
    }
    let old = shape[i];
    let hi = prev_old.min(bound[i]).min(old + left);
    if hi < old {
        return;
    }
    for v in old..=hi {
        shape[i] = v;
        grow(i + 1, left - (v - old), old, shape, bound, emit);
    }
    shape[i] = old;
}

/// The semistandard Young tableaux of shape `lambda` and weight `weight`, each
/// as its list of rows, top row first.
///
/// The list [`kostka`] counts. Everything the count gets to merge, this has to
/// keep apart, so the cost is `K_{λμ}` rather than the number of shapes inside
/// λ. The two functions are the same walk read at different resolutions. Asking
/// for the tableaux when the count will do is the expensive mistake.
///
/// The empty vector when `|λ| ≠ |weight|` or `K_{λμ} = 0`; a single empty
/// tableau `[[]]` when both are empty.
///
/// **The weight is a composition, not a partition**, and the difference is not
/// cosmetic: entry `i` of `weight` is how many `i + 1`s the tableau carries, so
/// a zero is a value that goes unused and shifts every later label. Weight
/// `(2, 0, 1)` gives `1 1 3` where `(2, 1)` gives `1 1 2` — the same count,
/// different tableaux. [`kostka`] may sort its weight because `K_{λμ}` is
/// symmetric in μ; this may not. Sage reaches exactly this case, since
/// `SemistandardTableaux(λ)` iterates over every content vector of `|λ|` and
/// most of them are not weakly decreasing.
///
/// **Order**: increasing lexicographic in the row-major reading word — row 1
/// left to right, then row 2, and so on. This is Symmetrica's `kostka_tab`
/// order, which Sage's `SemistandardTableaux(λ, μ)` doctests print verbatim, so
/// it is part of the contract rather than an artifact of the traversal. The
/// chain walk below does not produce it, so the result is sorted; `kostka_tab`
/// in `scripts/check_backend.py` checks the order against Symmetrica itself
/// over every pair up to degree 9.
///
/// # Examples
///
/// The two SSYT of shape `(3,1)` and weight `(2,1,1)`, in that order — the 2 is
/// placed in the first row before it is placed in the second.
///
/// ```
/// use symfn::{semistandard_tableaux, Partition};
/// let ts = semistandard_tableaux(&Partition::new([3, 1]), &[2, 1, 1]);
/// assert_eq!(ts, vec![vec![vec![1, 1, 2], vec![3]], vec![vec![1, 1, 3], vec![2]]]);
/// ```
///
/// The same shape at weight `(2, 0, 1, 1)` — as many tableaux, relabeled,
/// because the value `2` is now unused.
///
/// ```
/// use symfn::{semistandard_tableaux, Partition};
/// let ts = semistandard_tableaux(&Partition::new([3, 1]), &[2, 0, 1, 1]);
/// assert_eq!(ts, vec![vec![vec![1, 1, 3], vec![4]], vec![vec![1, 1, 4], vec![3]]]);
/// ```
pub fn semistandard_tableaux(lambda: &Partition, weight: &[u32]) -> Vec<Vec<Vec<u32>>> {
    // Widened, because a caller's weight is arbitrary data and its parts need
    // not fit `u32` in sum the way a partition's do.
    if u64::from(lambda.size()) != weight.iter().map(|&x| u64::from(x)).sum::<u64>() {
        return Vec::new();
    }
    if lambda.is_empty() {
        return vec![Vec::new()];
    }
    // The chain of intermediate shapes, padded to ℓ(λ) so a row that is still
    // empty is a 0 rather than a missing entry.
    let mut chain: Vec<Vec<u32>> = vec![vec![0u32; lambda.len()]];
    let mut out = Vec::new();
    chains(lambda, weight, &mut chain, &mut out);
    // The walk groups by value — all chains sharing where the 1s went come out
    // together — and the reading word orders by position instead, so the two
    // disagree from the first shape with three rows onwards. Sorting is a log
    // factor on an enumeration that already costs `K_{λμ}·|λ|`.
    out.sort_by_key(|a| reading_word(a));
    out
}

/// The row-major reading word: row 1 left to right, then row 2, and so on.
fn reading_word(t: &[Vec<u32>]) -> Vec<u32> {
    t.concat()
}

/// Extend `chain` by one horizontal strip per remaining part of μ, recording a
/// tableau at every completed chain.
fn chains(
    lambda: &Partition,
    left: &[u32],
    chain: &mut Vec<Vec<u32>>,
    out: &mut Vec<Vec<Vec<u32>>>,
) {
    let Some((&r, rest)) = left.split_first() else {
        out.push(tableau_of(chain));
        return;
    };
    let mut buf = chain
        .last()
        .expect("the chain starts at the empty shape")
        .clone();
    grow(0, r, u32::MAX, &mut buf, lambda.parts(), &mut |grown| {
        // Prune what the remaining strips can no longer carry up to λ: with k
        // strips left that needs λ_{j+k} ≤ ν_j, the k-strip form of "s_λ
        // vanishes in fewer than ℓ(λ) variables".
        let k = rest.len();
        for j in 0..lambda.len() {
            let need = if j + k < lambda.len() {
                lambda.part(j + k)
            } else {
                0
            };
            if grown[j] < need {
                return;
            }
        }
        chain.push(grown.to_vec());
        chains(lambda, rest, chain, out);
        chain.pop();
    });
}

/// The tableau of a chain: the cells of `λⁱ/λⁱ⁻¹` hold `i`.
///
/// Reading a row left to right meets the values in increasing order of `i`, so
/// rows come out weakly increasing for free; column strictness is what the
/// horizontal-strip condition on each step already guaranteed.
fn tableau_of(chain: &[Vec<u32>]) -> Vec<Vec<u32>> {
    let rows = chain.last().map_or(0, Vec::len);
    let mut out = vec![Vec::new(); rows];
    for (step, pair) in chain.windows(2).enumerate() {
        let value = step as u32 + 1;
        for row in 0..rows {
            for _ in pair[0][row]..pair[1][row] {
                out[row].push(value);
            }
        }
    }
    out
}

/// The whole Kostka table of degree `n`, as `table[i][j] = K_{λⁱ λʲ}`, with
/// rows and columns both indexed by
/// `memo::partitions_cached`.
///
/// Computes `p(n)²` values in `p(n)` chain sweeps. Dropping [`kostka`]'s
/// bound on λ leaves the final layer of μ's chain holding the whole column.
/// K_{λμ} depends on μ only as a multiset, so consuming parts in descending
/// order lets every μ with a common prefix share that initial segment.
///
/// Symmetrica's entry point is `kostka_tafel`, which
/// `scripts/compare_symmetrica.py` drives; [`kostka`] names the single-pair
/// equivalents.
///
/// **Range.** Entries pass `u128` near n ≈ 58 — the largest is
/// K_{λ,1ⁿ} = f^λ ≈ √(n!) — but the table is `p(n)²` values, 1.1 GB at
/// n = 32, so memory walls about twenty degrees earlier.
/// [`kostka_table_in`] takes an arbitrary ring and does not move that.
pub fn kostka_table(n: u32) -> Vec<Vec<u128>> {
    // `i128` internally, then cast: Kostka numbers are non-negative, and the
    // half-bit given up sits past the memory wall in the range above.
    kostka_table_in::<i128>(n)
        .into_iter()
        .map(|row| row.into_iter().map(|v| v as u128).collect())
        .collect()
}

/// [`kostka_table`] over an arbitrary coefficient ring.
///
/// Not a widening: the entries are the same counts, and the memory wall on
/// [`kostka_table`] arrives first whatever `C` is. It exists so the sweep can
/// carry a different accumulator — Kostka–Foulkes refines K_{λμ} by charge, so
/// K_{λμ}(t) is this same walk over horizontal strips carrying a polynomial
/// rather than a count.
pub fn kostka_table_in<C: Ring>(n: u32) -> Vec<Vec<C>> {
    let parts = crate::memo::partitions_cached(n);
    let mut table = vec![vec![C::zero(); parts.len()]; parts.len()];
    if n == 0 {
        table[0][0] = C::one();
        return table;
    }
    let l = n as usize;
    if l <= MASK_LIMIT {
        // The h → s trie on β-masks, whose leaf for μ is h_μ in the Schur
        // basis — mask ↦ K_{λμ} — written down as μ's column. Layer counts are
        // `i128` for the reason `convert::expand_multiplicative` gives: K_{λμ}
        // ≤ f^λ ≤ √(n!), inside `i128` through the mask width, so `C` is
        // touched once per entry. `partitions_cached` is in descending
        // lexicographic order, which keeps common prefixes contiguous. 1.5x
        // over the `Vec`-keyed sweep below at n = 20 and 1.3x at 24
        // (`bench_ops`, `kostka_table_n20/24`; `docs/record/transitions.md`).
        let index: Map<u64, usize> = parts
            .iter()
            .enumerate()
            .map(|(i, p)| (beta_mask(p, l), i))
            .collect();
        let leaves: Vec<(&Partition, usize)> = parts.iter().zip(0..).collect();
        pieri_trie(
            &leaves,
            0,
            &unit_mask_layer(l),
            false,
            &mut |_, &col, layer| {
                crate::interrupt::poll();
                for (mask, &v) in layer {
                    table[index[mask]][col] = C::from_i128(v);
                }
            },
        );
        return table;
    }
    table_on_partitions(n, &parts, table)
}

/// [`kostka_table_in`] on partition-keyed layers: the same trie, with no wall
/// before the table's own memory wall. What runs past the β-mask width.
fn table_on_partitions<C: Ring>(
    n: u32,
    parts: &[Partition],
    mut table: Vec<Vec<C>>,
) -> Vec<Vec<C>> {
    let index: HashMap<&[u32], usize> = parts
        .iter()
        .enumerate()
        .map(|(i, p)| (p.parts(), i))
        .collect();
    // Descending part order, so the longest common prefixes are shared.
    let mut order: Vec<usize> = (0..parts.len()).collect();
    order.sort_by(|&a, &b| parts[a].parts().cmp(parts[b].parts()));
    let mut root: HashMap<Vec<u32>, C> = HashMap::new();
    root.insert(Vec::new(), C::one());
    table_sweep(n, parts, &order, 0, &root, &index, &mut table);
    table
}

fn table_sweep<C: Ring>(
    n: u32,
    parts: &[Partition],
    group: &[usize],
    depth: usize,
    layer: &HashMap<Vec<u32>, C>,
    index: &HashMap<&[u32], usize>,
    table: &mut [Vec<C>],
) {
    let mut i = 0;
    // Columns that end here: the layer is exactly this μ's column.
    while i < group.len() && parts[group[i]].len() == depth {
        crate::interrupt::poll();
        let col = group[i];
        for (shape, ways) in layer {
            if let Some(&row) = index.get(shape.as_slice()) {
                table[row][col] = ways.clone();
            }
        }
        i += 1;
    }
    // The rest are grouped by their next part, each group sharing one step.
    while i < group.len() {
        crate::interrupt::poll();
        let r = parts[group[i]].part(depth);
        let start = i;
        while i < group.len() && parts[group[i]].part(depth) == r {
            i += 1;
        }
        // After `depth` strips a shape has at most `depth` rows, and one strip
        // adds at most one new row: a strip needs shape_{i-1} ≥ new_i, so the
        // first empty row can grow but the next is pinned to 0.
        let bound = vec![n; depth + 1];
        let mut next: HashMap<Vec<u32>, C> = HashMap::new();
        let mut buf: Vec<u32> = Vec::new();
        for (shape, ways) in layer {
            buf.clear();
            buf.extend_from_slice(shape);
            buf.resize(bound.len(), 0);
            grow(0, r, u32::MAX, &mut buf, &bound, &mut |grown: &[u32]| {
                let end = grown.iter().rposition(|&x| x > 0).map_or(0, |i| i + 1);
                next.entry(grown[..end].to_vec())
                    .or_insert_with(C::zero)
                    .add_assign(ways);
            });
        }
        table_sweep(n, parts, &group[start..i], depth + 1, &next, index, table);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    /// The swept table must equal the per-pair `kostka` at every entry.
    ///
    /// Two things could silently go wrong and neither shows up in a spot check:
    /// the trie can misattribute a column if the prefix grouping is off by one,
    /// and dropping the λ-bound changes which shapes the layer reaches, so a
    /// row could go missing rather than wrong.
    #[test]
    fn swept_table_matches_per_pair_kostka() {
        for n in 0..=12u32 {
            let parts = crate::memo::partitions_cached(n);
            let table = kostka_table(n);
            assert_eq!(table.len(), parts.len(), "table is p({n}) square");
            for (i, lambda) in parts.iter().enumerate() {
                for (j, mu) in parts.iter().enumerate() {
                    assert_eq!(
                        table[i][j],
                        kostka(lambda, mu),
                        "K_{{{lambda},{mu}}} at degree {n}"
                    );
                }
            }
        }
    }

    #[test]
    fn small_kostka_values() {
        assert_eq!(kostka(&p(&[2]), &p(&[2])), 1);
        assert_eq!(kostka(&p(&[2]), &p(&[1, 1])), 1); // tableau "1 2"
        assert_eq!(kostka(&p(&[1, 1]), &p(&[2])), 0); // column can't repeat
        assert_eq!(kostka(&p(&[1, 1]), &p(&[1, 1])), 1);
        assert_eq!(kostka(&p(&[2, 1]), &p(&[1, 1, 1])), 2); // standard tableaux of (2,1)
        assert_eq!(kostka(&p(&[3, 2, 1]), &p(&[1, 1, 1, 1, 1, 1])), 16); // f^{321}
    }

    /// The old implementation stored tableau entries in `u8` and wrapped its
    /// loop bound, so `ℓ(μ) = 256` silently returned 0 instead of the true
    /// count. The chain DP stores no entries at all, so that boundary cannot
    /// exist — and the case is now cheap enough to assert through the public
    /// API, where before it was exponential and needed a white-box test.
    ///
    /// K_{1ⁿ,1ⁿ} = 1: a single column admits exactly one SSYT.
    #[test]
    fn value_range_extends_past_the_u8_boundary() {
        for len in [200usize, 255, 256, 300] {
            let tall = Partition::new(std::iter::repeat_n(1, len));
            assert_eq!(kostka(&tall, &tall), 1, "K_{{1^{len},1^{len}}}");
        }
    }

    /// `K_{λ,1ⁿ}` counts standard tableaux, so the hook-length formula
    /// `n!/∏h(i,j)` gives it independently — sharing no code with the chain DP.
    ///
    /// This reaches degrees the Sage fixture does not, which is the point: a
    /// correctness check covering only what the exponential predecessor could
    /// afford would never exercise the DP where it now operates.
    #[test]
    fn standard_tableaux_count_matches_the_hook_formula() {
        for parts in [
            &[3u32, 2, 1][..],
            &[4, 3, 2, 1],
            &[5, 4, 3, 2, 1],
            &[6, 5, 4, 3, 2],
            &[7, 6, 5, 4, 3],
        ] {
            let lam = p(parts);
            let n = lam.size();
            let ones = Partition::new(std::iter::repeat_n(1, n as usize));

            // n! / ∏ hooks, both exact in u128 at these degrees.
            let conj = lam.conjugate();
            let mut hooks: u128 = 1;
            for (i, &row) in lam.parts().iter().enumerate() {
                for j in 0..row as usize {
                    hooks *= (row as u128 - j as u128) + (conj.part(j) as u128 - i as u128) - 1;
                }
            }
            let factorial: u128 = (1..=u128::from(n)).product();
            assert_eq!(
                factorial % hooks,
                0,
                "hook product must divide n! for {lam}"
            );
            assert_eq!(kostka(&lam, &ones), factorial / hooks, "K_{{{lam},1^{n}}}");
        }
    }

    /// The partition-keyed sweep is what runs past the β-mask width, where
    /// no test can afford the table; it is pinned here at the degrees the
    /// mask route serves, against that route.
    #[test]
    fn partition_keyed_table_matches_the_mask_route() {
        for n in 1..=11u32 {
            let parts = crate::memo::partitions_cached(n);
            let blank = vec![vec![0i128; parts.len()]; parts.len()];
            let via_parts = table_on_partitions(n, &parts, blank);
            let via_masks: Vec<Vec<i128>> = kostka_table_in(n);
            assert_eq!(via_parts, via_masks, "Kostka table at degree {n}");
        }
    }

    /// The generic sweep must agree with the `u128` one entry for entry, and
    /// must run over a ring that is not an integer type at all — which is the
    /// case Kostka–Foulkes will be.
    ///
    /// `QtPoly<i64>` carrying constants gives back the same counts, so the
    /// layer is genuinely ring-agnostic rather than accidentally working for
    /// things that look like integers.
    #[test]
    fn generic_table_agrees_and_accepts_a_polynomial_ring() {
        use crate::qt::QtPoly;
        for n in 0..=9u32 {
            let plain = kostka_table(n);
            let wide: Vec<Vec<i128>> = kostka_table_in(n);
            let poly: Vec<Vec<QtPoly<i64>>> = kostka_table_in(n);
            for i in 0..plain.len() {
                for j in 0..plain.len() {
                    assert_eq!(plain[i][j], wide[i][j] as u128, "i128 at ({i},{j}) deg {n}");
                    assert_eq!(
                        poly[i][j].coeff(0, 0),
                        plain[i][j] as i64,
                        "QtPoly at ({i},{j}) deg {n}"
                    );
                    // A count lands entirely in the constant term.
                    assert!(poly[i][j].len() <= 1);
                }
            }
        }
    }

    /// The enumeration and the count are the same walk at different
    /// resolutions, so `#SSYT(λ, μ) = K_{λμ}` must hold at every pair — and the
    /// chain DP shares no code with `chains`, which reconstructs tableaux.
    #[test]
    fn tableau_count_matches_the_kostka_number() {
        for n in 0..=8u32 {
            for lambda in crate::memo::partitions_cached(n).iter() {
                for mu in crate::memo::partitions_cached(n).iter() {
                    let ts = semistandard_tableaux(lambda, mu.parts());
                    assert_eq!(
                        ts.len() as u128,
                        kostka(lambda, mu),
                        "#SSYT({lambda}, {mu})"
                    );
                }
            }
        }
    }

    /// Everything returned is a semistandard tableau of the right shape and
    /// weight, and no tableau is returned twice.
    ///
    /// The count test above cannot see any of this on its own: emitting the
    /// same tableau twice while dropping another would keep the total right.
    #[test]
    fn every_tableau_is_semistandard_of_the_given_shape_and_weight() {
        for n in 0..=7u32 {
            for lambda in crate::memo::partitions_cached(n).iter() {
                for mu in crate::memo::partitions_cached(n).iter() {
                    let ts = semistandard_tableaux(lambda, mu.parts());
                    let mut seen = std::collections::HashSet::new();
                    for t in &ts {
                        assert!(seen.insert(t.clone()), "{t:?} twice for {lambda}, {mu}");
                        let shape: Vec<u32> = t.iter().map(|r| r.len() as u32).collect();
                        assert_eq!(shape, lambda.parts(), "shape of {t:?}");
                        let mut weight = vec![0u32; mu.len()];
                        for row in t {
                            for &v in row {
                                weight[v as usize - 1] += 1;
                            }
                        }
                        assert_eq!(weight, mu.parts(), "weight of {t:?}");
                        for row in t {
                            assert!(row.windows(2).all(|w| w[0] <= w[1]), "row of {t:?}");
                        }
                        for (i, row) in t.iter().enumerate().skip(1) {
                            for (j, &v) in row.iter().enumerate() {
                                assert!(t[i - 1][j] < v, "column {j} of {t:?}");
                            }
                        }
                    }
                }
            }
        }
    }

    /// The order is part of the contract, because Sage's
    /// `SemistandardTableaux(λ, μ)` doctests print these lists verbatim. These
    /// four are copied from `sage/combinat/tableau.py`; the exhaustive check
    /// against Symmetrica itself is `scripts/check_backend.py`.
    #[test]
    fn tableau_order_matches_symmetricas() {
        assert_eq!(
            semistandard_tableaux(&p(&[3, 1]), &[2, 1, 1]),
            vec![vec![vec![1, 1, 2], vec![3]], vec![vec![1, 1, 3], vec![2]]]
        );
        assert_eq!(
            semistandard_tableaux(&p(&[3, 2, 1]), &[2, 2, 2]),
            vec![
                vec![vec![1, 1, 2], vec![2, 3], vec![3]],
                vec![vec![1, 1, 3], vec![2, 2], vec![3]],
            ]
        );
        assert_eq!(
            semistandard_tableaux(&p(&[2, 2]), &[2, 1, 1]),
            vec![vec![vec![1, 1], vec![2, 3]]]
        );
        assert_eq!(
            semistandard_tableaux(&p(&[2, 2, 2]), &[2, 2, 1, 1]),
            vec![vec![vec![1, 1], vec![2, 2], vec![3, 4]]]
        );
    }

    /// The degenerate pairs, which are convention choices rather than corner
    /// cases: one empty tableau for the empty shape, nothing when the sizes
    /// disagree. Both match `kostka_tab`, checked against Symmetrica directly.
    #[test]
    fn empty_shape_gives_one_tableau_and_mismatched_sizes_give_none() {
        assert_eq!(
            semistandard_tableaux(&p(&[]), &[]),
            vec![Vec::<Vec<u32>>::new()]
        );
        assert!(semistandard_tableaux(&p(&[2]), &[1, 1, 1]).is_empty());
        assert!(semistandard_tableaux(&p(&[1, 1]), &[2]).is_empty());
    }

    #[test]
    fn diagonal_is_one() {
        for parts in [&[2, 1][..], &[3, 1], &[2, 2], &[4]] {
            let lam = p(parts);
            assert_eq!(kostka(&lam, &lam), 1, "K_{{{lam},{lam}}}");
        }
    }
}
