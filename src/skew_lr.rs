//! Single-traversal Littlewood–Richardson: expand a whole skew Schur function
//! in one enumeration.
//!
//! Both existing backends answer the question "given λ, μ, ν, what is c^λ_{μν}?"
//! and then *repeat the whole computation* for every candidate output shape.
//! [`NaiveLr::schur_product`](crate::lr::NaiveLr) loops over all p(n) partitions
//! of n and runs a full backtrack for each; [`skew_schur`](crate::hopf::skew_schur)
//! does the same. That is the dominant inefficiency, and it is structural rather
//! than constant-factor.
//!
//! This backend inverts the question. It enumerates the LR skew tableaux of a
//! shape **once** and bins each one by the content it happens to have, yielding
//!
//! ```text
//!   s_{outer/inner} = Σ_ν c^{outer}_{inner,ν} s_ν
//! ```
//!
//! in a single pass — every nonzero coefficient, from one traversal.
//!
//! Three things make the traversal itself cheap:
//!
//! 1. **Content doubles as the ballot check.** Instead of tracking a target
//!    content, a cap, and a separate running ballot count, we keep one vector
//!    `cont` of value multiplicities. A cell may take value v exactly when
//!    `cont[v] < cont[v-1]`, which simultaneously enforces the lattice word
//!    condition and keeps `cont` a partition. At the end of a branch `cont`
//!    *is* the content — no separate bookkeeping, no post-fill scan.
//!
//! 2. **Neighbours are precomputed as indices.** The row and column constraints
//!    need the cell to the right and the cell above. Those are resolved once
//!    into a flat cell array, so the inner loop is array indexing rather
//!    than repeated partition lookups and bounds arithmetic.
//!
//! 3. **The value bound is per-row.** In an LR tableau the entry in row r
//!    (0-indexed) is at most r+1, since row 0 can hold only 1s and columns
//!    increase strictly. Rows near the top branch very little.
//!
//! Products come from the same primitive. The skew shape
//!
//! ```text
//!     . . . │ ν₁          (μ₁ columns of "inner", then ν)
//!     . . . │ ν₂
//!     μ₁                  (μ, flush left)
//!     μ₂
//! ```
//!
//! is **disconnected** — the ν block and the μ block share no row and no column
//! — and a skew Schur function of a disconnected shape is the product of its
//! pieces. So s_μ·s_ν is one call to [`expand_skew`] on that shape, and the
//! entire product falls out of a single traversal.
//!
//! Finally, the shape is normalized before counting: since (λ/μ)' = λ'/μ' and
//! conjugating the shape conjugates every output partition, we transpose
//! whenever that reduces the row count, because rows are what drive branching.
//!
//! Correctness is pinned by exhaustive agreement with
//! [`NaiveLr`](crate::lr::NaiveLr) and, independently, against the `lrcalc`
//! binary (see `scripts/oracle_lrcalc.sh`).
//!
//! *Provenance:* the algorithmic ideas here — one traversal binned by content,
//! the content-as-ballot test, shape normalization — are the standard ones,
//! used by Anders Buch's `lrcalc` among others. This is an independent Rust
//! implementation, not a translation; symfn shares no code with `lrcalc` and
//! stays MIT/Apache-2.0.

use std::collections::HashMap;

use crate::lr::LrBackend;
use crate::memo::skew_cached;
use crate::partition::Partition;

/// A skew cell, with its constraint neighbours already resolved to indices into
/// the cell array. `NONE` means "no such neighbour, so no bound from it".
#[derive(Clone, Copy)]
struct Cell {
    /// Cell to the right (row: weakly increasing left→right).
    right: u32,
    /// Cell above (column: strictly increasing top→bottom).
    above: u32,
    /// Largest value allowed here: row index + 1.
    max: u32,
}

const NONE: u32 = u32::MAX;

/// Row-level DP backend replacement: expand a skew Schur function in one pass.
#[derive(Clone, Copy, Debug, Default)]
pub struct SkewLr;

/// Expand s_{outer/inner} = Σ_ν c^{outer}_{inner,ν} s_ν, sorted by ν.
///
/// Returns empty when `inner ⊄ outer`. The degree-0 case (outer == inner) is the
/// empty partition with coefficient 1, i.e. s_{λ/λ} = 1.
pub fn expand_skew(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
    if !outer.contains(inner) {
        return Vec::new();
    }
    (*skew_cached(outer, inner, || expand_skew_uncached(outer, inner))).clone()
}

fn expand_skew_uncached(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
    // Branching is driven by the number of rows (row r admits values 1..=r+1),
    // so count in the transposed orientation when it is shorter. Conjugating the
    // shape conjugates every output partition, so we simply undo it afterwards.
    if outer.len() > outer.part(0) as usize {
        let mut out: Vec<(Partition, u128)> = count_shape(&outer.conjugate(), &inner.conjugate())
            .into_iter()
            .map(|(nu, c)| (nu.conjugate(), c))
            .collect();
        out.sort_by(|a, b| a.0.cmp(&b.0));
        return out;
    }
    let mut out = count_shape(outer, inner);
    out.sort_by(|a, b| a.0.cmp(&b.0));
    out
}

/// The traversal proper: enumerate LR tableaux of `outer/inner`, bin by content.
fn count_shape(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)> {
    let rows = outer.len();
    let width = outer.part(0) as usize;

    // Cells in reading-word order: top→bottom, right→left. In this order a
    // cell's right neighbour and its neighbour above are both already filled,
    // so every constraint prunes at the moment of placement.
    let mut pos = vec![NONE; rows * width];
    let mut coords: Vec<(usize, usize)> = Vec::new();
    for r in 0..rows {
        let start = inner.part(r) as usize;
        let end = outer.part(r) as usize;
        for c in (start..end).rev() {
            pos[r * width + c] = coords.len() as u32;
            coords.push((r, c));
        }
    }

    let cells: Vec<Cell> = coords
        .iter()
        .map(|&(r, c)| Cell {
            right: if c + 1 < outer.part(r) as usize {
                pos[r * width + c + 1]
            } else {
                NONE
            },
            above: if r > 0
                && c >= inner.part(r - 1) as usize
                && c < outer.part(r - 1) as usize
            {
                pos[(r - 1) * width + c]
            } else {
                NONE
            },
            max: r as u32 + 1,
        })
        .collect();

    // cont[v] = how many v's placed so far. cont[0] is a sentinel fixed at
    // "infinity" so the ballot test `cont[v] < cont[v-1]` needs no special case
    // for v = 1 (a 1 is always placeable).
    let mut cont = vec![0u32; rows + 2];
    cont[0] = u32::MAX;
    let mut vals = vec![0u32; cells.len()];
    let mut acc: HashMap<Vec<u32>, u128> = HashMap::new();

    rec(0, &cells, &mut vals, &mut cont, &mut acc);

    acc.into_iter()
        .map(|(parts, c)| (Partition::from_sorted(parts), c))
        .collect()
}

fn rec(
    i: usize,
    cells: &[Cell],
    vals: &mut [u32],
    cont: &mut [u32],
    acc: &mut HashMap<Vec<u32>, u128>,
) {
    if i == cells.len() {
        // `cont` is weakly decreasing by the ballot test, so it is already a
        // partition; trim the trailing zeros and record it. Hashing the slice
        // avoids allocating for contents we have already seen.
        let len = cont[1..].iter().rposition(|&x| x != 0).map_or(0, |k| k + 1);
        let key = &cont[1..=len];
        if let Some(n) = acc.get_mut(key) {
            *n += 1;
        } else {
            acc.insert(key.to_vec(), 1);
        }
        return;
    }

    let cell = cells[i];
    // Row weakly increasing (we fill right→left), capped by the per-row bound.
    let hi = if cell.right != NONE {
        vals[cell.right as usize].min(cell.max)
    } else {
        cell.max
    };
    // Column strictly increasing.
    let lo = if cell.above != NONE {
        vals[cell.above as usize] + 1
    } else {
        1
    };

    for v in lo..=hi {
        // Ballot / lattice-word condition, which is exactly "cont stays a
        // partition": we may add a v only while v-1 is strictly ahead.
        if cont[v as usize] >= cont[v as usize - 1] {
            continue;
        }
        vals[i] = v;
        cont[v as usize] += 1;
        rec(i + 1, cells, vals, cont, acc);
        cont[v as usize] -= 1;
    }
}

/// The disconnected skew shape whose skew Schur function is s_μ · s_ν.
///
/// ν sits in the top ℓ(ν) rows, pushed right past μ₁ empty columns; μ sits
/// flush left in the rows below. The two blocks share no column, so the shape
/// is a disjoint union and its skew Schur function factors as the product.
fn product_shape(mu: &Partition, nu: &Partition) -> (Partition, Partition) {
    let w = mu.part(0);
    let mut outer = Vec::with_capacity(nu.len() + mu.len());
    let mut inner = Vec::with_capacity(nu.len());
    for i in 0..nu.len() {
        outer.push(w + nu.part(i));
        inner.push(w);
    }
    outer.extend_from_slice(mu.parts());
    // outer is weakly decreasing: w+ν₁ ≥ … ≥ w+ν_last ≥ w = μ₁ ≥ μ₂ ≥ …
    // inner is (w,…,w) which is weakly decreasing and contained in outer.
    (
        Partition::from_sorted(outer),
        Partition::new(inner.into_iter()),
    )
}

impl SkewLr {
    /// The full expansion of s_μ · s_ν, from a single traversal.
    pub fn product(&self, mu: &Partition, nu: &Partition) -> Vec<(Partition, u128)> {
        let (outer, inner) = product_shape(mu, nu);
        expand_skew(&outer, &inner)
    }
}

impl LrBackend for SkewLr {
    fn lr_coeff(&self, lambda: &Partition, mu: &Partition, nu: &Partition) -> u128 {
        if lambda.size() != mu.size() + nu.size() || !lambda.contains(mu) {
            return 0;
        }
        // One traversal of λ/μ yields every ν at once, and it is cached on
        // (λ, μ), so a sweep over ν costs a single enumeration in total.
        expand_skew(lambda, mu)
            .into_iter()
            .find(|(n, _)| n == nu)
            .map(|(_, c)| c)
            .unwrap_or(0)
    }

    fn schur_product(&self, mu: &Partition, nu: &Partition) -> Vec<(Partition, u128)> {
        self.product(mu, nu)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lr::NaiveLr;
    use crate::partition::partitions_of;
    use crate::strip_lr::StripLr;

    fn p(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    #[test]
    fn known_products() {
        let got = SkewLr.product(&p(&[1]), &p(&[1]));
        let shapes: Vec<(Vec<u32>, u128)> =
            got.iter().map(|(l, c)| (l.parts().to_vec(), *c)).collect();
        assert_eq!(shapes, vec![(vec![1, 1], 1), (vec![2], 1)]);

        assert_eq!(SkewLr.lr_coeff(&p(&[3, 2, 1]), &p(&[2, 1]), &p(&[2, 1])), 2);
    }

    #[test]
    fn skew_expansion_matches_lrcalc_documented_example() {
        // `lrcalc skew 3 2 1 / 2 1` -> 2 (2,1); 1 (1,1,1); 1 (3)
        let got = expand_skew(&p(&[3, 2, 1]), &p(&[2, 1]));
        let shapes: Vec<(Vec<u32>, u128)> =
            got.iter().map(|(l, c)| (l.parts().to_vec(), *c)).collect();
        assert_eq!(
            shapes,
            vec![(vec![1, 1, 1], 1), (vec![2, 1], 2), (vec![3], 1)]
        );
    }

    #[test]
    fn empty_skew_is_the_unit() {
        // s_{λ/λ} = 1, i.e. the empty partition with coefficient 1.
        let got = expand_skew(&p(&[3, 2]), &p(&[3, 2]));
        assert_eq!(got, vec![(Partition::default(), 1)]);
        // inner ⊄ outer is the zero expansion.
        assert!(expand_skew(&p(&[2, 1]), &p(&[3])).is_empty());
    }

    #[test]
    fn agrees_with_other_backends_exhaustively() {
        // Every pair with |mu| + |nu| <= 7, full expansions compared.
        let mut checked = 0;
        for a in 0..=7u32 {
            for b in 0..=(7 - a) {
                for mu in partitions_of(a) {
                    for nu in partitions_of(b) {
                        let fast = SkewLr.schur_product(&mu, &nu);
                        assert_eq!(fast, NaiveLr.schur_product(&mu, &nu), "s{mu}*s{nu}");
                        assert_eq!(fast, StripLr.schur_product(&mu, &nu), "s{mu}*s{nu}");
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 200, "expected a real sweep, got {checked}");
    }

    #[test]
    fn skew_expansion_agrees_with_naive_exhaustively() {
        // Directly exercise the skew primitive (not just products) against the
        // coefficient-at-a-time backend, including disconnected inner shapes.
        let mut checked = 0;
        for n in 0..=8u32 {
            for outer in partitions_of(n) {
                for m in 0..=n {
                    for inner in partitions_of(m) {
                        if !outer.contains(&inner) {
                            continue;
                        }
                        let fast = expand_skew(&outer, &inner);
                        let mut slow: Vec<(Partition, u128)> = partitions_of(n - m)
                            .into_iter()
                            .filter_map(|nu| {
                                let c = NaiveLr.lr_coeff(&outer, &inner, &nu);
                                (c != 0).then_some((nu, c))
                            })
                            .collect();
                        slow.sort_by(|a, b| a.0.cmp(&b.0));
                        assert_eq!(fast, slow, "s_{{{outer}/{inner}}}");
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 200, "expected a real sweep, got {checked}");
    }

    #[test]
    fn conjugate_symmetry_of_the_shape_normalization() {
        // The transpose path and the direct path must agree: expanding a tall
        // shape exercises the conjugation branch, a wide one does not.
        let tall = expand_skew(&p(&[2, 2, 2, 1, 1]), &p(&[1]));
        let wide = expand_skew(&p(&[5, 3]), &p(&[1]));
        // c^{λ'}_{μ'ν'} = c^λ_{μν}: conjugating shape and contents is an involution.
        let back: Vec<(Partition, u128)> = {
            let mut v: Vec<(Partition, u128)> = tall
                .iter()
                .map(|(n, c)| (n.conjugate(), *c))
                .collect();
            v.sort_by(|a, b| a.0.cmp(&b.0));
            v
        };
        let direct = {
            let mut v: Vec<(Partition, u128)> = expand_skew(
                &p(&[2, 2, 2, 1, 1]).conjugate(),
                &p(&[1]).conjugate(),
            );
            v.sort_by(|a, b| a.0.cmp(&b.0));
            v
        };
        assert_eq!(back, direct);
        assert!(!wide.is_empty());
    }
}
