//! An optimized Littlewood–Richardson backend: a row-level dynamic program.
//!
//! [`NaiveLr`](crate::lr::NaiveLr) enumerates LR tableaux **cell by cell**. This
//! backend instead builds the chain
//!
//! ```text
//!   μ = λ⁰ ⊆ λ¹ ⊆ … ⊆ λ^ℓ(ν) = λ,   λⁱ/λⁱ⁻¹ a horizontal strip of size νᵢ
//! ```
//!
//! one **row-strip at a time**, which is a genuine algorithmic improvement for
//! two reasons: strips are enumerated at row granularity rather than per cell,
//! and distinct tableaux that reach the same state collapse into a single DP
//! entry carrying a multiplicity, so work is shared instead of repeated. It
//! also computes the *entire* product s_μ·s_ν in one pass, rather than testing
//! candidate λ one at a time.
//!
//! The lattice (Yamanouchi) condition becomes a condition on consecutive strips:
//! writing θⁱ_j for the cells added in row j at step i,
//!
//! ```text
//!   for all i ≥ 2, j ≥ 1:   Σ_{k ≤ j} θⁱ_k  ≤  Σ_{k ≤ j-1} θⁱ⁻¹_k
//! ```
//!
//! i.e. the number of i's in the first j rows never exceeds the number of
//! (i−1)'s in the first j−1 rows. (Taking j = 1 recovers the familiar fact that
//! only 1's may appear in the first row.)
//!
//! Correctness is pinned by exhaustive agreement with `NaiveLr` and by the Sage
//! oracle — exactly the cross-check the [`LrBackend`] trait was designed for.

use std::collections::HashMap;
use std::sync::Arc;

use crate::lr::LrBackend;
use crate::memo::product_cached;
use crate::partition::Partition;

/// Row-level DP backend for Littlewood–Richardson coefficients.
#[derive(Clone, Copy, Debug, Default)]
pub struct StripLr;

/// Every partition obtained from `lambda` by adding a horizontal strip of size
/// `r`, paired with the per-row counts θ of what was added.
///
/// μ ⊇ λ with μ/λ a horizontal strip means the interlacing
/// μ₁ ≥ λ₁ ≥ μ₂ ≥ λ₂ ≥ …, so at most one new cell lands in any column.
fn horizontal_strips(lambda: &[u32], r: u32) -> Vec<(Vec<u32>, Vec<u32>)> {
    let rows = lambda.len() + 1; // a strip may open one new row
    let mut out = Vec::new();
    let mut cur = vec![0u32; rows];

    fn rec(
        j: usize,
        rows: usize,
        left: u32,
        lambda: &[u32],
        cur: &mut Vec<u32>,
        out: &mut Vec<(Vec<u32>, Vec<u32>)>,
    ) {
        if j == rows {
            if left == 0 {
                let mut parts: Vec<u32> = cur.iter().copied().filter(|&x| x > 0).collect();
                let theta: Vec<u32> = (0..rows)
                    .map(|k| cur[k] - lambda.get(k).copied().unwrap_or(0))
                    .collect();
                parts.truncate(cur.iter().filter(|&&x| x > 0).count());
                out.push((parts, theta));
            }
            return;
        }
        let lam_j = lambda.get(j).copied().unwrap_or(0);
        // Row j may grow up to the previous row's *old* value (interlacing),
        // and is unbounded above only for j = 0.
        let upper = if j == 0 {
            lam_j + left
        } else {
            let prev_old = lambda.get(j - 1).copied().unwrap_or(0);
            prev_old.min(lam_j + left)
        };
        for v in lam_j..=upper {
            let added = v - lam_j;
            if added > left {
                break;
            }
            cur[j] = v;
            rec(j + 1, rows, left - added, lambda, cur, out);
        }
        cur[j] = 0;
    }

    rec(0, rows, r, lambda, &mut cur, &mut out);
    out
}

/// The lattice condition between consecutive strips (see module docs).
fn lattice_ok(theta: &[u32], prev: &[u32]) -> bool {
    let n = theta.len().max(prev.len());
    let mut run_theta = 0u32;
    let mut run_prev = 0u32;
    for j in 0..n {
        run_theta += theta.get(j).copied().unwrap_or(0);
        // Σ_{k ≤ j} θⁱ_k ≤ Σ_{k ≤ j-1} θⁱ⁻¹_k
        if run_theta > run_prev {
            return false;
        }
        run_prev += prev.get(j).copied().unwrap_or(0);
    }
    true
}

impl StripLr {
    /// The full expansion of s_μ · s_ν, computed in a single DP pass.
    ///
    /// Returns only the nonzero terms, sorted by λ.
    pub fn product(&self, mu: &Partition, nu: &Partition) -> Vec<(Partition, u128)> {
        (*self.product_shared(mu, nu)).clone()
    }

    /// [`product`](Self::product) without the copy: the memoized expansion
    /// itself. See [`expand_skew_shared`](crate::skew_lr::expand_skew_shared) —
    /// the deep clone is one allocation per term, for data the cache already
    /// holds.
    pub fn product_shared(&self, mu: &Partition, nu: &Partition) -> Arc<Vec<(Partition, u128)>> {
        product_cached(mu, nu, || self.product_uncached(mu, nu))
    }

    fn product_uncached(&self, mu: &Partition, nu: &Partition) -> Vec<(Partition, u128)> {
        // State: (current partition, strip added at the previous step) -> count.
        let mut states: HashMap<(Vec<u32>, Vec<u32>), u128> = HashMap::new();
        states.insert((mu.parts().to_vec(), Vec::new()), 1);

        for (i, &r) in nu.parts().iter().enumerate() {
            let mut next: HashMap<(Vec<u32>, Vec<u32>), u128> = HashMap::new();
            for ((cur, prev), mult) in &states {
                for (grown, theta) in horizontal_strips(cur, r) {
                    // The first strip is unconstrained; later ones must satisfy
                    // the ballot condition against the strip before them.
                    if i > 0 && !lattice_ok(&theta, prev) {
                        continue;
                    }
                    *next.entry((grown, theta)).or_insert(0) += mult;
                }
            }
            states = next;
        }

        // Collapse the auxiliary strip component; sum multiplicities per shape.
        let mut totals: HashMap<Vec<u32>, u128> = HashMap::new();
        for ((shape, _), mult) in states {
            *totals.entry(shape).or_insert(0) += mult;
        }
        let mut out: Vec<(Partition, u128)> = totals
            .into_iter()
            .map(|(parts, c)| (Partition::from_sorted(parts), c))
            .collect();
        out.sort_by(|a, b| a.0.cmp(&b.0));
        out
    }
}

impl LrBackend for StripLr {
    fn lr_coeff(&self, lambda: &Partition, mu: &Partition, nu: &Partition) -> u128 {
        if lambda.size() != mu.size() + nu.size() || !lambda.contains(mu) {
            return 0;
        }
        self.product_shared(mu, nu)
            .iter()
            .find(|(l, _)| l == lambda)
            .map(|(_, c)| *c)
            .unwrap_or(0)
    }

    fn schur_product(&self, mu: &Partition, nu: &Partition) -> Vec<(Partition, u128)> {
        self.product(mu, nu)
    }

    fn schur_product_shared(&self, mu: &Partition, nu: &Partition) -> Arc<Vec<(Partition, u128)>> {
        self.product_shared(mu, nu)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lr::NaiveLr;
    use crate::partition::partitions_of;

    fn p(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    #[test]
    fn known_products() {
        // s_1·s_1 = s_2 + s_{11}
        let got = StripLr.product(&p(&[1]), &p(&[1]));
        let shapes: Vec<(Vec<u32>, u128)> =
            got.iter().map(|(l, c)| (l.parts().to_vec(), *c)).collect();
        assert_eq!(shapes, vec![(vec![1, 1], 1), (vec![2], 1)]);

        // the multiplicity-2 case
        assert_eq!(
            StripLr.lr_coeff(&p(&[3, 2, 1]), &p(&[2, 1]), &p(&[2, 1])),
            2
        );
    }

    #[test]
    fn agrees_with_naive_backend_exhaustively() {
        // Every pair with |mu| + |nu| <= 7, full expansions compared.
        let mut checked = 0;
        for a in 0..=7u32 {
            for b in 0..=(7 - a) {
                for mu in partitions_of(a) {
                    for nu in partitions_of(b) {
                        let fast = StripLr.schur_product(&mu, &nu);
                        let slow = NaiveLr.schur_product(&mu, &nu);
                        assert_eq!(fast, slow, "s{mu}*s{nu}");
                        checked += 1;
                    }
                }
            }
        }
        assert!(checked > 200, "expected a real sweep, got {checked}");
    }

    /// The shared form of a product is the owned one without the copy, on
    /// every backend — the provided default and each override.
    #[test]
    fn schur_product_shared_agrees_with_schur_product_on_every_backend() {
        use crate::skew_lr::SkewLr;
        for a in 0..=5u32 {
            for b in 0..=(5 - a) {
                for mu in partitions_of(a) {
                    for nu in partitions_of(b) {
                        let want = NaiveLr.schur_product(&mu, &nu);
                        assert_eq!(
                            *NaiveLr.schur_product_shared(&mu, &nu),
                            want,
                            "Naive s{mu}*s{nu}"
                        );
                        assert_eq!(
                            *StripLr.schur_product_shared(&mu, &nu),
                            want,
                            "Strip s{mu}*s{nu}"
                        );
                        assert_eq!(
                            *SkewLr.schur_product_shared(&mu, &nu),
                            want,
                            "Skew s{mu}*s{nu}"
                        );
                        assert_eq!(
                            *AutoLr.schur_product_shared(&mu, &nu),
                            want,
                            "Auto s{mu}*s{nu}"
                        );
                    }
                }
            }
        }
    }

    /// Whatever route computes a product — the rectangle form, two-row or
    /// three-row counting, or the layer — `AutoLr` stores it under the one
    /// entry `SkewLr` keeps for that product, in either argument order. So a
    /// repeat is a lookup, and `lr_coeff`, which peeks that entry, answers a
    /// sweep of coefficients from it, zeros included.
    ///
    /// The pairs are chosen to fire each route (asserted, so a recalibration
    /// of the predicates cannot silently retarget the test) and to be used by
    /// no other test, since the table is process-wide.
    #[test]
    fn every_route_of_auto_lr_is_memoized_under_one_entry() {
        use crate::skew_lr::SkewLr;
        let rect = (p(&[3, 3, 3]), p(&[2, 2]));
        let two = (p(&[15, 13, 11]), p(&[19, 17]));
        let three = (p(&[10, 8, 6]), p(&[9, 8, 7]));
        let general = (p(&[5, 3, 2, 1]), p(&[4, 2, 1]));
        assert!(crate::rect::okada_product(&rect.0, &rect.1).is_some());
        assert!(crate::two_row::prefer_counting(&two.0, &two.1));
        assert!(crate::two_row::two_row_product(&two.0, &two.1).is_some());
        assert!(crate::three_row::prefer_counting(&three.0, &three.1));
        assert!(crate::three_row::three_row_product(&three.0, &three.1).is_some());
        assert!(!crate::two_row::prefer_counting(&general.0, &general.1));
        assert!(!crate::three_row::prefer_counting(&general.0, &general.1));

        for (mu, nu) in [rect, two, three, general] {
            let first = AutoLr.schur_product_shared(&mu, &nu);
            let again = AutoLr.schur_product_shared(&mu, &nu);
            assert!(Arc::ptr_eq(&first, &again), "s{mu}·s{nu} recomputed");
            let swapped = AutoLr.schur_product_shared(&nu, &mu);
            assert!(Arc::ptr_eq(&first, &swapped), "s{nu}·s{mu} stored apart");
            let layer = SkewLr.schur_product_shared(&mu, &nu);
            assert!(
                Arc::ptr_eq(&first, &layer),
                "SkewLr keeps s{mu}·s{nu} apart"
            );

            // Every λ of the degree where that is affordable; past that, the
            // product's own terms and two shapes that contain μ, have the
            // right size, and are zero by theorem — too many rows, and ν's
            // cells all in one row — so the peek is asked for zeros too.
            let n = mu.size() + nu.size();
            let sweep: Vec<Partition> = if n <= 20 {
                partitions_of(n)
            } else {
                let mut v: Vec<Partition> = first.iter().map(|(l, _)| l.clone()).collect();
                let mut tall = mu.parts().to_vec();
                tall[0] += nu.part(0);
                tall.extend(std::iter::repeat_n(1, (nu.size() - nu.part(0)) as usize));
                v.push(Partition::new(tall));
                let mut wide = mu.parts().to_vec();
                wide[0] += nu.size();
                v.push(Partition::new(wide));
                v
            };
            let mut seen = 0;
            for lambda in &sweep {
                let want = first
                    .binary_search_by(|(l, _)| l.cmp(lambda))
                    .map_or(0, |i| first[i].1);
                assert_eq!(
                    AutoLr.lr_coeff(lambda, &mu, &nu),
                    want,
                    "c^{lambda}_{{{mu},{nu}}}"
                );
                seen += usize::from(want != 0);
            }
            assert_eq!(seen, first.len(), "s{mu}·s{nu}: the sweep missed terms");
            assert!(seen < sweep.len(), "s{mu}·s{nu}: no zero was asked for");
        }
    }
}

/// The backend the library uses by default: a closed form when both factors
/// are rectangles, a counting route on some two- and three-row products, and
/// [`SkewLr`](crate::skew_lr::SkewLr) otherwise.
///
/// [`two_row::prefer_counting`](crate::two_row::prefer_counting) and
/// [`three_row::prefer_counting`](crate::three_row::prefer_counting) select the
/// counting routes. Only
/// [`schur_product`](crate::lr::LrBackend::schur_product) takes them;
/// `lr_coeff` chooses between the rectangle form and `SkewLr`.
///
/// Every product is memoized under one key whatever route computed it
/// (`skew_lr::memoized_product`), so a repeated product costs a lookup and
/// `lr_coeff`'s peek at the product cache answers a sweep of coefficients
/// off any of the routes.
///
/// ## The rectangle path
///
/// `s_{(aᵖ)}·s_{(bᑫ)}` is multiplicity-free with an explicitly describable
/// support (Okada 1998), so [`crate::rect`] *generates* the answer instead of
/// searching for it — and this is precisely the shape the general engine is
/// worst at, since `skew_lr`'s conjugate-dispatch heuristic records rectangles
/// as a known loss. The gain grows with the product, and a *single*
/// coefficient gains most of all, since the predicate is O(ℓ(λ)) and replaces
/// the whole search (`docs/record/littlewood-richardson.md`).
///
/// ## The general path
///
/// [`SkewLr`](crate::skew_lr::SkewLr) beats both
/// [`NaiveLr`](crate::lr::NaiveLr) and [`StripLr`] at every size measured,
/// down to every pair with |μ|+|ν| ≤ 12, so there is no crossover to dispatch
/// on and no small-input regime where the layer map's overhead loses
/// (`examples/bench_lr.rs`, `docs/record/littlewood-richardson.md`).
///
/// `AutoLr` is the one place dispatch lives. All backends are verified
/// equivalent — against each other exhaustively, and against both Sage and
/// lrcalc — so the choice is unobservable except in timing; [`crate::rect`] is
/// pinned to `SkewLr` over every small rectangle pair, both argument orders,
/// and every λ including the zeros.
#[derive(Clone, Copy, Debug, Default)]
pub struct AutoLr;

/// Total degree |μ|+|ν| at or above which [`StripLr`] beats
/// [`NaiveLr`](crate::lr::NaiveLr).
///
/// Not a dispatch threshold: [`SkewLr`](crate::skew_lr::SkewLr) beats both at
/// every size, so nothing chooses between these two. Kept as the documented
/// crossover.
pub const STRIP_THRESHOLD: u32 = 36;

impl LrBackend for AutoLr {
    fn lr_coeff(&self, lambda: &Partition, mu: &Partition, nu: &Partition) -> u128 {
        if let Some(c) = crate::rect::okada_coeff(lambda, mu, nu) {
            return c;
        }
        crate::skew_lr::SkewLr.lr_coeff(lambda, mu, nu)
    }

    fn schur_product(&self, mu: &Partition, nu: &Partition) -> Vec<(Partition, u128)> {
        (*self.schur_product_shared(mu, nu)).clone()
    }

    fn schur_product_shared(&self, mu: &Partition, nu: &Partition) -> Arc<Vec<(Partition, u128)>> {
        crate::skew_lr::memoized_product(mu, nu, || {
            // Rectangles first: a closed form beats a fiber count.
            if let Some(v) = crate::rect::okada_product(mu, nu) {
                return Some(v);
            }
            if crate::two_row::prefer_counting(mu, nu) {
                if let Some(v) = crate::two_row::two_row_product(mu, nu) {
                    return Some(v);
                }
            }
            if crate::three_row::prefer_counting(mu, nu) {
                if let Some(v) = crate::three_row::three_row_product(mu, nu) {
                    return Some(v);
                }
            }
            None
        })
    }
}
