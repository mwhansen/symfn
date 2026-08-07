//! `s_μ · s_ν` when one factor has exactly two rows, by counting fibres per
//! output instead of accumulating over tableaux.
//!
//! ## The idea
//!
//! With ν = (ν₁, ν₂) the LR chain is μ ⊆ λ¹ ⊆ λ, where λ¹/μ is a horizontal
//! strip of size ν₁ and λ/λ¹ one of size ν₂, subject to the lattice condition.
//! So `c^λ_{μν}` is the *number of admissible λ¹* — and iterating over λ to
//! count λ¹ directly builds no chain at all.
//!
//! Interlacing alone would pin each part independently,
//! `λ¹ᵢ ∈ [max(μᵢ, λᵢ₊₁), min(μᵢ₋₁, λᵢ)]`, making the fibre a box on the
//! hyperplane Σλ¹ = |μ|+ν₁ — a closed-form composition count. The lattice
//! condition spoils that: it says #2's in rows 1..j ≤ #1's in rows 1..j−1,
//! which in prefix sums (Λⱼ, Lⱼ, Mⱼ for λ, λ¹, μ) reads
//!
//! ```text
//!   Λⱼ − Lⱼ ≤ Lⱼ₋₁ − Mⱼ₋₁      ⟺      λ¹ⱼ ≥ Λⱼ + Mⱼ₋₁ − 2·Lⱼ₋₁
//! ```
//!
//! a constraint on *prefix sums*, so the fibre is an order polytope rather than
//! a box. What saves the method is that the admissible range for λ¹ⱼ **stays
//! contiguous** once Lⱼ₋₁ is fixed. So a DP over rows keyed on the running
//! prefix sum works, and every step is a range-add on a difference array — never
//! an enumeration. Cost is O(rows × span) per output term, span = |μ|+ν₁.
//!
//! (Taking j = 1 gives λ¹₁ ≥ λ₁, and interlacing gives λ¹₁ ≤ λ₁, so λ¹₁ = λ₁ —
//! the familiar fact that only 1's may appear in the first row.)
//!
//! ## When it wins
//!
//! Counting is O(terms × rows × span); [`SkewLr`](crate::skew_lr::SkewLr) is
//! O(tableaux). Tableaux outgrow terms, so counting wins asymptotically and
//! the gap widens with the size of the product — per-term cost grows roughly
//! linearly here where `SkewLr`'s grows faster. The crossover is measured
//! rather than derived, over products from 8 thousand to 500 thousand terms
//! (`docs/record/littlewood-richardson.md`). [`prefer_counting`] decides when
//! the dispatch turns this on; below that `SkewLr` stays in charge.

// The two *value* narrowings in this module carry their own checks at the sites
// below. Everything else is DP index arithmetic, bounded by the shape.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use crate::partition::Partition;

/// The two factors ordered so the second has exactly two rows.
///
/// `s_μ·s_ν = s_ν·s_μ`, and the chain formulation walks ν's rows, so the
/// two-row factor has to be the one supplying the strips.
fn orient<'a>(a: &'a Partition, b: &'a Partition) -> Option<(&'a Partition, &'a Partition)> {
    if b.len() == 2 {
        Some((a, b))
    } else if a.len() == 2 {
        Some((b, a))
    } else {
        None
    }
}

/// Whether [`two_row_product`] would handle this pair.
pub fn applies(a: &Partition, b: &Partition) -> bool {
    orient(a, b).is_some()
}

/// Whether counting is expected to *beat*
/// [`SkewLr`](crate::skew_lr::SkewLr) here.
///
/// Counting costs O(candidates × rows × span) and `SkewLr` O(tableaux), so the
/// crossover depends on how large the fibres are, which is not known before
/// computing them. These bounds are therefore empirical, fitted to
/// `examples/calibrate_two_row.rs` and **deliberately conservative**: every
/// measured loss is excluded, at the cost of declining two measured wins
/// (`s[110,66]²`, ruled out by `rows ≥ 3`, and
/// `s[20,18,16,14,12,10,8,6]·s[20,16]`, ruled out by `rows ≤ 6`). The worst
/// case that still dispatches did not regress.
///
/// Each clause corresponds to a failure mode that was measured, not guessed,
/// and names the shape that exhibits it:
///
/// * `rows ≥ 3` — one-row μ is Pieri, which `SkewLr` already does cheaply
///   (`s[160]·s[80,50]`).
/// * `rows ≤ 6` — many rows is `SkewLr`'s best regime, where its
///   compression compounds (`s[14,…,5]·s[14,11]`, ℓ = 10).
/// * `4·ν₂ ≥ ν₁` — a lopsided ν makes most candidates vanish, so the candidate
///   sweep stops paying for itself (`s[30,24,18]·s[40,2]`).
/// * `3·|ν| ≥ |μ|` — span is driven by |μ|, so a small output cannot amortize
///   it (`s[30,24,18]·s[6,5]`).
/// * `2·|μ| ≥ |ν|` — the mirror case, few tableaux for `SkewLr` to walk
///   (`s[12,9,6]·s[60,48]`).
/// * `|μ|+|ν| ≥ 75` — below this the whole product is small enough that
///   `SkewLr`'s constant factors win (`s[10,8,6]·s[10,8]`).
///
/// Fitted to ~24 measured pairs, so treat it as a starting point rather than a
/// law; widening it needs new measurements, not reasoning.
pub fn prefer_counting(a: &Partition, b: &Partition) -> bool {
    let Some((mu, nu)) = orient(a, b) else {
        return false;
    };
    let rows = mu.len();
    let (n1, n2) = (nu.part(0), nu.part(1));
    let (m, v) = (mu.size(), nu.size());
    (3..=6).contains(&rows) && 4 * n2 >= n1 && 3 * v >= m && 2 * m >= v && m + v >= 75
}

/// `c^λ_{μν}` when one factor has exactly two rows, else `None`.
///
/// This is the fibre count for a single λ: O(rows × span), no enumeration and
/// no product expansion.
///
/// Returns `Some(0)` when |λ| ≠ |μ|+|ν| or μ ⊄ λ, where the coefficient is
/// zero by theorem. `None` means only that neither factor has two rows.
///
/// Not currently wired into [`AutoLr`](crate::strip_lr::AutoLr), whose
/// `lr_coeff` first peeks the product cache so that a caller sweeping many λ
/// against one (μ, ν) pays for the expansion once.
///
/// # Panics
///
/// Panics if an intermediate fibre count in the row DP does not fit `i128`.
pub fn two_row_coeff(lambda: &Partition, a: &Partition, b: &Partition) -> Option<u128> {
    let (mu, nu) = orient(a, b)?;
    if lambda.size() != mu.size() + nu.size() || !lambda.contains(mu) {
        return Some(0);
    }
    let mu_v: Vec<u32> = mu.parts().to_vec();
    let mu_size: u32 = mu_v.iter().sum();
    let n1 = nu.part(0);
    let span = (mu_size + n1) as usize + 1;
    let mut st = Fibre {
        mu: &mu_v,
        target: i64::from(mu_size) + i64::from(n1),
        cur: vec![0u128; span + 1],
        nxt: vec![0i128; span + 2],
    };
    Some(st.count(lambda.parts()))
}

/// `s_a · s_b` when one factor has exactly two rows, else `None`.
///
/// Output matches [`LrBackend::schur_product`](crate::lr::LrBackend::schur_product):
/// sorted by partition, zero coefficients omitted.
///
/// # Panics
///
/// Panics if an intermediate fibre count in the row DP does not fit `i128`.
pub fn two_row_product(a: &Partition, b: &Partition) -> Option<Vec<(Partition, u128)>> {
    let (mu, nu) = orient(a, b)?;
    let mu: Vec<u32> = mu.parts().to_vec();
    let (n1, n2) = (nu.part(0), nu.part(1));

    let mu_size: u32 = mu.iter().sum();
    let span = (mu_size + n1) as usize + 1;
    let mut st = Fibre {
        mu: &mu,
        target: i64::from(mu_size) + i64::from(n1),
        cur: vec![0u128; span + 1],
        nxt: vec![0i128; span + 2],
    };

    // λ ⊇ μ with |λ| = |μ|+|ν|, at most two new rows, and λⱼ ≤ μⱼ₋₂ — the last
    // because λⱼ ≤ λ¹ⱼ₋₁ ≤ μⱼ₋₂ through the two strips. These bounds are tight
    // enough that almost every candidate has a non-zero coefficient.
    let rows = mu.len() + 2;
    let mut cand = vec![0u32; rows];
    let mut out = Vec::new();
    walk(0, rows, n1 + n2, u32::MAX, &mut cand, &mut st, &mut out);
    out.sort_by(|x: &(Partition, u128), y| x.0.cmp(&y.0));
    Some(out)
}

/// Scratch shared by every fibre count in one product.
struct Fibre<'a> {
    mu: &'a [u32],
    /// Σλ¹, fixed by the first strip's size.
    target: i64,
    cur: Vec<u128>,
    nxt: Vec<i128>,
}

impl Fibre<'_> {
    fn mu_at(&self, i: usize) -> u32 {
        self.mu.get(i).copied().unwrap_or(0)
    }

    /// `c^λ_{μν}` — the number of admissible λ¹ — by a DP over rows whose state
    /// is the running prefix sum Lⱼ₋₁.
    fn count(&mut self, lam: &[u32]) -> u128 {
        let target = self.target;
        let span = target as usize + 1;
        self.cur[..=span].fill(0);
        self.cur[0] = 1; // L₀ = 0

        // Reachable prefix sums stay in a shrinking window; tracking it keeps
        // the scan proportional to what is live rather than to `span`.
        let (mut win_lo, mut win_hi) = (0usize, 0usize);
        let (mut big_lam, mut big_mu) = (0i64, 0i64); // Λⱼ and Mⱼ₋₁

        for (i, &lam_i) in lam.iter().enumerate() {
            big_lam += i64::from(lam_i);
            let lo_i =
                i64::from(self.mu_at(i)).max(i64::from(lam.get(i + 1).copied().unwrap_or(0)));
            // λ¹ᵢ ≤ μᵢ₋₁; past μ's length that bound is 0, correctly forcing
            // λ¹ᵢ = 0 there.
            let hi_i = if i == 0 {
                i64::from(lam[0])
            } else {
                i64::from(self.mu_at(i - 1)).min(i64::from(lam_i))
            };
            if hi_i < lo_i {
                return 0;
            }

            self.nxt[..=span + 1].fill(0);
            let (mut new_lo, mut new_hi) = (usize::MAX, 0usize);
            for l_prev in win_lo..=win_hi {
                let ways = self.cur[l_prev];
                if ways == 0 {
                    continue;
                }
                let need = big_lam + big_mu - 2 * l_prev as i64;
                let lo = lo_i.max(need);
                if lo > hi_i {
                    continue;
                }
                let s = l_prev as i64 + lo;
                if s > target {
                    continue;
                }
                let e = (l_prev as i64 + hi_i).min(target);
                // The one *value* narrowing in this DP: `ways` is a tableau
                // count with no a-priori bound, and the difference array is
                // signed because it subtracts. One compare per entry, in a loop
                // dominated by the two writes below it.
                let ways = i128::try_from(ways).unwrap_or_else(|_| {
                    panic!("a two-row multiplicity of {ways} does not fit i128")
                });
                self.nxt[s as usize] += ways;
                self.nxt[e as usize + 1] -= ways;
                new_lo = new_lo.min(s as usize);
                new_hi = new_hi.max(e as usize);
            }
            if new_lo == usize::MAX {
                return 0; // nothing survived this row
            }

            let mut run: i128 = 0;
            for l in 0..=new_hi {
                run += self.nxt[l];
                // A prefix sum of the difference array, so it is a count of
                // fillings and non-negative by construction; the cast is exact
                // for every value that invariant admits.
                debug_assert!(run >= 0, "a tableau count went negative at {l}");
                self.cur[l] = run as u128;
            }
            win_lo = new_lo;
            win_hi = new_hi;
            big_mu += i64::from(self.mu_at(i));
        }
        if win_lo <= target as usize && target as usize <= win_hi {
            self.cur[target as usize]
        } else {
            0
        }
    }
}

/// Enumerate candidate λ and count each fibre.
fn walk(
    j: usize,
    rows: usize,
    left: u32,
    prev: u32,
    cand: &mut Vec<u32>,
    st: &mut Fibre,
    out: &mut Vec<(Partition, u128)>,
) {
    if j == rows {
        if left != 0 {
            return;
        }
        let end = cand.iter().rposition(|&x| x > 0).map_or(0, |i| i + 1);
        let c = st.count(&cand[..end]);
        if c > 0 {
            out.push((Partition::new(cand[..end].iter().copied()), c));
        }
        return;
    }
    let lo = st.mu_at(j);
    let hi = if j < 2 {
        prev.min(left + lo)
    } else {
        st.mu_at(j - 2).min(prev).min(left + lo)
    };
    if hi < lo {
        return;
    }
    for v in lo..=hi {
        cand[j] = v;
        walk(j + 1, rows, left - (v - lo), v, cand, st, out);
    }
    cand[j] = 0;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lr::LrBackend;
    use crate::skew_lr::SkewLr;

    fn p(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    /// Agrees with the general engine wherever it applies, in both argument
    /// orders (the orientation swap must not be observable).
    #[test]
    fn agrees_with_skew_lr_on_two_row_factors() {
        let mut cases = 0;
        for n1 in 1..=5u32 {
            for n2 in 1..=n1 {
                let nu = p(&[n1, n2]);
                for mu_parts in [
                    &[1u32][..],
                    &[2],
                    &[2, 1],
                    &[3, 1],
                    &[2, 2],
                    &[3, 2, 1],
                    &[4, 2, 1],
                    &[5, 3, 2],
                    &[3, 3, 3],
                    &[4, 4, 2, 1],
                ] {
                    let mu = p(mu_parts);
                    let want = SkewLr.schur_product(&mu, &nu);
                    assert_eq!(
                        two_row_product(&mu, &nu).as_ref(),
                        Some(&want),
                        "s{mu} · s{nu}"
                    );
                    assert_eq!(
                        two_row_product(&nu, &mu).as_ref(),
                        Some(&SkewLr.schur_product(&nu, &mu)),
                        "s{nu} · s{mu} (swapped)"
                    );
                    cases += 1;
                }
            }
        }
        assert!(cases >= 100, "expected a real sweep, got {cases}");
    }

    /// Both factors two rows — the orientation is then ambiguous, and either
    /// choice must give the same answer.
    #[test]
    fn two_row_times_two_row_is_orientation_independent() {
        for a in 1..=6u32 {
            for b in 1..=a {
                for c in 1..=6u32 {
                    for d in 1..=c {
                        let (x, y) = (p(&[a, b]), p(&[c, d]));
                        let want = SkewLr.schur_product(&x, &y);
                        assert_eq!(two_row_product(&x, &y).as_ref(), Some(&want), "s{x}·s{y}");
                        assert_eq!(
                            two_row_product(&y, &x).as_ref(),
                            Some(&SkewLr.schur_product(&y, &x))
                        );
                    }
                }
            }
        }
    }

    /// Degenerate factors still have to come out right: an empty μ makes the
    /// product s_ν, and a single-row μ is Pieri.
    #[test]
    fn degenerate_factors() {
        let nu = p(&[3, 2]);
        assert_eq!(
            two_row_product(&Partition::default(), &nu),
            Some(vec![(nu.clone(), 1)])
        );
        for mu in [p(&[1]), p(&[4]), p(&[7])] {
            assert_eq!(
                two_row_product(&mu, &nu).as_ref(),
                Some(&SkewLr.schur_product(&mu, &nu)),
                "s{mu} · s{nu}"
            );
        }
    }

    /// Must decline rather than answer wrongly when neither factor has two rows.
    #[test]
    fn declines_without_a_two_row_factor() {
        for (a, b) in [
            (p(&[3, 2, 1]), p(&[2, 1, 1])),
            (p(&[4]), p(&[3])),
            (p(&[2, 2, 2]), p(&[1])),
            (Partition::default(), Partition::default()),
        ] {
            assert_eq!(two_row_product(&a, &b), None, "s{a} · s{b}");
            assert!(!applies(&a, &b));
        }
    }

    /// The predicate must match what was measured — these are the exact shapes
    /// from the calibration sweep, so a change to the bounds that silently
    /// re-admits a known loss fails here.
    #[test]
    fn dispatch_predicate_matches_the_calibration() {
        let cases: &[(&[u32], &[u32], bool)] = &[
            (&[20, 16, 12], &[20, 16], true),                         // 1.28x
            (&[40, 32, 24], &[40, 32], true),                         // 1.66x
            (&[24, 20, 16, 12], &[24, 20], true),                     // 1.61x
            (&[30, 24, 18], &[30, 24], true),                         // 1.16x
            (&[10, 8, 6], &[10, 8], false),                           // 0.92x — too small
            (&[160], &[80, 50], false),                               // 0.17x — one-row μ
            (&[30, 24, 18], &[40, 2], false),                         // 0.08x — lopsided ν
            (&[30, 24, 18], &[6, 5], false),                          // 0.66x — ν tiny vs μ
            (&[12, 9, 6], &[60, 48], false),                          // 0.18x — μ tiny vs ν
            (&[14, 13, 12, 11, 10, 9, 8, 7, 6, 5], &[14, 11], false), // 0.55x — tall μ
            (&[70, 42], &[70, 42], false),                            // 0.91x — two-row μ
        ];
        for (m, n, want) in cases {
            let (mu, nu) = (p(m), p(n));
            assert_eq!(prefer_counting(&mu, &nu), *want, "s{mu} · s{nu}");
            // Commutative: the predicate must not depend on argument order.
            assert_eq!(prefer_counting(&nu, &mu), *want, "s{nu} · s{mu}");
        }
    }

    /// End-to-end: on a shape the dispatch actually fires for, `AutoLr` must
    /// still agree with the general engine. The unit tests above run far below
    /// the threshold, so without this nothing checks the live path.
    #[test]
    fn dispatched_path_agrees_with_the_engine() {
        let (mu, nu) = (p(&[16, 13, 10, 7]), p(&[16, 13]));
        assert!(prefer_counting(&mu, &nu), "shape must actually dispatch");
        assert_eq!(
            crate::strip_lr::AutoLr.schur_product(&mu, &nu),
            SkewLr.schur_product(&mu, &nu)
        );
    }

    /// Coefficients above 1 are where a fibre count can go wrong; pin a case
    /// with real multiplicity against the engine, coefficient by coefficient.
    #[test]
    fn multiplicities_match() {
        let (mu, nu) = (p(&[4, 3, 2]), p(&[4, 3]));
        let got = two_row_product(&mu, &nu).expect("two-row");
        let want = SkewLr.schur_product(&mu, &nu);
        assert_eq!(got, want);
        assert!(
            got.iter().any(|(_, c)| *c >= 3),
            "expected a coefficient ≥ 3 to make this test meaningful"
        );
    }
}
