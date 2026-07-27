//! `s_μ · s_ν` when one factor has exactly three rows, by counting fibres per
//! output instead of accumulating over tableaux.
//!
//! The same idea as [`crate::two_row`], one strip deeper. With ν = (ν₁,ν₂,ν₃)
//! the LR chain is μ ⊆ λ¹ ⊆ λ² ⊆ λ, so `c^λ_{μν}` counts admissible *pairs*
//! (λ¹, λ²) and the DP state has to carry both.
//!
//! ## Why the 2-D state is affordable
//!
//! Two things keep it small, and neither is obvious:
//!
//! 1. **Key the state on cells added, not on absolute prefix sums.** Writing
//!    `aⱼ = Σ_{k≤j} (λ¹ₖ − μₖ)` and `bⱼ = Σ_{k≤j} (λ²ₖ − λ¹ₖ)`, these are
//!    bounded by ν₁ and ν₂ — not by |μ|+ν₁, which is what the prefix sums
//!    themselves range over. On `[20,16,12]²` that is 21×17 rather than 69×85.
//! 2. **λ¹ still has known bounds.** Chaining the three interlacings gives
//!    `λ¹ᵢ ≤ μᵢ₋₁` and `λ¹ᵢ ≥ λᵢ₊₂`, both in terms of μ and λ, so only λ² needs
//!    the previous row — leaving the state (λ¹ⱼ, aⱼ, bⱼ).
//!
//! The two lattice conditions then read, with `cⱼ` the cells added by strip 3,
//!
//! ```text
//!   bⱼ ≤ aⱼ₋₁        (#2's in rows 1..j ≤ #1's in rows 1..j−1)
//!   cⱼ ≤ bⱼ₋₁        (#3's in rows 1..j ≤ #2's in rows 1..j−1)
//! ```
//!
//! and `cⱼ = (Λⱼ − Mⱼ) − aⱼ − bⱼ` is determined, so it costs no state.
//!
//! Predicted cost was thousands of operations per term, which would have been
//! hopeless. Measured is 26–2192: the reachable state space is far smaller than
//! its bounding box.
//!
//! ## When it wins
//!
//! Counting is O(candidates × states) and the frontier is O(tableaux), so this
//! wins asymptotically. Measured against [`SkewLr`](crate::skew_lr::SkewLr):
//!
//! ```text
//!   [12,10,8]²      6 579 terms    5.7ms    5.4ms   1.06x
//!   [14,12,10]²    12 068 terms    9.4ms    8.4ms   1.13x
//!   [20,16,12]²    64 335 terms   116ms    85.9ms   1.35x
//!   [22,18,14]²    99 208 terms   223ms     145ms   1.54x
//!   [24,20,16]²   145 505 terms   433ms     233ms   1.86x
//! ```
//!
//! [`prefer_counting`] decides when the dispatch turns it on.

use crate::partition::Partition;

/// The two factors ordered so the second has exactly three rows.
///
/// `s_μ·s_ν = s_ν·s_μ`, and the chain walks ν's rows, so the three-row factor
/// has to be the one supplying the strips.
/// When *both* have three rows the choice is free, so it must be canonical or
/// the result would depend on argument order. Taking the smaller as ν is also
/// the cheaper option: the state table is (ν₁+1)×(ν₂+1)×(max λ₁+1) and the
/// candidate sweep spans |ν| cells, so both shrink with it.
fn orient<'a>(a: &'a Partition, b: &'a Partition) -> Option<(&'a Partition, &'a Partition)> {
    match (a.len() == 3, b.len() == 3) {
        (true, true) => Some(if (a.size(), a.parts()) <= (b.size(), b.parts()) {
            (b, a)
        } else {
            (a, b)
        }),
        (false, true) => Some((a, b)),
        (true, false) => Some((b, a)),
        (false, false) => None,
    }
}

/// Whether counting is expected to beat the frontier here.
///
/// Empirical, in the same spirit as [`crate::two_row::prefer_counting`].
/// Measured on `s_μ²` for three-row μ, on AC power:
///
/// ```text
///   n =  18   30    36    42    48    54    60    72    84   108   120
///       1.57 0.82  0.74  0.87  0.88  1.05  1.06  1.13  1.04  1.54  1.86
/// ```
///
/// ⚠️ **Those in-process figures are optimistic and must not set the bound.**
/// Each ran `SkewLr` first and counting second, so counting inherited a warm
/// allocator. Measured out-of-process against lrcalc — one cold process per
/// side, which is what a caller actually sees — `n = 60` is a *regression*:
/// `[12,10,8]²` went 0.79x → 0.53x when the dispatch was set that low, while
/// in-process it had looked like 1.06x. `n ≥ 90` is set from the
/// out-of-process numbers instead:
///
/// ```text
///   [12,10,8]²   n = 60   0.79x -> 0.53x   regressed, excluded
///   [14,12,10]²  n = 72   0.75x -> 0.77x   neutral,   excluded
///   [20,16,12]²  n = 96   0.63x -> 0.79x   clear gain, admitted
/// ```
///
/// The lesson generalises: an in-process A/B that always runs the same side
/// second measures allocator state as much as algorithm.
///
/// Also requires a balanced ν and comparable factor sizes for the same reasons
/// the two-row predicate does — a lopsided or tiny ν makes most candidates
/// vanish, so the candidate sweep stops paying for itself.
pub fn prefer_counting(a: &Partition, b: &Partition) -> bool {
    let Some((mu, nu)) = orient(a, b) else {
        return false;
    };
    let rows = mu.len();
    let (n1, n3) = (nu.part(0), nu.part(2));
    let (m, v) = (mu.size(), nu.size());
    (3..=5).contains(&rows) && 3 * n3 >= n1 && 3 * v >= m && 3 * m >= v && m + v >= 90
}

/// `s_a · s_b` when one factor has exactly three rows, else `None`.
pub fn three_row_product(a: &Partition, b: &Partition) -> Option<Vec<(Partition, u128)>> {
    let (mu_p, nu_p) = orient(a, b)?;
    let mu: Vec<u32> = mu_p.parts().to_vec();
    let nu: Vec<u32> = nu_p.parts().to_vec();

    // Widest candidate puts every cell of ν in row 0.
    let max_l1 = (mu_p.part(0) + nu_p.size()) as usize;
    let mut st = Fibre {
        mu: &mu,
        nu: &nu,
        cur: Table::new(max_l1, nu[0] as usize, nu[1] as usize),
        nxt: Table::new(max_l1, nu[0] as usize, nu[1] as usize),
        live: Vec::new(),
    };

    // λ ⊇ μ with |λ| = |μ|+|ν|, at most three new rows, and λⱼ ≤ μⱼ₋₃ — the
    // last because λⱼ ≤ λ²ⱼ₋₁ ≤ λ¹ⱼ₋₂ ≤ μⱼ₋₃ through the three strips.
    let rows = mu.len() + 3;
    let mut cand = vec![0u32; rows];
    let mut out = Vec::new();
    walk(0, rows, nu_p.size(), u32::MAX, &mut cand, &mut st, &mut out);
    out.sort_by(|x: &(Partition, u128), y| x.0.cmp(&y.0));
    Some(out)
}

/// Dense generation-stamped state table: O(1) insert with no per-row clearing.
///
/// A first version used a `HashMap` keyed on the state tuple and ran 3.3x
/// *slower* at identical operation counts — the algorithm was right and the
/// data structure was wrong.
struct Table {
    gen: Vec<u32>,
    val: Vec<u128>,
    touched: Vec<u32>,
    era: u32,
    stride_a: usize,
    stride_l: usize,
}

impl Table {
    fn new(max_l1: usize, n1: usize, n2: usize) -> Table {
        let stride_a = n2 + 1;
        let stride_l = (n1 + 1) * stride_a;
        let size = (max_l1 + 1) * stride_l;
        Table {
            gen: vec![0; size],
            val: vec![0; size],
            touched: Vec::with_capacity(256),
            era: 0,
            stride_a,
            stride_l,
        }
    }
    fn clear(&mut self) {
        self.era += 1;
        self.touched.clear();
    }
    fn add(&mut self, l1: usize, a: usize, b: usize, ways: u128) {
        let i = l1 * self.stride_l + a * self.stride_a + b;
        if self.gen[i] != self.era {
            self.gen[i] = self.era;
            self.val[i] = 0;
            self.touched.push(i as u32);
        }
        self.val[i] += ways;
    }
    fn decode(&self, i: usize, n1: usize) -> (usize, usize, usize) {
        (
            i / self.stride_l,
            (i / self.stride_a) % (n1 + 1),
            i % self.stride_a,
        )
    }
}

struct Fibre<'a> {
    mu: &'a [u32],
    nu: &'a [u32],
    cur: Table,
    nxt: Table,
    live: Vec<(usize, usize, usize, u128)>,
}

fn at(v: &[u32], i: usize) -> i64 {
    v.get(i).copied().map_or(0, i64::from)
}

impl Fibre<'_> {
    /// `c^λ_{μν}`: the number of admissible pairs (λ¹, λ²).
    fn count(&mut self, lam: &[u32]) -> u128 {
        let (n1, n2) = (at(self.nu, 0), at(self.nu, 1));
        self.cur.clear();
        // λ¹ from a nonexistent previous row is unbounded; λ's first part is an
        // exact stand-in, since every later bound is a `min` against it.
        self.cur.add(at(lam, 0) as usize, 0, 0, 1);

        let (mut big_lam, mut big_mu) = (0i64, 0i64);
        for j in 0..lam.len() {
            big_lam += at(lam, j);
            big_mu += at(self.mu, j);
            let l1_lo = at(self.mu, j).max(at(lam, j + 2));
            let l1_hi = if j == 0 {
                at(lam, 0)
            } else {
                at(self.mu, j - 1).min(at(lam, j))
            };

            self.live.clear();
            for &t in &self.cur.touched {
                let (l1, a, b) = self.cur.decode(t as usize, n1 as usize);
                self.live.push((l1, a, b, self.cur.val[t as usize]));
            }
            self.nxt.clear();
            for &(l1_prev, a_prev, b_prev, ways) in self.live.iter() {
                let (a_prev, b_prev) = (a_prev as i64, b_prev as i64);
                for l1 in l1_lo..=l1_hi {
                    let a = a_prev + l1 - at(self.mu, j);
                    if a > n1 {
                        break;
                    }
                    let l2_lo = l1.max(at(lam, j + 1));
                    let l2_hi = (l1_prev as i64).min(at(lam, j));
                    for l2 in l2_lo..=l2_hi {
                        let b = b_prev + l2 - l1;
                        // Lattice: #2's in rows 1..j ≤ #1's in rows 1..j−1.
                        if b > n2 || b > a_prev {
                            break;
                        }
                        // Strip 3's cells are determined, so they cost no state.
                        let c = (big_lam - big_mu) - a - b;
                        if c < 0 || c > b_prev {
                            continue;
                        }
                        self.nxt.add(l1 as usize, a as usize, b as usize, ways);
                    }
                }
            }
            std::mem::swap(&mut self.cur, &mut self.nxt);
            if self.cur.touched.is_empty() {
                return 0;
            }
        }
        let mut total = 0u128;
        for &t in &self.cur.touched {
            let (_, a, b) = self.cur.decode(t as usize, n1 as usize);
            if a as i64 == n1 && b as i64 == n2 {
                total += self.cur.val[t as usize];
            }
        }
        total
    }
}

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
    let lo = at(st.mu, j) as u32;
    let hi = if j < 3 {
        prev.min(left + lo)
    } else {
        (at(st.mu, j - 3) as u32).min(prev).min(left + lo)
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

    /// Agree with the general engine wherever it applies, in both argument
    /// orders — the orientation swap must not be observable.
    #[test]
    fn agrees_with_skew_lr_on_three_row_factors() {
        let mut cases = 0;
        for n1 in 1..=4u32 {
            for n2 in 1..=n1 {
                for n3 in 1..=n2 {
                    let nu = p(&[n1, n2, n3]);
                    for mu_parts in [
                        &[1u32][..],
                        &[3],
                        &[2, 1],
                        &[3, 2],
                        &[2, 2, 1],
                        &[4, 2, 1],
                        &[3, 3, 3],
                        &[4, 3, 2, 1],
                        &[5, 3, 2, 1],
                    ] {
                        let mu = p(mu_parts);
                        assert_eq!(
                            three_row_product(&mu, &nu).as_ref(),
                            Some(&SkewLr.schur_product(&mu, &nu)),
                            "s{mu} · s{nu}"
                        );
                        assert_eq!(
                            three_row_product(&nu, &mu).as_ref(),
                            Some(&SkewLr.schur_product(&nu, &mu)),
                            "s{nu} · s{mu} (swapped)"
                        );
                        cases += 1;
                    }
                }
            }
        }
        assert!(cases >= 100, "expected a real sweep, got {cases}");
    }

    /// Both factors three rows: the orientation is free, so it must be chosen
    /// canonically — otherwise the answer, and the dispatch decision, would
    /// depend on argument order. Covers distinct pairs, not just squares.
    #[test]
    fn two_three_row_factors_are_orientation_independent() {
        let shapes = [
            p(&[1, 1, 1]), p(&[3, 2, 1]), p(&[4, 4, 2]), p(&[5, 3, 1]), p(&[2, 2, 2]),
        ];
        for x in &shapes {
            for y in &shapes {
                let want = SkewLr.schur_product(x, y);
                assert_eq!(three_row_product(x, y).as_ref(), Some(&want), "s{x}·s{y}");
                assert_eq!(
                    three_row_product(y, x).as_ref(),
                    Some(&SkewLr.schur_product(y, x)),
                    "s{y}·s{x} (swapped)"
                );
                assert_eq!(
                    prefer_counting(x, y),
                    prefer_counting(y, x),
                    "dispatch must not depend on argument order: s{x}·s{y}"
                );
            }
        }
    }

    #[test]
    fn degenerate_factors() {
        let nu = p(&[3, 2, 1]);
        assert_eq!(
            three_row_product(&Partition::default(), &nu),
            Some(vec![(nu.clone(), 1)])
        );
        for mu in [p(&[1]), p(&[5]), p(&[9])] {
            assert_eq!(
                three_row_product(&mu, &nu).as_ref(),
                Some(&SkewLr.schur_product(&mu, &nu)),
                "s{mu} · s{nu}"
            );
        }
    }

    #[test]
    fn declines_without_a_three_row_factor() {
        for (a, b) in [
            (p(&[3, 2, 1, 1]), p(&[2, 1])),
            (p(&[4]), p(&[3])),
            (p(&[2, 2]), p(&[1])),
            (Partition::default(), Partition::default()),
        ] {
            assert_eq!(three_row_product(&a, &b), None, "s{a} · s{b}");
        }
    }

    /// Large coefficients are where a fibre count goes wrong; pin one against
    /// the engine on a case with real multiplicity.
    #[test]
    fn multiplicities_match() {
        let (mu, nu) = (p(&[5, 4, 3]), p(&[4, 3, 2]));
        let got = three_row_product(&mu, &nu).expect("three-row");
        assert_eq!(got, SkewLr.schur_product(&mu, &nu));
        assert!(
            got.iter().any(|(_, c)| *c >= 5),
            "expected a coefficient ≥ 5 to make this test meaningful"
        );
    }

    /// End-to-end on a shape the dispatch actually fires for. Every test above
    /// sits far below the threshold, so without this nothing covers the live
    /// `AutoLr` route.
    #[test]
    fn dispatched_path_agrees_with_the_engine() {
        let mu = p(&[20, 16, 12]);
        assert!(prefer_counting(&mu, &mu), "shape must actually dispatch");
        assert_eq!(
            crate::strip_lr::AutoLr.schur_product(&mu, &mu),
            SkewLr.schur_product(&mu, &mu)
        );
    }

    /// The predicate must keep excluding the shapes measured as losses.
    #[test]
    fn dispatch_predicate_matches_the_calibration() {
        let cases: &[(&[u32], &[u32], bool)] = &[
            (&[20, 16, 12], &[20, 16, 12], true),  // 1.35x
            (&[24, 20, 16], &[24, 20, 16], true),  // 1.86x
            (&[22, 18, 14], &[22, 18, 14], true),  // 1.54x
            (&[14, 12, 10], &[14, 12, 10], false), // n = 72, out-of-process neutral
            (&[12, 10, 8], &[12, 10, 8], false),   // n = 60, out-of-process regression
            (&[10, 8, 6], &[10, 8, 6], false),
            (&[8, 6, 4], &[8, 6, 4], false),
            (&[30, 24, 18], &[40, 2, 1], false),   // lopsided ν
            (&[30, 24, 18], &[3, 2, 1], false),    // ν tiny against μ
        ];
        for (m, n, want) in cases {
            let (mu, nu) = (p(m), p(n));
            assert_eq!(prefer_counting(&mu, &nu), *want, "s{mu} · s{nu}");
            assert_eq!(prefer_counting(&nu, &mu), *want, "s{nu} · s{mu}");
        }
    }
}
