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
//! Measured cost is 26–2192 operations per term: the reachable state space is
//! far smaller than its bounding box (`docs/record/littlewood-richardson.md`).
//!
//! ## When it wins
//!
//! Counting is O(candidates × states) and `SkewLr` is O(tableaux), so this
//! wins asymptotically — and the packed state, with the candidates counted in
//! parallel across workers (`candidates.rs`), makes the constants
//! competitive from n ≈ 48 up, which is where [`prefer_counting`] turns the
//! dispatch on. Calibrated against [`SkewLr`](crate::skew_lr::SkewLr) — itself
//! parallel over each row's states, so wall time against wall time — over
//! products from 3 thousand to 1.5 million terms
//! (`docs/record/littlewood-richardson.md`).

// Every `as` here is DP index arithmetic — window bounds, row lengths, and the
// `i64` offsets that let a bound go negative before it is clamped. All are
// bounded by the shape, whose parts are `u32`. The counts themselves are
// `u128` and are never cast.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

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

/// Whether counting is expected to beat
/// [`SkewLr`](crate::skew_lr::SkewLr) here.
///
/// Empirical, in the same spirit as [`crate::two_row::prefer_counting`], and
/// fitted the same way: out of process, interleaved `lr_cli` runs with counting
/// forced on and off, min of 5, on AC, one cold process per run — because an
/// in-process A/B hands whichever side runs second a warm allocator. Both
/// routes are parallel — the candidates over workers (`candidates.rs`),
/// the layer over each row's states — so the comparison is wall time against
/// wall time. The sweep is `s_μ²` for three-row μ plus asymmetric pairs with
/// three- to twelve-row μ. At n = 36 counting and `SkewLr` tie at the process
/// floor, and that shape stays below the bound; the margin grows from there,
/// 1.1x at n = 48 to 3.7x at n = 144. The ladder, the absolute times and the
/// recalibration story live in `docs/record/littlewood-richardson.md`.
///
/// `rows ≥ 3` because a one-row μ is Pieri and a two-row μ is the two-row
/// route's; there is no upper bound, since six- to twelve-row μ measured
/// 1.2–10x (`s[12,11,10,9,8,7]·s[6,5,4]` to `s[12,…,1]·s[9,7,5]`). The
/// balance clauses are the two-row predicate's reasons: a lopsided or tiny ν
/// makes most candidates vanish, so the candidate sweep stops paying for
/// itself. `8·|ν| ≥ |μ|` sits past the measured win at |μ|/|ν| = 6.3
/// (`s[14,…,5]·s[6,5,4]`, 2.15x) and short of the tie at ratio 12
/// (`s[30,24,18]·s[3,2,1]`); nothing between them is measured.
pub fn prefer_counting(a: &Partition, b: &Partition) -> bool {
    let Some((mu, nu)) = orient(a, b) else {
        return false;
    };
    let rows = mu.len();
    let (n1, n3) = (nu.part(0), nu.part(2));
    let (m, v) = (mu.size(), nu.size());
    rows >= 3 && 3 * n3 >= n1 && 8 * v >= m && 3 * m >= v && m + v >= 48
}

/// `s_a · s_b` when one factor has exactly three rows, else `None`.
///
/// Output matches
/// [`LrBackend::schur_product`](crate::lr::LrBackend::schur_product): sorted by
/// partition, zero coefficients omitted.
///
/// Also declines — same `None`, and the caller's fallback engine answers —
/// when the three-row factor is wider than 1023 or a candidate first row could
/// exceed 4095, the widths the packed state representation carries.
pub fn three_row_product(a: &Partition, b: &Partition) -> Option<Vec<(Partition, u128)>> {
    let (mu_p, nu_p) = orient(a, b)?;
    let mu: Vec<u32> = mu_p.parts().to_vec();
    let nu: Vec<u32> = nu_p.parts().to_vec();

    // Widest candidate puts every cell of ν in row 0.
    let max_l1 = (mu_p.part(0) + nu_p.size()) as usize;
    // States live in one packed u32 (λ¹ in 12 bits, a and b in 10 each), so
    // decoding costs shifts rather than the divisions a flat index needs.
    // Shapes past those widths decline; the caller falls back to the general
    // engine, which owns that regime anyway — a product this size has more
    // pressing costs than dispatch.
    if nu[0] > 1023 || max_l1 > 4095 {
        return None;
    }
    // The candidates — λ ⊇ μ with |λ| = |μ|+|ν|, at most three new rows, and
    // λⱼ ≤ μⱼ₋₃ through the three strips — are `candidates::walk` with three
    // strips; each is counted by a `Fibre`, one per worker.
    let threads = std::thread::available_parallelism().map_or(1, |n| n.get());
    Some(crate::candidates::count_all(
        &mu,
        3,
        nu_p.size(),
        threads,
        || {
            let mut st = Fibre::new(&mu, &nu, max_l1);
            move |lam: &[u32]| st.count(lam)
        },
    ))
}

/// Dense generation-stamped state table: O(1) insert with no per-row clearing.
///
/// `touched` holds each live state as a packed `(λ¹ << 20) | (a << 10) | b`
/// rather than its flat cell index: unpacking is then three shift-masks where a
/// flat index costs two integer divisions by run-time strides
/// (`examples/calibrate_three_row.rs`; the measurement is in
/// `docs/record/littlewood-richardson.md`). The 12/10/10 split is why
/// [`three_row_product`] declines ν₁ ≥ 1024 or first-row candidates ≥ 4096.
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
    #[inline]
    fn add(&mut self, l1: usize, a: usize, b: usize, ways: u128) {
        let i = l1 * self.stride_l + a * self.stride_a + b;
        if self.gen[i] != self.era {
            self.gen[i] = self.era;
            self.val[i] = 0;
            self.touched.push(((l1 << 20) | (a << 10) | b) as u32);
        }
        self.val[i] += ways;
    }
    #[inline]
    fn get(&self, l1: usize, a: usize, b: usize) -> u128 {
        self.val[l1 * self.stride_l + a * self.stride_a + b]
    }
}

#[inline]
fn unpack(t: u32) -> (i64, i64, i64) {
    (
        i64::from(t >> 20),
        i64::from((t >> 10) & 0x3ff),
        i64::from(t & 0x3ff),
    )
}

struct Fibre<'a> {
    mu: &'a [u32],
    nu: &'a [u32],
    cur: Table,
    nxt: Table,
}

impl<'a> Fibre<'a> {
    fn new(mu: &'a [u32], nu: &'a [u32], max_l1: usize) -> Fibre<'a> {
        Fibre {
            mu,
            nu,
            cur: Table::new(max_l1, nu[0] as usize, nu[1] as usize),
            nxt: Table::new(max_l1, nu[0] as usize, nu[1] as usize),
        }
    }
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
            let mu_j = at(self.mu, j);
            let l1_lo = mu_j.max(at(lam, j + 2));
            let l1_hi = if j == 0 {
                at(lam, 0)
            } else {
                at(self.mu, j - 1).min(at(lam, j))
            };
            let (lam_j, lam_j1) = (at(lam, j), at(lam, j + 1));
            let d = big_lam - big_mu;

            self.nxt.clear();
            let (cur, nxt) = (&self.cur, &mut self.nxt);
            for &t in &cur.touched {
                let (l1_prev, a_prev, b_prev) = unpack(t);
                let ways = cur.get(l1_prev as usize, a_prev as usize, b_prev as usize);
                let l2_hi = l1_prev.min(lam_j);
                if l2_hi < lam_j1 {
                    continue; // no admissible λ²ⱼ for any λ¹ⱼ: λ²ⱼ ≥ λⱼ₊₁ already fails
                }
                // a ≤ ν₁ and l2_lo = max(l1, λⱼ₊₁) ≤ l2_hi, both monotone in l1.
                let l1_top = l1_hi.min(mu_j + n1 - a_prev).min(l2_hi);
                for l1 in l1_lo..=l1_top {
                    let a = a_prev + l1 - mu_j;
                    // b = b_prev + λ²ⱼ − λ¹ⱼ sweeps an interval as λ²ⱼ does; the
                    // lattice conditions (b ≤ aⱼ₋₁, c = D−a−b ∈ [0, bⱼ₋₁]) and
                    // b ≤ ν₂ clamp that interval rather than puncture it.
                    let b_lo = (b_prev + l1.max(lam_j1) - l1).max(d - a - b_prev);
                    let b_hi = (b_prev + l2_hi - l1).min(n2).min(a_prev).min(d - a);
                    for b in b_lo..=b_hi {
                        nxt.add(l1 as usize, a as usize, b as usize, ways);
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
            let (l1, a, b) = unpack(t);
            if a == n1 && b == n2 {
                total += self.cur.get(l1 as usize, a as usize, b as usize);
            }
        }
        total
    }
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
            p(&[1, 1, 1]),
            p(&[3, 2, 1]),
            p(&[4, 4, 2]),
            p(&[5, 3, 1]),
            p(&[2, 2, 2]),
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

    /// The packed state carries λ¹ in 12 bits and a, b in 10 each, so shapes
    /// past those widths must decline rather than truncate: a wide ν, and a μ
    /// whose first row plus |ν| overflows the λ¹ field. Just inside the bound
    /// still answers.
    #[test]
    fn declines_shapes_wider_than_the_packed_state() {
        // A two-row μ forces the wide factor into the ν role.
        assert_eq!(three_row_product(&p(&[5, 4]), &p(&[1030, 2, 1])), None);
        assert_eq!(three_row_product(&p(&[4090, 8, 4]), &p(&[8, 6, 4])), None);
        assert!(three_row_product(&p(&[4070, 8, 4]), &p(&[8, 6, 4])).is_some());
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

    /// The predicate must keep admitting the shapes measured as wins and
    /// excluding the ties and losses.
    #[test]
    fn dispatch_predicate_matches_the_calibration() {
        let cases: &[(&[u32], &[u32], bool)] = &[
            (&[20, 16, 12], &[20, 16, 12], true),    // 2.08x out-of-process
            (&[24, 20, 16], &[24, 20, 16], true),    // 2.31x
            (&[22, 18, 14], &[22, 18, 14], true),    // 1.99x
            (&[14, 12, 10], &[14, 12, 10], true),    // 1.73x
            (&[12, 10, 8], &[12, 10, 8], true),      // 1.44x
            (&[10, 8, 6], &[10, 8, 6], true),        // 1.12x
            (&[16, 13, 10, 7], &[8, 6, 4], true),    // 1.48x, four-row μ
            (&[14, 12, 10, 8, 6], &[7, 5, 3], true), // 1.37x, five-row μ
            (&[12, 11, 10, 9, 8, 7], &[6, 5, 4], true), // 1.24x, six-row μ
            (&[16, 14, 12, 10, 8, 6, 4], &[10, 8, 6], true), // 9.96x, seven-row μ
            (&[12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1], &[9, 7, 5], true), // 7.68x, twelve-row μ
            (&[14, 13, 12, 11, 10, 9, 8, 7, 6, 5], &[6, 5, 4], true), // 2.15x, |μ|/|ν| = 6.3
            (&[8, 6, 4], &[8, 6, 4], false),         // n = 36, tie at the floor
            (&[30, 24, 18], &[40, 2, 1], false),     // lopsided ν
            (&[30, 24, 18], &[3, 2, 1], false),      // ν tiny against μ, 0.98x
        ];
        for (m, n, want) in cases {
            let (mu, nu) = (p(m), p(n));
            assert_eq!(prefer_counting(&mu, &nu), *want, "s{mu} · s{nu}");
            assert_eq!(prefer_counting(&nu, &mu), *want, "s{nu} · s{mu}");
        }
    }
}
