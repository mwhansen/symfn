//! LR with a two-row ν: counting per output vs the frontier DP.
//!
//! This is the smallest LR case that still carries the lattice condition, and
//! it settles whether counting fibres per output beats accumulating over
//! tableaux. It does, asymptotically — crossover near `[40,32,24]·[40,32]`,
//! 3.7x by `[60,48,36]·[60,48]` and widening.
//!
//! Note the baseline that matters is `AutoLr`, not the naive chain enumerator
//! also implemented here: the frontier DP is far better than chain enumeration,
//! so measuring against chains flatters this method by an order of magnitude.
//!
//! Chain: μ ⊆ λ¹ ⊆ λ, horizontal strips of sizes ν₁ and ν₂, plus lattice
//! (#2's in rows 1..j ≤ #1's in rows 1..j−1).
//!
//!   (A) enumerate both strips, filter by lattice, accumulate.
//!   (B) for each candidate λ, count the λ¹ with a DP over rows keyed on the
//!       running prefix sum L. Interlacing gives λ¹ᵢ ∈ [max(μᵢ,λᵢ₊₁),
//!       min(μᵢ₋₁,λᵢ)] and lattice gives λ¹ⱼ ≥ Λⱼ + Mⱼ₋₁ − 2Lⱼ₋₁, so the
//!       admissible range stays contiguous and each DP step is a range-add
//!       (difference array), never an enumeration.
//!
//! Both are checked against the library.

use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

use symfn::{AutoLr, LrBackend, Partition};

macro_rules! p {
    ($($a:tt)*) => {{ println!($($a)*); std::io::stdout().flush().ok(); }};
}

fn h_strips(lambda: &[u32], r: u32) -> Vec<Vec<u32>> {
    let rows = lambda.len() + 1;
    let mut out = Vec::new();
    let mut cur = vec![0u32; rows];
    fn rec(j: usize, rows: usize, left: u32, lam: &[u32], cur: &mut Vec<u32>, out: &mut Vec<Vec<u32>>) {
        if j == rows {
            if left == 0 {
                let mut v = cur.clone();
                while v.last() == Some(&0) {
                    v.pop();
                }
                out.push(v);
            }
            return;
        }
        let lam_j = lam.get(j).copied().unwrap_or(0);
        let upper = if j == 0 {
            lam_j + left
        } else {
            lam.get(j - 1).copied().unwrap_or(0).min(lam_j + left)
        };
        for v in lam_j..=upper {
            cur[j] = v;
            rec(j + 1, rows, left - (v - lam_j), lam, cur, out);
        }
        cur[j] = 0;
    }
    rec(0, rows, r, lambda, &mut cur, &mut out);
    out
}

/// (A) Enumerate chains, keep the lattice-admissible ones.
fn by_chains(mu: &[u32], nu: &[u32]) -> (Vec<(Vec<u32>, u128)>, u64) {
    let mut acc: HashMap<Vec<u32>, u128> = HashMap::new();
    let mut work = 0u64;
    for l1 in h_strips(mu, nu[0]) {
        for l2 in h_strips(&l1, nu[1]) {
            work += 1;
            // lattice: for every j, Σ_{k≤j} θ²_k ≤ Σ_{k≤j-1} θ¹_k
            let n = l2.len();
            let (mut t2, mut t1_prev, mut ok) = (0i64, 0i64, true);
            for j in 0..n {
                let a = i64::from(l1.get(j).copied().unwrap_or(0));
                let b = i64::from(l2[j]);
                let m = i64::from(mu.get(j).copied().unwrap_or(0));
                t2 += b - a;
                if t2 > t1_prev {
                    ok = false;
                    break;
                }
                t1_prev += a - m;
            }
            if ok {
                *acc.entry(l2).or_insert(0) += 1;
            }
        }
    }
    let mut v: Vec<_> = acc.into_iter().collect();
    v.sort();
    (v, work)
}

/// Count the λ¹ for one λ: DP over rows, state = running prefix sum of λ¹.
/// `cur`/`nxt` are caller-owned scratch so the hot path never allocates.
fn count_fibre(mu: &[u32], lam: &[u32], nu1: u32, cur: &mut [u128], nxt: &mut [i128]) -> u128 {
    let mu_size: i64 = mu.iter().map(|&x| i64::from(x)).sum();
    let target = mu_size + i64::from(nu1); // Σλ¹
    let span = target as usize + 1;
    cur[..span].fill(0);
    cur[0] = 1; // L₀ = 0

    let rows = lam.len();
    let (mut big_lam, mut big_mu) = (0i64, 0i64); // Λ_j and M_{j-1}
    for i in 0..rows {
        big_lam += i64::from(lam[i]);
        let lo_i = i64::from(mu.get(i).copied().unwrap_or(0))
            .max(i64::from(lam.get(i + 1).copied().unwrap_or(0)));
        // λ¹ᵢ ≤ μᵢ₋₁ (horizontal strip); past μ's length that bound is 0, which
        // correctly forces λ¹ᵢ = 0 there.
        let hi_i = if i == 0 {
            i64::from(lam[0])
        } else {
            i64::from(mu.get(i - 1).copied().unwrap_or(0)).min(i64::from(lam[i]))
        };
        nxt[..=span].fill(0);
        if hi_i >= lo_i {
            for l_prev in 0..span {
                let ways = cur[l_prev];
                if ways == 0 {
                    continue;
                }
                // lattice: λ¹ᵢ ≥ Λ_j + M_{j-1} − 2·L_{j-1}
                let need = big_lam + big_mu - 2 * l_prev as i64;
                let lo = lo_i.max(need);
                if lo > hi_i {
                    continue;
                }
                // Range-add over L' = l_prev + v for v ∈ [lo, hi_i].
                let s = l_prev as i64 + lo;
                let e = l_prev as i64 + hi_i;
                if s > target {
                    continue;
                }
                let e = e.min(target);
                nxt[s as usize] += ways as i128;
                nxt[e as usize + 1] -= ways as i128;
            }
        }
        // Integrate the difference array back into `cur`.
        let mut run: i128 = 0;
        for l in 0..span {
            run += nxt[l];
            cur[l] = run.max(0) as u128;
        }
        big_mu += i64::from(mu.get(i).copied().unwrap_or(0));
    }
    cur[target as usize]
}

/// (B) Iterate candidate λ, count each fibre. No chain is built.
fn by_counting(mu: &[u32], nu: &[u32]) -> (Vec<(Vec<u32>, u128)>, u64) {
    let mu_size: u32 = mu.iter().sum();
    let total = mu_size + nu[0] + nu[1];
    let rows = mu.len() + 2;
    let mut out = Vec::new();
    let mut work = 0u64;
    let mut cand = vec![0u32; rows];
    let span = (mu_size + nu[0]) as usize + 2;
    let mut cur = vec![0u128; span];
    let mut nxt = vec![0i128; span + 1];

    #[allow(clippy::too_many_arguments)]
    fn rec(
        j: usize,
        rows: usize,
        left: u32,
        prev: u32,
        cand: &mut Vec<u32>,
        mu: &[u32],
        nu: &[u32],
        out: &mut Vec<(Vec<u32>, u128)>,
        work: &mut u64,
        cur: &mut [u128],
        nxt: &mut [i128],
    ) {
        let at = |i: usize| mu.get(i).copied().unwrap_or(0);
        if j == rows {
            if left != 0 {
                return;
            }
            *work += 1;
            let mut lam: Vec<u32> = cand.clone();
            while lam.last() == Some(&0) {
                lam.pop();
            }
            let c = count_fibre(mu, &lam, nu[0], cur, nxt);
            if c > 0 {
                out.push((lam, c));
            }
            return;
        }
        let lo = at(j);
        let hi = if j < 2 {
            prev.min(left + lo)
        } else {
            at(j - 2).min(prev).min(left + lo)
        };
        if hi < lo {
            return;
        }
        for v in lo..=hi {
            cand[j] = v;
            rec(j + 1, rows, left - (v - lo), v, cand, mu, nu, out, work, cur, nxt);
        }
        cand[j] = 0;
    }
    rec(0, rows, total - mu_size, u32::MAX, &mut cand, mu, nu, &mut out, &mut work, &mut cur, &mut nxt);
    out.sort();
    (out, work)
}

fn main() {
    for (mu_s, nu_s) in [
        ("6,4,2", "6,4"),
        ("20,16,12", "20,16"),
        ("30,24,18", "30,24"),
        ("40,32,24", "40,32"),
        ("60,48,36", "60,48"),
    ] {
        let mu: Vec<u32> = mu_s.split(',').map(|x| x.parse().unwrap()).collect();
        let nu: Vec<u32> = nu_s.split(',').map(|x| x.parse().unwrap()).collect();

        symfn::clear_caches();
        let t_lib = Instant::now();
        let _ = AutoLr.schur_product(
            &Partition::new(mu.iter().copied()),
            &Partition::new(nu.iter().copied()),
        );
        let t_lib = t_lib.elapsed();

        let truth: Vec<(Vec<u32>, u128)> = {
            let mut v: Vec<_> = AutoLr
                .schur_product(&Partition::new(mu.iter().copied()), &Partition::new(nu.iter().copied()))
                .into_iter()
                .map(|(l, c)| (l.parts().to_vec(), c))
                .collect();
            v.sort();
            v
        };

        let big = mu[0] >= 60;
        let t = Instant::now();
        let (ca, wa) = if big { (Vec::new(), 0) } else { by_chains(&mu, &nu) };
        let ta = t.elapsed();

        let t = Instant::now();
        let (cb, wb) = by_counting(&mu, &nu);
        let tb = t.elapsed();

        p!("s[{mu_s}]·s[{nu_s}]:  {} terms", truth.len());
        p!("   AutoLr   {:>12} {:>19.4?}", "(frontier)", t_lib);
        if !big {
            p!("   chains   {:>12} built   {:>11.4?}   {}", wa, ta, if ca == truth { "ok" } else { "MISMATCH" });
        }
        p!(
            "   counting {:>12} tested  {:>11.4?}   {}   -> {:>7}, {:.2}x vs AutoLr, {:.2}µs/term",
            wb, tb,
            if cb == truth { "ok" } else { "MISMATCH" },
            if big { "--".to_string() } else { format!("{:.0}x work", wa as f64 / wb.max(1) as f64) },
            t_lib.as_secs_f64() / tb.as_secs_f64().max(1e-12),
            tb.as_secs_f64() * 1e6 / truth.len().max(1) as f64
        );
        if cb != truth {
            p!("      got {} terms, want {}", cb.len(), truth.len());
            for (k, v) in cb.iter().take(3) {
                p!("      got {k:?} -> {v}");
            }
            for (k, v) in truth.iter().take(3) {
                p!("      want {k:?} -> {v}");
            }
        }
    }
}
