//! Model problem for "count per output instead of enumerating chains".
//!
//! Computes s_μ·h_a·h_b two ways and compares work:
//!
//!   (A) chain enumeration — for every λ¹ = μ + horizontal a-strip, for every
//!       λ² = λ¹ + horizontal b-strip, accumulate. This is what the frontier DP
//!       does, minus the lattice condition.
//!
//!   (B) per-output counting — for each candidate λ², count the λ¹ directly.
//!       Interlacing pins λ¹ᵢ to [max(μᵢ, λ²ᵢ₊₁), min(μᵢ₋₁, λ²ᵢ)] independently,
//!       with Σλ¹ fixed, so the coefficient is the number of lattice points in a
//!       box on a hyperplane — a bounded-composition count, closed form by
//!       inclusion–exclusion. No chain is ever built.
//!
//! Both are checked against the library. h_r = s_(r), so the truth is two
//! Schur products.

use std::collections::HashMap;
use std::io::Write;
use std::time::Instant;

use symfn::{AutoLr, LrBackend, Partition};

macro_rules! p {
    ($($a:tt)*) => {{ println!($($a)*); std::io::stdout().flush().ok(); }};
}

/// Every λ' with λ'/λ a horizontal strip of size `r`.
fn h_strips(lambda: &[u32], r: u32) -> Vec<Vec<u32>> {
    let rows = lambda.len() + 1;
    let mut out = Vec::new();
    let mut cur = vec![0u32; rows];
    fn rec(
        j: usize,
        rows: usize,
        left: u32,
        lam: &[u32],
        cur: &mut Vec<u32>,
        out: &mut Vec<Vec<u32>>,
    ) {
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
        // Interlacing: row j may not exceed the previous row's *old* value.
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

/// (A) Build every chain.
fn by_chains(mu: &[u32], a: u32, b: u32) -> (Vec<(Vec<u32>, u128)>, u64) {
    let mut acc: HashMap<Vec<u32>, u128> = HashMap::new();
    let mut work = 0u64;
    for l1 in h_strips(mu, a) {
        for l2 in h_strips(&l1, b) {
            work += 1;
            *acc.entry(l2).or_insert(0) += 1;
        }
    }
    let mut v: Vec<_> = acc.into_iter().collect();
    v.sort();
    (v, work)
}

/// Number of integer x with 0 ≤ xᵢ ≤ capᵢ and Σxᵢ = m, by inclusion–exclusion
/// over which coordinates break their cap.
fn bounded_compositions(cap: &[u32], m: i64) -> u128 {
    if m < 0 {
        return 0;
    }
    let n = cap.len();
    if n == 0 {
        return u128::from(m == 0);
    }
    let mut total: i128 = 0;
    for mask in 0u32..(1 << n) {
        let mut rem = m;
        for i in 0..n {
            if mask >> i & 1 == 1 {
                rem -= i64::from(cap[i]) + 1;
            }
        }
        if rem < 0 {
            continue;
        }
        // C(rem + n - 1, n - 1)
        let mut c: i128 = 1;
        for k in 1..n {
            c = c * (rem as i128 + k as i128) / k as i128;
        }
        if mask.count_ones() % 2 == 0 {
            total += c;
        } else {
            total -= c;
        }
    }
    total.max(0) as u128
}

/// (B) For each candidate λ², count λ¹ in closed form.
fn by_counting(mu: &[u32], a: u32, b: u32) -> (Vec<(Vec<u32>, u128)>, u64) {
    let total: u32 = mu.iter().sum::<u32>() + a + b;
    let rows = mu.len() + 2; // two strips can open at most two new rows
    let at = |i: usize| mu.get(i).copied().unwrap_or(0); // 0-based μ
    let mut out = Vec::new();
    let mut work = 0u64;
    let mut cand = vec![0u32; rows];

    // λ² ⊇ μ, and λ²ⱼ ≤ μⱼ₋₂ because λ²ⱼ ≤ λ¹ⱼ₋₁ ≤ μⱼ₋₂ (0-based: ≤ at(j-2)).
    #[allow(clippy::too_many_arguments)]
    fn rec(
        j: usize,
        rows: usize,
        left: u32,
        prev: u32,
        cand: &mut Vec<u32>,
        mu: &[u32],
        a: u32,
        out: &mut Vec<(Vec<u32>, u128)>,
        work: &mut u64,
    ) {
        let at = |i: usize| mu.get(i).copied().unwrap_or(0);
        if j == rows {
            if left != 0 {
                return;
            }
            *work += 1;
            let mut l2: Vec<u32> = cand.clone();
            while l2.last() == Some(&0) {
                l2.pop();
            }
            // λ¹ᵢ ∈ [max(μᵢ, λ²ᵢ₊₁), min(μᵢ₋₁, λ²ᵢ)], independently.
            let n1 = l2.len().max(mu.len());
            let (mut lo_sum, mut cap) = (0i64, Vec::with_capacity(n1));
            for i in 0..n1 {
                let l2i = l2.get(i).copied().unwrap_or(0);
                let l2i1 = l2.get(i + 1).copied().unwrap_or(0);
                let lo = at(i).max(l2i1);
                let hi = if i == 0 { l2i } else { at(i - 1).min(l2i) };
                if hi < lo {
                    return; // empty box: coefficient 0
                }
                lo_sum += i64::from(lo);
                cap.push(hi - lo);
            }
            let target = mu.iter().sum::<u32>() as i64 + i64::from(a) - lo_sum;
            let c = bounded_compositions(&cap, target);
            if c > 0 {
                out.push((l2, c));
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
            rec(j + 1, rows, left - (v - lo), v, cand, mu, a, out, work);
        }
        cand[j] = 0;
    }
    let mu_size: u32 = mu.iter().sum();
    rec(
        0,
        rows,
        total - mu_size,
        u32::MAX,
        &mut cand,
        mu,
        a,
        &mut out,
        &mut work,
    );
    let _ = at;
    out.sort();
    (out, work)
}

fn main() {
    for (mu_s, a, b) in [
        ("6,4,2", 6u32, 4u32),
        ("20,16,12", 20, 16),
        ("30,24,18", 30, 24),
    ] {
        let mu: Vec<u32> = mu_s.split(',').map(|x| x.parse().unwrap()).collect();
        let mup = Partition::new(mu.iter().copied());

        // Ground truth: h_r = s_(r), so two Schur products.
        let step1 = AutoLr.schur_product(&mup, &Partition::new([a]));
        let mut truth: HashMap<Vec<u32>, u128> = HashMap::new();
        for (l1, c1) in &step1 {
            for (l2, c2) in AutoLr.schur_product(l1, &Partition::new([b])) {
                *truth.entry(l2.parts().to_vec()).or_insert(0) += c1 * c2;
            }
        }
        let mut truth: Vec<_> = truth.into_iter().collect();
        truth.sort();

        let t = Instant::now();
        let (ca, wa) = by_chains(&mu, a, b);
        let ta = t.elapsed();

        let t = Instant::now();
        let (cb, wb) = by_counting(&mu, a, b);
        let tb = t.elapsed();

        p!("s[{mu_s}]·h_{a}·h_{b}:  {} terms", truth.len());
        p!(
            "   chains   {:>12} built   {:>10.4?}   {}",
            wa,
            ta,
            if ca == truth { "ok" } else { "MISMATCH" }
        );
        p!(
            "   counting {:>12} tested  {:>10.4?}   {}   -> {:.1}x work, {:.2}x time",
            wb,
            tb,
            if cb == truth { "ok" } else { "MISMATCH" },
            wa as f64 / wb.max(1) as f64,
            ta.as_secs_f64() / tb.as_secs_f64().max(1e-12)
        );
    }
}
