//! An independent check on a whole Schur expansion, including shapes no
//! external program can reach.
//!
//! `[24,20,16,12]²` is the one case in the comparison sweep with **no external
//! oracle** — lrcalc cannot finish it, so its 5.3M terms have only ever been
//! checked against our own conjugate orientation. That is a real consistency
//! check but not an independent one: a bug in the shared frontier code would
//! reproduce itself in both orientations.
//!
//! Principal specialization gives an independent one. Evaluating
//! `s_μ · s_ν = Σ_λ c^λ_{μν} s_λ` at `x = (1,1,…,1)` (n ones) turns the identity
//! into a scalar equation:
//!
//! ```text
//!   s_μ(1ⁿ) · s_ν(1ⁿ)  =  Σ_λ c^λ_{μν} · s_λ(1ⁿ)
//! ```
//!
//! and each `s_λ(1ⁿ)` follows from the hook-content formula
//! `∏_{(i,j)∈λ} (n + j − i) / h(i,j)`, which shares no code with the LR
//! machinery — it needs only the diagram and its conjugate. The left side never
//! touches a coefficient at all. So agreement tests every coefficient at once
//! against something the frontier had no hand in.
//!
//! It is a *weighted checksum*, not a proof: distinct wrong expansions could
//! collide. Each extra n is another independent equation, and the weights
//! `s_λ(1ⁿ)` vary wildly with λ, so a handful of n makes a coincidental pass
//! vanishingly unlikely — while a single wrong coefficient fails almost surely.
//!
//! Arithmetic is mod a prime so the checksum stays in `u64`: the true values
//! here have thousands of digits.
//!
//!   cargo run --release --example verify_specialization            # small shapes
//!   cargo run --release --example verify_specialization 24,20,16,12  # the big one

use std::io::Write;
use std::time::Instant;

use symfn::{AutoLr, LrBackend, Partition};

macro_rules! p {
    ($($a:tt)*) => {{ println!($($a)*); std::io::stdout().flush().ok(); }};
}

const P: u64 = (1 << 61) - 1; // Mersenne prime, so products fit in u128

fn mul(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % P as u128) as u64
}

fn pow(mut a: u64, mut e: u64) -> u64 {
    let mut r = 1;
    while e > 0 {
        if e & 1 == 1 {
            r = mul(r, a);
        }
        a = mul(a, a);
        e >>= 1;
    }
    r
}

/// `s_λ(1ⁿ) mod P` by the hook-content formula.
///
/// Zero when `ℓ(λ) > n`, which the formula produces on its own: the cell
/// `(n, 1)` then has content `1 − n − 1`, making the numerator `n + j − i = 0`.
fn schur_at_ones(lambda: &Partition, n: u64) -> u64 {
    let conj = lambda.conjugate();
    let mut num = 1u64;
    let mut den = 1u64;
    for (i, &row) in lambda.parts().iter().enumerate() {
        for j in 0..row as usize {
            // Hook: arm + leg + 1, using the conjugate for the column length.
            let hook = (row as usize - j) + (conj.part(j) as usize - i) - 1;
            // Content j − i, offset by n. Both are small, so this stays exact
            // before reduction; a true zero (ℓ(λ) > n) survives as zero.
            let content = n as i64 + j as i64 - i as i64;
            if content <= 0 {
                return 0;
            }
            num = mul(num, content as u64 % P);
            den = mul(den, hook as u64 % P);
        }
    }
    mul(num, pow(den, P - 2))
}

/// A checksum that cannot fail proves nothing. Perturb one coefficient in an
/// otherwise correct expansion and confirm the check rejects it — for every
/// term in turn, so this also shows no position is a blind spot.
fn negative_control() {
    let mu = Partition::new([4u32, 3, 2, 1]);
    let prod = AutoLr.schur_product(&mu, &mu);
    let base = 2 * mu.len() as u64;
    let ns: Vec<u64> = [base + 1, base + 3, base + 8, base + 21, base + 55].into();
    let lhs: Vec<u64> = ns.iter().map(|&n| mul(schur_at_ones(&mu, n), schur_at_ones(&mu, n))).collect();

    let (mut caught, mut total) = (0usize, 0usize);
    for k in 0..prod.len() {
        for delta in [1i128, -1] {
            if prod[k].1 as i128 + delta < 0 {
                continue;
            }
            total += 1;
            let detected = ns.iter().zip(&lhs).any(|(&n, &want)| {
                let mut rhs = 0u64;
                for (i, (lambda, c)) in prod.iter().enumerate() {
                    let c = if i == k { (*c as i128 + delta) as u128 } else { *c };
                    rhs = (rhs + mul((c % P as u128) as u64, schur_at_ones(lambda, n))) % P;
                }
                rhs != want
            });
            if detected {
                caught += 1;
            }
        }
    }
    p!(
        "negative control: {caught}/{total} single-coefficient perturbations detected  {}",
        if caught == total { "OK" } else { "<-- BLIND SPOT" }
    );
}

fn main() {
    let arg = std::env::args().nth(1);
    if arg.is_none() {
        negative_control();
    }
    let shapes: Vec<Vec<u32>> = match &arg {
        Some(s) => vec![s.split(',').map(|x| x.trim().parse().expect("part")).collect()],
        None => vec![
            vec![2, 1],
            vec![3, 2, 1],
            vec![5, 4, 3, 2, 1],
            vec![6, 5, 4, 3, 2],
            vec![8, 6, 4],
            vec![12, 10, 8],
            vec![7, 7, 7, 7],
        ],
    };

    for sh in shapes {
        let mu = Partition::new(sh.iter().copied());
        let t = Instant::now();
        let prod = AutoLr.schur_product(&mu, &mu);
        let elapsed = t.elapsed();

        // n must exceed the longest λ or every weight is zero and the check is
        // vacuous; ℓ(λ) ≤ 2·ℓ(μ) here.
        let base = 2 * mu.len() as u64;
        let mut all_ok = true;
        let mut checked = 0;
        for n in [base + 1, base + 3, base + 8, base + 21, base + 55] {
            let lhs = mul(schur_at_ones(&mu, n), schur_at_ones(&mu, n));
            let mut rhs = 0u64;
            for (lambda, c) in &prod {
                rhs = (rhs + mul((c % P as u128) as u64, schur_at_ones(lambda, n))) % P;
            }
            if lhs == 0 && rhs == 0 {
                continue; // vacuous, do not count it as evidence
            }
            checked += 1;
            if lhs != rhs {
                all_ok = false;
                p!("  n={n}: MISMATCH lhs={lhs} rhs={rhs}");
            }
        }
        p!(
            "s{mu}²  {:>9} terms  product {:>9.3?}  {} independent n  {}",
            prod.len(),
            elapsed,
            checked,
            if all_ok && checked > 0 { "OK" } else { "FAILED" }
        );
    }
}
