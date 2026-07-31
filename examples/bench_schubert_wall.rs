//! Where does E2 stop? The standing requirement: "Report where it stops. Find this
//! engine's wall (S₁₄? S₁₅?) and put it in docs/record/schubert.md next to the
//! others."
//!
//! The incumbent ladder ended at S₁₃ because that is where the *incumbents* stopped —
//! Sage/Symmetrica walls at S₁₀–S₁₁ and the C `schubmult` at S₁₃–S₁₄. This
//! walks past both. Deterministic xorshift so the pairs are reproducible
//! without a rand dependency.
//!
//! Also re-checks coefficient growth, which an early note had at "≤16" on the strength
//! of Symmetrica's easy rows and which is already known to be wrong (c = 130
//! at S₁₃).
//!
//! `CAP=<seconds>` bounds each case; `PAIRS=<n>` sets pairs per size.

use std::time::Instant;
use symfn::permutation::Perm;
use symfn::schubert::{schubert_monomial_mass_of, Schubert};

fn main() {
    let cap: f64 = std::env::var("CAP")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(120.0);
    let pairs: usize = std::env::var("PAIRS")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let mut seed: u64 = 0x9E37_79B9_7F4A_7C15;
    let mut rng = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    let mut rand_perm = |n: u32, r: &mut dyn FnMut() -> u64| {
        let mut v: Vec<u32> = (1..=n).collect();
        for i in (1..n as usize).rev() {
            let j = (r() % (i as u64 + 1)) as usize;
            v.swap(i, j);
        }
        Perm::new(v).unwrap()
    };

    println!(
        "{:<20} {:>10} {:>12} {:>10} {:>9}",
        "case", "mass", "terms", "max|c|", "time"
    );
    for n in [14u32, 15, 16, 17] {
        for i in 0..pairs {
            let u = rand_perm(n, &mut rng);
            let v = rand_perm(n, &mut rng);
            let label = format!("S_{n}.{i} l={},{}", u.length(), v.length());
            let mass = schubert_monomial_mass_of(&u, &v);
            let a: Schubert<i128> = Schubert::monomial(u, 1);
            let b: Schubert<i128> = Schubert::monomial(v, 1);
            let t = Instant::now();
            let p = a.mul_e2(&b);
            let dt = t.elapsed().as_secs_f64();
            let mx = p.terms().values().map(|c| c.abs()).max().unwrap_or(0);
            println!(
                "{label:<20} {mass:>10.2e} {:>12} {mx:>10} {dt:>8.3}s",
                p.terms().len()
            );
            if dt > cap {
                println!("  (over cap {cap}s — this size is the wall)");
                return;
            }
        }
    }
}
