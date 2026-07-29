//! Check both versions of the Delta conjecture as far as the enumeration runs.
//!
//! ```text
//!   cargo run --release --example delta_conjecture -- 9        # whole symmetric function
//!   cargo run --release --example delta_conjecture -- 10 h1n   # the <·, h_1^n> slice only
//! ```
//!
//! `Δ'_{e_k} e_n` comes from [`symfn::delta_prime_e`]; the two right-hand sides
//! come from [`symfn::dyck`]. The **rise** version is a theorem, so a mismatch
//! there is a bug in this crate. The **valley** version is open — a mismatch
//! there, with the rise version still agreeing, would be a counterexample and is
//! printed as one rather than swallowed.
//!
//! One enumeration serves the whole `k` ladder ([`symfn::ladder`]), so the time
//! reported per `n` covers **all** of it, not one `(n,k)`.

use std::time::Instant;

use symfn::coeff::Rational;
use symfn::convert::ToSchur;
use symfn::sym::SymFn;
use symfn::{Partition, QtPoly, Ring, Side};

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);
    let h1n = std::env::args().nth(2).is_some_and(|s| s == "h1n");

    println!(
        "Delta conjecture, {} — Delta'_{{e_k}} e_n against labelled Dyck paths",
        if h1n {
            "coefficient of x_1...x_n"
        } else {
            "whole symmetric function"
        }
    );
    println!(
        "{:>3} {:>10} {:>11} {:>11}   {}",
        "n", "operator", "rise(s)", "valley(s)", "verdict (all k)"
    );

    let mut bad = 0;
    for n in 1..=top {
        let t0 = Instant::now();
        let full: Vec<_> = (0..n)
            .map(|k| symfn::delta_prime_e::<Rational>(k, n))
            .collect();
        // <f, h_1^n> = Σ_λ coeff · f^λ, the coefficient of x_1⋯x_n.
        let slice: Vec<QtPoly<Rational>> = full
            .iter()
            .map(|f| {
                let mut acc: QtPoly<Rational> = QtPoly::zero();
                for (lambda, c) in f.terms() {
                    let d = symfn::dimension(lambda).expect("a partition has a dimension");
                    for (&(a, b), v) in c.terms() {
                        acc.add_term(a, b, v.mul(&Rational::from_u128(d as u128)));
                    }
                }
                acc
            })
            .collect();
        let op = t0.elapsed().as_secs_f64();

        let mut secs = [0.0f64; 2];
        let mut good = [0usize; 2];
        for (slot, which) in [Side::Rise, Side::Valley].into_iter().enumerate() {
            let t0 = Instant::now();
            let (got_slice, got_full) = if h1n {
                let ones = Partition::new(std::iter::repeat(1).take(n as usize));
                (
                    symfn::ladder_at_content::<Rational>(&ones, which),
                    Vec::new(),
                )
            } else {
                (Vec::new(), symfn::ladder::<Rational>(n, which))
            };
            secs[slot] = t0.elapsed().as_secs_f64();
            for k in 0..n as usize {
                let ok = if h1n {
                    got_slice[k] == slice[k]
                } else {
                    got_full[k].to_schur() == full[k]
                };
                if ok {
                    good[slot] += 1;
                } else {
                    bad += 1;
                    println!("    MISMATCH n={n} k={k} {which:?}");
                }
            }
        }
        println!(
            "{n:>3} {op:>10.3} {:>11.3} {:>11.3}   rise {}/{n}, valley {}/{n}",
            secs[0], secs[1], good[0], good[1]
        );
    }

    println!();
    if bad == 0 {
        println!("both versions agree with the operator for every k checked.");
        println!("the rise version is a theorem; the valley version is open, and this");
        println!("is evidence for it at every (n,k) above.");
    } else {
        println!("{bad} MISMATCH(ES) — if only the valley column disagrees, check Val(P)");
        println!("before believing it (see src/dyck.rs module docs), then report it.");
    }
}
