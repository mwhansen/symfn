//! The Macdonald operator ladder: ∇e_n, Δ'_{e_k}e_n, Θ.
//!
//! ```text
//!   cargo run --release --example bench_deltaop -- 12
//! ```
//!
//! The Sage side of the same ladder is in `docs/spec-macdonald-operators.md` §2
//! (29.165s at n=10, 81.506s at n=11, 234.631s at n=12).

use std::time::Instant;

use symfn::coeff::Rational;
use symfn::sym::SymFn;
use symfn::{Partition, QtPoly, Schur};

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(11);

    println!("== nabla e_n: closed form over ℤ, over ℚ, and the general pairing ==");
    println!(
        "{:>3} {:>6} {:>12} {:>12} {:>12} {:>8}",
        "n", "p(n)", "i128(s)", "Rational(s)", "pairing(s)", "terms"
    );
    for n in 1..=top {
        symfn::clear_caches();
        let t0 = Instant::now();
        let z = symfn::nabla_e::<i128>(n);
        let integral = t0.elapsed().as_secs_f64();

        symfn::clear_caches();
        let t0 = Instant::now();
        let a = symfn::nabla_e::<Rational>(n);
        let closed = t0.elapsed().as_secs_f64();

        // the two widths must agree, or `i128` wrapped
        for (lambda, c) in z.terms() {
            for (&k, &v) in c.terms() {
                assert_eq!(
                    Rational::from_int(v),
                    a.coeff(lambda).coeff(k.0, k.1),
                    "i128 disagrees with Rational at degree {n}, {lambda}"
                );
            }
        }

        symfn::clear_caches();
        let e: Schur<QtPoly<Rational>> = Schur::monomial(
            Partition::new(std::iter::repeat(1).take(n as usize)),
            <QtPoly<Rational> as symfn::Ring>::one(),
        );
        let t0 = Instant::now();
        let b = symfn::nabla(&e);
        let pairing = t0.elapsed().as_secs_f64();

        assert_eq!(a, b, "the two routes must agree at degree {n}");
        println!(
            "{n:>3} {:>6} {integral:>12.4} {closed:>12.4} {pairing:>12.4} {:>8}",
            symfn::partitions_of(n).len(),
            a.terms().len()
        );
    }

    println!();
    println!("== Delta'_{{e_k}} e_n, the Delta conjecture ladder ==");
    println!("{:>3} {:>4} {:>12} {:>8}", "n", "k", "sec", "terms");
    for n in 1..=top {
        for k in [0, n / 2, n.saturating_sub(1)] {
            if k >= n {
                continue;
            }
            symfn::clear_caches();
            let t0 = Instant::now();
            let f = symfn::delta_prime_e::<i128>(k, n);
            let dt = t0.elapsed().as_secs_f64();
            println!("{n:>3} {k:>4} {dt:>12.4} {:>8}", f.terms().len());
        }
    }

    println!();
    println!("== Theta_{{e_k}} nabla e_{{n-k}}  (degree n, expansion at n) ==");
    println!("{:>3} {:>4} {:>12} {:>8}", "n", "k", "sec", "terms");
    for n in 2..=top.min(10) {
        for k in [1, n / 2] {
            if k >= n {
                continue;
            }
            symfn::clear_caches();
            let inner = symfn::nabla_e::<Rational>(n - k);
            let t0 = Instant::now();
            let f = symfn::theta(&symfn::deltaop::elementary(k), &inner);
            let dt = t0.elapsed().as_secs_f64();
            println!("{n:>3} {k:>4} {dt:>12.4} {:>8}", f.terms().len());
        }
    }
}
