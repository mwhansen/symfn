//! The two routes to a whole (q,t)-Kostka table, on the same work.
//!
//! ```text
//!   cargo run --release --example bench_qtk_routes -- 10
//! ```
//!
//! Like for like, unlike `bench_llm.rs`: both produce the same `p(n)²`
//! polynomials in the same basis, and a test asserts they agree.

use std::time::Instant;

use symfn::coeff::Rational;
use symfn::QtPoly;

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);
    println!(
        "{:>3} {:>6} {:>12} {:>12} {:>9}",
        "n", "p(n)", "branching", "operator", "ratio"
    );
    for n in 1..=top {
        symfn::clear_caches();
        let t0 = Instant::now();
        let a: Vec<Vec<QtPoly<Rational>>> = symfn::qt_kostka_table(n);
        let branching = t0.elapsed().as_secs_f64();

        symfn::clear_caches();
        let t0 = Instant::now();
        let b: Vec<Vec<QtPoly<Rational>>> = symfn::qt_kostka_table_via_operator(n);
        let operator = t0.elapsed().as_secs_f64();

        assert_eq!(a, b, "the two routes must agree at degree {n}");
        let ratio = if operator > 0.0 { branching / operator } else { 0.0 };
        println!(
            "{n:>3} {:>6} {branching:>11.4}s {operator:>11.4}s {ratio:>8.1}x",
            a.len()
        );
    }
}
