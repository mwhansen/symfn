//! The two routes to `J_λ`, on the same work.
//!
//! ```text
//!   cargo run --release --example bench_llm -- 11
//! ```
//!
//! **Not yet like for like.** The branching formula returns `J_λ` in the
//! monomial basis; the eigenvector returns it in `S_μ[X^{tq}]`, and the crossing
//! between them is not written yet. So this compares the part that is comparable
//! — the work that produces the coefficients — and the crossing will be charged
//! against the eigenvector column when it exists. Read it as an upper bound on
//! how good the new route can be, not as a result.
//!
//! Both are timed per degree over *every* shape, which is the unit a
//! (q,t)-Kostka table needs.

use std::time::Instant;

use symfn::coeff::Rational;
use symfn::sym::SymFn;
use symfn::{Frac, Monomial, QtPoly};

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!(
        "{:>3} {:>6} {:>12} {:>12} {:>9}",
        "n", "p(n)", "branching", "eigenvector", "ratio"
    );
    for n in 1..=top {
        symfn::clear_caches();
        let start = Instant::now();
        let mut sink = 0usize;
        for lambda in symfn::partitions_of(n) {
            let j: Monomial<Frac<Rational>> = symfn::macdonald_j(&lambda);
            sink += j.terms().len();
        }
        let branching = start.elapsed().as_secs_f64();

        symfn::clear_caches();
        let start = Instant::now();
        let all: Vec<(_, Vec<QtPoly<Rational>>, _)> = symfn::eigenvectors(n);
        for (_, b, _) in &all {
            sink += b.iter().filter(|p| !p.is_empty()).count();
        }
        let eigen = start.elapsed().as_secs_f64();

        let ratio = if eigen > 0.0 { branching / eigen } else { 0.0 };
        println!(
            "{n:>3} {:>6} {branching:>11.4}s {eigen:>11.4}s {ratio:>8.1}x  ({sink})",
            all.len()
        );
    }
}
