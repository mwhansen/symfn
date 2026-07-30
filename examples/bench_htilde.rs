//! The cost of the whole modified-Macdonald table for one degree.
//!
//! ```text
//!   cargo run --release --example bench_htilde -- 12
//! ```
//!
//! `bh::htilde_table` is the input every Macdonald *operator* needs: ∇, Δ_f and
//! Θ_f are all diagonal on `{H̃_μ}`, so a whole degree is the unit of work. The
//! comparable Sage cost is inside `Ht(e[n])` — see
//! `docs/record/macdonald-operators-spec.md` §2, where both sides are
//! tabulated.

use std::time::Instant;

use symfn::coeff::Rational;
use symfn::sym::SymFn;

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(12);
    println!(
        "{:>3} {:>6} {:>12} {:>10}",
        "n", "p(n)", "htilde(s)", "terms"
    );
    for n in 1..=top {
        symfn::clear_caches();
        let t0 = Instant::now();
        let table = symfn::bh::htilde_table::<Rational>(n);
        let dt = t0.elapsed().as_secs_f64();
        let terms: usize = table.iter().map(|(_, s)| s.terms().len()).sum();
        println!("{n:>3} {:>6} {dt:>12.4} {terms:>10}", table.len());
    }
}
