//! How wide do the eigenvector coefficients get, and where does a fixed width
//! stop being safe?
//!
//! ```text
//!   cargo run --release --example llm_coeff_sizes -- 12
//! ```
//!
//! `impl_ring_for_int` multiplies with a plain `*`, so a wrapped `i64` or `i128`
//! is **silent**. The only honest way to find the ceiling is to compute the same
//! thing in two widths and look for the first disagreement — the same argument
//! `mac_coeff_sizes.rs` makes for the branching formula.
//!
//! These are larger than that route's coefficients and grow faster (~6 bits per
//! degree against ~3.5) because `solve` clears denominators: the returned `b_κ`
//! is `a_κ · v`, not `a_κ`.

use symfn::QtPoly;

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(11);
    println!(
        "{:>3}  {:>20}  {:>5}  i64 agrees with i128?",
        "n", "widest |coeff|", "bits"
    );
    for n in 1..=top {
        let a: Vec<(_, Vec<QtPoly<i64>>, QtPoly<i64>)> = symfn::macop::eigenvectors(n);
        let b: Vec<(_, Vec<QtPoly<i128>>, QtPoly<i128>)> = symfn::macop::eigenvectors(n);
        let widest = b
            .iter()
            .flat_map(|(_, v, _)| v.iter())
            .flat_map(|p| p.terms().map(|(_, c)| c.abs()))
            .max()
            .unwrap_or(0);
        let agree = a.iter().zip(b.iter()).all(|((_, x, _), (_, y, _))| {
            x.iter().zip(y.iter()).all(|(p, q)| {
                p.terms()
                    .map(|(k, c)| (*k, *c as i128))
                    .eq(q.terms().map(|(k, c)| (*k, *c)))
            })
        });
        println!(
            "{n:>3}  {widest:>20}  {:>5}  {}",
            128 - widest.leading_zeros(),
            if agree { "yes" } else { "NO - i64 has wrapped" }
        );
    }
}
