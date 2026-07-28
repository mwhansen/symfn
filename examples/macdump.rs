//! Dump `P_λ(x;q,t)` in the monomial basis for `scripts/check_macdonald.py`.
//!
//! One line per (λ, μ) with a nonzero coefficient:
//!
//! ```text
//!   lambda | mu | numerator | denominator
//! ```
//!
//! Numerator and denominator are `a,b:c` terms for `c·qᵃtᵇ`, semicolon
//! separated. The denominator is **expanded** here purely so the comparison is
//! simple; the library keeps it factored.

use symfn::sym::SymFn;
use symfn::{macdonald_p, Frac, Monomial, Ring};

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    for n in 1..=top {
        for lambda in symfn::partitions_of(n) {
            let p: Monomial<Frac<i128>> = macdonald_p(&lambda);
            for (mu, c) in p.terms() {
                let (num, _) = c.parts();
                println!(
                    "{} | {} | {} | {}",
                    join(lambda.parts()),
                    join(mu.parts()),
                    poly(num),
                    poly(&c.denominator())
                );
            }
        }
    }
}

fn poly(p: &symfn::QtPoly<i128>) -> String {
    if p.is_zero() {
        return "0,0:0".into();
    }
    p.terms()
        .map(|((a, b), c)| format!("{a},{b}:{c}"))
        .collect::<Vec<_>>()
        .join(";")
}

fn join(p: &[u32]) -> String {
    p.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}
