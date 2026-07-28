//! Dump `Q'_λ` in the Schur basis for every λ up to a degree, for `scripts/check_hl.py`.
//!
//! One line per (λ, μ) pair with a nonzero coefficient:
//!
//! ```text
//!   lambda | mu | e0:c0,e1:c1,...     (exponent of t : coefficient)
//! ```

use symfn::sym::SymFn;
use symfn::{hall_littlewood_p_table, hall_littlewood_table, QtPoly, Schur};

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);
    // "P" selects the other Hall–Littlewood basis; the default is Q'.
    let p_basis = std::env::args().nth(2).is_some_and(|a| a == "P");
    for n in 1..=top {
        let table = if p_basis {
            hall_littlewood_p_table::<i64>(n)
        } else {
            hall_littlewood_table::<i64>(n)
        };
        for (lambda, hl) in table {
            let hl: Schur<QtPoly<i64>> = hl;
            for (mu, c) in hl.terms() {
                let body: Vec<String> = c
                    .terms()
                    .map(|((a, b), v)| {
                        assert_eq!(*a, 0, "Hall-Littlewood must not involve q");
                        format!("{b}:{v}")
                    })
                    .collect();
                println!(
                    "{} | {} | {}",
                    join(lambda.parts()),
                    join(mu.parts()),
                    body.join(",")
                );
            }
        }
    }
}

fn join(p: &[u32]) -> String {
    p.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}
