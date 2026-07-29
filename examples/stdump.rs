//! Dump character-basis data for `scripts/check_st.py`.
//!
//! Three sections, so that a wrong answer and a right answer in the wrong slot
//! are distinguishable — the transitions can be checked even where the product
//! cannot, and vice versa.
//!
//! ```text
//!   cargo run --release --example stdump -- 6 > /tmp/st.txt
//!   sage -python scripts/check_st.py /tmp/st.txt
//! ```
//!
//! Line formats:
//!
//! ```text
//!   st2s | lambda | nu:c,nu:c,...        s̃_λ in the Schur basis
//!   s2st | nu     | lambda:c,...         s_ν in the s̃ basis
//!   prod | lam;mu | nu:c,...             s̃_λ · s̃_μ  (reduced Kronecker)
//! ```

use symfn::sym::SymFn;
use symfn::{partitions_of, reduced_kronecker_product, Partition, Schur, St};
use symfn::{FromSchur, ToSchur};

fn join(parts: &[u32]) -> String {
    if parts.is_empty() {
        return "0".to_string();
    }
    parts
        .iter()
        .map(|p| p.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn body(terms: impl Iterator<Item = (Partition, i128)>) -> String {
    terms
        .map(|(p, c)| format!("{}:{}", join(p.parts()), c))
        .collect::<Vec<_>>()
        .join(" ")
}

fn main() {
    // Degree for the transitions; products run to `top_prod` on each side, since
    // a product of two degree-n inputs reaches degree 2n.
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    let top_prod: u32 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(4);

    for n in 0..=top {
        for lambda in partitions_of(n) {
            let x: St<i128> = St::monomial(lambda.clone(), 1);
            let s = x.to_schur();
            println!(
                "st2s | {} | {}",
                join(lambda.parts()),
                body(s.terms().iter().map(|(p, c)| (p.clone(), *c)))
            );

            let y: Schur<i128> = Schur::monomial(lambda.clone(), 1);
            let st: St<i128> = St::from_schur(&y);
            println!(
                "s2st | {} | {}",
                join(lambda.parts()),
                body(st.terms().iter().map(|(p, c)| (p.clone(), *c)))
            );
        }
    }

    for a in 0..=top_prod {
        for lambda in partitions_of(a) {
            for b in 0..=top_prod {
                for mu in partitions_of(b) {
                    let prod: St<i128> = reduced_kronecker_product(&lambda, &mu);
                    println!(
                        "prod | {};{} | {}",
                        join(lambda.parts()),
                        join(mu.parts()),
                        body(prod.terms().iter().map(|(p, c)| (p.clone(), *c)))
                    );
                }
            }
        }
    }

    // The cases the spec is actually about: two-row shapes at and past Sage's
    // wall. Emitted unconditionally so the comparison script can decide which
    // of them Sage is willing to answer.
    for (lambda, mu) in [
        (vec![3u32, 2], vec![3u32, 2]),
        (vec![4, 2], vec![4, 2]),
        (vec![4, 3], vec![4, 3]),
        (vec![5, 3], vec![5, 3]),
        (vec![3, 2, 1], vec![3, 2, 1]),
    ] {
        let l = Partition::new(lambda);
        let m = Partition::new(mu);
        let prod: St<i128> = reduced_kronecker_product(&l, &m);
        println!(
            "prod | {};{} | {}",
            join(l.parts()),
            join(m.parts()),
            body(prod.terms().iter().map(|(p, c)| (p.clone(), *c)))
        );
    }
}
