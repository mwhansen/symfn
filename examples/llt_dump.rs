//! Dump the LLT families for `scripts/check_llt.py`.
//!
//! ```text
//!   cargo run --release --example llt_dump -- 5 > /tmp/llt.txt
//! ```
//!
//! One record per line, `kind | key | ... | body`, with the body a
//! `exponent:coefficient` list in the single variable each family uses:
//!
//! ```text
//!   hspin   | k | mu     | nu   | e:c,...      H^(k)_mu, m-basis        [LLT] (28)
//!   hcospin | k | mu     | nu   | e:c,...      H~^(k)_mu, m-basis       [LLT] (27)
//!   cospin  | k | lambda | nu   | e:c,...      G~^(k)_lambda, m-basis   [LLT] (26)
//!   tuple   | shapes    | nu    | e:c,...      G_nu, m-basis, RAW inv grading
//!   floor   | shapes    | m                    min inv of the same tuple
//!   quot    | k | lambda | shapes              the k-quotient, runner order
//! ```
//!
//! `tuple` is deliberately the **unfloored** grading: Sage floors its tuple
//! entry point, so the dictionary `sage.cospin(tuple) = q^{−min inv} G_ν` is
//! only testable if both halves cross the boundary separately. That is what the
//! `floor` records are for.

use symfn::llt::{llt_g, llt_gtilde, llt_h, llt_h_tilde, llt_min_inv, SkewTuple};
use symfn::qt::QtPoly;
use symfn::sym::{Monomial, SymFn};
use symfn::{partitions_of, Partition};

fn join(p: &[u32]) -> String {
    p.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

/// A `q`-only coefficient as `e:c` pairs. Asserts the `t` slot is unused, so a
/// stray `t` becomes a failure here rather than a silently dropped exponent.
fn body(c: &QtPoly<i64>) -> String {
    c.terms()
        .map(|((a, b), v)| {
            assert_eq!(*b, 0, "the LLT families live in q alone");
            format!("{a}:{v}")
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn emit(kind: &str, key: &str, f: &Monomial<QtPoly<i64>>) {
    for (nu, c) in f.terms() {
        println!("{kind} | {key} | {} | {}", join(nu.parts()), body(c));
    }
}

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(4);

    for k in 1..=4u32 {
        for n in 1..=top {
            for mu in partitions_of(n) {
                let key = format!("{k} | {}", join(mu.parts()));
                emit("hspin", &key, &llt_h::<i64>(&mu, k));
                emit("hcospin", &key, &llt_h_tilde::<i64>(&mu, k));
            }
        }
    }

    // Cospin on plain shapes, every empty-core λ of the small degrees.
    for k in 2..=4u32 {
        for r in 1..=top {
            for lambda in partitions_of(k * r) {
                if !lambda.has_empty_k_core(k) {
                    continue;
                }
                let key = format!("{k} | {}", join(lambda.parts()));
                emit("cospin", &key, &llt_gtilde::<i64>(&lambda, k));
                let quot: Vec<String> = lambda
                    .k_quotient(k)
                    .iter()
                    .map(|p| join(p.parts()))
                    .collect();
                println!("quot | {key} | {}", quot.join(";"));
            }
        }
    }

    // Tuples, straight and skew, in the raw grading plus their floor.
    let straight: &[&[&[u32]]] = &[
        &[&[1], &[1]],
        &[&[2], &[1]],
        &[&[1, 1], &[1]],
        &[&[2], &[2]],
        &[&[2, 1], &[1]],
        &[&[1], &[1], &[1]],
        &[&[2], &[1], &[1]],
        &[&[2, 1], &[2], &[1]],
        &[&[2, 2], &[2, 1]],
        &[&[1], &[1], &[1], &[1]],
        &[&[2, 1], &[1, 1], &[2]],
    ];
    for shapes in straight {
        let ps: Vec<Partition> = shapes
            .iter()
            .map(|s| Partition::new(s.iter().copied()))
            .collect();
        let nu = SkewTuple::from_partitions(&ps, &vec![0i32; ps.len()]);
        let key: Vec<String> = ps.iter().map(|p| join(p.parts())).collect();
        let key = key.join(";");
        emit("tuple", &key, &llt_g::<i64>(&nu));
        println!("floor | {key} | {}", llt_min_inv(&nu));
    }
}
