//! Time the inverse expansions `s → H̃` and `s → J`, one degree per process.
//!
//! ```text
//!   cargo run --release --example bench_inverse -- 6
//! ```
//!
//! Prints `<workload> <seconds>` for one degree and exits. The degree is the
//! unit and the process is the boundary: both this crate and Sage memoize
//! across degrees, and `docs/record/qt-kostka.md` measures what a sequential
//! run does to the numbers. `scripts/bench_inverse.py` drives one process per
//! degree and per arm.
//!
//! The workload is **every** λ of the degree, which is what a caller asking
//! whether some family of functions is `H̃`- or `J`-positive does.
//! `j_table` is the same answer from the whole-degree entry point, which is
//! where a caller with every shape of one degree should go.

use std::time::Instant;

use symfn::{
    clear_caches, partitions_of, schur_in_j_table, schur_to_macdonald_ht, schur_to_macdonald_j,
    QtPoly, Rational, Schur, SymFn,
};

fn s_lambda(lambda: &symfn::Partition) -> Schur<QtPoly<Rational>> {
    Schur::monomial(lambda.clone(), QtPoly::term(0, 0, Rational::from_int(1)))
}

fn main() {
    let n: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);

    let jobs: Vec<(&str, fn(u32) -> usize)> = vec![
        ("ht", |n| {
            partitions_of(n)
                .iter()
                .map(|l| schur_to_macdonald_ht(&s_lambda(l)).len())
                .sum()
        }),
        ("j", |n| {
            partitions_of(n)
                .iter()
                .map(|l| schur_to_macdonald_j(&s_lambda(l)).len())
                .sum()
        }),
        ("j_table", |n| {
            schur_in_j_table::<Rational>(n)
                .iter()
                .flatten()
                .filter(|c| !symfn::Ring::is_zero(*c))
                .count()
        }),
    ];

    for (name, f) in jobs {
        // Cold, so no workload is charged for another's tables.
        clear_caches();
        let t0 = Instant::now();
        let terms = f(n);
        println!("{name} {:.4} {terms}", t0.elapsed().as_secs_f64());
    }
}
