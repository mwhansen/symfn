//! Time the inverse expansions `s → H̃`, `s → J`, `m → P`, `m → Q` and the
//! three Jack ones, one degree per process.
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
//! where a caller with every shape of one degree should go. `p`, `q` and the
//! three `jack_*` workloads take the monomial basis rather than the Schur
//! one, because that is what those expansions are written in.

use std::time::Instant;

use symfn::{
    clear_caches, monomial_to_jack_j, monomial_to_jack_p, monomial_to_jack_q,
    monomial_to_macdonald_p, monomial_to_macdonald_q, partitions_of, schur_in_j_table,
    schur_to_macdonald_ht, schur_to_macdonald_j, AFrac, Frac, Monomial, QtPoly, Rational, Ring,
    Schur, SymFn,
};

fn s_lambda(lambda: &symfn::Partition) -> Schur<QtPoly<Rational>> {
    Schur::monomial(lambda.clone(), QtPoly::term(0, 0, Rational::from_int(1)))
}

fn m_lambda(lambda: &symfn::Partition) -> Monomial<Frac<Rational>> {
    Monomial::monomial(lambda.clone(), <Frac<Rational> as Ring>::one())
}

fn m_alpha(lambda: &symfn::Partition) -> Monomial<AFrac<Rational>> {
    Monomial::monomial(lambda.clone(), <AFrac<Rational> as Ring>::one())
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
                .filter(|c| !Ring::is_zero(*c))
                .count()
        }),
        ("p", |n| {
            partitions_of(n)
                .iter()
                .map(|l| monomial_to_macdonald_p(&m_lambda(l)).len())
                .sum()
        }),
        ("q", |n| {
            partitions_of(n)
                .iter()
                .map(|l| monomial_to_macdonald_q(&m_lambda(l)).len())
                .sum()
        }),
        ("jack_p", |n| {
            partitions_of(n)
                .iter()
                .map(|l| monomial_to_jack_p(&m_alpha(l)).len())
                .sum()
        }),
        ("jack_q", |n| {
            partitions_of(n)
                .iter()
                .map(|l| monomial_to_jack_q(&m_alpha(l)).len())
                .sum()
        }),
        ("jack_j", |n| {
            partitions_of(n)
                .iter()
                .map(|l| monomial_to_jack_j(&m_alpha(l)).len())
                .sum()
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
