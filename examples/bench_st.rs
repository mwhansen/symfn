//! The reduced Kronecker product, against the wall it exists to break.
//!
//! Unlike `bench_ops`, this one *does* have an external reference: the same
//! cases were timed in Sage on 2026-07-29 and are recorded in
//! `docs/record/kronecker.md`. Three of them Sage cannot finish at all,
//! so the comparison is deliberately asymmetric — for those rows the only
//! meaningful statement is "completes, in this long".
//!
//! ```text
//!   cargo run --release --example bench_st
//! ```
//!
//! Every case clears the caches first, so nothing is warmed by an earlier row;
//! the memoized behaviour is the last section rather than a contaminant. The
//! bit-width column answers the open question: the reduced
//! Kronecker coefficients are non-negative but have no `√(n!)`-style bound of
//! the kind that justified `i128` for `K̃`, so their size is measured and not
//! assumed. `peak_bits` covers the *answer*; the intermediate rationals inside
//! Γ⁻¹ are a separate question that a coefficient survey cannot see.

use std::time::Instant;

use symfn::{clear_caches, reduced_kronecker_product, Partition, St, SymFn};

fn p(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

fn bits(n: i128) -> u32 {
    128 - n.unsigned_abs().leading_zeros()
}

fn row(lambda: &[u32], mu: &[u32], sage: &str) {
    clear_caches();
    let t = Instant::now();
    let out: St<i128> = reduced_kronecker_product(&p(lambda), &p(mu));
    let dt = t.elapsed().as_secs_f64();
    let peak = out.terms().values().map(|c| bits(*c)).max().unwrap_or(0);
    let biggest = out.terms().values().copied().max().unwrap_or(0);
    println!(
        "  {:<24} {:>7} {:>10.4} {:>12} {:>10} {:>9}",
        format!("st{lambda:?} · st{mu:?}"),
        out.terms().len(),
        dt,
        sage,
        biggest,
        peak
    );
}

fn main() {
    std::panic::set_hook(Box::new(|_| {})); // the wall section expects a panic
    println!("reduced Kronecker product — cold caches per row");
    println!(
        "  {:<24} {:>7} {:>10} {:>12} {:>10} {:>9}",
        "case", "terms", "sec", "sage(sec)", "max coeff", "bits"
    );
    row(&[2, 1], &[2, 1], "0.0221");
    row(&[3, 1], &[3, 1], "0.0577");
    row(&[3, 2], &[3, 2], "1.4418");
    row(&[4, 2], &[4, 2], "8.6438");
    row(&[4, 3], &[4, 3], "30.2569");
    row(&[5, 3], &[5, 3], ">90");
    row(&[6, 4], &[6, 4], ">90");
    row(&[8, 5], &[7, 4], ">90");
    row(&[2, 1, 1], &[2, 1, 1], "0.0108");
    row(&[3, 2, 1], &[3, 2, 1], "0.1467");
    row(&[4, 3, 2], &[4, 3, 2], "-");
    row(&[1, 1, 1, 1, 1], &[1, 1, 1, 1, 1], "-");

    println!();
    println!("degree ladder in |λ|+|μ| — the curve, not the point");
    println!(
        "  {:<24} {:>7} {:>10} {:>12} {:>10} {:>9}",
        "case", "terms", "sec", "sage(sec)", "max coeff", "bits"
    );
    for n in 2..=8u32 {
        row(&[n, n / 2], &[n, n / 2], "-");
    }

    println!();
    println!("the fixed-width wall (build with --features bignum to pass it)");
    println!(
        "  {:<24} {:>7} {:>10} {:>12} {:>10} {:>9}",
        "case", "terms", "sec", "sage(sec)", "max coeff", "bits"
    );
    for (a, b) in [
        (vec![7u32, 4], vec![7u32, 4]),
        (vec![8, 5], vec![7, 4]),
        (vec![8, 5], vec![8, 5]),
        (vec![9, 6], vec![9, 6]),
    ] {
        // Without `bignum` the escalation has nowhere to go and panics, which is
        // the designed behaviour: a wrapped intermediate must never be returned.
        let caught = std::panic::catch_unwind(|| {
            clear_caches();
            let t = Instant::now();
            let out: St<i128> = reduced_kronecker_product(&p(&a), &p(&b));
            (
                out.terms().len(),
                t.elapsed().as_secs_f64(),
                out.terms().values().copied().max().unwrap_or(0),
            )
        });
        match caught {
            Ok((n, dt, mx)) => println!(
                "  {:<24} {:>7} {:>10.4} {:>12} {:>10} {:>9}",
                format!("st{a:?} · st{b:?}"),
                n,
                dt,
                ">90",
                mx,
                bits(mx)
            ),
            Err(_) => println!(
                "  {:<24} {:>7} {:>10} {:>12} {:>10} {:>9}",
                format!("st{a:?} · st{b:?}"),
                "-",
                "overflow",
                ">90",
                "-",
                "-"
            ),
        }
    }

    println!();
    println!("warm: a whole table of pairs, caches shared");
    for n in [4u32, 5, 6] {
        clear_caches();
        let parts = symfn::partitions_of(n);
        let t = Instant::now();
        let mut total = 0usize;
        for a in &parts {
            for b in &parts {
                let out: St<i128> = reduced_kronecker_product(a, b);
                total += out.terms().len();
            }
        }
        println!(
            "  degree {n}: {:>4} pairs {:>10.4}s  {:>7} terms",
            parts.len() * parts.len(),
            t.elapsed().as_secs_f64(),
            total
        );
    }
}
