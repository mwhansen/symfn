//! The LLT ladder: whole-degree ribbon tables, tuples, the by-path shuffle
//! refinement, and Kazhdan–Lusztig columns.
//!
//! ```text
//!   cargo run --release --example bench_llt -- 12
//! ```
//!
//! The Sage side is `docs/spec-llt.md` §2.2, measured on this machine. Whole
//! degree `HSp_k[μ] → s`, fresh parent per (k, n):
//!
//! ```text
//!   k \ n      7        8        9        10       11
//!    2       0.65 s   2.22 s  10.48 s   51.16 s   >120 s
//!    3       2.11 s  13.01 s  90.01 s   >120 s      —
//!    4      18.72 s   >120 s     —         —        —
//! ```
//!
//! and 3-tuples: `((2,2),(2,2),(2,1))` (n = 11) 19.69 s,
//! `((3,2),(2,2),(2,1))` (n = 12) >120 s.
//!
//! ⚠️ **Those Sage rows are battery numbers** and this machine drifts about
//! 1.8× at low charge, so a ratio below ~2× means nothing. Run both sides on
//! mains before quoting anything.
//!
//! Every ribbon row runs the two-width ladder — `i64` against `i128` — and
//! asserts they agree term for term, so a wrapped fixed-width intermediate is a
//! failure rather than a quietly different number.

use std::time::Instant;

use symfn::convert::ToSchur;
use symfn::llt::{
    llt_g, llt_gtilde_table, llt_h, llt_h_table, llt_kl_column, nabla_e_by_path, SkewTuple,
};
use symfn::qt::QtPoly;
use symfn::sym::{Monomial, Schur, SymFn};
use symfn::Partition;

fn part(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

/// The two coefficient widths must agree, or `i64` wrapped somewhere.
fn assert_widths_agree(
    narrow: &[(Partition, Monomial<QtPoly<i64>>)],
    wide: &[(Partition, Monomial<QtPoly<i128>>)],
    what: &str,
) {
    assert_eq!(narrow.len(), wide.len(), "{what}: table sizes");
    for ((la, a), (lb, b)) in narrow.iter().zip(wide.iter()) {
        assert_eq!(la, lb, "{what}: shape order");
        assert_eq!(a.terms().len(), b.terms().len(), "{what}: {la} support");
        for (mu, p) in a.terms() {
            let w = b.coeff(mu);
            for (&(x, y), c) in p.terms() {
                assert_eq!(
                    i128::from(*c),
                    w.coeff(x, y),
                    "{what}: i64 and i128 disagree at {la} / {mu}"
                );
            }
        }
    }
}

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(11);

    println!("== whole-degree H^(k) tables: the unit Sage's walls are measured in ==");
    println!(
        "{:>3} {:>3} {:>6} {:>11} {:>11} {:>9} {:>10}",
        "k", "n", "p(n)", "i64(s)", "i128(s)", "terms", "growth"
    );
    for k in 2..=4u32 {
        let mut prev = 0f64;
        for n in 1..=top {
            symfn::clear_caches();
            let t0 = Instant::now();
            let narrow: Vec<(Partition, Monomial<QtPoly<i64>>)> = llt_h_table(n, k);
            let secs = t0.elapsed().as_secs_f64();

            let t0 = Instant::now();
            let wide: Vec<(Partition, Monomial<QtPoly<i128>>)> = llt_h_table(n, k);
            let wsecs = t0.elapsed().as_secs_f64();
            assert_widths_agree(&narrow, &wide, &format!("H^({k}) table at n={n}"));

            let terms: usize = narrow.iter().map(|(_, f)| f.terms().len()).sum();
            let growth = if prev > 0.0 {
                format!("{:>10.2}", secs / prev)
            } else {
                format!("{:>10}", "—")
            };
            prev = secs;
            println!(
                "{k:>3} {n:>3} {:>6} {secs:>11.4} {wsecs:>11.4} {terms:>9} {growth}",
                narrow.len()
            );
        }
        println!();
    }

    println!("== the single shape μ = (n), 94-100% of a whole degree for Sage ==");
    println!(
        "{:>3} {:>3} {:>11} {:>9}",
        "k", "n", "H^(k)_(n) (s)", "terms"
    );
    for k in 2..=4u32 {
        for n in 1..=top {
            symfn::clear_caches();
            let t0 = Instant::now();
            let f: Monomial<QtPoly<i64>> = llt_h(&part(&[n]), k);
            println!(
                "{k:>3} {n:>3} {:>13.5} {:>9}",
                t0.elapsed().as_secs_f64(),
                f.terms().len()
            );
        }
    }

    println!();
    println!("== every G̃^(k)_λ of a degree from one walk (no Sage entry point) ==");
    println!("{:>3} {:>3} {:>9} {:>12}", "k", "n", "shapes", "table(s)");
    for k in 2..=3u32 {
        for n in 1..=top.min(9) {
            symfn::clear_caches();
            let t0 = Instant::now();
            let table: Vec<(Partition, Monomial<QtPoly<i64>>)> = llt_gtilde_table(n, k);
            println!(
                "{k:>3} {n:>3} {:>9} {:>12.4}",
                table.len(),
                t0.elapsed().as_secs_f64()
            );
        }
    }

    println!();
    println!("== tuples via R1 (SYT + descent buckets); Sage dies on the last row ==");
    println!(
        "{:>26} {:>4} {:>12} {:>9}",
        "tuple", "n", "G_nu(s)", "terms"
    );
    let tuples: &[&[&[u32]]] = &[
        &[&[2, 2], &[2, 1], &[2]],
        &[&[2, 2], &[2, 2], &[2]],
        &[&[3, 2], &[2, 2], &[1]],
        &[&[2, 2], &[2, 2], &[2, 1]],
        &[&[3, 2], &[2, 2], &[2, 1]],
    ];
    for shapes in tuples {
        let ps: Vec<Partition> = shapes.iter().map(|s| part(s)).collect();
        let nu = SkewTuple::from_partitions(&ps, &vec![0i32; ps.len()]);
        let t0 = Instant::now();
        let f: Monomial<QtPoly<i64>> = llt_g(&nu);
        let label: Vec<String> = ps.iter().map(|p| p.to_string()).collect();
        println!(
            "{:>26} {:>4} {:>12.4} {:>9}",
            label.join(""),
            nu.size(),
            t0.elapsed().as_secs_f64(),
            f.terms().len()
        );
    }

    println!();
    println!("== the by-path shuffle refinement: ∇e_n as Σ_D t^area G_D, per path ==");
    println!(
        "{:>3} {:>8} {:>12} {:>12}  cross-check",
        "n", "paths", "by-path(s)", "to Schur(s)"
    );
    for n in 1..=top.min(10) {
        let t0 = Instant::now();
        let pieces = nabla_e_by_path::<i128>(n);
        let build = t0.elapsed().as_secs_f64();
        let t0 = Instant::now();
        let mut total: Schur<QtPoly<i128>> = Schur::zero();
        let mut negative = 0u64;
        for (_, g) in &pieces {
            let s = g.to_schur();
            for (_, p) in s.terms() {
                for (_, c) in p.terms() {
                    // Schur positivity per path: a theorem for tuples of
                    // partitions, an unpublished preprint in general. A negative
                    // coefficient is a result to report, not a bug to fix.
                    if *c < 0 {
                        negative += 1;
                    }
                }
            }
            total = total.add(&s);
        }
        let conv = t0.elapsed().as_secs_f64();
        let want = symfn::deltaop::nabla_e::<i128>(n);
        let verdict = if total != want {
            "*** DISAGREES WITH deltaop::nabla_e ***".to_string()
        } else if negative > 0 {
            format!("*** {negative} NEGATIVE SCHUR COEFFS -- REPORT ***")
        } else {
            "matches nabla_e, all Schur-positive".to_string()
        };
        println!(
            "{n:>3} {:>8} {build:>12.4} {conv:>12.4}  {verdict}",
            pieces.len()
        );
    }

    println!();
    println!("== [LLT] Conj 6.4: is H^(k+1)_mu − H^(k)_mu Schur-positive?  OPEN ==");
    println!(
        "{:>3} {:>3} {:>6} {:>10} {:>10} {:>12}  verdict",
        "k", "n", "p(n)", "sweep(s)", "coeffs", "negatives"
    );
    for k in 1..=4u32 {
        for n in 1..=top {
            let t0 = Instant::now();
            let lo = llt_h_table::<i128>(n, k);
            let hi = llt_h_table::<i128>(n, k + 1);
            let (mut coeffs, mut negative) = (0u64, 0u64);
            let mut witness = String::new();
            for ((mu, a), (mu2, b)) in lo.iter().zip(hi.iter()) {
                assert_eq!(mu, mu2, "the two tables must be in shape order");
                let diff = b.to_schur().sub(&a.to_schur());
                for (nu, p) in diff.terms() {
                    for (&(e, _), c) in p.terms() {
                        coeffs += 1;
                        // OPEN CONJECTURE. A violation is a result to report,
                        // not a bug to fix.
                        if *c < 0 {
                            negative += 1;
                            if witness.is_empty() {
                                witness = format!(" first: mu={mu} s_{nu} q^{e} -> {c}");
                            }
                        }
                    }
                }
            }
            let verdict = if negative == 0 {
                "Schur-positive".to_string()
            } else {
                format!("*** COUNTEREXAMPLE -- REPORT ***{witness}")
            };
            println!(
                "{k:>3} {n:>3} {:>6} {:>10.4} {coeffs:>10} {negative:>12}  {verdict}",
                lo.len(),
                t0.elapsed().as_secs_f64()
            );
        }
        println!();
    }

    println!("== Kazhdan–Lusztig columns via [KMS] straightening (no incumbent) ==");
    println!(
        "{:>3} {:>10} {:>6} {:>12}",
        "k", "lambda", "shapes", "column(s)"
    );
    for k in 2..=3u32 {
        for d in 1..=top.min(6) {
            let t0 = Instant::now();
            let mut entries = 0usize;
            for lambda in symfn::partitions_of(d) {
                entries += llt_kl_column::<i64>(&lambda, k).len();
            }
            println!(
                "{k:>3} {:>10} {entries:>6} {:>12.4}",
                format!("all |-{d}"),
                t0.elapsed().as_secs_f64()
            );
        }
    }
}
