//! Baselines for every operation that is **not** Littlewood–Richardson.
//!
//! LR has three harnesses (`bench_lr`, `bench_shapes`, `compare_lrcalc.py`);
//! everything else in the library had none. That asymmetry is not harmless:
//! closing the equivalent hole in the LR comparison immediately exposed a 21x
//! defect in `lr_coeff` that no existing benchmark could have shown. This gives
//! the rest of the library the same visibility.
//!
//! There is no external oracle here the way lrcalc is for LR, so these are
//! *self-relative*: they answer "did this commit make plethysm slower?", not
//! "are we fast?". Correctness lives in `tests/sage_oracle.rs`.
//!
//! ```text
//!   cargo build --release --example bench_ops
//!   cp target/release/examples/bench_ops /tmp/after     # and likewise before
//!   for i in 1 2 3; do /tmp/before before; /tmp/after after; done
//! ```
//!
//! **Interleave builds and take the min per (build, case).** Two lessons paid
//! for elsewhere in this project: an A/B that always runs one side second
//! measures allocator warmth as much as algorithm, and a laptop on battery
//! throttles progressively, so later rows of a long run are pessimistic. Every
//! case clears the caches first, so none is warmed by an earlier one — the
//! memoized behavior is a separate row rather than a contaminant.
//!
//! Columns: tag, operation, seconds, and a work unit (what was computed), so a
//! time can be read against how much it did.

use std::time::Instant;

use symfn::{
    antipode, character, clear_caches, convert, coproduct, hall, kostka, omega, partitions_of,
    plethysm, skew_schur, Elementary, Homogeneous, Monomial, Partition, PowerSum, Rational, Schur,
    SymFn,
};

fn p(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

/// Time `f` once from cold caches, printing `tag  name  seconds  work`.
fn bench<T>(tag: &str, name: &str, work: impl Fn(&T) -> String, f: impl FnOnce() -> T) {
    clear_caches();
    let t = Instant::now();
    let out = f();
    let dt = t.elapsed();
    println!("{tag}\t{name}\t{:.6}\t{}", dt.as_secs_f64(), work(&out));
}

fn count<T>(n: usize) -> impl Fn(&T) -> String {
    move |_| format!("{n} values")
}

fn main() {
    let tag = std::env::args().nth(1).unwrap_or_else(|| "bench".into());

    // --- Kostka numbers: the s → m transition -------------------------------
    for n in [12u32, 16, 20] {
        let parts = partitions_of(n);
        let k = parts.len();
        bench(
            &tag,
            &format!("kostka_all_pairs_n{n}"),
            count(k * k),
            || {
                let mut acc = 0u128;
                for a in &parts {
                    for b in &parts {
                        acc = acc.wrapping_add(kostka(a, b));
                    }
                }
                acc
            },
        );
    }

    // The whole table in one sweep, against the pairwise rows above: same
    // p(n)² values, one trie of Pieri steps instead of p(n)² chain DPs.
    for n in [20u32, 24] {
        let k = partitions_of(n).len();
        bench(&tag, &format!("kostka_table_n{n}"), count(k * k), || {
            symfn::kostka::kostka_table(n)
                .iter()
                .flatten()
                .fold(0u128, |a, &x| a.wrapping_add(x))
        });
    }

    // --- Characters: Murnaghan–Nakayama -------------------------------------
    for n in [16u32, 18] {
        let parts = partitions_of(n);
        let k = parts.len();
        bench(&tag, &format!("character_table_n{n}"), count(k * k), || {
            let mut acc = 0i128;
            for a in &parts {
                for b in &parts {
                    acc = acc.wrapping_add(character(a, b));
                }
            }
            acc
        });
    }

    // The whole table is a different engine from the pairwise loop above: one
    // β-mask sweep through `p_expand`, not p(n)² independent MN recursions.
    for n in [24u32, 28] {
        let k = partitions_of(n).len();
        bench(
            &tag,
            &format!("character_beta_sweep_n{n}"),
            count(k * k),
            || {
                clear_caches();
                symfn::character::character_table(n)
                    .iter()
                    .flatten()
                    .fold(0i128, |a, &x| a.wrapping_add(x))
            },
        );
    }

    // --- Basis conversions ---------------------------------------------------
    // Integral bases stay over ℤ; anything through p needs ℚ, which is the
    // expensive family and the reason it gets its own rows.
    // Degree 20. This was degree 12 while `kostka` was exponential — at 20 a
    // single `convert_s_to_m` took 400 seconds. It now takes 12 ms, so the
    // sizes here move up to stay informative.
    let lam = p(&[6, 5, 4, 3, 2]);
    let s_i: Schur<i128> = Schur::monomial(lam.clone(), 1);
    let s_q: Schur<Rational> = Schur::monomial(lam.clone(), Rational::new(1, 1));

    bench(
        &tag,
        "convert_s_to_m",
        |m: &Monomial<i128>| format!("{} terms", m.terms().len()),
        || convert::<i128, _, Monomial<i128>>(&s_i),
    );
    bench(
        &tag,
        "convert_s_to_e",
        |e: &Elementary<i128>| format!("{} terms", e.terms().len()),
        || convert::<i128, _, Elementary<i128>>(&s_i),
    );
    bench(
        &tag,
        "convert_s_to_h",
        |h: &Homogeneous<i128>| format!("{} terms", h.terms().len()),
        || convert::<i128, _, Homogeneous<i128>>(&s_i),
    );
    bench(
        &tag,
        "convert_s_to_p",
        |v: &PowerSum<Rational>| format!("{} terms", v.terms().len()),
        || convert::<Rational, _, PowerSum<Rational>>(&s_q),
    );

    // The reverse directions: m → s inverts the Kostka matrix, p → s uses
    // characters, so these are the ones that can blow up.
    let m_i: Monomial<i128> = convert::<i128, _, Monomial<i128>>(&s_i);
    let m_in = m_i.terms().len();
    bench(
        &tag,
        "convert_m_to_s",
        move |x: &Schur<i128>| format!("{m_in} m-terms in, {} out", x.terms().len()),
        || convert::<i128, _, Schur<i128>>(&m_i),
    );
    let p_q: PowerSum<Rational> = convert::<Rational, _, PowerSum<Rational>>(&s_q);
    let p_in = p_q.terms().len();
    bench(
        &tag,
        "convert_p_to_s",
        move |x: &Schur<Rational>| format!("{p_in} p-terms in, {} out", x.terms().len()),
        || convert::<Rational, _, Schur<Rational>>(&p_q),
    );

    // --- ω and the Hall inner product ---------------------------------------
    bench(
        &tag,
        "omega_on_h",
        |h: &Homogeneous<i128>| format!("{} terms", h.terms().len()),
        || {
            let h: Homogeneous<i128> = convert::<i128, _, Homogeneous<i128>>(&s_i);
            omega::<i128, Homogeneous<i128>>(&h)
        },
    );
    let h_in = p_q.terms().len();
    bench(
        &tag,
        "hall_s_p_degree17",
        move |_: &Rational| format!("1 pairing over {h_in} p-terms"),
        || hall::<Rational, _, _>(&s_q, &p_q),
    );

    // --- Hopf structure ------------------------------------------------------
    bench(
        &tag,
        "skew_schur_[9,8,7,6,5]/[3,2,1]",
        |x: &Schur<i128>| format!("{} terms", x.terms().len()),
        || skew_schur::<i128>(&p(&[9, 8, 7, 6, 5]), &p(&[3, 2, 1])),
    );
    let s_big: Schur<i128> = Schur::monomial(p(&[8, 7, 6, 5, 4]), 1);
    bench(
        &tag,
        "coproduct_[8,7,6,5,4]",
        |t: &symfn::SymTensor<i128>| format!("{} terms", t.terms().len()),
        || coproduct(&s_big),
    );
    bench(
        &tag,
        "antipode_[8,7,6,5,4]",
        |x: &Schur<i128>| format!("{} terms", x.terms().len()),
        || antipode(&s_big),
    );

    // --- Plethysm ------------------------------------------------------------
    // Routed through the power-sum basis, so it pays rational arithmetic and
    // is the most expensive single operation in the library.
    for (f, g) in [
        (&[2, 1][..], &[2, 1][..]),
        (&[4][..], &[3][..]),
        (&[3, 2][..], &[2, 1][..]),
    ] {
        let (fs, gs) = (p(f), p(g));
        let name = format!("plethysm_{fs}[{gs}]");
        bench(
            &tag,
            &name,
            |x: &Schur<Rational>| format!("{} terms", x.terms().len()),
            || {
                let a: Schur<Rational> = Schur::monomial(fs.clone(), Rational::new(1, 1));
                let b: Schur<Rational> = Schur::monomial(gs.clone(), Rational::new(1, 1));
                plethysm(&a, &b)
            },
        );
    }

    // --- Where s→m dies ------------------------------------------------------
    // `kostka` enumerates SSYT one at a time, so a single value is exponential
    // in the shape. This row exists to keep that visible and to show any future
    // fix immediately.
    // Explicit multi-row shapes: a single-row λ has K_{λμ} = 1 for every μ, so
    // a shape list that degenerates to one row measures nothing at all.
    for parts in [&[5u32, 4, 3, 2, 1][..], &[6, 5, 4, 3, 2], &[8, 7, 6, 5, 4]] {
        let lam = p(parts);
        let mus = partitions_of(lam.size());
        let n = mus.len();
        bench(&tag, &format!("kostka_row_{lam}"), count(n), || {
            mus.iter()
                .map(|m| kostka(&lam, m))
                .fold(0u128, u128::wrapping_add)
        });
    }

    // --- Memoization ---------------------------------------------------------
    // The warm path is the real Sage usage pattern (many queries against one
    // object) and was never measured. Deliberately does NOT clear caches
    // between the repeats, which is the entire point of the row.
    clear_caches();
    let parts = partitions_of(12);
    let _ = parts.iter().map(|a| kostka(a, a)).sum::<u128>(); // warm
    let t = Instant::now();
    let mut acc = 0u128;
    for _ in 0..50 {
        for a in &parts {
            acc = acc.wrapping_add(kostka(a, a));
        }
    }
    println!(
        "{tag}\tkostka_warm_50x\t{:.6}\t{} lookups ({acc})",
        t.elapsed().as_secs_f64(),
        50 * parts.len()
    );
}
