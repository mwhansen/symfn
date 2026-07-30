//! The §6 ladder of `docs/spec-schubert.md`: E3 against the two incumbents on
//! the exact cases of §2.
//!
//! Run with `--release`. Single-threaded, per the standing rule — the C
//! `schubmult` this is measured against is single-threaded, and a parallel
//! number set against it would conflate an algorithmic claim with a hardware
//! one.
//!
//! The reference columns are recorded, not re-run here:
//!   - **schubmult C** (Buch, lrcalc, in Sage's own tree) — re-measured
//!     2026-07-30 and reproducible to within 5% of the §2 values. ⚠️ Each of
//!     its timings carries ~3ms of process startup, so its rows at or below
//!     ~0.03s are largely measuring `exec` and are not real targets.
//!   - **Sage/Symmetrica** — what Sage ships today. `>120` means it did not
//!     finish inside the spec's per-case alarm.
//!
//! `states` is E3's cost model (distinct peel-DAG states of the peeled
//! factor); `pipe_dreams` is E1's and Symmetrica's. The gap between those two
//! columns is the whole thesis, and this harness exists to find out how much
//! of it converts into time.

use std::time::Instant;
use symfn::permutation::Perm;
use symfn::schubert::{dimension, peel_states, Schubert};

/// (label, u, v, schubmult C seconds, Sage/Symmetrica seconds) — `None` = did
/// not finish in 120s.
type Case = (&'static str, Vec<u32>, Vec<u32>, Option<f64>, Option<f64>);

fn stair(k: u32) -> Vec<u32> {
    (1..=k)
        .map(|i| 2 * i)
        .chain((1..=k).map(|i| 2 * i - 1))
        .collect()
}

fn cases() -> Vec<Case> {
    let mut v: Vec<Case> = vec![
        (
            "stair4^2 (S_8)",
            stair(4),
            stair(4),
            Some(0.005),
            Some(0.078),
        ),
        (
            "stair5^2 (S_10)",
            stair(5),
            stair(5),
            Some(0.016),
            Some(6.85),
        ),
        ("stair6^2 (S_12)", stair(6), stair(6), Some(0.355), None),
        ("stair7^2 (S_14)", stair(7), stair(7), Some(20.6), None),
    ];
    // seed-1 random pairs, identical to spec_schubert_walls2.py / _schubmult.py
    let rand: Vec<Case> = vec![
        (
            "S_10.0 l=21,22",
            vec![6, 3, 7, 4, 1, 8, 10, 9, 5, 2],
            vec![5, 4, 3, 9, 2, 7, 8, 10, 6, 1],
            Some(0.066),
            Some(47.2),
        ),
        (
            "S_10.1 l=19,18",
            vec![10, 4, 5, 1, 6, 2, 7, 8, 3, 9],
            vec![9, 6, 3, 1, 4, 7, 2, 5, 8, 10],
            Some(0.0043),
            Some(0.0004),
        ),
        (
            "S_10.2 l=28,20",
            vec![8, 9, 3, 10, 4, 2, 5, 7, 1, 6],
            vec![8, 5, 2, 3, 6, 7, 4, 10, 1, 9],
            Some(0.0045),
            Some(0.0067),
        ),
        (
            "S_11.0 l=32,27",
            vec![1, 10, 7, 11, 9, 2, 4, 8, 5, 6, 3],
            vec![10, 4, 7, 5, 6, 2, 1, 9, 3, 11, 8],
            Some(0.222),
            None,
        ),
        (
            "S_11.1 l=35,26",
            vec![4, 8, 5, 11, 10, 6, 9, 1, 7, 3, 2],
            vec![7, 5, 3, 2, 9, 6, 11, 10, 1, 4, 8],
            Some(0.748),
            None,
        ),
        (
            "S_11.2 l=31,21",
            vec![8, 11, 5, 7, 1, 2, 4, 10, 6, 9, 3],
            vec![3, 1, 6, 8, 11, 5, 7, 4, 2, 10, 9],
            Some(3.86),
            None,
        ),
        (
            "S_12.0 l=37,24",
            vec![12, 1, 4, 6, 11, 3, 10, 8, 9, 5, 7, 2],
            vec![3, 2, 5, 4, 12, 6, 9, 11, 8, 7, 1, 10],
            Some(21.0),
            None,
        ),
        (
            "S_12.1 l=27,21",
            vec![1, 2, 4, 10, 8, 11, 12, 7, 9, 3, 5, 6],
            vec![4, 2, 8, 3, 5, 7, 10, 9, 12, 1, 6, 11],
            Some(1.61),
            None,
        ),
        (
            "S_12.2 l=43,32",
            vec![7, 12, 10, 3, 9, 5, 2, 8, 4, 11, 6, 1],
            vec![3, 4, 12, 2, 10, 11, 6, 5, 8, 1, 9, 7],
            Some(19.0),
            None,
        ),
        (
            "S_13.0 l=25,36",
            vec![1, 8, 2, 12, 9, 3, 4, 6, 11, 7, 5, 13, 10],
            vec![1, 12, 7, 3, 2, 11, 8, 6, 13, 10, 5, 9, 4],
            None,
            None,
        ),
        (
            "S_13.1 l=36,41",
            vec![9, 13, 3, 7, 2, 1, 4, 10, 11, 5, 12, 8, 6],
            vec![9, 5, 8, 10, 2, 12, 11, 4, 3, 6, 1, 13, 7],
            Some(20.9),
            None,
        ),
        (
            "S_13.2 l=40,41",
            vec![6, 7, 10, 1, 12, 11, 8, 3, 9, 2, 5, 13, 4],
            vec![2, 12, 9, 11, 4, 10, 6, 3, 1, 8, 7, 5, 13],
            Some(105.8),
            None,
        ),
    ];
    v.extend(rand);
    v
}

fn fmt(t: Option<f64>) -> String {
    match t {
        Some(x) => format!("{x:.4}s"),
        None => ">120s".to_string(),
    }
}

fn main() {
    let budget: f64 = std::env::var("BUDGET")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(120.0);

    println!(
        "{:<18} {:>7} {:>13} {:>10} {:>10} {:>9} {:>8}",
        "case", "terms", "pipe_dreams", "states", "symfn E3", "schubmult", "vs C"
    );
    for (label, u, v, c_sec, sage_sec) in cases() {
        let pu = Perm::new(u).unwrap();
        let pv = Perm::new(v).unwrap();
        let a: Schubert<i128> = Schubert::monomial(pu.clone(), 1);
        let b: Schubert<i128> = Schubert::monomial(pv.clone(), 1);

        // the two cost models, for the factor each engine would actually peel
        let pd = dimension(&pu).min(dimension(&pv));
        let st = peel_states(&pu).min(peel_states(&pv));

        let t = Instant::now();
        let prod = a.mul(&b);
        let dt = t.elapsed().as_secs_f64();

        let ratio = match c_sec {
            Some(c) => format!("{:.1}x", c / dt),
            None => "wins".to_string(),
        };
        println!(
            "{:<18} {:>7} {:>13} {:>10} {:>9.4}s {:>9} {:>8}",
            label,
            prod.terms().len(),
            pd,
            st,
            dt,
            fmt(c_sec),
            ratio
        );
        let _ = sage_sec;
        if dt > budget {
            println!("  (over budget {budget}s — stopping)");
            break;
        }
    }
}
