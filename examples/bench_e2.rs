//! E2 (memoized transition) against E3 (merged peel DAG), on the incumbent cases.
//!
//! Both engines cost (nodes) × (size of the running element), so the question
//! this answers is narrow and factual: **is the transition tree smaller than
//! the peel DAG, and does that show up in time?** E3 is ~4× behind the C
//! `schubmult` on staircases but 51× behind on `S_11.1`, whose output is
//! 118 822 terms — so the interesting row is the large-output one, not the
//! staircases.
//!
//! `E2 nodes` / `E2 passes` are counted; `E3 nodes` is `peel_states`. A pass
//! is one `mul_variable` over a whole element, which is the unit of work both
//! engines are really spending.
//!
//! Run with `--release`. `SKIP_E3=1` drops the E3 column for rows where it is
//! known to be slow.

use std::time::Instant;
use symfn::permutation::Perm;
use symfn::schubert::{peel_states, Schubert};

fn stair(k: u32) -> Vec<u32> {
    (1..=k)
        .map(|i| 2 * i)
        .chain((1..=k).map(|i| 2 * i - 1))
        .collect()
}

/// (label, u, v, schubmult C seconds)
fn cases() -> Vec<(String, Vec<u32>, Vec<u32>, Option<f64>)> {
    let mut v: Vec<(String, Vec<u32>, Vec<u32>, Option<f64>)> = (4..=7)
        .map(|k| (format!("stair{k}^2"), stair(k), stair(k), None))
        .collect();
    v[0].3 = Some(0.005);
    v[1].3 = Some(0.016);
    v[2].3 = Some(0.355);
    v[3].3 = Some(20.6);
    for (label, u, w, c) in [
        (
            "S_10.0 l=21,22",
            vec![6, 3, 7, 4, 1, 8, 10, 9, 5, 2],
            vec![5, 4, 3, 9, 2, 7, 8, 10, 6, 1],
            0.066,
        ),
        (
            "S_11.0 l=32,27",
            vec![1, 10, 7, 11, 9, 2, 4, 8, 5, 6, 3],
            vec![10, 4, 7, 5, 6, 2, 1, 9, 3, 11, 8],
            0.222,
        ),
        (
            "S_11.1 l=35,26",
            vec![4, 8, 5, 11, 10, 6, 9, 1, 7, 3, 2],
            vec![7, 5, 3, 2, 9, 6, 11, 10, 1, 4, 8],
            0.748,
        ),
        (
            "S_11.2 l=31,21",
            vec![8, 11, 5, 7, 1, 2, 4, 10, 6, 9, 3],
            vec![3, 1, 6, 8, 11, 5, 7, 4, 2, 10, 9],
            3.86,
        ),
        (
            "S_12.1 l=27,21",
            vec![1, 2, 4, 10, 8, 11, 12, 7, 9, 3, 5, 6],
            vec![4, 2, 8, 3, 5, 7, 10, 9, 12, 1, 6, 11],
            1.61,
        ),
        (
            "S_12.0 l=37,24",
            vec![12, 1, 4, 6, 11, 3, 10, 8, 9, 5, 7, 2],
            vec![3, 2, 5, 4, 12, 6, 9, 11, 8, 7, 1, 10],
            21.0,
        ),
        (
            "S_12.2 l=43,32",
            vec![7, 12, 10, 3, 9, 5, 2, 8, 4, 11, 6, 1],
            vec![3, 4, 12, 2, 10, 11, 6, 5, 8, 1, 9, 7],
            19.0,
        ),
        (
            "S_13.1 l=36,41",
            vec![9, 13, 3, 7, 2, 1, 4, 10, 11, 5, 12, 8, 6],
            vec![9, 5, 8, 10, 2, 12, 11, 4, 3, 6, 1, 13, 7],
            20.9,
        ),
        (
            "S_13.2 l=40,41",
            vec![6, 7, 10, 1, 12, 11, 8, 3, 9, 2, 5, 13, 4],
            vec![2, 12, 9, 11, 4, 10, 6, 3, 1, 8, 7, 5, 13],
            105.8,
        ),
    ] {
        v.push((label.to_string(), u, w, Some(c)));
    }
    v.push((
        "S_13.0 l=25,36 *".to_string(),
        vec![1, 8, 2, 12, 9, 3, 4, 6, 11, 7, 5, 13, 10],
        vec![1, 12, 7, 3, 2, 11, 8, 6, 13, 10, 5, 9, 4],
        None, // schubmult C: >120s
    ));
    v.push(("stair8^2".to_string(), stair(8), stair(8), None));
    v
}

fn main() {
    let skip_e3 = std::env::var("SKIP_E3").is_ok();
    println!(
        "{:<16} {:>7} {:>9} {:>10} {:>10} {:>10} {:>10} {:>8}",
        "case", "terms", "E3 nodes", "E2 nodes", "E2 passes", "E3 time", "E2 time", "E2 vs C"
    );
    let only = std::env::var("ONLY").unwrap_or_default();
    for (label, u, v, c_sec) in cases() {
        if !only.is_empty() && !label.contains(&only) {
            continue;
        }
        let pu = Perm::new(u).unwrap();
        let pv = Perm::new(v).unwrap();
        let a: Schubert<i128> = Schubert::monomial(pu, 1);
        let b: Schubert<i128> = Schubert::monomial(pv, 1);
        let e3n = peel_states(&pu).min(peel_states(&pv));

        // E2 recurses on the second factor; try both orders and keep the
        // better, since c^w_{uv} = c^w_{vu} makes that free and the trees
        // differ (the `lr_coeff` lesson: peeling the wrong side cost 21x).
        let t = Instant::now();
        let (p1, n1, k1) = a.mul_e2_depth(&b, 4096);
        let d1 = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let (p2, n2, k2) = b.mul_e2_depth(&a, 4096);
        let d2 = t.elapsed().as_secs_f64();
        assert_eq!(p1, p2, "{label}: E2 disagreed with itself across orders");
        let (e2t, e2n, e2k) = if d1 <= d2 { (d1, n1, k1) } else { (d2, n2, k2) };

        let e3t = if skip_e3 {
            f64::NAN
        } else {
            let t = Instant::now();
            let p3 = a.mul_e3(&b);
            assert_eq!(p1, p3, "{label}: E2 vs E3 mismatch");
            t.elapsed().as_secs_f64()
        };

        println!(
            "{:<16} {:>7} {:>9} {:>10} {:>10} {:>9.4}s {:>9.4}s {:>8}",
            label,
            p1.terms().len(),
            e3n,
            e2n,
            e2k,
            e3t,
            e2t,
            match c_sec {
                Some(c) => format!("{:.2}x", c / e2t),
                None => "-".to_string(),
            }
        );
    }
}
