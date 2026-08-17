//! Compare the LR backends head-to-head (run with --release).
//!
//! ```text
//!   cargo run --release --example bench_lr
//! ```
//!
//! **The baseline is [`StripLr`], not `NaiveLr`.** The naive backend is the
//! reference oracle — one independent backtrack per candidate content — and on
//! the shapes worth timing it does not finish: the harness ran past ten
//! minutes on `[8,7,6,5,4,3]²` and had to be killed, which makes it useless
//! for the thing a benchmark is for. Its job is to disagree with a production
//! route, and that job is done in `tests/lrcalc_oracle.rs` and
//! `tests/lr_specialization.rs`, against a shape size chosen for an oracle
//! rather than for a stopwatch. A benchmark that cannot be run teaches
//! nothing, and correctness asserted here would be asserted in the one place
//! it is least affordable.
//!
//! So both columns are production backends and the comparison is the one a
//! change can actually move. The equality assertions stay: two backends
//! disagreeing is still a bug, and checking it costs nothing once both have
//! run.
use std::time::Instant;
use symfn::{clear_caches, LrBackend, Partition, SkewLr, StripLr};

fn time(f: impl FnOnce() -> Vec<(Partition, u128)>) -> (f64, Vec<(Partition, u128)>) {
    clear_caches();
    let t = Instant::now();
    let v = f();
    (t.elapsed().as_secs_f64(), v)
}

fn main() {
    let shapes: Vec<Vec<u32>> = vec![
        vec![5, 4, 3, 2, 1],
        vec![6, 5, 4, 3, 2],
        vec![6, 5, 4, 3, 2, 1],
        vec![7, 6, 5, 4, 3],
        vec![8, 7, 6, 5, 4, 3],
    ];
    println!(
        "{:<18} {:>11} {:>11} {:>9}",
        "shape^2", "StripLr", "SkewLr", "speedup"
    );
    for sh in &shapes {
        let p = Partition::new(sh.iter().copied());
        let (ts, a) = time(|| StripLr.schur_product(&p, &p));
        let (tk, b) = time(|| SkewLr.schur_product(&p, &p));
        assert_eq!(a, b, "StripLr vs SkewLr on {p}");
        println!(
            "{:<18} {:>9.4}s {:>9.4}s {:>8.1}x",
            format!("{p}"),
            ts,
            tk,
            ts / tk
        );
    }

    // The skew primitive, where the win is structural rather than constant:
    // expanding s_{λ/μ} costs one query per candidate content on a backend that
    // answers one coefficient at a time, and one traversal on this one. The
    // per-nu side runs through StripLr for the same reason as above — the
    // structural factor is p(n) either way, and the naive backend only makes
    // the constant in front of it unaffordable.
    println!("\nskew expansion s_(lambda/mu)");
    println!(
        "{:<26} {:>11} {:>11} {:>9}",
        "lambda/mu", "per-nu", "SkewLr", "speedup"
    );
    for (o, i) in [
        (vec![6u32, 5, 4, 3, 2], vec![2u32, 1]),
        (vec![7, 6, 5, 4, 3], vec![3, 2, 1]),
        (vec![8, 7, 6, 5, 4], vec![3, 2, 1]),
        // Wide shapes with a small inner part: many cells free to move, so the
        // per-nu sweep has both a large p(n) and expensive individual searches.
        (vec![9, 9, 8, 8, 7], vec![2, 1]),
        (vec![10, 9, 8, 7, 6, 5], vec![4, 3, 2, 1]),
    ] {
        let outer = Partition::new(o.iter().copied());
        let inner = Partition::new(i.iter().copied());
        let n = outer.size() - inner.size();
        let (told, a) = time(|| {
            let mut v: Vec<(Partition, u128)> = symfn::partitions_of(n)
                .into_iter()
                .filter_map(|nu| {
                    let c = StripLr.lr_coeff(&outer, &inner, &nu);
                    (c != 0).then_some((nu, c))
                })
                .collect();
            v.sort_by(|x, y| x.0.cmp(&y.0));
            v
        });
        let (tnew, b) = time(|| symfn::expand_skew(&outer, &inner));
        assert_eq!(a, b, "skew {outer}/{inner}");
        println!(
            "{:<26} {:>9.4}s {:>9.4}s {:>8.1}x",
            format!("{outer}/{inner}"),
            told,
            tnew,
            told / tnew
        );
    }
}
