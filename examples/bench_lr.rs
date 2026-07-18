//! Compare the LR backends head-to-head (run with --release).
use std::time::Instant;
use symfn::{clear_caches, LrBackend, NaiveLr, Partition, SkewLr, StripLr};

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
        "{:<18} {:>11} {:>11} {:>11} {:>9}",
        "shape^2", "NaiveLr", "StripLr", "SkewLr", "vs best"
    );
    for sh in &shapes {
        let p = Partition::new(sh.iter().copied());
        let (tn, a) = time(|| NaiveLr.schur_product(&p, &p));
        let (ts, b) = time(|| StripLr.schur_product(&p, &p));
        let (tk, c) = time(|| SkewLr.schur_product(&p, &p));
        assert_eq!(a, b, "NaiveLr vs StripLr on {p}");
        assert_eq!(a, c, "NaiveLr vs SkewLr on {p}");
        println!(
            "{:<18} {:>9.4}s {:>9.4}s {:>9.4}s {:>8.1}x",
            format!("{p}"),
            tn,
            ts,
            tk,
            tn.min(ts) / tk
        );
    }

    // The skew primitive, where the win is structural rather than constant:
    // expanding s_{λ/μ} used to run one full backtrack per candidate content.
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
                    let c = NaiveLr.lr_coeff(&outer, &inner, &nu);
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
