//! Time `s_shape^2` through [`SkewLr`] alone, one tab-separated line per case.
//!
//! `bench_lr` compares the backends against each other and `compare_lrcalc.py`
//! compares against lrcalc; neither answers "did this commit make SkewLr
//! faster?". This does, and it is built to be *interleaved*: keep a binary from
//! each side of a change and alternate them, because run-to-run noise here is
//! large enough (tens of percent) that two consecutive runs of different builds
//! prove nothing. Take the min per (build, case) and only believe differences
//! that survive several rounds.
//!
//!   cargo build --release --example bench_shapes
//!   cp target/release/examples/bench_shapes /tmp/after      # and likewise before
//!   for i in 1 2 3; do /tmp/before before; /tmp/after after; done
//!
//! The shape list leans on wide shapes (few rows, large parts) on purpose:
//! those stress the frontier hardest, and the other two harnesses barely cover
//! them. Each case clears the caches first, so none is warmed by an earlier one.
use std::time::Instant;
use symfn::{clear_caches, LrBackend, Partition, SkewLr};

fn main() {
    let cases: Vec<Vec<u32>> = vec![
        // Wide: few rows, large parts. The hard regime.
        vec![12, 10, 8],
        vec![20, 16, 12],
        vec![16, 13, 10, 7],
        // Staircases: the regime the frontier merges best.
        vec![8, 7, 6, 5, 4, 3],
        vec![9, 8, 7, 6, 5],
        vec![8, 7, 6, 5, 4],
        vec![6, 5, 4, 3, 2, 1],
        // Rectangle and tall, as controls.
        vec![7, 7, 7, 7, 7],
        vec![3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ];
    let tag = std::env::args().nth(1).unwrap_or_else(|| "?".into());
    for sh in cases {
        let p = Partition::new(sh.iter().copied());
        clear_caches();
        let t = Instant::now();
        let v = SkewLr.schur_product(&p, &p);
        let dt = t.elapsed().as_secs_f64();
        println!("{tag}\t{p}\t{dt:.4}\t{}", v.len());
    }
}
