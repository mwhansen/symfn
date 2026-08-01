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
//! Columns: tag, shape, seconds, terms, peak live layer states. The peak
//! state count is the memory axis — bytes-per-state times that number bounds
//! the layer's residency, and unlike RSS it is allocator-independent.
//!
//! Extra arguments after the tag select shapes explicitly (`24,20,16,12`),
//! replacing the built-in list — that is how to put one case alone in a
//! process so `/usr/bin/time -l` attributes peak RSS to it. `ORIENT=conj`
//! runs the same product on the conjugate diagram (the coefficients are equal
//! by c^λ_{μν} = c^{λ'}_{μ'ν'}), to compare the two orientations' layers.
//!
//! The shape list leans on wide shapes (few rows, large parts) on purpose:
//! those stress the layer hardest, and the other two harnesses barely cover
//! them. Each case clears the caches first, so none is warmed by an earlier one.
use std::time::Instant;
use symfn::skew_lr::take_peak_layer_states;
use symfn::{clear_caches, LrBackend, Partition, SkewLr};

/// s_p · s_p by expanding the *conjugate* juxtaposed shape.
///
/// Builds the disconnected diagram for p' ⊔ p' and expands that; by
/// c^λ_{μν} = c^{λ'}_{μ'ν'} the multiset of coefficients is the same, so
/// timing and layer size are comparable case-for-case with the direct
/// orientation (terms are reported conjugated back, as a correctness check).
fn product_conjugate(p: &Partition) -> Vec<(Partition, u128)> {
    let q = p.conjugate();
    let shift = q.part(0);
    let mut outer: Vec<u32> = q.parts().iter().map(|&x| x + shift).collect();
    outer.extend_from_slice(q.parts());
    let inner = vec![shift; if shift == 0 { 0 } else { q.len() }];
    let expansion = symfn::expand_skew(&Partition::new(outer), &Partition::new(inner));
    expansion
        .into_iter()
        .map(|(l, c)| (l.conjugate(), c))
        .collect()
}

fn main() {
    let default_cases: Vec<Vec<u32>> = vec![
        // Wide: few rows, large parts. The hard regime.
        vec![12, 10, 8],
        vec![20, 16, 12],
        vec![16, 13, 10, 7],
        // Staircases: the regime the layer merges best.
        vec![8, 7, 6, 5, 4, 3],
        vec![9, 8, 7, 6, 5],
        vec![8, 7, 6, 5, 4],
        vec![6, 5, 4, 3, 2, 1],
        // Rectangle and tall, as controls.
        vec![7, 7, 7, 7, 7],
        vec![3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3],
    ];
    let mut args = std::env::args().skip(1);
    let tag = args.next().unwrap_or_else(|| "?".into());
    let explicit: Vec<Vec<u32>> = args
        .map(|a| {
            a.split(',')
                .map(|x| x.parse().expect("shape part"))
                .collect()
        })
        .collect();
    let cases = if explicit.is_empty() {
        default_cases
    } else {
        explicit
    };
    let conj = std::env::var("ORIENT").as_deref() == Ok("conj");

    for sh in cases {
        let p = Partition::new(sh.iter().copied());
        clear_caches();
        take_peak_layer_states();
        let t = Instant::now();
        let v = if conj {
            product_conjugate(&p)
        } else {
            SkewLr.schur_product(&p, &p)
        };
        let dt = t.elapsed().as_secs_f64();
        let peak = take_peak_layer_states();
        println!("{tag}\t{p}\t{dt:.4}\t{}\t{peak}", v.len());
    }
}
