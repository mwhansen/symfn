//! Peak *heap* accounting for one skew expansion, via the shared counting
//! allocator in `symfn::measure`.
//!
//! RSS conflates live bytes with allocator retention; this separates them, and
//! divides by the peak frontier size to give bytes-per-state — the number that
//! says whether the *representation* is the problem. `heapstat` is the general
//! version; this stays because it is parameterised by shape and reports the
//! frontier counter alongside.
#[global_allocator]
static ALLOC: symfn::measure::Counting = symfn::measure::Counting::new();

use symfn::{
    clear_caches, measure,
    skew_lr::{expand_skew, take_peak_frontier_states},
    Partition,
};

fn p(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

fn main() {
    let which = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "0".into())
        .parse::<usize>()
        .unwrap();
    let shapes: [&[u32]; 4] = [
        &[10, 8, 6, 4],
        &[12, 10, 8, 6],
        &[8, 7, 6, 5, 4, 3],
        &[16, 13, 10, 7],
    ];
    let mu = p(shapes[which]);
    let w = mu.part(0);
    let outer: Vec<u32> = mu
        .parts()
        .iter()
        .map(|x| x + w)
        .chain(mu.parts().iter().copied())
        .collect();
    let inner: Vec<u32> = mu
        .parts()
        .iter()
        .map(|_| w)
        .chain(std::iter::repeat(0).take(mu.len()))
        .collect();
    clear_caches();
    let _ = take_peak_frontier_states();
    measure::reset();
    let r = expand_skew(&p(&outer), &Partition::new(inner.into_iter()));
    let states = take_peak_frontier_states();
    let peak = measure::snapshot().peak;
    eprintln!(
        "{mu}^2  {} terms  peak states {states}  peak heap {:.1} MB  = {:.0} bytes/state",
        r.len(),
        peak as f64 / 1048576.0,
        peak as f64 / states as f64
    );
}
