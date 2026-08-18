//! Peak *layer* size for one skew expansion of `mu^2`, with no allocator
//! instrumentation in the way.
//!
//! `lrheap` reports the same expansion's heap; this one leaves the global
//! allocator alone, so it is the build to time and to sample-profile. The shape
//! is given on the command line — `lrmem 20 16 12 8` — because the shapes that
//! expose a scaling effect are larger than the ones a preset list would carry.

use symfn::{
    skew_lr::{expand_skew_shared, take_peak_layer_states},
    Partition,
};

fn main() {
    let parts: Vec<u32> = std::env::args()
        .skip(1)
        .map(|a| a.parse::<u32>().expect("integer part"))
        .collect();
    let mu = Partition::new(if parts.is_empty() {
        vec![10, 8, 6, 4]
    } else {
        parts
    });
    // mu^2 is the square skew shape whose expansion is the product s_mu * s_mu:
    // mu shifted right by its own width, sitting above an unshifted copy.
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
        .chain(std::iter::repeat_n(0, mu.len()))
        .collect();
    let r = expand_skew_shared(&Partition::new(outer), &Partition::new(inner));
    eprintln!(
        "{mu}^2\t{} terms\tpeak_states {}",
        r.len(),
        take_peak_layer_states()
    );
}
