//! A single, attributable Schubert workload for the sampling profiler.
//!
//! ```text
//!   cargo build --profile profiling --example profile_schubert
//!   ./target/profiling/examples/profile_schubert e3 6 &
//!   sample $! 20 -f /tmp/schubert_e3.txt
//! ```
//!
//! `bench_schubert` runs the whole ladder and its profile mixes engines and
//! sizes. This one does exactly one product with one engine.
//!
//! - `e3 k` — E3 on `stair_k²` (the limiting case; k=6 is ~2.5s, k=7 ~140s)
//! - `e1 k` — E1 on the same, for contrast
//! - `monk k` — `mul_variable` alone, hammered on a large element, since that
//!   is the operation every E3 state performs

use std::time::Instant;
use symfn::permutation::Perm;
use symfn::schubert::{peel_states, Schubert};

fn stair(k: u32) -> Vec<u32> {
    (1..=k)
        .map(|i| 2 * i)
        .chain((1..=k).map(|i| 2 * i - 1))
        .collect()
}

fn main() {
    let route = std::env::args().nth(1).unwrap_or_else(|| "e3".to_string());
    let k: u32 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);
    let w = Perm::new(stair(k)).unwrap();
    let a: Schubert<i128> = Schubert::monomial(w, 1);
    eprintln!("stair{k}: states={} ", peel_states(&w));

    let t = Instant::now();
    let n = match route.as_str() {
        "e3" => a.mul(&a).terms().len(),
        "e1" => a.mul_naive(&a).terms().len(),
        // E2 is fast enough that one product is too short to sample; loop it
        // until there is enough signal, and report per-iteration time.
        "e2" => {
            let mut n = 0;
            let mut iters = 0;
            while t.elapsed().as_secs_f64() < 25.0 {
                n = a.mul_e2(&a).terms().len();
                iters += 1;
            }
            eprintln!(
                "e2 iters={iters} per-iter={:.4}s",
                t.elapsed().as_secs_f64() / iters as f64
            );
            n
        }
        "monk" => {
            // one big element, then many single-variable multiplications --
            // the inner loop of every E3 state
            let mut f = a.mul(&a);
            for _ in 0..200 {
                f = f.mul_variable(1 + (f.terms().len() as u32 % 3));
            }
            f.terms().len()
        }
        other => panic!("unknown route {other}"),
    };
    eprintln!(
        "{route} stair{k}: {n} terms in {:.4}s",
        t.elapsed().as_secs_f64()
    );
}
