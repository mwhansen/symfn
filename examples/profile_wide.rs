//! Profiling harness for the one regime lrcalc still beats us on: three-row
//! wide shapes between roughly 0.02s and 0.4s.
//!
//! A single `[20,16,12]²` takes ~0.3s, which is far too short for a sampling
//! profiler to say anything, so this repeats the product until a deadline and
//! prints the per-iteration time. Attach with the profiler macOS already ships:
//!
//! ```text
//!   cargo build --profile profiling --example profile_wide
//!   ./target/profiling/examples/profile_wide 20,16,12 25 &
//!   sample $! 20 -mayDie -f /tmp/wide.txt
//! ```
//!
//! The caches are cleared every iteration. Without that the second and every
//! later product is a memo hit and the profile is all hashing — measuring the
//! cache rather than the algorithm it is supposed to be hiding.
//!
//! Compare the printed per-iteration time against `scripts/compare_lrcalc.py`
//! for the same shape: if they disagree, the profile is of something other than
//! what the comparison measures and the attribution below it means nothing.

use std::time::{Duration, Instant};

use symfn::{AutoLr, LrBackend, Partition};

fn main() {
    let mut args = std::env::args().skip(1);
    let shape = args.next().unwrap_or_else(|| "20,16,12".into());
    let secs: u64 = args
        .next()
        .and_then(|s| s.parse().ok())
        .unwrap_or(25);

    let mu = Partition::new(shape.split(',').map(|x| x.trim().parse::<u32>().expect("part")));
    eprintln!("profiling s{mu}² for {secs}s (pid {})", std::process::id());

    let deadline = Instant::now() + Duration::from_secs(secs);
    let (mut iters, mut terms) = (0u32, 0usize);
    let start = Instant::now();
    while Instant::now() < deadline {
        symfn::clear_caches();
        terms = AutoLr.schur_product(&mu, &mu).len();
        iters += 1;
    }
    let per = start.elapsed().as_secs_f64() / f64::from(iters);
    eprintln!("{iters} iterations, {terms} terms, {per:.4}s per product");
}
