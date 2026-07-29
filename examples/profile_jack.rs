//! A single, attributable workload for the sampling profiler.
//!
//! ```text
//!   cargo build --profile profiling --example profile_jack
//!   ./target/profiling/examples/profile_jack 20 &
//!   sample $! 20 -f /tmp/jack.txt
//! ```
//!
//! `bench_jack` measures several different things and its profile mixes them.
//! This one does exactly one: the whole-degree `P → m` table, over `i128`, on
//! repeat, so every sample lands in the engine.
//!
//! Second argument picks the route, so the two engines can be profiled
//! separately — they share no mathematics and there is no reason to expect
//! their hot spots to coincide.

use symfn::afrac::AFrac;
use symfn::sym::Monomial;

fn main() {
    let n: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(18);
    let branching = std::env::args().nth(2).is_some_and(|s| s == "branch");

    let mut terms = 0usize;
    loop {
        symfn::clear_caches();
        for lambda in symfn::partitions_of(n) {
            let p: Monomial<AFrac<i128>> = if branching {
                symfn::jack_p_branching(&lambda)
            } else {
                symfn::jack_p_lb(&lambda)
            };
            // Keep the result alive so nothing is optimised away.
            terms += <Monomial<AFrac<i128>> as symfn::sym::SymFn<_>>::terms(&p).len();
        }
        if terms == usize::MAX {
            break; // unreachable; defeats dead-code elimination
        }
    }
}
