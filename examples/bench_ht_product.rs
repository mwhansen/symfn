//! The reduced-Kronecker route through the `h̃` basis, timed and profiled.
//!
//! `reduced_kronecker_via_ht` runs `ht_product_terms` once per pair of `h̃`
//! rows, and that function's leaf enumerates matrices under a budget of
//! `HT_PRODUCT_BUDGET` — so its per-leaf cost is the whole cost of the route.
//! `bench_htilde` measures a different thing (the modified-Macdonald table) and
//! nothing here had a harness.
//!
//! ```text
//!   cargo build --release --example bench_ht_product
//!   cp target/release/examples/bench_ht_product /tmp/after    # and likewise before
//!   for i in 1 2 3; do /tmp/before before; /tmp/after after; done
//!
//!   cargo build --profile profiling --example bench_ht_product
//!   ./target/profiling/examples/bench_ht_product loop 7 &
//!   sample $! 8 -mayDie -f /tmp/ht.txt
//! ```
//!
//! Interleave and take the min per (build, case): as everywhere else in this
//! tree, two consecutive runs of different builds prove nothing. Columns: tag,
//! degree, seconds, terms. Caches are cleared per case.

use std::time::Instant;

use symfn::{clear_caches, partitions_of, reduced_kronecker_via_ht, Partition, St, SymFn};

fn sweep(n: u32) -> (f64, usize) {
    let parts: Vec<Partition> = partitions_of(n).to_vec();
    clear_caches();
    let t = Instant::now();
    let mut terms = 0usize;
    for a in &parts {
        for b in &parts {
            if let Some(r) = reduced_kronecker_via_ht::<i128>(a, b) {
                terms += <St<i128> as SymFn<i128>>::terms(&r).len();
            }
        }
    }
    (t.elapsed().as_secs_f64(), terms)
}

fn main() {
    let mut args = std::env::args().skip(1);
    let first = args.next().unwrap_or_else(|| "bench".into());
    if first == "loop" {
        // Repeat one degree for a sampling profiler; caches cleared each pass,
        // or every pass after the first measures the memo instead.
        let n: u32 = args.next().and_then(|s| s.parse().ok()).unwrap_or(7);
        let mut keep = 0usize;
        loop {
            keep += sweep(n).1;
            if keep == usize::MAX {
                break;
            }
        }
    }
    for n in [5u32, 6, 7] {
        let (dt, terms) = sweep(n);
        println!("{first}\tht_kronecker_n{n}\t{dt:.6}\t{terms} terms");
    }
}
