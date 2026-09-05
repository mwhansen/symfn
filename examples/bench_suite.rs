//! symfn against its own history: one timing per workload in the catalog.
//!
//! ```text
//!   cargo run --release --example bench_suite                 # every workload, min of 3
//!   cargo run --release --example bench_suite -- 5 hl jack    # 5 repetitions, two workloads
//!   cargo run --release --example bench_suite > after.tsv
//!   python3 scripts/bench_compare.py docs/record/bench_suite.tsv after.tsv
//! ```
//!
//! The `bench_*` examples each measure one subsystem against an external
//! baseline and go as deep as that comparison needs. This one is shallow and
//! wide: the twelve-plus workloads of [`symfn::measure::workloads`], which
//! already span the subsystems and already have a size chosen where the work
//! lives, timed with the caches cleared before every repetition and the
//! minimum kept. The output is tab-separated so that two runs can be diffed
//! by `scripts/bench_compare.py`, which is what turns a change into a ratio
//! per subsystem. The committed run is `docs/record/bench_suite.tsv`, with the
//! machine and power state in its header; a comparison is only meaningful on
//! the machine that produced the file it is compared against.
//!
//! Not a test. Wall time is not assertable (`docs/record/memory.md`, "Why
//! memory can be a test when time cannot"), which is why the workload catalog
//! carries memory budgets and this carries none. Run it without the counting
//! allocator — `heapstat` is the harness that installs it — because the four
//! atomic adds per allocation are visible in these numbers.

use std::time::Instant;

use symfn::measure::workloads::{self, WORKLOADS};

fn main() {
    let mut args = std::env::args().skip(1);
    let reps: u32 = args.next().and_then(|s| s.parse().ok()).unwrap_or(3);
    let names: Vec<String> = args.collect();
    let chosen: Vec<&workloads::Workload> = if names.is_empty() {
        WORKLOADS.iter().collect()
    } else {
        names
            .iter()
            .map(|n| workloads::find(n).unwrap_or_else(|| panic!("no workload named {n}")))
            .collect()
    };

    println!(
        "# bench_suite\tmin of {reps}\tsymfn {}",
        env!("CARGO_PKG_VERSION")
    );
    println!("workload\tseconds\tnote");
    for w in chosen {
        let mut best = f64::INFINITY;
        let mut note = String::new();
        for _ in 0..reps {
            symfn::clear_caches();
            let t0 = Instant::now();
            note = (w.run)();
            best = best.min(t0.elapsed().as_secs_f64());
        }
        println!("{}\t{best:.4}\t{note}", w.name);
    }
}
