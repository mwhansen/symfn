//! Memory budgets — the regression half of the harness in `symfn::measure`.
//!
//! ```text
//!   cargo test --release --test memory -- --nocapture
//! ```
//!
//! Every workload in `symfn::measure::workloads` is run once and checked against
//! its budget. This is possible as a *test* — where a wall-clock benchmark is
//! not — because peak live bytes and allocation count are reproducible: a
//! single-threaded workload gives identical numbers run to run, and the threaded
//! ones move by about 0.1% (allocations) and 1% (peak) with the shard split.
//! RSS, by contrast, drifts with allocation history and is not assertable.
//!
//! **Release only.** A debug build allocates differently, and the numbers here
//! would be meaningless; the test asserts nothing under `debug_assertions`.
//!
//! Budgets are ceilings with headroom, not golden values, so ordinary churn does
//! not touch them. A failure prints the measured numbers and the `Budget` line
//! to paste back — a legitimate increase should be a visible, deliberate edit
//! rather than a silent drift.
//!
//! All workloads run inside **one** `#[test]`, sequentially. The counters are
//! process-wide (they must be: the LR expansion fans out across threads, and a
//! thread-local counter would miss the workers), so two measurements must never
//! overlap — and the test harness runs `#[test]` functions in parallel.

#[global_allocator]
static ALLOC: symfn::measure::Counting = symfn::measure::Counting::new();

use symfn::measure::{measure, workloads};

#[test]
fn workloads_stay_within_their_memory_budgets() {
    if cfg!(debug_assertions) {
        eprintln!("skipped: memory budgets are calibrated for --release");
        return;
    }

    // Freeing memory that predates the measurement drives the live counter
    // below zero — legal, and the reason it is signed. Unsigned it wrapped, and
    // the high-water mark latched ~2^64 bytes; run before the budgets so a
    // regression here is not read as a workload growing.
    let held: Vec<u8> = vec![7; 4 << 20];
    let ((), stats) = measure(move || drop(held));
    assert!(
        stats.peak < 1 << 20,
        "freeing pre-existing memory reported a peak of {} bytes",
        stats.peak
    );

    let mut failures = Vec::new();
    for w in workloads::WORKLOADS {
        let (note, stats) = measure(w.run);
        println!("{:<16} {stats}   {note}", w.name);
        if let Err(over) = w.budget.check(&stats) {
            failures.push(over);
        }
    }

    assert!(
        failures.is_empty(),
        "{} of {} workloads over budget\n\n{}",
        failures.len(),
        workloads::WORKLOADS.len(),
        failures.join("\n\n")
    );
}
