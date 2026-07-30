//! Heap accounting per subsystem — the exploratory half of the memory harness.
//!
//! ```text
//!   cargo run --release --example heapstat            # list the workloads
//!   cargo run --release --example heapstat -- htilde  # measure one
//!   HEAPSTAT_HIST=1 cargo run --release --example heapstat -- htilde
//! ```
//!
//! Reports peak live bytes, total bytes, and allocation count — the three
//! numbers `docs/roadmap/memory.md` is about, and unlike RSS the first two
//! repeat exactly on a single-threaded workload. `HEAPSTAT_HIST=1` adds the
//! size-class histogram, which is what attributes churn to a specific buffer.
//!
//! **One workload per process**: the memo tables are global, so a second
//! workload in the same run would measure a warm cache.
//!
//! Workloads live in `symfn::measure::workloads`, shared with `tests/memory.rs`
//! so a new entry gets a report *and* a regression budget. This binary counts
//! allocations, so it perturbs timing — do not benchmark wall clock with it.

#[global_allocator]
static ALLOC: symfn::measure::Counting = symfn::measure::Counting::new();

use symfn::measure::{self, workloads};

fn main() {
    let Some(which) = std::env::args().nth(1) else {
        println!("workloads: {}", workloads::names().join(" "));
        return;
    };
    let Some(w) = workloads::find(&which) else {
        println!(
            "unknown workload `{which}`; try: {}",
            workloads::names().join(" ")
        );
        std::process::exit(1);
    };

    let (note, stats) = measure::measure(w.run);
    println!("{:<16} {stats}   {note}", w.name);

    if let Err(over) = w.budget.check(&stats) {
        println!("\n{over}");
    }

    if std::env::var("HEAPSTAT_HIST").is_ok() {
        println!("\n  size class        allocs         bytes");
        for (size, n, bytes) in measure::histogram() {
            println!(
                "  {size:>9}   {n:>11}   {:>9.1} MB",
                bytes as f64 / 1048576.0
            );
        }
    }
}
