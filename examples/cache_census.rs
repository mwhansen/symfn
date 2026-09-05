//! The caches over a session-shaped sweep: what each table holds at each
//! degree, and what a byte budget costs in wall time.
//!
//! ```text
//!   cargo run --release --example cache_census
//!   cargo run --release --example cache_census -- 16 12 8 5   # the four ceilings
//! ```
//!
//! The sweep is the one `measure::workloads` runs as `session`, sized by four
//! degree ceilings: characters, Kostka numbers, products and skews, and the
//! `s → J` transition. Unbounded first, printing every table's bytes after
//! each degree — the census `docs/record/memory.md` Rule 4 asks for — then
//! the same sweep cold under a ladder of budgets, printing wall time and what
//! was held at the end. The wheel's default budget is argued from this
//! table in the record; here it is only measured.
//!
//! AC power, and say so when recording (`CLAUDE.md`, "Measurement
//! discipline").

use std::time::Instant;

use symfn::measure::workloads::session_sweep;
use symfn::{cache_stats, clear_caches, set_cache_budget};

fn held() -> usize {
    cache_stats().iter().map(|r| r.bytes).sum()
}

fn mb(bytes: usize) -> f64 {
    bytes as f64 / 1_048_576.0
}

fn main() {
    let args: Vec<u32> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("a degree ceiling"))
        .collect();
    let [chars, kostka, products, transition] = match args.as_slice() {
        [] => [16, 12, 8, 5],
        [a, b, c, d] => [*a, *b, *c, *d],
        _ => panic!("usage: cache_census [chars kostka products transition]"),
    };
    let ceilings = (chars, kostka, products, transition);

    println!("census, unbounded: bytes per table after each degree");
    clear_caches();
    set_cache_budget(None);
    let mut last = cache_stats();
    let t0 = Instant::now();
    session_sweep(ceilings, |degree| {
        let now = cache_stats();
        let grew: Vec<String> = now
            .iter()
            .zip(last.iter())
            .filter(|(a, b)| a.bytes != b.bytes)
            .map(|(a, _)| format!("{} {:.2} MB", a.name, mb(a.bytes)))
            .collect();
        println!(
            "  degree {degree:>2}  held {:>8.2} MB  {:>7.2}s   {}",
            mb(held()),
            t0.elapsed().as_secs_f64(),
            grew.join(", ")
        );
        last = now;
    });
    println!(
        "  unbounded, with the census reads: {:.2}s, held {:.2} MB at the end",
        t0.elapsed().as_secs_f64(),
        mb(held())
    );

    // The ladder is cold at every rung and the unbounded rung is its own
    // baseline, so the census reads above are not in the ratio.
    println!("\nthe same sweep, cold, under a budget");
    println!(
        "  {:>12} {:>9} {:>9} {:>12}",
        "budget", "time", "vs none", "held at end"
    );
    let ladder: Vec<Option<usize>> = std::iter::once(None)
        .chain([1 << 30, 256 << 20, 64 << 20, 16 << 20, 4 << 20, 1 << 20].map(Some))
        .collect();
    let mut baseline = 0.0f64;
    for budget in ladder {
        clear_caches();
        set_cache_budget(budget);
        let t = Instant::now();
        session_sweep(ceilings, |_| {});
        let dt = t.elapsed().as_secs_f64();
        if budget.is_none() {
            baseline = dt;
        }
        println!(
            "  {:>12} {:>8.2}s {:>8.2}x {:>9.2} MB",
            budget.map_or("none".to_string(), |b| format!("{} MB", b >> 20)),
            dt,
            dt / baseline,
            mb(held())
        );
    }
    set_cache_budget(None);
    clear_caches();
}
