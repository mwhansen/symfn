//! Hall–Littlewood: shared-suffix table vs one call per shape.
//!
//! Symmetrica's `hall_littlewood` recurses from scratch on every call. The
//! recursion peels the largest part, so its subproblems are the *suffixes* of λ,
//! and the partitions of n share those heavily — computing a whole degree should
//! therefore cost far less than p(n) independent calls. This measures that.
//!
//! Passes are **rotated**, not merely interleaved: an earlier benchmark in this
//! crate reported a 36% difference that turned out to be entirely position — the
//! pass that ran after a much larger one paid for a cold cache.

use std::time::Instant;

use symfn::sym::SymFn;
use symfn::{hall_littlewood, hall_littlewood_table, partitions_of, QtPoly, Schur};

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(16);

    println!(
        "{:>4} {:>6} {:>11} {:>11} {:>8}",
        "n", "p(n)", "table", "one-by-one", "speedup"
    );
    for n in 4..=top {
        let parts = partitions_of(n);
        let mut table_s = f64::MAX;
        let mut solo_s = f64::MAX;
        // Two rounds with the order swapped; keep the best of each so whichever
        // ran second cannot be charged for the other's cache footprint.
        for round in 0..2 {
            if round == 0 {
                table_s = table_s.min(time(|| {
                    let t: Vec<(_, Schur<QtPoly<i64>>)> = hall_littlewood_table(n);
                    t.len()
                }));
                solo_s = solo_s.min(time(|| solo(&parts)));
            } else {
                solo_s = solo_s.min(time(|| solo(&parts)));
                table_s = table_s.min(time(|| {
                    let t: Vec<(_, Schur<QtPoly<i64>>)> = hall_littlewood_table(n);
                    t.len()
                }));
            }
        }
        println!(
            "{n:>4} {:>6} {table_s:>11.4} {solo_s:>11.4} {:>7.2}x",
            parts.len(),
            solo_s / table_s
        );
    }
}

fn solo(parts: &[symfn::Partition]) -> usize {
    parts
        .iter()
        .map(|l| {
            let f: Schur<QtPoly<i64>> = hall_littlewood(l);
            f.terms().len()
        })
        .sum()
}

fn time<T>(f: impl FnOnce() -> T) -> f64 {
    let start = Instant::now();
    let out = f();
    let secs = start.elapsed().as_secs_f64();
    std::hint::black_box(out);
    secs
}
