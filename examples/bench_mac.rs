//! Time the Macdonald and Hall–Littlewood-P layers.
//!
//! ```text
//!   cargo run --release --example bench_mac -- 8
//! ```
//!
//! Prints one line per (workload, degree). The order the workloads run in is
//! **rotated** per degree rather than fixed, so a warm allocator or a warm memo
//! cannot systematically favor whichever one happens to go first.

use std::time::Instant;

use symfn::sym::SymFn;
use symfn::{macdonald_j, macdonald_p, macdonald_q, Frac, Monomial};

fn mac_all(n: u32, basis: u8) -> usize {
    let mut terms = 0;
    for lambda in symfn::partitions_of(n) {
        let f: Monomial<Frac<i128>> = match basis {
            0 => macdonald_p(&lambda),
            1 => macdonald_q(&lambda),
            _ => macdonald_j(&lambda),
        };
        terms += f.terms().len();
    }
    terms
}

fn hl_p(n: u32) -> usize {
    symfn::hall_littlewood_p_table::<i64>(n)
        .iter()
        .map(|(_, f)| f.terms().len())
        .sum()
}

fn hl_qp(n: u32) -> usize {
    symfn::hall_littlewood_table::<i64>(n)
        .iter()
        .map(|(_, f)| f.terms().len())
        .sum()
}

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    let jobs: Vec<(&str, fn(u32) -> usize)> = vec![
        ("mac P  ", |n| mac_all(n, 0)),
        ("mac Q  ", |n| mac_all(n, 1)),
        ("mac J  ", |n| mac_all(n, 2)),
        ("hl P   ", hl_p),
        ("hl Q'  ", hl_qp),
    ];

    for n in 1..=top {
        // Rotate the starting job with the degree.
        let off = (n as usize) % jobs.len();
        let mut line: Vec<(usize, String)> = Vec::new();
        for k in 0..jobs.len() {
            let i = (k + off) % jobs.len();
            let (name, f) = jobs[i];
            let t0 = Instant::now();
            let terms = f(n);
            let dt = t0.elapsed().as_secs_f64();
            line.push((i, format!("  {name} n={n:2}  {dt:8.4}s  {terms:6} terms")));
        }
        line.sort_by_key(|(i, _)| *i);
        for (_, s) in line {
            println!("{s}");
        }
        println!();
    }
}
