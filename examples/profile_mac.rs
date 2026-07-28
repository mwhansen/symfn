//! A long single-workload run for a sampling profiler.
//!
//! ```text
//!   cargo run --release --example profile_mac -- 10 P &
//!   sample $! 10 -f /tmp/mac.txt
//! ```

use symfn::sym::SymFn;
use symfn::{macdonald_j, macdonald_p, macdonald_q, Frac, Monomial};

fn main() {
    let n: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(10);
    let basis = std::env::args().nth(2).unwrap_or_else(|| "P".into());
    let reps: u32 = std::env::args()
        .nth(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(1);
    let mut sink = 0usize;
    for _ in 0..reps {
        for lambda in symfn::partitions_of(n) {
            let f: Monomial<Frac<i128>> = match basis.as_str() {
                "P" => macdonald_p(&lambda),
                "Q" => macdonald_q(&lambda),
                "J" => macdonald_j(&lambda),
                other => panic!("unknown basis {other}"),
            };
            sink += f.terms().len();
        }
    }
    println!("{sink}");
}
