//! A single, attributable plethysm workload for the sampling profiler.
//!
//! ```text
//!   cargo build --profile profiling --example profile_plethysm
//!   ./target/profiling/examples/profile_plethysm ladder 7 7 &
//!   sample $! 20 -f /tmp/plethysm.txt
//! ```
//!
//! `bench_ops` carries plethysm cases, but they are sub-millisecond and their
//! profile is marshalling and setup. These run one case, long enough to
//! sample, on one route:
//!
//! - `ladder n m` — `s_n[s_m]` through the Schur-basis recursion, the route a
//!   one-row inner takes. `7 7` is ~18s, `6 6` ~0.3s.
//! - `power n mu..` — `s_n[s_mu]` through the power-sum basis, by giving an
//!   inner with more than one row so the ladder declines it. `5 4 2` is ~1s
//!   in release and is the route whose cost is the p → s conversion.
//!
//! The two are worth profiling separately because they share nothing: one is
//! Littlewood–Richardson products over an abacus rule, the other is z_μ,
//! characters and a β-mask sweep.

use std::time::Instant;
use symfn::{Partition, Rational, Ring, Schur, SymFn};

fn s(v: &[u32]) -> Schur<Rational> {
    Schur::monomial(Partition::new(v.iter().copied()), Rational::one())
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let route = args.first().map(String::as_str).unwrap_or("ladder");
    let nums: Vec<u32> = args[1..].iter().filter_map(|a| a.parse().ok()).collect();

    let (outer, inner) = match route {
        "ladder" => {
            let n = *nums.first().unwrap_or(&7);
            let m = *nums.get(1).unwrap_or(&7);
            (s(&[n]), s(&[m]))
        }
        "power" => {
            let n = *nums.first().unwrap_or(&5);
            let mu: Vec<u32> = nums[1..].to_vec();
            let mu = if mu.is_empty() { vec![4, 2] } else { mu };
            (s(&[n]), s(&mu))
        }
        other => {
            eprintln!("unknown route {other:?}; expected `ladder` or `power`");
            std::process::exit(2);
        }
    };

    let t = Instant::now();
    let r = symfn::plethysm(&outer, &inner);
    println!(
        "{route}: {:?}, {} terms, degree {:?}",
        t.elapsed(),
        r.terms().len(),
        r.degree()
    );
}
