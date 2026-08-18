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
//! - `general n mu..` — `s_n[s_mu]` with a multi-row inner, which is what
//!   drives the k-quotient sweep in `adams_schur`. `6 3 2` is ~2s.
//! - `power n mu..` — kept as a name for the same shape; the power-sum route
//!   is now the tests' oracle rather than a shipped path.
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
        "general" => {
            // A multi-row inner, which is what exercises the k-quotient sweep
            // in `adams_schur`; a one-row inner takes the closed form and
            // never enters it.
            let n = *nums.first().unwrap_or(&6);
            let mu: Vec<u32> = nums[1..].to_vec();
            let mu = if mu.is_empty() { vec![3, 2] } else { mu };
            (s(&[n]), s(&mu))
        }
        "power" => {
            let n = *nums.first().unwrap_or(&5);
            let mu: Vec<u32> = nums[1..].to_vec();
            let mu = if mu.is_empty() { vec![4, 2] } else { mu };
            (s(&[n]), s(&mu))
        }
        other => {
            eprintln!("unknown route {other:?}; expected `ladder`, `general` or `power`");
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
