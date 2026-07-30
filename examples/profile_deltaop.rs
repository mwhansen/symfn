//! A single, attributable workload for the sampling profiler.
//!
//! ```text
//!   cargo build --profile profiling --example profile_deltaop
//!   ./target/profiling/examples/profile_deltaop nabla 11 &
//!   sample $! 15 -f /tmp/deltaop.txt
//! ```
//!
//! `bench_deltaop` measures several things and its profile mixes them. This one
//! does exactly one, chosen by the first argument:
//!
//! - `nabla` — `deltaop::nabla_e(n)`, the closed form. This is the crate's
//!   answer for ∇ totals and what beats Sage's `e[n].nabla()`.
//! - `delta` — `delta_prime_e(k, n)` over the whole k ladder.
//! - `htilde` — `bh::htilde_table(n)`, the Macdonald table both of the above
//!   are built on, so its profile is the shared substrate.

use symfn::qt::QtPoly;
use symfn::sym::{Schur, SymFn};

fn main() {
    let route = std::env::args()
        .nth(1)
        .unwrap_or_else(|| "nabla".to_string());
    let n: u32 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(11);

    let mut keep = 0usize;
    loop {
        symfn::clear_caches();
        keep += match route.as_str() {
            "nabla" => {
                let f: Schur<QtPoly<i128>> = symfn::nabla_e(n);
                f.terms().len()
            }
            "delta" => (0..n)
                .map(|k| {
                    let f: Schur<QtPoly<i128>> = symfn::delta_prime_e(k, n);
                    f.terms().len()
                })
                .sum::<usize>(),
            "htilde" => symfn::bh::htilde_table::<i128>(n).len(),
            other => panic!("unknown route {other}; want nabla, delta or htilde"),
        };
        if keep == usize::MAX {
            break;
        }
    }
}
