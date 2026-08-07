//! A single, attributable workload for the sampling profiler.
//!
//! ```text
//!   cargo build --profile profiling --example profile_llt
//!   ./target/profiling/examples/profile_llt r1 9 &
//!   sample $! 20 -f /tmp/llt_r1.txt
//! ```
//!
//! `bench_llt` measures six different things and its profile mixes all of
//! them. This one does exactly one, chosen by the first argument, because the
//! routes share no mathematics and there is no reason to expect their hot spots
//! to coincide:
//!
//! - `r1` — the by-path shuffle refinement, i.e. `#SYT` walks over Dyck-path
//!   tuples. This is where `nabla_e_by_path(10)`'s 80 s lives.
//! - `r2` — whole-degree `H^(k)` tables, i.e. the abacus strip walk.
//! - `r3` — Kazhdan–Lusztig columns, i.e. Fock straightening.

use symfn::llt::{llt_h_table, llt_kl_column, nabla_e_by_path};
use symfn::qt::QtPoly;
use symfn::sym::{Monomial, SymFn};
use symfn::Partition;

fn main() {
    let route = std::env::args().nth(1).unwrap_or_else(|| "r1".to_string());
    let n: u32 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);
    let k: u32 = std::env::args()
        .nth(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);

    let mut keep = 0usize;
    loop {
        symfn::clear_caches();
        keep += match route.as_str() {
            "r1" => nabla_e_by_path::<i128>(n)
                .iter()
                .map(|(_, g)| g.terms().len())
                .sum::<usize>(),
            "r2" => {
                let t: Vec<(Partition, Monomial<QtPoly<i64>>)> = llt_h_table(n, k);
                t.iter().map(|(_, f)| f.terms().len()).sum::<usize>()
            }
            "r3" => symfn::partitions_of(n)
                .iter()
                .map(|l| llt_kl_column::<i64>(l, k).len())
                .sum::<usize>(),
            other => panic!("unknown route {other}; want r1, r2 or r3"),
        };
        // Keep the result alive so nothing is optimized away.
        if keep == usize::MAX {
            break;
        }
    }
}
