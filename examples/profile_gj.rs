//! Sampling targets for the two things that consume `J → p`.
//!
//! ```text
//!   cargo build --profile profiling --example profile_gj
//!   ./target/profiling/examples/profile_gj gj 9      & sample $! 20 -f /tmp/gj.txt
//!   ./target/profiling/examples/profile_gj stanley 6 & sample $! 20 -f /tmp/st.txt
//! ```
//!
//! `J → p` is our slowest unit against Sage (1102x where `P → m` is 8353x), so
//! it looks like the thing to optimise. Whether it *is* depends entirely on
//! what fraction of its two consumers it accounts for, and that is a
//! measurement, not an inference.

use symfn::afrac::AFrac;

fn main() {
    let what = std::env::args().nth(1).unwrap_or_else(|| "gj".into());
    let n: u32 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    loop {
        symfn::clear_caches();
        match what.as_str() {
            "stanley" => {
                let parts = symfn::partitions_of(n);
                let big = symfn::partitions_of(2 * n);
                let mut live = 0usize;
                for la in &parts {
                    for mu in &parts {
                        for nu in &big {
                            let g: AFrac<i128> = symfn::jack_structure_constant(la, mu, nu);
                            live += usize::from(!<AFrac<i128> as symfn::Ring>::is_zero(&g));
                        }
                    }
                }
                if live == usize::MAX {
                    break;
                }
            }
            "jp" => {
                let t = symfn::jack_powersum_table::<i128>(n);
                if t.is_empty() {
                    break;
                }
            }
            "modular" => {
                let t = symfn::gj_connection_tables_modular(n);
                if !t.laws_hold() {
                    break;
                }
            }
            _ => {
                let t = symfn::gj_connection_tables(n);
                if !t.laws_hold() {
                    break;
                }
            }
        }
    }
}
