//! Is `c^λ_{μν}(1)` the double coset connection coefficient, and in what
//! normalization?
//!
//! ```text
//!   cargo run --release --example gj_b1 -- 5
//! ```
//!
//! [GJ] specialize at `b = 1` to the double coset algebra of the hyperoctahedral
//! group, as `b = 0` gives the class algebra of `S_n`. The `b = 0` pin is
//! already in place. This finds the constant relating the two sides at `b = 1`
//! before it gets written into a test as though it were known.

use symfn::gj::{double_coset_coefficient, gj_connection_tables};
use symfn::Partition;

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(4);

    for n in 1..=top {
        let t = gj_connection_tables(n);
        let parts = symfn::partitions_of(n);
        println!("== n = {n} ==");
        let mut ratios: Vec<(String, f64)> = Vec::new();
        for la in &parts {
            for mu in &parts {
                for nu in &parts {
                    let key = (la.clone(), mu.clone(), nu.clone());
                    let c = t.c.get(&key).map_or((0i128, 1u128), |p| {
                        // c(1) = (Σ num) / den
                        (p.num.iter().sum::<i128>(), p.den)
                    });
                    let b = double_coset_coefficient(la, mu, nu);
                    if c.0 == 0 && b == 0 {
                        continue;
                    }
                    let cval = c.0 as f64 / c.1 as f64;
                    let ratio = if b == 0 { f64::NAN } else { cval / b as f64 };
                    ratios.push((format!("{la} {mu} {nu}"), ratio));
                    if n <= 4 {
                        println!(
                            "  {:<8} {:<8} {:<8}  c(1) = {:<10} b = {:<8} ratio = {ratio}",
                            format!("{la}"),
                            format!("{mu}"),
                            format!("{nu}"),
                            format!("{}/{}", c.0, c.1),
                            b
                        );
                    }
                }
            }
        }
        // Is the ratio constant? That is the only question here.
        let finite: Vec<f64> = ratios
            .iter()
            .map(|(_, r)| *r)
            .filter(|r| r.is_finite())
            .collect();
        let (lo, hi) = finite
            .iter()
            .fold((f64::INFINITY, f64::NEG_INFINITY), |(a, b), &r| {
                (a.min(r), b.max(r))
            });
        let nan = ratios.len() - finite.len();
        println!(
            "  {} live triples, ratio in [{lo}, {hi}], {nan} with b = 0 but c(1) != 0",
            ratios.len()
        );
        if nan > 0 {
            println!("  ^^ a nonzero c(1) where the count is zero kills a pure scalar relation");
            for (k, r) in ratios.iter().filter(|(_, r)| !r.is_finite()).take(5) {
                println!("     {k}  ratio {r}");
            }
        }
        println!();
    }
}
