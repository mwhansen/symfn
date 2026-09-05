//! Whole-element `s ↔ s̃` conversion, cold and warm, at full support.
//!
//! Full support per degree is the shape the Sage character-basis peel hands
//! this route once the peel is intercepted into whole-element conversions
//! (`docs/record/kronecker.md`, "`s → s̃` is now the whole cost of the
//! character bases in Sage"). The cold `s → s̃` column is the one the
//! per-degree row engine was measured with; `st → s` rides along because the
//! two directions share nothing past the caches and drift separately.
//!
//! ```text
//!   cargo run --release --example bench_s2st [degrees...]
//! ```
//!
//! With no arguments, degrees 10–16. Caches are cleared before each cold
//! measurement; the warm repeat is the memoized read path.

use std::time::Instant;

use symfn::{clear_caches, partitions_of, FromSchur, Schur, St, SymFn, ToSchur};

fn main() {
    let degrees: Vec<u32> = {
        let args: Vec<u32> = std::env::args()
            .skip(1)
            .map(|a| a.parse().expect("degrees are u32"))
            .collect();
        if args.is_empty() {
            vec![10, 12, 14, 16]
        } else {
            args
        }
    };
    println!("full-support conversions, cold caches per measurement");
    println!(
        "  {:>3} {:>8} {:>12} {:>12} {:>12}",
        "n", "shapes", "st->s cold", "s->st cold", "s->st warm"
    );
    for n in degrees {
        let mut s: Schur<i128> = Schur::zero();
        let mut st_full: St<i128> = St::zero();
        for lam in partitions_of(n) {
            s.add_term(lam.clone(), 1);
            st_full.add_term(lam, 1);
        }
        clear_caches();
        let t = Instant::now();
        let _back: Schur<i128> = st_full.to_schur();
        let st2s = t.elapsed().as_secs_f64();
        clear_caches();
        let t = Instant::now();
        let cold_out: St<i128> = St::from_schur(&s);
        let cold = t.elapsed().as_secs_f64();
        let t = Instant::now();
        let warm_out: St<i128> = St::from_schur(&s);
        let warm = t.elapsed().as_secs_f64();
        assert_eq!(cold_out.terms(), warm_out.terms());
        println!(
            "  {:>3} {:>8} {:>11.4}s {:>11.4}s {:>11.6}s",
            n,
            partitions_of(n).len(),
            st2s,
            cold,
            warm
        );
    }
}
