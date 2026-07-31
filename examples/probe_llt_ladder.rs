//! The Delta-conjecture **rise** ladder: the per-path LLT route against the
//! labelled-path walk it replaced.
//!
//! ```text
//!   cargo run --release --example probe_llt_ladder -- 8
//! ```
//!
//! `dyck::ladder(n, Side::Rise)` now goes through the per-path vertical-strip
//! LLT polynomials — `Rise` and its weights `t^{−a_i}` are functions of the
//! area sequence alone, so the `z`-extraction factors out of the labelling sum
//! (see `dyck.rs`'s module docs and `docs/record/llt.md`). The
//! labelled walk survives as `dyck::ladder_at_content`, which is both the
//! correctness oracle (`the_rise_ladder_via_llt_agrees_with_the_labelled_walk`)
//! and the cheaper route for a single coarse content.
//!
//! This measures what the swap bought, and re-checks that the two agree at every
//! `k` while it does — a redistribution between `k`s is exactly what an aggregate
//! check would miss.
//!
//! ⚠️ Same-run A/B ratios, so they hold on battery; the absolute seconds do not
//! compare to `docs/record/llt.md`'s mains table.

use std::time::Instant;

use symfn::qt::QtPoly;
use symfn::sym::{Monomial, SymFn};
use symfn::{ladder, ladder_at_content, Side};

/// The pre-LLT route: one labelled-path walk per content.
fn via_labellings(n: u32) -> Vec<Monomial<QtPoly<i128>>> {
    let mut out = vec![Monomial::zero(); n as usize];
    for mu in symfn::partitions_of(n) {
        for (k, c) in ladder_at_content::<i128>(&mu, Side::Rise)
            .into_iter()
            .enumerate()
        {
            if !c.is_empty() {
                out[k].add_term(mu.clone(), c);
            }
        }
    }
    out
}

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(7);

    println!("== the rise ladder: per-path LLT vs the labelled walk ==");
    println!(
        "{:>3} {:>12} {:>14} {:>9}  agree at every k?",
        "n", "llt(s)", "labellings(s)", "speedup"
    );
    for n in 1..=top {
        let t0 = Instant::now();
        let fast = ladder::<i128>(n, Side::Rise);
        let a = t0.elapsed().as_secs_f64();

        let t0 = Instant::now();
        let slow = via_labellings(n);
        let b = t0.elapsed().as_secs_f64();

        let bad: Vec<usize> = (0..n as usize).filter(|&k| fast[k] != slow[k]).collect();
        let verdict = if bad.is_empty() {
            "yes".to_string()
        } else {
            format!("*** DISAGREE at k = {bad:?} ***")
        };
        println!(
            "{n:>3} {a:>12.4} {b:>14.4} {:>8.1}x  {verdict}",
            b / a.max(1e-9)
        );
    }
}
