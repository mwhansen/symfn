//! The Goulden–Jackson tables `c^λ_{μν}(b)` and `h^λ_{μν}(b)`, degree by
//! degree, and where the engine stops.
//!
//! ```text
//!   cargo run --release --example gj_tables -- 9        # the ladder
//!   cargo run --release --example gj_tables -- 6 show   # and print degree 6
//! ```
//!
//! Two open conjectures live on these coefficients — Matchings-Jack on `c`,
//! the b-conjecture on `h` — and no package computes either table.
//! ℚ[b]-polynomiality is a theorem ([DF]) and `c`'s integrality is a theorem
//! ([BD]), so both are asserted; **positivity is the open question and is only
//! observed**, with any negative coefficient printed as a finding.
//!
//! Every degree is also pinned at `b = 0` against the class algebra of `S_n`,
//! computed from characters alone — the check that says [GJ]'s (1) and (5)
//! were transcribed correctly.

use std::time::Instant;

use symfn::gj::{gj_connection_tables, BPoly};
use symfn::Partition;

fn shape(p: &Partition) -> String {
    format!("{p}")
}

fn poly(b: &BPoly) -> String {
    let body = b
        .num
        .iter()
        .enumerate()
        .filter(|(_, &c)| c != 0)
        .map(|(k, c)| match k {
            0 => format!("{c}"),
            1 => format!("{c}b"),
            _ => format!("{c}b^{k}"),
        })
        .collect::<Vec<_>>()
        .join(" + ");
    let body = if body.is_empty() {
        "0".to_string()
    } else {
        body
    };
    if b.den == 1 {
        body
    } else {
        format!("({body})/{}", b.den)
    }
}

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);
    let show = std::env::args().nth(2).is_some_and(|s| s == "show");

    println!("== Goulden-Jackson connection tables ==");
    println!(
        "{:>3} {:>6} {:>10} {:>9} {:>9} {:>6} {:>5}  laws",
        "n", "p(n)", "build(s)", "c terms", "h terms", "b deg", "bits"
    );

    let mut findings = Vec::new();
    for n in 1..=top {
        symfn::clear_caches();
        let t0 = Instant::now();
        let t = gj_connection_tables(n);
        let secs = t0.elapsed().as_secs_f64();

        // b = 0 must be the class algebra — the transcription check.
        let parts = symfn::partitions_of(n);
        let mut mismatch = 0;
        for la in &parts {
            for mu in &parts {
                for nu in &parts {
                    let key = (la.clone(), mu.clone(), nu.clone());
                    let got = t.c.get(&key).map_or((0, 1), BPoly::at_zero);
                    let want = symfn::class_algebra_coefficient(la, mu, nu);
                    if got != (want, 1) {
                        mismatch += 1;
                    }
                }
            }
        }

        let bdeg =
            t.c.values()
                .chain(t.h.values())
                .map(|p| p.num.len().saturating_sub(1))
                .max()
                .unwrap_or(0);
        let laws = if !t.laws_hold() {
            "BROKEN"
        } else if mismatch > 0 {
            "b=0 MISMATCH"
        } else if t.negative.is_empty() {
            "ok"
        } else {
            "NEGATIVE"
        };
        println!(
            "{n:>3} {:>6} {secs:>10.3} {:>9} {:>9} {bdeg:>6} {:>5}  {laws}",
            parts.len(),
            t.c.len(),
            t.h.len(),
            t.peak_bits
        );

        for (tag, key, p) in &t.negative {
            findings.push(format!(
                "{tag}^{}_{{{},{}}} = {}",
                shape(&key.0),
                shape(&key.1),
                shape(&key.2),
                poly(p)
            ));
        }
        assert!(
            t.laws_hold(),
            "a PROVEN law failed at n = {n} — this is a bug"
        );
        assert_eq!(mismatch, 0, "b = 0 is not the class algebra at n = {n}");

        if show && n == top {
            println!();
            println!("-- degree {n}, the nonzero c^lambda_{{mu,nu}}(b) --");
            for (key, p) in &t.c {
                println!(
                    "  c^{:<10} {:<10} {:<10} = {}",
                    shape(&key.0),
                    shape(&key.1),
                    shape(&key.2),
                    poly(p)
                );
            }
            println!();
            println!("-- degree {n}, the nonzero h^lambda_{{mu,nu}}(b) --");
            for (key, p) in &t.h {
                println!(
                    "  h^{:<10} {:<10} {:<10} = {}",
                    shape(&key.0),
                    shape(&key.1),
                    shape(&key.2),
                    poly(p)
                );
            }
        }
    }

    println!();
    if findings.is_empty() {
        println!("Matchings-Jack and the b-conjecture: every coefficient computed lies");
        println!("in N[b]. Both are OPEN; this is evidence, at every degree above.");
    } else {
        println!("*** NEGATIVE COEFFICIENTS — A RESULT, NOT A BUG ***");
        for f in &findings {
            println!("    {f}");
        }
        println!("Check the normalization against docs/spec-jack.md 3.6 first,");
        println!("then report it.");
    }
}
