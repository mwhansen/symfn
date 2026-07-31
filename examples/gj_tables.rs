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
//! Every degree is pinned at **both** known specializations: `b = 0` against
//! the class algebra of `S_n` from characters alone, and `b = 1` against the
//! double coset algebra of `(S_2n, H_n)` by enumerating matchings. Neither
//! touches a Jack polynomial, so together they say [GJ]'s (1) and (5) were
//! transcribed correctly.
//!
//! The closing table reports **how much of the output is not already a
//! theorem**, because most triples fall in cases Goulden–Jackson or
//! Kanunnikov–Promyslov–Vassilieva already proved, and a bulk count of those is
//! not evidence for anything.

use std::time::Instant;

use symfn::gj::{
    double_coset_table, gj_connection_tables, matchings_jack_coverage, BPoly, Coverage,
};
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
        "{:>3} {:>6} {:>10} {:>10} {:>9} {:>9} {:>6} {:>5} {:>7}  laws",
        "n", "p(n)", "modular(s)", "exact(s)", "c terms", "h terms", "b deg", "bits", "primes"
    );

    // Past this degree the exact engine is not run: it is the cross-check, not
    // the workhorse, and its cost is what the modular engine exists to escape.
    let exact_top: u32 = std::env::args()
        .nth(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    // The `b = 1` pin enumerates `(2n−1)!!` matchings per λ: 10,395 at n = 6 and
    // 135,135 at n = 7, so it stops being free quickly.
    let b1_top: u32 = std::env::args()
        .nth(4)
        .and_then(|s| s.parse().ok())
        .unwrap_or(6);

    let mut findings = Vec::new();
    let mut coverage: Vec<(u32, usize, usize, usize)> = Vec::new();
    for n in 1..=top {
        symfn::clear_caches();
        let t0 = Instant::now();
        let t = symfn::gj_connection_tables_modular(n);
        let secs = t0.elapsed().as_secs_f64();

        // While both are affordable, require them to agree entry for entry.
        // They share no arithmetic: one works in Q(alpha) with factored linear
        // denominators, the other never forms a rational function at all.
        let exact = if n <= exact_top {
            symfn::clear_caches();
            let t1 = Instant::now();
            let e = gj_connection_tables(n);
            let es = t1.elapsed().as_secs_f64();
            assert_eq!(e.c, t.c, "engines disagree on c at n = {n}");
            assert_eq!(e.h, t.h, "engines disagree on h at n = {n}");
            format!("{es:>10.3}")
        } else {
            format!("{:>10}", "—")
        };

        // b = 0 must be the class algebra and b = 1 the double coset algebra —
        // the two transcription checks. b = 1 costs (2n−1)!! per λ, so it is run
        // while that is affordable and skipped after.
        let parts = symfn::partitions_of(n);
        let dc = if n <= b1_top {
            Some(double_coset_table(n))
        } else {
            None
        };
        let mut mismatch = 0;
        let mut open_live = 0usize;
        let mut gj_live = 0usize;
        let mut variant_live = 0usize;
        for la in &parts {
            for mu in &parts {
                for nu in &parts {
                    let key = (la.clone(), mu.clone(), nu.clone());
                    let got = t.c.get(&key).map_or((0, 1), BPoly::at_zero);
                    let want = symfn::class_algebra_coefficient(la, mu, nu);
                    if got != (want, 1) {
                        mismatch += 1;
                    }
                    if let Some(dc) = &dc {
                        let at_one =
                            t.c.get(&key)
                                .map_or((0i128, 1u128), |p| (p.num.iter().sum::<i128>(), p.den));
                        if at_one != (*dc.get(&key).unwrap_or(&0) as i128, 1) {
                            mismatch += 1;
                        }
                    }
                    // Coverage is counted on the entries that exist, since a
                    // zero coefficient is not evidence for anything.
                    if t.c.contains_key(&key) {
                        match matchings_jack_coverage(la, mu, nu) {
                            Coverage::Open => open_live += 1,
                            Coverage::ProvedByGj => gj_live += 1,
                            Coverage::SinglePartVariant => variant_live += 1,
                        }
                    }
                }
            }
        }
        coverage.push((n, open_live, gj_live, variant_live));

        let bdeg =
            t.c.values()
                .chain(t.h.values())
                .map(|p| p.num.len().saturating_sub(1))
                .max()
                .unwrap_or(0);
        let laws = if !t.laws_hold() {
            "BROKEN"
        } else if mismatch > 0 {
            "SPECIALIZATION MISMATCH"
        } else if t.negative.is_empty() {
            "ok"
        } else {
            "NEGATIVE"
        };
        println!(
            "{n:>3} {:>6} {secs:>10.3} {exact} {:>9} {:>9} {bdeg:>6} {:>5} {:>7}  {laws}",
            parts.len(),
            t.c.len(),
            t.h.len(),
            t.peak_bits,
            t.primes_used
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
        assert_eq!(
            mismatch, 0,
            "b = 0 is not the class algebra, or b = 1 not the double coset \
             algebra, at n = {n}"
        );

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
    println!("== how much of this is not already a theorem ==");
    println!(
        "{:>3} {:>10} {:>12} {:>12} {:>9}",
        "n", "open", "proved [GJ]", "(n)-variant", "open %"
    );
    for &(n, open, gj, var) in &coverage {
        let total = open + gj + var;
        println!(
            "{n:>3} {open:>10} {gj:>12} {var:>12} {:>8.1}%",
            100.0 * open as f64 / total.max(1) as f64
        );
    }
    println!();
    println!("Coverage, so a count is not mistaken for evidence:");
    println!("  proved [GJ]   lambda = [1^n] or [2,1^(n-2)] -- Goulden and Jackson");
    println!("                built the statistic and proved these outright.");
    println!("  (n)-variant   one of the three partitions is (n). Kanunnikov-Vassilieva");
    println!("                proved mu = nu = (n); with Promyslov, any one of the");
    println!("                three -- but for a VARIATION with labelled matchings,");
    println!("                so this is weaker than the column to its left.");
    println!("  open          everything else. The smallest is lambda = mu = nu = (2,2)");
    println!("                at n = 4. Two- and three-part lambda with neither mu nor");
    println!("                nu equal to (n) is the uncovered region.");
    println!();
    println!("Both specializations are pinned against objects with independent");
    println!("definitions: b = 0 against the class algebra of S_n from characters,");
    println!("b = 1 against the double coset algebra of (S_2n, H_n) by counting");
    println!("matchings. Neither touches a Jack polynomial.");
    println!();

    if findings.is_empty() {
        println!("Matchings-Jack and the b-conjecture: every coefficient computed lies");
        println!("in N[b]. Both are OPEN, and positivity is the open half -- so this is");
        println!("evidence only on the `open` column above.");
        println!();
        println!("A bulk sign check is NOT the interesting output here. Positivity and");
        println!("integrality are already theorems ([DF], [BD]), the degree bound is");
        println!("characterized, and Ben Dali's marginal sums are already known");
        println!("b-positive -- so a counterexample has to hide inside a marginal sum");
        println!("with its siblings cancelling it. The interesting computation is the");
        println!("STATISTIC wt_lambda and how rigid it is, on the open triples, at");
        println!("n <= 9; see docs/research-gaps.md.");
    } else {
        println!("*** NEGATIVE COEFFICIENTS — A RESULT, NOT A BUG ***");
        for f in &findings {
            println!("    {f}");
        }
        println!("Check the normalization against docs/record/jack.md first,");
        println!("then report it.");
    }
}
