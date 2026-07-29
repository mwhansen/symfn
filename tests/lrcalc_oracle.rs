//! Cross-validation against **lrcalc** as a second independent oracle.
//!
//! Sage already validates symfn on small degrees (`tests/sage_oracle.rs`).
//! lrcalc earns a separate fixture by reaching much larger shapes, which is
//! exactly where [`SkewLr`] departs from the coefficient-at-a-time backends:
//! it enumerates a skew shape once and bins by content, so a bug there would
//! show up as a wrong *set* of terms, not merely a wrong single coefficient.
//! Small-degree agreement would not catch that; these shapes do.
//!
//! The fixture is committed, so these tests need no lrcalc installed.
//! Regenerate with:
//!   LRCALC=/path/to/lrcalc python3 scripts/gen_lrcalc_oracle.py \
//!     > tests/fixtures/lrcalc_oracle.txt
//!
//! Only lrcalc's output is used; symfn shares no code with it.

use symfn::{expand_skew, AutoLr, LrBackend, NaiveLr, Partition, SkewLr, StripLr};

const FIXTURE: &str = include_str!("fixtures/lrcalc_oracle.txt");

fn parse_partition(s: &str) -> Partition {
    if s.is_empty() {
        return Partition::default();
    }
    Partition::new(s.split(',').map(|x| x.parse::<u32>().expect("part")))
}

/// Parse `PART:COEFF ...` into the sorted expansion symfn returns.
fn parse_expansion(rest: &str) -> Vec<(Partition, u128)> {
    let mut v: Vec<(Partition, u128)> = rest
        .split_whitespace()
        .map(|tok| {
            let (part, coeff) = tok.rsplit_once(':').expect("PART:COEFF");
            (parse_partition(part), coeff.parse::<u128>().expect("coeff"))
        })
        .collect();
    v.sort_by(|a, b| a.0.cmp(&b.0));
    v
}

/// Iterate `(tag, lhs, rhs, expansion)` over fixture lines.
fn lines() -> impl Iterator<Item = (&'static str, Partition, Partition, Vec<(Partition, u128)>)> {
    FIXTURE.lines().filter(|l| !l.trim().is_empty()).map(|l| {
        let mut it = l.splitn(3, ' ');
        let tag = it.next().expect("tag");
        let arg = it.next().expect("LHS|RHS");
        let rest = it.next().unwrap_or("");
        let (a, b) = arg.split_once('|').expect("LHS|RHS");
        (
            tag,
            parse_partition(a),
            parse_partition(b),
            parse_expansion(rest),
        )
    })
}

#[test]
fn fixture_is_present_and_covers_both_operations() {
    let (mut muls, mut skews) = (0, 0);
    for (tag, ..) in lines() {
        match tag {
            "smul" => muls += 1,
            "skew" => skews += 1,
            other => panic!("unknown fixture tag {other:?}"),
        }
    }
    assert!(muls >= 10, "expected a real product sweep, got {muls}");
    assert!(skews >= 8, "expected a real skew sweep, got {skews}");
}

#[test]
fn schur_products_match_lrcalc() {
    for (tag, mu, nu, expected) in lines() {
        if tag != "smul" {
            continue;
        }
        assert_eq!(
            SkewLr.schur_product(&mu, &nu),
            expected,
            "SkewLr: s{mu} * s{nu}"
        );
        assert_eq!(
            AutoLr.schur_product(&mu, &nu),
            expected,
            "AutoLr: s{mu} * s{nu}"
        );
    }
}

#[test]
fn skew_expansions_match_lrcalc() {
    for (tag, lambda, mu, expected) in lines() {
        if tag != "skew" {
            continue;
        }
        assert_eq!(
            expand_skew(&lambda, &mu),
            expected,
            "expand_skew: s_{{{lambda}/{mu}}}"
        );
    }
}

/// The slower backends are checked on the fixture entries they can finish, so
/// agreement is pinned across implementations *and* against lrcalc, not just
/// between the fast path and lrcalc.
#[test]
fn legacy_backends_match_lrcalc_on_small_entries() {
    let mut checked = 0;
    for (tag, mu, nu, expected) in lines() {
        if tag != "smul" || mu.size() + nu.size() > 16 {
            continue;
        }
        assert_eq!(
            NaiveLr.schur_product(&mu, &nu),
            expected,
            "NaiveLr s{mu}*s{nu}"
        );
        assert_eq!(
            StripLr.schur_product(&mu, &nu),
            expected,
            "StripLr s{mu}*s{nu}"
        );
        checked += 1;
    }
    assert!(
        checked >= 6,
        "expected several small entries, got {checked}"
    );
}

/// Individual coefficients read out of the fixture must match `lr_coeff`,
/// including the zeros: every λ of the right size that lrcalc did *not* list
/// must come back as 0.
#[test]
fn individual_coefficients_and_zeros_match_lrcalc() {
    for (tag, mu, nu, expected) in lines() {
        if tag != "smul" || mu.size() + nu.size() > 16 {
            continue;
        }
        for (lambda, c) in &expected {
            assert_eq!(
                SkewLr.lr_coeff(lambda, &mu, &nu),
                *c,
                "c^{lambda}_{{{mu},{nu}}}"
            );
        }
        for lambda in symfn::partitions_of(mu.size() + nu.size()) {
            let want = expected
                .iter()
                .find(|(l, _)| *l == lambda)
                .map_or(0, |(_, c)| *c);
            assert_eq!(
                SkewLr.lr_coeff(&lambda, &mu, &nu),
                want,
                "c^{lambda}_{{{mu},{nu}}}"
            );
        }
    }
}
