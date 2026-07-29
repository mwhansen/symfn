//! Oracle-style integration tests.
//!
//! Two kinds of checks:
//!  1. **Known-value table** — LR products whose expansions are textbook /
//!     Sage-verified. This is where a `proptest` harness will later assert
//!     agreement with live Sage output across randomized inputs; for now the
//!     table is a hand-picked, high-signal subset (including a multiplicity-2
//!     coefficient).
//!  2. **Algebraic laws** — commutativity, associativity, and degree additivity
//!     of the Schur product on small inputs. These need no external oracle and
//!     catch structural regressions cheaply.

use symfn::{partitions_of, Partition, Schur, SymFn};

fn s(parts: &[u32]) -> Schur<i64> {
    Schur::monomial(Partition::new(parts.iter().copied()), 1)
}

fn term(parts: &[u32], c: i64) -> (Partition, i64) {
    (Partition::new(parts.iter().copied()), c)
}

/// Assert that `prod` equals exactly the given list of (partition, coeff) terms.
fn assert_expansion(prod: &Schur<i64>, expected: &[(Partition, i64)]) {
    for (p, c) in expected {
        assert_eq!(prod.coeff(p), *c, "coeff of {p}");
    }
    assert_eq!(
        prod.terms().len(),
        expected.iter().filter(|(_, c)| *c != 0).count(),
        "unexpected extra terms in {prod}"
    );
}

#[test]
fn known_schur_products() {
    // s_1 · s_1 = s_2 + s_{11}
    assert_expansion(&s(&[1]).mul(&s(&[1])), &[term(&[2], 1), term(&[1, 1], 1)]);

    // s_2 · s_{11} = s_{31} + s_{211}
    assert_expansion(
        &s(&[2]).mul(&s(&[1, 1])),
        &[term(&[3, 1], 1), term(&[2, 1, 1], 1)],
    );

    // s_{21} · s_{21} = s_{42}+s_{411}+s_{33}+2 s_{321}+s_{3111}+s_{222}+s_{2211}
    assert_expansion(
        &s(&[2, 1]).mul(&s(&[2, 1])),
        &[
            term(&[4, 2], 1),
            term(&[4, 1, 1], 1),
            term(&[3, 3], 1),
            term(&[3, 2, 1], 2),
            term(&[3, 1, 1, 1], 1),
            term(&[2, 2, 2], 1),
            term(&[2, 2, 1, 1], 1),
        ],
    );
}

#[test]
fn schur_product_commutes_on_all_small_pairs() {
    let basis: Vec<Partition> = (0..=3).flat_map(partitions_of).collect();
    for a in &basis {
        for b in &basis {
            let sa = Schur::<i64>::monomial(a.clone(), 1);
            let sb = Schur::<i64>::monomial(b.clone(), 1);
            assert_eq!(sa.mul(&sb), sb.mul(&sa), "s{a} · s{b} not commutative");
        }
    }
}

#[test]
fn schur_product_is_associative() {
    let a = s(&[2]);
    let b = s(&[1]);
    let c = s(&[1]);
    assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)));
}

#[test]
fn degree_adds_under_product() {
    for m in 0..=3u32 {
        for n in 0..=3u32 {
            for mu in partitions_of(m) {
                for nu in partitions_of(n) {
                    let prod =
                        Schur::<i64>::monomial(mu.clone(), 1).mul(&Schur::monomial(nu.clone(), 1));
                    if let Some(d) = prod.degree() {
                        assert_eq!(d, m + n, "deg(s{mu} · s{nu})");
                    }
                }
            }
        }
    }
}
