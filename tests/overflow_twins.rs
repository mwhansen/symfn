//! Every entry point that panics when a value leaves its fixed width has a
//! twin that returns `None` there instead, and the two agree everywhere else.
//!
//! The pairs: `Partition::z` / `try_z`, `Ring::from_u128` / `try_from_u128`
//! and the `i128` pair, `character` / `try_character`, `kostka` / `try_kostka`,
//! `class_algebra_coefficient` / `try_class_algebra_coefficient`, and under
//! `bignum` `kronecker_coeff` / `try_kronecker_coeff` (`tests/bignum.rs`).
//! `principal_specialization_q` and `schubert::dimension` changed shape
//! instead and return `Option` like their siblings. The decision is recorded
//! in `docs/record/failure-and-overflow.md`.

use symfn::guard::overflow_count;
use symfn::permutation::Perm;
use symfn::schubert::{dimension, schubert_monomial_mass_of};
use symfn::{
    class_algebra_coefficient, kostka, partitions_of, principal_specialization,
    principal_specialization_q, try_class_algebra_coefficient, try_kostka, Guarded, GuardedRat,
    Partition, Rational, Ring,
};

/// 34! is the last factorial in `u128`; z_{1^n} = n! is where `z` walls.
#[test]
fn try_z_declines_exactly_where_z_panics() {
    let last: u128 = (1..=34u128).product();
    assert_eq!(Partition::new([1; 34]).try_z(), Some(last));
    assert_eq!(Partition::new([1; 34]).z(), last);
    assert_eq!(Partition::new([1; 35]).try_z(), None);
    assert_eq!(Partition::new([2, 1]).try_z(), Some(2));
    assert_eq!(Partition::default().try_z(), Some(1));
    for lam in (0..=8).flat_map(partitions_of) {
        assert_eq!(lam.try_z(), Some(lam.z()), "z at {lam}");
    }
}

#[test]
#[should_panic(expected = "does not fit u128")]
fn z_panics_at_the_first_factorial_past_u128() {
    let _ = Partition::new([1; 35]).z();
}

#[test]
fn try_from_u128_declines_exactly_where_from_u128_panics() {
    let past_i64 = u128::try_from(i64::MAX).unwrap() + 1;
    let past_i128 = u128::try_from(i128::MAX).unwrap() + 1;
    assert_eq!(<i64 as Ring>::try_from_u128(past_i64 - 1), Some(i64::MAX));
    assert_eq!(<i64 as Ring>::try_from_u128(past_i64), None);
    assert_eq!(
        <i128 as Ring>::try_from_u128(past_i64),
        Some(i128::from(i64::MAX) + 1)
    );
    assert_eq!(<i128 as Ring>::try_from_u128(past_i128), None);
    assert_eq!(
        <Rational as Ring>::try_from_u128(past_i64),
        Some(Rational::new(i128::from(i64::MAX) + 1, 1))
    );
    assert_eq!(<Rational as Ring>::try_from_u128(past_i128), None);
    assert_eq!(
        <i64 as Ring>::try_from_i128(i128::from(i64::MIN)),
        Some(i64::MIN)
    );
    assert_eq!(<i64 as Ring>::try_from_i128(i128::from(i64::MIN) - 1), None);
    assert_eq!(<i128 as Ring>::try_from_i128(i128::MIN), Some(i128::MIN));
}

#[test]
#[should_panic(expected = "does not fit i64")]
fn from_u128_on_i64_panics_naming_the_ring() {
    let _ = <i64 as Ring>::from_u128(u128::try_from(i64::MAX).unwrap() + 1);
}

#[test]
#[should_panic(expected = "does not fit i128")]
fn from_u128_on_i128_panics_naming_the_ring() {
    let _ = <i128 as Ring>::from_u128(u128::MAX);
}

/// A `try_` asks rather than acts, so the guard's counter stays where it was.
#[test]
fn the_guarded_twins_decline_without_reporting() {
    let before = overflow_count();
    assert_eq!(<Guarded as Ring>::try_from_u128(u128::MAX), None);
    assert_eq!(<GuardedRat as Ring>::try_from_u128(u128::MAX), None);
    assert_eq!(<GuardedRat as Ring>::try_from_i128(i128::MIN), None);
    assert_eq!(
        <Guarded as Ring>::try_from_i128(i128::MIN),
        Some(Guarded::from_i128(i128::MIN))
    );
    assert_eq!(
        overflow_count(),
        before,
        "a try_ moved the overflow counter"
    );
}

/// K_{λμ} through degree 6, both forms; the `None` side needs n ≈ 58 and is
/// not a value this suite can reach.
#[test]
fn try_kostka_agrees_with_kostka_below_the_wall() {
    for n in 0..=6 {
        for lam in partitions_of(n) {
            for mu in partitions_of(n) {
                assert_eq!(
                    try_kostka(&lam, &mu),
                    Some(kostka(&lam, &mu)),
                    "K at {lam},{mu}"
                );
            }
        }
    }
    assert_eq!(
        try_kostka(&Partition::new([2]), &Partition::new([1])),
        Some(0)
    );
}

/// The leading `n!` leaves `i128` at n = 34, before any character is formed.
#[test]
fn try_class_algebra_coefficient_declines_at_the_factorial_wall() {
    let t = Partition::new([2, 1]);
    assert_eq!(
        try_class_algebra_coefficient(&Partition::new([1, 1, 1]), &t, &t),
        Some(3)
    );
    assert_eq!(
        try_class_algebra_coefficient(&Partition::new([1, 1, 1]), &t, &t),
        Some(class_algebra_coefficient(
            &Partition::new([1, 1, 1]),
            &t,
            &t
        ))
    );
    assert_eq!(
        try_class_algebra_coefficient(&Partition::new([2]), &t, &t),
        Some(0)
    );
    let n34 = Partition::new([34]);
    assert_eq!(try_class_algebra_coefficient(&n34, &n34, &n34), None);
}

#[test]
#[should_panic(expected = "leaves i128")]
fn class_algebra_coefficient_panics_at_the_factorial_wall() {
    let n34 = Partition::new([34]);
    let _ = class_algebra_coefficient(&n34, &n34, &n34);
}

/// The q-analogue declines where its inversion leaves `i128`, and answers
/// everywhere its content-formula sibling does.
#[test]
fn the_q_specialization_declines_past_i128_like_its_sibling() {
    let s21 = Partition::new([2, 1]);
    assert_eq!(
        principal_specialization_q(&s21, 3),
        Some(vec![0, 1, 2, 2, 2, 1])
    );
    assert_eq!(principal_specialization_q(&s21, 1), Some(Vec::new()));
    assert_eq!(
        principal_specialization_q(&Partition::default(), 0),
        Some(vec![1])
    );
    for lam in (1..=5).flat_map(partitions_of) {
        for n in 0..=4 {
            let total: i128 = principal_specialization_q(&lam, n).unwrap().iter().sum();
            assert_eq!(
                u128::try_from(total).ok(),
                principal_specialization(&lam, n),
                "q = 1 at {lam}, n = {n}"
            );
        }
    }
    let row = Partition::new([80]);
    assert_eq!(principal_specialization_q(&row, 80), None);
}

#[test]
fn schubert_dimension_declines_rather_than_saturating() {
    let w = Perm::from_code(&[2, 1, 0]);
    assert_eq!(dimension(&w), Some(1));
    assert_eq!(dimension(&Perm::identity()), Some(1));
    let u = Perm::from_code(&[1, 2, 0]);
    assert_eq!(
        schubert_monomial_mass_of(&u, &u),
        dimension(&u).unwrap().pow(2)
    );
}
