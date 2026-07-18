//! GMP-backed coefficient tests (`--features gmp`).
//!
//! These demonstrate the payoff of parameterizing on the coefficient ring: the
//! *same* basis, conversion, and Hopf code runs over `rug::Integer` /
//! `rug::Rational` with no changes, and stays exact in ranges where the `i64`
//! scaffold would silently wrap.

#![cfg(feature = "gmp")]

use rug::{Integer, Rational as RugRational};
use symfn::{
    convert, FromSchur, Homogeneous, Monomial, Partition, PowerSum, Ring, Schur, SymFn, ToSchur,
};

fn p(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

#[test]
fn schur_algebra_works_over_gmp_integers() {
    let a: Schur<Integer> = Schur::monomial(p(&[2, 1]), Integer::from(1));
    let prod = a.mul(&a);
    // Same textbook expansion as the i64 path, multiplicity 2 included.
    assert_eq!(prod.coeff(&p(&[3, 2, 1])), Integer::from(2));
    assert_eq!(prod.coeff(&p(&[4, 2])), Integer::from(1));
    assert_eq!(prod.terms().len(), 7);
}

#[test]
fn conversions_round_trip_over_gmp() {
    for parts in [&[2, 1][..], &[3, 1], &[3, 2, 1]] {
        let s: Schur<Integer> = Schur::monomial(p(parts), Integer::from(1));
        let via_h: Schur<Integer> = convert(&Homogeneous::from_schur(&s));
        let via_m: Schur<Integer> = convert(&Monomial::from_schur(&s));
        assert_eq!(via_h, s, "s→h→s over GMP at {:?}", parts);
        assert_eq!(via_m, s, "s→m→s over GMP at {:?}", parts);
    }
}

#[test]
fn power_sums_round_trip_over_gmp_rationals() {
    for parts in [&[2, 1][..], &[3], &[2, 2]] {
        let s: Schur<RugRational> = Schur::monomial(p(parts), RugRational::from(1));
        let ps: PowerSum<RugRational> = PowerSum::from_schur(&s);
        assert_eq!(ps.to_schur(), s, "s→p→s over GMP-ℚ at {:?}", parts);
    }
}

#[test]
fn huge_coefficients_stay_exact() {
    // (2^64)² per factor — squaring this overflows i64 by a wide margin, so an
    // exact result here is real evidence the bignum path carries through.
    let big = Integer::from(u64::MAX) * Integer::from(u64::MAX);
    let a: Schur<Integer> = Schur::monomial(p(&[1]), big.clone());
    let prod = a.mul(&a); // s_1·s_1 = s_2 + s_{11}
    let want = big.clone() * big;
    assert_eq!(prod.coeff(&p(&[2])), want);
    assert_eq!(prod.coeff(&p(&[1, 1])), want);
}

#[test]
fn from_u128_is_exact_beyond_i64_range() {
    let big: u128 = u64::MAX as u128 + 12_345;
    let exact = <Integer as Ring>::from_u128(big);
    assert_eq!(exact.to_string(), big.to_string());
    // The fixed-width path necessarily loses it — which is why the seam exists.
    assert_ne!(<i64 as Ring>::from_u128(big) as u128, big);
}

/// Characters past the `i128` ceiling are exact over GMP.
///
/// |χ^λ(μ)| ≤ d_λ and max d_λ ≈ √(n!), so the fixed-width path tops out near
/// n = 58. Above that `try_character` reports overflow and `character_in`
/// re-runs the recursion in the coefficient ring itself — which is the whole
/// point of the ring being a parameter. Ground truth is the hook-length
/// formula, computed independently in `rug::Integer`.
#[test]
fn characters_beyond_i128_are_exact_over_gmp() {
    use symfn::{character_in, try_character};

    // d_λ = m!/∏hooks, in exact integers — an independent computation, not the
    // Murnaghan–Nakayama recursion being tested.
    fn dimension_by_hooks(lam: &[u32]) -> Integer {
        let m: u32 = lam.iter().sum();
        let mut conj = vec![0i64; lam[0] as usize];
        for &part in lam {
            for c in conj.iter_mut().take(part as usize) {
                *c += 1;
            }
        }
        let mut num = Integer::from(1);
        for k in 2..=m {
            num *= k;
        }
        let mut den = Integer::from(1);
        for (i, &row) in lam.iter().enumerate() {
            for j in 0..row as usize {
                den *= (row as i64 - j as i64) + (conj[j] - i as i64) - 1;
            }
        }
        num / den
    }

    // A shape whose dimension comfortably exceeds i128 (max ≈ 1.7e38).
    let parts: Vec<u32> = vec![18, 16, 14, 12, 10, 8];
    let n: u32 = parts.iter().sum();
    let lam = Partition::new(parts.iter().copied());
    let ones = Partition::new(std::iter::repeat(1).take(n as usize));

    let want = dimension_by_hooks(&parts);
    assert!(
        want > Integer::from(i128::MAX),
        "test is pointless unless it exceeds i128; got {want}"
    );
    // The fixed-width path must *decline*, not wrap.
    assert_eq!(try_character(&lam, &ones), None, "should report overflow");
    // The ring-generic path must be exact.
    let got: Integer = character_in(&lam, &ones);
    assert_eq!(got, want, "χ^{lam}(1^{n}) over GMP");
}
