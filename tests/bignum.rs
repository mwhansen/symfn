//! Arbitrary-precision coefficient tests (`--features bignum`).
//!
//! These demonstrate the payoff of parameterizing on the coefficient ring: the
//! *same* basis, conversion, and Hopf code runs over `BigInt` / `BigRational`
//! with no changes, and stays exact in ranges where the fixed-width scaffold
//! would silently wrap. That is what makes the escalation path possible: a
//! computation that overflows `i128` can simply be run again here.

#![cfg(feature = "bignum")]

use num_bigint::BigInt;
use num_rational::BigRational;
use symfn::{
    convert, FromSchur, Homogeneous, Monomial, Partition, PowerSum, Ring, Schur, SymFn, ToSchur,
};

fn p(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

#[test]
fn schur_algebra_works_over_bignum_integers() {
    let a: Schur<BigInt> = Schur::monomial(p(&[2, 1]), BigInt::from(1));
    let prod = a.mul(&a);
    // Same textbook expansion as the i64 path, multiplicity 2 included.
    assert_eq!(prod.coeff(&p(&[3, 2, 1])), BigInt::from(2));
    assert_eq!(prod.coeff(&p(&[4, 2])), BigInt::from(1));
    assert_eq!(prod.terms().len(), 7);
}

#[test]
fn conversions_round_trip_over_bignum() {
    for parts in [&[2, 1][..], &[3, 1], &[3, 2, 1]] {
        let s: Schur<BigInt> = Schur::monomial(p(parts), BigInt::from(1));
        let via_h: Schur<BigInt> = convert(&Homogeneous::from_schur(&s));
        let via_m: Schur<BigInt> = convert(&Monomial::from_schur(&s));
        assert_eq!(via_h, s, "s→h→s over bignums at {:?}", parts);
        assert_eq!(via_m, s, "s→m→s over bignums at {:?}", parts);
    }
}

#[test]
fn power_sums_round_trip_over_bignum_rationals() {
    for parts in [&[2, 1][..], &[3], &[2, 2]] {
        let s: Schur<BigRational> = Schur::monomial(p(parts), BigRational::from(BigInt::from(1)));
        let ps: PowerSum<BigRational> = PowerSum::from_schur(&s);
        assert_eq!(ps.to_schur(), s, "s→p→s over bignums-ℚ at {:?}", parts);
    }
}

#[test]
fn huge_coefficients_stay_exact() {
    // (2^64)² per factor — squaring this overflows i64 by a wide margin, so an
    // exact result here is real evidence the bignum path carries through.
    let big = BigInt::from(u64::MAX) * BigInt::from(u64::MAX);
    let a: Schur<BigInt> = Schur::monomial(p(&[1]), big.clone());
    let prod = a.mul(&a); // s_1·s_1 = s_2 + s_{11}
    let want = big.clone() * big;
    assert_eq!(prod.coeff(&p(&[2])), want);
    assert_eq!(prod.coeff(&p(&[1, 1])), want);
}

#[test]
fn from_u128_is_exact_beyond_i64_range() {
    let big: u128 = u64::MAX as u128 + 12_345;
    let exact = <BigInt as Ring>::from_u128(big);
    assert_eq!(exact.to_string(), big.to_string());
    // The fixed-width path necessarily loses it — which is why the seam exists.
    assert_ne!(<i64 as Ring>::from_u128(big) as u128, big);
}

/// Characters past the `i128` ceiling are exact over bignums.
///
/// |χ^λ(μ)| ≤ d_λ and max d_λ ≈ √(n!), so the fixed-width path tops out near
/// n = 58. Above that `try_character` reports overflow and `character_in`
/// re-runs the recursion in the coefficient ring itself — which is the whole
/// point of the ring being a parameter. Ground truth is the hook-length
/// formula, computed independently in `BigInt`.
#[test]
fn characters_beyond_i128_are_exact_over_bignum() {
    use symfn::{character_in, try_character};

    // d_λ = m!/∏hooks, in exact integers — an independent computation, not the
    // Murnaghan–Nakayama recursion being tested.
    fn dimension_by_hooks(lam: &[u32]) -> BigInt {
        let m: u32 = lam.iter().sum();
        let mut conj = vec![0i64; lam[0] as usize];
        for &part in lam {
            for c in conj.iter_mut().take(part as usize) {
                *c += 1;
            }
        }
        let mut num = BigInt::from(1);
        for k in 2..=m {
            num *= k;
        }
        let mut den = BigInt::from(1);
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
        want > BigInt::from(i128::MAX),
        "test is pointless unless it exceeds i128; got {want}"
    );
    // The fixed-width path must *decline*, not wrap.
    assert_eq!(try_character(&lam, &ones), None, "should report overflow");
    // The ring-generic path must be exact.
    let got: BigInt = character_in(&lam, &ones);
    assert_eq!(got, want, "χ^{lam}(1^{n}) over bignums");
}

/// The reduced Kronecker engine past its fixed-width wall.
///
/// `st[8,5] · st[8,5]` (total degree 26) cannot be computed over `i128`: the
/// intermediate rationals carry `z_γ`, which passes 10²⁶ around there, and the
/// guarded path reports that rather than wrapping. With this feature on, the
/// same call re-runs over `BigRational` and answers.
///
/// The answer is checked two ways that do not involve the engine. Its
/// top-degree part must be the Littlewood–Richardson expansion — `ḡ^ν_{λμ} =
/// c^ν_{λμ}` when `|ν| = |λ|+|μ|` — and every coefficient must be non-negative,
/// since reduced Kronecker coefficients count multiplicities. No oracle
/// available anywhere can check this size directly: Sage exceeds 90s at total
/// degree 16.
///
/// **Ignored by default, and the only ignored test in the suite.** Escalation
/// re-runs the whole computation over `BigRational`, which costs ~17s in release
/// and ~220s in the debug profile `cargo test` uses — and there is no cheaper
/// case, because nothing escalates below total degree 26 (measured: `[2,1^6]²`
/// at degree 16 stays fixed-width and takes 0.02s). Run it deliberately:
///
/// ```text
///   cargo test --release --features bignum -- --ignored
/// ```
#[test]
#[ignore = "escalated path: ~17s in release, ~220s in debug"]
fn reduced_kronecker_escalates_past_the_fixed_width_wall() {
    use symfn::{reduced_kronecker_product, St, SymAlgebra};

    let lambda = p(&[8, 5]);
    let mu = p(&[8, 5]);
    let red: St<i128> = reduced_kronecker_product(&lambda, &mu);
    assert!(red.terms().len() > 4000, "got {} terms", red.terms().len());

    for (nu, c) in red.terms() {
        assert!(
            *c >= 0,
            "reduced Kronecker coefficient at {nu} is negative: {c}"
        );
    }

    let sl: Schur<i128> = Schur::monomial(lambda.clone(), 1);
    let sm: Schur<i128> = Schur::monomial(mu.clone(), 1);
    for (nu, c) in sl.mul(&sm).terms() {
        assert_eq!(red.coeff(nu), *c, "top degree at {nu}");
    }

    // The unit still behaves at this size, which exercises the escalated path
    // twice more rather than reusing the cached column.
    let x: St<i128> = St::monomial(lambda.clone(), 1);
    assert_eq!(x.mul(&St::unit()), x);
}
