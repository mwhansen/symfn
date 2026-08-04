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
}

/// The other side of the seam. This used to assert that the fixed-width path
/// *truncates* — "necessarily loses it, which is why the seam exists" — and
/// truncating a structure constant is exactly what R8 forbids
/// (`docs/policies/failure.md`): the value that comes back is a wrong LR
/// coefficient or a wrong `z_λ` inside an otherwise exact computation. It now
/// refuses, and the escalation the seam exists for is the caller's answer.
#[test]
#[should_panic(expected = "does not fit i64")]
fn from_u128_past_i64_refuses_rather_than_truncating() {
    let big: u128 = u64::MAX as u128 + 12_345;
    let _ = <i64 as Ring>::from_u128(big);
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
    let ones = Partition::new(std::iter::repeat_n(1, n as usize));

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

/// `AFrac` over the arbitrary-precision rings.
///
/// This test exists because of the warning in `coeff.rs`: a trait bound is only
/// checked where it is *instantiated*, so `cargo build` can be green while a
/// coefficient type nothing in the crate constructs is broken. `AFrac<C>` leans
/// on [`Ring::div_exact`] for every cancellation, and `BigInt` implements it
/// exactly while `BigRational` implements it as a field — two different
/// meanings, both correct here, and neither exercised by the default suite.
///
/// The Python boundary escalates to `AFrac<BigInt>`, so a silent failure here
/// would be a silent failure there.
#[test]
fn jack_runs_over_bignum_coefficients() {
    use symfn::afrac::AFrac;

    // The convention gate, at both widths: J_(2) = (α+1)m_2 + 2m_11.
    let two = p(&[2]);
    let big: Monomial<AFrac<BigInt>> = symfn::jack_j(&two);
    assert_eq!(big.coeff(&two), AFrac::linear(1, 1));
    assert_eq!(big.coeff(&p(&[1, 1])), <AFrac<BigInt> as Ring>::from_i64(2));

    // Whole degrees must agree term for term with the fixed-width run — the
    // two-width ladder, at the widths that cannot wrap.
    for n in 0..=6u32 {
        for lambda in symfn::partitions_of(n) {
            let narrow: Monomial<AFrac<i128>> = symfn::jack_p(&lambda);
            let wide: Monomial<AFrac<BigInt>> = symfn::jack_p(&lambda);
            let exact: Monomial<AFrac<BigRational>> = symfn::jack_p(&lambda);
            assert_eq!(narrow.terms().len(), wide.terms().len(), "at {lambda}");
            assert_eq!(narrow.terms().len(), exact.terms().len(), "at {lambda}");
            for (mu, c) in narrow.terms() {
                // Compare through the numerator/denominator shape rather than
                // the values, since the three rings are different types: the
                // atoms are canonical, so they must match exactly.
                let (n_num, n_den, n_scale) = c.parts();
                let w = wide.coeff(mu);
                let (w_num, w_den, w_scale) = w.parts();
                assert_eq!(n_num.len(), w_num.len(), "degree at {lambda}/{mu}");
                assert_eq!(
                    n_den.map(|(k, m)| (*k, *m)).collect::<Vec<_>>(),
                    w_den.map(|(k, m)| (*k, *m)).collect::<Vec<_>>(),
                    "atoms at {lambda}/{mu}"
                );
                assert_eq!(n_scale, w_scale, "scalar at {lambda}/{mu}");
                for (a, b) in n_num.iter().zip(w_num.iter()) {
                    assert_eq!(BigInt::from(*a), *b, "coefficient at {lambda}/{mu}");
                }
            }
        }
    }

    // And [KS] Thm 1.1 still holds over BigInt: J is integral, so the
    // denominator must clear entirely.
    for lambda in symfn::partitions_of(6) {
        let j: Monomial<AFrac<BigInt>> = symfn::jack_j(&lambda);
        for (mu, c) in j.terms() {
            assert!(
                c.clone().into_poly().is_some(),
                "J_{lambda} at {mu} kept a denominator over BigInt"
            );
        }
    }
}

/// The one `(q,t)` wall a caller reaches casually, and the ladder that crosses
/// it.
///
/// `llt_h` at μ = 1ⁿ, k = 3 leaves `i128` at n = 87 — in about a second, which
/// is what makes it worth a mechanism rather than a documented refusal
/// (`docs/record/failure-and-overflow.md`). This pins the two halves of the
/// contract: the fixed-width pass **reports** rather than answering, and the
/// wide pass produces a coefficient that provably could not have fitted.
///
/// Release-only, the way `tests/memory.rs` is: the same walk costs 78 s in a
/// debug build against 9 s here, and CI's release lane is where it runs.
#[test]
fn llt_h_reports_at_the_wall_and_escalates_past_it() {
    use num_traits::Signed;
    use symfn::{guarded, Guarded, Partition, SymFn};

    if cfg!(debug_assertions) {
        eprintln!("skipped: the n = 87 wall costs 78 s in a debug build");
        return;
    }

    let mu = Partition::new(std::iter::repeat_n(1, 87));

    // The fast pass must not hand back a `Some`: at this degree an answer
    // coefficient passes i128::MAX, and a `Some` here would be the wrapped
    // value the whole guard exists to prevent.
    assert!(
        guarded(|| symfn::llt::llt_h::<Guarded>(&mu, 3)).is_none(),
        "the fixed-width pass must report at n = 87"
    );

    // And the wide pass answers, with something that did not fit.
    let wide = symfn::llt::llt_h::<BigInt>(&mu, 3);
    let ceiling = BigInt::from(i128::MAX);
    let widest = wide
        .terms()
        .values()
        .flat_map(|p| p.terms().map(|(_, c)| c.abs()))
        .max()
        .expect("H^(3)_{1^87} is not empty");
    assert!(
        widest > ceiling,
        "n = 87 should exceed i128::MAX; widest was {widest}"
    );
}

/// R6: a guarded scope reports at the z_μ wall instead of panicking through it.
///
/// `powersum_scalar` needs z_μ, and used to form it with `Partition::z`, which
/// accumulates in native `u128`. Inside `guarded` that is the one failure the
/// escalation ladder cannot handle: past |μ| = 34 the native multiply overflows
/// and **panics**, while the ladder is watching for `None`. The window is
/// sharp — at |μ| = 34, z_μ = 34! still fits `u128` and only the injection into
/// `i128` refuses, so it reported correctly and the bug began one degree later.
///
/// Reachable from `jack_scalar` and `jack_structure_constant`, both of which
/// wrap their fast pass in `guarded`, and cheaply from the Rust API — this test
/// runs in milliseconds. Fixed by accumulating z_μ in the coefficient ring
/// (`Partition::z_in`), which `z`'s own docs already name as the escape for a
/// caller multiplying *by* z_λ.
#[test]
fn the_z_wall_reports_inside_a_guarded_scope_rather_than_panicking() {
    use symfn::guard::{guarded, Guarded};
    use symfn::{AFrac, PowerSum};

    let factorial = |n: u32| -> BigInt { (1..=n).map(BigInt::from).product() };

    for n in [33u32, 34, 35, 40] {
        let mu = p(&vec![1u32; n as usize]);

        let mut fast: PowerSum<AFrac<Guarded>> = PowerSum::zero();
        fast.add_term(mu.clone(), <AFrac<Guarded> as Ring>::one());
        // The assertion is that this *returns* at all past n = 34. Before the
        // fix it panicked, which no `escalate` can catch.
        let narrow = guarded(|| symfn::powersum_scalar(&fast, &fast));

        let mut wide: PowerSum<AFrac<BigInt>> = PowerSum::zero();
        wide.add_term(mu.clone(), <AFrac<BigInt> as Ring>::one());
        let exact = symfn::powersum_scalar(&wide, &wide);

        // ⟨p_{1ⁿ}, p_{1ⁿ}⟩ = z_{1ⁿ} = n!, which the wide pass must carry exactly.
        assert!(
            format!("{exact:?}").contains(&factorial(n).to_string()),
            "the wide pass lost z_{{1^{n}}} = {n}!"
        );

        if n <= 33 {
            assert!(
                narrow.is_some(),
                "z_{{1^{n}}} = {n}! fits, so the fast pass should answer"
            );
        } else {
            assert!(
                narrow.is_none(),
                "z_{{1^{n}}} = {n}! is past i128, so the fast pass must report and escalate"
            );
        }
    }
}

/// The `s → J` table is the same over `BigRational` as over `Rational`.
///
/// `schur_in_macdonald_j` runs the fixed-width ring first and re-runs here when
/// the guard reports, so this is the branch nothing at a reachable degree
/// exercises. Compared as numerator terms and denominator factors rather than
/// through a common ring, because there is no common ring — which is the point.
#[test]
fn the_macdonald_j_inverse_agrees_over_bignum_rationals() {
    for n in 0..=5u32 {
        let small = symfn::schur_in_j_table::<symfn::Rational>(n);
        let big = symfn::schur_in_j_table::<BigRational>(n);
        for (i, row) in small.iter().enumerate() {
            for (j, c) in row.iter().enumerate() {
                let (num, den) = c.parts();
                let (wide_num, wide_den) = big[i][j].parts();
                assert_eq!(
                    den.map(|(&e, &m)| (e, m)).collect::<Vec<_>>(),
                    wide_den.map(|(&e, &m)| (e, m)).collect::<Vec<_>>(),
                    "the denominator of w_[{i}][{j}] at degree {n}"
                );
                let narrow: Vec<_> = num
                    .terms()
                    .map(|(&e, v)| (e, BigInt::from(v.numer()), BigInt::from(v.denom())))
                    .collect();
                let wide: Vec<_> = wide_num
                    .terms()
                    .map(|(&e, v)| (e, v.numer().clone(), v.denom().clone()))
                    .collect();
                assert_eq!(narrow, wide, "the numerator of w_[{i}][{j}] at degree {n}");
            }
        }
    }
}
