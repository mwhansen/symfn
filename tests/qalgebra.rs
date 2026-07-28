//! The library over a coefficient ring that is **not a field**.
//!
//! This is the test the `Field` → [`QAlgebra`] split exists for. Before it,
//! `s → p`, the internal product and plethysm were bounded on `Field`, which
//! demanded far more than the mathematics does — every division in the library
//! is by `z_μ`, a positive integer. ℚ[t] is not a field, so those operations
//! were simply unavailable over it, and they are the operations a Macdonald or
//! Hall–Littlewood coefficient ring needs most.
//!
//! `Poly` below is ℚ[t]. It deliberately does **not** implement `Field`, so if
//! any of these paths regressed to that bound this file would stop compiling —
//! which is the real assertion here, more than any `assert_eq!` in it.
//!
//! Expected values are Sage's, quoted in each test.

use symfn::{
    convert::{FromSchur, ToSchur},
    Partition, Plethystic, PowerSum, QAlgebra, Rational, Ring, Schur, SymFn,
};

// --- ℚ[t] --------------------------------------------------------------------

/// A dense polynomial in t over ℚ, coefficients lowest-degree first, always
/// trimmed of trailing zeros so `PartialEq` is the mathematical one.
#[derive(Clone, PartialEq, Debug)]
struct Poly(Vec<Rational>);

impl Poly {
    fn trim(mut v: Vec<Rational>) -> Self {
        while v.last().map_or(false, |c| c.is_zero()) {
            v.pop();
        }
        Poly(v)
    }
    /// c·t^k
    fn term(c: i128, k: usize) -> Self {
        let mut v = vec![Rational::zero(); k];
        v.push(Rational::from_int(c));
        Poly::trim(v)
    }
}

impl Ring for Poly {
    fn zero() -> Self {
        Poly(Vec::new())
    }
    fn one() -> Self {
        Poly(vec![Rational::one()])
    }
    fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    fn add_assign(&mut self, other: &Self) {
        if other.0.len() > self.0.len() {
            self.0.resize(other.0.len(), Rational::zero());
        }
        for (a, b) in self.0.iter_mut().zip(other.0.iter()) {
            a.add_assign(b);
        }
        *self = Poly::trim(std::mem::take(&mut self.0));
    }
    fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Poly::zero();
        }
        let mut v = vec![Rational::zero(); self.0.len() + other.0.len() - 1];
        for (i, a) in self.0.iter().enumerate() {
            for (j, b) in other.0.iter().enumerate() {
                let p = a.mul(b);
                v[i + j].add_assign(&p);
            }
        }
        Poly::trim(v)
    }
    fn neg(&self) -> Self {
        Poly(self.0.iter().map(|c| c.neg()).collect())
    }
    fn from_i64(n: i64) -> Self {
        Poly::trim(vec![Rational::from_int(n as i128)])
    }
    fn from_u128(n: u128) -> Self {
        Poly::trim(vec![Rational::from_int(n as i128)])
    }
    fn from_i128(n: i128) -> Self {
        Poly::trim(vec![Rational::from_int(n)])
    }
}

impl QAlgebra for Poly {
    /// ℚ[t] contains ℚ, so dividing by an integer is coefficientwise and exact
    /// — no inverse of t is ever needed, which is precisely why `Field` was the
    /// wrong requirement.
    fn div_u128(&self, n: u128) -> Self {
        Poly(self.0.iter().map(|c| c.div_u128(n)).collect())
    }
}

impl Plethystic for Poly {
    /// t ↦ t^n, i.e. spread the coefficients out by a factor of n. This is
    /// Sage's default convention (`p[2](t*p[1]) == t^2*p[2]`).
    fn frobenius(&self, n: u32) -> Self {
        if self.is_zero() {
            return Poly::zero();
        }
        let mut v = vec![Rational::zero(); (self.0.len() - 1) * n as usize + 1];
        for (k, c) in self.0.iter().enumerate() {
            v[k * n as usize] = *c;
        }
        Poly::trim(v)
    }
}

// Note the absence: `impl Field for Poly` does not exist and cannot, since t
// has no inverse here. Every operation exercised below therefore proves its
// bound is no stronger than `QAlgebra`/`Plethystic`.

fn part(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

fn s_poly(v: &[u32], c: Poly) -> Schur<Poly> {
    Schur::monomial(part(v), c)
}

// --- the tests ---------------------------------------------------------------

/// `s → p` over ℚ[t]: the conversion that carries z_μ⁻¹, and the one that was
/// unavailable over any non-field.
///
/// Checked against the ℚ answer: coefficients are polynomial multiples of the
/// rational ones, since the conversion is linear over the coefficient ring and
/// t is inert to it.
#[test]
fn s_to_p_works_over_a_non_field() {
    for deg in 1..=6u32 {
        for lambda in symfn::partitions_of(deg) {
            // (3t² + 1)·s_λ
            let c = Poly::trim(vec![Rational::from_int(1), Rational::zero(), Rational::from_int(3)]);
            let over_qt: PowerSum<Poly> = PowerSum::from_schur(&s_poly(lambda.parts(), c.clone()));

            let over_q: PowerSum<Rational> =
                PowerSum::from_schur(&Schur::monomial(lambda.clone(), Rational::one()));

            assert_eq!(
                over_qt.terms().len(),
                over_q.terms().len(),
                "support differs for {lambda}"
            );
            for (mu, pc) in over_qt.terms() {
                let want = Poly::trim(vec![Rational::one()]).mul(&c).mul(&Poly(vec![over_q.coeff(mu)]));
                assert_eq!(pc, &want, "s_{lambda} -> p at {mu}");
            }
        }
    }
}

/// A full round trip s → p → s over ℚ[t]. The forward leg divides by z_μ and
/// the back leg does not, so this exercises both bounds at once.
#[test]
fn s_to_p_to_s_round_trips_over_a_non_field() {
    for deg in 1..=6u32 {
        for lambda in symfn::partitions_of(deg) {
            let f = s_poly(lambda.parts(), Poly::term(5, 3));
            let p: PowerSum<Poly> = PowerSum::from_schur(&f);
            assert_eq!(p.to_schur(), f, "round trip at {lambda}");
        }
    }
}

/// The internal (Kronecker) product over ℚ[t]. It routes through the power-sum
/// basis in both directions, so it inherits the division twice over.
#[test]
fn kronecker_works_over_a_non_field() {
    // s_[2,1] * s_[2,1] = s_[3] + s_[2,1] + s_[1,1,1]  (Sage: itensor)
    let a = s_poly(&[2, 1], Poly::term(1, 1)); // t·s_{21}
    let r = symfn::internal(&a, &a);
    let t2 = Poly::term(1, 2);
    assert_eq!(r.coeff(&part(&[3])), t2, "t²·s_3");
    assert_eq!(r.coeff(&part(&[2, 1])), t2);
    assert_eq!(r.coeff(&part(&[1, 1, 1])), t2);
    assert_eq!(r.terms().len(), 3);
}

/// Plethysm over ℚ[t], where the coefficient ring is **not** inert.
///
/// `p_n` substitutes into the alphabet and t belongs to it, so the coefficients
/// are raised too. Every value here is Sage's:
///
/// ```text
///   sage: R.<t> = QQ[]; Sym = SymmetricFunctions(R); s = Sym.schur()
///   sage: s[2](t*s[1])            ->  t^2*s[2]
///   sage: s[1,1](t*s[1])          ->  t^2*s[1,1]
///   sage: s[2](t*s[2])            ->  t^2*s[4] + t^2*s[2,2]
/// ```
///
/// The exponent is the tell: a naive implementation that carries coefficients
/// through unchanged returns `t*s[2]` for the first of these, and over ℚ — the
/// only ring available before this change — it would have been right.
#[test]
fn plethysm_raises_coefficient_variables() {
    let t = Poly::term(1, 1);
    let t2 = Poly::term(1, 2);

    let inner = s_poly(&[1], t.clone()); // t·s_1
    let r = symfn::plethysm(&s_poly(&[2], Poly::one()), &inner);
    assert_eq!(r.coeff(&part(&[2])), t2, "s_2[t·s_1] = t²·s_2");
    assert_eq!(r.terms().len(), 1);

    let r = symfn::plethysm(&s_poly(&[1, 1], Poly::one()), &inner);
    assert_eq!(r.coeff(&part(&[1, 1])), t2, "s_11[t·s_1] = t²·s_11");
    assert_eq!(r.terms().len(), 1);

    let r = symfn::plethysm(&s_poly(&[2], Poly::one()), &s_poly(&[2], t));
    assert_eq!(r.coeff(&part(&[4])), t2, "s_2[t·s_2] = t²(s_4 + s_22)");
    assert_eq!(r.coeff(&part(&[2, 2])), t2);
    assert_eq!(r.terms().len(), 2);
}

/// Degree three, where the Frobenius is visible as t ↦ t³ rather than a square
/// that a doubling bug could also produce.
///
/// ```text
///   sage: s[3](t*s[1])   ->  t^3*s[3]
///   sage: s[2,1](t*s[1]) ->  t^3*s[2,1]
/// ```
#[test]
fn plethysm_frobenius_is_the_right_power() {
    let inner = s_poly(&[1], Poly::term(1, 1));
    for outer in [&[3u32][..], &[2, 1], &[1, 1, 1]] {
        let r = symfn::plethysm(&s_poly(outer, Poly::one()), &inner);
        assert_eq!(r.terms().len(), 1, "s_{outer:?}[t·s_1]");
        assert_eq!(
            r.coeff(&part(outer)),
            Poly::term(1, 3),
            "s_{outer:?}[t·s_1] should be t³·s_{outer:?}"
        );
    }
}

/// The integral bases never divide, so they must work over rings that are not
/// ℚ-algebras at all. Nothing here may quietly acquire a `QAlgebra` bound.
#[test]
fn integral_paths_need_no_division() {
    use symfn::{Elementary, Homogeneous, Monomial};
    // i64 is a Ring and neither a Field nor a QAlgebra.
    let f: Schur<i64> = Schur::monomial(part(&[3, 2, 1]), 2);
    let h: Homogeneous<i64> = Homogeneous::from_schur(&f);
    let e: Elementary<i64> = Elementary::from_schur(&f);
    let m: Monomial<i64> = Monomial::from_schur(&f);
    assert_eq!(h.to_schur(), f);
    assert_eq!(e.to_schur(), f);
    assert_eq!(m.to_schur(), f);
    // p → s is integral too, even though s → p is not.
    let p: PowerSum<i64> = PowerSum::monomial(part(&[3, 2, 1]), 1);
    assert!(!p.to_schur().is_zero());
}
