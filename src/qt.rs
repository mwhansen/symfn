//! Bivariate `(q, t)` polynomials as a coefficient ring.
//!
//! This is the coefficient type the non-classical families need: Hall–Littlewood
//! lives over ℤ[t], Macdonald over ℚ(q,t), Jack over ℚ(α). [`QtPoly`] covers the
//! polynomial half of that — the fraction field is a separate layer on top, and
//! is only needed once Macdonald arrives.
//!
//! ## Why this can exist at all
//!
//! ℚ[q,t] is **not a field**, and until the dividing paths were re-bounded on
//! [`QAlgebra`](crate::coeff::QAlgebra) rather than
//! [`Field`](crate::coeff::Field) they were unavailable over it. Every division
//! in this library is by z_μ — an integer — so a ring containing ℚ suffices.
//! That is what makes `s → p`, the internal product and plethysm work here.
//!
//! ## Sparse, and generic over the coefficients
//!
//! Terms are held as a map from exponent pair to coefficient, with zeros
//! removed, so equality is structural and a polynomial costs what it uses.
//! Hall–Littlewood and Macdonald expansions are sparse in (q, t) — Kostka–
//! Foulkes polynomials in particular have few terms relative to their degree —
//! so a dense representation would mostly store zeros.
//!
//! The coefficient ring is a parameter for the same reason it is everywhere else
//! here: `QtPoly<i64>` is ℤ[q,t] for exact small work, `QtPoly<Rational>` is
//! ℚ[q,t], and `QtPoly<BigInt>` (under `bignum`) has no ceiling. Kostka–Foulkes
//! coefficients are integers, Macdonald's are not, and neither should force the
//! other's representation.
//!
//! ## Plethysm acts on q and t
//!
//! [`Plethystic::frobenius`] raises the *variables*, so `p_n` sends q^a t^b to
//! q^{an} t^{bn}. This is Sage's default convention, and getting it wrong is
//! invisible over ℚ — see [`crate::plethysm`].

use core::fmt;
use std::collections::BTreeMap;

use crate::coeff::{Plethystic, QAlgebra, Ring};

/// A polynomial in `q` and `t`, sparse in the exponent pair `(a, b)` for
/// `q^a t^b`, with coefficients in `C`.
///
/// Never stores a zero coefficient, so `PartialEq` is mathematical equality.
#[derive(Clone, PartialEq, Eq, Debug, Default)]
pub struct QtPoly<C: Ring>(BTreeMap<(u32, u32), C>);

impl<C: Ring> QtPoly<C> {
    /// `c · q^a t^b`.
    pub fn term(a: u32, b: u32, c: C) -> Self {
        let mut m = BTreeMap::new();
        if !c.is_zero() {
            m.insert((a, b), c);
        }
        QtPoly(m)
    }

    /// The variable `q`.
    pub fn q() -> Self {
        Self::term(1, 0, C::one())
    }

    /// The variable `t`.
    pub fn t() -> Self {
        Self::term(0, 1, C::one())
    }

    /// Coefficient of `q^a t^b`.
    pub fn coeff(&self, a: u32, b: u32) -> C {
        self.0.get(&(a, b)).cloned().unwrap_or_else(C::zero)
    }

    /// The terms, ascending by exponent pair.
    pub fn terms(&self) -> impl Iterator<Item = (&(u32, u32), &C)> {
        self.0.iter()
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// `self += c · q^a t^b`, dropping the term if it cancels.
    pub fn add_term(&mut self, a: u32, b: u32, c: C) {
        if c.is_zero() {
            return;
        }
        let key = (a, b);
        let now_zero = {
            let e = self.0.entry(key).or_insert_with(C::zero);
            e.add_assign(&c);
            e.is_zero()
        };
        if now_zero {
            self.0.remove(&key);
        }
    }

    /// Substitute numbers for `q` and `t`.
    ///
    /// The specialisations that matter are exactly this: Hall–Littlewood at
    /// t = 0 is Schur and at t = 1 is monomial, which is how the family gets
    /// checked against bases that already have oracles.
    pub fn eval(&self, q: &C, t: &C) -> C {
        let mut total = C::zero();
        for ((a, b), c) in &self.0 {
            let mut term = c.clone();
            for _ in 0..*a {
                term = term.mul(q);
            }
            for _ in 0..*b {
                term = term.mul(t);
            }
            total.add_assign(&term);
        }
        total
    }

    /// Highest `q` and `t` exponents present, or `None` when zero.
    pub fn degrees(&self) -> Option<(u32, u32)> {
        let a = self.0.keys().map(|k| k.0).max()?;
        let b = self.0.keys().map(|k| k.1).max()?;
        Some((a, b))
    }
}

impl<C: Ring> Ring for QtPoly<C> {
    fn zero() -> Self {
        QtPoly(BTreeMap::new())
    }
    fn one() -> Self {
        Self::term(0, 0, C::one())
    }
    fn is_zero(&self) -> bool {
        self.0.is_empty()
    }
    fn add_assign(&mut self, other: &Self) {
        for (k, c) in &other.0 {
            self.add_term(k.0, k.1, c.clone());
        }
    }
    fn mul(&self, other: &Self) -> Self {
        let mut out = Self::zero();
        for ((a1, b1), c1) in &self.0 {
            for ((a2, b2), c2) in &other.0 {
                out.add_term(a1 + a2, b1 + b2, c1.mul(c2));
            }
        }
        out
    }
    fn neg(&self) -> Self {
        QtPoly(self.0.iter().map(|(k, c)| (*k, c.neg())).collect())
    }
    fn from_i64(n: i64) -> Self {
        Self::term(0, 0, C::from_i64(n))
    }
    fn from_u128(n: u128) -> Self {
        Self::term(0, 0, C::from_u128(n))
    }
    fn from_i128(n: i128) -> Self {
        Self::term(0, 0, C::from_i128(n))
    }

    // `as_ratio` / `from_ratio` are deliberately left declining. They would have
    // to answer in `i128`, which can only represent a *constant* polynomial, and
    // a batch mixing constants with genuine polynomials would have to fall back
    // anyway. Declining keeps `convert::integral_sweep` on its generic path,
    // which is always correct.
}

impl<C: QAlgebra> QAlgebra for QtPoly<C> {
    /// Coefficientwise, and exact: ℚ[q,t] contains ℚ, so dividing by an integer
    /// never needs an inverse of q or t. This is the whole reason the bound is
    /// `QAlgebra` and not `Field` — see the module docs.
    fn div_u128(&self, n: u128) -> Self {
        QtPoly(self.0.iter().map(|(k, c)| (*k, c.div_u128(n))).collect())
    }
}

impl<C: Plethystic> Plethystic for QtPoly<C> {
    /// `p_n` raises the variables: q^a t^b ↦ q^{an} t^{bn}, and the coefficients
    /// are pushed through their own Frobenius.
    ///
    /// Exponents only grow, so the map is injective on monomials and no two
    /// terms can collide — the result is built directly rather than accumulated.
    fn frobenius(&self, n: u32) -> Self {
        QtPoly(
            self.0
                .iter()
                .map(|((a, b), c)| ((a * n, b * n), c.frobenius(n)))
                .collect(),
        )
    }
}

impl<C: Ring> fmt::Display for QtPoly<C> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.0.is_empty() {
            return f.write_str("0");
        }
        let mut parts = Vec::new();
        for ((a, b), c) in &self.0 {
            let mut s = format!("{c:?}");
            if *a > 0 {
                s.push_str(&if *a == 1 { "*q".into() } else { format!("*q^{a}") });
            }
            if *b > 0 {
                s.push_str(&if *b == 1 { "*t".into() } else { format!("*t^{b}") });
            }
            parts.push(s);
        }
        f.write_str(&parts.join(" + "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::convert::{FromSchur, ToSchur};
    use crate::partition::Partition;
    use crate::sym::{PowerSum, Schur, SymFn};

    type P = QtPoly<i64>;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    #[test]
    fn ring_axioms_hold() {
        let q: P = QtPoly::q();
        let t: P = QtPoly::t();
        let one = <P as Ring>::one();
        // (q + t)(q - t) = q^2 - t^2
        let sum = q.add_ring(&t);
        let diff = q.add_ring(&t.neg());
        let prod = sum.mul(&diff);
        assert_eq!(prod.coeff(2, 0), 1);
        assert_eq!(prod.coeff(0, 2), -1);
        assert_eq!(prod.len(), 2, "cross terms must cancel: {prod}");
        // distributivity and identity
        assert_eq!(q.mul(&one), q);
        assert_eq!(q.mul(&sum), q.mul(&q).add_ring(&q.mul(&t)));
        // no explicit zeros are ever stored
        assert!(q.add_ring(&q.neg()).is_zero());
    }

    /// Helper: `Ring` has `add_assign` but no `add`, so tests build one.
    trait AddRing: Ring {
        fn add_ring(&self, other: &Self) -> Self {
            let mut x = self.clone();
            x.add_assign(other);
            x
        }
    }
    impl<T: Ring> AddRing for T {}

    #[test]
    fn division_by_integers_is_exact_and_needs_no_field() {
        let x: QtPoly<Rational> = QtPoly::term(2, 3, Rational::from_int(6));
        let half = x.div_u128(4);
        assert_eq!(half.coeff(2, 3), Rational::new(3, 2));
        // and it is genuinely a division, not a truncation
        assert_eq!(half.mul(&<QtPoly<Rational> as Ring>::from_i64(4)), x);
    }

    /// `frobenius` must be a ring homomorphism and the identity at n = 1 —
    /// the contract `Plethystic` states, and what plethysm relies on.
    ///
    /// Over ℚ[q,t] rather than ℤ[q,t] because `Plethystic: QAlgebra`, and that
    /// is deliberate: plethysm routes through the power-sum basis and so
    /// carries z_μ⁻¹. ℤ[q,t] is a perfectly good ring for *holding*
    /// Hall–Littlewood coefficients and cannot support plethysm, which the
    /// bound says out loud.
    #[test]
    fn frobenius_is_a_ring_homomorphism() {
        type R = QtPoly<Rational>;
        let r = Rational::from_int;
        let a: R = QtPoly::term(1, 2, r(3)).add_ring(&QtPoly::term(0, 1, r(-1)));
        let b: R = QtPoly::term(2, 0, r(5)).add_ring(&<R as Ring>::one());
        for n in 1..=4u32 {
            assert_eq!(
                a.mul(&b).frobenius(n),
                a.frobenius(n).mul(&b.frobenius(n)),
                "multiplicative at n = {n}"
            );
            assert_eq!(
                a.add_ring(&b).frobenius(n),
                a.frobenius(n).add_ring(&b.frobenius(n)),
                "additive at n = {n}"
            );
        }
        assert_eq!(a.frobenius(1), a, "identity at n = 1");
        assert_eq!(a.frobenius(3).coeff(3, 6), r(3), "q t^2 -> q^3 t^6");
    }

    /// Sage's convention, checked on the two-variable case that ℚ[t] alone
    /// cannot distinguish: `s_2[q·t·s_1] = q²t²·s_2`.
    #[test]
    fn plethysm_raises_both_variables() {
        let qt: QtPoly<Rational> = QtPoly::term(1, 1, Rational::from_int(1));
        let inner: Schur<QtPoly<Rational>> = Schur::monomial(part(&[1]), qt);
        let outer: Schur<QtPoly<Rational>> =
            Schur::monomial(part(&[2]), <QtPoly<Rational> as Ring>::one());
        let r = crate::plethysm::plethysm(&outer, &inner);
        assert_eq!(r.terms().len(), 1);
        assert_eq!(r.coeff(&part(&[2])).coeff(2, 2), Rational::from_int(1));
        assert_eq!(r.coeff(&part(&[2])).len(), 1, "no other (q,t) term");
    }

    /// The case a single monomial cannot distinguish: a **sum** in the
    /// coefficient, where the Frobenius has to act on each term and the answer
    /// mixes. Sage:
    ///
    /// ```text
    ///   sage: R.<q,t> = QQ[]; s = SymmetricFunctions(R).schur()
    ///   sage: s[2]((q+t)*s[1])
    ///   q*t*s[1, 1] + (q^2+q*t+t^2)*s[2]
    /// ```
    ///
    /// A frobenius that scaled only the leading term, or that scaled the whole
    /// polynomial by q^n t^n rather than raising each variable, still passes
    /// `plethysm_raises_both_variables` and fails here.
    #[test]
    fn plethysm_of_a_polynomial_coefficient_matches_sage() {
        type R = QtPoly<Rational>;
        let r = Rational::from_int;
        let q_plus_t: R = QtPoly::term(1, 0, r(1)).add_ring(&QtPoly::term(0, 1, r(1)));
        let inner: Schur<R> = Schur::monomial(part(&[1]), q_plus_t);
        let outer: Schur<R> = Schur::monomial(part(&[2]), <R as Ring>::one());
        let res = crate::plethysm::plethysm(&outer, &inner);

        assert_eq!(res.terms().len(), 2);
        let c11 = res.coeff(&part(&[1, 1]));
        assert_eq!(c11.len(), 1);
        assert_eq!(c11.coeff(1, 1), r(1), "q*t on s_11");

        let c2 = res.coeff(&part(&[2]));
        assert_eq!(c2.len(), 3, "q^2 + q t + t^2, got {c2}");
        assert_eq!(c2.coeff(2, 0), r(1));
        assert_eq!(c2.coeff(1, 1), r(1));
        assert_eq!(c2.coeff(0, 2), r(1));
    }

    /// The library's dividing paths must work over ℚ[q,t], which is the point
    /// of the type. `s → p` carries z_μ⁻¹ and would need a `Field` if the bound
    /// had not been fixed.
    #[test]
    fn s_to_p_round_trips_over_qt() {
        for n in 1..=6u32 {
            for lambda in crate::partitions_of(n) {
                let c: QtPoly<Rational> = QtPoly::term(2, 1, Rational::from_int(3));
                let f: Schur<QtPoly<Rational>> = Schur::monomial(lambda.clone(), c);
                let p: PowerSum<QtPoly<Rational>> = PowerSum::from_schur(&f);
                assert_eq!(p.to_schur(), f, "s -> p -> s at {lambda}");
            }
        }
    }

    /// Evaluation is what makes the specialisations checkable: Hall–Littlewood
    /// at t = 0 is Schur and at t = 1 is monomial, so this is the harness those
    /// tests will use.
    #[test]
    fn evaluation_specialises() {
        // 1 + t + t^2, at t = 0 and t = 1
        let f: P = QtPoly::term(0, 0, 1)
            .add_ring(&QtPoly::term(0, 1, 1))
            .add_ring(&QtPoly::term(0, 2, 1));
        assert_eq!(f.eval(&0, &0), 1);
        assert_eq!(f.eval(&0, &1), 3);
        assert_eq!(f.eval(&0, &2), 7);
        // q and t are independent
        let g: P = QtPoly::term(1, 1, 1);
        assert_eq!(g.eval(&5, &7), 35);
        assert_eq!(g.eval(&5, &0), 0);
    }
}
