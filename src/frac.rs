//! ℚ(q,t) — but only the part of it Macdonald actually inhabits.
//!
//! Macdonald's `P_λ(x; q, t)` has coefficients in the fraction field of ℤ[q,t],
//! which the library has had no way to express: [`QtPoly`] is a ring and the
//! dividing paths ask for [`QAlgebra`](crate::coeff::QAlgebra), a ring
//! *containing* ℚ, which ℚ[q,t] is and ℚ(q,t) needs more than.
//!
//! ## Why this is not a general fraction field
//!
//! A general `Frac<R>` needs a gcd in `R` to stay reduced, and bivariate
//! polynomial gcd is a real algorithm — content and primitive parts over
//! ℤ[q][t], with the coefficient swell that implies. Building it would be the
//! bulk of the work and none of the point.
//!
//! It is also unnecessary. Every denominator Macdonald produces is a product of
//! **binomials `1 − qᵃtᵇ`**: the hook coefficient is
//!
//! ```text
//!   b_λ(s) = (1 − q^{a(s)}   t^{l(s)+1}) / (1 − q^{a(s)+1} t^{l(s)})
//! ```
//!
//! and `ψ_{λ/μ}` is a product of ratios of those. That class is closed under
//! what a fraction field does: the product of two such denominators is one, and
//! so is their **lcm**, which is all addition needs. So the denominator is held
//! *factored* — a multiset of exponent pairs — and never expanded.
//!
//! Cancellation is then trial division of the numerator by a binomial, which is
//! a short exact loop rather than a gcd (see [`divide_by_factor`]).
//!
//! ## Equality is cross-multiplied, deliberately
//!
//! The factored form is **not canonical**, because `1 − qᵃtᵇ` need not be
//! irreducible: `(1 + q)/(1 − q²)` and `1/(1 − q)` are the same element, both
//! fully reduced against whole binomial factors, and structurally different.
//! `PartialEq` therefore cross-multiplies instead of comparing representations.
//! A derived `PartialEq` would silently call equal things unequal.

use std::collections::BTreeMap;

use crate::coeff::{Field, QAlgebra, Ring};
use crate::qt::QtPoly;

/// An element of ℚ(q,t) whose denominator is a product of binomials `1 − qᵃtᵇ`.
///
/// The denominator is a multiset of exponent pairs; `(a, b) ↦ m` means the
/// factor `(1 − qᵃtᵇ)` appears `m` times. The pair `(0, 0)` never appears — that
/// binomial is zero.
#[derive(Clone, Debug)]
pub struct Frac<C: Ring> {
    num: QtPoly<C>,
    den: BTreeMap<(u32, u32), u32>,
}

impl<C: Ring> Frac<C> {
    /// A polynomial, as a fraction with denominator 1.
    pub fn from_poly(num: QtPoly<C>) -> Self {
        let mut f = Frac {
            num,
            den: BTreeMap::new(),
        };
        f.reduce();
        f
    }

    /// `1 / (1 − qᵃtᵇ)`. Panics on `(0, 0)`, which would be `1/0`.
    pub fn inv_factor(a: u32, b: u32) -> Self {
        assert!(a > 0 || b > 0, "1 - q^0 t^0 is zero");
        let mut den = BTreeMap::new();
        den.insert((a, b), 1);
        Frac {
            num: <QtPoly<C> as Ring>::one(),
            den,
        }
    }

    /// The polynomial `1 − qᵃtᵇ`.
    pub fn factor(a: u32, b: u32) -> Self {
        Frac::from_poly(binomial(a, b))
    }

    /// A product of binomial powers, `∏ (1 − qᵃtᵇ)^{m}`, with negative `m`
    /// meaning a denominator factor.
    ///
    /// The fast path for anything built as a product of ratios of binomials —
    /// Macdonald's ψ is exactly that. Multiplying the ratios one at a time with
    /// [`Ring::mul`] expands and re-reduces the numerator at every step; here
    /// the exponents are summed first, so factors appearing on both sides
    /// cancel before anything is expanded at all.
    pub fn from_factors(factors: &BTreeMap<(u32, u32), i32>) -> Self {
        let mut num = <QtPoly<C> as Ring>::one();
        let mut den = BTreeMap::new();
        for (&(a, b), &m) in factors {
            debug_assert!(a > 0 || b > 0, "1 - q^0 t^0 is zero");
            match m.cmp(&0) {
                core::cmp::Ordering::Greater => {
                    for _ in 0..m {
                        num = num.mul(&binomial(a, b));
                    }
                }
                core::cmp::Ordering::Less => {
                    den.insert((a, b), (-m) as u32);
                }
                core::cmp::Ordering::Equal => {}
            }
        }
        let mut f = Frac { num, den };
        f.reduce();
        f
    }

    /// The numerator, and the denominator's factors with their multiplicities.
    pub fn parts(&self) -> (&QtPoly<C>, impl Iterator<Item = (&(u32, u32), &u32)>) {
        (&self.num, self.den.iter())
    }

    /// The denominator, expanded. Only for display and testing — the whole point
    /// of the factored form is not to do this.
    pub fn denominator(&self) -> QtPoly<C> {
        let mut d = <QtPoly<C> as Ring>::one();
        for (&(a, b), &m) in &self.den {
            for _ in 0..m {
                d = d.mul(&binomial(a, b));
            }
        }
        d
    }

    /// Divide out every denominator factor that also divides the numerator.
    ///
    /// Public because [`Ring::add_assign`] deliberately does *not* do it: see
    /// there. A sum accumulated over many terms should be reduced once at the
    /// end, not after every addition.
    pub fn reduce(&mut self) {
        if self.num.is_zero() {
            self.den.clear();
            return;
        }
        self.den.retain(|&(a, b), m| {
            while *m > 0 {
                match divide_by_factor(&self.num, a, b) {
                    Some(q) => {
                        self.num = q;
                        *m -= 1;
                    }
                    None => break,
                }
            }
            *m > 0
        });
    }

    /// `self` rewritten over `target`, which must be a multiple of `self.den`.
    fn lift(&self, target: &BTreeMap<(u32, u32), u32>) -> QtPoly<C> {
        let mut num = self.num.clone();
        for (&(a, b), &m) in target {
            let extra = m - self.den.get(&(a, b)).copied().unwrap_or(0);
            for _ in 0..extra {
                num = num.mul(&binomial(a, b));
            }
        }
        num
    }
}

impl<C: Field> Frac<C> {
    /// Substitute values for `q` and `t`; `None` if the denominator vanishes.
    ///
    /// This is how the specialisations get checked — `q = t` should give the
    /// Schur function, `q = 0` Hall–Littlewood — and the `None` is real: those
    /// are exactly the points where individual coefficients blow up even though
    /// the whole family stays finite.
    pub fn eval(&self, q: &C, t: &C) -> Option<C> {
        let mut d = C::one();
        for (&(a, b), &m) in &self.den {
            let mut f = C::one();
            f.add_assign(&pow(q, a).mul(&pow(t, b)).neg());
            if f.is_zero() {
                return None;
            }
            for _ in 0..m {
                d = d.mul(&f);
            }
        }
        Some(self.num.eval(q, t).div(&d))
    }
}

fn pow<C: Ring>(x: &C, n: u32) -> C {
    let mut acc = C::one();
    for _ in 0..n {
        acc = acc.mul(x);
    }
    acc
}

/// The polynomial `1 − qᵃtᵇ`.
fn binomial<C: Ring>(a: u32, b: u32) -> QtPoly<C> {
    let mut p = <QtPoly<C> as Ring>::one();
    p.add_term(a, b, C::one().neg());
    p
}

/// Exact division by `1 − qᵃtᵇ`, or `None` if it does not divide.
///
/// Multiplying by `qᵃtᵇ` strictly increases the lexicographic key, and [`QtPoly`]
/// keeps its terms sorted by that key, so the lex-least term of the remainder
/// must be a term of the quotient. Peel it off, subtract its multiple of the
/// divisor, repeat: each step raises the least degree, so this terminates.
///
/// The bound is what makes non-divisibility detectable rather than an infinite
/// loop — a quotient term cannot have degree above `deg(N) − (a + b)`.
fn divide_by_factor<C: Ring>(n: &QtPoly<C>, a: u32, b: u32) -> Option<QtPoly<C>> {
    if n.is_zero() {
        return Some(QtPoly::zero());
    }
    let bound = n.terms().map(|(k, _)| k.0 + k.1).max().unwrap();
    // A `BTreeMap` rather than the sorted `Vec` [`QtPoly`] uses: this loop pops
    // the least key and inserts a larger one every step, which on a `Vec` shifts
    // the whole tail twice per term.
    let mut rem: BTreeMap<(u32, u32), C> = n.terms().map(|(k, c)| (*k, c.clone())).collect();
    let mut quot = QtPoly::zero();
    while let Some((&(x, y), c)) = rem.iter().next() {
        let c = c.clone();
        if x + y + a + b > bound {
            return None;
        }
        // The quotient's keys come out ascending, so this only ever appends.
        quot.add_term(x, y, c.clone());
        rem.remove(&(x, y));
        // rem += c·q^{x+a}t^{y+b}
        let slot = rem.entry((x + a, y + b)).or_insert_with(C::zero);
        slot.add_assign(&c);
        if slot.is_zero() {
            rem.remove(&(x + a, y + b));
        }
    }
    Some(quot)
}

/// Cross-multiplied — see the module docs on why the representation is not
/// canonical.
impl<C: Ring> PartialEq for Frac<C> {
    fn eq(&self, other: &Self) -> bool {
        if self.den == other.den {
            return self.num == other.num;
        }
        let mut lcm = self.den.clone();
        for (k, &m) in &other.den {
            let e = lcm.entry(*k).or_insert(0);
            *e = (*e).max(m);
        }
        self.lift(&lcm) == other.lift(&lcm)
    }
}

impl<C: Ring> Eq for Frac<C> {}

impl<C: Ring> Ring for Frac<C> {
    fn zero() -> Self {
        Frac {
            num: QtPoly::zero(),
            den: BTreeMap::new(),
        }
    }
    fn one() -> Self {
        Frac {
            num: <QtPoly<C> as Ring>::one(),
            den: BTreeMap::new(),
        }
    }
    fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
    fn add_assign(&mut self, other: &Self) {
        if other.is_zero() {
            return;
        }
        if self.is_zero() {
            *self = other.clone();
            return;
        }
        // Over the lcm of the two denominators, which stays a product of
        // binomials — the property that makes the factored form usable at all.
        let mut lcm = self.den.clone();
        for (k, &m) in &other.den {
            let e = lcm.entry(*k).or_insert(0);
            *e = (*e).max(m);
        }
        let mut num = self.lift(&lcm);
        num.add_assign(&other.lift(&lcm));
        self.num = num;
        self.den = lcm;
        // Deliberately not reduced. Trial division is the expensive operation
        // here -- 825 of ~2500 profile samples on Macdonald -- and a running
        // sum reduced after every addition pays it once per term for a
        // cancellation that can only be decided once the sum is complete.
        // Callers that accumulate should call `reduce` at the end; correctness
        // does not depend on it, since `is_zero` reads the numerator and
        // equality cross-multiplies.
    }
    fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let mut den = self.den.clone();
        for (k, &m) in &other.den {
            *den.entry(*k).or_insert(0) += m;
        }
        let mut f = Frac {
            num: self.num.mul(&other.num),
            den,
        };
        f.reduce();
        f
    }
    fn neg(&self) -> Self {
        Frac {
            num: self.num.neg(),
            den: self.den.clone(),
        }
    }
    fn from_i64(n: i64) -> Self {
        Frac::from_poly(<QtPoly<C> as Ring>::from_i64(n))
    }
    fn from_u128(n: u128) -> Self {
        Frac::from_poly(<QtPoly<C> as Ring>::from_u128(n))
    }
    fn from_i128(n: i128) -> Self {
        Frac::from_poly(<QtPoly<C> as Ring>::from_i128(n))
    }
}

impl<C: QAlgebra> QAlgebra for Frac<C> {
    fn div_u128(&self, n: u128) -> Self {
        Frac {
            num: self.num.div_u128(n),
            den: self.den.clone(),
        }
    }
}

/// There is deliberately no [`Field`] impl. Inverting a general numerator would
/// put an arbitrary polynomial in the denominator and destroy the closed class
/// the whole design rests on. Macdonald never asks for that — it divides by
/// `b_λ(s)`, a *ratio of binomials*, and [`Frac::ratio`] is that operation.
impl<C: Ring> Frac<C> {
    /// `(1 − q^{a₁}t^{b₁}) / (1 − q^{a₂}t^{b₂})`.
    ///
    /// `(0, 0)` in the numerator means the factor is absent, i.e. a numerator of
    /// 1 — which is how `b_λ(s)` behaves for a cell outside the diagram. In the
    /// denominator it would be division by zero and panics.
    pub fn ratio(a1: u32, b1: u32, a2: u32, b2: u32) -> Self {
        assert!(a2 > 0 || b2 > 0, "1 - q^0 t^0 is zero");
        let num = if a1 == 0 && b1 == 0 {
            <QtPoly<C> as Ring>::one()
        } else {
            binomial(a1, b1)
        };
        let mut den = BTreeMap::new();
        den.insert((a2, b2), 1);
        let mut f = Frac { num, den };
        f.reduce();
        f
    }
}

impl<C: Ring> core::fmt::Display for Frac<C> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "({})", self.num)?;
        for (&(a, b), &m) in &self.den {
            write!(f, "/(1-q^{a}t^{b})")?;
            if m > 1 {
                write!(f, "^{m}")?;
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;

    type F = Frac<Rational>;

    fn r(n: i128) -> Rational {
        Rational::from_int(n)
    }

    #[test]
    fn exact_division_by_a_binomial() {
        // (1 - q^2) = (1 - q)(1 + q)
        let one_minus_q2: QtPoly<Rational> = binomial(2, 0);
        let q = divide_by_factor(&one_minus_q2, 1, 0).expect("(1-q) divides (1-q^2)");
        // quotient is 1 + q
        assert_eq!(q.coeff(0, 0), r(1));
        assert_eq!(q.coeff(1, 0), r(1));
        assert_eq!(q.len(), 2, "{q}");
        // (1 - t) does not divide (1 - q^2)
        assert!(divide_by_factor(&one_minus_q2, 0, 1).is_none());
        // and division is genuinely exact
        assert_eq!(q.mul(&binomial(1, 0)), one_minus_q2);
    }

    /// The case the module docs call out: representations differ, values do not.
    #[test]
    fn equality_is_not_structural() {
        let mut one_plus_q: QtPoly<Rational> = <QtPoly<Rational> as Ring>::one();
        one_plus_q.add_term(1, 0, r(1));
        // (1 + q) / (1 - q^2)
        let a = F::from_poly(one_plus_q).mul(&F::inv_factor(2, 0));
        // 1 / (1 - q)
        let b = F::inv_factor(1, 0);
        assert_eq!(a, b, "{a} must equal {b}");
        // and they really are stored differently, or this test proves nothing
        assert_ne!(a.den, b.den, "the representations must differ");
    }

    #[test]
    fn ring_axioms_hold() {
        let a = F::inv_factor(1, 0); // 1/(1-q)
        let b = F::inv_factor(0, 1); // 1/(1-t)
        let one = <F as Ring>::one();
        assert_eq!(a.mul(&one), a);
        assert_eq!(a.mul(&b), b.mul(&a));
        // a + (-a) = 0, and zero clears the denominator
        let mut z = a.clone();
        z.add_assign(&a.neg());
        assert!(z.is_zero());
        assert_eq!(z, <F as Ring>::zero());
        // distributivity
        let mut sum = a.clone();
        sum.add_assign(&b);
        let mut want = a.mul(&b);
        want.add_assign(&b.mul(&b));
        assert_eq!(sum.mul(&b), want);
    }

    /// 1/(1-q) + 1/(1-t) = (2 - q - t)/((1-q)(1-t)), and the denominator must
    /// be the lcm rather than the product.
    #[test]
    fn addition_uses_the_lcm_not_the_product() {
        let mut s = F::inv_factor(1, 0);
        s.add_assign(&F::inv_factor(1, 0));
        // 1/(1-q) + 1/(1-q) = 2/(1-q): the denominator must not have squared
        let (num, den) = s.parts();
        assert_eq!(num.coeff(0, 0), r(2));
        assert_eq!(num.len(), 1, "{num}");
        let den: Vec<_> = den.collect();
        assert_eq!(den, vec![(&(1, 0), &1)], "denominator must stay (1-q)^1");
    }

    #[test]
    fn reduction_cancels_common_factors() {
        // (1 - q) / (1 - q) = 1
        let f = F::ratio(1, 0, 1, 0);
        assert_eq!(f, <F as Ring>::one());
        let (_, den) = f.parts();
        assert_eq!(den.count(), 0, "the denominator must be gone, not merely equal");
        // (1 - q^2)/(1 - q) = 1 + q, a polynomial
        let g = F::ratio(2, 0, 1, 0);
        let (num, den) = g.parts();
        assert_eq!(den.count(), 0);
        assert_eq!(num.coeff(0, 0), r(1));
        assert_eq!(num.coeff(1, 0), r(1));
    }

    #[test]
    fn evaluation_and_its_poles() {
        let f = F::inv_factor(1, 1); // 1/(1 - q t)
        assert_eq!(f.eval(&r(2), &r(3)), Some(Rational::new(-1, 5)));
        assert_eq!(f.eval(&r(1), &r(1)), None, "q t = 1 is a pole");
        // a pole that cancels is not a pole
        let g = F::ratio(1, 1, 1, 1); // (1-qt)/(1-qt) = 1
        assert_eq!(g.eval(&r(1), &r(1)), Some(r(1)));
    }

    #[test]
    fn division_by_integers_is_exact() {
        let f = F::from_i64(6).mul(&F::inv_factor(1, 0));
        let h = f.div_u128(4);
        assert_eq!(h.eval(&r(0), &r(0)), Some(Rational::new(3, 2)));
    }
}
