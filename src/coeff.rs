//! The coefficient-ring abstraction: [`Ring`], then [`QAlgebra`] and
//! [`Plethystic`] for the two things a ring may additionally have to supply.
//!
//! Every symmetric-function type is generic over `C: Ring`, so the *same* basis
//! code works over machine integers today and over GMP bignums (`rug::Integer`)
//! or a `(q,t)`-polynomial ring tomorrow — the coefficient ring is a parameter,
//! never baked in.
//!
//! The layering is deliberate and load-bearing, and each step up is demanded by
//! exactly one thing:
//!
//! | bound | what needs it | why |
//! |---|---|---|
//! | [`Ring`] | s, h, e, m, f and their conversions; products; skewing | exact over ℤ |
//! | [`QAlgebra`] | `s → p`, internal product, plethysm | z_μ⁻¹ — division by an **integer** |
//! | [`Plethystic`] | plethysm | `p_n` acts on the coefficients too |
//!
//! [`Field`] appears in none of those rows, which is the point. It was the
//! bound on the dividing paths, and it was too strong: they divide only by z_μ,
//! so a ring containing ℚ suffices and need not invert its own elements. ℚ[t]
//! and ℚ[q,t] are the cases that matter — neither is a field, both are fine —
//! and they are precisely the coefficient rings Hall–Littlewood and Macdonald
//! need. `Field` is kept because it is a real thing to name and [`Rational`] is
//! one, but nothing in the library requires it.
//!
//! Implementors: `i64`/`i128` are rings only; [`Rational`] and (under the `gmp`
//! feature) `rug::Rational` implement all of them, `rug::Integer` is a ring.
//! `tests/qalgebra.rs` carries a ℚ[t] implementing the two upper traits and
//! deliberately **not** `Field` — which is what keeps the dividing paths from
//! quietly drifting back to the stronger bound.
//!
//! A trait bound is only checked where it is instantiated, so adding a bound
//! can pass `cargo build` and still break a coefficient type nothing in the
//! crate constructs. `cargo test --features gmp` is what catches that; it is
//! not part of the default test run.

/// A commutative ring usable as a symmetric-function coefficient.
///
/// Intentionally small: enough for additive structure, multiplication, and the
/// subtraction that determinant / Newton-identity code needs.
pub trait Ring: Clone + PartialEq + core::fmt::Debug {
    fn zero() -> Self;
    fn one() -> Self;
    fn is_zero(&self) -> bool;
    /// `self += other`
    fn add_assign(&mut self, other: &Self);
    /// `self * other`
    fn mul(&self, other: &Self) -> Self;
    /// The additive inverse `-self`.
    fn neg(&self) -> Self;
    /// Injection of a (small) integer, used for structure constants.
    fn from_i64(n: i64) -> Self;

    /// Injection of an unsigned integer. Structure constants — Littlewood–
    /// Richardson coefficients, Kostka numbers, z_λ — are naturally `u128`, and
    /// routing them through [`Ring::from_i64`] would silently truncate the large
    /// ones. Fixed-width coefficient types necessarily lose range here, but a
    /// bignum type (`rug::Integer` under the `gmp` feature) overrides this to be
    /// exact — which is the whole point of the seam.
    fn from_u128(n: u128) -> Self {
        Self::from_i64(n as i64)
    }

    /// This value as an exact ratio of `i128`s, if it is one.
    ///
    /// The seam for putting a batch of coefficients over a **common
    /// denominator**, so that a sum of (coefficient × integer) terms can be
    /// accumulated in `i128` and converted back once at the end instead of
    /// doing rational arithmetic per term. `p → s` does exactly that shape of
    /// sum, and profiling put 55% of plethysm's runtime in the gcds it implies.
    ///
    /// Default `None`: a ring that cannot answer simply keeps the generic path,
    /// which is always correct. Returning `Some` is a promise that
    /// [`Ring::from_ratio`] inverts it exactly.
    fn as_ratio(&self) -> Option<(i128, i128)> {
        None
    }

    /// Exact `num/den`, inverting [`Ring::as_ratio`]. `None` if not
    /// representable — the caller then falls back rather than rounding.
    fn from_ratio(_num: i128, _den: i128) -> Option<Self> {
        None
    }

    /// Injection of a *signed* wide integer. Symmetric-group characters are the
    /// motivating case: |χ^λ(μ)| ≤ d_λ and max d_λ ≈ √(n!), which passes `i64`
    /// at n ≈ 35, so routing them through [`Ring::from_i64`] would silently
    /// truncate. Same seam as [`Ring::from_u128`], for values that can be
    /// negative.
    fn from_i128(n: i128) -> Self {
        Self::from_i64(n as i64)
    }

    /// `self -= other`, provided via [`Ring::neg`].
    fn sub_assign(&mut self, other: &Self) {
        let neg = other.neg();
        self.add_assign(&neg);
    }
}

/// A ring in which every nonzero element is invertible.
///
/// Note this is **not** what the library's dividing paths require — see
/// [`QAlgebra`], which is weaker and is what they are actually bounded by. A
/// field is still a useful thing to name, and [`Rational`] is one.
pub trait Field: Ring {
    /// The multiplicative inverse `self⁻¹`. Panics if `self` is zero.
    fn inv(&self) -> Self;

    /// `self / other`, provided via [`Field::inv`].
    fn div(&self, other: &Self) -> Self {
        self.mul(&other.inv())
    }
}

/// A ring containing ℚ — equivalently, one in which every nonzero *integer* is
/// invertible.
///
/// **This, and not [`Field`], is what dividing in this library actually needs**,
/// and the difference is the whole point of the trait existing. Every division
/// the crate performs is by `z_μ`, a positive integer: `s → p` carries z_μ⁻¹,
/// and the internal product and plethysm inherit it by routing through the
/// power-sum basis. Nothing ever divides by a general ring element.
///
/// Bounding those paths on `Field` therefore demanded far more than the
/// mathematics does, and it excluded exactly the rings this library most wants
/// to serve. ℚ[t] is not a field, nor is ℚ[q,t] — but z_μ⁻¹ lives in both, so
/// `s → p` over them is perfectly well defined and was simply unavailable. The
/// same applies to any ℚ-algebra a caller brings across the Sage boundary, and
/// to the (q,t)-coefficient rings the Macdonald and Hall–Littlewood work needs.
///
/// A field is a ℚ-algebra as soon as it has characteristic 0, but the
/// implication is deliberately *not* written as a blanket impl: that would
/// occupy the impl for every downstream type, and a polynomial ring — the case
/// this exists for — could then never implement it. Implementors state both.
pub trait QAlgebra: Ring {
    /// `self / n` for a positive integer `n`.
    ///
    /// The contract is exactness: `n` is invertible by assumption, so this
    /// neither rounds nor fails. Panics if `n` is zero, matching
    /// [`Field::inv`].
    fn div_u128(&self, n: u128) -> Self;
}

/// A coefficient ring that knows how plethysm acts on **its own elements**.
///
/// Plethysm is the one operation in this library that cannot treat the
/// coefficient ring as inert. Everything else — products, basis changes,
/// skewing — is linear with structure constants that are plain integers, so a
/// coefficient is only ever multiplied and added. `p_n[·]` is different: it
/// substitutes into the *alphabet*, and a coefficient ring with variables in it
/// is part of that alphabet. The standard convention raises them:
///
/// ```text
///   p_n[t · p_1] = t^n · p_n            (over ℚ[t])
/// ```
///
/// which is what Sage computes, and what `exclude=[t]` there turns off.
///
/// So a ring is usable for plethysm only once it says which map that is, and
/// this trait is that statement. For ℚ there is nothing to raise and it is the
/// identity; for ℚ[t] it is t ↦ t^n; for a ring whose variables should be held
/// *constant* it is the identity again, which is a real convention choice the
/// implementor makes rather than a default this crate can pick.
///
/// Splitting it from [`QAlgebra`] keeps the bound honest: `s → p` divides but
/// never substitutes, so it must not demand this.
pub trait Plethystic: QAlgebra {
    /// The nth plethystic Frobenius: raise every variable of `self` to the nth
    /// power, fixing the constants. Must be a ring homomorphism, and must be
    /// the identity when `n == 1`.
    fn frobenius(&self, n: u32) -> Self;
}

macro_rules! impl_ring_for_int {
    ($($t:ty),*) => {$(
        impl Ring for $t {
            #[inline] fn zero() -> Self { 0 }
            #[inline] fn one() -> Self { 1 }
            #[inline] fn is_zero(&self) -> bool { *self == 0 }
            #[inline] fn add_assign(&mut self, other: &Self) { *self += *other; }
            #[inline] fn mul(&self, other: &Self) -> Self { *self * *other }
            #[inline] fn neg(&self) -> Self { -*self }
            #[inline] fn from_i64(n: i64) -> Self { n as $t }
            #[inline] fn from_u128(n: u128) -> Self { n as $t }
            #[inline] fn from_i128(n: i128) -> Self { n as $t }
        }
    )*};
}

impl_ring_for_int!(i64, i128);

// ----------------------------------------------------------------------------
// Rational — an exact field over i128, always kept in lowest terms with den > 0.
// ----------------------------------------------------------------------------

/// An exact rational number `num/den`, normalized so `den > 0` and
/// `gcd(|num|, den) == 1`. The field used wherever symmetric-function
/// coefficients leave ℤ (power sums, inner products).
#[derive(Clone, Copy, PartialEq, Eq)]
pub struct Rational {
    num: i128,
    den: i128,
}

fn gcd(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

impl Rational {
    /// Construct `num/den` in lowest terms. Panics if `den == 0`.
    pub fn new(num: i128, den: i128) -> Self {
        assert!(den != 0, "Rational with zero denominator");
        let mut n = num;
        let mut d = den;
        if d < 0 {
            n = -n;
            d = -d;
        }
        if n == 0 {
            return Rational { num: 0, den: 1 };
        }
        let g = gcd(n, d);
        Rational {
            num: n / g,
            den: d / g,
        }
    }

    /// The integer `n` as `n/1`.
    pub fn from_int(n: i128) -> Self {
        Rational { num: n, den: 1 }
    }

    pub fn numer(&self) -> i128 {
        self.num
    }
    pub fn denom(&self) -> i128 {
        self.den
    }
}

impl core::fmt::Debug for Rational {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if self.den == 1 {
            write!(f, "{}", self.num)
        } else {
            write!(f, "{}/{}", self.num, self.den)
        }
    }
}

impl Ring for Rational {
    fn zero() -> Self {
        Rational { num: 0, den: 1 }
    }
    fn one() -> Self {
        Rational { num: 1, den: 1 }
    }
    fn is_zero(&self) -> bool {
        self.num == 0
    }
    fn add_assign(&mut self, other: &Self) {
        // a/b + c/d = (ad + cb) / bd, then normalize.
        let num = self.num * other.den + other.num * self.den;
        let den = self.den * other.den;
        *self = Rational::new(num, den);
    }
    fn mul(&self, other: &Self) -> Self {
        Rational::new(self.num * other.num, self.den * other.den)
    }
    fn neg(&self) -> Self {
        Rational {
            num: -self.num,
            den: self.den,
        }
    }
    fn from_i64(n: i64) -> Self {
        Rational::from_int(n as i128)
    }
    fn from_u128(n: u128) -> Self {
        Rational::from_int(n as i128)
    }
    fn from_i128(n: i128) -> Self {
        Rational::from_int(n)
    }
    fn as_ratio(&self) -> Option<(i128, i128)> {
        // Always in lowest terms with a positive denominator, by construction.
        Some((self.num, self.den))
    }
    fn from_ratio(num: i128, den: i128) -> Option<Self> {
        (den != 0).then(|| Rational::new(num, den))
    }
}

impl Field for Rational {
    fn inv(&self) -> Self {
        assert!(self.num != 0, "inverse of zero Rational");
        Rational::new(self.den, self.num)
    }
}

impl Plethystic for Rational {
    /// ℚ has no variables to raise, so p_n fixes every scalar.
    fn frobenius(&self, _n: u32) -> Self {
        *self
    }
}

impl QAlgebra for Rational {
    fn div_u128(&self, n: u128) -> Self {
        assert!(n != 0, "division of Rational by zero");
        let n = n as i128;
        // Cancel against the numerator *before* multiplying the denominator.
        // The divisor here is z_μ, which reaches |μ|! — so `den * n` overflows
        // i128 far sooner than the reduced form does, and these two share
        // factors constantly in the formulas that call this.
        let g = gcd(self.num, n);
        Rational::new(self.num / g, self.den * (n / g))
    }
}

// ----------------------------------------------------------------------------
// GMP-backed coefficients (feature = "gmp").
//
// This is the payoff of making the coefficient ring a parameter: swapping in
// `rug::Integer` / `rug::Rational` gives GMP's Karatsuba/Toom-Cook/FFT
// multiplication for the large coefficients that plethysm and high-degree
// expansions produce, with no change to any basis or conversion code.
//
// Note `from_u128` is exact here — the whole reason that seam exists.
// ----------------------------------------------------------------------------

#[cfg(feature = "gmp")]
mod gmp_impls {
    use super::{Field, Plethystic, QAlgebra, Ring};
    use rug::{Integer, Rational as RugRational};

    impl Ring for Integer {
        fn zero() -> Self {
            Integer::new()
        }
        fn one() -> Self {
            Integer::from(1)
        }
        fn is_zero(&self) -> bool {
            *self == 0
        }
        fn add_assign(&mut self, other: &Self) {
            *self += other;
        }
        fn mul(&self, other: &Self) -> Self {
            Integer::from(self * other)
        }
        fn neg(&self) -> Self {
            Integer::from(-self)
        }
        fn from_i64(n: i64) -> Self {
            Integer::from(n)
        }
        fn from_u128(n: u128) -> Self {
            Integer::from(n) // exact, no truncation
        }
        fn from_i128(n: i128) -> Self {
            Integer::from(n) // exact, no truncation
        }
    }

    impl Ring for RugRational {
        fn zero() -> Self {
            RugRational::new()
        }
        fn one() -> Self {
            RugRational::from(1)
        }
        fn is_zero(&self) -> bool {
            *self == 0
        }
        fn add_assign(&mut self, other: &Self) {
            *self += other;
        }
        fn mul(&self, other: &Self) -> Self {
            RugRational::from(self * other)
        }
        fn neg(&self) -> Self {
            RugRational::from(-self)
        }
        fn from_i64(n: i64) -> Self {
            RugRational::from(n)
        }
        fn from_u128(n: u128) -> Self {
            RugRational::from(Integer::from(n)) // exact
        }
        fn from_i128(n: i128) -> Self {
            RugRational::from(Integer::from(n)) // exact
        }
    }

    impl Field for RugRational {
        fn inv(&self) -> Self {
            assert!(*self != 0, "inverse of zero Rational");
            RugRational::from(self.clone().recip())
        }
    }

    impl QAlgebra for RugRational {
        fn div_u128(&self, n: u128) -> Self {
            assert!(n != 0, "division of Rational by zero");
            RugRational::from(self / Integer::from(n))
        }
    }

    impl Plethystic for RugRational {
        /// ℚ, arbitrary precision: still no variables to raise.
        fn frobenius(&self, _n: u32) -> Self {
            self.clone()
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rational_arithmetic() {
        let mut a = Rational::new(1, 2);
        a.add_assign(&Rational::new(1, 3)); // 5/6
        assert_eq!(a, Rational::new(5, 6));
        assert_eq!(a.mul(&Rational::new(6, 5)), Rational::one());
        assert_eq!(Rational::new(2, 4), Rational::new(1, 2)); // normalization
        assert_eq!(Rational::new(3, -6), Rational::new(-1, 2)); // sign to numerator
    }

    #[test]
    fn rational_field_ops() {
        let a = Rational::new(3, 7);
        assert_eq!(a.mul(&a.inv()), Rational::one());
        assert_eq!(Rational::new(1, 2).div(&Rational::new(3, 4)), Rational::new(2, 3));
        let mut b = Rational::from_int(5);
        b.sub_assign(&Rational::new(1, 2)); // 9/2
        assert_eq!(b, Rational::new(9, 2));
    }
}
