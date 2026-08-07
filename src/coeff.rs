//! The coefficient-ring abstraction: [`Ring`], then [`QAlgebra`] and
//! [`Plethystic`] for the two things a ring may additionally have to supply.
//!
//! Every symmetric-function type is generic over `C: Ring`, so the *same* basis
//! code works over machine integers today and over arbitrary-precision integers
//! (`BigInt`) or a `(q,t)`-polynomial ring tomorrow — the coefficient ring is a
//! parameter, never baked in.
//!
//! The layering is deliberate, and each step up is demanded by exactly one
//! thing:
//!
//! | bound | what needs it | why |
//! |---|---|---|
//! | [`Ring`] | s, h, e, m, f and their conversions; products; skewing | exact over ℤ |
//! | [`QAlgebra`] | `s → p`, internal product, plethysm | z_μ⁻¹ — division by an **integer** |
//! | [`Plethystic`] | plethysm | `p_n` acts on the coefficients too |
//!
//! [`Field`] appears in none of those rows, which is the point: the dividing
//! paths divide only by z_μ, so a ring containing ℚ suffices and need not
//! invert its own elements. `ℚ[t]` and `ℚ[q,t]` are the cases that matter —
//! neither is a field, both are fine — and they are precisely the coefficient
//! rings Hall–Littlewood and Macdonald need. `Field` is kept because it is a
//! real thing to name and [`Rational`] is one, but nothing in the library
//! requires it.
//!
//! Implementors: `i64`/`i128` are rings only; [`Rational`] and (under the
//! `bignum` feature) `BigRational` implement all of them, `BigInt` is a ring.
//! `tests/qalgebra.rs` carries a `ℚ[t]` implementing the two upper traits and
//! deliberately **not** `Field` — which is what keeps the dividing paths from
//! quietly drifting back to the stronger bound.
//!
//! A trait bound is only checked where it is instantiated, so adding a bound
//! can pass `cargo build` and still break a coefficient type nothing in the
//! crate constructs. `cargo test --features bignum` is what catches that; it is
//! not part of the default test run.

/// A commutative ring usable as a symmetric-function coefficient.
///
/// Intentionally small: enough for additive structure, multiplication, and the
/// subtraction that determinant / Newton-identity code needs.
///
/// # The seams
///
/// Several methods below are referred to throughout this crate and the record
/// as **seams** — places where behavior changes by swapping the coefficient
/// type, with no call site edited (the sense of Feathers, *Working Effectively
/// with Legacy Code*, ch. 4). What is swapped here is exactness rather than
/// testability: [`from_u128`](Ring::from_u128) refuses on `i64` and is exact on
/// `BigInt`, and the generic basis code above it is identical either way.
/// [`from_u128`](Ring::from_u128), [`from_i128`](Ring::from_i128),
/// [`div_exact`](Ring::div_exact) and the
/// [`as_ratio`](Ring::as_ratio)/[`from_ratio`](Ring::from_ratio) pair are the
/// seams that matter, because structure constants and division are where a
/// fixed-width ring runs out of range — see `docs/policies/failure.md` (R8).
pub trait Ring: Clone + PartialEq + core::fmt::Debug {
    /// The additive identity.
    fn zero() -> Self;
    /// The multiplicative identity.
    fn one() -> Self;
    /// Whether this is the additive identity.
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
    /// routing them through [`Ring::from_i64`] would silently truncate the
    /// large ones. A fixed-width coefficient type necessarily runs out of range
    /// here and says so; a bignum type (`BigInt` under the `bignum` feature)
    /// overrides this to be exact — which is what the seam is for.
    ///
    /// An implementor that narrows here must **not** truncate: check and panic
    /// naming the constant and the ring, or report and escalate the way
    /// [`Guarded`](crate::guard::Guarded) does
    /// (`docs/policies/failure.md`, R8). The default holds itself to that.
    ///
    /// # Panics
    ///
    /// Panics if `n` is past `i64::MAX` and the implementor has not overridden
    /// this.
    fn from_u128(n: u128) -> Self {
        Self::from_i64(i64::try_from(n).unwrap_or_else(|_| {
            panic!("the structure constant {n} does not fit i64; this ring must override from_u128")
        }))
    }

    /// This value as an exact ratio of `i128`s, if it is one.
    ///
    /// The seam for putting a batch of coefficients over a **common
    /// denominator**, so that a sum of (coefficient × integer) terms can be
    /// accumulated in `i128` and converted back once at the end instead of
    /// doing rational arithmetic per term. `p → s` does exactly that shape of
    /// sum, where the gcds it implies are otherwise most of plethysm's runtime
    /// (`docs/record/plethysm.md`).
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
    /// negative, and held to the same rule.
    ///
    /// # Panics
    ///
    /// Panics if `n` is outside `i64` and the implementor has not overridden
    /// this.
    fn from_i128(n: i128) -> Self {
        Self::from_i64(i64::try_from(n).unwrap_or_else(|_| {
            panic!("the character value {n} does not fit i64; this ring must override from_i128")
        }))
    }

    /// `self -= other`, provided via [`Ring::neg`].
    fn sub_assign(&mut self, other: &Self) {
        let neg = other.neg();
        self.add_assign(&neg);
    }

    /// `self / other` when the quotient exists **in this ring**; `None`
    /// otherwise, including when `other` is zero.
    ///
    /// The seam for exact polynomial division — see
    /// [`QtPoly::divide_exact`](crate::qt::QtPoly::divide_exact), which needs
    /// to divide coefficients as it eliminates leading terms. A general ring
    /// has no division at all, and this asks for much less than one: not an
    /// inverse, only the answer in the cases where one exists. ℤ has that and
    /// is not a field, which is why this is separate from
    /// [`Field::inv`].
    ///
    /// Default `None`, the same shape as [`Ring::as_ratio`]: a ring that cannot
    /// answer declines, and every caller must already handle `None` because
    /// *not divisible* is an ordinary outcome. Declining is therefore
    /// conservative rather than wrong — a division that could have succeeded is
    /// reported as a failure, never the reverse.
    fn div_exact(&self, _other: &Self) -> Option<Self> {
        None
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
/// **This, and not [`Field`], is what dividing in this library actually
/// needs**, and that difference is why the trait exists. Every
/// division the crate performs is by `z_μ`, a positive integer: `s → p` carries
/// z_μ⁻¹, and the internal product and plethysm inherit it by routing through
/// the power-sum basis. Nothing ever divides by a general ring element.
///
/// Bounding those paths on [`Field`] would demand far more than the
/// mathematics does, and would exclude exactly the rings this library most
/// wants to serve. `ℚ[t]` is not a field, nor is `ℚ[q,t]` — but z_μ⁻¹ lives in
/// both, so `s → p` over them is well defined. The same applies to any
/// ℚ-algebra a caller brings across the Sage boundary, and to the
/// (q,t)-coefficient rings the Macdonald and Hall–Littlewood work needs.
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
///   p_n[t · p_1] = t^n · p_n            (over `ℚ[t]`)
/// ```
///
/// which is what Sage computes, and what `exclude=[t]` there turns off.
///
/// So a ring is usable for plethysm only once it says which map that is, and
/// this trait is that statement. For ℚ there is nothing to raise and it is the
/// identity; for `ℚ[t]` it is t ↦ t^n; for a ring whose variables should be
/// held *constant* it is the identity again, which is a real convention choice
/// the implementor makes rather than a default this crate can pick.
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
            /// Widening for `i128`, the identity for `i64`: exact either way.
            #[inline] fn from_i64(n: i64) -> Self { n as $t }
            /// # Panics
            ///
            /// Panics if the constant does not fit. This is the seam large
            /// values enter through — LR coefficients, Kostka numbers, `z_λ`,
            /// characters — so truncating here would put a wrong structure
            /// constant into an otherwise exact computation. Injection is never
            /// a hot loop, so the check costs nothing that matters
            /// (`docs/policies/failure.md`, R8).
            #[inline] fn from_u128(n: u128) -> Self {
                <$t>::try_from(n).unwrap_or_else(|_| panic!(
                    "the structure constant {n} does not fit {}; use the bignum ring",
                    stringify!($t)
                ))
            }
            /// # Panics
            ///
            /// Panics if the constant does not fit — see [`Ring::from_u128`].
            #[inline] fn from_i128(n: i128) -> Self {
                <$t>::try_from(n).unwrap_or_else(|_| panic!(
                    "the structure constant {n} does not fit {}; use the bignum ring",
                    stringify!($t)
                ))
            }
            /// Exact in ℤ: divides only when the remainder is zero.
            #[inline] fn div_exact(&self, other: &Self) -> Option<Self> {
                (*other != 0 && *self % *other == 0).then(|| *self / *other)
            }
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

/// The message every `Rational` site that needs a magnitude or a sign flip
/// prints: `i128::MIN` has no negation inside the width, so `-MIN`, `MIN.abs()`
/// and `gcd(MIN, ·)` are overflow rather than arithmetic. Reachable only by
/// constructing one directly or by an arithmetic result landing exactly on
/// `MIN`, and a panic naming the requirement is the documented wall — the
/// escalating entry points run [`GuardedRat`](crate::guard::GuardedRat), which
/// reports instead and re-runs over `BigRational` (`docs/policies/failure.md`,
/// R5).
const NO_NEGATION: &str =
    "a Rational part of i128::MIN has no negation in i128; use the bignum ring";

fn gcd(mut a: i128, mut b: i128) -> i128 {
    assert!(a != i128::MIN && b != i128::MIN, "{NO_NEGATION}");
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
    /// Construct `num/den` in lowest terms.
    ///
    /// # Panics
    ///
    /// Panics if `den == 0`, or if either part is `i128::MIN`, which has no
    /// negation inside the width and so cannot be normalized.
    pub fn new(num: i128, den: i128) -> Self {
        assert!(den != 0, "Rational with zero denominator");
        assert!(num != i128::MIN && den != i128::MIN, "{NO_NEGATION}");
        let mut n = num;
        let mut d = den;
        if d < 0 {
            n = -n;
            d = -d;
        }
        if n == 0 {
            return Rational { num: 0, den: 1 };
        }
        if d == 1 {
            return Rational { num: n, den: 1 };
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

    /// The numerator, in lowest terms. The sign of the rational lives here.
    pub fn numer(&self) -> i128 {
        self.num
    }
    /// The denominator, in lowest terms and always positive; 1 for an integer,
    /// and 1 for zero.
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
        // Integers stay integers, and `gcd(n, 1) == 1` needs no Euclid to
        // discover. Most `Rational` arithmetic in this library never leaves ℤ —
        // a Macdonald `J` over ℚ(q,t) is integral throughout, and the fractions
        // only appear at `s → p`, where `z_ν⁻¹` enters. Without this guard every
        // one of those integer additions paid a 128-bit gcd: `u128_div_rem` was
        // 29% of a (q,t)-Kostka profile and `Rational::add_assign` another 21%.
        if self.den == 1 && other.den == 1 {
            self.num += other.num;
            return;
        }
        // a/b + c/d = (ad + cb) / bd, then normalize.
        let num = self.num * other.den + other.num * self.den;
        let den = self.den * other.den;
        *self = Rational::new(num, den);
    }
    fn mul(&self, other: &Self) -> Self {
        if self.den == 1 && other.den == 1 {
            return Rational {
                num: self.num * other.num,
                den: 1,
            };
        }
        Rational::new(self.num * other.num, self.den * other.den)
    }
    /// # Panics
    ///
    /// Panics if the numerator is `i128::MIN` — see [`Rational::new`]. The
    /// arithmetic fast paths below construct the fields directly, so this
    /// cannot be ruled out by construction and is checked here instead.
    fn neg(&self) -> Self {
        assert!(self.num != i128::MIN, "{NO_NEGATION}");
        Rational {
            num: -self.num,
            den: self.den,
        }
    }
    fn from_i64(n: i64) -> Self {
        Rational::from_int(n as i128)
    }
    /// # Panics
    ///
    /// Panics if the constant is past `i128::MAX` — `z_λ` reaches `|λ|!` and
    /// passes `i128` at λ ⊢ 34, so this is a wall a caller can reach, and
    /// truncating it would put a wrong `z_λ` under an otherwise exact division
    /// (`docs/policies/failure.md`, R8).
    fn from_u128(n: u128) -> Self {
        Rational::from_int(i128::try_from(n).unwrap_or_else(|_| {
            panic!("the structure constant {n} does not fit i128; use the bignum ring")
        }))
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
    /// A field: every nonzero divisor works.
    fn div_exact(&self, other: &Self) -> Option<Self> {
        (!other.is_zero()).then(|| self.mul(&other.inv()))
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
    /// # Panics
    ///
    /// Panics if `n == 0`, or if `n` is past `i128::MAX` — the divisor here is
    /// `z_μ`, which reaches `|μ|!`, so this is a wall a caller can reach rather
    /// than a contract violation, and it is the one
    /// [`GuardedRat`](crate::guard::GuardedRat) reports instead of panicking.
    fn div_u128(&self, n: u128) -> Self {
        assert!(n != 0, "division of Rational by zero");
        let n = i128::try_from(n)
            .unwrap_or_else(|_| panic!("a divisor of {n} is past i128::MAX; use the bignum ring"));
        // Cancel against the numerator *before* multiplying the denominator.
        // The divisor here is z_μ, which reaches |μ|! — so `den * n` overflows
        // i128 far sooner than the reduced form does, and these two share
        // factors constantly in the formulas that call this.
        let g = gcd(self.num, n);
        Rational::new(self.num / g, self.den * (n / g))
    }
}

// ----------------------------------------------------------------------------
// Arbitrary-precision coefficients (feature = "bignum").
//
// This is the payoff of making the coefficient ring a parameter: the same basis
// and conversion code runs over `BigInt` / `BigRational` with no change, which
// is what lets a fixed-width computation that overflows be re-run exactly
// rather than returning a wrapped answer.
//
// `num-bigint` and not GMP, on measurement rather than taste: coefficients past
// `i128` are 2-5 limbs (`examples/coeff_sizes.rs`), where every library runs
// schoolbook and GMP's asymptotically-fast paths never engage, and coefficient
// arithmetic is ~4% of runtime. Exactness is the requirement; speed is not. The
// license settles the rest -- `rug` is LGPL-3.0+, this crate is MIT OR
// Apache-2.0, and the wheel has to be distributable under the latter.
//
// `from_u128` and `from_i128` are exact here -- the reason those
// seams exist.
// ----------------------------------------------------------------------------

#[cfg(feature = "bignum")]
mod bignum_impls {
    use super::{Field, Plethystic, QAlgebra, Ring};
    use num_bigint::BigInt;
    use num_rational::BigRational;
    use num_traits::{One, Signed, Zero};

    impl Ring for BigInt {
        fn zero() -> Self {
            <BigInt as Zero>::zero()
        }
        fn one() -> Self {
            <BigInt as One>::one()
        }
        fn is_zero(&self) -> bool {
            Zero::is_zero(self)
        }
        fn add_assign(&mut self, other: &Self) {
            *self += other;
        }
        fn mul(&self, other: &Self) -> Self {
            self * other
        }
        fn neg(&self) -> Self {
            -self
        }
        fn from_i64(n: i64) -> Self {
            BigInt::from(n)
        }
        fn from_u128(n: u128) -> Self {
            BigInt::from(n) // exact
        }
        fn from_i128(n: i128) -> Self {
            BigInt::from(n) // exact
        }
        /// Exact in ℤ: divides only when the remainder is zero.
        fn div_exact(&self, other: &Self) -> Option<Self> {
            (!Zero::is_zero(other) && Zero::is_zero(&(self % other))).then(|| self / other)
        }
    }

    impl Ring for BigRational {
        fn zero() -> Self {
            <BigRational as Zero>::zero()
        }
        fn one() -> Self {
            <BigRational as One>::one()
        }
        fn is_zero(&self) -> bool {
            Zero::is_zero(self)
        }
        fn add_assign(&mut self, other: &Self) {
            *self += other;
        }
        fn mul(&self, other: &Self) -> Self {
            self * other
        }
        fn neg(&self) -> Self {
            -self
        }
        fn from_i64(n: i64) -> Self {
            BigRational::from(BigInt::from(n))
        }
        fn from_u128(n: u128) -> Self {
            BigRational::from(BigInt::from(n)) // exact
        }
        fn from_i128(n: i128) -> Self {
            BigRational::from(BigInt::from(n)) // exact
        }
        fn div_exact(&self, other: &Self) -> Option<Self> {
            (!Zero::is_zero(other)).then(|| self / other)
        }

        // Deliberately *not* implemented: `as_ratio` / `from_ratio` would have
        // to answer in `i128`, and a value that needed this ring in the first
        // place is exactly one that does not fit. Declining keeps
        // `convert::integral_sweep` on its generic path, which is always
        // correct, instead of narrowing and losing digits.
    }

    impl Field for BigRational {
        fn inv(&self) -> Self {
            assert!(!Zero::is_zero(self), "inverse of zero BigRational");
            self.recip()
        }
    }

    impl QAlgebra for BigRational {
        fn div_u128(&self, n: u128) -> Self {
            assert!(n != 0, "division of BigRational by zero");
            self / BigRational::from(BigInt::from(n))
        }
    }

    impl Plethystic for BigRational {
        /// Arbitrary-precision ℚ: still no variables to raise.
        fn frobenius(&self, _n: u32) -> Self {
            self.clone()
        }
    }

    /// `BigInt` is a ring but not a ℚ-algebra, so it cannot carry `s → p`.
    /// This is the seam the escalation path uses to decide which of the two
    /// bignum types a given operation needs.
    pub fn _assert_bounds() {
        fn ring<T: Ring>() {}
        fn qalg<T: QAlgebra>() {}
        ring::<BigInt>();
        ring::<BigRational>();
        qalg::<BigRational>();
    }

    /// Whether `v` is representable in `i128`, used to decide if a bignum
    /// answer can be handed back through a fixed-width channel.
    pub fn fits_i128(v: &BigInt) -> bool {
        v.abs() <= BigInt::from(i128::MAX)
    }
}

#[cfg(feature = "bignum")]
pub use bignum_impls::fits_i128;

#[cfg(test)]
mod tests {
    use super::*;

    /// Structure constants — LR coefficients, Kostka numbers, `z_λ`,
    /// characters — enter through these two seams, so a truncation here is a
    /// wrong constant inside an otherwise exact computation, and the loudest
    /// available failure is the cheapest one: injection is never a hot loop.
    #[test]
    fn fixed_width_injections_refuse_rather_than_truncate() {
        // Inside the width, both are exact.
        assert_eq!(
            <i64 as Ring>::from_u128(u64::MAX as u128 / 2),
            9223372036854775807
        );
        assert_eq!(
            <i128 as Ring>::from_u128(u128::MAX / 2),
            170141183460469231731687303715884105727
        );
        assert_eq!(
            <Rational as Ring>::from_u128(1 << 100),
            Rational::from_int(1 << 100)
        );

        for (name, f) in [
            (
                "i64 from_u128",
                (|| {
                    <i64 as Ring>::from_u128(u64::MAX as u128 + 1);
                }) as fn(),
            ),
            ("i64 from_i128", || {
                <i64 as Ring>::from_i128(i64::MAX as i128 + 1);
            }),
            ("i128 from_u128", || {
                <i128 as Ring>::from_u128(i128::MAX as u128 + 1);
            }),
            ("Rational from_u128", || {
                <Rational as Ring>::from_u128(i128::MAX as u128 + 1);
            }),
        ] {
            let out = std::panic::catch_unwind(f);
            assert!(out.is_err(), "{name} must refuse a constant past its width");
        }
    }

    /// `Rational` is the fixed-width field, so its `i128::MIN` corner is a
    /// documented wall rather than a report: the escalating entry points run
    /// `GuardedRat`, which reports it and re-runs over `BigRational`.
    #[test]
    #[should_panic(expected = "no negation in i128")]
    fn rational_refuses_the_width_minimum_as_a_numerator() {
        let _ = Rational::new(i128::MIN, 3);
    }

    #[test]
    #[should_panic(expected = "no negation in i128")]
    fn negating_a_width_minimum_rational_refuses() {
        // Past the constructor: the arithmetic fast paths build the fields
        // directly, so `neg` cannot assume normalization ruled this out.
        let min = Rational::from_int(i128::MIN);
        let _ = min.neg();
    }

    /// `z_μ` reaches `|μ|!` and is passed as `u128`, so a divisor past
    /// `i128::MAX` is a wall a caller can reach — it used to narrow with `as`,
    /// which turns a large positive divisor into a negative one.
    #[test]
    #[should_panic(expected = "past i128::MAX")]
    fn dividing_by_a_divisor_past_the_width_refuses() {
        let _ = Rational::one().div_u128(u128::MAX);
    }

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
        assert_eq!(
            Rational::new(1, 2).div(&Rational::new(3, 4)),
            Rational::new(2, 3)
        );
        let mut b = Rational::from_int(5);
        b.sub_assign(&Rational::new(1, 2)); // 9/2
        assert_eq!(b, Rational::new(9, 2));
    }
}
