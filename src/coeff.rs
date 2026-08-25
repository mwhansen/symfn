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
    /// The multiplicative inverse `self⁻¹`.
    ///
    /// # Panics
    ///
    /// Panics if `self` is zero.
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
    /// never rounds.
    ///
    /// # Panics
    ///
    /// Panics if `n` is zero. [`Rational`] additionally panics on a divisor
    /// past `i128::MAX`.
    fn div_u128(&self, n: u128) -> Self;
}

/// A ring with a gcd, on top of [`Ring::div_exact`]: the *integer* rings this
/// crate puts in the numerator and denominator of a fraction type.
///
/// [`AFrac`](crate::afrac::AFrac)'s tail is what needs it. A polynomial gcd
/// over a field is plain Euclid, but running it on rational coefficients
/// makes them grow
/// multiplicatively — the dense form of a Jack coefficient overflows `i128`
/// by degree 6 that way (`docs/record/jack.md`). The primitive-part algorithm
/// keeps the coefficients the size the integer ring already holds, and what it
/// needs beyond division is exactly this: the content of a polynomial, which
/// is the gcd of its coefficients.
///
/// A field may implement this with `one()`, which is correct — every nonzero
/// element is a unit — and reduces the algorithm to ordinary Euclid, growth
/// included. That is why the bound is a separate trait rather than a default
/// on [`Ring`]: a ring that cannot control the growth should not silently
/// look like one that can.
pub trait Integral: Ring {
    /// A greatest common divisor, never negative by
    /// [`is_negative`](Integral::is_negative)'s reckoning, and zero only when
    /// both arguments are.
    fn gcd(&self, other: &Self) -> Self;

    /// Whether this is negative, so that a sign can be normalized.
    fn is_negative(&self) -> bool;
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
    ///
    /// `n = 0` is outside the contract — `p_0` is not a raising of variables —
    /// and an implementation with variables to raise refuses it
    /// ([`QtPoly`](crate::qt::QtPoly), its `# Panics`). The constant rings
    /// accept it vacuously, since fixing every scalar asks nothing of `n`.
    fn frobenius(&self, n: u32) -> Self;
}

macro_rules! impl_ring_for_int {
    ($($t:ty => $div_exact:ident),*) => {$(
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
                $div_exact(*self, *other)
            }
        }
    )*};
}

impl_ring_for_int!(i64 => div_exact_i64, i128 => div_exact_i128);

macro_rules! impl_integral_for_int {
    ($($t:ty),*) => {$(
        impl Integral for $t {
            /// Euclid on the magnitudes. Never negative, so the sign
            /// normalization a primitive part does has a fixed target.
            ///
            /// # Panics
            ///
            /// Panics if either operand is `MIN` and the other is zero, where
            /// the gcd is the magnitude of `MIN` and does not fit. Every call
            /// site here is a polynomial's content, whose coefficients came
            /// through this ring's own checked injection, so this is a
            /// contract violation rather than a wall (`docs/policies/
            /// failure.md`, R8).
            #[inline]
            fn gcd(&self, other: &Self) -> Self {
                let g = gcd_u128(self.unsigned_abs() as u128, other.unsigned_abs() as u128);
                <$t>::try_from(g).unwrap_or_else(|_| {
                    panic!("the gcd {g} does not fit {}", stringify!($t))
                })
            }
            #[inline]
            fn is_negative(&self) -> bool {
                *self < 0
            }
        }
    )*};
}

impl_integral_for_int!(i64, i128);

#[inline]
fn div_exact_i64(a: i64, b: i64) -> Option<i64> {
    (b != 0 && a % b == 0).then(|| a / b)
}

/// The `i128` exact quotient, done in 64-bit arithmetic when both operands
/// fit — which in this library is nearly always. On aarch64 (and x86-64) a
/// 128-bit `%` or `/` is a call into `compiler_builtins`, some 9 ns each, and
/// a 64-bit one is an instruction; measured on Jack's `AFrac<i128>` reduction,
/// where this runs per coefficient, the narrowing alone is 1.12-1.15x on
/// `jack_table` (`docs/record/coefficient-arithmetic.md`).
///
/// `b == -1` is separated because `i64::MIN / -1` overflows `i64` while the
/// `i128` quotient is fine.
#[inline]
pub(crate) fn div_exact_i128(a: i128, b: i128) -> Option<i128> {
    if b == 0 {
        return None;
    }
    if let (Ok(a64), Ok(b64)) = (i64::try_from(a), i64::try_from(b)) {
        if b64 == -1 {
            return Some(-a);
        }
        return (a64 % b64 == 0).then(|| i128::from(a64 / b64));
    }
    (a % b == 0).then(|| a / b)
}

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
/// prints: `i128::MIN` has no negation inside the width, so `-MIN` and
/// `MIN.abs()` are overflow rather than arithmetic. Reachable only by
/// constructing one directly or by an arithmetic result landing exactly on
/// `MIN`, and a panic naming the requirement is the documented wall — the
/// escalating entry points run [`GuardedRat`](crate::guard::GuardedRat), which
/// reports instead and re-runs over `BigRational` (`docs/policies/failure.md`,
/// R5).
const NO_NEGATION: &str =
    "a Rational part of i128::MIN has no negation in i128; use the bignum ring";

/// The gcd every fixed-width rational in the crate reduces by, in 64-bit
/// arithmetic when both operands fit.
///
/// They nearly always do: instrumented over `s → p`, plethysm and the Jack and
/// (q,t)-Kostka routes, more than 99% of the operand pairs were below 2³² and
/// took two or three Euclid steps (`docs/record/coefficient-arithmetic.md`).
/// What the narrowing removes is not arithmetic but the call: a 128-bit `%`
/// is a `compiler_builtins` routine on every target this ships to, and a
/// 64-bit one is an instruction. Plain Euclid at 64 bits measured the same as
/// a Euclid-then-binary hybrid in situ, so it is the simple form; the wide
/// fallback is binary because it never divides.
pub(crate) fn gcd_u128(a: u128, b: u128) -> u128 {
    if (a | b) >> 64 == 0 {
        // Bounded by the check just above.
        #[allow(clippy::cast_possible_truncation)]
        let (mut a, mut b) = (a as u64, b as u64);
        while b != 0 {
            let t = a % b;
            a = b;
            b = t;
        }
        return u128::from(a);
    }
    gcd_u128_wide(a, b)
}

/// Binary (Stein) gcd on wide operands: shifts and subtractions only.
fn gcd_u128_wide(mut a: u128, mut b: u128) -> u128 {
    if a == 0 {
        return b;
    }
    if b == 0 {
        return a;
    }
    let shift = (a | b).trailing_zeros();
    a >>= a.trailing_zeros();
    loop {
        b >>= b.trailing_zeros();
        if a > b {
            core::mem::swap(&mut a, &mut b);
        }
        b -= a;
        if b == 0 {
            break;
        }
    }
    a << shift
}

/// `a / g` for a `g` that divides `a`, in 64-bit arithmetic when both fit —
/// the same call-versus-instruction saving as [`gcd_u128`], and skipped
/// outright for `g == 1`, which after a gcd is the common case.
#[inline]
pub(crate) fn quo(a: i128, g: i128) -> i128 {
    if g == 1 {
        return a;
    }
    if let (Ok(a64), Ok(g64)) = (i64::try_from(a), i64::try_from(g)) {
        // `g` is a gcd, so `g ≥ 1` at every call site, and the quotient's
        // magnitude is at most `|a|`: it fits where `a` fits.
        return i128::from(a64 / g64);
    }
    a / g
}

/// What the fixed-width rational arithmetic does with an integer operation
/// that leaves `i128` — the one thing [`Rational`] and
/// [`GuardedRat`](crate::guard::GuardedRat) differ on. `Rational` runs native
/// arithmetic, which panics under `overflow-checks` in every profile
/// (`docs/policies/failure.md`, R3); `GuardedRat` reports and continues on a
/// `0`, and its constructor refuses the pair afterwards. Everything else the
/// two share — the integer fast paths, Henrici's addition, cross-cancelled
/// multiplication, the `z_μ` division and normalization — is the `rat_*`
/// functions below, written once, so that a fast path can no longer reach one
/// ring and not the other (`docs/record/coefficient-arithmetic.md`).
pub(crate) trait Overflow {
    fn add(a: i128, b: i128) -> i128;
    fn mul(a: i128, b: i128) -> i128;
}

/// [`Rational`]'s policy: native arithmetic.
pub(crate) struct Panics;

impl Overflow for Panics {
    #[inline]
    fn add(a: i128, b: i128) -> i128 {
        a + b
    }
    #[inline]
    fn mul(a: i128, b: i128) -> i128 {
        a * b
    }
}

/// `gcd(|a|, |b|)` as an `i128`. At every call site below at least one operand
/// is a positive denominator or a positive divisor, so the gcd is at most that
/// operand and fits — even when the other is `i128::MIN`, whose magnitude
/// `unsigned_abs` still has.
#[inline]
#[allow(clippy::cast_possible_wrap)]
pub(crate) fn gcd_i128(a: i128, b: i128) -> i128 {
    gcd_u128(a.unsigned_abs(), b.unsigned_abs()) as i128
}

/// `num/den` in lowest terms with `den > 0`, for `den ≠ 0` and neither part
/// `i128::MIN` — the caller has applied its own wall to those.
#[inline]
pub(crate) fn rat_normalize(num: i128, den: i128) -> (i128, i128) {
    let (mut n, mut d) = (num, den);
    if d < 0 {
        n = -n;
        d = -d;
    }
    if n == 0 {
        return (0, 1);
    }
    if d == 1 {
        return (n, 1);
    }
    let g = gcd_i128(n, d);
    (quo(n, g), quo(d, g))
}

/// `a/b + c/d` in lowest terms, both inputs in lowest terms with positive
/// denominators.
///
/// Integers stay integers, and `gcd(n, 1) == 1` needs no Euclid to discover;
/// most rational arithmetic in this library never leaves ℤ. Then Henrici's
/// addition (Knuth, TAOCP 4.5.1) with `g = gcd(b, d)`: equal denominators
/// need one gcd, of the sum against `b`. If `g == 1`, `(ad + cb)/(bd)` is
/// already in lowest terms — a prime dividing `b` and `ad + cb` divides `ad`,
/// hence `a` — so no gcd at all. Otherwise `t = a(d/g) + c(b/g)` and
/// `g' = gcd(t, g)` give `(t/g') / ((b/g)(d/g'))`. Every gcd runs on operands
/// no wider than the denominators, and the intermediates are smaller than the
/// cross products by a factor of `g`, which is also what keeps them inside
/// `i128` longer.
#[inline]
pub(crate) fn rat_add<P: Overflow>(a: i128, b: i128, c: i128, d: i128) -> (i128, i128) {
    if b == 1 && d == 1 {
        return (P::add(a, c), 1);
    }
    if b == d {
        let t = P::add(a, c);
        if t == 0 {
            return (0, 1);
        }
        let g = gcd_i128(t, b);
        return (quo(t, g), quo(b, g));
    }
    let g = gcd_i128(b, d);
    if g == 1 {
        return (P::add(P::mul(a, d), P::mul(c, b)), P::mul(b, d));
    }
    let (bg, dg) = (quo(b, g), quo(d, g));
    let t = P::add(P::mul(a, dg), P::mul(c, bg));
    if t == 0 {
        return (0, 1);
    }
    let g2 = gcd_i128(t, g);
    (quo(t, g2), P::mul(bg, quo(d, g2)))
}

/// `(a/b)(c/d)` in lowest terms, both inputs in lowest terms with positive
/// denominators: the integer fast path, then cross-cancellation — `gcd(a, d)`
/// and `gcd(c, b)` first, so the products are formed from reduced factors and
/// are in lowest terms without a gcd of the products (Knuth, TAOCP 4.5.1).
#[inline]
pub(crate) fn rat_mul<P: Overflow>(a: i128, b: i128, c: i128, d: i128) -> (i128, i128) {
    if b == 1 && d == 1 {
        return (P::mul(a, c), 1);
    }
    if a == 0 || c == 0 {
        return (0, 1);
    }
    let g1 = gcd_i128(a, d);
    let g2 = gcd_i128(c, b);
    (
        P::mul(quo(a, g1), quo(c, g2)),
        P::mul(quo(b, g2), quo(d, g1)),
    )
}

/// `(a/b) / n` for a positive `n`, in lowest terms.
///
/// Cancels against the numerator *before* multiplying the denominator: the
/// divisor here is `z_μ`, which reaches `|μ|!`, so `b · n` overflows `i128`
/// far sooner than the reduced form does, and the two share factors constantly
/// in the formulas that call this. The result is already in lowest terms —
/// `gcd(a, b) = 1` and `gcd(a/g, n/g) = 1` give `gcd(a/g, b·(n/g)) = 1` — and
/// normalizing again would repeat the gcd to learn nothing; this is the one
/// division `s → p` does per term.
#[inline]
pub(crate) fn rat_div<P: Overflow>(a: i128, b: i128, n: i128) -> (i128, i128) {
    let g = gcd_i128(a, n);
    (quo(a, g), P::mul(b, quo(n, g)))
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
        let (num, den) = rat_normalize(num, den);
        Rational { num, den }
    }

    /// The integer `n` as `n/1`.
    ///
    /// Accepts `i128::MIN`, which [`Rational::new`] refuses. Negating that
    /// value then panics: `MIN` has no negation inside `i128`.
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
        let (num, den) = rat_add::<Panics>(self.num, self.den, other.num, other.den);
        *self = Rational { num, den };
    }
    fn mul(&self, other: &Self) -> Self {
        let (num, den) = rat_mul::<Panics>(self.num, self.den, other.num, other.den);
        Rational { num, den }
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

/// Every nonzero element of a field is a unit, so a gcd carries no
/// information — the primitive part of a polynomial over a field is the
/// polynomial, and the pseudo-division that [`Integral`] exists to control
/// degenerates to ordinary Euclid. Correct, and the growth this bound is meant
/// to bound comes back; the rings that matter for that are the integer ones.
impl Integral for Rational {
    fn gcd(&self, other: &Self) -> Self {
        if self.is_zero() && other.is_zero() {
            <Self as Ring>::zero()
        } else {
            <Self as Ring>::one()
        }
    }
    fn is_negative(&self) -> bool {
        self.num < 0
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
        let (num, den) = rat_div::<Panics>(self.num, self.den, n);
        Rational { num, den }
    }
}

// ----------------------------------------------------------------------------
// Arbitrary-precision coefficients (feature = "bignum").
//
// These cost nothing beyond the impls below, because the coefficient ring is a
// parameter: the same basis and conversion code runs over `BigInt` /
// `BigRational` with no change, which is what lets a fixed-width computation
// that overflows be re-run exactly rather than returning a wrapped answer.
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
    use super::{Field, Integral, Plethystic, QAlgebra, Ring};
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

    impl Integral for BigInt {
        /// Euclid on the magnitudes, written out rather than pulled from
        /// `num-integer`: this crate's dependency list is deliberately short
        /// (`Cargo.toml`), and the loop is three lines.
        fn gcd(&self, other: &Self) -> Self {
            let (mut a, mut b) = (self.abs(), other.abs());
            while !Zero::is_zero(&b) {
                let t = a % &b;
                a = b;
                b = t;
            }
            a
        }
        fn is_negative(&self) -> bool {
            Signed::is_negative(self)
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

    /// A field: see [`Integral`] for `Rational`.
    impl Integral for BigRational {
        fn gcd(&self, other: &Self) -> Self {
            if Zero::is_zero(self) && Zero::is_zero(other) {
                <Self as Ring>::zero()
            } else {
                <Self as Ring>::one()
            }
        }
        fn is_negative(&self) -> bool {
            Signed::is_negative(self)
        }
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

    /// Every branch of the addition and multiplication shortcuts — equal
    /// denominators, coprime ones, ones sharing a factor, integer fast paths,
    /// zeros, both signs — must give what "form the cross product and reduce
    /// it" gives. The reference here does exactly that, in `i128` with its own
    /// Euclid, so it shares no code with the shortcuts.
    #[test]
    fn shortcut_arithmetic_matches_reduce_after_the_fact() {
        fn plain_gcd(mut a: i128, mut b: i128) -> i128 {
            (a, b) = (a.abs(), b.abs());
            while b != 0 {
                (a, b) = (b, a % b);
            }
            a
        }
        fn reduce(n: i128, d: i128) -> (i128, i128) {
            if n == 0 {
                return (0, 1);
            }
            let g = plain_gcd(n, d);
            let (n, d) = (n / g, d / g);
            if d < 0 {
                (-n, -d)
            } else {
                (n, d)
            }
        }
        let nums = [-7i128, -6, -1, 0, 1, 2, 3, 5, 6, 12, 35, 1 << 40];
        let dens = [1i128, 2, 3, 4, 6, 7, 12, 30, 1 << 40, (1 << 40) + 1];
        for &a in &nums {
            for &b in &dens {
                for &c in &nums {
                    for &d in &dens {
                        let x = Rational::new(a, b);
                        let y = Rational::new(c, d);
                        let mut sum = x;
                        sum.add_assign(&y);
                        assert_eq!(
                            (sum.num, sum.den),
                            reduce(a * d + c * b, b * d),
                            "{a}/{b} + {c}/{d}"
                        );
                        let prod = x.mul(&y);
                        assert_eq!(
                            (prod.num, prod.den),
                            reduce(a * c, b * d),
                            "{a}/{b} * {c}/{d}"
                        );
                        for n in [1u128, 2, 6, 7, 1 << 40] {
                            let q = x.div_u128(n);
                            assert_eq!(
                                (q.num, q.den),
                                reduce(a, b * i128::try_from(n).unwrap()),
                                "{a}/{b} / {n}"
                            );
                        }
                    }
                }
            }
        }
    }

    /// The narrowed gcd and exact quotient must agree with the wide forms on
    /// operands that straddle 64 bits, and the wide gcd is binary — checked
    /// against Euclid on products of large primes, where a slip in the shift
    /// bookkeeping would show.
    #[test]
    fn narrowed_gcd_and_quotient_agree_with_the_wide_forms() {
        fn euclid(mut a: u128, mut b: u128) -> u128 {
            while b != 0 {
                (a, b) = (b, a % b);
            }
            a
        }
        let p = 1_000_000_007u128;
        let q = 998_244_353u128;
        let r = (1u128 << 61) - 1;
        let cases = [
            (0u128, 0u128),
            (0, 5),
            (12, 18),
            (u64::MAX as u128, 6),
            (u64::MAX as u128 + 1, 6),
            (p * q, q * r),
            (p * r * 8, q * r * 12),
            (r * r, r * p),
            (u128::MAX / 3, u128::MAX / 5),
        ];
        for &(a, b) in &cases {
            assert_eq!(gcd_u128(a, b), euclid(a, b), "gcd({a}, {b})");
            assert_eq!(gcd_u128(b, a), euclid(a, b), "gcd({b}, {a})");
        }
        for &(a, g) in &[
            (36i128, 6i128),
            (-36, 6),
            (i128::from(i64::MAX) * 4, 4),
            (i128::from(i64::MIN) * 3, 3),
            (i128::from(i64::MIN), 1),
            (i128::MAX, i128::MAX),
        ] {
            assert_eq!(quo(a, g), a / g, "{a} / {g}");
        }
        assert_eq!(
            div_exact_i128(i128::from(i64::MIN), -1),
            Some(-i128::from(i64::MIN))
        );
        assert_eq!(div_exact_i128(7, 0), None);
        assert_eq!(div_exact_i128(7, 2), None);
        assert_eq!(div_exact_i128(i128::MAX - 1, 2), Some((i128::MAX - 1) / 2));
        assert_eq!(
            div_exact_i128(-(1i128 << 100), 1 << 30),
            Some(-(1i128 << 70))
        );
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
