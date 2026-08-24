//! ℚ(α) with no restriction on the denominator.
//!
//! [`AFrac`] holds ℚ(α) in the shape the Jack engines produce it: a numerator
//! over an integer times a product of *linear* forms `uα + v`. Every
//! denominator those engines build is a product of hooks, so the class is
//! closed and the factored form is what keeps the arithmetic cheap — matching
//! factors cancel before anything is expanded.
//!
//! Plethysm leaves that class. The nth plethystic Frobenius raises the
//! variables of the coefficient ring, so over ℚ(α) it is α ↦ α^n, and a
//! denominator `α + 1` becomes `α² + 1` at `n = 2` — irreducible over ℚ, and
//! not a product of linear forms. `docs/record/jack.md` records the Sage
//! values showing that is the mathematics and not the encoding.
//!
//! ## Why the factored form cannot simply be widened
//!
//! The obvious repair is to let an atom be `u·α^k + v` rather than `uα + v`,
//! since that family *is* closed under the Frobenius. It does not work,
//! because such an atom can factor: `α³ + 8 = (α + 2)(α² − 2α + 4)`. A
//! denominator holding `α³ + 8` and one holding `α + 2` times `α² − 2α + 4`
//! would then be two spellings of one element, and the representation would
//! stop being canonical. The same objection defeats every partly-factored
//! form: only a full gcd reduction decides equality by structure.
//!
//! So this type keeps numerator and denominator **dense, integral, and
//! coprime**, with the denominator primitive and positively led. That is a
//! normal form — for a nonzero element there is exactly one such pair — so
//! [`PartialEq`] is structural, unlike [`AFrac`]'s, which cross-multiplies
//! because its integer content is not canonical.
//!
//! ## Integral, not monic, and that is a measured choice
//!
//! The textbook normal form makes the denominator monic, which needs only a
//! field. It was written that way first and it does not work here: dividing by
//! the leading coefficient puts fractions in both parts, they multiply up, and
//! the dense form of a Jack coefficient **overflows `i128` by degree 6** —
//! where [`AFrac`] reaches degree 7 comfortably. The integral form keeps every
//! coefficient the size the integer ring already holds, at the cost of a
//! primitive-part gcd rather than a monic one (\[GCL\] ch. 2), and of the
//! [`Integral`] bound that gcd needs.
//!
//! The rings this is instantiated over are therefore the *integer* pair
//! [`Guarded`](crate::guard::Guarded) and `BigInt`, the same pair
//! [`AFrac`] crosses the Python boundary over.

use crate::afrac::{AFrac, Alpha, FromAFrac, Linears};
use crate::coeff::{Field, Integral, Plethystic, Ring};

/// An element of ℚ(α), numerator and denominator dense in α.
///
/// The invariants, together the normal form the module doc describes:
///
/// - `num` has no trailing zeros; empty is 0.
/// - `den` is never empty, is primitive, and its leading coefficient is
///   positive. When `num` is empty, `den` is `[1]`.
/// - `num` and `den` are coprime in `ℚ[α]`.
///
/// ```
/// use symfn::{ARat, Plethystic, Ring};
///
/// // 1/(α + 1)
/// let a = ARat::<i128>::from_parts(vec![1], vec![1, 1]);
/// // and the same element written with a common factor left in
/// let same = ARat::<i128>::from_parts(vec![2, 2], vec![2, 4, 2]);
/// assert_eq!(a, same);
///
/// assert_eq!(a.frobenius(2), ARat::<i128>::from_parts(vec![1], vec![1, 0, 1]));
/// ```
///
/// `p_2` sends `1/(α + 1)` to `1/(α² + 1)`, **not** to `1/(α + 1)²`: it raises
/// the variable rather than the whole form. `α² + 1` is irreducible over ℚ,
/// which is why this cannot be an [`AFrac`].
#[derive(Clone, Debug)]
pub struct ARat<C: Integral> {
    num: Vec<C>,
    den: Vec<C>,
}

fn trim<C: Ring>(v: &mut Vec<C>) {
    while v.last().is_some_and(C::is_zero) {
        v.pop();
    }
}

/// `a + b`, dense.
fn padd<C: Integral>(a: &[C], b: &[C]) -> Vec<C> {
    let mut out = Vec::with_capacity(a.len().max(b.len()));
    for k in 0..a.len().max(b.len()) {
        let mut c = a.get(k).cloned().unwrap_or_else(C::zero);
        if let Some(y) = b.get(k) {
            c.add_assign(y);
        }
        out.push(c);
    }
    trim(&mut out);
    out
}

/// `a · b`, dense. Schoolbook: the α-degrees here stay small next to the
/// partition counts around them.
fn pmul<C: Integral>(a: &[C], b: &[C]) -> Vec<C> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let mut out = vec![C::zero(); a.len() + b.len() - 1];
    for (i, x) in a.iter().enumerate() {
        if x.is_zero() {
            continue;
        }
        for (j, y) in b.iter().enumerate() {
            let t = x.mul(y);
            out[i + j].add_assign(&t);
        }
    }
    trim(&mut out);
    out
}

/// `a · c` for a scalar `c`.
fn pscale<C: Integral>(a: &[C], c: &C) -> Vec<C> {
    if c.is_zero() {
        return Vec::new();
    }
    a.iter().map(|x| x.mul(c)).collect()
}

/// `a / c` for a scalar that divides every coefficient.
///
/// # Panics
///
/// Panics if it does not, which is a broken invariant rather than a wall: this
/// is only ever called with a content or a gcd of the very coefficients it
/// divides.
fn pdiv_exact<C: Integral>(a: &[C], c: &C) -> Vec<C> {
    a.iter()
        .map(|x| {
            x.div_exact(c)
                .expect("a content divides the coefficients it was taken from")
        })
        .collect()
}

/// The gcd of the coefficients, non-negative, and zero only for the zero
/// polynomial.
fn content<C: Integral>(a: &[C]) -> C {
    let mut g = C::zero();
    for c in a {
        g = g.gcd(c);
    }
    g
}

/// `a` divided by its content, with a positive leading coefficient.
fn primitive<C: Integral>(a: &[C]) -> Vec<C> {
    if a.is_empty() {
        return Vec::new();
    }
    let g = content(a);
    let mut out = pdiv_exact(a, &g);
    if out[out.len() - 1].is_negative() {
        out = out.iter().map(C::neg).collect();
    }
    out
}

/// The pseudo-remainder of `a` by `b`, **up to content**.
///
/// Pseudo-division is how a division algorithm runs over a ring that does not
/// invert its leading coefficients: scaling by `lc(b)` before each subtraction
/// is what makes the step exact (\[GCL\] ch. 2). The textbook form leaves the
/// whole `lc(b)^{d+1}` in, and that is the inflation this takes back out —
/// the primitive part is taken after *every* step rather than once at the end.
///
/// ⚠️ **The result is therefore not the pseudo-remainder**, only an integer
/// multiple of it, which is all a gcd needs: content is divided out anyway.
/// Left in, `lc(b)^{d+1}` overflows `i128` on a degree-7 Jack coefficient,
/// where the answer itself does not.
///
/// # Panics
///
/// Panics if `b` is zero.
fn prem<C: Integral>(a: &[C], b: &[C]) -> Vec<C> {
    assert!(!b.is_empty(), "pseudo-division by the zero polynomial");
    if a.len() < b.len() {
        return a.to_vec();
    }
    let lc = b[b.len() - 1].clone();
    let mut r = a.to_vec();
    while r.len() >= b.len() && !r.is_empty() {
        let shift = r.len() - b.len();
        let factor = r[r.len() - 1].clone();
        r = pscale(&r, &lc);
        for (i, y) in b.iter().enumerate() {
            let t = y.mul(&factor).neg();
            r[shift + i].add_assign(&t);
        }
        trim(&mut r);
        if !r.is_empty() {
            r = primitive(&r);
        }
    }
    r
}

/// The gcd of `a` and `b` in `ℤ[α]`, primitive and positively led.
///
/// The primitive-part algorithm: divide out the content, run pseudo-division,
/// and take the primitive part of every remainder so the pseudo-division's
/// inflation does not accumulate. `gcd(0, 0)` is `[1]`, which is what the
/// callers here want — a zero numerator carries no denominator at all.
fn pgcd<C: Integral>(a: &[C], b: &[C]) -> Vec<C> {
    if a.is_empty() && b.is_empty() {
        return vec![C::one()];
    }
    if a.is_empty() {
        return primitive(b);
    }
    if b.is_empty() {
        return primitive(a);
    }
    let (mut x, mut y) = (primitive(a), primitive(b));
    if x.len() < y.len() {
        core::mem::swap(&mut x, &mut y);
    }
    while !y.is_empty() {
        let r = prem(&x, &y);
        x = y;
        y = if r.is_empty() { r } else { primitive(&r) };
    }
    primitive(&x)
}

impl<C: Integral> ARat<C> {
    /// A polynomial in α, given densely by its coefficients.
    pub fn from_coeffs(mut num: Vec<C>) -> Self {
        trim(&mut num);
        ARat {
            num,
            den: vec![C::one()],
        }
    }

    /// `num / den`, put in the normal form.
    ///
    /// # Panics
    ///
    /// Panics if `den` is zero.
    pub fn from_parts(mut num: Vec<C>, mut den: Vec<C>) -> Self {
        trim(&mut num);
        trim(&mut den);
        assert!(!den.is_empty(), "a zero denominator in Q(alpha)");
        let mut out = ARat { num, den };
        out.normalize();
        out
    }

    /// Restore the normal form: divide out the polynomial gcd, then the
    /// integer content the two parts share, then fix the denominator's sign.
    fn normalize(&mut self) {
        if self.num.is_empty() {
            self.den = vec![C::one()];
            return;
        }
        let g = pgcd(&self.num, &self.den);
        if g.len() > 1 {
            self.num = divide_exact(&self.num, &g);
            self.den = divide_exact(&self.den, &g);
        }
        let shared = content(&self.num).gcd(&content(&self.den));
        if !shared.is_zero() && shared != C::one() {
            self.num = pdiv_exact(&self.num, &shared);
            self.den = pdiv_exact(&self.den, &shared);
        }
        if self.den[self.den.len() - 1].is_negative() {
            self.num = self.num.iter().map(C::neg).collect();
            self.den = self.den.iter().map(C::neg).collect();
        }
    }

    /// The numerator and denominator, both dense in α and in the normal form.
    pub fn parts(&self) -> (&[C], &[C]) {
        (&self.num, &self.den)
    }

    /// The α-degree of the numerator, `None` for zero.
    pub fn degree(&self) -> Option<usize> {
        self.num.len().checked_sub(1)
    }

    /// This element as a polynomial in α, or `None` if it is not one.
    pub fn into_poly(self) -> Option<Vec<C>> {
        (self.den.len() == 1 && self.den[0] == C::one()).then_some(self.num)
    }

    /// The polynomial `uα + v`.
    pub fn linear(u: u32, v: u32) -> Self {
        Self::from_coeffs(vec![C::from_u128(v as u128), C::from_u128(u as u128)])
    }

    /// `1 / (uα + v)`.
    ///
    /// # Panics
    ///
    /// Panics on `0α + 0`, which would be `1/0`.
    pub fn inv_linear(u: u32, v: u32) -> Self {
        assert!(u > 0 || v > 0, "0*alpha + 0 is zero");
        Self::from_parts(
            vec![C::one()],
            vec![C::from_u128(v as u128), C::from_u128(u as u128)],
        )
    }

    /// Multiply by `∏ (uα + v)^m`, negative `m` meaning a denominator factor.
    ///
    /// The `Q` and `J` normalizers are products of hooks and arrive in exactly
    /// this shape, so this is the same seam
    /// [`AFrac::mul_factors`](crate::afrac::AFrac::mul_factors) offers. Here it
    /// is only a convenience over [`Ring::mul`], since nothing stays factored:
    /// the hooks are accumulated into one numerator and one denominator before
    /// a single reduction, rather than reduced against the value one at a time.
    ///
    /// # Panics
    ///
    /// Panics if a key with a nonzero multiplicity is `(0, 0)`, the zero form.
    pub fn mul_factors(&self, factors: &Linears) -> Self {
        let mut up = vec![C::one()];
        let mut down = vec![C::one()];
        for (&(u, v), &m) in factors {
            if m == 0 {
                continue;
            }
            assert!(u > 0 || v > 0, "0*alpha + 0 is zero");
            let linear = vec![C::from_u128(v as u128), C::from_u128(u as u128)];
            let side = if m > 0 { &mut up } else { &mut down };
            for _ in 0..m.unsigned_abs() {
                *side = pmul(side, &linear);
            }
        }
        Self::from_parts(pmul(&self.num, &up), pmul(&self.den, &down))
    }

    /// The same element, read out of the factored form the Jack engines build.
    ///
    /// The one direction that always works: every [`AFrac`] is an element of
    /// ℚ(α), and expanding its atoms loses only the factorization. The reverse
    /// does not exist, which is the whole reason this type is here.
    pub fn from_afrac(f: &AFrac<C>) -> Self {
        let (num, atoms, scale) = f.parts();
        let mut den = vec![C::from_u128(scale)];
        for (&(u, v), &m) in atoms {
            let linear = vec![C::from_u128(v as u128), C::from_u128(u as u128)];
            for _ in 0..m {
                den = pmul(&den, &linear);
            }
        }
        Self::from_parts(num.to_vec(), den)
    }
}

/// `a / b`, exact by assumption — `b` is a gcd of `a` and something else.
///
/// Written as a pseudo-division with the inflation divided back out, so no
/// coefficient inversion is needed.
///
/// # Panics
///
/// Panics if the division is not exact, which is a broken invariant.
fn divide_exact<C: Integral>(a: &[C], b: &[C]) -> Vec<C> {
    assert!(!b.is_empty(), "exact division by the zero polynomial");
    if a.is_empty() {
        return Vec::new();
    }
    let lc = b[b.len() - 1].clone();
    let mut r = a.to_vec();
    let mut q = vec![C::zero(); a.len() - b.len() + 1];
    while r.len() >= b.len() && !r.is_empty() {
        let shift = r.len() - b.len();
        let factor = r[r.len() - 1]
            .div_exact(&lc)
            .expect("a gcd's leading coefficient divides the quotient's");
        q[shift] = factor.clone();
        for (i, y) in b.iter().enumerate() {
            let t = y.mul(&factor).neg();
            r[shift + i].add_assign(&t);
        }
        trim(&mut r);
    }
    assert!(r.is_empty(), "an exact division left a remainder");
    trim(&mut q);
    q
}

impl<C: Integral + Field> ARat<C> {
    /// The value at `alpha`, or `None` where the denominator vanishes.
    ///
    /// The `None` is real rather than defensive: a coefficient of a Jack
    /// polynomial genuinely blows up at the α its hooks vanish at, even where
    /// the whole family stays finite.
    pub fn eval(&self, alpha: &C) -> Option<C> {
        let at = |p: &[C]| {
            let mut acc = C::zero();
            for c in p.iter().rev() {
                acc = acc.mul(alpha);
                acc.add_assign(c);
            }
            acc
        };
        let d = at(&self.den);
        (!d.is_zero()).then(|| at(&self.num).mul(&d.inv()))
    }
}

impl<C: Integral> Alpha for ARat<C> {
    fn mul_linears(&self, factors: &Linears) -> Self {
        self.mul_factors(factors)
    }
    /// Nothing: every operation here lands in the normal form already, since
    /// the denominator is dense and a sum that left it unreduced would grow in
    /// degree rather than in a multiset of atoms.
    fn settle(&mut self) {}
}

impl<C: Integral> FromAFrac<C> for ARat<C> {
    fn lift(f: &AFrac<C>) -> Self {
        Self::from_afrac(f)
    }
}

/// Structural, and sound because the representation is canonical: the
/// denominator is primitive and positively led, the numerator coprime to it,
/// and such a pair is unique. This is the property [`AFrac`] gives up in
/// exchange for the factored form, and the reason its own `PartialEq`
/// cross-multiplies.
impl<C: Integral> PartialEq for ARat<C> {
    fn eq(&self, other: &Self) -> bool {
        self.num == other.num && self.den == other.den
    }
}

impl<C: Integral> Eq for ARat<C> {}

impl<C: Integral> Ring for ARat<C> {
    fn zero() -> Self {
        ARat {
            num: Vec::new(),
            den: vec![C::one()],
        }
    }
    fn one() -> Self {
        ARat {
            num: vec![C::one()],
            den: vec![C::one()],
        }
    }
    fn is_zero(&self) -> bool {
        self.num.is_empty()
    }
    /// `(a·d + c·b)/(b·d)`, reduced. Unlike [`AFrac`], which leaves a running
    /// sum unreduced because a cancellation can only be decided once the sum is
    /// complete, this reduces every time: the denominator is dense, so an
    /// unreduced sum grows in degree rather than in a multiset of atoms.
    fn add_assign(&mut self, other: &Self) {
        if other.is_zero() {
            return;
        }
        if self.is_zero() {
            *self = other.clone();
            return;
        }
        // Over the lcm, not the product. The two are the same element either
        // way, but the product's degree is the sum of the degrees and its
        // coefficients the product of the coefficients — a sum over the p(n)
        // shapes of one degree reaches it, and a degree-7 Jack coefficient
        // overflows `i128` there while its lcm form does not. The gcd would
        // take the same factor back out afterwards; the point is not to form
        // it (Knuth, TAOCP 4.5.1, for the integer case).
        let g = pgcd(&self.den, &other.den);
        let (mine, theirs) = if g.len() > 1 {
            (divide_exact(&other.den, &g), divide_exact(&self.den, &g))
        } else {
            (other.den.clone(), self.den.clone())
        };
        let num = padd(&pmul(&self.num, &mine), &pmul(&other.num, &theirs));
        let den = pmul(&self.den, &mine);
        *self = ARat::from_parts(num, den);
    }
    /// Cross-cancelled before the products are formed, for the reason
    /// [`Ring::add_assign`] sums over the lcm.
    fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let (mut an, mut bd) = (self.num.clone(), other.den.clone());
        let g1 = pgcd(&an, &bd);
        if g1.len() > 1 {
            an = divide_exact(&an, &g1);
            bd = divide_exact(&bd, &g1);
        }
        let (mut cn, mut dd) = (other.num.clone(), self.den.clone());
        let g2 = pgcd(&cn, &dd);
        if g2.len() > 1 {
            cn = divide_exact(&cn, &g2);
            dd = divide_exact(&dd, &g2);
        }
        ARat::from_parts(pmul(&an, &cn), pmul(&bd, &dd))
    }
    fn neg(&self) -> Self {
        ARat {
            num: self.num.iter().map(C::neg).collect(),
            den: self.den.clone(),
        }
    }
    fn from_i64(n: i64) -> Self {
        Self::from_coeffs(vec![C::from_i64(n)])
    }
    fn from_u128(n: u128) -> Self {
        Self::from_coeffs(vec![C::from_u128(n)])
    }
    fn from_i128(n: i128) -> Self {
        Self::from_coeffs(vec![C::from_i128(n)])
    }
}

/// ℚ(α) contains ℚ even when `C` is only the integers: dividing by an integer
/// multiplies the denominator by it, which is exact and asks nothing of `C`.
/// This is the same trick [`AFrac`]'s integer `scale` plays, and it is why the
/// Jack engines can run over `i128` and still be handed to `s → p`.
impl<C: Integral> crate::coeff::QAlgebra for ARat<C> {
    fn div_u128(&self, n: u128) -> Self {
        assert!(n != 0, "division of an element of Q(alpha) by zero");
        ARat::from_parts(self.num.clone(), pscale(&self.den, &C::from_u128(n)))
    }
}

/// ℚ(α) is a field however integral `C` is, since inverting only swaps the two
/// parts.
impl<C: Integral> Field for ARat<C> {
    /// # Panics
    ///
    /// Panics on zero.
    fn inv(&self) -> Self {
        assert!(!self.is_zero(), "inverse of zero in Q(alpha)");
        ARat::from_parts(self.den.clone(), self.num.clone())
    }
}

impl<C: Integral> Plethystic for ARat<C> {
    /// `p_n` raises the variable: α ↦ α^n, which spreads a dense polynomial's
    /// coefficients over every nth slot. This is the map [`AFrac`] cannot
    /// carry, and the reason this type exists.
    ///
    /// The normal form survives untouched: spreading changes no coefficient, so
    /// both parts stay primitive and positively led, and `p(α^n)` and `q(α^n)`
    /// are coprime whenever `p` and `q` are — the resultant of the raised pair
    /// is a power of the resultant of the original, and so nonzero.
    ///
    /// The coefficients are not pushed through a Frobenius of their own, and
    /// nothing is missed by that: [`Integral`] is an *integer* ring, so its
    /// elements are constants and α is the only variable there is to raise.
    ///
    /// # Panics
    ///
    /// Panics if `n == 0`, which is not a raising of variables: every exponent
    /// would land on zero, an evaluation at `α = 1` rather than a Frobenius.
    fn frobenius(&self, n: u32) -> Self {
        assert!(
            n > 0,
            "the plethystic Frobenius needs n >= 1: p_0 does not raise variables"
        );
        let spread = |p: &[C]| {
            if p.is_empty() {
                return Vec::new();
            }
            let mut out = vec![C::zero(); (p.len() - 1) * n as usize + 1];
            for (k, c) in p.iter().enumerate() {
                out[k * n as usize] = c.clone();
            }
            out
        };
        ARat {
            num: spread(&self.num),
            den: spread(&self.den),
        }
    }
}

impl<C: Integral> core::fmt::Display for ARat<C> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        fn poly<C: Integral>(f: &mut core::fmt::Formatter<'_>, p: &[C]) -> core::fmt::Result {
            if p.is_empty() {
                return f.write_str("0");
            }
            let mut first = true;
            for (k, c) in p.iter().enumerate().rev() {
                if c.is_zero() {
                    continue;
                }
                if !first {
                    f.write_str(" + ")?;
                }
                first = false;
                match k {
                    0 => write!(f, "{c:?}")?,
                    1 => write!(f, "{c:?}*alpha")?,
                    _ => write!(f, "{c:?}*alpha^{k}")?,
                }
            }
            if first {
                f.write_str("0")?;
            }
            Ok(())
        }
        f.write_str("(")?;
        poly(f, &self.num)?;
        f.write_str(")")?;
        if self.den.len() > 1 || !self.den.first().is_some_and(|c| *c == C::one()) {
            f.write_str("/(")?;
            poly(f, &self.den)?;
            f.write_str(")")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    type R = ARat<i128>;

    /// `1/(α+1) + 1/(α+2)` is `(2α+3)/(α²+3α+2)`, and the sum must arrive in
    /// the normal form rather than over the product of the two denominators
    /// with a common factor left in.
    #[test]
    fn arithmetic_lands_in_the_normal_form() {
        let a = R::from_parts(vec![1], vec![1, 1]);
        let b = R::from_parts(vec![1], vec![2, 1]);
        let mut sum = a.clone();
        sum.add_assign(&b);
        assert_eq!(sum, R::from_parts(vec![3, 2], vec![2, 3, 1]));
        // (α+1)/(α+1) is 1, not a pair with a common factor.
        assert_eq!(R::from_parts(vec![1, 1], vec![1, 1]), <R as Ring>::one());
        // 2/(2α+2) reduces to 1/(α+1): the shared content goes too.
        assert_eq!(R::from_parts(vec![2], vec![2, 2]), a);
        // A negative denominator is carried into the numerator.
        assert_eq!(R::from_parts(vec![1], vec![-1, -1]), a.neg());
    }

    /// Equality is structural here, which is only sound because the form is
    /// canonical — the same element written four ways.
    #[test]
    fn equal_elements_are_structurally_equal() {
        let want = R::from_parts(vec![1], vec![1, 1]);
        for (num, den) in [
            (vec![2], vec![2, 2]),
            (vec![-3], vec![-3, -3]),
            (vec![1, 1], vec![1, 2, 1]),
        ] {
            assert_eq!(R::from_parts(num, den), want, "1/(alpha+1) rewritten");
        }
    }

    /// The Frobenius is what [`AFrac`] cannot carry: `α + 1` must become
    /// `α² + 1`, which is irreducible, and not `(α + 1)²`.
    #[test]
    fn frobenius_raises_alpha_out_of_the_linear_class() {
        let a = R::from_parts(vec![1], vec![1, 1]); // 1/(α+1)
        assert_eq!(
            a.frobenius(2),
            R::from_parts(vec![1], vec![1, 0, 1]),
            "1/(alpha+1) at n = 2 is 1/(alpha^2+1)"
        );
        assert_eq!(a.frobenius(1), a, "n = 1 is the identity");
        let b = R::from_parts(vec![0, 1], vec![2, 1]); // α/(α+2)
        for n in 1..4 {
            assert_eq!(
                a.mul(&b).frobenius(n),
                a.frobenius(n).mul(&b.frobenius(n)),
                "multiplicative at n = {n}"
            );
            let mut sum = a.clone();
            sum.add_assign(&b);
            let mut raised = a.frobenius(n);
            raised.add_assign(&b.frobenius(n));
            assert_eq!(sum.frobenius(n), raised, "additive at n = {n}");
        }
    }

    /// Every [`AFrac`] is an element of ℚ(α), and reading it in must agree with
    /// the arithmetic done here — including the integer scale, which `AFrac`
    /// keeps outside the atoms.
    #[test]
    fn the_factored_form_reads_in() {
        // 6/(2·(α+1)(α+2)) = 3/((α+1)(α+2))
        let f = AFrac::<i128>::from_coeffs(vec![6])
            .div_linear(1, 1)
            .div_linear(1, 2)
            .div_int(2);
        assert_eq!(R::from_afrac(&f), R::from_parts(vec![3], vec![2, 3, 1]));
        assert_eq!(
            R::from_afrac(&<AFrac<i128> as Ring>::zero()),
            <R as Ring>::zero()
        );
    }

    /// The ring axioms, and that division by an integer is exact.
    #[test]
    fn ring_axioms_hold() {
        use crate::coeff::QAlgebra;
        let a = R::from_parts(vec![1, 2], vec![1, 1]);
        let b = R::from_parts(vec![3], vec![0, 1]);
        let c = R::from_parts(vec![1, 0, 1], vec![5, 1]);
        let mut lhs = a.mul(&b);
        lhs.add_assign(&a.mul(&c));
        let mut sum = b.clone();
        sum.add_assign(&c);
        assert_eq!(lhs, a.mul(&sum), "distributive");
        assert_eq!(a.mul(&b), b.mul(&a), "commutative");
        assert_eq!(a.mul(&a.inv()), <R as Ring>::one(), "a field");
        let mut back = a.div_u128(7);
        for _ in 0..6 {
            let more = a.div_u128(7);
            back.add_assign(&more);
        }
        assert_eq!(back, a, "seven sevenths");
    }

    /// The gcd is the primitive-part one, so a pair with a nontrivial common
    /// factor of positive degree reduces rather than only losing its content.
    #[test]
    fn the_polynomial_gcd_cancels_a_shared_factor() {
        // (α² − 1)/(α² + 2α + 1) = (α − 1)/(α + 1)
        assert_eq!(
            R::from_parts(vec![-1, 0, 1], vec![1, 2, 1]),
            R::from_parts(vec![-1, 1], vec![1, 1])
        );
        // and one that shares nothing stays as it is
        let coprime = R::from_parts(vec![1, 0, 1], vec![1, 1]);
        assert_eq!(coprime.parts(), (&[1i128, 0, 1][..], &[1i128, 1][..]));
    }
}
