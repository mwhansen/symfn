//! ℚ(α) with no restriction on the denominator.
//!
//! [`AFrac`] holds ℚ(α) in the shape the Jack engines
//! produce it: a numerator over an integer times a product of *linear* forms
//! `uα + v`. Every denominator those engines build is a product of hooks, so
//! the class is closed and the factored form is what keeps the arithmetic
//! cheap — matching factors cancel before anything is expanded.
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
//! So this type keeps the denominator **dense and monic, coprime to the
//! numerator**. That is a normal form: for a nonzero element there is exactly
//! one such pair, so [`PartialEq`] is structural — unlike
//! [`AFrac`], whose integer content is not canonical and
//! which therefore cross-multiplies.
//!
//! ## The bound is [`Field`], and that is what the gcd needs
//!
//! Euclid's algorithm divides by leading coefficients, so the coefficient ring
//! must invert them. [`AFrac`] runs over integer rings
//! because dividing by a *linear* form is a recurrence with an exact integer
//! division at each step; a general gcd has no such route. The rings this is
//! instantiated over are [`GuardedRat`](crate::guard::GuardedRat) and
//! `BigRational`, the same pair every other rational-coefficient boundary
//! escalates through.
//!
//! ## References
//!
//! - \[GCL\] Geddes, Czapor, Labahn, *Algorithms for Computer Algebra*, ch. 2.

use crate::afrac::{AFrac, Linears};
use crate::coeff::{Field, Plethystic, QAlgebra, Ring};

/// An element of ℚ(α), numerator and denominator dense in α.
///
/// The invariants, together the normal form the module doc describes:
///
/// - `num` has no trailing zeros; empty is 0.
/// - `den` is monic and never empty. When `num` is empty, `den` is `[1]`.
/// - `num` and `den` are coprime in `ℚ[α]`.
///
/// ```
/// use symfn::{ARat, Plethystic, Rational, Ring};
///
/// let r = |n| Rational::from_int(n);
/// // 1/(α + 1)
/// let a = ARat::<Rational>::from_parts(vec![r(1)], vec![r(1), r(1)]);
/// // and the same element written with a common factor left in
/// let same = ARat::<Rational>::from_parts(vec![r(2), r(2)], vec![r(2), r(4), r(2)]);
/// assert_eq!(a, same);
///
/// assert_eq!(
///     a.frobenius(2),
///     ARat::<Rational>::from_parts(vec![r(1)], vec![r(1), r(0), r(1)]),
/// );
/// ```
///
/// `p_2` sends `1/(α + 1)` to `1/(α² + 1)`, **not** to `1/(α + 1)²`: it raises
/// the variable rather than the whole form. `α² + 1` is irreducible over ℚ,
/// which is why this cannot be an [`AFrac`].
#[derive(Clone, Debug)]
pub struct ARat<C: Field> {
    num: Vec<C>,
    den: Vec<C>,
}

fn trim<C: Ring>(v: &mut Vec<C>) {
    while v.last().is_some_and(C::is_zero) {
        v.pop();
    }
}

/// `a + b`, dense.
fn padd<C: Field>(a: &[C], b: &[C]) -> Vec<C> {
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

/// `a · b`, dense. Schoolbook: the degrees here are the α-degrees of Jack
/// coefficients, which stay small next to the partition counts around them.
fn pmul<C: Field>(a: &[C], b: &[C]) -> Vec<C> {
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

fn pneg<C: Field>(a: &[C]) -> Vec<C> {
    a.iter().map(C::neg).collect()
}

/// `a · c` for a scalar `c`.
fn pscale<C: Field>(a: &[C], c: &C) -> Vec<C> {
    if c.is_zero() {
        return Vec::new();
    }
    a.iter().map(|x| x.mul(c)).collect()
}

/// `(quotient, remainder)` of `a` by `b`.
///
/// # Panics
///
/// Panics if `b` is zero.
fn pdivrem<C: Field>(a: &[C], b: &[C]) -> (Vec<C>, Vec<C>) {
    assert!(!b.is_empty(), "division of a polynomial in alpha by zero");
    if a.len() < b.len() {
        return (Vec::new(), a.to_vec());
    }
    let inv = b[b.len() - 1].inv();
    let mut r = a.to_vec();
    let mut q = vec![C::zero(); a.len() - b.len() + 1];
    while r.len() >= b.len() && !r.is_empty() {
        let shift = r.len() - b.len();
        let factor = r[r.len() - 1].mul(&inv);
        q[shift] = factor.clone();
        for (i, y) in b.iter().enumerate() {
            let t = y.mul(&factor).neg();
            r[shift + i].add_assign(&t);
        }
        trim(&mut r);
    }
    trim(&mut q);
    (q, r)
}

/// The monic gcd of `a` and `b`, `[1]` when they are coprime.
///
/// Euclid over a field, which is the simple case — the coefficient growth that
/// makes a fraction-free variant worth having (\[GCL\] ch. 2) is bounded here
/// by the α-degree, which is at most the degree of the symmetric function.
///
/// `gcd(0, 0)` is `[1]`, which is what the callers here want: a zero numerator
/// carries no denominator at all.
fn pgcd<C: Field>(a: &[C], b: &[C]) -> Vec<C> {
    let mut x = a.to_vec();
    let mut y = b.to_vec();
    while !y.is_empty() {
        let (_, r) = pdivrem(&x, &y);
        x = y;
        y = r;
    }
    if x.is_empty() {
        return vec![C::one()];
    }
    let inv = x[x.len() - 1].inv();
    pscale(&x, &inv)
}

impl<C: Field> ARat<C> {
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

    /// Restore the normal form: divide out the gcd and make the denominator
    /// monic.
    fn normalize(&mut self) {
        if self.num.is_empty() {
            self.den = vec![C::one()];
            return;
        }
        let g = pgcd(&self.num, &self.den);
        if g.len() > 1 {
            self.num = pdivrem(&self.num, &g).0;
            self.den = pdivrem(&self.den, &g).0;
        }
        let inv = self.den[self.den.len() - 1].inv();
        self.num = pscale(&self.num, &inv);
        self.den = pscale(&self.den, &inv);
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
        (self.den.len() == 1).then_some(self.num)
    }

    /// The polynomial `uα + v`.
    pub fn linear(u: u32, v: u32) -> Self {
        Self::from_coeffs(vec![C::from_u128(v as u128), C::from_u128(u as u128)])
    }

    /// Multiply by `∏ (uα + v)^m`, negative `m` meaning a denominator factor.
    ///
    /// The `Q` and `J` normalizers are products of hooks and arrive in exactly
    /// this shape, so this is the same seam
    /// [`AFrac::mul_factors`](crate::afrac::AFrac::mul_factors) offers — but
    /// here it is only a convenience over [`Ring::mul`], since nothing stays
    /// factored.
    ///
    /// # Panics
    ///
    /// Panics if a key with a nonzero multiplicity is `(0, 0)`, the zero form.
    pub fn mul_factors(&self, factors: &Linears) -> Self {
        let mut out = self.clone();
        for (&(u, v), &m) in factors {
            if m == 0 {
                continue;
            }
            assert!(u > 0 || v > 0, "0*alpha + 0 is zero");
            let f = Self::linear(u, v);
            for _ in 0..m.unsigned_abs() {
                out = if m > 0 {
                    out.mul(&f)
                } else {
                    out.mul(&f.inv())
                };
            }
        }
        out
    }

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

impl<C: Field> ARat<C> {
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

/// Structural, and sound because the representation is canonical: the
/// denominator is monic and coprime to the numerator, and such a pair is
/// unique. This is the property [`AFrac`] gives up in exchange for the
/// factored form, and the reason its own `PartialEq` cross-multiplies.
impl<C: Field> PartialEq for ARat<C> {
    fn eq(&self, other: &Self) -> bool {
        self.num == other.num && self.den == other.den
    }
}

impl<C: Field> Eq for ARat<C> {}

impl<C: Field> Ring for ARat<C> {
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
    /// complete, this reduces every time: the denominator here is dense, so an
    /// unreduced sum grows in degree rather than in a multiset of atoms.
    fn add_assign(&mut self, other: &Self) {
        if other.is_zero() {
            return;
        }
        if self.is_zero() {
            *self = other.clone();
            return;
        }
        let num = padd(&pmul(&self.num, &other.den), &pmul(&other.num, &self.den));
        let den = pmul(&self.den, &other.den);
        *self = ARat::from_parts(num, den);
    }
    fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        ARat::from_parts(pmul(&self.num, &other.num), pmul(&self.den, &other.den))
    }
    fn neg(&self) -> Self {
        ARat {
            num: pneg(&self.num),
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

/// ℚ(α) is a field, so this divides by inverting the integer — exactly, since
/// the integer is a nonzero constant polynomial.
impl<C: Field> QAlgebra for ARat<C> {
    fn div_u128(&self, n: u128) -> Self {
        assert!(n != 0, "division of an element of Q(alpha) by zero");
        ARat::from_parts(self.num.clone(), pscale(&self.den, &C::from_u128(n)))
    }
}

impl<C: Field> Field for ARat<C> {
    /// # Panics
    ///
    /// Panics on zero.
    fn inv(&self) -> Self {
        assert!(!self.is_zero(), "inverse of zero in Q(alpha)");
        ARat::from_parts(self.den.clone(), self.num.clone())
    }
}

impl<C: Field + Plethystic> Plethystic for ARat<C> {
    /// `p_n` raises the variable: α ↦ α^n, which spreads a dense polynomial's
    /// coefficients over every nth slot. This is the map
    /// [`AFrac`] cannot carry, and the reason this type
    /// exists.
    ///
    /// The normal form survives: raising leaves the denominator monic, and
    /// `p(α^n)` and `q(α^n)` are coprime whenever `p` and `q` are, since the
    /// resultant of the raised pair is a power of the resultant of the
    /// original and so nonzero.
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
                out[k * n as usize] = c.frobenius(n);
            }
            out
        };
        ARat {
            num: spread(&self.num),
            den: spread(&self.den),
        }
    }
}

impl<C: Field> core::fmt::Display for ARat<C> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        fn poly<C: Field>(f: &mut core::fmt::Formatter<'_>, p: &[C]) -> core::fmt::Result {
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
    use crate::coeff::Rational;

    type R = ARat<Rational>;

    fn r(n: i128) -> Rational {
        Rational::from_int(n)
    }

    /// `1/(α+1) + 1/(α+2)` is `(2α+3)/(α²+3α+2)`, and the sum must arrive in
    /// the normal form rather than over the product of the two denominators
    /// with a common factor left in.
    #[test]
    fn arithmetic_lands_in_the_normal_form() {
        let a = R::from_parts(vec![r(1)], vec![r(1), r(1)]);
        let b = R::from_parts(vec![r(1)], vec![r(2), r(1)]);
        let mut sum = a.clone();
        sum.add_assign(&b);
        assert_eq!(sum, R::from_parts(vec![r(3), r(2)], vec![r(2), r(3), r(1)]));
        // (α+1)/(α+1) is 1, not a pair with a common factor.
        let one = R::from_parts(vec![r(1), r(1)], vec![r(1), r(1)]);
        assert_eq!(one, <R as Ring>::one());
        // 2/(2α+2) reduces to 1/(α+1): the denominator is monic.
        let scaled = R::from_parts(vec![r(2)], vec![r(2), r(2)]);
        assert_eq!(scaled, a);
    }

    /// Equality is structural here, which is only sound because the form is
    /// canonical — the same two elements written four ways.
    #[test]
    fn equal_elements_are_structurally_equal() {
        let want = R::from_parts(vec![r(1)], vec![r(1), r(1)]);
        for (num, den) in [
            (vec![r(2)], vec![r(2), r(2)]),
            (vec![r(-3)], vec![r(-3), r(-3)]),
            (vec![r(1), r(1)], vec![r(1), r(2), r(1)]),
        ] {
            assert_eq!(R::from_parts(num, den), want, "1/(alpha+1) rewritten");
        }
    }

    /// The Frobenius is what `AFrac` cannot carry: `α + 1` must become
    /// `α² + 1`, which is irreducible, and not `(α + 1)²`.
    #[test]
    fn frobenius_raises_alpha_out_of_the_linear_class() {
        let a = R::from_parts(vec![r(1)], vec![r(1), r(1)]); // 1/(α+1)
        assert_eq!(
            a.frobenius(2),
            R::from_parts(vec![r(1)], vec![r(1), r(0), r(1)]),
            "1/(alpha+1) at n = 2 is 1/(alpha^2+1)"
        );
        assert_eq!(a.frobenius(1), a, "n = 1 is the identity");
        let b = R::from_parts(vec![r(0), r(1)], vec![r(2), r(1)]); // α/(α+2)
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

    /// Every `AFrac` is an element of ℚ(α), and reading it in must agree with
    /// the arithmetic done here — including the integer scale, which `AFrac`
    /// keeps outside the atoms.
    #[test]
    fn the_factored_form_reads_in() {
        // 6/(2·(α+1)(α+2)) = 3/((α+1)(α+2))
        let f = AFrac::<Rational>::from_coeffs(vec![r(6)])
            .div_linear(1, 1)
            .div_linear(1, 2)
            .div_int(2);
        let want = R::from_parts(vec![r(3)], vec![r(2), r(3), r(1)]);
        assert_eq!(R::from_afrac(&f), want);
        assert_eq!(
            R::from_afrac(&<AFrac<Rational> as Ring>::zero()),
            <R as Ring>::zero()
        );
    }

    /// The ring axioms, and that division by an integer is exact.
    #[test]
    fn ring_axioms_hold() {
        let a = R::from_parts(vec![r(1), r(2)], vec![r(1), r(1)]);
        let b = R::from_parts(vec![r(3)], vec![r(0), r(1)]);
        let c = R::from_parts(vec![r(1), r(0), r(1)], vec![r(5), r(1)]);
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
        assert_eq!(a.eval(&r(1)), Some(Rational::new(3, 2)), "at alpha = 1");
        assert_eq!(b.eval(&r(0)), None, "a pole is reported, not rounded");
    }
}
