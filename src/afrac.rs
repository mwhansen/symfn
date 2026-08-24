//! ℚ(α), with denominators that factor into linear forms `uα + v` and a
//! general factor beside them for the one operation that leaves that class.
//!
//! Every scalar Jack polynomials produce is a ratio of **integer-linear forms
//! `uα + v`**: the two hooks
//!
//! ```text
//!   h_low(s) = α·a(s)     + l(s) + 1          h_up(s) = α·(a(s)+1) + l(s)
//! ```
//!
//! are linear, ψ is a product of ratios of them, and the Laplace–Beltrami
//! denominator `E(κ) − E(λ) = (n(κ')−n(λ'))·α + (n(λ)−n(κ))` is a single one.
//! So the representation [`Frac`](crate::frac::Frac) uses for ℚ(q,t) works
//! here too, on an atom family that behaves better in every way below.
//!
//! ## Why this family is tamer than the binomial one
//!
//! Normalize each atom to be **primitive** — `gcd(u, v) = 1`, the content
//! pulled out into an integer — and three things become true that are false for
//! `1 − qᵃtᵇ`:
//!
//! 1. **Distinct primitive linear forms are irreducible and pairwise coprime**
//!    in `ℚ[α]`. `Frac`'s docs explain why *its* factored form is not canonical
//!    (`1 − q²` is reducible, so `(1+q)/(1−q²)` and `1/(1−q)` are the same
//!    element stored differently). Nothing of the sort happens here.
//! 2. **Taking the larger exponent of each atom gives the exact lcm** of two
//!    denominators, not merely a common multiple, so addition grows the
//!    denominator no more than it has to — unlike the `1 − qᵃtᵇ` family, where
//!    [`Ratio`](crate::deltaop::Ratio) addition settles for a common multiple
//!    because `q² − t²` factors.
//! 3. **A failed cancellation is detected on its first step.** Synthetic
//!    division by `uα + v` starts at the top coefficient and needs `u` to
//!    divide it; by Gauss's lemma that is *necessary* for divisibility in
//!    `ℚ[α]`, so the usual failure exits immediately. The (q,t) engine has the
//!    opposite problem — `divide_exact` runs its failures to completion — and
//!    needs a bespoke necessary-condition pre-pass (`frac::diff_may_divide`).
//!
//! ⚠️ **Skip the primitive part and the answers leave the ring.** Eigenvalue
//! differences are genuinely non-primitive — κ = (2,2), λ = (1,1,1,1) gives
//! `E(κ)−E(λ) = 2α+4` — and dividing `ℤ[α]` by a non-primitive linear form
//! leaves `ℤ[α]`.
//! Splitting off the content is what makes Gauss's lemma apply and the quotient
//! integral.
//!
//! ## The representation
//!
//! ```text
//!   value = num(α) / (scale · ∏ (uα + v)^m · tail(α))
//! ```
//!
//! with `num` a **dense** `Vec<C>` (a `QtPoly` with a dead `t` would be a
//! sparse two-variable key for a dense univariate object), the atoms a
//! `BTreeMap`, `scale` a positive integer, and `tail` empty for the polynomial
//! 1. The integer denominator is what lets `C = i128` stay integral: a value
//! like `1/(2(α+2))` has no home in `ℤ[α]` otherwise.
//!
//! ## The tail, and why it is beside the atoms rather than instead of them
//!
//! One operation leaves the linear class: the plethystic Frobenius raises the
//! variable, α ↦ α^n, so an atom `α + 1` becomes `α² + 1` — irreducible over
//! ℚ. `tail` is where that goes. It is **empty for every value the Jack
//! engines build**, since their denominators are products of hooks, so nothing
//! below costs anything until a plethysm runs.
//!
//! Making the *whole* denominator dense instead was tried and measured
//! (`docs/record/jack.md`). It does not work: the factored form never runs a
//! polynomial gcd — addition takes the lcm of two multisets and cancellation
//! is an exact division by a linear form — and a dense denominator replaces
//! that with a gcd whose pseudo-remainders leave `i128` at degree 7, where the
//! answers are still 13 bits wide. Keeping the atoms keeps that route for the
//! values that never leave the class, which is nearly all of them.
//!
//! **The tail carries no rational root**, which is what keeps the split
//! well-defined: [`AFrac::frobenius`] extracts the linear factors of
//! `uα^n + v` into the atoms as it produces them, a product of root-free
//! polynomials is root-free, and so is what a cancellation leaves behind.
//! Without that, `α³ + 8 = (α + 2)(α² − 2α + 4)` would be storable two ways.
//! Cancelling a tail against the numerator is the one place a polynomial gcd
//! runs here, and it needs [`Integral`] — a gcd on `C` — which is why the
//! bound is that rather than [`Ring`].
//!
//! `C` must implement [`Ring::div_exact`] and [`Integral::gcd`] faithfully —
//! ℤ-like (`i128`, `BigInt`) or a field ([`Rational`](crate::coeff::Rational),
//! `BigRational`, where every nonzero element is a unit and the gcd is 1).
//! A ring that declines every division would make [`AFrac::reduce`] a no-op and
//! the denominators grow without bound; nothing would be *wrong*, but nothing
//! would cancel either.
//!
//! ## Equality cross-multiplies anyway
//!
//! The *atom* half of the representation is canonical, which is the property
//! `docs/record/jack.md` identifies. The **integer content** half is
//! not, because cancelling it needs a gcd inside `C` that the [`Ring`] trait
//! does not offer. [`AFrac::reduce`] gets it by trial-dividing by the prime
//! factors of `scale` — which are always tiny, since `scale` is only ever built
//! from hook contents and eigenvalue gcds — but the search gives up on a large
//! prime residue rather than factoring it. Missing a cancellation is a size
//! inefficiency and never an error, so the give-up is safe; it does mean two
//! representations of one element can survive, and cross-multiplying is what
//! makes that harmless. One extra scalar multiply is cheaper than demanding a
//! gcd from `C`.

// Polynomial degrees and the `i32` atom offsets, bounded by the degree.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::BTreeMap;

use crate::coeff::{Field, Integral, Plethystic, QAlgebra, Ring};

/// A **primitive** integer-linear atom `uα + v`: `gcd(u, v) = 1` and not both
/// zero. `(1, 0)` is `α` itself, which does occur — `h_up` at a cell with no
/// leg is `α·(a+1)`.
pub type Atom = (u32, u32);

/// Linear forms with signed multiplicities: `(u, v) ↦ m` means `(uα + v)^m`,
/// negative `m` meaning a denominator factor.
///
/// The forms here need **not** be primitive — hooks arrive as `αa + (l+1)` with
/// whatever content they have, and [`AFrac::from_factors`] splits it off. This
/// is the α-analog of the `Factors` multiset `macdonald.rs` passes around, and
/// it is used the same way: accumulate exponents first so matching factors
/// cancel before anything is expanded.
pub type Linears = BTreeMap<(u32, u32), i32>;

/// An element of ℚ(α) whose denominator is an integer times a product of
/// primitive linear forms `uα + v`.
#[derive(Clone, Debug)]
pub struct AFrac<C: Integral> {
    /// Dense in α: `num[k]` multiplies `α^k`. No trailing zeros; empty is 0.
    num: Vec<C>,
    /// The primitive atoms with their multiplicities.
    den: BTreeMap<Atom, u32>,
    /// A further denominator factor that is **not** a product of linear forms,
    /// dense in α, primitive and positively led. Empty means the polynomial 1,
    /// which is what it is for every value the Jack engines build; only the
    /// plethystic Frobenius puts anything here.
    ///
    /// It carries no rational root, and cannot acquire one: the Frobenius
    /// extracts the linear factors of `uα^n + v` into `den` as it produces it,
    /// a product of root-free polynomials is root-free, and so is a factor of
    /// one that a cancellation leaves behind.
    tail: Vec<C>,
    /// A positive integer denominator. Always ≥ 1.
    scale: u128,
}

fn gcd32(mut a: u32, mut b: u32) -> u32 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// Split `uα + v` into an integer content and a primitive atom.
///
/// `None` for the atom means the form was a constant (`u == 0`), so the whole
/// thing *is* the content — `h_low` at a cell with no arm is exactly that.
fn split(u: u32, v: u32) -> (u128, Option<Atom>) {
    assert!(u > 0 || v > 0, "0·α + 0 is zero");
    let g = gcd32(u, v);
    if u == 0 {
        (g as u128, None) // gcd(0, v) = v
    } else {
        (g as u128, Some((u / g, v / g)))
    }
}

/// Trial-division bound for the scalar denominator.
///
/// Every prime that actually appears in an `AFrac` scalar comes from a hook
/// content (`≤ 2n`) or a `gcd` of two eigenvalue differences (`≤ n²`), so this
/// is never reached in practice. It exists so a pathological scalar costs a
/// bounded number of divisions instead of a factorization.
const TRIAL_BOUND: u128 = 1 << 10;

// ---------------------------------------------------------------------- tail
//
// The general-denominator half. Everything below runs only when a `tail` is
// present, which is only after a plethysm: the Jack engines build denominators
// out of hooks, which are linear, and never reach it.

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

/// The product of two tails, either of which may be empty for 1.
fn mul_tails<C: Integral>(a: &[C], b: &[C]) -> Vec<C> {
    if a.is_empty() {
        return b.to_vec();
    }
    if b.is_empty() {
        return a.to_vec();
    }
    pmul(a, b)
}

/// `(what a picks up, what b picks up, their lcm)` for two tails.
///
/// The three cheap cases are the ones that happen: both empty (every value the
/// Jack engines build), one empty, or the same tail on both sides — which is
/// what a sum of coefficients from one plethysm looks like. Only two genuinely
/// different tails pay for a polynomial gcd.
fn tail_lcm<C: Integral>(a: &[C], b: &[C]) -> (Vec<C>, Vec<C>, Vec<C>) {
    if a == b {
        return (Vec::new(), Vec::new(), a.to_vec());
    }
    if a.is_empty() {
        return (b.to_vec(), Vec::new(), b.to_vec());
    }
    if b.is_empty() {
        return (Vec::new(), a.to_vec(), a.to_vec());
    }
    let g = pgcd(a, b);
    let for_a = divide_exact(b, &g);
    let for_b = divide_exact(a, &g);
    let lcm = pmul(a, &for_a);
    (for_a, for_b, lcm)
}

/// Drop trailing zero coefficients, so the degree is honest and `is_zero` is a
/// length test.
fn trim<C: Ring>(v: &mut Vec<C>) {
    while v.last().is_some_and(C::is_zero) {
        v.pop();
    }
}

/// `p · (uα + v)`, dense and **in place**.
///
/// `result[k] = p[k]·v + p[k−1]·u`, so walking *downward* reads `p[k−1]`
/// before it is overwritten and the whole product needs one `push` rather than
/// a fresh buffer. [`AFrac::lift`] runs this once per denominator atom on every
/// addition, so the allocating version costs one heap allocation per atom per
/// addition.
fn mul_linear_in_place<C: Ring>(p: &mut Vec<C>, u: u32, v: u32) {
    if p.is_empty() {
        return;
    }
    let (cu, cv) = (C::from_u128(u as u128), C::from_u128(v as u128));
    p.push(C::zero());
    for k in (0..p.len()).rev() {
        let hi = if k > 0 { p[k - 1].mul(&cu) } else { C::zero() };
        let mut lo = p[k].mul(&cv);
        lo.add_assign(&hi);
        p[k] = lo;
    }
    trim(p);
}

/// `p · (uα + v)`, allocating — the form the tests and [`AFrac::mul_factors`]
/// read most naturally.
fn mul_linear<C: Ring>(p: &[C], u: u32, v: u32) -> Vec<C> {
    let mut out = p.to_vec();
    mul_linear_in_place(&mut out, u, v);
    out
}

/// Exact division of `p` by the **primitive** form `uα + v`, or `None`.
///
/// Synthetic division from the top:
///
/// ```text
///   q[d−1] = p[d]/u,   p[i−1] −= q[i−1]·v,   … ,   remainder p[0]
/// ```
///
/// Two things make this the whole divisibility test and not just the quotient.
/// By Gauss's lemma, a primitive `uα + v` dividing `p` in `ℚ[α]` divides it in
/// `C[α]` — so `u ∤ p[d]` is already a proof of non-divisibility, and it is the
/// *first* thing checked. And the leftover `p[0]` is the remainder, which must
/// vanish.
///
/// `docs/record/jack.md` states the test as the integer root
/// evaluation `Σ_k p_k (−v)^k u^{d−k} = 0`. That is the same predicate at the
/// same asymptotic cost, but it forms `u^d`, which for a degree-12 numerator
/// with `v ≈ 200` (reachable around n = 30) is a 90-bit intermediate on top of
/// the coefficient. Synthetic division never builds one, and returns the
/// quotient in the same pass rather than needing a second.
// Used by the tests, which is why it survives `dead_code`: it is the
// allocating form of `divide_in_place`, and the doc above it is where the
// reason synthetic division beats the record's stated predicate is written
// down (`docs/record/jack.md`).
#[allow(dead_code)]
pub(crate) fn divide_by_linear<C: Ring>(p: &[C], u: u32, v: u32) -> Option<Vec<C>> {
    let mut num = p.to_vec();
    divide_in_place(&mut num, u, v).then_some(num)
}

/// Whether `uα + v` divides `p`, **allocating nothing**.
///
/// The recurrence `p = q·(uα + v)` reads
///
/// ```text
///   q[d−1] = p[d]/u,   q[i−1] = (p[i] − q[i]·v)/u,   p[0] = q[0]·v
/// ```
///
/// and each step needs only the *previous* `q`, never the whole array. So the
/// whole predicate is one running scalar — no quotient is materialized, and a
/// failure costs nothing but the walk.
///
/// That matters because `reduce` trial-divides by every denominator atom and
/// most of those fail: building the quotient before knowing the division
/// succeeds spends two heap allocations on each failure
/// (`docs/record/jack.md`). Same trap as `deltaop`'s `divide_exact` — a cheap
/// failure test that was not actually cheap — reached from the other
/// direction.
fn divides_by_linear<C: Ring>(p: &[C], u: u32, v: u32) -> bool {
    debug_assert!(gcd32(u, v) == 1 && u > 0, "divisor must be primitive in α");
    let d = p.len() - 1;
    let cu = C::from_u128(u as u128);
    let cv = C::from_u128(v as u128);
    // Gauss's lemma: a primitive `uα + v` dividing `p` in ℚ[α] divides it in
    // `C[α]`, so `u ∤ p[d]` refutes it outright — and that is the first step.
    let Some(mut q) = p[d].div_exact(&cu) else {
        return false;
    };
    for i in (1..d).rev() {
        let mut r = p[i].clone();
        r.sub_assign(&q.mul(&cv));
        match r.div_exact(&cu) {
            Some(next) => q = next,
            None => return false,
        }
    }
    // The remainder: p[0] must be exactly q[0]·v.
    let mut r = p[0].clone();
    r.sub_assign(&q.mul(&cv));
    r.is_zero()
}

/// Replace `num` by `num / (uα + v)` if it divides, reporting whether it did.
///
/// Tests first with [`divides_by_linear`], which allocates nothing, and only
/// then rewrites — **in place**, with no second buffer. The rewrite stores
/// `q[i−1]` at index `i`, which never clobbers a coefficient still to be read,
/// and one `remove(0)` shifts the answer down. `num` is untouched on failure.
fn divide_in_place<C: Ring>(num: &mut Vec<C>, u: u32, v: u32) -> bool {
    if num.is_empty() {
        return true; // 0 is divisible by everything
    }
    if num.len() == 1 {
        return false; // a nonzero constant is divisible by no linear form
    }
    if !divides_by_linear(num, u, v) {
        return false;
    }
    let cu = C::from_u128(u as u128);
    let cv = C::from_u128(v as u128);
    let d = num.len() - 1;
    num[d] = num[d]
        .div_exact(&cu)
        .expect("divisibility was just established");
    for i in (1..d).rev() {
        let back = num[i + 1].mul(&cv);
        let mut r = num[i].clone();
        r.sub_assign(&back);
        num[i] = r.div_exact(&cu).expect("divisibility was just established");
    }
    num.remove(0);
    trim(num);
    true
}

impl<C: Integral> AFrac<C> {
    /// A polynomial in α, given densely by its coefficients.
    pub fn from_coeffs(mut num: Vec<C>) -> Self {
        trim(&mut num);
        AFrac {
            num,
            den: BTreeMap::new(),
            tail: Vec::new(),
            scale: 1,
        }
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
        let (content, atom) = split(u, v);
        let mut den = BTreeMap::new();
        if let Some(a) = atom {
            den.insert(a, 1);
        }
        AFrac {
            num: vec![C::one()],
            den,
            tail: Vec::new(),
            scale: content,
        }
    }

    /// `∏ (uα + v)^m` from a signed multiplicity map — the α-analog of
    /// [`Frac::from_factors`](crate::frac::Frac::from_factors), and used the
    /// same way: ψ is built as a multiset of linear forms so that matching
    /// factors cancel *before* any expansion happens.
    ///
    /// Not reduced, on the same policy as [`Ring::add_assign`]: this is called
    /// once per tableau, and reduction is a per-*coefficient* operation.
    ///
    /// # Panics
    ///
    /// Panics if a key with a nonzero multiplicity is `(0, 0)`, the zero form.
    /// The panic is [`mul_factors`](Self::mul_factors)', which this delegates
    /// to.
    pub fn from_factors(factors: &Linears) -> Self {
        <Self as Ring>::one().mul_factors(factors)
    }

    /// Multiply by `∏ (uα + v)^m`, negative `m` meaning a denominator factor.
    ///
    /// The `Q` and `J` normalizers are exactly this shape — products of hooks —
    /// and reaching them through [`Ring::mul`] would expand `H_λ` into a
    /// polynomial first. Here every step is one linear multiply or one exact
    /// linear division.
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
            let (content, atom) = split(u, v);
            let reps = m.unsigned_abs();
            if m > 0 {
                for _ in 0..reps {
                    out.scale_content(content);
                    if let Some((au, av)) = atom {
                        out.num = mul_linear(&out.num, au, av);
                    }
                }
            } else {
                for _ in 0..reps {
                    out.scale *= content;
                    if let Some(a) = atom {
                        *out.den.entry(a).or_insert(0) += 1;
                    }
                }
            }
        }
        out
    }

    /// Divide by the linear form `uα + v`, cancelling immediately if it goes.
    ///
    /// This is the Laplace–Beltrami step: `E(κ) − E(λ)` is one such form, and
    /// the recursion performs exactly one of these per coefficient.
    ///
    /// # Panics
    ///
    /// Panics on `0α + 0`, the zero form.
    pub fn div_linear(&self, u: u32, v: u32) -> Self {
        let (content, atom) = split(u, v);
        let mut out = self.clone();
        out.scale *= content;
        if let Some(a) = atom {
            *out.den.entry(a).or_insert(0) += 1;
            out.reduce_at(a);
        }
        out.content_reduce();
        out
    }

    /// Multiply by an integer, without going through a polynomial product.
    pub fn scale_int(&self, k: i64) -> Self {
        if k == 0 {
            return <Self as Ring>::zero();
        }
        let c = C::from_i64(k);
        let mut out = self.clone();
        for x in &mut out.num {
            *x = x.mul(&c);
        }
        trim(&mut out.num);
        out.content_reduce();
        out
    }

    /// Divide by a positive integer, cancelling it against the numerator's
    /// content where it goes and keeping the rest in the scalar denominator.
    ///
    /// This is the only way to put a `scale` *back*, which is what a
    /// coefficient arriving from outside the crate needs:
    /// [`parts`](Self::parts) hands out `num / (scale · ∏ atoms)`, and
    /// [`from_coeffs`](Self::from_coeffs) and [`mul_factors`](Self::mul_factors)
    /// rebuild everything but the scale.
    ///
    /// # Panics
    ///
    /// Panics on `k = 0`.
    pub fn div_int(&self, k: u128) -> Self {
        assert!(k != 0, "division by zero");
        let mut out = self.clone();
        out.scale *= k;
        out.content_reduce();
        out
    }

    /// Multiply the *value* by the positive integer `k`, cancelling it against
    /// the scalar denominator first.
    fn scale_content(&mut self, k: u128) {
        if k == 1 {
            return;
        }
        let g = crate::coeff::gcd_u128(self.scale, k);
        self.scale /= g;
        let rest = k / g;
        if rest != 1 {
            let c = C::from_u128(rest);
            for x in &mut self.num {
                *x = x.mul(&c);
            }
            trim(&mut self.num);
        }
    }

    /// Divide out every denominator atom that also divides the numerator, and
    /// cancel the integer scalar against the numerator's content.
    ///
    /// Public because [`Ring::add_assign`] deliberately does not do it — the
    /// same policy [`Frac`](crate::frac::Frac) documents. The Jack engines call
    /// it once per coefficient.
    pub fn reduce(&mut self) {
        if self.num.is_empty() {
            self.den.clear();
            self.tail.clear();
            self.scale = 1;
            return;
        }
        let atoms: Vec<Atom> = self.den.keys().copied().collect();
        for a in atoms {
            self.reduce_at(a);
        }
        self.reduce_tail();
        self.content_reduce();
    }

    /// Cancel the tail against the numerator, and normalize what is left.
    ///
    /// The one place a **polynomial** gcd runs. It costs nothing when there is
    /// no tail, which is every value the Jack engines build — only a plethysm
    /// puts one there. See `docs/record/jack.md` for the measurement that put
    /// the general denominator here rather than in place of the atoms: the
    /// dense form's gcd leaves `i128` at degree 7, where the answers are still
    /// 13 bits wide.
    fn reduce_tail(&mut self) {
        if self.tail.is_empty() {
            return;
        }
        let g = pgcd(&self.num, &self.tail);
        if g.len() > 1 {
            self.num = divide_exact(&self.num, &g);
            self.tail = divide_exact(&self.tail, &g);
        }
        // A primitive tail divided by a primitive factor stays primitive
        // (Gauss), so only the sign can need fixing.
        if self.tail[self.tail.len() - 1].is_negative() {
            self.tail = self.tail.iter().map(C::neg).collect();
            self.num = self.num.iter().map(C::neg).collect();
        }
        if self.tail.len() == 1 {
            // A primitive constant is ±1, and the sign has just been fixed.
            self.tail.clear();
        }
    }

    /// Cancel one atom as far as it goes.
    ///
    /// The whole loop allocates nothing: [`divide_in_place`] tests without a
    /// buffer and rewrites without one. See [`divides_by_linear`] for the two
    /// heap allocations a failure would otherwise cost.
    fn reduce_at(&mut self, atom: Atom) {
        let (u, v) = atom;
        while let Some(&m) = self.den.get(&atom) {
            if m == 0 || self.num.is_empty() || !divide_in_place(&mut self.num, u, v) {
                break;
            }
            if m == 1 {
                self.den.remove(&atom);
            } else {
                self.den.insert(atom, m - 1);
            }
        }
    }

    /// Cancel `scale` against the numerator's content, one prime at a time.
    ///
    /// [`Ring`] has no gcd, so the content cannot be *read* — but it can be
    /// *tested*, and only divisors of `scale` are worth testing. `scale` is
    /// built solely from hook contents and eigenvalue-difference gcds, so its
    /// prime factors are bounded by roughly `n²`. The trial-division loop below
    /// finds all of them in a few dozen steps.
    ///
    /// A large prime residue is tried once as a lump and then abandoned rather
    /// than factored. That can miss a cancellation; it cannot
    /// produce a wrong value, which is why [`PartialEq`] cross-multiplies
    /// instead of comparing representations.
    ///
    /// For `C = i128` this is exact integer content extraction. For a field it
    /// always succeeds and drives `scale` to 1, moving the content into the
    /// coefficients — also correct, just a different place to keep it.
    fn content_reduce(&mut self) {
        if self.scale == 1 || self.num.is_empty() {
            return;
        }
        // `residue` is the part of the scale not yet identified as a small
        // prime, and it shrinks whether or not the cancellation succeeds —
        // that separation is the point. Reading the loop bound off
        // `self.scale` instead makes the final step try the *uncancelled*
        // scalar as one lump, so `202 = 2·101` against a numerator of content
        // 101 keeps its 101 forever. Nothing is collected into a `Vec`: this
        // runs on every division, and an allocation here is pure overhead.
        let mut residue = self.scale;
        let mut p = 2u128;
        while p < TRIAL_BOUND && p * p <= residue {
            if residue.is_multiple_of(p) {
                while residue.is_multiple_of(p) {
                    residue /= p;
                }
                self.cancel_scalar(p);
            }
            p += if p == 2 { 1 } else { 2 };
        }
        // Whatever is left is prime or a product of large primes; try it once
        // as a lump rather than factoring it. Missing a cancellation is a size
        // inefficiency, never a wrong value — which is why `PartialEq`
        // cross-multiplies.
        if residue > 1 {
            self.cancel_scalar(residue);
        }
    }

    /// Divide numerator and `scale` by `p` for as long as both allow.
    fn cancel_scalar(&mut self, p: u128) {
        let d = C::from_u128(p);
        while self.scale.is_multiple_of(p) {
            let q: Option<Vec<C>> = self.num.iter().map(|x| x.div_exact(&d)).collect();
            match q {
                Some(q) => {
                    self.num = q;
                    self.scale /= p;
                }
                None => return,
            }
        }
    }

    /// The numerator coefficients — dense in α, `num[k]` multiplying `α^k`,
    /// with no trailing zeros — the denominator's atoms in ascending order, and
    /// the positive integer scalar.
    pub fn parts(&self) -> (&[C], impl Iterator<Item = (&Atom, &u32)>, u128) {
        (&self.num, self.den.iter(), self.scale)
    }

    /// The denominator's general factor, dense in α — empty for the polynomial
    /// 1, which is what it is unless a plethysm put something there.
    ///
    /// Separate from [`parts`](Self::parts) because it is empty almost always,
    /// and because a caller that cannot hold it should have to ask.
    pub fn tail(&self) -> &[C] {
        &self.tail
    }

    /// The value as a polynomial in α, if the denominator cancels away — `None`
    /// if this element genuinely is not one.
    ///
    /// Reduces first, so it answers about the *element*. Callers that know on
    /// mathematical grounds that the answer must be a polynomial — every
    /// coefficient of `J_λ`, by \[KS\] Thm 1.1 — should `expect` it and let a
    /// `None` be the loud failure it is.
    pub fn into_poly(mut self) -> Option<Vec<C>> {
        self.reduce();
        if !self.den.is_empty() || !self.tail.is_empty() {
            return None;
        }
        if self.scale == 1 {
            return Some(self.num);
        }
        // A ℚ-algebra can still absorb a leftover scalar; ℤ cannot, and then
        // this really is not a polynomial over `C`.
        let d = C::from_u128(self.scale);
        self.num.iter().map(|x| x.div_exact(&d)).collect()
    }

    /// The value as `(polynomial in α, positive integer denominator)`, if the
    /// *atoms* cancel away — `None` if a genuine linear form survives.
    ///
    /// The weaker sibling of [`AFrac::into_poly`], and the one the \[GJ\]
    /// pipeline needs: \[DF\] proves `c` and `h` are polynomials in `b` over ℚ,
    /// and \[BD\] proves only `c`'s are integral. Collapsing to `ℚ[b]` and
    /// *reporting* the denominator keeps those two claims distinguishable
    /// instead of failing on the weaker one.
    pub fn into_rational_poly(mut self) -> Option<(Vec<C>, u128)> {
        self.reduce();
        (self.den.is_empty() && self.tail.is_empty()).then_some((self.num, self.scale))
    }

    /// Substitute `α ↦ 1/α`, exactly, staying inside the family.
    ///
    /// ```text
    ///   num(1/α) = α^{−D}·rev(num),      (uα+v)|_{1/α} = (vα+u)/α
    /// ```
    ///
    /// so with `D` the numerator degree and `M` the total atom multiplicity the
    /// value picks up `α^{M−D}` and every atom `(u,v)` becomes `(v,u)` — which
    /// is still primitive, since `gcd` is symmetric. `α` itself, the `(1,0)`
    /// atom, becomes the constant 1 and disappears into the scalar; that is the
    /// only case where the atom count changes, and it is why this cannot be
    /// done by swapping the pairs alone.
    ///
    /// Needed for the `ω_α`-duality law `ω_α P_λ^{(α)} = Q_{λ'}^{(1/α)}`, which
    /// is the one specialization in `docs/record/jack.md` that no
    /// other test reaches — it is the only statement relating `P` to `Q`,
    /// conjugation, and the parameter inversion at once.
    ///
    /// # Panics
    ///
    /// Panics if the value carries a [`tail`](Self::tail). The rewriting above
    /// is exact only because every denominator factor is linear, and a general
    /// factor `T(α)` would need `α^{deg T}·T(1/α)` — reversing its
    /// coefficients, which can make it reducible and so leave the normal form.
    /// Only a plethysm produces a tail, and duality is asked of Jack
    /// polynomials rather than of plethysms, so this is a contract violation
    /// rather than a wall (`docs/policies/failure.md`, R2).
    pub fn invert_alpha(&self) -> Self {
        assert!(
            self.tail.is_empty(),
            "alpha-inversion is written for the factored denominator only"
        );
        if self.num.is_empty() {
            return <Self as Ring>::zero();
        }
        let degree = self.num.len() - 1;
        let mult: u32 = self.den.values().sum();
        let mut num: Vec<C> = self.num.iter().rev().cloned().collect();
        let mut den: BTreeMap<Atom, u32> = BTreeMap::new();
        let mut scale = self.scale;
        for (&(u, v), &m) in &self.den {
            // (u, v) primitive ⟹ (v, u) primitive; but (1,0) ↦ (0,1), the
            // constant 1, which `split` folds into the content.
            let (content, atom) = split(v, u);
            match atom {
                Some(a) => *den.entry(a).or_insert(0) += m,
                None => {
                    debug_assert_eq!(content, 1, "a primitive atom's swap has unit content");
                }
            }
        }
        // The leftover α^{M−D}: a numerator shift one way, an `α` atom the other.
        let shift = mult as i64 - degree as i64;
        match shift.cmp(&0) {
            core::cmp::Ordering::Greater => {
                let mut shifted = vec![C::zero(); shift as usize];
                shifted.append(&mut num);
                num = shifted;
            }
            core::cmp::Ordering::Less => {
                *den.entry((1, 0)).or_insert(0) += (-shift) as u32;
            }
            core::cmp::Ordering::Equal => {}
        }
        trim(&mut num);
        let mut out = AFrac {
            num,
            den,
            tail: Vec::new(),
            scale: 1,
        };
        out.scale = core::mem::replace(&mut scale, 1);
        out.reduce();
        out
    }

    /// Degree in α of the numerator; `None` for zero.
    pub fn degree(&self) -> Option<usize> {
        (!self.num.is_empty()).then(|| self.num.len() - 1)
    }

    /// `self` rewritten over `target`, a multiple of `self.den`, and over the
    /// integer `scale` — the numerator only.
    fn lift(&self, target: &BTreeMap<Atom, u32>, scale: u128, tail_extra: &[C]) -> Vec<C> {
        let mut num = self.num.clone();
        for (&(u, v), &m) in target {
            let extra = m - self.den.get(&(u, v)).copied().unwrap_or(0);
            for _ in 0..extra {
                mul_linear_in_place(&mut num, u, v);
            }
        }
        let k = scale / self.scale;
        if k != 1 {
            let c = C::from_u128(k);
            for x in &mut num {
                *x = x.mul(&c);
            }
            trim(&mut num);
        }
        if !tail_extra.is_empty() {
            num = pmul(&num, tail_extra);
        }
        num
    }
}

impl AFrac<i128> {
    /// Evaluate at a rational α, over ℚ — the integral analogue of
    /// [`AFrac::eval`], which needs `C` to be a field and `i128` is not.
    ///
    /// `None` only if the denominator vanishes, which for `α > 0` cannot
    /// happen: every atom is `uα + v` with `u, v ≥ 0` and not both zero.
    pub fn eval_i128(&self, alpha: &crate::coeff::Rational) -> Option<crate::coeff::Rational> {
        use crate::coeff::{Field, Rational};
        let (num, den, scale) = self.parts();
        if num.is_empty() {
            return Some(Rational::zero());
        }
        let mut acc = Rational::zero();
        for c in num.iter().rev() {
            acc = acc.mul(alpha);
            acc.add_assign(&Rational::from_i128(*c));
        }
        let mut d = Rational::from_u128(scale);
        for (&(u, v), &m) in den {
            let mut f = Rational::from_u128(u as u128).mul(alpha);
            f.add_assign(&Rational::from_u128(v as u128));
            if f.is_zero() {
                return None;
            }
            for _ in 0..m {
                d = d.mul(&f);
            }
        }
        if !self.tail.is_empty() {
            let mut t = Rational::zero();
            for c in self.tail.iter().rev() {
                t = t.mul(alpha);
                t.add_assign(&Rational::from_i128(*c));
            }
            if t.is_zero() {
                return None;
            }
            d = d.mul(&t);
        }
        Some(acc.div(&d))
    }

    /// The value as a polynomial in α with **non-negative integer**
    /// coefficients, if it is one.
    ///
    /// \[KS\] Thm 1.1 (every `J` coefficient, divided by `u_μ`, lies in `ℕ[α]`)
    /// and Stanley's open positivity conjecture are both exactly this question,
    /// and both are laws the engine is *held to* rather than facts it arranges
    /// — see `jack.rs`. Concrete over `i128` because "non-negative" needs an
    /// order that [`Ring`] does not carry.
    pub fn into_natural_poly(self) -> Option<Vec<i128>> {
        let coeffs = self.into_poly()?;
        coeffs.iter().all(|&c| c >= 0).then_some(coeffs)
    }
}

impl<C: Integral + Field> AFrac<C> {
    /// Substitute a value for α; `None` if the denominator vanishes there.
    pub fn eval(&self, alpha: &C) -> Option<C> {
        if self.num.is_empty() {
            return Some(C::zero());
        }
        let mut d = C::from_u128(self.scale);
        for (&(u, v), &m) in &self.den {
            let mut f = C::from_u128(u as u128).mul(alpha);
            f.add_assign(&C::from_u128(v as u128));
            if f.is_zero() {
                return None;
            }
            for _ in 0..m {
                d = d.mul(&f);
            }
        }
        // Horner, from the top.
        let horner = |p: &[C]| {
            let mut acc = C::zero();
            for c in p.iter().rev() {
                acc = acc.mul(alpha);
                acc.add_assign(c);
            }
            acc
        };
        if !self.tail.is_empty() {
            let t = horner(&self.tail);
            if t.is_zero() {
                return None;
            }
            d = d.mul(&t);
        }
        Some(horner(&self.num).div(&d))
    }
}

/// Cross-multiplied — see the module docs. The atoms are canonical; the integer
/// content is not always, so structural comparison would call equal things
/// unequal for `C = Rational`.
impl<C: Integral> PartialEq for AFrac<C> {
    fn eq(&self, other: &Self) -> bool {
        if self.num.is_empty() || other.num.is_empty() {
            return self.num.is_empty() && other.num.is_empty();
        }
        if self.den == other.den && self.scale == other.scale && self.tail == other.tail {
            return self.num == other.num;
        }
        let mut lcm = self.den.clone();
        for (k, &m) in &other.den {
            let e = lcm.entry(*k).or_insert(0);
            *e = (*e).max(m);
        }
        let g = crate::coeff::gcd_u128(self.scale, other.scale);
        let s = self.scale / g * other.scale;
        // Each side picks up the other's tail. That is a common multiple
        // rather than the lcm, which is all an equality test needs: both sides
        // gain the same factor, so it cancels out of the comparison.
        self.lift(&lcm, s, &other.tail) == other.lift(&lcm, s, &self.tail)
    }
}

impl<C: Integral> Eq for AFrac<C> {}

impl<C: Integral> Ring for AFrac<C> {
    fn zero() -> Self {
        AFrac {
            num: Vec::new(),
            den: BTreeMap::new(),
            tail: Vec::new(),
            scale: 1,
        }
    }
    fn one() -> Self {
        AFrac {
            num: vec![C::one()],
            den: BTreeMap::new(),
            tail: Vec::new(),
            scale: 1,
        }
    }
    fn is_zero(&self) -> bool {
        self.num.is_empty()
    }
    fn add_assign(&mut self, other: &Self) {
        if other.is_zero() {
            return;
        }
        if self.is_zero() {
            *self = other.clone();
            return;
        }
        // Over the lcm — and here it really is the lcm, because distinct
        // primitive linear forms are coprime. `Frac` can only manage a common
        // multiple.
        //
        // `self` is the accumulator in every caller that matters, so check
        // before rewriting it: once its denominator has stopped growing, the
        // lcm *is* its denominator and lifting would re-walk the numerator to
        // multiply it by nothing.
        let needs_atoms = other
            .den
            .iter()
            .any(|(k, &m)| self.den.get(k).copied().unwrap_or(0) < m);
        let needs_scale = !self.scale.is_multiple_of(other.scale);
        // The tails' lcm, and what each numerator must pick up to reach it.
        // Both empty is the overwhelmingly common case — every value the Jack
        // engines build — and costs two length tests.
        let (mine, theirs, common) = tail_lcm(&self.tail, &other.tail);
        if needs_atoms || needs_scale || !mine.is_empty() {
            let mut lcm = self.den.clone();
            for (k, &m) in &other.den {
                let e = lcm.entry(*k).or_insert(0);
                *e = (*e).max(m);
            }
            let g = crate::coeff::gcd_u128(self.scale, other.scale);
            let s = self.scale / g * other.scale;
            self.num = self.lift(&lcm, s, &mine);
            self.den = lcm;
            self.scale = s;
            self.tail = common;
        }
        let lifted = other.lift(&self.den, self.scale, &theirs);
        if self.num.len() < lifted.len() {
            self.num.resize(lifted.len(), C::zero());
        }
        for (x, y) in self.num.iter_mut().zip(lifted.iter()) {
            x.add_assign(y);
        }
        trim(&mut self.num);
        // Deliberately not reduced — the `Frac::add_assign` policy, for the
        // same reason: a running sum should be reduced once at the end, not
        // after every term. Cancellation cannot even be decided until the sum
        // is complete.
    }
    fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let mut num = vec![C::zero(); self.num.len() + other.num.len() - 1];
        for (i, a) in self.num.iter().enumerate() {
            if a.is_zero() {
                continue;
            }
            for (j, b) in other.num.iter().enumerate() {
                let v = a.mul(b);
                num[i + j].add_assign(&v);
            }
        }
        trim(&mut num);
        let mut den = self.den.clone();
        for (k, &m) in &other.den {
            *den.entry(*k).or_insert(0) += m;
        }
        let mut f = AFrac {
            num,
            den,
            tail: mul_tails(&self.tail, &other.tail),
            scale: self.scale * other.scale,
        };
        // Reduced, unlike `add_assign` — the same asymmetry `Frac::mul`
        // documents and for the same measured reason: this is where a
        // denominator first gets cut down, before the value is used many times
        // by additions that would each lift it.
        f.reduce();
        f
    }
    fn neg(&self) -> Self {
        AFrac {
            num: self.num.iter().map(C::neg).collect(),
            den: self.den.clone(),
            tail: self.tail.clone(),
            scale: self.scale,
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

/// **A ℚ-algebra whatever `C` is** — including `C = i128`, which is not one.
///
/// Dividing by a positive integer here is multiplying `scale`, which is exact
/// and needs nothing from `C`. That is why the engines can run
/// over `AFrac<i128>` and still be handed to `s → p`, which asks for
/// [`QAlgebra`] because it divides by `z_μ`.
impl<C: Integral> QAlgebra for AFrac<C> {
    fn div_u128(&self, n: u128) -> Self {
        assert!(n != 0, "division of AFrac by zero");
        let mut out = self.clone();
        out.scale *= n;
        out.content_reduce();
        out
    }
}

/// The integer nth root of `k`, or `None` if `k` is not a perfect nth power.
fn nth_root(k: u32, n: u32) -> Option<u32> {
    if n == 1 {
        return Some(k);
    }
    let mut r = 0u32;
    while r.pow(n) < k {
        r += 1;
    }
    (r.pow(n) == k).then_some(r)
}

/// `(uα + v)` raised: `uα^n + v`, split into the linear factors it has and the
/// factor that has none.
///
/// **A primitive `uα + v` with `v ≥ 1` has at most one linear factor after
/// raising, and it is there only when `n` is odd and `u` and `v` are both
/// perfect nth powers.** A rational root `−s/t` in lowest terms of `uα^n + v`
/// satisfies `u·s^n = ±v·t^n`; with `gcd(u, v) = 1` and `gcd(s, t) = 1` that
/// forces `t^n = u` and `s^n = v`, and the sign forces `n` odd, since `u` and
/// `v` are both positive. There is at most one such root because `u x^n + v`
/// has exactly one real root for odd `n`, so the cofactor is root-free — which
/// is the invariant `tail` carries.
///
/// `v = 0` is the exception and is handled first: a primitive `(u, 0)` is
/// `(1, 0)`, the atom α, and `α^n` is `n` copies of it.
fn raise_atom<C: Integral>(u: u32, v: u32, n: u32) -> (Vec<(Atom, u32)>, Vec<C>) {
    if v == 0 {
        debug_assert_eq!(u, 1, "a primitive atom with no constant term is alpha");
        return (vec![((1, 0), n)], Vec::new());
    }
    let mut raised = vec![C::zero(); n as usize + 1];
    raised[0] = C::from_u128(v as u128);
    raised[n as usize] = C::from_u128(u as u128);
    if n % 2 == 0 {
        return (Vec::new(), raised);
    }
    let (Some(t), Some(s)) = (nth_root(u, n), nth_root(v, n)) else {
        return (Vec::new(), raised);
    };
    let cofactor = divide_by_linear(&raised, t, s)
        .expect("t*alpha + s divides u*alpha^n + v when t^n = u and s^n = v");
    (vec![((t, s), 1)], cofactor)
}

/// **The one operation that leaves the factored class**, and the reason
/// [`AFrac`] carries a `tail` at all.
///
/// `p_n` raises the variable, so over ℚ(α) it is α ↦ α^n. The numerator's
/// coefficients spread over every nth slot; a denominator atom `uα + v`
/// becomes `uα^n + v`, which is a linear form only when `n = 1`. What linear
/// factors it does have go back into the atoms through [`raise_atom`], and the
/// root-free rest joins the tail.
///
/// The coefficients are not pushed through a Frobenius of their own, and
/// nothing is missed by that: [`Integral`] is an integer ring, so its elements
/// are constants and α is the only variable there is to raise.
///
/// # Panics
///
/// Panics if `n == 0`, which is not a raising of variables: every exponent
/// would land on zero, an evaluation at `α = 1` rather than a Frobenius.
impl<C: Integral> Plethystic for AFrac<C> {
    fn frobenius(&self, n: u32) -> Self {
        assert!(
            n > 0,
            "the plethystic Frobenius needs n >= 1: p_0 does not raise variables"
        );
        if n == 1 || self.is_zero() {
            return self.clone();
        }
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
        let mut den: BTreeMap<Atom, u32> = BTreeMap::new();
        let mut tail = spread(&self.tail);
        for (&(u, v), &m) in &self.den {
            let (atoms, rest) = raise_atom::<C>(u, v, n);
            for (a, k) in atoms {
                *den.entry(a).or_insert(0) += k * m;
            }
            for _ in 0..m {
                tail = mul_tails(&tail, &rest);
            }
        }
        let mut out = AFrac {
            num: spread(&self.num),
            den,
            tail,
            scale: self.scale,
        };
        out.reduce();
        out
    }
}

impl<C: Integral> core::fmt::Display for AFrac<C> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if self.num.is_empty() {
            return f.write_str("0");
        }
        let mut first = true;
        f.write_str("(")?;
        for (k, c) in self.num.iter().enumerate() {
            if c.is_zero() {
                continue;
            }
            if !first {
                f.write_str(" + ")?;
            }
            first = false;
            match k {
                0 => write!(f, "{c:?}")?,
                1 => write!(f, "{c:?}a")?,
                _ => write!(f, "{c:?}a^{k}")?,
            }
        }
        f.write_str(")")?;
        if self.scale != 1 {
            write!(f, "/{}", self.scale)?;
        }
        for (&(u, v), &m) in &self.den {
            match (u, v) {
                (1, 0) => f.write_str("/(a)")?,
                (1, _) => write!(f, "/(a+{v})")?,
                (_, 0) => write!(f, "/({u}a)")?,
                _ => write!(f, "/({u}a+{v})")?,
            }
            if m > 1 {
                write!(f, "^{m}")?;
            }
        }
        if !self.tail.is_empty() {
            f.write_str("/(")?;
            let mut first = true;
            for (k, c) in self.tail.iter().enumerate() {
                if c.is_zero() {
                    continue;
                }
                if !first {
                    f.write_str("+")?;
                }
                first = false;
                match k {
                    0 => write!(f, "{c:?}")?,
                    1 => write!(f, "{c:?}a")?,
                    _ => write!(f, "{c:?}a^{k}")?,
                }
            }
            f.write_str(")")?;
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;

    fn r(n: i128) -> Rational {
        Rational::from_int(n)
    }

    /// The plethystic Frobenius is the one operation that leaves the factored
    /// class, and the tail is where what leaves it goes.
    ///
    /// `1/(α+1)` at `n = 2` must become `1/(α²+1)` — irreducible, so no atom
    /// survives — and **not** `1/(α+1)²`, which is the reading that confuses
    /// raising the variable with raising the form.
    #[test]
    fn frobenius_raises_alpha_and_the_tail_catches_what_leaves() {
        type F = AFrac<i128>;
        let a = F::inv_linear(1, 1); // 1/(α+1)
        let two = a.frobenius(2);
        assert_eq!(two.parts().1.count(), 0, "alpha^2 + 1 has no linear factor");
        assert_eq!(two.tail(), &[1, 0, 1], "1/(alpha^2 + 1)");
        assert_eq!(a.frobenius(1), a, "n = 1 is the identity");

        // α³ + 1 = (α + 1)(α² − α + 1): the linear factor goes back to the
        // atoms and only the root-free cofactor stays in the tail. Leaving it
        // whole would give one element two spellings.
        let three = a.frobenius(3);
        let atoms: Vec<_> = three.parts().1.map(|(a, &m)| (*a, m)).collect();
        assert_eq!(atoms, vec![((1, 1), 1)]);
        assert_eq!(three.tail(), &[1, -1, 1]);

        // α itself is the one atom whose raising is all linear.
        let over_alpha = F::inv_linear(1, 0);
        let cubed = over_alpha.frobenius(3);
        let atoms: Vec<_> = cubed.parts().1.map(|(a, &m)| (*a, m)).collect();
        assert_eq!(atoms, vec![((1, 0), 3)], "1/alpha^3");
        assert!(cubed.tail().is_empty());
    }

    /// A Frobenius must be a ring homomorphism, and the tail is where that is
    /// easiest to get wrong: two values with different tails have to reach a
    /// common denominator before they can be added.
    #[test]
    fn frobenius_is_a_ring_homomorphism() {
        type F = AFrac<i128>;
        let a = F::inv_linear(1, 1);
        let b = F::from_coeffs(vec![0, 1]).div_linear(1, 2); // α/(α+2)
        for n in 1..5 {
            assert_eq!(
                a.mul(&b).frobenius(n),
                a.frobenius(n).mul(&b.frobenius(n)),
                "multiplicative at n = {n}"
            );
            let mut sum = a.clone();
            sum.add_assign(&b);
            sum.reduce();
            let mut raised = a.frobenius(n);
            raised.add_assign(&b.frobenius(n));
            raised.reduce();
            assert_eq!(sum.frobenius(n), raised, "additive at n = {n}");
        }
    }

    /// Values with tails must still compare, add and cancel as elements, not
    /// as representations — the property the whole type rests on.
    #[test]
    fn a_tail_cancels_against_the_numerator() {
        type F = AFrac<i128>;
        let tailed = F::inv_linear(1, 1).frobenius(2); // 1/(α²+1)
        let mut whole = tailed.mul(&F::from_coeffs(vec![1, 0, 1])); // ·(α²+1)
        whole.reduce();
        assert_eq!(whole, <F as Ring>::one(), "the tail divides out");
        assert!(whole.tail().is_empty());

        // 1/(α²+1) + 1/(α²+1) is 2/(α²+1), over one tail rather than its square.
        let mut doubled = tailed.clone();
        doubled.add_assign(&tailed);
        doubled.reduce();
        assert_eq!(doubled.tail(), &[1, 0, 1]);
        assert_eq!(doubled, tailed.scale_int(2));
        assert_eq!(
            tailed.eval_i128(&r(1)),
            Some(Rational::new(1, 2)),
            "at alpha = 1"
        );
    }

    /// The atom normalization, including the two edge shapes that actually
    /// occur: a constant hook (`h_low` with no arm) and `α` itself (`h_up` with
    /// no leg).
    #[test]
    fn atoms_normalize_to_primitive_forms() {
        assert_eq!(split(0, 5), (5, None), "0·α + 5 is the constant 5");
        assert_eq!(split(3, 0), (3, Some((1, 0))), "3α = 3 · α");
        assert_eq!(split(2, 4), (2, Some((1, 2))), "2α + 4 = 2 · (α + 2)");
        assert_eq!(split(3, 2), (1, Some((3, 2))), "already primitive");
    }

    /// ⚠️ The spec's own example of why primitivity matters: κ = (2,2),
    /// λ = (1,1,1,1) gives `E(κ) − E(λ) = 2α + 4`, and dividing `ℤ[α]` by that
    /// leaves `ℤ[α]` unless the 2 is split off first.
    #[test]
    fn dividing_by_a_non_primitive_form_stays_integral() {
        // (2α + 4) divides (2α + 4) — trivially, but the *representation* must
        // put the 2 in `scale` and α+2 in the atoms, or the numerator leaves ℤ.
        let f: AFrac<i128> = AFrac::linear(2, 4);
        let q = f.div_linear(2, 4);
        assert_eq!(q, <AFrac<i128> as Ring>::one(), "{q}");
        // And a genuinely fractional value keeps its 2 in the scalar rather
        // than putting 1/2 into an integer numerator.
        let half: AFrac<i128> = <AFrac<i128> as Ring>::one().div_linear(2, 4);
        let (num, _, scale) = half.parts();
        assert_eq!(scale, 2, "the content lives in the scalar: {half}");
        assert_eq!(num, &[1i128]);
    }

    #[test]
    fn exact_division_by_a_linear_form() {
        // (α + 1)(α + 2) = α² + 3α + 2
        let p = vec![2i128, 3, 1];
        let q = divide_by_linear(&p, 1, 1).expect("(α+1) divides it");
        assert_eq!(q, vec![2i128, 1], "α + 2");
        assert!(divide_by_linear(&p, 1, 3).is_none(), "(α+3) does not");
        // The Gauss's-lemma early exit: (2α+1) is primitive and 2 ∤ 1.
        assert!(divide_by_linear(&p, 2, 1).is_none());
    }

    /// Division by `α` itself — the `(1, 0)` atom, which `h_up` produces at
    /// every cell in the last row.
    #[test]
    fn division_by_alpha() {
        let p = vec![0i128, 3, 1]; // α² + 3α
        assert_eq!(divide_by_linear(&p, 1, 0), Some(vec![3i128, 1]));
        let q = vec![1i128, 3, 1]; // α² + 3α + 1
        assert_eq!(divide_by_linear(&q, 1, 0), None);
    }

    #[test]
    fn round_trip_through_every_atom() {
        for p in [vec![1i128], vec![2, 3, 1], vec![0, 0, 5], vec![-7, 4]] {
            for (u, v) in [(1u32, 0u32), (1, 1), (1, 5), (2, 1), (3, 2)] {
                let prod = mul_linear(&p, u, v);
                let back = divide_by_linear(&prod, u, v)
                    .unwrap_or_else(|| panic!("({u}α+{v}) must divide its own multiple"));
                assert_eq!(back, p, "round trip through {u}α+{v}");
            }
        }
    }

    #[test]
    fn ring_axioms_hold() {
        type F = AFrac<Rational>;
        let a = F::inv_linear(1, 1); // 1/(α+1)
        let b = F::inv_linear(1, 2); // 1/(α+2)
        let one = <F as Ring>::one();
        assert_eq!(a.mul(&one), a);
        assert_eq!(a.mul(&b), b.mul(&a));
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

    /// Addition must use the lcm, and for *this* family the lcm is exact —
    /// adding `1/(α+1)` to itself must not square the denominator.
    #[test]
    fn addition_uses_the_exact_lcm() {
        type F = AFrac<Rational>;
        let mut s = F::inv_linear(1, 1);
        s.add_assign(&F::inv_linear(1, 1));
        let (num, den, scale) = s.parts();
        assert_eq!(num, &[r(2)]);
        assert_eq!(scale, 1);
        assert_eq!(den.collect::<Vec<_>>(), vec![(&(1, 1), &1)], "{s}");
    }

    /// 1/(α+1) + 1/(α+2) = (2α+3)/((α+1)(α+2)) — two coprime atoms, so the lcm
    /// is the product and nothing cancels.
    #[test]
    fn addition_over_two_atoms() {
        type F = AFrac<Rational>;
        let mut s = F::inv_linear(1, 1);
        s.add_assign(&F::inv_linear(1, 2));
        assert_eq!(s.eval(&r(3)), Some(Rational::new(9, 20)), "{s}");
        let (num, _, _) = s.parts();
        assert_eq!(num, &[r(3), r(2)]);
    }

    #[test]
    fn reduction_cancels_and_into_poly_answers_about_the_element() {
        type F = AFrac<i128>;
        // (α² + 3α + 2)/(α + 1) = α + 2
        let f = F::from_coeffs(vec![2, 3, 1]).div_linear(1, 1);
        assert_eq!(f.clone().into_poly(), Some(vec![2i128, 1]), "{f}");
        // 1/(α+1) is genuinely not a polynomial
        assert_eq!(F::inv_linear(1, 1).into_poly(), None);
        // and over ℤ a surviving 1/2 is not one either
        assert_eq!(F::from_coeffs(vec![1]).div_linear(0, 2).into_poly(), None);
        // over ℚ it is
        let g: AFrac<Rational> = AFrac::from_coeffs(vec![r(1)]).div_linear(0, 2);
        assert_eq!(g.into_poly(), Some(vec![Rational::new(1, 2)]));
    }

    /// Cross-multiplied equality: the same element stored with the scalar in
    /// two different places.
    #[test]
    fn equality_survives_an_unnormalized_content() {
        type F = AFrac<Rational>;
        let a = F::from_coeffs(vec![r(1)]).div_linear(0, 2); // absorbed: [1/2]
        let b = F::from_coeffs(vec![r(3)]).div_linear(0, 6); // absorbed: [1/2]
        assert_eq!(a, b, "{a} vs {b}");
        // over ℤ neither absorbs, and they are still equal
        let c: AFrac<i128> = AFrac::from_coeffs(vec![1]).div_linear(0, 2);
        let d: AFrac<i128> = AFrac::from_coeffs(vec![3]).div_linear(0, 6);
        assert_eq!(c, d, "{c} vs {d}");
    }

    #[test]
    fn factors_multiply_and_cancel_before_expanding() {
        type F = AFrac<Rational>;
        let mut f: Linears = Linears::new();
        f.insert((1, 2), 1); // (α + 2)
        f.insert((0, 3), 1); // the constant 3
        f.insert((2, 1), -1); // / (2α + 1)
        let v = F::from_factors(&f);
        // at α = 1: 3·3/3 = 3
        assert_eq!(v.eval(&r(1)), Some(r(3)), "{v}");
        // multiplying by the inverse multiset gives 1
        let inv: Linears = f.iter().map(|(&k, &m)| (k, -m)).collect();
        let one = v.mul_factors(&inv);
        assert_eq!(one, <F as Ring>::one(), "{one}");
    }

    #[test]
    fn division_by_integers_is_exact_over_a_non_qalgebra() {
        // The point of the `QAlgebra for AFrac<C>` impl: `i128` is not a
        // ℚ-algebra, but `AFrac<i128>` is.
        type F = AFrac<i128>;
        let f = <F as Ring>::from_i64(6).div_u128(4);
        assert_eq!(f.eval_rational(), (vec![3i128], 2));
    }

    impl AFrac<i128> {
        /// Test helper: the reduced numerator and scalar of a denominator-free
        /// value.
        fn eval_rational(mut self) -> (Vec<i128>, u128) {
            self.reduce();
            assert!(self.den.is_empty());
            (self.num, self.scale)
        }
    }

    /// `α ↦ 1/α` must be an exact involution that stays in the family, and it
    /// must agree with evaluating at the reciprocal.
    #[test]
    fn alpha_inversion_is_an_involution() {
        type F = AFrac<Rational>;
        let cases = [
            F::from_coeffs(vec![r(2), r(3), r(1)]),      // α² + 3α + 2
            F::inv_linear(1, 1),                         // 1/(α+1)
            F::inv_linear(1, 0),                         // 1/α — the atom that vanishes
            F::linear(3, 2).div_linear(2, 5),            // (3α+2)/(2α+5)
            F::from_coeffs(vec![r(7)]).div_linear(0, 3), // 7/3, no α at all
        ];
        for f in &cases {
            let back = f.invert_alpha().invert_alpha();
            assert_eq!(
                &back,
                f,
                "involution: {f} -> {} -> {back}",
                f.invert_alpha()
            );
            // and it really is the substitution: value at 1/x, evaluated at x
            for &x in &[2i128, 3, 5] {
                let want = f.eval(&Rational::new(1, x));
                let got = f.invert_alpha().eval(&r(x));
                assert_eq!(got, want, "{f} at α = 1/{x}");
            }
        }
    }

    /// `1/α ↦ α` is the case where an atom disappears into the scalar, which is
    /// the reason the swap cannot be done pairwise.
    #[test]
    fn inverting_alpha_can_remove_an_atom() {
        type F = AFrac<Rational>;
        let f = F::inv_linear(1, 0); // 1/α
        let g = f.invert_alpha(); // α
        assert_eq!(g, F::linear(1, 0), "{g}");
        let (_, den, _) = g.parts();
        assert_eq!(den.count(), 0, "the α atom must be gone, not merely equal");
    }

    #[test]
    fn evaluation_and_its_poles() {
        type F = AFrac<Rational>;
        let f = F::inv_linear(1, 1);
        assert_eq!(f.eval(&r(3)), Some(Rational::new(1, 4)));
        assert_eq!(f.eval(&r(-1)), None, "α = −1 is a pole");
        // a pole that cancels is not a pole
        let g = F::linear(1, 1).div_linear(1, 1);
        assert_eq!(g.eval(&r(-1)), Some(r(1)));
    }
}
