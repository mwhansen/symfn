//! ℚ(α), restricted to denominators that factor into linear forms `uα + v`.
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
//!    [`Ratio::add_mul`](crate::deltaop::Ratio) settles for a common multiple
//!    because `q² − t²` factors.
//! 3. **A failed cancellation is detected on its first step.** Synthetic
//!     division by `uα + v` starts at the top coefficient and needs `u` to
//!     divide it; by Gauss's lemma that is *necessary* for divisibility in
//!     `ℚ[α]`, so the usual failure exits immediately. The (q,t) engine has the
//!     opposite problem — `divide_exact` runs its failures to completion — and
//!     needs a bespoke necessary-condition pre-pass
//!     (`deltaop::diff_may_divide`).
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
//!   value = num(α) / (scale · ∏ (uα + v)^m)
//! ```
//!
//! with `num` a **dense** `Vec<C>` (a `QtPoly` with a dead `t` would be a
//! sparse two-variable key for a dense univariate object), the atoms a
//! `BTreeMap`, and `scale` a positive integer. The integer denominator is what
//! lets `C = i128` stay integral: a value like `1/(2(α+2))` has no home in
//! `ℤ[α]` otherwise.
//!
//! `C` must implement [`Ring::div_exact`] faithfully — ℤ-like (`i128`,
//! `BigInt`) or a field ([`Rational`](crate::coeff::Rational), `BigRational`).
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

use crate::coeff::{Field, QAlgebra, Ring};

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
pub struct AFrac<C: Ring> {
    /// Dense in α: `num[k]` multiplies `α^k`. No trailing zeros; empty is 0.
    num: Vec<C>,
    /// The primitive atoms with their multiplicities.
    den: BTreeMap<Atom, u32>,
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

fn gcd128(mut a: u128, mut b: u128) -> u128 {
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

impl<C: Ring> AFrac<C> {
    /// A polynomial in α, given densely by its coefficients.
    pub fn from_coeffs(mut num: Vec<C>) -> Self {
        trim(&mut num);
        AFrac {
            num,
            den: BTreeMap::new(),
            scale: 1,
        }
    }

    /// The polynomial `uα + v`.
    pub fn linear(u: u32, v: u32) -> Self {
        Self::from_coeffs(vec![C::from_u128(v as u128), C::from_u128(u as u128)])
    }

    /// `1 / (uα + v)`. Panics on `0α + 0`.
    pub fn inv_linear(u: u32, v: u32) -> Self {
        let (content, atom) = split(u, v);
        let mut den = BTreeMap::new();
        if let Some(a) = atom {
            den.insert(a, 1);
        }
        AFrac {
            num: vec![C::one()],
            den,
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
    pub fn from_factors(factors: &Linears) -> Self {
        <Self as Ring>::one().mul_factors(factors)
    }

    /// Multiply by `∏ (uα + v)^m`, negative `m` meaning a denominator factor.
    ///
    /// The `Q` and `J` normalizers are exactly this shape — products of hooks —
    /// and reaching them through [`Ring::mul`] would expand `H_λ` into a
    /// polynomial first. Here every step is one linear multiply or one exact
    /// linear division.
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

    /// Multiply the *value* by the positive integer `k`, cancelling it against
    /// the scalar denominator first.
    fn scale_content(&mut self, k: u128) {
        if k == 1 {
            return;
        }
        let g = gcd128(self.scale, k);
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
            self.scale = 1;
            return;
        }
        let atoms: Vec<Atom> = self.den.keys().copied().collect();
        for a in atoms {
            self.reduce_at(a);
        }
        self.content_reduce();
    }

    /// Cancel one atom as far as it goes.
    ///
    /// The whole loop allocates nothing: [`divide_in_place`] tests without a
    /// buffer and rewrites without one. This is the hot path of the entire
    /// module — see [`divides_by_linear`] for what it used to cost.
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
    /// prime factors are bounded by roughly `n²` and [`coarse_factors`] finds
    /// all of them by trial division in a few dozen steps.
    ///
    /// A large prime residue is tried once as a lump and then abandoned rather
    /// than factored. That can leave a cancellation on the table; it cannot
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

    /// The numerator coefficients, the denominator's atoms, and the scalar.
    pub fn parts(&self) -> (&[C], impl Iterator<Item = (&Atom, &u32)>, u128) {
        (&self.num, self.den.iter(), self.scale)
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
        if !self.den.is_empty() {
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
        self.den.is_empty().then_some((self.num, self.scale))
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
    pub fn invert_alpha(&self) -> Self {
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
        let mut out = AFrac { num, den, scale: 1 };
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
    fn lift(&self, target: &BTreeMap<Atom, u32>, scale: u128) -> Vec<C> {
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

impl<C: Field> AFrac<C> {
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
        let mut acc = C::zero();
        for c in self.num.iter().rev() {
            acc = acc.mul(alpha);
            acc.add_assign(c);
        }
        Some(acc.div(&d))
    }
}

/// Cross-multiplied — see the module docs. The atoms are canonical; the integer
/// content is not always, so structural comparison would call equal things
/// unequal for `C = Rational`.
impl<C: Ring> PartialEq for AFrac<C> {
    fn eq(&self, other: &Self) -> bool {
        if self.num.is_empty() || other.num.is_empty() {
            return self.num.is_empty() && other.num.is_empty();
        }
        if self.den == other.den && self.scale == other.scale {
            return self.num == other.num;
        }
        let mut lcm = self.den.clone();
        for (k, &m) in &other.den {
            let e = lcm.entry(*k).or_insert(0);
            *e = (*e).max(m);
        }
        let g = gcd128(self.scale, other.scale);
        let s = self.scale / g * other.scale;
        self.lift(&lcm, s) == other.lift(&lcm, s)
    }
}

impl<C: Ring> Eq for AFrac<C> {}

impl<C: Ring> Ring for AFrac<C> {
    fn zero() -> Self {
        AFrac {
            num: Vec::new(),
            den: BTreeMap::new(),
            scale: 1,
        }
    }
    fn one() -> Self {
        AFrac {
            num: vec![C::one()],
            den: BTreeMap::new(),
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
        if needs_atoms || needs_scale {
            let mut lcm = self.den.clone();
            for (k, &m) in &other.den {
                let e = lcm.entry(*k).or_insert(0);
                *e = (*e).max(m);
            }
            let g = gcd128(self.scale, other.scale);
            let s = self.scale / g * other.scale;
            self.num = self.lift(&lcm, s);
            self.den = lcm;
            self.scale = s;
        }
        let lifted = other.lift(&self.den, self.scale);
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
impl<C: Ring> QAlgebra for AFrac<C> {
    fn div_u128(&self, n: u128) -> Self {
        assert!(n != 0, "division of AFrac by zero");
        let mut out = self.clone();
        out.scale *= n;
        out.content_reduce();
        out
    }
}

impl<C: Ring> core::fmt::Display for AFrac<C> {
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
