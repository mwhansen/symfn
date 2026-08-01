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

// Exponent arithmetic on `1 − qᵃtᵇ` atoms, bounded by the degree.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

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

    /// `1 / (1 − qᵃtᵇ)`.
    ///
    /// # Panics
    ///
    /// On `(0, 0)`, which would be `1/0`.
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
    ///
    /// **Not reduced**, on the same policy as [`Ring::add_assign`]: this is
    /// called once per tableau and reduction is a per-*coefficient* operation.
    /// Calling [`reduce`](Self::reduce) here cost 4× — 72% of its trial
    /// divisions fail — and changed nothing, because the one `reduce` at the end
    /// of a coefficient reaches the same form. The Macdonald dumps are
    /// byte-identical with it and without it.
    /// # Panics
    ///
    /// If any key is `(0, 0)`. That factor is `1 − q⁰t⁰ = 0`: with a positive
    /// exponent the whole product is zero, and with a negative one it is a zero
    /// *denominator* factor, which nothing downstream would catch.
    pub fn from_factors(factors: &BTreeMap<(u32, u32), i32>) -> Self {
        let mut num = <QtPoly<C> as Ring>::one();
        let mut den = BTreeMap::new();
        for (&(a, b), &m) in factors {
            assert!(a > 0 || b > 0, "1 - q^0 t^0 is zero");
            match m.cmp(&0) {
                core::cmp::Ordering::Greater => {
                    for _ in 0..m {
                        num = num.mul_binomial(a, b);
                    }
                }
                core::cmp::Ordering::Less => {
                    den.insert((a, b), (-m) as u32);
                }
                core::cmp::Ordering::Equal => {}
            }
        }
        Frac { num, den }
    }

    /// Multiply by `∏ (1 − qᵃtᵇ)^{m}`, negative `m` meaning a denominator
    /// factor — the same encoding [`from_factors`](Self::from_factors) reads.
    ///
    /// The scalars `Q` and `J` apply to `P` are exactly of this shape, and
    /// reaching them through [`Ring::mul`] would expand `b_λ` or `c_λ` into a
    /// polynomial first and then run the general product against it. Applying
    /// the factors one at a time keeps every multiplication a
    /// [`QtPoly::mul_binomial`].
    /// # Panics
    ///
    /// If any key is `(0, 0)`; see [`from_factors`](Self::from_factors).
    pub fn mul_factors(&self, factors: &BTreeMap<(u32, u32), i32>) -> Self {
        let mut num = self.num.clone();
        let mut den = self.den.clone();
        for (&(a, b), &m) in factors {
            assert!(a > 0 || b > 0, "1 - q^0 t^0 is zero");
            match m.cmp(&0) {
                core::cmp::Ordering::Greater => {
                    for _ in 0..m {
                        num = num.mul_binomial(a, b);
                    }
                }
                core::cmp::Ordering::Less => {
                    *den.entry((a, b)).or_insert(0) += (-m) as u32;
                }
                core::cmp::Ordering::Equal => {}
            }
        }
        Frac { num, den }
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
                d = d.mul_binomial(a, b);
            }
        }
        d
    }

    /// The numerator, if the denominator cancels away entirely — `None` if this
    /// element genuinely is not a polynomial.
    ///
    /// Reduces first, so it answers about the *element* and not about the
    /// representation, which is not canonical. Callers that know on
    /// mathematical grounds that the result must be a polynomial — Macdonald's
    /// `J`, the (q,t)-Kostka polynomials — should unwrap it and let a `None`
    /// be the loud failure it is.
    pub fn into_poly(mut self) -> Option<QtPoly<C>> {
        self.reduce();
        self.den.is_empty().then_some(self.num)
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
                num = num.mul_binomial(a, b);
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
/// From `N = Q·(1 − qᵃtᵇ)`, matching coefficients gives
///
/// ```text
///   Q[k] = N[k] + Q[k − δ]        δ = (a, b)
/// ```
///
/// so `Q` along a **chain** `k, k+δ, k+2δ, …` is just the running sum of `N`
/// along it, and the chains are independent. Divisibility is the statement that
/// every chain sums to zero — nothing else is needed, and it falls out of the
/// same walk that builds the quotient.
///
/// Terms arrive lex-ascending and `+δ` is monotone for that order, so a chain's
/// members are encountered in order and the first *unconsumed* term reached is
/// always the start of its chain: its predecessor `k − δ` is lex-smaller, so if
/// it were a term of `N` its own walk would have consumed this one.
///
/// This replaced a `BTreeMap` remainder that popped the least key and inserted a
/// larger one per step. That was 782 samples of a 3300-sample profile with
/// another ~500 in the B-tree itself, and **72% of the calls fail** — `reduce`
/// trial-divides by every denominator factor and only 28% divide — so the
/// failures were most of the cost. Here a failure is detected by a chain sum
/// that will not vanish, at the same price as the success.
///
/// The two exits both rest on `deg(Q) ≤ deg(N) − (a + b)`, for the total degree:
/// if `M` is a maximal-degree term of `Q` then `Q[M + δ] = 0`, so
/// `N[M + δ] = −Q[M] ≠ 0` and `M + δ` is a term of `N`.
///
/// * A nonzero running sum at `p` with `deg(p) + a + b > bound` cannot be a
///   quotient term, so the division is inexact.
/// * A zero running sum at `p` with `deg(p) > bound` ends the chain: every term
///   of `N` has degree at most `bound`, so none can remain further along it.
pub(crate) fn divide_by_factor<C: Ring>(n: &QtPoly<C>, a: u32, b: u32) -> Option<QtPoly<C>> {
    let terms = n.raw();
    // Emptiness and the degree bound are one question: no terms, no bound, and
    // the zero polynomial divides.
    let Some(bound) = terms.iter().map(|(k, _)| k.0 + k.1).max() else {
        return Some(QtPoly::zero());
    };
    let mut consumed = vec![false; terms.len()];
    let mut out: Vec<((u32, u32), C)> = Vec::with_capacity(terms.len());

    for i in 0..terms.len() {
        if consumed[i] {
            continue;
        }
        let mut sum = C::zero();
        let mut p = terms[i].0;
        // Chain members sit at increasing indices, so the search window only
        // ever shrinks from the left.
        let mut lo = i;
        loop {
            match terms[lo..].binary_search_by_key(&p, |e| e.0) {
                Ok(off) => {
                    let j = lo + off;
                    sum.add_assign(&terms[j].1);
                    consumed[j] = true;
                    lo = j + 1;
                }
                Err(off) => lo += off,
            }
            if sum.is_zero() {
                if p.0 + p.1 > bound {
                    break;
                }
            } else {
                if p.0 + p.1 + a + b > bound {
                    return None;
                }
                out.push((p, sum.clone()));
            }
            p = (p.0 + a, p.1 + b);
        }
    }
    // Chains interleave, so the pieces come out sorted individually but not
    // together. Only the 28% of calls that divide ever reach this.
    //
    // `sort_unstable` and not `sort`: [`divide_by_diff`] below records the
    // measurement, but the short version is that these runs are many and short,
    // so the stable sort's run merging has nothing to exploit and buys 0%.
    out.sort_unstable_by_key(|x| x.0);
    Some(QtPoly::from_sorted(out))
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
        //
        // `self` is the accumulator in every caller that matters, so it is the
        // large side. Its denominator stops growing after the first few terms,
        // and from then on the lcm *is* its denominator: check before rewriting
        // it, or every addition clones and re-walks the whole numerator to
        // multiply it by nothing.
        if other
            .den
            .iter()
            .any(|(k, &m)| self.den.get(k).copied().unwrap_or(0) < m)
        {
            let mut lcm = self.den.clone();
            for (k, &m) in &other.den {
                let e = lcm.entry(*k).or_insert(0);
                *e = (*e).max(m);
            }
            self.num = self.lift(&lcm);
            self.den = lcm;
        }
        let lifted = other.lift(&self.den);
        self.num.add_assign(&lifted);
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
        // Reduced, unlike `add_assign` and `from_factors`. That looks like an
        // inconsistency and was removed once on the strength of a profile —
        // `divide_by_factor` reached from here was 38% of a (q,t)-Kostka run,
        // and the multiplier in `PowerSum::to_schur` is a symmetric-group
        // *character*, a constant, which over ℚ is a unit and cannot make a
        // binomial newly divide anything. Every trial division that ran was
        // doomed before it started.
        //
        // Removing it made that run 1.7× **slower**. The reasoning about this
        // product was right and the conclusion was wrong: the coefficient
        // arriving here has never been reduced by anything else, so this call
        // was where a denominator first got cut down — and it is cut down once,
        // before `add_assign` uses it p(n) times and lifts every other term to
        // the lcm. Cheap here, quadratic to skip.
        //
        // So the policy is not "never reduce eagerly", it is "reduce once,
        // before the value is used many times". See `qtkostka::invert_s_basis`,
        // which reduces for the same reason and makes this one redundant.
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
    /// 1 — which is how `b_λ(s)` behaves for a cell outside the diagram.
    ///
    /// # Panics
    ///
    /// If `(a₂, b₂) == (0, 0)`. The same pair is a *legal* numerator and an
    /// illegal denominator, which is why only one side is checked.
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

// --- the `qᵃ − tᵇ` atom family -----------------------------------------------
//
// Lives here rather than in `deltaop`, which is where it was written, because
// `bh`'s Bergeron–Haiman recursion divides by the same atoms and was reaching
// for the generic `QtPoly::divide_exact` instead — 35% of an `htilde_table`
// profile. Two callers, one home, next to `divide_by_factor` which is the same
// job for the `1 − qᵃtᵇ` family.

/// A one-pass **necessary** condition for `qᵃ − tᵇ` to divide `n`.
///
/// This is the single most important line in the module for speed, and it is
/// there for the reason [`Frac`](crate::Frac)'s own division notes give: in
/// [`Ratio::reduce`] most trial divisions *fail*, so the cost of the failures is
/// the cost of the reduction.
///
/// [`divide_by_factor`](crate::frac) detects a failure early because a factor
/// `1 − qᵃtᵇ` has its leading term at `(0,0)` and the chain sums run out
/// quickly. [`QtPoly::divide_exact`](crate::qt::QtPoly) does not: the leading
/// term of `qᵃ − tᵇ` is `qᵃ` (lex, and `a ≥ 1` for every `Diff` atom), so the
/// "leading monomial is not a multiple" exit never fires on the `t` exponent and
/// a doomed division still runs the **whole** elimination — building a
/// `BTreeMap` of the numerator and eliminating every term — only to find a
/// nonempty remainder at the end. Measured before this filter existed: `∇e_11`
/// took 44.4s, essentially all of it here.
///
/// The test is exact and needs no arithmetic beyond addition. Write
/// `d = gcd(a,b)`, `a' = a/d`, `b' = b/d`; the substitution `q ↦ s^{b'}`,
/// `t ↦ s^{a'}` sends `q^{a'} − t^{b'}` to zero, so it annihilates every
/// multiple of it. Since `(qᵃ − tᵇ) | n` implies `(q^{a'} − t^{b'}) | n`, a
/// nonzero image proves non-divisibility. The image is one univariate
/// polynomial in `s`, accumulated by bucketing each term of `n` at degree
/// `b'·x + a'·y`.
///
/// It is deliberately **not** a sufficient condition — when `gcd(a,b) > 1` the
/// substitution only kills one irreducible factor — so a `true` still runs the
/// real division. That keeps it a filter and not a second answer.
pub(crate) fn diff_may_divide<C: Ring>(n: &QtPoly<C>, a: u32, b: u32) -> bool {
    let mut d = (a, b);
    while d.1 != 0 {
        d = (d.1, d.0 % d.1);
    }
    let (ap, bp) = (a / d.0, b / d.0);
    let mut acc: Vec<C> = Vec::new();
    for (&(x, y), c) in n.terms() {
        let deg = (x * bp + y * ap) as usize;
        if acc.len() <= deg {
            acc.resize(deg + 1, C::zero());
        }
        acc[deg].add_assign(c);
    }
    acc.iter().all(C::is_zero)
}

/// Exact division by `qᵃ − tᵇ` (both exponents ≥ 1), or `None`.
///
/// The counterpart of [`divide_by_factor`](crate::frac) for the other atom
/// family, and it exists for the reason that one records: routing this through
/// [`QtPoly::divide_exact`](crate::qt::QtPoly) keeps the remainder in a
/// `BTreeMap` and pays a node rebalance per elimination step. Sampled at degree
/// 12 (`sample`, the workflow `Cargo.toml` documents), that put **83% of the
/// whole profile** inside `Atom::divide`, nearly all of it in B-tree
/// `remove_kv_tracking` / `bulk_steal_left` / `memmove`. `frac.rs` found the
/// same thing for Macdonald `P` and fixed it the same way.
///
/// ## The chain
///
/// Eliminating leading terms against `qᵃ − tᵇ` is a *flow*. Its lex-leading
/// monomial is `qᵃ` (any `Diff` atom has `a ≥ 1`), so the step at the current
/// maximum key `(x,y)` emits the quotient term `(x−a, y)` and moves that
/// coefficient — added, since the other monomial of the divisor is `−tᵇ` — to
/// `(x−a, y+b)`. Writing `σ(x,y) = (x−a, y+b)`, the whole division is
///
/// ```text
///   running := 0
///   for p along a σ-chain:   running += N[p];  Q[(p.0−a, p.1)] := running
/// ```
///
/// so `Q` is a **running sum of `N` along the chain**, exactly as
/// `divide_by_factor` is a running sum along `k + δ`. The chains are
/// independent, and `σ` strictly decreases the lexicographic key, so walking
/// `N`'s terms in **descending** order means the first unconsumed term reached
/// is always the head of its chain: its predecessor `(x+a, y−b)` is lex-greater,
/// so had it been a term of `N` its own walk would have consumed this one.
///
/// ## Both exits
///
/// `deg_t(N) = deg_t(Q) + b`, because `q^a·Q` and `t^b·Q` cannot cancel at the
/// top `t`-degree, so no quotient term may have `t`-exponent above `max_t(N)`:
///
/// * a nonzero running sum at `p` with `p.0 < a` or `p.1 > max_t` cannot be a
///   quotient term, so the division is inexact;
/// * a zero running sum with `p.1 > max_t` ends the chain, since no term of `N`
///   can lie further along it.
pub(crate) fn divide_by_diff<C: Ring>(n: &QtPoly<C>, a: u32, b: u32) -> Option<QtPoly<C>> {
    debug_assert!(a > 0 && b > 0, "divide_by_diff is for genuine Diff atoms");
    // Borrowed, not collected: this runs on every trial division, and cloning
    // the term list first is the same mistake in miniature as lifting before
    // multiplying above.
    let terms = n.raw();
    if terms.is_empty() {
        return Some(QtPoly::zero());
    }
    let max_t = terms
        .iter()
        .map(|(k, _)| k.1)
        .max()
        .expect("the empty case returned three lines above, so `terms` has a maximum");
    let mut consumed = vec![false; terms.len()];
    let mut out: Vec<((u32, u32), C)> = Vec::with_capacity(terms.len());

    for i in (0..terms.len()).rev() {
        if consumed[i] {
            continue;
        }
        let mut running = C::zero();
        let mut p = terms[i].0;
        // Chain positions descend lexicographically while `terms` ascends, so
        // the search window only ever shrinks from the right — the mirror of the
        // shrink-from-the-left `divide_by_factor` uses.
        let mut hi = terms.len();
        loop {
            match terms[..hi].binary_search_by_key(&p, |e| e.0) {
                Ok(j) => {
                    if !consumed[j] {
                        running.add_assign(&terms[j].1);
                        consumed[j] = true;
                    }
                    hi = j;
                }
                Err(j) => hi = j,
            }
            if !running.is_zero() {
                if p.0 < a || p.1 > max_t {
                    return None;
                }
                out.push(((p.0 - a, p.1), running.clone()));
            } else if p.1 > max_t || hi == 0 {
                // Nothing can change again: either every remaining position on
                // the chain is past the degree bound, or the window is empty.
                //
                // ⚠️ Breaking on `running == 0` *alone* is also correct — the
                // flow onward is zero, so the terms ahead are exactly a fresh
                // chain and the outer scan reaches them later — and it is
                // **slower**, measured 12.3s → 14.6s on `∇e_12`. Each restart
                // resets the search window to the full term list, and losing the
                // shrink costs more than the walking it saves. Recorded because
                // it is the obvious optimisation and it loses.
                break;
            }
            if p.0 < a {
                break;
            }
            p = (p.0 - a, p.1 + b);
        }
    }
    // Chains interleave, so the pieces are each sorted but not jointly: `out` is
    // a concatenation of descending runs, one per chain.
    //
    // ⚠️ **This sort is 15.6% of a `∇e_11` profile** — more than any arithmetic
    // in the module — and the obvious fix does not work. A concatenation of
    // ordered runs is exactly the shape `slice::sort`'s merge detects and
    // `sort_unstable`'s pattern-defeating quicksort does not, so the stable sort
    // ought to win outright. Measured: `∇e_12` went 7.95s → 7.85s, i.e. nothing.
    // The runs are **many and short** rather than few and long, so there is no
    // run structure to exploit and this is honest sorting work. Recorded because
    // it is the obvious optimisation and it loses — the same posture the
    // `running == 0` note above takes.
    //
    // Untried, and the reason it stays untried: a counting sort on the first
    // coordinate would work (`x` strictly decreases along a chain, so each chain
    // contributes at most one term per `x` and the `x`-groups are small), but it
    // needs a scratch buffer sized by the q-degree, and the ceiling on the whole
    // idea is ~1.15× on `∇e_n`.
    out.sort_unstable_by_key(|x| x.0);
    Some(QtPoly::from_sorted(out))
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

    /// The chain formulation has to handle a quotient term at a key the
    /// *numerator* does not have: `(1 − q³)/(1 − q) = 1 + q + q²` visits `q` and
    /// `q²`, neither of which is a term of `1 − q³`. A walk that only stepped
    /// between existing terms would skip them.
    #[test]
    fn division_visits_keys_the_numerator_lacks() {
        let n: QtPoly<Rational> = binomial(3, 0);
        let q = divide_by_factor(&n, 1, 0).expect("(1-q) divides (1-q^3)");
        for a in 0..3 {
            assert_eq!(q.coeff(a, 0), r(1), "coefficient of q^{a}");
        }
        assert_eq!(q.len(), 3, "{q}");
    }

    /// Several independent chains at once, with a divisor that moves in both
    /// variables — the case where the chains genuinely interleave in lex order
    /// and the pieces have to be re-sorted before they are a `QtPoly`.
    #[test]
    fn division_handles_interleaved_chains() {
        // (1 + q + t)·(1 − q t) — three chains under δ = (1,1)
        let mut f: QtPoly<Rational> = <QtPoly<Rational> as Ring>::one();
        f.add_term(1, 0, r(1));
        f.add_term(0, 1, r(1));
        let n = f.mul_binomial(1, 1);
        let got = divide_by_factor(&n, 1, 1).expect("(1-qt) divides its own multiple");
        assert_eq!(got, f, "got {got}, want {f}");
    }

    /// Non-divisibility must be *detected*, not looped on, including when the
    /// chains all start below the degree bound.
    #[test]
    fn non_divisibility_is_detected() {
        let mut f: QtPoly<Rational> = <QtPoly<Rational> as Ring>::one();
        f.add_term(2, 3, r(1)); // 1 + q^2 t^3, divisible by no 1 - q^a t^b
        for (a, b) in [(1, 0), (0, 1), (1, 1), (2, 3), (1, 2)] {
            assert!(
                divide_by_factor(&f, a, b).is_none(),
                "(1 - q^{a} t^{b}) must not divide {f}"
            );
        }
    }

    /// Round trip on a spread of shapes: multiplying by a binomial and dividing
    /// it back out must return the original exactly.
    #[test]
    fn division_inverts_multiplication() {
        let mut dense: QtPoly<Rational> = QtPoly::zero();
        for a in 0..4 {
            for b in 0..3 {
                dense.add_term(a, b, r(i128::from(a) - i128::from(b) + 1));
            }
        }
        let mut sparse: QtPoly<Rational> = QtPoly::zero();
        sparse.add_term(0, 0, r(1));
        sparse.add_term(4, 1, r(-3));
        sparse.add_term(1, 6, r(5));
        for f in [dense, sparse] {
            for (a, b) in [(1, 0), (0, 1), (1, 1), (2, 1), (3, 3)] {
                let prod = f.mul_binomial(a, b);
                let back = divide_by_factor(&prod, a, b)
                    .unwrap_or_else(|| panic!("(1 - q^{a} t^{b}) must divide its own multiple"));
                assert_eq!(back, f, "round trip through (1 - q^{a} t^{b})");
            }
        }
    }

    /// The two exact divisions must agree wherever both apply.
    ///
    /// [`divide_by_factor`] is the chain walk specialised to `1 − qᵃtᵇ`;
    /// [`QtPoly::divide_exact`] is leading-term elimination against an arbitrary
    /// divisor. They share no code and no idea — one runs along arithmetic
    /// progressions of exponents, the other down a monomial order — so agreement
    /// between them is evidence and not tautology. This is the check that lets
    /// the general routine be trusted where the specialised one cannot reach.
    ///
    /// Both the divisible and the non-divisible cases: a general divider that
    /// silently returned a truncated quotient would pass a round-trip test.
    #[test]
    fn the_two_exact_divisions_agree() {
        let mut dense: QtPoly<Rational> = QtPoly::zero();
        for a in 0..4 {
            for b in 0..3 {
                dense.add_term(a, b, r(i128::from(a) - i128::from(b) + 1));
            }
        }
        let mut sparse: QtPoly<Rational> = QtPoly::zero();
        sparse.add_term(0, 0, r(1));
        sparse.add_term(4, 1, r(-3));
        sparse.add_term(1, 6, r(5));
        for f in [dense, sparse] {
            for (a, b) in [(1, 0), (0, 1), (1, 1), (2, 1), (3, 3)] {
                let d = binomial::<Rational>(a, b);
                // Divisible: same quotient, not merely both succeeding.
                let prod = f.mul_binomial(a, b);
                assert_eq!(
                    prod.divide_exact(&d),
                    divide_by_factor(&prod, a, b),
                    "(1 - q^{a} t^{b}) dividing its own multiple"
                );
                // Not divisible: both must decline. `f` itself is a multiple of
                // no binomial here — the round trip above is what says the
                // divisible case is not vacuous.
                assert_eq!(
                    f.divide_exact(&d),
                    divide_by_factor(&f, a, b),
                    "(1 - q^{a} t^{b}) against a non-multiple"
                );
            }
        }
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
        assert_eq!(
            den.count(),
            0,
            "the denominator must be gone, not merely equal"
        );
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
