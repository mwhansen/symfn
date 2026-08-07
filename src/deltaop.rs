//! The Macdonald operator algebra: ∇, Δ_f, Δ'_f, Π and Θ_f.
//!
//! ## What they are
//!
//! All of them are **diagonal on the modified Macdonald basis** `{H̃_μ}`, which
//! [`bh::htilde_table`](crate::bh) already produces for a whole degree:
//!
//! ```text
//!   ∇ H̃_μ    = T_μ H̃_μ                T_μ = ∏_{c∈μ} q^{a'(c)} t^{l'(c)}
//!   Δ_f H̃_μ  = f[B_μ] H̃_μ             B_μ = Σ_{c∈μ} q^{a'(c)} t^{l'(c)}
//!   Δ'_f H̃_μ = f[B_μ − 1] H̃_μ
//!   Π H̃_μ    = Π_μ H̃_μ                Π_μ = ∏_{c≠(0,0)} (1 − q^{a'} t^{l'})
//!   Θ_f F    = Π f* Π⁻¹ F              f* = f[X/M],  M = (1−q)(1−t)
//! ```
//!
//! `a'` and `l'` are the **co**arm and **co**leg — the column and row indices,
//! both 0-based. `w_μ` below is the only quantity built from the ordinary arm
//! and leg, and swapping the two families is silent: `T_μ` built from arms and
//! legs is still a monomial, still distinct across μ, and still yields a
//! plausible answer. `t_mu_is_the_hook_monomial` is what pins it.
//!
//! `Δ'` needs **no virtual alphabet**. `B_μ − 1` looks like a formal alphabet
//! with a negative letter, which [`Ratio`] cannot hold; it does not have to,
//! because the cell `(0,0)` contributes exactly `q⁰t⁰ = 1` to `B_μ`, so
//! `B_μ − 1` is `B_μ` with that cell deleted — the same cell set `Π_μ` skips.
//!
//! ## The expansion into `H̃`, and why there is no linear solve
//!
//! Applying any of these to an arbitrary `F` means writing `F` in the `H̃`
//! basis. The obvious route inverts the `K̃` matrix — `p(n) × p(n)` over
//! ℚ(q,t), 77×77 at degree 12. It is not needed. `H̃` is **orthogonal** for the
//! star scalar product [IR], so the expansion is diagonal:
//!
//! ```text
//!   ⟨p_ρ, p_σ⟩_* = δ_{ρσ} · z_ρ · (−1)^{|ρ|−ℓ(ρ)} · ∏_i (1−q^{ρ_i})(1−t^{ρ_i})
//!   ⟨H̃_μ, H̃_ν⟩_* = δ_{μν} · w_μ
//!   F = Σ_μ ( ⟨F, H̃_μ⟩_* / w_μ ) H̃_μ
//! ```
//!
//! and the `z_ρ` cancels out of the pairing, which is what keeps it integral:
//! pairing `F` against `s_κ` rather than against `p_ρ` leaves
//! `Σ_ρ F_ρ χ^κ_ρ ε_ρ W_ρ`, with no `z_ρ⁻¹` in it.
//!
//! So `H̃_μ` is never converted to the power sums at all. `F` is converted
//! once; everything after that is `⟨F,s_κ⟩_*` (integer-scaled additions)
//! followed by `Σ_κ K̃_{κμ}⟨F,s_κ⟩_*`. Converting all `p(n)` of the `H̃_μ`
//! instead is `p(n)³` scaled additions on polynomials that are never needed in
//! that basis.
//!
//! ## Arithmetic
//!
//! Exactly two families of denominator arise and [`Atom`] is their union:
//! `1 − qᵃtᵇ` from `M`, `Π_μ`, `f[X/M]` and the star weights, and `qᵃ − tᵇ`
//! from `w_μ`. [`Frac`](crate::Frac) closes over the first and
//! [`bh::Rat`](crate::bh) over the second; this is the first thing in the crate
//! that needs both at once, so [`Ratio`] holds a denominator as a multiset over
//! the union and never expands it.
//!
//! The two families **overlap**, and normalizing the overlap away is
//! necessary, not tidiness: `q^a − 1` is `−(1 − q^a)` and `q⁰ − t^b` is
//! `1 − t^b`, both of which `w_μ` produces, and both of which the star weights
//! produce independently. Left as distinct atoms they would never cancel
//! against each other. [`Atom::diff`] is where that normalization happens.
//!
//! ## Range
//!
//! Over a fixed-width `C` these operators refuse rather than wrapping past
//! their wall (`docs/policies/failure.md`, R3), and runtime arrives first by a
//! wide margin: at `i128` both [`nabla_e`] and [`delta_prime_e`] gain ~2.8 bits
//! per degree, reaching 127 bits near n ≈ 51–52, while the computations stop
//! finishing around n = 16. Measurements and the harness are in
//! `docs/record/failure-and-overflow.md` (`examples/probe_qt_walls.rs`).
//!
//! ## Reduction policy, which is the opposite of `Frac`'s
//!
//! [`Frac::add_assign`](crate::Frac) deliberately does *not* reduce: a running
//! sum should be reduced once at the end, because trial division is the
//! expensive operation and a cancellation can only be decided once the sum is
//! complete. **Here the accumulator reduces after every term**, and the
//! difference is not a matter of taste.
//!
//! The sum being accumulated is `Σ_μ (n_μ / w_μ) K̃_{λμ}` over all `p(n)`
//! shapes, and the `w_μ` are largely coprime. Without reduction the denominator
//! grows to their lcm and every earlier numerator is lifted against all of it.
//! Reducing every step holds the peak numerator to 1393 terms and the
//! denominator to 15 atoms at degree 9 (`docs/record/macdonald-operators.md`),
//! and the denominator cancels to nothing at the end — which it must, since
//! the answer is a polynomial. That is [`macop::Coeff`](crate::macop)'s policy,
//! for [`macop`](crate::macop)'s reason.
//!
//! ## References
//!
//! - **[DM]** M. D'Adderio, A. Mellit, *A proof of the compositional Delta
//!   conjecture*, [arXiv:2011.11467](https://arxiv.org/abs/2011.11467) —
//!   (5)–(8) `M, B_μ, T_μ, Π_μ`; (10) `f*`; (11) ∇; (12) Δ_f and Δ'_f;
//!   (21) **Π**; (22) Θ_f.
//! - **[DIV]** M. D'Adderio, A. Iraci, A. Vanden Wyngaerd, *Theta operators,
//!   refined Delta conjectures, and coinvariants*,
//!   [arXiv:1906.02623](https://arxiv.org/abs/1906.02623) — (9) `w_μ`.
//! - **[IR]** A. Iraci, M. Romero, *Delta and Theta operator expansions*,
//!   [arXiv:2203.10342](https://arxiv.org/abs/2203.10342) — the star scalar
//!   product `⟨F,G⟩_* = ⟨F,(ωG)[MX]⟩`.
//!
//! `docs/record/macdonald-operators.md` has the measured record, and
//! `scripts/verify_deltaop_formulas.py` verified every formula below against
//! Sage before any of it was written. Sage's entry point into this family is
//! the `nabla` method on a symmetric function, and it is the only one: Δ_f,
//! Δ'_f, Π, Θ_f and the star scalar product have no Sage equivalent, which that
//! record establishes by measurement rather than by assumption.
//!
//! [DIV]: https://arxiv.org/abs/1906.02623
//! [DM]: https://arxiv.org/abs/2011.11467
//! [IR]: https://arxiv.org/abs/2203.10342

// Shape indices; the coefficients are `Ratio<C>` and are never cast.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::BTreeMap;

use crate::coeff::{QAlgebra, Ring};
use crate::convert::FromSchur;
use crate::partition::Partition;
use crate::qt::QtPoly;
use crate::sym::{PowerSum, Schur, SymFn};

// ---------------------------------------------------------------------------
// Atoms
// ---------------------------------------------------------------------------

/// A denominator factor: one of the two binomial families this module divides
/// by, in a normalized form that lets them cancel against each other.
///
/// [`Frac`](crate::Frac) exists for the first family and [`bh::Rat`](crate::bh)
/// for the second; see the module docs on why both are needed here and why the
/// normalization matters.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub enum Atom {
    /// `1 − qᵃtᵇ`, with `(a,b) ≠ (0,0)`.
    Unit(u32, u32),
    /// `qᵃ − tᵇ`, with `a ≥ 1` and `b ≥ 1` — the cases where it is *not* one of
    /// the above. See [`Atom::diff`].
    Diff(u32, u32),
}

impl Atom {
    /// `1 − qᵃtᵇ`.
    ///
    /// # Panics
    ///
    /// Panics on `(0, 0)`, which is zero.
    pub fn unit(a: u32, b: u32) -> Self {
        assert!(a > 0 || b > 0, "1 - q^0 t^0 is zero");
        Atom::Unit(a, b)
    }

    /// `qᵃ − tᵇ` as `(atom, negated)`: the value is `−atom` when `negated`.
    ///
    /// The two boundary cases are why this is a constructor rather than a
    /// variant. `q⁰ − tᵇ` *is* `1 − tᵇ` and `qᵃ − t⁰` is `−(1 − qᵃ)`, so
    /// both belong to the `Unit` family — and `w_μ` produces both (a cell with
    /// zero arm gives the first, zero leg the second) while the star weights
    /// produce `1 − q^k` and `1 − t^k` directly. Keeping them as separate atoms
    /// would leave two names for one polynomial and no cancellation between
    /// them.
    ///
    /// # Panics
    ///
    /// Panics on `(0, 0)`: `q⁰ − t⁰` is zero, and unlike the two boundary cases
    /// above it belongs to no atom family.
    pub fn diff(a: u32, b: u32) -> (Self, bool) {
        assert!(a > 0 || b > 0, "q^0 - t^0 is zero");
        if a == 0 {
            (Atom::Unit(0, b), false)
        } else if b == 0 {
            (Atom::Unit(a, 0), true)
        } else {
            (Atom::Diff(a, b), false)
        }
    }

    /// The atom as a polynomial. Only for display, testing, and the general
    /// division path — the point of the factored form is not to do this.
    pub fn poly<C: Ring>(self) -> QtPoly<C> {
        match self {
            Atom::Unit(a, b) => {
                let mut p = <QtPoly<C> as Ring>::one();
                p.add_term(a, b, C::one().neg());
                p
            }
            Atom::Diff(a, b) => {
                let mut p = QtPoly::term(a, 0, C::one());
                p.add_term(0, b, C::one().neg());
                p
            }
        }
    }

    /// `n · atom`, through the specialized two-run merge for each family —
    /// never [`Ring::mul`], which would sort a concatenation of two already
    /// sorted runs. See [`QtPoly::mul_diff`](crate::qt::QtPoly::mul_diff) for
    /// what that cost when measured.
    fn mul_into<C: Ring>(self, n: &QtPoly<C>) -> QtPoly<C> {
        match self {
            Atom::Unit(a, b) => n.mul_binomial(a, b),
            Atom::Diff(a, b) => n.mul_diff(a, b),
        }
    }

    /// `n / atom`, or `None` if it does not divide exactly.
    ///
    /// Two different algorithms, because the families are genuinely different
    /// and each already has the right one in the crate.
    /// [`frac::divide_by_factor`](crate::frac) walks arithmetic progressions of
    /// exponents and is specialized to `1 − qᵃtᵇ`; a factor `qᵃ − tᵇ` is not of
    /// that shape and goes through
    /// [`QtPoly::divide_exact`](crate::qt::QtPoly::divide_exact), the general
    /// leading-term elimination that [`bh`](crate::bh) uses for the same
    /// reason.
    fn divide<C: Ring>(self, n: &QtPoly<C>) -> Option<QtPoly<C>> {
        match self {
            Atom::Unit(a, b) => crate::frac::divide_by_factor(n, a, b),
            Atom::Diff(a, b) => {
                if !crate::frac::diff_may_divide(n, a, b) {
                    return None;
                }
                crate::frac::divide_by_diff(n, a, b)
            }
        }
    }
}

impl core::fmt::Display for Atom {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Atom::Unit(a, b) => write!(f, "(1-q^{a}t^{b})"),
            Atom::Diff(a, b) => write!(f, "(q^{a}-t^{b})"),
        }
    }
}

/// A multiset of atoms, `atom ↦ multiplicity`.
type Atoms = BTreeMap<Atom, u32>;

fn push(map: &mut Atoms, atom: Atom, times: u32) {
    if times > 0 {
        *map.entry(atom).or_insert(0) += times;
    }
}

// ---------------------------------------------------------------------------
// Ratio
// ---------------------------------------------------------------------------

/// An element of ℚ(q,t) whose denominator is a product of [`Atom`]s, held
/// factored and never expanded.
///
/// The same design as [`Frac`](crate::Frac) over a wider family — see the
/// module docs — including the reason equality must cross-multiply: the
/// representation is not canonical, because `1 − q²` and `(1−q)(1+q)` are the
/// same polynomial with different atom multisets.
#[derive(Clone, Debug)]
pub struct Ratio<C: Ring> {
    num: QtPoly<C>,
    den: Atoms,
}

impl<C: Ring> Ratio<C> {
    /// A polynomial, as a fraction with denominator 1.
    pub fn from_poly(num: QtPoly<C>) -> Self {
        Ratio {
            num,
            den: Atoms::new(),
        }
    }

    /// `num / ∏ atoms`, unreduced.
    fn over(num: QtPoly<C>, den: Atoms) -> Self {
        Ratio { num, den }
    }

    /// The numerator and the denominator's factors with their multiplicities.
    ///
    /// The factors come in ascending [`Atom`] order.
    pub fn parts(&self) -> (&QtPoly<C>, impl Iterator<Item = (&Atom, &u32)>) {
        (&self.num, self.den.iter())
    }

    /// Multiply by `∏ atoms` (into the numerator).
    pub fn mul_atoms(&self, atoms: &Atoms) -> Self {
        let mut num = self.num.clone();
        for (&a, &m) in atoms {
            for _ in 0..m {
                num = a.mul_into(&num);
            }
        }
        Ratio {
            num,
            den: self.den.clone(),
        }
    }

    /// Divide by `∏ atoms` (into the denominator).
    pub fn div_atoms(&self, atoms: &Atoms) -> Self {
        let mut den = self.den.clone();
        for (&a, &m) in atoms {
            push(&mut den, a, m);
        }
        Ratio {
            num: self.num.clone(),
            den,
        }
    }

    /// `self` rewritten over `target`, which must be a multiple of `self.den`.
    fn lift(&self, target: &Atoms) -> QtPoly<C> {
        let mut num = self.num.clone();
        for (&a, &m) in target {
            for _ in 0..(m - self.den.get(&a).copied().unwrap_or(0)) {
                num = a.mul_into(&num);
            }
        }
        num
    }

    /// `self += other · p`, over the lcm of the two denominators.
    ///
    /// The lcm is a *common multiple* and not always the least one, because the
    /// atoms are not pairwise coprime — `q² − t²` is `(q−t)(q+t)`. Correctness
    /// does not need minimality (only a common multiple is required to add) and
    /// [`reduce`](Self::reduce) takes the excess back out. `bh::Rat` lives with
    /// exactly this. **Multiply before lifting.** `p` is a `K̃` entry — tens of
    /// terms — and `other.num` is comparable, while the lifted form is the size
    /// of the whole running answer. `lift(other) · p` and `lift(other · p)` are
    /// the same element, but the first runs a general product with a
    /// several-thousand-term operand and quicksorts the result, and the second
    /// runs it on the small pair and then makes a handful of linear merging
    /// passes.
    ///
    /// The two orders are timed against each other in
    /// `docs/record/macdonald-operators.md`.
    fn add_mul(&mut self, other: &Self, p: &QtPoly<C>) {
        if other.num.is_empty() || p.is_empty() {
            return;
        }
        let scaled = other.num.mul(p);
        if self.num.is_empty() && self.den.is_empty() {
            self.num = scaled;
            self.den = other.den.clone();
            return;
        }
        if other.den != self.den {
            let mut lcm = self.den.clone();
            for (&a, &m) in &other.den {
                let e = lcm.entry(a).or_insert(0);
                *e = (*e).max(m);
            }
            self.num = self.lift(&lcm);
            self.den = lcm;
        }
        let mut lifted = scaled;
        for (&a, &m) in &self.den {
            for _ in 0..(m - other.den.get(&a).copied().unwrap_or(0)) {
                lifted = a.mul_into(&lifted);
            }
        }
        self.num.add_assign(&lifted);
    }

    /// Divide out every denominator atom that also divides the numerator.
    pub fn reduce(&mut self) {
        if self.num.is_empty() {
            self.den.clear();
            return;
        }
        self.den.retain(|&a, m| {
            while *m > 0 {
                match a.divide(&self.num) {
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

    /// The numerator, if the denominator cancels away entirely.
    ///
    /// Reduces first, so it answers about the *element* rather than the
    /// representation. Callers that know on mathematical grounds that the
    /// answer is a polynomial — ∇, Δ, Δ', Θ of an integral input — should
    /// unwrap and let a `None` be the loud failure it is.
    pub fn into_poly(mut self) -> Option<QtPoly<C>> {
        self.reduce();
        self.den.is_empty().then_some(self.num)
    }
}

impl<C: Ring> PartialEq for Ratio<C> {
    fn eq(&self, other: &Self) -> bool {
        if self.den == other.den {
            return self.num == other.num;
        }
        let mut lcm = self.den.clone();
        for (&a, &m) in &other.den {
            let e = lcm.entry(a).or_insert(0);
            *e = (*e).max(m);
        }
        self.lift(&lcm) == other.lift(&lcm)
    }
}

impl<C: Ring> Eq for Ratio<C> {}

impl<C: Ring> Ring for Ratio<C> {
    fn zero() -> Self {
        Ratio {
            num: QtPoly::zero(),
            den: Atoms::new(),
        }
    }
    fn one() -> Self {
        Ratio {
            num: <QtPoly<C> as Ring>::one(),
            den: Atoms::new(),
        }
    }
    fn is_zero(&self) -> bool {
        self.num.is_zero()
    }
    fn add_assign(&mut self, other: &Self) {
        let one = <QtPoly<C> as Ring>::one();
        self.add_mul(other, &one);
    }
    fn mul(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Self::zero();
        }
        let mut den = self.den.clone();
        for (&a, &m) in &other.den {
            push(&mut den, a, m);
        }
        let mut r = Ratio {
            num: self.num.mul(&other.num),
            den,
        };
        // Reduced here and not in `add_assign`, the same split `Frac` documents
        // and for the same measured reason: a value produced by a product is
        // about to be used many times, and cutting its denominator down once is
        // cheaper than lifting against it repeatedly.
        r.reduce();
        r
    }
    fn neg(&self) -> Self {
        Ratio {
            num: self.num.neg(),
            den: self.den.clone(),
        }
    }
    fn from_i64(n: i64) -> Self {
        Ratio::from_poly(<QtPoly<C> as Ring>::from_i64(n))
    }
    fn from_u128(n: u128) -> Self {
        Ratio::from_poly(<QtPoly<C> as Ring>::from_u128(n))
    }
    fn from_i128(n: i128) -> Self {
        Ratio::from_poly(<QtPoly<C> as Ring>::from_i128(n))
    }
}

impl<C: QAlgebra> QAlgebra for Ratio<C> {
    fn div_u128(&self, n: u128) -> Self {
        Ratio {
            num: self.num.div_u128(n),
            den: self.den.clone(),
        }
    }
}

impl<C: Ring> core::fmt::Display for Ratio<C> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "({})", self.num)?;
        for (a, &m) in &self.den {
            write!(f, "/{a}")?;
            if m > 1 {
                write!(f, "^{m}")?;
            }
        }
        Ok(())
    }
}

// ---------------------------------------------------------------------------
// Cell statistics
// ---------------------------------------------------------------------------

/// Every cell of μ as its `(coarm, coleg) = (column, row)`, both 0-based.
///
/// These are the statistics `B_μ`, `T_μ` and `Π_μ` are built from — *not* the
/// arm and leg, which only `w_μ` uses.
fn coarms(mu: &Partition) -> Vec<(u32, u32)> {
    let mut out = Vec::with_capacity(mu.size() as usize);
    for (i, &row) in mu.parts().iter().enumerate() {
        for j in 0..row {
            out.push((j, i as u32));
        }
    }
    out
}

/// `B_μ = Σ_{c∈μ} q^{a'(c)} t^{l'(c)}`, \[DM\] (6).
///
/// Every cell gives a distinct monomial, so every coefficient is 1.
fn b_mu<C: Ring>(cells: &[(u32, u32)]) -> QtPoly<C> {
    let mut out = QtPoly::zero();
    for &(a, l) in cells {
        out.add_term(a, l, C::one());
    }
    out
}

/// `T_μ = ∏_{c∈μ} q^{a'(c)} t^{l'(c)} = q^{n(μ')} t^{n(μ)}`, \[DM\] (7).
fn t_mu<C: Ring>(cells: &[(u32, u32)]) -> QtPoly<C> {
    let (mut a, mut l) = (0, 0);
    for &(x, y) in cells {
        a += x;
        l += y;
    }
    QtPoly::term(a, l, C::one())
}

/// `Π_μ = ∏_{c ≠ (0,0)} (1 − q^{a'} t^{l'})`, \[DM\] (8), as atoms.
fn pi_atoms(cells: &[(u32, u32)]) -> Atoms {
    let mut out = Atoms::new();
    for &(a, l) in cells {
        if (a, l) != (0, 0) {
            push(&mut out, Atom::unit(a, l), 1);
        }
    }
    out
}

/// `w_μ = ∏_{c∈μ} (q^{a} − t^{l+1})(t^{l} − q^{a+1})`, \[DIV\] (9), as
/// `(atoms, negated)`.
///
/// The second factor is `−(q^{a+1} − t^{l})`, so the product carries a global
/// `(−1)^{|μ|}` before the per-atom normalization of [`Atom::diff`] contributes
/// its own signs. **Arm and leg here**, unlike everything else in this module.
fn w_atoms(mu: &Partition) -> (Atoms, bool) {
    let shape = mu.parts();
    let mut out = Atoms::new();
    let mut negated = false;
    for i in 0..shape.len() {
        for j in 0..shape[i] as usize {
            let a = crate::macdonald::arm(shape, i, j);
            let l = crate::macdonald::leg(shape, i, j);
            for (x, y) in [(a, l + 1), (a + 1, l)] {
                let (atom, neg) = Atom::diff(x, y);
                push(&mut out, atom, 1);
                negated ^= neg;
            }
            // the `(t^l − q^{a+1}) = −(q^{a+1} − t^l)` sign, once per cell
            negated = !negated;
        }
    }
    (out, negated)
}

/// `M = (1 − q)(1 − t)`, as atoms.
fn m_atoms() -> Atoms {
    let mut out = Atoms::new();
    push(&mut out, Atom::unit(1, 0), 1);
    push(&mut out, Atom::unit(0, 1), 1);
    out
}

// ---------------------------------------------------------------------------
// The star scalar product
// ---------------------------------------------------------------------------

/// `⟨F, s_κ⟩_*` for every κ of the degree, indexed against
/// `memo::partitions_cached`.
///
/// ```text
///   ⟨F, s_κ⟩_* = Σ_ρ F_ρ · χ^κ_ρ · (−1)^{|ρ|−ℓ(ρ)} · ∏_i (1−q^{ρ_i})(1−t^{ρ_i})
/// ```
///
/// **`z_ρ` does not appear.** `⟨p_ρ,p_ρ⟩_*` carries a `z_ρ` and `s_κ`'s
/// `p_ρ`-coefficient carries a `z_ρ⁻¹`, and they cancel. That is what keeps
/// this step free of division and is why the pairing is taken against the Schur
/// basis rather than shape by shape against `H̃_μ`.
fn star_against_schur<C: QAlgebra>(f: &Schur<Ratio<C>>, n: u32) -> Vec<Ratio<C>> {
    let parts = crate::memo::partitions_cached(n);
    let fp: PowerSum<Ratio<C>> = PowerSum::from_schur(f);

    // v_ρ = F_ρ · ε_ρ · W_ρ, once per ρ.
    let weighted: Vec<(Partition, Ratio<C>)> = fp
        .terms()
        .iter()
        .map(|(rho, c)| {
            let mut w = Atoms::new();
            for &part in rho.parts() {
                push(&mut w, Atom::unit(part, 0), 1);
                push(&mut w, Atom::unit(0, part), 1);
            }
            let mut v = c.mul_atoms(&w);
            if (rho.size() - rho.len() as u32) % 2 == 1 {
                v = v.neg();
            }
            // Reduced once, before p(n) uses of it below — the policy `Frac`'s
            // `mul` documents.
            v.reduce();
            (rho.clone(), v)
        })
        .collect();

    parts
        .iter()
        .map(|kappa| {
            let mut acc = <Ratio<C> as Ring>::zero();
            for (rho, v) in &weighted {
                let chi = crate::character::character(kappa, rho);
                if chi != 0 {
                    acc.add_mul(v, &QtPoly::term(0, 0, C::from_i128(chi)));
                }
            }
            acc.reduce();
            acc
        })
        .collect()
}

/// The `H̃`-basis coefficients of `f`: `c_μ = ⟨f,H̃_μ⟩_* / w_μ`, in the order
/// of `memo::partitions_cached`.
fn coefficients<C: QAlgebra>(
    f: &Schur<Ratio<C>>,
    n: u32,
    htilde: &[(Partition, Schur<QtPoly<C>>)],
) -> Vec<Ratio<C>> {
    let u = star_against_schur(f, n);
    let parts = crate::memo::partitions_cached(n);
    htilde
        .iter()
        .map(|(mu, ht)| {
            // ⟨f, H̃_μ⟩_* = Σ_κ K̃_{κμ} ⟨f, s_κ⟩_*
            let mut acc = <Ratio<C> as Ring>::zero();
            for (k, kappa) in parts.iter().enumerate() {
                let kt = ht.coeff(kappa);
                if !kt.is_empty() {
                    acc.add_mul(&u[k], &kt);
                }
            }
            acc.reduce();
            let (w, negated) = w_atoms(mu);
            let mut c = acc.div_atoms(&w);
            if negated {
                c = c.neg();
            }
            c.reduce();
            c
        })
        .collect()
}

/// `Σ_μ c_μ · H̃_μ`, back in the Schur basis.
///
/// One accumulator per λ, each reduced after every μ. See the module docs on
/// why that is the opposite of `Frac`'s policy and why it is required.
fn combine<C: Ring>(
    coeffs: &[Ratio<C>],
    n: u32,
    htilde: &[(Partition, Schur<QtPoly<C>>)],
) -> Schur<Ratio<C>> {
    let parts = crate::memo::partitions_cached(n);
    let mut out = Schur::zero();
    for lambda in parts.iter() {
        let items: Vec<Ratio<C>> = coeffs
            .iter()
            .zip(htilde.iter())
            .filter(|(c, _)| !c.is_zero())
            .filter_map(|(c, (_, ht))| {
                let kt = ht.coeff(lambda);
                (!kt.is_empty()).then(|| {
                    let mut r = Ratio::over(c.num.mul(&kt), c.den.clone());
                    r.reduce();
                    r
                })
            })
            .collect();
        out.add_term(lambda.clone(), sum_tree(items));
    }
    out
}

/// Add a list of fractions **pairwise in a balanced tree**, not as a running
/// sum.
///
/// Both orders do `p(n) − 1` additions, and the difference is entirely in what
/// each addition costs. A running sum reaches the size of the whole answer
/// after a handful of terms and then pays a lift *and a full trial-division
/// sweep* at that size for every one of the remaining `p(n)` steps, against a
/// denominator carrying atoms accumulated from every `w_μ` seen so far. A tree
/// does most of its work near the leaves, where both the numerators and the
/// atom multisets are small, and only the last few merges are full size.
///
/// The reduce sweeps are quadratic in `p(n)` under a running sum and
/// `O(p(n) log p(n))` here; the two orders are timed against each other in
/// `docs/record/macdonald-operators.md`.
fn sum_tree<C: Ring>(mut items: Vec<Ratio<C>>) -> Ratio<C> {
    let one = <QtPoly<C> as Ring>::one();
    while items.len() > 1 {
        let mut next: Vec<Ratio<C>> = Vec::with_capacity(items.len().div_ceil(2));
        let mut it = items.into_iter();
        while let Some(mut a) = it.next() {
            if let Some(b) = it.next() {
                a.add_mul(&b, &one);
                a.reduce();
            }
            next.push(a);
        }
        items = next;
    }
    items.pop().unwrap_or_else(<Ratio<C> as Ring>::zero)
}

// ---------------------------------------------------------------------------
// The generic operator
// ---------------------------------------------------------------------------

/// The degree of a homogeneous element, or a panic naming the caller's error.
fn degree_of<C: Ring, S: SymFn<C>>(f: &S, what: &str) -> Option<u32> {
    let mut it = f.terms().keys().map(Partition::size);
    let first = it.next()?;
    assert!(
        it.all(|d| d == first),
        "{what} requires a homogeneous argument"
    );
    Some(first)
}

/// Expand in `H̃`, scale the μ coefficient by `eigen(μ)`, sum back.
fn diagonal<C: QAlgebra>(
    f: &Schur<Ratio<C>>,
    eigen: impl Fn(&Partition, &[(u32, u32)]) -> Ratio<C>,
) -> Schur<Ratio<C>> {
    let n = match degree_of(f, "a Macdonald operator") {
        Some(n) => n,
        None => return Schur::zero(),
    };
    let htilde = crate::bh::htilde_table::<C>(n);
    let mut coeffs = coefficients(f, n, &htilde);
    for (c, (mu, _)) in coeffs.iter_mut().zip(htilde.iter()) {
        let e = eigen(mu, &coarms(mu));
        *c = c.mul(&e);
    }
    combine(&coeffs, n, &htilde)
}

fn lift_in<C: Ring>(f: &Schur<QtPoly<C>>) -> Schur<Ratio<C>> {
    Schur::from_terms(
        f.terms()
            .iter()
            .map(|(l, c)| (l.clone(), Ratio::from_poly(c.clone())))
            .collect(),
    )
}

/// Bring a result back to polynomials, panicking on a surviving denominator —
/// which is a bug in the mathematics above, not a case to handle.
fn lift_out<C: Ring>(f: Schur<Ratio<C>>, what: &str) -> Schur<QtPoly<C>> {
    let mut out = Schur::zero();
    for (lambda, c) in f.terms() {
        let p = c.clone().into_poly().unwrap_or_else(|| {
            panic!("{what} is not a polynomial at {lambda}: {c}");
        });
        out.add_term(lambda.clone(), p);
    }
    out
}

// ---------------------------------------------------------------------------
// Eigenvalues
// ---------------------------------------------------------------------------

/// `f[B]` where `B` is the multiset of cell monomials `q^{a'} t^{l'}`.
///
/// Through the power sums, where the substitution is one line:
/// `p_k[B] = Σ_c q^{k a'} t^{k l'}`, and `f[B] = Σ_ρ c_ρ ∏_i p_{ρ_i}[B]`.
///
/// `f` carries integer coefficients — see the module docs and
/// `docs/record/macdonald-operators.md` on why the signature refuses
/// ℚ(q,t) there: `f[·]` is a plethysm at a formal alphabet, `q` and `t` are
/// letters of that alphabet, and for an `f` with `(q,t)` coefficients the two
/// readings are different operators. For constant coefficients the distinction
/// is vacuous.
fn plethystic_eval<C: QAlgebra>(f: &Schur<i128>, cells: &[(u32, u32)]) -> QtPoly<C> {
    let lifted: Schur<C> = Schur::from_terms(
        f.terms()
            .iter()
            .map(|(l, &c)| (l.clone(), C::from_i128(c)))
            .collect(),
    );
    let fp: PowerSum<C> = PowerSum::from_schur(&lifted);
    let mut out = QtPoly::zero();
    for (rho, c) in fp.terms() {
        let mut term = QtPoly::term(0, 0, c.clone());
        for &k in rho.parts() {
            let mut pk = QtPoly::zero();
            for &(a, l) in cells {
                pk.add_term(a * k, l * k, C::one());
            }
            term = term.mul(&pk);
        }
        out.add_assign(&term);
    }
    out
}

/// The cells `Δ'` and `Π` use: every cell but `(0,0)`.
fn without_corner(cells: &[(u32, u32)]) -> Vec<(u32, u32)> {
    cells.iter().copied().filter(|&c| c != (0, 0)).collect()
}

/// `e_k[B]` for `B` the multiset of cell monomials — **over ℤ**, with no
/// power-sum detour.
///
/// [`plethystic_eval`] handles an arbitrary `f` and pays for it: `s → p`
/// divides by `z_ρ`, so it needs a [`QAlgebra`]. `e_k` is the subscript the
/// Delta-conjecture literature uses, and for it the answer is one
/// coefficient of `∏_c (1 + z·m_c)` — a `k`-term convolution over the cells
/// with no division anywhere. That is what lets [`nabla_e`] and
/// [`delta_prime_e`] run over `QtPoly<i128>`, which is the same move
/// [`qt_kostka_table_via_bh`](crate::qtkostka) makes and for the same measured
/// reason (see [`nabla_e`]).
fn elementary_eval<C: Ring>(cells: &[(u32, u32)], k: u32) -> QtPoly<C> {
    let k = k as usize;
    let mut acc: Vec<QtPoly<C>> = vec![QtPoly::zero(); k + 1];
    acc[0] = <QtPoly<C> as Ring>::one();
    for &(a, l) in cells {
        let m = QtPoly::term(a, l, C::one());
        for j in (1..=k).rev() {
            let add = acc[j - 1].mul(&m);
            acc[j].add_assign(&add);
        }
    }
    acc.swap_remove(k)
}

// ---------------------------------------------------------------------------
// Public operators
// ---------------------------------------------------------------------------

/// `∇F`, \[DM\] (11).
///
/// Diagonal on `H̃` with eigenvalue `T_μ`, a monomial — so ∇ is the cheapest of
/// the family and everything it costs is the change of basis.
///
/// # Panics
///
/// Panics if `f` is not homogeneous.
pub fn nabla<C: QAlgebra>(f: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>> {
    let got = diagonal(&lift_in(f), |_, cells| Ratio::from_poly(t_mu::<C>(cells)));
    lift_out(got, "nabla")
}

/// `∇^r F`, ∇ applied `r` times.
///
/// Applying [`nabla`] `r` times would redo the change of basis every time; ∇ is
/// diagonal, so the `r`th power is the `r`th power of the eigenvalue and one
/// expansion suffices.
///
/// # Panics
///
/// Panics if `f` is not homogeneous.
pub fn nabla_power<C: QAlgebra>(f: &Schur<QtPoly<C>>, r: u32) -> Schur<QtPoly<C>> {
    let got = diagonal(&lift_in(f), |_, cells| {
        let t = t_mu::<C>(cells);
        let mut acc = <QtPoly<C> as Ring>::one();
        for _ in 0..r {
            acc = acc.mul(&t);
        }
        Ratio::from_poly(acc)
    });
    lift_out(got, "nabla_power")
}

/// `Δ_f F`, \[DM\] (12) — eigenvalue `f[B_μ]`.
///
/// # Panics
///
/// Panics if `x` is not homogeneous.
pub fn delta<C: QAlgebra>(f: &Schur<i128>, x: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>> {
    let got = diagonal(&lift_in(x), |_, cells| {
        Ratio::from_poly(plethystic_eval::<C>(f, cells))
    });
    lift_out(got, "delta")
}

/// `Δ'_f F`, \[DM\] (12) — eigenvalue `f[B_μ − 1]`, i.e. `f[·]` over every cell
/// but `(0,0)`. See the module docs on why no virtual alphabet is needed.
///
/// # Panics
///
/// Panics if `x` is not homogeneous.
pub fn delta_prime<C: QAlgebra>(f: &Schur<i128>, x: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>> {
    let got = diagonal(&lift_in(x), |_, cells| {
        Ratio::from_poly(plethystic_eval::<C>(f, &without_corner(cells)))
    });
    lift_out(got, "delta_prime")
}

/// `ΠF`, \[DM\] (21) — eigenvalue `Π_μ`, a polynomial.
///
/// # Panics
///
/// Panics if `x` is not homogeneous.
pub fn big_pi<C: QAlgebra>(x: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>> {
    let got = diagonal(&lift_in(x), |_, cells| {
        <Ratio<C> as Ring>::one().mul_atoms(&pi_atoms(cells))
    });
    lift_out(got, "big_pi")
}

/// `Π⁻¹F`.
///
/// Returns [`Ratio`] coefficients and not polynomials, because it genuinely is
/// not one: `Π⁻¹` divides by `Π_μ` and nothing puts it back. Only the composite
/// [`theta`] is polynomial.
///
/// # Panics
///
/// Panics if `x` is not homogeneous.
pub fn big_pi_inverse<C: QAlgebra>(x: &Schur<QtPoly<C>>) -> Schur<Ratio<C>> {
    diagonal(&lift_in(x), |_, cells| {
        <Ratio<C> as Ring>::one().div_atoms(&pi_atoms(cells))
    })
}

/// `f* = f[X/M]`, \[DM\] (10): `p_k ↦ p_k / ((1−q^k)(1−t^k))`.
///
/// **Linear on the alphabet, not a plethysm on the coefficients** — the
/// distinction `qtkostka.rs` documents at length for `φ_t`, and the reason `f`
/// is restricted to integer coefficients (see [`plethystic_eval`]).
fn star_substitute<C: QAlgebra>(f: &Schur<i128>) -> Schur<Ratio<C>> {
    use crate::convert::ToSchur;
    let lifted: Schur<Ratio<C>> = Schur::from_terms(
        f.terms()
            .iter()
            .map(|(l, &c)| (l.clone(), <Ratio<C> as Ring>::from_i128(c)))
            .collect(),
    );
    let fp: PowerSum<Ratio<C>> = PowerSum::from_schur(&lifted);
    let mut out: PowerSum<Ratio<C>> = PowerSum::zero();
    for (rho, c) in fp.terms() {
        let mut d = Atoms::new();
        for &k in rho.parts() {
            push(&mut d, Atom::unit(k, 0), 1);
            push(&mut d, Atom::unit(0, k), 1);
        }
        let mut v = c.div_atoms(&d);
        v.reduce();
        out.add_term(rho.clone(), v);
    }
    out.to_schur()
}

/// `Θ_f F = Π f* Π⁻¹ F`, \[DM\] (22).
///
/// The only operator here that is not a scalar per μ: `f*` raises the degree by
/// `deg f`, so the middle step is an ordinary Schur product and the second
/// expansion happens one degree band higher than the first.
///
/// The degree-0 cases are \[DM\]'s and are a genuine special case rather than
/// an accident of the formula — the general route would divide by `Π_∅`. When
/// `x` has degree 0 the answer is `f · x` if `f` is constant too, and zero
/// otherwise.
///
/// # Panics
///
/// Panics if `f` or `x` is not homogeneous.
pub fn theta<C: QAlgebra>(f: &Schur<i128>, x: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>> {
    let k = match degree_of(f, "theta's subscript") {
        Some(k) => k,
        None => return Schur::zero(),
    };
    let n = match degree_of(x, "theta") {
        Some(n) => n,
        None => return Schur::zero(),
    };
    if n == 0 {
        // Θ_f F = 0 for deg f ≥ 1, and f·F when both are constants.
        return if k == 0 {
            let c = x.coeff(&Partition::new([]));
            let scale = QtPoly::term(0, 0, C::from_i128(f.coeff(&Partition::new([]))));
            Schur::monomial(Partition::new([]), c.mul(&scale))
        } else {
            Schur::zero()
        };
    }

    let inner = big_pi_inverse(x);
    let product = star_substitute::<C>(f).mul(&inner);
    let htilde = crate::bh::htilde_table::<C>(n + k);
    let mut coeffs = coefficients(&product, n + k, &htilde);
    for (c, (mu, _)) in coeffs.iter_mut().zip(htilde.iter()) {
        *c = c.mul_atoms(&pi_atoms(&coarms(mu)));
        c.reduce();
    }
    lift_out(combine(&coeffs, n + k, &htilde), "theta")
}

// ---------------------------------------------------------------------------
// The closed forms
// ---------------------------------------------------------------------------

/// The `H̃`-coefficients of `e_n` in closed form: `M B_μ Π_μ / w_μ`.
///
/// Verified against the pairing route for n ≤ 6 before being written down, and
/// again by `the_closed_form_agrees_with_the_pairing` below. It skips
/// [`star_against_schur`] entirely, which is the whole cost of the general
/// path.
fn e_coefficients<C: Ring>(n: u32) -> Vec<Ratio<C>> {
    crate::memo::partitions_cached(n)
        .iter()
        .map(|mu| {
            let cells = coarms(mu);
            let mut num = b_mu::<C>(&cells);
            for (&a, &m) in m_atoms().iter().chain(pi_atoms(&cells).iter()) {
                for _ in 0..m {
                    num = a.mul_into(&num);
                }
            }
            let (w, negated) = w_atoms(mu);
            let mut c = Ratio::over(if negated { num.neg() } else { num }, w);
            c.reduce();
            c
        })
        .collect()
}

/// `∇e_n`, by the closed form rather than the general pairing.
///
/// At `n = 0` the answer is `s_∅` with coefficient 1.
///
/// Bounded on [`Ring`] and **not** [`QAlgebra`], unlike [`nabla`]. The closed
/// form never divides by an integer, so this runs over `QtPoly<i128>` where the
/// general path needs ℚ. That is the same trade
/// [`qt_kostka_table_via_bh`](crate::qtkostka) makes.
///
/// `i128` is measurably faster than `Rational` here, and by less than it looks
/// like it should be: `Rational::add_assign` already short-circuits when both
/// operands are integers, which they always are here
/// (`docs/record/macdonald-operators.md`).
///
/// ⚠️ `i128` refuses rather than wraps past its width
/// (`docs/policies/failure.md`, R3). A degree that exceeds it panics, and
/// [`guard`](crate::guard) is the escape hatch that reports and re-runs wide
/// instead.
pub fn nabla_e<C: Ring>(n: u32) -> Schur<QtPoly<C>> {
    closed_form(n, |cells| t_mu::<C>(cells))
}

/// `Δ'_{e_k} e_n`, the Delta conjecture's object, by the closed form.
///
/// At `n = 0` the answer is `s_∅` with coefficient 1, for every `k`.
///
/// Same [`Ring`] bound and the same reason as [`nabla_e`]: computing the
/// eigenvalue never divides.
pub fn delta_prime_e<C: Ring>(k: u32, n: u32) -> Schur<QtPoly<C>> {
    closed_form(n, move |cells| {
        elementary_eval::<C>(&without_corner(cells), k)
    })
}

fn closed_form<C: Ring>(n: u32, eigen: impl Fn(&[(u32, u32)]) -> QtPoly<C>) -> Schur<QtPoly<C>> {
    if n == 0 {
        return Schur::monomial(Partition::new([]), <QtPoly<C> as Ring>::one());
    }
    let htilde = crate::bh::htilde_table::<C>(n);
    let mut coeffs = e_coefficients::<C>(n);
    for (c, (mu, _)) in coeffs.iter_mut().zip(htilde.iter()) {
        *c = c.mul(&Ratio::from_poly(eigen(&coarms(mu))));
    }
    lift_out(combine(&coeffs, n, &htilde), "the closed form for e_n")
}

/// `e_k = s_{1^k}` in the Schur basis, the subscript everything here uses.
pub fn elementary(k: u32) -> Schur<i128> {
    Schur::monomial(Partition::new(std::iter::repeat_n(1, k as usize)), 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;

    type Q = QtPoly<Rational>;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn e_schur(n: u32) -> Schur<Q> {
        Schur::monomial(part(&vec![1; n as usize]), <Q as Ring>::one())
    }

    fn h_schur(n: u32) -> Schur<Q> {
        Schur::monomial(part(&[n]), <Q as Ring>::one())
    }

    fn s_schur(lambda: &[u32]) -> Schur<Q> {
        Schur::monomial(part(lambda), <Q as Ring>::one())
    }

    /// `T_μ = q^{n(μ')} t^{n(μ)}`, with `n(μ) = Σ (i−1)μ_i`.
    ///
    /// The coarm/leg confusion this module's docs warn about is silent
    /// otherwise: arms and legs also give a monomial, also distinct across μ,
    /// and also produce a triangular-looking answer.
    #[test]
    fn t_mu_is_the_hook_monomial() {
        for n in 0..=7u32 {
            for mu in crate::partitions_of(n) {
                let n_mu: u32 = mu
                    .parts()
                    .iter()
                    .enumerate()
                    .map(|(i, &p)| i as u32 * p)
                    .sum();
                let n_conj: u32 = mu
                    .conjugate()
                    .parts()
                    .iter()
                    .enumerate()
                    .map(|(i, &p)| i as u32 * p)
                    .sum();
                let got: Q = t_mu(&coarms(&mu));
                assert_eq!(
                    got,
                    QtPoly::term(n_conj, n_mu, Rational::from_int(1)),
                    "{mu}"
                );
            }
        }
    }

    /// `⟨H̃_μ, H̃_ν⟩_* = δ_{μν} w_μ`, which pins `w_μ`, the star weights and
    /// the modified Macdonald polynomials against each other at once.
    ///
    /// The off-diagonal half is the sharp one: a wrong weight or a wrong `w`
    /// would still make the diagonal look self-consistent.
    #[test]
    fn htilde_is_star_orthogonal_with_norm_w() {
        for n in 0..=6u32 {
            let parts = crate::partitions_of(n);
            let htilde = crate::bh::htilde_table::<Rational>(n);
            for (mu, ht) in &htilde {
                let lifted = lift_in(ht);
                let u = star_against_schur(&lifted, n);
                for (j, (nu, other)) in htilde.iter().enumerate() {
                    let _ = j;
                    let mut acc = <Ratio<Rational> as Ring>::zero();
                    for (k, kappa) in parts.iter().enumerate() {
                        let kt = other.coeff(kappa);
                        if !kt.is_empty() {
                            acc.add_mul(&u[k], &kt);
                        }
                    }
                    acc.reduce();
                    let got = acc.into_poly().expect("the pairing is a polynomial");
                    if mu == nu {
                        let (w, negated) = w_atoms(mu);
                        let mut want = <Q as Ring>::one();
                        for (&a, &m) in &w {
                            for _ in 0..m {
                                want = a.mul_into(&want);
                            }
                        }
                        if negated {
                            want = want.neg();
                        }
                        assert_eq!(got, want, "<H~_{mu},H~_{mu}>_* must be w_{mu}");
                    } else {
                        assert!(got.is_empty(), "<H~_{mu},H~_{nu}>_* must vanish");
                    }
                }
            }
        }
    }

    /// Expanding into `H̃` and summing back is the identity.
    #[test]
    fn the_htilde_expansion_round_trips() {
        for n in 0..=6u32 {
            let htilde = crate::bh::htilde_table::<Rational>(n);
            for lambda in crate::partitions_of(n) {
                let f = s_schur(lambda.parts());
                let lifted = lift_in(&f);
                let coeffs = coefficients(&lifted, n, &htilde);
                let back = lift_out(combine(&coeffs, n, &htilde), "round trip");
                assert_eq!(back, f, "round trip at {lambda}");
            }
        }
    }

    /// `∇e_2 = s_2 + (q+t) s_{11}`, hand-computed.
    #[test]
    fn nabla_of_e_two_is_the_hand_computation() {
        let got = nabla(&e_schur(2));
        assert_eq!(got.coeff(&part(&[2])), <Q as Ring>::one());
        let mut want: Q = QtPoly::term(1, 0, Rational::from_int(1));
        want.add_term(0, 1, Rational::from_int(1));
        assert_eq!(got.coeff(&part(&[1, 1])), want);
        assert_eq!(got.terms().len(), 2, "{got}");
    }

    /// `⟨∇e_n, h_1^n⟩` at q = t = 1 is `dim DH_n = (n+1)^{n−1}`.
    ///
    /// Pairing with `h_1^n` is `Σ_λ c_λ f^λ`, so this reads every coefficient
    /// rather than one slice of one — a wrong `T_μ` cannot survive it.
    #[test]
    fn nabla_e_counts_the_diagonal_harmonics() {
        // Haiman's count (n+1)^{n-1} is stated for n >= 1; nabla e_0 = 1 is
        // covered by the sweeps that start at 0.
        for n in 1..=7u32 {
            let got = nabla_e::<Rational>(n);
            let one = Rational::from_int(1);
            let mut total = Rational::from_int(0);
            for (lambda, c) in got.terms() {
                let d = crate::dimension(lambda).expect("a partition has a dimension");
                let v = c.eval(&one, &one).mul(&Rational::from_u128(d));
                total.add_assign(&v);
            }
            let want = Rational::from_u128((n as u128 + 1).pow(n - 1));
            assert_eq!(total, want, "dim DH_{n}");
        }
    }

    /// `∇e_n` is Schur positive with coefficients in `ℕ[q,t]` — the shuffle
    /// theorem's statement, and not something this code arranges: the answer
    /// arrives over ℚ(q,t) and every denominator has to cancel first.
    #[test]
    fn nabla_e_is_schur_positive() {
        for n in 0..=7u32 {
            for (lambda, c) in nabla_e::<Rational>(n).terms() {
                for (_, v) in c.terms() {
                    assert_eq!(v.denom(), 1, "nabla e_{n} at {lambda} has {v:?}");
                    assert!(v.numer() > 0, "nabla e_{n} at {lambda} has {v:?}");
                }
            }
        }
    }

    /// The closed form and the general pairing must agree.
    ///
    /// They share `combine` and nothing else: one reads `M B_μ Π_μ / w_μ` off
    /// the shape, the other runs `e_n` through the power sums, the character
    /// table and `⟨,⟩_*`.
    #[test]
    fn the_closed_form_agrees_with_the_pairing() {
        for n in 0..=6u32 {
            assert_eq!(nabla_e::<Rational>(n), nabla(&e_schur(n)), "nabla e_{n}");
            for k in 0..n {
                assert_eq!(
                    delta_prime_e::<Rational>(k, n),
                    delta_prime(&elementary(k), &e_schur(n)),
                    "Delta'_e{k} e_{n}"
                );
            }
        }
    }

    /// `Δ_{e_n} = ∇` on degree n, since `e_n[B_μ] = T_μ`.
    #[test]
    fn delta_of_e_n_is_nabla() {
        for n in 0..=5u32 {
            for f in [e_schur(n), h_schur(n)] {
                assert_eq!(delta(&elementary(n), &f), nabla(&f), "degree {n}");
            }
        }
    }

    /// `Δ'_{e_{n−1}} e_n = ∇e_n` — the k = n−1 end of the Delta conjecture,
    /// where it degenerates to the shuffle theorem. Catches an off-by-one in
    /// the cell deletion.
    #[test]
    fn delta_prime_at_the_top_is_nabla() {
        // Delta'_{e_{n-1}} needs e_{n-1} to exist, so the identity starts at n = 1.
        for n in 1..=6u32 {
            assert_eq!(
                delta_prime_e::<Rational>(n - 1, n),
                nabla_e::<Rational>(n),
                "Delta'_e{} e_{n}",
                n - 1
            );
        }
    }

    /// **`Θ_{e_k} ∇ e_{n−k} = Δ'_{e_{n−k−1}} e_n`**, a published theorem
    /// relating all three families.
    ///
    /// The sharpest test available, and the only one that exercises Θ at all:
    /// it runs the product step, both directions of Π, the Δ' cell deletion,
    /// ∇'s monomial, and the expansion at two different degrees. There is no
    /// external oracle for Δ, Δ' or Θ anywhere, so this carries the weight that
    /// a Sage comparison carries elsewhere in the crate.
    #[test]
    fn theta_composed_with_nabla_is_delta_prime() {
        for n in 2..=5u32 {
            for k in 1..n {
                let lhs = theta(&elementary(k), &nabla_e::<Rational>(n - k));
                let rhs = delta_prime_e::<Rational>(n - k - 1, n);
                assert_eq!(lhs, rhs, "Theta_e{k} nabla e_{} vs Delta'", n - k);
            }
        }
    }

    /// Π and Π⁻¹ are inverse, and Π's answer is a polynomial while Π⁻¹'s is
    /// not.
    #[test]
    fn big_pi_and_its_inverse_are_inverse() {
        for n in 0..=5u32 {
            for lambda in crate::partitions_of(n) {
                let f = s_schur(lambda.parts());
                let back = lift_out(
                    diagonal(&big_pi_inverse(&f), |_, cells| {
                        <Ratio<Rational> as Ring>::one().mul_atoms(&pi_atoms(cells))
                    }),
                    "pi round trip",
                );
                assert_eq!(back, f, "Pi Pi^-1 at {lambda}");
                let there = big_pi(&f);
                let and_back = big_pi_inverse(&there);
                assert_eq!(
                    lift_out(and_back, "pi round trip"),
                    f,
                    "Pi^-1 Pi at {lambda}"
                );
            }
        }
    }

    /// ∇ is linear, and `∇^r` agrees with r applications.
    #[test]
    fn nabla_is_linear_and_powers_compose() {
        let f = s_schur(&[3, 1]).add(&s_schur(&[2, 2]).scale(&QtPoly::term(
            0,
            0,
            Rational::from_int(3),
        )));
        let a = nabla(&s_schur(&[3, 1]));
        let b = nabla(&s_schur(&[2, 2]));
        let want = a.add(&b.scale(&QtPoly::term(0, 0, Rational::from_int(3))));
        assert_eq!(nabla(&f), want, "linearity");

        let mut iterated = f.clone();
        for _ in 0..3 {
            iterated = nabla(&iterated);
        }
        assert_eq!(nabla_power(&f, 3), iterated, "nabla^3");
    }

    /// \[QZ\] 2026: `(−1)^{|μ|−ℓ(μ)} ∇^r m_μ` is Schur positive.
    ///
    /// A 2026 theorem about an object no other package computes, checked here
    /// on the range that runs quickly. `m_μ` has to be built by conversion,
    /// which also exercises ∇ on an input that is not a single Schur function.
    #[test]
    fn nabla_of_a_monomial_is_signed_schur_positive() {
        use crate::convert::ToSchur;
        for n in 0..=5u32 {
            for mu in crate::partitions_of(n) {
                let m: crate::Monomial<Q> =
                    crate::Monomial::monomial(mu.clone(), <Q as Ring>::one());
                let f = m.to_schur();
                for r in 1..=2u32 {
                    let got = nabla_power(&f, r);
                    let sign = (mu.size() - mu.len() as u32) % 2 == 1;
                    for (lambda, c) in got.terms() {
                        for (_, v) in c.terms() {
                            let ok = if sign { v.numer() < 0 } else { v.numer() > 0 };
                            assert!(ok, "nabla^{r} m_{mu} at {lambda} has {v:?}");
                        }
                    }
                }
            }
        }
    }

    /// Atom normalization must actually merge the two families, or the
    /// cancellation the whole design rests on never happens.
    #[test]
    fn atoms_normalize_across_the_two_families() {
        // q^a - 1 = -(1 - q^a) and q^0 - t^b = 1 - t^b
        assert_eq!(Atom::diff(3, 0), (Atom::Unit(3, 0), true));
        assert_eq!(Atom::diff(0, 3), (Atom::Unit(0, 3), false));
        assert_eq!(Atom::diff(2, 3), (Atom::Diff(2, 3), false));
        // and the normalized forms really are the polynomials they claim to be
        let (a, neg) = Atom::diff(3, 0);
        let mut want: Q = QtPoly::term(3, 0, Rational::from_int(1));
        want.add_term(0, 0, Rational::from_int(-1));
        let got: Q = a.poly();
        assert_eq!(if neg { got.neg() } else { got }, want);
    }

    /// [`crate::frac::divide_by_diff`] must agree with [`QtPoly::divide_exact`]
    /// everywhere, on multiples **and** on non-multiples.
    ///
    /// The two share no code and no idea — one runs a running sum along
    /// arithmetic progressions of exponents, the other eliminates leading terms
    /// down a monomial order through a `BTreeMap` — so agreement is evidence
    /// rather than tautology. This is the same check `frac.rs` keeps for
    /// `divide_by_factor`, and it is what licenses replacing the general
    /// routine with the specialized one in the hot path.
    ///
    /// The non-divisible half matters as much as the divisible half: a
    /// specialized divider that silently returned a truncated quotient would
    /// pass a round-trip test on its own.
    #[test]
    fn the_specialized_diff_division_agrees_with_the_general_one() {
        let mut dense: Q = QtPoly::zero();
        for a in 0..5 {
            for b in 0..4 {
                dense.add_term(
                    a,
                    b,
                    Rational::from_int(i128::from(a) - 2 * i128::from(b) + 1),
                );
            }
        }
        let mut sparse: Q = QtPoly::zero();
        sparse.add_term(0, 0, Rational::from_int(1));
        sparse.add_term(4, 1, Rational::from_int(-3));
        sparse.add_term(1, 6, Rational::from_int(5));
        sparse.add_term(3, 3, Rational::from_int(2));

        for f in [dense, sparse, <Q as Ring>::one(), QtPoly::zero()] {
            for (a, b) in [(1u32, 1u32), (1, 2), (2, 1), (2, 3), (3, 2), (1, 5), (4, 4)] {
                let d: Q = Atom::Diff(a, b).poly();
                // A genuine multiple: same quotient, not merely both succeeding.
                let prod = f.mul_diff(a, b);
                assert_eq!(
                    crate::frac::divide_by_diff(&prod, a, b),
                    prod.divide_exact(&d),
                    "(q^{a} - t^{b}) dividing its own multiple"
                );
                assert_eq!(
                    crate::frac::divide_by_diff(&prod, a, b),
                    Some(f.clone()),
                    "(q^{a} - t^{b}) must recover the cofactor"
                );
                // And a non-multiple (for the generic f above it is one).
                assert_eq!(
                    crate::frac::divide_by_diff(&f, a, b),
                    f.divide_exact(&d),
                    "(q^{a} - t^{b}) against {f}"
                );
            }
        }
    }

    /// The filter may never reject a genuine multiple.
    ///
    /// It is only a *necessary* condition, so a `true` on a non-multiple is
    /// fine; a `false` on a multiple would silently corrupt every reduction.
    #[test]
    fn the_divisibility_filter_never_rejects_a_multiple() {
        let mut f: Q = QtPoly::zero();
        for a in 0..5 {
            for b in 0..4 {
                f.add_term(a, b, Rational::from_int(i128::from(a) + i128::from(b) + 1));
            }
        }
        for (a, b) in [(1u32, 1u32), (2, 2), (2, 4), (4, 2), (3, 3), (2, 3), (6, 4)] {
            assert!(
                crate::frac::diff_may_divide(&f.mul_diff(a, b), a, b),
                "the filter rejected a multiple of (q^{a} - t^{b})"
            );
        }
    }

    /// Both division paths must be exact and must agree with multiplication.
    #[test]
    fn atom_division_inverts_multiplication() {
        let mut f: Q = QtPoly::zero();
        for a in 0..4 {
            for b in 0..3 {
                f.add_term(a, b, Rational::from_int(i128::from(a) - i128::from(b) + 1));
            }
        }
        for atom in [
            Atom::Unit(1, 0),
            Atom::Unit(0, 1),
            Atom::Unit(2, 3),
            Atom::Diff(1, 1),
            Atom::Diff(3, 2),
        ] {
            let prod = atom.mul_into(&f);
            assert_eq!(atom.divide(&prod), Some(f.clone()), "{atom}");
        }
    }
}
