//! The Orellana–Zabrocki character bases, and a reduced-Kronecker engine.
//!
//! Two inhomogeneous bases live here: [`St`] (`s̃_λ`, the irreducible character
//! basis) and [`Ht`] (`h̃_λ`, the induced trivial character basis). The source
//! is Orellana–Zabrocki, *Symmetric group characters as symmetric functions*,
//! [arXiv:1605.06672](https://arxiv.org/abs/1605.06672); equation and theorem
//! numbers below are that paper's. `docs/record/kronecker.md` has the measured
//! record, including the measurements that chose between routes.
//!
//! The reason to want them is Theorem 7:
//!
//! ```text
//!     s̃_λ · s̃_μ = Σ_{|ν| ≤ |λ|+|μ|} ḡ^ν_{λμ} · s̃_ν
//! ```
//!
//! — an *ordinary* product of symmetric functions whose structure constants are
//! the **reduced (stable) Kronecker coefficients**. Nothing here is ever
//! indexed by a partition of the large n, so this sidesteps the p(n)×p(n)
//! character table that makes [`ops::internal`](crate::ops::internal) run out
//! of memory (1.1 GB at n = 32, per `docs/record/kronecker.md`) long before the
//! interesting cases.
//!
//! ## How it computes
//!
//! Theorem 14 gives the whole engine, once read the right way. The linear map Γ
//! defined by `Γ(s_λ) = s̃_λ` satisfies
//!
//! ```text
//!     Γ(p_γ) = 𝐩_γ = Π_i 𝐩_{i^{m_i(γ)}}                                  (24)
//!     𝐩_{i^r} = Σ_{k=0}^{r} (−1)^{r−k} i^k C(r,k) · (P_i)_k
//!     P_i     = (1/i) Σ_{d|i} μ(i/d) p_d      ⟺      p_i = Σ_{d|i} d·P_d
//! ```
//!
//! where `(x)_k` is the falling factorial *taken in the ring*. The consequence
//! the paper does not spell out, and that this module is built on: `𝐩_{i^r}` is
//! a univariate polynomial in `P_i` of degree r with leading coefficient `i^r`,
//! so the change of basis between monomials in the `P_i` and the `𝐩_γ` is a
//! **tensor product of univariate triangular matrices** — one per part size,
//! and each invertible on its own. Γ⁻¹ therefore costs a per-variable back
//! substitution rather than a p(n)×p(n) matrix inversion.
//!
//! That makes the product cheap in a way the Schur route is not:
//!
//! ```text
//!     s̃_λ · s̃_μ = Γ⁻¹( Γ(s_λ) · Γ(s_μ) ),  read in the Schur basis
//! ```
//!
//! and the multiplication in the middle happens in the **power-sum basis**,
//! where a product is a multiset union of indices. No Littlewood–Richardson
//! coefficient is computed anywhere in a reduced Kronecker calculation.
//!
//! ## What is checked, and against what
//!
//! Sage is an oracle here and was never read — the same clean-room rule as
//! `docs/cleanroom-spec-skew-lr.md`. Before any of this was written, the Γ/Γ⁻¹
//! route was prototyped and compared against Sage's `st` basis: exact agreement
//! on the transitions and on products through `st[4,2]·st[4,2]`, which is 186
//! terms.

// Shape indices. The two structure-constant narrowings check at their sites.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::BTreeMap;
use std::sync::Arc;

use crate::coeff::{QAlgebra, Rational, Ring};
use crate::convert::{FromSchur, ToSchur};
use crate::guard::{guarded, overflow_count, GuardedRat};
use crate::kostka::kostka;
use crate::memo::{
    bold_p_peek, bold_p_store, ht_to_st_cached, reduced_kronecker_cached, schur_to_st_cached,
    st_to_ht_cached, st_to_schur_cached,
};
use crate::partition::Partition;
use crate::sym::{Ht, PowerSum, Schur, St, SymAlgebra, SymFn};

// --- small number theory ----------------------------------------------------

fn divisors(n: u32) -> Vec<u32> {
    (1..=n).filter(|d| n.is_multiple_of(*d)).collect()
}

/// The number-theoretic Möbius function μ(n).
fn mobius(n: u32) -> i128 {
    let mut n = n;
    let mut primes = 0;
    let mut d = 2;
    while d * d <= n {
        if n.is_multiple_of(d) {
            n /= d;
            if n.is_multiple_of(d) {
                return 0; // a squared factor
            }
            primes += 1;
        }
        d += 1;
    }
    if n > 1 {
        primes += 1;
    }
    if primes % 2 == 0 {
        1
    } else {
        -1
    }
}

/// The parts of γ as (part, multiplicity) pairs, ascending by part.
fn multiplicities(gamma: &Partition) -> Vec<(u32, usize)> {
    let mut out: Vec<(u32, usize)> = Vec::new();
    for &part in gamma.parts().iter().rev() {
        match out.last_mut() {
            Some((p, m)) if *p == part => *m += 1,
            _ => out.push((part, 1)),
        }
    }
    out
}

// --- the univariate transition (24), and its inverse ------------------------

/// Coefficients of `𝐩_{i^r}` as a polynomial in `x = P_i`, indexed by power.
///
/// `B[m]` is the coefficient of `x^m`; `B[r] = i^r` is the leading one, which
/// is what makes the family a basis of `ℚ[x]` and the inverse below exist.
///
/// The entries stay small for a reason worth stating, since `i^k` looks
/// alarming: r is only ever a multiplicity `m_i(γ)` of a part i in a partition
/// of n, so `i·r ≤ n`. That makes `i^r` maximized around i = 2, r = n/2, where
/// it is ~10³ at n = 20 rather than 20²⁰. Callers that do not respect
/// `i·r ≤ n` are outside the
/// contract.
fn bold_uni(i: u32, r: usize) -> Vec<i128> {
    // Falling factorials (x)_k as coefficient vectors, built up by multiplying
    // by (x − (k−1)).
    let mut falls: Vec<Vec<i128>> = vec![vec![1]];
    for k in 1..=r {
        let prev = &falls[k - 1];
        let mut next = vec![0i128; k + 1];
        for (m, &c) in prev.iter().enumerate() {
            next[m + 1] += c; // times x
            next[m] -= c * (k as i128 - 1); // times −(k−1)
        }
        falls.push(next);
    }

    let mut out = vec![0i128; r + 1];
    let mut binom = 1i128; // C(r, k)
    for k in 0..=r {
        let sign = if (r - k).is_multiple_of(2) { 1 } else { -1 };
        let scale = sign * (i as i128).pow(k as u32) * binom;
        for (m, &c) in falls[k].iter().enumerate() {
            out[m] += scale * c;
        }
        binom = binom * (r as i128 - k as i128) / (k as i128 + 1);
    }
    out
}

/// A rational coefficient ring the engine can run over, and can recognize
/// integers in on the way out.
///
/// Two instantiations that matter, and the pair is the whole overflow story:
/// [`GuardedRat`], which is fixed-width and *reports* leaving it, and (under
/// the `bignum` feature) `BigRational`, which cannot leave it. The engine runs
/// over the first and re-runs over the second when the first says it lost.
trait RatLike: QAlgebra {
    /// This value as an integer, or `None` if it is not one.
    fn as_integer(&self) -> Option<i128>;
}

impl RatLike for Rational {
    fn as_integer(&self) -> Option<i128> {
        (self.denom() == 1).then(|| self.numer())
    }
}

impl RatLike for GuardedRat {
    fn as_integer(&self) -> Option<i128> {
        (self.denom() == 1).then(|| self.numer())
    }
}

#[cfg(feature = "bignum")]
impl RatLike for num_rational::BigRational {
    fn as_integer(&self) -> Option<i128> {
        use num_traits::ToPrimitive;
        self.is_integer()
            .then(|| self.to_integer().to_i128())
            .flatten()
    }
}

/// The row `x^m = Σ_r N[r] · 𝐩_{i^r}`, by back substitution against
/// [`bold_uni`].
///
/// This is the piece that makes Γ⁻¹ cheap: the full change of basis is a tensor
/// product over part sizes, so inverting it never means inverting anything
/// bigger than an (m+1)×(m+1) triangular matrix.
fn bold_uni_inverse<R: RatLike>(i: u32, m: usize) -> Vec<R> {
    let rows: Vec<Vec<i128>> = (0..=m).map(|r| bold_uni(i, r)).collect();
    let mut coeffs = vec![R::zero(); m + 1];
    // Remainder in the monomial basis, starting at x^m.
    let mut rem = vec![R::zero(); m + 1];
    rem[m] = R::one();
    for r in (0..=m).rev() {
        if rem[r].is_zero() {
            continue;
        }
        // The diagonal entry is rows[r][r] = i^r, so the division is by a
        // positive integer — `QAlgebra`, never `Field`.
        let c = rem[r].div_u128((i as u128).pow(r as u32));
        coeffs[r] = c.clone();
        for (j, &b) in rows[r].iter().enumerate().take(r + 1) {
            let term = c.mul(&R::from_i128(b));
            rem[j].sub_assign(&term);
        }
    }
    coeffs
}

/// `P_i = (1/i) Σ_{d|i} μ(i/d) p_d`, in the power-sum basis.
fn p_variable<R: RatLike>(i: u32) -> PowerSum<R> {
    let mut out = PowerSum::zero();
    for d in divisors(i) {
        let m = mobius(i / d);
        if m != 0 {
            out.add_term(Partition::new([d]), R::from_i128(m).div_u128(i as u128));
        }
    }
    out
}

/// `𝐩_γ` in the power-sum basis (OZ Eq 24).
fn bold_p_in<R: RatLike>(gamma: &Partition) -> PowerSum<R> {
    let mut acc: PowerSum<R> = PowerSum::unit();
    for (i, r) in multiplicities(gamma) {
        let coeffs = bold_uni(i, r);
        let pi = p_variable::<R>(i);
        // Powers of P_i, accumulated rather than recomputed.
        let mut power: PowerSum<R> = PowerSum::unit();
        let mut factor = PowerSum::zero();
        for (m, &c) in coeffs.iter().enumerate() {
            if m > 0 {
                power = power.mul(&pi);
            }
            if c != 0 {
                factor = factor.add(&power.scale(&R::from_i128(c)));
            }
        }
        acc = acc.mul(&factor);
    }
    acc
}

/// `𝐩_γ` over [`GuardedRat`], memoized — but **only when the computation stayed
/// inside the fixed width**.
///
/// A value produced by a call that overflowed is garbage, and caching it would
/// convert a detected overflow into an undetected one: the counter is compared
/// around the call that computes a value, so a later reader of a poisoned entry
/// would see a clean counter and accept the answer. Recomputing is the cheap
/// side of that trade — past the wall every call is about to escalate anyway.
fn bold_guarded(gamma: &Partition) -> PowerSum<GuardedRat> {
    if let Some(v) = bold_p_peek(gamma) {
        return (*v).clone();
    }
    let before = overflow_count();
    let value = bold_p_in::<GuardedRat>(gamma);
    if overflow_count() == before {
        bold_p_store(gamma, value.clone());
    }
    value
}

// --- Γ and Γ⁻¹ ---------------------------------------------------------------

/// Apply Γ to a power-sum element: `p_γ ↦ 𝐩_γ` (OZ Thm 14).
///
/// Because `Γ(s_λ) = s̃_λ`, this is how an `st` element becomes an ordinary
/// one: write it over the Schur basis formally, convert to power sums, apply
/// this.
fn gamma<R: RatLike>(f: &PowerSum<R>, bold: fn(&Partition) -> PowerSum<R>) -> PowerSum<R> {
    let mut out = PowerSum::zero();
    for (g, c) in f.terms() {
        crate::interrupt::poll();
        for (idx, v) in bold(g).terms() {
            out.add_term(idx.clone(), v.mul(c));
        }
    }
    out
}

/// Apply Γ⁻¹ to a power-sum element.
///
/// Two steps, both of them tensor products over part sizes:
///
/// 1. `p → P`-monomials, by `p_i = Σ_{d|i} d·P_d` (the Möbius inverse of the
///    definition of `P_i`, so this direction has no denominators at all);
/// 2. `P`-monomials → the `𝐩` basis, by [`bold_uni_inverse`] one variable at a
///    time.
///
/// The result is read back with its indices interpreted as `p_γ` rather than
/// `𝐩_γ`, which is exactly Γ⁻¹.
fn gamma_inverse<R: RatLike>(f: &PowerSum<R>) -> PowerSum<R> {
    // Step 1: expand into monomials in the P variables.
    let mut monomials: BTreeMap<Partition, R> = BTreeMap::new();
    for (g, c) in f.terms() {
        crate::interrupt::poll();
        let mut acc: BTreeMap<Vec<u32>, R> = BTreeMap::new();
        acc.insert(Vec::new(), c.clone());
        for &part in g.parts() {
            let mut next: BTreeMap<Vec<u32>, R> = BTreeMap::new();
            for (mono, cc) in &acc {
                for d in divisors(part) {
                    let mut key = mono.clone();
                    key.push(d);
                    key.sort_unstable_by(|a, b| b.cmp(a));
                    let v = cc.mul(&R::from_i128(d as i128));
                    next.entry(key)
                        .and_modify(|e| e.add_assign(&v))
                        .or_insert(v);
                }
            }
            acc = next;
        }
        for (mono, v) in acc {
            if v.is_zero() {
                continue;
            }
            let key = Partition::new(mono);
            monomials
                .entry(key)
                .and_modify(|e| e.add_assign(&v))
                .or_insert(v);
        }
    }

    // Step 2: each P-monomial into the 𝐩 basis, one variable at a time.
    let mut out = PowerSum::zero();
    for (mono, c) in monomials {
        crate::interrupt::poll();
        if c.is_zero() {
            continue;
        }
        let mut acc: BTreeMap<Vec<u32>, R> = BTreeMap::new();
        acc.insert(Vec::new(), c);
        for (i, m) in multiplicities(&mono) {
            let inv = bold_uni_inverse::<R>(i, m);
            let mut next: BTreeMap<Vec<u32>, R> = BTreeMap::new();
            for (idx, cc) in &acc {
                for (r, n) in inv.iter().enumerate() {
                    if n.is_zero() {
                        continue;
                    }
                    let mut key = idx.clone();
                    key.extend(std::iter::repeat_n(i, r));
                    key.sort_unstable_by(|a, b| b.cmp(a));
                    let v = cc.mul(n);
                    next.entry(key)
                        .and_modify(|e| e.add_assign(&v))
                        .or_insert(v);
                }
            }
            acc = next;
        }
        for (idx, v) in acc {
            out.add_term(Partition::new(idx), v);
        }
    }
    out
}

/// A Schur element over ℚ whose coefficients are known to be integers, as
/// integers.
///
/// Both transitions are integral — Γ⁻¹ by OZ Thm 1(2), where the coefficients
/// are multiplicities in the restriction of a GLₙ-module, and Γ because it
/// inverts one. The rationals are an artefact of routing through power sums,
/// which divides by `z_γ`. Panicking on a surviving denominator is the same
/// discipline `Rat::into_poly` applies on the Macdonald side: the theorem is
/// enforced where it is used, not assumed.
fn integral_row<R: RatLike>(f: &Schur<R>) -> Option<Vec<(Partition, i128)>> {
    f.terms()
        .iter()
        .map(|(p, c)| c.as_integer().map(|n| (p.clone(), n)))
        .collect()
}

/// `s̃_λ` in the power-sum basis.
fn st_in_power_sum<R: RatLike>(
    lambda: &Partition,
    bold: fn(&Partition) -> PowerSum<R>,
) -> PowerSum<R> {
    let s: Schur<R> = Schur::monomial(lambda.clone(), R::one());
    gamma(&PowerSum::from_schur(&s), bold)
}

/// Run a row computation over fixed-width coefficients, and re-run it exactly
/// if anything left the width.
///
/// The measured wall is `|λ|+|μ| = 24`: `st[8,5]·st[7,4]` completes and
/// `st[8,5]·st[8,5]` does not. The overflow is entirely in the **intermediate**
/// rationals — the answers there are 16-bit — because routing through power
/// sums divides by `z_γ`, and `z_γ` alone reaches 10²⁶ by degree 26, before the
/// two sides are multiplied together.
///
/// Without the `bignum` feature this panics rather than returning something
/// wrong, which is the only acceptable behavior: an intermediate that left the
/// width has no exact continuation here. A wrapped intermediate can land on a
/// denominator of 1 and be accepted as an integer answer. That is the failure
/// mode [`guarded`] exists to remove.
fn escalating(
    what: &str,
    fast: impl FnOnce() -> Option<Vec<(Partition, i128)>>,
    exact: impl FnOnce() -> Option<Vec<(Partition, i128)>>,
) -> Vec<(Partition, i128)> {
    // Two ways the fast path can lose, and both must be non-fatal: the guard
    // counter moves, or a wrapped intermediate lands on a value that is not an
    // integer. Panicking on the second — which an earlier version did, inside
    // the closure — makes the escalation unreachable, since the panic escapes
    // `guarded` before it can report `None`.
    if let Some(Some(v)) = guarded(fast) {
        return v;
    }
    #[cfg(feature = "bignum")]
    {
        let _ = what;
        exact().unwrap_or_else(|| {
            panic!(
                "{what}: a coefficient is not an integer over BigRational, which is a bug \
                    rather than an overflow"
            )
        })
    }
    #[cfg(not(feature = "bignum"))]
    {
        let _ = exact;
        panic!(
            "{what}: an intermediate coefficient left i128. The answer itself is \
             almost certainly small — this is the z_γ in the power-sum route, not \
             the reduced Kronecker coefficients. Rebuild with --features bignum."
        )
    }
}

/// `s̃_λ` in the Schur basis, memoized.
fn st_to_schur_row(lambda: &Partition) -> Arc<Vec<(Partition, i128)>> {
    st_to_schur_cached(lambda, || {
        escalating(
            &format!("s̃_{lambda} → s"),
            || integral_row(&st_in_power_sum(lambda, bold_guarded).to_schur()),
            || {
                #[cfg(feature = "bignum")]
                {
                    integral_row(
                        &st_in_power_sum::<num_rational::BigRational>(lambda, bold_p_in).to_schur(),
                    )
                }
                #[cfg(not(feature = "bignum"))]
                None
            },
        )
    })
}

/// `s_ν` in the `s̃` basis, memoized — the coefficients `r_{νμ}` of OZ Thm
/// 1(2).
fn schur_to_st_row(nu: &Partition) -> Arc<Vec<(Partition, i128)>> {
    schur_to_st_cached(nu, || {
        fn run<R: RatLike>(nu: &Partition) -> Option<Vec<(Partition, i128)>> {
            let s: Schur<R> = Schur::monomial(nu.clone(), R::one());
            integral_row(&gamma_inverse(&PowerSum::from_schur(&s)).to_schur())
        }
        escalating(
            &format!("s_{nu} → s̃"),
            || run::<GuardedRat>(nu),
            || {
                #[cfg(feature = "bignum")]
                {
                    run::<num_rational::BigRational>(nu)
                }
                #[cfg(not(feature = "bignum"))]
                None
            },
        )
    })
}

/// One whole column of reduced Kronecker coefficients, `s̃_λ · s̃_μ`, memoized.
///
/// The multiplication is a multiset union of power-sum indices; the two Γ's
/// either side are what cost anything.
fn reduced_kronecker_row(lambda: &Partition, mu: &Partition) -> Arc<Vec<(Partition, i128)>> {
    reduced_kronecker_cached(lambda, mu, || {
        fn run<R: RatLike>(
            lambda: &Partition,
            mu: &Partition,
            bold: fn(&Partition) -> PowerSum<R>,
        ) -> Option<Vec<(Partition, i128)>> {
            let prod = st_in_power_sum(lambda, bold).mul(&st_in_power_sum(mu, bold));
            integral_row(&gamma_inverse(&prod).to_schur())
        }
        escalating(
            &format!("s̃_{lambda} · s̃_{mu}"),
            || run(lambda, mu, bold_guarded),
            || {
                #[cfg(feature = "bignum")]
                {
                    run::<num_rational::BigRational>(lambda, mu, bold_p_in)
                }
                #[cfg(not(feature = "bignum"))]
                None
            },
        )
    })
}

// --- the public surface ------------------------------------------------------

impl<C: Ring> ToSchur<C> for St<C> {
    fn to_schur(&self) -> Schur<C> {
        let mut out = Schur::zero();
        for (lambda, c) in self.terms() {
            for (nu, k) in st_to_schur_row(lambda).iter() {
                out.add_term(nu.clone(), C::from_i128(*k).mul(c));
            }
        }
        out
    }
}

impl<C: Ring> FromSchur<C> for St<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        let mut out = St::zero();
        for (nu, c) in s.terms() {
            for (lambda, k) in schur_to_st_row(nu).iter() {
                out.add_term(lambda.clone(), C::from_i128(*k).mul(c));
            }
        }
        out
    }
}

impl<C: Ring> St<C> {
    /// The product `s̃_λ · s̃_μ`, whose structure constants are the reduced
    /// Kronecker coefficients (OZ Thm 7).
    ///
    /// Bilinear over one memoized column per (λ, μ), the same shape as
    /// [`Schur::mul_with`] over an LR backend — so a memoized column is shared
    /// by every element that mentions the pair, and the coefficient ring never
    /// enters the engine.
    ///
    /// # Panics
    ///
    /// Panics where [`reduced_kronecker_product`] does, whose columns this
    /// reads.
    pub fn mul(&self, other: &Self) -> Self {
        let mut out = Self::zero();
        for (lambda, cl) in self.terms() {
            for (mu, cm) in other.terms() {
                let scale = cl.mul(cm);
                for (nu, k) in reduced_kronecker_row(lambda, mu).iter() {
                    out.add_term(nu.clone(), C::from_i128(*k).mul(&scale));
                }
            }
        }
        out
    }
}

impl<C: Ring> SymAlgebra<C> for St<C> {
    fn times(&self, o: &Self) -> Self {
        self.mul(o)
    }
}

/// The reduced (stable) Kronecker product `s̃_λ · s̃_μ = Σ_ν ḡ^ν_{λμ} s̃_ν`.
///
/// The whole column at once, because that is the engine's unit of work.
/// `s̃_∅` is the unit, so an empty λ returns `s̃_μ`.
///
/// # Panics
///
/// Panics without the `bignum` feature, past the measured wall at `|λ|+|μ| =
/// 24`: `s̃_{(8,5)}·s̃_{(7,4)}` completes and `s̃_{(8,5)}·s̃_{(8,5)}` does not.
/// The wall is `z_γ` in the intermediate rationals, not the answers, which stay
/// under 20 bits — a fact about `i128`, and with `bignum` the same call
/// escalates and returns exactly (`docs/record/kronecker.md`).
pub fn reduced_kronecker_product<C: Ring>(lambda: &Partition, mu: &Partition) -> St<C> {
    let mut out = St::zero();
    for (nu, k) in reduced_kronecker_row(lambda, mu).iter() {
        out.add_term(nu.clone(), C::from_i128(*k));
    }
    out
}

/// A single reduced Kronecker coefficient `ḡ^ν_{λμ}`.
///
/// Returns zero for a ν the column does not mention.
///
/// Convenience over [`reduced_kronecker_product`], and honest about it: asking
/// for one coefficient costs what the whole column costs, exactly as
/// [`ops::kronecker`](crate::ops::kronecker) does for the unreduced case.
///
/// # Panics
///
/// Panics where [`reduced_kronecker_product`] does, whose column this reads.
pub fn reduced_kronecker<C: Ring>(lambda: &Partition, mu: &Partition, nu: &Partition) -> C {
    reduced_kronecker_row(lambda, mu)
        .iter()
        .find(|(p, _)| p == nu)
        .map(|(_, k)| C::from_i128(*k))
        .unwrap_or_else(C::zero)
}

// --- the induced trivial basis, and the second route ------------------------

/// Partitions γ ⊢ |λ| + k with γ/λ a horizontal strip.
///
/// The interlacing condition γ₁ ≥ λ₁ ≥ γ₂ ≥ λ₂ ≥ … is what "at most one cell
/// per column" comes to, and it bounds ℓ(γ) by ℓ(λ)+1.
fn horizontal_strips(lambda: &Partition, k: u32) -> Vec<Partition> {
    let l = lambda.parts();
    let mut out = Vec::new();
    let mut current: Vec<u32> = Vec::with_capacity(l.len() + 1);

    fn rec(i: usize, left: u32, l: &[u32], current: &mut Vec<u32>, out: &mut Vec<Partition>) {
        if i == l.len() + 1 {
            if left == 0 {
                out.push(Partition::new(current.clone()));
            }
            return;
        }
        // γ_i ranges over [λ_i, λ_{i−1}] (no upper bound for i = 0).
        let lo = if i < l.len() { l[i] } else { 0 };
        let hi = if i == 0 {
            lo + left
        } else {
            l[i - 1].min(lo + left)
        };
        for v in lo..=hi {
            current.push(v);
            rec(i + 1, left - (v - lo), l, current, out);
            current.pop();
        }
    }

    rec(0, k, l, &mut current, &mut out);
    out
}

/// The coefficient of `s̃_λ` in `h̃_μ` (OZ Eq 7).
///
/// The paper states this with *stable* Kostka numbers
/// `K_{(n−|λ|,λ)(n−|μ|,μ)}` for n ≥ 2|μ|, and then gives the n-free form used
/// here: sum the ordinary `K_{γμ}` over γ ⊢ |μ| with γ/λ a horizontal strip.
/// Same number, and no partition of 2|μ| is ever enumerated to get it.
fn ht_to_st_coeff(lambda: &Partition, mu: &Partition) -> u128 {
    if lambda.size() > mu.size() {
        return 0;
    }
    horizontal_strips(lambda, mu.size() - lambda.size())
        .iter()
        .map(|g| kostka(g, mu))
        .sum()
}

/// `h̃_μ` expanded in the `s̃` basis.
fn ht_to_st_row(mu: &Partition) -> Arc<Vec<(Partition, i128)>> {
    ht_to_st_cached(mu, || {
        let mut out = Vec::new();
        for size in 0..=mu.size() {
            for lambda in crate::memo::partitions_cached(size).iter() {
                let c = ht_to_st_coeff(lambda, mu);
                if c != 0 {
                    // A structure constant crossing into the signed ring, so it
                    // checks rather than proving: `ht_to_st_coeff` sums ordinary
                    // Kostka numbers, which have no a-priori `i128` bound here.
                    let c = i128::try_from(c).unwrap_or_else(|_| {
                        panic!("the h̃ → s̃ coefficient at ({lambda}, {mu}) does not fit i128")
                    });
                    out.push((lambda.clone(), c));
                }
            }
        }
        out
    })
}

/// Order on partitions in which the `h̃ → s̃` matrix is triangular with unit
/// diagonal: by size, then by **reverse** lexicographic order.
///
/// Both halves are forced. `ht_to_st_coeff` vanishes unless |λ| ≤ |μ|, and at
/// equal size it is the Kostka number `K_{λμ}`, which vanishes unless λ ⊵ μ —
/// and dominance refines to lex. So every off-diagonal entry sits strictly
/// earlier in this order than its column, which is what back substitution
/// needs.
fn triangular_key(p: &Partition) -> (u32, std::cmp::Reverse<Vec<u32>>) {
    (p.size(), std::cmp::Reverse(p.parts().to_vec()))
}

/// `s̃_λ` expanded in the `h̃` basis, by back substitution against
/// [`ht_to_st_row`].
fn st_to_ht_row(lambda: &Partition) -> Arc<Vec<(Partition, i128)>> {
    st_to_ht_cached(lambda, || {
        let mut rem: BTreeMap<Partition, i128> = BTreeMap::new();
        rem.insert(lambda.clone(), 1);
        let mut out: Vec<(Partition, i128)> = Vec::new();

        // One lookup yields both the pivot and its coefficient, so neither the
        // emptiness test nor the re-fetch can disagree with it.
        while let Some((pivot, c)) = rem
            .iter()
            .max_by_key(|(p, _)| triangular_key(p))
            .map(|(p, &c)| (p.clone(), c))
        {
            rem.remove(&pivot);
            if c == 0 {
                continue;
            }
            out.push((pivot.clone(), c));
            for (nu, k) in ht_to_st_row(&pivot).iter() {
                if *nu == pivot {
                    continue; // the unit diagonal, already consumed
                }
                *rem.entry(nu.clone()).or_insert(0) -= c * k;
            }
            rem.retain(|_, v| *v != 0);
        }
        out
    })
}

impl<C: Ring> ToSchur<C> for Ht<C> {
    fn to_schur(&self) -> Schur<C> {
        let mut st: St<C> = St::zero();
        for (mu, c) in self.terms() {
            for (lambda, k) in ht_to_st_row(mu).iter() {
                st.add_term(lambda.clone(), C::from_i128(*k).mul(c));
            }
        }
        st.to_schur()
    }
}

impl<C: Ring> FromSchur<C> for Ht<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        let st: St<C> = St::from_schur(s);
        let mut out = Ht::zero();
        for (lambda, c) in st.terms() {
            for (mu, k) in st_to_ht_row(lambda).iter() {
                out.add_term(mu.clone(), C::from_i128(*k).mul(c));
            }
        }
        out
    }
}

/// How many matrices [`ht_product_terms`] will look at before giving up.
///
/// The enumeration is cheap exactly where the partitions are short — λ = μ =
/// (6,4) is a 2×2 free block, at most 1225 matrices — and hopeless where they
/// are long: λ = μ = (1¹⁰) is a 10×10 block with row sums 1, which is 11¹⁰.
/// This route is a cross-check, so refusing is a perfectly good answer.
const HT_PRODUCT_BUDGET: u64 = 4_000_000;

/// `h̃_λ · h̃_μ` by the matrix rule, or `None` if the enumeration is too large.
///
/// The terms come back in ascending [`Partition`] order. Every multiplicity is
/// positive.
///
/// **Derived rather than read**, and the derivation is one paragraph. `h̃_λ` is
/// the character of the permutation module `M^{(n−|λ|,λ)}`. A tensor product of
/// permutation modules is the permutation module on the product of the two
/// coset spaces; its orbits are the double cosets, indexed by non-negative
/// integer
/// matrices with row sums `(n−|λ|, λ₁, λ₂, …)` and column sums
/// `(n−|μ|, μ₁, μ₂, …)`; and the stabilizer of an orbit is the Young subgroup on
/// the entries. So
///
/// ```text
///     h̃_λ · h̃_μ = Σ_A h̃_{entries of A, minus the (0,0) corner}
/// ```
///
/// Only the corner grows with n, so dropping it leaves a finite sum that never
/// mentions n. Checked against Sage on 9 pairs before it was implemented,
/// including `ht[3,1]·ht[2,2]`.
pub fn ht_product_terms(lambda: &Partition, mu: &Partition) -> Option<Vec<(Partition, u128)>> {
    let rows = lambda.parts();
    let cols = mu.parts();
    let mut acc: BTreeMap<Partition, u128> = BTreeMap::new();
    let mut cap: Vec<u32> = cols.to_vec();
    let mut entries: Vec<u32> = Vec::new();
    let mut scratch: Vec<u32> = Vec::new();
    let mut budget = HT_PRODUCT_BUDGET;

    // Choose one row of the free block at a time; `cap` carries the column
    // capacity left, and `left` the row capacity left.
    #[allow(clippy::too_many_arguments)]
    fn row(
        i: usize,
        j: usize,
        left: u32,
        rows: &[u32],
        cap: &mut Vec<u32>,
        entries: &mut Vec<u32>,
        scratch: &mut Vec<u32>,
        acc: &mut BTreeMap<Partition, u128>,
        budget: &mut u64,
    ) -> bool {
        if *budget == 0 {
            return false;
        }
        if j == cap.len() {
            entries.push(left); // the row slack, A[i][0]
            let ok = block(i + 1, rows, cap, entries, scratch, acc, budget);
            entries.pop();
            return ok;
        }
        for v in 0..=left.min(cap[j]) {
            cap[j] -= v;
            entries.push(v);
            let ok = row(i, j + 1, left - v, rows, cap, entries, scratch, acc, budget);
            entries.pop();
            cap[j] += v;
            if !ok {
                return false;
            }
        }
        true
    }

    fn block(
        i: usize,
        rows: &[u32],
        cap: &mut Vec<u32>,
        entries: &mut Vec<u32>,
        scratch: &mut Vec<u32>,
        acc: &mut BTreeMap<Partition, u128>,
        budget: &mut u64,
    ) -> bool {
        if i == rows.len() {
            if *budget == 0 {
                return false;
            }
            *budget -= 1;
            // The column slacks A[0][j] complete the matrix.
            //
            // Built in a reused buffer, and sorted and stripped of zeros in
            // place, so the leaf allocates exactly once — for the key the map
            // has to own. It used to allocate twice, cloning `entries` and then
            // letting `Partition::new` collect a second vector out of the
            // filter, on every one of up to `HT_PRODUCT_BUDGET` leaves.
            scratch.clear();
            scratch.extend_from_slice(entries);
            scratch.extend_from_slice(cap);
            scratch.retain(|&x| x != 0);
            scratch.sort_unstable_by(|a, b| b.cmp(a));
            *acc.entry(Partition::from_sorted(scratch.clone()))
                .or_insert(0) += 1;
            return true;
        }
        row(i, 0, rows[i], rows, cap, entries, scratch, acc, budget)
    }

    if !block(
        0,
        rows,
        &mut cap,
        &mut entries,
        &mut scratch,
        &mut acc,
        &mut budget,
    ) {
        return None;
    }
    Some(acc.into_iter().collect())
}

impl<C: Ring> Ht<C> {
    /// The product in the `h̃` basis, by the matrix rule.
    ///
    /// # Panics
    ///
    /// Panics if the enumeration exceeds its budget. This is a capacity wall,
    /// not a violated precondition: [`ht_product_terms`] returns `None` on the
    /// same input, and a caller who needs to handle the case rather than die on
    /// it should go through that.
    pub fn mul(&self, other: &Self) -> Self {
        let mut out = Self::zero();
        for (lambda, cl) in self.terms() {
            for (mu, cm) in other.terms() {
                let scale = cl.mul(cm);
                let terms = ht_product_terms(lambda, mu)
                    .expect("ht product enumeration exceeded its budget");
                for (nu, k) in terms {
                    out.add_term(nu, C::from_u128(k).mul(&scale));
                }
            }
        }
        out
    }
}

impl<C: Ring> SymAlgebra<C> for Ht<C> {
    fn times(&self, o: &Self) -> Self {
        self.mul(o)
    }
}

/// `s̃_λ · s̃_μ` computed the other way: through the `h̃` basis and its matrix
/// rule, sharing no code with [`reduced_kronecker_product`].
///
/// `None` when the matrix enumeration would be too large (long partitions).
/// This exists to be disagreed with — it is the crate's standing pattern of
/// holding a fast engine to an independent one. Here it matters more than
/// usual: past `st[4,3]·st[4,3]` there is no third-party package left to ask
/// (`docs/record/kronecker.md` surveys them).
///
/// # Panics
///
/// Panics if an `h̃` product multiplicity exceeds `i128`. It counts double
/// cosets and is unbounded in principle. The narrowing therefore checks as it
/// crosses into the signed ring rather than being absorbed into the sum (R5).
/// The enumeration budget is the *other* limit and is not a panic — that is
/// the `None`.
pub fn reduced_kronecker_via_ht<C: Ring>(lambda: &Partition, mu: &Partition) -> Option<St<C>> {
    let mut acc: BTreeMap<Partition, i128> = BTreeMap::new();
    for (a, ca) in st_to_ht_row(lambda).iter() {
        for (b, cb) in st_to_ht_row(mu).iter() {
            for (nu, k) in ht_product_terms(a, b)? {
                // As above: `k` counts double cosets and is unbounded in
                // principle, so the narrowing is checked at the seam rather
                // than absorbed into the sum.
                let k = i128::try_from(k).unwrap_or_else(|_| {
                    panic!("the h̃ product multiplicity at {nu} does not fit i128")
                });
                *acc.entry(nu).or_insert(0) += ca * cb * k;
            }
        }
    }
    let mut out = St::zero();
    for (nu, c) in acc {
        if c == 0 {
            continue;
        }
        for (lam, k) in ht_to_st_row(&nu).iter() {
            out.add_term(lam.clone(), C::from_i128(c * k));
        }
    }
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sym::{Elementary, Homogeneous};

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn st(v: &[u32]) -> St<i128> {
        St::monomial(part(v), 1)
    }

    /// OZ Thm 1(3): `s̃_{1^r} = Σ_{i=0}^{r} (−1)^i e_{r−i}`.
    ///
    /// Hand-checkable, and it pins the *embedding* rather than any structure
    /// constant — the paper uses exactly this as the initial condition that
    /// makes the basis unique, so getting it right is not implied by getting
    /// the product right.
    #[test]
    fn column_shapes_are_alternating_sums_of_elementaries() {
        for r in 1..=6u32 {
            let col = st(&vec![1; r as usize]);
            let e: Elementary<i128> = Elementary::from_schur(&col.to_schur());
            for i in 0..=r {
                let want = if i % 2 == 0 { 1 } else { -1 };
                assert_eq!(e.coeff(&part(&[r - i])), want, "e_{} in s̃_(1^{r})", r - i);
            }
            assert_eq!(e.terms().len(), r as usize + 1);
        }
    }

    /// `s̃_1 = s_1 − 1`, the smallest non-trivial case, and the one where the
    /// inhomogeneity first shows up.
    #[test]
    fn first_character_is_the_defect_of_the_trivial_one() {
        let s = st(&[1]).to_schur();
        assert_eq!(s.coeff(&part(&[1])), 1);
        assert_eq!(s.coeff(&Partition::default()), -1);
        assert_eq!(s.terms().len(), 2);
    }

    /// OZ Eq (20): `h_{21} = s̃_3 + s̃_{21} + 4s̃_2 + 3s̃_{11} + 7s̃_1 +
    /// 4s̃_∅`.
    ///
    /// Printed in the paper, so this is a check against a *published* value
    /// rather than against another routine in this crate.
    #[test]
    fn published_expansion_of_h_21() {
        let h: Homogeneous<i128> = Homogeneous::monomial(part(&[2, 1]), 1);
        let got: St<i128> = St::from_schur(&h.to_schur());
        for (p, want) in [
            (vec![3], 1),
            (vec![2, 1], 1),
            (vec![2], 4),
            (vec![1, 1], 3),
            (vec![1], 7),
            (vec![], 4),
        ] {
            assert_eq!(got.coeff(&part(&p)), want, "coefficient of s̃_{p:?}");
        }
        assert_eq!(got.terms().len(), 6);
    }

    /// OZ Eq (21): the decomposition of `V^{⊗4}` as an Sₙ-module, also printed
    /// in the paper.
    #[test]
    fn published_expansion_of_h_1111() {
        let h: Homogeneous<i128> = Homogeneous::monomial(part(&[1, 1, 1, 1]), 1);
        let got: St<i128> = St::from_schur(&h.to_schur());
        for (p, want) in [
            (vec![], 15),
            (vec![1], 37),
            (vec![1, 1], 31),
            (vec![1, 1, 1], 10),
            (vec![1, 1, 1, 1], 1),
            (vec![2], 31),
            (vec![2, 1], 20),
            (vec![2, 1, 1], 3),
            (vec![2, 2], 2),
            (vec![3], 10),
            (vec![3, 1], 3),
            (vec![4], 1),
        ] {
            assert_eq!(got.coeff(&part(&p)), want, "coefficient of s̃_{p:?}");
        }
        assert_eq!(got.terms().len(), 12);
    }

    /// The two transitions invert each other, on every partition up to size 8.
    #[test]
    fn st_and_schur_round_trip() {
        for n in 0..=8u32 {
            for lambda in crate::partition::partitions_of(n) {
                let x = St::<i128>::monomial(lambda.clone(), 1);
                let back: St<i128> = St::from_schur(&x.to_schur());
                assert_eq!(back, x, "round trip at {lambda}");
            }
        }
    }

    /// Same for the `h̃` basis, which goes through `s̃` in both directions.
    #[test]
    fn ht_and_schur_round_trip() {
        for n in 0..=6u32 {
            for lambda in crate::partition::partitions_of(n) {
                let x = Ht::<i128>::monomial(lambda.clone(), 1);
                let back: Ht<i128> = Ht::from_schur(&x.to_schur());
                assert_eq!(back, x, "round trip at {lambda}");
            }
        }
    }

    /// `ḡ^ν_{λμ} = c^ν_{λμ}` when |ν| = |λ| + |μ| — the top-degree part of the
    /// reduced Kronecker product is the Littlewood–Richardson rule.
    ///
    /// This ties the new engine to the oldest one in the crate, and it is the
    /// check that would catch a systematic misreading of the paper: the two
    /// computations have no step in common.
    #[test]
    fn top_degree_is_littlewood_richardson() {
        for lambda in crate::partition::partitions_of(4) {
            for mu in crate::partition::partitions_of(3) {
                let red: St<i128> = reduced_kronecker_product(&lambda, &mu);
                let sl: Schur<i128> = Schur::monomial(lambda.clone(), 1);
                let sm: Schur<i128> = Schur::monomial(mu.clone(), 1);
                let lr = sl.mul(&sm);
                for (nu, c) in lr.terms() {
                    assert_eq!(
                        red.coeff(nu),
                        *c,
                        "top degree of s̃_{lambda} · s̃_{mu} at {nu}"
                    );
                }
            }
        }
    }

    /// The product agrees with the independent `h̃`-matrix route.
    #[test]
    fn two_product_routes_agree() {
        for lambda in crate::partition::partitions_of(3) {
            for mu in crate::partition::partitions_of(3) {
                let a: St<i128> = reduced_kronecker_product(&lambda, &mu);
                let b: St<i128> =
                    reduced_kronecker_via_ht(&lambda, &mu).expect("small enough to enumerate");
                assert_eq!(a, b, "routes disagree on s̃_{lambda} · s̃_{mu}");
            }
        }
    }

    /// `s̃_∅ = 1` is the unit, and the product is commutative.
    #[test]
    fn product_has_a_unit_and_commutes() {
        let a = st(&[2, 1]);
        let b = st(&[2]).add(&st(&[1, 1]));
        assert_eq!(a.mul(&St::unit()), a);
        assert_eq!(a.mul(&b), b.mul(&a));
    }

    /// Associativity, which is where a wrong Γ⁻¹ would show up even if every
    /// individual product looked plausible.
    #[test]
    fn product_is_associative() {
        let a = st(&[2]);
        let b = st(&[1, 1]);
        let c = st(&[2, 1]);
        assert_eq!(a.mul(&b).mul(&c), a.mul(&b.mul(&c)));
    }

    /// The product matches multiplying the two elements as ordinary symmetric
    /// functions — i.e. OZ Thm 7 really is a statement about the outer product,
    /// and the power-sum engine computes that and not something else.
    #[test]
    fn product_agrees_with_the_schur_route() {
        for lambda in crate::partition::partitions_of(3) {
            for mu in crate::partition::partitions_of(2) {
                let fast: St<i128> = reduced_kronecker_product(&lambda, &mu);
                let a: St<i128> = St::monomial(lambda.clone(), 1);
                let b: St<i128> = St::monomial(mu.clone(), 1);
                let slow: St<i128> = St::from_schur(&a.to_schur().mul(&b.to_schur()));
                assert_eq!(fast, slow, "s̃_{lambda} · s̃_{mu}");
            }
        }
    }

    /// Reduced Kronecker coefficients are the *stable* value of the ordinary
    /// ones: `ḡ^ν_{λμ} = g((n−|λ|,λ), (n−|μ|,μ), (n−|ν|,ν))` for n large.
    ///
    /// The only test that reaches the existing Kronecker implementation, and so
    /// the only one that could catch both of them being wrong together — which
    /// it cannot, because they share nothing but `Partition`.
    #[test]
    fn agrees_with_the_ordinary_kronecker_once_stable() {
        use crate::coeff::Rational;
        let lambda = part(&[2]);
        let mu = part(&[1, 1]);
        let red: St<i128> = reduced_kronecker_product(&lambda, &mu);
        let n = 9u32; // ≥ 2(|λ| + |μ|), comfortably inside stability
        let big = |p: &Partition| {
            let mut v = vec![n - p.size()];
            v.extend_from_slice(p.parts());
            Partition::new(v)
        };
        for nu in crate::partition::partitions_of(3) {
            let g: Rational = crate::ops::kronecker(&big(&lambda), &big(&mu), &big(&nu));
            assert_eq!(Rational::from_int(red.coeff(&nu)), g, "ḡ vs g at ν = {nu}");
        }
    }
}
