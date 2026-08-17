//! Plethysm `f[g]` of two symmetric functions in the Schur basis.
//!
//! Sage spells the same operation `s[a](s[b])` (`scripts/compare_sage.py`).
//!
//! There are **two routes**, and which one runs is decided by the *inner*
//! argument alone.
//!
//! ## The general route: through the power-sum basis
//!
//! Any g at all. Plethysm is almost trivial in p:
//!
//! - `p_n[p_m] = p_{nm}`, so `p_n[g]` is just g's p-expansion with every part
//!   multiplied by n (rational coefficients are fixed by `p_n`);
//! - plethysm is additive and *multiplicative* in the outer argument, so
//!   `f[g] = Σ_λ c_λ ∏_i p_{λ_i}[g]` once `f = Σ_λ c_λ p_λ`;
//! - products in the p-basis are multiset unions — the cheapest product in
//!   the crate.
//!
//! So the cost is one p→s conversion at degree d·e, which runs on memoized
//! Murnaghan–Nakayama characters and is essentially the whole runtime
//! (`docs/record/plethysm.md`).
//!
//! ## The one-row route: never leaving the Schur basis
//!
//! When g is `s_m` — one row, coefficient one — that conversion is avoidable
//! entirely, and avoiding it is worth two orders of magnitude: `s_6[s_6]` in
//! 0.25s against 28s. The Newton recursion
//!
//! ```text
//!   n·h_n[g] = Σ_{k=1..n} p_k[g] · h_{n−k}[g]
//! ```
//!
//! builds a ladder of `h_n[g]` using only Littlewood–Richardson products, and
//! the outer argument is carried onto it by its h-expansion, since
//! `h_μ[g] = ∏_i h_{μ_i}[g]`. What makes the recursion usable is that
//! [`adams_one_row`] gives each `p_k[s_m]` combinatorially, off an abacus,
//! with no conversion behind it.
//!
//! Both routes require a [`Plethystic`] ring: `z_μ⁻¹` and the division by n
//! both need division by an integer, and `p_n` acts on the coefficients as
//! well as the parts.

use crate::coeff::{Plethystic, Ring};
use crate::convert::{FromSchur, ToSchur};
use crate::partition::Partition;
use crate::sym::{Homogeneous, PowerSum, Schur, SymFn};

/// `p_k[h_m]` in the Schur basis, as a signed sum with no conversion behind it.
///
/// The power-sum route computes this by expanding in p and converting back at
/// degree k·m, which is the expensive half of every plethysm. It is avoidable
/// here: when every part of the cycle type is divisible by k, χ^λ vanishes
/// unless λ has empty k-core, and otherwise factors through the k-quotient. For
/// h_m the surviving quotients are exactly the k-tuples of **one-row**
/// partitions summing to m, so
///
/// ```text
///   p_k[h_m] = Σ ±s_λ   over (a_0, …, a_{k−1}) with Σa_i = m
/// ```
///
/// — `C(m+k−1, k−1)` terms, 462 at k = m = 6 against p(36) = 17,977.
///
/// λ is read off the abacus with **one bead per runner**: β_i = k·a_i + i,
/// sorted descending, λ_j = β_j − (k−1−j). One bead is enough because these λ
/// never have more than k parts, and it is what keeps the sign O(k²) rather
/// than O((km)²).
///
/// ## The sign is measured against the empty configuration
///
/// The sign is the parity of the permutation sorting the β list — *plus* the
/// parity of the all-zeros configuration, `(−1)^{k(k−1)/2}`, which must come
/// out as λ = ∅ with sign +1. Dropping that normalization leaves the whole sum
/// multiplied by a constant depending on k, which is a **right support with
/// flipped signs**: the answer stays plausible and is wrong for k = 2, 3, 6 and
/// right for k = 4, 5. It went wrong that way three times while this was being
/// derived (`docs/record/plethysm.md`), which is why the doctest below is a k
/// where the constant is −1.
///
/// `p_2[h_2] = s_4 − s_{3,1} + s_{2,2}`. Every sign here is what the
/// normalization decides — an unnormalized k = 2 returns the same three
/// partitions with all three signs flipped, which is why the values and not
/// merely the support are asserted.
///
/// ```
/// use symfn::{Partition, Rational, Ring, Schur, SymFn};
/// let r: Schur<Rational> = symfn::plethysm::adams_one_row(2, 2);
/// let at = |v: &[u32]| r.coeff(&Partition::new(v.iter().copied()));
/// assert_eq!(at(&[4]), Rational::one());
/// assert_eq!(at(&[3, 1]), Rational::one().neg());
/// assert_eq!(at(&[2, 2]), Rational::one());
/// assert_eq!(r.terms().len(), 3);
/// ```
///
/// `m = 0` is `h_0 = 1`, and `p_k[1] = 1`: the only composition is all zeros,
/// which reads back as λ = ∅ with sign +1.
///
/// # Panics
///
/// Panics if `k` is zero: p_0 is not an operation this expresses.
pub fn adams_one_row<C: Ring>(k: u32, m: u32) -> Schur<C> {
    assert!(k > 0, "p_k[h_m] needs k ≥ 1");
    // Runner indices stay `u32` — the width β lives in — so the arithmetic
    // below needs no narrowing cast to state (`docs/policies/failure.md`, R5).
    // `slots` is the same count as a length, and widening to `usize` is exact.
    let slots = k as usize;
    let mut out = Schur::zero();
    // The all-zeros configuration has β ascending 0..k−1, so C(k,2) inversions.
    let base_even = (slots * (slots - 1) / 2).is_multiple_of(2);
    let mut rows = Vec::with_capacity(slots);
    compositions(m, slots, &mut rows, &mut |a: &[u32]| {
        let betas: Vec<u32> = (0..k).zip(a).map(|(i, &x)| k * x + i).collect();
        let inversions = (0..slots)
            .flat_map(|x| (x + 1..slots).map(move |y| (x, y)))
            .filter(|&(x, y)| betas[x] < betas[y])
            .count();
        let mut sorted = betas;
        sorted.sort_unstable_by(|a, b| b.cmp(a));
        // β distinct and descending, so β_j ≥ k−1−j and no part is negative.
        let parts: Vec<u32> = sorted
            .iter()
            .zip(0..k)
            .map(|(&b, j)| b - (k - 1 - j))
            .filter(|&p| p > 0)
            .collect();
        let positive = inversions.is_multiple_of(2) == base_even;
        let c = if positive { C::one() } else { C::one().neg() };
        out.add_term(Partition::from_sorted(parts), c);
    });
    out
}

/// Every `(a_0, …, a_{k−1})` of nonnegative integers summing to `m`.
fn compositions(m: u32, k: usize, cur: &mut Vec<u32>, emit: &mut impl FnMut(&[u32])) {
    if cur.len() + 1 == k {
        cur.push(m);
        emit(cur);
        cur.pop();
        return;
    }
    for a in 0..=m {
        cur.push(a);
        compositions(m - a, k, cur, emit);
        cur.pop();
    }
}

/// `h_0[h_m] … h_upto[h_m]` in the Schur basis, by the Newton recursion
///
/// ```text
///   n·h_n[g] = Σ_{k=1..n} p_k[g] · h_{n−k}[g]
/// ```
///
/// Every product here is an ordinary Littlewood–Richardson product, so the
/// whole ladder is built without ever leaving the Schur basis — which is the
/// point, since the p-basis route's cost is the single conversion at degree
/// n·m that this replaces.
///
/// The division by n is exact: h_n[g] has integer coefficients, and
/// [`QAlgebra::div_u128`](crate::coeff::QAlgebra::div_u128) is what a
/// [`Plethystic`] ring already promises for z_μ⁻¹.
///
/// The ladder is the reusable object rather than its last rung: a caller
/// wanting `s_λ[h_m]` for several λ pays for one ladder
/// (`docs/record/plethysm.md`).
fn h_ladder<C: Plethystic>(m: u32, upto: u32) -> Vec<Schur<C>> {
    // Once each, not once per rung. Every rung consumes every earlier p_k, so
    // building them inside the loop recomputes each one up to `upto` times —
    // it cost 2.1x on s_7[s_7] before this line moved out
    // (`docs/record/plethysm.md`).
    let adams: Vec<Schur<C>> = (1..=upto).map(|k| adams_one_row(k, m)).collect();
    let mut h: Vec<Schur<C>> = vec![Schur::monomial(Partition::default(), C::one())];
    for n in 1..=upto {
        crate::interrupt::poll();
        let mut acc: Schur<C> = Schur::zero();
        for k in 1..=n {
            acc = acc.add(&adams[(k - 1) as usize].mul(&h[(n - k) as usize]));
        }
        let mut hn: Schur<C> = Schur::zero();
        for (lambda, c) in acc.terms() {
            hn.add_term(lambda.clone(), c.div_u128(u128::from(n)));
        }
        h.push(hn);
    }
    h
}

/// `f[s_m]` for a one-row inner argument, without converting at degree d·m.
///
/// The outer argument is expanded in the **h basis** once — plethysm is a ring
/// homomorphism in it, so `h_μ[g] = ∏_i h_{μ_i}[g]` — and every factor is a
/// rung of the `h_ladder`. That replaces the Jacobi–Trudi determinant the same
/// identity suggests, which would be factorial in the number of rows.
///
/// Measured against the power-sum route on this laptop: `s_6[s_6]` in 0.26s
/// where the conversion route took 27s, agreeing on all 2002 terms
/// (`docs/record/plethysm.md`).
///
/// # Panics
///
/// Panics if `m` is zero. Panics where the ring does: a fixed-width `C` whose
/// values leave its width panics rather than wrapping, and `GuardedRat`
/// reports instead.
pub fn plethysm_by_one_row<C: Plethystic>(f: &Schur<C>, m: u32) -> Schur<C> {
    assert!(m > 0, "f[s_m] needs m ≥ 1");
    let hf: Homogeneous<C> = Homogeneous::from_schur(f);
    // No parts at all is degree 0, not nothing: a constant outer argument is
    // h_∅ = 1, whose plethysm is itself. Reading the absent maximum as "return
    // zero" made `s_∅[s_1]` vanish, which the Sage fixture caught.
    let upto = hf
        .terms()
        .keys()
        .flat_map(Partition::parts)
        .copied()
        .max()
        .unwrap_or(0);
    let ladder = h_ladder::<C>(m, upto);

    let mut out = Schur::zero();
    for (mu, c) in hf.terms() {
        crate::interrupt::poll();
        let mut prod: Schur<C> = Schur::monomial(Partition::default(), C::one());
        for &part in mu.parts() {
            prod = prod.mul(&ladder[part as usize]);
        }
        out = out.add(&prod.scale(c));
    }
    out
}

/// The one-row inner argument `s_m`, if `g` is exactly that.
///
/// Exactly: one term, one row, coefficient one. A coefficient other than one
/// would need `p_k[c·g] = c^k p_k[g]` threaded through the ladder, which is
/// expressible but is not what the measured case needs.
fn one_row_inner<C: Ring>(g: &Schur<C>) -> Option<u32> {
    let terms = g.terms();
    if terms.len() != 1 {
        return None;
    }
    let (mu, c) = terms.iter().next()?;
    (mu.len() == 1 && *c == C::one()).then(|| mu.part(0))
}

/// `p_n[g]`: substitute `p_k ↦ p_{nk}` throughout g's power-sum expansion, and
/// apply the same substitution to the **coefficients**.
///
/// The coefficient half is easy to miss. `p_n` substitutes into the alphabet,
/// and the variables of a coefficient ring are part of that alphabet, so
/// `p_n[t·p_1] = t^n·p_n`. Over ℚ there is nothing to raise and
/// [`Plethystic::frobenius`] is the identity, so only a ring like `ℚ[t]` makes
/// the coefficient half visible.
fn scale_parts<C: Plethystic>(g: &PowerSum<C>, n: u32) -> PowerSum<C> {
    let mut out = PowerSum::zero();
    for (mu, d) in g.terms() {
        let scaled = Partition::new(mu.parts().iter().map(|&x| x * n));
        out.add_term(scaled, d.frobenius(n));
    }
    out
}

/// The plethysm `f[g]` of two symmetric functions in the Schur basis.
///
/// If f and g are homogeneous of degrees d and e, the result is homogeneous of
/// degree d·e.
///
/// # Panics
///
/// Panics if a value leaves the width of `C`. A `bignum` ring has no wall.
/// [`GuardedRat`](crate::GuardedRat) reports the overflow instead of panicking.
pub fn plethysm<C: Plethystic>(f: &Schur<C>, g: &Schur<C>) -> Schur<C> {
    // A one-row inner argument takes the ladder, which never converts at
    // degree d·m and is the difference between 0.26s and 27s on s_6[s_6]
    // (`docs/record/plethysm.md`). Every other g takes the power-sum route
    // below, which is general.
    if let Some(m) = one_row_inner(g) {
        return plethysm_by_one_row(f, m);
    }
    plethysm_via_power_sum(f, g)
}

/// The general route: expand both arguments in p, where plethysm is part
/// scaling and a multiset union, and convert the answer back.
///
/// Kept reachable by name because it is what the ladder is checked against —
/// two routes with no shared machinery, agreeing term for term
/// (`the_two_routes_are_one_operation`).
fn plethysm_via_power_sum<C: Plethystic>(f: &Schur<C>, g: &Schur<C>) -> Schur<C> {
    let pf: PowerSum<C> = PowerSum::from_schur(f);
    let pg: PowerSum<C> = PowerSum::from_schur(g);

    let mut acc = PowerSum::zero();
    for (lambda, c) in pf.terms() {
        // p_λ[g] = ∏_i p_{λ_i}[g]
        let mut prod: PowerSum<C> = PowerSum::monomial(Partition::default(), C::one());
        for &part in lambda.parts() {
            crate::interrupt::poll();
            prod = prod.mul(&scale_parts(&pg, part));
        }
        acc = acc.add(&prod.scale(c));
    }
    acc.to_schur()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::{Rational, Ring};

    fn s(v: &[u32]) -> Schur<Rational> {
        Schur::monomial(Partition::new(v.iter().copied()), Rational::one())
    }

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    #[test]
    fn classic_plethysms() {
        // s_2[s_2] = s_4 + s_{22}
        let r = plethysm(&s(&[2]), &s(&[2]));
        assert_eq!(r.coeff(&part(&[4])), Rational::one());
        assert_eq!(r.coeff(&part(&[2, 2])), Rational::one());
        assert_eq!(r.terms().len(), 2);

        // s_{11}[s_2] = s_{31}
        let r = plethysm(&s(&[1, 1]), &s(&[2]));
        assert_eq!(r.coeff(&part(&[3, 1])), Rational::one());
        assert_eq!(r.terms().len(), 1);
    }

    /// A constant outer argument is degree 0, not nothing.
    ///
    /// `s_∅ = 1` and `1[g] = 1`, but the ladder route reads the largest part of
    /// the h-expansion to size itself, and ∅ has no largest part. Treating that
    /// absence as "no work" returned zero. The Sage fixture caught it on
    /// `s_∅[s_1]`; this states it where the reason is visible.
    #[test]
    fn a_constant_outer_argument_survives_the_ladder() {
        let one: Schur<Rational> = Schur::monomial(Partition::default(), Rational::one());
        for m in 1..=3u32 {
            assert_eq!(plethysm_by_one_row(&one, m), one, "1[s_{m}] is not 1");
        }
        let zero: Schur<Rational> = Schur::zero();
        assert_eq!(plethysm_by_one_row(&zero, 2), zero, "0[s_2] is not 0");
    }

    #[test]
    fn plethysm_identities() {
        // s_1[g] = g and f[s_1] = f
        let g = s(&[2, 1]);
        assert_eq!(plethysm(&s(&[1]), &g), g);
        assert_eq!(plethysm(&g, &s(&[1])), g);
    }

    /// The ladder and the power-sum route share no machinery — one is
    /// Littlewood–Richardson products of an abacus rule, the other is z_μ,
    /// characters and a β-mask sweep — so agreeing on every term is the
    /// strongest statement available here, and the only one at degrees Sage
    /// cannot reach to check.
    ///
    /// The outer partitions deliberately include more than one row: those
    /// exercise the h-expansion of the outer argument, where a wrong Kostka
    /// inverse would show, while one-row outers read a single rung and would
    /// not.
    #[test]
    fn the_two_routes_are_one_operation() {
        for m in 1..=4u32 {
            for parts in [
                &[1u32][..],
                &[2],
                &[3],
                &[1, 1],
                &[2, 1],
                &[2, 2],
                &[3, 1],
                &[2, 1, 1],
            ] {
                let f = s(parts);
                let ladder = plethysm_by_one_row(&f, m);
                let via_p = plethysm_via_power_sum(&f, &s(&[m]));
                assert_eq!(ladder, via_p, "s{parts:?}[s_{m}] differs by route");
            }
        }
    }

    /// A one-row *outer* is the case the ladder answers with no product at
    /// all, and the classical values pin the convention rather than merely the
    /// agreement above.
    #[test]
    fn the_ladder_reproduces_the_classical_plethysms() {
        // s_2[s_2] = s_4 + s_{22}
        let r = plethysm_by_one_row(&s(&[2]), 2);
        assert_eq!(r.coeff(&part(&[4])), Rational::one());
        assert_eq!(r.coeff(&part(&[2, 2])), Rational::one());
        assert_eq!(r.terms().len(), 2);

        // The two orders, which is the pair that pins the convention: they have
        // the same degree and different answers, so a route that silently
        // transposed its arguments would pass either one alone.
        // s_2[s_3] = s_6 + s_{4,2}, and s_3[s_2] = s_6 + s_{4,2} + s_{2,2,2}
        // (both checked against Sage).
        let r = plethysm_by_one_row(&s(&[2]), 3);
        for lambda in [&[6u32][..], &[4, 2]] {
            assert_eq!(
                r.coeff(&part(lambda)),
                Rational::one(),
                "s_2[s_3] misses s{lambda:?}"
            );
        }
        assert_eq!(r.terms().len(), 2, "s_2[s_3] is not two terms");

        let r = plethysm_by_one_row(&s(&[3]), 2);
        for lambda in [&[6u32][..], &[4, 2], &[2, 2, 2]] {
            assert_eq!(
                r.coeff(&part(lambda)),
                Rational::one(),
                "s_3[s_2] misses s{lambda:?}"
            );
        }
        assert_eq!(r.terms().len(), 3, "s_3[s_2] is not three terms");
    }

    #[test]
    fn degree_is_multiplicative() {
        for f in [&[2][..], &[1, 1], &[2, 1]] {
            for g in [&[2][..], &[1, 1]] {
                let r = plethysm(&s(f), &s(g));
                let want = part(f).size() * part(g).size();
                assert_eq!(r.degree(), Some(want), "deg s{:?}[s{:?}]", f, g);
            }
        }
    }
}
