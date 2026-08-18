//! Plethysm `f[g]` of two symmetric functions in the Schur basis.
//!
//! Sage spells the same operation `s[a](s[b])` (`scripts/compare_sage.py`).
//!
//! There are **two routes**, and which one runs is decided by both shapes:
//! see `takes_the_ladder`. Neither dominates, and the spread between them at
//! one degree is three orders of magnitude in each direction.
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
//! ## The ladder route: never leaving the Schur basis
//!
//! That conversion is avoidable entirely, and avoiding it is worth two orders
//! of magnitude where it is the cost: `s_6[s_6]` in 0.25s against 28s, and
//! `s_3[s_{10,10}]` at degree 60 in 0.75s where the conversion route does not
//! finish at all. The Newton recursion
//!
//! ```text
//!   n·h_n[g] = Σ_{k=1..n} p_k[g] · h_{n−k}[g]
//! ```
//!
//! builds a ladder of `h_n[g]` using only Littlewood–Richardson products, and
//! the outer argument is carried onto it by its h-expansion, since
//! `h_μ[g] = ∏_i h_{μ_i}[g]`. What makes the recursion usable is that each
//! `p_k[g]` is available combinatorially off an abacus with no conversion
//! behind it — [`adams_one_row`] in closed form when g is a single row, and
//! [`adams`] by a k-quotient sweep otherwise.
//!
//! **The ladder is not always cheaper**, and that is why the choice is
//! measured rather than assumed. Its cost follows the *outer's* h-depth, since
//! it builds one rung per level and each rung needs a `p_k[g]`; the conversion
//! route's follows the total degree. A general inner makes each `p_k[g]` a
//! Littlewood–Richardson sweep rather than a closed form, so a deep outer over
//! a small inner is the case the ladder loses badly
//! (`docs/record/plethysm.md`).
//!
//! Both routes require a [`Plethystic`] ring: `z_μ⁻¹` and the division by n
//! both need division by an integer, and `p_n` acts on the coefficients as
//! well as the parts.

use crate::coeff::{Plethystic, Ring};
use crate::convert::FromSchur;
use crate::partition::Partition;
use crate::sym::{Homogeneous, Schur, SymFn};

use crate::convert::ToSchur;
use crate::sym::PowerSum;

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

/// `p_k[s_ν]` in the Schur basis for **any** ν, by the same k-quotient rule.
///
/// [`adams_one_row`] is the case where the inner product below is 1 exactly
/// when every quotient is a row. In general it is a
/// multi-Littlewood–Richardson coefficient:
///
/// ```text
///   p_k[s_ν] = Σ ⟨s_ν, ∏_i s_{ν^(i)}⟩ · ±s_λ
/// ```
///
/// over k-tuples with `Σ|ν^(i)| = |ν|`, λ having empty k-core and that
/// k-quotient.
///
/// **The tuples are swept with a shared prefix, not enumerated and filtered.**
/// Every tuple agreeing in its first slots shares the same partial product, so
/// the sweep extends one slot at a time and reads each tuple's coefficient off
/// the layer it lands in — the same trick `p → s` batching uses. Partial
/// products are pruned by containment in ν: later slots only add cells, so a
/// term already outgrowing ν in some row can never come back. Enumerating
/// instead wastes 1.6x to 9.3x on tuples that contribute nothing, and the
/// pruning is what makes those never exist (`docs/record/plethysm.md`).
///
/// Coefficients are `i128` because they are multi-LR coefficients — counts —
/// and the ring only enters when the caller scales by them.
fn adams_schur(k: u32, nu: &Partition) -> Vec<(Partition, i128)> {
    let n = nu.size();
    if n == 0 {
        return vec![(Partition::default(), 1)];
    }
    let beads = n + 1;
    let base_positive = quotient_sign(&vec![Partition::default(); k as usize], k, beads).1;

    let cx = Sweep {
        k,
        nu,
        beads,
        base_positive,
    };
    let mut out: Vec<(Partition, i128)> = Vec::new();
    let mut quot: Vec<Partition> = Vec::with_capacity(k as usize);
    let root: Vec<(Partition, i128)> = vec![(Partition::default(), 1)];
    sweep_slots(&cx, n, 0, &mut quot, &root, &mut out);
    out
}

/// One slot of [`adams_schur`]'s sweep: choose this runner's partition, extend
/// the partial product, and recurse on what survives the containment prune.
struct Sweep<'a> {
    k: u32,
    nu: &'a Partition,
    beads: u32,
    base_positive: bool,
}

fn sweep_slots(
    cx: &Sweep,
    left: u32,
    slot: u32,
    quot: &mut Vec<Partition>,
    layer: &[(Partition, i128)],
    out: &mut Vec<(Partition, i128)>,
) {
    let (k, nu, beads, base_positive) = (cx.k, cx.nu, cx.beads, cx.base_positive);
    if slot == k {
        if left == 0 {
            if let Some(&(_, c)) = layer.iter().find(|(lambda, _)| lambda == nu) {
                let (lambda, positive) = quotient_sign(quot, k, beads);
                let signed = if positive == base_positive { c } else { -c };
                out.push((lambda, signed));
            }
        }
        return;
    }
    crate::interrupt::poll();
    let remaining_slots = k - slot - 1;
    for size in 0..=left {
        // The last slot must take everything left; an earlier one may not take
        // so much that the rest cannot be filled (they can each take zero, so
        // only the last is constrained).
        if remaining_slots == 0 && size != left {
            continue;
        }
        for q in crate::partitions_of(size) {
            let next = extend_layer(layer, &q, nu);
            if next.is_empty() {
                continue;
            }
            quot.push(q);
            sweep_slots(cx, left - size, slot + 1, quot, &next, out);
            quot.pop();
        }
    }
}

/// Multiply a partial product by `s_q` and keep only what can still reach ν.
///
/// Containment is the whole prune: `s_α · s_β` produces only λ ⊇ α, and every
/// later slot adds cells, so a term outside ν is dead rather than merely
/// unpromising.
fn extend_layer(
    layer: &[(Partition, i128)],
    q: &Partition,
    nu: &Partition,
) -> Vec<(Partition, i128)> {
    if q.is_empty() {
        return layer.to_vec();
    }
    let mut acc: std::collections::BTreeMap<Partition, i128> = std::collections::BTreeMap::new();
    for (lambda, c) in layer {
        let a: Schur<i128> = Schur::monomial(lambda.clone(), *c);
        let b: Schur<i128> = Schur::monomial(q.clone(), 1);
        for (mu, d) in a.mul(&b).terms() {
            if fits_inside(mu, nu) {
                *acc.entry(mu.clone()).or_insert(0) += *d;
            }
        }
    }
    acc.into_iter().filter(|&(_, c)| c != 0).collect()
}

/// Whether λ ⊆ ν as diagrams.
fn fits_inside(lambda: &Partition, nu: &Partition) -> bool {
    lambda.len() <= nu.len() && (0..lambda.len()).all(|i| lambda.part(i) <= nu.part(i))
}

/// λ with empty k-core and k-quotient `quot`, and whether its sort is even.
///
/// `beads` is passed in and **fixed across the whole sum**. Deriving it from
/// the tuple compares terms computed in different abacuses, which is a wrong
/// sign on a right support — the defect this construction produced four
/// separate times (`docs/record/plethysm.md`). Fixing it at the caller is what
/// makes that unavailable rather than merely discouraged.
fn quotient_sign(quot: &[Partition], k: u32, beads: u32) -> (Partition, bool) {
    // Bead and runner indices stay `u32`, the width β lives in, so nothing
    // here narrows (`docs/policies/failure.md`, R5).
    let mut betas: Vec<u32> = Vec::with_capacity(quot.len() * beads as usize);
    for (i, q) in (0..).zip(quot) {
        for j in 0..beads {
            let part = if (j as usize) < q.len() {
                q.part(j as usize)
            } else {
                0
            };
            betas.push(k * (part + beads - 1 - j) + i);
        }
    }
    let inversions = (0..betas.len())
        .flat_map(|x| (x + 1..betas.len()).map(move |y| (x, y)))
        .filter(|&(x, y)| betas[x] < betas[y])
        .count();
    let total = k * beads;
    let mut sorted = betas;
    sorted.sort_unstable_by(|a, b| b.cmp(a));
    let parts: Vec<u32> = sorted
        .iter()
        .zip(0..total)
        .map(|(&b, j)| b - (total - 1 - j))
        .filter(|&p| p > 0)
        .collect();
    (Partition::from_sorted(parts), inversions.is_multiple_of(2))
}

/// `p_k[g]` in the Schur basis, for any g.
///
/// The Adams operation ψ^k is a ring homomorphism, so it is **additive** in its
/// argument and the rule above extends to a whole element term by term. That
/// additivity is what removes the one-row restriction from the ladder: the
/// recursion never cared what g was, only that `p_k[g]` was available without
/// a conversion.
///
/// **`p_k` acts on the coefficients too**, through
/// [`Plethystic::frobenius`]: `p_k` is a
/// substitution on the alphabet and a coefficient ring's variables are part of
/// that alphabet, so `p_k[t·s_ν] = t^k · p_k[s_ν]`. Invisible over ℚ, and the
/// whole answer over `ℚ[q,t]`.
///
/// # Panics
///
/// Panics if `k` is zero.
pub fn adams<C: Plethystic>(k: u32, g: &Schur<C>) -> Schur<C> {
    assert!(k > 0, "p_k[g] needs k ≥ 1");
    let mut out = Schur::zero();
    for (nu, c) in g.terms() {
        crate::interrupt::poll();
        // One row is the case with a closed form over compositions, and it is
        // the one the ladder hits most, so it keeps its own path.
        let terms: Vec<(Partition, i128)> = if nu.len() == 1 {
            adams_one_row::<i128>(k, nu.part(0))
                .terms()
                .iter()
                .map(|(l, v)| (l.clone(), *v))
                .collect()
        } else {
            adams_schur(k, nu)
        };
        let scale = c.frobenius(k);
        for (lambda, v) in terms {
            out.add_term(lambda, C::from_i128(v).mul(&scale));
        }
    }
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
fn h_ladder<C: Plethystic>(g: &Schur<C>, upto: u32) -> Vec<Schur<C>> {
    // Once each, not once per rung. Every rung consumes every earlier p_k, so
    // building them inside the loop recomputes each one up to `upto` times —
    // it cost 2.1x on s_7[s_7] before this line moved out
    // (`docs/record/plethysm.md`).
    let adams: Vec<Schur<C>> = (1..=upto).map(|k| adams(k, g)).collect();
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
pub fn plethysm_via_ladder<C: Plethystic>(f: &Schur<C>, g: &Schur<C>) -> Schur<C> {
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
    let ladder = h_ladder::<C>(g, upto);

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
    if takes_the_ladder(f, g) {
        return plethysm_via_ladder(f, g);
    }
    plethysm_via_power_sum(f, g)
}

/// Which route is cheaper, decided from the shapes rather than tried.
///
/// The two costs are driven by different things, which is what makes the
/// choice predictable rather than a guess. The ladder builds one rung per
/// level of the **outer's h-depth** and each rung needs `p_k[g]`; the
/// power-sum route pays for one conversion at degree d·e whatever the shapes.
///
/// So a one-row inner always takes the ladder: [`adams_one_row`] is a closed
/// form over compositions and stays cheap however deep the ladder goes, and
/// the measured win runs from 20x to 1,500x.
///
/// A general inner makes each `p_k[g]` a k-fold Littlewood–Richardson sweep
/// whose cost climbs with k, so a deep ladder stops being worth it. Measured
/// at degree 36, `s_3[s_{6,6}]` is **1062x** and `s_4[s_{4,4}]` **34x** for
/// the ladder, while `s_9[s_{2,2}]` is 0.07x and `s_{12}[s_{2,1}]` 0.02x —
/// three orders of magnitude of spread at one degree, ordered by depth against
/// |ν| and nothing else.
///
/// **The threshold is a fit, not a theorem.** `depth ≤ |ν|` separates every
/// case measured (`docs/record/plethysm.md`) and is stated in the units the
/// two costs actually scale in, but six points do not make a law; a case that
/// straddles it belongs in the record rather than in a silently moved
/// constant.
fn takes_the_ladder<C: Plethystic>(f: &Schur<C>, g: &Schur<C>) -> bool {
    let inner_is_one_row = g.terms().keys().all(|nu| nu.len() <= 1);
    if inner_is_one_row {
        return true;
    }
    let depth = Homogeneous::<C>::from_schur(f)
        .terms()
        .keys()
        .flat_map(Partition::parts)
        .copied()
        .max()
        .unwrap_or(0);
    let inner_size = g.terms().keys().map(Partition::size).max().unwrap_or(0);
    depth <= inner_size
}

/// The general route: expand both arguments in p, where plethysm is part
/// scaling and a multiset union, and convert the answer back.
///
/// Kept reachable by name because it is what the ladder is checked against —
/// two routes with no shared machinery, agreeing term for term
/// (`the_two_routes_are_one_operation`). It is both the oracle and the route
/// a deep outer over a general inner still takes.
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
            assert_eq!(
                plethysm_via_ladder(&one, &s(&[m])),
                one,
                "1[s_{m}] is not 1"
            );
        }
        let zero: Schur<Rational> = Schur::zero();
        assert_eq!(
            plethysm_via_ladder(&zero, &s(&[2])),
            zero,
            "0[s_2] is not 0"
        );
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
        let outers: [&[u32]; 8] = [
            &[1],
            &[2],
            &[3],
            &[1, 1],
            &[2, 1],
            &[2, 2],
            &[3, 1],
            &[2, 1, 1],
        ];
        // Inners with more than one row are the case the general Adams
        // operation exists for; a one-row inner takes the closed form over
        // compositions and would not exercise the k-quotient sweep at all.
        let inners: [&[u32]; 7] = [&[1], &[2], &[3], &[4], &[2, 1], &[2, 2], &[3, 1]];
        for inner in inners {
            for outer in outers {
                let (f, g) = (s(outer), s(inner));
                let ladder = plethysm_via_ladder(&f, &g);
                let via_p = plethysm_via_power_sum(&f, &g);
                assert_eq!(ladder, via_p, "s{outer:?}[s{inner:?}] differs by route");
            }
        }
    }

    /// A multi-term inner is where additivity of the Adams operation is doing
    /// the work, and where a route that only handled single Schur inners would
    /// still pass everything above.
    #[test]
    fn a_sum_as_the_inner_argument_agrees_by_route() {
        let mut g: Schur<Rational> = Schur::zero();
        g.add_term(part(&[2]), Rational::one());
        g.add_term(part(&[1, 1]), Rational::one());
        for outer in [&[2u32][..], &[1, 1], &[3], &[2, 1]] {
            let f = s(outer);
            assert_eq!(
                plethysm_via_ladder(&f, &g),
                plethysm_via_power_sum(&f, &g),
                "s{outer:?}[s_2 + s_{{1,1}}] differs by route"
            );
        }
        // And a negative coefficient, which additivity has to carry too.
        let mut d: Schur<Rational> = Schur::zero();
        d.add_term(part(&[2]), Rational::one());
        d.add_term(part(&[1, 1]), Rational::one().neg());
        for outer in [&[2u32][..], &[3], &[2, 1]] {
            let f = s(outer);
            assert_eq!(
                plethysm_via_ladder(&f, &d),
                plethysm_via_power_sum(&f, &d),
                "s{outer:?}[s_2 − s_{{1,1}}] differs by route"
            );
        }
    }

    /// A one-row *outer* is the case the ladder answers with no product at
    /// all, and the classical values pin the convention rather than merely the
    /// agreement above.
    #[test]
    fn the_ladder_reproduces_the_classical_plethysms() {
        // s_2[s_2] = s_4 + s_{22}
        let r = plethysm_via_ladder(&s(&[2]), &s(&[2]));
        assert_eq!(r.coeff(&part(&[4])), Rational::one());
        assert_eq!(r.coeff(&part(&[2, 2])), Rational::one());
        assert_eq!(r.terms().len(), 2);

        // The two orders, which is the pair that pins the convention: they have
        // the same degree and different answers, so a route that silently
        // transposed its arguments would pass either one alone.
        // s_2[s_3] = s_6 + s_{4,2}, and s_3[s_2] = s_6 + s_{4,2} + s_{2,2,2}
        // (both checked against Sage).
        let r = plethysm_via_ladder(&s(&[2]), &s(&[3]));
        for lambda in [&[6u32][..], &[4, 2]] {
            assert_eq!(
                r.coeff(&part(lambda)),
                Rational::one(),
                "s_2[s_3] misses s{lambda:?}"
            );
        }
        assert_eq!(r.terms().len(), 2, "s_2[s_3] is not two terms");

        let r = plethysm_via_ladder(&s(&[3]), &s(&[2]));
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
