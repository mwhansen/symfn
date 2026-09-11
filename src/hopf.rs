//! The Hopf-algebra structure of Sym: skewing, comultiplication, counit, and
//! antipode.
//!
//! - skew Schur:   s_{λ/μ} = Σ_ν c^λ_{μν} s_ν
//! - skewing:      g^⊥, the adjoint of multiplication by g ([`SkewBy`])
//! - coproduct:    Δ(s_λ)  = Σ_{μ,ν} c^λ_{μν} s_μ ⊗ s_ν
//! - counit:       ε(s_λ)  = δ_{λ,∅}
//! - antipode:     S(s_λ)  = (−1)^{|λ|} s_{λ'}
//!
//! ## References
//!
//! - **\[M\]** I. G. Macdonald, *Symmetric Functions and Hall Polynomials*,
//!   2nd ed., Oxford, 1995 — I (5.3) for the skew expansion; I.5
//!   Example 3 for `g^⊥` and `s_μ^⊥ s_λ = s_{λ/μ}`; I.5 Example 25(a)–(c)
//!   for the coproduct, the counit, and `Δ s_λ = Σ_μ s_{λ/μ} ⊗ s_μ`, which
//!   with (5.3) is the coproduct formula above.

use crate::coeff::Ring;
use crate::convert::ToSchur;
use crate::memo::partitions_cached;
use crate::partition::Partition;
use crate::sym::{
    impl_linear_ops, Elementary, Forgotten, Homogeneous, InPlaceArith, Monomial, PowerSum, Schur,
    SymFn,
};
use std::collections::BTreeMap;

/// An element of Sym ⊗ Sym in the Schur basis: a formal `C`-combination of
/// pairs of partitions (s_μ ⊗ s_ν). The codomain of the coproduct.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SymTensor<C: Ring>(BTreeMap<(Partition, Partition), C>);

impl<C: Ring> SymTensor<C> {
    /// The zero tensor.
    pub fn zero() -> Self {
        SymTensor(BTreeMap::new())
    }

    /// The terms, keyed `(μ, ν)` for s_μ ⊗ s_ν. Explicit zeros are never
    /// stored, so the map is empty exactly when the tensor is.
    pub fn terms(&self) -> &BTreeMap<(Partition, Partition), C> {
        &self.0
    }

    /// Whether this is 0.
    pub fn is_zero(&self) -> bool {
        self.0.is_empty()
    }

    /// Coefficient of s_μ ⊗ s_ν.
    pub fn coeff(&self, mu: &Partition, nu: &Partition) -> C {
        self.0
            .get(&(mu.clone(), nu.clone()))
            .cloned()
            .unwrap_or_else(C::zero)
    }

    /// `self += c · (s_μ ⊗ s_ν)`, keeping the map free of explicit zeros.
    pub fn add_term(&mut self, key: (Partition, Partition), c: C) {
        if c.is_zero() {
            return;
        }
        let now_zero = {
            let e = self.0.entry(key.clone()).or_insert_with(C::zero);
            e.add_assign(&c);
            e.is_zero()
        };
        if now_zero {
            self.0.remove(&key);
        }
    }

    /// The sum, with terms that cancel dropped rather than stored as zero.
    pub fn add(&self, other: &Self) -> Self {
        let mut out = self.clone();
        for (k, c) in &other.0 {
            out.add_term(k.clone(), c.clone());
        }
        out
    }
}

impl<C: Ring> InPlaceArith<C> for SymTensor<C> {
    fn add_from(&mut self, other: &Self) {
        for (k, c) in &other.0 {
            self.add_term(k.clone(), c.clone());
        }
    }

    fn add_owned(&mut self, mut other: Self) {
        // Addition commutes, so the side with more terms keeps its map.
        if self.0.len() < other.0.len() {
            std::mem::swap(self, &mut other);
        }
        for (k, c) in std::mem::take(&mut other.0) {
            self.add_term(k, c);
        }
    }

    fn sub_from(&mut self, other: &Self) {
        for (k, c) in &other.0 {
            self.add_term(k.clone(), c.neg());
        }
    }

    fn negate(&mut self) {
        // -c is zero only when c is, so no term can cancel here.
        for c in self.0.values_mut() {
            *c = c.neg();
        }
    }

    fn scale_by(&mut self, c: &C) {
        if c.is_zero() {
            self.0.clear();
            return;
        }
        for v in self.0.values_mut() {
            *v = v.mul(c);
        }
        // A ring with zero divisors can send a nonzero product to zero.
        self.0.retain(|_, v| !v.is_zero());
    }
}

impl_linear_ops!(SymTensor);

/// The skew Schur function s_{λ/μ} = Σ_ν c^λ_{μν} s_ν.
///
/// Returns zero when μ ⊄ λ.
pub fn skew_schur<C: Ring>(lambda: &Partition, mu: &Partition) -> Schur<C> {
    let mut out = Schur::zero();
    if !lambda.contains(mu) {
        return out;
    }
    // One traversal of the shape yields every ν with a nonzero coefficient.
    // (This used to sweep all p(n) partitions, running a full LR backtrack per
    // candidate — the same answer for orders of magnitude more work.)
    for (nu, c) in crate::skew_lr::expand_skew_shared(lambda, mu).iter() {
        crate::interrupt::poll();
        out.add_term(nu.clone(), C::from_u128(*c));
    }
    out
}

// --- skewing by an arbitrary symmetric function ------------------------------

/// Skewing: `f.skew_by(g)` is g^⊥ f, the **adjoint of multiplication by g**
/// under the Hall inner product,
///
/// ```text
///   ⟨g^⊥ f, h⟩ = ⟨f, g·h⟩   for every h.
/// ```
///
/// [`skew_schur`] is the special case g = s_μ, f = s_λ, which is why the
/// classical notation for it is a quotient: s_μ^⊥ s_λ = s_{λ/μ}.
///
/// The generic parameter is the *basis of g*, and that is the whole design.
/// The map g ↦ g^⊥ is linear, so any g could be expanded into Schur and handed
/// to the Littlewood–Richardson route. But three of the six bases have adjoints
/// with far cheaper direct rules, each of them a Pieri or Murnaghan–Nakayama
/// rule read backwards:
///
/// | g in basis | g^⊥ on s_λ                            | machinery |
/// |------------|---------------------------------------|-----------|
/// | h          | remove a horizontal strip             | Pieri     |
/// | e          | remove a vertical strip               | dual Pieri|
/// | p          | remove a rim hook, signed by height    | M–N       |
/// | s, m, f    | Σ_μ d_μ s_{λ/μ}                        | LR        |
///
/// The p case is the one that most repays a native path: it needs **no
/// division**, while routing a power sum through the Schur basis would require
/// ℚ for the s → p direction. So `Schur<i64>` can be skewed by a power sum, and
/// the LR route could not have offered that at all.
///
/// Because (fg)^⊥ = f^⊥ g^⊥, a multi-part index just iterates: h_μ^⊥ is
/// h_{μ_1}^⊥ ∘ h_{μ_2}^⊥ ∘ ….
pub trait SkewBy<C: Ring, G> {
    /// `g^⊥` applied to `self` — the adjoint of multiplication by `g` under the
    /// Hall inner product.
    fn skew_by(&self, g: &G) -> Self;
}

/// g = Σ_μ d_μ s_μ, so g^⊥ f = Σ_μ d_μ s_{f/μ}. The general route, and the
/// fallback the m and f bases use after expanding.
impl<C: Ring> SkewBy<C, Schur<C>> for Schur<C> {
    fn skew_by(&self, g: &Schur<C>) -> Self {
        let mut out = Schur::zero();
        for (lambda, c) in self.terms() {
            for (mu, d) in g.terms() {
                if !lambda.contains(mu) {
                    continue;
                }
                let cd = c.mul(d);
                if cd.is_zero() {
                    continue;
                }
                for (nu, k) in crate::skew_lr::expand_skew_shared(lambda, mu).iter() {
                    out.add_term(nu.clone(), C::from_u128(*k).mul(&cd));
                }
            }
        }
        out
    }
}

/// h_r^⊥ s_λ = Σ s_ν over ν with λ/ν a horizontal r-strip.
///
/// Pieri says h_r·s_ν = Σ s_λ over λ with λ/ν a horizontal r-strip, all
/// coefficients 1; the adjoint therefore *removes* what Pieri adds, and the
/// multiplicity-free-ness of Pieri is why no coefficient appears here either.
impl<C: Ring> SkewBy<C, Homogeneous<C>> for Schur<C> {
    fn skew_by(&self, g: &Homogeneous<C>) -> Self {
        let mut out = Schur::zero();
        for (mu, d) in g.terms() {
            let mut acc = self.clone();
            for &r in mu.parts() {
                acc = strip_off(&acc, r);
                if acc.is_zero() {
                    break;
                }
            }
            out = out.add(&acc.scale(d));
        }
        out
    }
}

/// e_r^⊥ s_λ = Σ s_ν over ν with λ/ν a *vertical* r-strip.
///
/// Implemented as ω ∘ h_μ^⊥ ∘ ω rather than by a second strip enumerator. ω is
/// an isometry with ω(h_r) = e_r, so
///
/// ```text
///   ⟨e_r^⊥ s_λ, s_ν⟩ = ⟨s_λ, e_r s_ν⟩ = ⟨s_{λ'}, h_r s_{ν'}⟩ = ⟨h_r^⊥ s_{λ'}, s_{ν'}⟩
/// ```
///
/// — an identity, not an approximation, and it costs one transpose per term
/// against a duplicate of the enumerator and its own separate bugs.
impl<C: Ring> SkewBy<C, Elementary<C>> for Schur<C> {
    fn skew_by(&self, g: &Elementary<C>) -> Self {
        let as_h = Homogeneous::from_terms(g.terms().clone());
        SkewBy::skew_by(&self.omega(), &as_h).omega()
    }
}

/// p_r^⊥ s_λ = Σ (−1)^{ht(λ/ν)} s_ν over ν obtained by **removing a rim hook**
/// of size r from λ.
///
/// Murnaghan–Nakayama read backwards: p_r·s_ν = Σ (−1)^{ht} s_λ over λ ⊇ ν with
/// λ/ν a rim hook of size r, so the adjoint removes them with the same sign.
/// The enumeration is the same β-number bit arithmetic the character table
/// runs on.
///
/// Note the ring bound: only [`Ring`], not [`Field`](crate::coeff::Field). The
/// generic route would have to expand p_μ into Schur, which needs division by
/// z_μ; this one never leaves ℤ.
impl<C: Ring> SkewBy<C, PowerSum<C>> for Schur<C> {
    fn skew_by(&self, g: &PowerSum<C>) -> Self {
        let mut out = Schur::zero();
        for (mu, d) in g.terms() {
            let mut acc = self.clone();
            for &r in mu.parts() {
                let mut next = Schur::zero();
                for (lambda, c) in acc.terms() {
                    for (nu, height) in crate::character::border_strips(lambda, r) {
                        next.add_term(nu, if height % 2 == 1 { c.neg() } else { c.clone() });
                    }
                }
                acc = next;
                if acc.is_zero() {
                    break;
                }
            }
            out = out.add(&acc.scale(d));
        }
        out
    }
}

impl<C: Ring> SkewBy<C, Monomial<C>> for Schur<C> {
    fn skew_by(&self, g: &Monomial<C>) -> Self {
        SkewBy::skew_by(self, &g.to_schur())
    }
}

impl<C: Ring> SkewBy<C, Forgotten<C>> for Schur<C> {
    fn skew_by(&self, g: &Forgotten<C>) -> Self {
        SkewBy::skew_by(self, &g.to_schur())
    }
}

/// Remove a horizontal `r`-strip from every term.
///
/// Exposed inside the crate because this *is* `h_r^⊥` on a Schur expansion, and
/// the Hall–Littlewood recursion calls it directly for a single `r` rather than
/// going through [`SkewBy`] with a one-part `h`.
pub(crate) fn strip_off<C: Ring>(f: &Schur<C>, r: u32) -> Schur<C> {
    let mut out = Schur::zero();
    let mut found: Vec<Partition> = Vec::new();
    for (lambda, c) in f.terms() {
        crate::interrupt::poll();
        found.clear();
        let mut cur = lambda.parts().to_vec();
        remove_horizontal(lambda.parts(), 0, r, &mut cur, &mut |nu| {
            found.push(Partition::from_sorted(
                nu.iter().copied().filter(|&v| v > 0).collect(),
            ));
        });
        for nu in found.drain(..) {
            out.add_term(nu, c.clone());
        }
    }
    out
}

/// Every ν ⊆ λ with λ/ν a horizontal strip of exactly `left` cells.
///
/// λ/ν is a horizontal strip iff λ_1 ≥ ν_1 ≥ λ_2 ≥ ν_2 ≥ …, so row i is free in
/// [λ_{i+1}, λ_i] independently of the others — the only coupling is the cell
/// budget. Rows i onward can give up Σ_{j≥i}(λ_j − λ_{j+1}) = λ_i cells in
/// total, which telescopes, so the capacity prune is a single comparison.
fn remove_horizontal(
    lam: &[u32],
    i: usize,
    left: u32,
    cur: &mut Vec<u32>,
    emit: &mut impl FnMut(&[u32]),
) {
    if i == lam.len() {
        if left == 0 {
            emit(cur);
        }
        return;
    }
    if left > lam[i] {
        return; // more cells than every remaining row can supply
    }
    let floor = if i + 1 < lam.len() { lam[i + 1] } else { 0 };
    let lo = floor.max(lam[i] - left);
    for v in lo..=lam[i] {
        cur[i] = v;
        remove_horizontal(lam, i + 1, left - (lam[i] - v), cur, emit);
    }
    cur[i] = lam[i];
}

/// The coproduct Δ(f) = Σ_{μ,ν} c^λ_{μν} s_μ ⊗ s_ν, extended linearly.
///
/// Computed as Δ(s_λ) = Σ_{μ⊆λ} s_μ ⊗ s_{λ/μ}, so one traversal per μ yields
/// every ν that occurs.
pub fn coproduct<C: Ring>(f: &Schur<C>) -> SymTensor<C> {
    let mut out = SymTensor::zero();
    for (lambda, c) in f.terms() {
        crate::interrupt::poll();
        for k in 0..=lambda.size() {
            for mu in partitions_cached(k).iter() {
                // Only μ ⊆ λ contribute; the rest have s_{λ/μ} = 0.
                if !lambda.contains(mu) {
                    continue;
                }
                for (nu, coeff) in crate::skew_lr::expand_skew_shared(lambda, mu).iter() {
                    out.add_term((mu.clone(), nu.clone()), C::from_u128(*coeff).mul(c));
                }
            }
        }
    }
    out
}

/// The counit ε(f): project onto the degree-0 component, ε(s_λ) = δ_{λ,∅}.
pub fn counit<C: Ring>(f: &Schur<C>) -> C {
    f.coeff(&Partition::default())
}

/// The antipode S(f), with S(s_λ) = (−1)^{|λ|} s_{λ'}.
pub fn antipode<C: Ring>(f: &Schur<C>) -> Schur<C> {
    let mut out = Schur::zero();
    for (lambda, c) in f.terms() {
        crate::interrupt::poll();
        let coeff = if lambda.size() % 2 == 1 {
            c.neg()
        } else {
            c.clone()
        };
        out.add_term(lambda.conjugate(), coeff);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn s(v: &[u32]) -> Schur<i64> {
        Schur::monomial(part(v), 1)
    }

    #[test]
    fn skew_schur_known_value() {
        // s_{(2,1)/(1)} = s_2 + s_{11} (two disconnected cells ≅ s_1·s_1).
        let sk: Schur<i64> = skew_schur(&part(&[2, 1]), &part(&[1]));
        assert_eq!(sk.coeff(&part(&[2])), 1);
        assert_eq!(sk.coeff(&part(&[1, 1])), 1);
        assert_eq!(sk.terms().len(), 2);
    }

    /// The **defining** property: ⟨g^⊥ f, h⟩ = ⟨f, g·h⟩, for every f, g, h in a
    /// range of degrees.
    ///
    /// This is the right test for skewing because it never mentions how any of
    /// it is computed — it pins the adjointness that *is* the definition, and
    /// the two sides share no code: the left runs the skew traversal and a
    /// coefficient dot product, the right runs Littlewood–Richardson.
    #[test]
    fn skewing_is_adjoint_to_multiplication() {
        for df in 0..=6u32 {
            for dg in 0..=df {
                for lambda in partitions_cached(df).iter() {
                    let f: Schur<i64> = Schur::monomial(lambda.clone(), 1);
                    for mu in partitions_cached(dg).iter() {
                        let g: Schur<i64> = Schur::monomial(mu.clone(), 1);
                        let skewed = SkewBy::skew_by(&f, &g);
                        // h ranges over a basis of the relevant degree, which is
                        // enough: both sides are linear in h.
                        for nu in partitions_cached(df - dg).iter() {
                            let h: Schur<i64> = Schur::monomial(nu.clone(), 1);
                            let lhs = crate::ops::hall::<i64, _, _>(&skewed, &h);
                            let rhs = crate::ops::hall::<i64, _, _>(&f, &g.mul(&h));
                            assert_eq!(lhs, rhs, "<{mu}^perp s_{lambda}, s_{nu}>");
                        }
                    }
                }
            }
        }
    }

    /// Each native rule — Pieri for h, dual Pieri for e, Murnaghan–Nakayama for
    /// p — must agree with expanding g into Schur and running the LR route.
    ///
    /// The point of the fast paths is that they avoid LR entirely, so this is
    /// the check that they are shortcuts and not different operations. Note the
    /// coefficient ring: the *comparison* has to be rational, because expanding
    /// p_μ into Schur divides by z_μ — while the p rule under test needs only a
    /// ring, and is exercised over `i64` in
    /// `skewing_by_power_sums_stays_integral`.
    #[test]
    fn native_skew_rules_match_the_lr_route() {
        for df in 0..=6u32 {
            for dg in 0..=df {
                for lambda in partitions_cached(df).iter() {
                    let f: Schur<Rational> = Schur::monomial(lambda.clone(), Rational::from_int(1));
                    for mu in partitions_cached(dg).iter() {
                        let one = Rational::from_int(1);

                        let gh: Homogeneous<Rational> = Homogeneous::monomial(mu.clone(), one);
                        assert_eq!(
                            SkewBy::skew_by(&f, &gh),
                            SkewBy::skew_by(&f, &gh.to_schur()),
                            "h_{mu}^perp s_{lambda}"
                        );

                        let ge: Elementary<Rational> = Elementary::monomial(mu.clone(), one);
                        assert_eq!(
                            SkewBy::skew_by(&f, &ge),
                            SkewBy::skew_by(&f, &ge.to_schur()),
                            "e_{mu}^perp s_{lambda}"
                        );

                        let gp: PowerSum<Rational> = PowerSum::monomial(mu.clone(), one);
                        assert_eq!(
                            SkewBy::skew_by(&f, &gp),
                            SkewBy::skew_by(&f, &gp.to_schur()),
                            "p_{mu}^perp s_{lambda}"
                        );

                        let gm: Monomial<Rational> = Monomial::monomial(mu.clone(), one);
                        assert_eq!(
                            SkewBy::skew_by(&f, &gm),
                            SkewBy::skew_by(&f, &gm.to_schur()),
                            "m_{mu}^perp s_{lambda}"
                        );
                    }
                }
            }
        }
    }

    /// Skewing by a power sum never divides, so it works over ℤ — a thing the
    /// generic route cannot do at all, since expanding p_μ into Schur needs
    /// 1/z_μ. Checked against the rational route's values.
    #[test]
    fn skewing_by_power_sums_stays_integral() {
        for df in 1..=6u32 {
            for dg in 1..=df {
                for lambda in partitions_cached(df).iter() {
                    for mu in partitions_cached(dg).iter() {
                        let f: Schur<i64> = Schur::monomial(lambda.clone(), 1);
                        let g: PowerSum<i64> = PowerSum::monomial(mu.clone(), 1);
                        let got = SkewBy::skew_by(&f, &g);

                        let fr: Schur<Rational> =
                            Schur::monomial(lambda.clone(), Rational::from_int(1));
                        let gr: PowerSum<Rational> =
                            PowerSum::monomial(mu.clone(), Rational::from_int(1));
                        let want = SkewBy::skew_by(&fr, &gr.to_schur());

                        assert_eq!(got.terms().len(), want.terms().len(), "{mu} on {lambda}");
                        for (nu, c) in got.terms() {
                            assert_eq!(
                                Rational::from_int((*c).into()),
                                want.coeff(nu),
                                "p_{mu}^perp s_{lambda} at {nu}"
                            );
                        }
                    }
                }
            }
        }
    }

    /// Skewing by s_μ must reproduce `skew_schur`, and h_1^⊥ must remove
    /// exactly one corner box.
    #[test]
    fn skew_by_reproduces_the_classical_special_cases() {
        let lam = part(&[3, 2, 1]);
        let f: Schur<i64> = Schur::monomial(lam.clone(), 1);
        for mu in partitions_cached(2).iter() {
            let g: Schur<i64> = Schur::monomial(mu.clone(), 1);
            assert_eq!(SkewBy::skew_by(&f, &g), skew_schur(&lam, mu), "s_{mu}^perp");
        }

        // λ = (3,2,1) has corners in all three rows.
        let h1: Homogeneous<i64> = Homogeneous::monomial(part(&[1]), 1);
        let one_off = SkewBy::skew_by(&f, &h1);
        assert_eq!(one_off.coeff(&part(&[2, 2, 1])), 1);
        assert_eq!(one_off.coeff(&part(&[3, 1, 1])), 1);
        assert_eq!(one_off.coeff(&part(&[3, 2])), 1);
        assert_eq!(one_off.terms().len(), 3);
    }

    /// Skewing lowers degree by |g|, and by exactly that: g^⊥ kills anything of
    /// lower degree.
    #[test]
    fn skewing_is_degree_lowering() {
        let f: Schur<i64> = Schur::monomial(part(&[2, 1]), 1);
        let big: Homogeneous<i64> = Homogeneous::monomial(part(&[4]), 1);
        assert!(SkewBy::skew_by(&f, &big).is_zero());
        let exact: Homogeneous<i64> = Homogeneous::monomial(part(&[2, 1]), 1);
        let r = SkewBy::skew_by(&f, &exact);
        assert_eq!(r.degree(), Some(0));
    }

    #[test]
    fn coproduct_known_values() {
        // Δ(s_1) = s_∅⊗s_1 + s_1⊗s_∅.
        let d1 = coproduct(&s(&[1]));
        assert_eq!(d1.coeff(&part(&[]), &part(&[1])), 1);
        assert_eq!(d1.coeff(&part(&[1]), &part(&[])), 1);
        assert_eq!(d1.terms().len(), 2);

        // Δ(s_2) = s_∅⊗s_2 + s_1⊗s_1 + s_2⊗s_∅.
        let d2 = coproduct(&s(&[2]));
        assert_eq!(d2.coeff(&part(&[]), &part(&[2])), 1);
        assert_eq!(d2.coeff(&part(&[1]), &part(&[1])), 1);
        assert_eq!(d2.coeff(&part(&[2]), &part(&[])), 1);
        assert_eq!(d2.terms().len(), 3);
    }

    #[test]
    fn antipode_and_counit() {
        assert_eq!(antipode(&s(&[2])).coeff(&part(&[1, 1])), 1); // (−1)² s_{11}
        assert_eq!(antipode(&s(&[1])).coeff(&part(&[1])), -1); // (−1)¹ s_1
        assert_eq!(counit(&s(&[])), 1);
        assert_eq!(counit(&s(&[2])), 0);
    }

    /// m ∘ (S ⊗ id) ∘ Δ, which the Hopf axiom says equals ε(·)·1.
    fn antipode_axiom(f: &Schur<i64>) -> Schur<i64> {
        let d = coproduct(f);
        let mut acc = Schur::zero();
        for ((mu, nu), c) in d.terms() {
            let left = antipode(&Schur::monomial(mu.clone(), 1));
            let prod = left.mul(&Schur::monomial(nu.clone(), 1)).scale(c);
            acc = acc.add(&prod);
        }
        acc
    }

    #[test]
    fn hopf_antipode_axiom_holds() {
        // For |λ| ≥ 1 the axiom gives 0; for λ = ∅ it gives the unit s_∅.
        for parts in [&[1][..], &[2], &[2, 1], &[3, 1]] {
            assert!(antipode_axiom(&s(parts)).is_zero(), "axiom at {:?}", parts);
        }
        let unit = antipode_axiom(&s(&[]));
        assert_eq!(unit.coeff(&part(&[])), 1);
        assert_eq!(unit.terms().len(), 1);
    }
}
