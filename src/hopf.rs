//! The Hopf-algebra structure of Sym: skew Schur functions, comultiplication,
//! counit, and antipode.
//!
//! All of it reduces to the Littlewood–Richardson coefficients already computed
//! by [`NaiveLr`]:
//!
//! - skew Schur:   s_{λ/μ} = Σ_ν c^λ_{μν} s_ν
//! - coproduct:    Δ(s_λ)  = Σ_{μ,ν} c^λ_{μν} s_μ ⊗ s_ν
//! - counit:       ε(s_λ)  = δ_{λ,∅}
//! - antipode:     S(s_λ)  = (−1)^{|λ|} s_{λ'}

use crate::coeff::Ring;
use crate::lr::{LrBackend, NaiveLr};
use crate::memo::partitions_cached;
use crate::partition::Partition;
use crate::sym::{Schur, SymFn};
use std::collections::BTreeMap;

/// An element of Sym ⊗ Sym in the Schur basis: a formal `C`-combination of pairs
/// of partitions (s_μ ⊗ s_ν). The codomain of the coproduct.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SymTensor<C: Ring>(BTreeMap<(Partition, Partition), C>);

impl<C: Ring> SymTensor<C> {
    pub fn zero() -> Self {
        SymTensor(BTreeMap::new())
    }

    pub fn terms(&self) -> &BTreeMap<(Partition, Partition), C> {
        &self.0
    }

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

    pub fn add(&self, other: &Self) -> Self {
        let mut out = self.clone();
        for (k, c) in &other.0 {
            out.add_term(k.clone(), c.clone());
        }
        out
    }
}

/// The skew Schur function s_{λ/μ} = Σ_ν c^λ_{μν} s_ν.
pub fn skew_schur<C: Ring>(lambda: &Partition, mu: &Partition) -> Schur<C> {
    let mut out = Schur::zero();
    if !lambda.contains(mu) {
        return out;
    }
    // One traversal of the shape yields every ν with a nonzero coefficient.
    // (This used to sweep all p(n) partitions, running a full LR backtrack per
    // candidate — the same answer for orders of magnitude more work.)
    for (nu, c) in crate::skew_lr::expand_skew(lambda, mu) {
        out.add_term(nu, C::from_u128(c));
    }
    out
}

/// The coproduct Δ(f) = Σ_{μ,ν} c^λ_{μν} s_μ ⊗ s_ν, extended linearly.
pub fn coproduct<C: Ring>(f: &Schur<C>) -> SymTensor<C> {
    let mut out = SymTensor::zero();
    for (lambda, c) in f.terms() {
        let n = lambda.size();
        for k in 0..=n {
            for mu in partitions_cached(k).iter() {
                if !lambda.contains(mu) {
                    continue;
                }
                for nu in partitions_cached(n - k).iter() {
                    let coeff = NaiveLr.lr_coeff(lambda, mu, nu);
                    if coeff != 0 {
                        out.add_term(
                            (mu.clone(), nu.clone()),
                            C::from_u128(coeff).mul(c),
                        );
                    }
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
        let coeff = if lambda.size() % 2 == 1 { c.neg() } else { c.clone() };
        out.add_term(lambda.conjugate(), coeff);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

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
