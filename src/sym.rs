//! Symmetric functions as **real types**, one per basis, unified by a trait.
//!
//! This is the deliberate departure from Symmetrica's single untyped `OP`
//! object (and from a single tagged `enum`): `Schur<C>`, `PowerSum<C>`, and
//! `Monomial<C>` are *distinct types*. You cannot add a Schur element to a
//! power-sum element, or feed one into the other's multiplication, because the
//! partitions mean different things in different bases and the compiler now
//! enforces that. Multiplication is basis-specific and dispatched statically:
//! `Schur::mul` goes through the Littlewood–Richardson backend, while
//! `PowerSum::mul` is a trivial multiset union — different rules, chosen at
//! compile time, no runtime tag.
//!
//! The shared *linear* structure (add, scale, coefficients, degree) lives on the
//! [`SymFn`] trait so it is written once and reused by every basis.

use crate::coeff::Ring;
use crate::lr::LrBackend;
use crate::partition::Partition;
use crate::strip_lr::AutoLr;
use std::collections::BTreeMap;

/// The linear structure common to every basis of the ring of symmetric
/// functions: a finite formal `C`-combination of partitions. `BTreeMap` keeps a
/// deterministic term order (nice for display and tests); a `HashMap` is the
/// obvious swap for raw throughput later.
pub trait SymFn<C: Ring>: Sized {
    /// Single-character basis symbol used when printing (`s`, `p`, `m`, …).
    const SYMBOL: &'static str;

    /// The terms, keyed by partition. Explicit zeros are never stored, so the
    /// map is empty exactly when the element is.
    fn terms(&self) -> &BTreeMap<Partition, C>;
    /// The terms, mutably. Inserting a zero coefficient through this breaks the
    /// invariant every other method reads — [`add_term`](Self::add_term) is the
    /// accumulating form that maintains it.
    fn terms_mut(&mut self) -> &mut BTreeMap<Partition, C>;
    /// Wrap a term map that is already free of explicit zeros.
    fn from_terms(terms: BTreeMap<Partition, C>) -> Self;

    /// The additive identity 0.
    fn zero() -> Self {
        Self::from_terms(BTreeMap::new())
    }

    /// A single basis element with coefficient `c` (dropped if `c` is zero).
    fn monomial(p: Partition, c: C) -> Self {
        let mut t = BTreeMap::new();
        if !c.is_zero() {
            t.insert(p, c);
        }
        Self::from_terms(t)
    }

    /// Coefficient of a given partition (zero if absent).
    fn coeff(&self, p: &Partition) -> C {
        self.terms().get(p).cloned().unwrap_or_else(C::zero)
    }

    /// Whether this is 0.
    fn is_zero(&self) -> bool {
        self.terms().is_empty()
    }

    /// The degree |λ| of the top term, or `None` for 0.
    fn degree(&self) -> Option<u32> {
        self.terms().keys().map(Partition::size).max()
    }

    /// `self += c · [p]`, keeping the map free of explicit zeros.
    fn add_term(&mut self, p: Partition, c: C) {
        if c.is_zero() {
            return;
        }
        let map = self.terms_mut();
        let now_zero = {
            let e = map.entry(p.clone()).or_insert_with(C::zero);
            e.add_assign(&c);
            e.is_zero()
        };
        if now_zero {
            map.remove(&p);
        }
    }

    /// Formal sum `self + other`.
    fn add(&self, other: &Self) -> Self {
        let mut out = Self::from_terms(self.terms().clone());
        for (p, c) in other.terms() {
            out.add_term(p.clone(), c.clone());
        }
        out
    }

    /// Scalar multiple `c · self`.
    fn scale(&self, c: &C) -> Self {
        let mut t = BTreeMap::new();
        if !c.is_zero() {
            for (p, cc) in self.terms() {
                let v = cc.mul(c);
                if !v.is_zero() {
                    t.insert(p.clone(), v);
                }
            }
        }
        Self::from_terms(t)
    }

    /// Additive inverse −self.
    fn neg(&self) -> Self {
        self.scale(&C::from_i64(-1))
    }

    /// Difference self − other.
    fn sub(&self, other: &Self) -> Self {
        self.add(&other.neg())
    }

    /// Human-readable form, e.g. `s[2] + s[1,1]` or `2*s[3,2,1]`.
    fn format(&self) -> String {
        if self.is_zero() {
            return "0".to_string();
        }
        let one = C::one();
        let mut out = Vec::new();
        for (p, c) in self.terms() {
            let prefix = if *c == one {
                String::new()
            } else {
                format!("{c:?}*")
            };
            out.push(format!("{prefix}{}{p}", Self::SYMBOL));
        }
        out.join(" + ")
    }
}

/// Define a basis newtype `$name` (symbol `$sym`) plus its `SymFn` and `Display`
/// impls. Keeps the per-basis boilerplate in one place.
macro_rules! basis {
    ($(#[$m:meta])* $name:ident, $sym:literal) => {
        $(#[$m])*
        #[derive(Clone, PartialEq, Eq, Debug)]
        pub struct $name<C: Ring>(BTreeMap<Partition, C>);

        impl<C: Ring> SymFn<C> for $name<C> {
            const SYMBOL: &'static str = $sym;
            fn terms(&self) -> &BTreeMap<Partition, C> { &self.0 }
            fn terms_mut(&mut self) -> &mut BTreeMap<Partition, C> { &mut self.0 }
            fn from_terms(terms: BTreeMap<Partition, C>) -> Self { $name(terms) }
        }

        impl<C: Ring> core::fmt::Display for $name<C> {
            fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
                f.write_str(&self.format())
            }
        }
    };
}

basis!(
    /// Symmetric functions in the **Schur** basis {s_λ}.
    Schur, "s"
);
basis!(
    /// Symmetric functions in the **power-sum** basis {p_λ}.
    PowerSum, "p"
);
basis!(
    /// Symmetric functions in the **monomial** basis {m_λ}.
    Monomial, "m"
);
basis!(
    /// Symmetric functions in the **elementary** basis {e_λ}.
    Elementary, "e"
);
basis!(
    /// Symmetric functions in the **complete homogeneous** basis {h_λ}.
    Homogeneous, "h"
);
basis!(
    /// Symmetric functions in the **forgotten** basis {f_λ}, defined by
    /// f_λ = ω(m_λ).
    ///
    /// This is the sixth and last of Macdonald's classical bases (I.2), and the
    /// only one with no independent combinatorial description — it is *defined*
    /// as the image of the monomial basis under ω, which is where the name comes
    /// from. Because ω is an involutive algebra automorphism, {f_λ} is a basis
    /// the moment {m_λ} is, and every conversion involving it is a conversion
    /// involving m with an ω on one side; see `convert`.
    ///
    /// Note it is deliberately **not** a [`SymAlgebra`]: unlike p, e, and h it is
    /// not multiplicative, so f_μ · f_ν is not f_{μ ∪ ν} and there is no
    /// index-concatenation product to expose.
    Forgotten, "f"
);
basis!(
    /// Symmetric functions in the Orellana–Zabrocki **irreducible character**
    /// basis {s̃_λ} — the `st` basis.
    ///
    /// `s̃_λ` is the unique symmetric function whose evaluation at the
    /// eigenvalues of a permutation matrix of cycle type γ ⊢ n is the *symmetric
    /// group* character `χ^{(n−|λ|,λ)}(γ)`, for every n ≥ |λ| + λ₁ (OZ Thm 1).
    /// Schur functions are the characters of GLₙ; these are the characters of
    /// Sₙ, living in the same ring. The long first row is implicit, which is why
    /// λ indexes a shape of *any* large size.
    ///
    /// Two consequences shape the type. It is **inhomogeneous** — `s̃_λ` has
    /// components in every degree from 0 to |λ| — which is why nothing in
    /// [`SymFn`] may assume a single degree; [`SymFn::degree`] returning the max
    /// over terms is already the right answer. And its structure constants are
    /// the **reduced (stable) Kronecker coefficients** (OZ Thm 7), so
    /// [`St::mul`](crate::character_basis) is a Kronecker engine wearing an
    /// ordinary product's clothes. See `character_basis` for both.
    St, "st"
);
basis!(
    /// Symmetric functions in the Orellana–Zabrocki **induced trivial
    /// character** basis {h̃_λ} — the `ht` basis.
    ///
    /// `h̃_λ` evaluates to the character of the trivial representation induced
    /// from a Young subgroup, i.e. of the permutation module `M^{(n−|λ|,λ)}`
    /// (OZ Def 4). It is to [`St`] what [`Homogeneous`] is to [`Schur`], and it
    /// earns its place the same way: the transition between the two is a Kostka
    /// matrix, and its own product is a sum over integer matrices. That is an
    /// independent second route to the reduced Kronecker coefficients, and it
    /// is how the first one gets checked at sizes no other package can reach.
    Ht, "ht"
);

/// The product shared by every *multiplicative* basis: since p_λ, e_λ, h_λ are
/// each defined as a product of one-part generators, x_λ · x_μ = x_{λ ∪ μ} (the
/// multiset union of parts). One implementation, reused by p, e, and h.
fn multiplicative_product<C: Ring, S: SymFn<C>>(a: &S, b: &S) -> S {
    let mut out = S::zero();
    for (mu, cmu) in a.terms() {
        for (nu, cnu) in b.terms() {
            let mut parts = mu.parts().to_vec();
            parts.extend_from_slice(nu.parts());
            out.add_term(Partition::new(parts), cmu.mul(cnu));
        }
    }
    out
}

impl<C: Ring> Schur<C> {
    /// Product in the Schur basis via an explicit LR backend.
    ///
    /// Structure constants arrive as `u128` and are injected through
    /// [`Ring::from_u128`], so a bignum coefficient type carries them exactly;
    /// only fixed-width types lose range, and that is inherent to them.
    pub fn mul_with<B: LrBackend>(&self, other: &Self, backend: &B) -> Self {
        let mut out = Self::zero();
        for (mu, cmu) in self.terms() {
            for (nu, cnu) in other.terms() {
                let cprod = cmu.mul(cnu);
                for (lambda, k) in backend.schur_product(mu, nu) {
                    let term = C::from_u128(k).mul(&cprod);
                    out.add_term(lambda, term);
                }
            }
        }
        out
    }

    /// Product in the Schur basis using the default backend,
    /// [`crate::strip_lr::AutoLr`] — currently a single traversal of the
    /// product's skew shape, binned by content.
    /// `NaiveLr` remains available as the reference oracle.
    pub fn mul(&self, other: &Self) -> Self {
        self.mul_with(other, &AutoLr)
    }
}

impl<C: Ring> PowerSum<C> {
    /// Product in the power-sum basis: p_μ · p_ν = p_{μ ∪ ν}. Note how different
    /// — and how much simpler — the rule is from the Schur case, yet it is
    /// selected statically by the type.
    pub fn mul(&self, other: &Self) -> Self {
        multiplicative_product(self, other)
    }
}

impl<C: Ring> Elementary<C> {
    /// Product in the elementary basis: e_μ · e_ν = e_{μ ∪ ν}.
    pub fn mul(&self, other: &Self) -> Self {
        multiplicative_product(self, other)
    }
}

impl<C: Ring> Homogeneous<C> {
    /// Product in the complete-homogeneous basis: h_μ · h_ν = h_{μ ∪ ν}.
    pub fn mul(&self, other: &Self) -> Self {
        multiplicative_product(self, other)
    }
}

impl<C: Ring> Monomial<C> {
    /// Product in the monomial basis: `m_μ · m_ν = Σ_λ c_λ m_λ`, where `c_λ` is
    /// the number of ways to write the exponent vector of λ as `α + β` with α a
    /// distinct rearrangement of μ and β one of ν.
    ///
    /// Not a multiset union — that rule belongs to the *multiplicative* bases
    /// `p`, `e`, `h`, and `m` is not one of them: `m_1 · m_1 = 2·m_{11} + m_2`,
    /// where a union would give `m_{11}` alone. The coefficients are
    /// non-negative but not always 0 or 1, which is the fact the doctest pins.
    ///
    /// Cost is bounded by the pairs of rearrangements the recursion visits,
    /// which is far below `R(μ)·R(ν)`. The slots are chosen jointly and the
    /// weakly-decreasing constraint prunes on the way down, so a pair whose sum
    /// is not a partition is abandoned at the first slot that proves it rather
    /// than after both rearrangements are complete.
    ///
    /// This is the operation behind Sage's product in the `m` basis, whose
    /// backend is Symmetrica's `mult_monomial_monomial`.
    ///
    /// # Examples
    ///
    /// ```
    /// use symfn::{Monomial, Partition, SymFn};
    /// let m1: Monomial<i64> = Monomial::monomial(Partition::new([1]), 1);
    /// let sq = m1.mul(&m1);
    /// assert_eq!(sq.coeff(&Partition::new([1, 1])), 2);
    /// assert_eq!(sq.coeff(&Partition::new([2])), 1);
    /// ```
    pub fn mul(&self, other: &Self) -> Self {
        let mut out = Self::zero();
        for (mu, cmu) in self.terms() {
            for (nu, cnu) in other.terms() {
                let c = cmu.mul(cnu);
                if c.is_zero() {
                    continue;
                }
                overlays(mu, nu, &mut |lambda: &[u32]| {
                    out.add_term(Partition::new(lambda.iter().copied()), c.clone());
                });
            }
        }
        out
    }
}

/// Every λ obtained as the slotwise sum of a distinct rearrangement of μ and one
/// of ν, emitted once per *pair* — so the multiplicity with which λ arrives is
/// the structure constant.
///
/// `ℓ(μ) + ℓ(ν)` slots are enough and never too few: a nonzero slot needs a part
/// from at least one side, so no λ in the product has more rows than that. The
/// bound matters because the answer is read off a *fixed* exponent vector — the
/// weakly decreasing one — and that is only legitimate while every λ that
/// occurs still fits in the alphabet.
fn overlays(mu: &Partition, nu: &Partition, emit: &mut impl FnMut(&[u32])) {
    let n = mu.len() + nu.len();
    let mut a = crate::eval::multiplicities(mu);
    let mut b = crate::eval::multiplicities(nu);
    let mut sum = vec![0u32; n];
    walk(
        &mut sum,
        0,
        u32::MAX,
        &mut a,
        &mut b,
        mu.len(),
        nu.len(),
        emit,
    );
}

/// Fill slot `slot` with `α_slot + β_slot`, keeping the running vector weakly
/// decreasing, and recurse.
#[allow(clippy::too_many_arguments)]
fn walk(
    sum: &mut [u32],
    slot: usize,
    prev: u32,
    a: &mut [(u32, u32)],
    b: &mut [(u32, u32)],
    left_a: usize,
    left_b: usize,
    emit: &mut impl FnMut(&[u32]),
) {
    if left_a == 0 && left_b == 0 {
        // Every later slot would be 0, so the vector ends here.
        emit(&sum[..slot]);
        return;
    }
    // Each side still needs a slot per unplaced part.
    if sum.len() - slot < left_a.max(left_b) {
        return;
    }
    // Index `len` means "this side contributes nothing to this slot". Both
    // sides declining would make the slot 0 with parts still to place, and a 0
    // forces every later slot to 0 — so that branch can never complete and is
    // skipped rather than explored and abandoned.
    let (na, nb) = (a.len(), b.len());
    for i in 0..=na {
        let va = if i == na {
            0
        } else if a[i].1 == 0 {
            continue;
        } else {
            a[i].0
        };
        for j in 0..=nb {
            let vb = if j == nb {
                0
            } else if b[j].1 == 0 {
                continue;
            } else {
                b[j].0
            };
            if i == na && j == nb {
                continue;
            }
            if va + vb > prev {
                continue;
            }
            if i < na {
                a[i].1 -= 1;
            }
            if j < nb {
                b[j].1 -= 1;
            }
            sum[slot] = va + vb;
            walk(
                sum,
                slot + 1,
                va + vb,
                a,
                b,
                left_a - usize::from(i < na),
                left_b - usize::from(j < nb),
                emit,
            );
            if i < na {
                a[i].1 += 1;
            }
            if j < nb {
                b[j].1 += 1;
            }
        }
    }
    sum[slot] = 0;
}

/// A basis that is also a *ring* under its own multiplication, with a unit (the
/// empty partition, = 1). This lets generic code — determinants in particular —
/// multiply basis elements without knowing which basis it holds.
pub trait SymAlgebra<C: Ring>: SymFn<C> {
    /// The multiplicative unit 1 = x_∅.
    fn unit() -> Self {
        Self::monomial(Partition::default(), C::one())
    }
    /// The basis product (delegates to each basis's inherent `mul`).
    fn times(&self, other: &Self) -> Self;
}

impl<C: Ring> SymAlgebra<C> for Homogeneous<C> {
    fn times(&self, o: &Self) -> Self {
        self.mul(o)
    }
}
impl<C: Ring> SymAlgebra<C> for Elementary<C> {
    fn times(&self, o: &Self) -> Self {
        self.mul(o)
    }
}
impl<C: Ring> SymAlgebra<C> for PowerSum<C> {
    fn times(&self, o: &Self) -> Self {
        self.mul(o)
    }
}
impl<C: Ring> SymAlgebra<C> for Schur<C> {
    fn times(&self, o: &Self) -> Self {
        self.mul(o)
    }
}
impl<C: Ring> SymAlgebra<C> for Monomial<C> {
    fn times(&self, o: &Self) -> Self {
        self.mul(o)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn s(parts: &[u32], c: i64) -> Schur<i64> {
        Schur::monomial(Partition::new(parts.iter().copied()), c)
    }

    /// The direct overlay rule and the `m → s → m` route must agree.
    ///
    /// They share no step: one counts pairs of rearrangements slot by slot, the
    /// other inverts the Kostka matrix, runs Littlewood–Richardson, and applies
    /// the Kostka matrix again. Agreement across a whole degree is the check
    /// that the joint recursion's pruning drops only branches that could not
    /// have completed.
    #[test]
    fn monomial_product_agrees_with_the_schur_route() {
        use crate::memo::partitions_cached;
        for a in 0..=4u32 {
            for b in 0..=4u32 {
                for mu in partitions_cached(a).iter() {
                    for nu in partitions_cached(b).iter() {
                        let x: Monomial<i64> = Monomial::monomial(mu.clone(), 1);
                        let y: Monomial<i64> = Monomial::monomial(nu.clone(), 1);
                        assert_eq!(x.mul(&y), x.mul_via_schur(&y), "m_{mu} · m_{nu}");
                    }
                }
            }
        }
    }

    /// The product is what evaluation says it is: `(fg)(x) = f(x)·g(x)` at a
    /// concrete alphabet, for any alphabet long enough to separate the terms.
    ///
    /// This pins the structure constants against something outside the basis
    /// entirely, so a rule that were self-consistently wrong — double-counting
    /// equal parts, say — would still fail here.
    #[test]
    fn monomial_product_evaluates_as_a_product() {
        use crate::memo::partitions_cached;
        let xs = [2i64, -1, 3, 1, -2, 4];
        for a in 0..=4u32 {
            for b in 0..=4u32 {
                for mu in partitions_cached(a).iter() {
                    for nu in partitions_cached(b).iter() {
                        let x: Monomial<i64> = Monomial::monomial(mu.clone(), 1);
                        let y: Monomial<i64> = Monomial::monomial(nu.clone(), 1);
                        assert_eq!(
                            x.mul(&y).eval(&xs),
                            x.eval(&xs) * y.eval(&xs),
                            "m_{mu} · m_{nu} at the alphabet"
                        );
                    }
                }
            }
        }
    }

    /// `m_1 · m_1 = 2·m_{11} + m_2`, the smallest value that distinguishes the
    /// monomial product from the multiset union p, e and h use.
    #[test]
    fn monomial_product_is_not_a_multiset_union() {
        let m1: Monomial<i64> = Monomial::monomial(Partition::new([1]), 1);
        let sq = m1.mul(&m1);
        assert_eq!(sq.coeff(&Partition::new([1, 1])), 2);
        assert_eq!(sq.coeff(&Partition::new([2])), 1);
        assert_eq!(sq.terms().len(), 2);
    }

    #[test]
    fn schur_product_matches_known_expansion() {
        // s_2 · s_1 = s_3 + s_{21}
        let prod = s(&[2], 1).mul(&s(&[1], 1));
        assert_eq!(prod.coeff(&Partition::new([3])), 1);
        assert_eq!(prod.coeff(&Partition::new([2, 1])), 1);
        assert_eq!(prod.terms().len(), 2);
    }

    #[test]
    fn schur_product_is_commutative() {
        let a = s(&[2, 1], 1);
        let b = s(&[2], 1).add(&s(&[1, 1], 1));
        assert_eq!(a.mul(&b), b.mul(&a));
    }

    #[test]
    fn schur_multiplicity_two_survives_scaling() {
        // 3·s_{21} · s_{21} has coefficient 6 on s_{321}
        let prod = s(&[2, 1], 3).mul(&s(&[2, 1], 1));
        assert_eq!(prod.coeff(&Partition::new([3, 2, 1])), 6);
    }

    #[test]
    fn power_sum_product_is_multiset_union() {
        // p_2 · p_{31} = p_{321}
        let a: PowerSum<i64> = PowerSum::monomial(Partition::new([2]), 1);
        let b: PowerSum<i64> = PowerSum::monomial(Partition::new([3, 1]), 1);
        let prod = a.mul(&b);
        assert_eq!(prod.coeff(&Partition::new([3, 2, 1])), 1);
        assert_eq!(prod.terms().len(), 1);
    }

    #[test]
    fn elementary_and_homogeneous_products_are_multiplicative() {
        // e_2 · e_1 = e_{21}, h_{11} · h_2 = h_{211}
        let e: Elementary<i64> = Elementary::monomial(Partition::new([2]), 1);
        let e1: Elementary<i64> = Elementary::monomial(Partition::new([1]), 1);
        assert_eq!(e.mul(&e1).coeff(&Partition::new([2, 1])), 1);

        let h: Homogeneous<i64> = Homogeneous::monomial(Partition::new([1, 1]), 1);
        let h2: Homogeneous<i64> = Homogeneous::monomial(Partition::new([2]), 1);
        let prod = h.mul(&h2);
        assert_eq!(prod.coeff(&Partition::new([2, 1, 1])), 1);
        assert_eq!(prod.terms().len(), 1);
    }

    #[test]
    fn display_omits_unit_coefficient() {
        let f = s(&[2], 1).add(&s(&[1, 1], 3));
        // BTreeMap order: [1,1] < [2] lexicographically? parts vecs: [1,1] vs [2] -> [1,1] < [2]
        assert_eq!(f.to_string(), "3*s[1,1] + s[2]");
    }
}
