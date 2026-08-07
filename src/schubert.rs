//! Schubert polynomials `S_w`, indexed by [`Perm`].
//!
//! `{S_w : w ∈ S_∞}` is a ℤ-basis of the *whole* polynomial ring
//! `ℤ[x₁, x₂, …]`, not of Sym — which is why this type deliberately does not
//! implement [`SymFn`](crate::sym::SymFn): that trait is `Partition`-indexed
//! and everything generic over it assumes symmetric functions.
//! [`Forgotten`](crate::sym::Forgotten) is the precedent for declining a trait
//! rather than half-satisfying it. The handful of shared helpers are small
//! enough to reimplement here; if a third basis-indexed-by-something-else ever
//! arrives, *then* extract the common shape.
//!
//! Everything in this module is 1-based, per [`crate::permutation`].
//!
//! # What is implemented here
//!
//! The reference layer: the verified primitives and the E1 engine that exists
//! to be the in-house oracle — the role [`NaiveLr`](crate::lr::NaiveLr) plays for LR. E1 is
//! Symmetrica's own route: expand one factor to monomials, then push the other
//! through a Monk chain per monomial. Its cost is the number of **pipe
//! dreams**, which is why it walls at S₁₀–S₁₁ and why E3 exists. E1
//! is kept forever regardless of what wins.

// Every `as` here is a variable index or a permutation position, bounded by the
// permutation's length. Coefficients are the generic `C` and are never cast.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use crate::coeff::Ring;
use crate::partition::Partition;
use crate::permutation::Perm;
use crate::sym::Schur;
use std::collections::{BTreeMap, HashMap};

/// A monomial exponent vector: index `j-1` holds the exponent of `x_j`,
/// trailing zeros stripped so equal monomials compare equal.
pub type Expo = Vec<u32>;

/// A finite formal `C`-combination of permutations, i.e. an element of
/// `ℤ[x₁, x₂, …]` written in the Schubert basis.
///
/// `BTreeMap` for the same reason [`crate::sym`] uses one — deterministic
/// iteration — plus one this module actually depends on: [`Schubert::
/// from_polynomial`]'s greedy peel needs the lex-minimal monomial, and an
/// ordered map hands it over for free.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Schubert<C: Ring> {
    terms: BTreeMap<Perm, C>,
}

impl<C: Ring> Default for Schubert<C> {
    fn default() -> Self {
        Schubert::zero()
    }
}

fn strip(mut e: Expo) -> Expo {
    while e.last() == Some(&0) {
        e.pop();
    }
    e
}

impl<C: Ring> Schubert<C> {
    /// The zero element.
    pub fn zero() -> Self {
        Schubert {
            terms: BTreeMap::new(),
        }
    }

    /// `S_w` scaled by `c` (dropped if `c` is zero).
    pub fn monomial(w: Perm, c: C) -> Self {
        let mut t = BTreeMap::new();
        if !c.is_zero() {
            t.insert(w, c);
        }
        Schubert { terms: t }
    }

    /// `1 = S_id`.
    pub fn one() -> Self {
        Schubert::monomial(Perm::identity(), C::one())
    }

    /// The terms, keyed by permutation. Explicit zeros are never stored, so
    /// the map is empty exactly when the element is.
    pub fn terms(&self) -> &BTreeMap<Perm, C> {
        &self.terms
    }

    /// Whether this is 0.
    pub fn is_zero(&self) -> bool {
        self.terms.is_empty()
    }

    /// The coefficient of `S_w`, zero if `w` does not occur.
    pub fn coeff(&self, w: &Perm) -> C {
        self.terms.get(w).cloned().unwrap_or_else(C::zero)
    }

    /// `self += c·S_w`, removing the entry if it cancels to zero.
    pub fn add_term(&mut self, w: Perm, c: &C) {
        if c.is_zero() {
            return;
        }
        match self.terms.get_mut(&w) {
            Some(slot) => {
                slot.add_assign(c);
                if slot.is_zero() {
                    self.terms.remove(&w);
                }
            }
            None => {
                self.terms.insert(w, c.clone());
            }
        }
    }

    /// The sum, with terms that cancel dropped rather than stored as zero.
    pub fn add(&self, other: &Self) -> Self {
        let mut out = self.clone();
        out.add_assign(other);
        out
    }

    /// `self += other`, in place.
    ///
    /// The in-place form exists because E3's BRANCH step accumulates over
    /// covers, and `acc = acc.add(&part)` copies the whole accumulator once
    /// per cover — which a profile of `stair7²` showed as the dominant cost
    /// once the key allocations were gone.
    pub fn add_assign(&mut self, other: &Self) {
        for (w, c) in &other.terms {
            self.add_term(*w, c);
        }
    }

    /// The difference, with terms that cancel dropped rather than stored as
    /// zero.
    pub fn sub(&self, other: &Self) -> Self {
        let mut out = self.clone();
        for (w, c) in &other.terms {
            out.add_term(*w, &c.neg());
        }
        out
    }

    /// `self` scaled by `c`. A zero `c` gives the zero element rather than an
    /// element carrying zero coefficients.
    pub fn scale(&self, c: &C) -> Self {
        let mut out = Schubert::zero();
        for (w, a) in &self.terms {
            out.add_term(*w, &a.mul(c));
        }
        out
    }

    /// `S_w = s_λ(x₁, …, x_k)` for the Grassmannian `w` of descent `k` and
    /// shape `λ`.
    ///
    /// The convention — `code(w) = (λ_k, …, λ_1)`, i.e. `w(i) = λ_{k+1−i} + i`
    /// for `i ≤ k` with the unused values following in increasing order — was
    /// pinned against Sage over 128 shapes on 2026-07-30, together with the
    /// flag-truncation identity that makes leaf dispatch legal.
    ///
    /// `λ` must have at most `k` parts; otherwise `s_λ(x₁..x_k) = 0` and there
    /// is no such permutation, so this returns `None` rather than inventing
    /// one.
    pub fn grassmannian(lambda: &Partition, k: u32) -> Option<Self> {
        if lambda.len() as u32 > k {
            return None;
        }
        Some(Schubert::monomial(grassmannian_perm(lambda, k)?, C::one()))
    }

    /// `xᵢ · f`, the **signed Monk rule** (1-based), verified against Sage
    /// exhaustively on `S₄ × i ∈ 1..4`:
    ///
    /// ```text
    ///     xᵢ · S_w = Σ_{j>i, cover} S_{w t_{ij}} − Σ_{j<i, cover} S_{w t_{ji}}
    /// ```
    ///
    /// Both sums run over Bruhat covers only. The left sum is why
    /// intermediates are signed even though structure constants are not, and
    /// therefore why `guard`'s checked arithmetic is doing real work here.
    ///
    /// # Panics
    ///
    /// Panics if `i == 0` — variable indices are 1-based — or if a cover scan
    /// reaches past [`MAX_SUPPORT`](crate::permutation::MAX_SUPPORT).
    pub fn mul_variable(&self, i: u32) -> Self {
        assert!(i >= 1, "variable indices are 1-based");
        let mut out = Schubert::zero();
        for (w, c) in &self.terms {
            for (_, u) in w.covers_right(i) {
                out.add_term(u, c);
            }
            let neg = c.neg();
            for (_, u) in w.covers_left(i) {
                out.add_term(u, &neg);
            }
        }
        out
    }

    /// `∂ᵢ f` on the Schubert basis: `∂ᵢ S_w = S_{w sᵢ}` when `w(i) > w(i+1)`,
    /// and `0` otherwise. 1-based.
    ///
    /// A second, entirely independent implementation exists on *expanded
    /// polynomials* — `(f − sᵢf)/(xᵢ − xᵢ₊₁)` by exact division — in the tests,
    /// so the two can disagree loudly. They share no machinery, which is the
    /// point.
    ///
    /// # Panics
    ///
    /// Panics if `i == 0` — variable indices are 1-based — or if `i + 1`
    /// reaches past [`MAX_SUPPORT`](crate::permutation::MAX_SUPPORT).
    pub fn divided_difference(&self, i: u32) -> Self {
        assert!(i >= 1, "variable indices are 1-based");
        let mut out = Schubert::zero();
        for (w, c) in &self.terms {
            if w.at(i) > w.at(i + 1) {
                out.add_term(w.transpose(i, i + 1), c);
            }
        }
        out
    }

    /// `∂_w`, composing one-letter steps along a reduced word of `w`, applied
    /// left to right.
    ///
    /// With that convention the basis action is `S_v ↦ S_{v·w}` when
    /// `ℓ(v·w) = ℓ(v) − ℓ(w)`, and `0` otherwise. Independent of *which*
    /// reduced word, since `∂` satisfies the braid relations — asserted in the
    /// tests rather than assumed.
    pub fn divided_difference_perm(&self, w: &Perm) -> Self {
        let mut out = self.clone();
        for i in w.reduced_word() {
            out = out.divided_difference(i);
            if out.is_zero() {
                break;
            }
        }
        out
    }

    /// Expand into monomials of `ℤ[x₁, x₂, …]`.
    ///
    /// This is Symmetrica's `algorithmus2` peel (sb.c:323) with the one thing
    /// it lacks: memoization. The recursion tree's leaf count is
    /// `S_w(1,…,1)` — the pipe-dream count, which is super-exponential — while
    /// its distinct states are far fewer (measured 2026-07-30: 2 097 152
    /// leaves against 3 052 states on `stair7`). Here the memo is keyed on
    /// `(perm, stufe)`, the granularity that needs **no** bookkeeping to be
    /// sound.
    ///
    /// The sharper key is `perm` alone — `stufe` is read only at the DESCEND
    /// branch, so it moves the subtree value by a power of a single variable —
    /// which buys a further 1.3–1.6×. That refinement belongs with E3, not
    /// with the reference expansion, because it is exactly the kind of shift
    /// bookkeeping that is wrong in a way tests notice late.
    pub fn expand(&self) -> Vec<(Expo, C)> {
        let mut acc: BTreeMap<Expo, C> = BTreeMap::new();
        let mut memo = PeelMemo::default();
        for (w, c) in &self.terms {
            for (e, k) in memo.peel_perm(w) {
                let term = c.mul(&C::from_u128(k));
                match acc.get_mut(&e) {
                    Some(slot) => {
                        slot.add_assign(&term);
                        if slot.is_zero() {
                            acc.remove(&e);
                        }
                    }
                    None => {
                        if !term.is_zero() {
                            acc.insert(e, term);
                        }
                    }
                }
            }
        }
        acc.into_iter().collect()
    }

    /// Write a polynomial in the Schubert basis: the greedy peel of
    /// `t_POLYNOM_SCHUBERT` (sb.c:232), minus its recompute-everything defect.
    ///
    /// Correct because of **lex-triangularity**: `x^{code(w)}` is the
    /// lex-*minimal* monomial of `S_w` (verified against Sage on 40 random
    /// `S₇`). So the lex-minimal monomial of the remaining polynomial names
    /// its own permutation, and subtracting `c·S_w` strictly shrinks the
    /// problem. Symmetrica relies on this silently; here it is the documented
    /// reason the loop terminates.
    pub fn from_polynomial(terms: &[(Expo, C)]) -> Self {
        let mut rem: BTreeMap<Expo, C> = BTreeMap::new();
        for (e, c) in terms {
            if c.is_zero() {
                continue;
            }
            let e = strip(e.clone());
            match rem.get_mut(&e) {
                Some(slot) => {
                    slot.add_assign(c);
                    if slot.is_zero() {
                        rem.remove(&e);
                    }
                }
                None => {
                    rem.insert(e, c.clone());
                }
            }
        }

        let mut out = Schubert::zero();
        let mut memo = PeelMemo::default();
        while let Some((alpha, c)) = rem.iter().next().map(|(a, c)| (a.clone(), c.clone())) {
            let w = Perm::from_code(&alpha);
            out.add_term(w, &c);
            for (e, k) in memo.peel_perm(&w) {
                let sub = c.mul(&C::from_u128(k)).neg();
                match rem.get_mut(&e) {
                    Some(slot) => {
                        slot.add_assign(&sub);
                        if slot.is_zero() {
                            rem.remove(&e);
                        }
                    }
                    None => {
                        if !sub.is_zero() {
                            rem.insert(e, sub);
                        }
                    }
                }
            }
            debug_assert!(
                !rem.contains_key(&alpha),
                "peel failed to clear its own lex-minimal monomial — triangularity violated"
            );
        }
        out
    }

    /// **E1**, the reference product and permanent oracle: expand the factor
    /// with the smaller pipe-dream count, then run a Monk chain per monomial.
    ///
    /// This is Symmetrica's route, with its factor-choice heuristic fixed —
    /// `mult_schubert_schubert` compares only the *head term's* vector length,
    /// which is not even a proxy for the cost it is trying to avoid. Here the
    /// choice uses [`dimension`], the actual leaf count.
    ///
    /// Cost is `#pipe dreams × Monk chains`, the two compounding blowups
    /// behind the S₁₀ wall. That is intended: this engine's job is to be
    /// obviously correct.
    pub fn mul_naive(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Schubert::zero();
        }
        let (expand_side, keep_side) = if total_dimension(self) <= total_dimension(other) {
            (self, other)
        } else {
            (other, self)
        };
        let mut out = Schubert::zero();
        for (alpha, c) in expand_side.expand() {
            let mut t = keep_side.clone();
            for (j, &e) in alpha.iter().enumerate() {
                for _ in 0..e {
                    t = t.mul_variable((j + 1) as u32);
                }
            }
            for (w, a) in &t.terms {
                out.add_term(*w, &a.mul(&c));
            }
        }
        out
    }

    /// The product `S_u · S_v`, by whichever engine is currently the best one.
    ///
    /// **That is [`mul_e2`](Self::mul_e2), the memoized transition** — measured,
    /// on `docs/record/schubert.md`'s ladder, at 1.9–27.6× the C `schubmult`
    /// where that finishes and completing rows it does not. The engine this
    /// used to call, [`mul_e3`](Self::mul_e3), is 4–97× *slower* and is kept as
    /// the historical comparison rather than as a route anyone should take.
    ///
    /// Callers wanting a specific engine name it: [`mul_naive`](Self::mul_naive)
    /// is E1, the permanent oracle; `mul_e2` and `mul_e3` are the two engines.
    /// This function is free to change which one it delegates to, and has.
    pub fn mul(&self, other: &Self) -> Self {
        self.mul_e2(other)
    }

    /// **E3**: evaluate the merged peel DAG of one factor once, carrying the
    /// partial product at every state.
    ///
    /// E1 walks a *tree* whose leaves are pipe dreams and re-runs a Monk chain
    /// per leaf. E3 walks the same recursion as a **DAG**, merging states, and
    /// stores `V(state)·S_u` at each. The transitions become operators:
    /// DESCEND is `stufe` applications of [`Schubert::mul_variable`], BRANCH
    /// is an addition, and the leaf is `S_u` or `x_level·S_u`.
    ///
    /// **Superseded by [`mul_e2`](Self::mul_e2), and the reason is the whole
    /// lesson of this module.** A pre-implementation measurement of *state
    /// compression* — `pipe dreams ÷ states`, 687× on `stair7` and 125 599× on
    /// an S₁₃ element — pointed hard at this engine, so it was built first. It
    /// duly beat E1 by 11.2× on `stair6²`, and then lost to the C `schubmult`
    /// by 4–51×, and then lost to E2 by up to 97×. Both engines cost (nodes) ×
    /// (size of the running element); **the metric counted only nodes**. On
    /// `S_11.1` E2 uses *more* nodes than E3 and is 97× faster, because E1 and
    /// E3 expand a factor into monomials so the running element inflates to
    /// answer-size early, while the transition recursion never expands.
    ///
    /// A cost model that omits a factor will rank engines confidently and
    /// wrongly. Kept, and kept tested, so the comparison stays reproducible.
    ///
    /// Keyed on `(perm, level, stufe)` — the granularity that is sound with no
    /// bookkeeping. Merging on `perm` alone is worth a further 1.3–1.6× and was
    /// never done: it needs a shift by a power of `x_level`, hence more Monk
    /// passes, and E2 overtook the engine before the trade was worth measuring.
    pub fn mul_e3(&self, other: &Self) -> Self {
        if self.is_zero() || other.is_zero() {
            return Schubert::zero();
        }
        // Choose by *state* count, not pipe dreams: states are what E3 pays.
        // (E1's choice by `dimension` is right for E1 and wrong here — the two
        // orders differ, since compression varies by permutation shape.)
        let (peel_side, keep_side) = if total_states(self) <= total_states(other) {
            (self, other)
        } else {
            (other, self)
        };
        let mut out = Schubert::zero();
        for (w, c) in peel_side.terms() {
            let part = E3 {
                keep: keep_side,
                memo: HashMap::new(),
            }
            .eval_perm(w);
            for (v, a) in &part.terms {
                out.add_term(*v, &a.mul(c));
            }
        }
        out
    }

    /// **E2**, the memoized transition recursion. Recurse on the *second*
    /// factor with Lascoux–Schützenberger:
    ///
    /// ```text
    ///     S_u · S_v = x_r · (S_u · S_v')  +  Σ (S_u · S_v'')
    /// ```
    ///
    /// where `(r, v', {v''}) = v.transition()`, `ℓ(v') = ℓ(v) − 1`, and each
    /// `ℓ(v'') = ℓ(v)`. Base case: `v` dominant, where `S_v` is the single
    /// monomial `x^{code(v)}` and the product is `ℓ(v)` Monk passes.
    ///
    /// **Why this exists even though E3 already works.** E1 and E3 both expand
    /// one factor into monomials, so both pay (number of nodes) × (size of the
    /// running element). On a large-output case that second factor is the size
    /// of the answer, which is where E3 lands 51× behind the C `schubmult`
    /// while being only ~4× behind on staircases. E2 has the same
    /// *shape* of cost, so it is not automatically better; what differs is the
    /// node count, and whether the transition tree is smaller than the peel
    /// DAG is a measurement, not an argument. `examples/bench_e2.rs` makes it.
    ///
    /// ⚠️ Termination is not obvious and is not proved here: the `v''` have the
    /// *same* length as `v`, so the recursion does not descend on `ℓ`. It is
    /// the recursion `newtrans` (mss.c:44) has relied on for thirty years, and
    /// Macdonald's notes prove the tree finite. [`Schubert::mul_e2_depth`]
    /// carries an explicit depth cap so a counterexample is a clean panic
    /// rather than a hang.
    pub fn mul_e2(&self, other: &Self) -> Self {
        self.mul_e2_depth(other, 4096).0
    }

    /// [`Schubert::mul_e2`] with an explicit recursion cap, also returning the
    /// two cost counters — nodes visited and Monk passes performed — because
    /// those are what the E2-vs-E3 comparison actually turns on.
    pub fn mul_e2_depth(&self, other: &Self, cap: u32) -> (Self, u64, u64) {
        if self.is_zero() || other.is_zero() {
            return (Schubert::zero(), 0, 0);
        }
        let mut out = Schubert::zero();
        let (mut nodes, mut passes) = (0u64, 0u64);
        for (w, c) in other.terms() {
            let mut refs = HashMap::new();
            count_transition_parents(w, &mut refs);
            let mut e2 = E2 {
                keep: self,
                target: None,
                memo: HashMap::new(),
                refs,
                cap,
                passes: 0,
                nodes: 0,
                peak_live: 0,
            };
            let part = e2.eval(w, 0);
            nodes += e2.nodes;
            passes += e2.passes;
            for (v, a) in &part.terms {
                out.add_term(*v, &a.mul(c));
            }
        }
        (out, nodes, passes)
    }

    /// The Poincaré pairing on `H*(Fl(n))`: `⟨S_u, S_v⟩ = δ_{v, w₀u}`, i.e.
    /// the coefficient of `S_{w₀⁽ⁿ⁾}` in the product.
    ///
    /// **`n` is an explicit argument**, unlike the incumbent.
    /// `scalarproduct_schubert` (sb.c:1840) reads `n` off however long the
    /// stored vectors happen to be, so the same mathematical inputs give
    /// different answers depending on prior padding — semantics a caller
    /// cannot predict. It also computes the full product and then runs
    /// `n(n−1)/2` divided-difference passes; reading one coefficient is the
    /// same number.
    ///
    /// # Panics
    ///
    /// Panics if `n` exceeds [`MAX_SUPPORT`](crate::permutation::MAX_SUPPORT),
    /// which is the only way `w₀⁽ⁿ⁾` fails to be representable — the reversal
    /// of `1..=n` is a permutation for every other `n`, including `n = 0`.
    pub fn pairing(&self, other: &Self, n: u32) -> C {
        let w0 = Perm::new((1..=n).rev())
            .expect("the reversal of 1..=n is a permutation unless n > MAX_SUPPORT");
        self.mul(other).coeff(&w0)
    }
}

/// The Grassmannian permutation of descent `k` and shape `λ` (at most `k`
/// parts): `w(i) = λ_{k+1−i} + i` for `i ≤ k`, unused values in order after.
pub fn grassmannian_perm(lambda: &Partition, k: u32) -> Option<Perm> {
    if lambda.len() as u32 > k {
        return None;
    }
    let k = k as usize;
    let mut head: Vec<u32> = Vec::with_capacity(k);
    for i in 1..=k {
        head.push(lambda.part(k - i) + i as u32);
    }
    let n = head.iter().copied().max().unwrap_or(0);
    let used: std::collections::BTreeSet<u32> = head.iter().copied().collect();
    let mut v = head;
    v.extend((1..=n).filter(|x| !used.contains(x)));
    Perm::new(v).ok()
}

/// `1 × w`: the shift `w'(1) = 1`, `w'(i+1) = w(i) + 1`.
///
/// Stanley symmetric functions are invariant under it, which is what makes it
/// a legal move rather than a different question.
fn shift_up(w: &Perm) -> Perm {
    let m = w.support_len();
    let mut v = Vec::with_capacity(m as usize + 1);
    v.push(1);
    for i in 1..=m {
        v.push(w.at(i) + 1);
    }
    Perm::new(v).expect("the shift 1 x w is a permutation unless it exceeds MAX_SUPPORT")
}

/// The **Stanley symmetric function** `F_w`, expanded in Schur functions.
///
/// The transition recursion with the `x_r·S_v` term dropped:
/// `F_w = Σ F_{v t_{qr}}` over the same left covers, with the Grassmannian
/// base case `F_w = s_λ` and the `F_w = F_{1×w}` shift when the left sum is
/// empty. Iterative with an explicit stack, per the incumbent, minus its
/// static `char[1000]`.
///
/// ⚠️ **This is not "Schubert → Schur".** Symmetrica exports it as
/// `t_SCHUBERT_SCHUR` and Sage inherits that name, but `F_w = S_w` only in the
/// stable range: `newtrans([2,1,4,3]) = s₂ + s₁₁` while `S_{2143}` is not even
/// symmetric. The name here says what it computes.
///
/// This is the module's one bridge into `Sym`, and the only reason `Schur`
/// appears in this file.
///
/// # Panics
///
/// Panics if the transition tree takes more than `2²⁴` steps. The cap is a
/// backstop against a non-terminating recursion, which would otherwise hang
/// rather than fail; it is **not** a measured capacity wall, and where in `S_n`
/// a legitimate expansion first reaches it is unmeasured.
///
/// Also if the shift `1 × w` would exceed
/// [`MAX_SUPPORT`](crate::permutation::MAX_SUPPORT).
pub fn stanley<C: Ring>(w: &Perm) -> Schur<C> {
    use crate::sym::SymFn;
    let mut out: Schur<C> = Schur::zero();
    let mut stack = vec![*w];
    let mut steps = 0u64;
    while let Some(p) = stack.pop() {
        steps += 1;
        assert!(
            steps < 1 << 24,
            "stanley: transition tree failed to terminate"
        );
        // last_descent() is the same test as is_grassmannian() for this
        // purpose and allocates nothing: no descent, or one, means vexillary
        // with a single Schur term.
        let Some(r) = p.last_descent() else {
            // identity: F = s_empty
            let lambda = Partition::default();
            let c = out.coeff(&lambda);
            let mut c = c;
            c.add_assign(&C::one());
            out.terms_mut().insert(lambda, c);
            continue;
        };
        if p.is_grassmannian().is_some() {
            // Grassmannian is vexillary with shape = its code sorted; for a
            // single descent the code is already weakly increasing.
            let lambda = Partition::new(p.code());
            let c = out.coeff(&lambda);
            let mut c = c;
            c.add_assign(&C::one());
            out.terms_mut().insert(lambda, c);
            continue;
        }
        // transition, inlined so nothing is allocated: the ups go straight
        // onto the work stack instead of through two intermediate `Vec`s.
        let wr = p.at(r);
        let s = (r + 1..=p.support_len())
            .rfind(|&s| p.at(s) < wr)
            .expect("a descent has something smaller to its right");
        let v = p.transpose(r, s);
        let before = stack.len();
        v.for_each_cover_left(r, |_, u| stack.push(u));
        if stack.len() == before {
            // no room to the left of r; shifting makes room and F is invariant
            stack.push(shift_up(&p));
        }
    }
    out
}

/// `S_u(1,…,1)·S_v(1,…,1)`: the product's total monomial mass, in
/// microseconds. The one cheap quantity that flags an out-of-family pair, and
/// deliberately not enforced anywhere — whether to refuse such a pair or
/// attempt it and die is still open (`docs/record/schubert.md`).
///
/// # Range
///
/// **Saturates** at `u128::MAX` rather than overflowing, which is why R4
/// (`docs/policies/failure.md`) permits it here: this is a magnitude used to
/// decide "too big to attempt", and a pair whose true mass exceeds `u128` is
/// one the saturated value classifies identically. A caller reading it as the
/// exact mass rather than as a cost signal is reading it wrong — at that
/// magnitude the product is unattemptable regardless.
pub fn schubert_monomial_mass_of(u: &Perm, v: &Perm) -> u128 {
    dimension(u).saturating_mul(dimension(v))
}

/// A single structure constant `c^w_{uv}`, without building the whole product.
///
/// The motivating case is the one no engine can complete: `S_13 ℓ=25,36` has a
/// monomial mass of 4.3×10¹⁶, so its full expansion cannot be materialized on
/// any machine. One coefficient of it is still a perfectly reasonable
/// question, and asking it is what positivity searches and rule-hunting
/// actually do.
///
/// Runs E2 with every term not `≤ w` in Bruhat order discarded at each node.
/// The pruning is sound because signed Monk moves strictly up the Bruhat order
/// (`monk_covers_move_up_the_bruhat_order`), so nothing below the cut can
/// reach `w`. It is also *cheap* to check — the tableau criterion on an inline
/// `Perm` — which matters because it runs on every surviving term.
///
/// Returns zero immediately unless `ℓ(w) = ℓ(u) + ℓ(v)` and `u ≤ w`, `v ≤ w`;
/// those are necessary conditions
/// (`product_support_lies_above_both_factors`), and they alone answer most
/// queries with no work at all.
pub fn schubert_coeff<C: Ring>(u: &Perm, v: &Perm, w: &Perm) -> C {
    if w.length() != u.length() + v.length() || !u.bruhat_le(w) || !v.bruhat_le(w) {
        return C::zero();
    }
    // c^w_{uv} = c^w_{vu}, so recurse on whichever side has the smaller
    // transition tree. `transition_tree` builds no elements, so asking is
    // microseconds. This is the `lr_coeff` lesson applied before the fact
    // rather than after: peeling the wrong side cost 21x there.
    let (recurse, keep_perm) = if transition_tree(v).0 <= transition_tree(u).0 {
        (v, u)
    } else {
        (u, v)
    };
    let keep: Schubert<C> = Schubert::monomial(*keep_perm, C::one());
    let mut refs = HashMap::new();
    count_transition_parents(recurse, &mut refs);
    let mut e2 = E2 {
        keep: &keep,
        target: Some(*w),
        memo: HashMap::new(),
        refs,
        cap: 4096,
        passes: 0,
        nodes: 0,
        peak_live: 0,
    };
    let out = e2.eval(recurse, 0);
    out.coeff(w)
}

/// `S_w(1,…,1)` — the number of pipe dreams of `w`, and exactly the leaf count
/// of the peel recursion.
///
/// Not re-exported at the crate root: [`crate::eval::dimension`] is a
/// different function on partitions, and two `dimension`s in one namespace is
/// how a caller gets a plausible wrong answer.
pub fn dimension(w: &Perm) -> u128 {
    PeelMemo::default().count(w)
}

/// Number of distinct `(perm, level, stufe)` states in the peel DAG of `w` —
/// E3's cost, the way [`dimension`] is E1's.
///
/// Cheap: it walks the DAG the product would walk anyway, without building any
/// element at the states.
pub fn peel_states(w: &Perm) -> u64 {
    let m = w.support_len();
    if m == 0 {
        return 1;
    }
    let n = m;
    let mut seen = std::collections::HashSet::new();
    let mut stack = vec![(w.padded(m), m.saturating_sub(1))];
    while let Some((p, stufe)) = stack.pop() {
        let len = p.len() as u32;
        if !seen.insert((p.clone(), n - len + 1, stufe)) {
            continue;
        }
        if len <= 2 {
            continue;
        }
        if p[0] == len {
            stack.push((p[1..].to_vec(), len.saturating_sub(2)));
        } else {
            for q in covers_at_1(&p) {
                stack.push((q, stufe.saturating_sub(1)));
            }
        }
    }
    seen.len() as u64
}

/// Saturating by R4's estimate clause: the sum is only ever compared against
/// another `total_states` to pick which side to expand, so saturation changes
/// the answer only when both sides saturate — where the two are equally
/// hopeless and either choice is as good. An exact width here would buy a
/// panic on a quantity that is not part of any result.
fn total_states<C: Ring>(f: &Schubert<C>) -> u64 {
    f.terms()
        .keys()
        .map(peel_states)
        .fold(0u64, |a, b| a.saturating_add(b))
}

/// E2's evaluator: the transition tree with `S_u·S_v` cached at each node
/// **only until its last parent has consumed it**.
///
/// The memo has to be evicted, not merely populated. A node's cached value is
/// a whole product — on `S_13.2` the answer alone is 3.2M terms — so retaining
/// every node for the length of the run is gigabytes. Before this existed,
/// `S_13 ℓ=25,36` died out of memory (`docs/record/schubert.md`): the failure
/// predicted for E3 in as many words, which E2 then shipped anyway.
///
/// So a pre-pass counts how many parents each node has, and `eval` drops a
/// child's value the moment the last of them is done with it. The pre-pass is
/// cheap in the way that matters: it walks the same tree using only `Perm`
/// operations, computing no products.
struct E2<'a, C: Ring> {
    keep: &'a Schubert<C>,
    /// When set, drop every term not `≤ target` in Bruhat order. Sound
    /// because Monk only moves up (both facts are tested), so a dropped term
    /// can reach the target neither directly nor by cancelling against
    /// something that does.
    target: Option<Perm>,
    memo: HashMap<Perm, std::rc::Rc<Schubert<C>>>,
    /// Parents not yet served, per node. At zero the memo entry is dropped.
    refs: HashMap<Perm, u32>,
    cap: u32,
    passes: u64,
    /// Distinct nodes evaluated. Not `memo.len()` any more — that is now the
    /// *live* set, which is the point.
    nodes: u64,
    /// High-water mark of the live set.
    peak_live: usize,
}

/// Size of E2's transition tree for `w`: `(distinct nodes, total edges)`.
///
/// **Computes no products.** It walks exactly the tree `mul_e2` walks using
/// only `Perm` operations, so it costs microseconds on inputs whose actual
/// product does not finish — which makes it the probe for asking *why* a case
/// is hard before paying to find out. `edges − nodes + 1` is the amount of
/// sharing the memo captures.
pub fn transition_tree(w: &Perm) -> (u64, u64) {
    let mut refs = HashMap::new();
    count_transition_parents(w, &mut refs);
    let nodes = refs.len() as u64;
    let edges: u64 = refs.values().map(|&n| n as u64).sum();
    (nodes, edges)
}

/// Count parent references over the transition tree, expanding each distinct
/// node once. The root gets one reference, from the caller.
fn count_transition_parents(root: &Perm, refs: &mut HashMap<Perm, u32>) {
    let mut stack = vec![*root];
    let mut expanded = std::collections::HashSet::new();
    while let Some(p) = stack.pop() {
        *refs.entry(p).or_insert(0) += 1;
        if !expanded.insert(p) || p.is_dominant() {
            continue;
        }
        let (_, vp, ups) = p.transition().expect("non-dominant implies a descent");
        stack.push(vp);
        stack.extend(ups);
    }
}

impl<C: Ring> E2<'_, C> {
    /// Note the child's value has been consumed; drop it if no parent remains.
    fn release(&mut self, child: &Perm) {
        if let Some(n) = self.refs.get_mut(child) {
            *n = n.saturating_sub(1);
            if *n == 0 {
                self.memo.remove(child);
            }
        }
    }

    fn eval(&mut self, v: &Perm, depth: u32) -> std::rc::Rc<Schubert<C>> {
        if let Some(r) = self.memo.get(v) {
            return std::rc::Rc::clone(r);
        }
        assert!(
            depth < self.cap,
            "E2 recursion exceeded depth {} at {v} — the transition tree was \
             supposed to be finite",
            self.cap
        );
        let out = if v.is_dominant() {
            // S_v is the single monomial x^code(v); the identity lands here
            // too, with an empty code and no work.
            let mut acc = self.keep.clone();
            for (i, &e) in v.code().iter().enumerate() {
                for _ in 0..e {
                    acc = acc.mul_variable(i as u32 + 1);
                    self.passes += 1;
                }
            }
            acc
        } else {
            let (r, vp, ups) = v.transition().expect("non-dominant implies a descent");
            let child = self.eval(&vp, depth + 1);
            let mut acc = child.mul_variable(r);
            self.passes += 1;
            drop(child);
            self.release(&vp);
            for u in &ups {
                let part = self.eval(u, depth + 1);
                acc.add_assign(&part);
                drop(part);
                self.release(u);
            }
            acc
        };
        self.nodes += 1;
        let mut out = out;
        if let Some(t) = &self.target {
            out.terms.retain(|z, _| z.bruhat_le(t));
        }
        let out = std::rc::Rc::new(out);
        self.memo.insert(*v, std::rc::Rc::clone(&out));
        self.peak_live = self.peak_live.max(self.memo.len());
        out
    }
}

/// E3's evaluator: the peel DAG with `V(state)·S_u` cached at every state.
struct E3<'a, C: Ring> {
    keep: &'a Schubert<C>,
    /// `Rc`, not the element itself: a state's value is read once per parent
    /// and the DAG is what makes parents plural, so a by-value memo copies
    /// whole elements for no reason.
    memo: HashMap<(Vec<u32>, u32, u32), std::rc::Rc<Schubert<C>>>,
}

impl<C: Ring> E3<'_, C> {
    fn eval_perm(mut self, w: &Perm) -> Schubert<C> {
        let m = w.support_len();
        if m == 0 {
            return self.keep.clone();
        }
        let v = w.padded(m);
        let r = self.eval(&v, m, m.saturating_sub(1));
        std::rc::Rc::try_unwrap(r).unwrap_or_else(|rc| (*rc).clone())
    }

    fn eval(&mut self, p: &[u32], n: u32, stufe: u32) -> std::rc::Rc<Schubert<C>> {
        let m = p.len() as u32;
        let level = n - m + 1;
        let key = (p.to_vec(), level, stufe);
        if let Some(v) = self.memo.get(&key) {
            return std::rc::Rc::clone(v);
        }
        let out = if m <= 2 {
            // LEAF: V = x_level for the transposition, 1 for the identity.
            if m == 2 && p[0] == 2 {
                self.keep.mul_variable(level)
            } else {
                self.keep.clone()
            }
        } else if p[0] == m {
            // DESCEND: V = x_level^stufe · (the subtree one alphabet along)
            let child = self.eval(&p[1..], n, m.saturating_sub(2));
            if stufe == 0 {
                (*child).clone()
            } else {
                let mut acc = child.mul_variable(level);
                for _ in 1..stufe {
                    acc = acc.mul_variable(level);
                }
                acc
            }
        } else {
            // BRANCH: sum over the covers raising p[0]
            let mut acc = Schubert::zero();
            for q in covers_at_1(p) {
                let part = self.eval(&q, n, stufe.saturating_sub(1));
                acc.add_assign(&part);
            }
            acc
        };
        let out = std::rc::Rc::new(out);
        self.memo.insert(key, std::rc::Rc::clone(&out));
        out
    }
}

/// Saturating for the same reason as [`total_states`]: a dispatch magnitude
/// compared only against its own counterpart, never returned as a result.
fn total_dimension<C: Ring>(f: &Schubert<C>) -> u128 {
    f.terms()
        .keys()
        .map(dimension)
        .fold(0u128, |a, b| a.saturating_add(b))
}

// ---------------------------------------------------------------------------
// the peel recursion (Symmetrica `algorithmus2`, sb.c:323)
// ---------------------------------------------------------------------------

/// Memo for the peel, keyed on `(one-line prefix, level, stufe)`.
///
/// The C carries `(vorfaktor, alphabetindex, stufe, perm)` and merges nothing.
/// One of those is genuinely derivable — the accumulated prefix `vorfaktor`
/// belongs outside the subtree, so the subtree value is a function of the
/// rest.
///
/// **`level` must stay in the key even though `len(perm)` pins it.** It pins
/// it only *within one top-level call*: `level = n − len(p) + 1`, and `n`
/// differs between permutations. A memo shared across `S_{132}` (n = 3) and
/// `S_{1423}` (n = 4) would hand a length-3 state computed in `x₁, x₂` back to
/// a caller expecting `x₂, x₃`. That is not hypothetical — it is what this key
/// looked like an hour ago, and the single-permutation round trip passed
/// exhaustively through S₆ while only [`Schubert::from_polynomial`], which
/// shares one memo across permutations of different sizes, saw it.
#[derive(Default)]
struct PeelMemo {
    poly: HashMap<(Vec<u32>, u32, u32), Vec<(Expo, u128)>>,
    /// Leaf counts need no `level` or `stufe`: neither can change how many
    /// leaves a subtree has, only which monomials they carry. This is the
    /// `states_p` granularity of the compression table.
    count: HashMap<Vec<u32>, u128>,
}

impl PeelMemo {
    fn peel_perm(&mut self, w: &Perm) -> Vec<(Expo, u128)> {
        let m = w.support_len();
        if m == 0 {
            return vec![(Vec::new(), 1)];
        }
        let v = w.padded(m);
        let n = v.len() as u32;
        self.peel(&v, n, m.saturating_sub(1))
    }

    /// `V(p, level, stufe)` in absolute variable indices, where
    /// `level = n − len(p) + 1`.
    fn peel(&mut self, p: &[u32], n: u32, stufe: u32) -> Vec<(Expo, u128)> {
        let m = p.len() as u32;
        let level = n - m + 1;
        let key = (p.to_vec(), level, stufe);
        if let Some(v) = self.poly.get(&key) {
            return v.clone();
        }
        let out: Vec<(Expo, u128)> = if m <= 2 {
            // LEAF: the identity contributes 1, the transposition x_level.
            let mut e = Vec::new();
            if m == 2 && p[0] == 2 {
                e = vec![0; level as usize];
                e[(level - 1) as usize] = 1;
            }
            vec![(strip(e), 1)]
        } else if p[0] == m {
            // DESCEND: x_level^stufe · (the rest, one alphabet along)
            let sub = self.peel(&p[1..], n, m.saturating_sub(2));
            let mut acc: BTreeMap<Expo, u128> = BTreeMap::new();
            for (e, k) in sub {
                let mut e2 = e;
                if stufe > 0 {
                    if e2.len() < level as usize {
                        e2.resize(level as usize, 0);
                    }
                    e2[(level - 1) as usize] += stufe;
                }
                *acc.entry(strip(e2)).or_insert(0) += k;
            }
            acc.into_iter().collect()
        } else {
            // BRANCH: the covers raising p[0], stufe one lower
            let mut acc: BTreeMap<Expo, u128> = BTreeMap::new();
            for q in covers_at_1(p) {
                for (e, k) in self.peel(&q, n, stufe.saturating_sub(1)) {
                    *acc.entry(e).or_insert(0) += k;
                }
            }
            acc.into_iter().collect()
        };
        self.poly.insert(key, out.clone());
        out
    }

    /// Leaf count only — the same DAG without building polynomials.
    fn count(&mut self, w: &Perm) -> u128 {
        let m = w.support_len();
        if m == 0 {
            return 1;
        }
        let v = w.padded(m);
        self.count_rec(&v)
    }

    fn count_rec(&mut self, p: &[u32]) -> u128 {
        if let Some(&v) = self.count.get(p) {
            return v;
        }
        let m = p.len() as u32;
        let out = if m <= 2 {
            1
        } else if p[0] == m {
            self.count_rec(&p[1..])
        } else {
            covers_at_1(p)
                .into_iter()
                .map(|q| self.count_rec(&q))
                .fold(0u128, |a, b| a.saturating_add(b))
        };
        self.count.insert(p.to_vec(), out);
        out
    }
}

/// The Bruhat covers `p ⋖ p·t_{1,i}` raising `p[0]`, in Symmetrica's
/// running-minimum scan order.
///
/// This is [`Perm::covers_right`] at `i = 1` specialized to a raw padded
/// vector: the peel works on fixed-length vectors of a known `n`, and routing
/// through `Perm` would strip fixed points and lose the level bookkeeping.
fn covers_at_1(p: &[u32]) -> Vec<Vec<u32>> {
    let mut out = Vec::new();
    let mut running_min = u32::MAX;
    for i in 1..p.len() {
        if p[i] > p[0] && p[i] < running_min {
            running_min = p[i];
            let mut q = p.to_vec();
            q.swap(0, i);
            out.push(q);
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(v: &[u32]) -> Perm {
        Perm::new(v.iter().copied()).unwrap()
    }

    fn sch(v: &[u32]) -> Schubert<i64> {
        Schubert::monomial(p(v), 1)
    }

    fn expo(e: &[u32]) -> Expo {
        strip(e.to_vec())
    }

    // --- layer 1: hand values ----------------------------------------------

    #[test]
    fn hand_values() {
        // S_132 = x1 + x2
        assert_eq!(
            sch(&[1, 3, 2]).expand(),
            vec![(expo(&[0, 1]), 1), (expo(&[1]), 1)]
        );
        // S_321 = x1^2 x2   (dominant, one monomial)
        assert_eq!(sch(&[3, 2, 1]).expand(), vec![(expo(&[2, 1]), 1)]);
        // S_1423 = x1^2 + x1x2 + x2^2   (Knutson p.6, the pipe-dream example)
        assert_eq!(
            sch(&[1, 4, 2, 3]).expand(),
            vec![(expo(&[0, 2]), 1), (expo(&[1, 1]), 1), (expo(&[2]), 1)]
        );
        // S_id = 1
        assert_eq!(Schubert::<i64>::one().expand(), vec![(expo(&[]), 1)]);
    }

    #[test]
    fn dominant_is_a_single_monomial() {
        for n in 0..=6u32 {
            for w in crate::permutation::tests::all_perms(n) {
                if w.is_dominant() {
                    let e = sch(&w.padded(n)).expand();
                    assert_eq!(e.len().max(1), 1, "{w} expanded to {e:?}");
                    assert_eq!(e[0].0, strip(w.code()), "{w}");
                    assert_eq!(dimension(&w), 1, "{w}");
                }
            }
        }
    }

    #[test]
    fn grassmannian_matches_schur_in_k_variables() {
        // S_w = s_lambda(x_1..x_k); check the two smallest by hand, plus the
        // shape/permutation convention.
        let w = grassmannian_perm(&Partition::new([1]), 2).unwrap();
        assert_eq!(w, p(&[1, 3, 2]));
        let w = grassmannian_perm(&Partition::new([2, 1]), 2).unwrap();
        assert_eq!(w, p(&[2, 4, 1, 3]));
        assert_eq!(w.is_grassmannian(), Some(2));
        // s_{2,1}(x1,x2) = x1^2 x2 + x1 x2^2
        assert_eq!(
            sch(&[2, 4, 1, 3]).expand(),
            vec![(expo(&[1, 2]), 1), (expo(&[2, 1]), 1)]
        );
        // more parts than variables has no Grassmannian permutation
        assert!(grassmannian_perm(&Partition::new([1, 1, 1]), 2).is_none());
    }

    // --- layer 2: structural laws ------------------------------------------

    #[test]
    fn degree_equals_length() {
        for n in 0..=6u32 {
            for w in crate::permutation::tests::all_perms(n) {
                for (e, _) in sch(&w.padded(n)).expand() {
                    assert_eq!(e.iter().sum::<u32>(), w.length(), "{w} term {e:?}");
                }
            }
        }
    }

    #[test]
    fn divided_difference_is_nilpotent_and_braided() {
        for n in 2..=5u32 {
            for w in crate::permutation::tests::all_perms(n) {
                let f = sch(&w.padded(n));
                for i in 1..n {
                    // d_i d_i = 0
                    assert!(
                        f.divided_difference(i).divided_difference(i).is_zero(),
                        "{w} {i}"
                    );
                    // commuting braid relation
                    for j in (i + 2)..n {
                        let a = f.divided_difference(i).divided_difference(j);
                        let b = f.divided_difference(j).divided_difference(i);
                        assert_eq!(a, b, "{w} commute {i},{j}");
                    }
                    // long braid relation
                    if i + 1 < n {
                        let j = i + 1;
                        let a = f
                            .divided_difference(i)
                            .divided_difference(j)
                            .divided_difference(i);
                        let b = f
                            .divided_difference(j)
                            .divided_difference(i)
                            .divided_difference(j);
                        assert_eq!(a, b, "{w} braid {i},{j}");
                    }
                }
            }
        }
    }

    /// `∂` computed on the basis, against `∂` computed by exact division on
    /// the expanded polynomial. The two share no machinery.
    #[test]
    fn divided_difference_agrees_with_exact_division() {
        for n in 2..=5u32 {
            for w in crate::permutation::tests::all_perms(n) {
                let f = sch(&w.padded(n));
                for i in 1..=(n + 1) {
                    let via_basis = f.divided_difference(i).expand();
                    let via_poly = poly_divided_difference(&f.expand(), i);
                    assert_eq!(via_basis, via_poly, "{w} at i={i}");
                }
            }
        }
    }

    /// `(f − s_i f)/(x_i − x_{i+1})`, by exact division, on exponent vectors.
    fn poly_divided_difference(terms: &[(Expo, i64)], i: u32) -> Vec<(Expo, i64)> {
        let i = (i - 1) as usize;
        let mut num: BTreeMap<Expo, i64> = BTreeMap::new();
        for (e, c) in terms {
            let mut a = e.clone();
            if a.len() <= i + 1 {
                a.resize(i + 2, 0);
            }
            *num.entry(strip(a.clone())).or_insert(0) += c;
            a.swap(i, i + 1);
            *num.entry(strip(a)).or_insert(0) -= c;
        }
        num.retain(|_, c| *c != 0);
        // divide by (x_i - x_{i+1}) with a monomial-order long division
        let mut quot: BTreeMap<Expo, i64> = BTreeMap::new();
        while let Some((lead, c)) = num.iter().next_back().map(|(e, c)| (e.clone(), *c)) {
            let mut q = lead.clone();
            if q.len() <= i {
                q.resize(i + 1, 0);
            }
            assert!(q[i] >= 1, "not divisible by (x_i - x_{{i+1}})");
            q[i] -= 1;
            let q = strip(q);
            *quot.entry(q.clone()).or_insert(0) += c;
            // subtract c·q·(x_i − x_{i+1})
            for (slot, sign) in [(i, 1i64), (i + 1, -1i64)] {
                let mut t = q.clone();
                if t.len() <= slot {
                    t.resize(slot + 1, 0);
                }
                t[slot] += 1;
                let t = strip(t);
                *num.entry(t.clone()).or_insert(0) -= c * sign;
                if num[&t] == 0 {
                    num.remove(&t);
                }
            }
        }
        quot.retain(|_, c| *c != 0);
        quot.into_iter().collect()
    }

    #[test]
    fn divided_difference_perm_lands_where_it_should() {
        for n in 2..=5u32 {
            for w in crate::permutation::tests::all_perms(n) {
                for v in crate::permutation::tests::all_perms(n) {
                    let got = sch(&v.padded(n)).divided_difference_perm(&w);
                    // S_v |-> S_{v w} when the length drops by exactly l(w)
                    let mut vw = v.padded(n);
                    let wp = w.padded(n);
                    let composed: Vec<u32> =
                        (0..n as usize).map(|i| vw[(wp[i] - 1) as usize]).collect();
                    vw = composed;
                    let vw = Perm::new(vw).unwrap();
                    let want = if vw.length() + w.length() == v.length() {
                        sch(&vw.padded(n))
                    } else {
                        Schubert::zero()
                    };
                    assert_eq!(got, want, "d_{w} on S_{v}");
                }
            }
        }
    }

    // --- layer 3: the defining recursion ------------------------------------

    /// `S_{w₀⁽ⁿ⁾} = x^δ`, and the downward `∂` recursion reproduces `expand()`
    /// for all of S₅.
    #[test]
    fn defining_recursion_reproduces_expand() {
        for n in 2..=5u32 {
            let w0: Vec<u32> = (1..=n).rev().collect();
            let top = sch(&w0);
            let delta: Expo = strip((1..=n).rev().map(|x| x - 1).collect::<Vec<u32>>());
            assert_eq!(top.expand(), vec![(delta, 1)], "w0 in S_{n}");

            for w in crate::permutation::tests::all_perms(n) {
                // w0 * w^{-1}-ish: reach S_w from S_{w0} by divided differences
                let u = compose(&Perm::new(w0.iter().copied()).unwrap(), &w, n);
                let got = top.divided_difference_perm(&u);
                assert_eq!(got, sch(&w.padded(n)), "reaching {w} in S_{n}");
            }
        }
    }

    /// `u` with `w₀ · u = w`, i.e. `u = w₀⁻¹ w`.
    fn compose(w0: &Perm, w: &Perm, n: u32) -> Perm {
        let a = w0.inverse().padded(n);
        let b = w.padded(n);
        // (a then b): position i -> a[b[i]-1]? we want w0 * u = w  =>  u = w0^{-1} w
        let v: Vec<u32> = (0..n as usize).map(|i| a[(b[i] - 1) as usize]).collect();
        // that is w composed with w0^{-1} in the one-line convention used by
        // transpose(); verify by the length identity the caller asserts.
        Perm::new(v).unwrap()
    }

    // --- layer 4: engines and round trips ----------------------------------

    #[test]
    fn from_polynomial_inverts_expand() {
        for n in 0..=6u32 {
            for w in crate::permutation::tests::all_perms(n) {
                let f = sch(&w.padded(n));
                assert_eq!(Schubert::from_polynomial(&f.expand()), f, "{w}");
            }
        }
    }

    #[test]
    fn from_polynomial_on_sums() {
        let f = sch(&[1, 3, 2])
            .add(&sch(&[3, 2, 1]))
            .add(&sch(&[1, 4, 2, 3]).scale(&3));
        assert_eq!(Schubert::from_polynomial(&f.expand()), f);
    }

    #[test]
    fn dimension_is_the_expansion_size() {
        for n in 0..=6u32 {
            for w in crate::permutation::tests::all_perms(n) {
                let total: i64 = sch(&w.padded(n)).expand().iter().map(|(_, c)| c).sum();
                assert_eq!(dimension(&w) as i64, total, "{w}");
            }
        }
    }

    #[test]
    fn monk_matches_expanded_multiplication() {
        for n in 0..=5u32 {
            for w in crate::permutation::tests::all_perms(n) {
                let f = sch(&w.padded(n));
                for i in 1..=(n + 1) {
                    let via_monk = f.mul_variable(i).expand();
                    let mut want: BTreeMap<Expo, i64> = BTreeMap::new();
                    for (e, c) in f.expand() {
                        let mut e2 = e;
                        if e2.len() < i as usize {
                            e2.resize(i as usize, 0);
                        }
                        e2[(i - 1) as usize] += 1;
                        *want.entry(strip(e2)).or_insert(0) += c;
                    }
                    let want: Vec<(Expo, i64)> = want.into_iter().collect();
                    assert_eq!(via_monk, want, "x_{i} · S_{w}");
                }
            }
        }
    }

    #[test]
    fn product_is_multiplication_of_polynomials() {
        for n in 0..=4u32 {
            for u in crate::permutation::tests::all_perms(n) {
                for v in crate::permutation::tests::all_perms(n) {
                    let a = sch(&u.padded(n));
                    let b = sch(&v.padded(n));
                    let got = a.mul(&b).expand();
                    let want = poly_mul(&a.expand(), &b.expand());
                    assert_eq!(got, want, "S_{u} · S_{v}");
                }
            }
        }
    }

    fn poly_mul(a: &[(Expo, i64)], b: &[(Expo, i64)]) -> Vec<(Expo, i64)> {
        let mut acc: BTreeMap<Expo, i64> = BTreeMap::new();
        for (e, c) in a {
            for (f, d) in b {
                let n = e.len().max(f.len());
                let mut g = vec![0u32; n];
                for (i, slot) in g.iter_mut().enumerate() {
                    *slot = e.get(i).copied().unwrap_or(0) + f.get(i).copied().unwrap_or(0);
                }
                *acc.entry(strip(g)).or_insert(0) += c * d;
            }
        }
        acc.retain(|_, c| *c != 0);
        acc.into_iter().collect()
    }

    /// The engines agree. E3 merges states that E1 visits separately, so
    /// a merge that is subtly wrong shows up here and nowhere else.
    #[test]
    fn e1_and_e3_agree_exhaustively() {
        for n in 0..=5u32 {
            for u in crate::permutation::tests::all_perms(n) {
                for v in crate::permutation::tests::all_perms(n) {
                    let a = sch(&u.padded(n));
                    let b = sch(&v.padded(n));
                    assert_eq!(a.mul_e3(&b), a.mul_naive(&b), "S_{u} · S_{v}");
                }
            }
        }
    }

    /// The verified example, and the point of the name: `F_w` is the
    /// *Stanley* symmetric function, not `S_w` rewritten.
    #[test]
    fn stanley_hand_values() {
        use crate::sym::SymFn;
        let f: Schur<i64> = stanley(&p(&[2, 1, 4, 3]));
        assert_eq!(f.coeff(&Partition::new([2])), 1);
        assert_eq!(f.coeff(&Partition::new([1, 1])), 1);
        assert_eq!(f.terms().len(), 2);
        // ...and S_2143 is not even symmetric, so F_w != S_w here
        assert_ne!(sch(&[2, 1, 4, 3]).expand().len(), 1);
        // identity
        let e: Schur<i64> = stanley(&Perm::identity());
        assert_eq!(e.coeff(&Partition::new([] as [u32; 0])), 1);
    }

    /// Grassmannian `w` is vexillary: a single Schur term, the shape being the
    /// same one `grassmannian_perm` was built from. Ties the bridge into `Sym`
    /// back to the convention pinned against Sage.
    #[test]
    fn stanley_of_grassmannian_is_one_schur() {
        use crate::sym::SymFn;
        for k in 1..=4u32 {
            for size in 0..=7u32 {
                for lam in crate::partition::partitions_of(size) {
                    if lam.len() as u32 > k {
                        continue;
                    }
                    let w = grassmannian_perm(&lam, k).unwrap();
                    let f: Schur<i64> = stanley(&w);
                    assert_eq!(f.terms().len().max(1), 1, "{w} -> {f:?}");
                    assert_eq!(f.coeff(&lam), 1, "{w} shape {lam:?}");
                }
            }
        }
    }

    /// `F_w` is Schur-positive, and invariant under the `1×w` shift the
    /// recursion relies on — a law test for the one move that is not a
    /// transition step.
    #[test]
    fn stanley_is_positive_and_shift_invariant() {
        use crate::sym::SymFn;
        for n in 0..=6u32 {
            for w in crate::permutation::tests::all_perms(n) {
                let f: Schur<i64> = stanley(&w);
                for (lam, c) in f.terms() {
                    assert!(*c > 0, "{w}: negative coefficient at {lam:?}");
                    assert_eq!(lam.size(), w.length(), "{w}: wrong degree");
                }
                let g: Schur<i64> = stanley(&shift_up(&w));
                assert_eq!(f, g, "shift invariance at {w}");
            }
        }
    }

    /// The pruned query must agree with reading the coefficient off the whole
    /// product — exhaustively, including all the zero answers, since pruning
    /// bugs show up as spurious zeros.
    #[test]
    fn schubert_coeff_agrees_with_the_full_product() {
        for n in 0..=4u32 {
            for u in crate::permutation::tests::all_perms(n) {
                for v in crate::permutation::tests::all_perms(n) {
                    let full = sch(&u.padded(n)).mul_naive(&sch(&v.padded(n)));
                    for w in crate::permutation::tests::all_perms(n + 2) {
                        let got: i64 = schubert_coeff(&u, &v, &w);
                        assert_eq!(got, full.coeff(&w), "c^{w}_{{{u},{v}}}");
                    }
                }
            }
        }
    }

    /// The pruning rule behind a single-coefficient query: every `w` in the
    /// support of `S_u · S_v` satisfies `u ≤ w` and `v ≤ w` in Bruhat order,
    /// and every Monk cover moves strictly up. Both are checked here because
    /// the query's soundness rests on them and neither is verified elsewhere.
    #[test]
    fn product_support_lies_above_both_factors() {
        for n in 0..=5u32 {
            for u in crate::permutation::tests::all_perms(n) {
                for v in crate::permutation::tests::all_perms(n) {
                    for w in sch(&u.padded(n))
                        .mul_naive(&sch(&v.padded(n)))
                        .terms()
                        .keys()
                    {
                        assert!(u.bruhat_le(w), "{u} !<= {w} in S_{u}·S_{v}");
                        assert!(v.bruhat_le(w), "{v} !<= {w} in S_{u}·S_{v}");
                    }
                }
            }
        }
    }

    #[test]
    fn monk_covers_move_up_the_bruhat_order() {
        for n in 0..=5u32 {
            for w in crate::permutation::tests::all_perms(n) {
                for i in 1..=(n + 1) {
                    for (_, u) in w.covers_right(i) {
                        assert!(w.bruhat_le(&u), "{w} !<= {u}");
                    }
                    for (_, u) in w.covers_left(i) {
                        assert!(w.bruhat_le(&u), "{w} !<= {u} (left)");
                    }
                }
            }
        }
    }

    #[test]
    fn bruhat_is_a_partial_order_with_known_bounds() {
        for n in 0..=5u32 {
            let w0 = p(&(1..=n).rev().collect::<Vec<u32>>());
            for a in crate::permutation::tests::all_perms(n) {
                assert!(a.bruhat_le(&a), "reflexive {a}");
                assert!(Perm::identity().bruhat_le(&a), "id <= {a}");
                assert!(a.bruhat_le(&w0), "{a} <= w0");
                for b in crate::permutation::tests::all_perms(n) {
                    if a.bruhat_le(&b) && b.bruhat_le(&a) {
                        assert_eq!(a, b, "antisymmetry {a} {b}");
                    }
                    if a.bruhat_le(&b) {
                        assert!(a.length() <= b.length(), "{a} <= {b} but longer");
                    }
                }
            }
        }
    }

    /// Engine agreement again, for the third. E2 recurses on a completely different
    /// decomposition (transition, not the peel), so agreement here is a real
    /// cross-check rather than a restatement.
    #[test]
    fn e2_agrees_with_e1_exhaustively() {
        for n in 0..=5u32 {
            for u in crate::permutation::tests::all_perms(n) {
                for v in crate::permutation::tests::all_perms(n) {
                    let a = sch(&u.padded(n));
                    let b = sch(&v.padded(n));
                    assert_eq!(a.mul_e2(&b), a.mul_naive(&b), "S_{u} · S_{v}");
                }
            }
        }
    }

    /// E2 recurses on the *second* factor while E1/E3 pick a side, so the two
    /// argument orders exercise different trees and must still agree.
    #[test]
    fn e2_is_order_independent() {
        for n in 2..=5u32 {
            for u in crate::permutation::tests::all_perms(n) {
                for v in crate::permutation::tests::all_perms(n) {
                    let a = sch(&u.padded(n));
                    let b = sch(&v.padded(n));
                    assert_eq!(a.mul_e2(&b), b.mul_e2(&a), "S_{u} · S_{v}");
                }
            }
        }
    }

    #[test]
    fn e2_agrees_on_staircases() {
        for k in 2..=4u32 {
            let w: Vec<u32> = (1..=k)
                .map(|i| 2 * i)
                .chain((1..=k).map(|i| 2 * i - 1))
                .collect();
            let a = sch(&w);
            assert_eq!(a.mul_e2(&a), a.mul_naive(&a), "stair{k}^2");
        }
    }

    /// E1 and E3 agree on *sums*, and on mixed sizes — the shape that caught
    /// the peel memo's missing level key.
    #[test]
    fn e1_and_e3_agree_on_sums() {
        let f = sch(&[1, 3, 2])
            .add(&sch(&[3, 2, 1]))
            .add(&sch(&[1, 4, 2, 3]).scale(&3));
        let g = sch(&[2, 4, 1, 3]).add(&sch(&[2, 1]).scale(&2));
        assert_eq!(f.mul_e3(&g), f.mul_naive(&g));
        assert_eq!(g.mul_e3(&f), f.mul_naive(&g));
    }

    /// Bigger than the exhaustive sweep reaches, where the DAG actually
    /// compresses (`stair4` has 64 pipe dreams against 114 states).
    #[test]
    fn e1_and_e3_agree_on_staircases() {
        for k in 2..=4u32 {
            let w: Vec<u32> = (1..=k)
                .map(|i| 2 * i)
                .chain((1..=k).map(|i| 2 * i - 1))
                .collect();
            let a = sch(&w);
            assert_eq!(a.mul_e3(&a), a.mul_naive(&a), "stair{k}^2");
        }
    }

    #[test]
    fn peel_states_counts_the_dag() {
        // matches scripts/spec_schubert_peel.py's states_ps column (the corrected
        // one -- the first version of that script undercounted it at 114/33)
        let stair4: Vec<u32> = vec![2, 4, 6, 8, 1, 3, 5, 7];
        assert_eq!(peel_states(&p(&stair4)), 158);
        let stair3: Vec<u32> = vec![2, 4, 6, 1, 3, 5];
        assert_eq!(peel_states(&p(&stair3)), 38);
        assert_eq!(dimension(&p(&stair4)), 64);
        assert_eq!(dimension(&p(&stair3)), 8);
    }

    #[test]
    fn product_is_commutative_and_unital() {
        for n in 0..=4u32 {
            for u in crate::permutation::tests::all_perms(n) {
                let a = sch(&u.padded(n));
                assert_eq!(a.mul(&Schubert::one()), a, "{u} · 1");
                for v in crate::permutation::tests::all_perms(n) {
                    let b = sch(&v.padded(n));
                    assert_eq!(a.mul(&b), b.mul(&a), "{u} · {v}");
                }
            }
        }
    }

    /// Structure constants are non-negative — the geometric fact, and the
    /// cheapest possible smoke test that the signed Monk intermediates cancel
    /// the way they must.
    #[test]
    fn structure_constants_are_non_negative() {
        for n in 0..=5u32 {
            for u in crate::permutation::tests::all_perms(n) {
                for v in crate::permutation::tests::all_perms(n) {
                    for (w, c) in sch(&u.padded(n)).mul(&sch(&v.padded(n))).terms() {
                        assert!(*c > 0, "c^{w}_{{{u},{v}}} = {c}");
                    }
                }
            }
        }
    }

    #[test]
    fn pairing_is_the_top_coefficient() {
        // <S_u, S_v> = 1 iff v = w0 u, on Fl(3)
        let n = 3u32;
        let w0 = p(&[3, 2, 1]);
        for u in crate::permutation::tests::all_perms(n) {
            for v in crate::permutation::tests::all_perms(n) {
                let got = sch(&u.padded(n)).pairing(&sch(&v.padded(n)), n);
                let w0u: Vec<u32> = (0..n as usize)
                    .map(|i| w0.padded(n)[(u.padded(n)[i] - 1) as usize])
                    .collect();
                let want = i64::from(v == Perm::new(w0u).unwrap());
                assert_eq!(got, want, "<S_{u}, S_{v}>");
            }
        }
    }

    /// Stability: padding the input permutation changes nothing.
    #[test]
    fn stability_under_padding() {
        let a = sch(&[3, 1, 4, 2]);
        let b = sch(&[3, 1, 4, 2, 5, 6]);
        assert_eq!(a, b);
        assert_eq!(a.expand(), b.expand());
        assert_eq!(a.mul(&a), b.mul(&b));
    }
}
