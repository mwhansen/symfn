//! The basis-change engine: conversions between all five classical bases,
//! routed through the Schur hub.
//!
//! Each basis implements [`ToSchur`] (expand into the Schur basis) and
//! [`FromSchur`] (contract a Schur element into this basis); [`convert`] then
//! composes them to reach any ordered pair. The per-basis algorithms are the
//! classical ones:
//!
//! | conversion            | method                                  | ring    |
//! |-----------------------|-----------------------------------------|---------|
//! | h → s, e → s          | products of one-row / one-column Schurs | ℤ       |
//! | s → h                 | Jacobi–Trudi, by signed permutations    | ℤ       |
//! | s → e                 | dual Jacobi–Trudi, likewise             | ℤ       |
//! | p → s                 | iterated Murnaghan–Nakayama             | ℤ       |
//! | s → p                 | s_λ = Σ z_μ⁻¹ χ^λ(μ) p_μ                 | **ℚ**   |
//! | s → m                 | Kostka numbers                          | ℤ       |
//! | m → s                 | Muir's rule                             | ℤ       |
//!
//! Only s → p needs a [`Field`]; every other path stays exact over ℤ.
//!
//! Three of these were rewritten after a degree ladder against Sage
//! (`scripts/compare_sage.py`) showed them *scaling* badly rather than merely
//! being slow. That distinction is the reason the ladder exists: at a single
//! size each looked like an acceptable constant factor, and s → e and m → s were
//! both **faster than Sage at degree 8** while losing badly by degree 20.

use std::collections::HashMap;

use crate::character::character_in;
use crate::coeff::{Field, Ring};
use crate::fasthash::Map;
use crate::kostka::kostka;
use crate::memo::{inverse_kostka_row_cached, lex_parts_cached, partitions_cached};
use crate::partition::Partition;
use crate::sym::{Elementary, Homogeneous, Monomial, PowerSum, Schur, SymAlgebra, SymFn};

/// Expand `self` into the Schur basis.
pub trait ToSchur<C: Ring> {
    fn to_schur(&self) -> Schur<C>;
}

/// Contract a Schur element into `Self`'s basis.
pub trait FromSchur<C: Ring>: Sized {
    fn from_schur(s: &Schur<C>) -> Self;
}

/// Convert between any two bases by composing through Schur.
pub fn convert<C, A, B>(a: &A) -> B
where
    C: Ring,
    A: ToSchur<C>,
    B: FromSchur<C>,
{
    B::from_schur(&a.to_schur())
}

// --- Schur: the hub, identity both ways -------------------------------------

impl<C: Ring> ToSchur<C> for Schur<C> {
    fn to_schur(&self) -> Schur<C> {
        self.clone()
    }
}
impl<C: Ring> FromSchur<C> for Schur<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        s.clone()
    }
}

// --- the Jacobi–Trudi determinants ------------------------------------------

/// The determinant det(y_{c_i − i + j})_{i,j} expanded over a **free
/// multiplicative** basis, as a signed list of partitions.
///
/// Both Jacobi–Trudi determinants have this shape, and both are taken over a
/// basis (h or e) whose elements are indexed by partitions and multiply by
/// *concatenating* indices: y_a · y_b = y_{a ∪ b}. So a single permutation term
///
/// ```text
///   sgn(w) · ∏_i y_{c_i − i + w(i)}
/// ```
///
/// is one signed **monomial** y_μ, with μ the sorted multiset {c_i − i + w(i)}.
/// The whole determinant is therefore a signed count of permutations grouped by
/// that multiset — no polynomial arithmetic anywhere.
///
/// That is the entire fix. This used to be a generic Laplace expansion over the
/// symmetric-function algebra, which cloned an (n−1)×(n−1) matrix of
/// *polynomials* at every node and did a full polynomial add and multiply per
/// term. Its cost was factorial in the matrix size — and the matrix size is
/// ℓ(λ) for s → h but **λ₁** for s → e, since that one is built from the
/// conjugate. Same helper, opposite behaviour: s → h stayed fast on the wide-
/// but-shallow shapes a degree ladder produces while s → e crossed over and
/// fell behind Sage past degree 16.
///
/// Enumerating permutations directly also prunes where the determinant is
/// sparse: a negative index means the entry is zero, so that whole subtree is
/// skipped rather than multiplied out.
fn jt_terms(c: &[u32]) -> Vec<(Partition, i64)> {
    if c.is_empty() {
        return vec![(Partition::default(), 1)];
    }
    let mut acc: HashMap<Partition, i64> = HashMap::new();
    let mut used = vec![false; c.len()];
    let mut idx: Vec<u32> = Vec::with_capacity(c.len());
    jt_rec(c, 0, &mut used, &mut idx, 1, &mut acc);
    acc.into_iter().filter(|(_, s)| *s != 0).collect()
}

fn jt_rec(
    c: &[u32],
    i: usize,
    used: &mut [bool],
    idx: &mut Vec<u32>,
    sign: i64,
    acc: &mut HashMap<Partition, i64>,
) {
    if i == c.len() {
        // `Partition::new` drops the zero indices (y_0 is the unit) and sorts.
        *acc.entry(Partition::new(idx.iter().copied())).or_insert(0) += sign;
        return;
    }
    // Choosing w(i) = j inverts against every still-unused column below j —
    // those are all assigned to rows after i. Counting them as we scan gives
    // the permutation sign incrementally, with no cycle decomposition.
    let mut below = 0i64;
    for j in 0..c.len() {
        if used[j] {
            continue;
        }
        let k = c[i] as i64 - i as i64 + j as i64;
        if k >= 0 {
            used[j] = true;
            idx.push(k as u32);
            let s = if below % 2 == 0 { sign } else { -sign };
            jt_rec(c, i + 1, used, idx, s, acc);
            idx.pop();
            used[j] = false;
        }
        below += 1;
    }
}

/// Assemble a signed partition list into a basis element.
fn from_jt<C: Ring, S: SymAlgebra<C>>(terms: Vec<(Partition, i64)>) -> S {
    let mut out = S::zero();
    for (mu, s) in terms {
        out.add_term(mu, C::from_i64(s));
    }
    out
}

// --- Homogeneous <-> Schur --------------------------------------------------

impl<C: Ring> ToSchur<C> for Homogeneous<C> {
    fn to_schur(&self) -> Schur<C> {
        // h_λ = ∏_i h_{λ_i} = ∏_i s_{(λ_i)}.
        let mut out = Schur::zero();
        for (lambda, c) in self.terms() {
            let mut prod = Schur::unit();
            for &part in lambda.parts() {
                prod = prod.mul(&Schur::monomial(Partition::new([part]), C::one()));
            }
            out = out.add(&prod.scale(c));
        }
        out
    }
}

impl<C: Ring> FromSchur<C> for Homogeneous<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = det(h_{λ_i − i + j}) (Jacobi–Trudi).
        let mut out = Homogeneous::zero();
        for (lambda, c) in s.terms() {
            out = out.add(&jacobi_trudi::<C>(lambda).scale(c));
        }
        out
    }
}

/// s_λ = det(h_{λ_i − i + j}); the matrix is ℓ(λ)×ℓ(λ).
fn jacobi_trudi<C: Ring>(lambda: &Partition) -> Homogeneous<C> {
    from_jt(jt_terms(lambda.parts()))
}

// --- Elementary <-> Schur ---------------------------------------------------

impl<C: Ring> ToSchur<C> for Elementary<C> {
    fn to_schur(&self) -> Schur<C> {
        // e_λ = ∏_i e_{λ_i} = ∏_i s_{(1^{λ_i})} (single columns).
        let mut out = Schur::zero();
        for (lambda, c) in self.terms() {
            let mut prod = Schur::unit();
            for &part in lambda.parts() {
                let column = Partition::new(std::iter::repeat(1).take(part as usize));
                prod = prod.mul(&Schur::monomial(column, C::one()));
            }
            out = out.add(&prod.scale(c));
        }
        out
    }
}

impl<C: Ring> FromSchur<C> for Elementary<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = det(e_{λ'_i − i + j}) (dual Jacobi–Trudi).
        let mut out = Elementary::zero();
        for (lambda, c) in s.terms() {
            out = out.add(&dual_jacobi_trudi::<C>(lambda).scale(c));
        }
        out
    }
}

/// s_λ = det(e_{λ'_i − i + j}); the matrix is λ₁×λ₁, since it is built from the
/// conjugate. That is why this direction, and not [`jacobi_trudi`], was the one
/// that fell behind on wide shapes.
fn dual_jacobi_trudi<C: Ring>(lambda: &Partition) -> Elementary<C> {
    from_jt(jt_terms(lambda.conjugate().parts()))
}

// --- PowerSum <-> Schur -----------------------------------------------------

impl<C: Ring> ToSchur<C> for PowerSum<C> {
    fn to_schur(&self) -> Schur<C> {
        let mut out = Schur::zero();
        // Terms are batched by degree, because p_μ and p_ν share Murnaghan–
        // Nakayama work whenever they share parts — and expanding each one
        // separately from ∅ throws all of that away. The β-mask width is |μ|,
        // so only equal degrees can share a sweep.
        let mut by_degree: HashMap<u32, Vec<(&Partition, &C)>> = HashMap::new();
        for (mu, c) in self.terms() {
            let n = mu.size();
            if n as usize <= MASK_LIMIT {
                by_degree.entry(n).or_default().push((mu, c));
            } else {
                // Degree past the β-mask width: fall back to characters, which
                // are exact in `C` and so stay correct for bignum rings.
                for lambda in partitions_cached(n).iter() {
                    let chi = character_in::<C>(lambda, mu);
                    if !chi.is_zero() {
                        out.add_term(lambda.clone(), chi.mul(c));
                    }
                }
            }
        }
        let mut degrees: Vec<u32> = by_degree.keys().copied().collect();
        degrees.sort_unstable();
        for n in degrees {
            let mut items = by_degree.remove(&n).unwrap();
            // Descending part order, so the longest common prefixes are shared.
            items.sort_by(|a, b| a.0.parts().cmp(b.0.parts()));
            let l = n as usize;
            let mut root: Map<u64, i128> = Map::default();
            if l == 0 {
                for (_, c) in &items {
                    out.add_term(Partition::default(), (*c).clone());
                }
                continue;
            }
            root.insert((1u64 << l) - 1, 1);
            // Accumulate on the β-mask, not on a Partition. Every leaf of the
            // traversal touches the whole frontier, so keying by partition
            // allocated and sorted a fresh Vec — and hashed a heap key — once
            // per (μ, mask) pair to produce a few dozen distinct terms. Masks
            // are u64, and the partitions get built once at the end.
            let mut acc: Map<u64, C> = Map::default();
            p_expand_shared(&items, 0, &root, &mut acc);
            for (mask, c) in acc {
                if !c.is_zero() {
                    out.add_term(mask_to_partition(mask, l), c);
                }
            }
        }
        out
    }
}

/// Expand a batch of p_μ of one degree in a single traversal, sharing the
/// Murnaghan–Nakayama sweep across every common prefix.
///
/// `items` is sorted by part sequence, so partitions agreeing in their first
/// `depth` parts are contiguous; each such run continues from *one* frontier
/// instead of rebuilding it. Plethysm is the case that motivates this: it
/// finishes by converting a p-element of degree d·e with dozens of terms, and
/// that conversion was ~99% of its runtime.
fn p_expand_shared<C: Ring>(
    items: &[(&Partition, &C)],
    depth: usize,
    frontier: &Map<u64, i128>,
    out: &mut Map<u64, C>,
) {
    let mut i = 0;
    // Partitions that end here: emit the frontier scaled by their coefficient.
    while i < items.len() && items[i].0.len() == depth {
        for (&mask, &chi) in frontier {
            if chi != 0 {
                let term = C::from_i128(chi).mul(items[i].1);
                out.entry(mask).or_insert_with(C::zero).add_assign(&term);
            }
        }
        i += 1;
    }
    // The rest are grouped by their next part, each group sharing one step.
    while i < items.len() {
        let k = items[i].0.part(depth);
        let start = i;
        while i < items.len() && items[i].0.part(depth) == k {
            i += 1;
        }
        let next = p_step(frontier, k);
        p_expand_shared(&items[start..i], depth + 1, &next, out);
    }
}

/// p_μ in the Schur basis, by **iterated Murnaghan–Nakayama** rather than by
/// evaluating characters.
///
/// `p_μ = Σ_λ χ^λ(μ) s_λ`, and the obvious implementation asks for χ^λ(μ) once
/// per λ — p(n) independent recursions for one p_μ. But MN is itself a
/// multiplication rule,
///
/// ```text
///   p_k · s_λ = Σ (−1)^{ht(ξ)} s_{λ ∪ ξ},   ξ a k-rim-hook added to λ,
/// ```
///
/// so multiplying successively by p_{μ₁}, p_{μ₂}, … builds the entire expansion
/// in ℓ(μ) passes and never computes a character at all. Same "produce the whole
/// answer in one sweep rather than query it entry by entry" shape as the Kostka
/// and coproduct fixes.
///
/// Rim hooks are handled in β-numbers (first-column hook lengths, strictly
/// decreasing): adding a k-rim-hook is replacing some β by β+k when that value
/// is free, and the height is how many β lie strictly between them. Using a
/// fixed β-length of n keeps every intermediate comparable — a partition of n
/// has at most n rows, so no representable shape is lost.
///
/// Accumulation is in **i128**, not the caller's ring, and that is a deliberate
/// and safe narrowing rather than the usual truncation hazard. Every value here
/// is a character χ^λ(μ), and |χ^λ(μ)| ≤ f^λ ≤ √(n!). This function already
/// declines when l > 32 (the β-mask does not fit), and √(32!) ≈ 1.6·10¹⁸ — so on
/// every input it *accepts*, i128 cannot overflow. Callers past that width take
/// the character fallback, which stays exact in `C` for bignum rings.
///
/// The narrowing matters because the ring is the hot loop. Plethysm runs this
/// over ℚ, where each rim hook cost a rational add — a gcd — plus a temporary
/// from negating the coefficient, to combine two integers. That made p → s
/// ~99% of plethysm's runtime.
///
/// Retained as the **reference form**: production now takes the batched
/// [`p_expand_shared`] path, which shares this sweep across every p_μ of a
/// degree, and `batched_p_expansion_matches_per_term` holds the two to
/// agreement. The prefix grouping is the part of that batching that can quietly
/// go wrong, so it gets an oracle rather than trust.
#[cfg(test)]
fn p_expand(mu: &Partition) -> Option<Vec<(Partition, i128)>> {
    let l = mu.size() as usize;
    if l == 0 {
        return Some(vec![(Partition::default(), 1)]);
    }
    if l > MASK_LIMIT {
        return None;
    }
    // β-numbers of ∅ with l slots: {0, 1, …, l−1}.
    let mut cur: Map<u64, i128> = Map::default();
    cur.insert((1u64 << l) - 1, 1);
    for &k in mu.parts() {
        cur = p_step(&cur, k);
    }
    Some(
        cur.into_iter()
            .filter(|&(_, c)| c != 0)
            .map(|(mask, c)| (mask_to_partition(mask, l), c))
            .collect(),
    )
}

/// β values run from 0 to at most (l−1) + max part < 2l, so a 64-bit mask holds
/// the whole set for l ≤ 32. Past that the mask — which is what makes any of
/// this worth doing — no longer fits, and callers fall back.
const MASK_LIMIT: usize = 32;

/// One Murnaghan–Nakayama step: multiply a frontier of β-masks by p_k.
fn p_step(cur: &Map<u64, i128>, k: u32) -> Map<u64, i128> {
    // The frontier grows monotonically through a sweep, so a default-capacity
    // map rehashes several times per step. Sizing to the input is a floor on
    // the output, not a guess.
    let mut next: Map<u64, i128> =
        Map::with_capacity_and_hasher(cur.len() * 2, Default::default());
    for (&mask, &c) in cur {
        let mut rest = mask;
        while rest != 0 {
            let b = rest.trailing_zeros();
            rest &= rest - 1;
            let nb = b + k;
            if mask >> nb & 1 == 1 {
                continue; // that β is taken: no such rim hook
            }
            // Height = how many β lie strictly between b and b+k.
            let between = mask & (((1u64 << nb) - 1) ^ ((1u64 << (b + 1)) - 1));
            let m = (mask & !(1u64 << b)) | (1u64 << nb);
            let slot = next.entry(m).or_insert(0);
            if between.count_ones() % 2 == 0 {
                *slot += c;
            } else {
                *slot -= c;
            }
        }
    }
    next
}

/// Bits high-to-low are β₀ > β₁ > …, and λ_i = β_i − (l−1−i).
fn mask_to_partition(mask: u64, l: usize) -> Partition {
    let mut parts = Vec::with_capacity(l);
    let mut rest = mask;
    let mut i = 0usize;
    while rest != 0 {
        let b = 63 - rest.leading_zeros() as usize;
        rest &= !(1u64 << b);
        let part = b - (l - 1 - i);
        if part > 0 {
            parts.push(part as u32);
        }
        i += 1;
    }
    Partition::new(parts)
}

impl<C: Field> FromSchur<C> for PowerSum<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = Σ_μ z_μ⁻¹ χ^λ(μ) p_μ  (needs division by z_μ).
        let mut out = PowerSum::zero();
        for (lambda, c) in s.terms() {
            for mu in partitions_cached(lambda.size()).iter() {
                let chi = character_in::<C>(lambda, mu);
                if !chi.is_zero() {
                    let z_inv = C::from_u128(mu.z()).inv();
                    let coeff = c.mul(&chi).mul(&z_inv);
                    out.add_term(mu.clone(), coeff);
                }
            }
        }
        out
    }
}

// --- Monomial <-> Schur -----------------------------------------------------

impl<C: Ring> FromSchur<C> for Monomial<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = Σ_μ K_{λμ} m_μ.
        let mut out = Monomial::zero();
        for (lambda, c) in s.terms() {
            for mu in partitions_cached(lambda.size()).iter() {
                let k = kostka(lambda, mu);
                if k != 0 {
                    out.add_term(mu.clone(), C::from_u128(k).mul(c));
                }
            }
        }
        out
    }
}

impl<C: Ring> ToSchur<C> for Monomial<C> {
    fn to_schur(&self) -> Schur<C> {
        let mut out = Schur::zero();
        for (mu, c) in self.terms() {
            match muir_expand(mu) {
                Some(terms) => {
                    for (lambda, v) in terms {
                        out.add_term(lambda, C::from_i128(v).mul(c));
                    }
                }
                // Past the β-mask width: fall back to the row solve, which is
                // slow but has no degree ceiling.
                None => {
                    let parts = lex_parts_cached(mu.size());
                    let row = inverse_kostka_row_cached(mu, || inverse_kostka_row(&parts, mu));
                    for (i, lambda) in parts.iter().enumerate() {
                        if row[i] != 0 {
                            // `row` is i128 because that is what the solve is
                            // built in; inject at that width, not through i64.
                            out.add_term(lambda.clone(), C::from_i128(row[i]).mul(c));
                        }
                    }
                }
            }
        }
        out
    }
}

/// m_μ in the Schur basis, by **Muir's rule** — no Kostka numbers and no linear
/// solve.
///
/// In n variables the bialternant gives s_λ = a_{λ+δ}/a_δ, and multiplying by
/// m_μ = Σ_α x^α (over the distinct rearrangements α of μ) just shifts exponents:
///
/// ```text
///   m_μ · a_δ = Σ_α a_{α+δ},   so   m_μ = Σ_α ± s_{sort(α+δ) − δ}
/// ```
///
/// with `a_β = 0` when β repeats and `± a_{sorted β}` otherwise. Taking n = |μ|
/// slots loses nothing: (K⁻¹)_{μλ} ≠ 0 forces μ ⊵ λ, so ℓ(λ) ≤ |μ|.
///
/// Working in β-numbers β = α + δ, each part of μ is *added to a distinct slot*
/// of the initial set {0, 1, …, l−1}, which makes this the same machinery as
/// [`p_expand`] — a u64 mask, values moved one at a time, sign flipped by the
/// number of occupied values jumped over. The one difference from
/// Murnaghan–Nakayama is that MN may move the same value repeatedly while here
/// every slot receives at most one part.
///
/// **Slots are processed in decreasing initial value, and that is what makes it
/// fast.** At the moment slot v is handled, every slot above it is final and
/// every slot below still sits at its own value, which is < v < v+k. So a
/// collision at v+k can only be with an already-final value: it can never be
/// resolved later, and the branch dies immediately. Processing upward instead
/// would make the same test unsound, since the occupant might yet move. That
/// prune is the difference between exploring the C(l, ℓ(μ)) · (rearrangements)
/// placements and exploring only the surviving ones — the placements chain
/// together exactly as rim hooks do, so the survivors are few.
///
/// This replaces a forward solve of the unitriangular system K⁻¹K = I, which
/// needed O(p(n)²) Kostka numbers to read p(n) of them. That solve was the
/// single largest deficit in the library: 244× slower than Sage at degree 20 and
/// widening, because the improvements before it attacked the constant and left
/// the complexity alone.
fn muir_expand(mu: &Partition) -> Option<Vec<(Partition, i128)>> {
    let l = mu.size() as usize;
    if l == 0 {
        return Some(vec![(Partition::default(), 1)]);
    }
    // β values reach (l−1) + μ₁ ≤ 2l−1, so a 64-bit mask holds them for l ≤ 32.
    if l > 32 {
        return None;
    }
    // Parts grouped by value: choosing "a part equal to k" once is what makes
    // the rearrangements *distinct*, so equal parts are never double-counted.
    let mut avail: Vec<(u32, u32)> = Vec::new();
    for &k in mu.parts() {
        match avail.last_mut() {
            Some((v, n)) if *v == k => *n += 1,
            _ => avail.push((k, 1)),
        }
    }
    let mut acc: HashMap<Partition, i128> = HashMap::new();
    muir_rec(l, l as i32 - 1, (1u64 << l) - 1, &mut avail, mu.len(), 1, &mut acc);
    Some(acc.into_iter().filter(|(_, v)| *v != 0).collect())
}

fn muir_rec(
    l: usize,
    v: i32,
    mask: u64,
    avail: &mut [(u32, u32)],
    left: usize,
    sign: i128,
    acc: &mut HashMap<Partition, i128>,
) {
    if v < 0 {
        if left == 0 {
            // β sorted descending, then λ_i = β_i − (l − i).
            let mut lam = Vec::with_capacity(l);
            let mut i = 0usize;
            for b in (0..64u32).rev() {
                if mask >> b & 1 == 1 {
                    lam.push(b - (l - 1 - i) as u32);
                    i += 1;
                }
            }
            *acc.entry(Partition::new(lam)).or_insert(0) += sign;
        }
        return;
    }
    // Every remaining part needs its own slot, and v+1 slots remain.
    if left > v as usize + 1 {
        return;
    }
    // Leave slot v where it is. Safe without a collision test: every already
    // final value is > v.
    muir_rec(l, v - 1, mask, avail, left, sign, acc);

    if left == 0 {
        return;
    }
    for i in 0..avail.len() {
        if avail[i].1 == 0 {
            continue;
        }
        let nb = v as u32 + avail[i].0;
        if mask >> nb & 1 == 1 {
            continue; // permanent collision — see the note above
        }
        // Jumping an occupied value transposes the two, flipping the sign.
        let between = mask & (((1u64 << nb) - 1) ^ ((1u64 << (v + 1)) - 1));
        let s = if between.count_ones() % 2 == 0 { sign } else { -sign };
        avail[i].1 -= 1;
        muir_rec(l, v - 1, (mask & !(1 << v)) | (1 << nb), avail, left - 1, s, acc);
        avail[i].1 += 1;
    }
}

/// Row μ of the inverse Kostka matrix, indexed against `parts` (decreasing-lex).
///
/// In that order the Kostka matrix K is upper-unitriangular — `K[i][j] ≠ 0`
/// needs λᵢ ⊵ λⱼ, and dominance implies lex — so `K⁻¹K = I` restricted to row μ
/// solves forward:
///
/// ```text
///   w[j] = 1,   w[jj] = −Σ_{m=j}^{jj−1} w[m]·K[m][jj]   for jj > j
/// ```
///
/// Only this one row is ever needed: m_μ = Σ_λ (K⁻¹)_{μλ} s_λ. Building the
/// whole matrix and inverting it, as this used to, computed p(n)² Kostka numbers
/// to read p(n) of them — 2.0 s for a single degree-20 conversion, of which the
/// matrix was ~97%.
///
/// The `w[m] == 0` skip is the part that matters: a zero coefficient makes its
/// Kostka number irrelevant, so the call is never made rather than made and
/// multiplied by zero.
fn inverse_kostka_row(parts: &[Partition], mu: &Partition) -> Vec<i128> {
    let size = parts.len();
    let j = parts
        .iter()
        .position(|p| p == mu)
        .expect("μ ∈ partitions(|μ|)");
    let mut w = vec![0i128; size];
    w[j] = 1;
    for jj in (j + 1)..size {
        let mut acc = 0i128;
        for m in j..jj {
            if w[m] == 0 {
                continue;
            }
            acc += w[m] * kostka(&parts[m], &parts[jj]) as i128;
        }
        w[jj] = -acc;
    }
    w
}

// --- Monomial multiplication (routed through Schur) -------------------------

impl<C: Ring> Monomial<C> {
    /// Product in the monomial basis. Unlike p/e/h this is not multiplicative,
    /// so it is computed by expanding into Schur, multiplying, and contracting
    /// back — entirely over ℤ.
    pub fn mul(&self, other: &Self) -> Self {
        let prod = self.to_schur().mul(&other.to_schur());
        Monomial::from_schur(&prod)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    #[test]
    fn h_and_e_expand_into_schur() {
        // h_{11} = s_2 + s_{11}; e_2 = s_{11}.
        let h11: Homogeneous<i64> = Homogeneous::monomial(part(&[1, 1]), 1);
        let s = h11.to_schur();
        assert_eq!(s.coeff(&part(&[2])), 1);
        assert_eq!(s.coeff(&part(&[1, 1])), 1);

        let e2: Elementary<i64> = Elementary::monomial(part(&[2]), 1);
        let se = e2.to_schur();
        assert_eq!(se.coeff(&part(&[1, 1])), 1);
        assert_eq!(se.terms().len(), 1);
    }

    #[test]
    fn jacobi_trudi_inverts_expansion() {
        // s_{11} = h_{11} − h_2 (Jacobi–Trudi).
        let s11: Schur<i64> = Schur::monomial(part(&[1, 1]), 1);
        let h = Homogeneous::from_schur(&s11);
        assert_eq!(h.coeff(&part(&[1, 1])), 1);
        assert_eq!(h.coeff(&part(&[2])), -1);
    }

    /// The batched p → s must equal expanding each p_μ on its own and summing.
    /// Batching groups terms by shared prefix and continues one frontier per
    /// group; an off-by-one in that grouping would silently attribute a term to
    /// the wrong μ, which no round-trip test would catch.
    #[test]
    fn batched_p_expansion_matches_per_term() {
        for n in 1..=10u32 {
            let mus = partitions_cached(n);
            // A single element carrying *every* p_μ of the degree at once, which
            // is the shape plethysm produces and the case batching exists for.
            let mut elt: PowerSum<i128> = PowerSum::zero();
            for (i, mu) in mus.iter().enumerate() {
                elt.add_term(mu.clone(), (i as i128) + 1);
            }
            let got: Schur<i128> = elt.to_schur();

            let mut want: Schur<i128> = Schur::zero();
            for (i, mu) in mus.iter().enumerate() {
                for (lambda, chi) in p_expand(mu).expect("small degree") {
                    want.add_term(lambda, chi * ((i as i128) + 1));
                }
            }
            assert_eq!(got, want, "batched p → s at degree {n}");
        }
    }

    /// Muir's rule against the linear solve it replaced, on **every** μ up to
    /// degree 12 — the two are independent routes to the same row of K⁻¹, and
    /// the solve was the shipped implementation, so this is a real oracle rather
    /// than a self-consistency check.
    #[test]
    fn muir_agrees_with_the_inverse_kostka_solve() {
        for n in 0..=12u32 {
            let parts = lex_parts_cached(n);
            for mu in parts.iter() {
                let row = inverse_kostka_row(&parts, mu);
                let mut want: Vec<(Partition, i128)> = parts
                    .iter()
                    .zip(&row)
                    .filter(|(_, &v)| v != 0)
                    .map(|(l, &v)| (l.clone(), v))
                    .collect();
                let mut got = muir_expand(mu).expect("degree 12 is inside the mask");
                want.sort();
                got.sort();
                assert_eq!(got, want, "m_{mu} in the Schur basis");
            }
        }
    }

    /// The two shapes worked by hand when the rule was derived. `m_{21}` is the
    /// one that catches a wrong triangularity assumption: (K⁻¹)_{μλ} is nonzero
    /// for μ ⊵ λ, so λ may be *longer* than μ, and taking only ℓ(μ) slots would
    /// silently drop the s_{111} term and return the 2-variable answer.
    #[test]
    fn muir_matches_hand_computation() {
        let got = muir_expand(&part(&[2, 1])).unwrap();
        let mut got: Vec<_> = got.iter().map(|(l, c)| (l.parts().to_vec(), *c)).collect();
        got.sort();
        assert_eq!(got, vec![(vec![1, 1, 1], -2), (vec![2, 1], 1)]);

        // m_{(n)} = Σ_i (−1)^i s_{(n−i, 1^i)}: the hooks, alternating.
        let got = muir_expand(&part(&[5])).unwrap();
        let mut got: Vec<_> = got.iter().map(|(l, c)| (l.parts().to_vec(), *c)).collect();
        got.sort_by_key(|(l, _)| l.len());
        let want = vec![
            (vec![5], 1),
            (vec![4, 1], -1),
            (vec![3, 1, 1], 1),
            (vec![2, 1, 1, 1], -1),
            (vec![1, 1, 1, 1, 1], 1),
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn round_trips_are_identity() {
        // s → h → s and s → e → s and s → m → s recover the original.
        for parts in [&[2, 1][..], &[3, 1], &[2, 2], &[3, 2, 1]] {
            let s: Schur<i64> = Schur::monomial(part(parts), 1);
            let back_h: Schur<i64> = convert(&Homogeneous::from_schur(&s));
            let back_e: Schur<i64> = convert(&Elementary::from_schur(&s));
            let back_m: Schur<i64> = convert(&Monomial::from_schur(&s));
            assert_eq!(back_h, s, "s→h→s at {:?}", parts);
            assert_eq!(back_e, s, "s→e→s at {:?}", parts);
            assert_eq!(back_m, s, "s→m→s at {:?}", parts);
        }
    }

    #[test]
    fn power_sum_round_trip_over_rationals() {
        // s → p → s over ℚ is the identity.
        for parts in [&[2, 1][..], &[3], &[2, 2], &[3, 1]] {
            let s: Schur<Rational> =
                Schur::monomial(part(parts), Rational::one());
            let p: PowerSum<Rational> = PowerSum::from_schur(&s);
            let back: Schur<Rational> = p.to_schur();
            assert_eq!(back, s, "s→p→s at {:?}", parts);
        }
    }

    #[test]
    fn power_sum_two_is_s2_minus_s11() {
        // p_2 = s_2 − s_{11}.
        let p2: PowerSum<i64> = PowerSum::monomial(part(&[2]), 1);
        let s = p2.to_schur();
        assert_eq!(s.coeff(&part(&[2])), 1);
        assert_eq!(s.coeff(&part(&[1, 1])), -1);
    }

    #[test]
    fn monomial_multiplication_via_schur() {
        // m_1 · m_1 = 2 m_{11} + m_2.
        let m1: Monomial<i64> = Monomial::monomial(part(&[1]), 1);
        let prod = m1.mul(&m1);
        assert_eq!(prod.coeff(&part(&[1, 1])), 2);
        assert_eq!(prod.coeff(&part(&[2])), 1);
        assert_eq!(prod.terms().len(), 2);
    }
}
