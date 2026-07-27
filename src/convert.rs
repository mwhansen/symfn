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
//! | s → h                 | Jacobi–Trudi determinant                | ℤ       |
//! | s → e                 | dual Jacobi–Trudi determinant           | ℤ       |
//! | p → s                 | Murnaghan–Nakayama characters           | ℤ       |
//! | s → p                 | s_λ = Σ z_μ⁻¹ χ^λ(μ) p_μ                 | **ℚ**   |
//! | s → m                 | Kostka numbers                          | ℤ       |
//! | m → s                 | inverse (unitriangular) Kostka          | ℤ       |
//!
//! Only s → p needs a [`Field`]; every other path stays exact over ℤ.

use std::collections::HashMap;

use crate::character::character_in;
use crate::coeff::{Field, Ring};
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

// --- determinant over a symmetric-function algebra (for Jacobi–Trudi) -------

/// Determinant of a small square matrix whose entries live in a symmetric-
/// function algebra, by Laplace expansion. Exponential in size, but the size is
/// ℓ(λ) (few parts), so fine as a baseline.
fn det<C: Ring, S: SymAlgebra<C> + Clone>(m: &[Vec<S>]) -> S {
    let n = m.len();
    match n {
        0 => S::unit(),
        1 => m[0][0].clone(),
        _ => {
            let mut acc = S::zero();
            for j in 0..n {
                let minor: Vec<Vec<S>> = m[1..]
                    .iter()
                    .map(|row| {
                        row.iter()
                            .enumerate()
                            .filter(|(c, _)| *c != j)
                            .map(|(_, x)| x.clone())
                            .collect()
                    })
                    .collect();
                let term = m[0][j].times(&det(&minor));
                if j % 2 == 0 {
                    acc = acc.add(&term);
                } else {
                    acc = acc.sub(&term);
                }
            }
            acc
        }
    }
}

/// h_k as a homogeneous element: 0 for k<0, the unit for k=0, h_{(k)} for k>0.
fn h_entry<C: Ring>(k: i64) -> Homogeneous<C> {
    match k {
        k if k < 0 => Homogeneous::zero(),
        0 => Homogeneous::unit(),
        k => Homogeneous::monomial(Partition::new([k as u32]), C::one()),
    }
}

/// e_k as an elementary element (same convention as [`h_entry`]).
fn e_entry<C: Ring>(k: i64) -> Elementary<C> {
    match k {
        k if k < 0 => Elementary::zero(),
        0 => Elementary::unit(),
        k => Elementary::monomial(Partition::new([k as u32]), C::one()),
    }
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

fn jacobi_trudi<C: Ring>(lambda: &Partition) -> Homogeneous<C> {
    let l = lambda.len();
    if l == 0 {
        return Homogeneous::unit();
    }
    let mat: Vec<Vec<Homogeneous<C>>> = (0..l)
        .map(|i| {
            (0..l)
                .map(|j| h_entry::<C>(lambda.part(i) as i64 - i as i64 + j as i64))
                .collect()
        })
        .collect();
    det(&mat)
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

fn dual_jacobi_trudi<C: Ring>(lambda: &Partition) -> Elementary<C> {
    let conj = lambda.conjugate();
    let m = conj.len();
    if m == 0 {
        return Elementary::unit();
    }
    let mat: Vec<Vec<Elementary<C>>> = (0..m)
        .map(|i| {
            (0..m)
                .map(|j| e_entry::<C>(conj.part(i) as i64 - i as i64 + j as i64))
                .collect()
        })
        .collect();
    det(&mat)
}

// --- PowerSum <-> Schur -----------------------------------------------------

impl<C: Ring> ToSchur<C> for PowerSum<C> {
    fn to_schur(&self) -> Schur<C> {
        let mut out = Schur::zero();
        for (mu, c) in self.terms() {
            match p_expand::<C>(mu) {
                Some(terms) => {
                    for (lambda, chi) in terms {
                        out.add_term(lambda, chi.mul(c));
                    }
                }
                // Degree past the β-mask width: fall back to characters, which
                // are exact in `C` and so stay correct for bignum rings.
                None => {
                    for lambda in partitions_cached(mu.size()).iter() {
                        let chi = character_in::<C>(lambda, mu);
                        if !chi.is_zero() {
                            out.add_term(lambda.clone(), chi.mul(c));
                        }
                    }
                }
            }
        }
        out
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
/// Accumulation is in `C`, not `i128`, so a bignum coefficient ring stays exact
/// past the i128 character ceiling (n ≈ 58) exactly as the character-based
/// version did.
fn p_expand<C: Ring>(mu: &Partition) -> Option<Vec<(Partition, C)>> {
    let l = mu.size() as usize;
    if l == 0 {
        return Some(vec![(Partition::default(), C::one())]);
    }
    // β values run from 0 to at most (l−1) + max part < 2l, so a 64-bit mask
    // holds the whole set for l ≤ 32. Beyond that, decline and let the caller
    // fall back — the mask is what makes this worth doing at all.
    if l > 32 {
        return None;
    }
    // β-numbers of ∅ with l slots: {0, 1, …, l−1}.
    let mut cur: HashMap<u64, C> = HashMap::new();
    cur.insert((1u64 << l) - 1, C::one());

    for &k in mu.parts() {
        let k = k as u32;
        let mut next: HashMap<u64, C> = HashMap::new();
        for (&mask, c) in &cur {
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
                let slot = next.entry(m).or_insert_with(C::zero);
                if between.count_ones() % 2 == 0 {
                    slot.add_assign(c);
                } else {
                    slot.add_assign(&c.neg());
                }
            }
        }
        cur = next;
    }

    Some(
        cur.into_iter()
            .filter(|(_, c)| !c.is_zero())
            .map(|(mask, c)| {
                // Bits high-to-low are β₀ > β₁ > …, and λ_i = β_i − (l−1−i).
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
                (Partition::new(parts), c)
            })
            .collect(),
    )
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
        // m_μ = Σ_λ (K⁻¹)_{μλ} s_λ, from the unitriangular inverse of Kostka.
        // The inverse-Kostka data is cached globally per degree.
        let mut out = Schur::zero();
        for (mu, c) in self.terms() {
            let parts = lex_parts_cached(mu.size());
            let row = inverse_kostka_row_cached(mu, || inverse_kostka_row(&parts, mu));
            for (i, lambda) in parts.iter().enumerate() {
                let v = row[i];
                if v != 0 {
                    // `v` is i128 because that is what the matrix is built in;
                    // inject at that width rather than narrowing through i64.
                    out.add_term(lambda.clone(), C::from_i128(v).mul(c));
                }
            }
        }
        out
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
