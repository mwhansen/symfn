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

use crate::character::character;
use crate::coeff::{Field, Ring};
use crate::kostka::kostka;
use crate::memo::{inverse_kostka_cached, partitions_cached};
use crate::partition::{partitions_of, Partition};
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
        // p_μ = Σ_λ χ^λ(μ) s_λ  (integer characters).
        let mut out = Schur::zero();
        for (mu, c) in self.terms() {
            for lambda in partitions_cached(mu.size()).iter() {
                let chi = character(lambda, mu);
                if chi != 0 {
                    out.add_term(lambda.clone(), C::from_i64(chi).mul(c));
                }
            }
        }
        out
    }
}

impl<C: Field> FromSchur<C> for PowerSum<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = Σ_μ z_μ⁻¹ χ^λ(μ) p_μ  (needs division by z_μ).
        let mut out = PowerSum::zero();
        for (lambda, c) in s.terms() {
            for mu in partitions_cached(lambda.size()).iter() {
                let chi = character(lambda, mu);
                if chi != 0 {
                    let z_inv = C::from_u128(mu.z()).inv();
                    let coeff = c.mul(&C::from_i64(chi)).mul(&z_inv);
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
            let n = mu.size();
            let data = inverse_kostka_cached(n, || inverse_kostka(n));
            let (parts, kinv) = (&data.0, &data.1);
            let j = parts.iter().position(|p| p == mu).expect("μ ∈ partitions(n)");
            for (i, lambda) in parts.iter().enumerate() {
                let v = kinv[j][i];
                if v != 0 {
                    out.add_term(lambda.clone(), C::from_i64(v as i64).mul(c));
                }
            }
        }
        out
    }
}

/// Partitions of `n` in decreasing-lex order (a linear extension of dominance),
/// together with the integer inverse of the Kostka matrix in that order.
///
/// In decreasing-lex order the Kostka matrix is upper-unitriangular, so its
/// inverse is exact over ℤ and computed by back-substitution.
fn inverse_kostka(n: u32) -> (Vec<Partition>, Vec<Vec<i128>>) {
    let mut parts = partitions_of(n);
    parts.sort_by(|a, b| b.parts().cmp(a.parts())); // decreasing lex
    let size = parts.len();

    // K[i][j] = kostka(parts[i], parts[j]); upper-unitriangular.
    let k: Vec<Vec<i128>> = (0..size)
        .map(|i| (0..size).map(|j| kostka(&parts[i], &parts[j]) as i128).collect())
        .collect();

    // Invert the upper-unitriangular matrix: V with K·V = I.
    // V[i][i] = 1; for i < j, V[i][j] = −Σ_{k=i+1}^{j} K[i][k]·V[k][j].
    let mut v = vec![vec![0i128; size]; size];
    for i in (0..size).rev() {
        v[i][i] = 1;
        for j in (i + 1)..size {
            let mut acc = 0i128;
            for kk in (i + 1)..=j {
                acc += k[i][kk] * v[kk][j];
            }
            v[i][j] = -acc;
        }
    }
    (parts, v)
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
