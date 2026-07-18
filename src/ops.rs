//! Standard operations on symmetric functions: the ω involution and the Hall
//! inner product.
//!
//! Both have cheap native forms in a preferred basis (ω on Schur/power sums, the
//! inner product via Schur orthonormality) and a generic form for any basis,
//! obtained by routing through the Schur hub.

use crate::coeff::Ring;
use crate::convert::{FromSchur, ToSchur};
use crate::sym::{PowerSum, Schur, SymFn};

impl<C: Ring> Schur<C> {
    /// The ω involution in the Schur basis: ω(s_λ) = s_{λ'} (conjugate shape).
    pub fn omega(&self) -> Self {
        let mut out = Schur::zero();
        for (lambda, c) in self.terms() {
            out.add_term(lambda.conjugate(), c.clone());
        }
        out
    }
}

impl<C: Ring> PowerSum<C> {
    /// The ω involution in the power-sum basis: ω(p_λ) = ε_λ p_λ, where
    /// ε_λ = (−1)^{|λ|−ℓ(λ)}.
    pub fn omega(&self) -> Self {
        let mut out = PowerSum::zero();
        for (lambda, c) in self.terms() {
            let odd = (lambda.size() - lambda.len() as u32) % 2 == 1;
            out.add_term(lambda.clone(), if odd { c.neg() } else { c.clone() });
        }
        out
    }
}

/// The ω involution for any basis, via Schur: ω(x) = (from Schur)(ω(to Schur x)).
/// For example ω of an `h`-element lands back in the `h`-basis (ω(h_λ) = e_λ,
/// re-expressed in h).
pub fn omega<C, B>(x: &B) -> B
where
    C: Ring,
    B: ToSchur<C> + FromSchur<C>,
{
    B::from_schur(&x.to_schur().omega())
}

/// The Hall inner product ⟨a, b⟩, computed via Schur orthonormality
/// ⟨s_λ, s_μ⟩ = δ_{λμ}: expand both into Schur and take the coefficient dot
/// product. Works across bases (e.g. ⟨h_λ, m_μ⟩ = δ_{λμ}).
pub fn hall<C, A, B>(a: &A, b: &B) -> C
where
    C: Ring,
    A: ToSchur<C>,
    B: ToSchur<C>,
{
    let sa = a.to_schur();
    let sb = b.to_schur();
    let mut acc = C::zero();
    for (lambda, ca) in sa.terms() {
        if let Some(cb) = sb.terms().get(lambda) {
            acc.add_assign(&ca.mul(cb));
        }
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::partition::Partition;
    use crate::sym::{Homogeneous, Monomial};

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    #[test]
    fn omega_on_schur_conjugates() {
        let s3: Schur<i64> = Schur::monomial(part(&[3]), 1);
        assert_eq!(s3.omega().coeff(&part(&[1, 1, 1])), 1);
        // ω² = id
        assert_eq!(s3.omega().omega(), s3);
        // s_{21} is self-conjugate
        let s21: Schur<i64> = Schur::monomial(part(&[2, 1]), 1);
        assert_eq!(s21.omega(), s21);
    }

    #[test]
    fn omega_on_power_sums_signs() {
        // ω(p_2) = −p_2 (ε = (−1)^{2−1}); ω(p_{11}) = +p_{11} (ε = (−1)^{2−2}).
        let p2: PowerSum<i64> = PowerSum::monomial(part(&[2]), 1);
        assert_eq!(p2.omega().coeff(&part(&[2])), -1);
        let p11: PowerSum<i64> = PowerSum::monomial(part(&[1, 1]), 1);
        assert_eq!(p11.omega().coeff(&part(&[1, 1])), 1);
    }

    #[test]
    fn omega_via_schur_maps_h_to_e_expression() {
        // ω(h_2) = e_2 = h_{11} − h_2.
        let h2: Homogeneous<i64> = Homogeneous::monomial(part(&[2]), 1);
        let w = omega(&h2);
        assert_eq!(w.coeff(&part(&[1, 1])), 1);
        assert_eq!(w.coeff(&part(&[2])), -1);
    }

    #[test]
    fn hall_products_match_known_pairings() {
        // Schur orthonormality.
        let s21: Schur<i64> = Schur::monomial(part(&[2, 1]), 1);
        let s3: Schur<i64> = Schur::monomial(part(&[3]), 1);
        assert_eq!(hall(&s21, &s21), 1);
        assert_eq!(hall(&s21, &s3), 0);

        // ⟨p_λ, p_λ⟩ = z_λ.
        let p2: PowerSum<i64> = PowerSum::monomial(part(&[2]), 1);
        assert_eq!(hall(&p2, &p2), part(&[2]).z() as i64); // = 2

        // ⟨h_λ, m_μ⟩ = δ_{λμ}: h and m are dual bases.
        let h21: Homogeneous<i64> = Homogeneous::monomial(part(&[2, 1]), 1);
        let m21: Monomial<i64> = Monomial::monomial(part(&[2, 1]), 1);
        let m3: Monomial<i64> = Monomial::monomial(part(&[3]), 1);
        assert_eq!(hall(&h21, &m21), 1);
        assert_eq!(hall(&h21, &m3), 0);
    }
}
