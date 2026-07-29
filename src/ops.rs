//! Standard operations on symmetric functions: the ω involution and the Hall
//! inner product.
//!
//! Both have cheap native forms in a preferred basis (ω on Schur/power sums, the
//! inner product via Schur orthonormality) and a generic form for any basis,
//! obtained by routing through the Schur hub.

use crate::coeff::{QAlgebra, Ring};
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

/// The **internal (Kronecker) product** `a * b`.
///
/// The third product on symmetric functions, after the ordinary one and
/// plethysm (Macdonald I.7). Under the Frobenius characteristic it is the
/// tensor product of S_n representations, so
///
/// ```text
///   s_λ * s_μ = Σ_ν g^ν_{λμ} s_ν,     g^ν_{λμ} = ⟨χ^λ χ^μ, χ^ν⟩
/// ```
///
/// with the `g` the Kronecker coefficients — famously harder than
/// Littlewood–Richardson, and with no known positive combinatorial rule.
///
/// Computing them here is nonetheless almost free, because the internal product
/// is **diagonal in the power-sum basis**:
///
/// ```text
///   p_λ * p_μ = δ_{λμ} · z_λ · p_λ
/// ```
///
/// So it is s → p on both sides, a coefficientwise multiply weighted by z_λ, and
/// p → s back. Nothing enumerates anything. That the hard object falls out of a
/// diagonal basis is the whole point of keeping the power-sum route fast.
///
/// Requires a [`Field`] for the z_λ⁻¹ that s → p introduces; the result of two
/// Schur inputs is integral regardless.
///
/// Degrees need no special handling. A term survives only when the same λ occurs
/// on both sides, which forces |a| = |b| — so the product of elements of
/// different degrees is zero, exactly as the grading demands.
pub fn internal<C: QAlgebra>(a: &Schur<C>, b: &Schur<C>) -> Schur<C> {
    let pa: PowerSum<C> = PowerSum::from_schur(a);
    let pb: PowerSum<C> = PowerSum::from_schur(b);
    // Iterate the smaller side; the intersection is what contributes.
    let (small, large) = if pa.terms().len() <= pb.terms().len() {
        (&pa, &pb)
    } else {
        (&pb, &pa)
    };
    let mut acc = PowerSum::zero();
    for (lambda, cs) in small.terms() {
        let cl = match large.terms().get(lambda) {
            Some(c) => c,
            None => continue,
        };
        acc.add_term(lambda.clone(), cs.mul(cl).mul(&C::from_u128(lambda.z())));
    }
    acc.to_schur()
}

/// A single Kronecker coefficient g^ν_{λμ}.
///
/// Convenience over [`internal`]; computing one costs the same as computing the
/// whole product, since the power-sum route produces every ν at once.
pub fn kronecker<C: QAlgebra>(
    lambda: &crate::partition::Partition,
    mu: &crate::partition::Partition,
    nu: &crate::partition::Partition,
) -> C {
    let sl: Schur<C> = Schur::monomial(lambda.clone(), C::one());
    let sm: Schur<C> = Schur::monomial(mu.clone(), C::one());
    internal(&sl, &sm).coeff(nu)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::memo::partitions_cached;
    use crate::partition::Partition;
    use crate::sym::{Homogeneous, Monomial};

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn q(n: i128) -> Rational {
        Rational::new(n, 1)
    }

    fn sch(v: &[u32]) -> Schur<Rational> {
        Schur::monomial(part(v), q(1))
    }

    /// The S_3 character table, worked by hand: s_3 is the trivial character,
    /// s_{111} the sign, s_{21} the 2-dimensional standard one.
    #[test]
    fn kronecker_matches_s3_representation_theory() {
        // trivial ⊗ anything = anything
        assert_eq!(internal(&sch(&[3]), &sch(&[2, 1])), sch(&[2, 1]));
        // sign ⊗ sign = trivial
        assert_eq!(internal(&sch(&[1, 1, 1]), &sch(&[1, 1, 1])), sch(&[3]));
        // sign ⊗ standard = standard
        assert_eq!(internal(&sch(&[1, 1, 1]), &sch(&[2, 1])), sch(&[2, 1]));
        // standard ⊗ standard = trivial + standard + sign  (dim 2·2 = 1+2+1)
        let r = internal(&sch(&[2, 1]), &sch(&[2, 1]));
        assert_eq!(r.coeff(&part(&[3])), q(1));
        assert_eq!(r.coeff(&part(&[2, 1])), q(1));
        assert_eq!(r.coeff(&part(&[1, 1, 1])), q(1));
        assert_eq!(r.terms().len(), 3);
    }

    /// Against the character-theoretic definition, which shares no code with the
    /// power-sum route: g^ν_{λμ} = Σ_ρ χ^λ(ρ)·χ^μ(ρ)·χ^ν(ρ) / z_ρ.
    ///
    /// This is the real oracle. The implementation leans entirely on
    /// `p_λ * p_μ = δ z_λ p_λ` plus s ↔ p; if either the identity or a
    /// conversion were wrong, a self-consistent but false answer would follow.
    #[test]
    fn kronecker_matches_the_character_formula() {
        for n in 1..=6u32 {
            let parts = partitions_cached(n);
            for lambda in parts.iter() {
                for mu in parts.iter() {
                    let got = internal(
                        &Schur::monomial(lambda.clone(), q(1)),
                        &Schur::monomial(mu.clone(), q(1)),
                    );
                    for nu in parts.iter() {
                        let mut want = Rational::new(0, 1);
                        for rho in parts.iter() {
                            let t = crate::character::character(lambda, rho)
                                * crate::character::character(mu, rho)
                                * crate::character::character(nu, rho);
                            want.add_assign(&Rational::new(t, rho.z() as i128));
                        }
                        assert_eq!(got.coeff(nu), want, "g^{nu}_{{{lambda},{mu}}}");
                    }
                }
            }
        }
    }

    /// Structural facts the coefficients must satisfy: symmetry in all three
    /// indices, invariance under conjugating any two, and Σ_ν g·dim(ν) =
    /// dim(λ)·dim(μ) — the dimension of the tensor product.
    #[test]
    fn kronecker_symmetries_and_dimension_count() {
        for n in 1..=6u32 {
            let parts = partitions_cached(n);
            let ones = Partition::new(std::iter::repeat(1).take(n as usize));
            for lambda in parts.iter() {
                for mu in parts.iter() {
                    let g = internal(
                        &Schur::monomial(lambda.clone(), q(1)),
                        &Schur::monomial(mu.clone(), q(1)),
                    );
                    let swapped = internal(
                        &Schur::monomial(mu.clone(), q(1)),
                        &Schur::monomial(lambda.clone(), q(1)),
                    );
                    assert_eq!(g, swapped, "g symmetric in λ, μ");

                    // Conjugating two of the three indices leaves g unchanged
                    // (tensoring both factors by the sign representation).
                    let conj = internal(
                        &Schur::monomial(lambda.conjugate(), q(1)),
                        &Schur::monomial(mu.clone(), q(1)),
                    );
                    for nu in parts.iter() {
                        assert_eq!(
                            g.coeff(nu),
                            conj.coeff(&nu.conjugate()),
                            "g^{{ν\'}}_{{λ\',μ}} = g^ν_{{λ,μ}}"
                        );
                    }

                    let dim = |p: &Partition| crate::character::character(p, &ones);
                    let total: i128 = parts
                        .iter()
                        .map(|nu| {
                            let c = g.coeff(nu);
                            assert_eq!(c.denom(), 1, "Kronecker coefficients are integers");
                            c.numer() * dim(nu)
                        })
                        .sum();
                    assert_eq!(total, dim(lambda) * dim(mu), "dim of the tensor product");
                }
            }
        }
    }

    #[test]
    fn internal_product_across_degrees_is_zero() {
        // The internal product is defined degree-wise; Sym_m * Sym_n = 0.
        let r = internal(&sch(&[2]), &sch(&[2, 1]));
        assert!(r.is_zero(), "different degrees must annihilate");
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
