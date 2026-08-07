//! Plethysm `f[g]` of two symmetric functions in the Schur basis.
//!
//! The algorithm runs through the power-sum basis, where plethysm becomes
//! almost trivial:
//!
//! - `p_n[p_m] = p_{nm}`, so `p_n[g]` is just g's p-expansion with every part
//!   multiplied by n (rational coefficients are fixed by `p_n`);
//! - plethysm is additive and *multiplicative* in the outer argument, so
//!   `f[g] = Σ_λ c_λ ∏_i p_{λ_i}[g]` once `f = Σ_λ c_λ p_λ`;
//! - products in the p-basis are multiset unions — the cheapest product in
//!    the crate.
//!
//! So the cost is dominated by the p→s conversion at the end, which runs on
//! memoized Murnaghan–Nakayama characters (`docs/record/plethysm.md`). Requires
//! a [`Plethystic`] ring: `z_μ⁻¹` needs division by an integer, and `p_n` acts
//! on the coefficients as well as the parts.

use crate::coeff::Plethystic;
use crate::convert::{FromSchur, ToSchur};
use crate::partition::Partition;
use crate::sym::{PowerSum, Schur, SymFn};

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
pub fn plethysm<C: Plethystic>(f: &Schur<C>, g: &Schur<C>) -> Schur<C> {
    let pf: PowerSum<C> = PowerSum::from_schur(f);
    let pg: PowerSum<C> = PowerSum::from_schur(g);

    let mut acc = PowerSum::zero();
    for (lambda, c) in pf.terms() {
        // p_λ[g] = ∏_i p_{λ_i}[g]
        let mut prod: PowerSum<C> = PowerSum::monomial(Partition::default(), C::one());
        for &part in lambda.parts() {
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

    #[test]
    fn plethysm_identities() {
        // s_1[g] = g and f[s_1] = f
        let g = s(&[2, 1]);
        assert_eq!(plethysm(&s(&[1]), &g), g);
        assert_eq!(plethysm(&g, &s(&[1])), g);
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
