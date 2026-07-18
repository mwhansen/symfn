//! Cross-cutting algebraic-law tests over a sweep of partitions.
//!
//! With no live Sage oracle available offline, these assert *structural* laws
//! that independently pin correctness: conversions are ring homomorphisms,
//! round-trips are the identity, ω is an involutive algebra map, the Hall inner
//! product realizes the known dual pairings, and Δ is an algebra map. The
//! proptest-vs-Sage harness (Phase 5, online) layers on top of these.

use symfn::{
    convert, coproduct, hall, partitions_of, Elementary, FromSchur, Homogeneous, Monomial,
    Partition, PowerSum, Rational, Ring, Schur, SymFn, SymTensor, ToSchur,
};

fn sweep() -> Vec<Partition> {
    (0..=5).flat_map(partitions_of).collect()
}

fn schur(p: &Partition) -> Schur<i64> {
    Schur::monomial(p.clone(), 1)
}

#[test]
fn schur_to_basis_and_back_is_identity() {
    for lam in sweep() {
        let s = schur(&lam);
        let via_h: Schur<i64> = convert(&Homogeneous::from_schur(&s));
        let via_e: Schur<i64> = convert(&Elementary::from_schur(&s));
        let via_m: Schur<i64> = convert(&Monomial::from_schur(&s));
        assert_eq!(via_h, s, "s→h→s at {lam}");
        assert_eq!(via_e, s, "s→e→s at {lam}");
        assert_eq!(via_m, s, "s→m→s at {lam}");
    }
}

#[test]
fn power_sum_round_trip_over_rationals() {
    for lam in sweep() {
        let s: Schur<Rational> = Schur::monomial(lam.clone(), Rational::one());
        let p: PowerSum<Rational> = PowerSum::from_schur(&s);
        assert_eq!(p.to_schur(), s, "s→p→s at {lam}");
    }
}

#[test]
fn to_schur_is_a_ring_homomorphism() {
    // (h_λ · h_μ) → s  ==  (h_λ → s) · (h_μ → s), and likewise for e.
    let small: Vec<Partition> = (0..=3).flat_map(partitions_of).collect();
    for a in &small {
        for b in &small {
            let ha: Homogeneous<i64> = Homogeneous::monomial(a.clone(), 1);
            let hb: Homogeneous<i64> = Homogeneous::monomial(b.clone(), 1);
            assert_eq!(
                ha.mul(&hb).to_schur(),
                ha.to_schur().mul(&hb.to_schur()),
                "h hom at {a},{b}"
            );

            let ea: Elementary<i64> = Elementary::monomial(a.clone(), 1);
            let eb: Elementary<i64> = Elementary::monomial(b.clone(), 1);
            assert_eq!(
                ea.mul(&eb).to_schur(),
                ea.to_schur().mul(&eb.to_schur()),
                "e hom at {a},{b}"
            );
        }
    }
}

#[test]
fn omega_is_an_involutive_algebra_map() {
    let small: Vec<Partition> = (0..=3).flat_map(partitions_of).collect();
    for a in &small {
        // ω² = id
        let s = schur(a);
        assert_eq!(s.omega().omega(), s, "ω² at {a}");
        // ω(s_λ) = s_{λ'}
        assert_eq!(s.omega().coeff(&a.conjugate()), 1);

        for b in &small {
            // ω(f·g) = ω(f)·ω(g)
            let sb = schur(b);
            assert_eq!(
                s.mul(&sb).omega(),
                s.omega().mul(&sb.omega()),
                "ω hom at {a},{b}"
            );
        }
    }
}

#[test]
fn hall_pairings_are_correct() {
    for a in sweep() {
        for b in sweep() {
            if a.size() != b.size() {
                continue;
            }
            let delta = if a == b { 1 } else { 0 };

            // ⟨s_λ, s_μ⟩ = δ
            let sa = schur(&a);
            let sb = schur(&b);
            assert_eq!(hall(&sa, &sb), delta, "⟨s{a},s{b}⟩");

            // ⟨h_λ, m_μ⟩ = δ (dual bases)
            let ha: Homogeneous<i64> = Homogeneous::monomial(a.clone(), 1);
            let mb: Monomial<i64> = Monomial::monomial(b.clone(), 1);
            assert_eq!(hall(&ha, &mb), delta, "⟨h{a},m{b}⟩");

            // ⟨p_λ, p_μ⟩ = z_λ δ
            let pa: PowerSum<i64> = PowerSum::monomial(a.clone(), 1);
            let pb: PowerSum<i64> = PowerSum::monomial(b.clone(), 1);
            let expected = if a == b { a.z() as i64 } else { 0 };
            assert_eq!(hall(&pa, &pb), expected, "⟨p{a},p{b}⟩");
        }
    }
}

/// Multiply two Sym⊗Sym elements: (s_a⊗s_b)(s_c⊗s_d) = (s_a·s_c)⊗(s_b·s_d).
fn tensor_mul(x: &SymTensor<i64>, y: &SymTensor<i64>) -> SymTensor<i64> {
    let mut out = SymTensor::zero();
    for ((a, b), cx) in x.terms() {
        for ((c, d), cy) in y.terms() {
            let left = Schur::monomial(a.clone(), 1).mul(&Schur::monomial(c.clone(), 1));
            let right = Schur::monomial(b.clone(), 1).mul(&Schur::monomial(d.clone(), 1));
            let scale = cx * cy;
            for (lp, lc) in left.terms() {
                for (rp, rc) in right.terms() {
                    out.add_term((lp.clone(), rp.clone()), lc * rc * scale);
                }
            }
        }
    }
    out
}

#[test]
fn coproduct_is_an_algebra_map() {
    // Δ(f·g) = Δ(f)·Δ(g) — the bialgebra compatibility, on small products.
    let small: Vec<Partition> = (0..=2).flat_map(partitions_of).collect();
    for a in &small {
        for b in &small {
            let sa = schur(a);
            let sb = schur(b);
            let lhs = coproduct(&sa.mul(&sb));
            let rhs = tensor_mul(&coproduct(&sa), &coproduct(&sb));
            assert_eq!(lhs, rhs, "Δ algebra map at {a},{b}");
        }
    }
}
