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

/// R10 in this file: both sides share `i64`, and the bound is what makes that legal.
///
/// Every law above computes the two sides it compares over the *same* fixed
/// width, which is the shape R10 names — two sides sharing one width wrap
/// identically and agree on the same wrong answer. What makes it legitimate
/// here is R10's own escape clause, an a-priori bound, and until this test the
/// bound was assumed rather than stated.
///
/// Measured: the widest value anywhere in these sweeps is **120**, and it is
/// `z_{1⁵} = 5!` from the Hall pairing — not a structure constant at all.
/// Against `i64::MAX ≈ 9.2·10¹⁸` that is 56 bits of headroom.
///
/// The bound is a property of the **degree caps** (5 for `sweep`, 3 for the
/// products, 2 for the coproduct), not of the laws, so raising a cap spends the
/// headroom silently. This test is what makes that spending visible.
///
/// Two things it does not do. It re-walks the sweeps rather than observing the
/// laws themselves, so a law added above without a line here is not covered —
/// the drift is real and the fix is to extend both together. And it is not the
/// last line of defense: since `overflow-checks` went into the release profile,
/// no profile wraps `i64` silently, so an overflow here would panic rather than
/// produce the agreeing-wrong-answer R10 is about. That backstop covers native
/// integers only — `as` casts, `wrapping_*`, and `Guarded` all still sit
/// outside it, which is why the rule is not retired.
#[test]
fn the_laws_run_far_below_the_width_both_sides_share() {
    /// The measured maximum. A change here is not a failure — it is a request
    /// to re-derive the headroom before accepting the new number.
    const MEASURED_MAX: i64 = 120;

    let mut worst: (i64, String) = (0, String::new());
    let mut see = |v: i64, what: String| {
        if v.abs() > worst.0 {
            worst = (v.abs(), what);
        }
    };

    for lam in sweep() {
        let s = schur(&lam);
        for c in Homogeneous::from_schur(&s).terms().values() {
            see(*c, format!("s→h at {lam}"));
        }
        for c in Elementary::from_schur(&s).terms().values() {
            see(*c, format!("s→e at {lam}"));
        }
        for c in Monomial::from_schur(&s).terms().values() {
            see(*c, format!("s→m at {lam}"));
        }
        see(lam.z() as i64, format!("z at {lam}"));
    }

    let small: Vec<Partition> = (0..=3).flat_map(partitions_of).collect();
    for a in &small {
        for b in &small {
            let ha: Homogeneous<i64> = Homogeneous::monomial(a.clone(), 1);
            let hb: Homogeneous<i64> = Homogeneous::monomial(b.clone(), 1);
            for c in ha.mul(&hb).to_schur().terms().values() {
                see(*c, format!("h·h→s at {a},{b}"));
            }
            let ea: Elementary<i64> = Elementary::monomial(a.clone(), 1);
            let eb: Elementary<i64> = Elementary::monomial(b.clone(), 1);
            for c in ea.mul(&eb).to_schur().terms().values() {
                see(*c, format!("e·e→s at {a},{b}"));
            }
        }
    }

    let tiny: Vec<Partition> = (0..=2).flat_map(partitions_of).collect();
    for a in &tiny {
        for b in &tiny {
            let t = tensor_mul(&coproduct(&schur(a)), &coproduct(&schur(b)));
            for c in t.terms().values() {
                see(*c, format!("Δ⊗Δ at {a},{b}"));
            }
        }
    }

    assert_eq!(
        worst.0,
        MEASURED_MAX,
        "the laws' widest value moved to {} (at {}), from {MEASURED_MAX}.\n  \
         These tests compare two sides that share `i64`, so the bound is what \
         keeps them honest (R10).\n  \
         Re-derive the headroom against `i64::MAX` = {}, then set MEASURED_MAX \
         to the new value.",
        worst.0,
        worst.1,
        i64::MAX,
    );
}
