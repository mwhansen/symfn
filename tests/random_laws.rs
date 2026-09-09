//! The algebraic laws of `tests/algebra_laws.rs`, over inputs nobody chose.
//!
//! That file sweeps every partition of degree at most 5, and its product laws
//! stop at degree 3: a fixed list, so a defect that needs a shape outside it
//! is never asked about. Here the shapes come from a deterministic generator
//! instead, at degrees the enumeration cannot afford, and each law is asked
//! about a thousand of them. The reasoning is the one `check_backend.py`
//! confirmed on the Sage side (`docs/record/python-and-sage-interop.md`): an
//! input the author did not pick finds what the ladder never generated.
//!
//! No `proptest`. The default build has no dependencies and `cargo test` on
//! the published tarball is meant to run offline, so the generator is a
//! xorshift with a fixed seed, the way
//! `examples/bench_schubert_wall.rs` draws its permutations. What that gives
//! up is shrinking; what it keeps is that every assertion names the input that
//! failed and the seed that produced it, so a failure reproduces by rerunning.
//!
//! `SYMFN_LAWS_SEED` and `SYMFN_LAWS_CASES` override the seed and the number
//! of cases per law, for a longer run by hand.
//!
//! Both sides of every comparison share `i128`. That is the shape
//! `docs/policies/failure.md` R10 warns about, and the escape clause is an
//! a-priori bound: `widest_value_stays_far_below_the_width` measures the
//! largest coefficient any law here produces and holds it under `i64::MAX`,
//! which leaves 64 bits between the values and the width they are computed in.

use symfn::{
    convert, coproduct, hall, Elementary, FromSchur, Homogeneous, Monomial, Partition, PowerSum,
    Rational, Ring, Schur, SymFn, SymTensor, ToSchur,
};

/// xorshift64*, seeded once per test so each law sees the same shapes on
/// every run.
struct Rng(u64);

impl Rng {
    fn new(law: &str) -> Rng {
        let seed: u64 = std::env::var("SYMFN_LAWS_SEED")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(0x9e37_79b9_7f4a_7c15);
        // Mix the law's name in, so two laws do not walk the same sequence.
        let mut h = seed;
        for b in law.bytes() {
            h = (h ^ u64::from(b)).wrapping_mul(0x100_0000_01b3);
        }
        Rng(h | 1)
    }

    fn next(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.0 = x;
        x.wrapping_mul(0x2545_f491_4f6c_dd1d)
    }

    fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }

    /// A partition of `n`, drawn by choosing parts from the largest down:
    /// each part is uniform between 1 and the previous part, capped by what
    /// is left. Not uniform over partitions of `n`, and not meant to be — it
    /// reaches long thin shapes and short wide ones in the same run, which is
    /// what a law that fails on one shape family needs.
    fn partition_of(&mut self, n: u32) -> Partition {
        let mut parts = Vec::new();
        let mut left = n;
        let mut cap = n;
        while left > 0 {
            let part = 1 + self.below(u64::from(cap.min(left))) as u32;
            parts.push(part);
            left -= part;
            cap = part;
        }
        Partition::new(parts)
    }

    fn partition_upto(&mut self, max: u32) -> Partition {
        let n = self.below(u64::from(max) + 1) as u32;
        self.partition_of(n)
    }
}

fn cases() -> usize {
    std::env::var("SYMFN_LAWS_CASES")
        .ok()
        .and_then(|s| s.parse().ok())
        .unwrap_or(1000)
}

fn schur(p: &Partition) -> Schur<i128> {
    Schur::monomial(p.clone(), 1)
}

#[test]
fn schur_to_basis_and_back_is_identity_at_random_degrees_up_to_10() {
    let mut rng = Rng::new("round trips");
    for _ in 0..cases() {
        let lam = rng.partition_upto(10);
        let s = schur(&lam);
        let via_h: Schur<i128> = convert(&Homogeneous::from_schur(&s));
        let via_e: Schur<i128> = convert(&Elementary::from_schur(&s));
        let via_m: Schur<i128> = convert(&Monomial::from_schur(&s));
        assert_eq!(via_h, s, "s→h→s at {lam}");
        assert_eq!(via_e, s, "s→e→s at {lam}");
        assert_eq!(via_m, s, "s→m→s at {lam}");
    }
}

#[test]
fn power_sum_round_trip_over_rationals_at_random_degrees_up_to_9() {
    let mut rng = Rng::new("power sums");
    for _ in 0..cases() {
        let lam = rng.partition_upto(9);
        let s: Schur<Rational> = Schur::monomial(lam.clone(), Rational::one());
        let p: PowerSum<Rational> = PowerSum::from_schur(&s);
        assert_eq!(p.to_schur(), s, "s→p→s at {lam}");
    }
}

#[test]
fn to_schur_is_a_ring_homomorphism_on_random_pairs_up_to_degree_7_each() {
    let mut rng = Rng::new("homomorphism");
    for _ in 0..cases() {
        let a = rng.partition_upto(7);
        let b = rng.partition_upto(7);
        let ha: Homogeneous<i128> = Homogeneous::monomial(a.clone(), 1);
        let hb: Homogeneous<i128> = Homogeneous::monomial(b.clone(), 1);
        assert_eq!(
            (&ha * &hb).to_schur(),
            ha.to_schur() * hb.to_schur(),
            "h hom at {a},{b}"
        );
        let ea: Elementary<i128> = Elementary::monomial(a.clone(), 1);
        let eb: Elementary<i128> = Elementary::monomial(b.clone(), 1);
        assert_eq!(
            (&ea * &eb).to_schur(),
            ea.to_schur() * eb.to_schur(),
            "e hom at {a},{b}"
        );
    }
}

#[test]
fn omega_is_an_involutive_algebra_map_on_random_pairs_up_to_degree_7_each() {
    let mut rng = Rng::new("omega");
    for _ in 0..cases() {
        let a = rng.partition_upto(7);
        let b = rng.partition_upto(7);
        let sa = schur(&a);
        let sb = schur(&b);
        assert_eq!(sa.omega().omega(), sa, "ω² at {a}");
        assert_eq!(sa.omega().coeff(&a.conjugate()), 1, "ω(s_λ) = s_λ' at {a}");
        assert_eq!(
            (&sa * &sb).omega(),
            sa.omega() * sb.omega(),
            "ω hom at {a},{b}"
        );
    }
}

#[test]
fn hall_pairings_are_correct_on_random_pairs_up_to_degree_10() {
    let mut rng = Rng::new("hall");
    for _ in 0..cases() {
        let n = rng.below(11) as u32;
        // Same degree by construction, and equal half the time, so the δ = 1
        // branch is asked about as often as the δ = 0 one.
        let a = rng.partition_of(n);
        let b = if rng.below(2) == 0 {
            a.clone()
        } else {
            rng.partition_of(n)
        };
        let delta: i128 = if a == b { 1 } else { 0 };

        assert_eq!(hall(&schur(&a), &schur(&b)), delta, "⟨s{a},s{b}⟩");

        let ha: Homogeneous<i128> = Homogeneous::monomial(a.clone(), 1);
        let mb: Monomial<i128> = Monomial::monomial(b.clone(), 1);
        assert_eq!(hall(&ha, &mb), delta, "⟨h{a},m{b}⟩");

        let pa: PowerSum<i128> = PowerSum::monomial(a.clone(), 1);
        let pb: PowerSum<i128> = PowerSum::monomial(b.clone(), 1);
        let expected = if a == b {
            i128::try_from(a.z()).expect("z fits")
        } else {
            0
        };
        assert_eq!(hall(&pa, &pb), expected, "⟨p{a},p{b}⟩");
    }
}

/// `(s_a⊗s_b)(s_c⊗s_d) = (s_a·s_c)⊗(s_b·s_d)`, the product on `Sym⊗Sym`.
fn tensor_mul(x: &SymTensor<i128>, y: &SymTensor<i128>) -> SymTensor<i128> {
    let mut out = SymTensor::zero();
    for ((a, b), cx) in x.terms() {
        for ((c, d), cy) in y.terms() {
            let left = schur(a) * schur(c);
            let right = schur(b) * schur(d);
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
fn coproduct_is_an_algebra_map_on_random_pairs_up_to_degree_4_each() {
    // The tensor product above is quadratic in the number of terms on each
    // side, which is why this law's degrees stay lowest.
    let mut rng = Rng::new("coproduct");
    for _ in 0..cases() / 4 {
        let a = rng.partition_upto(4);
        let b = rng.partition_upto(4);
        let lhs = coproduct(&(schur(&a) * schur(&b)));
        let rhs = tensor_mul(&coproduct(&schur(&a)), &coproduct(&schur(&b)));
        assert_eq!(lhs, rhs, "Δ algebra map at {a},{b}");
    }
}

#[test]
fn widest_value_stays_far_below_the_width() {
    // Re-walks the widest sweeps above rather than observing the laws, as
    // `tests/algebra_laws.rs` does for its own; a law added above without a
    // line here is not bounded by this test.
    let mut worst: (i128, String) = (0, String::new());
    let mut see = |v: i128, what: String| {
        if v.abs() > worst.0 {
            worst = (v.abs(), what);
        }
    };
    let mut rng = Rng::new("width");
    for _ in 0..cases() {
        let lam = rng.partition_upto(10);
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
        see(
            i128::try_from(lam.z()).expect("z fits"),
            format!("z at {lam}"),
        );

        let a = rng.partition_upto(7);
        let b = rng.partition_upto(7);
        let ha: Homogeneous<i128> = Homogeneous::monomial(a.clone(), 1);
        let hb: Homogeneous<i128> = Homogeneous::monomial(b.clone(), 1);
        for c in (&ha * &hb).to_schur().terms().values() {
            see(*c, format!("h·h→s at {a},{b}"));
        }
        let ea: Elementary<i128> = Elementary::monomial(a.clone(), 1);
        let eb: Elementary<i128> = Elementary::monomial(b.clone(), 1);
        for c in (&ea * &eb).to_schur().terms().values() {
            see(*c, format!("e·e→s at {a},{b}"));
        }
    }
    assert!(
        worst.0 < i128::from(i64::MAX),
        "the widest value reached {} (at {}), within 64 bits of i128's width; \
         lower the degree caps or re-derive the headroom (failure.md, R10)",
        worst.0,
        worst.1
    );
}
