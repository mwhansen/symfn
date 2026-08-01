//! Standard operations on symmetric functions: the ω involution and the Hall
//! inner product.
//!
//! Both have cheap native forms in a preferred basis (ω on Schur/power sums, the
//! inner product via Schur orthonormality) and a generic form for any basis,
//! obtained by routing through the Schur hub.

// A partition length.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use crate::character::character_in;
use crate::coeff::{QAlgebra, Ring};
use crate::convert::{FromSchur, ToSchur};
#[cfg(feature = "bignum")]
use crate::guard::{guarded, GuardedRat};
use crate::partition::Partition;
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
        // `z_in`, not `from_u128(z())`: the latter caps the internal product at
        // degree 34, since z_{1^35} has no `u128` representation.
        acc.add_term(lambda.clone(), cs.mul(cl).mul(&lambda.z_in::<C>()));
    }
    acc.to_schur()
}

/// A single Kronecker coefficient g^ν_{λμ}.
///
/// Convenience over [`internal`]; computing one costs the same as computing the
/// whole product, since the power-sum route produces every ν at once. When only
/// one ν is wanted at a degree where the whole product does not fit, use
/// [`kronecker_via_characters`] instead.
pub fn kronecker<C: QAlgebra>(
    lambda: &crate::partition::Partition,
    mu: &crate::partition::Partition,
    nu: &crate::partition::Partition,
) -> C {
    let sl: Schur<C> = Schur::monomial(lambda.clone(), C::one());
    let sm: Schur<C> = Schur::monomial(mu.clone(), C::one());
    internal(&sl, &sm).coeff(nu)
}

/// A single Kronecker coefficient g^ν_{λμ}, **without forming the product**.
///
/// The orthogonality formula, which is the definition unwound:
///
/// ```text
///   g^ν_{λμ} = ⟨χ^λ χ^μ, χ^ν⟩ = Σ_{ρ ⊢ n} χ^λ(ρ) χ^μ(ρ) χ^ν(ρ) / z_ρ
/// ```
///
/// This is the same identity `tests/` already checks [`internal`] against — the
/// cross-check exists because the product route rests entirely on the p-basis
/// diagonality plus s ↔ p, so a wrong identity there would be self-consistently
/// wrong. Read the other way it is an algorithm, and a different one: three
/// character *rows* and a weighted dot product.
///
/// **What it drops.** [`internal`] is s → p, a coefficientwise multiply, then
/// p → s back. The first two steps are cheap; the third is
/// [`PowerSum::to_schur`](crate::convert::ToSchur::to_schur), which expands
/// every p_ρ into every λ ⊢ n. That is the p(n) × p(n) work, and the memory
/// ceiling `docs/record/kronecker.md` records (1.1 GB at n = 32). Here the
/// back-transition is replaced by a third character row, so the cost is
/// **3·p(n) Murnaghan–Nakayama evaluations** — heavily shared, since
/// [`try_character`](crate::character::try_character) memoizes and the
/// recursion re-enters itself — and the memory is O(p(n)).
///
/// Past n = 32 the gap widens for a second reason: `p → s` batches its work
/// behind a β-mask of width [`MASK_LIMIT`](crate::convert), and above that falls
/// back to per-character evaluation, losing the sharing. This route never enters
/// that code.
///
/// It is *not* an asymptotic improvement. p(n) grows like exp(c√n), so this is
/// subexponential, not polynomial; computing Kronecker coefficients is #P-hard
/// and nothing here changes that. What it buys is the constant and the memory,
/// which is the difference between "does not finish" and "one coefficient" in
/// the n = 32–50 range. The polynomial-time bounded-parameter algorithms
/// (Christandl–Doran–Walter lattice-point counting; Panova, arXiv:2502.20253)
/// are a different axis — bounded *rows*, unbounded n — and are unimplemented
/// here; this routine is the intended oracle for them.
///
/// # Exactness
///
/// Divides by z_ρ **one small factor at a time** rather than forming z_ρ and
/// dividing once. That is not a micro-optimisation: [`Partition::z`] returns
/// `u128`, and z_{1^n} = n!, which leaves `u128` at n = 35 — precisely the range
/// this routine exists to reach. Every divisor used here is a part of ρ or a
/// multiplicity of one, so all of them are ≤ n.
///
/// The characters go through [`character_in`](crate::character::character_in),
/// so a bignum `C` is exact past the i128 character ceiling at n ≈ 58.
///
/// # Panics
///
/// A fixed-width `C` is the binding constraint, and it binds **much earlier than
/// the characters do: measured, the wall is at n ≈ 26** — pinned by
/// `unguarded_fixed_width_refuses_where_the_guarded_path_escalates`. The reason
/// is the same one `docs/record/kronecker.md` records for the `st` basis — the
/// running sum is a rational whose denominator divides lcm(z_ρ) even though the
/// answer is a small integer, so the *intermediates* leave i128 while the result
/// would fit comfortably. Over [`Rational`](crate::coeff::Rational) that
/// overflow now panics in every profile (`docs/policies/failure.md`, R3); before
/// the release profile carried `overflow-checks` it wrapped, and a wrapped
/// intermediate can land on a denominator of 1 and be accepted as an integer.
/// So **this generic form should not be called over `Rational` at n ≳ 26**: use
/// [`kronecker_coeff`], which runs the guarded ring and escalates.
///
/// [`Partition::z`]: crate::partition::Partition::z
pub fn kronecker_via_characters<C: QAlgebra>(
    lambda: &Partition,
    mu: &Partition,
    nu: &Partition,
) -> C {
    let n = lambda.size();
    // The internal product is defined degree-wise, so unequal degrees pair to
    // zero — the same convention [`internal`] reaches by having no shared λ.
    if mu.size() != n || nu.size() != n {
        return C::zero();
    }
    if n == 0 {
        return C::one();
    }

    let mut acc = C::zero();
    for rho in crate::memo::partitions_cached(n).iter() {
        // Characters vanish often, and each factor tested before the next is
        // computed saves the two Murnaghan-Nakayama sweeps behind it.
        let a: C = character_in(lambda, rho);
        if a.is_zero() {
            continue;
        }
        let b: C = character_in(mu, rho);
        if b.is_zero() {
            continue;
        }
        let c: C = character_in(nu, rho);
        if c.is_zero() {
            continue;
        }

        // term = χ^λ(ρ)·χ^μ(ρ)·χ^ν(ρ) / z_ρ, with z_ρ never formed — see the
        // exactness note above.
        acc.add_assign(&rho.div_by_z(&a.mul(&b).mul(&c)));
    }
    acc
}

/// A single Kronecker coefficient g^ν_{λμ}, exact, with overflow escalation.
///
/// The entry point [`kronecker_via_characters`] should be reached through unless
/// you are supplying your own coefficient ring. It runs the character sum over
/// [`GuardedRat`], which *reports* leaving the fixed width rather than wrapping,
/// and re-runs over `BigRational` when it does.
///
/// # Requires `bignum`
///
/// This function exists only under the `bignum` feature, rather than existing
/// everywhere and panicking without it. Measured, the fixed-width path stops
/// being trustworthy at **n ≈ 26** — far below the n ≈ 58 character ceiling,
/// because the partial sums are rationals over lcm(z_ρ) even though the answer
/// is a small integer. A single-coefficient query is wanted precisely at the
/// degrees where the whole product does not fit, so a build that cannot escalate
/// could serve almost none of its intended range; the honest form of that is an
/// absent function rather than one that panics on most of its inputs.
///
/// The default build keeps its zero dependencies and keeps [`kronecker`], which
/// is the faster route below the crossover over a fixed-width ring anyway.
///
/// ```
/// use symfn::ops::kronecker_coeff;
/// use symfn::partition::Partition;
///
/// // standard ⊗ standard = trivial + standard + sign, in S_3.
/// let std = Partition::new(vec![2, 1]);
/// assert_eq!(kronecker_coeff(&std, &std, &Partition::new(vec![3])), 1);
/// assert_eq!(kronecker_coeff(&std, &std, &std), 1);
/// ```
///
/// # Panics
///
/// If `g^ν_{λμ}` does not fit `i128` — a capacity wall of this signature, not
/// of the computation, which runs over `BigRational`. `symfn.kronecker_coeff`
/// on the Python side returns the same quantity unbounded.
///
/// Also if the exact value is not an integer, which is a bug in this crate
/// rather than an overflow. The two are separate messages on purpose:
/// collapsing them sends the reader after the wrong bug.
#[cfg(feature = "bignum")]
pub fn kronecker_coeff(lambda: &Partition, mu: &Partition, nu: &Partition) -> i128 {
    use num_traits::ToPrimitive;

    // Two ways the fast path can lose, and neither may be fatal: the guard
    // counter moves, or a wrapped intermediate lands on a non-integer. The
    // integrality test must therefore return `None`, not panic — a panic inside
    // the closure escapes `guarded` before it can report the loss.
    let fast = || {
        let v: GuardedRat = kronecker_via_characters(lambda, mu, nu);
        (v.denom() == 1).then(|| v.numer())
    };
    if let Some(Some(v)) = guarded(fast) {
        return v;
    }
    let v: num_rational::BigRational = kronecker_via_characters(lambda, mu, nu);
    // Two distinct failures, and collapsing them into one message sends the
    // reader after the wrong bug: a non-integer means the mathematics is wrong,
    // while a value past i128 means only that this signature is too narrow for
    // it. The Python boundary returns the same quantity unbounded.
    assert!(
        v.is_integer(),
        "g^{nu}_{{{lambda},{mu}}} is not an integer over BigRational, \
         which is a bug rather than an overflow"
    );
    let g = v.to_integer();
    g.to_i128()
        .unwrap_or_else(|| panic!("g^{nu}_{{{lambda},{mu}}} = {g} does not fit i128"))
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

    /// The single-coefficient route against the product route, over every
    /// (λ, μ, ν) triple through degree 7.
    ///
    /// The two share `character_in` and nothing else: `internal` reaches its
    /// answer through s → p, a diagonal multiply and p → s back, and this one
    /// never builds a symmetric function at all. Degree 7 is one past the range
    /// the hand-rolled formula test above covers, so the three checks overlap
    /// rather than merely chain.
    #[test]
    fn kronecker_via_characters_agrees_with_the_product() {
        for n in 1..=7u32 {
            let parts = partitions_cached(n);
            for lambda in parts.iter() {
                for mu in parts.iter() {
                    let prod = internal(
                        &Schur::monomial(lambda.clone(), q(1)),
                        &Schur::monomial(mu.clone(), q(1)),
                    );
                    for nu in parts.iter() {
                        let got: Rational = kronecker_via_characters(lambda, mu, nu);
                        assert_eq!(got, prod.coeff(nu), "g^{nu}_{{{lambda},{mu}}}");
                    }
                }
            }
        }
    }

    /// The grading and the empty case. `Sym_m * Sym_n = 0` for m ≠ n, which the
    /// product route gets from having no shared λ and this one has to state; and
    /// g^∅_{∅∅} = 1, the degree-zero product s_∅ · s_∅ = s_∅.
    #[test]
    fn kronecker_via_characters_respects_the_grading() {
        let z: Rational = kronecker_via_characters(&part(&[2]), &part(&[2, 1]), &part(&[2, 1]));
        assert_eq!(z, q(0));
        let z: Rational = kronecker_via_characters(&part(&[2, 1]), &part(&[2, 1]), &part(&[3, 1]));
        assert_eq!(z, q(0));
        let e: Rational = kronecker_via_characters(&part(&[]), &part(&[]), &part(&[]));
        assert_eq!(e, q(1));
    }

    /// The divisor schedule never forms z_ρ, so the ceiling z_ρ itself has does
    /// not apply to it: z_{1^n} = n!, and 34! is the *last* one inside `u128`
    /// (2.95e38 against a ceiling of 3.40e38; 35! is 1.03e40).
    ///
    /// Pinned rather than asserted by inequality, because the boundary is one
    /// off from the obvious guess — an earlier version of this test claimed 34!
    /// wrapped and was wrong.
    #[test]
    fn factorial_ceiling_for_z_is_at_thirty_five() {
        assert_eq!(
            part(&[1; 34]).z(),
            295_232_799_039_604_140_847_618_609_643_520_000_000
        );
    }

    /// Past n = 34 there is **no product route to compare against** — `internal`
    /// reaches z_μ⁻¹ through `s → p`, which forms z_μ as a `u128` and so is
    /// itself capped by the ceiling pinned above. So the checks here are
    /// identities rather than oracles, which is the point: they hold at degrees
    /// where nothing else in the crate can produce the answer.
    ///
    /// - g^ν_{λ,(n)} = δ_{λν}, tensoring with the trivial character.
    /// - g is symmetric in its three indices.
    ///
    /// Both exercise the divisor schedule at a ρ whose z_ρ is far outside
    /// `u128` — 40! ≈ 8.2e47 — which a routine that formed z_ρ and divided once
    /// could not do at all.
    #[test]
    #[cfg(feature = "bignum")]
    fn kronecker_coeff_holds_identities_past_every_fixed_width_ceiling() {
        let lambda = part(&[38, 2]);
        let trivial = part(&[40]);
        let other = part(&[37, 3]);
        assert_eq!(kronecker_coeff(&lambda, &trivial, &lambda), 1);
        assert_eq!(kronecker_coeff(&lambda, &trivial, &other), 0);

        let (a, b, c) = (part(&[36, 3, 1]), part(&[35, 5]), part(&[37, 2, 1]));
        let g = kronecker_coeff(&a, &b, &c);
        assert_eq!(g, kronecker_coeff(&b, &a, &c));
        assert_eq!(g, kronecker_coeff(&c, &b, &a));
        assert_eq!(g, kronecker_coeff(&a, &c, &b));
    }

    /// The fixed-width path must **report** rather than wrap. Measured, plain
    /// `Rational` at n = 40 returns confident nonsense — a fraction, where the
    /// answer is 1 — while [`kronecker_coeff`] refuses and escalates.
    ///
    /// Pinning the bad value's badness is deliberate: it is the exact failure
    /// the guarded ring exists to remove, and a future "`Rational` is fine
    /// here" would otherwise pass unnoticed.
    ///
    /// This test used to run only under `--ignored`, and to assert the *wrong
    /// answer*: with no `[profile.release]`, `Rational` wrapped in the profile
    /// users ship, and at n = 40 returned a fraction where the answer is 1. The
    /// premise inverted when `overflow-checks = true` landed
    /// (`docs/policies/failure.md`, R3) — the same call now panics in every
    /// profile, so this runs in the ordinary suite and pins the panic instead.
    /// It is a canary for the flag as much as a fact about `Rational`.
    #[test]
    #[cfg(feature = "bignum")]
    #[should_panic(expected = "overflow")]
    fn unguarded_fixed_width_refuses_where_the_guarded_path_escalates() {
        let lambda = part(&[38, 2]);
        let trivial = part(&[40]);
        let _: Rational = kronecker_via_characters(&lambda, &trivial, &lambda);
    }

    /// The other half of the pair above: the same input the unguarded ring
    /// cannot survive, answered exactly by the escalating entry point.
    #[test]
    #[cfg(feature = "bignum")]
    fn guarded_path_escalates_where_fixed_width_refuses() {
        let lambda = part(&[38, 2]);
        let trivial = part(&[40]);
        assert_eq!(kronecker_coeff(&lambda, &trivial, &lambda), 1);
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
