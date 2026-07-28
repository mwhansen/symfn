//! The (q,t)-Kostka polynomials `K_{λμ}(q,t)`.
//!
//! ```text
//!   J_μ(x; q, t) = Σ_λ K_{λμ}(q,t) S_λ(x; t),      S_λ(x; t) = s_λ[X(1 − t)]
//! ```
//!
//! Macdonald, *Symmetric Functions and Hall Polynomials*, 2nd ed., VI (8.11).
//! `K_{λμ}(0,t)` is the Kostka–Foulkes polynomial `K_{λμ}(t)` and `K_{λμ}(0,1)`
//! the ordinary Kostka number, so this is the two-variable top of the tower
//! [`kostka`](crate::kostka) and [`kf`](crate::kf) already occupy — and both of
//! them are available here as independent checks, which is most of why the test
//! module below is worth its length.
//!
//! Unlike the Kostka and Kostka–Foulkes tables this one is **dense**: `K_{λμ}`
//! is generally nonzero even when λ does not dominate μ. `K_{(21),(3)} = q² + q`
//! is the smallest witness. Nor is the diagonal 1: `K_{(21),(21)} = 1 + qt`.
//! Both the triangularity and the unit diagonal are `q = 0` phenomena, not
//! properties of these — which is worth knowing before writing a test that
//! assumes the Kostka–Foulkes shape carries over. Two of the ones below did.
//!
//! ## Inverting the S basis
//!
//! Reading `K` off means expanding `J_μ` in `{S_λ}` rather than in `{s_λ}`, so
//! what is needed is the inverse of `f ↦ f[X(1−t)]`. Call it
//!
//! ```text
//!   φ_t : p_n ↦ p_n / (1 − t^n),   extended ℚ(q,t)-linearly.
//! ```
//!
//! `S_λ = s_λ[X(1−t)]` has `p_ν`-coefficient `z_ν⁻¹ χ^λ_ν ∏_i (1 − t^{ν_i})`, so
//! `φ_t(S_λ) = s_λ` and therefore `φ_t(J_μ) = Σ_λ K_{λμ}(q,t) s_λ`. The whole
//! computation is: expand `J_μ` in power sums, divide the `p_ν` coefficient by
//! `∏_i (1 − t^{ν_i})`, expand back into Schur.
//!
//! ## It is not a plethysm, and that matters
//!
//! The plan this crate was carrying said the step needed a
//! [`Plethystic`](crate::coeff::Plethystic) impl for [`Frac`]. It does not, and
//! writing one and using it here would give the wrong answer.
//!
//! `φ_t` is **ℚ(q,t)-linear**: it acts on the alphabet `X` and holds the
//! coefficients fixed. A plethysm does the opposite — [`Plethystic::frobenius`]
//! exists precisely because `p_n` substitutes into the coefficient ring too, so
//! `p_n[t·p_1] = t^n·p_n`. Routing `J_μ` through the genuine plethysm would
//! raise the `q` and `t` already sitting in `J`'s coefficients, and `K` would
//! come out in the wrong variables.
//!
//! The two are easy to conflate because both are written `f[X/(1−t)]`. The
//! `(1 − t^n)` appearing in the divisor *is* the Frobenius image of `1 − t` —
//! but it comes from `S_λ`'s definition, where the coefficients are rational
//! constants and there is nothing else to raise. Nothing in `J` gets raised.
//! [`scale_parts`](crate::plethysm) is the operation that would have been
//! wrong here; the loop in [`invert_s_basis`] is the one that is right.
//!
//! ## Why the coefficient ring is a ℚ-algebra
//!
//! `s → p` divides by `z_ν` and `φ_t` divides by `1 − t^n`, so the route runs
//! through ℚ(q,t) even though `K_{λμ}(q,t) ∈ ℤ[q,t]` at the end (Haiman's
//! theorem gives non-negative integers, which the tests check but do not rely
//! on). The denominators must cancel exactly, and
//! [`Frac::into_poly`] is where that is enforced: a leftover factor is a bug,
//! not a fallback.

use std::collections::BTreeMap;

use crate::coeff::{QAlgebra, Ring};
use crate::convert::{FromSchur, ToSchur};
use crate::frac::Frac;
use crate::partition::Partition;
use crate::qt::QtPoly;
use crate::sym::{Monomial, PowerSum, Schur, SymFn};

/// `K_{λμ}(q,t)`.
///
/// This computes the whole of `J_μ` and reads one coefficient out of it, exactly
/// as [`kostka_foulkes`](crate::kf::kostka_foulkes) does. For more than one λ at
/// a fixed μ use [`qt_kostka_column`], and for a whole degree
/// [`qt_kostka_table`].
pub fn qt_kostka<C: QAlgebra>(lambda: &Partition, mu: &Partition) -> QtPoly<C> {
    if lambda.size() != mu.size() {
        return QtPoly::zero();
    }
    column_expansion::<C>(mu).coeff(lambda)
}

/// Every `K_{λμ}(q,t)` for a fixed μ, paired with its λ.
///
/// The natural unit of work: one `J_μ` is the whole column. Every λ of the
/// degree appears — see the note on density in the module docs.
pub fn qt_kostka_column<C: QAlgebra>(mu: &Partition) -> Vec<(Partition, QtPoly<C>)> {
    column_expansion::<C>(mu)
        .terms()
        .iter()
        .map(|(lambda, k)| (lambda.clone(), k.clone()))
        .collect()
}

/// The whole table of degree `n`, as `table[i][j] = K_{λⁱ λʲ}(q,t)` indexed
/// against [`partitions_cached`](crate::memo::partitions_cached) — the same
/// orientation as [`kostka_table`](crate::kostka::kostka_table) and
/// [`kostka_foulkes_table`](crate::kf::kostka_foulkes_table).
///
/// Each column costs its own `J_μ`; unlike the Hall–Littlewood recursion there
/// is nothing shared between columns to exploit, so this is p(n) independent
/// computations and not a cheaper joint one.
pub fn qt_kostka_table<C: QAlgebra>(n: u32) -> Vec<Vec<QtPoly<C>>> {
    let parts = crate::memo::partitions_cached(n);
    let mut table = vec![vec![QtPoly::zero(); parts.len()]; parts.len()];
    for (j, mu) in parts.iter().enumerate() {
        for (lambda, k) in column_expansion::<C>(mu).terms() {
            let i = parts
                .iter()
                .position(|p| p == lambda)
                .expect("K_{λμ} is supported on partitions of |μ|");
            table[i][j] = k.clone();
        }
    }
    table
}

/// `H̃_μ(x; q, t) = Σ_λ K̃_{λμ}(q,t) s_λ`, the modified Macdonald polynomial in
/// the Schur basis.
///
/// `K̃_{λμ}(q,t) = t^{n(μ)} K_{λμ}(q, 1/t)` with `n(μ) = Σ (i−1)μ_i`, so this is
/// [`qt_kostka_column`] with the `t`-exponents reflected — bookkeeping, not
/// arithmetic. The whole computation is above; this is the form the modern
/// literature states, and the one in which Haiman's positivity theorem reads
/// "non-negative integers" without a normalising power in the way.
///
/// Reflection is exact: `K_{λμ}` has `t`-degree at most `n(μ)`, which is what
/// makes `K̃` a polynomial at all, and the loop asserts it rather than assuming.
///
/// `H̃_{(2)} = s_2 + q·s_{11}` and `H̃_{(11)} = s_2 + t·s_{11}` are the smallest
/// pair, and show the `q ↔ t` symmetry under conjugating μ that the twisted form
/// has and `K` does not.
pub fn macdonald_ht<C: QAlgebra>(mu: &Partition) -> Schur<QtPoly<C>> {
    let n_mu: u32 = mu
        .parts()
        .iter()
        .enumerate()
        .map(|(i, &p)| i as u32 * p)
        .sum();
    let mut out = Schur::zero();
    for (lambda, k) in column_expansion::<C>(mu).terms() {
        let mut flipped = QtPoly::zero();
        for (&(a, b), c) in k.terms() {
            assert!(
                b <= n_mu,
                "K_{{{lambda},{mu}}} has t-degree {b} > n(mu) = {n_mu}"
            );
            flipped.add_term(a, n_mu - b, c.clone());
        }
        out.add_term(lambda.clone(), flipped);
    }
    out
}

/// `K̃_{λμ}(q,t)`, the modified (q,t)-Kostka polynomial — one coefficient of
/// [`macdonald_ht`], which is what it computes.
pub fn modified_qt_kostka<C: QAlgebra>(lambda: &Partition, mu: &Partition) -> QtPoly<C> {
    if lambda.size() != mu.size() {
        return QtPoly::zero();
    }
    macdonald_ht::<C>(mu).coeff(lambda)
}

/// `φ_t(J_μ)` in the Schur basis — the column, before it is taken apart.
fn column_expansion<C: QAlgebra>(mu: &Partition) -> Schur<QtPoly<C>> {
    let j: Monomial<Frac<C>> = crate::macdonald_j(mu);
    let p: PowerSum<Frac<C>> = PowerSum::from_schur(&j.to_schur());
    let expanded: Schur<Frac<C>> = invert_s_basis(&p).to_schur();

    let mut out = Schur::zero();
    for (lambda, c) in expanded.terms() {
        let k = c.clone().into_poly().unwrap_or_else(|| {
            panic!("K_{{{lambda},{mu}}} is not a polynomial: {c}");
        });
        out.add_term(lambda.clone(), k);
    }
    out
}

/// `φ_t`: divide the coefficient of `p_ν` by `∏_i (1 − t^{ν_i})`.
///
/// Applied through [`Frac::mul_factors`], so the divisor is never expanded —
/// `(0, ν_i) ↦ −1` is a denominator factor in the encoding [`Frac`] already
/// speaks, and factors shared with the numerator cancel without a product being
/// formed. A partition has no zero parts, so `(0,0)` cannot arise.
///
/// **Reduced here**, which is where this route's whole cost turned out to live.
/// `PowerSum::to_schur` uses each of these p(n) times — once per λ — and every
/// use is an `add_assign` that lifts the running sum to the lcm of the two
/// denominators. Cutting a denominator down once, here, is p(n) lifts it does
/// not have to widen. Doing it costs 0.04s at degree 9 and takes `p → s` from
/// 0.88s to 0.10s: the phase split went 1.39s → 0.66s on that line alone.
///
/// It is the same shape as the reduce in [`Ring::mul`](Frac::mul) — reduce once,
/// before a value is used many times — and it makes that one redundant for this
/// path, since a coefficient now arrives already reduced.
fn invert_s_basis<C: QAlgebra>(p: &PowerSum<Frac<C>>) -> PowerSum<Frac<C>> {
    let mut out = PowerSum::zero();
    for (nu, c) in p.terms() {
        let mut factors: BTreeMap<(u32, u32), i32> = BTreeMap::new();
        for &part in nu.parts() {
            *factors.entry((0, part)).or_insert(0) -= 1;
        }
        let mut v = c.mul_factors(&factors);
        v.reduce();
        out.add_term(nu.clone(), v);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::{Rational, Ring};

    type Q = QtPoly<Rational>;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn k(lambda: &[u32], mu: &[u32]) -> Q {
        qt_kostka(&part(lambda), &part(mu))
    }

    /// Hand-computed in the module's own terms, and the case that pins the
    /// orientation: `K_{(11),(2)} = q` while `K_{(2),(11)} = t`. Swapping the
    /// two indices survives every symmetric check further down.
    ///
    /// `J_(11) = (1−t)(1−t²) e_2`, and `φ_t(e_2) = (t·s_2 + s_11)/((1−t)(1−t²))`,
    /// so `K_{(2),(11)} = t` and `K_{(11),(11)} = 1`.
    #[test]
    fn degree_two_is_the_hand_computation() {
        assert_eq!(k(&[2], &[2]), <Q as Ring>::one());
        assert_eq!(k(&[1, 1], &[2]), QtPoly::term(1, 0, Rational::from_int(1)));
        assert_eq!(k(&[2], &[1, 1]), QtPoly::term(0, 1, Rational::from_int(1)));
        assert_eq!(k(&[1, 1], &[1, 1]), <Q as Ring>::one());
    }

    /// `K_{λμ}(0,t)` is the Kostka–Foulkes polynomial.
    ///
    /// The strongest check available in-crate, and cheap to state: this module
    /// runs `J_μ` through the power-sum basis over ℚ(q,t), while
    /// [`kostka_foulkes`](crate::kf::kostka_foulkes) reads its answer off the
    /// Morris recursion in ℤ[t]. No shared code below `Partition`.
    #[test]
    fn at_q_zero_it_is_kostka_foulkes() {
        for n in 1..=6u32 {
            for mu in crate::partitions_of(n) {
                for (lambda, kqt) in qt_kostka_column::<Rational>(&mu) {
                    let want: QtPoly<Rational> = crate::kostka_foulkes(&lambda, &mu);
                    let got = set_q_to_zero(&kqt);
                    assert_eq!(got, want, "K_{{{lambda},{mu}}} at q = 0");
                }
            }
        }
    }

    /// `K_{λμ}(q,t) = K_{λ'μ'}(t,q)` — Macdonald VI (8.15).
    ///
    /// Independent of the `q = 0` check in a way that matters: that one only
    /// ever looks at the `q⁰` slice, so it says nothing about how the two
    /// variables are stored relative to each other. This one is entirely about
    /// that.
    #[test]
    fn conjugating_both_shapes_swaps_q_and_t() {
        for n in 1..=6u32 {
            for mu in crate::partitions_of(n) {
                let conj = qt_kostka_column::<Rational>(&mu.conjugate());
                for (lambda, kqt) in qt_kostka_column::<Rational>(&mu) {
                    let other = conj
                        .iter()
                        .find(|(l, _)| l == &lambda.conjugate())
                        .map(|(_, k)| swap_variables(k))
                        .unwrap_or_else(QtPoly::zero);
                    assert_eq!(kqt, other, "K_{{{lambda},{mu}}} under conjugation");
                }
            }
        }
    }

    /// Haiman's theorem: the coefficients are non-negative integers.
    ///
    /// Not something this code arranges — it arrives over ℚ(q,t) and every
    /// denominator has to cancel — so it is a real check on the whole route and
    /// not a restatement of it. `into_poly` already refuses a surviving
    /// denominator; this catches a surviving `1/2`.
    #[test]
    fn coefficients_are_non_negative_integers() {
        for n in 1..=6u32 {
            for mu in crate::partitions_of(n) {
                for (lambda, kqt) in qt_kostka_column::<Rational>(&mu) {
                    for (_, c) in kqt.terms() {
                        assert_eq!(c.denom(), 1, "K_{{{lambda},{mu}}} has {c:?} in it");
                        assert!(c.numer() > 0, "K_{{{lambda},{mu}}} has {c:?} in it");
                    }
                }
            }
        }
    }

    /// `K_{λμ}(1,1) = f^λ`, the number of standard tableaux of shape λ — the
    /// **same value for every μ**.
    ///
    /// The Garsia–Haiman module `R_μ` is a graded version of the regular
    /// representation of `S_n`, and these are its graded multiplicities, so
    /// specialising the grading away leaves `dim S^λ` no matter which μ was
    /// chosen. That independence is what makes this a sharp test: `q = 0` only
    /// ever inspects one slice of each polynomial, and this reads all of it.
    ///
    /// It is also where the density claim in the module docs is checked, since
    /// `f^λ > 0` always. `K_{(21),(3)} = q² + q` is nonzero even though (21)
    /// does not dominate (3) — an implementation that imposed the dominance
    /// triangularity of the *ordinary* Kostka numbers would pass the `q = 0`
    /// test and fail here.
    #[test]
    fn at_q_and_t_one_it_counts_standard_tableaux() {
        for n in 1..=6u32 {
            for mu in crate::partitions_of(n) {
                let column = qt_kostka_column::<Rational>(&mu);
                assert_eq!(
                    column.len(),
                    crate::partitions_of(n).len(),
                    "the column of {mu} should be dense"
                );
                for (lambda, kqt) in column {
                    let got: i128 = kqt.terms().map(|(_, c)| c.numer()).sum();
                    let want = crate::dimension(&lambda).unwrap() as i128;
                    assert_eq!(got, want, "K_{{{lambda},{mu}}} at q = t = 1");
                }
            }
        }
    }

    /// The table must be the columns, in the documented orientation.
    #[test]
    fn the_table_agrees_with_the_columns() {
        for n in 1..=5u32 {
            let parts = crate::partitions_of(n);
            let table = qt_kostka_table::<Rational>(n);
            for (j, mu) in parts.iter().enumerate() {
                for (i, lambda) in parts.iter().enumerate() {
                    assert_eq!(table[i][j], qt_kostka(lambda, mu), "[{i}][{j}]");
                }
            }
        }
        // ...and the orientation check only has teeth on an asymmetric table.
        let t = qt_kostka_table::<Rational>(3);
        assert!(
            (0..t.len()).any(|i| (0..t.len()).any(|j| t[i][j] != t[j][i])),
            "a symmetric table would make the orientation untestable"
        );
    }

    /// `H̃_μ` is symmetric under conjugating μ together with swapping q and t —
    /// `K̃_{λμ}(q,t) = K̃_{λ'μ'}(t,q)`... with λ **not** conjugated, unlike the
    /// relation `K` satisfies.
    ///
    /// That difference is the whole reason to have this form, and it is what a
    /// wrong `n(μ)` would break: the reflection is by a shape-dependent power,
    /// so getting `n(μ) = Σ(i−1)μ_i` confused with `n(μ') = Σ binom(μ_i, 2)`
    /// still yields polynomials and still passes an integrality check.
    #[test]
    fn the_modified_form_is_symmetric_in_q_and_t() {
        for n in 1..=6u32 {
            for mu in crate::partitions_of(n) {
                let conj = macdonald_ht::<Rational>(&mu.conjugate());
                for (lambda, kt) in macdonald_ht::<Rational>(&mu).terms() {
                    assert_eq!(
                        kt,
                        &swap_variables(&conj.coeff(lambda)),
                        "K~_{{{lambda},{mu}}}"
                    );
                }
            }
        }
    }

    /// The two smallest, which pin the direction of the reflection: an inverted
    /// `n(μ)` sends `H̃_{(11)}` to `s_2 + t·s_{11}` as well, and the pair would
    /// still look symmetric.
    #[test]
    fn h_tilde_of_degree_two_is_the_known_pair() {
        let one = Rational::from_int(1);
        let ht = macdonald_ht::<Rational>(&part(&[2]));
        assert_eq!(ht.coeff(&part(&[2])), <Q as Ring>::one());
        assert_eq!(ht.coeff(&part(&[1, 1])), QtPoly::term(1, 0, one));

        let ht = macdonald_ht::<Rational>(&part(&[1, 1]));
        assert_eq!(ht.coeff(&part(&[2])), <Q as Ring>::one());
        assert_eq!(ht.coeff(&part(&[1, 1])), QtPoly::term(0, 1, one));
    }

    fn set_q_to_zero(p: &QtPoly<Rational>) -> QtPoly<Rational> {
        let mut out = QtPoly::zero();
        for (&(a, b), c) in p.terms() {
            if a == 0 {
                out.add_term(0, b, c.clone());
            }
        }
        out
    }

    fn swap_variables(p: &QtPoly<Rational>) -> QtPoly<Rational> {
        let mut out = QtPoly::zero();
        for (&(a, b), c) in p.terms() {
            out.add_term(b, a, c.clone());
        }
        out
    }
}
