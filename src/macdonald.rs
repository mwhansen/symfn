//! Macdonald's `P_λ(x; q, t)`, by the branching formula.
//!
//! ```text
//!   P_λ(x; q, t) = Σ_T ψ_T(q, t) x^T
//! ```
//!
//! over semistandard tableaux `T` of shape λ. A tableau of content μ is the
//! same thing as a chain `∅ = ν⁰ ⊂ ν¹ ⊂ … ⊂ ν^k = λ` of horizontal strips with
//! `|νⁱ/νⁱ⁻¹| = μ_i`, and `ψ_T` is the product of a factor per strip — so the
//! coefficient of `m_μ` is a sum over exactly the chains
//! [`charge`](mod@crate::charge) already enumerates for Kostka–Foulkes. The two
//! share `build`/`strips` for that reason.
//!
//! Macdonald, *Symmetric Functions and Hall Polynomials*, 2nd ed., Chapter VI,
//! (6.24) and (7.13'). Symmetrica has no Macdonald polynomials at all, so
//! unlike Hall–Littlewood there is no C implementation to compare against —
//! Sage is the only external oracle here.
//!
//! ## Range
//!
//! Over a fixed-width `C` this family refuses rather than wrapping past its
//! wall (`docs/policies/failure.md`, R3). Through the Python boundary all
//! three of `P`, `Q` and `J` escalate — the fixed-width pass reports and the
//! same generic code re-runs over `BigInt` — so the walls below are what a
//! *Rust* caller at `C = i128` meets.
//!
//! **The extremal shape is the single row `λ = (n)`**, at every degree
//! measured, and the walls there are reachable in about a minute per call:
//! `P` gives out at n = 30 (122 bits at n = 29), `Q` and `J` at n = 26. The
//! jump from 101 bits to overflow in one degree says where the wall is — in
//! the **intermediates** of the `Frac` arithmetic rather than in the answers,
//! which is the shape this crate keeps meeting.
//!
//! Shape dominates degree here as it does for Hall–Littlewood, but in the
//! opposite direction: at λ = 1ⁿ the same `J` gains ~0.3 bits per degree and
//! is still 23 bits at n = 96, so 1ⁿ is nowhere near extremal for `J` even
//! though it is exactly extremal for `Q'`.
//!
//! Denominators carry no such number: `Frac` keeps them factored as a multiset
//! of binomials `1 − qᵃtᵇ` and never expands one.
//!
//! Measurements and the harness are in
//! `docs/record/failure-and-overflow.md` (`examples/probe_qt_walls.rs`).
//!
//! ## The coefficients
//!
//! For a cell `s` of λ with arm `a` and leg `l`,
//!
//! ```text
//!   b_λ(s) = (1 − q^a t^{l+1}) / (1 − q^{a+1} t^l)
//! ```
//!
//! and for a horizontal strip λ/μ,
//!
//! ```text
//!   ψ_{λ/μ} = ∏ b_μ(s) / b_λ(s)
//! ```
//!
//! over the cells `s` of μ lying in a row that meets λ/μ but **not** in a
//! column that meets it. Both conditions matter and they are easy to swap; the
//! tests pin them on λ = (2), where `ψ` is `(1−t)(1+q)/(1−qt)` and getting
//! either condition wrong gives something else.
//!
//! Every factor is a ratio of binomials `1 − qᵃtᵇ`, which is the entire reason
//! [`Frac`] can avoid a general gcd — see its module docs. Neither exponent
//! pair can be `(0,0)`: the numerator of `b` has `t`-exponent `l+1 ≥ 1` and the
//! denominator has `q`-exponent `a+1 ≥ 1`.

// Shape indices; coefficients are `Frac<C>`.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::{BTreeMap, HashMap};

use crate::coeff::Ring;
use crate::frac::Frac;
use crate::partition::Partition;
use crate::sym::{Monomial, SymFn};

/// `P_λ(x; q, t)` in the monomial basis.
///
/// Monic and triangular: the coefficient of `m_λ` is 1 and every other `m_μ`
/// that appears has μ strictly below λ in dominance order.
pub fn macdonald_p<C: Ring>(lambda: &Partition) -> Monomial<Frac<C>> {
    let mut out = Monomial::zero();
    if lambda.is_empty() {
        out.add_term(Partition::new([]), <Frac<C> as Ring>::one());
        return out;
    }
    let target: Vec<u32> = lambda.parts().to_vec();
    // ψ of a strip depends only on the two shapes, and the same strip recurs
    // across chains and across contents — a tall shape has few distinct strips
    // and very many tableaux using them.
    let mut cache: HashMap<(Vec<u32>, Vec<u32>), Factors> = HashMap::new();
    let mut acc: Factors = BTreeMap::new();

    for mu in crate::partitions_of(lambda.size()) {
        let mut total = <Frac<C> as Ring>::zero();
        let mut chain: Vec<Vec<u32>> = vec![vec![0; target.len()]];
        crate::charge::build(&target, mu.parts(), 0, &mut chain, &mut |ch| {
            // Accumulate the whole tableau's exponents before building anything.
            acc.clear();
            for k in 1..ch.len() {
                let key = (ch[k].clone(), ch[k - 1].clone());
                let f = cache
                    .entry(key)
                    .or_insert_with(|| psi_factors(&ch[k], &ch[k - 1]));
                for (&e, &m) in f.iter() {
                    *acc.entry(e).or_insert(0) += m;
                }
            }
            acc.retain(|_, m| *m != 0);
            total.add_assign(&Frac::from_factors(&acc));
        });
        if !total.is_zero() {
            // The one reduction for this coefficient — see `Frac::add_assign`.
            total.reduce();
            out.add_term(mu, total);
        }
    }
    out
}

/// `Q_λ(x; q, t) = b_λ(q,t) · P_λ`, the other normalization of the same basis.
///
/// `b_λ = ∏_{s∈λ} b_λ(s)`, so this is `P` scaled by a single ratio of binomial
/// products — no re-enumeration.
pub fn macdonald_q<C: Ring>(lambda: &Partition) -> Monomial<Frac<C>> {
    scale(macdonald_p(lambda), &b_factors(lambda.parts()))
}

/// `J_λ(x; q, t) = c_λ(q,t) · P_λ`, the integral form.
///
/// `c_λ = ∏_{s∈λ} (1 − q^{a(s)} t^{l(s)+1})` — a *polynomial*, and exactly the
/// numerator of `b_λ`, which is what clears `P`'s denominators and makes `J`
/// the form with coefficients in `ℤ[q,t]`.
pub fn macdonald_j<C: Ring>(lambda: &Partition) -> Monomial<Frac<C>> {
    scale(macdonald_p(lambda), &c_factors(lambda.parts()))
}

/// Multiply every coefficient by a product of binomial powers.
///
/// The scalar stays *factored* all the way through — see [`Frac::mul_factors`].
fn scale<C: Ring>(f: Monomial<Frac<C>>, by: &Factors) -> Monomial<Frac<C>> {
    let mut out = Monomial::zero();
    for (mu, c) in f.terms() {
        let mut v = c.mul_factors(by);
        v.reduce();
        out.add_term(mu.clone(), v);
    }
    out
}

/// `b_λ = ∏_{s∈λ} (1 − q^{a} t^{l+1}) / (1 − q^{a+1} t^{l})`.
fn b_factors(lambda: &[u32]) -> Factors {
    let mut f = Factors::new();
    for_each_cell(lambda, &mut |a, l| {
        *f.entry((a, l + 1)).or_insert(0) += 1;
        *f.entry((a + 1, l)).or_insert(0) -= 1;
    });
    f.retain(|_, m| *m != 0);
    f
}

/// `c_λ = ∏_{s∈λ} (1 − q^{a} t^{l+1})` — the **numerator** of `b_λ`.
///
/// Not `(1 − q^{a+1} t^{l})`, which is `c'_λ`, the denominator. Writing that
/// one gives a `J_(2)` whose leading coefficient is `(1−q²)(1−q)` where Sage
/// has `(1−t)(1−qt)`; the polynomiality test below is what caught it.
pub(crate) fn c_factors(lambda: &[u32]) -> Factors {
    let mut f = Factors::new();
    for_each_cell(lambda, &mut |a, l| {
        *f.entry((a, l + 1)).or_insert(0) += 1;
    });
    f.retain(|_, m| *m != 0);
    f
}

/// `c'_λ = ∏_{s∈λ}(1 − q^{a(s)+1} t^{l(s)})` — the **denominator** of `b_λ`,
/// and the scalar Lapointe–Lascoux–Morse call `c_{λ'}(t,q)` in their 3.15.
///
/// Conjugating a shape swaps arms and legs, and swapping `q` with `t` then puts
/// `c_{λ'}(t,q) = ∏_{s∈λ}(1 − q^{a(s)+1} t^{l(s)})` — this. It is what the
/// eigenvector route has to multiply by, since that normalizes the leading
/// coefficient to 1 where `J` has this.
pub(crate) fn c_prime_factors(lambda: &[u32]) -> BTreeMap<(u32, u32), u32> {
    let mut f = BTreeMap::new();
    for_each_cell(lambda, &mut |a, l| {
        *f.entry((a + 1, l)).or_insert(0) += 1;
    });
    f
}

/// Every cell of λ, as its (arm, leg).
fn for_each_cell(lambda: &[u32], visit: &mut impl FnMut(u32, u32)) {
    for i in 0..lambda.len() {
        for j in 0..lambda[i] as usize {
            visit(arm(lambda, i, j), leg(lambda, i, j));
        }
    }
}

/// Binomial exponents with signed multiplicities: `(a, b) ↦ m` is
/// `(1 − qᵃtᵇ)^m`, negative meaning a denominator factor.
type Factors = BTreeMap<(u32, u32), i32>;

/// The factors of `ψ_{λ/μ}` for a horizontal strip, both shapes padded to one
/// length.
///
/// Returned unexpanded. Nothing here multiplies a polynomial: the whole of ψ is
/// a product of ratios of binomials, so it is entirely described by counting
/// them, and the counting lets matching factors cancel before any expansion.
fn psi_factors(lam: &[u32], mu: &[u32]) -> Factors {
    let width = lam.iter().copied().max().unwrap_or(0) as usize;
    // Column j meets the strip when λ has a cell there that μ does not.
    let met: Vec<bool> = (0..width)
        .map(|j| count_above(lam, j) > count_above(mu, j))
        .collect();

    let mut f = Factors::new();
    for i in 0..mu.len() {
        if lam[i] == mu[i] {
            continue; // row i does not meet the strip
        }
        for j in 0..mu[i] as usize {
            if met[j] {
                continue; // column j does meet it
            }
            let (am, lm) = (arm(mu, i, j), leg(mu, i, j));
            let (al, ll) = (arm(lam, i, j), leg(lam, i, j));
            // b_μ(s) = (1 − q^{am} t^{lm+1}) / (1 − q^{am+1} t^{lm})
            *f.entry((am, lm + 1)).or_insert(0) += 1;
            *f.entry((am + 1, lm)).or_insert(0) -= 1;
            // 1 / b_λ(s)
            *f.entry((al + 1, ll)).or_insert(0) += 1;
            *f.entry((al, ll + 1)).or_insert(0) -= 1;
        }
    }
    f.retain(|_, m| *m != 0);
    f
}

/// Cells strictly to the right of `(i, j)` in row `i`.
pub(crate) fn arm(shape: &[u32], i: usize, j: usize) -> u32 {
    shape[i] - j as u32 - 1
}

/// Cells strictly below `(i, j)` in column `j`.
pub(crate) fn leg(shape: &[u32], i: usize, j: usize) -> u32 {
    count_above(shape, j) as u32 - 1 - i as u32
}

/// The conjugate part `shape'_j` — how many rows reach past column `j`.
pub(crate) fn count_above(shape: &[u32], j: usize) -> usize {
    shape.iter().filter(|&&v| v as usize > j).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::qt::QtPoly;

    type F = Frac<Rational>;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn r(n: i128) -> Rational {
        Rational::from_int(n)
    }

    /// Hand-computed, and the case that pins both halves of the row/column
    /// condition in ψ.
    ///
    /// P_(2) = m_(2) + [(1−t)(1+q)/(1−qt)] m_(1,1). Sage:
    ///
    /// ```text
    ///   sage: P = SymmetricFunctions(QQ['q','t'].fraction_field()).macdonald().P()
    ///   sage: SymmetricFunctions(...).monomial()(P[2])
    ///   m[1, 1] + ... m[2]
    /// ```
    #[test]
    fn p_of_two_matches_the_hand_computation() {
        let p: Monomial<F> = macdonald_p(&part(&[2]));
        assert_eq!(p.coeff(&part(&[2])), <F as Ring>::one(), "monic in m_λ");

        let c = p.coeff(&part(&[1, 1]));
        // (1 − t)(1 + q) / (1 − q t), built independently of the ψ machinery
        let mut one_plus_q: QtPoly<Rational> = <QtPoly<Rational> as Ring>::one();
        one_plus_q.add_term(1, 0, r(1));
        let want = F::factor(0, 1) // 1 − t
            .mul(&F::from_poly(one_plus_q))
            .mul(&F::inv_factor(1, 1)); // / (1 − q t)
        assert_eq!(c, want, "got {c}, want {want}");
    }

    /// q = t collapses Macdonald to Schur. Checked in the monomial basis
    /// against the library's own s → m, which is tested by other means.
    #[test]
    fn at_q_equals_t_it_is_the_schur_function() {
        use crate::convert::FromSchur;
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = macdonald_p(&lambda);
                let s: Monomial<Rational> =
                    Monomial::from_schur(&crate::Schur::monomial(lambda.clone(), r(1)));
                // q = t = 2 is a generic point: no denominator 1 − qᵃtᵇ vanishes
                // there, unlike q = t = 1.
                for (mu, c) in p.terms() {
                    let got = c.eval(&r(2), &r(2)).expect("no pole at q = t = 2");
                    assert_eq!(got, s.coeff(mu), "P_{lambda} at q=t, term {mu}");
                }
                for (mu, want) in s.terms() {
                    if p.coeff(mu) == <F as Ring>::zero() {
                        assert!(want.is_zero(), "missing term {mu} of s_{lambda}");
                    }
                }
            }
        }
    }

    /// t = 1 collapses P to m_λ: the transition becomes the identity.
    #[test]
    fn at_t_equals_one_it_is_the_monomial() {
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = macdonald_p(&lambda);
                for (mu, c) in p.terms() {
                    let got = c.eval(&r(3), &r(1)).expect("no pole at q = 3, t = 1");
                    assert_eq!(
                        got,
                        Rational::from_int(i128::from(mu == &lambda)),
                        "P_{lambda} at t = 1, term {mu}"
                    );
                }
            }
        }
    }

    /// **q = 0 is Hall–Littlewood P.** The check this whole layer was built to
    /// make possible, and the strongest evidence available for either side:
    /// `macdonald_p` is a branching formula over rational functions, while
    /// `hall_littlewood_p` inverts the Kostka–Foulkes matrix produced by the
    /// Morris recursion. They share no code and no algorithm.
    ///
    /// Compared at several values of t rather than symbolically: `Frac::eval`
    /// substitutes numbers, and setting q = 0 while keeping t formal would need
    /// a separate exact division in `ℤ[t]`. Small values keep the exact
    /// rationals well inside i128 — the symbolic comparison is Sage's job.
    #[test]
    fn at_q_zero_it_is_hall_littlewood_p() {
        use crate::convert::FromSchur;
        use crate::qt::QtPoly;
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let mac: Monomial<F> = macdonald_p(&lambda);
                // Hall–Littlewood P arrives in the Schur basis; move it to m.
                let hl: crate::Schur<QtPoly<Rational>> = crate::hall_littlewood_p(&lambda);
                let hl: Monomial<QtPoly<Rational>> = Monomial::from_schur(&hl);
                for &tv in &[2i128, 3, 5] {
                    let t = r(tv);
                    for (mu, c) in mac.terms() {
                        let got = c.eval(&r(0), &t).expect("q = 0 is not a pole");
                        assert_eq!(
                            got,
                            hl.coeff(mu).eval(&Rational::from_int(0), &t),
                            "P_{lambda} at q=0, t={tv}, term {mu}"
                        );
                    }
                    for (mu, c) in hl.terms() {
                        if mac.coeff(mu) == <F as Ring>::zero() {
                            assert!(
                                c.eval(&Rational::from_int(0), &t).is_zero(),
                                "Hall–Littlewood has {mu} where Macdonald does not"
                            );
                        }
                    }
                }
            }
        }
    }

    /// Q = b_λ · P and J = c_λ · P, so both must agree with P after dividing
    /// the scalar back out — and J must be a *polynomial*, which is its whole
    /// point.
    #[test]
    fn q_and_j_are_scalar_multiples_of_p() {
        for n in 0..=5u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = macdonald_p(&lambda);
                let q: Monomial<F> = macdonald_q(&lambda);
                let j: Monomial<F> = macdonald_j(&lambda);
                // The scalar is read off the leading term, where P is 1.
                let b = q.coeff(&lambda);
                let c = j.coeff(&lambda);
                for (mu, pc) in p.terms() {
                    assert_eq!(q.coeff(mu), pc.mul(&b), "Q at {mu}");
                    assert_eq!(j.coeff(mu), pc.mul(&c), "J at {mu}");
                }
                // J has no denominator: c_λ clears exactly what P carries.
                for (mu, jc) in j.terms() {
                    let (_, den) = jc.parts();
                    assert_eq!(den.count(), 0, "J_{lambda} coefficient at {mu} is {jc}");
                }
            }
        }
    }

    /// Monic and strictly triangular in dominance order.
    #[test]
    fn expansion_is_unitriangular_in_dominance() {
        for n in 0..=7u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = macdonald_p(&lambda);
                assert_eq!(p.coeff(&lambda), <F as Ring>::one(), "monic at {lambda}");
                for mu in p.terms().keys() {
                    assert!(
                        crate::kostka::kostka(&lambda, mu) > 0,
                        "{mu} appears in P_{lambda} but λ does not dominate it"
                    );
                }
            }
        }
    }
}
