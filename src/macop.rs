//! The first Macdonald operator `M₁`, as a matrix on modified Schur functions.
//!
//! ## Reference
//!
//! L. Lapointe, A. Lascoux, J. Morse, *Determinantal expressions for Macdonald
//! polynomials*, International Mathematics Research Notices **1998** no. 18,
//! 957–978 (arXiv:math/9808050). Cited throughout as **[LLM]**, by their own
//! numbering — every bare `3.7` or `Theorem 2.1` below is theirs.
//!
//! The result this module exists for is [LLM] 3.6–3.7: on the basis
//! `S_μ[X^{tq}]` with `X^{tq} = X(t−1)/(q−1)`,
//!
//! ```text
//!   M₁ S_μ[X^{tq}] = Σ_α ε(μ,α) [|α|] S^α[X^{tq}],
//!   [|α|] = q^{α₁}t^{n−1} + q^{α₂}t^{n−2} + … + q^{α_n},
//! ```
//!
//! where `ε(μ,α)` is the sign of the Jacobi–Trudi permutation term and `S^α` is
//! the product of complete functions `h_{α₁}···h_{α_n}`. **α is indexed by the
//! column of the determinant, not the row** — see [`eigenvalue`], where getting
//! that wrong is both the easy reading and a silent wrong answer. The action is
//! **triangular** with `[|μ|]` on the diagonal, and the eigenvalues are
//! distinct, so `J_λ` is recoverable as the eigenvector for `[|λ|]` — with no
//! tableau enumeration anywhere. That is the point: the branching formula's cost
//! grows 4.7× per degree (see `docs/record/qt-kostka.md`), and this does not.
//!
//! ## Why the composition and not the multiset
//!
//! `S^α = h_{α₁}···h_{α_n}` does not care about the order of α — but `[|α|]`
//! does, since it pairs `α_i` with `t^{n−i}`. So the expansion has to be visited
//! per permutation term, which is what
//! [`jt_compositions`](crate::convert::jt_compositions) is for; the aggregating
//! `jt_terms` the classical conversions use would have summed away exactly the
//! information the eigenvalue needs.
//!
//! ## Padding, and why the answer does not depend on it
//!
//! `[|α|]` names an `n`, and every partition of the degree is padded to a common
//! length. LLM note (p.10) that the matrix entries are independent of that
//! choice because each is a *difference* of two eigenvalues, whose last
//! `n − ℓ(λ)` components agree and cancel. The entries produced here are the
//! eigenvalues themselves rather than differences, so they do depend on the
//! padding — the degree is used throughout, and the subtraction that makes it
//! irrelevant happens in the eigenvector solve.
//!
//! Padding a partition with zeros is also what makes the enumeration cheap. A
//! row with `c[row] = 0` admits only `j ≥ row`, so the zero rows of a short
//! partition force most of the permutation and the tree collapses: `μ = (n)`
//! padded to length `n` yields exactly one term, not `n!`.

use std::collections::BTreeMap;

use crate::coeff::Ring;
use crate::convert::jt_compositions;
use crate::partition::Partition;
use crate::qt::QtPoly;

/// `[|α|] = Σ_i q^{α_i} t^{n−i}`, the eigenvalue symbol, for α laid out by
/// **column** of the Jacobi–Trudi determinant.
///
/// [LLM] 3.2 (see the module docs for the citation). The subtlety is entirely
/// in what indexes α, and
/// it cost a wrong answer before a hand computation at `λ = (1,1)` found it.
///
/// [LLM] 3.5 writes `M₁` through the formal-operator notation of their 2.4: it
/// adds the alphabet `X^t` to **one column** of `det(S_{μ_i−i+j}[X^{tq}])` and
/// sums over which column, weighted `t^{n−i}`. Expanding that determinant, the
/// factor contributed by row `j` picks up `q^{α_j}` exactly when `σ(j)` is the
/// chosen column — so the power of `t` travels with the *column* `σ(j)`, and the
/// entry sitting at column `i` is the one from row `σ⁻¹(i)`. Their `α` is a
/// rearrangement `σ(μ+ρ)−ρ`, which is that column indexing.
///
/// Laying α out by row instead is the natural reading and is wrong. It is also
/// nearly undetectable: it still yields distinct eigenvalues, still yields a
/// dominance-triangular matrix, and still satisfies `M₁b = [|λ|]b` — that
/// equation only says the solve agrees with the matrix it was handed. It fails
/// at the first comparison with a real Macdonald polynomial: at `λ = (1,1)` the
/// coefficient ratio must be `−(q−t)/(1−qt)` and the row indexing forces 1,
/// whatever power of `t` is paired with whatever position.
///
/// Every monomial here is distinct — the `t`-exponent alone separates them — so
/// terms cannot collide and every coefficient is 1. That is what keeps
/// `[|λ|] − [|μ|]` unit-coefficient, which is what lets the route run over ℤ
/// rather than ℚ (see [`QtPoly::divide_exact`](crate::qt::QtPoly::divide_exact)).
pub(crate) fn eigenvalue<C: Ring>(alpha: &[u32]) -> QtPoly<C> {
    let n = alpha.len();
    let mut out = QtPoly::zero();
    for (i, &a) in alpha.iter().enumerate() {
        out.add_term(a, (n - 1 - i) as u32, C::one());
    }
    out
}

/// `[|λ|]` for a partition, padded to length `n`.
pub(crate) fn eigenvalue_of<C: Ring>(lambda: &Partition, n: usize) -> QtPoly<C> {
    let padded: Vec<u32> = (0..n).map(|i| lambda.part(i)).collect();
    eigenvalue(&padded)
}

/// The matrix of `M₁` on `{S_μ[X^{tq}]}` for one degree, as
/// `a[i][j]` = the coefficient of `S_{λⁱ}` in `M₁ S_{λʲ}`, indexed against
/// [`partitions_cached`](crate::memo::partitions_cached).
///
/// Column `j` is built in the `h` basis, where the expansion lands naturally,
/// and carried into Schur by the **Kostka matrix**: `h_ν = Σ_κ K_{κν} s_κ`.
/// `h_α = h_{sort α}`, so the composition collapses to a partition on the way
/// into the accumulator — after its eigenvalue has been read off it.
///
/// Going through [`Homogeneous::to_schur`](crate::convert) instead is the
/// obvious spelling. That routine expands `h_ν` as a product `∏ s_{(ν_i)}` of
/// Schur functions — a chain of Littlewood–Richardson products, run over
/// `QtPoly` coefficients, once per (column, term) pair, to recompute a
/// transition that depends on nothing but the degree. The Kostka table *is* that
/// transition, it is integral, and the crate already memoises it.
///
/// Worth **2.8×** on this step (0.0063s against 0.0178s at degree 10), which is
/// less than it sounds like it should be and is recorded because the first
/// version of this comment guessed "a factor of 30" without measuring. It is
/// also, for now, worth nothing at all: the matrix is 0.01% of the route, and
/// [`eigenvector`]'s solve is the other 99.99%.
pub fn operator_matrix<C: Ring>(n: u32) -> Vec<Vec<QtPoly<C>>> {
    let parts = crate::memo::partitions_cached(n);
    let index: std::collections::HashMap<&Partition, usize> =
        parts.iter().enumerate().map(|(i, p)| (p, i)).collect();
    let width = n as usize;
    let kostka = crate::kostka::kostka_table(n);
    let mut out = vec![vec![QtPoly::zero(); parts.len()]; parts.len()];

    let mut d: Vec<QtPoly<C>> = vec![QtPoly::zero(); parts.len()];
    for (j, mu) in parts.iter().enumerate() {
        for e in d.iter_mut() {
            *e = QtPoly::zero();
        }
        let padded: Vec<u32> = (0..width).map(|i| mu.part(i)).collect();
        jt_compositions(&padded, &mut |alpha, sign| {
            let ev: QtPoly<C> = eigenvalue(alpha);
            let nu = Partition::new(alpha.iter().copied());
            let slot = &mut d[index[&nu]];
            if sign < 0 {
                slot.sub_assign(&ev);
            } else {
                slot.add_assign(&ev);
            }
        });
        for (nu, dv) in d.iter().enumerate() {
            if dv.is_empty() {
                continue;
            }
            for kappa in 0..parts.len() {
                let k = kostka[kappa][nu];
                if k == 0 {
                    continue;
                }
                let scaled = dv.mul(&QtPoly::term(0, 0, C::from_u128(k)));
                out[kappa][j].add_assign(&scaled);
            }
        }
    }
    out
}

/// The eigenvector of `M₁` for the eigenvalue `[|λ|]`, cleared of denominators.
///
/// Returns `(b, v)` indexed against
/// [`partitions_cached`](crate::memo::partitions_cached), where `v` is the
/// common denominator and `Σ_κ (b_κ / v) S_κ[X^{tq}]` is `J_λ` up to the scalar
/// LLM call `c_{λ'}(t,q)` — this normalises the `S_λ` coefficient to 1 rather
/// than to theirs.
///
/// ## The solve
///
/// Triangularity ([LLM] Theorem 3.7, reproduced as a test below) turns the
/// eigenvector equation into back
/// substitution. Row κ of `(M₁ − [|λ|])a = 0` reads
///
/// ```text
///   ([|κ|] − [|λ|]) a_κ  =  − Σ_{κ ▷ μ ⊵ λ} A_{κμ} a_μ
/// ```
///
/// because `A_{κκ} = [|κ|]` and `A_{κμ}` vanishes unless κ ⊵ μ. Taking κ in any
/// linear extension of dominance — lexicographic ascending is one, since
/// λ ⊵ μ implies λ ≥ μ lexicographically — every `a_μ` on the right is already
/// known.
///
/// ## Why this stays in ℤ[q,t]
///
/// `a_κ` is a genuine rational function, so the recursion is run on
/// `b_κ = a_κ · v` with `v = ∏_{κ ▷ λ}([|κ|] − [|λ|])`. Each `a_κ`'s denominator
/// divides a sub-product of `v`, so every `b_κ` is a polynomial and every step
/// is an exact division — [`QtPoly::divide_exact`](crate::qt::QtPoly), never a
/// gcd and never a field. A `None` from it is a bug in this reasoning and is
/// raised as one, not swallowed.
pub fn eigenvector<C: Ring>(lambda: &Partition) -> (Vec<QtPoly<C>>, QtPoly<C>) {
    let n = lambda.size();
    let a: Vec<Vec<QtPoly<C>>> = operator_matrix(n);
    solve(n, &a, lambda)
}

/// Every `J_λ` of one degree, sharing the operator matrix.
///
/// The natural unit of work: `M₁` depends only on the degree, so building it
/// once and solving p(n) times is what a whole table costs — asking
/// [`eigenvector`] p(n) times rebuilds the matrix every time.
pub fn eigenvectors<C: Ring>(n: u32) -> Vec<(Partition, Vec<QtPoly<C>>, QtPoly<C>)> {
    let a: Vec<Vec<QtPoly<C>>> = operator_matrix(n);
    crate::memo::partitions_cached(n)
        .iter()
        .map(|lambda| {
            let (b, v) = solve(n, &a, lambda);
            (lambda.clone(), b, v)
        })
        .collect()
}

/// One eigenvector coefficient: a numerator over a **factored** denominator,
/// the factors drawn from the fixed family `gap[k] = [|κ_k|] − [|λ|]`.
///
/// This is [`Frac`](crate::Frac)'s design over a different family. `Frac` holds
/// denominators as a multiset of binomials `1 − qᵃtᵇ` because that class is
/// closed under products and lcms, which is all a sum needs, and so never
/// expands one. The same is true here: the divisors are p(n) known polynomials,
/// enumerated before the solve starts, so a denominator is a multiset of indices
/// into `gap` and no gcd is required to combine two of them.
///
/// It exists because the first version of this solve cleared every denominator
/// at once — `b_κ = a_κ · ∏_{κ ▷ λ} gap_κ` — which put 48,419 terms in `v` at
/// degree 10 against 5,630 in the entire operator matrix, and ran every one of
/// the p(n)³ products in the solve at that size. See
/// `docs/record/qt-kostka.md`.
#[derive(Clone)]
struct Coeff<C: Ring> {
    num: QtPoly<C>,
    den: BTreeMap<usize, u32>,
}

impl<C: Ring> Coeff<C> {
    fn zero() -> Self {
        Coeff {
            num: QtPoly::zero(),
            den: BTreeMap::new(),
        }
    }

    fn one() -> Self {
        Coeff {
            num: <QtPoly<C> as Ring>::one(),
            den: BTreeMap::new(),
        }
    }

    fn is_zero(&self) -> bool {
        self.num.is_empty()
    }

    /// `self` rewritten over `target`, which must be a multiple of `self.den`.
    fn lift(&self, target: &BTreeMap<usize, u32>, gap: &[QtPoly<C>]) -> QtPoly<C> {
        let mut num = self.num.clone();
        for (&k, &m) in target {
            for _ in 0..(m - self.den.get(&k).copied().unwrap_or(0)) {
                num = num.mul(&gap[k]);
            }
        }
        num
    }

    /// `self += other · p`, over the lcm of the two denominators.
    fn add_mul(&mut self, other: &Self, p: &QtPoly<C>, gap: &[QtPoly<C>]) {
        if other.is_zero() || p.is_empty() {
            return;
        }
        if self.is_zero() {
            self.num = other.num.mul(p);
            self.den = other.den.clone();
            return;
        }
        if other.den != self.den {
            let mut lcm = self.den.clone();
            for (&k, &m) in &other.den {
                let e = lcm.entry(k).or_insert(0);
                *e = (*e).max(m);
            }
            self.num = self.lift(&lcm, gap);
            self.den = lcm;
        }
        let lifted = other.lift(&self.den, gap);
        self.num.add_assign(&lifted.mul(p));
    }

    /// Divide out every denominator factor that also divides the numerator.
    ///
    /// The same trial division [`Frac::reduce`](crate::Frac::reduce) performs,
    /// and load-bearing for the same reason plus one more: without it the
    /// denominators only grow, every later row lifts against them, and the swell
    /// this type exists to avoid comes back.
    fn reduce(&mut self, gap: &[QtPoly<C>]) {
        if self.num.is_empty() {
            self.den.clear();
            return;
        }
        self.den.retain(|&k, m| {
            while *m > 0 {
                match self.num.divide_exact(&gap[k]) {
                    Some(q) => {
                        self.num = q;
                        *m -= 1;
                    }
                    None => break,
                }
            }
            *m > 0
        });
    }
}

fn solve<C: Ring>(n: u32, a: &[Vec<QtPoly<C>>], lambda: &Partition) -> (Vec<QtPoly<C>>, QtPoly<C>) {
    let parts = crate::memo::partitions_cached(n);
    let width = n as usize;
    let evs: Vec<QtPoly<C>> = parts.iter().map(|p| eigenvalue_of(p, width)).collect();
    let li = parts
        .iter()
        .position(|p| p == lambda)
        .expect("lambda must be a partition of its own size");

    // `gap[k]` is `[|κ|] − [|λ|]`, the divisor for row κ. Nonzero off the
    // diagonal because the eigenvalues are distinct.
    let gap: Vec<QtPoly<C>> = evs
        .iter()
        .map(|e| {
            let mut d = e.clone();
            d.sub_assign(&evs[li]);
            d
        })
        .collect();

    let above: Vec<usize> = (0..parts.len())
        .filter(|&k| k != li && crate::kostka::dominates(&parts[k], lambda))
        .collect();

    let mut coeff: Vec<Coeff<C>> = vec![Coeff::zero(); parts.len()];
    coeff[li] = Coeff::one();

    // Lexicographic ascending refines dominance, so every μ needed by row κ is
    // solved before κ is reached.
    let mut order = above.clone();
    order.sort_by(|&x, &y| parts[x].parts().cmp(parts[y].parts()));

    for &k in &order {
        let mut sum = Coeff::zero();
        for m in 0..parts.len() {
            if m == k {
                continue;
            }
            let prev = coeff[m].clone();
            sum.add_mul(&prev, &a[k][m], &gap);
        }
        sum.num = sum.num.neg();
        // a_κ = −(…)/gap_κ. The factor goes into the denominator and `reduce`
        // takes it straight back out whenever it divides, which is the common
        // case and is why this stays small.
        *sum.den.entry(k).or_insert(0) += 1;
        sum.reduce(&gap);
        coeff[k] = sum;
    }

    // Back to a common denominator for the caller — over the lcm of what
    // survived reduction, which is the point: the full product is never formed.
    let mut lcm: BTreeMap<usize, u32> = BTreeMap::new();
    for c in &coeff {
        for (&k, &m) in &c.den {
            let e = lcm.entry(k).or_insert(0);
            *e = (*e).max(m);
        }
    }
    let mut v = <QtPoly<C> as Ring>::one();
    for (&k, &m) in &lcm {
        for _ in 0..m {
            v = v.mul(&gap[k]);
        }
    }
    let b = coeff.iter().map(|c| c.lift(&lcm, &gap)).collect();
    (b, v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::sym::{Homogeneous, Schur, SymFn};

    /// [LLM] Theorem 3.7: the action is triangular with `[|μ|]` on the diagonal.
    ///
    /// This is the whole of step two under test, and it is a sharp check on both
    /// halves at once. The diagonal pins the eigenvalue convention — which way
    /// `t^{n−i}` runs against the parts, an off-by-one that would still produce a
    /// plausible triangular matrix. The zeros pin `ε(μ,α)`: a sign error or a
    /// dropped permutation term shows up as a nonzero below the diagonal, since
    /// the cancellations that produce those zeros are exactly the ones the signs
    /// are responsible for.
    #[test]
    fn the_action_is_triangular_with_the_eigenvalues_on_the_diagonal() {
        for n in 1..=7u32 {
            let parts = crate::partitions_of(n);
            let a: Vec<Vec<QtPoly<Rational>>> = operator_matrix(n);
            for (j, mu) in parts.iter().enumerate() {
                for (i, kappa) in parts.iter().enumerate() {
                    if i == j {
                        assert_eq!(a[i][j], eigenvalue_of(mu, n as usize), "diagonal at {mu}");
                    } else if !crate::kostka::dominates(kappa, mu) {
                        assert!(
                            a[i][j].is_empty(),
                            "M1 S_{mu} has a {kappa} term, but {kappa} does not dominate {mu}"
                        );
                    }
                }
            }
        }
    }

    /// The vector the solve returns must actually be an eigenvector — checked
    /// against the matrix, exactly, over ℤ[q,t].
    ///
    /// `M₁ b = [|λ|] b` in every row, including the rows below λ where both
    /// sides are zero. This is the milestone for the whole route: it says the
    /// enumeration, the eigenvalue convention, the dominance ordering of the
    /// back substitution and every exact division agree, with no oracle and no
    /// reference implementation involved. Whether the eigenvector is *`J_λ`* is
    /// a separate question about the basis, and a separate test.
    #[test]
    fn the_solution_really_is_an_eigenvector() {
        for n in 1..=7u32 {
            let parts = crate::partitions_of(n);
            let a: Vec<Vec<QtPoly<i64>>> = operator_matrix(n);
            for lambda in &parts {
                let (b, v) = eigenvector::<i64>(lambda);
                assert!(!v.is_empty(), "the common denominator must be nonzero");
                let ev = eigenvalue_of::<i64>(lambda, n as usize);
                for k in 0..parts.len() {
                    let mut lhs = QtPoly::zero();
                    for m in 0..parts.len() {
                        lhs.add_assign(&a[k][m].mul(&b[m]));
                    }
                    let rhs = ev.mul(&b[k]);
                    assert_eq!(lhs, rhs, "row {} of M1 b = [|{lambda}|] b", parts[k]);
                }
            }
        }
    }

    /// The eigenvector is supported exactly where [LLM] 3.15 says it is: on μ ⊵ λ,
    /// with λ itself present.
    ///
    /// Support is not implied by the eigenvector equation — the zero vector
    /// satisfies that too, and so would a solution that had collapsed onto the
    /// wrong shape. `b_λ = v ≠ 0` is what says the normalisation survived.
    #[test]
    fn the_eigenvector_is_supported_above_lambda() {
        for n in 1..=7u32 {
            let parts = crate::partitions_of(n);
            for (li, lambda) in parts.iter().enumerate() {
                let (b, v) = eigenvector::<i64>(lambda);
                assert_eq!(b[li], v, "the {lambda} coefficient is the normalisation");
                for (k, kappa) in parts.iter().enumerate() {
                    if !crate::kostka::dominates(kappa, lambda) {
                        assert!(b[k].is_empty(), "J_{lambda} should have no {kappa} term");
                    }
                }
            }
        }
    }

    /// The payoff: the eigenvector really is Macdonald's `J_λ`, checked against
    /// the branching formula.
    ///
    /// Two routes with nothing in common — one enumerates semistandard tableaux
    /// and multiplies out ψ, the other solves a linear system built from
    /// permutation signs — so agreement is the strongest evidence available for
    /// either.
    ///
    /// The basis has to be crossed first. `S_κ[X^{tq}]` with
    /// `X^{tq} = X(t−1)/(q−1)` is an ordinary Schur function of a modified
    /// alphabet, and in the power sums that modification is diagonal:
    /// `p_k ↦ p_k (1−t^k)/(1−q^k)`. Both factors are binomials `1 − qᵃtᵇ`, so
    /// [`Frac`] carries them without expanding anything — the same closure
    /// property the whole ℚ(q,t) design rests on.
    ///
    /// Compared by **cross-multiplication** rather than by normalising. An
    /// eigenvector is only defined up to a scalar, so requiring a particular one
    /// would be testing a convention rather than the mathematics; and dividing
    /// by `v` is not available anyway, since `v` is a product of
    /// `[|κ|] − [|λ|]` and leaves the class of denominators `Frac` can hold.
    #[test]
    fn the_eigenvector_is_the_integral_form() {
        use crate::convert::{FromSchur, ToSchur};
        use crate::frac::Frac;
        use crate::sym::{Monomial, PowerSum};
        use std::collections::BTreeMap;

        /// `f[X(t−1)/(q−1)]`: multiply the `p_ρ` coefficient by
        /// `∏ (1−t^{ρ_i})/(1−q^{ρ_i})`.
        fn modify(f: &Schur<Frac<Rational>>) -> Schur<Frac<Rational>> {
            let p: PowerSum<Frac<Rational>> = PowerSum::from_schur(f);
            let mut out = PowerSum::zero();
            for (rho, c) in p.terms() {
                let mut factors: BTreeMap<(u32, u32), i32> = BTreeMap::new();
                for &k in rho.parts() {
                    *factors.entry((0, k)).or_insert(0) += 1; // 1 − t^k
                    *factors.entry((k, 0)).or_insert(0) -= 1; // / (1 − q^k)
                }
                let mut v = c.mul_factors(&factors);
                v.reduce();
                out.add_term(rho.clone(), v);
            }
            out.to_schur()
        }

        for n in 1..=5u32 {
            let parts = crate::partitions_of(n);
            for lambda in &parts {
                let (b, v) = eigenvector::<Rational>(lambda);

                // Σ_κ b_κ S_κ[X^{tq}], carried across into ordinary Schur.
                let mut formal: Schur<Frac<Rational>> = Schur::zero();
                for (k, kappa) in parts.iter().enumerate() {
                    formal.add_term(kappa.clone(), Frac::from_poly(b[k].clone()));
                }
                let got = modify(&formal);

                // v · J_λ, from the branching formula.
                let j: Monomial<Frac<Rational>> = crate::macdonald_j(lambda);
                let j: Schur<Frac<Rational>> = j.to_schur();
                let scale = Frac::from_poly(v.clone());
                let want: Schur<Frac<Rational>> = Schur::from_terms(
                    j.terms()
                        .iter()
                        .map(|(k, c)| (k.clone(), c.mul(&scale)))
                        .collect(),
                );

                // Proportional, cross-multiplied against the λ coefficient.
                let (gl, wl) = (got.coeff(lambda), want.coeff(lambda));
                assert!(
                    !gl.is_zero() && !wl.is_zero(),
                    "{lambda} coefficient vanished"
                );
                let keys: std::collections::BTreeSet<_> =
                    got.terms().keys().chain(want.terms().keys()).collect();
                for kappa in keys {
                    assert_eq!(
                        got.coeff(kappa).mul(&wl),
                        want.coeff(kappa).mul(&gl),
                        "J_{lambda} at {kappa}"
                    );
                }
            }
        }
    }

    /// The solve is exact in fixed width, checked by running it twice.
    ///
    /// `impl_ring_for_int` multiplies with a plain `*`, so a wrap is silent —
    /// the only way to know is to compute the same thing in two widths. `i64`
    /// and `i128` agreeing means the answer fits comfortably in the narrower
    /// one, which is a much stronger statement than either alone.
    /// `examples/llm_coeff_sizes.rs` carries the measured widths: ~6 bits per
    /// degree, 50 bits at degree 12, so `i64` holds to about degree 14 and
    /// `i128` well past where the enumeration is feasible.
    #[test]
    fn the_solve_is_exact_in_fixed_width() {
        for n in 1..=8u32 {
            let narrow: Vec<(_, Vec<QtPoly<i64>>, QtPoly<i64>)> = eigenvectors(n);
            let wide: Vec<(_, Vec<QtPoly<i128>>, QtPoly<i128>)> = eigenvectors(n);
            for ((lambda, x, _), (_, y, _)) in narrow.iter().zip(wide.iter()) {
                for (p, q) in x.iter().zip(y.iter()) {
                    let widened: Vec<_> = p.terms().map(|(k, c)| (*k, *c as i128)).collect();
                    let got: Vec<_> = q.terms().map(|(k, c)| (*k, *c)).collect();
                    assert_eq!(widened, got, "J_{lambda} differs between i64 and i128");
                }
            }
        }
    }

    /// The eigenvalues must be **distinct**, or the eigenvector solve has no
    /// unique answer and the divisions in it are by zero.
    ///
    /// [LLM] assert this in passing after their 3.7; it is cheap to check and the
    /// entire recursion rests on it.
    #[test]
    fn the_eigenvalues_are_distinct() {
        for n in 1..=8u32 {
            let parts = crate::partitions_of(n);
            let evs: Vec<QtPoly<Rational>> =
                parts.iter().map(|p| eigenvalue_of(p, n as usize)).collect();
            for i in 0..evs.len() {
                for j in (i + 1)..evs.len() {
                    assert_ne!(evs[i], evs[j], "{} and {}", parts[i], parts[j]);
                }
            }
        }
    }

    /// The expansion must reproduce `S_μ` itself when the eigenvalue is dropped.
    ///
    /// `Σ_α ε(μ,α) h_{sort α}` is the ordinary `s → h` transition, so running
    /// `jt_compositions` with the eigenvalue replaced by 1 has to give back
    /// exactly `s_μ`. That isolates the new enumeration from the new eigenvalue:
    /// if this passes and the triangularity test fails, the fault is in
    /// `[|α|]`, and if this fails the fault is in the permutation walk.
    #[test]
    fn the_composition_walk_reproduces_the_s_to_h_transition() {
        use crate::convert::{FromSchur, ToSchur};
        for n in 1..=7u32 {
            for mu in crate::partitions_of(n) {
                let padded: Vec<u32> = (0..n as usize).map(|i| mu.part(i)).collect();
                let mut got: Homogeneous<Rational> = Homogeneous::zero();
                jt_compositions(&padded, &mut |alpha, sign| {
                    got.add_term(
                        Partition::new(alpha.iter().copied()),
                        Rational::from_int(i128::from(sign)),
                    );
                });
                let want: Homogeneous<Rational> =
                    Homogeneous::from_schur(&Schur::monomial(mu.clone(), Rational::from_int(1)));
                assert_eq!(got, want, "s -> h at {mu}");
            }
        }
    }
}
