//! Evaluation at a finite alphabet, and the principal specializations.
//!
//! Everything else in this crate computes *with* symmetric functions as formal
//! objects. This module is the bridge back to concrete numbers: given
//! x_1, …, x_n, what is f(x_1, …, x_n)?
//!
//! Two quite different things live here, and the distinction matters:
//!
//! * **[`Schur::eval`] and friends** — evaluation at an *arbitrary* alphabet,
//!   generic over [`Ring`]. Each basis uses the algorithm that is natural for
//!   it, which is why this is per-basis rather than "convert to X, then
//!   evaluate": p, e, and h are products of one-variable generators and cost
//!   almost nothing, while Schur uses the branching rule.
//! * **[`principal_specialization`] and [`principal_specialization_q`]** — the
//!   *closed forms* for the specific alphabets (1, 1, …, 1) and
//!   (1, q, q², …, q^{n−1}). These do not enumerate anything: they are products
//!   over the cells of λ. Where they apply they are incomparably faster than
//!   evaluation, which is the whole reason to special-case them.
//!
//! Note the deliberate asymmetry in the bialternant's absence. s_λ = a_{λ+δ}/a_δ
//! is the textbook formula and would be an O(n³) determinant, but it needs
//! *division* and it is 0/0 whenever two of the x_i coincide — so it is neither
//! generic over `Ring` nor correct on repeated alphabets. The branching rule
//! below is slower on generic input and always right.

use crate::coeff::Ring;
use crate::fasthash::Map;
use crate::partition::Partition;
use crate::sym::{Elementary, Forgotten, Homogeneous, Monomial, PowerSum, Schur, SymFn};
use crate::convert::ToSchur;

/// `x^k`, by repeated multiplication. `Ring` has no `pow`, and adding one for
/// this would be adding an operation every coefficient type must implement to
/// serve a single caller. Exponents here are bounded by |λ|.
fn powers<C: Ring>(x: &C, upto: u32) -> Vec<C> {
    let mut out = Vec::with_capacity(upto as usize + 1);
    out.push(C::one());
    for k in 1..=upto as usize {
        out.push(out[k - 1].mul(x));
    }
    out
}

// --- the multiplicative bases: one generator at a time ----------------------

impl<C: Ring> PowerSum<C> {
    /// p_λ(x_1, …, x_n) = ∏_i (Σ_j x_j^{λ_i}).
    pub fn eval(&self, xs: &[C]) -> C {
        if self.is_zero() {
            return C::zero();
        }
        let top = self.terms().keys().map(|l| l.part(0)).max().unwrap_or(0);
        // p_k for every k that appears, computed once across all terms.
        let pk: Vec<C> = {
            let mut acc = vec![C::zero(); top as usize + 1];
            for x in xs {
                for (k, p) in powers(x, top).into_iter().enumerate() {
                    acc[k].add_assign(&p);
                }
            }
            acc
        };
        let mut total = C::zero();
        for (lambda, c) in self.terms() {
            let mut term = c.clone();
            for &part in lambda.parts() {
                term = term.mul(&pk[part as usize]);
            }
            total.add_assign(&term);
        }
        total
    }
}

/// The elementary symmetric polynomials e_0, …, e_{upto} of `xs`, by the
/// one-variable-at-a-time recurrence e_k(x_1..x_j) = e_k(x_1..x_{j−1}) +
/// x_j·e_{k−1}(x_1..x_{j−1}) — i.e. "does x_j appear in this product?".
fn elementary_of<C: Ring>(xs: &[C], upto: u32) -> Vec<C> {
    let mut e = vec![C::zero(); upto as usize + 1];
    e[0] = C::one();
    for x in xs {
        // Descending, so each e[k−1] read is still the previous variable's.
        for k in (1..=upto as usize).rev() {
            let add = e[k - 1].mul(x);
            e[k].add_assign(&add);
        }
    }
    e
}

/// The complete homogeneous polynomials h_0, …, h_{upto} of `xs`. Same shape of
/// recurrence, ascending rather than descending because h_k(x_1..x_j) =
/// h_k(x_1..x_{j−1}) + x_j·h_{k−1}(x_1..x_j) lets x_j repeat.
fn homogeneous_of<C: Ring>(xs: &[C], upto: u32) -> Vec<C> {
    let mut h = vec![C::zero(); upto as usize + 1];
    h[0] = C::one();
    for x in xs {
        for k in 1..=upto as usize {
            let add = h[k - 1].mul(x);
            h[k].add_assign(&add);
        }
    }
    h
}

/// Evaluate a basis whose elements are products of one-row generators.
fn eval_multiplicative<C: Ring, S: SymFn<C>>(f: &S, gens: &[C]) -> C {
    let mut total = C::zero();
    for (lambda, c) in f.terms() {
        let mut term = c.clone();
        for &part in lambda.parts() {
            term = term.mul(&gens[part as usize]);
        }
        total.add_assign(&term);
    }
    total
}

impl<C: Ring> Elementary<C> {
    /// e_λ(x_1, …, x_n) = ∏_i e_{λ_i}(x).
    ///
    /// Note e_k = 0 for k > n, so this is 0 as soon as λ has a part exceeding
    /// the alphabet size — the recurrence produces that on its own.
    pub fn eval(&self, xs: &[C]) -> C {
        if self.is_zero() {
            return C::zero();
        }
        let top = self.terms().keys().map(|l| l.part(0)).max().unwrap_or(0);
        eval_multiplicative(self, &elementary_of(xs, top))
    }
}

impl<C: Ring> Homogeneous<C> {
    /// h_λ(x_1, …, x_n) = ∏_i h_{λ_i}(x).
    pub fn eval(&self, xs: &[C]) -> C {
        if self.is_zero() {
            return C::zero();
        }
        let top = self.terms().keys().map(|l| l.part(0)).max().unwrap_or(0);
        eval_multiplicative(self, &homogeneous_of(xs, top))
    }
}

// --- the monomial basis: its own definition ---------------------------------

impl<C: Ring> Monomial<C> {
    /// m_λ(x_1, …, x_n) = Σ_α x^α over the **distinct** rearrangements α of λ
    /// into n slots.
    ///
    /// Generated rather than filtered: parts are grouped by value and each slot
    /// picks a value from the multiset, so equal parts never produce the same α
    /// twice. That is the same device as `convert::muir_expand`, for the same
    /// reason — "distinct rearrangements" is a statement about a multiset, and
    /// enumerating permutations of a list and deduplicating afterwards would do
    /// ℓ! work to emit ℓ!/∏mᵢ! terms.
    pub fn eval(&self, xs: &[C]) -> C {
        let mut total = C::zero();
        for (lambda, c) in self.terms() {
            if lambda.len() > xs.len() {
                continue; // no room: m_λ vanishes in fewer than ℓ(λ) variables
            }
            let pows: Vec<Vec<C>> = xs.iter().map(|x| powers(x, lambda.part(0))).collect();
            let mut avail: Vec<(u32, u32)> = Vec::new();
            for &k in lambda.parts() {
                match avail.last_mut() {
                    Some((v, n)) if *v == k => *n += 1,
                    _ => avail.push((k, 1)),
                }
            }
            let mut acc = C::zero();
            rearrange(&pows, 0, &mut avail, lambda.len(), &C::one(), &mut acc);
            total.add_assign(&c.mul(&acc));
        }
        total
    }
}

fn rearrange<C: Ring>(
    pows: &[Vec<C>],
    slot: usize,
    avail: &mut [(u32, u32)],
    left: usize,
    run: &C,
    acc: &mut C,
) {
    if left == 0 {
        acc.add_assign(run);
        return;
    }
    // Every remaining part needs a slot of its own.
    if pows.len() - slot < left {
        return;
    }
    // This slot gets exponent 0 …
    rearrange(pows, slot + 1, avail, left, run, acc);
    // … or one of the still-unplaced distinct part values.
    for i in 0..avail.len() {
        if avail[i].1 == 0 {
            continue;
        }
        avail[i].1 -= 1;
        let next = run.mul(&pows[slot][avail[i].0 as usize]);
        rearrange(pows, slot + 1, avail, left - 1, &next, acc);
        avail[i].1 += 1;
    }
}

// --- the Schur basis: the branching rule ------------------------------------

impl<C: Ring> Schur<C> {
    /// s_λ(x_1, …, x_n) = Σ over semistandard tableaux of shape λ with entries
    /// in 1..n, of x^{content}.
    ///
    /// Computed by the **branching rule** rather than by summing over tableaux:
    /// a tableau is exactly a chain ∅ = ν⁰ ⊆ ν¹ ⊆ … ⊆ νⁿ = λ whose successive
    /// differences are horizontal strips (νᵏ is the set of cells holding an
    /// entry ≤ k), and the weight ∏ x_k^{|νᵏ/νᵏ⁻¹|} depends only on the strip
    /// sizes. So sweeping variable by variable and keeping a frontier of
    /// *shapes* collapses all tableaux sharing a prefix into one number, and
    /// the cost is the number of shapes inside λ rather than the number of
    /// tableaux — which is the difference between polynomial and exponential.
    ///
    /// This is the same chain DP as [`crate::kostka`], carrying ring elements
    /// instead of counts; K_{λμ} is what it degenerates to when the x_i are
    /// formal.
    pub fn eval(&self, xs: &[C]) -> C {
        let mut total = C::zero();
        for (lambda, c) in self.terms() {
            let v = eval_schur_at(lambda, xs);
            if !v.is_zero() {
                total.add_assign(&c.mul(&v));
            }
        }
        total
    }
}

fn eval_schur_at<C: Ring>(lambda: &Partition, xs: &[C]) -> C {
    let l = lambda.len();
    if l == 0 {
        return C::one();
    }
    // A chain of k horizontal strips reaches λ from ν iff λ_{j+k} ≤ ν_j for
    // every j. With ν = ∅ and k = n that is ℓ(λ) ≤ n, the familiar statement
    // that s_λ vanishes in fewer variables than λ has rows.
    if l > xs.len() {
        return C::zero();
    }
    let lam: Vec<u32> = lambda.parts().to_vec();
    let n = xs.len();

    let mut cur: Map<Vec<u32>, C> = Map::default();
    cur.insert(vec![0u32; l], C::one());
    let mut scratch = vec![0u32; l];

    for (step, x) in xs.iter().enumerate() {
        let pw = powers(x, lambda.size());
        // Strips still to come after this one.
        let remaining = n - step - 1;
        let mut next: Map<Vec<u32>, C> = Map::default();
        for (mu, c) in cur.drain() {
            grow_strip(&lam, &mu, 0, u32::MAX, &mut scratch, 0, &mut |nu, added| {
                // Prune anything the remaining variables can no longer complete.
                if !reachable(&lam, nu, remaining) {
                    return;
                }
                let v = c.mul(&pw[added as usize]);
                if v.is_zero() {
                    return;
                }
                match next.get_mut(nu) {
                    Some(e) => e.add_assign(&v),
                    None => {
                        next.insert(nu.to_vec(), v);
                    }
                }
            });
        }
        cur = next;
    }
    cur.remove(&lam).unwrap_or_else(C::zero)
}

/// Can `k` further horizontal strips carry `nu` up to `lam`? Exactly when
/// λ_{j+k} ≤ ν_j for every j — the k-strip analogue of ℓ(λ) ≤ n.
fn reachable(lam: &[u32], nu: &[u32], k: usize) -> bool {
    for j in 0..lam.len() {
        let need = if j + k < lam.len() { lam[j + k] } else { 0 };
        if nu[j] < need {
            return false;
        }
    }
    true
}

/// Every ν with μ ⊆ ν ⊆ λ and ν/μ a horizontal strip, with the number of cells
/// added.
///
/// The interlacing ν_1 ≥ μ_1 ≥ ν_2 ≥ μ_2 ≥ … makes the rows *independent* once
/// scanned top-down: row i is free in [μ_i, min(μ_{i−1}, λ_i)] whatever the
/// other rows do. So this is a plain product of ranges — no capacity
/// bookkeeping, unlike `kostka::grow`, which additionally has to place a fixed
/// number of cells.
fn grow_strip(
    lam: &[u32],
    mu: &[u32],
    i: usize,
    prev: u32,
    cur: &mut Vec<u32>,
    added: u32,
    emit: &mut impl FnMut(&[u32], u32),
) {
    if i == lam.len() {
        emit(cur, added);
        return;
    }
    let hi = prev.min(lam[i]);
    for v in mu[i]..=hi {
        cur[i] = v;
        grow_strip(lam, mu, i + 1, mu[i], cur, added + (v - mu[i]), emit);
    }
}

impl<C: Ring> Forgotten<C> {
    /// f_λ(x_1, …, x_n), via the Schur expansion.
    ///
    /// ω is an automorphism of the ring of symmetric functions in *infinitely*
    /// many variables and does not commute with restricting the alphabet, so
    /// there is no shortcut of the form "evaluate m_λ and fix it up": the
    /// expansion has to happen first.
    pub fn eval(&self, xs: &[C]) -> C {
        self.to_schur().eval(xs)
    }
}

// --- closed-form specializations --------------------------------------------

/// Hook lengths of λ, row-major, in the same order as the cells are visited by
/// the specialization products.
fn hooks(lambda: &Partition) -> Vec<u32> {
    let conj = lambda.conjugate();
    let mut out = Vec::with_capacity(lambda.size() as usize);
    for (i, &row) in lambda.parts().iter().enumerate() {
        for j in 0..row as usize {
            out.push((row - j as u32) + (conj.part(j) - i as u32) - 1);
        }
    }
    out
}

/// f^λ, the number of standard Young tableaux of shape λ — equivalently the
/// dimension of the irreducible S_{|λ|} representation indexed by λ.
///
/// The hook length formula, |λ|! / ∏_u h(u). The division is exact, and it is
/// performed *incrementally* rather than by forming |λ|! first: multiplying by
/// k and dividing by whichever hooks then divide keeps the running value near
/// the answer instead of near |λ|!, which overflows far sooner. `None` if even
/// the answer exceeds `u128`.
pub fn dimension(lambda: &Partition) -> Option<u128> {
    let n = lambda.size();
    if n == 0 {
        return Some(1);
    }
    let mut hs = hooks(lambda);
    hs.sort_unstable_by(|a, b| b.cmp(a));
    let mut acc: u128 = 1;
    for k in 1..=u128::from(n) {
        acc = acc.checked_mul(k)?;
        // Cancel any hook that now divides, largest first.
        for h in hs.iter_mut() {
            if *h != 0 && acc % u128::from(*h) == 0 {
                acc /= u128::from(*h);
                *h = 0;
            }
        }
    }
    // Any hook left un-cancelled must still divide the product.
    for h in hs {
        if h != 0 {
            acc /= u128::from(h);
        }
    }
    Some(acc)
}

/// s_λ(1, 1, …, 1) with `n` ones — the number of semistandard tableaux of shape
/// λ with entries in 1..n, and the dimension of the GL_n irreducible.
///
/// Stanley's content formula, ∏_u (n + c(u)) / h(u) where c(u) = j − i. Zero
/// exactly when ℓ(λ) > n, which the product reports on its own: the cell at the
/// bottom of the first column then has content −ℓ(λ) + 1 ≤ −n.
///
/// `None` on `u128` overflow of the *numerator product*, which can happen well
/// before the answer would — the divisions are interleaved but not perfectly.
pub fn principal_specialization(lambda: &Partition, n: u32) -> Option<u128> {
    if lambda.len() as u32 > n {
        return Some(0);
    }
    let mut hs = hooks(lambda);
    hs.sort_unstable_by(|a, b| b.cmp(a));
    let mut acc: u128 = 1;
    for (i, &row) in lambda.parts().iter().enumerate() {
        for j in 0..row as usize {
            // n + c(u) > 0 here, since ℓ(λ) ≤ n was checked above.
            let f = u128::from(n) + j as u128 - i as u128;
            acc = acc.checked_mul(f)?;
            for h in hs.iter_mut() {
                if *h != 0 && acc % u128::from(*h) == 0 {
                    acc /= u128::from(*h);
                    *h = 0;
                }
            }
        }
    }
    for h in hs {
        if h != 0 {
            acc /= u128::from(h);
        }
    }
    Some(acc)
}

/// s_λ(1, q, q², …, q^{n−1}) as the coefficient list of a polynomial in q,
/// lowest degree first.
///
/// The q-analogue of [`principal_specialization`], from the same product:
///
/// ```text
///   s_λ(1, q, …, q^{n−1}) = q^{n(λ)} · ∏_u (1 − q^{n + c(u)}) / (1 − q^{h(u)})
/// ```
///
/// with n(λ) = Σ (i−1)λ_i. Setting q = 1 recovers the content formula, and the
/// coefficient of q^k counts semistandard tableaux of shape λ with entries in
/// 1..n and Σ(entries − 1) = k.
///
/// Neither factor divides the other cell-by-cell, so the quotient is taken once
/// at the end. Both products have constant term 1, which makes the division a
/// truncated power-series inversion — no leading-coefficient case analysis, and
/// exact in ℤ because the quotient is known in advance to be a polynomial.
pub fn principal_specialization_q(lambda: &Partition, n: u32) -> Vec<i128> {
    if lambda.len() as u32 > n {
        return Vec::new();
    }
    if lambda.is_empty() {
        return vec![1];
    }
    let mut num = vec![1i128];
    for (i, &row) in lambda.parts().iter().enumerate() {
        for j in 0..row as usize {
            let a = (n as i64 + j as i64 - i as i64) as usize;
            num = mul_one_minus_q(&num, a);
        }
    }
    let mut den = vec![1i128];
    for h in hooks(lambda) {
        den = mul_one_minus_q(&den, h as usize);
    }

    let qdeg = (num.len() - 1) - (den.len() - 1);
    let mut q = vec![0i128; qdeg + 1];
    for k in 0..=qdeg {
        let mut v = *num.get(k).unwrap_or(&0);
        for (j, d) in den.iter().enumerate().skip(1) {
            if j > k {
                break;
            }
            v -= d * q[k - j];
        }
        q[k] = v; // den[0] == 1, so no division
    }

    // The q^{n(λ)} prefactor.
    let shift: usize = lambda
        .parts()
        .iter()
        .enumerate()
        .map(|(i, &r)| i * r as usize)
        .sum();
    let mut out = vec![0i128; shift];
    out.extend_from_slice(&q);
    out
}

/// Multiply a polynomial by (1 − q^a).
fn mul_one_minus_q(p: &[i128], a: usize) -> Vec<i128> {
    let mut out = vec![0i128; p.len() + a];
    out[..p.len()].copy_from_slice(p);
    for (i, c) in p.iter().enumerate() {
        out[i + a] -= c;
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::convert::FromSchur;
    use crate::memo::partitions_cached;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn ints(v: &[i64]) -> Vec<i64> {
        v.to_vec()
    }

    /// Every basis must agree with every other at the same alphabet, since they
    /// are five descriptions of one function. This is the check that catches a
    /// wrong recurrence in any single one of them: the per-basis algorithms
    /// share no code, so they can only agree by all being right.
    #[test]
    fn all_bases_agree_at_the_same_alphabet() {
        let xs = ints(&[2, -1, 3, 1, -2]);
        for n in 1..=6u32 {
            for lambda in partitions_cached(n).iter() {
                let s: Schur<i64> = Schur::monomial(lambda.clone(), 1);
                let want = s.eval(&xs);

                let m: Monomial<i64> = Monomial::from_schur(&s);
                assert_eq!(m.eval(&xs), want, "m for {lambda}");
                let e: Elementary<i64> = Elementary::from_schur(&s);
                assert_eq!(e.eval(&xs), want, "e for {lambda}");
                let h: Homogeneous<i64> = Homogeneous::from_schur(&s);
                assert_eq!(h.eval(&xs), want, "h for {lambda}");
                let f: Forgotten<i64> = Forgotten::from_schur(&s);
                assert_eq!(f.eval(&xs), want, "f for {lambda}");

                let sr: Schur<Rational> = Schur::monomial(lambda.clone(), Rational::from_int(1));
                let p: PowerSum<Rational> = PowerSum::from_schur(&sr);
                let xr: Vec<Rational> = xs.iter().map(|&v| Rational::from_int(v.into())).collect();
                assert_eq!(p.eval(&xr), Rational::from_int(want.into()), "p for {lambda}");
            }
        }
    }

    /// A repeated alphabet is the case the bialternant cannot do: a_δ vanishes
    /// when two variables coincide, so s_λ = a_{λ+δ}/a_δ is 0/0 there. The
    /// branching rule has no such hole, and this pins that it does not.
    #[test]
    fn repeated_and_zero_alphabet_entries_are_fine() {
        let xs = ints(&[2, 2, 2, 0, 5]);
        for n in 1..=6u32 {
            for lambda in partitions_cached(n).iter() {
                let s: Schur<i64> = Schur::monomial(lambda.clone(), 1);
                let h: Homogeneous<i64> = Homogeneous::from_schur(&s);
                assert_eq!(s.eval(&xs), h.eval(&xs), "{lambda}");
            }
        }
    }

    /// s_λ vanishes in fewer than ℓ(λ) variables, and e_k in fewer than k.
    #[test]
    fn vanishing_below_the_row_count() {
        let xs = ints(&[3, 1, 4]);
        for n in 1..=7u32 {
            for lambda in partitions_cached(n).iter() {
                let s: Schur<i64> = Schur::monomial(lambda.clone(), 1);
                if lambda.len() > xs.len() {
                    assert_eq!(s.eval(&xs), 0, "s_{lambda} in 3 variables");
                }
            }
        }
        let e4: Elementary<i64> = Elementary::monomial(part(&[4]), 1);
        assert_eq!(e4.eval(&xs), 0);
    }

    /// Evaluating at n ones must reproduce the content formula, and the two
    /// arrive by completely different routes — a chain DP over shapes versus a
    /// product over cells.
    #[test]
    fn evaluation_at_ones_matches_the_content_formula() {
        for deg in 1..=8u32 {
            for lambda in partitions_cached(deg).iter() {
                for n in 1..=6usize {
                    let ones = vec![1i64; n];
                    let s: Schur<i64> = Schur::monomial(lambda.clone(), 1);
                    let want = principal_specialization(lambda, n as u32).unwrap();
                    assert_eq!(s.eval(&ones) as u128, want, "s_{lambda}(1^{n})");
                }
            }
        }
    }

    /// The q-specialization, checked two ways: term by term against a genuine
    /// evaluation at (1, q, …) for a concrete q, and at q = 1 against the
    /// content formula.
    #[test]
    fn q_specialization_matches_evaluation_and_collapses_at_one() {
        for deg in 1..=7u32 {
            for lambda in partitions_cached(deg).iter() {
                for n in 1..=5usize {
                    let coeffs = principal_specialization_q(lambda, n as u32);

                    // q = 1: the sum of coefficients is s_λ(1^n).
                    let at_one: i128 = coeffs.iter().sum();
                    assert_eq!(
                        at_one as u128,
                        principal_specialization(lambda, n as u32).unwrap(),
                        "q=1 for {lambda}, n={n}"
                    );

                    // q = 2: evaluate the alphabet (1, 2, 4, …) directly and
                    // compare against the polynomial's own value there.
                    let alphabet: Vec<i64> = (0..n).map(|k| 1i64 << k).collect();
                    let s: Schur<i64> = Schur::monomial(lambda.clone(), 1);
                    let direct = s.eval(&alphabet) as i128;
                    let from_poly: i128 =
                        coeffs.iter().enumerate().map(|(k, c)| c << k).sum();
                    assert_eq!(from_poly, direct, "q=2 for {lambda}, n={n}");
                }
            }
        }
    }

    /// The q-specialization's coefficients count tableaux, so they are
    /// non-negative — a property no step of the construction makes obvious,
    /// since both the numerator and the denominator have negative coefficients
    /// and the division is a power-series inversion.
    #[test]
    fn q_specialization_coefficients_are_non_negative() {
        for deg in 1..=8u32 {
            for lambda in partitions_cached(deg).iter() {
                for n in 1..=5u32 {
                    for (k, c) in principal_specialization_q(lambda, n).iter().enumerate() {
                        assert!(*c >= 0, "q^{k} of s_{lambda}(1,q,..q^{n}) is {c}");
                    }
                }
            }
        }
    }

    /// dim λ against the Kostka number K_{λ,1^n}, which counts the same standard
    /// tableaux by an unrelated algorithm.
    #[test]
    fn dimension_matches_the_kostka_count() {
        for deg in 1..=10u32 {
            for lambda in partitions_cached(deg).iter() {
                let ones = Partition::new(std::iter::repeat(1).take(deg as usize));
                assert_eq!(
                    dimension(lambda).unwrap(),
                    crate::kostka::kostka(lambda, &ones),
                    "dim {lambda}"
                );
            }
        }
    }

    /// Σ_λ (dim λ)² = n!, the decomposition of the regular representation.
    #[test]
    fn dimensions_square_sum_to_the_factorial() {
        for n in 1..=12u32 {
            let total: u128 = partitions_cached(n)
                .iter()
                .map(|l| {
                    let d = dimension(l).unwrap();
                    d * d
                })
                .sum();
            let factorial: u128 = (1..=u128::from(n)).product();
            assert_eq!(total, factorial, "n = {n}");
        }
    }

    /// The incremental cancellation must not overflow where the answer fits.
    ///
    /// For the staircase λ = (10,9,…,1), f^λ has 35 digits — inside `u128`,
    /// which holds 39 — while 55! has **74**. Forming the factorial first would
    /// overflow by thirty-five orders of magnitude, so this case only computes
    /// at all because the divisions are interleaved with the multiplications.
    #[test]
    fn dimension_survives_past_the_factorial_ceiling() {
        let lam = part(&[10, 9, 8, 7, 6, 5, 4, 3, 2, 1]);
        assert_eq!(lam.size(), 55);
        assert_eq!(
            dimension(&lam),
            Some(44_261_486_084_874_072_183_645_699_204_710_400)
        );
    }
}
