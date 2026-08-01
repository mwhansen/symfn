//! Hall–Littlewood polynomials, by the Morris (1963) recursion.
//!
//! Returns `Q'_λ(x; t)` expanded in the Schur basis with coefficients in
//! [`QtPoly`] — that is, `Q'_λ = Σ_μ K_{μλ}(t) s_μ`, so the Kostka–Foulkes
//! polynomials fall out of the expansion rather than being computed separately.
//!
//! ## The recursion
//!
//! Peel the **largest** part, recurse, then rebuild:
//!
//! ```text
//!   HL(∅)  = s_∅
//!   HL(λ)  = Σ_{i=0}^{|λ⁻|} t^i · straighten( (h_i^⊥ HL(λ⁻)) with a part λ_1 + i )
//! ```
//!
//! where λ⁻ is λ with its largest part removed. Each term of `h_i^⊥ HL(λ⁻)` has
//! size `|λ⁻| − i`, and the new part contributes `λ_1 + i`, so every branch lands
//! back on degree `|λ|` — the loop bound `i ≤ |λ⁻|` is exactly where `h_i^⊥`
//! starts annihilating.
//!
//! Two pieces of this were already built and already fast: `h_i^⊥` on a Schur
//! expansion is [`strip_off`](crate::hopf) — Pieri run backwards, horizontal
//! strips, no Littlewood–Richardson — and the coefficient ring is [`QtPoly`].
//!
//! ## Straightening
//!
//! The new part is inserted at the **top** row, where it need not be the largest,
//! so the result is a Schur function of a non-partition sequence and has to be
//! straightened. This is β-numbers: `β_j = c_j + j` on the sequence written
//! *ascending*, then sort. Equal β's mean the term is zero; each transposition
//! flips the sign; a negative part after sorting also kills the term.
//!
//! ## Provenance
//!
//! The recursion is Symmetrica's `hall_littlewood` in `sr.c` (public domain, see
//! `NOTICE.md`), which follows A. O. Morris, *The characters of the group
//! GL(n,q)*, Math. Zeitschr. **81** (1963) 112–123. Symmetrica works with
//! partitions stored ascending and appends the new part at the end of that
//! vector; here partitions are descending, so the same step prepends. Its
//! `reorder_hall_littlewood` performs the straightening by repeated adjacent
//! swaps; sorting the β-numbers once is the same map, and the sign is the parity
//! of the sort.
//!
//! Symmetrica has no Kostka–Foulkes entry point and its `hall_littlewood` uses
//! no charge statistic, so [`crate::charge::kostka_foulkes_by_charge`] is a genuinely
//! independent check on everything here.
//!
//! ## Reach
//!
//! Over a fixed-width `C` this family is exact until a coefficient leaves the
//! width, and then it refuses rather than wrapping (`docs/policies/failure.md`,
//! R3). Through the Python boundary [`hall_littlewood`] escalates — the
//! fixed-width pass reports and the same generic code re-runs over `BigInt` —
//! so the wall below is what a *Rust* caller at `C = i128` meets.
//! **The two entry points have different reach, and the difference is the
//! point:**
//!
//! * [`hall_littlewood_table`] is p(n) polynomials, and stops finishing long
//!   before it stops fitting: at `i128` its widest coefficient gains ~2.2 bits
//!   per degree and would reach 127 bits near n ≈ 67, which no machine reaches.
//! * [`hall_littlewood`] on one λ reaches much further in degree, so the wall
//!   is real: at λ = 1ⁿ — the worst case, see below — coefficients are 87 bits
//!   at n = 47 and project to 127 bits near n ≈ 63, which is hours rather than
//!   never.
//!
//! **The shape matters more than the degree.** `Q'_{1ⁿ}` is extremal because
//! its coefficients are the Kostka–Foulkes ones at content `1ⁿ`, which sum to
//! `K_{μ,1ⁿ} = f^μ`, and `Σ_μ (f^μ)² = n!` caps every coefficient at `√(n!)` —
//! the same bound that puts the character ceiling near n ≈ 58, and the measured
//! widths track it about 9 bits below. At the other extreme `Q'_{(n-1,1)}` has
//! 1-bit coefficients at *every* degree tested and no wall at all. A reach
//! statement in n alone is therefore wrong for one of them.
//!
//! [`hall_littlewood_p`] has no such bound — `P` comes from inverting the
//! Kostka–Foulkes matrix, so its coefficients are signed and are not counts —
//! and is measured only: 14 bits at n = 24 for the table, ~0.9 bits per degree.
//!
//! Degrees, slopes and the harness are in `docs/record/failure-and-overflow.md`
//! (`examples/probe_qt_walls.rs`).

// β-numbers and part indices, bounded by |λ|; coefficients are `QtPoly<C>`.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::HashMap;
use std::rc::Rc;

use crate::coeff::Ring;
use crate::partition::Partition;
use crate::qt::QtPoly;
use crate::sym::{Schur, SymFn};

/// `Q'_λ(x; t)` in the Schur basis: `Σ_μ K_{μλ}(t) s_μ`.
pub fn hall_littlewood<C: Ring>(lambda: &Partition) -> Schur<QtPoly<C>> {
    let mut memo = HashMap::new();
    Rc::try_unwrap(hl_suffix(lambda.parts(), &mut memo)).unwrap_or_else(|rc| (*rc).clone())
}

/// `Q'_λ` for **every** λ ⊢ n, in the order of [`partitions_of`](crate::partitions_of).
///
/// The recursion peels from the front, so its subproblems are the *suffixes* of
/// λ, and partitions of n share those heavily — this computes each distinct
/// suffix once. Same pattern as [`kostka_table`](crate::kostka) and the
/// character table: a table is cheaper than p(n) separate calls.
pub fn hall_littlewood_table<C: Ring>(n: u32) -> Vec<(Partition, Schur<QtPoly<C>>)> {
    let mut memo = HashMap::new();
    crate::partitions_of(n)
        .into_iter()
        .map(|lambda| {
            let v = hl_suffix(lambda.parts(), &mut memo);
            (lambda, (*v).clone())
        })
        .collect()
}

/// `P_λ(x; t)` in the Schur basis, for every λ ⊢ n.
///
/// The *other* Hall–Littlewood basis. [`hall_littlewood`] gives `Q'`, and the
/// two are related by
///
/// ```text
///   s_μ = Σ_λ K_{μλ}(t) P_λ        while       Q'_λ = Σ_μ K_{μλ}(t) s_μ
/// ```
///
/// — the same Kostka–Foulkes matrix, used in the two directions. So `P` is what
/// comes out of **inverting** it, and no new enumeration is needed.
///
/// The inverse stays in ℤ[t]: `K` is unitriangular in dominance order, so
/// solving `s_μ = P_μ + Σ_{λ ◁ μ} K_{μλ} P_λ` for `P_μ` never divides. The
/// partitions are visited lex-ascending, which is a linear extension of
/// dominance (λ ⊵ μ implies λ ≥ μ lexicographically), so every `P_λ` the sum
/// needs is already known.
pub fn hall_littlewood_p_table<C: Ring>(n: u32) -> Vec<(Partition, Schur<QtPoly<C>>)> {
    let parts = crate::memo::partitions_cached(n);
    let k = crate::kf::kostka_foulkes_table::<C>(n);
    let mut out: Vec<Schur<QtPoly<C>>> = vec![Schur::zero(); parts.len()];
    // `partitions_of` is lex-descending, so counting down visits lex-ascending.
    for j in (0..parts.len()).rev() {
        let mut acc = Schur::monomial(parts[j].clone(), <QtPoly<C> as Ring>::one());
        for l in (j + 1)..parts.len() {
            let coeff = &k[j][l];
            if coeff.is_zero() {
                continue;
            }
            let done = std::mem::replace(&mut out[l], Schur::zero());
            for (nu, c) in done.terms() {
                acc.add_term(nu.clone(), coeff.mul(c).neg());
            }
            out[l] = done;
        }
        out[j] = acc;
    }
    parts.iter().cloned().zip(out).collect()
}

/// `P_λ(x; t)` in the Schur basis.
///
/// Computes the whole degree: the inversion needs every dominance-smaller `P`
/// anyway, so a single shape costs what the table costs — use
/// [`hall_littlewood_p_table`] when more than one is wanted.
pub fn hall_littlewood_p<C: Ring>(lambda: &Partition) -> Schur<QtPoly<C>> {
    hall_littlewood_p_table(lambda.size())
        .into_iter()
        .find(|(mu, _)| mu == lambda)
        .map(|(_, f)| f)
        .unwrap_or_else(|| Schur::monomial(Partition::new([]), <QtPoly<C> as Ring>::one()))
}

type Memo<C> = HashMap<Vec<u32>, Rc<Schur<QtPoly<C>>>>;

/// `HL` of a descending part list, memoised on the list itself.
fn hl_suffix<C: Ring>(parts: &[u32], memo: &mut Memo<C>) -> Rc<Schur<QtPoly<C>>> {
    if parts.is_empty() {
        return Rc::new(Schur::monomial(
            Partition::new([]),
            <QtPoly<C> as Ring>::one(),
        ));
    }
    if let Some(hit) = memo.get(parts) {
        return Rc::clone(hit);
    }

    let first = parts[0] as i64;
    let prev = hl_suffix(&parts[1..], memo);

    // Every (ν, i) at once.
    //
    // Written as the recursion reads — `for i in 0..=|λ⁻| { h_i^⊥ prev }` — this
    // built a whole intermediate `Schur<QtPoly>` per i, walked it once, and
    // dropped it. A sampling profile put 660 samples in map insertion, ~750 in
    // malloc/free and 44 in the strip enumeration itself: the cost was the
    // temporaries, not the combinatorics.
    //
    // Removing a horizontal strip of *any* size from ν is one interlacing walk
    // (ν_1 ≥ μ_1 ≥ ν_2 ≥ μ_2 ≥ …) with `i = |ν| − |μ|` falling out at the leaf,
    // so the sizes come for free from a single pass and nothing is materialised
    // between `prev` and `out`.
    let mut out = Schur::zero();
    let mut v: Vec<i64> = Vec::with_capacity(parts.len() + 1);
    // `straighten` works in place, so it gets its own buffer: `v` is the live
    // stack of the strip walk and must survive the call unchanged.
    let mut beta: Vec<i64> = Vec::with_capacity(parts.len() + 1);
    for (nu, c) in prev.terms() {
        let nu = nu.parts();
        removals(nu, nu.len(), &mut v, 0, &mut |mu, i| {
            beta.clear();
            beta.extend_from_slice(mu);
            beta.push(first + i as i64);
            if let Some((negate, mu)) = straighten(&mut beta) {
                // `out.add_term(mu, c.shift_t(i))` reads better and was 464 of
                // ~1800 profile samples: `shift_t` allocates a whole polynomial
                // so that `add_assign` can immediately walk it and drop it. The
                // shift is a rename of the exponents, so the terms are merged
                // straight into the slot instead.
                out.terms_mut()
                    .entry(mu)
                    .or_insert_with(<QtPoly<C> as Ring>::zero)
                    .add_shifted(c, (0, i), negate);
            }
        });
    }

    // Inserting into the slot directly bypasses `add_term`'s zero check, so a
    // coefficient that cancelled to zero can survive as an explicit zero. Every
    // consumer — equality, `terms().len()`, the Python boundary — assumes the
    // map holds no zeros, so prune once here rather than checking per insert.
    out.terms_mut().retain(|_, c| !c.is_zero());

    let out = Rc::new(out);
    memo.insert(parts.to_vec(), Rc::clone(&out));
    out
}

/// Every μ with ν/μ a horizontal strip, of every size at once.
///
/// Emits μ **ascending** — rows are chosen bottom-up, which is already the order
/// [`straighten`] wants — together with the number of cells removed. Zero parts
/// are kept: in ascending form they lead, and a leading zero shifts every
/// β-number and every index by the same amount, so it changes neither the sort,
/// the collisions, nor the sign.
fn removals(
    nu: &[u32],
    row: usize,
    cur: &mut Vec<i64>,
    removed: u32,
    emit: &mut impl FnMut(&mut Vec<i64>, u32),
) {
    if row == 0 {
        emit(cur, removed);
        return;
    }
    let r = row - 1;
    let floor = if r + 1 < nu.len() { nu[r + 1] } else { 0 };
    for m in floor..=nu[r] {
        cur.push(m as i64);
        removals(nu, r, cur, removed + (nu[r] - m), emit);
        cur.pop();
    }
}

/// Sort an **ascending** integer sequence into a partition by β-numbers.
///
/// Returns `(negate, μ)`, or `None` when the term vanishes. `v` is scratch and
/// is left in an unspecified state.
fn straighten(v: &mut [i64]) -> Option<(bool, Partition)> {
    let l = v.len();
    // β_j = v_j + j. Ascending v with no forbidden pattern gives strictly
    // increasing β, so sorting β is the whole of the straightening rule.
    for (j, x) in v.iter_mut().enumerate() {
        *x += j as i64;
    }
    // Insertion sort: the sequence is nearly sorted (only the pushed part is out
    // of place), so this is a short walk, and it counts its own transpositions.
    let mut swaps = 0usize;
    for i in 1..l {
        let mut j = i;
        while j > 0 && v[j - 1] > v[j] {
            v.swap(j - 1, j);
            swaps += 1;
            j -= 1;
        }
    }
    for i in 1..l {
        if v[i] == v[i - 1] {
            return None; // repeated β-number: the Schur function is zero
        }
    }
    if l > 0 && v[0] < 0 {
        return None; // a negative part survives the sort
    }
    // Back to parts, then descending, dropping the zeros.
    let parts: Vec<u32> = (0..l)
        .rev()
        .map(|j| (v[j] - j as i64) as u32)
        .filter(|&x| x > 0)
        .collect();
    Some((swaps % 2 == 1, Partition::from_sorted(parts)))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::convert::ToSchur;
    use crate::sym::{Homogeneous, Monomial};

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    type Q = QtPoly<i64>;

    /// t = 0 collapses Q'_λ to s_λ — only the i = 0 branch survives.
    #[test]
    fn at_t_zero_it_is_the_schur_function() {
        for n in 0..=8u32 {
            for lambda in crate::partitions_of(n) {
                let hl: Schur<Q> = hall_littlewood(&lambda);
                for (mu, c) in hl.terms() {
                    let want = i64::from(*mu == lambda);
                    assert_eq!(c.eval(&0, &0), want, "Q'_{lambda} at t=0, term {mu}");
                }
            }
        }
    }

    /// t = 1 gives h_λ. This is the identity that pins *which* Hall–Littlewood
    /// this is: Q'_λ(x;1) = h_λ, while P_λ(x;1) = m_λ.
    #[test]
    fn at_t_one_it_is_the_complete_homogeneous() {
        for n in 0..=8u32 {
            for lambda in crate::partitions_of(n) {
                let hl: Schur<Q> = hall_littlewood(&lambda);
                let h: Schur<i64> = Homogeneous::monomial(lambda.clone(), 1).to_schur();
                for (mu, c) in hl.terms() {
                    assert_eq!(c.eval(&0, &1), h.coeff(mu), "Q'_{lambda} at t=1, {mu}");
                }
                assert_eq!(
                    hl.terms().len(),
                    h.terms().len(),
                    "same support at t = 1 for {lambda}"
                );
            }
        }
    }

    /// The coefficients *are* the Kostka–Foulkes polynomials — checked against
    /// the charge enumeration, which shares no code with this recursion.
    #[test]
    fn coefficients_are_kostka_foulkes() {
        for n in 0..=8u32 {
            let parts = crate::partitions_of(n);
            for lambda in &parts {
                let hl: Schur<Q> = hall_littlewood(lambda);
                for mu in &parts {
                    let want = crate::charge::kostka_foulkes_by_charge::<i64>(mu, lambda);
                    assert_eq!(hl.coeff(mu), want, "K_{mu}{lambda}(t)");
                }
            }
        }
    }

    /// Triangularity: K_{μλ}(t) = 0 unless μ ⊵ λ, and K_{λλ}(t) = 1.
    #[test]
    fn expansion_is_unitriangular() {
        for n in 0..=8u32 {
            for lambda in crate::partitions_of(n) {
                let hl: Schur<Q> = hall_littlewood(&lambda);
                assert_eq!(hl.coeff(&lambda), <Q as Ring>::one(), "K_{lambda}{lambda}");
                for mu in hl.terms().keys() {
                    assert!(
                        crate::kostka::kostka(mu, &lambda) > 0,
                        "{mu} appears in Q'_{lambda} but is outside the Kostka support"
                    );
                }
            }
        }
    }

    /// P is the other Hall–Littlewood basis, and the relation that defines it
    /// must hold against the Kostka–Foulkes matrix it was inverted from:
    /// `s_μ = Σ_λ K_{μλ}(t) P_λ`.
    #[test]
    fn p_inverts_the_kostka_foulkes_matrix() {
        for n in 0..=8u32 {
            let parts = crate::partitions_of(n);
            let table: Vec<(Partition, Schur<Q>)> = hall_littlewood_p_table(n);
            for mu in &parts {
                let mut acc: Schur<Q> = Schur::zero();
                for (lambda, p) in &table {
                    let k = crate::kostka_foulkes::<i64>(mu, lambda);
                    if k.is_zero() {
                        continue;
                    }
                    for (nu, c) in p.terms() {
                        acc.add_term(nu.clone(), k.mul(c));
                    }
                }
                let want: Schur<Q> = Schur::monomial(mu.clone(), <Q as Ring>::one());
                assert_eq!(acc, want, "s_{mu} = sum_lambda K P_lambda");
            }
        }
    }

    /// t = 0 gives s_λ, t = 1 gives m_λ — the two specialisations of P, and the
    /// pair that distinguishes it from Q' (which gives s_λ and h_λ).
    #[test]
    fn p_specialises_to_schur_and_monomial() {
        use crate::sym::Monomial;
        for n in 0..=8u32 {
            for lambda in crate::partitions_of(n) {
                let p: Schur<Q> = hall_littlewood_p(&lambda);
                for (mu, c) in p.terms() {
                    assert_eq!(
                        c.eval(&0, &0),
                        i64::from(*mu == lambda),
                        "P_{lambda} at t = 0 must be s_{lambda}, term {mu}"
                    );
                }
                let m: Schur<i64> = Monomial::monomial(lambda.clone(), 1).to_schur();
                for (mu, c) in p.terms() {
                    assert_eq!(c.eval(&0, &1), m.coeff(mu), "P_{lambda} at t = 1, {mu}");
                }
                // Both directions, but *not* equal supports: a P coefficient may
                // be a nonzero polynomial that happens to vanish at t = 1, since
                // unlike Q' its coefficients are not sign-definite. Comparing
                // `terms().len()` fails here for a correct answer.
                for (mu, want) in m.terms() {
                    assert_eq!(p.coeff(mu).eval(&0, &1), *want, "missing {mu} at t = 1");
                }
            }
        }
    }

    /// The table must agree term for term with the one-at-a-time calls; sharing
    /// suffixes is an optimisation, not a different algorithm.
    #[test]
    fn table_agrees_with_individual_calls() {
        for n in 0..=9u32 {
            for (lambda, got) in hall_littlewood_table::<i64>(n) {
                assert_eq!(got, hall_littlewood(&lambda), "table at {lambda}");
            }
        }
    }

    /// Straightening in isolation, on the cases the recursion actually hits.
    #[test]
    fn straightening_follows_the_beta_rules() {
        // already a partition, ascending (1,2) -> (2,1)
        assert_eq!(straighten(&mut [1, 2]), Some((false, part(&[2, 1]))));
        // s_{(1,2)} written descending as a descent: ascending (2,1) has
        // β = (2,2) -> repeated, zero
        assert_eq!(straighten(&mut [2, 1]), None);
        // ascending (0,3): β = (0,4) -> parts (0,3) -> (3)
        assert_eq!(straighten(&mut [0, 3]), Some((false, part(&[3]))));
        // ascending (3,1): β = (3,2) -> sort to (2,3), one swap -> parts (2,2)
        assert_eq!(straighten(&mut [3, 1]), Some((true, part(&[2, 2]))));
        // a negative part that survives: ascending (-1, 5) is already sorted in β
        assert_eq!(straighten(&mut [-1, 5]), None);
        assert_eq!(straighten(&mut []), Some((false, part(&[]))));
    }

    /// Q'_λ is not P_λ, and the two are easy to conflate. m_λ is what P would
    /// give at t = 1; check we do *not* give it, on the smallest case where the
    /// two differ.
    #[test]
    fn it_is_q_prime_not_p() {
        let hl: Schur<Q> = hall_littlewood(&part(&[2, 1]));
        let m: Schur<i64> = Monomial::monomial(part(&[2, 1]), 1).to_schur();
        let at_one: Vec<(Partition, i64)> = hl
            .terms()
            .iter()
            .map(|(mu, c)| (mu.clone(), c.eval(&0, &1)))
            .collect();
        let as_m: Vec<(Partition, i64)> =
            m.terms().iter().map(|(mu, c)| (mu.clone(), *c)).collect();
        assert_ne!(at_one, as_m, "Q'_(2,1)(x;1) must be h_(2,1), not m_(2,1)");
    }
}
