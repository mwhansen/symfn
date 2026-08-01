//! Labelled Dyck paths, and the combinatorial side of the Delta conjecture.
//!
//! ## Reference
//!
//! J. Haglund, J. Remmel, A. Wilson, *The Delta Conjecture*,
//! [arXiv:1509.07058](https://arxiv.org/abs/1509.07058), cited as **[HRW]**.
//! Their Conjecture 1.1 has two combinatorial sides for the same symmetric
//! function:
//!
//! ```text
//!   Δ'_{e_k} e_n  =  Rise_{n,k}(x;q,t)  =  Valley_{n,k}(x;q,t)
//! ```
//!
//! The **rise** version is a theorem (D'Adderio–Mellit; Blasiak–Haiman–Morse–
//! Pun–Seelinger, by two independent routes). The **valley** version is still
//! open. [`deltaop`](crate::deltaop) computes the left-hand side; this module
//! computes the two right-hand sides, so the three can be held against each
//! other.
//!
//! That makes this module the point of the whole Macdonald-operator exercise
//! rather than a test fixture. `docs/record/dyck-paths.md` records the
//! measurement that redirected it: once `Δ'_{e_k} e_n` costs 0.1s at degree 8,
//! the operator is no longer what stops a search, and **this enumeration is**.
//!
//! ## The objects
//!
//! A Dyck path of size `n` is its **area sequence** `a_1 … a_n` with `a_1 = 0`
//! and `0 ≤ a_i ≤ a_{i−1} + 1`. A *labelled* Dyck path attaches a positive
//! integer `ℓ_i` to row `i`, strictly increasing up a rise:
//!
//! ```text
//!   a_i = a_{i−1} + 1   ⟹   ℓ_i > ℓ_{i−1}
//! ```
//!
//! and then, [HRW]:
//!
//! ```text
//!   area(P) = Σ_i a_i
//!   d_i(P)  = #{ j > i : a_i = a_j,     ℓ_i < ℓ_j }
//!           + #{ j > i : a_i = a_j + 1, ℓ_i > ℓ_j }
//!   dinv(P) = Σ_i d_i(P)
//!   Rise(P) = { i : a_i = a_{i−1} + 1 }
//!   Val(P)  = { i : a_i < a_{i−1} } ∪ { i : a_i = a_{i−1}, ℓ_i > ℓ_{i−1} }
//! ```
//!
//! `Val` is the definition that has to be read carefully: it is **not**
//! `a_i ≤ a_{i−1}`. The ties are included only when the label goes *up*, and
//! taking the easy reading instead makes the valley side disagree with the rise
//! side everywhere except `k = n−1`, where nothing is chosen and both collapse
//! to the shuffle theorem. That was the first thing tried here and it looked
//! like a counterexample to an open conjecture for about a minute.
//!
//! ## Getting the whole symmetric function
//!
//! `x^P = ∏_i x_{ℓ_i}`, so the coefficient of `m_μ` is the number-with-weights
//! of labelled paths whose labels have **content μ** — `⟨f, h_μ⟩`, since `h` and
//! `m` are dual. So the monomial expansion is one enumeration per partition of
//! `n`, and `μ = (1ⁿ)` (all labels distinct) is the `⟨·, h_1ⁿ⟩` coefficient on
//! its own — the cheapest useful check, and the one that counts
//! `(n+1)^{n−1}` paths at `k = n−1`.
//!
//! ## The `z` extraction is an elementary symmetric polynomial
//!
//! Both sides read off `z^{n−k−1}` from a product over the chosen statistic:
//!
//! ```text
//!   Rise:   ∏_{i ∈ Rise(P)} (1 + z / t^{a_i})
//!   Valley: ∏_{i ∈ Val(P)}  (1 + z / q^{d_i(P)+1})
//! ```
//!
//! which is `e_{n−k−1}` of those weights. Enumerating the subsets is
//! unnecessary — a two-line knapsack over the exponents gives the whole
//! selection polynomial at once, and it is what keeps the inner loop cheap.
//!
//! ## The rise side does not need the labels at all
//!
//! `Rise(P)` and the weights `t^{−a_i}` it selects over are functions of the
//! **area sequence alone** — no label appears in either. So the whole
//! `z`-extraction is a constant of the labelling sum and factors straight out of
//! it:
//!
//! ```text
//!   Rise_{n,k} = Σ_D [ Σ_{S ⊆ Rise(D), |S| = n−1−k} t^{area(D) − Σ_{i∈S} a_i} ] · G_D(x; q)
//! ```
//!
//! where `G_D(x;q) = Σ_labellings q^{dinv} x^ℓ` is the **vertical-strip LLT
//! polynomial** of the path. That is the whole rise ladder from `C_n` LLT
//! evaluations plus one knapsack each, in place of one labelled-path walk per
//! content — and [`crate::llt`] computes `G_D` from `#SYT` standard objects
//! rather than `#labellings`, by [HHL]'s standardization. Measured against the
//! labelled walk, identical at every `k` and **29× / 56×** faster at `n = 8, 9`
//! (~2× per degree); `docs/record/llt.md` has the table.
//!
//! `docs/record/dyck-paths.md` recorded this as the win left on the table and
//! named the obstruction: standardizing *labelled paths* has no
//! `dinv`-invariant tie-break. The way through is that it is the *tuple*
//! fillings that get standardized, where [HHL] (82) is an identity rather than
//! a convention.
//!
//! **The valley side does not factor, and that is the whole point of it.**
//! `Val(P)` reads the labels — its tie clause is `ℓ_i > ℓ_{i−1}` — and its
//! weights are `q^{d_i+1}`, per labelling. So [`Side::Valley`] keeps the honest
//! enumeration below, and since valley is the *open* side, this makes the rise
//! half of the comparison free rather than moving the conjecture.
//! [`ladder_at_content`] keeps the labelled walk for **both** sides: it is the
//! oracle the fast route is checked against, and for a coarse content like
//! `μ = (n)` — one labelling — it is also simply cheaper.
//!
//! ## Negative exponents, and the one offset
//!
//! The rise weights divide by `t^{a_i}` and the chosen exponents sum to at most
//! `area(P)`, so `t` stays non-negative. The valley weights divide by
//! `q^{d_i+1}` and `Σ_{i∈S}(d_i + 1)` can exceed `dinv(P)` — by at most `|S| ≤ n`
//! — so the valley side genuinely visits negative powers of `q` before the sum
//! closes up. [`QtPoly`] has unsigned exponents, so the accumulation carries a
//! uniform `q^n` and divides it out at the end; the division is exact, and a
//! remainder there would mean the offset was too small rather than that the
//! conjecture failed.

// Path indices, bounded by the path length.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use crate::coeff::Ring;
use crate::partition::Partition;
use crate::qt::QtPoly;
use crate::sym::{Monomial, SymFn};

/// Every area sequence of size `n`.
///
/// `a_1 = 0` and `a_i ≤ a_{i−1} + 1`, which is the standard bijection with Dyck
/// paths; there are `C_n` of them.
fn for_each_area(n: usize, visit: &mut impl FnMut(&[u32])) {
    let mut area = Vec::with_capacity(n);
    fn rec(n: usize, area: &mut Vec<u32>, visit: &mut impl FnMut(&[u32])) {
        if area.len() == n {
            visit(area);
            return;
        }
        let top = if area.is_empty() {
            0
        } else {
            area[area.len() - 1] + 1
        };
        for a in 0..=top {
            area.push(a);
            rec(n, area, visit);
            area.pop();
        }
    }
    rec(n, &mut area, visit);
    debug_assert!(area.is_empty());
}

/// Every labelling of `area` whose labels have content `counts`
/// (`counts[j]` copies of the label `j+1`), strictly increasing up a rise.
fn for_each_labelling(
    area: &[u32],
    counts: &mut [u32],
    labels: &mut Vec<u32>,
    visit: &mut impl FnMut(&[u32]),
) {
    let i = labels.len();
    if i == area.len() {
        visit(labels);
        return;
    }
    let rise = i > 0 && area[i] == area[i - 1] + 1;
    for j in 0..counts.len() {
        if counts[j] == 0 {
            continue;
        }
        let label = j as u32 + 1;
        if rise && label <= labels[i - 1] {
            continue;
        }
        counts[j] -= 1;
        labels.push(label);
        for_each_labelling(area, counts, labels, visit);
        labels.pop();
        counts[j] += 1;
    }
}

/// `d_i(P)`, the diagonal inversions beginning in row `i`.
fn d_row(area: &[u32], labels: &[u32], i: usize) -> u32 {
    let mut d = 0;
    for j in (i + 1)..area.len() {
        if (area[i] == area[j] && labels[i] < labels[j])
            || (area[i] == area[j] + 1 && labels[i] > labels[j])
        {
            d += 1;
        }
    }
    d
}

/// `e_j` of the weights `x^{-w}` for **every** `j ≤ top`, as
/// `j ↦ [(subtracted exponent, multiplicity)]`.
///
/// A knapsack rather than a walk over the `C(|w|, j)` subsets: `dp[j][s]` counts
/// the ways to choose `j` of the weights with exponent sum `s`.
///
/// Every `j` at once, because the enumeration around this is the expensive part
/// and the `j`s are what the `k` ladder ranges over. Asking for one `j` per pass
/// — which is what this did first — re-walks every labelled path `n` times to
/// produce `n` slices of a table that one walk already fills.
fn choose_all(weights: &[u32], top: usize) -> Vec<Vec<(u32, i128)>> {
    let top = top.min(weights.len());
    let span: u32 = weights.iter().sum();
    let width = span as usize + 1;
    let mut dp = vec![0i128; (top + 1) * width];
    dp[0] = 1;
    for &w in weights {
        for j in (0..top).rev() {
            for s in (0..width).rev() {
                let c = dp[j * width + s];
                if c != 0 {
                    dp[(j + 1) * width + s + w as usize] += c;
                }
            }
        }
    }
    (0..=top)
        .map(|j| {
            (0..width)
                .filter_map(|s| {
                    let c = dp[j * width + s];
                    (c != 0).then_some((s as u32, c))
                })
                .collect()
        })
        .collect()
}

/// Which combinatorial side to build.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Side {
    /// [HRW]'s rise version — a **theorem**, so a mismatch is our bug.
    Rise,
    /// [HRW]'s valley version — **open**, so a mismatch is a result and must be
    /// reported as one rather than debugged away.
    Valley,
}

/// `Rise_{n,k}` or `Valley_{n,k}` in the monomial basis, for **every** `k` at
/// once: the returned vector is indexed by `k`, `0 ≤ k < n`.
///
/// This is the unit of work, and the two sides reach it differently.
///
/// [`Side::Rise`] goes through the per-path LLT polynomials
/// ([`rise_ladder_via_llt`]) — the factorization in the module docs, which
/// replaces the labelled-path walk entirely. [`Side::Valley`] cannot factor and
/// so enumerates: one walk per content, with [`choose_all`] producing every
/// `k`'s slice from one knapsack, so asking for a single `k` costs the same as
/// asking for all of them.
///
/// Both routes are held against each other by
/// `the_rise_ladder_via_llt_agrees_with_the_labelled_walk`.
pub fn ladder<C: Ring>(n: u32, which: Side) -> Vec<Monomial<QtPoly<C>>> {
    if n == 0 {
        return Vec::new();
    }
    if which == Side::Rise {
        return rise_ladder_via_llt(n);
    }
    let mut out = vec![Monomial::zero(); n as usize];
    for mu in crate::partitions_of(n) {
        for (k, c) in ladder_at_content::<C>(&mu, which).into_iter().enumerate() {
            if !c.is_empty() {
                out[k].add_term(mu.clone(), c);
            }
        }
    }
    out
}

/// The rise ladder from the per-path LLT decomposition — see the module docs.
///
/// One [`llt::llt_g`](crate::llt::llt_g) per area sequence gives `G_D` in the
/// monomial basis for **all** contents at once, where the labelled walk pays one
/// enumeration per content; the `z`-extraction is then the same knapsack over
/// the rise weights, which depend only on the area sequence.
///
/// The `t` exponent `area(D) − Σ_{i∈S} a_i` is non-negative because `S` is a
/// subset of the rows and the weights are those rows' own `a_i`, which is the
/// same reason [`ladder_at_content`] needs no offset on this side.
fn rise_ladder_via_llt<C: Ring>(n: u32) -> Vec<Monomial<QtPoly<C>>> {
    let nn = n as usize;
    let mut out = vec![Monomial::zero(); nn];
    for_each_area(nn, &mut |area| {
        let g: Monomial<QtPoly<C>> = crate::llt::llt_g(&crate::llt::SkewTuple::from_area(area));
        let areasum: u32 = area.iter().sum();
        let rises: Vec<u32> = (1..nn)
            .filter(|&i| area[i] == area[i - 1] + 1)
            .map(|i| area[i])
            .collect();
        // Slot j of the knapsack is the `z^j` coefficient, i.e. k = n−1−j.
        for (j, row) in choose_all(&rises, nn - 1).into_iter().enumerate() {
            let k = nn - 1 - j;
            for (sub, mult) in row {
                let scale = C::from_i128(mult);
                for (mu, poly) in g.terms() {
                    out[k]
                        .terms_mut()
                        .entry(mu.clone())
                        .or_insert_with(<QtPoly<C> as Ring>::zero)
                        .add_scaled_shifted(poly, (0, areasum - sub), &scale);
                }
            }
        }
    });
    // The slots were filled through the raw entry API, which bypasses
    // `add_term`'s zero check, so an explicit zero can survive a cancellation.
    for slot in &mut out {
        slot.terms_mut().retain(|_, c| !c.is_zero());
    }
    out
}

/// `⟨Rise_{n,k}, h_μ⟩` or `⟨Valley_{n,k}, h_μ⟩` for every `k`, indexed by `k`.
///
/// `μ = (1ⁿ)` is the cheap slice — the coefficient of `x_1 ⋯ x_n`, which at
/// `k = n−1` counts the `(n+1)^{n−1}` parking functions.
///
/// # Panics
///
/// If `μ` is empty. The Delta conjecture is stated for `n > 0`, and the
/// returned vector is indexed by `k ∈ 0..n`, so `n = 0` has no answer to hold.
pub fn ladder_at_content<C: Ring>(mu: &Partition, which: Side) -> Vec<QtPoly<C>> {
    let n = mu.size() as usize;
    assert!(n > 0, "the Delta conjecture asks for n > 0");
    // The valley weights can subtract more powers of q than `dinv` supplies; see
    // the module docs. `n` is enough: at most `|S| ≤ n` excess.
    let offset = match which {
        Side::Rise => 0,
        Side::Valley => n as u32,
    };

    let mut counts: Vec<u32> = mu.parts().to_vec();
    let mut acc: Vec<QtPoly<C>> = vec![QtPoly::zero(); n];
    let mut labels: Vec<u32> = Vec::with_capacity(n);
    let mut d: Vec<u32> = vec![0; n];

    for_each_area(n, &mut |area| {
        let areasum: u32 = area.iter().sum();
        for_each_labelling(area, &mut counts, &mut labels, &mut |labels| {
            for (i, slot) in d.iter_mut().enumerate() {
                *slot = d_row(area, labels, i);
            }
            let dinv: u32 = d.iter().sum();

            let weights: Vec<u32> = match which {
                Side::Rise => (1..n)
                    .filter(|&i| area[i] == area[i - 1] + 1)
                    .map(|i| area[i])
                    .collect(),
                Side::Valley => (1..n)
                    .filter(|&i| {
                        area[i] < area[i - 1]
                            || (area[i] == area[i - 1] && labels[i] > labels[i - 1])
                    })
                    .map(|i| d[i] + 1)
                    .collect(),
            };

            // slot j of the knapsack is the `z^{j}` coefficient, i.e. k = n−1−j.
            for (j, row) in choose_all(&weights, n - 1).into_iter().enumerate() {
                let k = n - 1 - j;
                for (sub, mult) in row {
                    let (a, b) = match which {
                        Side::Rise => (dinv, areasum - sub),
                        Side::Valley => (dinv + offset - sub, areasum),
                    };
                    acc[k].add_term(a, b, C::from_i128(mult));
                }
            }
        });
    });

    if offset == 0 {
        return acc;
    }
    // Divide the uniform q^offset back out. Exact by construction — a failure
    // here means the offset was too small, not that the conjecture broke.
    let divisor: QtPoly<C> = QtPoly::term(offset, 0, C::one());
    acc.into_iter()
        .map(|f| {
            f.divide_exact(&divisor)
                .expect("the q-offset must divide out exactly")
        })
        .collect()
}

/// `Rise_{n,k}` or `Valley_{n,k}` for one `k`.
///
/// A slice of [`ladder`], and it costs the same as the whole ladder — see there.
///
/// # Panics
///
/// Unless `k < n`, which is the range the Delta conjecture is stated over.
pub fn side<C: Ring>(n: u32, k: u32, which: Side) -> Monomial<QtPoly<C>> {
    assert!(k < n, "the Delta conjecture asks for k < n");
    ladder::<C>(n, which).swap_remove(k as usize)
}

/// One coefficient of one `k`. A slice of [`ladder_at_content`], at the same
/// cost.
///
/// # Panics
///
/// Unless `k < |μ|`, and if `μ` is empty; see [`ladder_at_content`].
pub fn side_at_content<C: Ring>(mu: &Partition, k: u32, which: Side) -> QtPoly<C> {
    assert!(k < mu.size(), "the Delta conjecture asks for k < n");
    ladder_at_content::<C>(mu, which).swap_remove(k as usize)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::convert::ToSchur;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    /// The count of labelled Dyck paths with distinct labels is `(n+1)^{n−1}` —
    /// the dimension of the diagonal harmonics.
    ///
    /// This is the enumeration under test on its own, before any statistic is
    /// involved: `k = n−1` chooses nothing, so every path contributes 1.
    #[test]
    fn distinct_labellings_count_the_diagonal_harmonics() {
        for n in 1..=7u32 {
            let f: QtPoly<Rational> =
                side_at_content(&part(&vec![1; n as usize]), n - 1, Side::Rise);
            let total: i128 = f.terms().map(|(_, c)| c.numer()).sum();
            assert_eq!(
                total,
                (n as i128 + 1).pow(n - 1),
                "labelled Dyck paths of size {n}"
            );
        }
    }

    /// **The rise version, against the operator.** A theorem, so a mismatch here
    /// is a bug on one side or the other.
    ///
    /// Compared as whole symmetric functions, not just at `h_1ⁿ`: the monomial
    /// expansion is one enumeration per content, and a statistic that was wrong
    /// only on repeated labels would survive the `(1ⁿ)` check alone.
    #[test]
    fn the_rise_side_is_delta_prime() {
        for n in 1..=6u32 {
            for k in 0..n {
                let got = side::<Rational>(n, k, Side::Rise).to_schur();
                let want = crate::deltaop::delta_prime_e::<Rational>(k, n);
                assert_eq!(got, want, "Rise_{{{n},{k}}}");
            }
        }
    }

    /// **The valley version, against the operator — an open conjecture.**
    ///
    /// If this ever fails, the first thing to check is `Val(P)`: the tie case
    /// belongs in only when the label rises (see the module docs). If it fails
    /// with `Val` right, it is a counterexample and must be reported as one.
    #[test]
    fn the_valley_side_is_delta_prime() {
        for n in 1..=6u32 {
            for k in 0..n {
                let got = side::<Rational>(n, k, Side::Valley).to_schur();
                let want = crate::deltaop::delta_prime_e::<Rational>(k, n);
                assert_eq!(got, want, "Valley_{{{n},{k}}}");
            }
        }
    }

    /// The two sides must agree with each other, which is the conjecture stated
    /// without the operator in the way.
    #[test]
    fn the_two_sides_agree() {
        for n in 1..=6u32 {
            for k in 0..n {
                assert_eq!(
                    side::<Rational>(n, k, Side::Rise),
                    side::<Rational>(n, k, Side::Valley),
                    "the two sides at n={n}, k={k}"
                );
            }
        }
    }

    /// `k = n−1` is the shuffle theorem: both sides must be `∇e_n`.
    #[test]
    fn the_top_of_the_ladder_is_the_shuffle_theorem() {
        for n in 1..=6u32 {
            let want = crate::deltaop::nabla_e::<Rational>(n);
            for which in [Side::Rise, Side::Valley] {
                assert_eq!(
                    side::<Rational>(n, n - 1, which).to_schur(),
                    want,
                    "{which:?} at k = n-1, n={n}"
                );
            }
        }
    }

    /// **The rise route against the labelled walk it replaced.**
    ///
    /// [`ladder`] no longer enumerates labelled paths on the rise side, so the
    /// walk that used to be the implementation is now the oracle — and it has to
    /// be checked at every `k` and every content, not just in total: the
    /// factorization moves the `z`-extraction outside the labelling sum, and an
    /// error there would show up as a redistribution between `k`s that any
    /// aggregate check would miss.
    #[test]
    fn the_rise_ladder_via_llt_agrees_with_the_labelled_walk() {
        for n in 1..=6u32 {
            let fast = rise_ladder_via_llt::<Rational>(n);
            // The labelled walk, one content at a time — the pre-LLT route.
            let mut slow = vec![Monomial::zero(); n as usize];
            for mu in crate::partitions_of(n) {
                for (k, c) in ladder_at_content::<Rational>(&mu, Side::Rise)
                    .into_iter()
                    .enumerate()
                {
                    if !c.is_empty() {
                        slow[k].add_term(mu.clone(), c);
                    }
                }
            }
            for k in 0..n as usize {
                assert_eq!(fast[k], slow[k], "Rise_{{{n},{k}}} via LLT vs labellings");
            }
        }
    }

    /// The subset knapsack must agree with a direct walk over the subsets.
    #[test]
    fn the_selection_knapsack_agrees_with_enumerating_subsets() {
        let weights = [0u32, 1, 2, 2, 3, 1];
        for m in 0..=weights.len() {
            let mut want: std::collections::BTreeMap<u32, i128> = Default::default();
            for mask in 0u32..(1 << weights.len()) {
                if mask.count_ones() as usize != m {
                    continue;
                }
                let s: u32 = (0..weights.len())
                    .filter(|&i| mask >> i & 1 == 1)
                    .map(|i| weights[i])
                    .sum();
                *want.entry(s).or_insert(0) += 1;
            }
            let rows = choose_all(&weights, weights.len());
            let got: std::collections::BTreeMap<u32, i128> = rows[m].iter().copied().collect();
            assert_eq!(got, want, "choosing {m}");
        }
    }
}
