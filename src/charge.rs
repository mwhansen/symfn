//! The charge statistic, and Kostka–Foulkes polynomials as a **reference**.
//!
//! `K_{λμ}(t) = Σ_{T ∈ SSYT(λ, μ)} t^{charge(T)}` (Lascoux–Schützenberger).
//! This enumerates tableaux, so it is exponential and is not how the library
//! computes Hall–Littlewood — it is the obviously-correct implementation the
//! fast one gets held to, in the same role [`NaiveLr`](crate::NaiveLr) plays
//! for Littlewood–Richardson.
//!
//! Hall–Littlewood comes from a recursion over skewing and straightening (see
//! `docs/record/hall-littlewood.md`); this shares no code with that, so
//! agreement between them is evidence rather than tautology (V3,
//! `docs/policies/validation.md`). Symmetrica's
//! `hall_littlewood` uses no charge statistic, which is what makes the two
//! routes genuinely disjoint.
//!
//! ## Charge
//!
//! On a word whose content is a partition:
//!
//! * **standard** (every letter once) — `index(1) = 0`, and `index(i)` is
//!   `index(i−1) + 1` when `i` lies to the right of `i−1`, else `index(i−1)`.
//!   The charge is the sum of the indices.
//! * **general** — remove standard subwords by scanning right to left for a
//!   1, then continuing (wrapping at the left end) for a 2, and so on until no
//!   next letter is found; remove that subword and repeat. The charge is the
//!   sum over the subwords.
//!
//! The reading word of a tableau is taken **bottom row to top row, each row left
//! to right**. Conventions differ across the literature by reversal, and this is
//! the one that reproduces Sage's `KostkaFoulkesPolynomial`
//! (`sage.combinat.sf.kfpoly`): for λ = (2,1), μ = (1,1,1) the two tableaux
//! give charges 2 and 1, i.e. `t² + t`.
//!
//! \[M\] III.6 reads the tableau the other way — top row first, each row
//! right to left — which is this word reversed, and mirrors both rules to
//! match: the index rises when `i` lies to the *left* of `i−1`, and subwords
//! are extracted scanning from the left. Reversing the word and mirroring both
//! rules leaves every charge unchanged, so this module is \[M\]'s definition
//! read backwards. Reversing the word without mirroring the rules is the
//! plausible wrong answer.
//!
//! ## References
//!
//! - **\[M\]** I. G. Macdonald, *Symmetric Functions and Hall Polynomials*,
//!   2nd ed., Oxford, 1995 — III (6.5), the theorem, with the charge
//!   of a word and of a tableau defined in the three steps just before it.
//!   This is the normative source for the statistic.
//! - **\[LS\]** A. Lascoux, M.-P. Schützenberger, *Sur une conjecture de
//!   H. O. Foulkes*, C. R. Acad. Sci. Paris Sér. A-B **286** (1978) — the
//!   theorem this module implements, that `K_{λμ}(t)` is the charge generating
//!   function over `SSYT(λ, μ)`.
//! - **\[But\]** L. M. Butler, *Subgroup lattices and symmetric functions*,
//!   Mem. Amer. Math. Soc. **112** (1994), no. 539, Ch. 2.4 — a second
//!   account of the charge statistic and of that generating function.

// A word length.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use crate::coeff::Ring;
use crate::partition::Partition;
use crate::qt::QtPoly;

/// The charge of a word whose content is weakly decreasing.
///
/// Letters are 1-based. Returns 0 for the empty word.
pub fn charge(word: &[u32]) -> u32 {
    // The empty word and the `max` are one question, so they are asked once:
    // the `else` arm *is* the documented degenerate case.
    let Some(&top) = word.iter().max() else {
        return 0;
    };
    let mut counts = vec![0u32; top as usize + 1];
    for &c in word {
        counts[c as usize] += 1;
    }
    if counts[1..].iter().all(|&c| c <= 1) {
        return standard_charge(word);
    }
    // Peel standard subwords until nothing is left.
    let mut rest: Vec<Option<u32>> = word.iter().map(|&c| Some(c)).collect();
    let mut total = 0;
    while rest.iter().any(|c| c.is_some()) {
        let sub = peel_standard(&mut rest);
        if sub.is_empty() {
            break; // defensive: cannot happen for a partition content
        }
        total += standard_charge(&sub);
    }
    total
}

/// Charge of a word in which every letter occurs at most once.
fn standard_charge(word: &[u32]) -> u32 {
    let top = match word.iter().max() {
        Some(&m) => m as usize,
        None => return 0,
    };
    // position of each letter
    let mut at = vec![usize::MAX; top + 1];
    for (j, &c) in word.iter().enumerate() {
        at[c as usize] = j;
    }
    let mut index = 0u32;
    let mut total = 0u32;
    let mut prev = usize::MAX;
    for c in 1..=top {
        let pos = at[c];
        if pos == usize::MAX {
            continue;
        }
        if prev != usize::MAX && pos > prev {
            index += 1;
        }
        total += index;
        prev = pos;
    }
    total
}

/// Remove one standard subword, returning it in left-to-right order.
///
/// Scans right to left for a 1, then keeps going — wrapping past the left end —
/// for a 2, and so on. Wrapping is what makes the subword *standard* rather
/// than merely increasing, and is the step that is easy to omit.
fn peel_standard(rest: &mut [Option<u32>]) -> Vec<u32> {
    let n = rest.len();
    let mut picked: Vec<(usize, u32)> = Vec::new();
    let mut target = 1u32;
    let mut from = n; // scanning starts just past the right end
    loop {
        let mut found = None;
        for step in 0..n {
            let j = (from + n - 1 - step) % n;
            if rest[j] == Some(target) {
                found = Some(j);
                break;
            }
        }
        match found {
            Some(j) => {
                rest[j] = None;
                picked.push((j, target));
                target += 1;
                from = j;
            }
            None => break,
        }
    }
    picked.sort_unstable_by_key(|&(j, _)| j);
    picked.into_iter().map(|(_, c)| c).collect()
}

/// `K_{λμ}(t)`, by enumerating semistandard tableaux.
///
/// Returns the zero polynomial when |λ| ≠ |μ|, and `1` for the empty pair.
///
/// The **reference** implementation: exponential, and here to be disagreed with.
/// [`kostka_foulkes`](crate::kostka_foulkes) is the one to call — it reads the
/// polynomials off the Hall–Littlewood transition, shares no code with this, and
/// is what the `kf` tests hold to this one.
pub fn kostka_foulkes_by_charge<C: Ring>(lambda: &Partition, mu: &Partition) -> QtPoly<C> {
    let mut out = QtPoly::zero();
    if lambda.size() != mu.size() {
        return out;
    }
    if lambda.is_empty() {
        return <QtPoly<C> as Ring>::one();
    }
    let rows = lambda.len();
    let mut chain: Vec<Vec<u32>> = vec![vec![0; rows]];
    let target: Vec<u32> = lambda.parts().to_vec();
    build(&target, mu.parts(), 0, &mut chain, &mut |ch| {
        let w = reading_word(&target, ch);
        out.add_term(0, charge(&w), C::one());
    });
    out
}

/// Extend the chain by one horizontal strip per part of μ, in order.
///
/// Shared with [`crate::macdonald`], which needs the same chains of horizontal
/// strips — a semistandard tableau of shape λ and content μ *is* such a chain,
/// and both the charge statistic and Macdonald's ψ are functions of it.
pub(crate) fn build(
    target: &[u32],
    mu: &[u32],
    k: usize,
    chain: &mut Vec<Vec<u32>>,
    emit: &mut impl FnMut(&[Vec<u32>]),
) {
    if k == mu.len() {
        if chain
            .last()
            .expect("the chain starts at the empty shape")
            .as_slice()
            == target
        {
            emit(chain);
        }
        return;
    }
    let cur = chain
        .last()
        .expect("the chain starts at the empty shape")
        .clone();
    let mut next = cur.clone();
    strips(target, &cur, 0, u32::MAX, mu[k], &mut next, &mut |shape| {
        chain.push(shape.to_vec());
        build(target, mu, k + 1, chain, emit);
        chain.pop();
    });
}

/// Every ν with `cur ⊆ ν ⊆ target`, `ν/cur` a horizontal strip of `left` cells.
pub(crate) fn strips(
    target: &[u32],
    cur: &[u32],
    i: usize,
    prev: u32,
    left: u32,
    out: &mut Vec<u32>,
    emit: &mut impl FnMut(&[u32]),
) {
    if i == target.len() {
        if left == 0 {
            emit(out);
        }
        return;
    }
    // Interlacing: row i may reach the previous row's *old* value, and no more
    // than the target allows.
    let hi = prev.min(target[i]).min(cur[i] + left);
    for v in cur[i]..=hi {
        out[i] = v;
        strips(target, cur, i + 1, cur[i], left - (v - cur[i]), out, emit);
    }
    out[i] = cur[i];
}

/// The row reading word: bottom row to top row, left to right within a row.
fn reading_word(target: &[u32], chain: &[Vec<u32>]) -> Vec<u32> {
    let mut word = Vec::new();
    for i in (0..target.len()).rev() {
        for j in 0..target[i] as usize {
            // Cell (i, j) holds the least k whose shape already covers it.
            let letter = (1..chain.len())
                .find(|&k| chain[k][i] as usize > j)
                .unwrap_or(chain.len() - 1);
            word.push(letter as u32);
        }
    }
    word
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    /// Hand-computed, and the case that fixes the reading convention.
    ///
    /// Shape (2,1), content (1,1,1): the tableaux are [[1,2],[3]] and
    /// [[1,3],[2]], reading as 3,1,2 and 2,1,3, with charges 2 and 1.
    #[test]
    fn charge_on_standard_words_matches_hand_computation() {
        assert_eq!(charge(&[3, 1, 2]), 2);
        assert_eq!(charge(&[2, 1, 3]), 1);
        // (2,2) content 1^4: 3,4,1,2 and 2,4,1,3 give 4 and 2.
        assert_eq!(charge(&[3, 4, 1, 2]), 4);
        assert_eq!(charge(&[2, 4, 1, 3]), 2);
        // the identity word has charge 0 + 1 + ... + (n-1)
        assert_eq!(charge(&[1, 2, 3, 4]), 6);
        // fully decreasing: every letter left of its predecessor
        assert_eq!(charge(&[4, 3, 2, 1]), 0);
    }

    #[test]
    fn charge_of_a_repeated_word_splits_into_subwords() {
        // content (2,1) — one subword of length 2, one of length 1
        assert_eq!(charge(&[1, 1, 2]), charge(&[1, 2]) + charge(&[1]));
        assert_eq!(charge(&[]), 0);
        assert_eq!(charge(&[1, 1, 1]), 0);
    }

    /// K_{λλ}(t) = 1 and K_{λμ}(t) has K(1) equal to the ordinary Kostka number
    /// — the specialization that ties this to machinery already tested.
    #[test]
    fn kostka_foulkes_specializes_to_kostka_at_one() {
        for n in 0..=7u32 {
            let parts = crate::partitions_of(n);
            for lambda in &parts {
                for mu in &parts {
                    let k: QtPoly<i64> = kostka_foulkes_by_charge(lambda, mu);
                    assert_eq!(
                        k.eval(&0, &1) as u128,
                        crate::kostka::kostka(lambda, mu),
                        "K_{lambda}{mu}(1) must be the Kostka number"
                    );
                    if lambda == mu {
                        assert_eq!(k, <QtPoly<i64> as Ring>::one(), "K_{lambda}{lambda} = 1");
                    }
                }
            }
        }
    }

    #[test]
    fn known_kostka_foulkes_values() {
        // Sage: KostkaFoulkesPolynomial([2,1],[1,1,1],t) = t^2 + t
        let k: QtPoly<i64> = kostka_foulkes_by_charge(&part(&[2, 1]), &part(&[1, 1, 1]));
        assert_eq!(k.coeff(0, 1), 1);
        assert_eq!(k.coeff(0, 2), 1);
        assert_eq!(k.len(), 2, "{k}");
        // K_{[2,2],[1,1,1,1]} = t^4 + t^2
        let k: QtPoly<i64> = kostka_foulkes_by_charge(&part(&[2, 2]), &part(&[1, 1, 1, 1]));
        assert_eq!(k.coeff(0, 2), 1);
        assert_eq!(k.coeff(0, 4), 1);
        assert_eq!(k.len(), 2, "{k}");
        // K_{[3,2,1],[2,2,1,1]} = t^3 + 2t^2 + t
        let k: QtPoly<i64> = kostka_foulkes_by_charge(&part(&[3, 2, 1]), &part(&[2, 2, 1, 1]));
        assert_eq!(k.coeff(0, 1), 1);
        assert_eq!(k.coeff(0, 2), 2);
        assert_eq!(k.coeff(0, 3), 1);
        assert_eq!(k.len(), 3, "{k}");
    }

    /// Non-zero exactly when λ dominates μ, matching the Kostka support.
    #[test]
    fn support_matches_dominance() {
        for n in 0..=7u32 {
            let parts = crate::partitions_of(n);
            for lambda in &parts {
                for mu in &parts {
                    let k: QtPoly<i64> = kostka_foulkes_by_charge(lambda, mu);
                    let kostka = crate::kostka::kostka(lambda, mu);
                    assert_eq!(k.is_zero(), kostka == 0, "support at {lambda}, {mu}");
                }
            }
        }
    }
}
