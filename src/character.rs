//! Irreducible characters χ^λ(μ) of the symmetric group, via the
//! Murnaghan–Nakayama rule.
//!
//! These are the entries of the s ↔ p transition: p_μ = Σ_λ χ^λ(μ) s_λ, and
//! (over a field) s_λ = Σ_μ z_μ⁻¹ χ^λ(μ) p_μ. Characters are integers, so p → s
//! stays over ℤ; only s → p pulls in the z_μ⁻¹ denominators.

use crate::memo::character_cached;
use crate::partition::Partition;

/// χ^λ(μ): the value of the irreducible character indexed by λ on the conjugacy
/// class of cycle type μ. Requires |λ| = |μ|.
pub fn character(lambda: &Partition, mu: &Partition) -> i64 {
    if lambda.is_empty() && mu.is_empty() {
        return 1;
    }
    if mu.is_empty() || lambda.size() != mu.size() {
        return 0;
    }
    // Memoized: the recursion below re-enters `character`, so every shared
    // subproblem across the whole Murnaghan-Nakayama tree is cached too.
    character_cached(lambda, mu, || character_uncached(lambda, mu))
}

fn character_uncached(lambda: &Partition, mu: &Partition) -> i64 {
    // Peel off the first part of μ, summing signed contributions over every
    // border strip (rim hook) of that size removable from λ.
    let r = mu.part(0);
    let rest = Partition::from_sorted(mu.parts()[1..].to_vec());
    let mut total = 0i64;
    for (next, height) in border_strips(lambda, r) {
        let sign = if height % 2 == 0 { 1 } else { -1 };
        total += sign * character(&next, &rest);
    }
    total
}

/// Every way to remove a border strip of length `r` from λ, as
/// `(resulting partition, height)` where height = (#rows spanned) − 1.
///
/// Uses β-numbers: with L = ℓ(λ), set β_i = λ_i + (L−1−i) (strictly decreasing,
/// distinct). Removing a rim hook of length r ⇔ replacing some β_i by β_i − r,
/// provided it is ≥ 0 and not already present; the height is the number of β_j
/// strictly between the new and old values.
fn border_strips(lambda: &Partition, r: u32) -> Vec<(Partition, u32)> {
    let l = lambda.len();
    if l == 0 {
        return Vec::new();
    }
    let beta: Vec<i64> = (0..l)
        .map(|i| lambda.part(i) as i64 + (l as i64 - 1 - i as i64))
        .collect();
    let present: std::collections::HashSet<i64> = beta.iter().copied().collect();

    let mut out = Vec::new();
    for &b in &beta {
        let bp = b - r as i64;
        if bp < 0 || present.contains(&bp) {
            continue;
        }
        let height = beta.iter().filter(|&&x| x > bp && x < b).count() as u32;

        // New β-set: swap b → bp, re-sort descending, subtract the offsets back.
        let mut nb: Vec<i64> = beta.iter().map(|&x| if x == b { bp } else { x }).collect();
        nb.sort_unstable_by(|a, c| c.cmp(a));
        let mut parts = Vec::new();
        for (k, &x) in nb.iter().enumerate() {
            let val = x - (l as i64 - 1 - k as i64);
            if val > 0 {
                parts.push(val as u32);
            }
        }
        out.push((Partition::from_sorted(parts), height));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    #[test]
    fn s4_character_table_rows() {
        // Irrep (2,2) of S_4 on classes 1^4, 21^2, 2^2, 31, 4: 2, 0, 2, -1, 0.
        let lam = p(&[2, 2]);
        assert_eq!(character(&lam, &p(&[1, 1, 1, 1])), 2); // dimension
        assert_eq!(character(&lam, &p(&[2, 1, 1])), 0);
        assert_eq!(character(&lam, &p(&[2, 2])), 2);
        assert_eq!(character(&lam, &p(&[3, 1])), -1);
        assert_eq!(character(&lam, &p(&[4])), 0);
    }

    #[test]
    fn sign_and_trivial_characters() {
        // Trivial character (n) is all ones; sign character (1^n) is ε.
        for mu in [p(&[3]), p(&[2, 1]), p(&[1, 1, 1])] {
            assert_eq!(character(&p(&[3]), &mu), 1, "trivial at {mu}");
        }
        assert_eq!(character(&p(&[1, 1, 1]), &p(&[1, 1, 1])), 1); // ε(id) = 1
        assert_eq!(character(&p(&[1, 1, 1]), &p(&[2, 1])), -1); // ε(transposition)
        assert_eq!(character(&p(&[1, 1, 1]), &p(&[3])), 1); // ε(3-cycle)
    }
}
