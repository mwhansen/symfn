//! Irreducible characters χ^λ(μ) of the symmetric group, via the
//! Murnaghan–Nakayama rule.
//!
//! These are the entries of the s ↔ p transition: p_μ = Σ_λ χ^λ(μ) s_λ, and
//! (over a field) s_λ = Σ_μ z_μ⁻¹ χ^λ(μ) p_μ. Characters are integers, so p → s
//! stays over ℤ; only s → p pulls in the z_μ⁻¹ denominators.
//!
//! **Range.** |χ^λ(μ)| ≤ d_λ = χ^λ(1ⁿ), and Σ_λ d_λ² = n! gives max d_λ ≈ √(n!).
//! That passes `i64` at **n ≈ 35** and `i128` at **n ≈ 58**. Values are therefore
//! accumulated in `i128`, and every addition is checked: past the ceiling
//! [`try_character`] reports `None` rather than returning a wrapped answer.
//! This module previously returned `i64` and wrapped silently — χ^λ(1³⁶) came
//! back *negative*, which a dimension cannot be, and the corrupted value flowed
//! into the s ↔ p conversions unnoticed.

use crate::memo::character_cached;
use crate::partition::Partition;

/// χ^λ(μ): the value of the irreducible character indexed by λ on the conjugacy
/// class of cycle type μ. Requires |λ| = |μ|.
///
/// # Panics
///
/// If the value exceeds `i128` (around n ≈ 58 — see the module docs). Use
/// [`try_character`] to handle that case instead of panicking. Panicking is the
/// deliberate default: the alternative this replaced was a silently wrong
/// answer, and every caller in the crate feeds these into transition matrices
/// where a wrong value is indistinguishable from a right one.
pub fn character(lambda: &Partition, mu: &Partition) -> i128 {
    try_character(lambda, mu).unwrap_or_else(|| {
        panic!(
            "character χ^{lambda}({mu}) overflows i128 (|λ| = {}); \
             use try_character to handle this",
            lambda.size()
        )
    })
}

/// χ^λ(μ), returning `None` if the value does not fit in `i128`.
///
/// The fallible form of [`character`]. Overflow is a genuine possibility only
/// for |λ| ≳ 58; below that this never returns `None`.
pub fn try_character(lambda: &Partition, mu: &Partition) -> Option<i128> {
    if lambda.is_empty() && mu.is_empty() {
        return Some(1);
    }
    if mu.is_empty() || lambda.size() != mu.size() {
        return Some(0);
    }
    // Memoized: the recursion below re-enters `try_character`, so every shared
    // subproblem across the whole Murnaghan-Nakayama tree is cached too.
    // `None` (overflow) is deliberately *not* cached — it is not a value, and
    // the table would then have to distinguish "absent" from "unrepresentable".
    character_cached(lambda, mu, || character_uncached(lambda, mu))
}

/// `Some(v)` on success; `None` if any partial sum leaves `i128`.
fn character_uncached(lambda: &Partition, mu: &Partition) -> Option<i128> {
    // Peel off the first part of μ, summing signed contributions over every
    // border strip (rim hook) of that size removable from λ.
    let r = mu.part(0);
    let rest = Partition::from_sorted(mu.parts()[1..].to_vec());
    let mut total: i128 = 0;
    for (next, height) in border_strips(lambda, r) {
        let sub = try_character(&next, &rest)?;
        let term = if height % 2 == 0 { sub } else { -sub };
        total = total.checked_add(term)?;
    }
    Some(total)
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

    /// d_λ by the hook-length formula, with gcd cancellation so no intermediate
    /// factorial materializes. Independent of the Murnaghan–Nakayama recursion,
    /// which is the point: it is ground truth the recursion is checked against.
    fn dimension_by_hooks(lam: &[u32]) -> i128 {
        fn gcd(mut a: i128, mut b: i128) -> i128 {
            while b != 0 {
                let t = a % b;
                a = b;
                b = t;
            }
            a
        }
        if lam.is_empty() {
            return 1;
        }
        let m: usize = lam.iter().map(|&x| x as usize).sum();
        let mut conj = vec![0i64; lam[0] as usize];
        for &part in lam {
            for c in conj.iter_mut().take(part as usize) {
                *c += 1;
            }
        }
        // d_λ = m! / ∏ hooks; cancel each hook into the numerator factors, so
        // every hook reduces to 1 and what remains is exactly the integer d_λ.
        let mut num: Vec<i128> = (2..=m as i128).collect();
        for (i, &row) in lam.iter().enumerate() {
            for j in 0..row as usize {
                let mut h = (row as i64 - j as i64) + (conj[j] - i as i64) - 1;
                for ni in num.iter_mut() {
                    if h == 1 {
                        break;
                    }
                    let g = gcd(*ni, h as i128);
                    if g > 1 {
                        *ni /= g;
                        h /= g as i64;
                    }
                }
                assert_eq!(h, 1, "hook did not cancel — d_λ is not an integer");
            }
        }
        num.into_iter().product()
    }

    /// Regression: characters used to accumulate in `i64` and wrap silently.
    /// χ^λ(1³⁶) came back *negative* — impossible for a dimension — and fed the
    /// s ↔ p conversions as if it were correct. The first failing size is n=36,
    /// so the sweep straddles it.
    #[test]
    fn identity_class_matches_hook_formula_past_the_i64_ceiling() {
        for n in [20u32, 30, 34, 36, 40, 48] {
            // A staircase-ish λ ⊢ n: 1+2+3+… truncated to fit.
            let (mut parts, mut left, mut k) = (Vec::new(), n, 1);
            while left > 0 {
                let part = k.min(left);
                parts.push(part);
                left -= part;
                k += 1;
            }
            parts.sort_unstable_by(|a, b| b.cmp(a));
            let lam = Partition::new(parts.iter().copied());
            let ones = Partition::new(std::iter::repeat(1).take(n as usize));

            let want = dimension_by_hooks(&parts);
            assert_eq!(character(&lam, &ones), want, "χ^{lam}(1^{n})");
            assert!(want > 0, "a dimension must be positive");
            // The old i64 path could not even hold this above n = 35.
            if n >= 36 {
                assert!(want > i64::MAX as i128, "n={n} should exceed i64 range");
            }
        }
    }

    /// Past `i128` the fallible form reports overflow instead of wrapping.
    /// max d_λ ≈ √(n!) crosses i128::MAX around n ≈ 58.
    #[test]
    fn overflow_is_reported_not_wrapped() {
        let n = 70u32;
        let lam = Partition::new(std::iter::once(n)); // trivial character: χ = 1
        let ones = Partition::new(std::iter::repeat(1).take(n as usize));
        assert_eq!(try_character(&lam, &ones), Some(1), "χ^(n)(1ⁿ) = 1, no overflow");

        // A wide shape at the same size has an astronomically large dimension.
        let big = Partition::new([n / 2, n / 2 - 1, n / 2 - 2, n / 2 - 3].iter().copied());
        let m = big.size();
        let ones_m = Partition::new(std::iter::repeat(1).take(m as usize));
        match try_character(&big, &ones_m) {
            None => {}                              // reported overflow: correct
            Some(v) => assert!(v > 0, "if it fits, it must still be positive"),
        }
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
