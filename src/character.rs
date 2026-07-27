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

use std::collections::HashMap;

use crate::coeff::Ring;
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

/// χ^λ(μ) directly in the coefficient ring `C`, with no fixed-width ceiling.
///
/// The `i128` path above is the fast one and covers every size anyone computes
/// in practice (n ≲ 58), so it is tried first and its global memo does the
/// work. Only when it reports overflow does the recursion re-run in `C` itself,
/// which is **exact** for a bignum ring (`rug::Integer` under the `gmp`
/// feature) and merely wraps differently for a fixed-width one — a fixed-width
/// `C` cannot represent the value either way, and that limit is the caller's
/// choice of ring, not this module's.
///
/// This is the seam that lets the `gmp` feature lift the character ceiling
/// entirely, in the same spirit as [`Ring::from_u128`].
pub fn character_in<C: Ring>(lambda: &Partition, mu: &Partition) -> C {
    if let Some(v) = try_character(lambda, mu) {
        return C::from_i128(v);
    }
    // Past i128. The global memo is typed `i128` and cannot hold these, so the
    // fallback carries a memo of its own, local to this call — the recursion
    // still has heavily overlapping subproblems and is hopeless without one.
    let mut memo = HashMap::new();
    character_generic(lambda, mu, &mut memo)
}

fn character_generic<C: Ring>(
    lambda: &Partition,
    mu: &Partition,
    memo: &mut HashMap<(Partition, Partition), C>,
) -> C {
    if lambda.is_empty() && mu.is_empty() {
        return C::one();
    }
    if mu.is_empty() || lambda.size() != mu.size() {
        return C::zero();
    }
    let key = (lambda.clone(), mu.clone());
    if let Some(v) = memo.get(&key) {
        return v.clone();
    }
    let r = mu.part(0);
    let rest = Partition::from_sorted(mu.parts()[1..].to_vec());
    let mut total = C::zero();
    for (next, height) in border_strips(lambda, r) {
        let sub = character_generic(&next, &rest, memo);
        if height % 2 == 0 {
            total.add_assign(&sub);
        } else {
            total.sub_assign(&sub);
        }
    }
    memo.insert(key, total.clone());
    total
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
    // β₀ = λ₁ + ℓ − 1 is the largest β, so the whole set fits a u64 whenever it
    // is below 64 — true for every |λ| ≤ 32, which is the range anyone computes
    // characters in. Past that, fall through to the allocating form below.
    if lambda.part(0) as usize + l - 1 < 64 {
        return border_strips_masked(lambda, l, r);
    }
    border_strips_general(lambda, r)
}

/// The general form, with no width limit on β. Retained as the reference the
/// masked path is checked against, and as the fallback past |λ| = 32.
fn border_strips_general(lambda: &Partition, r: u32) -> Vec<(Partition, u32)> {
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

/// [`border_strips`] with the β-set held in a u64 instead of on the heap.
///
/// Same mathematics; the difference is entirely representation, and that
/// difference is most of the cost of a character. This runs at *every node* of
/// the Murnaghan–Nakayama recursion, and the general form above allocates a
/// `Vec` for β, a **`HashSet`** for membership, and then per strip another `Vec`
/// plus a sort — so a single χ^λ(μ) did thousands of heap allocations to do
/// arithmetic that fits in registers. Here membership is a bit test, "how many β
/// lie strictly between" is a masked `count_ones`, and the only remaining
/// allocation is the partition each strip has to return.
fn border_strips_masked(lambda: &Partition, l: usize, r: u32) -> Vec<(Partition, u32)> {
    let mut mask = 0u64;
    for i in 0..l {
        mask |= 1u64 << (lambda.part(i) as usize + l - 1 - i);
    }
    let mut out = Vec::new();
    let mut rest = mask;
    while rest != 0 {
        let b = rest.trailing_zeros();
        rest &= rest - 1;
        if b < r {
            continue; // β − r would be negative
        }
        let bp = b - r;
        if mask >> bp & 1 == 1 {
            continue; // that β is taken: no such rim hook
        }
        let between = mask & (((1u64 << b) - 1) ^ ((1u64 << (bp + 1)) - 1));
        let nm = (mask & !(1u64 << b)) | (1u64 << bp);
        out.push((mask_to_partition(nm, l), between.count_ones()));
    }
    out
}

/// Bits high-to-low are β₀ > β₁ > …, and λ_i = β_i − (ℓ−1−i).
fn mask_to_partition(mask: u64, l: usize) -> Partition {
    let mut parts = Vec::new();
    let mut rest = mask;
    let mut i = 0usize;
    while rest != 0 {
        let b = 63 - rest.leading_zeros() as usize;
        rest &= !(1u64 << b);
        let val = b - (l - 1 - i);
        if val > 0 {
            parts.push(val as u32);
        }
        i += 1;
    }
    Partition::from_sorted(parts)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    /// The masked and allocating forms of `border_strips` must yield the same
    /// strips with the same heights. Only the representation differs, and the
    /// general form is the one already checked against Sage.
    ///
    /// Compared as multisets: the masked form walks β upward and the general one
    /// downward, so the orders are reversed. That is invisible to every caller,
    /// which only sums the contributions — but the first version of this test
    /// asserted on order and failed, which is how the difference was noticed
    /// rather than assumed harmless.
    #[test]
    fn masked_border_strips_match_the_general_form() {
        for n in 1..=11u32 {
            for lambda in crate::memo::partitions_cached(n).iter() {
                for r in 1..=n {
                    let mut got = border_strips_masked(lambda, lambda.len(), r);
                    let mut want = border_strips_general(lambda, r);
                    got.sort();
                    want.sort();
                    assert_eq!(got, want, "strips of size {r} from {lambda}");
                }
            }
        }
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
