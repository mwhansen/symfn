//! Irreducible characters χ^λ(μ) of the symmetric group, via the
//! Murnaghan–Nakayama rule.
//!
//! These are the entries of the s ↔ p transition: p_μ = Σ_λ χ^λ(μ) s_λ, and
//! (over a field) s_λ = Σ_μ z_μ⁻¹ χ^λ(μ) p_μ. Characters are integers, so p → s
//! stays over ℤ; only s → p pulls in the z_μ⁻¹ denominators.
//!
//! **Range.** |χ^λ(μ)| ≤ d_λ = χ^λ(1ⁿ) and Σ_λ d_λ² = n!, so the largest value
//! is about √(n!): past `i64` at **n ≈ 35**, past `i128` at **n ≈ 58**.
//! Accumulation is in `i128` and every addition is checked — [`character`]
//! panics at the ceiling, [`try_character`] returns `None`, neither wraps.
//! [`character_in`] has no ceiling at all over a bignum ring.
//!
//! In Sage, one value is Symmetrica's `charvalue` and the whole table is
//! `chartafel`, both reached through `sage.libs.symmetrica.all`
//! (`scripts/compare_symmetrica.py`).

use std::collections::HashMap;

use crate::coeff::Ring;
use crate::convert::Beta;
use crate::memo::{character_cached, MaskKey};
use crate::partition::Partition;

/// χ^λ(μ): the value of the irreducible character indexed by λ on the conjugacy
/// class of cycle type μ. Zero when |λ| ≠ |μ|, and 1 for the empty pair.
///
/// # Panics
///
/// Panics if the value exceeds `i128`, around n ≈ 58. [`try_character`] returns
/// `None` there instead, and [`character_in`] has no ceiling over a bignum
/// ring.
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
/// `None` means that and only that — never an input error — and takes n ≈ 58
/// to reach. |λ| ≠ |μ| is `Some(0)`, and the empty pair is `Some(1)`.
/// [`character`] panics where this returns `None`; [`character_in`] has no
/// ceiling over a bignum ring.
pub fn try_character(lambda: &Partition, mu: &Partition) -> Option<i128> {
    if lambda.is_empty() && mu.is_empty() {
        return Some(1);
    }
    if mu.is_empty() || lambda.size() != mu.size() {
        return Some(0);
    }
    if fits_mask::<u64>(lambda) && fits_mask::<u64>(mu) {
        return character_masked::<u64>(lambda, mu);
    }
    if fits_mask::<u128>(lambda) && fits_mask::<u128>(mu) {
        return character_masked::<u128>(lambda, mu);
    }
    // Past both mask widths: the partition-keyed recursion, memoized the same
    // way — it re-enters `try_character`, so every shared subproblem across
    // the whole Murnaghan-Nakayama tree is cached too. `None` (overflow) is
    // deliberately *not* cached, on any path — it is not a value, and the
    // table would then have to distinguish "absent" from "unrepresentable".
    character_cached(lambda, mu, || character_uncached(lambda, mu))
}

/// Whether λ's β-set fits a mask of width `M`: β₀ = λ₁ + ℓ − 1 < `M::BITS`.
/// Since λ₁ + ℓ − 1 ≤ |λ|, a `u64` holds every partition of at most 63 —
/// already past the `i128` ceiling on the values, n ≈ 58 — and a `u128` every
/// partition of at most 127, where an `i128` character exists only for shapes
/// like a hook. The wide instantiation costs nothing beyond the trait, and
/// keeps the partition-keyed recursion for the sizes where nothing fixed-width
/// applies anyway.
fn fits_mask<M: Beta>(lambda: &Partition) -> bool {
    lambda.is_empty() || (lambda.part(0) as usize + lambda.len() - 1) < M::BITS as usize
}

/// The β-mask of λ with ℓ(λ) slots, which is the canonical one
/// (`Beta::canonical`).
// A bit position below `M::BITS`, which `fits_mask` has checked.
#[allow(clippy::cast_possible_truncation)]
fn beta_mask<M: Beta>(lambda: &Partition) -> M {
    let l = lambda.len();
    let mut mask = M::low_ones(0);
    for i in 0..l {
        mask = mask.with_bit((lambda.part(i) as usize + l - 1 - i) as u32);
    }
    mask
}

/// χ^λ(μ) with the whole recursion on β-masks and a word-keyed memo.
///
/// The recursion of [`character_uncached`] built a `Partition` for every strip
/// and for every μ-suffix, then cloned both into a `(Partition, Partition)`
/// key at every node to look it up — allocation and SipHash were most of a
/// character, and so most of `s → p` (`docs/record/transitions.md`). Here a
/// node is `(β-mask of λ', β-mask of the μ-suffix)`: the strips are bit
/// operations on the first, the suffix is the second with its top bit
/// removed, and the pair is the memo key. Nothing on the path allocates but
/// the local memo itself.
///
/// The recursion reads the shared memo under one guard held for the whole
/// call and writes only to a local map, merged in afterwards
/// (`memo::character_masks_read`).
fn character_masked<M: MaskKey>(lambda: &Partition, mu: &Partition) -> Option<i128> {
    let mut ctx = MaskedRecursion::<M> {
        mu: mu.parts(),
        local: Default::default(),
        shared: crate::memo::character_masks_read::<M>(),
    };
    let result = ctx.chi(beta_mask::<M>(lambda), beta_mask::<M>(mu), 0);
    let MaskedRecursion { local, shared, .. } = ctx;
    drop(shared);
    // Overflow (`None`) is never stored, as in `character_cached`; but every
    // value that was completed below the overflowing node is a character and
    // is kept.
    crate::memo::character_masks_store::<M>(local);
    result
}

/// One row of the character table: χ^λ(μ) for every μ ⊢ |λ|, in the order of
/// `memo::partitions_cached(|λ|)`. An entry is `None` exactly where
/// [`try_character`] is — the value passes `i128`, from |λ| ≈ 58.
///
/// One [`character_masked`] call per entry pays the memo round trip — the
/// shared-guard acquire, a fresh local map, the merge back — p(n) times per
/// row, and on `s → p` that traffic measured at 21% of the whole conversion
/// (`docs/record/coefficient-arithmetic.md`). Here the row runs under one
/// guard against one local map, merged once; a subproblem one μ completes is
/// read straight out of the local map by every later μ that reaches it.
pub(crate) fn character_row(lambda: &Partition) -> Vec<Option<i128>> {
    let n = lambda.size() as usize;
    let mus = crate::memo::partitions_cached(lambda.size());
    // β₀ = μ₁ + ℓ − 1 ≤ |μ|, so one width check on the degree covers every μ
    // of it (`fits_mask`) — and λ ⊢ n is covered by the same bound.
    if n < u64::BITS as usize {
        return character_row_masked::<u64>(lambda, &mus);
    }
    if n < u128::BITS as usize {
        return character_row_masked::<u128>(lambda, &mus);
    }
    mus.iter().map(|mu| try_character(lambda, mu)).collect()
}

/// [`character_row`] with the whole row in one [`MaskedRecursion`].
///
/// Correct for the same reason the shared table is: a memo key
/// `(canonical λ'-mask, canonical suffix mask)` determines its value with no
/// reference to which top-level μ the recursion entered through, so sibling
/// columns of a row may share one local map exactly as separate calls share
/// the global one.
fn character_row_masked<M: MaskKey>(lambda: &Partition, mus: &[Partition]) -> Vec<Option<i128>> {
    let shared = crate::memo::character_masks_read::<M>();
    // A cold row settles near p(n)·n/4 subproblems — 17,783 at the staircase
    // of 27, 151,776 at 36, 989,123 at 45 — and growing there from empty was
    // 14% of `s → p` in rehashes (`docs/record/coefficient-arithmetic.md`).
    // The estimate only has to land within a doubling, and what the shared
    // table already holds a warm row will not insert, so it comes off the top.
    let cold = mus.len() * (lambda.size() as usize / 4 + 1);
    let mut ctx = MaskedRecursion::<M> {
        mu: &[],
        local: crate::fasthash::Map::with_capacity_and_hasher(
            cold.saturating_sub(shared.len()),
            Default::default(),
        ),
        shared,
    };
    let lam = beta_mask::<M>(lambda);
    let row = mus
        .iter()
        .map(|mu| {
            ctx.mu = mu.parts();
            ctx.chi(lam, beta_mask::<M>(mu), 0)
        })
        .collect();
    let MaskedRecursion { local, shared, .. } = ctx;
    drop(shared);
    // As in `character_masked`: overflow (`None`) is never stored, and every
    // value completed below an overflowing node is a character and is kept.
    crate::memo::character_masks_store::<M>(local);
    row
}

/// The state of one [`character_masked`] call — or of one whole
/// [`character_row_masked`] row, which swaps `mu` between columns and keeps
/// the maps.
struct MaskedRecursion<'a, M: MaskKey> {
    mu: &'a [u32],
    local: crate::fasthash::Map<(M, M), i128>,
    shared: std::sync::RwLockReadGuard<'static, crate::fasthash::Map<(M, M), i128>>,
}

impl<M: MaskKey> MaskedRecursion<'_, M> {
    /// χ^{λ'}(μ[k..]) for the β-mask `mask` of λ' and the canonical β-mask
    /// `mu_mask` of the suffix μ[k..].
    ///
    /// `mask` keeps the same number of slots at every depth: removing a rim
    /// hook moves one bit down and never drops one, so a part that reaches
    /// zero is a low bit, not a missing one. `mu_mask` loses its top bit at
    /// each depth, which is exactly dropping the first part: the other β stay
    /// where they are, since their index and the length both fall by one.
    // `highest()` is a bit position of a mask, below `M::BITS`.
    #[allow(clippy::cast_possible_truncation)]
    fn chi(&mut self, mask: M, mu_mask: M, k: usize) -> Option<i128> {
        if k == self.mu.len() {
            // |λ'| = |μ[k..]| at every node, so λ' is empty here and χ = 1.
            return Some(1);
        }
        let key = (mask.canonical(), mu_mask);
        if let Some(&v) = self.local.get(&key) {
            return Some(v);
        }
        if let Some(&v) = self.shared.get(&key) {
            return Some(v);
        }
        crate::interrupt::poll();
        // As `border_strips_masked`: a rim hook of length r is a β moved down
        // by r onto a free position, with height the number of β passed over.
        let r = self.mu[k];
        let rest_mu = mu_mask.without_bit(mu_mask.highest() as u32);
        let mut total: i128 = 0;
        let mut rest = mask;
        while !rest.is_zero() {
            let b = rest.trailing_zeros();
            rest = rest.clear_lowest();
            if b < r {
                continue;
            }
            let bp = b - r;
            if mask.test(bp) {
                continue;
            }
            let height = mask.between(bp + 1, b).count_ones();
            let next = mask.without_bit(b).with_bit(bp);
            let sub = self.chi(next, rest_mu, k + 1)?;
            let term = if height.is_multiple_of(2) { sub } else { -sub };
            total = total.checked_add(term)?;
        }
        self.local.insert(key, total);
        Some(total)
    }
}

/// χ^λ(μ) directly in the coefficient ring `C`, with no fixed-width ceiling.
///
/// The memoized `i128` path runs first and answers everything below n ≈ 58.
/// Past that the recursion re-runs in `C`: **exact** for a bignum ring,
/// and no better than `C` itself for a fixed-width one.
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

pub(crate) fn character_generic<C: Ring>(
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
        crate::interrupt::poll();
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
        crate::interrupt::poll();
        let sub = try_character(&next, &rest)?;
        let term = if height % 2 == 0 { sub } else { -sub };
        total = total.checked_add(term)?;
    }
    Some(total)
}

/// The whole character table of S_n, as `table[i][j] = χ^{λⁱ}(λʲ)`, with rows
/// and columns both indexed by
/// `memo::partitions_cached`.
///
/// `p(n)²` values in `p(n)` Murnaghan–Nakayama sweeps: one sweep is a whole
/// column, and sweeps whose μ share a prefix share that initial segment.
///
/// **Range.** Entries pass `i128` near n ≈ 58, but the table is `p(n)²`
/// values — 1 GB at n = 32 — so memory walls about twenty degrees earlier.
/// [`character_table_in`] takes an arbitrary ring and does not move that.
pub fn character_table(n: u32) -> Vec<Vec<i128>> {
    character_table_in(n)
}

/// [`character_table`] over an arbitrary coefficient ring.
///
/// Not a widening: the entries are the same integers, and the memory wall on
/// [`character_table`] arrives first whatever `C` is. It exists so the sweep
/// can carry a ring whose elements are not integers at all — characters as
/// constants in a polynomial ring, the shape Kostka–Foulkes accumulates
/// through. For a single character past `i128`, [`character_in`] is the form
/// that escalates.
pub fn character_table_in<C: Ring>(n: u32) -> Vec<Vec<C>> {
    let parts = crate::memo::partitions_cached(n);
    let mut table = vec![vec![C::zero(); parts.len()]; parts.len()];
    let l = n as usize;
    if l == 0 {
        table[0][0] = C::one();
        return table;
    }
    if l > crate::convert::MASK_LIMIT {
        // Past the β-mask width; fall back to the per-entry recursion, which is
        // exact in `C` for a bignum ring.
        for (i, lambda) in parts.iter().enumerate() {
            for (j, mu) in parts.iter().enumerate() {
                crate::interrupt::poll();
                table[i][j] = character_in::<C>(lambda, mu);
            }
        }
        return table;
    }
    // β-mask of each λ, so a swept column can be indexed straight back to a row.
    // The sink below hits this once per emitted entry on a bare `u64`, which is
    // exactly the shape [`crate::fasthash`] is for; SipHash here cost 1.12x
    // (`character_beta_sweep_n28` in examples/bench_ops.rs).
    let index: crate::fasthash::Map<u64, usize> = parts
        .iter()
        .enumerate()
        .map(|(i, p)| {
            let mut mask = 0u64;
            for k in 0..l {
                mask |= 1u64 << (p.part(k) as usize + l - 1 - k);
            }
            (mask, i)
        })
        .collect();

    // Descending part order, so the longest common prefixes are shared.
    let mut order: Vec<usize> = (0..parts.len()).collect();
    order.sort_by(|&a, &b| parts[a].parts().cmp(parts[b].parts()));
    let items: Vec<(&Partition, usize)> = order.iter().map(|&i| (&parts[i], i)).collect();

    let mut root: crate::fasthash::Map<u64, C> = Default::default();
    root.insert((1u64 << l) - 1, C::one());
    crate::convert::p_expand_shared(&items, 0, &root, &mut |&col: &usize, mask, chi: &C| {
        if let Some(&row) = index.get(&mask) {
            table[row][col] = chi.clone();
        }
    });
    table
}

/// Every way to remove a border strip of length `r` from λ, as
/// `(resulting partition, height)` where height = (#rows spanned) − 1.
///
/// Uses β-numbers: with L = ℓ(λ), set β_i = λ_i + (L−1−i) (strictly decreasing,
/// distinct). Removing a rim hook of length r ⇔ replacing some β_i by β_i − r,
/// provided it is ≥ 0 and not already present; the height is the number of β_j
/// strictly between the new and old values.
pub(crate) fn border_strips(lambda: &Partition, r: u32) -> Vec<(Partition, u32)> {
    let l = lambda.len();
    if l == 0 {
        return Vec::new();
    }
    // β₀ = λ₁ + ℓ − 1 is the largest β, so the whole set fits a u64 whenever it
    // is below 64 — true for every |λ| ≤ 63, since λ₁ + ℓ − 1 ≤ |λ|. Past
    // that, fall through to the allocating form below.
    if lambda.part(0) as usize + l - 1 < 64 {
        return border_strips_masked(lambda, l, r);
    }
    border_strips_general(lambda, r)
}

/// The general form, with no width limit on β. Retained as the reference the
/// masked path is checked against, and as the fallback past |λ| = 63.
// Index arithmetic on a β-set: every part is `u32` by `Partition`'s own
// invariant and every offset is bounded by ℓ(λ), so each value here fits `u32`
// and no intermediate leaves `i64`. Nothing below carries a coefficient.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
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
/// Same mathematics; the difference is entirely representation
/// (`docs/record/oracles-and-comparisons.md`). The general form above
/// allocates a `Vec` for β, a **`HashSet`** for membership, and then per strip
/// another `Vec` plus a sort — heap allocations to do arithmetic that fits in
/// registers. Here membership is a bit test, "how many β lie strictly between"
/// is a masked `count_ones`, and the only remaining allocation is the
/// partition each strip has to return. `MaskedRecursion::chi` runs the same
/// bit operations inline and never builds the partition at all; this is the
/// form for a caller that wants the strips as partitions.
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
// `b` is a set bit of a `u64` β-mask and `l ≤ 32` (`MASK_LIMIT`), so
// `b - (l - 1 - i)` is a part of a partition: non-negative and far inside `u32`.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
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

    /// The swept table must equal the per-entry recursion at every position.
    ///
    /// Two independent things could go wrong quietly: the prefix trie can
    /// misattribute a column, and the β-mask index maps a swept shape back to a
    /// row — an off-by-one in the offset would drop rows rather than corrupt
    /// them, which a spot check would pass.
    /// The generic table must agree with the `i128` one and must accept a ring
    /// that is not an integer type — the Kostka–Foulkes shape.
    #[test]
    fn generic_character_table_agrees_and_accepts_a_polynomial_ring() {
        use crate::qt::QtPoly;
        for n in 0..=9u32 {
            let plain = character_table(n);
            let poly: Vec<Vec<QtPoly<i64>>> = character_table_in(n);
            for i in 0..plain.len() {
                for j in 0..plain.len() {
                    assert_eq!(
                        poly[i][j].coeff(0, 0),
                        plain[i][j] as i64,
                        "({i},{j}) deg {n}"
                    );
                    assert!(poly[i][j].len() <= 1, "a character is a constant");
                }
            }
        }
    }

    /// The β-mask recursion must agree with the partition-keyed one it
    /// replaced, at both widths. All pairs through degree 11, and a few
    /// wide-and-long pairs above that where a canonical-mask slip (a suffix
    /// mask off by a shift, say) would key two different subproblems together
    /// and pass every small case.
    #[test]
    fn masked_recursion_matches_the_partition_keyed_one() {
        let general = |lam: &Partition, mu: &Partition| -> Option<i128> {
            if lam.is_empty() && mu.is_empty() {
                return Some(1);
            }
            if mu.is_empty() || lam.size() != mu.size() {
                return Some(0);
            }
            character_cached(lam, mu, || character_uncached(lam, mu))
        };
        let check = |lambda: &Partition, mu: &Partition| {
            let want = general(lambda, mu);
            assert_eq!(
                character_masked::<u64>(lambda, mu),
                want,
                "χ^{lambda}({mu}), u64"
            );
            assert_eq!(
                character_masked::<u128>(lambda, mu),
                want,
                "χ^{lambda}({mu}), u128"
            );
        };
        for n in 0..=11u32 {
            for lambda in crate::memo::partitions_cached(n).iter() {
                for mu in crate::memo::partitions_cached(n).iter() {
                    check(lambda, mu);
                }
            }
        }
        let big = [
            (p(&[9, 8, 7, 6]), p(&[5, 5, 4, 4, 3, 3, 2, 2, 1, 1])),
            (p(&[12, 8, 4, 2]), p(&[7, 6, 5, 4, 3, 1])),
            (p(&[6, 6, 6, 6, 6]), p(&[6, 5, 4, 3, 3, 3, 2, 2, 1, 1])),
            (
                p(&[16, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]),
                p(&[3, 3, 3, 3, 3, 3, 3, 3, 3, 3]),
            ),
        ];
        for (lambda, mu) in &big {
            check(lambda, mu);
        }
        // Past the u64 width, where only the wide mask and the partition path
        // can run: the hook (33, 1³²) has β₀ = 65.
        let lambda = p(&std::iter::once(33)
            .chain(std::iter::repeat_n(1, 32))
            .collect::<Vec<_>>());
        let mu = p(&[8, 8, 8, 8, 8, 8, 8, 8, 1]);
        assert!(!fits_mask::<u64>(&lambda), "the hook must not fit a u64");
        assert_eq!(
            character_masked::<u128>(&lambda, &mu),
            general(&lambda, &mu),
            "χ^{lambda}({mu}), u128"
        );
    }

    /// The batched row must agree entry by entry with the ring-generic
    /// recursion, which touches neither β-mask table. That independence is the
    /// point of not checking against [`try_character`]: a row that stored a
    /// value under a wrong mask key would poison the shared table and then
    /// agree with every per-entry call that reads it. Both widths run, since
    /// [`character_row`] only ever picks `u64` below degree 64.
    #[test]
    fn character_row_matches_the_independent_recursion() {
        let mut memo = HashMap::new();
        for n in 0..=11u32 {
            let mus = crate::memo::partitions_cached(n);
            for lambda in mus.iter() {
                let row = character_row(lambda);
                let wide = character_row_masked::<u128>(lambda, &mus);
                assert_eq!(row.len(), mus.len(), "row of {lambda} is p({n}) long");
                for ((mu, got), w) in mus.iter().zip(&row).zip(&wide) {
                    let want: i128 = character_generic(lambda, mu, &mut memo);
                    assert_eq!(*got, Some(want), "χ^{lambda}({mu})");
                    assert_eq!(*w, Some(want), "χ^{lambda}({mu}), u128");
                }
            }
        }
    }

    #[test]
    fn swept_character_table_matches_per_entry() {
        for n in 0..=11u32 {
            let parts = crate::memo::partitions_cached(n);
            let table = character_table(n);
            assert_eq!(table.len(), parts.len(), "table is p({n}) square");
            for (i, lambda) in parts.iter().enumerate() {
                for (j, mu) in parts.iter().enumerate() {
                    assert_eq!(
                        table[i][j],
                        character(lambda, mu),
                        "χ^{lambda}({mu}) at degree {n}"
                    );
                }
            }
        }
    }

    /// The masked and allocating forms of `border_strips` must yield the same
    /// strips with the same heights. Only the representation differs, and the
    /// general form is the one already checked against Sage.
    ///
    /// Compared as multisets: the masked form walks β upward and the general
    /// one downward, so the orders are reversed. That is invisible to every
    /// caller, which only sums the contributions — but the first version of
    /// this test asserted on order and failed, which is how the difference was
    /// noticed rather than assumed harmless.
    #[test]
    fn masked_border_strips_match_the_general_form() {
        for n in 0..=11u32 {
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
            let ones = Partition::new(std::iter::repeat_n(1, n as usize));

            let want = dimension_by_hooks(&parts);
            assert_eq!(character(&lam, &ones), want, "χ^{lam}(1^{n})");
            assert!(want > 0, "a dimension must be positive");
            // The old i64 path could not even hold this above n = 35.
            if n >= 36 {
                assert!(want > i64::MAX as i128, "n={n} should exceed i64 range");
            }
        }
    }

    /// Past `i128` [`try_character`] reports overflow instead of wrapping.
    /// max d_λ ≈ √(n!) crosses i128::MAX around n ≈ 58.
    #[test]
    fn overflow_is_reported_not_wrapped() {
        let n = 70u32;
        let lam = Partition::new(std::iter::once(n)); // trivial character: χ = 1
        let ones = Partition::new(std::iter::repeat_n(1, n as usize));
        assert_eq!(
            try_character(&lam, &ones),
            Some(1),
            "χ^(n)(1ⁿ) = 1, no overflow"
        );

        // A wide shape at the same size has an astronomically large dimension.
        let big = Partition::new([n / 2, n / 2 - 1, n / 2 - 2, n / 2 - 3].iter().copied());
        let m = big.size();
        let ones_m = Partition::new(std::iter::repeat_n(1, m as usize));
        match try_character(&big, &ones_m) {
            None => {} // reported overflow: correct
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
