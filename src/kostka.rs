//! Kostka numbers K_{λμ} = number of semistandard Young tableaux of shape λ and
//! content μ (μ_i entries equal to i, rows weakly increasing, columns strictly
//! increasing).
//!
//! These are the s → m transition: s_λ = Σ_μ K_{λμ} m_μ. The reverse, m → s,
//! inverts the (unitriangular) Kostka matrix; see `convert`.
//!
//! ## How it is computed
//!
//! Not by enumerating tableaux. Grouping an SSYT's cells by value gives the
//! standard bijection with **chains of horizontal strips**
//!
//! ```text
//!   ∅ = λ⁰ ⊆ λ¹ ⊆ … ⊆ λ^ℓ(μ) = λ,   λⁱ/λⁱ⁻¹ a horizontal strip of size μᵢ
//! ```
//!
//! (the cells holding i form a horizontal strip, since two i's in one column
//! would break column-strictness). So K_{λμ} counts chains, and chains through
//! the same intermediate shape **merge** — which is the whole difference
//! between counting and enumerating. This is `strip_lr`'s dynamic program
//! without the lattice condition, which is exactly what makes it count SSYT
//! rather than LR tableaux.
//!
//! Every intermediate shape is pruned to λ, so the state space is bounded by the
//! partitions inside λ rather than by the tableaux of shape λ.
//!
//! The predecessor enumerated SSYT one cell at a time and was exponential: a
//! single K_{λμ} at degree 20 took 638 ms, which made `convert_s_to_m` on one
//! degree-20 Schur function take 400 seconds (`examples/bench_ops.rs`).

use std::collections::HashMap;

use crate::memo::kostka_cached;
use crate::partition::Partition;

/// K_{λμ}. Requires nothing of the arguments beyond being partitions; returns 0
/// unless |λ| = |μ|.
pub fn kostka(lambda: &Partition, mu: &Partition) -> u128 {
    if lambda.size() != mu.size() {
        return 0;
    }
    if lambda.is_empty() {
        return 1; // both empty: the empty tableau
    }
    kostka_cached(lambda, mu, || kostka_uncached(lambda, mu))
}

fn kostka_uncached(lambda: &Partition, mu: &Partition) -> u128 {
    if mu.is_empty() {
        return 0;
    }
    let bound = lambda.parts();
    // Frontier of the chain: intermediate shape -> number of ways to reach it.
    let mut cur: HashMap<Vec<u32>, u128> = HashMap::new();
    cur.insert(Vec::new(), 1);
    let mut next: HashMap<Vec<u32>, u128> = HashMap::new();
    let mut buf: Vec<u32> = Vec::new();

    for &r in mu.parts() {
        next.clear();
        for (shape, &ways) in cur.iter() {
            buf.clear();
            buf.extend_from_slice(shape);
            buf.resize(bound.len(), 0);
            grow(0, r, u32::MAX, &mut buf, bound, &mut |grown: &[u32]| {
                let end = grown.iter().rposition(|&x| x > 0).map_or(0, |i| i + 1);
                *next.entry(grown[..end].to_vec()).or_insert(0) += ways;
            });
        }
        std::mem::swap(&mut cur, &mut next);
        if cur.is_empty() {
            return 0;
        }
    }
    cur.get(bound).copied().unwrap_or(0)
}

/// Add a horizontal strip of `left` cells to `shape` in place, staying inside
/// `bound`, calling `emit` on each result.
///
/// A horizontal strip means the interlacing `shape_{i-1} ≥ new_i ≥ shape_i`, so
/// row `i` may grow only up to the *previous row's old* value — and never past
/// `bound_i`, which is what keeps the frontier proportional to the partitions
/// inside λ rather than to all partitions.
fn grow(
    i: usize,
    left: u32,
    prev_old: u32,
    shape: &mut Vec<u32>,
    bound: &[u32],
    emit: &mut impl FnMut(&[u32]),
) {
    if i == shape.len() {
        if left == 0 {
            emit(shape);
        }
        return;
    }
    // What all remaining rows can still absorb; bail when `left` cannot fit.
    let capacity: u32 = bound[i..]
        .iter()
        .zip(shape[i..].iter())
        .map(|(b, s)| b.saturating_sub(*s))
        .sum();
    if capacity < left {
        return;
    }
    let old = shape[i];
    let hi = prev_old.min(bound[i]).min(old + left);
    if hi < old {
        return;
    }
    for v in old..=hi {
        shape[i] = v;
        grow(i + 1, left - (v - old), old, shape, bound, emit);
    }
    shape[i] = old;
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    #[test]
    fn small_kostka_values() {
        assert_eq!(kostka(&p(&[2]), &p(&[2])), 1);
        assert_eq!(kostka(&p(&[2]), &p(&[1, 1])), 1); // tableau "1 2"
        assert_eq!(kostka(&p(&[1, 1]), &p(&[2])), 0); // column can't repeat
        assert_eq!(kostka(&p(&[1, 1]), &p(&[1, 1])), 1);
        assert_eq!(kostka(&p(&[2, 1]), &p(&[1, 1, 1])), 2); // standard tableaux of (2,1)
        assert_eq!(kostka(&p(&[3, 2, 1]), &p(&[1, 1, 1, 1, 1, 1])), 16); // f^{321}
    }

    /// The old implementation stored tableau entries in `u8` and wrapped its
    /// loop bound, so `ℓ(μ) = 256` silently returned 0 instead of the true
    /// count. The chain DP stores no entries at all, so that boundary cannot
    /// exist — and the case is now cheap enough to assert through the public
    /// API, where before it was exponential and needed a white-box test.
    ///
    /// K_{1ⁿ,1ⁿ} = 1: a single column admits exactly one SSYT.
    #[test]
    fn value_range_extends_past_the_u8_boundary() {
        for len in [200usize, 255, 256, 300] {
            let tall = Partition::new(std::iter::repeat(1).take(len));
            assert_eq!(kostka(&tall, &tall), 1, "K_{{1^{len},1^{len}}}");
        }
    }

    /// `K_{λ,1ⁿ}` counts standard tableaux, so the hook-length formula
    /// `n!/∏h(i,j)` gives it independently — sharing no code with the chain DP.
    ///
    /// This reaches degrees the Sage fixture does not: the rewrite made shapes
    /// like `[6,5,4,3,2]` go from 638 ms per value to microseconds, and a
    /// correctness check that only covers what the *old* implementation could
    /// afford would not exercise the new one where it now operates.
    #[test]
    fn standard_tableaux_count_matches_the_hook_formula() {
        for parts in [
            &[3u32, 2, 1][..],
            &[4, 3, 2, 1],
            &[5, 4, 3, 2, 1],
            &[6, 5, 4, 3, 2],
            &[7, 6, 5, 4, 3],
        ] {
            let lam = p(parts);
            let n = lam.size();
            let ones = Partition::new(std::iter::repeat(1).take(n as usize));

            // n! / ∏ hooks, both exact in u128 at these degrees.
            let conj = lam.conjugate();
            let mut hooks: u128 = 1;
            for (i, &row) in lam.parts().iter().enumerate() {
                for j in 0..row as usize {
                    hooks *= (row as u128 - j as u128)
                        + (conj.part(j) as u128 - i as u128)
                        - 1;
                }
            }
            let factorial: u128 = (1..=u128::from(n)).product();
            assert_eq!(factorial % hooks, 0, "hook product must divide n! for {lam}");
            assert_eq!(kostka(&lam, &ones), factorial / hooks, "K_{{{lam},1^{n}}}");
        }
    }

    #[test]
    fn diagonal_is_one() {
        for parts in [&[2, 1][..], &[3, 1], &[2, 2], &[4]] {
            let lam = p(parts);
            assert_eq!(kostka(&lam, &lam), 1, "K_{{{lam},{lam}}}");
        }
    }
}
