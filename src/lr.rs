//! Littlewood–Richardson coefficients, computed natively in Rust.
//!
//! LR is the central structure constant, and it is computed here rather than
//! by an external C library. The [`LrBackend`] trait keeps the
//! *implementation* an interchangeable detail: [`NaiveLr`] counts LR
//! (Yamanouchi) skew tableaux directly, pruning all three constraints —
//! row-weakly-increasing, column-strictly-increasing, and the ballot/lattice
//! condition — incrementally as it fills.
//!
//! Two faster backends plug in with no caller changes —
//! [`StripLr`](crate::strip_lr::StripLr), a row-strip DP, and
//! [`SkewLr`](crate::skew_lr::SkewLr), which expands a whole shape in one
//! traversal and is now the default via [`AutoLr`](crate::strip_lr::AutoLr).
//! `NaiveLr` keeps its place as the **reference implementation**: it is the
//! most obviously-correct of the three, and the faster ones are held to
//! exhaustive agreement with it. That is the oracle pattern, kept in-house —
//! the rule is V8 in `docs/policies/validation.md`.

// Shape indices.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use crate::memo::{lr_cached, partitions_cached};
use crate::partition::Partition;

/// A source of Littlewood–Richardson coefficients c^λ_{μν}.
pub trait LrBackend {
    /// The LR coefficient c^λ_{μν}: the multiplicity of s_λ in s_μ · s_ν, equal
    /// to the number of LR (Yamanouchi) skew tableaux of shape λ/μ and content ν.
    fn lr_coeff(&self, lambda: &Partition, mu: &Partition, nu: &Partition) -> u128;

    /// Expand the product s_μ · s_ν as Σ_λ c^λ_{μν} s_λ, returning only the
    /// nonzero terms, **sorted by λ**. Provided in terms of
    /// [`LrBackend::lr_coeff`]; a fast backend may override it to avoid scanning
    /// all partitions. Sorted order is part of the contract so that swapping
    /// backends is observationally transparent.
    fn schur_product(&self, mu: &Partition, nu: &Partition) -> Vec<(Partition, u128)> {
        let n = mu.size() + nu.size();
        let mut out = Vec::new();
        for lambda in partitions_cached(n).iter() {
            // Necessary shape bounds: μ ⊆ λ, at most ℓ(μ)+ℓ(ν) rows, and
            // λ₁ ≤ μ₁+ν₁. Cheap filters before the counting backtrack.
            if lambda.len() > mu.len() + nu.len()
                || lambda.part(0) > mu.part(0) + nu.part(0)
                || !lambda.contains(mu)
            {
                continue;
            }
            let c = self.lr_coeff(lambda, mu, nu);
            if c != 0 {
                out.push((lambda.clone(), c));
            }
        }
        out.sort_by(|a, b| a.0.cmp(&b.0));
        out
    }
}

/// Pure-Rust backend: counts Littlewood–Richardson skew tableaux with all
/// constraints pruned incrementally.
///
/// The reference implementation — the faster backends are validated against it.
/// For production use prefer [`AutoLr`](crate::strip_lr::AutoLr).
#[derive(Clone, Copy, Debug, Default)]
pub struct NaiveLr;

impl LrBackend for NaiveLr {
    fn lr_coeff(&self, lambda: &Partition, mu: &Partition, nu: &Partition) -> u128 {
        // Necessary conditions. A non-zero coefficient needs λ to contain both
        // factors, so checking ν too is a real early exit, not just symmetry.
        if lambda.size() != mu.size() + nu.size() || !lambda.contains(mu) || !lambda.contains(nu) {
            return 0;
        }
        // Filling λ/μ with content ν costs |ν| cells, and c^λ_{μν} = c^λ_{νμ},
        // so peel off the *larger* factor and fill the smaller shape. Without
        // this, c^λ_{μν} with |μ| = 15 and |ν| = 75 backtracks over 75 cells:
        // 99 ms, against 0.1 ms once the 15-cell side is chosen instead.
        let (inner, content) = if mu.size() >= nu.size() {
            (mu, nu)
        } else {
            (nu, mu)
        };
        lr_cached(lambda, inner, content, || {
            lr_coeff_uncached(lambda, inner, content)
        })
    }
}

fn lr_coeff_uncached(lambda: &Partition, mu: &Partition, nu: &Partition) -> u128 {
    {
        let maxval = nu.len();
        if maxval == 0 {
            // s_μ · s_∅ = s_μ.
            return if lambda == mu { 1 } else { 0 };
        }

        let rows = lambda.len();
        let width = lambda.part(0) as usize;
        // grid[r][c] = filled value (0 = unset). Only skew cells μ_r ≤ c < λ_r used.
        let mut grid = vec![vec![0u32; width]; rows];

        // Fill cells in *reading-word order*: top→bottom, right→left. In this
        // order the right neighbor (row constraint), the top neighbor (column
        // constraint), and the ballot prefix are all already determined, so
        // every constraint prunes immediately — no post-fill lattice scan.
        let mut cells = Vec::new();
        for r in 0..rows {
            let start = mu.part(r) as usize;
            let end = lambda.part(r) as usize;
            for c in (start..end).rev() {
                cells.push((r, c));
            }
        }

        let cap: Vec<u32> = nu.parts().to_vec();
        let mut used = vec![0u32; maxval]; // content used so far, per value
        let mut cnt = vec![0i64; maxval + 2]; // ballot running counts, in reading order
        let mut count: u128 = 0;

        let mut st = State {
            cells: &cells,
            grid: &mut grid,
            lambda,
            mu,
            maxval,
            cap: &cap,
            used: &mut used,
            cnt: &mut cnt,
            count: &mut count,
        };
        backtrack(0, &mut st);
        count
    }
}

/// Backtracking state, grouped so the recursion takes a single borrow.
struct State<'a> {
    cells: &'a [(usize, usize)],
    grid: &'a mut [Vec<u32>],
    lambda: &'a Partition,
    mu: &'a Partition,
    maxval: usize,
    cap: &'a [u32],
    used: &'a mut [u32],
    cnt: &'a mut [i64],
    count: &'a mut u128,
}

fn backtrack(idx: usize, st: &mut State) {
    if idx == st.cells.len() {
        *st.count += 1;
        return;
    }
    let (r, c) = st.cells[idx];
    let end_r = st.lambda.part(r) as usize;

    // Right neighbor (row weakly increasing, filled right→left): if present, the
    // entry here must be ≤ it. Absent ⇒ no upper bound (use maxval).
    let right = if c + 1 < end_r {
        st.grid[r][c + 1]
    } else {
        st.maxval as u32
    };
    // Top neighbor (column strictly increasing): the cell above, if it is a skew
    // cell. Absent ⇒ no lower bound (use 0).
    let top = if r > 0 && c >= st.mu.part(r - 1) as usize && c < st.lambda.part(r - 1) as usize {
        st.grid[r - 1][c]
    } else {
        0
    };

    for v in 1..=st.maxval as u32 {
        let vi = (v - 1) as usize;
        if st.used[vi] >= st.cap[vi] {
            continue; // content exhausted for this value
        }
        if v > right {
            continue; // row must be weakly increasing
        }
        if v <= top {
            continue; // column must be strictly increasing
        }
        // Ballot check: appending v keeps #(v-1) ≥ #v (only that pair can break).
        // `cnt` is indexed by value, so compare cnt[v] (after +1) against cnt[v-1].
        if v >= 2 && st.cnt[v as usize] + 1 > st.cnt[(v - 1) as usize] {
            continue;
        }

        st.grid[r][c] = v;
        st.used[vi] += 1;
        st.cnt[v as usize] += 1;
        backtrack(idx + 1, st);
        st.cnt[v as usize] -= 1;
        st.used[vi] -= 1;
        st.grid[r][c] = 0;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn p(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    #[test]
    fn pieri_s1_s1() {
        // s_1 · s_1 = s_2 + s_{11}
        let prod = NaiveLr.schur_product(&p(&[1]), &p(&[1]));
        let mut got: Vec<_> = prod.iter().map(|(l, c)| (l.parts().to_vec(), *c)).collect();
        got.sort();
        assert_eq!(got, vec![(vec![1, 1], 1), (vec![2], 1)]);
    }

    /// LR entries range over `1..=ℓ(ν)`, so the entry type must hold more than
    /// 255. It used to be `u8` with the loop bound written `st.maxval as u8`,
    /// which wrapped: at ℓ(ν) = 256 the bound became 0 and `lr_coeff` returned
    /// 0 instantly — a fast, silent wrong answer, while ℓ(ν) = 255 was correct.
    ///
    /// `s_∅ · s_ν = s_ν`, so every coefficient below must be 1. These stay cheap
    /// despite the 300 cells because the ballot condition forces the fill: in
    /// reading order the entry at each cell is pinned to one value, so the
    /// search never branches.
    #[test]
    fn value_range_extends_past_the_u8_boundary() {
        let empty = Partition::default();
        for len in [255usize, 256, 300] {
            let tall = Partition::new(std::iter::repeat_n(1, len));
            assert_eq!(
                NaiveLr.lr_coeff(&tall, &empty, &tall),
                1,
                "c^(1^{len})_(∅,1^{len}) with ℓ(ν) = {len}"
            );
        }
    }

    #[test]
    fn s21_squared_has_multiplicity_two() {
        // s_{21}·s_{21} = s_{42}+s_{411}+s_{33}+2·s_{321}+s_{3111}+s_{222}+s_{2211}
        let coeff = NaiveLr.lr_coeff(&p(&[3, 2, 1]), &p(&[2, 1]), &p(&[2, 1]));
        assert_eq!(coeff, 2, "c^{{321}}_{{21,21}} should be 2");

        let prod = NaiveLr.schur_product(&p(&[2, 1]), &p(&[2, 1]));
        let total: u128 = prod.iter().map(|(_, c)| *c).sum();
        assert_eq!(total, 8, "sum of LR coeffs = f^{{21}} squared count");
        // every term has degree 6
        assert!(prod.iter().all(|(l, _)| l.size() == 6));
    }
}
