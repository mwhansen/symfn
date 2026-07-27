//! Kostka numbers K_{λμ} = number of semistandard Young tableaux of shape λ and
//! content μ (μ_i entries equal to i, rows weakly increasing, columns strictly
//! increasing).
//!
//! These are the s → m transition: s_λ = Σ_μ K_{λμ} m_μ. The reverse, m → s,
//! inverts the (unitriangular) Kostka matrix; see `convert`.

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
    let maxval = mu.len();
    if maxval == 0 {
        return 0;
    }
    let rows = lambda.len();
    let width = lambda.part(0) as usize;
    let mut grid = vec![vec![0u32; width]; rows];

    // Row-major fill order (top→bottom, left→right): left and top SSYT
    // neighbours are always already placed.
    let cells: Vec<(usize, usize)> = (0..rows)
        .flat_map(|r| (0..lambda.part(r) as usize).map(move |c| (r, c)))
        .collect();

    let cap: Vec<u32> = mu.parts().to_vec();
    let mut used = vec![0u32; maxval];
    let mut count = 0u128;
    count_ssyt(0, &cells, &mut grid, maxval, &cap, &mut used, &mut count);
    count
}

#[allow(clippy::too_many_arguments)]
fn count_ssyt(
    idx: usize,
    cells: &[(usize, usize)],
    grid: &mut [Vec<u32>],
    maxval: usize,
    cap: &[u32],
    used: &mut [u32],
    count: &mut u128,
) {
    if idx == cells.len() {
        *count += 1;
        return;
    }
    let (r, c) = cells[idx];
    let left = if c > 0 { grid[r][c - 1] } else { 0 };
    // Because the shape is a partition, if (r,c) exists then so does (r-1,c).
    let top = if r > 0 { grid[r - 1][c] } else { 0 };

    for v in 1..=maxval as u32 {
        let vi = (v - 1) as usize;
        if used[vi] >= cap[vi] {
            continue;
        }
        if v < left {
            continue; // rows weakly increasing
        }
        if v <= top {
            continue; // columns strictly increasing (top == 0 ⇒ no constraint)
        }
        grid[r][c] = v;
        used[vi] += 1;
        count_ssyt(idx + 1, cells, grid, maxval, cap, used, count);
        used[vi] -= 1;
        grid[r][c] = 0;
    }
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

    /// Tableau entries are values in `1..=ℓ(μ)`, so the entry type must hold
    /// more than 255. It used to be `u8` with the loop bound written
    /// `maxval as u8`, which *wrapped*: at ℓ(μ) = 256 the bound became 0, the
    /// value loop ran empty, and `kostka` returned 0 instantly instead of the
    /// true count — silently, with no overflow anywhere to notice.
    ///
    /// This is white-box (it drives `count_ssyt` directly) because the public
    /// entry point cannot be pushed past 255 cheaply: μ with 256 parts forces
    /// |λ| ≥ 256, and unlike the LR fill there is no ballot condition to prune,
    /// so the honest count is exponential. Driving the counter with a sparse
    /// `cap` keeps the *shape* at two cells while putting the live values up at
    /// 299 and 300 — exactly the range the old bound could not reach.
    #[test]
    fn value_range_extends_past_the_u8_boundary() {
        let maxval = 300;
        // A single column of two cells: strictly increasing, so the only fill
        // using the two available values is (299, 300) and the count is 1.
        let cells = vec![(0usize, 0usize), (1, 0)];
        let mut grid = vec![vec![0u32; 1]; 2];
        let mut cap = vec![0u32; maxval];
        cap[298] = 1; // value 299
        cap[299] = 1; // value 300
        let mut used = vec![0u32; maxval];
        let mut count = 0u128;
        count_ssyt(0, &cells, &mut grid, maxval, &cap, &mut used, &mut count);
        assert_eq!(count, 1, "values above 255 must be reachable");
    }

    #[test]
    fn diagonal_is_one() {
        for parts in [&[2, 1][..], &[3, 1], &[2, 2], &[4]] {
            let lam = p(parts);
            assert_eq!(kostka(&lam, &lam), 1, "K_{{{lam},{lam}}}");
        }
    }
}
