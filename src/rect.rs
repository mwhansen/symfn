//! Okada's closed form for a product of two **rectangles**.
//!
//! When both factors are rectangles the Littlewood–Richardson expansion is
//! multiplicity-free and its support has an explicit description, so the whole
//! product can be *generated* rather than searched for. Okada (1998): writing
//! `(aᵖ)` for the partition with `p` parts equal to `a`, and taking `p ≥ q`,
//! the partitions occurring in `s_{(aᵖ)}·s_{(bᑫ)}` are exactly those λ with
//! ℓ(λ) ≤ p+q such that
//!
//! ```text
//!   λ_{q+1} = λ_{q+2} = … = λ_p = a
//!   λ_q ≥ max(a, b)
//!   λ_i + λ_{p+q+1-i} = a + b        for i = 1, …, q
//! ```
//!
//! all with coefficient 1.
//!
//! ## Why this enumerates without backtracking
//!
//! The conditions leave only λ₁ … λ_q free — the middle block is pinned to `a`
//! and the last `q` parts are pinned to `a+b−λ_i` — and *every* remaining
//! constraint collapses into a single chain:
//!
//! - the middle block needs λ_q ≥ a and a ≥ a+b−λ_q, i.e. λ_q ≥ b; when p = q
//!   the block is empty and the two neighbours ask 2λ_q ≥ a+b. All three follow
//!   from λ_q ≥ max(a,b), which is already required.
//! - the last block is `a+b−λ_q, …, a+b−λ_1`, so it decreases exactly when
//!   λ₁ ≥ … ≥ λ_q does. It stays non-negative exactly when λ₁ ≤ a+b.
//!
//! So λ occurs **iff** `a+b ≥ λ₁ ≥ λ₂ ≥ … ≥ λ_q ≥ max(a,b)`, and the support is
//! in bijection with the partitions in a `q × min(a,b)` box — `C(q+min(a,b), q)`
//! of them. Nothing is generated that must later be discarded, and |λ| = pa+qb
//! comes out automatically: the `q` paired conditions contribute `q(a+b)` and
//! the middle block `(p−q)a`.
//!
//! Verified against [`SkewLr`](crate::skew_lr::SkewLr) in the tests below, and
//! this is the case the general engine handles *worst* — the conjugate-dispatch
//! heuristic in `skew_lr` records rectangles as a known loss. Measured on
//! `s(12⁶)·s(12⁶)` (18564 terms): 37.0 ms via the DP, 3.8 ms here.

use crate::partition::Partition;

/// `Some((a, p))` if `p` is a non-empty partition with all parts equal to `a`.
fn as_rectangle(lambda: &Partition) -> Option<(u32, usize)> {
    let a = *lambda.parts().first()?;
    lambda
        .parts()
        .iter()
        .all(|&x| x == a)
        .then(|| (a, lambda.len()))
}

/// Both factors as rectangles, ordered so the first has at least as many rows.
fn both_rectangles(mu: &Partition, nu: &Partition) -> Option<(u32, usize, u32, usize)> {
    let (a, p) = as_rectangle(mu)?;
    let (b, q) = as_rectangle(nu)?;
    Some(if p >= q { (a, p, b, q) } else { (b, q, a, p) })
}

/// `s_μ · s_ν` when both μ and ν are rectangles, or `None` if either is not.
///
/// The result matches [`LrBackend::schur_product`](crate::lr::LrBackend::schur_product):
/// sorted by partition, coefficients non-zero. Every coefficient here is 1.
pub fn okada_product(mu: &Partition, nu: &Partition) -> Option<Vec<(Partition, u128)>> {
    let (a, p, b, q) = both_rectangles(mu, nu)?;
    let mut out = Vec::new();
    let mut head = vec![0u32; q];
    extend(0, a + b, &mut head, a, b, p, q, &mut out);
    out.sort_by(|x: &(Partition, u128), y| x.0.cmp(&y.0));
    Some(out)
}

/// Choose λ_{i+1} … λ_q weakly decreasing in `[max(a,b), hi]`, then emit.
#[allow(clippy::too_many_arguments)]
fn extend(
    i: usize,
    hi: u32,
    head: &mut Vec<u32>,
    a: u32,
    b: u32,
    p: usize,
    q: usize,
    out: &mut Vec<(Partition, u128)>,
) {
    if i == q {
        out.push((assemble(head, a, b, p, q), 1));
        return;
    }
    for v in (a.max(b)..=hi).rev() {
        head[i] = v;
        extend(i + 1, v, head, a, b, p, q, out);
    }
}

/// Build λ from its free first `q` parts: middle block `a`, tail `a+b−λ_i`.
fn assemble(head: &[u32], a: u32, b: u32, p: usize, q: usize) -> Partition {
    let mut full = vec![0u32; p + q];
    full[..q].copy_from_slice(head);
    full[q..p].fill(a);
    for (i, &v) in head.iter().enumerate() {
        full[p + q - 1 - i] = a + b - v;
    }
    // Weakly decreasing, so trailing zeros (if λ₁ = a+b) are a suffix.
    while full.last() == Some(&0) {
        full.pop();
    }
    Partition::new(full)
}

/// `c^λ_{μν}` when both μ and ν are rectangles, or `None` if either is not.
///
/// This tests Okada's conditions on λ directly — O(ℓ(λ)), no enumeration. The
/// conditions imply |λ| = |μ|+|ν|, so a λ of the wrong size is rejected by them
/// rather than needing a separate guard.
pub fn okada_coeff(lambda: &Partition, mu: &Partition, nu: &Partition) -> Option<u128> {
    let (a, p, b, q) = both_rectangles(mu, nu)?;
    if lambda.len() > p + q {
        return Some(0);
    }
    // λ_i for 1-based i, reading past the end as 0.
    let at = |i: usize| lambda.parts().get(i - 1).copied().unwrap_or(0);
    let occurs = (q + 1..=p).all(|i| at(i) == a)
        && at(q) >= a.max(b)
        && (1..=q).all(|i| at(i) + at(p + q + 1 - i) == a + b);
    Some(u128::from(occurs))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lr::LrBackend;
    use crate::skew_lr::SkewLr;

    fn rect(rows: usize, part: u32) -> Partition {
        Partition::new(std::iter::repeat(part).take(rows))
    }

    /// The whole point: agree with the general engine on every small rectangle
    /// pair, in both argument orders (the `p ≥ q` normalisation must not care).
    #[test]
    fn agrees_with_skew_lr_on_all_small_rectangle_pairs() {
        let mut cases = 0;
        for a in 1..=4u32 {
            for p in 1..=4usize {
                for b in 1..=4u32 {
                    for q in 1..=4usize {
                        let (mu, nu) = (rect(p, a), rect(q, b));
                        let want = SkewLr.schur_product(&mu, &nu);
                        assert_eq!(
                            okada_product(&mu, &nu).as_ref(),
                            Some(&want),
                            "s({a}^{p}) · s({b}^{q})"
                        );
                        assert_eq!(
                            okada_product(&nu, &mu).as_ref(),
                            Some(&SkewLr.schur_product(&nu, &mu)),
                            "s({b}^{q}) · s({a}^{p}) (swapped)"
                        );
                        cases += 1;
                    }
                }
            }
        }
        assert_eq!(cases, 256);
    }

    /// `okada_coeff` must agree with the engine on *every* λ of the right size,
    /// including the zeros — the enumeration test above only sees the support.
    #[test]
    fn coefficients_including_zeros_agree_with_skew_lr() {
        for a in 1..=3u32 {
            for p in 1..=3usize {
                for b in 1..=3u32 {
                    for q in 1..=3usize {
                        let (mu, nu) = (rect(p, a), rect(q, b));
                        let n = a * p as u32 + b * q as u32;
                        for lambda in crate::partition::partitions_of(n) {
                            assert_eq!(
                                okada_coeff(&lambda, &mu, &nu),
                                Some(SkewLr.lr_coeff(&lambda, &mu, &nu)),
                                "c^{lambda}_{{({a}^{p}),({b}^{q})}}"
                            );
                        }
                    }
                }
            }
        }
    }

    /// A λ of the wrong size must be rejected by the conditions themselves.
    #[test]
    fn wrong_size_lambda_is_zero() {
        let (mu, nu) = (rect(2, 3), rect(2, 3)); // |μ|+|ν| = 12
        for lambda in crate::partition::partitions_of(11) {
            assert_eq!(okada_coeff(&lambda, &mu, &nu), Some(0), "λ = {lambda}");
        }
        for lambda in crate::partition::partitions_of(13) {
            assert_eq!(okada_coeff(&lambda, &mu, &nu), Some(0), "λ = {lambda}");
        }
    }

    /// Non-rectangles (and the empty partition) must decline, not answer wrongly.
    #[test]
    fn declines_when_a_factor_is_not_a_rectangle() {
        let r = rect(2, 2);
        for other in [Partition::new([3u32, 1]), Partition::default()] {
            assert_eq!(okada_product(&r, &other), None, "{other}");
            assert_eq!(okada_product(&other, &r), None, "{other}");
            assert_eq!(okada_coeff(&Partition::new([4u32, 3]), &r, &other), None);
        }
    }

    /// Rows and columns are rectangles too, so the Pieri cases route through here.
    #[test]
    fn pieri_cases_are_rectangles() {
        // s_{(3)} · s_{(1^2)}: a horizontal times a vertical strip.
        let got = okada_product(&Partition::new([3u32]), &rect(2, 1)).expect("rectangles");
        assert_eq!(
            got,
            SkewLr.schur_product(&Partition::new([3u32]), &rect(2, 1))
        );
    }

    /// The support is the set of partitions in a `q × min(a,b)` box.
    #[test]
    fn support_size_is_the_expected_binomial() {
        for (a, p, b, q) in [
            (4u32, 4usize, 4u32, 4usize),
            (8, 5, 8, 5),
            (12, 6, 12, 6),
            (3, 5, 7, 2),
        ] {
            let n = okada_product(&rect(p, a), &rect(q, b))
                .expect("rectangles")
                .len();
            let (m, k) = (a.min(b) as u128, q.min(p) as u128);
            let want = (1..=k).fold(1u128, |acc, i| acc * (m + i) / i);
            assert_eq!(n as u128, want, "s({a}^{p})·s({b}^{q})");
        }
    }
}
