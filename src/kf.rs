//! Kostka–Foulkes polynomials `K_{λμ}(t)`, from the Hall–Littlewood transition.
//!
//! `Q'_μ = Σ_λ K_{λμ}(t) s_λ`, so the polynomials are already the coefficients
//! of [`hall_littlewood`](crate::hall_littlewood) and this module is the entry
//! point that says so. Symmetrica has no Kostka–Foulkes function at all — its
//! `hall_littlewood` is the only way to reach these, and the transition has to
//! be read off by hand (`docs/record/hall-littlewood.md`). Sage reaches the
//! same polynomials through `KostkaFoulkesPolynomial` in
//! `sage.combinat.sf.kfpoly`, and `scripts/check_kf.py` is the comparison.
//!
//! `t = 1` recovers the ordinary Kostka number, so [`kostka_foulkes_table`] is
//! the t-analog of [`kostka_table`](crate::kostka::kostka_table) and is
//! indexed the same way.
//!
//! ## Two routes, deliberately
//!
//! [`kostka_foulkes_by_charge`](crate::charge::kostka_foulkes_by_charge)
//! computes the same polynomials by enumerating semistandard tableaux and
//! summing `t^{charge}`. It is exponential and is not what anyone should call —
//! it is here because it shares no code with the recursion, which makes
//! agreement between the two evidence rather than tautology (V3,
//! `docs/policies/validation.md`). Same role [`NaiveLr`](crate::NaiveLr) plays
//! for Littlewood–Richardson.
//!
//! ## Cost
//!
//! A single `K_{λμ}(t)` costs a whole `Q'_μ`: the recursion produces the entire
//! column at once and there is no cheaper route to one entry of it. So asking
//! for a column, or for the table, is free relative to asking for one value —
//! and asking for p(n)² values one at a time would repeat each column p(n)
//! times, which is exactly the mistake [`kostka_table`](mod@crate::kostka)
//! documents.
//!
//! ## Range
//!
//! These polynomials *are* the coefficients
//! [`hall_littlewood`](crate::hall_littlewood) produces, so
//! the range is that module's: over a fixed-width `C` the call refuses rather
//! than wrapping past its wall, and at `i128` the whole-degree table gains
//! ~2.2 bits per degree, reaching 127 bits near n ≈ 67 — far past the degree
//! the table stops finishing at. Asking for one `K_{λμ}(t)` does not move that
//! wall, because one value costs a whole `Q'_μ` (above); what does move it is
//! the *shape*, and `μ = 1ⁿ` is the extremal case, where the coefficients are
//! the Kostka numbers `f^μ` and so are capped by `√(n!)`.
//!
//! Measurements and the harness are in
//! `docs/record/failure-and-overflow.md` (`examples/probe_qt_walls.rs`).
//!
//! ## References
//!
//! - **\[M\]** I. G. Macdonald, *Symmetric Functions and Hall Polynomials*,
//!   2nd ed., Oxford, 1995 — III.6, where `K(t) = M(s, P)` is defined by
//!   `s_λ = Σ_μ K_{λμ}(t) P_μ` (by way of III (2.6)), with `K(0)` the identity
//!   and `K(1)` the Kostka matrix; and III.5 Example 7(a), where `Q'_μ` is the
//!   basis dual to `P_μ` under the Hall inner product. Together these give
//!   `Q'_μ = Σ_λ K_{λμ}(t) s_λ`, the identity this module reads the
//!   polynomials off. III (6.5) is the charge formula
//!   [`kostka_foulkes_by_charge`](crate::charge::kostka_foulkes_by_charge)
//!   implements.

// A shape index.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use crate::coeff::Ring;
use crate::partition::Partition;
use crate::qt::QtPoly;
use crate::sym::SymFn;

/// `K_{λμ}(t)`.
///
/// Zero unless `|λ| = |μ|` and λ ⊵ μ; `K_{λλ}(t) = 1`.
///
/// This computes all of `Q'_μ` and reads one coefficient out of it. When more
/// than one λ is wanted for the same μ, take the column from
/// [`kostka_foulkes_column`] instead, and for a whole degree use
/// [`kostka_foulkes_table`].
pub fn kostka_foulkes<C: Ring>(lambda: &Partition, mu: &Partition) -> QtPoly<C> {
    if lambda.size() != mu.size() {
        return QtPoly::zero();
    }
    crate::hall_littlewood::<C>(mu).coeff(lambda)
}

/// Every `K_{λμ}(t)` for a fixed μ, paired with its λ, sorted by λ.
///
/// The natural unit of work: one `Q'_μ` *is* the column, so this costs what a
/// single [`kostka_foulkes`] costs. Only the λ with `K_{λμ} ≠ 0` appear, which
/// by triangularity are exactly those dominating μ.
pub fn kostka_foulkes_column<C: Ring>(mu: &Partition) -> Vec<(Partition, QtPoly<C>)> {
    crate::hall_littlewood::<C>(mu)
        .terms()
        .iter()
        .map(|(lambda, k)| (lambda.clone(), k.clone()))
        .collect()
}

/// The whole table of degree `n`, as `table[i][j] = K_{λⁱ λʲ}(t)` indexed
/// against `memo::partitions_cached` — the same
/// orientation as [`kostka_table`](crate::kostka::kostka_table).
///
/// Built from [`hall_littlewood_table`](crate::hall_littlewood_table), so the
/// recursion's shared suffixes are shared here too.
pub fn kostka_foulkes_table<C: Ring>(n: u32) -> Vec<Vec<QtPoly<C>>> {
    let parts = crate::memo::partitions_cached(n);
    let index: std::collections::HashMap<&Partition, usize> =
        parts.iter().enumerate().map(|(i, p)| (p, i)).collect();
    let mut table = vec![vec![QtPoly::zero(); parts.len()]; parts.len()];
    for (j, (mu, hl)) in crate::hall_littlewood_table::<C>(n).into_iter().enumerate() {
        crate::interrupt::poll();
        debug_assert_eq!(&mu, &parts[j], "hall_littlewood_table must share the order");
        for (lambda, k) in hl.terms() {
            table[index[lambda]][j] = k.clone();
        }
    }
    table
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::charge::kostka_foulkes_by_charge;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    type Q = QtPoly<i64>;

    /// Against the charge enumeration, which shares no code with the recursion.
    #[test]
    fn agrees_with_the_charge_enumeration() {
        for n in 0..=8u32 {
            let parts = crate::partitions_of(n);
            for lambda in &parts {
                for mu in &parts {
                    assert_eq!(
                        kostka_foulkes::<i64>(lambda, mu),
                        kostka_foulkes_by_charge(lambda, mu),
                        "K_{lambda}{mu}(t)"
                    );
                }
            }
        }
    }

    /// t = 1 is the ordinary Kostka number — the tie back to tested machinery.
    #[test]
    fn specializes_to_kostka_at_one() {
        for n in 0..=9u32 {
            let parts = crate::partitions_of(n);
            let table = kostka_foulkes_table::<i64>(n);
            let plain = crate::kostka::kostka_table(n);
            for i in 0..parts.len() {
                for j in 0..parts.len() {
                    assert_eq!(
                        table[i][j].eval(&0, &1) as u128,
                        plain[i][j],
                        "K_{}{}(1)",
                        parts[i],
                        parts[j]
                    );
                }
            }
        }
    }

    /// The table, the column and the single value must be the same numbers;
    /// only the amount of work differs.
    #[test]
    fn table_column_and_single_value_agree() {
        for n in 0..=8u32 {
            let parts = crate::partitions_of(n);
            let table = kostka_foulkes_table::<i64>(n);
            for (j, mu) in parts.iter().enumerate() {
                let column: std::collections::HashMap<_, _> =
                    kostka_foulkes_column::<i64>(mu).into_iter().collect();
                for (i, lambda) in parts.iter().enumerate() {
                    let want = &table[i][j];
                    assert_eq!(&kostka_foulkes::<i64>(lambda, mu), want);
                    let from_column = column.get(lambda).cloned().unwrap_or_else(QtPoly::zero);
                    assert_eq!(&from_column, want, "column at {lambda}, {mu}");
                }
            }
        }
    }

    #[test]
    fn known_values() {
        // Sage: KostkaFoulkesPolynomial([2,1],[1,1,1],t) = t^2 + t
        let k: Q = kostka_foulkes(&part(&[2, 1]), &part(&[1, 1, 1]));
        assert_eq!(k.coeff(0, 1), 1);
        assert_eq!(k.coeff(0, 2), 1);
        assert_eq!(k.len(), 2, "{k}");
        // K_{λλ} = 1, and mismatched sizes are zero rather than an error
        let l = part(&[3, 2, 1]);
        assert_eq!(kostka_foulkes::<i64>(&l, &l), <Q as Ring>::one());
        assert_eq!(kostka_foulkes::<i64>(&l, &part(&[2, 1])), QtPoly::zero());
    }

    /// K_{λμ}(0) = δ_{λμ}: the Hall–Littlewood basis degenerates to Schur,
    /// which is the same statement as unitriangularity of the transition.
    #[test]
    fn at_t_zero_the_transition_is_the_identity() {
        for n in 0..=8u32 {
            let parts = crate::partitions_of(n);
            let table = kostka_foulkes_table::<i64>(n);
            for i in 0..parts.len() {
                for j in 0..parts.len() {
                    assert_eq!(table[i][j].eval(&0, &0), i64::from(i == j), "at t = 0");
                }
            }
        }
    }
}
