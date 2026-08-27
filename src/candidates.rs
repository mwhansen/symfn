//! The candidate walk shared by the two counting routes, and its parallel
//! driver.
//!
//! [`crate::two_row`] and [`crate::three_row`] compute `s_μ · s_ν` by
//! enumerating every λ the chain of horizontal strips admits and counting each
//! one's fibre. The enumeration is the same walk with a different depth — a
//! chain of `strips` strips can open at most `strips` new rows and bounds row j
//! of λ by μ_{j−strips} — and the fibre counts are independent of one another,
//! so the walk is where the routes parallelize: the candidates are split by
//! their first two rows into work items, and each worker counts its items with
//! scratch of its own. Nothing here knows what a fibre is; the counter comes
//! from the caller.

use std::sync::atomic::{AtomicUsize, Ordering};

use crate::partition::Partition;

/// Work items per worker below which the walk stays on one thread. A worker
/// costs a spawn and its own scratch, tens of microseconds, and an item is the
/// subtree under one two-row prefix; a product with a few dozen items is a
/// few dozen fibre counts and finishes before a second thread would start.
const ITEMS_PER_WORKER: usize = 32;

/// Every prefix of a candidate λ, called at depth `stop` with the prefix and
/// the cells still to place.
///
/// The candidates are the λ ⊇ μ with |λ| = |μ| + `total` that a chain of
/// `strips` horizontal strips can reach: at most `strips` new rows, and
/// λⱼ ≤ μ_{j−strips} because λⱼ ≤ λ^{s−1}ⱼ₋₁ ≤ … ≤ μⱼ₋ₛ through the strips.
/// `cand` has `rows = ℓ(μ) + strips` slots; the walk fills `cand[..stop]` in
/// lexicographic order and leaves it zeroed on return. Called with `stop` equal
/// to `rows` it visits every complete candidate, and the leaf's `left` is zero
/// exactly when the candidate has the right size.
// The eight parameters are the recursion state the doc comment names; a
// struct would carry the same eight fields.
#[allow(clippy::too_many_arguments)]
pub(crate) fn walk(
    j: usize,
    stop: usize,
    left: u32,
    prev: u32,
    mu: &[u32],
    strips: usize,
    cand: &mut Vec<u32>,
    leaf: &mut impl FnMut(&[u32], u32),
) {
    if j == stop {
        leaf(&cand[..stop], left);
        return;
    }
    let at = |i: usize| mu.get(i).copied().unwrap_or(0);
    let lo = at(j);
    let hi = if j < strips {
        prev.min(left + lo)
    } else {
        at(j - strips).min(prev).min(left + lo)
    };
    if hi < lo {
        return;
    }
    for v in lo..=hi {
        cand[j] = v;
        walk(j + 1, stop, left - (v - lo), v, mu, strips, cand, leaf);
    }
    cand[j] = 0;
}

/// The nonzero terms of the product, sorted by λ: every candidate the chain
/// admits, counted by a counter from `make_counter`.
///
/// `make_counter` is called once per worker and returns that worker's fibre
/// counter, which owns whatever scratch it needs; a candidate is passed to it
/// as its parts without trailing zeros. `threads` is a cap: never more than
/// one worker per [`ITEMS_PER_WORKER`] work items, so a small product runs on
/// the calling thread alone.
/// The output does not depend on the thread count: candidates are counted
/// independently and sorted at the end.
///
/// A worker's panic is re-raised on the calling thread with its payload, so a
/// counter's own `# Panics` contract holds unchanged.
pub(crate) fn count_all<F>(
    mu: &[u32],
    strips: usize,
    total: u32,
    threads: usize,
    make_counter: impl Fn() -> F + Sync,
) -> Vec<(Partition, u128)>
where
    F: FnMut(&[u32]) -> u128,
{
    let rows = mu.len() + strips;
    let split = rows.min(2);
    // Work items: the first two rows of every candidate, with the cells left
    // to place below them. Lexicographic, so the deepest subtrees — small first
    // rows, most cells still to place — are claimed first.
    let mut items: Vec<([u32; 2], u32)> = Vec::new();
    {
        let mut cand = vec![0u32; rows];
        walk(
            0,
            split,
            total,
            u32::MAX,
            mu,
            strips,
            &mut cand,
            &mut |c, left| {
                items.push(([c[0], c.get(1).copied().unwrap_or(0)], left));
            },
        );
    }
    let threads = threads.min(items.len() / ITEMS_PER_WORKER).max(1);

    // Count everything below one work item.
    let run = |counter: &mut F,
               cand: &mut Vec<u32>,
               out: &mut Vec<(Partition, u128)>,
               item: ([u32; 2], u32)| {
        let ([v0, v1], left) = item;
        cand[..split].copy_from_slice(&[v0, v1][..split]);
        let prev = if split == 0 {
            u32::MAX
        } else {
            cand[split - 1]
        };
        walk(split, rows, left, prev, mu, strips, cand, &mut |c, left| {
            if left != 0 {
                return;
            }
            let end = c.iter().rposition(|&x| x > 0).map_or(0, |i| i + 1);
            let k = counter(&c[..end]);
            if k > 0 {
                out.push((Partition::new(c[..end].iter().copied()), k));
            }
        });
        // `walk` zeroes what it filled; the prefix is ours to clear.
        cand[..split].fill(0);
    };

    let mut out = if threads <= 1 {
        let mut counter = make_counter();
        let mut cand = vec![0u32; rows];
        let mut out = Vec::new();
        for &item in &items {
            run(&mut counter, &mut cand, &mut out, item);
        }
        out
    } else {
        let cursor = AtomicUsize::new(0);
        let parts: Vec<Vec<(Partition, u128)>> = std::thread::scope(|s| {
            let handles: Vec<_> = (0..threads)
                .map(|_| {
                    let (cursor, items, run, make_counter) = (&cursor, &items, &run, &make_counter);
                    s.spawn(move || {
                        let mut counter = make_counter();
                        let mut cand = vec![0u32; rows];
                        let mut out = Vec::new();
                        loop {
                            let i = cursor.fetch_add(1, Ordering::Relaxed);
                            if i >= items.len() {
                                break;
                            }
                            run(&mut counter, &mut cand, &mut out, items[i]);
                        }
                        out
                    })
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().unwrap_or_else(|e| std::panic::resume_unwind(e)))
                .collect()
        });
        parts.into_iter().flatten().collect()
    };
    out.sort_by(|x, y| x.0.cmp(&y.0));
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The walk visits exactly the partitions the bounds describe, once each,
    /// in lexicographic order — checked against a brute-force filter over
    /// every partition of the right size.
    #[test]
    fn walk_visits_each_admissible_candidate_once_in_order() {
        use crate::partition::partitions_of;
        for (mu, strips, total) in [
            (vec![3u32, 1], 2usize, 4u32),
            (vec![4, 2, 1], 3, 5),
            (vec![2, 2], 2, 6),
            (vec![], 3, 4),
            (vec![5], 2, 1),
        ] {
            let rows = mu.len() + strips;
            let mut seen = Vec::new();
            let mut cand = vec![0u32; rows];
            walk(
                0,
                rows,
                total,
                u32::MAX,
                &mu,
                strips,
                &mut cand,
                &mut |c, left| {
                    if left == 0 {
                        let end = c.iter().rposition(|&x| x > 0).map_or(0, |i| i + 1);
                        seen.push(c[..end].to_vec());
                    }
                },
            );
            assert!(cand.iter().all(|&x| x == 0), "walk left cand dirty");
            let n = mu.iter().sum::<u32>() + total;
            let at = |v: &[u32], i: usize| v.get(i).copied().unwrap_or(0);
            let mut want: Vec<Vec<u32>> = partitions_of(n)
                .into_iter()
                .map(|p| p.parts().to_vec())
                .filter(|lam| {
                    // Every row, including the ones λ leaves empty: containing
                    // μ is a condition there too.
                    lam.len() <= rows
                        && (0..rows).all(|j| {
                            at(lam, j) >= at(&mu, j)
                                && (j < strips || at(lam, j) <= at(&mu, j - strips))
                        })
                })
                .collect();
            want.sort();
            assert_eq!(seen, want, "μ = {mu:?}, {strips} strips, {total} cells");
        }
    }

    /// One thread or several, the same terms in the same order.
    #[test]
    fn count_all_does_not_depend_on_the_thread_count() {
        // A counter that is 1 on candidates with an even first part and the
        // number of parts otherwise: nonzero, distinguishable, deterministic.
        let make = || {
            |lam: &[u32]| {
                u128::from(if lam[0].is_multiple_of(2) {
                    1
                } else {
                    lam.len() as u32
                })
            }
        };
        let mu = vec![20u32, 16, 12];
        let one = count_all(&mu, 3, 48, 1, make);
        let many = count_all(&mu, 3, 48, 4, make);
        assert_eq!(one, many);
        assert!(one.len() > 1000, "expected a real sweep, got {}", one.len());
        assert!(
            one.windows(2).all(|w| w[0].0 < w[1].0),
            "terms must be sorted and distinct"
        );
    }

    /// A worker's panic surfaces on the caller with its own message.
    #[test]
    fn a_worker_panic_is_re_raised_with_its_payload() {
        let make = || {
            |lam: &[u32]| -> u128 {
                if lam[0] == 30 {
                    panic!("counter refused λ₁ = 30");
                }
                1
            }
        };
        let mu = vec![20u32, 16, 12];
        let err = std::panic::catch_unwind(|| count_all(&mu, 3, 48, 4, make)).unwrap_err();
        let msg = err
            .downcast_ref::<&str>()
            .map(|s| (*s).to_string())
            .or_else(|| err.downcast_ref::<String>().cloned())
            .unwrap_or_default();
        assert_eq!(msg, "counter refused λ₁ = 30");
    }
}
