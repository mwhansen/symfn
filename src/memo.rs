//! Memoization of the pure combinatorial kernels.
//!
//! Every table here caches the result of a **pure function** of its key, so the
//! library stays referentially transparent and thread-safe: two threads racing
//! on the same key compute the same value, and the worst case is duplicated
//! work, never a wrong answer. That is categorically different from Symmetrica's
//! mutable global scratch state — nothing here is observable in the results.
//!
//! Caches are unbounded, which suits interactive research (degrees stay modest
//! and reuse is high). [`clear_caches`] releases them if a long-running session
//! wants the memory back.

use std::collections::HashMap;
use std::hash::Hash;
use std::sync::{Arc, OnceLock, RwLock};

use crate::partition::{partitions_of, Partition};

type Table<K, V> = RwLock<HashMap<K, V>>;

/// Look up `key`, computing and inserting it on a miss.
///
/// The read guard is released before `compute` runs, so a `compute` that
/// recurses back into the same table cannot deadlock.
fn lookup<K, V>(table: &Table<K, V>, key: &K, compute: impl FnOnce() -> V) -> V
where
    K: Eq + Hash + Clone,
    V: Clone,
{
    if let Some(v) = table.read().unwrap().get(key) {
        return v.clone();
    }
    let v = compute();
    table.write().unwrap().insert(key.clone(), v.clone());
    v
}

macro_rules! table {
    ($name:ident, $key:ty, $val:ty) => {
        fn $name() -> &'static Table<$key, $val> {
            static T: OnceLock<Table<$key, $val>> = OnceLock::new();
            T.get_or_init(|| RwLock::new(HashMap::new()))
        }
    };
}

table!(partitions_table, u32, Arc<Vec<Partition>>);
table!(character_table, (Partition, Partition), i64);
table!(kostka_table, (Partition, Partition), u128);
table!(lr_table, (Partition, Partition, Partition), u128);
table!(inverse_kostka_table, u32, Arc<(Vec<Partition>, Vec<Vec<i128>>)>);
table!(product_table, (Partition, Partition), Arc<Vec<(Partition, u128)>>);
table!(skew_table, (Partition, Partition), Arc<Vec<(Partition, u128)>>);

/// The partitions of `n`, shared rather than regenerated.
///
/// This is the broadest win: conversions and the coproduct enumerate
/// `partitions_of` inside nested loops, so the same list was being rebuilt
/// thousands of times.
pub fn partitions_cached(n: u32) -> Arc<Vec<Partition>> {
    lookup(partitions_table(), &n, || Arc::new(partitions_of(n)))
}

/// Memoized χ^λ(μ). The Murnaghan–Nakayama recursion has heavily overlapping
/// subproblems, so this is the difference between exponential and near-linear
/// on repeated use.
pub fn character_cached(lambda: &Partition, mu: &Partition, compute: impl FnOnce() -> i64) -> i64 {
    lookup(character_table(), &(lambda.clone(), mu.clone()), compute)
}

/// Memoized Kostka number K_{λμ}.
pub fn kostka_cached(lambda: &Partition, mu: &Partition, compute: impl FnOnce() -> u128) -> u128 {
    lookup(kostka_table(), &(lambda.clone(), mu.clone()), compute)
}

/// Memoized Littlewood–Richardson coefficient c^λ_{μν}.
pub fn lr_cached(
    lambda: &Partition,
    mu: &Partition,
    nu: &Partition,
    compute: impl FnOnce() -> u128,
) -> u128 {
    lookup(
        lr_table(),
        &(lambda.clone(), mu.clone(), nu.clone()),
        compute,
    )
}

/// Memoized inverse-Kostka data for degree `n` (partition order + inverse
/// matrix). Building it is O(p(n)²) tableau counts, so caching per degree
/// matters a lot for repeated m → s conversions.
pub fn inverse_kostka_cached(
    n: u32,
    compute: impl FnOnce() -> (Vec<Partition>, Vec<Vec<i128>>),
) -> Arc<(Vec<Partition>, Vec<Vec<i128>>)> {
    lookup(inverse_kostka_table(), &n, || Arc::new(compute()))
}

/// Memoized full expansion of s_μ · s_ν. Caching the whole product (rather than
/// individual coefficients) is what makes repeated multiplication cheap.
pub fn product_cached(
    mu: &Partition,
    nu: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, u128)>,
) -> Arc<Vec<(Partition, u128)>> {
    lookup(product_table(), &(mu.clone(), nu.clone()), || {
        Arc::new(compute())
    })
}

/// Memoized full expansion of the skew Schur function s_{outer/inner}.
///
/// Keyed on the *shape* rather than on (λ, μ, ν), because one traversal of the
/// shape produces every ν at once. A caller sweeping over ν — `skew_schur`, the
/// coproduct, a single `lr_coeff` in a loop — therefore pays for exactly one
/// enumeration instead of p(n) of them.
pub fn skew_cached(
    outer: &Partition,
    inner: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, u128)>,
) -> Arc<Vec<(Partition, u128)>> {
    lookup(skew_table(), &(outer.clone(), inner.clone()), || {
        Arc::new(compute())
    })
}

/// Drop every cached table, releasing the memory.
pub fn clear_caches() {
    partitions_table().write().unwrap().clear();
    character_table().write().unwrap().clear();
    kostka_table().write().unwrap().clear();
    lr_table().write().unwrap().clear();
    inverse_kostka_table().write().unwrap().clear();
    product_table().write().unwrap().clear();
    skew_table().write().unwrap().clear();
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn partitions_are_shared_and_correct() {
        let a = partitions_cached(6);
        let b = partitions_cached(6);
        assert_eq!(a.len(), 11); // p(6) = 11
        assert!(Arc::ptr_eq(&a, &b), "second call should reuse the cached Arc");
    }

    #[test]
    fn lookup_computes_once_then_reuses() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        // A key no other test touches: `clear_caches` is global, and the test
        // harness runs tests in parallel, so isolation comes from the key.
        let calls = AtomicUsize::new(0);
        let lam = Partition::new([9, 7, 5, 3, 1]);
        let mu = Partition::new([5, 5, 5, 5, 5]);
        for _ in 0..5 {
            kostka_cached(&lam, &mu, || {
                calls.fetch_add(1, Ordering::SeqCst);
                12345
            });
        }
        assert_eq!(calls.load(Ordering::SeqCst), 1, "should compute exactly once");
    }
}
