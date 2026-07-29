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
use crate::bh::Rat;
use crate::guard::GuardedRat;
use crate::qt::QtPoly;
use crate::sym::{PowerSum, Schur};

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
table!(character_table, (Partition, Partition), i128);
table!(kostka_table, (Partition, Partition), u128);
table!(lr_table, (Partition, Partition, Partition), u128);
table!(lex_parts_table, u32, Arc<Vec<Partition>>);
table!(inverse_kostka_row_table, Partition, Arc<Vec<i128>>);
table!(product_table, (Partition, Partition), Arc<Vec<(Partition, u128)>>);
table!(skew_table, (Partition, Partition), Arc<Vec<(Partition, u128)>>);
table!(htilde_table, u32, Arc<Vec<(Partition, Schur<QtPoly<i128>>)>>);
table!(bh_pieri_table, (Partition, Partition), Rat<i128>);
table!(bh_ell_table, (Partition, Partition), Rat<i128>);
table!(bold_p_table, Partition, Arc<PowerSum<GuardedRat>>);
table!(st_to_schur_table, Partition, Arc<Vec<(Partition, i128)>>);
table!(schur_to_st_table, Partition, Arc<Vec<(Partition, i128)>>);
table!(
    reduced_kronecker_table,
    (Partition, Partition),
    Arc<Vec<(Partition, i128)>>
);

/// `H̃_μ` in the Schur basis for a whole degree — the modified (q,t)-Kostka
/// coefficients, cached at `i128` and converted by the caller.
///
/// Cached because the recursion's unit of work is the **degree**, while
/// `qt_kostka` and `qt_kostka_column` are asked for one value or one column.
/// Without this, taking every column of degree 9 one at a time costs 30× the
/// whole table — the mistake [`kostka_table`](crate::kostka) documents, arrived
/// at from the other direction.
///
/// `i128` rather than a type parameter, following
/// [`character_cached`](character_cached): a `static` cannot be generic, and the
/// values are integers. That is safe here by a bound and not just a measurement.
/// `K̃_{λμ}` has non-negative coefficients (Haiman) summing to `K̃_{λμ}(1,1) =
/// f^λ`, and `Σ_λ (f^λ)² = n!`, so no coefficient exceeds `√(n!)` — past `i128`
/// only around degree 57, which the enumeration never reaches. Measured, they
/// are 9 bits at degree 12 and growing about 1 per degree.
pub fn htilde_cached(
    n: u32,
    compute: impl FnOnce() -> Vec<(Partition, Schur<QtPoly<i128>>)>,
) -> Arc<Vec<(Partition, Schur<QtPoly<i128>>)>> {
    lookup(htilde_table(), &n, || Arc::new(compute()))
}

/// The Bergeron–Haiman Pieri coefficient `c⁽ʳ⁾_{μν}`, and `L_{μν} = ⟨H̃_μ, h_ν⟩`.
///
/// Cached **across degrees**, which is the point: computing degree `n` needs
/// both at every size below `n`, so a degree-12 run rebuilds most of what a
/// degree-11 run already knew. Sharing them is 1.5× on a walk up the degrees
/// (2.15s → 1.45s for 1..=12) for **no extra memory at all** — peak RSS moves
/// 156MB → 157MB, because the degree-12 call was building that cache inside
/// itself either way. All the sharing does is stop the smaller degrees
/// rebuilding it.
///
/// It buys nothing for a single cold degree, which is what
/// `bench_qt_kostka.py` measures, so the headline comparison against Sage is
/// unaffected either way.
///
/// `i128` for the same reason as [`htilde_cached`] — a `static` cannot be
/// generic — and safe by measurement rather than a bound here, since these are
/// intermediate rational functions and not the coefficients Haiman's theorem
/// constrains: 25 bits at degree 12, growing about 3 per degree, so `i128`
/// holds past degree 45.
pub fn bh_pieri_cached(
    key: &(Partition, Partition),
    compute: impl FnOnce() -> Rat<i128>,
) -> Rat<i128> {
    lookup(bh_pieri_table(), key, compute)
}

/// `L_{μν}`; see [`bh_pieri_cached`].
pub fn bh_ell_cached(
    key: &(Partition, Partition),
    compute: impl FnOnce() -> Rat<i128>,
) -> Rat<i128> {
    lookup(bh_ell_table(), key, compute)
}

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
///
/// `compute` returning `None` means the character overflows `i128`. That is
/// **not** cached: it is the absence of a representable value rather than a
/// value, and caching it would make the table unable to distinguish "not yet
/// computed" from "cannot be represented". Overflow is confined to |λ| ≳ 58, so
/// recomputing it is not a hot path.
pub fn character_cached(
    lambda: &Partition,
    mu: &Partition,
    compute: impl FnOnce() -> Option<i128>,
) -> Option<i128> {
    let key = (lambda.clone(), mu.clone());
    if let Some(&v) = character_table().read().unwrap().get(&key) {
        return Some(v);
    }
    let v = compute()?;
    character_table().write().unwrap().insert(key, v);
    Some(v)
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

/// The partitions of `n` in **decreasing-lex** order, the order in which the
/// Kostka matrix is upper-unitriangular.
pub fn lex_parts_cached(n: u32) -> Arc<Vec<Partition>> {
    lookup(lex_parts_table(), &n, || {
        let mut v = partitions_of(n);
        v.sort_by(|a, b| b.parts().cmp(a.parts()));
        Arc::new(v)
    })
}

/// Memoized **one row** of the inverse Kostka matrix, indexed against
/// [`lex_parts_cached`].
///
/// Keyed by μ rather than by degree because m → s needs a single row per term.
/// Caching the whole matrix instead meant one conversion paid for p(n)² tableau
/// counts to read p(n) of them.
pub fn inverse_kostka_row_cached(
    mu: &Partition,
    compute: impl FnOnce() -> Vec<i128>,
) -> Arc<Vec<i128>> {
    lookup(inverse_kostka_row_table(), mu, || Arc::new(compute()))
}

/// `𝐩_γ` (OZ Eq 24) in the power-sum basis — the image of `p_γ` under the map Γ
/// that carries `s_λ` to `s̃_λ`.
///
/// Keyed by γ alone and shared across every λ and every degree, because that is
/// what it is: a function of γ. One `s̃_λ` of degree n needs all p(n) of them,
/// and every other λ of that degree needs the same ones again.
///
/// **Peek and store are separate on purpose**, unlike every other table here.
/// The values are fixed-width rationals that *report* overflow rather than
/// wrapping, and a value computed by a call that overflowed is garbage. Caching
/// it would be worse than recomputing it: the overflow counter is checked around
/// the call that produced it, so a later reader of the poisoned entry would see
/// a clean counter and accept a wrong answer. The caller stores only what it has
/// confirmed clean; see `character_basis::bold_guarded`.
pub fn bold_p_peek(gamma: &Partition) -> Option<Arc<PowerSum<GuardedRat>>> {
    bold_p_table().read().unwrap().get(gamma).cloned()
}

/// See [`bold_p_peek`]. Storing an entry asserts it was computed without
/// overflow.
pub fn bold_p_store(gamma: &Partition, value: PowerSum<GuardedRat>) {
    bold_p_table()
        .write()
        .unwrap()
        .insert(gamma.clone(), Arc::new(value));
}

/// Memoized `s̃_λ` in the Schur basis, and its inverse `s_λ` in the `s̃` basis.
///
/// Cached at `i128` with the generic conversion at the edges — the same shape as
/// [`htilde_cached`], for the same reason: a `static` cannot be generic, and
/// both transitions are integral (OZ Thm 1(2) makes the `s → s̃` direction a
/// matrix of *non-negative* integers, being multiplicities in a restriction).
pub fn st_to_schur_cached(
    lambda: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, i128)>,
) -> Arc<Vec<(Partition, i128)>> {
    lookup(st_to_schur_table(), lambda, || Arc::new(compute()))
}

/// See [`st_to_schur_cached`]; this is the other direction.
pub fn schur_to_st_cached(
    nu: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, i128)>,
) -> Arc<Vec<(Partition, i128)>> {
    lookup(schur_to_st_table(), nu, || Arc::new(compute()))
}

/// Memoized `s̃_λ · s̃_μ` — one whole column of reduced Kronecker coefficients.
///
/// Whole expansion rather than single coefficients, following
/// [`product_cached`]: the engine produces every ν at once, so caching per
/// coefficient would recompute the column once per ν asked for.
pub fn reduced_kronecker_cached(
    lambda: &Partition,
    mu: &Partition,
    compute: impl FnOnce() -> Vec<(Partition, i128)>,
) -> Arc<Vec<(Partition, i128)>> {
    lookup(
        reduced_kronecker_table(),
        &(lambda.clone(), mu.clone()),
        || Arc::new(compute()),
    )
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

/// The cached expansion of s_{outer/inner}, if some earlier call computed it.
///
/// A read-only peek — never computes. This is what lets a coefficient query
/// answer from an expansion a previous caller already paid for (e.g. the
/// whole-product expansion behind c^λ_{μν}) instead of starting a fresh
/// traversal of its own.
pub fn skew_cache_peek(
    outer: &Partition,
    inner: &Partition,
) -> Option<Arc<Vec<(Partition, u128)>>> {
    skew_table()
        .read()
        .unwrap()
        .get(&(outer.clone(), inner.clone()))
        .cloned()
}

/// Drop every cached table, releasing the memory.
pub fn clear_caches() {
    htilde_table().write().unwrap().clear();
    bh_pieri_table().write().unwrap().clear();
    bh_ell_table().write().unwrap().clear();
    partitions_table().write().unwrap().clear();
    character_table().write().unwrap().clear();
    kostka_table().write().unwrap().clear();
    lr_table().write().unwrap().clear();
    lex_parts_table().write().unwrap().clear();
    inverse_kostka_row_table().write().unwrap().clear();
    product_table().write().unwrap().clear();
    skew_table().write().unwrap().clear();
    bold_p_table().write().unwrap().clear();
    st_to_schur_table().write().unwrap().clear();
    schur_to_st_table().write().unwrap().clear();
    reduced_kronecker_table().write().unwrap().clear();
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
