//! The (q,t)-Kostka polynomials `K_{λμ}(q,t)`.
//!
//! ```text
//!   J_μ(x; q, t) = Σ_λ K_{λμ}(q,t) S_λ(x; t),      S_λ(x; t) = s_λ[X(1 − t)]
//! ```
//!
//! Macdonald, *Symmetric Functions and Hall Polynomials*, 2nd ed., VI (8.11).
//! `K_{λμ}(0,t)` is the Kostka–Foulkes polynomial `K_{λμ}(t)` and `K_{λμ}(0,1)`
//! the ordinary Kostka number, so this is the two-variable top of the tower
//! [`kostka`](mod@crate::kostka) and [`kf`](mod@crate::kf) already occupy — and
//! both of them are available here as independent checks, which is most of why
//! the test module below is worth its length.
//!
//! Unlike the Kostka and Kostka–Foulkes tables this one is **dense**: `K_{λμ}`
//! is generally nonzero even when λ does not dominate μ. `K_{(21),(3)} = q² +
//! q` is the smallest witness. Nor is the diagonal 1: `K_{(21),(21)} = 1 + qt`.
//! Both the triangularity and the unit diagonal are `q = 0` phenomena, not
//! properties of these — which is worth knowing before writing a test that
//! assumes the Kostka–Foulkes shape carries over.
//!
//! ## Three routes, and which one runs
//!
//! - [`qt_kostka_table`] — the **Bergeron–Haiman** Pieri recursion
//!   ([`bh`](crate::bh)). This is what every entry point here reaches, and the
//!   reason they are all bounded on [`Ring`] rather than
//!   [`QAlgebra`]: neither the recursion nor the `m → s`
//!   transition ever divides by an integer, so ℚ is not needed and the Python
//!   layer runs the whole thing over `i128`.
//! - [`qt_kostka_table_via_branching`] — `J_μ` by the branching formula, then
//!   the `S`-basis inversion below. Slower, and kept for the reason
//!   [`NaiveLr`](crate::NaiveLr) and
//!   [`kostka_foulkes_by_charge`](crate::charge::kostka_foulkes_by_charge) are
//!   kept: it shares no code with the recursion, so agreement is evidence
//!   rather than tautology. It is also the route held to Sage every pair
//!   through degree 7, which is what the fast one inherits. Sage's equivalent
//!   is `sage.combinat.sf.macdonald.qt_kostka`, and [`macdonald_ht`]'s is the
//!   `Ht` basis; `scripts/check_qt_kostka.py` compares against both.
//! - [`qt_kostka_table_via_operator`] — Lapointe–Lascoux–Morse, via
//!   [`macop`](crate::macop). Slowest of the three, and kept on the same
//!   argument one step further: three algorithms sharing nothing above
//!   `Partition` is stronger than two, and this one is built from an
//!   eigenvector problem rather than a tableau sum.
//!
//! `examples/bench_qtk_routes.rs` asserts all three agree at every degree it
//! times, which is where that evidence is actually collected.
//!
//! ## Range
//!
//! Over a fixed-width `C` this family refuses rather than wrapping past its
//! wall (`docs/policies/failure.md`, R3), and **the wall is not reachable**:
//! at `i128` the widest coefficient of [`qt_kostka_table`] gains ~1.8 bits per
//! degree and would reach 127 bits near n ≈ 77, against a table that stops
//! finishing around n = 18. Asking for one value or one column does not change
//! that — both cost the whole table, as above — so the degree is the unit of
//! work whatever is asked for.
//!
//! `K̃` has a bound as well as a measurement: Haiman's positivity makes the
//! coefficients non-negative with `K̃_{λμ}(1,1) = f^λ`, and `Σ_λ (f^λ)² = n!`
//! caps them at `√(n!)`, which passes `i128` only near n ≈ 57.
//!
//! Escalation is deliberately absent, and would not work if added at the
//! boundary: [`crate::bh::htilde_table`] computes its cache at `i128` whatever
//! `C` is, so a `BigInt` instantiation walls in the same place.
//! `docs/policies/failure.md` owns that decision. Measurements in
//! `docs/record/failure-and-overflow.md` (`examples/probe_qt_walls.rs`).
//!
//! ## Inverting the S basis
//!
//! Reading `K` off means expanding `J_μ` in `{S_λ}` rather than in `{s_λ}`, so
//! what is needed is the inverse of `f ↦ f[X(1−t)]`. Call it
//!
//! ```text
//!   φ_t : p_n ↦ p_n / (1 − t^n),   extended ℚ(q,t)-linearly.
//! ```
//!
//! `S_λ = s_λ[X(1−t)]` has `p_ν`-coefficient `z_ν⁻¹ χ^λ_ν ∏_i (1 − t^{ν_i})`,
//! so `φ_t(S_λ) = s_λ` and therefore `φ_t(J_μ) = Σ_λ K_{λμ}(q,t) s_λ`. The
//! whole computation is: expand `J_μ` in power sums, divide the `p_ν`
//! coefficient by `∏_i (1 − t^{ν_i})`, expand back into Schur.
//!
//! ## It is not a plethysm, and that matters
//!
//! The step looks like it needs a [`Plethystic`](crate::coeff::Plethystic)
//! impl for [`Frac`]. It does not, and writing one and using it here would
//! give the wrong answer.
//!
//! `φ_t` is **ℚ(q,t)-linear**: it acts on the alphabet `X` and holds the
//! coefficients fixed. A plethysm does the opposite —
//! [`Plethystic::frobenius`](crate::coeff::Plethystic::frobenius) exists
//! precisely because `p_n` substitutes into the coefficient ring too, so
//! `p_n[t·p_1] = t^n·p_n`. Routing `J_μ` through the genuine plethysm would
//! raise the `q` and `t` already sitting in `J`'s coefficients, and `K` would
//! come out in the wrong variables.
//!
//! The two are easy to conflate because both are written `f[X/(1−t)]`. The
//! `(1 − t^n)` appearing in the divisor *is* the Frobenius image of `1 − t` —
//! but it comes from `S_λ`'s definition, where the coefficients are rational
//! constants and there is nothing else to raise. Nothing in `J` gets raised.
//! [`plethysm`](mod@crate::plethysm) is the operation that would have been
//! wrong here; dividing the `p_ν` coefficient by `∏_i (1 − t^{ν_i})`, as
//! above, is what is right.
//!
//! ## Why the coefficient ring is a ℚ-algebra
//!
//! `s → p` divides by `z_ν` and `φ_t` divides by `1 − t^n`, so the route runs
//! through ℚ(q,t) even though `K_{λμ}(q,t) ∈ ℤ[q,t]` at the end (Haiman's
//! theorem gives non-negative integers, which the tests check but do not rely
//! on). The denominators must cancel exactly, and
//! [`Frac::into_poly`] is where that is enforced: a leftover factor is a bug,
//! not a fallback.
//!
//! ## Going the other way
//!
//! [`schur_in_j_table`] is the `s → J` transition — the inverse of `J → s`,
//! computed as a projection rather than a solve — and
//! [`schur_to_macdonald_j`] applies it to one element, which is the form a
//! caller asking whether a function is `J`-positive wants.
//!
//! ## References
//!
//! - **\[Hai\]** M. Haiman, *Hilbert schemes, polygraphs and the Macdonald
//!   positivity conjecture*, J. Amer. Math. Soc. **14** (2001) — Theorem 2,
//!   which proves Conjecture 2.1.2, `K̃_{λμ}(q,t) ∈ ℕ[q,t]`: the positivity
//!   the bound above rests on.
//! - **\[M\]** I. G. Macdonald, *Symmetric Functions and Hall Polynomials*,
//!   2nd ed., Oxford, 1995 — VI (8.11), the definition above, and VI (8.16),
//!   `K_{λμ}(1,1) = f^λ`, which `K̃` shares and which is the other half of the
//!   bound.

// Degree and shape indices; coefficients are `QtPoly<C>`.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::BTreeMap;

use crate::coeff::{QAlgebra, Ring};
use crate::convert::{FromSchur, ToSchur};
use crate::frac::Frac;
use crate::partition::Partition;
use crate::qt::QtPoly;
use crate::sym::{add_at, by_degree, Monomial, PowerSum, Schur, SymFn};

/// `K_{λμ}(q,t)`.
///
/// Zero unless `|λ| = |μ|`.
///
/// This computes the whole table of degree `|μ|` and reads one coefficient out
/// of it, exactly as [`kostka_foulkes`](crate::kf::kostka_foulkes) does with
/// `Q'_μ`. For more than one λ at a fixed μ use [`qt_kostka_column`], and for a
/// whole degree [`qt_kostka_table`].
pub fn qt_kostka<C: Ring>(lambda: &Partition, mu: &Partition) -> QtPoly<C> {
    if lambda.size() != mu.size() {
        return QtPoly::zero();
    }
    qt_kostka_column::<C>(mu)
        .into_iter()
        .find(|(l, _)| l == lambda)
        .map(|(_, k)| k)
        .unwrap_or_else(QtPoly::zero)
}

/// Every `K_{λμ}(q,t)` for a fixed μ, paired with its λ.
///
/// Computed by taking a column of [`qt_kostka_table`], which is **not** the
/// waste it looks like. The Bergeron–Haiman recursion's unit of work is the
/// degree. Its whole table costs about what one column of the branching
/// formula costs (`docs/record/qt-kostka.md`).
///
/// Every λ of the degree appears — see the note on density in the module docs.
///
/// # Panics
///
/// Panics only on a bug in this crate, never on an input: μ is by construction
/// one of the partitions of `|μ|`, and the wall the table itself can hit is the
/// t-degree check in [`qt_kostka_table_via_bh`].
pub fn qt_kostka_column<C: Ring>(mu: &Partition) -> Vec<(Partition, QtPoly<C>)> {
    let parts = crate::memo::partitions_cached(mu.size());
    let j = parts
        .iter()
        .position(|p| p == mu)
        .expect("mu must be a partition of its own size");
    let table = qt_kostka_table::<C>(mu.size());
    parts
        .iter()
        .enumerate()
        .map(|(i, lambda)| (lambda.clone(), table[i][j].clone()))
        .collect()
}

/// The whole table of degree `n`, as `table[i][j] = K_{λⁱ λʲ}(q,t)` indexed
/// against `memo::partitions_cached` — the same
/// orientation as [`kostka_table`](crate::kostka::kostka_table) and
/// [`kostka_foulkes_table`](crate::kf::kostka_foulkes_table).
///
/// Computed by the Bergeron–Haiman recursion (see [`bh`](crate::bh)), which is
/// the fast route and shares its work across the whole degree — a value like
/// `c⁽¹⁾_{(3,1),(3)}` is reached from every μ and ν whose recursion passes
/// through it. Against the branching formula that produces
/// [`qt_kostka_column`] its margin grows with the degree
/// (`docs/record/qt-kostka.md`, "The three routes on one table").
///
/// [`qt_kostka_table_via_operator`] and [`qt_kostka_table_via_branching`] are
/// the same table by two other algorithms, kept because agreement between three
/// routes that share nothing above `Partition` is the evidence this rests on.
pub fn qt_kostka_table<C: Ring>(n: u32) -> Vec<Vec<QtPoly<C>>> {
    qt_kostka_table_via_bh(n)
}

/// The Schur functions of degree `n` in the Macdonald `J` basis:
/// `table[i][j]` is the coefficient of `J_{λʲ}` in `s_{λⁱ}`, indexed against
/// `memo::partitions_cached` as
/// [`qt_kostka_table`] is.
///
/// This is the **inverse** of the `J → s` transition, and it is not computed as
/// one. `{J_μ}` is orthogonal for the `(q,t)` scalar product, and
/// `⟨P_μ, Q_μ⟩ = 1` (Macdonald VI (6.19)) with `J_μ = c_μ P_μ = c'_μ Q_μ` makes
/// `⟨J_μ, J_μ⟩ = c_μ c'_μ`, so the expansion is a projection:
///
/// ```text
///   s_λ = Σ_μ ⟨s_λ, J_μ⟩_{q,t} / (c_μ c'_μ) · J_μ
/// ```
///
/// and `⟨f, g⟩_{q,t} = ⟨f, g[X(1−q)/(1−t)]⟩` against the Hall product, which
/// makes the numerator a coefficient read off a Schur expansion:
/// `⟨s_λ, J_μ⟩_{q,t} = [s_λ] H_μ[X(1−q)]`, since `H_μ = J_μ[X/(1−t)]` is the
/// `Σ_λ K_{λμ} s_λ` this module already computes. So a whole degree of the
/// inverse costs one [`qt_kostka_table`] plus `p(n)` power-sum round trips —
/// against the `O(p(n)³)` triangular solve over ℚ(q,t) that inverting the
/// matrix takes.
///
/// The numerators are polynomials: `H_μ[X(1−q)]` divides by nothing, and only
/// the hook products put anything under the line. They are the denominators
/// [`Frac`] keeps factored.
///
/// The table is triangular, so about half the entries are zero, and `n = 0`
/// gives the one-by-one table holding 1. Sage's equivalent is
/// `MacdonaldPolynomials_j._s_to_self_cache[n]`, which it fills by inverting
/// the other direction.
///
/// Unlike the four `K` entry points here this one cannot run over `i128`: the
/// power-sum round trip divides by `z_ν`, so `C` is a [`QAlgebra`] and the
/// range is that ring's. The Python boundary escalates, so there is no wall
/// there.
///
/// # Panics
///
/// Panics only on a bug in this crate: the table's λ come from a Schur
/// expansion of degree `n`, and every partition of `n` is indexed.
///
/// # Examples
///
/// The diagonal is `1/c_λ`, and at λ = (2) that is `(1−t)(1−qt)` — **not**
/// `c'_λ = (1−q)(1−q²)`, the other hook product, which would give `1/3` at the
/// point below instead of `1/10`:
///
/// ```
/// use symfn::{schur_in_j_table, Rational, Ring};
///
/// // Index 0 is (2) and index 1 is (11), as `partitions_of(2)` orders them.
/// let w = schur_in_j_table::<Rational>(2);
/// let (q, t) = (Rational::from_int(2), Rational::from_int(3));
///
/// assert_eq!(w[0][0].eval(&q, &t), Some(Rational::new(1, 10)));
/// // s_(2) reaches J_(11); s_(11) does not reach J_(2).
/// assert_eq!(w[0][1].eval(&q, &t), Some(Rational::new(-1, 80)));
/// assert!(w[1][0].is_zero());
/// ```
pub fn schur_in_j_table<C: QAlgebra + Send + Sync + 'static>(n: u32) -> Vec<Vec<Frac<C>>> {
    (*cached_table::<C>(n)).clone()
}

/// The same table, shared rather than copied out — what the element-wise
/// expansion reads, since it only borrows the rows it needs.
fn cached_table<C: QAlgebra + Send + Sync + 'static>(n: u32) -> std::sync::Arc<Vec<Vec<Frac<C>>>> {
    crate::memo::schur_in_j_cached(n, || schur_in_j_table_uncached::<C>(n))
}

fn schur_in_j_table_uncached<C: QAlgebra>(n: u32) -> Vec<Vec<Frac<C>>> {
    let parts = crate::memo::partitions_cached(n);
    let index: std::collections::HashMap<&Partition, usize> =
        parts.iter().enumerate().map(|(i, p)| (p, i)).collect();
    let k = qt_kostka_table::<C>(n);
    let mut table = vec![vec![<Frac<C> as Ring>::zero(); parts.len()]; parts.len()];

    for (j, mu) in parts.iter().enumerate() {
        crate::interrupt::poll();
        let mut h: Schur<QtPoly<C>> = Schur::zero();
        for (i, lambda) in parts.iter().enumerate() {
            if !k[i][j].is_zero() {
                h.add_term(lambda.clone(), k[i][j].clone());
            }
        }
        // `X ↦ X(1−q)` is ℚ(q,t)-linear on the alphabet, so it is the same
        // scaling of `p_ν` that `invert_s_basis` performs the other way round —
        // and for the same reason it is not the plethysm of that name.
        let p: PowerSum<QtPoly<C>> = PowerSum::from_schur(&h);
        let mut scaled: PowerSum<QtPoly<C>> = PowerSum::zero();
        for (nu, c) in p.terms() {
            let mut v = c.clone();
            for &part in nu.parts() {
                v = v.mul_binomial(part, 0);
            }
            scaled.add_term(nu.clone(), v);
        }

        let mut den: BTreeMap<(u32, u32), i32> = BTreeMap::new();
        for (e, m) in crate::macdonald::c_factors(mu.parts()) {
            *den.entry(e).or_insert(0) -= m;
        }
        for (e, m) in crate::macdonald::c_prime_factors(mu.parts()) {
            *den.entry(e).or_insert(0) -= m as i32;
        }

        let g: Schur<QtPoly<C>> = scaled.to_schur();
        for (lambda, c) in g.terms() {
            let mut w = Frac::from_poly(c.clone()).mul_factors(&den);
            w.reduce();
            let i = index[lambda];
            table[i][j] = w;
        }
    }
    table
}

/// `f`, given in the Schur basis, rewritten in the Macdonald `J` basis: the
/// `c_μ` of `f = Σ_μ c_μ J_μ(x; q, t)`.
///
/// The element-wise form of [`schur_in_j_table`], which holds a whole degree
/// at once. Mixed degrees are accepted and handled degree by degree; the zero
/// element gives the empty map. The result is in element order with no zeros,
/// and costs one [`schur_in_j_table`] per degree present in `f` — which is
/// what expanding a single `s_λ` costs anyway, so a caller with several shapes
/// of one degree should prefer this to a call per shape. Sage's equivalent is
/// `Sym.macdonald().J()(f)`.
///
/// The coefficients are [`Frac`]s: `J` is the integral form, so the `J → s`
/// direction has polynomial coefficients, but the inverse divides by the hook
/// products `c_μ c'_μ` and those do not cancel.
///
/// ```
/// use symfn::{schur_to_macdonald_j, Partition, QtPoly, Rational, Ring, Schur, SymFn};
///
/// let one = QtPoly::term(0, 0, Rational::from_int(1));
/// let s11: Schur<QtPoly<Rational>> = Schur::monomial(Partition::new([1, 1]), one);
/// let in_j = schur_to_macdonald_j(&s11);
///
/// assert_eq!(in_j.len(), 1);
/// let (num, den) = in_j[&Partition::new([1, 1])].parts();
/// assert_eq!(*num, <QtPoly<Rational> as Ring>::one());
/// assert_eq!(den.collect::<Vec<_>>(), vec![(&(0, 1), &1), (&(0, 2), &1)]);
/// ```
///
/// So `s_11 = J_11 / ((1−t)(1−t²))` and nothing else: the denominator is
/// `c_{(11)}`, **not** `c'_{(11)} = (1−q)(1−q²)`, which is the other hook
/// product and the twist to check. The table is triangular the other way from
/// the `J → s` one — `s_2` reaches both `J_2` and `J_11`, while `s_11` reaches
/// `J_11` alone.
pub fn schur_to_macdonald_j<C: QAlgebra + Send + Sync + 'static>(
    f: &Schur<QtPoly<C>>,
) -> BTreeMap<Partition, Frac<C>> {
    let mut out: BTreeMap<Partition, Frac<C>> = BTreeMap::new();
    for (n, terms) in by_degree(f) {
        let parts = crate::memo::partitions_cached(n);
        let index: std::collections::HashMap<&Partition, usize> =
            parts.iter().enumerate().map(|(i, p)| (p, i)).collect();
        let w = cached_table::<C>(n);
        for (lambda, c) in terms {
            crate::interrupt::poll();
            let scale = Frac::from_poly(c.clone());
            let row = &w[index[lambda]];
            for (mu, entry) in parts.iter().zip(row) {
                if !entry.is_zero() {
                    add_at(&mut out, mu, entry.mul(&scale));
                }
            }
        }
    }
    // Reduced once, at the end: `Frac::add_assign` deliberately leaves a
    // running sum over the lcm of the denominators, because a cancellation can
    // only be decided when the sum is complete (`src/frac.rs`).
    for v in out.values_mut() {
        v.reduce();
    }
    out
}

/// The table by the branching formula — the reference implementation.
///
/// Slower than [`qt_kostka_table`] and kept as the thing that route is checked
/// against: this one goes through `macdonald_j`, which is verified against Sage
/// independently.
///
/// # Panics
///
/// Panics if the column expansion produces a λ that is not a partition of `n`,
/// which is a bug in the branching formula rather than an input this rejects.
pub fn qt_kostka_table_via_branching<C: QAlgebra>(n: u32) -> Vec<Vec<QtPoly<C>>> {
    let parts = crate::memo::partitions_cached(n);
    let mut table = vec![vec![QtPoly::zero(); parts.len()]; parts.len()];
    for (j, mu) in parts.iter().enumerate() {
        crate::interrupt::poll();
        for (lambda, k) in column_expansion::<C>(mu).terms() {
            let i = parts
                .iter()
                .position(|p| p == lambda)
                .expect("K_{λμ} is supported on partitions of |μ|");
            table[i][j] = k.clone();
        }
    }
    table
}

/// `H̃_μ(x; q, t) = Σ_λ K̃_{λμ}(q,t) s_λ`, the modified Macdonald polynomial in
/// the Schur basis.
///
/// Read straight out of the Bergeron–Haiman recursion, which produces `K̃`
/// **natively** — `H̃` is what that recursion is about, and `K` is the
/// reflected one. So this is the cheaper of the two and [`qt_kostka_table`] is
/// the one paying for a reflection.
///
/// `K̃_{λμ}(q,t) = t^{n(μ)} K_{λμ}(q, 1/t)`, the form the modern literature
/// uses and the one in which Haiman's positivity reads "non-negative integers"
/// with no normalizing power in the way. `H̃_{(2)} = s_2 + q·s_{11}` and
/// `H̃_{(11)} = s_2 + t·s_{11}` are the smallest pair, and show the `q ↔ t`
/// symmetry under conjugating μ that the twisted form has and `K` does not.
///
/// # Panics
///
/// Panics only on a bug in this crate: `htilde_table(|μ|)` carries a row for
/// every partition of `|μ|`, and μ is one of them.
pub fn macdonald_ht<C: Ring>(mu: &Partition) -> Schur<QtPoly<C>> {
    crate::bh::htilde_table::<C>(mu.size())
        .into_iter()
        .find(|(m, _)| m == mu)
        .map(|(_, s)| s)
        .expect("mu must be a partition of its own size")
}

/// `K̃_{λμ}(q,t)`, the modified (q,t)-Kostka polynomial — one coefficient of
/// [`macdonald_ht`], which is what it computes.
///
/// Zero unless `|λ| = |μ|`.
pub fn modified_qt_kostka<C: Ring>(lambda: &Partition, mu: &Partition) -> QtPoly<C> {
    if lambda.size() != mu.size() {
        return QtPoly::zero();
    }
    macdonald_ht::<C>(mu).coeff(lambda)
}

/// The whole table through the Bergeron–Haiman recursion — the fast route.
///
/// `H̃` comes back in the Schur basis with coefficients `K̃_{λμ}`, and `K` is
/// the `t`-reversal of that: `K_{λμ}(q,t) = t^{n(μ)} K̃_{λμ}(q, 1/t)`, the same
/// involution [`macdonald_ht`] applies in the other direction. Bookkeeping, not
/// arithmetic.
///
/// Bounded on [`Ring`] and not [`QAlgebra`], unlike every other route here:
/// neither the recursion nor `m → s` ever divides by an integer, so this runs
/// over `QtPoly<i128>` where the others need ℚ.
///
/// # Panics
///
/// Panics if some `K̃_{λμ}` carries a `t`-degree above `n(μ)`. Macdonald theory
/// says it cannot, and the reflection to `K` subtracts that degree — so an
/// unchecked violation would silently produce negative exponents rather than
/// fail. It is asserted because it is the step the reflection rests on.
pub fn qt_kostka_table_via_bh<C: Ring>(n: u32) -> Vec<Vec<QtPoly<C>>> {
    let parts = crate::memo::partitions_cached(n);
    let index: std::collections::HashMap<&Partition, usize> =
        parts.iter().enumerate().map(|(i, p)| (p, i)).collect();
    let mut table = vec![vec![QtPoly::zero(); parts.len()]; parts.len()];
    for (j, (mu, s)) in crate::bh::htilde_table::<C>(n).into_iter().enumerate() {
        crate::interrupt::poll();
        let n_mu: u32 = mu
            .parts()
            .iter()
            .enumerate()
            .map(|(i, &p)| i as u32 * p)
            .sum();
        debug_assert_eq!(&mu, &parts[j], "htilde_table must share the order");
        for (lambda, kt) in s.terms() {
            let mut k = QtPoly::zero();
            for (&(a, b), c) in kt.terms() {
                assert!(
                    b <= n_mu,
                    "K~_{{{lambda},{mu}}} has t-degree {b} > n(mu) = {n_mu}"
                );
                k.add_term(a, n_mu - b, c.clone());
            }
            table[index[lambda]][j] = k;
        }
    }
    table
}

/// The same column, reached through the Macdonald operator instead of the
/// branching formula.
///
/// ## The two plethysms compose into one
///
/// [`macop::eigenvector`](crate::macop) returns `J_μ` in `S_κ[X^{tq}]` with
/// `X^{tq} = X(t−1)/(q−1)`, so getting to `K` looks like two substitutions —
/// cross to an ordinary alphabet, then apply `φ_t`. They collapse:
///
/// ```text
///   p_k ↦ p_k (t^k−1)/(q^k−1)   then   p_k ↦ p_k/(1−t^k)
///     = p_k (t^k−1) / ((q^k−1)(1−t^k))
///     = p_k / (1 − q^k)
/// ```
///
/// One pass through the power sums, and `1 − q^k` is a binomial [`Frac`]
/// already holds. The `t` half cancels outright — which is worth noticing,
/// because doing the two substitutions separately would build and then destroy
/// every `(1 − t^k)` in the expansion.
///
/// ## Normalization and the two divisions
///
/// The eigenvector fixes the `S_μ` coefficient at 1 where `J_μ` has
/// `c_{μ'}(t,q)` (\[LLM\] 3.15), which is
/// [`c_prime_factors`](crate::macdonald) — applied one binomial at a time, so
/// nothing expands it. And the solve returns `b_κ = a_κ · v`, so the last step
/// divides `v` back out with [`QtPoly::divide_exact`](crate::qt::QtPoly). Both
/// are exact by construction and neither is a gcd.
// Per-shape form of the operator route, used by the tests that hold the three
// routes to each other. The module doc names it, and the record keeps this
// route "twice over, since it shares no mathematics with either alternative"
// (`docs/record/qt-kostka.md`).
#[allow(dead_code)]
fn column_via_operator<C: QAlgebra>(mu: &Partition) -> Schur<QtPoly<C>> {
    let (b, v) = crate::macop::eigenvector::<C>(mu);
    kostka_from_eigenvector(mu, &b, &v)
}

/// The whole table through the operator route, sharing `M₁` across the degree.
///
/// The matrix depends only on the degree, so this is the unit of work that
/// route wants.
///
/// # Panics
///
/// Panics if a normalized `K_{λμ}` is not a polynomial, or if the solve's
/// common denominator does not divide it. Either is a bug in this crate rather
/// than an input it rejects.
pub fn qt_kostka_table_via_operator<C: QAlgebra>(n: u32) -> Vec<Vec<QtPoly<C>>> {
    let parts = crate::memo::partitions_cached(n);
    let index: std::collections::HashMap<&Partition, usize> =
        parts.iter().enumerate().map(|(i, p)| (p, i)).collect();
    let mut table = vec![vec![QtPoly::zero(); parts.len()]; parts.len()];
    for (j, (mu, b, v)) in crate::macop::eigenvectors::<C>(n).into_iter().enumerate() {
        crate::interrupt::poll();
        for (lambda, k) in kostka_from_eigenvector(&mu, &b, &v).terms() {
            table[index[lambda]][j] = k.clone();
        }
    }
    table
}

fn kostka_from_eigenvector<C: QAlgebra>(
    mu: &Partition,
    b: &[QtPoly<C>],
    v: &QtPoly<C>,
) -> Schur<QtPoly<C>> {
    let parts = crate::memo::partitions_cached(mu.size());

    let mut formal: Schur<Frac<C>> = Schur::zero();
    for (k, kappa) in parts.iter().enumerate() {
        crate::interrupt::poll();
        if !b[k].is_empty() {
            formal.add_term(kappa.clone(), Frac::from_poly(b[k].clone()));
        }
    }

    // Ψ: p_k ↦ p_k / (1 − q^k), the composite of the two substitutions.
    let p: PowerSum<Frac<C>> = PowerSum::from_schur(&formal);
    let mut scaled = PowerSum::zero();
    for (rho, c) in p.terms() {
        let mut factors: BTreeMap<(u32, u32), i32> = BTreeMap::new();
        for &k in rho.parts() {
            *factors.entry((k, 0)).or_insert(0) -= 1;
        }
        let mut w = c.mul_factors(&factors);
        w.reduce();
        scaled.add_term(rho.clone(), w);
    }
    let expanded: Schur<Frac<C>> = scaled.to_schur();

    // `c_{μ'}(t,q)` is applied **before** leaving `Frac`, not after. The
    // denominators `Ψ` leaves behind are canceled by it and not by anything
    // else: at μ = (1) the expansion is `s_1/(1−q)` and `c' = 1−q`, so asking
    // for a polynomial first fails on the smallest case there is.
    let cprime: BTreeMap<(u32, u32), i32> = crate::macdonald::c_prime_factors(mu.parts())
        .into_iter()
        .map(|(k, m)| (k, m as i32))
        .collect();
    let mut out = Schur::zero();
    for (lambda, c) in expanded.terms() {
        let scaled = c.mul_factors(&cprime);
        let num = scaled.clone().into_poly().unwrap_or_else(|| {
            panic!("K_{{{lambda},{mu}}} is not a polynomial after normalizing: {scaled}")
        });
        let k = num
            .divide_exact(v)
            .unwrap_or_else(|| panic!("v must divide the normalized K_{{{lambda},{mu}}}"));
        out.add_term(lambda.clone(), k);
    }
    out
}

/// `φ_t(J_μ)` in the Schur basis — the column, before it is taken apart.
fn column_expansion<C: QAlgebra>(mu: &Partition) -> Schur<QtPoly<C>> {
    let j: Monomial<Frac<C>> = crate::macdonald_j(mu);
    let p: PowerSum<Frac<C>> = PowerSum::from_schur(&j.to_schur());
    let expanded: Schur<Frac<C>> = invert_s_basis(&p).to_schur();

    let mut out = Schur::zero();
    for (lambda, c) in expanded.terms() {
        let k = c.clone().into_poly().unwrap_or_else(|| {
            panic!("K_{{{lambda},{mu}}} is not a polynomial: {c}");
        });
        out.add_term(lambda.clone(), k);
    }
    out
}

/// `φ_t`: divide the coefficient of `p_ν` by `∏_i (1 − t^{ν_i})`.
///
/// Applied through [`Frac::mul_factors`], so the divisor is never expanded —
/// `(0, ν_i) ↦ −1` is a denominator factor in the encoding [`Frac`] already
/// speaks, and factors shared with the numerator cancel without a product being
/// formed. A partition has no zero parts, so `(0,0)` cannot arise.
///
/// **Reduced here**, which is where this route's whole cost turned out to live.
/// `PowerSum::to_schur` uses each of these p(n) times — once per λ — and every
/// use is an `add_assign` that lifts the running sum to the lcm of the two
/// denominators. Cutting a denominator down once, here, is p(n) lifts it does
/// not have to widen, and it more than halves the whole route
/// (`docs/record/qt-kostka.md`).
///
/// It is the same shape as the reduce in [`Ring::mul`](Frac::mul) — reduce
/// once, before a value is used many times — and it makes that one redundant
/// for this path, since a coefficient now arrives already reduced.
fn invert_s_basis<C: QAlgebra>(p: &PowerSum<Frac<C>>) -> PowerSum<Frac<C>> {
    let mut out = PowerSum::zero();
    for (nu, c) in p.terms() {
        let mut factors: BTreeMap<(u32, u32), i32> = BTreeMap::new();
        for &part in nu.parts() {
            *factors.entry((0, part)).or_insert(0) -= 1;
        }
        let mut v = c.mul_factors(&factors);
        v.reduce();
        out.add_term(nu.clone(), v);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::{Rational, Ring};

    type Q = QtPoly<Rational>;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn k(lambda: &[u32], mu: &[u32]) -> Q {
        qt_kostka(&part(lambda), &part(mu))
    }

    /// Hand-computed in the module's own terms, and the case that pins the
    /// orientation: `K_{(11),(2)} = q` while `K_{(2),(11)} = t`. Swapping the
    /// two indices survives every symmetric check further down.
    ///
    /// `J_(11) = (1−t)(1−t²) e_2`, and `φ_t(e_2) = (t·s_2 +
    /// s_11)/((1−t)(1−t²))`, so `K_{(2),(11)} = t` and `K_{(11),(11)} = 1`.
    #[test]
    fn degree_two_is_the_hand_computation() {
        assert_eq!(k(&[2], &[2]), <Q as Ring>::one());
        assert_eq!(k(&[1, 1], &[2]), QtPoly::term(1, 0, Rational::from_int(1)));
        assert_eq!(k(&[2], &[1, 1]), QtPoly::term(0, 1, Rational::from_int(1)));
        assert_eq!(k(&[1, 1], &[1, 1]), <Q as Ring>::one());
    }

    /// `K_{λμ}(0,t)` is the Kostka–Foulkes polynomial.
    ///
    /// The strongest check available in-crate, and cheap to state: this module
    /// runs `J_μ` through the power-sum basis over ℚ(q,t), while
    /// [`kostka_foulkes`](crate::kf::kostka_foulkes) reads its answer off the
    /// Morris recursion in `ℤ[t]`. No shared code below `Partition`.
    #[test]
    fn at_q_zero_it_is_kostka_foulkes() {
        for n in 0..=6u32 {
            for mu in crate::partitions_of(n) {
                for (lambda, kqt) in qt_kostka_column::<Rational>(&mu) {
                    let want: QtPoly<Rational> = crate::kostka_foulkes(&lambda, &mu);
                    let got = set_q_to_zero(&kqt);
                    assert_eq!(got, want, "K_{{{lambda},{mu}}} at q = 0");
                }
            }
        }
    }

    /// `K_{λμ}(q,t) = K_{λ'μ'}(t,q)` — Macdonald VI (8.15).
    ///
    /// Independent of the `q = 0` check in a way that matters: that one only
    /// ever looks at the `q⁰` slice, so it says nothing about how the two
    /// variables are stored relative to each other. This one is entirely about
    /// that.
    #[test]
    fn conjugating_both_shapes_swaps_q_and_t() {
        for n in 0..=6u32 {
            for mu in crate::partitions_of(n) {
                let conj = qt_kostka_column::<Rational>(&mu.conjugate());
                for (lambda, kqt) in qt_kostka_column::<Rational>(&mu) {
                    let other = conj
                        .iter()
                        .find(|(l, _)| l == &lambda.conjugate())
                        .map(|(_, k)| swap_variables(k))
                        .unwrap_or_else(QtPoly::zero);
                    assert_eq!(kqt, other, "K_{{{lambda},{mu}}} under conjugation");
                }
            }
        }
    }

    /// Haiman's theorem: the coefficients are non-negative integers.
    ///
    /// Not something this code arranges — it arrives over ℚ(q,t) and every
    /// denominator has to cancel — so it is a real check on the whole route and
    /// not a restatement of it. `into_poly` already refuses a surviving
    /// denominator; this catches a surviving `1/2`.
    #[test]
    fn coefficients_are_non_negative_integers() {
        for n in 0..=6u32 {
            for mu in crate::partitions_of(n) {
                for (lambda, kqt) in qt_kostka_column::<Rational>(&mu) {
                    for (_, c) in kqt.terms() {
                        assert_eq!(c.denom(), 1, "K_{{{lambda},{mu}}} has {c:?} in it");
                        assert!(c.numer() > 0, "K_{{{lambda},{mu}}} has {c:?} in it");
                    }
                }
            }
        }
    }

    /// `K_{λμ}(1,1) = f^λ`, the number of standard tableaux of shape λ — the
    /// **same value for every μ**.
    ///
    /// The Garsia–Haiman module `R_μ` is a graded version of the regular
    /// representation of `S_n`, and these are its graded multiplicities, so
    /// specializing the grading away leaves `dim S^λ` no matter which μ was
    /// chosen. That independence is what makes this a sharp test: `q = 0` only
    /// ever inspects one slice of each polynomial, and this reads all of it.
    ///
    /// It is also where the density claim in the module docs is checked, since
    /// `f^λ > 0` always. `K_{(21),(3)} = q² + q` is nonzero even though (21)
    /// does not dominate (3) — an implementation that imposed the dominance
    /// triangularity of the *ordinary* Kostka numbers would pass the `q = 0`
    /// test and fail here.
    #[test]
    fn at_q_and_t_one_it_counts_standard_tableaux() {
        for n in 0..=6u32 {
            for mu in crate::partitions_of(n) {
                let column = qt_kostka_column::<Rational>(&mu);
                assert_eq!(
                    column.len(),
                    crate::partitions_of(n).len(),
                    "the column of {mu} should be dense"
                );
                for (lambda, kqt) in column {
                    let got: i128 = kqt.terms().map(|(_, c)| c.numer()).sum();
                    let want = crate::dimension(&lambda).unwrap() as i128;
                    assert_eq!(got, want, "K_{{{lambda},{mu}}} at q = t = 1");
                }
            }
        }
    }

    /// The table must be the columns, in the documented orientation.
    #[test]
    fn the_table_agrees_with_the_columns() {
        for n in 0..=5u32 {
            let parts = crate::partitions_of(n);
            let table = qt_kostka_table_via_branching::<Rational>(n);
            for (j, mu) in parts.iter().enumerate() {
                for (i, lambda) in parts.iter().enumerate() {
                    assert_eq!(table[i][j], qt_kostka(lambda, mu), "[{i}][{j}]");
                }
            }
        }
        // ...and the orientation check only distinguishes anything on an
        // asymmetric table.
        let t = qt_kostka_table::<Rational>(3);
        assert!(
            (0..t.len()).any(|i| (0..t.len()).any(|j| t[i][j] != t[j][i])),
            "a symmetric table would make the orientation untestable"
        );
    }

    /// `H̃_μ` is symmetric under conjugating μ together with swapping q and t —
    /// `K̃_{λμ}(q,t) = K̃_{λ'μ'}(t,q)`... with λ **not** conjugated, unlike the
    /// relation `K` satisfies.
    ///
    /// That difference is why this form exists, and it is what a
    /// wrong `n(μ)` would break: the reflection is by a shape-dependent power,
    /// so getting `n(μ) = Σ(i−1)μ_i` confused with `n(μ') = Σ binom(μ_i, 2)`
    /// still yields polynomials and still passes an integrality check.
    #[test]
    fn the_modified_form_is_symmetric_in_q_and_t() {
        for n in 0..=6u32 {
            for mu in crate::partitions_of(n) {
                let conj = macdonald_ht::<Rational>(&mu.conjugate());
                for (lambda, kt) in macdonald_ht::<Rational>(&mu).terms() {
                    assert_eq!(
                        kt,
                        &swap_variables(&conj.coeff(lambda)),
                        "K~_{{{lambda},{mu}}}"
                    );
                }
            }
        }
    }

    /// The two smallest, which pin the direction of the reflection: an inverted
    /// `n(μ)` sends `H̃_{(11)}` to `s_2 + t·s_{11}` as well, and the pair would
    /// still look symmetric.
    #[test]
    fn h_tilde_of_degree_two_is_the_known_pair() {
        let one = Rational::from_int(1);
        let ht = macdonald_ht::<Rational>(&part(&[2]));
        assert_eq!(ht.coeff(&part(&[2])), <Q as Ring>::one());
        assert_eq!(ht.coeff(&part(&[1, 1])), QtPoly::term(1, 0, one));

        let ht = macdonald_ht::<Rational>(&part(&[1, 1]));
        assert_eq!(ht.coeff(&part(&[2])), <Q as Ring>::one());
        assert_eq!(ht.coeff(&part(&[1, 1])), QtPoly::term(0, 1, one));
    }

    /// The two routes to the same column must agree, term for term.
    ///
    /// One enumerates semistandard tableaux and multiplies out ψ; the other
    /// solves an eigenvector problem for the Macdonald operator, crosses two
    /// alphabets that collapse into one, and divides out a normalization. They
    /// share `Partition`, `QtPoly` and nothing above that.
    ///
    /// The comparison is against [`qt_kostka_column`], which is itself checked
    /// against Sage every pair through degree 7 — so this inherits that, rather
    /// than being a second unverified thing agreeing with a first.
    #[test]
    fn the_operator_route_agrees_with_the_branching_route() {
        for n in 0..=6u32 {
            for mu in crate::partitions_of(n) {
                let want = column_expansion::<Rational>(&mu);
                let got = column_via_operator::<Rational>(&mu);
                for lambda in crate::partitions_of(n) {
                    assert_eq!(
                        got.coeff(&lambda),
                        want.coeff(&lambda),
                        "K_{{{lambda},{mu}}}"
                    );
                }
            }
        }
    }

    /// `schur_in_j_table` is the matrix inverse of the `J → s` expansion.
    ///
    /// The proposition the projection route has to earn, and the only check
    /// that touches both sides: `J_μ` comes off the branching formula in the
    /// monomial basis, while the table is built from the Bergeron–Haiman
    /// recursion and a scalar product. Nothing between them is shared, so a
    /// wrong hook product, a wrong plethysm or a transposed index all show up
    /// here as an off-diagonal entry that fails to cancel.
    #[test]
    fn the_schur_table_inverts_the_j_expansion() {
        for n in 0..=5u32 {
            let parts = crate::partitions_of(n);
            let w = schur_in_j_table::<Rational>(n);
            // b[j][i] = the coefficient of s_{λⁱ} in J_{λʲ}.
            let b: Vec<Schur<Frac<Rational>>> = parts
                .iter()
                .map(|mu| crate::macdonald_j::<Rational>(mu).to_schur())
                .collect();
            for (i, lambda) in parts.iter().enumerate() {
                for (l, nu) in parts.iter().enumerate() {
                    let mut sum = <Frac<Rational> as Ring>::zero();
                    for j in 0..parts.len() {
                        sum.add_assign(&w[i][j].mul(&b[j].coeff(nu)));
                    }
                    let want = if i == l {
                        <Frac<Rational> as Ring>::one()
                    } else {
                        <Frac<Rational> as Ring>::zero()
                    };
                    assert_eq!(sum, want, "row {lambda} against J of {nu}");
                }
            }
        }
    }

    /// The diagonal is `1/c_λ`, not `1/c'_λ`.
    ///
    /// `P_λ = s_λ + (lower in dominance)` and `J_λ = c_λ P_λ`, so this is
    /// forced — and it is the entry a swapped hook product survives everywhere
    /// else, since `c` and `c'` are exchanged by conjugating λ *and* swapping
    /// the variables, which most of the checks here are symmetric under.
    /// `λ = (2)` alone separates them: `c = (1−t)(1−qt)` against
    /// `c' = (1−q)(1−q²)`.
    #[test]
    fn the_diagonal_is_one_over_c_lambda() {
        for n in 0..=5u32 {
            let parts = crate::partitions_of(n);
            let w = schur_in_j_table::<Rational>(n);
            for (i, lambda) in parts.iter().enumerate() {
                let inv: BTreeMap<(u32, u32), i32> = crate::macdonald::c_factors(lambda.parts())
                    .into_iter()
                    .map(|(e, m)| (e, -m))
                    .collect();
                assert_eq!(w[i][i], Frac::from_factors(&inv), "the {lambda} diagonal");
            }
        }
    }

    fn set_q_to_zero(p: &QtPoly<Rational>) -> QtPoly<Rational> {
        let mut out = QtPoly::zero();
        for (&(a, b), c) in p.terms() {
            if a == 0 {
                out.add_term(0, b, *c);
            }
        }
        out
    }

    fn swap_variables(p: &QtPoly<Rational>) -> QtPoly<Rational> {
        let mut out = QtPoly::zero();
        for (&(a, b), c) in p.terms() {
            out.add_term(b, a, *c);
        }
        out
    }
    fn s_in_j(lambda: &[u32]) -> BTreeMap<Partition, Frac<Rational>> {
        let one = QtPoly::term(0, 0, Rational::from_int(1));
        schur_to_macdonald_j(&Schur::monomial(part(lambda), one))
    }

    /// `1 / ∏ (1 − qᵃtᵇ)^m` from the factors and their multiplicities.
    fn over(factors: &[((u32, u32), i32)]) -> Frac<Rational> {
        Frac::from_factors(&factors.iter().map(|&(e, m)| (e, -m)).collect())
    }

    /// The degree-2 expansions, confirmed against Sage
    /// (`SAGE_DISABLE_SYMFN=1`, `Sym.macdonald().J()(s[2])`).
    ///
    /// The denominators are the hook product `c_μ`, not `c'_μ`, which is the
    /// twist `the_diagonal_is_one_over_c_lambda` guards on the table and this
    /// guards on the element-wise form: at μ = (11) they are `(1−t)(1−t²)`
    /// against `(1−q)(1−q²)`. The off-diagonal entry is what a transposed
    /// index would move — `s_2` reaches `J_11`, and `s_11` does not reach
    /// `J_2`.
    #[test]
    fn s2_and_s11_in_j_are_the_hand_values() {
        let mut t_minus_q: Q = QtPoly::term(0, 1, Rational::from_int(1));
        t_minus_q.add_term(1, 0, Rational::from_int(-1));

        let in_j = s_in_j(&[2]);
        assert_eq!(in_j.len(), 2, "s_2 reaches both shapes");
        assert_eq!(
            in_j[&part(&[2])],
            over(&[((0, 1), 1), ((1, 1), 1)]),
            "s_2 at J_2"
        );
        assert_eq!(
            in_j[&part(&[1, 1])],
            over(&[((0, 1), 1), ((0, 2), 1), ((1, 1), 1)]).mul(&Frac::from_poly(t_minus_q)),
            "s_2 at J_11"
        );

        let in_j = s_in_j(&[1, 1]);
        assert_eq!(in_j.len(), 1, "s_11 reaches J_11 alone");
        assert_eq!(
            in_j[&part(&[1, 1])],
            over(&[((0, 1), 1), ((0, 2), 1)]),
            "s_11 at J_11"
        );
    }

    /// The memo hands back the table that was computed, and one degree's entry
    /// is per coefficient ring rather than shared across them.
    ///
    /// The sharing would be silent: `schur_in_macdonald_j` escalates from a
    /// guarded fixed-width pass to `BigRational`, and a cache that ignored the
    /// ring would hand the wide pass the narrow pass's values and make the
    /// escalation nominal (`docs/record/qt-kostka.md`).
    #[test]
    fn the_cached_table_is_the_computed_table_for_each_ring() {
        use crate::guard::GuardedRat;

        for n in 0..=4u32 {
            let want = schur_in_j_table_uncached::<Rational>(n);
            crate::clear_caches();
            assert_eq!(schur_in_j_table::<Rational>(n), want, "cold at degree {n}");
            assert_eq!(schur_in_j_table::<Rational>(n), want, "warm at degree {n}");

            let wide = schur_in_j_table_uncached::<GuardedRat>(n);
            assert_eq!(
                schur_in_j_table::<GuardedRat>(n),
                wide,
                "a second ring read the first ring's entry at degree {n}"
            );
        }
        crate::clear_caches();
    }

    /// The element-wise form is the table read by rows.
    ///
    /// The round trip is `the_schur_table_inverts_the_j_expansion`, which the
    /// table already carries; what only this can catch is the wrapper reading
    /// the table transposed, which would still be a plausible triangular
    /// answer.
    #[test]
    fn the_j_expansion_is_the_table_read_by_rows() {
        for n in 0..=5u32 {
            let parts = crate::partitions_of(n);
            let w = schur_in_j_table::<Rational>(n);
            for (i, lambda) in parts.iter().enumerate() {
                let got = s_in_j(lambda.parts());
                for (j, mu) in parts.iter().enumerate() {
                    let want = &w[i][j];
                    if want.is_zero() {
                        assert!(!got.contains_key(mu), "an explicit zero at {lambda}, {mu}");
                    } else {
                        assert_eq!(got[mu], *want, "{lambda} at {mu}");
                    }
                }
            }
        }
    }

    /// Mixed degrees are handled degree by degree, and the zero element gives
    /// the empty map — the contract every inverse expansion in the crate keeps.
    #[test]
    fn the_j_expansion_is_linear_and_takes_mixed_degrees() {
        let mut f: Schur<Q> = Schur::monomial(part(&[2]), <Q as Ring>::one());
        f.add_term(part(&[2, 1]), QtPoly::term(0, 0, Rational::from_int(3)));

        let got = schur_to_macdonald_j(&f);
        let low = s_in_j(&[2]);
        let high = s_in_j(&[2, 1]);
        let three = Frac::from_poly(QtPoly::term(0, 0, Rational::from_int(3)));

        assert_eq!(
            got.len(),
            low.len() + high.len(),
            "degrees 2 and 3 collided"
        );
        for (mu, c) in &low {
            assert_eq!(&got[mu], c, "degree 2 changed at {mu}");
        }
        for (mu, c) in &high {
            assert_eq!(got[mu], c.mul(&three), "degree 3 is not 3x at {mu}");
        }

        assert!(schur_to_macdonald_j(&Schur::<Q>::zero()).is_empty());
    }
}
