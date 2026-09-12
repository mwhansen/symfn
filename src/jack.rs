//! Jack polynomials `P_λ(x; α)`, and the two other normalizations `Q` and `J`.
//!
//! One deformation parameter instead of Macdonald's two, and everything the
//! (q,t) world does in binomials `1 − qᵃtᵇ` this world does in **linear forms
//! `uα + v`** — see [`AFrac`].
//!
//! ## Three engines, sharing nothing but `Partition` and `AFrac`
//!
//! | fn | route | what it is for |
//! |---|---|---|
//! | [`jack_p`] / [`jack_p_lb`] | the Laplace–Beltrami eigenoperator recursion | the engine |
//! | [`jack_p_branching`] | chains of horizontal strips, ψ^α | the cross-check |
//! | [`jack_j_tableaux`] | Knop–Sahi's tableau formula | the reference, and the *positive* route |
//!
//! The three routes worth taking are Knop–Sahi, the Lassalle recurrences and
//! the Laplace–Beltrami eigenoperator, rather than Gram–Schmidt. All three
//! are dealt with in `docs/record/jack.md`; the eigenoperator
//! one wins, and the reason is arithmetic rather than combinatorics — it
//! enumerates nothing at all, and its denominator at every step is a *single*
//! atom of [`AFrac`].
//!
//! E1 and E2 share no mathematics: one is an eigenvalue equation on the
//! monomial expansion, the other is branching/Pieri. That is the
//! shared-nothing cross-check standard `qtkostka.rs` sets with its three routes
//! to the (q,t)-Kostka polynomials. E3 is exponential and stays forever, in
//! [`NaiveLr`](crate::lr::NaiveLr)'s role.
//!
//! ## The three hooks, and the trap
//!
//! For a cell `s` with arm `a` and leg `l`:
//!
//! ```text
//!   lower hook   h_low(s) = α·a     + l + 1        H_λ  = ∏_s h_low(s)
//!   upper hook   h_up(s)  = α·(a+1) + l            H'_λ = ∏_s h_up(s)
//!   [KS] weight  d(s)     = α·(a+1) + l + 1
//! ```
//!
//! ⚠️ **Three near-identical linear hooks circulate and the swap is silent.**
//! `d` is equal to neither of the other two; all three are needed here
//! (`H`/`H'` normalize, `d` weights tableaux); any two agree on enough small
//! cells to pass a careless test; and the literature's `c_λ, c'_λ, j_λ`
//! notation packs them differently per paper.
//! `hooks_are_three_distinct_families` pins them, the way `deltaop` pins
//! coarm/coleg against arm/leg.
//!
//! Note `h_low` and `h_up` are exactly the α-limits of the two binomials in
//! Macdonald's `b_λ(s)`, under `(1 − qᵃtᵇ) ↦ aα + b`, which is why the
//! branching route's ψ factors ([`jack_p_branching`]) are Macdonald's with one
//! type changed. ⚠️
//! The limit is per **atom** and never per coefficient: a coefficient of `P` is
//! a *sum* of ψ-products, and pushing a finished Macdonald table through `q =
//! t^α, t → 1` would need L'Hôpital on every fraction. That dead end is
//! recorded in `docs/record/jack.md`.
//!
//! ## Normalizations
//!
//! ```text
//!   P_λ  monic in m_λ            Q_λ = (H_λ/H'_λ)·P_λ         J_λ = H_λ·P_λ
//!   ⟨P_λ,P_λ⟩ = H'_λ/H_λ         ⟨P_λ,Q_μ⟩ = δ_λμ             ⟨J_λ,J_λ⟩ = H_λH'_λ
//! ```
//!
//! The norms are *closed products of `2|λ|` linear factors* and never compute a
//! pairing: [`jack_norm_j`] returns the multiset, where Sage prices the same
//! table like a full expansion (`docs/record/jack.md`).
//!
//! Sage's equivalents are the `P`, `Q` and `J` bases of `Sym.jack()`
//! (`scripts/check_bindings.py`); the whole-degree unit [`jack_table`] has no
//! Sage entry point.
//!
//! ## References
//!
//! Dumitriu–Edelman–Shuman, *MOPS*,
//! [arXiv:math-ph/0409066](https://arxiv.org/abs/math-ph/0409066) (the LB
//! recursion); Knop–Sahi,
//! [arXiv:q-alg/9610016](https://arxiv.org/abs/q-alg/9610016) Thm 5.1 and
//! Thm 1.1 (tableaux, and the `u_μ`-divisibility law); Macdonald and Stanley
//! 1989 for the branching and normalization facts. Sage is an oracle here and
//! never a source.
//!
//! Macdonald, *Symmetric Functions and Hall Polynomials*, 2nd ed., VI.10:
//! `H_λ` and `H'_λ` above are his `c_λ(α)` and `c'_λ(α)` (10.21);
//! `Q_λ = (H_λ/H'_λ)·P_λ` is (10.16); `J_λ = H_λ·P_λ` is (10.22); `P_λ` is
//! characterized by (10.13) and (10.14); and the branching route's ψ factors
//! are the limits (10.10)–(10.12) of VI (6.24).
//!
//! Stanley, *Some combinatorial properties of Jack symmetric functions*, Adv.
//! Math. **77** (1989): Theorem 5.8, `⟨J_λ, J_λ⟩` as the product of upper and
//! lower hooks that [`jack_norm_j`] returns; Theorem 6.3, the tableau formula
//! for `J_{λ/μ}`, which is the branching route in the `J` normalization; and
//! Conjecture 8.3, that `⟨J_λ J_μ, J_ν⟩ ∈ ℕ[α]`, which [`stanley_table`]
//! tabulates and which his Theorem 6.1 proves for `ν = (n)`. Every formula
//! taken from Macdonald or Stanley was also verified numerically against
//! Sage before being written down, by `scripts/spec_jack_verify.py`.
//!
//! [KS]: https://arxiv.org/abs/q-alg/9610016

// Every `as` here is a row or column index of a shape, or the `i32` arm offset
// that `AFrac`'s atoms are indexed by — all bounded by |λ|. The coefficients
// are `AFrac<C>` and cannot be cast.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::{BTreeMap, HashMap};

use crate::afrac::{AFrac, Linears};
use crate::coeff::{Field, Integral, Rational, Ring};
use crate::convert::{FromSchur, ToSchur};
use crate::macdonald::{arm, count_above, leg};
use crate::partition::Partition;
use crate::sym::{add_at, by_degree, Monomial, PowerSum, SymFn};

// ---------------------------------------------------------------- hooks -----

/// Every cell of `shape`, as its `(arm, leg)`.
fn for_each_cell(shape: &[u32], visit: &mut impl FnMut(u32, u32)) {
    for i in 0..shape.len() {
        for j in 0..shape[i] as usize {
            visit(arm(shape, i, j), leg(shape, i, j));
        }
    }
}

/// `H_λ = ∏_s (α·a(s) + l(s) + 1)`, the **lower** hooks, as a factor multiset.
///
/// This is the scalar that turns `P` into the integral form `J`, and the
/// leading coefficient `[m_λ]J_λ`.
pub fn hook_lower(lambda: &Partition) -> Linears {
    let mut f = Linears::new();
    for_each_cell(lambda.parts(), &mut |a, l| {
        *f.entry((a, l + 1)).or_insert(0) += 1;
    });
    f
}

/// `H'_λ = ∏_s (α·(a(s)+1) + l(s))`, the **upper** hooks.
pub fn hook_upper(lambda: &Partition) -> Linears {
    let mut f = Linears::new();
    for_each_cell(lambda.parts(), &mut |a, l| {
        *f.entry((a + 1, l)).or_insert(0) += 1;
    });
    f
}

/// `⟨J_λ, J_λ⟩_α = H_λ · H'_λ` — a product of `2|λ|` linear factors, returned
/// factored and never expanded.
pub fn jack_norm_j(lambda: &Partition) -> Linears {
    let mut f = hook_lower(lambda);
    for (k, m) in hook_upper(lambda) {
        *f.entry(k).or_insert(0) += m;
    }
    f
}

/// `⟨P_λ, P_λ⟩_α = H'_λ / H_λ`.
pub fn jack_norm_p(lambda: &Partition) -> Linears {
    let mut f = hook_upper(lambda);
    for (k, m) in hook_lower(lambda) {
        *f.entry(k).or_insert(0) -= m;
    }
    f.retain(|_, m| *m != 0);
    f
}

// --------------------------------------------------- E1: Laplace–Beltrami ---

/// `n(λ) = Σ (i−1)·λ_i`, the statistic the eigenvalue is built from.
fn n_stat(parts: &[u32]) -> u32 {
    parts.iter().enumerate().map(|(i, &x)| i as u32 * x).sum()
}

/// Dominance: does `a` dominate `b` (both partitions of the same size)?
fn dominates(a: &[u32], b: &[u32]) -> bool {
    let (mut sa, mut sb) = (0u32, 0u32);
    for i in 0..a.len().max(b.len()) {
        sa += a.get(i).copied().unwrap_or(0);
        sb += b.get(i).copied().unwrap_or(0);
        if sa < sb {
            return false;
        }
    }
    true
}

/// `P_κ(x; α)` in the monomial basis, by the Laplace–Beltrami recursion.
///
/// With `E(ν) := α·n(ν') − n(ν)` — the \[MOPS\] eigenvalue `ρ^α_κ` cleared of
/// its `2/α` prefactor — the expansion `P_κ = Σ_λ c_{κλ} m_λ` satisfies
///
/// ```text
///   c_{κλ} = [ Σ_{(i,j,t)} (λ_i − λ_j + 2t) · c_{κμ} ] / (E(κ) − E(λ)),
/// ```
///
/// summing over **positions** `i < j` in λ and every `t ≥ 1` with `λ_j − t ≥
/// 0`, where `μ = sort(λ + t·e_i − t·e_j)` must satisfy `λ < μ ≤ κ` in
/// dominance. `c_{κκ} = 1`, and rows are filled in any linear extension of
/// reverse dominance — `n(λ)` ascending is one, since moving a box up strictly
/// decreases `n`.
///
/// Two details a prose reading of \[MOPS\] leaves ambiguous, both pinned
/// operationally by `scripts/spec_jack_verify.py` against Sage rather than
/// argued: the sum is over **positions**, so distinct `(i,j,t)` producing the
/// same μ each contribute; and `μ` is re-sorted, so a move can leave the
/// partition's row order.
///
/// The denominator is
/// `E(κ) − E(λ) = (n(κ')−n(λ'))·α + (n(λ)−n(κ))` with **both** integer
/// coefficients positive for `λ < κ` — \[MOPS\] Lemma 2.16 made visible, and a
/// single [`AFrac`] atom. Nothing here enumerates a tableau: one row is `p(n)`
/// coefficients, each a sum over `O(ℓ(λ)²·λ₁)` moves.
pub fn jack_p_lb<C: Integral>(kappa: &Partition) -> Monomial<AFrac<C>> {
    let mut out = Monomial::zero();
    let one = <AFrac<C> as Ring>::one();
    if kappa.is_empty() {
        out.add_term(Partition::default(), one);
        return out;
    }
    let n = kappa.size();
    let target: Vec<u32> = kappa.parts().to_vec();
    let n_kappa = n_stat(&target);
    let n_kappa_conj = n_stat(kappa.conjugate().parts());

    // `partitions_of` is not free and `jack_table` calls this p(n) times, so the
    // shared list comes from the cache; the two eigenvalue statistics are then
    // precomputed once per call rather than per (κ, λ) pair — `conjugate()` for
    // every λ inside the loop was p(n)² transposes per table.
    let mut order: Vec<(Partition, u32, u32)> = crate::memo::partitions_cached(n)
        .iter()
        .map(|p| {
            let conj = n_stat(p.conjugate().parts());
            (p.clone(), n_stat(p.parts()), conj)
        })
        .collect();
    order.sort_by_key(|(_, ns, _)| *ns);

    let mut coeffs: HashMap<Partition, AFrac<C>> = HashMap::new();
    coeffs.insert(kappa.clone(), one.clone());
    out.add_term(kappa.clone(), one);

    let mut moved: Vec<u32> = Vec::new();
    for (lambda, n_lambda, n_lambda_conj) in &order {
        crate::interrupt::poll();
        if lambda == kappa || !dominates(&target, lambda.parts()) {
            continue;
        }
        let lp = lambda.parts();
        let mut acc = <AFrac<C> as Ring>::zero();
        for i in 0..lp.len() {
            for j in i + 1..lp.len() {
                for t in 1..=lp[j] {
                    moved.clear();
                    moved.extend_from_slice(lp);
                    moved[i] += t;
                    moved[j] -= t;
                    let mu = Partition::new(moved.iter().copied());
                    if !dominates(&target, mu.parts()) {
                        continue;
                    }
                    if let Some(c) = coeffs.get(&mu) {
                        let w = i64::from(lp[i] - lp[j] + 2 * t);
                        acc.add_assign(&c.scale_int(w));
                        // Reduced per term, not per coefficient — the opposite
                        // of `Frac`'s policy, because here a failed
                        // cancellation exits on its first step and letting the
                        // atom multiset grow across the sum is the expensive
                        // side. This is the arithmetic `spec_jack_swell.py`
                        // measured, and the measurement is what licensed the
                        // design.
                        acc.reduce();
                    }
                }
            }
        }
        if acc.is_zero() {
            continue;
        }
        let u = n_kappa_conj - n_lambda_conj;
        let v = n_lambda - n_kappa;
        debug_assert!(u > 0 && v > 0, "E(kappa) - E(lambda) must be positive");
        let c = acc.div_linear(u, v);
        coeffs.insert(lambda.clone(), c.clone());
        out.add_term(lambda.clone(), c);
    }
    out
}

// ------------------------------------------------------- E2: branching ------

/// `ψ^α_{λ/μ}` as a multiset of linear forms, for a horizontal strip.
///
/// Cell for cell this is [`macdonald::psi_factors`](crate::macdonald) with the
/// exponent pair reinterpreted: `(a, b)` meant `1 − qᵃtᵇ` there and means
/// `aα + b` here. Both hooks appear, from both shapes:
///
/// ```text
///   ψ^α_{λ/μ} = ∏_{s ∈ rows(λ/μ) \ cols(λ/μ)}
///                 [h_low^μ(s) · h_up^λ(s)] / [h_up^μ(s) · h_low^λ(s)]
/// ```
///
/// over the cells of μ in a row meeting the strip but **not** in a column
/// meeting it. Both conditions matter and both are easy to swap; the hand case
/// `ψ^α_{(2)/(1)} = 2/(α+1)` pins them, and getting either half wrong gives
/// something else.
fn psi_alpha_factors(lam: &[u32], mu: &[u32]) -> Linears {
    let width = lam.iter().copied().max().unwrap_or(0) as usize;
    let met: Vec<bool> = (0..width)
        .map(|j| count_above(lam, j) > count_above(mu, j))
        .collect();

    let mut f = Linears::new();
    for i in 0..mu.len() {
        if lam[i] == mu[i] {
            continue; // row i does not meet the strip
        }
        for j in 0..mu[i] as usize {
            if met[j] {
                continue; // column j does meet it
            }
            let (am, lm) = (arm(mu, i, j), leg(mu, i, j));
            let (al, ll) = (arm(lam, i, j), leg(lam, i, j));
            // h_low^μ(s) / h_up^μ(s)
            *f.entry((am, lm + 1)).or_insert(0) += 1;
            *f.entry((am + 1, lm)).or_insert(0) -= 1;
            // h_up^λ(s) / h_low^λ(s)
            *f.entry((al + 1, ll)).or_insert(0) += 1;
            *f.entry((al, ll + 1)).or_insert(0) -= 1;
        }
    }
    f.retain(|_, m| *m != 0);
    f
}

/// `P_λ(x; α)` in the monomial basis, by the branching formula — the
/// **cross-check** for [`jack_p_lb`], sharing no mathematics with it.
///
/// A semistandard tableau of shape λ and content μ is a chain of horizontal
/// strips, so `[m_μ]P_λ` is a sum over exactly the chains
/// [`charge`](mod@crate::charge) already enumerates. Same chain enumeration,
/// same per-strip cache, and the same multiset encoding of the ψ factors as
/// [`macdonald_p`](crate::macdonald_p) — one function changed.
///
/// Enumeration-bound, unlike E1, so it is the slower route on a whole degree
/// and the natural one for a single coefficient.
pub fn jack_p_branching<C: Integral>(lambda: &Partition) -> Monomial<AFrac<C>> {
    let mut out = Monomial::zero();
    if lambda.is_empty() {
        out.add_term(Partition::default(), <AFrac<C> as Ring>::one());
        return out;
    }
    let target: Vec<u32> = lambda.parts().to_vec();
    let mut cache: HashMap<(Vec<u32>, Vec<u32>), Linears> = HashMap::new();
    let mut acc: Linears = BTreeMap::new();

    for mu in crate::partitions_of(lambda.size()) {
        crate::interrupt::poll();
        let mut total = <AFrac<C> as Ring>::zero();
        let mut chain: Vec<Vec<u32>> = vec![vec![0; target.len()]];
        crate::charge::build(&target, mu.parts(), 0, &mut chain, &mut |ch| {
            acc.clear();
            for k in 1..ch.len() {
                let key = (ch[k].clone(), ch[k - 1].clone());
                let f = cache
                    .entry(key)
                    .or_insert_with(|| psi_alpha_factors(&ch[k], &ch[k - 1]));
                for (&e, &m) in f.iter() {
                    *acc.entry(e).or_insert(0) += m;
                }
            }
            acc.retain(|_, m| *m != 0);
            total.add_assign(&AFrac::from_factors(&acc));
        });
        if !total.is_zero() {
            total.reduce();
            out.add_term(mu, total);
        }
    }
    out
}

// ------------------------------------------------------ E3: Knop–Sahi -------

/// `J_λ(x; α)` by \[KS\] Theorem 5.1 — the exponential reference route, and the
/// only **manifestly positive** one.
///
/// ```text
///   J_λ(x; α) = Σ_{T admissible} d_T(α) x^T,   d_T = ∏_{s critical} (α(a(s)+1) + l(s) + 1)
/// ```
///
/// `T` labels the cells of λ with `1..n`; admissible means `T(i,j) ≠ T(i',j)`
/// for `i' > i` and `T(i,j) ≠ T(i',j−1)` for `i' < i, j > 1`; a cell is
/// *critical* when `j > 1` and it repeats its left neighbor's label.
///
/// `n^{|λ|}` labelings before pruning, so this is `NaiveLr`'s role: the
/// reference implementation kept forever, and the route in which \[KS\] Thm 1.1
/// (`[m_μ]J_λ / u_μ ∈ ℕ[α]`) is manifest rather than a theorem about the
/// output.
pub fn jack_j_tableaux<C: Integral>(lambda: &Partition) -> Monomial<AFrac<C>> {
    let mut out = Monomial::zero();
    if lambda.is_empty() {
        out.add_term(Partition::default(), <AFrac<C> as Ring>::one());
        return out;
    }
    let shape: Vec<u32> = lambda.parts().to_vec();
    let n = lambda.size() as usize;
    let mut offset = Vec::with_capacity(shape.len());
    let mut cells: Vec<(usize, usize)> = Vec::with_capacity(n);
    for (i, &row) in shape.iter().enumerate() {
        offset.push(cells.len());
        for j in 0..row as usize {
            cells.push((i, j));
        }
    }

    let mut st = KsWalk {
        shape: &shape,
        offset: &offset,
        cells: &cells,
        n: n as u32,
        labels: vec![0; n],
        weight: Linears::new(),
        content: vec![0; n],
        totals: HashMap::new(),
    };
    st.walk(0);

    for mu in crate::partitions_of(lambda.size()) {
        let mut key = mu.parts().to_vec();
        key.resize(n, 0);
        if let Some(mut v) = st.totals.remove(&key) {
            v.reduce();
            out.add_term(mu, v);
        }
    }
    out
}

/// The backtracking state for [`jack_j_tableaux`], kept in a struct so the
/// recursion does not need nine parameters.
struct KsWalk<'a, C: Integral> {
    shape: &'a [u32],
    offset: &'a [usize],
    cells: &'a [(usize, usize)],
    n: u32,
    labels: Vec<u32>,
    weight: Linears,
    content: Vec<u32>,
    totals: HashMap<Vec<u32>, AFrac<C>>,
}

impl<C: Integral> KsWalk<'_, C> {
    fn walk(&mut self, k: usize) {
        if k == self.cells.len() {
            let e = self
                .totals
                .entry(self.content.clone())
                .or_insert_with(<AFrac<C> as Ring>::zero);
            e.add_assign(&AFrac::from_factors(&self.weight));
            return;
        }
        let (i, j) = self.cells[k];
        'value: for v in 0..self.n {
            // Cells above in this column, and above in the column to the left,
            // are exactly the already-assigned ones the two rules forbid.
            for ip in 0..i {
                if self.shape[ip] as usize > j && self.labels[self.offset[ip] + j] == v {
                    continue 'value;
                }
                if j > 0 && self.labels[self.offset[ip] + j - 1] == v {
                    continue 'value;
                }
            }
            self.labels[k] = v;
            self.content[v as usize] += 1;
            let critical = (j > 0 && self.labels[self.offset[i] + j - 1] == v).then(|| {
                let key = (arm(self.shape, i, j) + 1, leg(self.shape, i, j) + 1);
                *self.weight.entry(key).or_insert(0) += 1;
                key
            });
            self.walk(k + 1);
            if let Some(key) = critical {
                let e = self
                    .weight
                    .get_mut(&key)
                    .expect("this frame incremented the critical cell, and only it removes one");
                *e -= 1;
                if *e == 0 {
                    self.weight.remove(&key);
                }
            }
            self.content[v as usize] -= 1;
        }
    }
}

// ------------------------------------------------------- the public API -----

/// `P_λ(x; α)` in the monomial basis: monic in `m_λ` and triangular in
/// dominance order.
///
/// The empty partition gives the single term `∅ ↦ 1`.
///
/// Dispatches to [`jack_p_lb`]: the eigenoperator route wins the whole-degree
/// unit by a margin that grows with the degree (`docs/record/jack.md`, "The
/// two engines on a whole degree"), which is why [`jack_table`] calls it
/// rather than the branching formula.
pub fn jack_p<C: Integral>(lambda: &Partition) -> Monomial<AFrac<C>> {
    jack_p_lb(lambda)
}

/// Multiply every coefficient by a product of linear-form powers, keeping the
/// scalar *factored* the whole way — the [`Frac::mul_factors`] pattern.
///
/// [`Frac::mul_factors`]: crate::frac::Frac::mul_factors
fn scale_by<C: Integral>(f: Monomial<AFrac<C>>, by: &Linears) -> Monomial<AFrac<C>> {
    let mut out = Monomial::zero();
    for (mu, c) in f.terms() {
        let mut v = c.mul_factors(by);
        v.reduce();
        out.add_term(mu.clone(), v);
    }
    out
}

/// `Q_λ = (H_λ/H'_λ)·P_λ`, the basis dual to `P` under `⟨·,·⟩_α`.
pub fn jack_q<C: Integral>(lambda: &Partition) -> Monomial<AFrac<C>> {
    let mut f = hook_lower(lambda);
    for (k, m) in hook_upper(lambda) {
        *f.entry(k).or_insert(0) -= m;
    }
    f.retain(|_, m| *m != 0);
    scale_by(jack_p(lambda), &f)
}

/// `J_λ = H_λ·P_λ`, the integral form.
///
/// Every coefficient is a *polynomial* in α — with non-negative integer
/// coefficients, and divisible by `u_μ = ∏ m_i(μ)!` (\[KS\] Thm 1.1). None of
/// that is arranged here: the coefficients arrive through fraction arithmetic
/// and the denominators cancel, which is why the tests that check it are real
/// checks on the whole route.
pub fn jack_j<C: Integral>(lambda: &Partition) -> Monomial<AFrac<C>> {
    scale_by(jack_p(lambda), &hook_lower(lambda))
}

/// Every `P_λ` of degree `n` — the unit of work the walls in
/// `docs/record/jack.md` are measured in, and the one Sage has no entry
/// point for.
pub fn jack_table<C: Integral>(n: u32) -> Vec<(Partition, Monomial<AFrac<C>>)> {
    crate::partitions_of(n)
        .into_iter()
        .map(|l| {
            let p = jack_p(&l);
            (l, p)
        })
        .collect()
}

/// Every `J_λ` of degree `n`.
pub fn jack_j_table<C: Integral>(n: u32) -> Vec<(Partition, Monomial<AFrac<C>>)> {
    crate::partitions_of(n)
        .into_iter()
        .map(|l| {
            let j = jack_j(&l);
            (l, j)
        })
        .collect()
}

// ------------------------------------------------ the inverse direction -----

/// `f`, given in the monomial basis, rewritten in the Jack `P` basis: the
/// `c_λ` of `f = Σ_λ c_λ P_λ(x; α)`.
///
/// `P` is monic and dominance-unitriangular in the monomial basis, so
/// `m_λ = P_λ − Σ_{μ ◁ λ} c_{λμ} m_μ` solves downward and the whole `m → P`
/// transition is a back-substitution through [`jack_table`] — the same
/// coefficients `P → m` carries, resolved the other way, and no new
/// enumeration. This is [`monomial_to_macdonald_p`](crate::monomial_to_macdonald_p)
/// with `AFrac` in place of `Frac`.
///
/// `f` may mix degrees — each degree's table is applied to its own terms — and
/// the zero element gives the empty map. The result is keyed by partition in
/// the element order and holds no zeros. It is a plain map because the crate
/// has no `P`-basis type, and a [`Monomial`] holding `P`-coefficients would be
/// the confusion the basis types exist to prevent.
///
/// Costs one [`jack_table`] per degree present in `f`, memoized per degree and
/// ring. Sage's equivalent is `Sym.jack().P()(f)`.
///
/// ```
/// use symfn::{monomial_to_jack_p, AFrac, Monomial, Partition, Rational, Ring, SymFn};
///
/// type F = AFrac<Rational>;
/// let m2: Monomial<F> = Monomial::monomial(Partition::new([2]), <F as Ring>::one());
/// let in_p = monomial_to_jack_p(&m2);
///
/// assert_eq!(in_p[&Partition::new([2])], <F as Ring>::one());
/// assert_eq!(
///     in_p[&Partition::new([1, 1])],
///     <F as Ring>::from_i64(-2).div_linear(1, 1),
/// );
/// ```
///
/// So `m_2 = P_2 − [2/(α+1)] P_11`, which is `P_2` read backwards: the
/// coefficient `P → m` puts on the dominance-smaller shape comes back negated.
/// ⚠️ Sending `α → 1/α` — the direction Jack duality runs in, and the twist to
/// check for — would give `−2α/(α+1)` instead. At `α = 1` both are `−1`, which
/// is why the `P_λ(x; 1) = s_λ` check cannot see the difference on its own.
pub fn monomial_to_jack_p<C: Integral + Send + Sync + 'static>(
    f: &Monomial<AFrac<C>>,
) -> BTreeMap<Partition, AFrac<C>> {
    let mut out: BTreeMap<Partition, AFrac<C>> = BTreeMap::new();
    for (n, terms) in by_degree(f) {
        let parts = crate::memo::partitions_cached(n);
        let index: HashMap<&Partition, usize> =
            parts.iter().enumerate().map(|(i, p)| (p, i)).collect();
        let table = cached_in_p_table::<C>(n);
        for (mu, c) in terms {
            crate::interrupt::poll();
            for (lambda, v) in &table[index[mu]] {
                add_at(&mut out, lambda, v.mul(c));
            }
        }
    }
    // Reduced once, at the end: `AFrac::add_assign` deliberately leaves the
    // running sum unreduced, because a cancellation can only be decided when
    // the sum is complete (`src/afrac.rs`).
    for v in out.values_mut() {
        v.reduce();
    }
    out
}

/// `f`, given in the monomial basis, rewritten in the Jack `Q` basis: the
/// `c_λ` of `f = Σ_λ c_λ Q_λ(x; α)`.
///
/// `Q_λ = (H_λ/H'_λ)·P_λ`, so this is [`monomial_to_jack_p`] with each
/// coefficient multiplied by that shape's `H'_λ/H_λ` — which is
/// [`jack_norm_p`], the squared norm `⟨P_λ, P_λ⟩_α`, applied factored so no
/// second solve happens. Same contract as [`monomial_to_jack_p`]; Sage's
/// equivalent is `Sym.jack().Q()(f)`.
///
/// ```
/// use symfn::{monomial_to_jack_q, AFrac, Monomial, Partition, Rational, Ring, SymFn};
///
/// type F = AFrac<Rational>;
/// let m11: Monomial<F> = Monomial::monomial(Partition::new([1, 1]), <F as Ring>::one());
/// let in_q = monomial_to_jack_q(&m11);
///
/// assert_eq!(in_q.len(), 1);
/// // α(α+1)/2
/// let want = F::linear(1, 0).mul(&F::linear(1, 1)).div_int(2);
/// assert_eq!(in_q[&Partition::new([1, 1])], want);
/// ```
///
/// So `m_11 = [α(α+1)/2] Q_11`, where [`monomial_to_jack_p`] gives
/// `m_11 = P_11` outright — the smallest shape at which the two normalizations
/// differ. ⚠️ At `α = 1` this coefficient is `1`, exactly as `P`'s is, so a
/// test that only sets `α = 1` cannot tell `Q` from `P` either.
pub fn monomial_to_jack_q<C: Integral + Send + Sync + 'static>(
    f: &Monomial<AFrac<C>>,
) -> BTreeMap<Partition, AFrac<C>> {
    let mut out = monomial_to_jack_p(f);
    for (lambda, v) in &mut out {
        *v = v.mul_factors(&jack_norm_p(lambda));
        v.reduce();
    }
    out
}

/// `f`, given in the monomial basis, rewritten in the Jack `J` basis: the
/// `c_λ` of `f = Σ_λ c_λ J_λ(x; α)`.
///
/// `J_λ = H_λ·P_λ`, so this is [`monomial_to_jack_p`] with each coefficient
/// divided by that shape's lower hooks ([`hook_lower`]), applied factored. Same
/// contract as [`monomial_to_jack_p`]; Sage's equivalent is
/// `Sym.jack().J()(f)`.
///
/// The coefficients here are *not* the polynomials in α that `J → m` has:
/// dividing by `H_λ` puts the hooks in a denominator, and only the forward
/// direction is integral.
///
/// ```
/// use symfn::{monomial_to_jack_j, AFrac, Monomial, Partition, Rational, Ring, SymFn};
///
/// type F = AFrac<Rational>;
/// let m2: Monomial<F> = Monomial::monomial(Partition::new([2]), <F as Ring>::one());
/// let in_j = monomial_to_jack_j(&m2);
///
/// let over = F::inv_linear(1, 1); // 1/(α+1)
/// assert_eq!(in_j[&Partition::new([2])], over);
/// assert_eq!(in_j[&Partition::new([1, 1])], over.neg());
/// ```
///
/// So `m_2 = [J_2 − J_11]/(α+1)`, where `J_2 = (α+1)m_2 + 2m_11` and
/// `J_11 = 2m_11` — the convention gate read backwards.
pub fn monomial_to_jack_j<C: Integral + Send + Sync + 'static>(
    f: &Monomial<AFrac<C>>,
) -> BTreeMap<Partition, AFrac<C>> {
    let mut out = monomial_to_jack_p(f);
    for (lambda, v) in &mut out {
        let mut over_h = hook_lower(lambda);
        for m in over_h.values_mut() {
            *m = -*m;
        }
        *v = v.mul_factors(&over_h);
        v.reduce();
    }
    out
}

/// [`monomial_in_p_table`], memoized per ring and degree.
///
/// The solve is nearly all of what an `m → P` call costs and its unit is the
/// degree, while the entry point is asked for one element — so without this,
/// sweeping p(n) shapes rebuilds it p(n) times. The Macdonald pair measured
/// that (`docs/record/macdonald.md`); this one was built with the memoization
/// in place from the start, and `docs/record/jack.md` has what it is worth
/// here.
fn cached_in_p_table<C: Integral + Send + Sync + 'static>(
    n: u32,
) -> std::sync::Arc<Vec<BTreeMap<Partition, AFrac<C>>>> {
    crate::memo::jack_p_inverse_cached(n, || monomial_in_p_table::<C>(n))
}

/// The `m → P` transition for degree `n`: entry `j` is `m_{parts[j]}` written
/// in the `P` basis, keyed by partition.
///
/// [`jack_table`] inverted by back-substitution.
/// [`partitions_of`](crate::partitions_of) is lex-descending and λ ⊵ μ implies
/// λ ≥ μ lexicographically, so the dominance-smallest shape is the *last*
/// index and counting the index down solves every `m_μ` before the sums that
/// need it.
fn monomial_in_p_table<C: Integral>(n: u32) -> Vec<BTreeMap<Partition, AFrac<C>>> {
    let parts = crate::memo::partitions_cached(n);
    let table = jack_table::<C>(n);
    let mut out: Vec<BTreeMap<Partition, AFrac<C>>> = vec![BTreeMap::new(); parts.len()];
    for j in (0..parts.len()).rev() {
        crate::interrupt::poll();
        let mut acc = BTreeMap::new();
        acc.insert(parts[j].clone(), <AFrac<C> as Ring>::one());
        for l in (j + 1)..parts.len() {
            let c = table[j].1.coeff(&parts[l]);
            if c.is_zero() {
                continue;
            }
            for (nu, v) in &out[l] {
                add_at(&mut acc, nu, c.mul(v).neg());
            }
        }
        for v in acc.values_mut() {
            v.reduce();
        }
        out[j] = acc;
    }
    out
}

// ------------------------------------------ arithmetic in the basis ---------

/// `f + g`, both given as coefficients in one of the Jack bases.
///
/// Which basis is not asked and does not matter: addition is termwise in
/// whatever basis both are written in, and mixing two of them is the caller's
/// error to avoid. Coefficients are reduced, so the result is the
/// representation every other entry point returns — which is what the Python
/// layer's structural `==` on a sum depends on.
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::jack::jack_element_add;
/// use symfn::{AFrac, Partition, Rational, Ring};
///
/// type F = AFrac<Rational>;
/// let one: BTreeMap<Partition, F> =
///     [(Partition::new([2]), <F as Ring>::one())].into_iter().collect();
/// let minus: BTreeMap<Partition, F> =
///     [(Partition::new([2]), <F as Ring>::one().neg())].into_iter().collect();
///
/// assert!(jack_element_add(&one, &minus).is_empty());
/// assert_eq!(jack_element_add(&one, &BTreeMap::new()), one);
/// ```
///
/// A shape whose coefficients cancel leaves no entry at all.
pub fn jack_element_add<C: Integral>(
    f: &BTreeMap<Partition, AFrac<C>>,
    g: &BTreeMap<Partition, AFrac<C>>,
) -> BTreeMap<Partition, AFrac<C>> {
    let mut out = f.clone();
    for (mu, c) in g {
        add_at(&mut out, mu, c.clone());
    }
    for v in out.values_mut() {
        v.reduce();
    }
    out.retain(|_, v| !v.is_zero());
    out
}

/// `c·f`, `f` given as coefficients in one of the Jack bases.
///
/// Same basis-blindness as [`jack_element_add`]. Reducing matters in two ways
/// here and not only one: an atom `u·α + v` can cancel, and so can the
/// integer `scale` an `AFrac` carries, which no other type in the crate has.
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::jack::jack_element_scale;
/// use symfn::{AFrac, Partition, Rational, Ring};
///
/// type F = AFrac<Rational>;
/// let f: BTreeMap<Partition, F> =
///     [(Partition::new([2]), F::inv_linear(1, 1))].into_iter().collect();
/// let scaled = jack_element_scale(&f, &F::linear(1, 1));
///
/// assert_eq!(scaled[&Partition::new([2])], <F as Ring>::one());
/// ```
///
/// So `(α+1)·[1/(α+1)]` is 1 and not itself over itself.
pub fn jack_element_scale<C: Integral>(
    f: &BTreeMap<Partition, AFrac<C>>,
    c: &AFrac<C>,
) -> BTreeMap<Partition, AFrac<C>> {
    let mut out = BTreeMap::new();
    for (mu, v) in f {
        let mut w = v.mul(c);
        w.reduce();
        if !w.is_zero() {
            out.insert(mu.clone(), w);
        }
    }
    out
}

// ------------------------------------------- back to the monomial basis -----

/// `Σ_λ c_λ · one(λ)`, the expansion shared by the three normalizations.
///
/// Coefficients accumulate unreduced and are reduced once at the end, for the
/// reason [`monomial_to_jack_p`] gives: an `AFrac` cancellation can only be
/// decided when the sum is complete.
fn expand_jack<C: Integral>(
    f: &BTreeMap<Partition, AFrac<C>>,
    one: fn(&Partition) -> Monomial<AFrac<C>>,
) -> Monomial<AFrac<C>> {
    let mut acc: BTreeMap<Partition, AFrac<C>> = BTreeMap::new();
    for (lambda, c) in f {
        crate::interrupt::poll();
        for (mu, v) in one(lambda).terms() {
            add_at(&mut acc, mu, v.mul(c));
        }
    }
    let mut out = Monomial::zero();
    for (mu, mut v) in acc {
        v.reduce();
        out.add_term(mu, v);
    }
    out
}

/// The `P`-basis element `f = Σ_λ c_λ P_λ(x; α)`, expanded in the monomial
/// basis.
///
/// The inverse of [`monomial_to_jack_p`], and its input is that function's
/// output: a plain map from partition to coefficient, because the crate has no
/// `P`-basis type. Shapes of different degrees may be mixed and the empty map
/// gives zero.
///
/// One [`jack_p`] per shape present — not per shape of the degree, which is
/// what makes this the right route for an element with few terms and
/// [`jack_table`] the right one for a whole degree.
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::{jack_p_to_monomial, AFrac, Partition, Rational, Ring, SymFn};
///
/// type F = AFrac<Rational>;
/// let f: BTreeMap<Partition, F> =
///     [(Partition::new([2]), <F as Ring>::one())].into_iter().collect();
/// let m = jack_p_to_monomial(&f);
///
/// assert_eq!(m.coeff(&Partition::new([2])), <F as Ring>::one());
/// assert_eq!(
///     m.coeff(&Partition::new([1, 1])),
///     <F as Ring>::from_i64(2).div_linear(1, 1),
/// );
/// ```
///
/// So `P_2 = m_2 + [2/(α+1)] m_11`. ⚠️ Under `α → 1/α` — the direction Jack
/// duality runs in — the coefficient would be `2α/(α+1)`, and at `α = 1` both
/// are `1`, so a Schur specialization cannot tell them apart.
pub fn jack_p_to_monomial<C: Integral>(f: &BTreeMap<Partition, AFrac<C>>) -> Monomial<AFrac<C>> {
    expand_jack(f, jack_p::<C>)
}

/// The `Q`-basis element `f = Σ_λ c_λ Q_λ(x; α)`, expanded in the monomial
/// basis.
///
/// The inverse of [`monomial_to_jack_q`]; same contract as
/// [`jack_p_to_monomial`].
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::{jack_q_to_monomial, AFrac, Partition, Rational, Ring, SymFn};
///
/// type F = AFrac<Rational>;
/// let f: BTreeMap<Partition, F> =
///     [(Partition::new([1, 1]), <F as Ring>::one())].into_iter().collect();
/// let m = jack_q_to_monomial(&f);
///
/// // 2/(α·(α+1))
/// let want = <F as Ring>::from_i64(2)
///     .div_linear(1, 0)
///     .div_linear(1, 1);
/// assert_eq!(m.coeff(&Partition::new([1, 1])), want);
/// ```
///
/// So `Q_11 = [2/(α(α+1))] m_11`, where [`jack_p_to_monomial`] has
/// `P_11 = m_11` outright — the smallest shape at which the two
/// normalizations differ.
pub fn jack_q_to_monomial<C: Integral>(f: &BTreeMap<Partition, AFrac<C>>) -> Monomial<AFrac<C>> {
    expand_jack(f, jack_q::<C>)
}

/// The `J`-basis element `f = Σ_λ c_λ J_λ(x; α)`, expanded in the monomial
/// basis.
///
/// The inverse of [`monomial_to_jack_j`]; same contract as
/// [`jack_p_to_monomial`]. This is the one direction whose coefficients are
/// *polynomials* in α — \[KS\] Thm 1.1 — so a fraction surviving here is a
/// defect, not a normalization.
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::{jack_j_to_monomial, AFrac, Partition, Rational, Ring, SymFn};
///
/// type F = AFrac<Rational>;
/// let f: BTreeMap<Partition, F> =
///     [(Partition::new([2]), <F as Ring>::one())].into_iter().collect();
/// let m = jack_j_to_monomial(&f);
///
/// assert_eq!(m.coeff(&Partition::new([2])), F::linear(1, 1));
/// assert_eq!(m.coeff(&Partition::new([1, 1])), <F as Ring>::from_i64(2));
/// ```
///
/// So `J_2 = (α+1) m_2 + 2 m_11`, the convention gate this family is pinned
/// by, read forwards.
pub fn jack_j_to_monomial<C: Integral>(f: &BTreeMap<Partition, AFrac<C>>) -> Monomial<AFrac<C>> {
    expand_jack(f, jack_j::<C>)
}

/// `J_λ` in the **power-sum** basis — the Jack-character unit, and what the
/// \[GJ\] pipeline consumes.
///
/// Routed `m → s → p` through the [`convert`](mod@crate::convert) hub, which is
/// generic over the coefficient ring. The `z_ν` divisions in
/// `PowerSum::from_schur` ask for [`QAlgebra`](crate::coeff::QAlgebra), and
/// `AFrac<C>` is one for *any* `C` — including `i128`, which is not.
pub fn jack_j_powersum<C: Integral>(lambda: &Partition) -> PowerSum<AFrac<C>> {
    PowerSum::<AFrac<C>>::from_schur(&jack_j::<C>(lambda).to_schur())
}

/// Every `J_λ` of degree `n`, in the power-sum basis.
pub fn jack_powersum_table<C: Integral>(n: u32) -> Vec<(Partition, PowerSum<AFrac<C>>)> {
    crate::partitions_of(n)
        .into_iter()
        .map(|l| {
            let p = jack_j_powersum(&l);
            (l, p)
        })
        .collect()
}

/// `⟨f, g⟩_α` for two power-sum elements, where the form is diagonal:
/// `⟨p_λ, p_μ⟩_α = δ_λμ · z_λ · α^{ℓ(λ)}`.
pub fn powersum_scalar<C: Integral>(f: &PowerSum<AFrac<C>>, g: &PowerSum<AFrac<C>>) -> AFrac<C> {
    let mut out = <AFrac<C> as Ring>::zero();
    for (mu, a) in f.terms() {
        let Some(b) = g.terms().get(mu) else { continue };
        let mut alpha_pow = Linears::new();
        alpha_pow.insert((1, 0), mu.len() as i32); // α^{ℓ(μ)}
                                                   // `z_in`, not `from_u128(mu.z())`: this runs inside a `guarded` scope
                                                   // (`jack_escalate`), and `z` forms z_μ in native `u128`, which *panics*
                                                   // past |μ| = 34 instead of reporting — a wall the escalation ladder
                                                   // cannot catch, since it is watching for `None`. Accumulating in the
                                                   // ring instead makes every factor a `Guarded` multiply, so the same
                                                   // input reports and re-runs over `BigInt` (R6).
        let term = a.mul(b).mul(&mu.z_in::<AFrac<C>>());
        out.add_assign(&term.mul_factors(&alpha_pow));
    }
    out.reduce();
    out
}

/// `⟨f, g⟩_α` for arbitrary monomial-basis elements, through the power sums.
pub fn jack_scalar<C: Integral>(f: &Monomial<AFrac<C>>, g: &Monomial<AFrac<C>>) -> AFrac<C> {
    let fp = PowerSum::<AFrac<C>>::from_schur(&f.to_schur());
    let gp = PowerSum::<AFrac<C>>::from_schur(&g.to_schur());
    powersum_scalar(&fp, &gp)
}

/// `⟨J_λ J_μ, J_ν⟩_α` — **Stanley's object**, whose membership in `ℕ[α]` is his
/// 1989 conjecture and still open.
///
/// Computed in the power-sum basis, where the product is a multiset union and
/// the pairing is diagonal, so no basis change of the product is ever formed.
/// `J[3,2,1]²` is the first size at which the conjecture is interesting, and
/// is out of Sage's range (`docs/record/jack.md`).
///
/// A negative coefficient here is a **result to report, not a bug to fix** —
/// the `∇e_n`-positivity and valley-Delta rule, verbatim.
///
/// ⚠️ **For a whole table, use [`stanley_table`].** This recomputes all three
/// p-expansions on every call — the degree-12 table asks for 27 951
/// [`jack_j_powersum`] conversions of 99 distinct values. The single-shot form
/// is the honest primitive and is kept as one; the batch form is what a search
/// driver wants (`docs/record/jack.md`).
pub fn jack_structure_constant<C: Integral>(
    la: &Partition,
    mu: &Partition,
    nu: &Partition,
) -> AFrac<C> {
    if la.size() + mu.size() != nu.size() {
        return <AFrac<C> as Ring>::zero();
    }
    let a = jack_j_powersum::<C>(la);
    let b = jack_j_powersum::<C>(mu);
    let c = jack_j_powersum::<C>(nu);
    powersum_scalar(&a.mul(&b), &c)
}

/// Every `⟨J_λ J_μ, J_ν⟩_α` with `|λ| = |μ| = k` — **Stanley's whole table**,
/// with each p-expansion computed once.
///
/// Zero entries are omitted. The saving is not a constant factor: the table has
/// `p(k)²·p(2k)` entries and only `2p(k) + p(2k)` distinct expansions, so the
/// redundancy grows with `k` (282× at k = 6).
///
/// Deliberately **not** solved by memoizing [`jack_j_powersum`]. The Python
/// boundary runs over [`Guarded`](crate::guard::Guarded) precisely so an
/// overflowing intermediate is *detected*, and a cache filled at `i128` and
/// handed out to other widths would hide exactly that — the
/// `memo::bold_p` hazard, which is documented there for the same reason.
/// Hoisting the loop is the version with no correctness question in it.
///
/// Positivity is Stanley's 1989 conjecture and is **open**: this returns the
/// values, and asserts nothing about them.
pub fn stanley_table<C: Integral>(k: u32) -> Vec<(Partition, Partition, Partition, AFrac<C>)> {
    let small = crate::partitions_of(k);
    let large = crate::partitions_of(2 * k);
    let ps: Vec<PowerSum<AFrac<C>>> = small.iter().map(jack_j_powersum::<C>).collect();
    let qs: Vec<PowerSum<AFrac<C>>> = large.iter().map(jack_j_powersum::<C>).collect();

    let mut out = Vec::new();
    for (i, la) in small.iter().enumerate() {
        for (j, mu) in small.iter().enumerate() {
            // The product is a multiset union in the p basis, so it is formed
            // once per (λ, μ) and paired against every ν — never once per
            // triple, and never changed basis.
            let prod = ps[i].mul(&ps[j]);
            for (l, nu) in large.iter().enumerate() {
                let g = powersum_scalar(&prod, &qs[l]);
                if !g.is_zero() {
                    out.push((la.clone(), mu.clone(), nu.clone(), g));
                }
            }
        }
    }
    out
}

// ------------------------------------------------------------- zonal --------

/// `ω_α f`, where `ω_α := ω ∘ (p_r ↦ α·p_r)`.
///
/// On the power-sum basis `ω(p_μ) = ε_μ p_μ` with `ε_μ = (−1)^{|μ|−ℓ(μ)}`, and
/// the α-twist multiplies the `p_μ` coefficient by `α^{ℓ(μ)}` on top of that
/// sign.
///
/// ⚠️ **The twist is not optional and the failure is silent.** The duality law
/// `ω_α P_λ^{(α)} = Q_{λ'}^{(1/α)}` fails with plain `ω` already at λ = (1),
/// where `ωP_(1) = p_1` and `Q_(1)^{(1/α)} = α p_1`. Anything that checks the
/// law on symmetric shapes, or only up to a scalar, will not notice.
pub fn omega_alpha<C: Integral>(f: &PowerSum<AFrac<C>>) -> PowerSum<AFrac<C>> {
    let mut out = PowerSum::zero();
    for (mu, c) in f.terms() {
        let mut alpha_pow = Linears::new();
        alpha_pow.insert((1, 0), mu.len() as i32); // α^{ℓ(μ)}
        let mut v = c.mul_factors(&alpha_pow);
        if (mu.size() - mu.len() as u32) % 2 == 1 {
            v = v.neg();
        }
        v.reduce();
        out.add_term(mu.clone(), v);
    }
    out
}

/// Substitute a value for α in a whole expansion.
///
/// Returns `None` when α is a pole of one of the coefficients.
pub fn specialize<C: Integral + Field>(f: &Monomial<AFrac<C>>, alpha: &C) -> Option<Monomial<C>> {
    let mut out = Monomial::zero();
    for (mu, c) in f.terms() {
        out.add_term(mu.clone(), c.eval(alpha)?);
    }
    Some(out)
}

/// The zonal polynomial `Z_λ` in \[GJ\]'s normalization: `J_λ` at α = 2.
///
/// # Panics
///
/// Never, for any λ: the Jack hooks are products of `aα + b` with `a, b ≥ 0`
/// not both zero, so α = 2 is a pole of none of them. The `expect` inside is
/// that proof, not a wall — [`specialize`] returns `None` at a pole, for a
/// general α where one is possible.
pub fn zonal_j(lambda: &Partition) -> Monomial<Rational> {
    specialize(&jack_j::<Rational>(lambda), &Rational::from_int(2))
        .expect("alpha = 2 is not a pole of any Jack hook")
}

/// The zonal polynomial in **Sage's** normalization: `P_λ` at α = 2.
///
/// ⚠️ Sage's `zonal()` is `P^{(2)}`, not `J^{(2)}` — measured, not assumed
/// (`scripts/spec_jack_verify.py`), and the two differ by `H_λ(2)`. A fixture
/// that gets this backwards tests the wrong thing and still looks plausible,
/// which is why both are exposed under names that say which is which.
///
/// # Panics
///
/// Never; see [`zonal_j`].
pub fn zonal_p(lambda: &Partition) -> Monomial<Rational> {
    specialize(&jack_p::<Rational>(lambda), &Rational::from_int(2))
        .expect("alpha = 2 is not a pole of any Jack hook")
}

#[cfg(test)]
mod tests {
    use super::*;

    type F = AFrac<Rational>;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn r(n: i128) -> Rational {
        Rational::from_int(n)
    }

    /// The round trip: expanding `P_λ` in the monomial basis and reading it
    /// back gives `P_λ` again, for every shape through degree 7 — and the same
    /// for `Q` and `J`, because a solve that dropped the normalizer would
    /// still pass on `P` alone.
    #[test]
    fn every_jack_polynomial_comes_back_as_itself() {
        for n in 0..=7u32 {
            for lambda in crate::partitions_of(n) {
                let unit: BTreeMap<Partition, F> =
                    [(lambda.clone(), <F as Ring>::one())].into_iter().collect();
                let p: Monomial<F> = jack_p(&lambda);
                assert_eq!(monomial_to_jack_p(&p), unit, "m -> P of P_{lambda}");
                let q: Monomial<F> = jack_q(&lambda);
                assert_eq!(monomial_to_jack_q(&q), unit, "m -> Q of Q_{lambda}");
                let j: Monomial<F> = jack_j(&lambda);
                assert_eq!(monomial_to_jack_j(&j), unit, "m -> J of J_{lambda}");
            }
        }
    }

    /// **Expanding is linear.** The Jack counterpart of
    /// `adding_and_scaling_commute_with_expanding` in `src/macdonald.rs`, and
    /// the `AFrac` case matters on its own: its `reduce` also cancels the
    /// integer `scale`, which `Frac` has no analog of.
    #[test]
    fn adding_and_scaling_commute_with_expanding() {
        let a = Partition::new([2, 1]);
        let b = Partition::new([1, 1, 1]);
        let one = <F as Ring>::one();
        let f: BTreeMap<Partition, F> = [(a.clone(), one.clone())].into_iter().collect();
        let g: BTreeMap<Partition, F> = [(b.clone(), one.clone())].into_iter().collect();

        let mut want: Monomial<F> = jack_p_to_monomial(&f);
        for (mu, c) in jack_p_to_monomial(&g).terms() {
            want.add_term(mu.clone(), c.clone());
        }
        assert_eq!(
            jack_p_to_monomial(&jack_element_add(&f, &g)),
            want,
            "P_{a} + P_{b}"
        );

        let c = F::inv_linear(1, 1);
        let mut scaled: Monomial<F> = Monomial::zero();
        for (mu, v) in jack_p_to_monomial(&f).terms() {
            let mut w = v.mul(&c);
            w.reduce();
            scaled.add_term(mu.clone(), w);
        }
        assert_eq!(
            jack_p_to_monomial(&jack_element_scale(&f, &c)),
            scaled,
            "P_{a}/(alpha + 1)"
        );
    }

    /// The round trip the other way: solving `m_μ` into a normalization and
    /// expanding it back gives `m_μ`. Together with
    /// [`every_jack_polynomial_comes_back_as_itself`] this pins both
    /// composites, which a table that was inverted in only one direction
    /// would not survive.
    #[test]
    fn every_monomial_comes_back_as_itself() {
        for n in 0..=7u32 {
            for mu in crate::partitions_of(n) {
                let f: Monomial<F> = Monomial::monomial(mu.clone(), <F as Ring>::one());
                assert_eq!(
                    jack_p_to_monomial(&monomial_to_jack_p(&f)),
                    f,
                    "P at m_{mu}"
                );
                assert_eq!(
                    jack_q_to_monomial(&monomial_to_jack_q(&f)),
                    f,
                    "Q at m_{mu}"
                );
                assert_eq!(
                    jack_j_to_monomial(&monomial_to_jack_j(&f)),
                    f,
                    "J at m_{mu}"
                );
            }
        }
    }

    /// Mixed degrees expand term by term, each against its own degree's
    /// polynomials — so an element spanning two degrees must give what its
    /// terms give separately.
    #[test]
    fn expanding_is_linear_across_degrees() {
        let a = Partition::new([2, 1]);
        let b = Partition::new([2, 2]);
        let two = <F as Ring>::from_i64(2);
        let mixed: BTreeMap<Partition, F> = [(a.clone(), two.clone()), (b.clone(), two.clone())]
            .into_iter()
            .collect();

        let mut want: Monomial<F> = Monomial::zero();
        for lambda in [&a, &b] {
            for (mu, c) in jack_p::<Rational>(lambda).terms() {
                let mut v = c.mul(&two);
                v.reduce();
                want.add_term(mu.clone(), v);
            }
        }
        assert_eq!(jack_p_to_monomial(&mixed), want, "2·P_{a} + 2·P_{b}");
    }

    /// **The second engine.** `⟨P_λ, Q_μ⟩_α = δ_λμ` makes the `P`-coefficient
    /// of `f` the ratio `⟨f, P_λ⟩_α / ⟨P_λ, P_λ⟩_α`, and that route shares no
    /// mathematics with the back-substitution: it goes `m → s → p` and pairs
    /// diagonally, where the solve never leaves the monomial basis and never
    /// pairs anything (V3, `docs/policies/validation.md`).
    ///
    /// The norm is [`jack_norm_p`], a closed product of `2|λ|` linear forms,
    /// so the ratio costs one [`jack_scalar`] and one factored division.
    #[test]
    fn orthogonality_gives_the_same_coefficients_as_the_solve() {
        for n in 1..=5u32 {
            let parts = crate::partitions_of(n);
            for mu in &parts {
                let f: Monomial<F> = Monomial::monomial(mu.clone(), <F as Ring>::one());
                let solved = monomial_to_jack_p(&f);
                for lambda in &parts {
                    let mut paired = jack_scalar(&f, &jack_p::<Rational>(lambda));
                    let mut over = jack_norm_p(lambda);
                    for m in over.values_mut() {
                        *m = -*m;
                    }
                    paired = paired.mul_factors(&over);
                    paired.reduce();
                    let want = solved
                        .get(lambda)
                        .cloned()
                        .unwrap_or_else(<F as Ring>::zero);
                    assert_eq!(paired, want, "[P_{lambda}] m_{mu}");
                }
            }
        }
    }

    /// **`α = 1` is the Schur point.** `P_λ(x; 1) = s_λ`, so evaluating the
    /// `m → P` coefficients there must give the ordinary `m → s` transition —
    /// an outside check on the whole solve that costs nothing.
    ///
    /// ⚠️ It is blind to the `α → 1/α` twist, which fixes `α = 1`. That is
    /// what [`monomial_to_jack_p`]'s doctest is for.
    #[test]
    fn at_alpha_one_the_p_expansion_is_the_schur_expansion() {
        for n in 0..=6u32 {
            for mu in crate::partitions_of(n) {
                let f: Monomial<F> = Monomial::monomial(mu.clone(), <F as Ring>::one());
                let schur = Monomial::<Rational>::monomial(mu.clone(), r(1)).to_schur();
                let mut at_one: BTreeMap<Partition, Rational> = BTreeMap::new();
                for (lambda, c) in monomial_to_jack_p(&f) {
                    let v = c.eval(&r(1)).expect("α = 1 is not a pole of any Jack hook");
                    add_at(&mut at_one, &lambda, v);
                }
                // A coefficient nonzero in α may still vanish at α = 1, so the
                // supports match only after the zeros are dropped, which
                // `add_at` does.
                assert_eq!(at_one, *schur.terms(), "m_{mu} at α = 1");
            }
        }
    }

    /// The hand values, which is what the round trip cannot give: it is blind
    /// to any error the forward direction shares. `P_2 = m_2 + [2/(α+1)] m_11`
    /// and `P_11 = m_11`, so `m_2 = P_2 − [2/(α+1)] P_11`; `Q_11 =
    /// [2/(α(α+1))] P_11` and `J_λ = H_λ P_λ` give the other two.
    #[test]
    fn the_three_normalizations_have_their_hand_values_at_free_alpha() {
        let one = <F as Ring>::one();
        let m2: Monomial<F> = Monomial::monomial(part(&[2]), one.clone());
        let m11: Monomial<F> = Monomial::monomial(part(&[1, 1]), one.clone());

        let in_p = monomial_to_jack_p(&m2);
        assert_eq!(in_p[&part(&[2])], one, "m_2 is monic in P_2");
        assert_eq!(
            in_p[&part(&[1, 1])],
            <F as Ring>::from_i64(-2).div_linear(1, 1),
            "m_2 in P_11 is −2/(α+1), not −2α/(α+1)"
        );
        assert_eq!(
            monomial_to_jack_p(&m11),
            [(part(&[1, 1]), one.clone())].into_iter().collect(),
            "m_11 = P_11"
        );

        // α(α+1)/2, the reciprocal of Q_11's scalar 2/(α(α+1)).
        let in_q = monomial_to_jack_q(&m11);
        assert_eq!(in_q.len(), 1, "m_11 reaches Q_11 alone");
        assert_eq!(
            in_q[&part(&[1, 1])],
            F::linear(1, 0).mul(&F::linear(1, 1)).div_int(2),
            "m_11 in Q_11"
        );

        // J_2 = (α+1)m_2 + 2m_11 and J_11 = 2m_11, so m_2 = [J_2 − J_11]/(α+1).
        let in_j = monomial_to_jack_j(&m2);
        let over = F::inv_linear(1, 1);
        assert_eq!(in_j[&part(&[2])], over, "m_2 in J_2");
        assert_eq!(in_j[&part(&[1, 1])], over.neg(), "m_2 in J_11");
    }

    /// The memoized table is the computed one, at each ring separately: the
    /// key carries the ring because the Python boundary escalates, and a
    /// ring-blind cache would hand the wide pass the narrow pass's values.
    #[test]
    fn the_cached_jack_table_is_the_computed_table_for_each_ring() {
        use crate::guard::GuardedRat;

        for n in 0..=4u32 {
            let want = monomial_in_p_table::<Rational>(n);
            crate::clear_caches();
            assert_eq!(*cached_in_p_table::<Rational>(n), want, "cold at {n}");
            assert_eq!(*cached_in_p_table::<Rational>(n), want, "warm at {n}");

            let wide = monomial_in_p_table::<GuardedRat>(n);
            assert_eq!(
                *cached_in_p_table::<GuardedRat>(n),
                wide,
                "a second ring read the first ring's entry at degree {n}"
            );
        }
        crate::clear_caches();
    }

    /// All three transitions are linear and take an argument that mixes
    /// degrees; the zero element gives the empty map.
    #[test]
    fn the_jack_expansions_are_linear_and_take_mixed_degrees() {
        let two = <F as Ring>::from_i64(2);
        let three = <F as Ring>::from_i64(3);
        let mut f: Monomial<F> = Monomial::zero();
        f.add_term(part(&[1]), <F as Ring>::one());
        f.add_term(part(&[2]), two.clone());
        f.add_term(part(&[1, 1]), three.clone());

        for (name, to_basis) in [
            (
                "m -> P",
                monomial_to_jack_p as fn(&Monomial<F>) -> BTreeMap<Partition, F>,
            ),
            ("m -> Q", monomial_to_jack_q),
            ("m -> J", monomial_to_jack_j),
        ] {
            let mut want: BTreeMap<Partition, F> = BTreeMap::new();
            for (mu, scalar) in [
                (part(&[1]), <F as Ring>::one()),
                (part(&[2]), two.clone()),
                (part(&[1, 1]), three.clone()),
            ] {
                let piece = to_basis(&Monomial::monomial(mu, <F as Ring>::one()));
                for (lambda, c) in piece {
                    add_at(&mut want, &lambda, c.mul(&scalar));
                }
            }
            for v in want.values_mut() {
                v.reduce();
            }
            let mut got = to_basis(&f);
            for v in got.values_mut() {
                v.reduce();
            }
            assert_eq!(got, want, "{name} is linear across three degrees");
            assert!(to_basis(&Monomial::zero()).is_empty(), "{name} of 0");
        }
    }

    /// **The convention gate.** `J_(2) = (α+1)m_2 + 2m_11` is the line that
    /// dies first if the parameter convention drifts — Sage calls α "t", and
    /// the verify script's first assertion is exactly this.
    #[test]
    fn the_convention_gate() {
        let j: Monomial<F> = jack_j(&part(&[2]));
        assert_eq!(j.coeff(&part(&[2])), F::linear(1, 1), "{j}");
        assert_eq!(j.coeff(&part(&[1, 1])), <F as Ring>::from_i64(2), "{j}");
        assert_eq!(j.terms().len(), 2);
    }

    /// The hand values, including the branching ψ the whole of E2 rests on.
    #[test]
    fn hand_values() {
        assert_eq!(
            jack_j::<Rational>(&part(&[1])).coeff(&part(&[1])),
            <F as Ring>::one()
        );
        assert_eq!(
            jack_j::<Rational>(&part(&[1, 1])).coeff(&part(&[1, 1])),
            <F as Ring>::from_i64(2)
        );
        // P_(2) = m_2 + [2/(α+1)] m_11
        let p: Monomial<F> = jack_p(&part(&[2]));
        assert_eq!(p.coeff(&part(&[2])), <F as Ring>::one());
        let want = <F as Ring>::from_i64(2).div_linear(1, 1);
        assert_eq!(p.coeff(&part(&[1, 1])), want, "{p}");
        // ψ^α_{(2)/(1)} = 2/(α+1): the single cell (0,0) of μ = (1) gives
        // [1 · 2α] / [α · (α+1)].
        let psi = psi_alpha_factors(&[2], &[1]);
        assert_eq!(F::from_factors(&psi), want, "psi = {psi:?}");
    }

    /// ⚠️ `h_low`, `h_up` and the \[KS\] weight `d` are three *different*
    /// linear families, and any two of them agree on enough small cells to pass
    /// a careless test. λ = (2,1) separates all three.
    #[test]
    fn hooks_are_three_distinct_families() {
        let lam = part(&[2, 1]);
        // cells (0,0),(0,1),(1,0) with (arm,leg) = (1,1),(0,0),(0,0)
        // h_low = α·a + l + 1  →  α+2, 1, 1
        assert_eq!(
            hook_lower(&lam),
            Linears::from([((1, 2), 1), ((0, 1), 2)]),
            "H_(21) = (α+2)·1·1"
        );
        // h_up = α(a+1) + l  →  2α+1, α, α
        assert_eq!(
            hook_upper(&lam),
            Linears::from([((2, 1), 1), ((1, 0), 2)]),
            "H'_(21) = (2α+1)·α·α"
        );
        // d = α(a+1) + l + 1  →  2α+2, α+1, α+1 — equal to neither
        assert_ne!(hook_lower(&lam), hook_upper(&lam));
        let d: Linears = Linears::from([((2, 2), 1), ((1, 1), 2)]);
        assert_ne!(d, hook_lower(&lam));
        assert_ne!(d, hook_upper(&lam));
    }

    /// **E1 ≡ E2.** The eigenoperator recursion and the branching formula share
    /// `Partition` and `AFrac` and nothing else — no chains in one, no
    /// eigenvalues in the other — so agreement is evidence rather than
    /// tautology (V3, `docs/policies/validation.md`).
    #[test]
    fn the_two_engines_agree() {
        for n in 0..=7u32 {
            for lambda in crate::partitions_of(n) {
                let a: Monomial<F> = jack_p_lb(&lambda);
                let b: Monomial<F> = jack_p_branching(&lambda);
                assert_eq!(a, b, "P_{lambda}: LB {a} vs branching {b}");
            }
        }
    }

    /// **E3 ≡ both.** The Knop–Sahi tableau sum is exponential, so the range is
    /// short — but it is the only route in which positivity is manifest, and it
    /// shares nothing with either of the others.
    #[test]
    fn the_tableau_formula_agrees() {
        for n in 0..=5u32 {
            for lambda in crate::partitions_of(n) {
                let a: Monomial<F> = jack_j(&lambda);
                let b: Monomial<F> = jack_j_tableaux(&lambda);
                assert_eq!(a, b, "J_{lambda}: scaled-P {a} vs tableaux {b}");
            }
        }
    }

    /// Monic and strictly triangular in dominance order.
    #[test]
    fn expansion_is_unitriangular_in_dominance() {
        for n in 0..=8u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = jack_p(&lambda);
                assert_eq!(p.coeff(&lambda), <F as Ring>::one(), "monic at {lambda}");
                for mu in p.terms().keys() {
                    assert!(
                        dominates(lambda.parts(), mu.parts()),
                        "{mu} appears in P_{lambda} but λ does not dominate it"
                    );
                }
            }
        }
    }

    /// `[m_λ]J_λ = H_λ`, `⟨J,J⟩ = H·H'`, `⟨P,P⟩ = H'/H`, `⟨P_λ,Q_μ⟩ = δ`.
    ///
    /// The pairing goes through `m → s → p`, which divides by `z_ν` — so this
    /// also exercises `AFrac<C>` being a `QAlgebra` when `C` is not.
    #[test]
    fn norms_and_the_dual_basis() {
        for n in 0..=5u32 {
            let parts = crate::partitions_of(n);
            for lambda in &parts {
                let p: Monomial<F> = jack_p(lambda);
                let q: Monomial<F> = jack_q(lambda);
                let j: Monomial<F> = jack_j(lambda);
                assert_eq!(
                    j.coeff(lambda),
                    F::from_factors(&hook_lower(lambda)),
                    "[m_λ]J_{lambda} = H_λ"
                );
                assert_eq!(
                    jack_scalar(&j, &j),
                    F::from_factors(&jack_norm_j(lambda)),
                    "<J,J> at {lambda}"
                );
                assert_eq!(
                    jack_scalar(&p, &p),
                    F::from_factors(&jack_norm_p(lambda)),
                    "<P,P> at {lambda}"
                );
                for mu in &parts {
                    let got = jack_scalar(&p, &jack_q::<Rational>(mu));
                    let want = if lambda == mu {
                        <F as Ring>::one()
                    } else {
                        <F as Ring>::zero()
                    };
                    assert_eq!(got, want, "<P_{lambda}, Q_{mu}>");
                }
                let _ = q;
            }
        }
    }

    /// `⟨p_λ, p_μ⟩_α = δ · z_λ · α^{ℓ(λ)}` — the deformed form itself.
    #[test]
    fn the_deformed_power_sum_pairing() {
        for n in 0..=5u32 {
            let parts = crate::partitions_of(n);
            for la in &parts {
                for mu in &parts {
                    let a: PowerSum<F> = PowerSum::monomial(la.clone(), <F as Ring>::one());
                    let b: PowerSum<F> = PowerSum::monomial(mu.clone(), <F as Ring>::one());
                    let got = powersum_scalar(&a, &b);
                    if la == mu {
                        let mut f = Linears::new();
                        f.insert((1, 0), la.len() as i32);
                        let want = F::from_u128(la.z()).mul_factors(&f);
                        assert_eq!(got, want, "<p_{la}, p_{la}>");
                    } else {
                        assert!(got.is_zero(), "<p_{la}, p_{mu}> must vanish");
                    }
                }
            }
        }
    }

    /// **\[KS\] Thm 1.1, as a law.** Every `[m_μ]J_λ` clears its denominator,
    /// lies in `ℕ[α]`, **and** is divisible by `u_μ = ∏ m_i(μ)!`.
    ///
    /// The coefficients arrive through fraction arithmetic over `AFrac<i128>`,
    /// so every one of these is a check on the whole route rather than on a
    /// formula that already says so. Run over `i128` on purpose: integrality is
    /// then a property of the representation and not of a ℚ that could absorb
    /// anything.
    #[test]
    fn knop_sahi_positivity_and_u_divisibility() {
        for n in 0..=8u32 {
            for lambda in crate::partitions_of(n) {
                let j: Monomial<AFrac<i128>> = jack_j(&lambda);
                for (mu, c) in j.terms() {
                    let poly = c
                        .clone()
                        .into_natural_poly()
                        .unwrap_or_else(|| panic!("[m_{mu}]J_{lambda} = {c} not in N[alpha]"));
                    let mut u: i128 = 1;
                    let mut i = 0;
                    let p = mu.parts();
                    while i < p.len() {
                        let mut mult = 0i128;
                        let v = p[i];
                        while i < p.len() && p[i] == v {
                            mult += 1;
                            i += 1;
                        }
                        u *= (1..=mult).product::<i128>();
                    }
                    for co in &poly {
                        assert_eq!(co % u, 0, "u_{mu} = {u} must divide {co} in {c}");
                    }
                }
            }
        }
    }

    /// **α = 1 is the Schur function.** Checked against the crate's own hook
    /// machinery and `s → m`, which share no code with `jack.rs`:
    /// `J_λ = (∏ hooks) · s_λ`.
    #[test]
    fn at_alpha_one_it_is_the_schur_function() {
        use crate::Schur;
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let j = specialize(&jack_j::<Rational>(&lambda), &r(1)).expect("α = 1 is no pole");
                let mut hook = r(1);
                for_each_cell(lambda.parts(), &mut |a, l| {
                    hook = hook.mul(&r(i128::from(a + l + 1)));
                });
                let s: Monomial<Rational> =
                    Monomial::from_schur(&Schur::monomial(lambda.clone(), hook));
                assert_eq!(j, s, "J_{lambda} at α = 1");
            }
        }
    }

    /// **α = ∞ is the monomial basis, α = 0 is `e_{λ'}`** (\[KS\] p. 1). Taken
    /// as limits of `P`, not values: the leading behavior in α of each
    /// coefficient.
    ///
    /// At α → ∞ every off-diagonal coefficient must vanish, so `P_λ → m_λ`.
    /// Tested by evaluating at a large α and watching the off-diagonal
    /// coefficients shrink, which is the honest finite version of the claim.
    #[test]
    fn at_large_alpha_p_approaches_the_monomial() {
        for lambda in crate::partitions_of(5) {
            let p: Monomial<F> = jack_p(&lambda);
            for (mu, c) in p.terms() {
                if mu == &lambda {
                    continue;
                }
                let small = c.eval(&r(1000)).expect("no pole at α = 1000");
                let big = c.eval(&r(10)).expect("no pole at α = 10");
                assert!(
                    small.numer().abs() * big.denom().abs()
                        < big.numer().abs() * small.denom().abs(),
                    "[{mu}]P_{lambda} must shrink as α grows: {small:?} vs {big:?}"
                );
            }
        }
    }

    /// The principal specialization `J_λ(1^N; α) = ∏_{(i,j)∈λ} (N + α·j − i)`,
    /// against evaluating the expansion directly.
    ///
    /// 0-based coarm/coleg, and the product is over *cells*, so this is a
    /// genuinely different formula from anything the expansion used.
    #[test]
    fn principal_specialization() {
        for n in 0..=5u32 {
            for lambda in crate::partitions_of(n) {
                let j: Monomial<F> = jack_j(&lambda);
                for nvars in lambda.len() as u32..=4 {
                    // Σ_μ (number of monomials x^w with w a rearrangement of μ
                    // into N slots) · coefficient — i.e. m_μ(1^N).
                    let mut got = <F as Ring>::zero();
                    for (mu, c) in j.terms() {
                        let count = arrangements(mu, nvars);
                        got.add_assign(&c.mul(&F::from_u128(count)));
                    }
                    got.reduce();
                    let mut want: Linears = Linears::new();
                    for (i, &row) in lambda.parts().iter().enumerate() {
                        for jj in 0..row {
                            // N + α·j − i, with j the coarm and i the coleg
                            let v = i64::from(nvars) - i as i64;
                            assert!(v >= 0, "N ≥ ℓ(λ) keeps every factor positive");
                            *want.entry((jj, v as u32)).or_insert(0) += 1;
                        }
                    }
                    assert_eq!(got, F::from_factors(&want), "J_{lambda}(1^{nvars})");
                }
            }
        }
    }

    /// `m_μ(1^N)`: the number of distinct rearrangements of μ into `n` slots.
    fn arrangements(mu: &Partition, n: u32) -> u128 {
        let k = mu.len() as u32;
        if k > n {
            return 0;
        }
        // n! / ((n−k)! · ∏ m_i(μ)!)
        let mut out: u128 = ((n - k + 1)..=n).map(u128::from).product();
        let p = mu.parts();
        let mut i = 0;
        while i < p.len() {
            let mut mult = 0u128;
            let v = p[i];
            while i < p.len() && p[i] == v {
                mult += 1;
                i += 1;
            }
            out /= (1..=mult).product::<u128>();
        }
        out
    }

    /// Sage's `zonal()` is `P^{(2)}` and \[GJ\]'s `Z_λ` is `J^{(2)}`; the two
    /// differ by `H_λ(2)`. ⚠️ Recorded as a test because a fixture that gets
    /// this backwards still looks plausible.
    #[test]
    fn the_two_zonal_normalizations_differ_by_the_hook_product() {
        for n in 0..=5u32 {
            for lambda in crate::partitions_of(n) {
                let zj = zonal_j(&lambda);
                let zp = zonal_p(&lambda);
                let h = F::from_factors(&hook_lower(&lambda))
                    .eval(&r(2))
                    .expect("α = 2 is no pole");
                for (mu, c) in zp.terms() {
                    assert_eq!(zj.coeff(mu), c.mul(&h), "Z_{lambda} at {mu}");
                }
                assert_eq!(zp.coeff(&lambda), r(1), "P is monic");
            }
        }
    }

    /// **Stanley's conjecture, observed.** `⟨J_λ J_μ, J_ν⟩_α ∈ ℕ[α]` on every
    /// triple with `|λ|, |μ| ≤ 3`.
    ///
    /// This is an *open* conjecture: a negative coefficient here is a result to
    /// report, never something to fix. The assertion is what makes the engine
    /// say so loudly.
    #[test]
    fn stanley_structure_constants_are_observed_positive() {
        for (na, nb) in [(1u32, 2u32), (2, 2), (2, 3), (3, 3)] {
            for la in crate::partitions_of(na) {
                for mu in crate::partitions_of(nb) {
                    for nu in crate::partitions_of(na + nb) {
                        let g: AFrac<i128> = jack_structure_constant(&la, &mu, &nu);
                        assert!(
                            g.clone().into_natural_poly().is_some(),
                            "<J_{la} J_{mu}, J_{nu}> = {g} left N[alpha] — REPORT, do not fix"
                        );
                    }
                }
            }
        }
    }

    /// The batch table must agree with the single-shot primitive, entry for
    /// entry — including on the entries it *omits*, which is the half a
    /// spot-check would miss.
    ///
    /// `stanley_table` hoists the p-expansions out of the triple loop, and a
    /// hoist that quietly reused the wrong expansion would still produce a
    /// plausible, positive table.
    #[test]
    fn the_batch_stanley_table_agrees_with_the_primitive() {
        for k in 1..=3u32 {
            let table = stanley_table::<i128>(k);
            let seen: std::collections::HashMap<_, _> = table
                .iter()
                .map(|(a, b, c, g)| ((a.clone(), b.clone(), c.clone()), g.clone()))
                .collect();
            let mut checked = 0;
            for la in crate::partitions_of(k) {
                for mu in crate::partitions_of(k) {
                    for nu in crate::partitions_of(2 * k) {
                        let want: AFrac<i128> = jack_structure_constant(&la, &mu, &nu);
                        let key = (la.clone(), mu.clone(), nu.clone());
                        match seen.get(&key) {
                            Some(got) => assert_eq!(*got, want, "at {key:?}"),
                            None => {
                                assert!(want.is_zero(), "the table omits {key:?} but it is {want}")
                            }
                        }
                        checked += 1;
                    }
                }
            }
            assert!(checked > 0);
        }
    }

    /// **The `ω_α` duality:** `ω_α P_λ^{(α)} = Q_{λ'}^{(1/α)}`.
    ///
    /// The one law in `docs/record/jack.md` that no other test
    /// reaches — it is the only statement relating `P` to `Q`, conjugation, and
    /// the parameter inversion at once, so it independently pins `jack_q`'s
    /// normalization (which otherwise only appears in `⟨P,Q⟩ = δ`) and the
    /// α-limit conventions in the hooks.
    ///
    /// ⚠️ Both halves are traps. The α-twist is not optional: with plain `ω`
    /// this fails at λ = (1) already. And it must be checked on
    /// **non-self-conjugate** shapes, or λ' = λ hides the conjugation.
    #[test]
    fn the_omega_alpha_duality() {
        let mut asymmetric = 0;
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let p: PowerSum<F> = PowerSum::from_schur(&jack_p::<Rational>(&lambda).to_schur());
                let left = omega_alpha(&p);
                let conj = lambda.conjugate();
                if conj != lambda {
                    asymmetric += 1;
                }
                let right: PowerSum<F> =
                    PowerSum::from_schur(&jack_q::<Rational>(&conj).to_schur());
                // The left side is at α; the right side is at 1/α. Invert the
                // left rather than the right, so the comparison is exact and
                // never evaluates anything.
                for (mu, c) in left.terms() {
                    assert_eq!(
                        c.invert_alpha(),
                        right.coeff(mu),
                        "omega_alpha P_{lambda} vs Q_{conj} at p_{mu}"
                    );
                }
                for (mu, c) in right.terms() {
                    if left.coeff(mu).is_zero() {
                        assert!(
                            c.is_zero(),
                            "Q_{conj} has p_{mu} where omega_alpha P_{lambda} does not"
                        );
                    }
                }
            }
        }
        assert!(
            asymmetric > 0,
            "the sweep must include λ' ≠ λ, or it tests nothing"
        );
    }

    /// The α-twist really is needed: plain `ω` breaks the duality at the
    /// smallest possible shape.
    #[test]
    fn plain_omega_breaks_the_duality() {
        let one = part(&[1]);
        let p: PowerSum<F> = PowerSum::from_schur(&jack_p::<Rational>(&one).to_schur());
        let q: PowerSum<F> = PowerSum::from_schur(&jack_q::<Rational>(&one).to_schur());
        // ω P_(1) = p_1, but Q_(1)^{(1/α)} = α·p_1.
        let twisted = omega_alpha(&p).coeff(&one).invert_alpha();
        assert_eq!(twisted, q.coeff(&one), "with the twist it holds");
        let plain = p.coeff(&one); // ε = +1 at μ = (1), so plain ω is the identity here
        assert_ne!(
            plain.invert_alpha(),
            q.coeff(&one),
            "without the twist it must FAIL, or this test proves nothing"
        );
    }

    /// `J_λ` in the power-sum basis is the Jack character table; the sanity
    /// check available without an oracle is that it round-trips.
    #[test]
    fn the_power_sum_expansion_round_trips() {
        for n in 0..=5u32 {
            for lambda in crate::partitions_of(n) {
                let j: Monomial<F> = jack_j(&lambda);
                let p: PowerSum<F> = jack_j_powersum(&lambda);
                let back: Monomial<F> = Monomial::from_schur(&p.to_schur());
                assert_eq!(back, j, "J_{lambda} through the p basis");
            }
        }
    }
}
