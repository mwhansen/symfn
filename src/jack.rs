//! Jack polynomials `P_λ(x; α)`, and the two other normalizations `Q` and `J`.
//!
//! One deformation parameter instead of Macdonald's two, and everything the
//! (q,t) world does in binomials `1 − qᵃtᵇ` this world does in **linear forms
//! `uα + v`** — see [`AFrac`], which is the whole reason this is fast.
//!
//! ## Three engines, sharing nothing but `Partition` and `AFrac`
//!
//! | fn | route | what it is for |
//! |---|---|---|
//! | [`jack_p`] / [`jack_p_lb`] | the Laplace–Beltrami eigenoperator recursion | the engine |
//! | [`jack_p_branching`] | chains of horizontal strips, ψ^α | the cross-check |
//! | [`jack_j_tableaux`] | Knop–Sahi's tableau formula | the reference, and the *positive* route |
//!
//! `research-gaps.md` §2.6 asks for "Knop–Sahi, the Lassalle recurrences, or
//! the Laplace–Beltrami eigenoperator, rather than Gram–Schmidt". All three
//! named routes are dealt with in `docs/record/jack.md`; the eigenoperator
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
//! `d` is equal to neither of the other two; all three are needed here (`H`/`H'`
//! normalize, `d` weights tableaux); any two agree on enough small cells to pass
//! a careless test; and the literature's `c_λ, c'_λ, j_λ` notation packs them
//! differently per paper. `hooks_are_three_distinct_families` pins them, the way
//! `deltaop` pins coarm/coleg against arm/leg.
//!
//! Note `h_low` and `h_up` are exactly the α-limits of the two binomials in
//! Macdonald's `b_λ(s)`, under `(1 − qᵃtᵇ) ↦ aα + b`, which is why
//! [`psi_alpha_factors`] is [`macdonald::psi_factors`](crate::macdonald) with
//! one type changed. ⚠️ The limit is per **atom** and never per coefficient: a
//! coefficient of `P` is a *sum* of ψ-products, and pushing a finished Macdonald
//! table through `q = t^α, t → 1` would need L'Hôpital on every fraction. That
//! dead end is recorded in `docs/record/jack.md`.
//!
//! ## Normalizations
//!
//! ```text
//!   P_λ  monic in m_λ            Q_λ = (H_λ/H'_λ)·P_λ         J_λ = H_λ·P_λ
//!   ⟨P_λ,P_λ⟩ = H'_λ/H_λ         ⟨P_λ,Q_μ⟩ = δ_λμ             ⟨J_λ,J_λ⟩ = H_λH'_λ
//! ```
//!
//! The norms are *closed products of `2|λ|` linear factors* and never compute a
//! pairing: [`jack_norm_j`] returns the multiset. Sage prices the same table
//! like a full expansion and needs over 360 s at n = 12
//! (`docs/record/jack.md`).
//!
//! ## Sources
//!
//! Dumitriu–Edelman–Shuman, *MOPS*, [arXiv:math-ph/0409066] (the LB recursion);
//! Knop–Sahi, [arXiv:q-alg/9610016] Thm 5.1 and Thm 1.1 (tableaux, and the
//! `u_μ`-divisibility law); Macdonald VI.10 and Stanley 1989 for the branching
//! and normalization facts — those two were **not** fetched, and every formula
//! taken from that tradition was instead verified numerically against Sage
//! before being written down, by `scripts/spec_jack_verify.py`. Sage is an
//! oracle here and never a source.

use std::collections::{BTreeMap, HashMap};

use crate::afrac::{AFrac, Linears};
use crate::coeff::{Field, Rational, Ring};
use crate::convert::{FromSchur, ToSchur};
use crate::macdonald::{arm, count_above, leg};
use crate::partition::Partition;
use crate::sym::{Monomial, PowerSum, SymFn};

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
/// With `E(ν) := α·n(ν') − n(ν)` — the [MOPS] eigenvalue `ρ^α_κ` cleared of its
/// `2/α` prefactor — the expansion `P_κ = Σ_λ c_{κλ} m_λ` satisfies
///
/// ```text
///   c_{κλ} = [ Σ_{(i,j,t)} (λ_i − λ_j + 2t) · c_{κμ} ] / (E(κ) − E(λ)),
/// ```
///
/// summing over **positions** `i < j` in λ and every `t ≥ 1` with `λ_j − t ≥ 0`,
/// where `μ = sort(λ + t·e_i − t·e_j)` must satisfy `λ < μ ≤ κ` in dominance.
/// `c_{κκ} = 1`, and rows are filled in any linear extension of reverse
/// dominance — `n(λ)` ascending is one, since moving a box up strictly
/// decreases `n`.
///
/// Two details a prose reading of [MOPS] leaves ambiguous, both pinned
/// operationally by `scripts/spec_jack_verify.py` against Sage rather than
/// argued: the sum is over **positions**, so distinct `(i,j,t)` producing the
/// same μ each contribute; and `μ` is re-sorted, so a move can leave the
/// partition's row order.
///
/// The denominator is
/// `E(κ) − E(λ) = (n(κ')−n(λ'))·α + (n(λ)−n(κ))` with **both** integer
/// coefficients positive for `λ < κ` — [MOPS] Lemma 2.16 made visible, and a
/// single [`AFrac`] atom. Nothing here enumerates a tableau: one row is `p(n)`
/// coefficients, each a sum over `O(ℓ(λ)²·λ₁)` moves.
pub fn jack_p_lb<C: Ring>(kappa: &Partition) -> Monomial<AFrac<C>> {
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
/// [`charge`](crate::charge) already enumerates. Same `build`/`strips`, same
/// per-strip cache, same `Factors`-as-multiset encoding as
/// [`macdonald_p`](crate::macdonald_p) — one function changed.
///
/// Enumeration-bound, unlike E1, so it is the slower route on a whole degree
/// and the natural one for a single coefficient.
pub fn jack_p_branching<C: Ring>(lambda: &Partition) -> Monomial<AFrac<C>> {
    let mut out = Monomial::zero();
    if lambda.is_empty() {
        out.add_term(Partition::default(), <AFrac<C> as Ring>::one());
        return out;
    }
    let target: Vec<u32> = lambda.parts().to_vec();
    let mut cache: HashMap<(Vec<u32>, Vec<u32>), Linears> = HashMap::new();
    let mut acc: Linears = BTreeMap::new();

    for mu in crate::partitions_of(lambda.size()) {
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

/// `J_λ(x; α)` by [KS] Theorem 5.1 — the exponential reference route, and the
/// only **manifestly positive** one.
///
/// ```text
///   J_λ(x; α) = Σ_{T admissible} d_T(α) x^T,   d_T = ∏_{s critical} (α(a(s)+1) + l(s) + 1)
/// ```
///
/// `T` labels the cells of λ with `1..n`; admissible means `T(i,j) ≠ T(i',j)`
/// for `i' > i` and `T(i,j) ≠ T(i',j−1)` for `i' < i, j > 1`; a cell is
/// *critical* when `j > 1` and it repeats its left neighbour's label.
///
/// `n^{|λ|}` labelings before pruning, so this is `NaiveLr`'s role: the
/// reference implementation kept forever, and the route in which [KS] Thm 1.1
/// (`[m_μ]J_λ / u_μ ∈ ℕ[α]`) is manifest rather than a theorem about the
/// output.
pub fn jack_j_tableaux<C: Ring>(lambda: &Partition) -> Monomial<AFrac<C>> {
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
struct KsWalk<'a, C: Ring> {
    shape: &'a [u32],
    offset: &'a [usize],
    cells: &'a [(usize, usize)],
    n: u32,
    labels: Vec<u32>,
    weight: Linears,
    content: Vec<u32>,
    totals: HashMap<Vec<u32>, AFrac<C>>,
}

impl<C: Ring> KsWalk<'_, C> {
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
                let e = self.weight.get_mut(&key).expect("just inserted");
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
/// Dispatches to [`jack_p_lb`]: the eigenoperator route wins the whole-degree
/// unit by a growing margin (8.6× at n = 10), which is why [`jack_table`]
/// calls it rather than the branching formula.
pub fn jack_p<C: Ring>(lambda: &Partition) -> Monomial<AFrac<C>> {
    jack_p_lb(lambda)
}

/// Multiply every coefficient by a product of linear-form powers, keeping the
/// scalar *factored* the whole way — the [`Frac::mul_factors`] pattern.
///
/// [`Frac::mul_factors`]: crate::frac::Frac::mul_factors
fn scale_by<C: Ring>(f: Monomial<AFrac<C>>, by: &Linears) -> Monomial<AFrac<C>> {
    let mut out = Monomial::zero();
    for (mu, c) in f.terms() {
        let mut v = c.mul_factors(by);
        v.reduce();
        out.add_term(mu.clone(), v);
    }
    out
}

/// `Q_λ = (H_λ/H'_λ)·P_λ`, the basis dual to `P` under `⟨·,·⟩_α`.
pub fn jack_q<C: Ring>(lambda: &Partition) -> Monomial<AFrac<C>> {
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
/// coefficients, and divisible by `u_μ = ∏ m_i(μ)!` ([KS] Thm 1.1). None of
/// that is arranged here: the coefficients arrive through fraction arithmetic
/// and the denominators cancel, which is why the tests that check it are real
/// checks on the whole route.
pub fn jack_j<C: Ring>(lambda: &Partition) -> Monomial<AFrac<C>> {
    scale_by(jack_p(lambda), &hook_lower(lambda))
}

/// Every `P_λ` of degree `n` — the unit of work the walls in
/// `docs/record/jack.md` are measured in, and the one Sage has no entry
/// point for.
pub fn jack_table<C: Ring>(n: u32) -> Vec<(Partition, Monomial<AFrac<C>>)> {
    crate::partitions_of(n)
        .into_iter()
        .map(|l| {
            let p = jack_p(&l);
            (l, p)
        })
        .collect()
}

/// Every `J_λ` of degree `n`.
pub fn jack_j_table<C: Ring>(n: u32) -> Vec<(Partition, Monomial<AFrac<C>>)> {
    crate::partitions_of(n)
        .into_iter()
        .map(|l| {
            let j = jack_j(&l);
            (l, j)
        })
        .collect()
}

/// `J_λ` in the **power-sum** basis — the Jack-character unit, and what the
/// [GJ] pipeline consumes.
///
/// Routed `m → s → p` through the [`convert`](crate::convert) hub, which is
/// generic over the coefficient ring. The `z_ν` divisions in
/// `PowerSum::from_schur` ask for [`QAlgebra`], and `AFrac<C>` is one for
/// *any* `C` — including `i128`, which is not.
pub fn jack_j_powersum<C: Ring>(lambda: &Partition) -> PowerSum<AFrac<C>> {
    PowerSum::<AFrac<C>>::from_schur(&jack_j::<C>(lambda).to_schur())
}

/// Every `J_λ` of degree `n`, in the power-sum basis.
pub fn jack_powersum_table<C: Ring>(n: u32) -> Vec<(Partition, PowerSum<AFrac<C>>)> {
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
pub fn powersum_scalar<C: Ring>(f: &PowerSum<AFrac<C>>, g: &PowerSum<AFrac<C>>) -> AFrac<C> {
    let mut out = <AFrac<C> as Ring>::zero();
    for (mu, a) in f.terms() {
        let Some(b) = g.terms().get(mu) else { continue };
        let mut alpha_pow = Linears::new();
        alpha_pow.insert((1, 0), mu.len() as i32); // α^{ℓ(μ)}
        let term = a.mul(b).mul(&AFrac::from_u128(mu.z()));
        out.add_assign(&term.mul_factors(&alpha_pow));
    }
    out.reduce();
    out
}

/// `⟨f, g⟩_α` for arbitrary monomial-basis elements, through the power sums.
pub fn jack_scalar<C: Ring>(f: &Monomial<AFrac<C>>, g: &Monomial<AFrac<C>>) -> AFrac<C> {
    let fp = PowerSum::<AFrac<C>>::from_schur(&f.to_schur());
    let gp = PowerSum::<AFrac<C>>::from_schur(&g.to_schur());
    powersum_scalar(&fp, &gp)
}

/// `⟨J_λ J_μ, J_ν⟩_α` — **Stanley's object**, whose membership in ℕ[α] is his
/// 1989 conjecture and still open.
///
/// Computed in the power-sum basis, where the product is a multiset union and
/// the pairing is diagonal, so no basis change of the product is ever formed.
/// Sage cannot compute `J[3,2,1]²` at all inside 120 s
/// (`docs/record/jack.md`), which is the first size the conjecture is
/// interesting at.
///
/// A negative coefficient here is a **result to report, not a bug to fix** —
/// the `∇e_n`-positivity and valley-Delta posture, verbatim.
///
/// ⚠️ **For a whole table, use [`stanley_table`].** This recomputes all three
/// p-expansions on every call, and sampling the degree-12 table put **94.6% of
/// the runtime inside [`jack_j_powersum`]** — 27951 conversions for 99 distinct
/// values. The single-shot form is the honest primitive and is kept as one;
/// the batch form is 17× faster and is what a search driver wants.
pub fn jack_structure_constant<C: Ring>(
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
/// handed out to other widths would launder exactly that away — the
/// `memo::bold_p` hazard, which is documented there for the same reason.
/// Hoisting the loop is the version with no correctness question in it.
///
/// Positivity is Stanley's 1989 conjecture and is **open**: this returns the
/// values, and asserts nothing about them.
pub fn stanley_table<C: Ring>(k: u32) -> Vec<(Partition, Partition, Partition, AFrac<C>)> {
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
pub fn omega_alpha<C: Ring>(f: &PowerSum<AFrac<C>>) -> PowerSum<AFrac<C>> {
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
pub fn specialize<C: Field>(f: &Monomial<AFrac<C>>, alpha: &C) -> Option<Monomial<C>> {
    let mut out = Monomial::zero();
    for (mu, c) in f.terms() {
        out.add_term(mu.clone(), c.eval(alpha)?);
    }
    Some(out)
}

/// The zonal polynomial `Z_λ` in [GJ]'s normalization: `J_λ` at α = 2.
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

    /// ⚠️ `h_low`, `h_up` and the [KS] weight `d` are three *different* linear
    /// families, and any two of them agree on enough small cells to pass a
    /// careless test. λ = (2,1) separates all three.
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
    /// tautology.
    #[test]
    fn the_two_engines_agree() {
        for n in 1..=7u32 {
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
        for n in 1..=5u32 {
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
        for n in 1..=8u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = jack_p(&lambda);
                assert_eq!(p.coeff(&lambda), <F as Ring>::one(), "monic at {lambda}");
                for (mu, _) in p.terms() {
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
        for n in 1..=5u32 {
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
        for n in 1..=5u32 {
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

    /// **[KS] Thm 1.1, as a law.** Every `[m_μ]J_λ` clears its denominator,
    /// lies in ℕ[α], **and** is divisible by `u_μ = ∏ m_i(μ)!`.
    ///
    /// The coefficients arrive through fraction arithmetic over `AFrac<i128>`,
    /// so every one of these is a check on the whole route rather than on a
    /// formula that already says so. Run over `i128` on purpose: integrality is
    /// then a property of the representation and not of a ℚ that could absorb
    /// anything.
    #[test]
    fn knop_sahi_positivity_and_u_divisibility() {
        for n in 1..=8u32 {
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
        for n in 1..=6u32 {
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

    /// **α = ∞ is the monomial basis, α = 0 is `e_{λ'}`** ([KS] p. 1). Taken as
    /// limits of `P`, not values: the leading behaviour in α of each
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
        for n in 1..=5u32 {
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

    /// Sage's `zonal()` is `P^{(2)}` and [GJ]'s `Z_λ` is `J^{(2)}`; the two
    /// differ by `H_λ(2)`. ⚠️ Recorded as a test because a fixture that gets
    /// this backwards still looks plausible.
    #[test]
    fn the_two_zonal_normalizations_differ_by_the_hook_product() {
        for n in 1..=5u32 {
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
        for n in 1..=6u32 {
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

    /// The α-twist really is load-bearing: plain `ω` breaks the duality at the
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
        for n in 1..=5u32 {
            for lambda in crate::partitions_of(n) {
                let j: Monomial<F> = jack_j(&lambda);
                let p: PowerSum<F> = jack_j_powersum(&lambda);
                let back: Monomial<F> = Monomial::from_schur(&p.to_schur());
                assert_eq!(back, j, "J_{lambda} through the p basis");
            }
        }
    }
}
