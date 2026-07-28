//! The basis-change engine: conversions between all six classical bases,
//! routed through the Schur hub.
//!
//! Each basis implements [`ToSchur`] (expand into the Schur basis) and
//! [`FromSchur`] (contract a Schur element into this basis); [`convert`] then
//! composes them to reach any ordered pair. The per-basis algorithms are the
//! classical ones:
//!
//! | conversion            | method                                  | ring    |
//! |-----------------------|-----------------------------------------|---------|
//! | h → s, e → s          | products of one-row / one-column Schurs | ℤ       |
//! | s → h                 | Jacobi–Trudi, by signed permutations    | ℤ       |
//! | s → e                 | dual Jacobi–Trudi, likewise             | ℤ       |
//! | p → s                 | iterated Murnaghan–Nakayama             | ℤ       |
//! | s → p                 | s_λ = Σ z_μ⁻¹ χ^λ(μ) p_μ                 | **ℚ**   |
//! | s → m                 | Kostka numbers                          | ℤ       |
//! | m → s                 | Muir's rule                             | ℤ       |
//! | f ↔ s                 | the m conversion, composed with ω        | ℤ       |
//!
//! Only s → p divides, and what it needs is a [`QAlgebra`] — a ring containing
//! ℚ — not a [`Field`](crate::coeff::Field). The division is by z_μ, an
//! *integer*, so ℚ[t] and ℚ[q,t] qualify even though neither is a field. Every
//! other path stays exact over ℤ.
//!
//! Three of these were rewritten after a degree ladder against Sage
//! (`scripts/compare_sage.py`) showed them *scaling* badly rather than merely
//! being slow. That distinction is the reason the ladder exists: at a single
//! size each looked like an acceptable constant factor, and s → e and m → s were
//! both **faster than Sage at degree 8** while losing badly by degree 20.

use std::collections::{BTreeMap, HashMap};

use crate::character::character_in;
use crate::coeff::{QAlgebra, Ring};
use crate::fasthash::Map;
use crate::kostka::kostka;
use crate::memo::{inverse_kostka_row_cached, lex_parts_cached, partitions_cached};
use crate::partition::Partition;
use crate::sym::{Elementary, Forgotten, Homogeneous, Monomial, PowerSum, Schur, SymAlgebra, SymFn};

/// Expand `self` into the Schur basis.
pub trait ToSchur<C: Ring> {
    fn to_schur(&self) -> Schur<C>;
}

/// Contract a Schur element into `Self`'s basis.
pub trait FromSchur<C: Ring>: Sized {
    fn from_schur(s: &Schur<C>) -> Self;
}

/// Convert between any two bases by composing through Schur.
pub fn convert<C, A, B>(a: &A) -> B
where
    C: Ring,
    A: ToSchur<C>,
    B: FromSchur<C>,
{
    B::from_schur(&a.to_schur())
}

// --- Schur: the hub, identity both ways -------------------------------------

impl<C: Ring> ToSchur<C> for Schur<C> {
    fn to_schur(&self) -> Schur<C> {
        self.clone()
    }
}
impl<C: Ring> FromSchur<C> for Schur<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        s.clone()
    }
}

// --- the Jacobi–Trudi determinants ------------------------------------------

/// The determinant det(y_{c_i − i + j})_{i,j} expanded over a **free
/// multiplicative** basis, as a signed list of partitions.
///
/// Both Jacobi–Trudi determinants have this shape, and both are taken over a
/// basis (h or e) whose elements are indexed by partitions and multiply by
/// *concatenating* indices: y_a · y_b = y_{a ∪ b}. So a single permutation term
///
/// ```text
///   sgn(w) · ∏_i y_{c_i − i + w(i)}
/// ```
///
/// is one signed **monomial** y_μ, with μ the sorted multiset {c_i − i + w(i)}.
/// The whole determinant is therefore a signed count of permutations grouped by
/// that multiset — no polynomial arithmetic anywhere.
///
/// That is the entire fix. This used to be a generic Laplace expansion over the
/// symmetric-function algebra, which cloned an (n−1)×(n−1) matrix of
/// *polynomials* at every node and did a full polynomial add and multiply per
/// term. Its cost was factorial in the matrix size — and the matrix size is
/// ℓ(λ) for s → h but **λ₁** for s → e, since that one is built from the
/// conjugate. Same helper, opposite behaviour: s → h stayed fast on the wide-
/// but-shallow shapes a degree ladder produces while s → e crossed over and
/// fell behind Sage past degree 16.
///
/// Enumerating permutations directly also prunes where the determinant is
/// sparse: a negative index means the entry is zero, so that whole subtree is
/// skipped rather than multiplied out.
///
/// **Rows are assigned from the last to the first, and that is the whole
/// performance story.** Entry (i, j) vanishes when j < i - c_i, and c is weakly
/// decreasing, so i - c_i *increases* with i: the constraint tightens as the row
/// index grows. Walking rows forward therefore starts at the least constrained
/// row - row 0 accepts any column - and only meets the dead ends near the
/// leaves, long after the branching has happened. Walking them backwards puts
/// the tightest row first, so whole subtrees die at depth 1.
///
/// The difference is not a constant factor. For lambda = (14) this change alone
/// took s -> e from 1.5 seconds to microseconds. Symmetrica's `tsh_jt` builds
/// the transposed matrix (lambda_i + i - j, invalid when j > lambda_i + i),
/// which puts its tight row first for free: the same determinant and the same
/// algorithm, with the orientation accounting for the entire gap.
fn jt_terms(c: &[u32]) -> Vec<(Partition, i64)> {
    if c.is_empty() {
        return vec![(Partition::default(), 1)];
    }
    let mut acc: HashMap<Partition, i64> = HashMap::new();
    let mut used = vec![false; c.len()];
    let mut idx: Vec<u32> = Vec::with_capacity(c.len());
    jt_rec(c, c.len(), &mut used, &mut idx, 1, &mut acc);
    acc.into_iter().filter(|(_, s)| *s != 0).collect()
}

/// Assign row `i - 1`, rows `i..` being already placed.
fn jt_rec(
    c: &[u32],
    i: usize,
    used: &mut [bool],
    idx: &mut Vec<u32>,
    sign: i64,
    acc: &mut HashMap<Partition, i64>,
) {
    if i == 0 {
        // `Partition::new` drops the zero indices (y_0 is the unit) and sorts.
        *acc.entry(Partition::new(idx.iter().copied())).or_insert(0) += sign;
        return;
    }
    let row = i - 1;
    // Rows after this one are already placed, so choosing w(row) = j inverts
    // against every *used* column beneath j - the mirror of the forward walk,
    // which counted the unused ones.
    let mut below = 0i64;
    // Every column below this is a zero entry, so the scan can start there.
    let lo = (row as i64 - c[row] as i64).max(0) as usize;
    for &u in used.iter().take(lo) {
        if u {
            below += 1;
        }
    }
    for j in lo..c.len() {
        if used[j] {
            below += 1;
            continue;
        }
        let k = c[row] as i64 - row as i64 + j as i64;
        used[j] = true;
        idx.push(k as u32);
        let s = if below % 2 == 0 { sign } else { -sign };
        jt_rec(c, row, used, idx, s, acc);
        idx.pop();
        used[j] = false;
    }
}

/// Assemble a signed partition list into a basis element.
fn from_jt<C: Ring, S: SymAlgebra<C>>(terms: Vec<(Partition, i64)>) -> S {
    let mut out = S::zero();
    for (mu, s) in terms {
        out.add_term(mu, C::from_i64(s));
    }
    out
}

// --- Homogeneous <-> Schur --------------------------------------------------

impl<C: Ring> ToSchur<C> for Homogeneous<C> {
    fn to_schur(&self) -> Schur<C> {
        // h_λ = ∏_i h_{λ_i} = ∏_i s_{(λ_i)}.
        let mut out = Schur::zero();
        for (lambda, c) in self.terms() {
            let mut prod = Schur::unit();
            for &part in lambda.parts() {
                prod = prod.mul(&Schur::monomial(Partition::new([part]), C::one()));
            }
            out = out.add(&prod.scale(c));
        }
        out
    }
}

impl<C: Ring> FromSchur<C> for Homogeneous<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = det(h_{λ_i − i + j}); the matrix is ℓ(λ)×ℓ(λ).
        contract_multiplicative(s, false)
    }
}

// --- Elementary <-> Schur ---------------------------------------------------

impl<C: Ring> ToSchur<C> for Elementary<C> {
    fn to_schur(&self) -> Schur<C> {
        // e_λ = ∏_i e_{λ_i} = ∏_i s_{(1^{λ_i})} (single columns).
        let mut out = Schur::zero();
        for (lambda, c) in self.terms() {
            let mut prod = Schur::unit();
            for &part in lambda.parts() {
                let column = Partition::new(std::iter::repeat(1).take(part as usize));
                prod = prod.mul(&Schur::monomial(column, C::one()));
            }
            out = out.add(&prod.scale(c));
        }
        out
    }
}

impl<C: Ring> FromSchur<C> for Elementary<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = det(e_{λ'_i − i + j}); the matrix is λ₁×λ₁, being built from the
        // conjugate.
        contract_multiplicative(s, true)
    }
}

/// The h ↔ e transition, one generator at a time.
///
/// Newton's identity for the elementary/complete pair, Σ_{k=0}^{n} (−1)^k e_k
/// h_{n−k} = 0, rearranged:
///
/// ```text
///   h_n = Σ_{k=1}^{n} (−1)^{k−1} e_k · h_{n−k}
/// ```
///
/// and *symmetrically* with h and e exchanged — the identity is invariant under
/// the swap, so a single routine serves both directions. Each step multiplies by
/// one generator, which in a multiplicative basis is a multiset union, so this
/// is a linear recursion over cheap products rather than anything determinantal.
///
/// `gen[n]` is then the one-part generator of the *source* basis expanded in the
/// target, and a multi-part index is the product of those.
fn flip_generators<C: Ring, S: SymAlgebra<C>>(upto: u32) -> Vec<S> {
    // The table is the *same* for h→e and e→h, and its coefficients are
    // integers independent of `C`, so it is computed once in i64 and injected.
    // Recomputing it per call was the dominant cost of a flipped conversion:
    // the recursion touches O(n²) products over elements with up to p(n) terms,
    // which is cheap once and wasteful on every call.
    let table = flip_table(upto);
    table
        .iter()
        .map(|terms| {
            let mut x = S::zero();
            for (mu, c) in terms {
                x.add_term(mu.clone(), C::from_i64(*c));
            }
            x
        })
        .collect()
}

thread_local! {
    /// h_n in the e-basis (equivalently e_n in the h-basis) for n = 0.., as
    /// integer term lists. Grown monotonically, never invalidated: these are
    /// fixed integers, so a longer table subsumes a shorter one.
    static FLIP: std::cell::RefCell<Vec<Vec<(Partition, i64)>>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

fn flip_table(upto: u32) -> Vec<Vec<(Partition, i64)>> {
    FLIP.with(|cell| {
        let mut t = cell.borrow_mut();
        if t.is_empty() {
            t.push(vec![(Partition::default(), 1)]);
        }
        while t.len() <= upto as usize {
            let n = t.len();
            let mut acc: HashMap<Partition, i64> = HashMap::new();
            for k in 1..=n {
                let sign = if k % 2 == 1 { 1i64 } else { -1 };
                for (mu, c) in &t[n - k] {
                    let mut parts = mu.parts().to_vec();
                    parts.push(k as u32);
                    *acc.entry(Partition::new(parts)).or_insert(0) += sign * c;
                }
            }
            t.push(acc.into_iter().filter(|(_, c)| *c != 0).collect());
        }
        t[..=upto as usize].to_vec()
    })
}

/// The other of the two multiplicative bases a Schur element can contract into.
///
/// An associated type rather than a flag, so "compute in whichever Jacobi–Trudi
/// matrix is smaller and flip back" is expressible without either basis naming
/// the other at the call site — and so the recursion that does it is checked to
/// terminate by the compiler rather than by argument.
pub trait Dual<C: Ring>: SymAlgebra<C> {
    type Other: SymAlgebra<C> + Dual<C, Other = Self>;
}

impl<C: Ring> Dual<C> for Homogeneous<C> {
    type Other = Elementary<C>;
}

impl<C: Ring> Dual<C> for Elementary<C> {
    type Other = Homogeneous<C>;
}

/// Re-express an element of one multiplicative basis in the other, h ↔ e.
///
/// This is what makes a *wide* shape cheap for s → e and a *tall* one cheap for
/// s → h: whichever Jacobi–Trudi matrix is smaller can be used, and the result
/// flipped into the basis actually wanted. ℓ(λ) and λ₁ are the two matrix sizes
/// and they trade off against each other, so taking the minimum turns each
/// direction's worst family into the other's best.
fn flip_basis<C: Ring, A: SymFn<C>, B: SymAlgebra<C>>(x: &A) -> B {
    let upto = x.terms().keys().map(|p| p.part(0)).max().unwrap_or(0);
    let gens: Vec<B> = flip_generators::<C, B>(upto);
    let mut out = B::zero();
    for (mu, c) in x.terms() {
        let mut prod = B::unit();
        for &part in mu.parts() {
            prod = prod.times(&gens[part as usize]);
        }
        out = out.add(&prod.scale(c));
    }
    out
}

/// Jacobi–Trudi matrix size past which the determinant is abandoned for a Muir
/// sweep.
///
/// [`jt_terms`] enumerates permutations, so it is exponential in the matrix
/// size — and the matrix is ℓ(λ) for s → h but **λ₁** for s → e. Both
/// directions therefore have a family of shapes that ruins them, and they are
/// each other's mirror image: s → e dies on a single long row, s → h on a
/// single tall column.
///
/// Measured against Symmetrica at degree 14, s → e was 0.5–3x *faster* out to
/// λ₁ = 8 and then 2.5x, 42x, 148x slower at λ₁ = 10, 12, 13 — 1.5 seconds for
/// s_{(14)} against its 0.02.
///
/// **The degree ladder never saw it.** `scripts/compare_sage.py` builds its
/// shapes with several rows, so λ₁ stayed under 8 and this direction looked
/// healthy at every degree tested, while the module doc above had already
/// named wide shapes as the hazard. It surfaced the moment Sage was allowed to
/// pick the inputs (`scripts/check_backend.py`).
///
/// Taking the smaller of the two matrices (see [`flip_basis`]) removes most of
/// the problem on its own, since a shape that is bad for one direction is good
/// for the other. What survives is the family bad for *both* — hooks, where
/// ℓ(λ) and λ₁ are each about |λ|/2.
///
/// The limit is set high because **the sweep is nearly always the worse of the
/// two**, which was not the expectation. For the hook of degree 20 the
/// determinant takes 0.0048s against the sweep's 0.37s, and at degree 24,
/// 0.070s against 6.05s: p(n) Muir expansions cost more than a 10-to-12 wide
/// determinant, and p(n) grows steadily while the determinant is only bad once
/// the matrix is genuinely large. So the sweep is a backstop against a
/// pathological shape, not a fast path, and the threshold sits where the
/// determinant finally stops being sub-second.
const JT_LIMIT: usize = 14;

/// Matrix size at which flipping to the conjugate basis starts to pay.
///
/// Flipping is *not* free — see [`flip_basis`] — and taking the smaller matrix
/// whenever it is smaller at all made things worse in aggregate: over every
/// partition of degree 20 it cost 0.33x against Symmetrica, where never
/// flipping gave 2.20x. A shape like (5,5,5,5) has matrices of 4 and 5, and a
/// 5-wide determinant is far cheaper than expanding an h-element of degree 20
/// into e.
///
/// So the flip is reserved for the cases where the determinant is genuinely
/// exponential and the conjugate collapses it: a single row of degree 24 is
/// 1.66s direct and 0.001s flipped, because its h-expansion is the one term
/// h_24. Below this size the direct determinant always wins.
const FLIP_MIN: usize = 14;

/// s → h and s → e, which are the same computation on λ and on λ'.
///
/// Two routes, chosen per term by the matrix size:
///
/// * **the determinant**, for small matrices, via [`jt_terms`];
/// * **a Muir sweep**, for large ones, using
///
///   ```text
///     coefficient of h_μ in s_λ  =  coefficient of s_λ in m_μ
///   ```
///
///   which is just K⁻¹ read by column instead of by row: h_μ = Σ_λ K_{λμ} s_λ
///   gives s = (Kᵀ)⁻¹h, so the h-coefficients of s_λ are the λ-column of K⁻¹ —
///   and [`muir_expand`] already produces its rows, fast. Sweeping every μ ⊢ n
///   and keeping the λ entry costs p(n) expansions regardless of shape, which
///   is the same "a table is p(n) sweeps" trade as
///   [`kostka_table`](crate::kostka::kostka_table).
///
/// e is the conjugate case throughout: ω is an isometry with ω(h_μ) = e_μ, so
/// the e-coefficients of s_λ are the h-coefficients of s_{λ'}.
fn contract_multiplicative<C: Ring, S: Dual<C>>(s: &Schur<C>, dual: bool) -> S {
    let mut out = S::zero();
    // Terms needing the sweep, grouped by degree since a sweep serves a whole
    // degree at once.
    let mut swept: BTreeMap<u32, Vec<(Partition, C)>> = BTreeMap::new();
    // Terms cheaper in the *other* basis, to be flipped back at the end.
    let mut crossed: Schur<C> = Schur::zero();
    for (lambda, c) in s.terms() {
        let index = if dual { lambda.conjugate() } else { lambda.clone() };
        let other = index.conjugate();
        if index.len().min(other.len()) > JT_LIMIT {
            swept.entry(index.size()).or_default().push((index, c.clone()));
        } else if index.len() >= FLIP_MIN && other.len() < index.len() {
            // The conjugate determinant is smaller: compute there and flip.
            crossed.add_term(lambda.clone(), c.clone());
        } else {
            out = out.add(&from_jt::<C, S>(jt_terms(index.parts())).scale(c));
        }
    }
    if !crossed.is_zero() {
        // `dual` is inverted, so this lands in the opposite basis; `flip_basis`
        // brings it back. Both are multiplicative, so the flip is the Newton
        // recursion above and never a determinant.
        let mirror: S::Other = contract_multiplicative(&crossed, !dual);
        out = out.add(&flip_basis::<C, _, S>(&mirror));
    }
    for (n, targets) in swept {
        match muir_column(&targets, n) {
            Some(terms) => {
                for (mu, v) in terms {
                    out.add_term(mu, v);
                }
            }
            // Past the β-mask width, where `muir_expand` declines. The
            // determinant has no ceiling, only a cost.
            None => {
                for (index, c) in &targets {
                    out = out.add(&from_jt::<C, S>(jt_terms(index.parts())).scale(c));
                }
            }
        }
    }
    out
}

/// Σ_λ c_λ · (column λ of K⁻¹), by expanding every m_μ of degree `n` once.
fn muir_column<C: Ring>(targets: &[(Partition, C)], n: u32) -> Option<Vec<(Partition, C)>> {
    let want: HashMap<&Partition, &C> = targets.iter().map(|(p, c)| (p, c)).collect();
    let mut out = Vec::new();
    for mu in partitions_cached(n).iter() {
        let mut acc = C::zero();
        for (lambda, v) in muir_expand(mu)? {
            if let Some(c) = want.get(&lambda) {
                acc.add_assign(&C::from_i128(v).mul(c));
            }
        }
        if !acc.is_zero() {
            out.push((mu.clone(), acc));
        }
    }
    Some(out)
}

// --- PowerSum <-> Schur -----------------------------------------------------

impl<C: Ring> ToSchur<C> for PowerSum<C> {
    fn to_schur(&self) -> Schur<C> {
        let mut out = Schur::zero();
        // Terms are batched by degree, because p_μ and p_ν share Murnaghan–
        // Nakayama work whenever they share parts — and expanding each one
        // separately from ∅ throws all of that away. The β-mask width is |μ|,
        // so only equal degrees can share a sweep.
        let mut by_degree: HashMap<u32, Vec<(&Partition, &C)>> = HashMap::new();
        for (mu, c) in self.terms() {
            let n = mu.size();
            if n as usize <= MASK_LIMIT {
                by_degree.entry(n).or_default().push((mu, c));
            } else {
                // Degree past the β-mask width: fall back to characters, which
                // are exact in `C` and so stay correct for bignum rings.
                for lambda in partitions_cached(n).iter() {
                    let chi = character_in::<C>(lambda, mu);
                    if !chi.is_zero() {
                        out.add_term(lambda.clone(), chi.mul(c));
                    }
                }
            }
        }
        let mut degrees: Vec<u32> = by_degree.keys().copied().collect();
        degrees.sort_unstable();
        for n in degrees {
            let mut items = by_degree.remove(&n).unwrap();
            // Descending part order, so the longest common prefixes are shared.
            items.sort_by(|a, b| a.0.parts().cmp(b.0.parts()));
            let l = n as usize;
            let mut root: Map<u64, i128> = Map::default();
            if l == 0 {
                for (_, c) in &items {
                    out.add_term(Partition::default(), (*c).clone());
                }
                continue;
            }
            root.insert((1u64 << l) - 1, 1);
            if integral_sweep(&items, l, &root, &mut out) {
                continue;
            }
            // Accumulate on the β-mask, not on a Partition. Every leaf of the
            // traversal touches the whole frontier, so keying by partition
            // allocated and sorted a fresh Vec — and hashed a heap key — once
            // per (μ, mask) pair to produce a few dozen distinct terms. Masks
            // are u64, and the partitions get built once at the end.
            let mut acc: Map<u64, C> = Map::default();
            p_expand_shared(&items, 0, &root, &mut |c: &&C, mask, chi| {
                let term = C::from_i128(chi).mul(c);
                acc.entry(mask).or_insert_with(C::zero).add_assign(&term);
            });
            for (mask, c) in acc {
                if !c.is_zero() {
                    out.add_term(mask_to_partition(mask, l), c);
                }
            }
        }
        out
    }
}

/// Expand a batch of p_μ of one degree in a single traversal, sharing the
/// Murnaghan–Nakayama sweep across every common prefix.
///
/// `items` is sorted by part sequence, so partitions agreeing in their first
/// `depth` parts are contiguous; each such run continues from *one* frontier
/// instead of rebuilding it. Plethysm is the case that motivates this: it
/// finishes by converting a p-element of degree d·e with dozens of terms, and
/// that conversion was ~99% of its runtime.
pub(crate) fn p_expand_shared<T, F>(
    items: &[(&Partition, T)],
    depth: usize,
    frontier: &Map<u64, i128>,
    emit: &mut F,
) where
    F: FnMut(&T, u64, i128),
{
    let mut i = 0;
    // Partitions that end here: emit the frontier against their coefficient.
    while i < items.len() && items[i].0.len() == depth {
        for (&mask, &chi) in frontier {
            if chi != 0 {
                emit(&items[i].1, mask, chi);
            }
        }
        i += 1;
    }
    // The rest are grouped by their next part, each group sharing one step.
    while i < items.len() {
        let k = items[i].0.part(depth);
        let start = i;
        while i < items.len() && items[i].0.part(depth) == k {
            i += 1;
        }
        let next = p_step(frontier, k);
        p_expand_shared(&items[start..i], depth + 1, &next, emit);
    }
}

/// Run one degree's sweep entirely in `i128`, over a common denominator.
///
/// Returns `false` if that is not possible — an unsupported ring, or an
/// overflow — in which case the caller takes the generic path and nothing has
/// been written to `out`.
///
/// The point is that `p → s` computes Σ_μ c_μ χ^λ(μ), a sum of
/// (coefficient × integer) terms. Done in ℚ that is a rational multiply and a
/// rational add per (μ, mask) leaf, each normalising by a gcd — and profiling
/// put **55%** of plethysm's runtime in those gcds and the i128 division
/// underneath them. Putting every c_μ over one denominator D makes the entire
/// accumulation integer, with a single conversion back per output term.
///
/// D is the lcm of the denominators, which measurement said is the right shape
/// for this: across the plethysms driving this work the lcm *equalled the
/// largest denominator* every time, and never exceeded ~5·10⁵. Overflow is
/// still checked at every step rather than argued away, because the guarantee
/// is only about the cases measured.
fn integral_sweep<C: Ring>(
    items: &[(&Partition, &C)],
    l: usize,
    root: &Map<u64, i128>,
    out: &mut Schur<C>,
) -> bool {
    let mut den: i128 = 1;
    let mut ratios = Vec::with_capacity(items.len());
    for (_, c) in items {
        let (n, d) = match c.as_ratio() {
            Some(r) => r,
            None => return false, // ring opts out
        };
        ratios.push((n, d));
        let g = gcd_i128(den, d);
        den = match den.checked_mul(d / g) {
            Some(v) => v,
            None => return false,
        };
    }
    // Numerators rescaled to the common denominator.
    let mut scaled = Vec::with_capacity(items.len());
    for &(n, d) in &ratios {
        match n.checked_mul(den / d) {
            Some(v) => scaled.push(v),
            None => return false,
        }
    }
    let paired: Vec<(&Partition, i128)> = items
        .iter()
        .zip(&scaled)
        .map(|((mu, _), &s)| (*mu, s))
        .collect();

    let mut acc: Map<u64, i128> = Map::default();
    let mut overflow = false;
    p_expand_shared(&paired, 0, root, &mut |&s: &i128, mask, chi| {
        if overflow {
            return;
        }
        let slot = acc.entry(mask).or_insert(0);
        match s.checked_mul(chi).and_then(|t| slot.checked_add(t)) {
            Some(v) => *slot = v,
            None => overflow = true,
        }
    });
    if overflow {
        return false;
    }
    // Convert back once per output term, not once per leaf.
    let mut built = Vec::with_capacity(acc.len());
    for (&mask, &v) in &acc {
        if v == 0 {
            continue;
        }
        match C::from_ratio(v, den) {
            Some(c) => built.push((mask, c)),
            None => return false,
        }
    }
    for (mask, c) in built {
        out.add_term(mask_to_partition(mask, l), c);
    }
    true
}

fn gcd_i128(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

/// p_μ in the Schur basis, by **iterated Murnaghan–Nakayama** rather than by
/// evaluating characters.
///
/// `p_μ = Σ_λ χ^λ(μ) s_λ`, and the obvious implementation asks for χ^λ(μ) once
/// per λ — p(n) independent recursions for one p_μ. But MN is itself a
/// multiplication rule,
///
/// ```text
///   p_k · s_λ = Σ (−1)^{ht(ξ)} s_{λ ∪ ξ},   ξ a k-rim-hook added to λ,
/// ```
///
/// so multiplying successively by p_{μ₁}, p_{μ₂}, … builds the entire expansion
/// in ℓ(μ) passes and never computes a character at all. Same "produce the whole
/// answer in one sweep rather than query it entry by entry" shape as the Kostka
/// and coproduct fixes.
///
/// Rim hooks are handled in β-numbers (first-column hook lengths, strictly
/// decreasing): adding a k-rim-hook is replacing some β by β+k when that value
/// is free, and the height is how many β lie strictly between them. Using a
/// fixed β-length of n keeps every intermediate comparable — a partition of n
/// has at most n rows, so no representable shape is lost.
///
/// Accumulation is in **i128**, not the caller's ring, and that is a deliberate
/// and safe narrowing rather than the usual truncation hazard. Every value here
/// is a character χ^λ(μ), and |χ^λ(μ)| ≤ f^λ ≤ √(n!). This function already
/// declines when l > 32 (the β-mask does not fit), and √(32!) ≈ 1.6·10¹⁸ — so on
/// every input it *accepts*, i128 cannot overflow. Callers past that width take
/// the character fallback, which stays exact in `C` for bignum rings.
///
/// The narrowing matters because the ring is the hot loop. Plethysm runs this
/// over ℚ, where each rim hook cost a rational add — a gcd — plus a temporary
/// from negating the coefficient, to combine two integers. That made p → s
/// ~99% of plethysm's runtime.
///
/// Retained as the **reference form**: production now takes the batched
/// [`p_expand_shared`] path, which shares this sweep across every p_μ of a
/// degree, and `batched_p_expansion_matches_per_term` holds the two to
/// agreement. The prefix grouping is the part of that batching that can quietly
/// go wrong, so it gets an oracle rather than trust.
#[cfg(test)]
fn p_expand(mu: &Partition) -> Option<Vec<(Partition, i128)>> {
    let l = mu.size() as usize;
    if l == 0 {
        return Some(vec![(Partition::default(), 1)]);
    }
    if l > MASK_LIMIT {
        return None;
    }
    // β-numbers of ∅ with l slots: {0, 1, …, l−1}.
    let mut cur: Map<u64, i128> = Map::default();
    cur.insert((1u64 << l) - 1, 1);
    for &k in mu.parts() {
        cur = p_step(&cur, k);
    }
    Some(
        cur.into_iter()
            .filter(|&(_, c)| c != 0)
            .map(|(mask, c)| (mask_to_partition(mask, l), c))
            .collect(),
    )
}

/// β values run from 0 to at most (l−1) + max part < 2l, so a 64-bit mask holds
/// the whole set for l ≤ 32. Past that the mask — which is what makes any of
/// this worth doing — no longer fits, and callers fall back.
pub(crate) const MASK_LIMIT: usize = 32;

/// One Murnaghan–Nakayama step: multiply a frontier of β-masks by p_k.
pub(crate) fn p_step(cur: &Map<u64, i128>, k: u32) -> Map<u64, i128> {
    // The frontier grows monotonically through a sweep, so a default-capacity
    // map rehashes several times per step. Sizing to the input is a floor on
    // the output, not a guess.
    let mut next: Map<u64, i128> =
        Map::with_capacity_and_hasher(cur.len() * 2, Default::default());
    for (&mask, &c) in cur {
        let mut rest = mask;
        while rest != 0 {
            let b = rest.trailing_zeros();
            rest &= rest - 1;
            let nb = b + k;
            if mask >> nb & 1 == 1 {
                continue; // that β is taken: no such rim hook
            }
            // Height = how many β lie strictly between b and b+k.
            let between = mask & (((1u64 << nb) - 1) ^ ((1u64 << (b + 1)) - 1));
            let m = (mask & !(1u64 << b)) | (1u64 << nb);
            let slot = next.entry(m).or_insert(0);
            if between.count_ones() % 2 == 0 {
                *slot += c;
            } else {
                *slot -= c;
            }
        }
    }
    next
}

/// Bits high-to-low are β₀ > β₁ > …, and λ_i = β_i − (l−1−i).
fn mask_to_partition(mask: u64, l: usize) -> Partition {
    let mut parts = Vec::with_capacity(l);
    let mut rest = mask;
    let mut i = 0usize;
    while rest != 0 {
        let b = 63 - rest.leading_zeros() as usize;
        rest &= !(1u64 << b);
        let part = b - (l - 1 - i);
        if part > 0 {
            parts.push(part as u32);
        }
        i += 1;
    }
    Partition::new(parts)
}

impl<C: QAlgebra> FromSchur<C> for PowerSum<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = Σ_μ z_μ⁻¹ χ^λ(μ) p_μ  (needs division by z_μ).
        let mut out = PowerSum::zero();
        for (lambda, c) in s.terms() {
            for mu in partitions_cached(lambda.size()).iter() {
                let chi = character_in::<C>(lambda, mu);
                if !chi.is_zero() {
                    // One division by an integer — never by a ring element.
                    // That is exactly the `QAlgebra` contract, and why this is
                    // not bounded on `Field`.
                    out.add_term(mu.clone(), c.mul(&chi).div_u128(mu.z()));
                }
            }
        }
        out
    }
}

// --- Monomial <-> Schur -----------------------------------------------------

impl<C: Ring> FromSchur<C> for Monomial<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = Σ_μ K_{λμ} m_μ.
        let mut out = Monomial::zero();
        for (lambda, c) in s.terms() {
            for mu in partitions_cached(lambda.size()).iter() {
                let k = kostka(lambda, mu);
                if k != 0 {
                    out.add_term(mu.clone(), C::from_u128(k).mul(c));
                }
            }
        }
        out
    }
}

impl<C: Ring> ToSchur<C> for Monomial<C> {
    fn to_schur(&self) -> Schur<C> {
        let mut out = Schur::zero();
        for (mu, c) in self.terms() {
            match muir_expand(mu) {
                Some(terms) => {
                    for (lambda, v) in terms {
                        out.add_term(lambda, C::from_i128(v).mul(c));
                    }
                }
                // Past the β-mask width: fall back to the row solve, which is
                // slow but has no degree ceiling.
                None => {
                    let parts = lex_parts_cached(mu.size());
                    let row = inverse_kostka_row_cached(mu, || inverse_kostka_row(&parts, mu));
                    for (i, lambda) in parts.iter().enumerate() {
                        if row[i] != 0 {
                            // `row` is i128 because that is what the solve is
                            // built in; inject at that width, not through i64.
                            out.add_term(lambda.clone(), C::from_i128(row[i]).mul(c));
                        }
                    }
                }
            }
        }
        out
    }
}

// --- Forgotten <-> Schur ----------------------------------------------------
//
// f_λ = ω(m_λ), and ω is an involution, so both directions are the monomial
// conversion with an ω applied on the Schur side. Nothing new is computed:
//
//   to_schur:    Σ c_λ f_λ = ω(Σ c_λ m_λ)        →  ω(monomial_to_schur(c))
//   from_schur:  x = Σ c_λ f_λ  ⟺  ω(x) = Σ c_λ m_λ  →  monomial coeffs of ω(x)
//
// ω on the Schur basis is just conjugation of every index, so the extra cost is
// one transpose per term. Both directions therefore inherit Muir's rule and the
// Kostka machinery — including their asymptotics — for free.

impl<C: Ring> ToSchur<C> for Forgotten<C> {
    fn to_schur(&self) -> Schur<C> {
        Monomial::from_terms(self.terms().clone()).to_schur().omega()
    }
}

impl<C: Ring> FromSchur<C> for Forgotten<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        let m: Monomial<C> = Monomial::from_schur(&s.omega());
        Forgotten::from_terms(m.terms().clone())
    }
}

/// m_μ in the Schur basis, by **Muir's rule** — no Kostka numbers and no linear
/// solve.
///
/// In n variables the bialternant gives s_λ = a_{λ+δ}/a_δ, and multiplying by
/// m_μ = Σ_α x^α (over the distinct rearrangements α of μ) just shifts exponents:
///
/// ```text
///   m_μ · a_δ = Σ_α a_{α+δ},   so   m_μ = Σ_α ± s_{sort(α+δ) − δ}
/// ```
///
/// with `a_β = 0` when β repeats and `± a_{sorted β}` otherwise. Taking n = |μ|
/// slots loses nothing: (K⁻¹)_{μλ} ≠ 0 forces μ ⊵ λ, so ℓ(λ) ≤ |μ|.
///
/// Working in β-numbers β = α + δ, each part of μ is *added to a distinct slot*
/// of the initial set {0, 1, …, l−1}, which makes this the same machinery as
/// [`p_expand`] — a u64 mask, values moved one at a time, sign flipped by the
/// number of occupied values jumped over. The one difference from
/// Murnaghan–Nakayama is that MN may move the same value repeatedly while here
/// every slot receives at most one part.
///
/// **Slots are processed in decreasing initial value, and that is what makes it
/// fast.** At the moment slot v is handled, every slot above it is final and
/// every slot below still sits at its own value, which is < v < v+k. So a
/// collision at v+k can only be with an already-final value: it can never be
/// resolved later, and the branch dies immediately. Processing upward instead
/// would make the same test unsound, since the occupant might yet move. That
/// prune is the difference between exploring the C(l, ℓ(μ)) · (rearrangements)
/// placements and exploring only the surviving ones — the placements chain
/// together exactly as rim hooks do, so the survivors are few.
///
/// This replaces a forward solve of the unitriangular system K⁻¹K = I, which
/// needed O(p(n)²) Kostka numbers to read p(n) of them. That solve was the
/// single largest deficit in the library: 244× slower than Sage at degree 20 and
/// widening, because the improvements before it attacked the constant and left
/// the complexity alone.
fn muir_expand(mu: &Partition) -> Option<Vec<(Partition, i128)>> {
    let l = mu.size() as usize;
    if l == 0 {
        return Some(vec![(Partition::default(), 1)]);
    }
    // β values reach (l−1) + μ₁ ≤ 2l−1, so a 64-bit mask holds them for l ≤ 32.
    if l > 32 {
        return None;
    }
    // Parts grouped by value: choosing "a part equal to k" once is what makes
    // the rearrangements *distinct*, so equal parts are never double-counted.
    let mut avail: Vec<(u32, u32)> = Vec::new();
    for &k in mu.parts() {
        match avail.last_mut() {
            Some((v, n)) if *v == k => *n += 1,
            _ => avail.push((k, 1)),
        }
    }
    let mut acc: HashMap<Partition, i128> = HashMap::new();
    muir_rec(l, l as i32 - 1, (1u64 << l) - 1, &mut avail, mu.len(), 1, &mut acc);
    Some(acc.into_iter().filter(|(_, v)| *v != 0).collect())
}

fn muir_rec(
    l: usize,
    v: i32,
    mask: u64,
    avail: &mut [(u32, u32)],
    left: usize,
    sign: i128,
    acc: &mut HashMap<Partition, i128>,
) {
    if v < 0 {
        if left == 0 {
            // β sorted descending, then λ_i = β_i − (l − i).
            let mut lam = Vec::with_capacity(l);
            let mut i = 0usize;
            for b in (0..64u32).rev() {
                if mask >> b & 1 == 1 {
                    lam.push(b - (l - 1 - i) as u32);
                    i += 1;
                }
            }
            *acc.entry(Partition::new(lam)).or_insert(0) += sign;
        }
        return;
    }
    // Every remaining part needs its own slot, and v+1 slots remain.
    if left > v as usize + 1 {
        return;
    }
    // Leave slot v where it is. Safe without a collision test: every already
    // final value is > v.
    muir_rec(l, v - 1, mask, avail, left, sign, acc);

    if left == 0 {
        return;
    }
    for i in 0..avail.len() {
        if avail[i].1 == 0 {
            continue;
        }
        let nb = v as u32 + avail[i].0;
        if mask >> nb & 1 == 1 {
            continue; // permanent collision — see the note above
        }
        // Jumping an occupied value transposes the two, flipping the sign.
        let between = mask & (((1u64 << nb) - 1) ^ ((1u64 << (v + 1)) - 1));
        let s = if between.count_ones() % 2 == 0 { sign } else { -sign };
        avail[i].1 -= 1;
        muir_rec(l, v - 1, (mask & !(1 << v)) | (1 << nb), avail, left - 1, s, acc);
        avail[i].1 += 1;
    }
}

/// Row μ of the inverse Kostka matrix, indexed against `parts` (decreasing-lex).
///
/// In that order the Kostka matrix K is upper-unitriangular — `K[i][j] ≠ 0`
/// needs λᵢ ⊵ λⱼ, and dominance implies lex — so `K⁻¹K = I` restricted to row μ
/// solves forward:
///
/// ```text
///   w[j] = 1,   w[jj] = −Σ_{m=j}^{jj−1} w[m]·K[m][jj]   for jj > j
/// ```
///
/// Only this one row is ever needed: m_μ = Σ_λ (K⁻¹)_{μλ} s_λ. Building the
/// whole matrix and inverting it, as this used to, computed p(n)² Kostka numbers
/// to read p(n) of them — 2.0 s for a single degree-20 conversion, of which the
/// matrix was ~97%.
///
/// The `w[m] == 0` skip is the part that matters: a zero coefficient makes its
/// Kostka number irrelevant, so the call is never made rather than made and
/// multiplied by zero.
fn inverse_kostka_row(parts: &[Partition], mu: &Partition) -> Vec<i128> {
    let size = parts.len();
    let j = parts
        .iter()
        .position(|p| p == mu)
        .expect("μ ∈ partitions(|μ|)");
    let mut w = vec![0i128; size];
    w[j] = 1;
    for jj in (j + 1)..size {
        let mut acc = 0i128;
        for m in j..jj {
            if w[m] == 0 {
                continue;
            }
            acc += w[m] * kostka(&parts[m], &parts[jj]) as i128;
        }
        w[jj] = -acc;
    }
    w
}

// --- Monomial multiplication (routed through Schur) -------------------------

impl<C: Ring> Monomial<C> {
    /// Product in the monomial basis. Unlike p/e/h this is not multiplicative,
    /// so it is computed by expanding into Schur, multiplying, and contracting
    /// back — entirely over ℤ.
    pub fn mul(&self, other: &Self) -> Self {
        let prod = self.to_schur().mul(&other.to_schur());
        Monomial::from_schur(&prod)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    /// f_{(n)} = (−1)^{n−1} p_n and f_{(1^n)} = h_n.
    ///
    /// Both are hand-derivable and neither mentions ω, which is the point: the
    /// implementation *is* "apply ω", so a test phrased in terms of ω would only
    /// restate it. These come from the two edge cases of the monomial basis,
    /// m_{(n)} = p_n and m_{(1^n)} = e_n, pushed through ω(p_n) = (−1)^{n−1} p_n
    /// and ω(e_n) = h_n — facts about the *other* bases.
    #[test]
    fn forgotten_endpoints_match_hand_computation() {
        for n in 1..=6u32 {
            let row: Forgotten<i64> = Forgotten::monomial(part(&[n]), 1);
            let want = if (n - 1) % 2 == 0 { 1 } else { -1 };
            let got: PowerSum<Rational> = PowerSum::from_schur(&to_rat(&row.to_schur()));
            assert_eq!(got.terms().len(), 1, "f_({n}) should be a single p term");
            assert_eq!(got.coeff(&part(&[n])), Rational::from_int(want));

            let col: Forgotten<i64> = Forgotten::monomial(part(&vec![1; n as usize]), 1);
            let h: Homogeneous<i64> = Homogeneous::from_schur(&col.to_schur());
            assert_eq!(h.coeff(&part(&[n])), 1);
            assert_eq!(h.terms().len(), 1, "f_(1^{n}) should be exactly h_{n}");
        }
    }

    /// {f_λ} is dual to {e_λ} under the Hall inner product: ⟨f_λ, e_μ⟩ = δ_{λμ}.
    ///
    /// This is the structural characterisation of the forgotten basis, and it
    /// reaches it through code the conversion never touches — `hall` and the
    /// e → s expansion. If `Forgotten` were wired to the wrong involution, or to
    /// conjugation on the *index* rather than on the Schur expansion, the
    /// duality would fail while a round trip still closed.
    #[test]
    fn forgotten_is_dual_to_elementary() {
        for n in 1..=7u32 {
            let parts = partitions_cached(n);
            for lambda in parts.iter() {
                let f: Forgotten<i64> = Forgotten::monomial(lambda.clone(), 1);
                for mu in parts.iter() {
                    let e: Elementary<i64> = Elementary::monomial(mu.clone(), 1);
                    let want = i64::from(lambda == mu);
                    assert_eq!(
                        crate::ops::hall::<i64, _, _>(&f, &e),
                        want,
                        "<f_{lambda}, e_{mu}>"
                    );
                }
            }
        }
    }

    /// Every basis pair round-trips through the Schur hub, forgotten included.
    #[test]
    fn forgotten_round_trips_through_every_basis() {
        for n in 1..=7u32 {
            for lambda in partitions_cached(n).iter() {
                let f: Forgotten<i64> = Forgotten::monomial(lambda.clone(), 1);
                let s = f.to_schur();
                assert_eq!(Forgotten::from_schur(&s), f, "f -> s -> f");

                // out through each other basis and back
                let m: Monomial<i64> = Monomial::from_schur(&s);
                assert_eq!(Forgotten::from_schur(&m.to_schur()), f, "via m");
                let e: Elementary<i64> = Elementary::from_schur(&s);
                assert_eq!(Forgotten::from_schur(&e.to_schur()), f, "via e");
                let h: Homogeneous<i64> = Homogeneous::from_schur(&s);
                assert_eq!(Forgotten::from_schur(&h.to_schur()), f, "via h");
            }
        }
    }

    fn to_rat(s: &Schur<i64>) -> Schur<Rational> {
        let mut out = Schur::zero();
        for (p, c) in s.terms() {
            out.add_term(p.clone(), Rational::from_int(*c as i128));
        }
        out
    }

    #[test]
    fn h_and_e_expand_into_schur() {
        // h_{11} = s_2 + s_{11}; e_2 = s_{11}.
        let h11: Homogeneous<i64> = Homogeneous::monomial(part(&[1, 1]), 1);
        let s = h11.to_schur();
        assert_eq!(s.coeff(&part(&[2])), 1);
        assert_eq!(s.coeff(&part(&[1, 1])), 1);

        let e2: Elementary<i64> = Elementary::monomial(part(&[2]), 1);
        let se = e2.to_schur();
        assert_eq!(se.coeff(&part(&[1, 1])), 1);
        assert_eq!(se.terms().len(), 1);
    }

    #[test]
    fn jacobi_trudi_inverts_expansion() {
        // s_{11} = h_{11} − h_2 (Jacobi–Trudi).
        let s11: Schur<i64> = Schur::monomial(part(&[1, 1]), 1);
        let h = Homogeneous::from_schur(&s11);
        assert_eq!(h.coeff(&part(&[1, 1])), 1);
        assert_eq!(h.coeff(&part(&[2])), -1);
    }

    /// The batched p → s must equal expanding each p_μ on its own and summing.
    /// Batching groups terms by shared prefix and continues one frontier per
    /// group; an off-by-one in that grouping would silently attribute a term to
    /// the wrong μ, which no round-trip test would catch.
    #[test]
    fn batched_p_expansion_matches_per_term() {
        for n in 1..=10u32 {
            let mus = partitions_cached(n);
            // A single element carrying *every* p_μ of the degree at once, which
            // is the shape plethysm produces and the case batching exists for.
            let mut elt: PowerSum<i128> = PowerSum::zero();
            for (i, mu) in mus.iter().enumerate() {
                elt.add_term(mu.clone(), (i as i128) + 1);
            }
            let got: Schur<i128> = elt.to_schur();

            let mut want: Schur<i128> = Schur::zero();
            for (i, mu) in mus.iter().enumerate() {
                for (lambda, chi) in p_expand(mu).expect("small degree") {
                    want.add_term(lambda, chi * ((i as i128) + 1));
                }
            }
            assert_eq!(got, want, "batched p → s at degree {n}");
        }
    }

    /// The common-denominator integer sweep must agree exactly with rational
    /// arithmetic done per term.
    ///
    /// `batched_p_expansion_matches_per_term` above cannot cover this: it uses
    /// `PowerSum<i128>`, and `i128` declines `as_ratio`, so it exercises the
    /// generic path only. Over ℚ the integral path is what actually runs, and
    /// it is the one with a common denominator, a rescale, and three overflow
    /// checks in it.
    ///
    /// Coefficients are z_μ⁻¹, which is not an arbitrary choice — that is
    /// exactly what `s → p` produces and what plethysm then feeds back through
    /// `p → s`.
    #[test]
    fn integral_sweep_matches_rational_arithmetic() {
        for n in 1..=10u32 {
            let mus = partitions_cached(n);
            let mut elt: PowerSum<Rational> = PowerSum::zero();
            for mu in mus.iter() {
                elt.add_term(mu.clone(), Rational::new(1, mu.z() as i128));
            }
            let got: Schur<Rational> = elt.to_schur();

            let mut want: Schur<Rational> = Schur::zero();
            for mu in mus.iter() {
                let c = Rational::new(1, mu.z() as i128);
                for (lambda, chi) in p_expand(mu).expect("small degree") {
                    want.add_term(lambda, Rational::from_i128(chi).mul(&c));
                }
            }
            assert_eq!(got, want, "integral p → s at degree {n}");
        }
    }

    /// A ring that declines `as_ratio` must still convert correctly — the
    /// integral sweep has to bail without having written anything.
    #[test]
    fn rings_that_decline_as_ratio_fall_back_cleanly() {
        for n in 1..=8u32 {
            let mus = partitions_cached(n);
            let mut elt: PowerSum<i128> = PowerSum::zero();
            for (i, mu) in mus.iter().enumerate() {
                elt.add_term(mu.clone(), (i as i128) - 3);
            }
            assert!(i128::from_ratio(1, 2).is_none(), "i128 should decline");
            let got: Schur<i128> = elt.to_schur();
            let mut want: Schur<i128> = Schur::zero();
            for (i, mu) in mus.iter().enumerate() {
                for (lambda, chi) in p_expand(mu).expect("small degree") {
                    want.add_term(lambda, chi * ((i as i128) - 3));
                }
            }
            assert_eq!(got, want, "fallback p → s at degree {n}");
        }
    }

    /// Muir's rule against the linear solve it replaced, on **every** μ up to
    /// degree 12 — the two are independent routes to the same row of K⁻¹, and
    /// the solve was the shipped implementation, so this is a real oracle rather
    /// than a self-consistency check.
    #[test]
    fn muir_agrees_with_the_inverse_kostka_solve() {
        for n in 0..=12u32 {
            let parts = lex_parts_cached(n);
            for mu in parts.iter() {
                let row = inverse_kostka_row(&parts, mu);
                let mut want: Vec<(Partition, i128)> = parts
                    .iter()
                    .zip(&row)
                    .filter(|(_, &v)| v != 0)
                    .map(|(l, &v)| (l.clone(), v))
                    .collect();
                let mut got = muir_expand(mu).expect("degree 12 is inside the mask");
                want.sort();
                got.sort();
                assert_eq!(got, want, "m_{mu} in the Schur basis");
            }
        }
    }

    /// The two shapes worked by hand when the rule was derived. `m_{21}` is the
    /// one that catches a wrong triangularity assumption: (K⁻¹)_{μλ} is nonzero
    /// for μ ⊵ λ, so λ may be *longer* than μ, and taking only ℓ(μ) slots would
    /// silently drop the s_{111} term and return the 2-variable answer.
    #[test]
    fn muir_matches_hand_computation() {
        let got = muir_expand(&part(&[2, 1])).unwrap();
        let mut got: Vec<_> = got.iter().map(|(l, c)| (l.parts().to_vec(), *c)).collect();
        got.sort();
        assert_eq!(got, vec![(vec![1, 1, 1], -2), (vec![2, 1], 1)]);

        // m_{(n)} = Σ_i (−1)^i s_{(n−i, 1^i)}: the hooks, alternating.
        let got = muir_expand(&part(&[5])).unwrap();
        let mut got: Vec<_> = got.iter().map(|(l, c)| (l.parts().to_vec(), *c)).collect();
        got.sort_by_key(|(l, _)| l.len());
        let want = vec![
            (vec![5], 1),
            (vec![4, 1], -1),
            (vec![3, 1, 1], 1),
            (vec![2, 1, 1, 1], -1),
            (vec![1, 1, 1, 1, 1], 1),
        ];
        assert_eq!(got, want);
    }

    #[test]
    fn round_trips_are_identity() {
        // s → h → s and s → e → s and s → m → s recover the original.
        for parts in [&[2, 1][..], &[3, 1], &[2, 2], &[3, 2, 1]] {
            let s: Schur<i64> = Schur::monomial(part(parts), 1);
            let back_h: Schur<i64> = convert(&Homogeneous::from_schur(&s));
            let back_e: Schur<i64> = convert(&Elementary::from_schur(&s));
            let back_m: Schur<i64> = convert(&Monomial::from_schur(&s));
            assert_eq!(back_h, s, "s→h→s at {:?}", parts);
            assert_eq!(back_e, s, "s→e→s at {:?}", parts);
            assert_eq!(back_m, s, "s→m→s at {:?}", parts);
        }
    }

    #[test]
    fn power_sum_round_trip_over_rationals() {
        // s → p → s over ℚ is the identity.
        for parts in [&[2, 1][..], &[3], &[2, 2], &[3, 1]] {
            let s: Schur<Rational> =
                Schur::monomial(part(parts), Rational::one());
            let p: PowerSum<Rational> = PowerSum::from_schur(&s);
            let back: Schur<Rational> = p.to_schur();
            assert_eq!(back, s, "s→p→s at {:?}", parts);
        }
    }

    #[test]
    fn power_sum_two_is_s2_minus_s11() {
        // p_2 = s_2 − s_{11}.
        let p2: PowerSum<i64> = PowerSum::monomial(part(&[2]), 1);
        let s = p2.to_schur();
        assert_eq!(s.coeff(&part(&[2])), 1);
        assert_eq!(s.coeff(&part(&[1, 1])), -1);
    }

    #[test]
    fn monomial_multiplication_via_schur() {
        // m_1 · m_1 = 2 m_{11} + m_2.
        let m1: Monomial<i64> = Monomial::monomial(part(&[1]), 1);
        let prod = m1.mul(&m1);
        assert_eq!(prod.coeff(&part(&[1, 1])), 2);
        assert_eq!(prod.coeff(&part(&[2])), 1);
        assert_eq!(prod.terms().len(), 2);
    }
}
