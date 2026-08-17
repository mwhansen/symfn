//! The basis-change engine: conversions between all six classical bases,
//! routed through the Schur hub except where a pair has a direct rule.
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
//! The hub is skipped for the six ordered pairs among h, e and p, which are
//! each free multiplicative bases and so need only their generators expanded in
//! each other ([`FromSchur::from_basis`]):
//!
//! | conversion            | method                                  | ring    |
//! |-----------------------|-----------------------------------------|---------|
//! | p → h, p → e          | Newton's identity, one generator at a time | ℤ    |
//! | h → p, e → p          | the two halves of Cauchy, p_μ/z_μ       | **ℚ**   |
//! | h → e, e → h          | Newton's identity for the h/e pair       | ℤ       |
//!
//! Going through Schur instead is not merely a longer road: it *inflates*. A
//! power-sum element with a handful of terms becomes a Schur element with p(n)
//! of them, and the contraction that follows is then p(n) determinants. That
//! cost is what the direct rules remove (`docs/record/transitions.md`).
//!
//! Only the conversions *into* p divide, and what they need is a [`QAlgebra`] —
//! a ring containing ℚ — not a [`Field`](crate::coeff::Field). The division is
//! by z_μ, an *integer*, so `ℚ[t]` and `ℚ[q,t]` qualify even though neither is
//! a field. Every other path stays exact over ℤ.
//!
//! In Sage these conversions are basis coercions —
//! `Sym = SymmetricFunctions(QQ)`, then `h(s[lam])`. The ordered pairs among
//! s, e, h, m and p dispatch into Symmetrica's C. Symmetrica has no forgotten
//! basis, so the f pairs fall back to Sage's own Python
//! (`scripts/compare_sage.py`, `docs/record/transitions.md`).

use std::collections::{BTreeMap, HashMap};

use crate::character::character_in;
use crate::coeff::{QAlgebra, Ring};
use crate::fasthash::Map;
use crate::interrupt;
use crate::kostka::kostka;
use crate::memo::{inverse_kostka_row_cached, jt_row_cached, lex_parts_cached, partitions_cached};
use crate::partition::Partition;
use crate::sym::{
    Elementary, Forgotten, Homogeneous, Monomial, PowerSum, Schur, SymAlgebra, SymFn,
};

/// Expand `self` into the Schur basis.
pub trait ToSchur<C: Ring> {
    /// The same element, written in the Schur basis.
    fn to_schur(&self) -> Schur<C>;
}

/// Contract a Schur element into `Self`'s basis.
pub trait FromSchur<C: Ring>: Sized {
    /// The same element, written in `Self`'s basis.
    fn from_schur(s: &Schur<C>) -> Self;

    /// A route from the basis whose [`SymFn::SYMBOL`] is `src` that reaches
    /// `Self` **without** passing through Schur, or `None` when there is none.
    ///
    /// The hub is a cost, not just a convention: `p → s → h` inflates a
    /// power-sum element with a handful of terms into a Schur element with
    /// p(n) of them, and only then contracts. Overriding this is how a pair
    /// with a direct rule skips that, and [`convert`] prefers it whenever it
    /// answers.
    ///
    /// It is keyed on the *source symbol* rather than on the source type
    /// because that is what makes the bound land where it belongs: h → p
    /// divides by z_μ and so is only available over a [`QAlgebra`], while
    /// p → h is integral. A generic `convert` cannot state either condition,
    /// and each impl states its own.
    ///
    /// # Examples
    ///
    /// p_2 = 2h_2 − h_{(1,1)}, whose signs and the factor 2 are what separate
    /// Newton's identity from the two rival expansions that agree on p_1:
    ///
    /// ```
    /// use symfn::convert::{convert, FromSchur, ToSchur};
    /// use symfn::partition::Partition;
    /// use symfn::sym::{Homogeneous, PowerSum, SymFn};
    ///
    /// let p2: PowerSum<i64> = PowerSum::monomial(Partition::new([2]), 1);
    /// let h: Homogeneous<i64> = convert(&p2);
    /// assert_eq!(h.coeff(&Partition::new([2])), 2);
    /// assert_eq!(h.coeff(&Partition::new([1, 1])), -1);
    ///
    /// // The same answer the Schur hub gives, which is what `convert` skipped.
    /// assert_eq!(h, Homogeneous::from_schur(&p2.to_schur()));
    /// ```
    fn from_basis(src: &'static str, terms: &BTreeMap<Partition, C>) -> Option<Self> {
        let _ = (src, terms);
        None
    }
}

/// Convert between any two bases: by the direct rule when the pair has one,
/// otherwise by composing through Schur.
///
/// Every pair stays exact over ℤ except the conversions into p. Those divide
/// by z_μ, so their [`FromSchur`] impls are bound on [`QAlgebra`] rather than
/// on [`Ring`].
pub fn convert<C, A, B>(a: &A) -> B
where
    C: Ring,
    A: SymFn<C> + ToSchur<C>,
    B: FromSchur<C>,
{
    match B::from_basis(A::SYMBOL, a.terms()) {
        Some(b) => b,
        None => B::from_schur(&a.to_schur()),
    }
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
/// That is what keeps it cheap. A generic Laplace expansion over the
/// symmetric-function algebra instead clones an (n−1)×(n−1) matrix of
/// *polynomials* at every node and does a full polynomial add and multiply per
/// term, at a cost factorial in the matrix size. The matrix size is ℓ(λ) for
/// s → h but **λ₁** for s → e, since that one is built from the conjugate.
/// Same helper, opposite behavior: s → h stayed fast on the wide-but-shallow
/// shapes a degree ladder produces. s → e crossed over and fell behind Sage
/// past degree 16 (`scripts/compare_sage.py`).
///
/// Enumerating permutations directly also prunes where the determinant is
/// sparse: a negative index means the entry is zero, so that whole subtree is
/// skipped rather than multiplied out.
///
/// **Rows are assigned from the last to the first, and that is the whole
/// performance story.** Entry (i, j) vanishes when j < i - c_i, and c is weakly
/// decreasing, so i - c_i *increases* with i: the constraint tightens as the
/// row index grows. Walking rows forward therefore starts at the least
/// constrained row - row 0 accepts any column - and only meets the dead ends
/// near the leaves, long after the branching has happened. Walking them
/// backwards puts the tightest row first, so whole subtrees die at depth 1.
///
/// The difference is not a constant factor. Symmetrica's `tsh_jt` builds the
/// transposed matrix (lambda_i + i - j, invalid when j > lambda_i + i), which
/// puts its tight row first for free: the same determinant and the same
/// algorithm, with the orientation accounting for the entire gap
/// (`docs/record/transitions.md`).
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

/// The same determinant, expanded into ordered **compositions** instead of
/// sorted multisets.
///
/// [`jt_terms`] can sort each permutation term because its basis multiplies by
/// concatenation — `h_a · h_b` does not care which came first. The Macdonald
/// operator does: its eigenvalue `[|α|] = Σ_i q^{α_i} t^{n−i}` reads α
/// *positionally*, so the same multiset in two orders is two different
/// polynomials. This visits each permutation term with α laid out **by column**
/// — the index the consumer needs; see the note at the assignment.
///
/// Deliberately a sibling rather than a refactor of [`jt_terms`] into a shared
/// callback. That one is on the `s → h` and `s → e` hot paths, where the whole
/// point is that no polynomial arithmetic happens and terms aggregate into a
/// `HashMap` as they are found. Threading a caller's closure through it would
/// put an indirect call in that inner loop. That path took `s → e` on λ=(14)
/// from 1.5 seconds to microseconds (`docs/record/transitions.md`). The pruning
/// argument in [`jt_terms`] — rows assigned last to first, so the tightest
/// constraint is met at depth 1 — is the part that matters and it is reproduced
/// here.
pub(crate) fn jt_compositions(c: &[u32], visit: &mut impl FnMut(&[u32], i64)) {
    if c.is_empty() {
        visit(&[], 1);
        return;
    }
    let mut used = vec![false; c.len()];
    let mut alpha = vec![0u32; c.len()];
    jtc_rec(c, c.len(), &mut used, &mut alpha, 1, visit);
}

// Jacobi–Trudi index arithmetic. `c[row]`, `row` and `j` are all bounded by
// ℓ(λ) and λ₁, both `u32` in `Partition`, and the guarded subtraction is what
// keeps the result non-negative — the value is a part, not a coefficient.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
fn jtc_rec(
    c: &[u32],
    i: usize,
    used: &mut [bool],
    alpha: &mut [u32],
    sign: i64,
    visit: &mut impl FnMut(&[u32], i64),
) {
    if i == 0 {
        visit(alpha, sign);
        return;
    }
    let row = i - 1;
    let mut below = 0i64;
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
        used[j] = true;
        // Indexed by **column**, not by row, and not pushed in visit order.
        // The consumer pairs position with a power of t, and in Lapointe-
        // Lascoux-Morse (IMRN 1998 no. 18, 957-978; arXiv:math/9808050) that
        // power comes from which column of the determinant the permutation
        // selected -- their 3.5, via the formal operators of their 2.4. See
        // `macop::eigenvalue`. Both other choices produce a plausible
        // triangular matrix with distinct eigenvalues and the wrong Macdonald
        // polynomials.
        alpha[j] = (c[row] as i64 - row as i64 + j as i64) as u32;
        let s = if below % 2 == 0 { sign } else { -sign };
        jtc_rec(c, row, used, alpha, s, visit);
        used[j] = false;
    }
}

/// Assign row `i - 1`, rows `i..` being already placed.
// As `jtc_rec`: subscripts of h, bounded by the shape, never coefficients.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
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

/// [`jt_terms`] behind the row cache, which is what makes a *repeated* small
/// conversion cheap. The determinant does not depend on the coefficient ring or
/// on the coefficient, so the row is shared by every caller that mentions the
/// index.
fn jt_row(index: &Partition) -> std::sync::Arc<Vec<(Partition, i64)>> {
    jt_row_cached(index, || jt_terms(index.parts()))
}

// --- Homogeneous <-> Schur --------------------------------------------------

impl<C: Ring> ToSchur<C> for Homogeneous<C> {
    fn to_schur(&self) -> Schur<C> {
        // h_λ = ∏_i h_{λ_i} = ∏_i s_{(λ_i)}, each factor a Pieri step.
        expand_multiplicative(self, false)
    }
}

impl<C: Ring> FromSchur<C> for Homogeneous<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = det(h_{λ_i − i + j}); the matrix is ℓ(λ)×ℓ(λ).
        contract_multiplicative(s, false)
    }

    fn from_basis(src: &'static str, terms: &BTreeMap<Partition, C>) -> Option<Self> {
        match src {
            "p" => Some(multiplicative_route(terms, |n| {
                power_generator::<C, Self>(n, false)
            })),
            "e" => Some(multiplicative_route(terms, flip_generator::<C, Self>)),
            _ => None,
        }
    }
}

// --- Elementary <-> Schur ---------------------------------------------------

impl<C: Ring> ToSchur<C> for Elementary<C> {
    fn to_schur(&self) -> Schur<C> {
        // e_λ = ∏_i e_{λ_i} = ∏_i s_{(1^{λ_i})} (single columns), each factor
        // a Pieri step on the conjugate side.
        expand_multiplicative(self, true)
    }
}

/// h_λ and e_λ in the Schur basis, by **Pieri steps shared across terms**.
///
/// Both are products of one-row or one-column Schur functions. Building them
/// term by term — `Schur::unit()`, then one `Schur::mul` per part, then
/// `out.add(&prod.scale(c))` — puts almost the whole cost on those two lines
/// and almost all of the self time in the allocator
/// (`examples/profile_convert.rs`, `loop e2s 20`). Two things are wrong with
/// it.
///
/// * **The multiply was the general Littlewood–Richardson engine.**
///   `Schur::mul` routes to `AutoLr`, which built and expanded a skew shape
///   for what is a Pieri step: multiplying by s_{(k)} adds a horizontal
///   k-strip and by s_{(1^k)} a vertical one, both a direct enumeration with
///   no LR machinery under them. That is where the time was.
/// * **The layer was rebuilt per term.** `terms()` is a `BTreeMap` keyed by
///   `Partition`, which orders lexicographically by parts, so partitions
///   sharing their first `depth` parts are *already contiguous* — no sort
///   needed, unlike [`p_expand_shared`], which is handed a `Vec`. Each such
///   run continues from one layer instead of rebuilding it from the unit.
///
/// The sharing is the smaller half and was measured before it was written.
/// What it removes are the short cheap prefixes, while the leaves — the
/// expensive steps — are exactly what no two terms share. So the saving
/// weighted by the degree each step multiplies into is much smaller than the
/// step count suggests (`docs/record/transitions.md`). It is kept because it
/// is nearly free once the traversal is written this way, not because it
/// carries the win.
/// The layer is a **β-mask**, not a `Schur`. With the Pieri step in place the
/// profile is mostly allocator and `memmove`, with actual strip enumeration a
/// small minority of it: a `Schur<C>` is a `BTreeMap<Partition, C>`, so every
/// shape a step emitted allocated a heap `Vec<u32>`, sorted it, and memmoved
/// its way into a B-tree. On the β-mask the same step is bit arithmetic on a
/// `u64` in a `Map`, and partitions are built once per *output* term rather
/// than once per emitted shape — the same trade `p_expand` and `muir_expand`
/// already make.
///
/// Layer coefficients are `i128`, not `C`. Pieri's structure constants are all
/// 1, so a layer coefficient is a plain multiplicity. For h_μ that multiplicity
/// is the Kostka number K_{λμ}, bounded by f^λ ≤ √(n!). This path only runs for
/// n ≤ [`MASK_LIMIT`] = 32, where √(32!) ≈ 1.6·10¹⁸ sits far inside `i128`. So
/// no ring arithmetic happens in the sweep at all; `C` is touched once per
/// output term. Same narrowing argument, and the same bound, as [`p_expand`].
///
/// Terms are batched by degree because the mask width is |λ|, exactly as
/// `PowerSum::to_schur` batches for the same reason. Degrees past the mask
/// width take the partition-keyed [`expand_shared`] below, which has no wall.
fn expand_multiplicative<C: Ring, S: SymFn<C>>(x: &S, vertical: bool) -> Schur<C> {
    let mut out = Schur::zero();
    // `terms()` is ordered by `Partition`, i.e. lexicographically by parts, and
    // pushing in that order keeps each degree's group in it — which is what
    // makes shared prefixes contiguous without a sort.
    let mut by_degree: BTreeMap<u32, Vec<(&Partition, &C)>> = BTreeMap::new();
    for (lambda, c) in x.terms() {
        by_degree
            .entry(lambda.size())
            .or_default()
            .push((lambda, c));
    }
    for (n, items) in by_degree {
        let l = n as usize;
        if l == 0 {
            for (_, c) in &items {
                out.add_term(Partition::default(), (*c).clone());
            }
        } else if l > MASK_LIMIT {
            expand_shared(&items, 0, &Schur::unit(), vertical, &mut out);
        } else {
            // β-numbers of ∅ with l slots: {0, 1, …, l−1}.
            let mut root: Map<u64, i128> = Map::default();
            root.insert((1u64 << l) - 1, 1);
            expand_shared_masks(&items, 0, &root, l, vertical, &mut out);
        }
    }
    out
}

/// [`expand_shared`] on the β-mask layer.
fn expand_shared_masks<C: Ring>(
    items: &[(&Partition, &C)],
    depth: usize,
    layer: &Map<u64, i128>,
    l: usize,
    vertical: bool,
    out: &mut Schur<C>,
) {
    let mut i = 0;
    while i < items.len() && items[i].0.len() == depth {
        for (&mask, &v) in layer {
            if v != 0 {
                out.add_term(mask_to_partition(mask, l), C::from_i128(v).mul(items[i].1));
            }
        }
        i += 1;
    }
    while i < items.len() {
        let k = items[i].0.part(depth);
        let start = i;
        while i < items.len() && items[i].0.part(depth) == k {
            i += 1;
        }
        let next = pieri_masks(layer, k, vertical);
        expand_shared_masks(&items[start..i], depth + 1, &next, l, vertical, out);
    }
}

/// One Pieri step on a layer of β-masks.
fn pieri_masks(cur: &Map<u64, i128>, k: u32, vertical: bool) -> Map<u64, i128> {
    // As `p_step`: the layer grows through a sweep, so sizing to the input is
    // a floor on the output rather than a guess.
    let mut next: Map<u64, i128> = Map::with_capacity_and_hasher(cur.len() * 2, Default::default());
    let mut shapes = Vec::new();
    for (&mask, &c) in cur {
        shapes.clear();
        strip_masks(mask, k, vertical, &mut shapes);
        for &m in &shapes {
            *next.entry(m).or_insert(0) += c;
        }
    }
    next
}

/// Every β-mask reachable from `mask` by adding a horizontal — or, if
/// `vertical`, a vertical — strip of `k` cells.
///
/// **The two strip conditions are different constraints on the β-set, and this
/// is the trap.** Writing one and reusing it for the other produces plausible
/// wrong answers, not errors, so both are derived here rather than shared.
///
/// With β_i = λ_i + (l−1−i), strictly decreasing:
///
/// * **Horizontal.** μ/λ is a horizontal strip iff μ_i ≥ λ_i and μ_i ≤ λ_{i−1},
///   which in β reads β^λ_{i−1} > β^μ_i ≥ β^λ_i. Each β therefore moves up
///   *within its own interval*, bounded by the **original** β above it, and the
///   intervals are disjoint — so the choices are independent and no collision
///   test is needed.
/// * **Vertical.** μ_i ∈ {λ_i, λ_i + 1}, so each β either stays or moves up by
///   one, with no interval bound at all. What constrains it is that β^μ must
///   stay strictly decreasing, which bites only when two β are adjacent: the
///   lower may not step onto a higher one that stayed. Rows are walked from the
///   **highest β down**, so everything above is already final and that test is
///   sound — the same argument [`muir_rec`] makes for its collision check.
fn strip_masks(mask: u64, k: u32, vertical: bool, out: &mut Vec<u64>) {
    let mut bits: Vec<u32> = Vec::with_capacity(mask.count_ones() as usize);
    let mut rest = mask;
    while rest != 0 {
        let b = 63 - rest.leading_zeros();
        bits.push(b);
        rest &= !(1u64 << b);
    }
    // Suffix capacity: the most the β from position i down can absorb between
    // them. The topmost is unbounded, so only the deeper levels can dead-end.
    let mut cap = vec![0u32; bits.len() + 1];
    for i in (1..bits.len()).rev() {
        let room = if vertical {
            1
        } else {
            bits[i - 1] - 1 - bits[i]
        };
        cap[i] = cap[i + 1].saturating_add(room);
    }
    cap[0] = u32::MAX;

    fn rec(
        bits: &[u32],
        cap: &[u32],
        i: usize,
        left: u32,
        vertical: bool,
        acc: u64,
        out: &mut Vec<u64>,
    ) {
        if i == bits.len() {
            if left == 0 {
                out.push(acc);
            }
            return;
        }
        if left > cap[i] {
            return;
        }
        let b = bits[i];
        // The ceiling is the *original* β above, which is what makes the
        // horizontal intervals disjoint; see the note above.
        let hi = if vertical {
            b + 1
        } else if i == 0 {
            b + left
        } else {
            (bits[i - 1] - 1).min(b + left)
        };
        for nb in b..=hi {
            let used = nb - b;
            if used > left {
                break;
            }
            // Only the vertical case can collide, and only with a β that stayed.
            if vertical && used > 0 && acc >> nb & 1 == 1 {
                continue;
            }
            rec(
                bits,
                cap,
                i + 1,
                left - used,
                vertical,
                acc | (1u64 << nb),
                out,
            );
        }
    }
    rec(&bits, &cap, 0, k, vertical, 0, out);
}

/// Walk one prefix group, mirroring [`p_expand_shared`]'s shape.
fn expand_shared<C: Ring>(
    items: &[(&Partition, &C)],
    depth: usize,
    layer: &Schur<C>,
    vertical: bool,
    out: &mut Schur<C>,
) {
    let mut i = 0;
    // Partitions that end here: the layer is their whole product. Emitted
    // term by term into `out` rather than through `out.add(&layer.scale(c))`,
    // which allocated a scaled copy and then a merged map per input term.
    while i < items.len() && items[i].0.len() == depth {
        for (lambda, v) in layer.terms() {
            out.add_term(lambda.clone(), v.mul(items[i].1));
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
        let next = pieri_step(layer, k, vertical);
        expand_shared(&items[start..i], depth + 1, &next, vertical, out);
    }
}

/// Multiply a Schur element by s_{(1^k)} (`vertical`) or s_{(k)}.
///
/// Pieri: s_λ · h_k = Σ s_μ over μ/λ a horizontal k-strip, and s_λ · e_k the
/// same over vertical strips. Coefficients are all 1, so this only ever moves
/// the caller's coefficient onto a new index — no ring multiplication happens.
fn pieri_step<C: Ring>(cur: &Schur<C>, k: u32, vertical: bool) -> Schur<C> {
    let mut next = Schur::zero();
    let mut shapes = Vec::new();
    for (lambda, c) in cur.terms() {
        shapes.clear();
        if vertical {
            vertical_strips(lambda.parts(), k, &mut shapes);
        } else {
            horizontal_strips(lambda.parts(), k, &mut shapes);
        }
        for mu in shapes.drain(..) {
            next.add_term(mu, c.clone());
        }
    }
    next
}

/// Every μ ⊇ λ with μ/λ a horizontal strip of `k` cells.
///
/// The strip condition is the interlacing μ₁ ≥ λ₁ ≥ μ₂ ≥ λ₂ ≥ … — at most one
/// new cell per column. Row `i` is therefore capped by λ_{i−1}, the *old*
/// value, and μ comes out weakly decreasing for free: μ_i ≤ λ_{i−1} ≤ μ_{i−1}.
fn horizontal_strips(lambda: &[u32], k: u32, out: &mut Vec<Partition>) {
    fn rec(i: usize, left: u32, lambda: &[u32], cur: &mut Vec<u32>, out: &mut Vec<Partition>) {
        // A strip can open at most one new row, so the walk ends one past λ.
        if i == lambda.len() + 1 {
            if left == 0 {
                out.push(Partition::new(cur.iter().copied()));
            }
            return;
        }
        let lam_i = lambda.get(i).copied().unwrap_or(0);
        let upper = if i == 0 {
            lam_i + left
        } else {
            lambda[i - 1].min(lam_i + left)
        };
        for v in lam_i..=upper {
            cur.push(v);
            rec(i + 1, left - (v - lam_i), lambda, cur, out);
            cur.pop();
        }
    }
    let mut cur = Vec::with_capacity(lambda.len() + 1);
    rec(0, k, lambda, &mut cur, out);
}

/// Every μ ⊇ λ with μ/λ a vertical strip of `k` cells.
///
/// No two cells in one row, so μ_i ∈ {λ_i, λ_i + 1} and the only other
/// constraint is that μ stay weakly decreasing. Unlike a horizontal strip this
/// may open up to `k` new rows, each holding one cell.
fn vertical_strips(lambda: &[u32], k: u32, out: &mut Vec<Partition>) {
    fn rec(
        i: usize,
        rows: usize,
        left: u32,
        prev: u32,
        lambda: &[u32],
        cur: &mut Vec<u32>,
        out: &mut Vec<Partition>,
    ) {
        if left == 0 {
            // Nothing left to place: the remaining rows keep their old values.
            let mut done = cur.clone();
            done.extend_from_slice(&lambda[i.min(lambda.len())..]);
            out.push(Partition::new(done));
            return;
        }
        // Each remaining row can absorb at most one cell.
        if i == rows || left as usize > rows - i {
            return;
        }
        let lam_i = lambda.get(i).copied().unwrap_or(0);
        if lam_i <= prev {
            cur.push(lam_i);
            rec(i + 1, rows, left, lam_i, lambda, cur, out);
            cur.pop();
        }
        if lam_i < prev {
            cur.push(lam_i + 1);
            rec(i + 1, rows, left - 1, lam_i + 1, lambda, cur, out);
            cur.pop();
        }
    }
    let rows = lambda.len() + k as usize;
    let mut cur = Vec::with_capacity(rows);
    rec(0, rows, k, u32::MAX, lambda, &mut cur, out);
}

impl<C: Ring> FromSchur<C> for Elementary<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = det(e_{λ'_i − i + j}); the matrix is λ₁×λ₁, being built from the
        // conjugate.
        contract_multiplicative(s, true)
    }

    fn from_basis(src: &'static str, terms: &BTreeMap<Partition, C>) -> Option<Self> {
        match src {
            "p" => Some(multiplicative_route(terms, |n| {
                power_generator::<C, Self>(n, true)
            })),
            "h" => Some(multiplicative_route(terms, flip_generator::<C, Self>)),
            _ => None,
        }
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
/// the swap, so a single routine serves both directions. Each step multiplies
/// by one generator, which in a multiplicative basis is a multiset union, so
/// this is a linear recursion over cheap products rather than anything
/// determinantal.
///
/// `gen[n]` is then the one-part generator of the *source* basis expanded in
/// the target, and a multi-part index is the product of those.
fn flip_generator<C: Ring, S: SymAlgebra<C>>(n: u32) -> S {
    // The table is the *same* for h→e and e→h, and its coefficients are
    // integers independent of `C`, so it is computed once in i64 and injected.
    // Recomputing it per call was the dominant cost of a flipped conversion:
    // the recursion touches O(n²) products over elements with up to p(n) terms,
    // which is cheap once and wasteful on every call.
    let table = flip_table(n);
    let mut x = S::zero();
    for (mu, c) in &table[n as usize] {
        x.add_term(mu.clone(), C::from_i64(*c));
    }
    x
}

thread_local! {
    /// h_n in the e-basis (equivalently e_n in the h-basis) for n = 0.., as
    /// integer term lists. Grown monotonically, never invalidated: these are
    /// fixed integers, so a longer table subsumes a shorter one.
    static FLIP: std::cell::RefCell<Vec<Vec<(Partition, i64)>>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

// `k` runs to `n`, the degree, which is a `u32` throughout the crate.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
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

/// Pairs each multiplicative basis with the other, h with e, so a contraction
/// can be computed in whichever Jacobi–Trudi matrix is smaller.
///
/// An associated type rather than a flag, so "compute in whichever Jacobi–Trudi
/// matrix is smaller and flip back" is expressible without either basis naming
/// the other at the call site — and so the recursion that does it is checked to
/// terminate by the compiler rather than by argument.
pub trait Dual<C: Ring>: SymAlgebra<C> {
    /// The other multiplicative basis: [`Elementary`] for [`Homogeneous`], and
    /// back again, so the pair is its own inverse.
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
    multiplicative_route(x.terms(), flip_generator::<C, B>)
}

/// Rewrite a **free multiplicative** basis in another, given the source's
/// generators already expanded in the target.
///
/// h, e and p are each free on one generator per positive integer, and each
/// multiplies by index concatenation: x_μ · x_ν = x_{μ ∪ ν}. So an index of
/// several parts is the product of its parts' expansions, and the whole
/// transition is fixed by `gens[n]`, the source's n-th generator written in the
/// target. No determinant, no Schur hub, and no term of the target is ever
/// touched that the answer does not mention.
///
/// `gens[0]` is the unit; a partition has no zero part, so it is only ever read
/// when the index is empty.
fn multiplicative_route<C: Ring, B: SymAlgebra<C>>(
    terms: &BTreeMap<Partition, C>,
    generator: impl Fn(u32) -> B,
) -> B {
    // Only the parts actually mentioned. A single h_(18) names one generator,
    // and building the whole ladder up to it is four times the work of the
    // conversion — the peel that motivates these routes asks for one index at
    // a time, so that ladder would be rebuilt on every call.
    let mut gens: BTreeMap<u32, B> = BTreeMap::new();
    for mu in terms.keys() {
        for &part in mu.parts() {
            gens.entry(part).or_insert_with(|| generator(part));
        }
    }
    let mut out = B::zero();
    for (mu, c) in terms {
        let mut prod = B::unit();
        for &part in mu.parts() {
            prod = prod.times(&gens[&part]);
        }
        out = out.add(&prod.scale(c));
    }
    out
}

/// p_n expanded in the h-basis, as an integer term list.
///
/// Newton's identity n·h_n = Σ_{i=1}^{n} p_i h_{n−i}, rearranged so that the
/// unknown is alone:
///
/// ```text
///   p_n = n·h_n − Σ_{i=1}^{n−1} h_{n−i} · p_i
/// ```
///
/// Each term on the right multiplies an already-known p_i by the single
/// generator h_{n−i}, which in a multiplicative basis appends one part. So the
/// whole table is built by appending parts to earlier rows — the same shape as
/// [`flip_table`], and for the same reason.
///
/// The **e-basis needs no second table**: log E(t) = Σ (−1)^{r−1} p_r t^r / r
/// against log H(t) = Σ p_r t^r / r differs only by that sign, so p_n in e is
/// p_n in h with every coefficient scaled by (−1)^{n−1} and the letter changed.
///
/// Coefficients are `i128` rather than `i64` because the largest of them is
/// n·(ℓ−1)!/∏ m_i!, which counts compositions and so grows like 2^n; i64 would
/// cap the table near degree 60 while every other path here keeps going.
#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn power_in_h_table(upto: u32) -> Vec<Vec<(Partition, i128)>> {
    POWER_IN_H.with(|cell| {
        let mut t = cell.borrow_mut();
        if t.is_empty() {
            t.push(vec![(Partition::default(), 1)]);
        }
        while t.len() <= upto as usize {
            let n = t.len();
            let mut acc: HashMap<Partition, i128> = HashMap::new();
            acc.insert(Partition::new([n as u32]), n as i128);
            for i in 1..n {
                for (mu, c) in &t[i] {
                    let mut parts = mu.parts().to_vec();
                    parts.push((n - i) as u32);
                    *acc.entry(Partition::new(parts)).or_insert(0) -= c;
                }
            }
            t.push(acc.into_iter().filter(|(_, c)| *c != 0).collect());
        }
        t[..=upto as usize].to_vec()
    })
}

thread_local! {
    /// p_n in the h-basis for n = 0.., as integer term lists. Grown
    /// monotonically and never invalidated: these are fixed integers, so a
    /// longer table subsumes a shorter one.
    static POWER_IN_H: std::cell::RefCell<Vec<Vec<(Partition, i128)>>> =
        const { std::cell::RefCell::new(Vec::new()) };
}

/// p_n in the h-basis (`dual` false) or the e-basis (`dual` true), injected
/// into the coefficient ring.
fn power_generator<C: Ring, S: SymAlgebra<C>>(n: u32, dual: bool) -> S {
    // (−1)^{n−1} for n ≥ 1; the n = 0 row is the unit and unsigned.
    let flip = dual && n.is_multiple_of(2) && n > 0;
    let table = power_in_h_table(n);
    let mut x = S::zero();
    for (mu, c) in &table[n as usize] {
        x.add_term(mu.clone(), C::from_i128(if flip { -c } else { *c }));
    }
    x
}

/// h_n (`dual` false) or e_n (`dual` true) in the power-sum basis.
///
/// ```text
///   h_n = Σ_{μ ⊢ n} p_μ / z_μ        e_n = Σ_{μ ⊢ n} (−1)^{n−ℓ(μ)} p_μ / z_μ
/// ```
///
/// the two halves of the Cauchy identity, and the only direction of this family
/// that divides. The division is by z_μ, an integer, which is what a
/// [`QAlgebra`] promises — and it goes through [`Partition::div_by_z`], which
/// divides by z_μ's factors one at a time rather than forming z_μ, so the
/// transition has no degree ceiling of its own.
fn multiplicative_in_power<C: QAlgebra>(n: u32, dual: bool) -> PowerSum<C> {
    let mut x = PowerSum::zero();
    for mu in partitions_cached(n).iter() {
        let sign = if dual && (n as usize - mu.len()) % 2 == 1 {
            -1
        } else {
            1
        };
        x.add_term(mu.clone(), mu.div_by_z(&C::from_i64(sign)));
    }
    x
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
/// The cost is therefore set by λ₁, not by the degree, and it goes
/// exponential once λ₁ is large (`docs/record/transitions.md`).
///
/// `scripts/compare_sage.py` builds its shapes with several rows, so λ₁ stays
/// under 8 and no wide shape is covered there. `scripts/check_backend.py` lets
/// Sage pick the inputs, and reaches them.
///
/// Taking the smaller of the two matrices (see [`flip_basis`]) removes most of
/// the problem on its own, since a shape that is bad for one direction is good
/// for the other. What survives is the family bad for *both* — hooks, where
/// ℓ(λ) and λ₁ are each about |λ|/2.
///
/// The limit is set high because **the sweep is nearly always the worse of the
/// two**, which was not the expectation: a Muir sweep is p(n) expansions, and
/// p(n) grows at every degree while the determinant is only bad once the
/// matrix is genuinely large. So the sweep is a backstop against a
/// pathological shape, not a fast path, and the threshold sits where the
/// determinant finally stops being cheap (`docs/record/transitions.md`).
const JT_LIMIT: usize = 14;

/// Matrix size at which flipping to the conjugate basis starts to pay.
///
/// Flipping is *not* free — see [`flip_basis`] — and taking the smaller matrix
/// whenever it is smaller at all was measured worse in aggregate over every
/// partition of degree 20 (`docs/record/transitions.md`). A shape like
/// (5,5,5,5) has matrices of 4 and 5, and a 5-wide determinant is far cheaper
/// than expanding an h-element of degree 20 into e.
///
/// So the flip is reserved for the cases where the determinant is genuinely
/// exponential and the conjugate collapses it: a single row of degree 24,
/// whose h-expansion is the one term h_24. Below this size the direct
/// determinant always wins.
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
        let index = if dual {
            lambda.conjugate()
        } else {
            lambda.clone()
        };
        let other = index.conjugate();
        if index.len().min(other.len()) > JT_LIMIT {
            swept
                .entry(index.size())
                .or_default()
                .push((index, c.clone()));
        } else if index.len() >= FLIP_MIN && other.len() < index.len() {
            // The conjugate determinant is smaller: compute there and flip.
            crossed.add_term(lambda.clone(), c.clone());
        } else {
            // Accumulated into `out` directly rather than built and merged:
            // one element per input term is an allocation and a second pass
            // over the row, and this loop runs once per term of the input.
            for (mu, v) in jt_row(&index).iter() {
                out.add_term(mu.clone(), C::from_i64(*v).mul(c));
            }
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
            // determinant has no wall, only a cost.
            None => {
                for (index, c) in &targets {
                    for (mu, v) in jt_row(index).iter() {
                        out.add_term(mu.clone(), C::from_i64(*v).mul(c));
                    }
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
            if n as usize <= WIDE_MASK_LIMIT {
                by_degree.entry(n).or_default().push((mu, c));
            } else {
                // Degree past what any β-mask holds: fall back to characters,
                // which are exact in `C` and so stay correct for bignum rings.
                // This is p(n) character recursions in the coefficient ring per
                // term, against one shared sweep below, and it is the slowest
                // reachable path in the crate — the reason the sweep is worth
                // widening at all (`docs/record/plethysm.md`).
                for lambda in partitions_cached(n).iter() {
                    interrupt::poll();
                    let chi = character_in::<C>(lambda, mu);
                    if !chi.is_zero() {
                        out.add_term(lambda.clone(), chi.mul(c));
                    }
                }
            }
        }
        // Drained in ascending degree. Taking the keys first and removing as we
        // go visits each bucket once, in order, without holding a borrow.
        let mut degrees: Vec<u32> = by_degree.keys().copied().collect();
        degrees.sort_unstable();
        for n in degrees {
            let Some(mut items) = by_degree.remove(&n) else {
                continue;
            };
            // Descending part order, so the longest common prefixes are shared.
            items.sort_by(|a, b| a.0.parts().cmp(b.0.parts()));
            let l = n as usize;
            if l == 0 {
                for (_, c) in &items {
                    out.add_term(Partition::default(), (*c).clone());
                }
                continue;
            }
            // The narrow mask carries every degree it can. Widening is not free
            // — `docs/record/plethysm.md` has what it costs at equal degree —
            // and degrees past 32 are the only ones that need it.
            if l <= MASK_LIMIT {
                sweep_degree::<u64, C>(&items, l, &mut out);
            } else {
                sweep_degree::<u128, C>(&items, l, &mut out);
            }
        }
        out
    }
}

/// Expand one degree's batch of p_μ, in whichever mask width holds it.
///
/// Split out of [`ToSchur::to_schur`] so the two widths are one body rather
/// than two, with the choice made once at the call above.
fn sweep_degree<M: Beta, C: Ring>(items: &[(&Partition, &C)], l: usize, out: &mut Schur<C>) {
    let mut root: Map<M, i128> = Map::default();
    root.insert(M::low_ones(l), 1);
    if integral_sweep(items, l, &root, out) {
        return;
    }
    // Accumulate on the β-mask, not on a Partition. Every leaf of the
    // traversal touches the whole layer, so keying by partition allocated and
    // sorted a fresh Vec — and hashed a heap key — once per (μ, mask) pair to
    // produce a few dozen distinct terms. Masks are words, and the partitions
    // get built once at the end.
    let mut acc: Map<M, C> = Map::default();
    p_expand_shared(items, 0, &root, &mut |c: &&C, mask, chi: &i128| {
        let term = C::from_i128(*chi).mul(c);
        acc.entry(mask).or_insert_with(C::zero).add_assign(&term);
    });
    for (mask, c) in acc {
        if !c.is_zero() {
            out.add_term(mask_to_partition(mask, l), c);
        }
    }
}

/// Expand a batch of p_μ of one degree in a single traversal, sharing the
/// Murnaghan–Nakayama sweep across every common prefix.
///
/// `items` is sorted by part sequence, so partitions agreeing in their first
/// `depth` parts are contiguous; each such run continues from *one* layer
/// instead of rebuilding it. Plethysm is the case that motivates this: it
/// finishes by converting a p-element of degree d·e with dozens of terms, and
/// that conversion is nearly all of its runtime
/// (`docs/record/plethysm.md`).
pub(crate) fn p_expand_shared<M: Beta, C: Ring, T, F>(
    items: &[(&Partition, T)],
    depth: usize,
    layer: &Map<M, C>,
    emit: &mut F,
) where
    F: FnMut(&T, M, &C),
{
    let mut i = 0;
    // Partitions that end here: emit the layer against their coefficient.
    while i < items.len() && items[i].0.len() == depth {
        interrupt::poll();
        for (&mask, chi) in layer {
            if !chi.is_zero() {
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
        let next = p_step(layer, k);
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
/// rational add per (μ, mask) leaf, each normalizing by a gcd — and those gcds
/// with the i128 division underneath them are most of plethysm's runtime
/// (`docs/record/plethysm.md`). Putting every c_μ over one denominator D makes
/// the entire accumulation integer, with a single conversion back per output
/// term.
///
/// D is the lcm of the denominators, which measurement said is the right shape
/// for this: across the plethysms driving this work the lcm *equalled the
/// largest denominator* every time, and never exceeded ~5·10⁵. Overflow is
/// still checked at every step rather than argued away, because the guarantee
/// is only about the cases measured.
fn integral_sweep<M: Beta, C: Ring>(
    items: &[(&Partition, &C)],
    l: usize,
    root: &Map<M, i128>,
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

    let mut acc: Map<M, i128> = Map::default();
    let mut overflow = false;
    p_expand_shared(&paired, 0, root, &mut |&s: &i128, mask, chi: &i128| {
        if overflow {
            return;
        }
        let slot = acc.entry(mask).or_insert(0);
        match s.checked_mul(*chi).and_then(|t| slot.checked_add(t)) {
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
/// in ℓ(μ) passes and never computes a character at all. Same "produce the
/// whole answer in one sweep rather than query it entry by entry" shape as the
/// Kostka and coproduct fixes.
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
/// declines when l > 32 (the β-mask does not fit), and √(32!) ≈ 1.6·10¹⁸ — so
/// on every input it *accepts*, i128 cannot overflow. Callers past that width
/// take the character fallback, which stays exact in `C` for bignum rings.
///
/// The narrowing matters because the ring is the hot loop. Plethysm runs this
/// over ℚ, where each rim hook cost a rational add — a gcd — plus a temporary
/// from negating the coefficient, to combine two integers. That makes p → s
/// nearly all of plethysm's runtime (`docs/record/plethysm.md`).
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

/// The largest degree a `u64` β-mask holds, 32.
///
/// β values run from 0 to at most (l−1) + max part < 2l, so a 64-bit mask holds
/// the whole set for l ≤ 32. Degrees past it are swept in [`u128`] instead, up
/// to [`WIDE_MASK_LIMIT`]; only past *that* does a caller fall back.
pub(crate) const MASK_LIMIT: usize = 32;

/// The largest degree the β-mask sweep accepts at all, 55 — and the ceiling is
/// the accumulator, not the mask.
///
/// A `u128` mask would hold l ≤ 64. What stops it sooner is that the layer
/// accumulates in `i128`: every value in it is a character, so |χ^λ(μ)| ≤ f^λ ≤
/// √(l!), and one slot takes at most l contributions before it settles, giving
/// a transient bound of l·√(l!). That is 6.2·10³⁷ at l = 55 and 1.5·10³⁹ at
/// l = 56, against an `i128` ceiling of 1.7·10³⁸ — so 55 is the last degree on
/// which the sweep provably cannot overflow. R3 is the backstop if this
/// reasoning is ever wrong: the accumulation is checked in every profile, so a
/// broken bound is a panic and not a wrapped character.
pub(crate) const WIDE_MASK_LIMIT: usize = 55;

/// The bit operations the Murnaghan–Nakayama sweep performs on a β-mask.
///
/// A type parameter rather than one wider mask everywhere, because the widening
/// is only ever needed above degree 32 and the sweep below it is the hot path
/// the whole batching design exists to serve (`docs/record/plethysm.md`). Both
/// instantiations monomorphize, so the `u64` sweep compiles to what it compiled
/// to before this trait existed; the measured cost of each is in
/// `docs/record/plethysm.md`.
pub(crate) trait Beta: Copy + Eq + std::hash::Hash {
    /// `(1 << n) - 1`: the β-set of the empty partition with n slots.
    fn low_ones(n: usize) -> Self;
    fn is_zero(self) -> bool;
    fn test(self, n: u32) -> bool;
    fn with_bit(self, n: u32) -> Self;
    fn without_bit(self, n: u32) -> Self;
    /// `self & (((1 << hi) - 1) ^ ((1 << lo) - 1))`: the β strictly between two
    /// positions, which is what a rim hook's height counts.
    fn between(self, lo: u32, hi: u32) -> Self;
    fn trailing_zeros(self) -> u32;
    fn count_ones(self) -> u32;
    /// `self & (self - 1)`: drop the lowest set bit.
    fn clear_lowest(self) -> Self;
    /// Index of the highest set bit; only called on a nonzero mask.
    fn highest(self) -> usize;
}

macro_rules! impl_beta {
    ($t:ty) => {
        impl Beta for $t {
            #[inline]
            fn low_ones(n: usize) -> Self {
                (1 as $t << n) - 1
            }
            #[inline]
            fn is_zero(self) -> bool {
                self == 0
            }
            #[inline]
            fn test(self, n: u32) -> bool {
                self >> n & 1 == 1
            }
            #[inline]
            fn with_bit(self, n: u32) -> Self {
                self | (1 as $t) << n
            }
            #[inline]
            fn without_bit(self, n: u32) -> Self {
                self & !((1 as $t) << n)
            }
            #[inline]
            fn between(self, lo: u32, hi: u32) -> Self {
                self & ((((1 as $t) << hi) - 1) ^ (((1 as $t) << lo) - 1))
            }
            #[inline]
            fn trailing_zeros(self) -> u32 {
                <$t>::trailing_zeros(self)
            }
            #[inline]
            fn count_ones(self) -> u32 {
                <$t>::count_ones(self)
            }
            #[inline]
            fn clear_lowest(self) -> Self {
                self & (self - 1)
            }
            #[inline]
            fn highest(self) -> usize {
                (<$t>::BITS - 1 - <$t>::leading_zeros(self)) as usize
            }
        }
    };
}

impl_beta!(u64);
impl_beta!(u128);

/// One Murnaghan–Nakayama step: multiply a layer of β-masks by p_k.
// β-mask bit positions, bounded by `MASK_LIMIT = 32`.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
pub(crate) fn p_step<M: Beta, C: Ring>(cur: &Map<M, C>, k: u32) -> Map<M, C> {
    // The layer grows monotonically through a sweep, so a default-capacity
    // map rehashes several times per step. Sizing to the input is a floor on
    // the output, not a guess.
    let mut next: Map<M, C> = Map::with_capacity_and_hasher(cur.len() * 2, Default::default());
    for (&mask, c) in cur {
        interrupt::poll();
        let mut rest = mask;
        while !rest.is_zero() {
            let b = rest.trailing_zeros();
            rest = rest.clear_lowest();
            let nb = b + k;
            if mask.test(nb) {
                continue; // that β is taken: no such rim hook
            }
            // Height = how many β lie strictly between b and b+k.
            let between = mask.between(b + 1, nb);
            let m = mask.without_bit(b).with_bit(nb);
            let slot = next.entry(m).or_insert_with(C::zero);
            if between.count_ones().is_multiple_of(2) {
                slot.add_assign(c);
            } else {
                slot.sub_assign(c);
            }
        }
    }
    next
}

/// Bits high-to-low are β₀ > β₁ > …, and λ_i = β_i − (l−1−i).
// As `character::mask_to_partition`: a β-mask bit minus an offset below 32.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
fn mask_to_partition<M: Beta>(mask: M, l: usize) -> Partition {
    let mut parts = Vec::with_capacity(l);
    let mut rest = mask;
    let mut i = 0usize;
    while !rest.is_zero() {
        let b = rest.highest();
        rest = rest.without_bit(b as u32);
        let part = b - (l - 1 - i);
        if part > 0 {
            parts.push(part as u32);
        }
        i += 1;
    }
    // Bits walk high to low and β is strictly decreasing, so `parts` comes out
    // weakly decreasing with the zeros already dropped -- `Partition::new`
    // would sort a sorted vector, which the profile charged for.
    Partition::from_sorted(parts)
}

impl<C: QAlgebra> FromSchur<C> for PowerSum<C> {
    fn from_schur(s: &Schur<C>) -> Self {
        // s_λ = Σ_μ z_μ⁻¹ χ^λ(μ) p_μ  (needs division by z_μ).
        let mut out = PowerSum::zero();
        for (lambda, c) in s.terms() {
            for mu in partitions_cached(lambda.size()).iter() {
                let chi = character_in::<C>(lambda, mu);
                if !chi.is_zero() {
                    // Divisions by an integer — never by a ring element. That
                    // is exactly the `QAlgebra` contract, and why this is not
                    // bounded on `Field`. `div_by_z` divides by z_μ's factors
                    // one at a time rather than forming z_μ, which would cap
                    // this conversion at degree 34 (z_{1^35} leaves `u128`).
                    out.add_term(mu.clone(), mu.div_by_z(&c.mul(&chi)));
                }
            }
        }
        out
    }

    fn from_basis(src: &'static str, terms: &BTreeMap<Partition, C>) -> Option<Self> {
        match src {
            "h" => Some(multiplicative_route(terms, |n| {
                multiplicative_in_power::<C>(n, false)
            })),
            "e" => Some(multiplicative_route(terms, |n| {
                multiplicative_in_power::<C>(n, true)
            })),
            _ => None,
        }
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
                // slow but has no degree wall.
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
        Monomial::from_terms(self.terms().clone())
            .to_schur()
            .omega()
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
/// m_μ = Σ_α x^α (over the distinct rearrangements α of μ) just shifts
/// exponents:
///
/// ```text
///   m_μ · a_δ = Σ_α a_{α+δ},   so   m_μ = Σ_α ± s_{sort(α+δ) − δ}
/// ```
///
/// with `a_β = 0` when β repeats and `± a_{sorted β}` otherwise. Taking n = |μ|
/// slots loses nothing: (K⁻¹)_{μλ} ≠ 0 forces μ ⊵ λ, so ℓ(λ) ≤ |μ|.
///
/// Working in β-numbers β = α + δ, each part of μ is *added to a distinct slot*
/// of the initial set {0, 1, …, l−1}. That makes this the same machinery as
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
/// needed O(p(n)²) Kostka numbers to read p(n) of them — the largest deficit
/// the library carried, and one no constant-factor work could close
/// (`docs/record/transitions.md`).
// `l = |μ|`, a `u32` degree, and the recursion's `v` walks down from it.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
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
    // Accumulate on the β-mask, not on a `Partition` — the same change
    // `p_expand_shared`'s caller makes, for the same reason and in the same
    // file. Every surviving placement reaches the leaf, and keying by partition
    // built, heap-allocated and SipHashed a fresh `Vec<u32>` there once per
    // leaf to produce a few hundred distinct terms. The mask is already the
    // recursion's state and is one word; the partitions get built once at the
    // end, p(n) of them rather than one per leaf.
    let mut acc: Map<u64, i128> = Map::default();
    muir_rec(
        l as i32 - 1,
        (1u64 << l) - 1,
        &mut avail,
        mu.len(),
        1,
        &mut acc,
    );
    Some(
        acc.into_iter()
            .filter(|&(_, v)| v != 0)
            .map(|(mask, v)| (mask_to_partition(mask, l), v))
            .collect(),
    )
}

// β-set slot indices, all bounded by 64 (the mask width). The signs and
// coefficients in `acc` are `i128` and never cast.
#[allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]
fn muir_rec(
    v: i32,
    mask: u64,
    avail: &mut [(u32, u32)],
    left: usize,
    sign: i128,
    acc: &mut Map<u64, i128>,
) {
    if v < 0 {
        if left == 0 {
            // The mask *is* the β-set; `mask_to_partition` reads λ_i = β_i −
            // (l − i) off it once per distinct term, not once per leaf.
            *acc.entry(mask).or_insert(0) += sign;
        }
        return;
    }
    // Every remaining part needs its own slot, and v+1 slots remain.
    if left > v as usize + 1 {
        return;
    }
    // Leave slot v where it is. Safe without a collision test: every already
    // final value is > v.
    muir_rec(v - 1, mask, avail, left, sign, acc);

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
        let s = if between.count_ones().is_multiple_of(2) {
            sign
        } else {
            -sign
        };
        avail[i].1 -= 1;
        muir_rec(
            v - 1,
            (mask & !(1 << v)) | (1 << nb),
            avail,
            left - 1,
            s,
            acc,
        );
        avail[i].1 += 1;
    }
}

/// Row μ of the inverse Kostka matrix, indexed against `parts`
/// (decreasing-lex).
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
/// whole matrix and inverting it computes p(n)² Kostka numbers to read p(n) of
/// them, and the matrix is then nearly all of the conversion's cost
/// (`docs/record/transitions.md`).
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
            // The one narrowing here that carries a *value* rather than an
            // index: `kostka` answers in `u128`, and past `i128::MAX` an `as`
            // would hand back a negative Kostka number, which the alternating
            // sum below would absorb without a trace (R5, and R8's rule for the
            // same constants at the `Ring` seam).
            let k = i128::try_from(kostka(&parts[m], &parts[jj])).unwrap_or_else(|_| {
                panic!(
                    "K_{{{},{}}} does not fit i128; this row solve is \
                     fixed-width whatever the coefficient ring",
                    parts[m], parts[jj]
                )
            });
            acc += w[m] * k;
        }
        w[jj] = -acc;
    }
    w
}

// --- Monomial multiplication (routed through Schur) -------------------------

impl<C: Ring> Monomial<C> {
    /// [`Monomial::mul`](crate::Monomial::mul) by the long way round: expand
    /// into Schur, multiply there, contract back.
    ///
    /// Kept as the **reference oracle**, the same role `NaiveLr` plays for the
    /// Littlewood–Richardson backends, and not as the default: `m → s` inverts
    /// the Kostka matrix, so this costs the whole degree — every partition of
    /// `|μ| + |ν|` participates — where the direct rule costs the answer. It
    /// also passes through `i128` in the inverse Kostka row solve, which the
    /// direct rule never needs, since the monomial structure constants are
    /// counts and no intermediate is larger than the result.
    ///
    /// `monomial_product_agrees_with_the_schur_route` in
    /// [`crate::sym`] is the agreement test.
    ///
    /// # Panics
    ///
    /// Panics when a structure constant does not fit the coefficient ring
    /// `C`. The constants are LR coefficients, Muir's-rule counts, and Kostka
    /// numbers.
    ///
    /// Panics when a Kostka number needed by the inverse Kostka row solve does
    /// not fit `i128`. That solve runs only on input terms of size above 32.
    /// It runs at `i128` whatever `C` is.
    pub fn mul_via_schur(&self, other: &Self) -> Self {
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
    /// implementation *is* "apply ω", so a test phrased in terms of ω would
    /// only restate it. These come from the two edge cases of the monomial
    /// basis, m_{(n)} = p_n and m_{(1^n)} = e_n, pushed through ω(p_n) =
    /// (−1)^{n−1} p_n and ω(e_n) = h_n — facts about the *other* bases.
    #[test]
    fn forgotten_endpoints_match_hand_computation() {
        // The endpoints are one row and one column of n >= 1 cells, and the sign
        // (-1)^{n-1} has no n = 0 case; the empty shape is covered by the
        // round-trip and duality tests, which sweep from 0.
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

    /// {f_λ} is dual to {e_λ} under the Hall inner product: ⟨f_λ, e_μ⟩ =
    /// δ_{λμ}.
    ///
    /// This is the structural characterization of the forgotten basis, and it
    /// reaches it through code the conversion never touches — `hall` and the e
    /// → s expansion. If `Forgotten` were wired to the wrong involution, or to
    /// conjugation on the *index* rather than on the Schur expansion, the
    /// duality would fail while a round trip still closed.
    #[test]
    fn forgotten_is_dual_to_elementary() {
        for n in 0..=7u32 {
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
        for n in 0..=7u32 {
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

    /// Every direct transition among h, e and p gives what the Schur hub gives,
    /// on every index of every degree through 9.
    ///
    /// The two routes share no step. The hub side is Pieri products, a
    /// Jacobi–Trudi determinant and — for p — a table of `S_n` characters
    /// divided by z_μ; the direct side is a linear recursion over index
    /// concatenation. A sign convention that were self-consistently wrong on
    /// one side would still fail here, which a round trip through the same rule
    /// would not.
    #[test]
    fn the_direct_multiplicative_routes_agree_with_the_hub() {
        for n in 0..=9u32 {
            for lambda in partitions_cached(n).iter() {
                let h: Homogeneous<i64> = Homogeneous::monomial(lambda.clone(), 1);
                let e: Elementary<i64> = Elementary::monomial(lambda.clone(), 1);
                let p: PowerSum<i64> = PowerSum::monomial(lambda.clone(), 1);

                assert_eq!(
                    convert::<i64, _, Elementary<i64>>(&h),
                    Elementary::from_schur(&h.to_schur()),
                    "h_{lambda} -> e"
                );
                assert_eq!(
                    convert::<i64, _, Homogeneous<i64>>(&e),
                    Homogeneous::from_schur(&e.to_schur()),
                    "e_{lambda} -> h"
                );
                assert_eq!(
                    convert::<i64, _, Homogeneous<i64>>(&p),
                    Homogeneous::from_schur(&p.to_schur()),
                    "p_{lambda} -> h"
                );
                assert_eq!(
                    convert::<i64, _, Elementary<i64>>(&p),
                    Elementary::from_schur(&p.to_schur()),
                    "p_{lambda} -> e"
                );

                // The two directions into p divide, so they need a ring
                // containing ℚ on both sides of the comparison.
                let hq: Homogeneous<Rational> =
                    Homogeneous::monomial(lambda.clone(), Rational::one());
                let eq: Elementary<Rational> =
                    Elementary::monomial(lambda.clone(), Rational::one());
                assert_eq!(
                    convert::<Rational, _, PowerSum<Rational>>(&hq),
                    PowerSum::from_schur(&hq.to_schur()),
                    "h_{lambda} -> p"
                );
                assert_eq!(
                    convert::<Rational, _, PowerSum<Rational>>(&eq),
                    PowerSum::from_schur(&eq.to_schur()),
                    "e_{lambda} -> p"
                );
            }
        }
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
    /// Batching groups terms by shared prefix and continues one layer per
    /// group; an off-by-one in that grouping would silently attribute a term to
    /// the wrong μ, which no round-trip test would catch.
    #[test]
    fn batched_p_expansion_matches_per_term() {
        for n in 0..=10u32 {
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

    /// The two mask widths are one algorithm, so they must not be two answers.
    ///
    /// Only degrees at most [`MASK_LIMIT`] can be run both ways, and those are
    /// exactly the degrees the wide path never sees in production — which is
    /// the point. Without this, `u128` would be exercised only where nothing
    /// else can check it, and a transcription slip in the wider `Beta` impl
    /// would surface as a wrong plethysm at a degree no oracle reaches.
    #[test]
    fn the_two_mask_widths_agree_where_both_apply() {
        // Small degrees only. Every partition of n is swept both ways here, so
        // the cost climbs with p(n), and what is being checked is a bit
        // transcription rather than anything that appears late: the widths
        // either agree on the whole layer algebra or they disagree at once.
        for n in 1..=14u32 {
            let mus = partitions_cached(n);
            let items: Vec<(&Partition, &i128)> = mus.iter().map(|mu| (mu, &1i128)).collect();

            let mut narrow: Schur<i128> = Schur::zero();
            sweep_degree::<u64, i128>(&items, n as usize, &mut narrow);
            let mut wide: Schur<i128> = Schur::zero();
            sweep_degree::<u128, i128>(&items, n as usize, &mut wide);

            assert_eq!(narrow, wide, "u64 and u128 sweeps disagree at degree {n}");
        }
    }

    /// p → s past the old 32 wall, against a closed form rather than a route.
    ///
    /// `p_n = Σ_{r<n} (−1)^r s_{(n−r, 1^r)}` — the alternating sum of hooks —
    /// is exact, independent of everything here, and cheap on both sides: a
    /// one-part μ is a single Murnaghan–Nakayama step from the root, so the
    /// layer never grows. A round trip through s → p would test the same wide
    /// sweep but costs p(n) terms to set up, which is minutes at these degrees
    /// in a debug build and does not belong in a suite that runs in seconds.
    ///
    /// Degrees 33 and 55 are the two that matter: the first is the one past
    /// the narrow mask, the second is [`WIDE_MASK_LIMIT`] itself, where the
    /// `i128` bound in its doc is tightest.
    #[test]
    fn the_wide_mask_expands_a_power_sum_to_its_hooks() {
        for n in [33u32, 40, WIDE_MASK_LIMIT as u32] {
            let p: PowerSum<i128> = PowerSum::monomial(Partition::new([n]), 1);

            let mut want: Schur<i128> = Schur::zero();
            for r in 0..n {
                let mut parts = vec![n - r];
                parts.extend(std::iter::repeat_n(1, r as usize));
                want.add_term(
                    Partition::from_sorted(parts),
                    if r % 2 == 0 { 1 } else { -1 },
                );
            }
            assert_eq!(p.to_schur(), want, "p_{n} → s is not the hook sum");
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
        for n in 0..=10u32 {
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
        for n in 0..=8u32 {
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
    /// the solve was the shipped implementation, so this is a real oracle
    /// rather than a self-consistency check.
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

    /// The β-mask of `lambda` in `l` slots: β_i = λ_i + (l − 1 − i).
    fn beta_mask(lambda: &Partition, l: usize) -> u64 {
        let mut mask = 0u64;
        for i in 0..l {
            mask |= 1 << (lambda.part(i) as usize + (l - 1 - i));
        }
        mask
    }

    fn strips_via_masks(lambda: &Partition, k: u32, vertical: bool) -> Vec<Partition> {
        let l = lambda.size() as usize + k as usize;
        let mut masks = Vec::new();
        strip_masks(beta_mask(lambda, l), k, vertical, &mut masks);
        let mut out: Vec<Partition> = masks.into_iter().map(|m| mask_to_partition(m, l)).collect();
        out.sort();
        out
    }

    /// The β-mask strip enumeration must agree with the partition-keyed one on
    /// every shape and every strip size through degree 10.
    ///
    /// The two are independent readings of the Pieri rule — one an interlacing
    /// on parts, the other on β-numbers — and the partition form is the one the
    /// Sage oracle fixtures already validate, so this is a real oracle rather
    /// than a self-consistency check. It also has to be, because the horizontal
    /// and vertical conditions are *different* constraints on the β-set and
    /// getting one of them wrong yields plausible partitions rather than
    /// errors.
    #[test]
    fn beta_mask_strips_agree_with_the_partition_enumeration() {
        for n in 0..=10u32 {
            for lambda in partitions_cached(n).iter() {
                for k in 0..=6u32 {
                    let mut want = Vec::new();
                    horizontal_strips(lambda.parts(), k, &mut want);
                    want.sort();
                    assert_eq!(
                        strips_via_masks(lambda, k, false),
                        want,
                        "horizontal {k}-strips on {lambda}"
                    );

                    let mut want = Vec::new();
                    vertical_strips(lambda.parts(), k, &mut want);
                    want.sort();
                    assert_eq!(
                        strips_via_masks(lambda, k, true),
                        want,
                        "vertical {k}-strips on {lambda}"
                    );
                }
            }
        }
    }

    /// A shape where the two strip types give *different* answers, worked by
    /// hand — the value that pins which is which.
    ///
    /// From λ = (2,1), the 2-strips that separate the two rules are (4,1) and
    /// (2,1,1,1). Horizontal means no two new cells in one **column**, vertical
    /// no two in one **row**:
    ///
    /// * (4,1)/λ puts both cells in row 0, columns 2 and 3 — horizontal only.
    /// * (2,1,1,1)/λ puts both in column 0, rows 2 and 3 — vertical only.
    /// * (3,2)/λ puts one at (0,2) and one at (1,1) — both rules accept it, and
    ///   testing it alone would not distinguish them.
    #[test]
    fn horizontal_and_vertical_strips_are_not_the_same_rule() {
        let lambda = part(&[2, 1]);
        let h = strips_via_masks(&lambda, 2, false);
        let v = strips_via_masks(&lambda, 2, true);
        assert!(h.contains(&part(&[4, 1])), "horizontal: two cells in row 0");
        assert!(!v.contains(&part(&[4, 1])), "vertical: never two in a row");
        assert!(
            v.contains(&part(&[2, 1, 1, 1])),
            "vertical: two cells in column 0"
        );
        assert!(
            !h.contains(&part(&[2, 1, 1, 1])),
            "horizontal: never two in a column"
        );
        assert!(
            h.contains(&part(&[3, 2])) && v.contains(&part(&[3, 2])),
            "both rules accept (3,2)"
        );
    }

    /// s_λ · h_k and s_λ · e_k through the mask layer must equal the general
    /// Littlewood–Richardson product, which is what this path replaced.
    #[test]
    fn pieri_steps_agree_with_the_lr_product() {
        for n in 0..=7u32 {
            for lambda in partitions_cached(n).iter() {
                for k in 1..=4u32 {
                    let s: Schur<i128> = Schur::monomial(lambda.clone(), 1);
                    for vertical in [false, true] {
                        let gen = if vertical {
                            Partition::new(std::iter::repeat_n(1, k as usize))
                        } else {
                            Partition::new([k])
                        };
                        let want = s.mul(&Schur::monomial(gen, 1));
                        let mut got: Schur<i128> = Schur::zero();
                        for mu in strips_via_masks(lambda, k, vertical) {
                            got.add_term(mu, 1);
                        }
                        assert_eq!(
                            got,
                            want,
                            "s_{lambda} * {}_{k}",
                            if vertical { "e" } else { "h" }
                        );
                    }
                }
            }
        }
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
            let s: Schur<Rational> = Schur::monomial(part(parts), Rational::one());
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
