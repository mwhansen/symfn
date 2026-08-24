//! Macdonald's `P_λ(x; q, t)`, by the branching formula.
//!
//! ```text
//!   P_λ(x; q, t) = Σ_T ψ_T(q, t) x^T
//! ```
//!
//! over semistandard tableaux `T` of shape λ. A tableau of content μ is the
//! same thing as a chain `∅ = ν⁰ ⊂ ν¹ ⊂ … ⊂ ν^k = λ` of horizontal strips with
//! `|νⁱ/νⁱ⁻¹| = μ_i`, and `ψ_T` is the product of a factor per strip — so the
//! coefficient of `m_μ` is a sum over exactly the chains
//! [`charge`](mod@crate::charge) already enumerates for Kostka–Foulkes. The two
//! share that chain enumeration for that reason.
//!
//! Macdonald, *Symmetric Functions and Hall Polynomials*, 2nd ed., Chapter VI,
//! (6.24) and (7.13'). Symmetrica has no Macdonald polynomials at all, so
//! unlike Hall–Littlewood there is no C implementation to compare against —
//! Sage is the only external oracle here. Sage's equivalents are
//! `Sym.macdonald().P()`, `.Q()` and `.J()`, which
//! `scripts/check_macdonald.py` compares against.
//!
//! ## Range
//!
//! Over a fixed-width `C` this family refuses rather than wrapping past its
//! wall (`docs/policies/failure.md`, R3). Through the Python boundary all
//! three of `P`, `Q` and `J` escalate — the fixed-width pass reports and the
//! same generic code re-runs over `BigInt` — so the walls below are what a
//! *Rust* caller at `C = i128` meets.
//!
//! **The extremal shape is the single row `λ = (n)`**, at every degree
//! measured, and the walls there are reachable: `P` gives out at n = 30
//! (122 bits at n = 29), `Q` and `J` at n = 26. The jump from 101 bits to
//! overflow in one degree says where the wall is — in the **intermediates**
//! of the `Frac` arithmetic rather than in the answers, which is the shape
//! this crate keeps meeting.
//!
//! Shape dominates degree here as it does for Hall–Littlewood, but in the
//! opposite direction: at λ = 1ⁿ the same `J` gains ~0.3 bits per degree and
//! is still 23 bits at n = 96, so 1ⁿ is nowhere near extremal for `J` even
//! though it is exactly extremal for `Q'`.
//!
//! Denominators carry no such number: `Frac` keeps them factored as a multiset
//! of binomials `1 − qᵃtᵇ` and never expands one.
//!
//! Measurements and the harness are in
//! `docs/record/failure-and-overflow.md` (`examples/probe_qt_walls.rs`).
//!
//! ## The coefficients
//!
//! For a cell `s` of λ with arm `a` and leg `l`,
//!
//! ```text
//!   b_λ(s) = (1 − q^a t^{l+1}) / (1 − q^{a+1} t^l)
//! ```
//!
//! and for a horizontal strip λ/μ,
//!
//! ```text
//!   ψ_{λ/μ} = ∏ b_μ(s) / b_λ(s)
//! ```
//!
//! over the cells `s` of μ lying in a row that meets λ/μ but **not** in a
//! column that meets it. Both conditions matter and they are easy to swap; the
//! tests pin them on λ = (2), where `ψ` is `(1−t)(1+q)/(1−qt)` and getting
//! either condition wrong gives something else.
//!
//! Every factor is a ratio of binomials `1 − qᵃtᵇ`, which is the entire reason
//! [`Frac`] can avoid a general gcd — see its module docs. Neither exponent
//! pair can be `(0,0)`: the numerator of `b` has `t`-exponent `l+1 ≥ 1` and the
//! denominator has `q`-exponent `a+1 ≥ 1`.

// Shape indices; coefficients are `Frac<C>`.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use std::collections::{BTreeMap, HashMap};

use crate::coeff::Ring;
use crate::frac::Frac;
use crate::partition::Partition;
use crate::sym::{add_at, by_degree, Monomial, SymFn};

/// `P_λ(x; q, t)` in the monomial basis.
///
/// Monic and triangular: the coefficient of `m_λ` is 1 and every other `m_μ`
/// that appears has μ strictly below λ in dominance order. At λ = ∅ the
/// expansion is the single term `m_∅` with coefficient 1.
pub fn macdonald_p<C: Ring>(lambda: &Partition) -> Monomial<Frac<C>> {
    p_with(lambda, &mut PsiCache::new())
}

/// `P_λ(x; q, t)` in the monomial basis, for every λ ⊢ n, in the order of
/// [`partitions_of`](crate::partitions_of).
///
/// One [`macdonald_p`] per shape, with the ψ cache shared across them: the
/// strips a tableau of shape λ is built from are the horizontal strips between
/// shapes contained in λ, and shapes of one degree contain many of the same
/// smaller ones. The sharing saves little — `docs/record/macdonald.md` has the
/// measurement — so this exists as the whole-degree unit rather than as a
/// faster route to it.
///
/// That unit is what [`monomial_to_macdonald_p`] needs: it back-substitutes
/// through every dominance-smaller `P`, so one shape costs what the table
/// costs.
pub fn macdonald_p_table<C: Ring>(n: u32) -> Vec<(Partition, Monomial<Frac<C>>)> {
    let mut cache = PsiCache::new();
    crate::partitions_of(n)
        .into_iter()
        .map(|lambda| {
            let f = p_with(&lambda, &mut cache);
            (lambda, f)
        })
        .collect()
}

/// ψ of a strip depends only on the two shapes, and the same strip recurs
/// across chains, across contents and across shapes of the same degree — a
/// tall shape has few distinct strips and very many tableaux using them.
type PsiCache = HashMap<(Vec<u32>, Vec<u32>), Factors>;

/// [`macdonald_p`] against a caller-owned ψ cache.
fn p_with<C: Ring>(lambda: &Partition, cache: &mut PsiCache) -> Monomial<Frac<C>> {
    let mut out = Monomial::zero();
    if lambda.is_empty() {
        out.add_term(Partition::new([]), <Frac<C> as Ring>::one());
        return out;
    }
    let target: Vec<u32> = lambda.parts().to_vec();
    let mut acc: Factors = BTreeMap::new();

    for mu in crate::partitions_of(lambda.size()) {
        crate::interrupt::poll();
        let mut total = <Frac<C> as Ring>::zero();
        let mut chain: Vec<Vec<u32>> = vec![vec![0; target.len()]];
        crate::charge::build(&target, mu.parts(), 0, &mut chain, &mut |ch| {
            // Accumulate the whole tableau's exponents before building anything.
            acc.clear();
            for k in 1..ch.len() {
                let key = (ch[k].clone(), ch[k - 1].clone());
                let f = cache
                    .entry(key)
                    .or_insert_with(|| psi_factors(&ch[k], &ch[k - 1]));
                for (&e, &m) in f.iter() {
                    *acc.entry(e).or_insert(0) += m;
                }
            }
            acc.retain(|_, m| *m != 0);
            total.add_assign(&Frac::from_factors(&acc));
        });
        if !total.is_zero() {
            // The one reduction for this coefficient — see `Frac::add_assign`.
            total.reduce();
            out.add_term(mu, total);
        }
    }
    out
}

/// `Q_λ(x; q, t) = b_λ(q,t) · P_λ`, the other normalization of the same basis.
///
/// `b_λ = ∏_{s∈λ} b_λ(s)`, so this is `P` scaled by a single ratio of binomial
/// products — no re-enumeration.
pub fn macdonald_q<C: Ring>(lambda: &Partition) -> Monomial<Frac<C>> {
    scale(macdonald_p(lambda), &b_factors(lambda.parts()))
}

/// `J_λ(x; q, t) = c_λ(q,t) · P_λ`, the integral form.
///
/// `c_λ = ∏_{s∈λ} (1 − q^{a(s)} t^{l(s)+1})` — a *polynomial*, and exactly the
/// numerator of `b_λ`, which is what clears `P`'s denominators and makes `J`
/// the form with coefficients in `ℤ[q,t]`.
pub fn macdonald_j<C: Ring>(lambda: &Partition) -> Monomial<Frac<C>> {
    scale(macdonald_p(lambda), &c_factors(lambda.parts()))
}

/// `f`, given in the monomial basis, rewritten in the Macdonald `P` basis: the
/// `c_λ` of `f = Σ_λ c_λ P_λ(x; q, t)`.
///
/// `P` is monic and dominance-unitriangular in the monomial basis, so
/// `m_λ = P_λ − Σ_{μ ◁ λ} c_{λμ} m_μ` solves downward and the whole `m → P`
/// transition is a back-substitution through [`macdonald_p_table`]. Nothing
/// here is a new enumeration; the coefficients are the same ones `P → m`
/// carries, resolved the other way.
///
/// `f` may mix degrees — each degree's table is applied to its own terms — and
/// the zero element gives the empty map. The result is keyed by partition in
/// the element order and holds no zeros. It is a plain map because the crate
/// has no `P`-basis type, and a [`Monomial`] holding `P`-coefficients would be
/// the confusion the basis types exist to prevent.
///
/// Costs one [`macdonald_p_table`] per degree present in `f`. Sage's
/// equivalent is `Sym.macdonald().P()(f)`.
///
/// ```
/// use symfn::{monomial_to_macdonald_p, Frac, Monomial, Partition, Rational, Ring, SymFn};
///
/// let one = <Frac<Rational> as Ring>::one();
/// let m2: Monomial<Frac<Rational>> = Monomial::monomial(Partition::new([2]), one.clone());
/// let in_p = monomial_to_macdonald_p(&m2);
///
/// assert_eq!(in_p[&Partition::new([2])], one);
/// let (num, den) = in_p[&Partition::new([1, 1])].parts();
/// assert_eq!(den.collect::<Vec<_>>(), vec![(&(1, 1), &1)]);
/// let r = Rational::from_int;
/// // −(1 − t)(1 + q), expanded
/// let want = [((0, 0), r(-1)), ((0, 1), r(1)), ((1, 0), r(-1)), ((1, 1), r(1))];
/// assert!(num.terms().map(|(&e, c)| (e, c.clone())).eq(want));
/// ```
///
/// So `m_2 = P_2 − [(1−t)(1+q)/(1−qt)] P_11`, which is `P_2` read backwards:
/// the coefficient `P → m` puts on the dominance-smaller shape comes back
/// negated, and the shape it sits on is `(1,1)` rather than `(2)` because the
/// triangularity runs the other way. Swapping `q` and `t` would give
/// `(1−q)(1+t)/(1−qt)` instead, which is the twist to check.
pub fn monomial_to_macdonald_p<C: Ring + Send + Sync + 'static>(
    f: &Monomial<Frac<C>>,
) -> BTreeMap<Partition, Frac<C>> {
    let mut out: BTreeMap<Partition, Frac<C>> = BTreeMap::new();
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
    // Reduced once, at the end: `Frac::add_assign` deliberately leaves a
    // running sum over the lcm of the denominators, because a cancellation can
    // only be decided when the sum is complete (`src/frac.rs`).
    for v in out.values_mut() {
        v.reduce();
    }
    out
}

/// `f`, given in the monomial basis, rewritten in the Macdonald `Q` basis: the
/// `c_λ` of `f = Σ_λ c_λ Q_λ(x; q, t)`.
///
/// `Q_λ = b_λ P_λ`, so this is [`monomial_to_macdonald_p`] with each
/// coefficient divided by that shape's `b_λ` — a product of binomials, applied
/// factored, so no second solve happens. Same contract: mixed degrees are
/// allowed, the zero element gives the empty map, and the result is in element
/// order with no zeros. Sage's equivalent is `Sym.macdonald().Q()(f)`.
///
/// ```
/// use symfn::{monomial_to_macdonald_q, Frac, Monomial, Partition, Rational, Ring, SymFn};
///
/// let one = <Frac<Rational> as Ring>::one();
/// let m11: Monomial<Frac<Rational>> = Monomial::monomial(Partition::new([1, 1]), one);
/// let in_q = monomial_to_macdonald_q(&m11);
///
/// assert_eq!(in_q.len(), 1);
/// let (num, den) = in_q[&Partition::new([1, 1])].parts();
/// assert_eq!(den.collect::<Vec<_>>(), vec![(&(0, 1), &1), (&(0, 2), &1)]);
/// let r = Rational::from_int;
/// // (1 − q t)(1 − q), expanded
/// let want = [((0, 0), r(1)), ((1, 0), r(-1)), ((1, 1), r(-1)), ((2, 1), r(1))];
/// assert!(num.terms().map(|(&e, c)| (e, c.clone())).eq(want));
/// ```
///
/// So `m_11 = [(1−qt)(1−q)/((1−t)(1−t²))] Q_11`, where
/// [`monomial_to_macdonald_p`] gives `m_11 = P_11` outright — the value that
/// separates the two normalizations at the smallest shape where they differ.
pub fn monomial_to_macdonald_q<C: Ring + Send + Sync + 'static>(
    f: &Monomial<Frac<C>>,
) -> BTreeMap<Partition, Frac<C>> {
    let mut out = monomial_to_macdonald_p(f);
    for (lambda, v) in &mut out {
        let mut over_b = b_factors(lambda.parts());
        for m in over_b.values_mut() {
            *m = -*m;
        }
        *v = v.mul_factors(&over_b);
        v.reduce();
    }
    out
}

/// [`monomial_in_p_table`], memoized per ring and degree.
///
/// The back-substitution is 98% of what an `m → P` call costs at degree 8, and
/// its unit is the degree while the entry point is asked for one element — so
/// without this, sweeping p(n) shapes rebuilds it p(n) times
/// (`docs/record/macdonald.md`). See
/// [`schur_in_j_cached`](crate::memo::schur_in_j_cached) for why the key
/// carries the ring and why the store is conditional.
fn cached_in_p_table<C: Ring + Send + Sync + 'static>(
    n: u32,
) -> std::sync::Arc<Vec<BTreeMap<Partition, Frac<C>>>> {
    crate::memo::mac_p_inverse_cached(n, || monomial_in_p_table::<C>(n))
}

/// The `m → P` transition for degree `n`: entry `j` is `m_{parts[j]}` written
/// in the `P` basis, keyed by partition.
///
/// [`macdonald_p_table`] inverted by back-substitution.
/// [`partitions_of`](crate::partitions_of) is lex-descending and λ ⊵ μ implies
/// λ ≥ μ lexicographically, so the dominance-smallest shape is the *last*
/// index and counting the index down solves every `m_μ` before the sums that
/// need it.
///
/// Unlike the Hall–Littlewood inversion this one divides — the entries are
/// rational functions of `q` and `t` rather than polynomials — so every step
/// reduces.
fn monomial_in_p_table<C: Ring>(n: u32) -> Vec<BTreeMap<Partition, Frac<C>>> {
    let parts = crate::memo::partitions_cached(n);
    let table = macdonald_p_table::<C>(n);
    let mut out: Vec<BTreeMap<Partition, Frac<C>>> = vec![BTreeMap::new(); parts.len()];
    for j in (0..parts.len()).rev() {
        crate::interrupt::poll();
        let mut acc = BTreeMap::new();
        acc.insert(parts[j].clone(), <Frac<C> as Ring>::one());
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

/// `f + g`, both given as coefficients in one of the Macdonald bases.
///
/// Which basis is not asked and does not matter: addition is termwise in
/// whatever basis both are written in, and mixing two of them is the caller's
/// error to avoid. Coefficients are reduced, so the result is the same
/// representation every other entry point returns — which is what makes the
/// Python layer's `==` on a sum meaningful, since that equality is structural.
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::{macdonald_element_add, Frac, Partition, Rational, Ring};
///
/// type F = Frac<Rational>;
/// let one: BTreeMap<Partition, F> =
///     [(Partition::new([2]), <F as Ring>::one())].into_iter().collect();
/// let minus: BTreeMap<Partition, F> =
///     [(Partition::new([2]), <F as Ring>::one().neg())].into_iter().collect();
///
/// assert!(macdonald_element_add(&one, &minus).is_empty());
/// assert_eq!(macdonald_element_add(&one, &BTreeMap::new()), one);
/// ```
///
/// A shape whose coefficients cancel leaves no entry at all, which is the
/// invariant the whole crate keeps: a map holds no explicit zeros.
pub fn macdonald_element_add<C: Ring>(
    f: &BTreeMap<Partition, Frac<C>>,
    g: &BTreeMap<Partition, Frac<C>>,
) -> BTreeMap<Partition, Frac<C>> {
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

/// `c·f`, `f` given as coefficients in one of the Macdonald bases.
///
/// Same basis-blindness as [`macdonald_element_add`], and the same reason for
/// reducing: multiplying by `1 − q·t` when that factor sits in a denominator
/// must cancel it, or the answer prints in a form nothing else produces.
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::{macdonald_element_scale, Frac, Partition, Rational, Ring};
///
/// type F = Frac<Rational>;
/// let f: BTreeMap<Partition, F> =
///     [(Partition::new([2]), F::inv_factor(1, 1))].into_iter().collect();
/// let scaled = macdonald_element_scale(&f, &F::factor(1, 1));
///
/// assert_eq!(scaled[&Partition::new([2])], <F as Ring>::one());
/// ```
///
/// So `(1 − q·t)·[1/(1 − q·t)]` is 1 and not itself over itself.
pub fn macdonald_element_scale<C: Ring>(
    f: &BTreeMap<Partition, Frac<C>>,
    c: &Frac<C>,
) -> BTreeMap<Partition, Frac<C>> {
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
/// reason [`monomial_to_macdonald_p`] gives: a `Frac` cancellation can only be
/// decided when the sum is complete.
fn expand_mac<C: Ring>(
    f: &BTreeMap<Partition, Frac<C>>,
    one: fn(&Partition) -> Monomial<Frac<C>>,
) -> Monomial<Frac<C>> {
    let mut acc: BTreeMap<Partition, Frac<C>> = BTreeMap::new();
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

/// The `P`-basis element `f = Σ_λ c_λ P_λ(x; q, t)`, expanded in the monomial
/// basis.
///
/// The inverse of [`monomial_to_macdonald_p`], and its input is that
/// function's output: a plain map from partition to coefficient, because the
/// crate has no `P`-basis type. Shapes of different degrees may be mixed and
/// the empty map gives zero.
///
/// One [`macdonald_p`] per shape present — not per shape of the degree, which
/// is what makes this the right route for an element with few terms and
/// [`macdonald_p_table`] the right one for a whole degree.
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::{macdonald_p_to_monomial, Frac, Partition, Rational, Ring, SymFn};
///
/// type F = Frac<Rational>;
/// let f: BTreeMap<Partition, F> =
///     [(Partition::new([2]), <F as Ring>::one())].into_iter().collect();
/// let m = macdonald_p_to_monomial(&f);
///
/// assert_eq!(m.coeff(&Partition::new([2])), <F as Ring>::one());
/// let c = m.coeff(&Partition::new([1, 1]));
/// let (num, den) = c.parts();
/// assert_eq!(den.collect::<Vec<_>>(), vec![(&(1, 1), &1)]);
/// let r = Rational::from_int;
/// // (1 − t)(1 + q), expanded
/// let want = [((0, 0), r(1)), ((0, 1), r(-1)), ((1, 0), r(1)), ((1, 1), r(-1))];
/// assert!(num.terms().map(|(&e, c)| (e, c.clone())).eq(want));
/// ```
///
/// So `P_2 = m_2 + [(1−t)(1+q)/(1−q·t)] m_11`. ⚠️ Under `q ↔ t` — the twist
/// most Macdonald conventions differ by — the numerator would be
/// `(1−q)(1+t)`, and at `q = t` both are the same, so a Schur specialization
/// cannot tell them apart.
pub fn macdonald_p_to_monomial<C: Ring>(f: &BTreeMap<Partition, Frac<C>>) -> Monomial<Frac<C>> {
    expand_mac(f, macdonald_p::<C>)
}

/// The `Q`-basis element `f = Σ_λ c_λ Q_λ(x; q, t)`, expanded in the monomial
/// basis.
///
/// The inverse of [`monomial_to_macdonald_q`]; same contract as
/// [`macdonald_p_to_monomial`].
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::{macdonald_q_to_monomial, Frac, Partition, Rational, Ring, SymFn};
///
/// type F = Frac<Rational>;
/// let f: BTreeMap<Partition, F> =
///     [(Partition::new([1, 1]), <F as Ring>::one())].into_iter().collect();
/// let m = macdonald_q_to_monomial(&f);
///
/// let c = m.coeff(&Partition::new([1, 1]));
/// let (num, den) = c.parts();
/// assert_eq!(den.collect::<Vec<_>>(), vec![(&(1, 0), &1), (&(1, 1), &1)]);
/// let r = Rational::from_int;
/// // (1 − t)(1 − t²), expanded
/// let want = [((0, 0), r(1)), ((0, 1), r(-1)), ((0, 2), r(-1)), ((0, 3), r(1))];
/// assert!(num.terms().map(|(&e, c)| (e, c.clone())).eq(want));
/// ```
///
/// So `Q_11 = [(1−t)(1−t²)/((1−q)(1−q·t))] m_11`, where
/// [`macdonald_p_to_monomial`] has `P_11 = m_11` outright — the smallest shape
/// at which the two normalizations differ.
pub fn macdonald_q_to_monomial<C: Ring>(f: &BTreeMap<Partition, Frac<C>>) -> Monomial<Frac<C>> {
    expand_mac(f, macdonald_q::<C>)
}

/// The `J`-basis element `f = Σ_λ c_λ J_λ(x; q, t)`, expanded in the monomial
/// basis.
///
/// The integral form, so the coefficients here are *polynomials* in `q` and
/// `t`; same contract as [`macdonald_p_to_monomial`]. Its own inverse takes
/// the Schur basis rather than this one — see `schur_in_macdonald_j`, on the
/// Python surface — because that is the direction the `J` triangularity runs
/// in.
///
/// ```
/// use std::collections::BTreeMap;
/// use symfn::{macdonald_j_to_monomial, Frac, Partition, Rational, Ring, SymFn};
///
/// type F = Frac<Rational>;
/// let f: BTreeMap<Partition, F> =
///     [(Partition::new([2]), <F as Ring>::one())].into_iter().collect();
/// let m = macdonald_j_to_monomial(&f);
///
/// let c = m.coeff(&Partition::new([1, 1]));
/// let (num, den) = c.parts();
/// assert_eq!(den.count(), 0);
/// let r = Rational::from_int;
/// // (1 + q)(1 − t)², expanded
/// let want = [
///     ((0, 0), r(1)),
///     ((0, 1), r(-2)),
///     ((0, 2), r(1)),
///     ((1, 0), r(1)),
///     ((1, 1), r(-2)),
///     ((1, 2), r(1)),
/// ];
/// assert!(num.terms().map(|(&e, c)| (e, c.clone())).eq(want));
/// ```
///
/// So the `m_11` coefficient of `J_2` is `(1+q)(1−t)²` with no denominator at
/// all, which is what "integral form" means: it is `c_λ = (1−t)(1−q·t)` times
/// the `P` coefficient above, and the `1−q·t` cancels.
pub fn macdonald_j_to_monomial<C: Ring>(f: &BTreeMap<Partition, Frac<C>>) -> Monomial<Frac<C>> {
    expand_mac(f, macdonald_j::<C>)
}

/// Multiply every coefficient by a product of binomial powers.
///
/// The scalar stays *factored* all the way through — see [`Frac::mul_factors`].
fn scale<C: Ring>(f: Monomial<Frac<C>>, by: &Factors) -> Monomial<Frac<C>> {
    let mut out = Monomial::zero();
    for (mu, c) in f.terms() {
        crate::interrupt::poll();
        let mut v = c.mul_factors(by);
        v.reduce();
        out.add_term(mu.clone(), v);
    }
    out
}

/// `b_λ = ∏_{s∈λ} (1 − q^{a} t^{l+1}) / (1 − q^{a+1} t^{l})`.
fn b_factors(lambda: &[u32]) -> Factors {
    let mut f = Factors::new();
    for_each_cell(lambda, &mut |a, l| {
        *f.entry((a, l + 1)).or_insert(0) += 1;
        *f.entry((a + 1, l)).or_insert(0) -= 1;
    });
    f.retain(|_, m| *m != 0);
    f
}

/// `c_λ = ∏_{s∈λ} (1 − q^{a} t^{l+1})` — the **numerator** of `b_λ`.
///
/// Not `(1 − q^{a+1} t^{l})`, which is `c'_λ`, the denominator. Writing that
/// one gives a `J_(2)` whose leading coefficient is `(1−q²)(1−q)` where Sage
/// has `(1−t)(1−qt)`; the polynomiality test below is what caught it.
pub(crate) fn c_factors(lambda: &[u32]) -> Factors {
    let mut f = Factors::new();
    for_each_cell(lambda, &mut |a, l| {
        *f.entry((a, l + 1)).or_insert(0) += 1;
    });
    f.retain(|_, m| *m != 0);
    f
}

/// `c'_λ = ∏_{s∈λ}(1 − q^{a(s)+1} t^{l(s)})` — the **denominator** of `b_λ`,
/// and the scalar Lapointe–Lascoux–Morse call `c_{λ'}(t,q)` in their 3.15.
///
/// Conjugating a shape swaps arms and legs, and swapping `q` with `t` then puts
/// `c_{λ'}(t,q) = ∏_{s∈λ}(1 − q^{a(s)+1} t^{l(s)})` — this. It is what the
/// eigenvector route has to multiply by, since that normalizes the leading
/// coefficient to 1 where `J` has this.
pub(crate) fn c_prime_factors(lambda: &[u32]) -> BTreeMap<(u32, u32), u32> {
    let mut f = BTreeMap::new();
    for_each_cell(lambda, &mut |a, l| {
        *f.entry((a + 1, l)).or_insert(0) += 1;
    });
    f
}

/// Every cell of λ, as its (arm, leg).
fn for_each_cell(lambda: &[u32], visit: &mut impl FnMut(u32, u32)) {
    for i in 0..lambda.len() {
        for j in 0..lambda[i] as usize {
            visit(arm(lambda, i, j), leg(lambda, i, j));
        }
    }
}

/// Binomial exponents with signed multiplicities: `(a, b) ↦ m` is
/// `(1 − qᵃtᵇ)^m`, negative meaning a denominator factor.
type Factors = BTreeMap<(u32, u32), i32>;

/// The factors of `ψ_{λ/μ}` for a horizontal strip, both shapes padded to one
/// length.
///
/// Returned unexpanded. Nothing here multiplies a polynomial: the whole of ψ is
/// a product of ratios of binomials, so it is entirely described by counting
/// them, and the counting lets matching factors cancel before any expansion.
fn psi_factors(lam: &[u32], mu: &[u32]) -> Factors {
    let width = lam.iter().copied().max().unwrap_or(0) as usize;
    // Column j meets the strip when λ has a cell there that μ does not.
    let met: Vec<bool> = (0..width)
        .map(|j| count_above(lam, j) > count_above(mu, j))
        .collect();

    let mut f = Factors::new();
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
            // b_μ(s) = (1 − q^{am} t^{lm+1}) / (1 − q^{am+1} t^{lm})
            *f.entry((am, lm + 1)).or_insert(0) += 1;
            *f.entry((am + 1, lm)).or_insert(0) -= 1;
            // 1 / b_λ(s)
            *f.entry((al + 1, ll)).or_insert(0) += 1;
            *f.entry((al, ll + 1)).or_insert(0) -= 1;
        }
    }
    f.retain(|_, m| *m != 0);
    f
}

/// Cells strictly to the right of `(i, j)` in row `i`.
pub(crate) fn arm(shape: &[u32], i: usize, j: usize) -> u32 {
    shape[i] - j as u32 - 1
}

/// Cells strictly below `(i, j)` in column `j`.
pub(crate) fn leg(shape: &[u32], i: usize, j: usize) -> u32 {
    count_above(shape, j) as u32 - 1 - i as u32
}

/// The conjugate part `shape'_j` — how many rows reach past column `j`.
pub(crate) fn count_above(shape: &[u32], j: usize) -> usize {
    shape.iter().filter(|&&v| v as usize > j).count()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coeff::Rational;
    use crate::qt::QtPoly;

    type F = Frac<Rational>;

    fn part(v: &[u32]) -> Partition {
        Partition::new(v.iter().copied())
    }

    fn r(n: i128) -> Rational {
        Rational::from_int(n)
    }

    /// Hand-computed, and the case that pins both halves of the row/column
    /// condition in ψ.
    ///
    /// P_(2) = m_(2) + [(1−t)(1+q)/(1−qt)] m_(1,1). Sage:
    ///
    /// ```text
    ///   sage: P = SymmetricFunctions(QQ['q','t'].fraction_field()).macdonald().P()
    ///   sage: SymmetricFunctions(...).monomial()(P[2])
    ///   m[1, 1] + ... m[2]
    /// ```
    #[test]
    fn p_of_two_matches_the_hand_computation() {
        let p: Monomial<F> = macdonald_p(&part(&[2]));
        assert_eq!(p.coeff(&part(&[2])), <F as Ring>::one(), "monic in m_λ");

        let c = p.coeff(&part(&[1, 1]));
        // (1 − t)(1 + q) / (1 − q t), built independently of the ψ machinery
        let mut one_plus_q: QtPoly<Rational> = <QtPoly<Rational> as Ring>::one();
        one_plus_q.add_term(1, 0, r(1));
        let want = F::factor(0, 1) // 1 − t
            .mul(&F::from_poly(one_plus_q))
            .mul(&F::inv_factor(1, 1)); // / (1 − q t)
        assert_eq!(c, want, "got {c}, want {want}");
    }

    /// q = t collapses Macdonald to Schur. Checked in the monomial basis
    /// against the library's own s → m, which is tested by other means.
    #[test]
    fn at_q_equals_t_it_is_the_schur_function() {
        use crate::convert::FromSchur;
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = macdonald_p(&lambda);
                let s: Monomial<Rational> =
                    Monomial::from_schur(&crate::Schur::monomial(lambda.clone(), r(1)));
                // q = t = 2 is a generic point: no denominator 1 − qᵃtᵇ vanishes
                // there, unlike q = t = 1.
                for (mu, c) in p.terms() {
                    let got = c.eval(&r(2), &r(2)).expect("no pole at q = t = 2");
                    assert_eq!(got, s.coeff(mu), "P_{lambda} at q=t, term {mu}");
                }
                for (mu, want) in s.terms() {
                    if p.coeff(mu) == <F as Ring>::zero() {
                        assert!(want.is_zero(), "missing term {mu} of s_{lambda}");
                    }
                }
            }
        }
    }

    /// t = 1 collapses P to m_λ: the transition becomes the identity.
    #[test]
    fn at_t_equals_one_it_is_the_monomial() {
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = macdonald_p(&lambda);
                for (mu, c) in p.terms() {
                    let got = c.eval(&r(3), &r(1)).expect("no pole at q = 3, t = 1");
                    assert_eq!(
                        got,
                        Rational::from_int(i128::from(mu == &lambda)),
                        "P_{lambda} at t = 1, term {mu}"
                    );
                }
            }
        }
    }

    /// **q = 0 is Hall–Littlewood P.** The check this whole layer was built to
    /// make possible, and the strongest evidence available for either side:
    /// `macdonald_p` is a branching formula over rational functions, while
    /// `hall_littlewood_p` inverts the Kostka–Foulkes matrix produced by the
    /// Morris recursion. They share no code and no algorithm.
    ///
    /// Compared at several values of t rather than symbolically: `Frac::eval`
    /// substitutes numbers, and setting q = 0 while keeping t formal would need
    /// a separate exact division in `ℤ[t]`. Small values keep the exact
    /// rationals well inside i128 — the symbolic comparison is Sage's job.
    #[test]
    fn at_q_zero_it_is_hall_littlewood_p() {
        use crate::convert::FromSchur;
        use crate::qt::QtPoly;
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let mac: Monomial<F> = macdonald_p(&lambda);
                // Hall–Littlewood P arrives in the Schur basis; move it to m.
                let hl: crate::Schur<QtPoly<Rational>> = crate::hall_littlewood_p(&lambda);
                let hl: Monomial<QtPoly<Rational>> = Monomial::from_schur(&hl);
                for &tv in &[2i128, 3, 5] {
                    let t = r(tv);
                    for (mu, c) in mac.terms() {
                        let got = c.eval(&r(0), &t).expect("q = 0 is not a pole");
                        assert_eq!(
                            got,
                            hl.coeff(mu).eval(&Rational::from_int(0), &t),
                            "P_{lambda} at q=0, t={tv}, term {mu}"
                        );
                    }
                    for (mu, c) in hl.terms() {
                        if mac.coeff(mu) == <F as Ring>::zero() {
                            assert!(
                                c.eval(&Rational::from_int(0), &t).is_zero(),
                                "Hall–Littlewood has {mu} where Macdonald does not"
                            );
                        }
                    }
                }
            }
        }
    }

    /// Q = b_λ · P and J = c_λ · P, so both must agree with P after dividing
    /// the scalar back out — and J must be a *polynomial*, which is the
    /// property J exists to have.
    #[test]
    fn q_and_j_are_scalar_multiples_of_p() {
        for n in 0..=5u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = macdonald_p(&lambda);
                let q: Monomial<F> = macdonald_q(&lambda);
                let j: Monomial<F> = macdonald_j(&lambda);
                // The scalar is read off the leading term, where P is 1.
                let b = q.coeff(&lambda);
                let c = j.coeff(&lambda);
                for (mu, pc) in p.terms() {
                    assert_eq!(q.coeff(mu), pc.mul(&b), "Q at {mu}");
                    assert_eq!(j.coeff(mu), pc.mul(&c), "J at {mu}");
                }
                // J has no denominator: c_λ clears exactly what P carries.
                for (mu, jc) in j.terms() {
                    let (_, den) = jc.parts();
                    assert_eq!(den.count(), 0, "J_{lambda} coefficient at {mu} is {jc}");
                }
            }
        }
    }

    /// The round trip: expanding `P_λ` in the monomial basis and reading it
    /// back gives `P_λ` again, for every shape through degree 8. Both
    /// normalizations, because a solve that dropped `b_λ` would still pass on
    /// `P` alone.
    #[test]
    fn every_macdonald_polynomial_comes_back_as_itself() {
        for n in 0..=8u32 {
            for lambda in crate::partitions_of(n) {
                let unit: BTreeMap<Partition, F> =
                    [(lambda.clone(), <F as Ring>::one())].into_iter().collect();
                let p: Monomial<F> = macdonald_p(&lambda);
                assert_eq!(monomial_to_macdonald_p(&p), unit, "m -> P of P_{lambda}");
                let q: Monomial<F> = macdonald_q(&lambda);
                assert_eq!(monomial_to_macdonald_q(&q), unit, "m -> Q of Q_{lambda}");
            }
        }
    }

    /// **Expanding is linear.** `(f + g)` and `c·f` expanded in the monomial
    /// basis must equal the expansions added and scaled there — which ties
    /// [`macdonald_element_add`] and [`macdonald_element_scale`] to a route
    /// that never touches them, since `Monomial` adds and scales through the
    /// ordinary `Ring` operations.
    #[test]
    fn adding_and_scaling_commute_with_expanding() {
        let a = part(&[2, 1]);
        let b = part(&[1, 1, 1]);
        let one = <F as Ring>::one();
        let f: BTreeMap<Partition, F> = [(a.clone(), one.clone())].into_iter().collect();
        let g: BTreeMap<Partition, F> = [(b.clone(), one.clone())].into_iter().collect();

        let mut want: Monomial<F> = macdonald_p_to_monomial(&f);
        for (mu, c) in macdonald_p_to_monomial(&g).terms() {
            want.add_term(mu.clone(), c.clone());
        }
        assert_eq!(
            macdonald_p_to_monomial(&macdonald_element_add(&f, &g)),
            want,
            "P_{a} + P_{b}"
        );

        let c = F::inv_factor(1, 1);
        let mut scaled: Monomial<F> = Monomial::zero();
        for (mu, v) in macdonald_p_to_monomial(&f).terms() {
            let mut w = v.mul(&c);
            w.reduce();
            scaled.add_term(mu.clone(), w);
        }
        assert_eq!(
            macdonald_p_to_monomial(&macdonald_element_scale(&f, &c)),
            scaled,
            "P_{a}/(1 - q t)"
        );
    }

    /// The round trip the other way: solving `m_μ` into a normalization and
    /// expanding it back gives `m_μ`. Together with
    /// [`every_macdonald_polynomial_comes_back_as_itself`] this pins both
    /// composites, which a table inverted in only one direction would not
    /// survive.
    #[test]
    fn every_monomial_comes_back_as_itself() {
        for n in 0..=7u32 {
            for mu in crate::partitions_of(n) {
                let f: Monomial<F> = Monomial::monomial(mu.clone(), <F as Ring>::one());
                assert_eq!(
                    macdonald_p_to_monomial(&monomial_to_macdonald_p(&f)),
                    f,
                    "P at m_{mu}"
                );
                assert_eq!(
                    macdonald_q_to_monomial(&monomial_to_macdonald_q(&f)),
                    f,
                    "Q at m_{mu}"
                );
            }
        }
    }

    /// `J` has no `m → J` solve to invert — its inverse takes the Schur basis
    /// — so what pins [`macdonald_j_to_monomial`] is `J_λ = c_λ·P_λ`: the
    /// expansion, solved back into the `P` basis, must be a *single* term, at
    /// λ itself. A `J` that expanded through the wrong shape's scalar would
    /// still be triangular, and this would catch it.
    #[test]
    fn expanding_j_gives_a_multiple_of_p() {
        for n in 0..=6u32 {
            for lambda in crate::partitions_of(n) {
                let unit: BTreeMap<Partition, F> =
                    [(lambda.clone(), <F as Ring>::one())].into_iter().collect();
                let in_p = monomial_to_macdonald_p(&macdonald_j_to_monomial(&unit));
                assert_eq!(in_p.len(), 1, "J_{lambda} written in P");
                assert!(in_p.contains_key(&lambda), "J_{lambda} written in P");
            }
        }
    }

    /// The hand values, which is what the round trip cannot give: it is blind
    /// to any error the forward direction shares.
    ///
    /// `P_2 = m_2 + [(1−t)(1+q)/(1−qt)] m_11` and `P_11 = m_11`, so
    /// `m_2 = P_2 − [(1−t)(1+q)/(1−qt)] P_11` and `m_11 = P_11`. Dividing by
    /// `b_11 = (1−t²)(1−t)/((1−qt)(1−q))` gives the `Q` value.
    #[test]
    fn m2_and_m11_in_p_and_q_are_the_hand_values() {
        let one = <F as Ring>::one();
        let m2: Monomial<F> = Monomial::monomial(part(&[2]), one.clone());
        let m11: Monomial<F> = Monomial::monomial(part(&[1, 1]), one.clone());

        // −(1 − t)(1 + q)/(1 − q t)
        let mut one_plus_q: QtPoly<Rational> = <QtPoly<Rational> as Ring>::one();
        one_plus_q.add_term(1, 0, r(1));
        let psi = F::factor(0, 1)
            .mul(&F::from_poly(one_plus_q))
            .mul(&F::inv_factor(1, 1));

        let in_p = monomial_to_macdonald_p(&m2);
        assert_eq!(in_p[&part(&[2])], one, "m_2 is monic in P_2");
        assert_eq!(in_p[&part(&[1, 1])], psi.neg(), "m_2 in P_11");

        assert_eq!(
            monomial_to_macdonald_p(&m11),
            [(part(&[1, 1]), one.clone())].into_iter().collect(),
            "m_11 = P_11"
        );

        // 1/b_11 = (1 − q t)(1 − q) / ((1 − t²)(1 − t))
        let over_b = F::factor(1, 1)
            .mul(&F::factor(1, 0))
            .mul(&F::inv_factor(0, 2))
            .mul(&F::inv_factor(0, 1));
        let in_q = monomial_to_macdonald_q(&m11);
        assert_eq!(in_q.len(), 1, "m_11 reaches Q_11 alone");
        assert_eq!(in_q[&part(&[1, 1])], over_b, "m_11 in Q_11");
    }

    /// The memoized table is the computed one, at each ring separately: the
    /// key carries the ring because the Python boundary escalates, and a
    /// ring-blind cache would hand the wide pass the narrow pass's values.
    #[test]
    fn the_cached_table_is_the_computed_table_for_each_ring() {
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

    /// Both transitions are linear and take an argument that mixes degrees;
    /// the zero element gives the empty map.
    #[test]
    fn the_expansions_are_linear_and_take_mixed_degrees() {
        let two = F::from_poly(QtPoly::term(0, 0, r(2)));
        let three = F::from_poly(QtPoly::term(0, 0, r(3)));
        let mut f: Monomial<F> = Monomial::zero();
        f.add_term(part(&[1]), <F as Ring>::one());
        f.add_term(part(&[2]), two.clone());
        f.add_term(part(&[1, 1]), three.clone());

        for (name, to_basis) in [
            (
                "m -> P",
                monomial_to_macdonald_p as fn(&Monomial<F>) -> BTreeMap<Partition, F>,
            ),
            ("m -> Q", monomial_to_macdonald_q),
        ] {
            let mut want: BTreeMap<Partition, F> = BTreeMap::new();
            for (mu, scalar) in [
                (part(&[1]), <F as Ring>::one()),
                (part(&[2]), two.clone()),
                (part(&[1, 1]), three.clone()),
            ] {
                let piece = to_basis(&Monomial::monomial(mu, <F as Ring>::one()));
                for (lambda, c) in piece {
                    crate::sym::add_at(&mut want, &lambda, c.mul(&scalar));
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

    /// Monic and strictly triangular in dominance order.
    #[test]
    fn expansion_is_unitriangular_in_dominance() {
        for n in 0..=7u32 {
            for lambda in crate::partitions_of(n) {
                let p: Monomial<F> = macdonald_p(&lambda);
                assert_eq!(p.coeff(&lambda), <F as Ring>::one(), "monic at {lambda}");
                for mu in p.terms().keys() {
                    assert!(
                        crate::kostka::kostka(&lambda, mu) > 0,
                        "{mu} appears in P_{lambda} but λ does not dominate it"
                    );
                }
            }
        }
    }
}
