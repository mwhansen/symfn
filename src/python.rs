//! PyO3 bridge (feature = "python"): a **coarse-grained** API importable into Sage.
//!
//! The design rule established up front: cross-language calls have overhead, so
//! every function here does a *whole-object* operation — multiply two complete
//! symmetric functions, convert an entire element between bases — never a
//! per-monomial call. Marshalling happens once at the boundary; the inner loops
//! stay in Rust. Exposing fine-grained accessors would reintroduce exactly the
//! per-call tax this crate exists to escape.
//!
//! Elements cross the boundary as lists of `(partition, coefficient)` pairs,
//! e.g. `[([2,1], 3), ([3], -1)]`, which maps directly onto Sage's
//! `.monomial_coefficients()` dicts.
//!
//! Build with: `maturin develop --features python`
//!
//! ## Coefficients have no ceiling, and are not slower for it
//!
//! Coefficients cross as **Python `int`s of arbitrary size**, in both
//! directions. PyO3's `num-bigint` conversion does that natively, so nothing is
//! encoded as a string and a caller never sees a mixed-type list.
//!
//! Internally each call runs twice at most, and almost always once:
//!
//! 1. over [`Guarded`] — `i128` that *reports* overflow instead of wrapping;
//! 2. if anything overflowed, again over `BigInt`, which cannot.
//!
//! This is [`character_in`](crate::character::character_in)'s pattern applied to
//! every coefficient, and it exists because the alternative was silent
//! corruption: the boundary used to be `i128` and `impl Ring for i128`
//! multiplies with a plain `*`, so a structure constant past ~1.7e38 came back
//! **wrapped, with no signal**. Refusing loudly would have been defensible;
//! returning a wrong number was not.
//!
//! The fast path costs 0–1% against unchecked arithmetic
//! (`examples/bench_guarded.rs`), and escalation is rare — measured coefficient
//! widths in these workloads are 1–2 limbs (`examples/coeff_sizes.rs`).

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::coeff::Ring;
use crate::convert::{FromSchur, ToSchur};
use crate::guard::{guarded, Guarded, GuardedRat};
use crate::hopf::{self, SkewBy};
use crate::lr::{LrBackend, NaiveLr};
use crate::ops;
use crate::partition::Partition;
use crate::sym::{Elementary, Forgotten, Homogeneous, Monomial, PowerSum, Schur, SymFn};

/// A coefficient crossing the boundary.
///
/// Python only ever sees an `int` — this distinction is invisible there. It
/// exists because `BigInt` is heap-allocated and almost every coefficient is
/// small: routing all of them through it cost **7.6%** on a term-heavy pass and
/// **22%** on `coproduct`, which is marshalling-dominated. PyO3 converts `i128`
/// with no allocation, so the common case now pays nothing and only genuinely
/// wide values allocate.
///
/// Note this is *not* a mixed-type list on the Python side. Both arms convert
/// to the same `int`; the enum never escapes Rust.
#[derive(Clone, Debug)]
enum Coeff {
    Small(i128),
    Big(BigInt),
}

impl<'py> IntoPyObject<'py> for Coeff {
    type Target = PyAny;
    type Output = Bound<'py, PyAny>;
    type Error = PyErr;
    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        Ok(match self {
            Coeff::Small(v) => v.into_pyobject(py)?.into_any(),
            Coeff::Big(v) => v.into_pyobject(py)?.into_any(),
        })
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for Coeff {
    type Error = PyErr;
    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, PyErr> {
        // Narrow first: a Python int that fits skips the bignum path entirely.
        match ob.extract::<i128>() {
            Ok(v) => Ok(Coeff::Small(v)),
            Err(_) => Ok(Coeff::Big(ob.extract::<BigInt>()?)),
        }
    }
}

impl Coeff {
    fn to_big(&self) -> BigInt {
        match self {
            Coeff::Small(v) => BigInt::from(*v),
            Coeff::Big(v) => v.clone(),
        }
    }
    fn as_i128(&self) -> Option<i128> {
        match self {
            Coeff::Small(v) => Some(*v),
            Coeff::Big(v) => i128::try_from(v).ok(),
        }
    }
}

type Terms = Vec<(Vec<u32>, Coeff)>;
type RatTerms = Vec<(Vec<u32>, (Coeff, Coeff))>;

fn part(p: &[u32]) -> Partition {
    Partition::new(p.iter().copied())
}

// --- the two coefficient rings the boundary runs over ------------------------

/// A ring that can carry a coefficient across the boundary.
///
/// Implemented by exactly two types, which are the two passes: [`Guarded`] (the
/// fixed-width attempt, which may decline an input that does not fit) and
/// `BigInt` (the fallback, which never declines).
trait Boundary: Ring + Sized {
    fn from_coeff(v: &Coeff) -> Option<Self>;
    fn to_coeff(&self) -> Coeff;
}

impl Boundary for Guarded {
    fn from_coeff(v: &Coeff) -> Option<Self> {
        v.as_i128().map(Guarded)
    }
    /// No allocation: this is the fast path's whole point.
    fn to_coeff(&self) -> Coeff {
        Coeff::Small(self.0)
    }
}

impl Boundary for BigInt {
    fn from_coeff(v: &Coeff) -> Option<Self> {
        Some(v.to_big())
    }
    fn to_coeff(&self) -> Coeff {
        Coeff::Big(self.clone())
    }
}

/// The same, for the rings that carry denominators.
trait BoundaryRat: Ring + Sized {
    fn from_coeff(v: &Coeff) -> Option<Self>;
    /// `(numerator, denominator)`, denominator positive and in lowest terms.
    fn split(&self) -> (Coeff, Coeff);
}

impl BoundaryRat for GuardedRat {
    fn from_coeff(v: &Coeff) -> Option<Self> {
        v.as_i128().map(<GuardedRat as Ring>::from_i128)
    }
    fn split(&self) -> (Coeff, Coeff) {
        (Coeff::Small(self.numer()), Coeff::Small(self.denom()))
    }
}

impl BoundaryRat for BigRational {
    fn from_coeff(v: &Coeff) -> Option<Self> {
        Some(BigRational::from(v.to_big()))
    }
    fn split(&self) -> (Coeff, Coeff) {
        (
            Coeff::Big(self.numer().clone()),
            Coeff::Big(self.denom().clone()),
        )
    }
}

fn build<C: Boundary, B: SymFn<C>>(terms: &Terms) -> Option<B> {
    let mut x = B::zero();
    for (p, c) in terms {
        x.add_term(part(p), C::from_coeff(c)?);
    }
    Some(x)
}

fn build_rat<C: BoundaryRat, B: SymFn<C>>(terms: &Terms) -> Option<B> {
    let mut x = B::zero();
    for (p, c) in terms {
        x.add_term(part(p), C::from_coeff(c)?);
    }
    Some(x)
}

fn dump<C: Boundary, S: SymFn<C>>(x: &S) -> Terms {
    x.terms()
        .iter()
        .map(|(p, c)| (p.parts().to_vec(), c.to_coeff()))
        .collect()
}

/// Dump a rational element, refusing any coefficient that is not an integer.
///
/// Plethysm and the internal product take Schur input to Schur output, so a
/// denominator would mean a bug in the power-sum route rather than a
/// representable answer. Reported, never truncated.
fn dump_integral<C: BoundaryRat, S: SymFn<C>>(x: &S, what: &str) -> PyResult<Terms> {
    let mut out = Vec::with_capacity(x.terms().len());
    for (p, c) in x.terms() {
        let (num, den) = c.split();
        if !den.to_big().is_one() {
            let (n, d) = (num.to_big(), den.to_big());
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "non-integral {what} coefficient {n}/{d}"
            )));
        }
        out.push((p.parts().to_vec(), num));
    }
    Ok(out)
}

/// Run the fixed-width attempt; fall back to arbitrary precision.
///
/// `fast` returns `None` for either reason that forces the fallback — an input
/// too wide to load, or an overflow during the computation — and the two are
/// deliberately not distinguished, because the response is the same.
fn escalate<T>(fast: impl FnOnce() -> Option<T>, slow: impl FnOnce() -> T) -> T {
    match fast() {
        Some(v) => v,
        None => slow(),
    }
}

// --- products ---------------------------------------------------------------

/// Multiply two Schur-basis elements (Littlewood–Richardson).
#[pyfunction]
fn schur_multiply(a: Terms, b: Terms) -> Terms {
    escalate(
        || {
            let (x, y): (Schur<Guarded>, Schur<Guarded>) = (build(&a)?, build(&b)?);
            Some(dump(&guarded(|| x.mul(&y))?))
        },
        || {
            let (x, y): (Schur<BigInt>, Schur<BigInt>) = (build(&a).unwrap(), build(&b).unwrap());
            dump(&x.mul(&y))
        },
    )
}

/// Drop every memo cache.
///
/// Exposed for benchmarking rather than for normal use: the caches are
/// referentially transparent, so clearing them cannot change a result, only a
/// timing. A comparison that reuses inputs measures the cache on the second
/// call and not the algorithm — which is exactly the trap
/// `scripts/compare_sage.py` documents on Sage's side.
#[pyfunction]
fn clear_caches() {
    crate::memo::clear_caches();
}

/// A single Littlewood–Richardson coefficient c^λ_{μν}.
///
/// Deliberately [`NaiveLr`] and not [`AutoLr`](crate::strip_lr::AutoLr), which
/// looks backwards — `NaiveLr` is the naive *reference* backend — but is what
/// the measurements say. A targeted backtrack costs roughly the coefficient's
/// own size, while `AutoLr` builds a whole expansion and indexes into it, so
/// for one coefficient:
///
/// ```text
///   c^[16,14,12,10,8,6]_{[8,7,6,5,4,3],[8,7,6,5,4,3]} = 1        5µs vs  771µs
///   c^[24,20,16,12]_{[12,10,8,6],[12,10,8,6]}         = 1        6µs vs  145µs
///   c^[13,12..2]_{[5,4,3,2,1],[12,11..3]}         = 14080     2446µs vs  347µs
/// ```
///
/// So the naive search wins by 10–100x whenever the coefficient is small, which
/// is the overwhelmingly common case, and loses only when it is large — because
/// then it enumerates that many tableaux. Whole *products* are a different
/// question and go through `AutoLr` (see `schur_multiply`).
#[pyfunction]
fn lr_coefficient(lambda: Vec<u32>, mu: Vec<u32>, nu: Vec<u32>) -> u128 {
    NaiveLr.lr_coeff(&part(&lambda), &part(&mu), &part(&nu))
}

// --- conversions out of Schur ----------------------------------------------

/// A conversion `Schur -> $basis`, run fixed-width first and re-run exactly if
/// that overflows.
macro_rules! out_of_schur {
    ($name:ident, $basis:ident) => {
        #[pyfunction]
        fn $name(a: Terms) -> Terms {
            escalate(
                || {
                    let s: Schur<Guarded> = build(&a)?;
                    Some(dump(&guarded(|| $basis::<Guarded>::from_schur(&s))?))
                },
                || {
                    let s: Schur<BigInt> = build(&a).unwrap();
                    dump(&$basis::<BigInt>::from_schur(&s))
                },
            )
        }
    };
}

out_of_schur!(schur_to_homogeneous, Homogeneous);
out_of_schur!(schur_to_elementary, Elementary);
out_of_schur!(schur_to_monomial, Monomial);
out_of_schur!(schur_to_forgotten, Forgotten);

/// s → p. Coefficients are rational, returned as `(numerator, denominator)`.
#[pyfunction]
fn schur_to_power(a: Terms) -> RatTerms {
    fn split<C: BoundaryRat>(p: &PowerSum<C>) -> RatTerms {
        p.terms()
            .iter()
            .map(|(part, c)| (part.parts().to_vec(), c.split()))
            .collect()
    }
    escalate(
        || {
            let s: Schur<GuardedRat> = build_rat(&a)?;
            Some(split(&guarded(|| PowerSum::from_schur(&s))?))
        },
        || {
            let s: Schur<BigRational> = build_rat(&a).unwrap();
            let p: PowerSum<BigRational> = PowerSum::from_schur(&s);
            split(&p)
        },
    )
}

// --- conversions into Schur -------------------------------------------------

macro_rules! into_schur {
    ($name:ident, $basis:ident) => {
        #[pyfunction]
        fn $name(a: Terms) -> Terms {
            escalate(
                || {
                    let x: $basis<Guarded> = build(&a)?;
                    Some(dump(&guarded(|| x.to_schur())?))
                },
                || {
                    let x: $basis<BigInt> = build(&a).unwrap();
                    dump(&x.to_schur())
                },
            )
        }
    };
}

into_schur!(homogeneous_to_schur, Homogeneous);
into_schur!(elementary_to_schur, Elementary);
into_schur!(monomial_to_schur, Monomial);
into_schur!(power_to_schur, PowerSum);
into_schur!(forgotten_to_schur, Forgotten);

/// Plethysm f[g] of two Schur-basis elements.
///
/// Computed through the power-sum basis (see `crate::plethysm`).
#[pyfunction]
fn plethysm(f: Terms, g: Terms) -> PyResult<Terms> {
    escalate(
        || {
            let (x, y): (Schur<GuardedRat>, Schur<GuardedRat>) = (build_rat(&f)?, build_rat(&g)?);
            let r = guarded(|| crate::plethysm::plethysm(&x, &y))?;
            Some(dump_integral(&r, "plethysm"))
        },
        || {
            let (x, y): (Schur<BigRational>, Schur<BigRational>) =
                (build_rat(&f).unwrap(), build_rat(&g).unwrap());
            dump_integral(&crate::plethysm::plethysm(&x, &y), "plethysm")
        },
    )
}

// --- classical quantities ---------------------------------------------------

/// Skew a Schur-basis element by `g`, given in `basis` — the adjoint of
/// multiplication by g under the Hall inner product.
///
/// `basis` selects which rule runs, not merely how `g` is read: `"h"`, `"e"`,
/// and `"p"` take the native Pieri / dual-Pieri / Murnaghan–Nakayama paths and
/// never touch Littlewood–Richardson, while `"s"`, `"m"`, and `"f"` go through
/// it. Passing the same function in a different basis gives the same answer by
/// a different algorithm, which is exactly what the oracle script checks.
#[pyfunction]
#[pyo3(signature = (f, g, basis = "s"))]
fn skew_by(f: Terms, g: Terms, basis: &str) -> PyResult<Terms> {
    fn run<C: Boundary>(f: &Terms, g: &Terms, basis: &str) -> Option<PyResult<Schur<C>>> {
        let sf: Schur<C> = build(f)?;
        Some(Ok(match basis {
            "s" => SkewBy::skew_by(&sf, &build::<C, Schur<C>>(g)?),
            "h" => SkewBy::skew_by(&sf, &build::<C, Homogeneous<C>>(g)?),
            "e" => SkewBy::skew_by(&sf, &build::<C, Elementary<C>>(g)?),
            "p" => SkewBy::skew_by(&sf, &build::<C, PowerSum<C>>(g)?),
            "m" => SkewBy::skew_by(&sf, &build::<C, Monomial<C>>(g)?),
            "f" => SkewBy::skew_by(&sf, &build::<C, Forgotten<C>>(g)?),
            other => {
                return Some(Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "unknown basis {other:?}; expected one of s, h, e, p, m, f"
                ))))
            }
        }))
    }
    escalate(
        || {
            let r = guarded(|| run::<Guarded>(&f, &g, basis))??;
            Some(r.map(|s| dump(&s)))
        },
        || run::<BigInt>(&f, &g, basis).unwrap().map(|s| dump(&s)),
    )
}

/// Evaluate a Schur-basis element at the alphabet `xs`.
///
/// Integer alphabet only: this is the bridge to concrete values, and a float
/// one would silently make an exact answer approximate. Rational alphabets are
/// the natural extension if a caller needs them.
#[pyfunction]
fn evaluate_schur(a: Terms, xs: Vec<Coeff>) -> Coeff {
    escalate(
        || {
            let s: Schur<Guarded> = build(&a)?;
            let alphabet: Option<Vec<Guarded>> =
                xs.iter().map(<Guarded as Boundary>::from_coeff).collect();
            let alphabet = alphabet?;
            Some(guarded(|| s.eval(&alphabet))?.to_coeff())
        },
        || {
            let s: Schur<BigInt> = build(&a).unwrap();
            let alphabet: Vec<BigInt> = xs.iter().map(Coeff::to_big).collect();
            Coeff::Big(s.eval(&alphabet))
        },
    )
}

/// f^λ — the number of standard Young tableaux of shape λ, i.e. the dimension
/// of the irreducible S_{|λ|} representation. `None` past `u128`.
#[pyfunction]
fn dimension(lambda: Vec<u32>) -> Option<u128> {
    crate::eval::dimension(&part(&lambda))
}

/// s_λ(1^n), the dimension of the GL_n irreducible. `None` on overflow.
#[pyfunction]
fn principal_specialization(lambda: Vec<u32>, n: u32) -> Option<u128> {
    crate::eval::principal_specialization(&part(&lambda), n)
}

/// s_λ(1, q, …, q^{n−1}) as a coefficient list in q, lowest degree first.
#[pyfunction]
fn principal_specialization_q(lambda: Vec<u32>, n: u32) -> Vec<i128> {
    crate::eval::principal_specialization_q(&part(&lambda), n)
}

/// Kostka number K_{λμ}.
#[pyfunction]
fn kostka_number(lambda: Vec<u32>, mu: Vec<u32>) -> u128 {
    crate::kostka::kostka(&part(&lambda), &part(&mu))
}

/// Symmetric-group character χ^λ(μ).
///
/// Exact at every size: `try_character` reports overflow rather than wrapping,
/// and the recursion then re-runs in `BigInt`. |χ^λ(μ)| ≤ √(|λ|!), which passes
/// `i128` around |λ| = 58 — reachable, so this is not hypothetical.
#[pyfunction]
fn character_value(lambda: Vec<u32>, mu: Vec<u32>) -> Coeff {
    let (l, m) = (part(&lambda), part(&mu));
    match crate::character::try_character(&l, &m) {
        Some(v) => Coeff::Small(v),
        None => Coeff::Big(crate::character::character_in::<BigInt>(&l, &m)),
    }
}

/// The internal (Kronecker) product of two Schur-basis elements.
#[pyfunction]
fn internal_product(a: Terms, b: Terms) -> PyResult<Terms> {
    escalate(
        || {
            let (x, y): (Schur<GuardedRat>, Schur<GuardedRat>) = (build_rat(&a)?, build_rat(&b)?);
            let r = guarded(|| ops::internal(&x, &y))?;
            Some(dump_integral(&r, "Kronecker"))
        },
        || {
            let (x, y): (Schur<BigRational>, Schur<BigRational>) =
                (build_rat(&a).unwrap(), build_rat(&b).unwrap());
            dump_integral(&ops::internal(&x, &y), "Kronecker")
        },
    )
}

/// A conversion whose output partitions are returned as **indices** rather than
/// as lists: `[(degree, index, coefficient), ...]`, where `index` is into
/// [`partitions`] of that degree.
///
/// The caller almost always has to turn each output partition into an object of
/// its own — a Sage `Partition`, say — and doing that per term dominates. A
/// list of parts must be copied, hashed and looked up before it can be mapped
/// to a cached object; an index is a direct array access. Measured on Sage's
/// conversion shim, that lookup was ~40% of everything outside symfn itself.
///
/// The order is `partitions(degree)`, which is exposed for exactly this reason,
/// so a caller can build its own table once per degree and never build another
/// partition object.
#[pyfunction]
fn convert_indexed(a: Terms, src: &str, dst: &str) -> PyResult<Vec<(u32, usize, Coeff)>> {
    let terms = if src == "Schur" {
        a
    } else {
        match src {
            "monomial" => monomial_to_schur(a),
            "homogeneous" => homogeneous_to_schur(a),
            "elementary" => elementary_to_schur(a),
            "powersum" => power_to_schur(a),
            "forgotten" => forgotten_to_schur(a),
            other => return Err(bad_basis(other)),
        }
    };
    let out = match dst {
        "Schur" => terms,
        "monomial" => schur_to_monomial(terms),
        "homogeneous" => schur_to_homogeneous(terms),
        "elementary" => schur_to_elementary(terms),
        "forgotten" => schur_to_forgotten(terms),
        other => return Err(bad_basis(other)),
    };
    Ok(out
        .into_iter()
        .map(|(p, c)| {
            let n: u32 = p.iter().sum();
            (n, index_of(n, &p), c)
        })
        .collect())
}

fn bad_basis(other: &str) -> PyErr {
    pyo3::exceptions::PyValueError::new_err(format!(
        "unknown basis {other:?}; expected one of Schur, monomial, homogeneous, elementary, powersum, forgotten"
    ))
}

thread_local! {
    /// Partition -> position in `partitions(n)`, per degree. Built on first use
    /// of a degree and reused; the ordering is fixed, so it never invalidates.
    static INDEX: std::cell::RefCell<std::collections::HashMap<u32, std::collections::HashMap<Vec<u32>, usize>>> =
        std::cell::RefCell::new(std::collections::HashMap::new());
}

fn index_of(n: u32, parts: &[u32]) -> usize {
    INDEX.with(|cell| {
        let mut m = cell.borrow_mut();
        let table = m.entry(n).or_insert_with(|| {
            crate::memo::partitions_cached(n)
                .iter()
                .enumerate()
                .map(|(i, p)| (p.parts().to_vec(), i))
                .collect()
        });
        table[parts]
    })
}

/// The partitions of `n`, in the order the table functions below index by.
///
/// Exposed so a caller can interpret [`character_table`] and [`kostka_table`]
/// without having to guess or replicate this crate's ordering.
#[pyfunction]
fn partitions(n: u32) -> Vec<Vec<u32>> {
    crate::memo::partitions_cached(n)
        .iter()
        .map(|p| p.parts().to_vec())
        .collect()
}

/// The full character table of S_n: `table[i][j]` = χ^{λⁱ}(λʲ).
///
/// A batched entry point, not a convenience wrapper. Built one value at a time
/// across the FFI boundary, a p(n)×p(n) table costs p(n)² calls — 393,129 at
/// n = 20 — and what that measures is Python dispatch, not the character
/// recursion. Symmetrica has had `chartafel` for the same reason.
///
/// Entries are `i128`-backed, so this is exact only while every χ^λ(μ) of
/// degree `n` fits — up to |λ| ≈ 58. Past that use [`character_value`], which
/// escalates per entry. The table sweep is built on an `i128` accumulator
/// throughout and cannot be widened by changing this signature alone.
#[pyfunction]
fn character_table(n: u32) -> Vec<Vec<i128>> {
    crate::character::character_table(n)
}

/// The full Kostka table of degree `n`: `table[i][j]` = K_{λⁱ λʲ}.
///
/// `u128`-backed, with the same caveat as [`character_table`].
#[pyfunction]
fn kostka_table(n: u32) -> Vec<Vec<u128>> {
    crate::kostka::kostka_table(n)
}

// --- operations -------------------------------------------------------------

/// The ω involution on a Schur-basis element.
///
/// Conjugates indices and copies coefficients, so it cannot overflow; it runs
/// once, over `BigInt`.
#[pyfunction]
fn omega(a: Terms) -> Terms {
    let s: Schur<BigInt> = build(&a).unwrap();
    dump(&s.omega())
}

/// The Hall inner product of two Schur-basis elements.
#[pyfunction]
fn hall_inner_product(a: Terms, b: Terms) -> Coeff {
    escalate(
        || {
            let (x, y): (Schur<Guarded>, Schur<Guarded>) = (build(&a)?, build(&b)?);
            Some(guarded(|| ops::hall::<Guarded, _, _>(&x, &y))?.to_coeff())
        },
        || {
            let (x, y): (Schur<BigInt>, Schur<BigInt>) = (build(&a).unwrap(), build(&b).unwrap());
            Coeff::Big(ops::hall::<BigInt, _, _>(&x, &y))
        },
    )
}

// --- Hopf structure ---------------------------------------------------------

/// The skew Schur function s_{λ/μ}.
#[pyfunction]
fn skew_schur(lambda: Vec<u32>, mu: Vec<u32>) -> Terms {
    let s: Schur<BigInt> = hopf::skew_schur(&part(&lambda), &part(&mu));
    dump(&s)
}

/// The coproduct Δ, as `[((mu, nu), coefficient), ...]`.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn coproduct(a: Terms) -> Vec<((Vec<u32>, Vec<u32>), Coeff)> {
    fn split<C: Boundary>(x: &Schur<C>) -> Vec<((Vec<u32>, Vec<u32>), Coeff)> {
        hopf::coproduct(x)
            .terms()
            .iter()
            .map(|((m, n), c)| ((m.parts().to_vec(), n.parts().to_vec()), c.to_coeff()))
            .collect()
    }
    escalate(
        || {
            let s: Schur<Guarded> = build(&a)?;
            guarded(|| split(&s))
        },
        || {
            let s: Schur<BigInt> = build(&a).unwrap();
            split(&s)
        },
    )
}

/// The antipode S.
///
/// Conjugates and negates, so like [`omega`] it cannot overflow.
#[pyfunction]
fn antipode(a: Terms) -> Terms {
    let s: Schur<BigInt> = build(&a).unwrap();
    dump(&hopf::antipode(&s))
}

// --- Hall–Littlewood --------------------------------------------------------

/// `Q'_λ(x; t) = Σ_μ K_{μλ}(t) s_μ`, as `[(mu, [(t_exponent, coefficient), ...])]`.
///
/// The coefficients are polynomials, so this cannot reuse [`Terms`]. Sparse in
/// the exponent, which is how [`QtPoly`](crate::QtPoly) already holds them.
#[pyfunction]
fn hall_littlewood(lambda: Vec<u32>) -> Vec<(Vec<u32>, Vec<(u32, Coeff)>)> {
    hl_rows(&crate::hall_littlewood::<i128>(&part(&lambda)))
}

/// Every `Q'_λ` for `λ ⊢ n`, sharing the recursion's suffixes across the degree.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn hall_littlewood_table(n: u32) -> Vec<(Vec<u32>, Vec<(Vec<u32>, Vec<(u32, Coeff)>)>)> {
    crate::hall_littlewood_table::<i128>(n)
        .into_iter()
        .map(|(lambda, hl)| (lambda.parts().to_vec(), hl_rows(&hl)))
        .collect()
}

/// `K_{λμ}(t)` as `[(t_exponent, coefficient), ...]`.
#[pyfunction]
fn kostka_foulkes(lambda: Vec<u32>, mu: Vec<u32>) -> Vec<(u32, Coeff)> {
    t_poly(&crate::kostka_foulkes::<i128>(&part(&lambda), &part(&mu)))
}

/// Every `K_{λμ}(t)` for a fixed μ, as `[(lambda, [(t_exponent, coefficient)])]`.
///
/// One `Q'_μ` *is* the column, so this costs what a single value costs — see
/// [`crate::kf`].
#[pyfunction]
fn kostka_foulkes_column(mu: Vec<u32>) -> Vec<(Vec<u32>, Vec<(u32, Coeff)>)> {
    crate::kostka_foulkes_column::<i128>(&part(&mu))
        .into_iter()
        .map(|(lambda, k)| (lambda.parts().to_vec(), t_poly(&k)))
        .collect()
}

/// `P_λ(x; t)` in the Schur basis — the other Hall–Littlewood normalisation.
///
/// Costs the whole degree: the inversion needs every dominance-smaller `P`, so
/// use [`hall_littlewood_p_table`] when more than one shape is wanted.
#[pyfunction]
fn hall_littlewood_p(lambda: Vec<u32>) -> Vec<(Vec<u32>, Vec<(u32, Coeff)>)> {
    hl_rows(&crate::hall_littlewood_p::<i128>(&part(&lambda)))
}

/// Every `P_λ` for `λ ⊢ n`, from one inversion of the Kostka–Foulkes matrix.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn hall_littlewood_p_table(n: u32) -> Vec<(Vec<u32>, Vec<(Vec<u32>, Vec<(u32, Coeff)>)>)> {
    crate::hall_littlewood_p_table::<i128>(n)
        .into_iter()
        .map(|(lambda, hl)| (lambda.parts().to_vec(), hl_rows(&hl)))
        .collect()
}

/// The whole `K_{λμ}(t)` matrix for degree `n`, indexed as `partitions(n)` is.
///
/// Same orientation as [`kostka_table`], of which this is the t-analogue:
/// `table[i][j]` is `K_{λⁱλʲ}(t)`, and `t = 1` recovers that table entry for
/// entry. Asking for the p(n)² values one at a time would recompute each column
/// p(n) times.
#[pyfunction]
fn kostka_foulkes_table(n: u32) -> Vec<Vec<Vec<(u32, Coeff)>>> {
    crate::kostka_foulkes_table::<i128>(n)
        .into_iter()
        .map(|row| row.iter().map(t_poly).collect())
        .collect()
}

/// A `Schur<QtPoly>` as `[(mu, [(t_exponent, coefficient), ...])]`.
fn hl_rows(hl: &Schur<crate::QtPoly<i128>>) -> Vec<(Vec<u32>, Vec<(u32, Coeff)>)> {
    hl.terms()
        .iter()
        .map(|(mu, c)| (mu.parts().to_vec(), t_poly(c)))
        .collect()
}

/// A `QtPoly` known not to involve q, as `[(t_exponent, coefficient)]`.
fn t_poly(p: &crate::QtPoly<i128>) -> Vec<(u32, Coeff)> {
    p.terms()
        .map(|((a, b), v)| {
            debug_assert_eq!(*a, 0, "Hall-Littlewood must not involve q");
            (*b, Coeff::Small(*v))
        })
        .collect()
}

// --- Macdonald --------------------------------------------------------------

/// One Macdonald expansion: per μ, the numerator's `(q_exp, t_exp, coeff)` terms
/// and the denominator's `(q_exp, t_exp, multiplicity)` **factors**.
///
/// The denominator is handed over factored rather than expanded, which is both
/// cheaper and what a caller wants: `prod((1 - q^a*t^b)^m)` builds the element
/// directly in a fraction field, where expanding here and re-factoring there
/// would be work done twice. See [`Frac`](crate::Frac) for why the factored form
/// is the representation and not an optimisation.
type MacTerms = Vec<(Vec<u32>, Vec<(u32, u32, Coeff)>, Vec<(u32, u32, u32)>)>;

fn mac_terms(f: &Monomial<crate::Frac<i128>>) -> MacTerms {
    f.terms()
        .iter()
        .map(|(mu, c)| {
            let (num, den) = c.parts();
            (
                mu.parts().to_vec(),
                num.terms()
                    .map(|((a, b), v)| (*a, *b, Coeff::Small(*v)))
                    .collect(),
                den.map(|(&(a, b), &m)| (a, b, m)).collect(),
            )
        })
        .collect()
}

/// Macdonald `P_λ(x; q, t)` in the monomial basis.
///
/// Coefficients run over `i128`, which is not the ceiling here: the widest
/// numerator coefficient through degree 10 is 31594374 — 25 bits against 127,
/// growing about 3.5 bits per degree (`examples/mac_coeff_sizes.rs`, which
/// checks the narrow run against a `BigInt` one). The enumeration becomes
/// impractical long before the width does.
#[pyfunction]
fn macdonald_p(lambda: Vec<u32>) -> MacTerms {
    mac_terms(&crate::macdonald_p::<i128>(&part(&lambda)))
}

/// Macdonald `Q_λ = b_λ · P_λ`.
#[pyfunction]
fn macdonald_q(lambda: Vec<u32>) -> MacTerms {
    mac_terms(&crate::macdonald_q::<i128>(&part(&lambda)))
}

/// Macdonald `J_λ = c_λ · P_λ`, the integral form — every coefficient is a
/// polynomial, so the denominator list comes back empty.
#[pyfunction]
fn macdonald_j(lambda: Vec<u32>) -> MacTerms {
    mac_terms(&crate::macdonald_j::<i128>(&part(&lambda)))
}

// --- Jack --------------------------------------------------------------------

/// One coefficient of a Jack expansion:
/// `(numerator, denominator atoms, integer scalar)`, meaning
///
/// ```text
///   (Σ_k num[k]·α^k) / (scale · ∏ (u·α + v)^mult)
/// ```
///
/// The numerator is **dense** — index is the α-exponent — because that is what
/// the object is; the denominator is handed over **factored**, for the same
/// reason [`MacTerms`] hands over its binomials factored: a caller rebuilding
/// this in ℚ(α) wants `prod(u*a + v)`, and expanding here to re-factor there is
/// work done twice.
///
/// The atoms are *primitive* (`gcd(u, v) = 1`), so the factorization is
/// canonical — unlike the (q,t) family, where `1 − q²` is reducible. See
/// [`AFrac`](crate::afrac::AFrac).
type JackCell = (Vec<Coeff>, Vec<(u32, u32, u32)>, u128);

/// One Jack expansion: per basis index μ, a [`JackCell`].
type JackTerms = Vec<(Vec<u32>, Vec<Coeff>, Vec<(u32, u32, u32)>, u128)>;

fn jack_cell<C: Boundary>(c: &crate::AFrac<C>) -> JackCell {
    let (num, den, scale) = c.parts();
    (
        num.iter().map(Boundary::to_coeff).collect(),
        den.map(|(&(u, v), &m)| (u, v, m)).collect(),
        scale,
    )
}

fn jack_terms<C: Boundary>(f: &Monomial<crate::AFrac<C>>) -> JackTerms {
    f.terms()
        .iter()
        .map(|(mu, c)| {
            let (n, d, s) = jack_cell(c);
            (mu.parts().to_vec(), n, d, s)
        })
        .collect()
}

fn jack_terms_p<C: Boundary>(f: &PowerSum<crate::AFrac<C>>) -> JackTerms {
    f.terms()
        .iter()
        .map(|(mu, c)| {
            let (n, d, s) = jack_cell(c);
            (mu.parts().to_vec(), n, d, s)
        })
        .collect()
}

/// Run a Jack computation over guarded `i128`, re-running over `BigInt` if
/// anything overflowed.
///
/// `AFrac<Guarded>` works only because [`Guarded`] implements
/// [`Ring::div_exact`] — without it the atom cancellation silently stops
/// happening and the denominators grow instead of the coefficients.
fn jack_escalate(
    fast: impl FnOnce() -> crate::AFrac<Guarded>,
    slow: impl FnOnce() -> crate::AFrac<BigInt>,
) -> JackCell {
    escalate(|| guarded(|| jack_cell(&fast())), || jack_cell(&slow()))
}

fn jack_escalate_m(
    fast: impl FnOnce() -> Monomial<crate::AFrac<Guarded>>,
    slow: impl FnOnce() -> Monomial<crate::AFrac<BigInt>>,
) -> JackTerms {
    escalate(|| guarded(|| jack_terms(&fast())), || jack_terms(&slow()))
}

/// Jack `P_λ(x; α)` in the monomial basis: monic in `m_λ`, dominance-triangular.
///
/// Computed by the Laplace–Beltrami eigenoperator recursion, which enumerates
/// no tableaux at all. Sage has no whole-degree entry point and walls at
/// n = 12; see [`jack_table`].
#[pyfunction]
fn jack_p(lambda: Vec<u32>) -> JackTerms {
    jack_escalate_m(
        || crate::jack_p(&part(&lambda)),
        || crate::jack_p(&part(&lambda)),
    )
}

/// Jack `Q_λ = (H_λ/H'_λ)·P_λ`, the basis dual to `P` under `⟨·,·⟩_α`.
#[pyfunction]
fn jack_q(lambda: Vec<u32>) -> JackTerms {
    jack_escalate_m(
        || crate::jack_q(&part(&lambda)),
        || crate::jack_q(&part(&lambda)),
    )
}

/// Jack `J_λ = H_λ·P_λ`, the integral form.
///
/// Every coefficient is a polynomial in α with non-negative integer
/// coefficients, divisible by `u_μ = ∏ m_i(μ)!` ([KS] Thm 1.1) — so the
/// denominator list comes back empty and `scale` comes back 1. None of that is
/// arranged: the coefficients arrive through fraction arithmetic and cancel.
#[pyfunction]
fn jack_j(lambda: Vec<u32>) -> JackTerms {
    jack_escalate_m(
        || crate::jack_j(&part(&lambda)),
        || crate::jack_j(&part(&lambda)),
    )
}

/// Every `P_λ` of degree `n` — the unit of work Sage has no entry point for,
/// and the one `docs/spec-jack.md` §2 measures the walls in.
#[pyfunction]
fn jack_table(n: u32) -> Vec<(Vec<u32>, JackTerms)> {
    crate::partitions_of(n)
        .into_iter()
        .map(|l| {
            let rows = jack_escalate_m(|| crate::jack_p(&l), || crate::jack_p(&l));
            (l.parts().to_vec(), rows)
        })
        .collect()
}

/// `J_λ` in the **power-sum** basis — the Jack character table, and the unit
/// the Goulden–Jackson pipeline consumes.
#[pyfunction]
fn jack_j_powersum(lambda: Vec<u32>) -> JackTerms {
    escalate(
        || guarded(|| jack_terms_p(&crate::jack_j_powersum::<Guarded>(&part(&lambda)))),
        || jack_terms_p(&crate::jack_j_powersum::<BigInt>(&part(&lambda))),
    )
}

/// `⟨J_λ, J_λ⟩_α = H_λ·H'_λ`, returned **factored** as `[(u, v, mult)]`.
///
/// A product of `2|λ|` linear forms and no pairing at all. Sage prices the same
/// table like a full expansion: over 360 s at n = 12.
#[pyfunction]
fn jack_norm_j(lambda: Vec<u32>) -> Vec<(u32, u32, u32)> {
    crate::jack_norm_j(&part(&lambda))
        .into_iter()
        .map(|((u, v), m)| (u, v, m as u32))
        .collect()
}

/// `⟨J_λ J_μ, J_ν⟩_α` — **Stanley's object**, whose membership in ℕ[α] is his
/// 1989 conjecture and still open.
///
/// A negative coefficient is a result to report, not a bug: nothing here
/// asserts positivity. Sage cannot compute `J[3,2,1]²` at all inside 120 s.
#[pyfunction]
fn jack_structure_constant(la: Vec<u32>, mu: Vec<u32>, nu: Vec<u32>) -> JackCell {
    let (a, b, c) = (part(&la), part(&mu), part(&nu));
    jack_escalate(
        || crate::jack_structure_constant(&a, &b, &c),
        || crate::jack_structure_constant(&a, &b, &c),
    )
}

/// Stanley's **whole table**: every `⟨J_λ J_μ, J_ν⟩_α` with `|λ| = |μ| = k`,
/// as `(lambda, mu, nu, numerator, denominator atoms, scalar)`.
///
/// Zero entries are omitted. Prefer this over looping
/// [`jack_structure_constant`]: it computes each p-expansion once, and
/// sampling the degree-12 table put 94.6% of the single-shot loop inside the
/// conversion it repeats. Measured 46.6 s → 1.44 s at k = 6, and k = 8
/// (degree 16, 111804 triples) is 74 s — sizes Sage cannot reach for even one
/// entry.
///
/// Positivity is Stanley's 1989 conjecture and is **open**. This returns the
/// values and asserts nothing about them.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn stanley_table(
    k: u32,
) -> Vec<(
    Vec<u32>,
    Vec<u32>,
    Vec<u32>,
    Vec<Coeff>,
    Vec<(u32, u32, u32)>,
    u128,
)> {
    escalate(
        || {
            guarded(|| {
                crate::stanley_table::<Guarded>(k)
                    .iter()
                    .map(|(la, mu, nu, g)| {
                        let (n, d, s) = jack_cell(g);
                        (
                            la.parts().to_vec(),
                            mu.parts().to_vec(),
                            nu.parts().to_vec(),
                            n,
                            d,
                            s,
                        )
                    })
                    .collect()
            })
        },
        || {
            crate::stanley_table::<BigInt>(k)
                .iter()
                .map(|(la, mu, nu, g)| {
                    let (n, d, s) = jack_cell(g);
                    (
                        la.parts().to_vec(),
                        mu.parts().to_vec(),
                        nu.parts().to_vec(),
                        n,
                        d,
                        s,
                    )
                })
                .collect()
        },
    )
}

/// `⟨f, g⟩_α` for two monomial-basis elements whose coefficients are
/// **integer** polynomials in α, given densely: `[(partition, [c0, c1, …])]`.
///
/// The restriction to integral coefficients is the honest boundary: a general
/// `AFrac` input would need the atoms marshalled in too, and every element a
/// caller actually pairs — `J_λ`, integer combinations of them — is already of
/// this shape. Clear denominators on the Python side first if yours is not.
#[pyfunction]
fn jack_scalar(f: Vec<(Vec<u32>, Vec<i128>)>, g: Vec<(Vec<u32>, Vec<i128>)>) -> JackCell {
    fn build<C: Ring>(rows: &[(Vec<u32>, Vec<i128>)]) -> Monomial<crate::AFrac<C>> {
        let mut out = Monomial::zero();
        for (mu, coeffs) in rows {
            let num: Vec<C> = coeffs.iter().map(|&v| C::from_i128(v)).collect();
            out.add_term(part(mu), crate::AFrac::from_coeffs(num));
        }
        out
    }
    jack_escalate(
        || crate::jack_scalar(&build::<Guarded>(&f), &build::<Guarded>(&g)),
        || crate::jack_scalar(&build::<BigInt>(&f), &build::<BigInt>(&g)),
    )
}

/// The zonal polynomial, in **both** circulating normalizations, as exact
/// `(numerator, denominator)` pairs.
///
/// ⚠️ Sage's `zonal()` is `P^{(2)}` and [GJ]'s `Z_λ` is `J^{(2)}`; the two
/// differ by `H_λ(2)`. Measured, not assumed. Both are returned rather than one
/// under an ambiguous name, because a caller that picks the wrong one still
/// gets plausible-looking output.
#[pyfunction]
fn zonal(lambda: Vec<u32>, integral_form: bool) -> Vec<(Vec<u32>, Coeff, Coeff)> {
    let l = part(&lambda);
    let f = if integral_form {
        crate::zonal_j(&l)
    } else {
        crate::zonal_p(&l)
    };
    f.terms()
        .iter()
        .map(|(mu, c)| {
            (
                mu.parts().to_vec(),
                Coeff::Small(c.numer()),
                Coeff::Small(c.denom()),
            )
        })
        .collect()
}

/// The Goulden–Jackson connection tables `c^λ_{μν}(b)` and `h^λ_{μν}(b)` at
/// degree `n`, as `(lambda, mu, nu, [b-coefficients], denominator)`.
///
/// Returns `(c, h)`. Two open conjectures live here — Matchings-Jack on `c`,
/// the b-conjecture on `h` — and no package computes either table.
/// ℚ[b]-polynomiality and `c`'s integrality are theorems and are enforced (a
/// failure raises); **positivity is the open question and is only observed**,
/// so a negative coefficient comes back as data rather than an exception.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn gj_connection_tables(
    n: u32,
) -> PyResult<(
    Vec<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<Coeff>, u128)>,
    Vec<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<Coeff>, u128)>,
)> {
    let t = crate::gj_connection_tables(n);
    if !t.laws_hold() {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "a PROVEN law failed at n = {n}: not polynomial at {:?}, c not integral at {:?}",
            t.not_polynomial, t.c_not_integral
        )));
    }
    let rows = |m: &std::collections::BTreeMap<crate::gj::Key, crate::BPoly>| {
        m.iter()
            .map(|((la, mu, nu), p)| {
                (
                    la.parts().to_vec(),
                    mu.parts().to_vec(),
                    nu.parts().to_vec(),
                    p.num.iter().map(|&v| Coeff::Small(v)).collect(),
                    p.den,
                )
            })
            .collect()
    };
    Ok((rows(&t.c), rows(&t.h)))
}

/// `a^λ_{μν}`, the class-algebra connection coefficient of `S_n`, from
/// characters alone.
///
/// The independent object the `b = 0` slice of [`gj_connection_tables`] is
/// pinned against — no Jack polynomial and no fraction field anywhere in it.
#[pyfunction]
fn class_algebra_coefficient(la: Vec<u32>, mu: Vec<u32>, nu: Vec<u32>) -> Coeff {
    Coeff::Small(crate::class_algebra_coefficient(
        &part(&la),
        &part(&mu),
        &part(&nu),
    ))
}

// --- (q,t)-Kostka -----------------------------------------------------------

/// A `QtPoly` over ℤ, as `[(q_exp, t_exp, coeff)]`.
///
/// No integrality assertion, and none is needed any more: these used to arrive
/// over ℚ from a route that divides by `z_ν`, where landing back in ℤ[q,t] was
/// Macdonald's theorem rather than anything the code arranged. The
/// Bergeron–Haiman recursion never divides by an integer, so this whole path now
/// runs over `i128` and a non-integral value is not representable rather than
/// merely unexpected. `Rat::into_poly` still refuses a surviving denominator,
/// and `divide_exact` still refuses an inexact division, which is where the
/// theorem is now enforced.
///
/// `i128` is not a ceiling: `K̃_{λμ}` has non-negative coefficients summing to
/// `f^λ`, and `Σ_λ (f^λ)² = n!`, so nothing here exceeds `√(n!)` — past `i128`
/// only around degree 57.
fn qt_poly(p: &crate::QtPoly<i128>) -> Vec<(u32, u32, Coeff)> {
    p.terms()
        .map(|(&(a, b), v)| (a, b, Coeff::Small(*v)))
        .collect()
}

/// The (q,t)-Kostka polynomial `K_{λμ}(q,t)`, from `J_μ = Σ_λ K_{λμ} S_λ(x;t)`.
///
/// Computes the whole of `J_μ`; use [`qt_kostka_column`] for more than one λ at
/// a fixed μ, and [`qt_kostka_table`] for a whole degree.
#[pyfunction]
fn qt_kostka(lambda: Vec<u32>, mu: Vec<u32>) -> Vec<(u32, u32, Coeff)> {
    qt_poly(&crate::qt_kostka::<i128>(&part(&lambda), &part(&mu)))
}

/// Every `K_{λμ}(q,t)` for a fixed μ — one `J_μ`, which is what a single
/// [`qt_kostka`] costs anyway.
///
/// Unlike a Kostka–Foulkes column this one is **dense**: every λ of the degree
/// appears, since `K_{λμ}` is generally nonzero without λ dominating μ.
#[pyfunction]
fn qt_kostka_column(mu: Vec<u32>) -> Vec<(Vec<u32>, Vec<(u32, u32, Coeff)>)> {
    crate::qt_kostka_column::<i128>(&part(&mu))
        .iter()
        .map(|(lambda, k)| (lambda.parts().to_vec(), qt_poly(k)))
        .collect()
}

/// The modified Macdonald polynomial `H̃_μ(x;q,t)` in the **Schur** basis, as
/// `[(lambda, [(q_exp, t_exp, coeff), ...])]`.
///
/// Its coefficients are the modified (q,t)-Kostka polynomials
/// `K̃_{λμ}(q,t) = t^{n(μ)} K_{λμ}(q, 1/t)` — the form the modern literature
/// uses, and where Haiman's positivity reads "non-negative integers" with no
/// normalising power in the way. `H̃_{(2)} = s_2 + q·s_{11}` and
/// `H̃_{(11)} = s_2 + t·s_{11}`.
#[pyfunction]
fn macdonald_ht(mu: Vec<u32>) -> Vec<(Vec<u32>, Vec<(u32, u32, Coeff)>)> {
    crate::macdonald_ht::<i128>(&part(&mu))
        .terms()
        .iter()
        .map(|(lambda, k)| (lambda.parts().to_vec(), qt_poly(k)))
        .collect()
}

/// The whole `K_{λμ}(q,t)` matrix for degree `n`, indexed as `partitions(n)` is
/// — the same orientation as [`kostka_table`] and [`kostka_foulkes_table`], of
/// which this is the two-variable analogue. `q = 0` recovers the latter.
#[pyfunction]
fn qt_kostka_table(n: u32) -> Vec<Vec<Vec<(u32, u32, Coeff)>>> {
    crate::qt_kostka_table::<i128>(n)
        .into_iter()
        .map(|row| row.iter().map(qt_poly).collect())
        .collect()
}

// --- the Macdonald operator algebra -----------------------------------------

/// A Schur element with `(q,t)`-polynomial coefficients, as
/// `[(lambda, [(q_exp, t_exp, coeff), ...]), ...]` — the same shape
/// [`macdonald_ht`] already returns, so an `H̃` row can be fed straight back in.
type QtSchur = Vec<(Vec<u32>, Vec<(u32, u32, Coeff)>)>;

fn qt_schur_out(f: &Schur<crate::QtPoly<i128>>) -> QtSchur {
    f.terms()
        .iter()
        .map(|(lambda, c)| (lambda.parts().to_vec(), qt_poly(c)))
        .collect()
}

/// The same, from the ℚ-bounded operators.
///
/// Every operator here maps `ℤ[q,t]`-Schur combinations to `ℤ[q,t]` ones, so a
/// surviving denominator is a bug and is raised rather than rounded. The general
/// path runs over `Rational` only because `s → p` divides by `z_ρ`; the answer
/// is integral by the time it reaches this boundary.
fn qt_schur_out_rat(f: &Schur<crate::QtPoly<crate::Rational>>, what: &str) -> PyResult<QtSchur> {
    let mut out = QtSchur::new();
    for (lambda, c) in f.terms() {
        let mut row = Vec::with_capacity(c.len());
        for (&(a, b), v) in c.terms() {
            if v.denom() != 1 {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "{what} is not integral at {lambda}: coefficient of q^{a}t^{b} is {v:?}"
                )));
            }
            row.push((a, b, Coeff::Small(v.numer())));
        }
        out.push((lambda.parts().to_vec(), row));
    }
    Ok(out)
}

fn qt_schur_in(rows: &QtSchur) -> PyResult<Schur<crate::QtPoly<crate::Rational>>> {
    let mut out = Schur::zero();
    for (lambda, terms) in rows {
        let mut c = crate::QtPoly::zero();
        for (a, b, v) in terms {
            let v = v.as_i128().ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err("coefficient does not fit in i128")
            })?;
            c.add_term(*a, *b, crate::Rational::from_int(v));
        }
        out.add_term(part(lambda), c);
    }
    Ok(out)
}

/// `∇e_n` in the Schur basis — the shuffle theorem's object.
///
/// The closed form, which never divides by an integer and so runs over ℤ.
/// Prefer this to `nabla(elementary)`: it skips the change of basis entirely.
#[pyfunction]
fn nabla_e(n: u32) -> QtSchur {
    qt_schur_out(&crate::nabla_e::<i128>(n))
}

/// `Δ'_{e_k} e_n` in the Schur basis — the Delta conjecture's object.
#[pyfunction]
fn delta_prime_e(k: u32, n: u32) -> QtSchur {
    qt_schur_out(&crate::delta_prime_e::<i128>(k, n))
}

/// `∇F` for an arbitrary homogeneous `F`, given in the Schur basis.
#[pyfunction]
fn nabla(f: QtSchur) -> PyResult<QtSchur> {
    qt_schur_out_rat(&crate::nabla(&qt_schur_in(&f)?), "nabla")
}

/// `∇^r F`, sharing one change of basis across the powers — the object
/// Qiu–Zhang's 2026 theorem is about.
#[pyfunction]
fn nabla_power(f: QtSchur, r: u32) -> PyResult<QtSchur> {
    qt_schur_out_rat(&crate::nabla_power(&qt_schur_in(&f)?, r), "nabla_power")
}

/// `Δ_{e_k} F`, with eigenvalue `e_k[B_μ]`.
#[pyfunction]
fn delta_ek(k: u32, f: QtSchur) -> PyResult<QtSchur> {
    let ek = crate::deltaop::elementary(k);
    qt_schur_out_rat(&crate::delta(&ek, &qt_schur_in(&f)?), "delta")
}

/// `Δ'_{e_k} F`, with eigenvalue `e_k[B_μ − 1]`.
#[pyfunction]
fn delta_prime_ek(k: u32, f: QtSchur) -> PyResult<QtSchur> {
    let ek = crate::deltaop::elementary(k);
    qt_schur_out_rat(&crate::delta_prime(&ek, &qt_schur_in(&f)?), "delta_prime")
}

/// `Θ_{e_k} F`, which raises the degree by `k`.
///
/// Note the cost: Θ expands at degree `n + k`, so it pays for the larger degree
/// and not the input's.
#[pyfunction]
fn theta_ek(k: u32, f: QtSchur) -> PyResult<QtSchur> {
    let ek = crate::deltaop::elementary(k);
    qt_schur_out_rat(&crate::theta(&ek, &qt_schur_in(&f)?), "theta")
}

/// `ΠF`, with eigenvalue `Π_μ`.
///
/// `Π⁻¹` is deliberately absent: it is genuinely not a polynomial, so it cannot
/// cross this boundary. Only the composite `Θ` can.
#[pyfunction]
fn big_pi(f: QtSchur) -> PyResult<QtSchur> {
    qt_schur_out_rat(&crate::big_pi(&qt_schur_in(&f)?), "big_pi")
}

/// The combinatorial side of the Delta conjecture, in the **monomial** basis,
/// for every `k` at once — entry `k` of the returned list.
///
/// `side` is `"rise"` (a theorem) or `"valley"` (open). One enumeration serves
/// the whole ladder, so asking for one `k` would cost the same.
///
/// The two sides no longer cost the same. `"rise"` factors through the per-path
/// LLT polynomials ([`crate::llt`], and `dyck.rs`'s module docs for why), which
/// is 29× faster at n = 8 and 56× at n = 9 — n = 9 costs about 1.3s. `"valley"`
/// keeps the `(n+1)^{n−1}`-ish labelled enumeration, because `Val` reads the
/// labels: ⚠️ about 75s at n = 9 and an order of magnitude more at n = 10.
#[pyfunction]
fn delta_conjecture_side(n: u32, side: &str) -> PyResult<Vec<QtSchur>> {
    let which = match side {
        "rise" => crate::Side::Rise,
        "valley" => crate::Side::Valley,
        other => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "side must be \"rise\" or \"valley\", got {other:?}"
            )))
        }
    };
    Ok(crate::ladder::<i128>(n, which)
        .iter()
        .map(|f| {
            f.terms()
                .iter()
                .map(|(mu, c)| (mu.parts().to_vec(), qt_poly(c)))
                .collect()
        })
        .collect())
}

// --- LLT polynomials --------------------------------------------------------

/// A **monomial**-basis element with `(q,t)`-polynomial coefficients, as
/// `[(mu, [(q_exp, t_exp, coeff), ...]), ...]`.
///
/// Structurally the same as [`QtSchur`] and deliberately a distinct alias: the
/// LLT families are naturally monomial-basis objects and the partitions index
/// *weights*, not Schur shapes. Feeding one to an operator expecting `QtSchur`
/// would type-check in Python and be wrong, so the names carry the warning that
/// the types cannot.
///
/// The `t` slot is zero for every LLT family proper — they live in `q` alone —
/// and carries the `area`/`maj` grading only in the assembly functions
/// ([`nabla_e_by_path`], [`htilde_by_llt`]).
type QtMon = Vec<(Vec<u32>, Vec<(u32, u32, Coeff)>)>;

/// `i128` for the same reason the (q,t)-Kostka family above uses it: every
/// coefficient here counts tableaux, so it is a non-negative integer bounded by
/// `n!` — 8.7e10 at n = 14, where `i128` holds 1.7e38. `bench_llt` runs the
/// whole ladder at `i64` *and* `i128` and asserts they agree term for term, so
/// the narrower width is checked rather than assumed.
fn qt_mon_out(f: &Monomial<crate::QtPoly<i128>>) -> QtMon {
    f.terms()
        .iter()
        .map(|(mu, c)| (mu.parts().to_vec(), qt_poly(c)))
        .collect()
}

/// A tuple of straight shapes with content offsets, as Python passes it.
fn skew_tuple(shapes: &[Vec<u32>], offsets: Option<Vec<i32>>) -> PyResult<crate::llt::SkewTuple> {
    let ps: Vec<Partition> = shapes.iter().map(|s| part(s)).collect();
    let offs = match offsets {
        None => vec![0i32; ps.len()],
        Some(o) if o.len() == ps.len() => o,
        Some(o) => {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "got {} offsets for {} components",
                o.len(),
                ps.len()
            )))
        }
    };
    Ok(crate::llt::SkewTuple::from_partitions(&ps, &offs))
}

fn decorated_graph(
    n: u32,
    weak: Vec<(u32, u32)>,
    strict: Vec<(u32, u32)>,
) -> PyResult<crate::llt::DecoratedGraph> {
    if let Some(&(u, v)) = strict.iter().find(|&&(u, v)| u >= v) {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "a strict edge is oriented u < v (the constraint is κ(u) < κ(v)); got ({u}, {v})"
        )));
    }
    if let Some(&(u, v)) = weak.iter().chain(&strict).find(|&&(u, v)| u >= n || v >= n) {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "edge ({u}, {v}) lands outside 0..{n}"
        )));
    }
    let un = |&(u, v): &(u32, u32)| (u.min(v), u.max(v));
    if let Some(e) = weak
        .iter()
        .map(un)
        .find(|e| strict.iter().map(un).any(|s| s == *e))
    {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "edge {e:?} is both weak (counted) and strict (constraining); it must be one or the other"
        )));
    }
    Ok(crate::llt::DecoratedGraph::new(n, &weak, &strict))
}

/// `G̃^(k)_λ(x;q)`, the **cospin** ribbon generating function of [LLT] (26), in
/// the monomial basis.
///
/// Empty when λ has no k-ribbon tableaux (nonempty k-core). Sage's
/// `llt(k).cospin(Partition(λ))` is the same object; `docs/spec-llt.md` §2.3 has
/// the mains-to-mains comparison.
#[pyfunction]
#[pyo3(signature = (lambda, k))]
fn llt_gtilde(lambda: Vec<u32>, k: u32) -> QtMon {
    qt_mon_out(&crate::llt::llt_gtilde::<i128>(&part(&lambda), k))
}

/// `H^(k)_μ(x;q) = Σ_R q^{s(R)} x^{w(R)}`, the **spin** family of [LLT] (28).
///
/// Takes a partition and a level, never a tuple, and that is a mathematical
/// constraint rather than an API choice: the k-quotient of a shape does not
/// determine `s*`, so there is no honest `H` of a bare tuple. Sage's
/// `llt(k).hspin()[μ]`.
#[pyfunction]
#[pyo3(signature = (mu, k))]
fn llt_h(mu: Vec<u32>, k: u32) -> QtMon {
    qt_mon_out(&crate::llt::llt_h::<i128>(&part(&mu), k))
}

/// `H̃^(k)_μ = G̃^(k)_{kμ}` ([LLT] (27)) — Sage's `llt(k).hcospin()[μ]`.
#[pyfunction]
#[pyo3(signature = (mu, k))]
fn llt_h_tilde(mu: Vec<u32>, k: u32) -> QtMon {
    qt_mon_out(&crate::llt::llt_h_tilde::<i128>(&part(&mu), k))
}

/// `Σ_R q^{2s(R)} x^{w(R)}`, the spin-generating grading of [LT] (43).
///
/// The rawest of the four normalizations, and the one [`llt_kl_column`] is
/// pinned against. Sage has no entry point for this grading.
#[pyfunction]
#[pyo3(signature = (lambda, k))]
fn llt_g_lt(lambda: Vec<u32>, k: u32) -> QtMon {
    qt_mon_out(&crate::llt::llt_g_lt::<i128>(&part(&lambda), k))
}

/// `H^(k)_μ` for **every** μ ⊢ n — the whole degree, which is the unit
/// `docs/spec-llt.md` §2 measures the walls in.
///
/// This is the entry point Sage lacks: there it is `p(n)` separate per-element
/// conversions, and the one-row shape alone is 94–100% of the cost.
#[pyfunction]
#[pyo3(signature = (n, k))]
fn llt_h_table(n: u32, k: u32) -> Vec<(Vec<u32>, QtMon)> {
    crate::llt::llt_h_table::<i128>(n, k)
        .iter()
        .map(|(mu, f)| (mu.parts().to_vec(), qt_mon_out(f)))
        .collect()
}

/// `G̃^(k)_λ` for **every** λ ⊢ k·n with empty k-core, from a single walk.
#[pyfunction]
#[pyo3(signature = (n, k))]
fn llt_gtilde_table(n: u32, k: u32) -> Vec<(Vec<u32>, QtMon)> {
    crate::llt::llt_gtilde_table::<i128>(n, k)
        .iter()
        .map(|(lambda, f)| (lambda.parts().to_vec(), qt_mon_out(f)))
        .collect()
}

/// `G̃^(k)_λ` in the **Schur** basis.
#[pyfunction]
#[pyo3(signature = (lambda, k))]
fn llt_schur(lambda: Vec<u32>, k: u32) -> QtSchur {
    qt_schur_out(&crate::llt::llt_schur::<i128>(&part(&lambda), k))
}

/// `G_ν(x;q)` for a tuple of shapes, in the monomial basis and the **raw** inv
/// grading.
///
/// `offsets` defaults to all zero. ⚠️ The floor is **not** divided out:
/// `min_T inv(T)` can be positive, and Sage's `llt(k).cospin(tuple)` returns
/// `q^{−min inv} G_ν` instead. Divide by `q^{llt_min_inv(...)}` to compare —
/// exposing the floor is deliberate, since it is real data about ν and hiding it
/// is how the quotient dictionary gets misread (`docs/spec-llt.md` §1.3(a)).
#[pyfunction]
#[pyo3(signature = (shapes, offsets=None))]
fn llt_g(shapes: Vec<Vec<u32>>, offsets: Option<Vec<i32>>) -> PyResult<QtMon> {
    Ok(qt_mon_out(&crate::llt::llt_g::<i128>(&skew_tuple(
        &shapes, offsets,
    )?)))
}

/// `min_T inv(T)` over the semistandard fillings of a tuple — the forced
/// `q`-floor that [`llt_g`] does not divide out.
#[pyfunction]
#[pyo3(signature = (shapes, offsets=None))]
fn llt_min_inv(shapes: Vec<Vec<u32>>, offsets: Option<Vec<i32>>) -> PyResult<u32> {
    Ok(crate::llt::llt_min_inv(&skew_tuple(&shapes, offsets)?))
}

/// The **fundamental quasisymmetric** expansion of `G_ν`, as
/// `[(composition, [(q_exp, t_exp, coeff), ...]), ...]`.
///
/// [HHL] (82)'s descent buckets read directly. No package ships this expansion,
/// and the crate has no QSym type — the compositions carry their own meaning and
/// nothing here multiplies them.
#[pyfunction]
#[pyo3(signature = (shapes, offsets=None))]
fn llt_fundamental(
    shapes: Vec<Vec<u32>>,
    offsets: Option<Vec<i32>>,
) -> PyResult<Vec<(Vec<u32>, Vec<(u32, u32, Coeff)>)>> {
    Ok(
        crate::llt::llt_fundamental::<i128>(&skew_tuple(&shapes, offsets)?)
            .iter()
            .map(|(comp, c)| (comp.clone(), qt_poly(c)))
            .collect(),
    )
}

/// The k-core and k-quotient of λ, as `(core, [component, ...])`.
///
/// The abacus primitives the ribbon model rests on. Component **order** (runner
/// 0 first) is load-bearing — `G_ν` is not symmetric in its components — and
/// agrees with Sage's `Partition(λ).quotient(k)`, which `scripts/check_llt.py`
/// checks.
#[pyfunction]
#[pyo3(signature = (lambda, k))]
fn k_core_quotient(lambda: Vec<u32>, k: u32) -> (Vec<u32>, Vec<Vec<u32>>) {
    let l = part(&lambda);
    (
        l.k_core(k).parts().to_vec(),
        l.k_quotient(k).iter().map(|p| p.parts().to_vec()).collect(),
    )
}

/// `∇e_n = Σ_D t^{area(D)} G_D(x;q)`, as `[(area_sequence, G_D), ...]`.
///
/// The by-path Schur-positive refinement of the shuffle theorem — `∇e_n` written
/// as a positive sum of positive pieces. No package emits this decomposition,
/// and it is what makes the rise side of the Delta conjecture cheap (see
/// [`delta_conjecture_side`]).
///
/// ⚠️ `C_n` pieces and `#SYT` work each: n = 10 is 16 796 pieces in about 19s.
/// Use [`nabla_e`] for the total, which is far cheaper.
#[pyfunction]
fn nabla_e_by_path(n: u32) -> Vec<(Vec<u32>, QtMon)> {
    crate::llt::nabla_e_by_path::<i128>(n)
        .iter()
        .map(|(area, g)| (area.clone(), qt_mon_out(g)))
        .collect()
}

/// One **column** of the Schur-expansion table: `c^λ_μ` for every shape μ ⊢ k|λ|,
/// in the [KMS] variable `v`, as `[(mu, [(v_exp, 0, coeff), ...]), ...]`.
///
/// These are parabolic affine Kazhdan–Lusztig polynomials ([LT] Thm 4.2),
/// computed by exact Fock-space straightening with no Hecke algebra in sight.
/// Sage has no entry point of this shape.
///
/// ⚠️ The variable is `v`, and the ribbon side's grading is recovered at
/// **`q = −v`**. Coefficients are signed for that reason.
#[pyfunction]
#[pyo3(signature = (lambda, k))]
fn llt_kl_column(lambda: Vec<u32>, k: u32) -> Vec<(Vec<u32>, Vec<(u32, u32, Coeff)>)> {
    crate::llt::llt_kl_column::<i128>(&part(&lambda), k)
        .iter()
        .map(|(mu, c)| (mu.parts().to_vec(), qt_poly(c)))
        .collect()
}

/// `G_Γ(x;q) = Σ_κ q^{asc(κ)} x^κ` over the colorings of a decorated graph.
///
/// Vertices are `0 … n−1`. `weak` edges are ordered pairs whose *ascents* are
/// counted (`κ(u) < κ(v)` scores a `q`); `strict` edges must run `u < v` and
/// *constrain* (`κ(u) < κ(v)` or the coloring does not count), never scoring.
/// The two sets must be disjoint as unordered pairs, or the statistic silently
/// gains a `q` per strict edge — which is why that is raised rather than
/// tolerated.
#[pyfunction]
#[pyo3(signature = (n, weak, strict))]
fn llt_graph(n: u32, weak: Vec<(u32, u32)>, strict: Vec<(u32, u32)>) -> PyResult<QtMon> {
    Ok(qt_mon_out(&crate::llt::llt_graph::<i128>(
        &decorated_graph(n, weak, strict)?,
    )))
}

/// The Shareshian–Wachs chromatic quasisymmetric function `X_Γ(x;q)` of Γ, from
/// its LLT polynomial by the `(q−1)`-plethysm of [CM] Prop 3.5.
///
/// Γ should carry the [CM] presentation — natural orientation, no strict edges —
/// for the answer to be the chromatic function of the graph rather than of a
/// decorated relative of it. **Isolated vertices are part of Γ and must be
/// counted in `n`**; dropping them is a well-trodden route to a plausible wrong
/// answer.
///
/// Runs over ℚ because the plethysm goes through the power-sum basis, which
/// divides by `z_ρ`; the answer is integral by the time it crosses, and a
/// surviving denominator is raised rather than rounded.
#[pyfunction]
#[pyo3(signature = (n, weak, strict))]
fn chromatic_from_llt(n: u32, weak: Vec<(u32, u32)>, strict: Vec<(u32, u32)>) -> PyResult<QtMon> {
    let g = decorated_graph(n, weak, strict)?;
    let f = crate::llt::chromatic_from_llt::<crate::Rational>(&g);
    let mut out = QtMon::new();
    for (mu, c) in f.terms() {
        let mut row = Vec::with_capacity(c.len());
        for (&(a, b), v) in c.terms() {
            if v.denom() != 1 {
                return Err(pyo3::exceptions::PyValueError::new_err(format!(
                    "the chromatic function is not integral at {mu}: \
                     coefficient of q^{a}t^{b} is {v:?}"
                )));
            }
            row.push((a, b, Coeff::Small(v.numer())));
        }
        out.push((mu.parts().to_vec(), row));
    }
    Ok(out)
}

/// The [AS] **e-expansion** of `Ĝ_Γ(x; q+1)`: `Σ_θ q^{asc(θ)} e_{λ(θ)}` over
/// orientations of the free edges, as `[(partition, poly), ...]`.
///
/// By [DA]'s theorem the coefficients are non-negative, so this is **certified
/// positive output** rather than a conjecture to check. Weak edges are read as
/// unordered pairs: [AS]'s formula orients them itself.
///
/// ⚠️ `2^{#free edges}` terms.
#[pyfunction]
#[pyo3(signature = (n, weak, strict))]
fn llt_e_expansion(
    n: u32,
    weak: Vec<(u32, u32)>,
    strict: Vec<(u32, u32)>,
) -> PyResult<Vec<(Vec<u32>, Vec<(u32, u32, Coeff)>)>> {
    let g = decorated_graph(n, weak, strict)?;
    Ok(crate::llt::llt_e_expansion::<i128>(&g)
        .iter()
        .map(|(lambda, c)| (lambda.parts().to_vec(), qt_poly(c)))
        .collect())
}

/// `H̃_μ(x;q,t) = Σ_D q^{−a(D)} t^{maj(D)} G_{ν(μ,D)}(x;q)` — the [HHL]
/// decomposition, in the monomial basis.
///
/// The fourth route to `H̃` in this crate and the only one positively graded at
/// every intermediate step: Macdonald positivity *is* LLT positivity, and this
/// is where that becomes a computation.
///
/// ⚠️ A reference route, not a fast one — `2^{|μ|−μ₁}` LLT evaluations. Measured,
/// it ties [`macdonald_ht`] at n = 8 and loses 3× at n = 9, growing. Use
/// [`macdonald_ht`] to *compute* `H̃`; use this to check it independently.
#[pyfunction]
fn htilde_by_llt(mu: Vec<u32>) -> QtMon {
    qt_mon_out(&crate::llt::htilde_by_llt::<i128>(&part(&mu)))
}

#[pymodule]
fn symfn(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hall_littlewood, m)?)?;
    m.add_function(wrap_pyfunction!(hall_littlewood_table, m)?)?;
    m.add_function(wrap_pyfunction!(kostka_foulkes, m)?)?;
    m.add_function(wrap_pyfunction!(kostka_foulkes_column, m)?)?;
    m.add_function(wrap_pyfunction!(kostka_foulkes_table, m)?)?;
    m.add_function(wrap_pyfunction!(hall_littlewood_p, m)?)?;
    m.add_function(wrap_pyfunction!(hall_littlewood_p_table, m)?)?;
    m.add_function(wrap_pyfunction!(macdonald_p, m)?)?;
    m.add_function(wrap_pyfunction!(macdonald_q, m)?)?;
    m.add_function(wrap_pyfunction!(macdonald_j, m)?)?;
    m.add_function(wrap_pyfunction!(jack_p, m)?)?;
    m.add_function(wrap_pyfunction!(jack_q, m)?)?;
    m.add_function(wrap_pyfunction!(jack_j, m)?)?;
    m.add_function(wrap_pyfunction!(jack_table, m)?)?;
    m.add_function(wrap_pyfunction!(jack_j_powersum, m)?)?;
    m.add_function(wrap_pyfunction!(jack_norm_j, m)?)?;
    m.add_function(wrap_pyfunction!(jack_scalar, m)?)?;
    m.add_function(wrap_pyfunction!(jack_structure_constant, m)?)?;
    m.add_function(wrap_pyfunction!(stanley_table, m)?)?;
    m.add_function(wrap_pyfunction!(zonal, m)?)?;
    m.add_function(wrap_pyfunction!(gj_connection_tables, m)?)?;
    m.add_function(wrap_pyfunction!(class_algebra_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(qt_kostka, m)?)?;
    m.add_function(wrap_pyfunction!(qt_kostka_column, m)?)?;
    m.add_function(wrap_pyfunction!(qt_kostka_table, m)?)?;
    m.add_function(wrap_pyfunction!(macdonald_ht, m)?)?;
    m.add_function(wrap_pyfunction!(nabla_e, m)?)?;
    m.add_function(wrap_pyfunction!(delta_prime_e, m)?)?;
    m.add_function(wrap_pyfunction!(nabla, m)?)?;
    m.add_function(wrap_pyfunction!(nabla_power, m)?)?;
    m.add_function(wrap_pyfunction!(delta_ek, m)?)?;
    m.add_function(wrap_pyfunction!(delta_prime_ek, m)?)?;
    m.add_function(wrap_pyfunction!(theta_ek, m)?)?;
    m.add_function(wrap_pyfunction!(big_pi, m)?)?;
    m.add_function(wrap_pyfunction!(delta_conjecture_side, m)?)?;
    m.add_function(wrap_pyfunction!(llt_gtilde, m)?)?;
    m.add_function(wrap_pyfunction!(llt_h, m)?)?;
    m.add_function(wrap_pyfunction!(llt_h_tilde, m)?)?;
    m.add_function(wrap_pyfunction!(llt_g_lt, m)?)?;
    m.add_function(wrap_pyfunction!(llt_h_table, m)?)?;
    m.add_function(wrap_pyfunction!(llt_gtilde_table, m)?)?;
    m.add_function(wrap_pyfunction!(llt_schur, m)?)?;
    m.add_function(wrap_pyfunction!(llt_g, m)?)?;
    m.add_function(wrap_pyfunction!(llt_min_inv, m)?)?;
    m.add_function(wrap_pyfunction!(llt_fundamental, m)?)?;
    m.add_function(wrap_pyfunction!(llt_kl_column, m)?)?;
    m.add_function(wrap_pyfunction!(llt_graph, m)?)?;
    m.add_function(wrap_pyfunction!(chromatic_from_llt, m)?)?;
    m.add_function(wrap_pyfunction!(llt_e_expansion, m)?)?;
    m.add_function(wrap_pyfunction!(htilde_by_llt, m)?)?;
    m.add_function(wrap_pyfunction!(nabla_e_by_path, m)?)?;
    m.add_function(wrap_pyfunction!(k_core_quotient, m)?)?;
    m.add_function(wrap_pyfunction!(clear_caches, m)?)?;
    m.add_function(wrap_pyfunction!(schur_multiply, m)?)?;
    m.add_function(wrap_pyfunction!(lr_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_homogeneous, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_elementary, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_monomial, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_forgotten, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_power, m)?)?;
    m.add_function(wrap_pyfunction!(homogeneous_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(elementary_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(monomial_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(forgotten_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(power_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(plethysm, m)?)?;
    m.add_function(wrap_pyfunction!(kostka_number, m)?)?;
    m.add_function(wrap_pyfunction!(evaluate_schur, m)?)?;
    m.add_function(wrap_pyfunction!(dimension, m)?)?;
    m.add_function(wrap_pyfunction!(principal_specialization, m)?)?;
    m.add_function(wrap_pyfunction!(principal_specialization_q, m)?)?;
    m.add_function(wrap_pyfunction!(character_value, m)?)?;
    m.add_function(wrap_pyfunction!(internal_product, m)?)?;
    m.add_function(wrap_pyfunction!(partitions, m)?)?;
    m.add_function(wrap_pyfunction!(convert_indexed, m)?)?;
    m.add_function(wrap_pyfunction!(character_table, m)?)?;
    m.add_function(wrap_pyfunction!(kostka_table, m)?)?;
    m.add_function(wrap_pyfunction!(omega, m)?)?;
    m.add_function(wrap_pyfunction!(hall_inner_product, m)?)?;
    m.add_function(wrap_pyfunction!(skew_schur, m)?)?;
    m.add_function(wrap_pyfunction!(skew_by, m)?)?;
    m.add_function(wrap_pyfunction!(coproduct, m)?)?;
    m.add_function(wrap_pyfunction!(antipode, m)?)?;
    Ok(())
}
