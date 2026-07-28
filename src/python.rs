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
            let (x, y): (Schur<GuardedRat>, Schur<GuardedRat>) =
                (build_rat(&f)?, build_rat(&g)?);
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
            let (x, y): (Schur<GuardedRat>, Schur<GuardedRat>) =
                (build_rat(&a)?, build_rat(&b)?);
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

#[pymodule]
fn symfn(m: &Bound<'_, PyModule>) -> PyResult<()> {
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
