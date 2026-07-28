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
//! Coefficients are `i128` on this boundary (rationals cross as `(num, den)`).
//! `i64` was too narrow: it silently truncated plethysm numerators and any
//! structure constant past ~9.2e18. Python integers are arbitrary precision, so
//! the ceiling here is symfn's, not Python's — see the note in `ROADMAP.md`.
//! Arbitrary-precision passthrough — carrying `rug::Integer` across as decimal
//! strings — is the natural follow-up once a workload needs it.

use pyo3::prelude::*;

use crate::convert::{FromSchur, ToSchur};
use crate::hopf;
use crate::lr::{LrBackend, NaiveLr};
use crate::ops;
use crate::partition::Partition;
use crate::sym::{Elementary, Homogeneous, Monomial, PowerSum, Schur, SymFn};
use crate::Rational;

type Terms = Vec<(Vec<u32>, i128)>;
type RatTerms = Vec<(Vec<u32>, (i128, i128))>;

fn part(p: &[u32]) -> Partition {
    Partition::new(p.iter().copied())
}

fn build_schur(terms: &Terms) -> Schur<i128> {
    let mut s = Schur::zero();
    for (p, c) in terms {
        s.add_term(part(p), *c);
    }
    s
}

fn dump<S: SymFn<i128>>(x: &S) -> Terms {
    x.terms()
        .iter()
        .map(|(p, c)| (p.parts().to_vec(), *c))
        .collect()
}

// --- products ---------------------------------------------------------------

/// Multiply two Schur-basis elements (Littlewood–Richardson).
#[pyfunction]
fn schur_multiply(a: Terms, b: Terms) -> Terms {
    dump(&build_schur(&a).mul(&build_schur(&b)))
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

#[pyfunction]
fn schur_to_homogeneous(a: Terms) -> Terms {
    dump(&Homogeneous::from_schur(&build_schur(&a)))
}

#[pyfunction]
fn schur_to_elementary(a: Terms) -> Terms {
    dump(&Elementary::from_schur(&build_schur(&a)))
}

#[pyfunction]
fn schur_to_monomial(a: Terms) -> Terms {
    dump(&Monomial::from_schur(&build_schur(&a)))
}

/// s → p. Coefficients are rational, returned as `(numerator, denominator)`.
#[pyfunction]
fn schur_to_power(a: Terms) -> RatTerms {
    let mut s: Schur<Rational> = Schur::zero();
    for (p, c) in &a {
        s.add_term(part(p), Rational::from_int(*c));
    }
    let p: PowerSum<Rational> = PowerSum::from_schur(&s);
    p.terms()
        .iter()
        .map(|(part, c)| (part.parts().to_vec(), (c.numer(), c.denom())))
        .collect()
}

// --- conversions into Schur -------------------------------------------------

macro_rules! into_schur {
    ($name:ident, $basis:ident) => {
        #[pyfunction]
        fn $name(a: Terms) -> Terms {
            let mut x = $basis::<i128>::zero();
            for (p, c) in &a {
                x.add_term(part(p), *c);
            }
            dump(&x.to_schur())
        }
    };
}

into_schur!(homogeneous_to_schur, Homogeneous);
into_schur!(elementary_to_schur, Elementary);
into_schur!(monomial_to_schur, Monomial);
into_schur!(power_to_schur, PowerSum);

/// Plethysm f[g] of two Schur-basis elements.
///
/// Computed through the power-sum basis (see `crate::plethysm`). Schur inputs
/// give integral output; a non-integral coefficient would signal a bug and is
/// reported as a Python `ValueError` rather than silently truncated.
#[pyfunction]
fn plethysm(f: Terms, g: Terms) -> PyResult<Terms> {
    let mut sf: Schur<Rational> = Schur::zero();
    for (pp, c) in &f {
        sf.add_term(part(pp), Rational::from_int(*c));
    }
    let mut sg: Schur<Rational> = Schur::zero();
    for (pp, c) in &g {
        sg.add_term(part(pp), Rational::from_int(*c));
    }
    let r = crate::plethysm::plethysm(&sf, &sg);
    let mut out = Vec::new();
    for (pp, c) in r.terms() {
        if c.denom() != 1 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "non-integral plethysm coefficient {c:?}"
            )));
        }
        out.push((pp.parts().to_vec(), c.numer()));
    }
    Ok(out)
}

// --- classical quantities ---------------------------------------------------

/// Kostka number K_{λμ}.
#[pyfunction]
fn kostka_number(lambda: Vec<u32>, mu: Vec<u32>) -> u128 {
    crate::kostka::kostka(&part(&lambda), &part(&mu))
}

/// Symmetric-group character χ^λ(μ).
#[pyfunction]
fn character_value(lambda: Vec<u32>, mu: Vec<u32>) -> i128 {
    crate::character::character(&part(&lambda), &part(&mu))
}

/// The internal (Kronecker) product of two Schur-basis elements.
///
/// Schur inputs give integral output; a non-integral coefficient would mean a
/// bug in the power-sum route rather than a representable answer, so it is
/// reported rather than truncated — same contract as `plethysm`.
#[pyfunction]
fn internal_product(a: Terms, b: Terms) -> PyResult<Terms> {
    let mut sa: Schur<Rational> = Schur::zero();
    for (p, c) in &a {
        sa.add_term(part(p), Rational::from_int(*c));
    }
    let mut sb: Schur<Rational> = Schur::zero();
    for (p, c) in &b {
        sb.add_term(part(p), Rational::from_int(*c));
    }
    let r = crate::ops::internal(&sa, &sb);
    let mut out = Vec::new();
    for (p, c) in r.terms() {
        if c.denom() != 1 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "non-integral Kronecker coefficient {c:?}"
            )));
        }
        out.push((p.parts().to_vec(), c.numer()));
    }
    Ok(out)
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
#[pyfunction]
fn character_table(n: u32) -> Vec<Vec<i128>> {
    crate::character::character_table(n)
}

/// The full Kostka table of degree `n`: `table[i][j]` = K_{λⁱ λʲ}.
#[pyfunction]
fn kostka_table(n: u32) -> Vec<Vec<u128>> {
    crate::kostka::kostka_table(n)
}

// --- operations -------------------------------------------------------------

/// The ω involution on a Schur-basis element.
#[pyfunction]
fn omega(a: Terms) -> Terms {
    dump(&build_schur(&a).omega())
}

/// The Hall inner product of two Schur-basis elements.
#[pyfunction]
fn hall_inner_product(a: Terms, b: Terms) -> i128 {
    ops::hall(&build_schur(&a), &build_schur(&b))
}

// --- Hopf structure ---------------------------------------------------------

/// The skew Schur function s_{λ/μ}.
#[pyfunction]
fn skew_schur(lambda: Vec<u32>, mu: Vec<u32>) -> Terms {
    let s: Schur<i128> = hopf::skew_schur(&part(&lambda), &part(&mu));
    dump(&s)
}

/// The coproduct Δ, as `[((mu, nu), coefficient), ...]`.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn coproduct(a: Terms) -> Vec<((Vec<u32>, Vec<u32>), i128)> {
    hopf::coproduct(&build_schur(&a))
        .terms()
        .iter()
        .map(|((m, n), c)| ((m.parts().to_vec(), n.parts().to_vec()), *c))
        .collect()
}

/// The antipode S.
#[pyfunction]
fn antipode(a: Terms) -> Terms {
    dump(&hopf::antipode(&build_schur(&a)))
}

#[pymodule]
fn symfn(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(clear_caches, m)?)?;
    m.add_function(wrap_pyfunction!(schur_multiply, m)?)?;
    m.add_function(wrap_pyfunction!(lr_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_homogeneous, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_elementary, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_monomial, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_power, m)?)?;
    m.add_function(wrap_pyfunction!(homogeneous_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(elementary_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(monomial_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(power_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(plethysm, m)?)?;
    m.add_function(wrap_pyfunction!(kostka_number, m)?)?;
    m.add_function(wrap_pyfunction!(character_value, m)?)?;
    m.add_function(wrap_pyfunction!(internal_product, m)?)?;
    m.add_function(wrap_pyfunction!(partitions, m)?)?;
    m.add_function(wrap_pyfunction!(character_table, m)?)?;
    m.add_function(wrap_pyfunction!(kostka_table, m)?)?;
    m.add_function(wrap_pyfunction!(omega, m)?)?;
    m.add_function(wrap_pyfunction!(hall_inner_product, m)?)?;
    m.add_function(wrap_pyfunction!(skew_schur, m)?)?;
    m.add_function(wrap_pyfunction!(coproduct, m)?)?;
    m.add_function(wrap_pyfunction!(antipode, m)?)?;
    Ok(())
}
