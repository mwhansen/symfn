//! PyO3 bridge (feature = "python"): a **coarse-grained** API importable into
//! Sage.
//!
//! The design rule established up front: cross-language calls have overhead, so
//! every function here does a *whole-object* operation — multiply two complete
//! symmetric functions, convert an entire element between bases — never a
//! per-monomial call. Marshalling happens once at the boundary; the inner loops
//! stay in Rust. Exposing fine-grained accessors would reintroduce exactly the
//! per-call tax this crate exists to escape.
//!
//! The standing policy for this surface — the three-layer boundary, where a
//! new entry point lands, and what each one owes its callers — is
//! `docs/policies/python.md`.
//!
//! Elements cross the boundary as lists of `(partition, coefficient)` pairs,
//! e.g. `[((2, 1), 3), ((3,), -1)]`, which maps directly onto Sage's
//! `.monomial_coefficients()` dicts.
//!
//! **A partition comes back as a tuple and goes in as any sequence.** Outbound
//! it is a tuple because every consumer uses it as a key, and a list would have
//! to be copied into one first, once per output term.
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
//! This is [`character_in`](crate::character::character_in)'s pattern applied
//! to every coefficient, and it exists because `impl Ring for i128` multiplies
//! with a plain `*`: without the escalation a structure constant past ~1.7e38
//! comes back **wrapped, with no signal**.
//!
//! Escalation is rare — measured coefficient widths in these workloads are 1–2
//! limbs (`examples/coeff_sizes.rs`) — and
//! `docs/record/python-and-sage-interop.md` carries the fast path's measured
//! cost against unchecked arithmetic (`examples/bench_guarded.rs`).
//!
//! ## Every entry point answers Ctrl-C
//!
//! **Any function here can raise `KeyboardInterrupt`**, and the individual
//! `# Raises` sections do not repeat it. Whether a call is interruptible is not
//! a property of its mathematics but of how long it runs, so the fact belongs
//! once, here (`docs/policies/python.md`, P8).
//!
//! Importing this module installs a checker into
//! [`interrupt`](crate::interrupt) that runs CPython's pending signal handlers,
//! and every body below is wrapped in `interruptible`. Without that, a long
//! call held the GIL for its whole duration and CPython could not run the
//! handler until it returned — a Ctrl-C during a multi-minute conversion did
//! nothing at all, which is how Sage users met it
//! (`docs/record/python-and-sage-interop.md`).
//!
//! The exception raised is the one the signal handler actually raised, not a
//! `KeyboardInterrupt` manufactured here: a handler is arbitrary Python and
//! need not raise that, or anything.

// A content offset, bounded by the tuple the caller passed.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::One;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyAny;

use crate::coeff::Ring;
use crate::convert::{convert, FromSchur, ToSchur};
use crate::guard::{guarded, Guarded, GuardedRat};
use crate::hopf::{self, SkewBy};
use crate::lr::{LrBackend, NaiveLr};
use crate::ops;
use crate::partition::Partition;
use crate::permutation::{Perm, MAX_SUPPORT};
use crate::schubert::Schubert;
use crate::sym::{Elementary, Forgotten, Homogeneous, Ht, Monomial, PowerSum, Schur, St, SymFn};

/// A coefficient crossing the boundary.
///
/// Python only ever sees an `int` — this distinction is invisible there. It
/// exists because `BigInt` is heap-allocated and almost every coefficient is
/// small: routing all of them through it puts an allocation on every term, and
/// a marshalling-dominated call such as `coproduct` is where that lands hardest
/// (`docs/record/python-and-sage-interop.md` has the measurement). PyO3
/// converts `i128` with no allocation, so the common case pays nothing and only
/// genuinely wide values allocate.
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

/// A combinatorial support crossing the boundary: a **tuple** on the Python
/// side. A partition for the symmetric-function families, a one-line
/// permutation word for the Schuberts.
///
/// Outbound it is a tuple and not a list because every consumer uses it as a
/// key. Sage's adapter builds `{Partition: coefficient}` dicts and interns on
/// the tuple of parts, and a bare CPython caller reaches for a `dict` or a
/// `set` just as fast; a list has to be copied into a tuple before either can
/// hold it, once per output term. That copy is the cost this avoids, and the
/// per-term loop is where this library's marshalling budget goes
/// (`docs/record/python-and-sage-interop.md`).
///
/// Inbound it accepts any sequence, so a caller may hand back what it received
/// or pass the list it already has. That asymmetry is deliberate: the boundary
/// is strict about what it promises and permissive about what it accepts.
#[derive(Clone, Debug, PartialEq, Eq)]
struct Key(Vec<u32>);

impl<'py> IntoPyObject<'py> for Key {
    type Target = pyo3::types::PyTuple;
    type Output = Bound<'py, pyo3::types::PyTuple>;
    type Error = PyErr;
    fn into_pyobject(self, py: Python<'py>) -> Result<Self::Output, Self::Error> {
        pyo3::types::PyTuple::new(py, self.0)
    }
}

impl<'a, 'py> FromPyObject<'a, 'py> for Key {
    type Error = PyErr;
    fn extract(ob: pyo3::Borrowed<'a, 'py, PyAny>) -> Result<Self, PyErr> {
        Ok(Key(ob.extract::<Vec<u32>>()?))
    }
}

impl From<Vec<u32>> for Key {
    fn from(v: Vec<u32>) -> Self {
        Key(v)
    }
}

/// Lets a `Key` stand in wherever the validators take a slice of parts, so
/// [`part_arg`] and [`perm_arg`] read the same on either side of the boundary.
impl std::ops::Deref for Key {
    type Target = [u32];
    fn deref(&self) -> &[u32] {
        &self.0
    }
}

type Terms = Vec<(Key, Coeff)>;
type RatTerms = Vec<(Key, (Coeff, Coeff))>;

/// A partition argument, **validated** rather than repaired.
///
/// Trailing zeros are padding and are dropped: Sage hands over fixed-width
/// lists, and that tolerance is the same one `Perm`'s normal form extends to
/// one-line words. Anything else that is not a partition — parts out of order,
/// a zero between positive parts — is a caller error and raises.
///
/// The distinction matters because the repair is invisible. `Partition::new`
/// sorts, so `[1, 3]` would otherwise reach the mathematics as `[3, 1]` and the
/// caller would get a well-formed answer to a question they had not asked — the
/// plausible wrong value `docs/policies/failure.md` ranks below a crash.
/// Sorting is right for a Rust caller who built the vector from a generator; it
/// is wrong for a foreign caller whose list is data.
fn part_arg(p: &[u32]) -> PyResult<Partition> {
    let end = p.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1);
    let body = &p[..end];
    if let Some(i) = body.iter().position(|&x| x == 0) {
        return Err(PyValueError::new_err(format!(
            "not a partition: {p:?} has a zero part at index {i}, before a positive one; \
             only trailing zeros are padding"
        )));
    }
    if let Some(i) = body.windows(2).position(|w| w[0] < w[1]) {
        return Err(PyValueError::new_err(format!(
            "not a partition: {p:?} increases at index {i} ({} < {}); \
             parts must be weakly decreasing",
            body[i],
            body[i + 1]
        )));
    }
    Ok(Partition::new(body.iter().copied()))
}

/// Every partition in a list argument, validated left to right.
fn parts_arg(ps: &[Vec<u32>]) -> PyResult<Vec<Partition>> {
    ps.iter().map(|p| part_arg(p)).collect()
}

/// Partitions the caller has asked to be compared, which must share a degree.
///
/// For the objects that use this, an off-degree argument is not a zero — it is
/// a question with no referent. `χ^λ(μ)` needs `μ` to index a conjugacy class
/// of `S_{|λ|}`, and `g^ν_{λμ}` needs all three in one `S_n`. The core
/// functions return `0` there as a documented *convention*
/// (`ops.rs`, "unequal degrees pair to zero"). That is the right total
/// behavior for a Rust caller composing them, and the wrong answer to give a
/// foreign caller who mistyped a partition (R11).
///
/// Deliberately **not** applied to `c^λ_{μν}`, `K_{λμ}`, `s_{λ/μ}` or
/// `s_λ(1^n)`, where the zero is a theorem rather than a convention and a
/// caller sweeping a range depends on getting it.
fn same_degree(named: &[(&str, &Partition)]) -> PyResult<()> {
    let (first_name, first) = named[0];
    for (name, p) in &named[1..] {
        if p.size() != first.size() {
            return Err(PyValueError::new_err(format!(
                "these are not all partitions of one integer: {first_name} = {first} has \
                 degree {}, but {name} = {p} has degree {}",
                first.size(),
                p.size()
            )));
        }
    }
    Ok(())
}

/// A ribbon level or quotient index, which the recursions require to be ≥ 1.
///
/// `k = 0` is not a degenerate case with an empty answer — `n % k` divides by
/// zero and the abacus has no runners — so the owning modules assert on it
/// (`a ribbon level needs k ≥ 1`). Caught here, it is a `ValueError` before any
/// of them runs.
fn level_arg(k: u32) -> PyResult<u32> {
    if k < 1 {
        return Err(PyValueError::new_err("a ribbon level needs k >= 1, got 0"));
    }
    Ok(k)
}

/// A shape and level the 128-bit abacus can actually hold.
///
/// A capacity wall rather than a violated precondition, so what R2 asks of it
/// is that a Sage caller meet it as an exception and not as a
/// `PanicException`. The range is [`crate::llt::abacus_reach`]'s expression and
/// not a copy of it — see its doc for why a second copy would drift.
fn abacus_arg(lambda: &Partition, k: u32) -> PyResult<()> {
    let reach = crate::llt::abacus_reach(lambda, k);
    if reach >= crate::llt::ABACUS_REACH_LIMIT {
        return Err(PyValueError::new_err(format!(
            "the shape {lambda} at level {k} needs beta-numbers up to {reach}, \
             past this abacus's {} bits",
            crate::llt::ABACUS_REACH_LIMIT
        )));
    }
    Ok(())
}

/// The same for the whole-degree walk, whose reach reflects every shape of the
/// degree rather than one.
fn abacus_table_arg(n: u32, k: u32) -> PyResult<()> {
    let reach = crate::llt::abacus_reach_table(n, k);
    if reach >= crate::llt::ABACUS_REACH_LIMIT {
        return Err(PyValueError::new_err(format!(
            "every shape of size {} at level {k} needs beta-numbers up to {reach}, \
             past this abacus's {} bits",
            k * n,
            crate::llt::ABACUS_REACH_LIMIT
        )));
    }
    Ok(())
}

/// A cell count a [`crate::llt::SkewTuple`] can hold.
fn cells_arg(cells: usize, what: &str) -> PyResult<()> {
    if cells > crate::llt::MAX_CELLS {
        return Err(PyValueError::new_err(format!(
            "{what} is {cells} cells, past the {} a SkewTuple holds; \
             this is a representation limit, not a mathematical one",
            crate::llt::MAX_CELLS
        )));
    }
    Ok(())
}

// --- the two coefficient rings the boundary runs over ------------------------

/// A ring that can carry a coefficient across the boundary.
///
/// Implemented by exactly two types, which are the two passes: [`Guarded`] (the
/// fixed-width attempt, which may decline an input that does not fit) and
/// `BigInt` (the fallback, which never declines).
trait Boundary: Ring + Sized + ToCoeff {
    fn from_coeff(v: &Coeff) -> Option<Self>;
}

/// The *outbound* half of the boundary on its own.
///
/// Separate from [`Boundary`], which means "one of the two escalation passes",
/// because the `(q,t)` families that carry no ladder still have to emit their
/// coefficients — and they run over plain `i128`, which is emphatically not a
/// pass: it panics at its wall rather than reporting. Keeping the two apart is
/// what stops `i128` from being accepted anywhere an escalating pass is meant.
trait ToCoeff {
    fn to_coeff(&self) -> Coeff;
}

impl ToCoeff for i128 {
    fn to_coeff(&self) -> Coeff {
        Coeff::Small(*self)
    }
}

impl ToCoeff for Guarded {
    /// No allocation, which is what the fast path is for.
    fn to_coeff(&self) -> Coeff {
        Coeff::Small(self.0)
    }
}

impl ToCoeff for BigInt {
    fn to_coeff(&self) -> Coeff {
        Coeff::Big(self.clone())
    }
}

impl Boundary for Guarded {
    fn from_coeff(v: &Coeff) -> Option<Self> {
        v.as_i128().map(Guarded)
    }
}

impl Boundary for BigInt {
    fn from_coeff(v: &Coeff) -> Option<Self> {
        Some(v.to_big())
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

/// A boundary ring that accepts **every** coefficient — the escalation target.
///
/// This is what the `Option` in [`build`] is about: only the fixed-width pass
/// can decline an input, and the slow pass by construction cannot. Saying that
/// in the type system rather than with an `unwrap` at every call site is not
/// tidiness. An `unwrap` on the slow path reads as "this cannot happen", which
/// is a claim nothing checks — and on [`build_schubert`] a malformed
/// permutation makes it false (`docs/policies/failure.md`, R2).
trait Wide: Boundary {
    fn from_coeff_wide(v: &Coeff) -> Self;
}

impl Wide for BigInt {
    fn from_coeff_wide(v: &Coeff) -> Self {
        v.to_big()
    }
}

/// The same guarantee for the rings that carry denominators.
trait WideRat: BoundaryRat {
    fn from_coeff_wide(v: &Coeff) -> Self;
}

impl WideRat for BigRational {
    fn from_coeff_wide(v: &Coeff) -> Self {
        BigRational::from(v.to_big())
    }
}

/// The terms of an element with every partition already validated, so the
/// builders below can decline for one reason only: a coefficient too wide for
/// the fixed-width pass.
type Parsed<'a> = Vec<(Partition, &'a Coeff)>;

/// Validate the partitions **before** either pass runs.
fn terms_arg(t: &Terms) -> PyResult<Parsed<'_>> {
    t.iter().map(|(p, c)| Ok((part_arg(p)?, c))).collect()
}

fn build<C: Boundary, B: SymFn<C>>(terms: &Parsed) -> Option<B> {
    let mut x = B::zero();
    for (p, c) in terms {
        x.add_term(p.clone(), C::from_coeff(c)?);
    }
    Some(x)
}

/// [`build`] over a ring that cannot decline, so there is nothing to unwrap.
fn build_wide<C: Wide, B: SymFn<C>>(terms: &Parsed) -> B {
    let mut x = B::zero();
    for (p, c) in terms {
        x.add_term(p.clone(), C::from_coeff_wide(c));
    }
    x
}

fn build_rat<C: BoundaryRat, B: SymFn<C>>(terms: &Parsed) -> Option<B> {
    let mut x = B::zero();
    for (p, c) in terms {
        x.add_term(p.clone(), C::from_coeff(c)?);
    }
    Some(x)
}

/// [`build_rat`] over a ring that cannot decline.
fn build_rat_wide<C: WideRat, B: SymFn<C>>(terms: &Parsed) -> B {
    let mut x = B::zero();
    for (p, c) in terms {
        x.add_term(p.clone(), C::from_coeff_wide(c));
    }
    x
}

fn dump<C: Boundary, S: SymFn<C>>(x: &S) -> Terms {
    x.terms()
        .iter()
        .map(|(p, c)| (p.parts().to_vec().into(), c.to_coeff()))
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
        out.push((p.parts().to_vec().into(), num));
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

// --- cancellation -----------------------------------------------------------

thread_local! {
    /// The error [`check_python_signals`] took off CPython, held between the
    /// poll that saw it and the entry point that reports it. Thread-local
    /// because the panic it accompanies unwinds one thread, and that is the
    /// thread that will read it back.
    static PENDING: std::cell::RefCell<Option<PyErr>> = const { std::cell::RefCell::new(None) };
}

/// The checker installed into [`crate::interrupt`]: run CPython's pending
/// signal handlers and say whether one of them raised.
///
/// `check_signals` is what a C extension is expected to call to stay
/// responsive, and calling it is the entire fix — CPython has already recorded
/// the signal and is only waiting for a bytecode boundary that a long Rust
/// call never reaches. Reattaching is cheap: an entry point runs holding the
/// GIL, so this re-enters an attachment rather than acquiring one.
///
/// The raised error is kept rather than discarded. A signal handler is
/// arbitrary Python and need not raise `KeyboardInterrupt` at all, so
/// reporting one here would replace whatever the program actually raised.
fn check_python_signals() -> bool {
    Python::attach(|py| match py.check_signals() {
        Ok(()) => false,
        Err(e) => {
            PENDING.with(|p| *p.borrow_mut() = Some(e));
            true
        }
    })
}

/// Run an entry point's body, reporting a cancellation as the exception that
/// caused it.
///
/// Every `#[pyfunction]` in this module wraps its body in this. Uniformity is
/// the point: poll sites are in the kernel, so which entry points can reach
/// one is a property of the mathematics rather than of this file, and an
/// unwrapped body would surface a cancellation as PyO3's `PanicException`
/// instead (`docs/policies/failure.md`, the cancellation row).
fn interruptible<T>(f: impl FnOnce() -> PyResult<T>) -> PyResult<T> {
    match crate::interrupt::catch_interrupt(f) {
        Ok(r) => r,
        Err(crate::interrupt::Interrupted) => Err(PENDING
            .with(|p| p.borrow_mut().take())
            .unwrap_or_else(|| pyo3::exceptions::PyKeyboardInterrupt::new_err("interrupted"))),
    }
}

/// Silence the default panic report for cancellations, and only for those.
///
/// A cancellation is an expected outcome, but it is still a panic, so the
/// default hook would print a Rust backtrace notice to stderr every time a
/// user pressed Ctrl-C. Every other payload goes to the hook this replaced,
/// so a violated contract still reports exactly as before.
fn install_quiet_cancellation_hook() {
    let previous = std::panic::take_hook();
    std::panic::set_hook(Box::new(move |info| {
        if !crate::interrupt::is_interrupt(info.payload()) {
            previous(info);
        }
    }));
}

// --- products ---------------------------------------------------------------

/// Multiply two Schur-basis elements (Littlewood–Richardson).
///
/// Returns the product's `(partition, coefficient)` pairs, ordered
/// lexicographically by parts. Each partition appears once, and no coefficient
/// is zero. An empty list is the zero element, and its products are empty.
///
/// The structure constants are the Littlewood–Richardson numbers `c^ν_{λμ}`,
/// which are nonnegative, so a negative coefficient in the answer came from a
/// negative one in an argument. Cost is driven by the shapes, not by
/// coefficient width, which has no ceiling here.
///
/// ```text
/// >>> symfn.schur_multiply([([2, 1], 1)], [([1], 1)])
/// [((2, 1, 1), 1), ((2, 2), 1), ((3, 1), 1)]
/// ```
///
/// Sage's equivalent is `s(la) * s(mu)`.
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a partition;
/// only trailing zeros are padding.
#[pyfunction]
fn schur_multiply(a: Terms, b: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let (a, b) = (terms_arg(&a)?, terms_arg(&b)?);
        Ok(escalate(
            || {
                let (x, y): (Schur<Guarded>, Schur<Guarded>) = (build(&a)?, build(&b)?);
                Some(dump(&guarded(|| x.mul(&y))?))
            },
            || {
                let (x, y): (Schur<BigInt>, Schur<BigInt>) = (build_wide(&a), build_wide(&b));
                dump(&x.mul(&y))
            },
        ))
    })
}

// --- the Orellana-Zabrocki character basis ----------------------------------

/// Multiply two `st`-basis elements — the Orellana–Zabrocki irreducible
/// character basis `s̃`, whose structure constants **are** the reduced (stable)
/// Kronecker coefficients.
///
/// This is the entry point with the largest measured gap to Sage, because it is
/// where Sage stops: `st[4,3]²` is the largest case Sage still answers, and
/// `st[5,3]²`, `st[6,4]²` and `st[8,5]·st[7,4]` all run here while Sage exceeds
/// 90 s (`docs/record/kronecker.md`).
///
/// ⚠️ `s̃_λ` is **inhomogeneous** — it has components in every degree from 0 to
/// `|λ|` — so unlike every other product on this surface the answer's degree is
/// not the sum of the inputs'. The long implicit first row is what λ omits, so
/// λ indexes a shape of any large size rather than a partition of one `n`.
/// The empty partition is the unit `s̃_∅ = 1` and is a legitimate support both
/// ways, which is what the first term below is.
///
/// ```text
/// >>> symfn.st_multiply([([1], 1)], [([1], 1)])
/// [((), 1), ((1,), 1), ((1, 1), 1), ((2,), 1)]
/// ```
///
/// Sage's equivalent is `SymmetricFunctions(QQ).st()`, and the four terms
/// above are where the reduced product differs from the ordinary Kronecker
/// one, which stays in degree 1.
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a partition.
#[pyfunction]
fn st_multiply(a: Terms, b: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let (a, b) = (terms_arg(&a)?, terms_arg(&b)?);
        Ok(escalate(
            || {
                let (x, y): (St<Guarded>, St<Guarded>) = (build(&a)?, build(&b)?);
                Some(dump(&guarded(|| x.mul(&y))?))
            },
            || {
                let (x, y): (St<BigInt>, St<BigInt>) = (build_wide(&a), build_wide(&b));
                dump(&x.mul(&y))
            },
        ))
    })
}

/// The reduced Kronecker product `s̃_λ · s̃_μ`, as one column.
///
/// The engine's unit of work, and cheaper than [`st_multiply`] on two single
/// terms only in that it skips the bilinear loop; the column is memoized either
/// way.
///
/// Returns `(ν, ḡ^ν_{λμ})` pairs in the element order, and the answer is
/// inhomogeneous for the reason [`st_multiply`] gives. Both arguments index
/// shapes with an implicit long first row, so neither is a partition of a
/// fixed `n` and the empty partition is the unit.
///
/// ```text
/// >>> symfn.reduced_kronecker_product([2], [1])
/// [((1,), 1), ((1, 1), 1), ((2,), 1), ((2, 1), 1), ((3,), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless both arguments are partitions.
#[pyfunction]
fn reduced_kronecker_product(la: Vec<u32>, mu: Vec<u32>) -> PyResult<Terms> {
    interruptible(move || {
        let (l, m) = (part_arg(&la)?, part_arg(&mu)?);
        Ok(escalate(
            || {
                let x: St<Guarded> =
                    guarded(|| crate::reduced_kronecker_product::<Guarded>(&l, &m))?;
                Some(dump(&x))
            },
            || dump(&crate::reduced_kronecker_product::<BigInt>(&l, &m)),
        ))
    })
}

/// One reduced Kronecker coefficient `ḡ^ν_{λμ}`.
///
/// Unlike [`kronecker_coefficient`] the three shapes need **not** share a
/// degree: they index shapes with an implicit long first row, so `ḡ^ν_{λμ}` is
/// defined for any three and off-degree is not a question with no referent.
///
/// Returns one `int`, which is nonnegative and may be zero.
///
/// ```text
/// >>> symfn.reduced_kronecker([2], [1, 1], [1])
/// 1
/// >>> symfn.reduced_kronecker([1], [1], [])
/// 1
/// ```
///
/// The second value is the one that separates this from
/// [`kronecker_coefficient`], which has no answer for shapes of unequal size.
///
/// # Raises
///
/// Raises `ValueError` unless all three arguments are partitions.
#[pyfunction]
fn reduced_kronecker(la: Vec<u32>, mu: Vec<u32>, nu: Vec<u32>) -> PyResult<Coeff> {
    interruptible(move || {
        let (l, m, n) = (part_arg(&la)?, part_arg(&mu)?, part_arg(&nu)?);
        Ok(escalate(
            || {
                let v: Guarded = guarded(|| crate::reduced_kronecker::<Guarded>(&l, &m, &n))?;
                Some(v.to_coeff())
            },
            || crate::reduced_kronecker::<BigInt>(&l, &m, &n).to_coeff(),
        ))
    })
}

/// `s → st`: rewrite a Schur-basis element in the character basis.
///
/// The image is inhomogeneous even when the argument is not, since `s̃_λ` is,
/// so the result carries supports of every size up to the argument's degree.
/// [`st_to_schur`] inverts this exactly.
///
/// ```text
/// >>> symfn.schur_to_st([([2], 1)])
/// [((), 2), ((1,), 2), ((2,), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition.
#[pyfunction]
fn schur_to_st(a: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let a = terms_arg(&a)?;
        Ok(escalate(
            || {
                let x: Schur<Guarded> = build(&a)?;
                Some(dump(&guarded(|| St::from_schur(&x))?))
            },
            || {
                let x: Schur<BigInt> = build_wide(&a);
                dump(&St::from_schur(&x))
            },
        ))
    })
}

/// `st → s`: rewrite a character-basis element in the Schur basis.
///
/// The inverse of [`schur_to_st`]. The value below pins the normalization:
/// `s̃_(2) = s_2 − 2·s_1`, an element of mixed degree rather than the
/// degree-2 one an index-shifted reading would give.
///
/// ```text
/// >>> symfn.st_to_schur([([2], 1)])
/// [((1,), -2), ((2,), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition.
#[pyfunction]
fn st_to_schur(a: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let a = terms_arg(&a)?;
        Ok(escalate(
            || {
                let x: St<Guarded> = build(&a)?;
                Some(dump(&guarded(|| x.to_schur())?))
            },
            || {
                let x: St<BigInt> = build_wide(&a);
                dump(&x.to_schur())
            },
        ))
    })
}

/// Multiply two `ht`-basis elements — the Orellana–Zabrocki **induced trivial**
/// character basis `h̃`, by the double-coset matrix rule.
///
/// Returns `None` rather than raising when the enumeration exceeds its budget.
/// That is a capacity wall and not a caller error. The matrix rule is cheap
/// exactly where the partitions are short: `h̃_{(6,4)}²` is a 2×2 free block,
/// at most 1225 matrices. It is hopeless where they are long, since `(1^10)²`
/// is a 10×10 block with row sums 1, i.e. `11^10`. A caller with another route
/// should be told, not raised at, so this reports the wall the way
/// `try_character` does.
///
/// ⚠️ `h̃_λ` is **inhomogeneous**, like `s̃_λ`.
///
/// ```text
/// >>> symfn.ht_multiply([([1], 1)], [([1], 1)])
/// [((1,), 1), ((1, 1), 1)]
/// >>> symfn.ht_multiply([([1] * 10, 1)], [([1] * 10, 1)]) is None
/// True
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a partition.
/// The capacity wall is not an error and returns `None`.
#[pyfunction]
fn ht_multiply(a: Terms, b: Terms) -> PyResult<Option<Terms>> {
    interruptible(move || {
        let (a, b) = (terms_arg(&a)?, terms_arg(&b)?);
        // Asked before multiplying, which is what makes `Ht::mul`'s panic
        // unreachable from here. The rows are memoized, so the second ask inside
        // the product is a lookup.
        for (l, _) in &a {
            for (m, _) in &b {
                if crate::ht_product_terms(l, m).is_none() {
                    return Ok(None);
                }
            }
        }
        Ok(Some(escalate(
            || {
                let (x, y): (Ht<Guarded>, Ht<Guarded>) = (build(&a)?, build(&b)?);
                Some(dump(&guarded(|| x.mul(&y))?))
            },
            || {
                let (x, y): (Ht<BigInt>, Ht<BigInt>) = (build_wide(&a), build_wide(&b));
                dump(&x.mul(&y))
            },
        )))
    })
}

/// `s → ht`: rewrite a Schur-basis element in the induced trivial character
/// basis.
///
/// Inhomogeneous for the same reason [`schur_to_st`] is, and inverted exactly
/// by [`ht_to_schur`].
///
/// ```text
/// >>> symfn.schur_to_ht([([2], 1)])
/// [((1,), 1), ((2,), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition.
#[pyfunction]
fn schur_to_ht(a: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let a = terms_arg(&a)?;
        Ok(escalate(
            || {
                let x: Schur<Guarded> = build(&a)?;
                Some(dump(&guarded(|| Ht::from_schur(&x))?))
            },
            || {
                let x: Schur<BigInt> = build_wide(&a);
                dump(&Ht::from_schur(&x))
            },
        ))
    })
}

/// `ht → s`: rewrite an induced trivial character element in the Schur basis.
///
/// The inverse of [`schur_to_ht`], and the value below pins the
/// normalization: `h̃_(2) = s_2 − s_1`.
///
/// ```text
/// >>> symfn.ht_to_schur([([2], 1)])
/// [((1,), -1), ((2,), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition.
#[pyfunction]
fn ht_to_schur(a: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let a = terms_arg(&a)?;
        Ok(escalate(
            || {
                let x: Ht<Guarded> = build(&a)?;
                Some(dump(&guarded(|| x.to_schur())?))
            },
            || {
                let x: Ht<BigInt> = build_wide(&a);
                dump(&x.to_schur())
            },
        ))
    })
}

// --- Schubert polynomials ----------------------------------------------------

/// Schubert elements cross the boundary as `(one-line permutation, coeff)`
/// pairs, permutations 1-based, exactly as partitions do — and normalized on
/// entry the same way [`part_arg`] normalizes trailing zeros, so a caller that
/// pads to a fixed `n` gets the same element as one that does not. That padding
/// tolerance is what `Perm`'s normal form is for, and it has to survive
/// the FFI, since Sage hands over fixed-width lists.
type SchubTerms = Vec<(Key, Coeff)>;

/// The same terms with every one-line word already checked to be a permutation,
/// so the builders below can decline for one reason only: a coefficient too
/// wide for the fixed-width pass.
type SchubParsed<'a> = Vec<(Perm, &'a Coeff)>;

/// Validate the one-line words **before** either pass runs.
///
/// A malformed word is a caller error and gets a `ValueError` naming it. An
/// `unwrap` on the escalation path would raise a `PanicException` in Sage
/// instead. That is by definition a bug report and never an interface
/// (`docs/policies/failure.md`, R2).
fn schub_terms(t: &SchubTerms) -> PyResult<SchubParsed<'_>> {
    t.iter().map(|(w, c)| Ok((perm_arg(w)?, c))).collect()
}

fn build_schubert<C: Boundary>(t: &SchubParsed) -> Option<Schubert<C>> {
    let mut out = Schubert::zero();
    for (p, c) in t {
        out.add_term(*p, &C::from_coeff(c)?);
    }
    Some(out)
}

/// [`build_schubert`] over a ring that cannot decline — see [`Wide`].
fn build_schubert_wide<C: Wide>(t: &SchubParsed) -> Schubert<C> {
    let mut out = Schubert::zero();
    for (p, c) in t {
        out.add_term(*p, &C::from_coeff_wide(c));
    }
    out
}

fn dump_schubert<C: Boundary>(f: &Schubert<C>) -> SchubTerms {
    f.terms()
        .iter()
        .map(|(w, c)| {
            (
                w.one_line()
                    .iter()
                    .map(|&x| x as u32)
                    .collect::<Vec<_>>()
                    .into(),
                c.to_coeff(),
            )
        })
        .collect()
}

fn perm_arg(w: &[u32]) -> PyResult<Perm> {
    Perm::new(w.iter().copied())
        .map_err(|e| PyValueError::new_err(format!("not a permutation: {e:?}")))
}

/// A **1-based** variable index, bounded by what a [`Perm`] can hold.
///
/// The operators that take one act by a transposition or a cover scan at
/// position `i`, so they reach `i + 1` points; [`MAX_SUPPORT`] is therefore the
/// ceiling on `i + 1` and not on `i`. Unchecked, `i = 32` asserts inside the
/// cover scan and `i = u32::MAX` overflows the `i + 1` itself — both
/// `PanicException`s from an ordinary integer argument.
fn variable_arg(i: u32) -> PyResult<u32> {
    if i < 1 {
        return Err(PyValueError::new_err("variable index is 1-based"));
    }
    if i as usize >= MAX_SUPPORT {
        return Err(PyValueError::new_err(format!(
            "variable index {i} is past this permutation representation's \
             ceiling of {} points",
            MAX_SUPPORT
        )));
    }
    Ok(i)
}

/// Every term of a Schubert element lies in `S_n`.
///
/// The pairing reads the coefficient of `w0(n)`, so a `w` outside `S_n` cannot
/// contribute and the answer comes back `0` — indistinguishable from the
/// honest zero of two classes whose degrees do not complement. That is exactly
/// the confusion this signature's explicit `n` exists to prevent: Symmetrica's
/// `scalarproduct_schubert` infers `n` from whatever padding it finds, so the
/// same inputs give different answers, and accepting an out-of-range `w` here
/// would reintroduce the same trap one level up.
fn in_flag(t: &SchubParsed, n: u32, which: &str) -> PyResult<()> {
    if let Some((w, _)) = t.iter().find(|(w, _)| w.support_len() > n) {
        return Err(PyValueError::new_err(format!(
            "{which} has a term {w} outside S_{n}: it moves {} points, so the \
             pairing on Fl({n}) is not defined for it",
            w.support_len()
        )));
    }
    Ok(())
}

/// The number of variables `n` in `H*(Fl(n))`, bounded the same way.
///
/// [`schubert_pairing`] builds `w0 = n, n−1, …, 1`, which needs `n` points.
fn rank_arg(n: u32) -> PyResult<u32> {
    if n as usize > MAX_SUPPORT {
        return Err(PyValueError::new_err(format!(
            "Fl({n}) needs {n} points, past this permutation representation's \
             ceiling of {}",
            MAX_SUPPORT
        )));
    }
    Ok(n)
}

/// An exponent vector that [`Schubert::from_polynomial`] will read as a Lehmer
/// code.
///
/// `Perm::from_code` sizes the permutation as `max(|c|, maxᵢ(i + 1 + cᵢ))`,
/// which is generally *more* than `|c|` — bounding by the length alone is the
/// mistake that function's own doc warns about, so the check here is the same
/// expression it uses.
fn code_arg(e: &[u32]) -> PyResult<()> {
    let end = e.iter().rposition(|&x| x != 0).map_or(0, |i| i + 1);
    let mut m = end as u32;
    for (i, &c) in e.iter().enumerate() {
        m = m.max(i as u32 + 1 + c);
    }
    if m as usize > MAX_SUPPORT {
        return Err(PyValueError::new_err(format!(
            "the exponent vector {e:?} is the code of a permutation of {m} points, \
             past this representation's ceiling of {}",
            MAX_SUPPORT
        )));
    }
    Ok(())
}

/// Multiply two Schubert polynomials.
///
/// This is the entry point that replaces Symmetrica's
/// `mult_schubert_schubert`, which is what Sage routes through today.
///
/// Goes through [`Schubert::mul`] rather than naming an engine, so the wheel
/// tracks whichever engine the crate considers best. An engine named here is
/// one this boundary can drift from, leaving a Rust caller and a Sage caller
/// on engines orders of magnitude apart (`docs/record/schubert.md`).
///
/// Arguments and result are `(one-line word, coefficient)` lists, the words
/// **1-based** and padded however the caller likes. The result is ordered
/// lexicographically by word, each word appears once, and no coefficient is
/// zero. Every structure constant is nonnegative.
///
/// ```text
/// >>> symfn.schubert_multiply([([1, 3, 2], 1)], [([2, 1, 3], 1)])
/// [((2, 3, 1), 1), ((3, 1, 2), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a permutation
/// in one-line notation.
#[pyfunction]
fn schubert_multiply(a: SchubTerms, b: SchubTerms) -> PyResult<SchubTerms> {
    interruptible(move || {
        let (a, b) = (schub_terms(&a)?, schub_terms(&b)?);
        Ok(escalate(
            || {
                let (x, y): (Schubert<Guarded>, Schubert<Guarded>) =
                    (build_schubert(&a)?, build_schubert(&b)?);
                Some(dump_schubert(&guarded(|| x.mul(&y))?))
            },
            || {
                let (x, y): (Schubert<BigInt>, Schubert<BigInt>) =
                    (build_schubert_wide(&a), build_schubert_wide(&b));
                dump_schubert(&x.mul(&y))
            },
        ))
    })
}

/// `x_i · f`, the signed Monk rule. **1-based**, unlike Symmetrica's
/// `mult_schubert_variable`, which is 0-based while its own
/// `divdiff_schubert` is 1-based. One convention, stated.
///
/// The value below is the one that separates the two: `i = 1` multiplies by
/// `x_1`, so `x_1 · S_{132}` is what comes back, not `x_2 · S_{132}`.
///
/// ```text
/// >>> symfn.schubert_multiply_variable([([1, 3, 2], 1)], 1)
/// [((2, 3, 1), 1), ((3, 1, 2), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` if a term is not a permutation, if `i` is zero, or if
/// `i` is past the permutation representation's ceiling.
#[pyfunction]
fn schubert_multiply_variable(a: SchubTerms, i: u32) -> PyResult<SchubTerms> {
    interruptible(move || {
        let (a, i) = (schub_terms(&a)?, variable_arg(i)?);
        Ok(escalate(
            || {
                let x: Schubert<Guarded> = build_schubert(&a)?;
                Some(dump_schubert(&guarded(|| x.mul_variable(i))?))
            },
            || {
                let x: Schubert<BigInt> = build_schubert_wide(&a);
                dump_schubert(&x.mul_variable(i))
            },
        ))
    })
}

/// `∂_i f` on the Schubert basis, 1-based.
///
/// `∂_i S_w = S_{w·s_i}` when `ℓ(w·s_i) < ℓ(w)`, and 0 otherwise, so the
/// answer is empty exactly when every term is killed. Terms come back
/// ordered lexicographically by word.
///
/// ```text
/// >>> symfn.schubert_divided_difference([([3, 1, 2], 1)], 1)
/// [((1, 3, 2), 1)]
/// >>> symfn.schubert_divided_difference([([3, 1, 2], 1)], 2)
/// []
/// ```
///
/// # Raises
///
/// Raises `ValueError` if a term is not a permutation, if `i` is zero, or if
/// `i` is past the permutation representation's ceiling.
#[pyfunction]
fn schubert_divided_difference(a: SchubTerms, i: u32) -> PyResult<SchubTerms> {
    interruptible(move || {
        let (a, i) = (schub_terms(&a)?, variable_arg(i)?);
        Ok(escalate(
            || {
                let x: Schubert<Guarded> = build_schubert(&a)?;
                Some(dump_schubert(&guarded(|| x.divided_difference(i))?))
            },
            || {
                let x: Schubert<BigInt> = build_schubert_wide(&a);
                dump_schubert(&x.divided_difference(i))
            },
        ))
    })
}

/// `∂_w f`, composing along a reduced word of `w`.
///
/// `w` is a permutation in one-line notation, not a word in the generators.
/// The result does not depend on which reduced word is chosen, because the
/// `∂_i` satisfy the braid relations. Ordered lexicographically by word.
///
/// ```text
/// >>> symfn.schubert_divided_difference_perm([([1, 4, 3, 2], 1)], [1, 3, 2])
/// [((1, 3, 4, 2), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term and `w` are permutations.
#[pyfunction]
fn schubert_divided_difference_perm(a: SchubTerms, w: Vec<u32>) -> PyResult<SchubTerms> {
    interruptible(move || {
        let (a, p) = (schub_terms(&a)?, perm_arg(&w)?);
        Ok(escalate(
            || {
                let x: Schubert<Guarded> = build_schubert(&a)?;
                Some(dump_schubert(&guarded(|| x.divided_difference_perm(&p))?))
            },
            || {
                let x: Schubert<BigInt> = build_schubert_wide(&a);
                dump_schubert(&x.divided_difference_perm(&p))
            },
        ))
    })
}

/// Expand into monomials: `(exponent vector, coefficient)` pairs.
///
/// ⚠️ The output is `S_w(1,…,1)` terms, which grows super-exponentially —
/// 84 084 monomials for one random S₁₂ element of length 33. Callers wanting
/// a size estimate first should ask [`schubert_dimension`], which is cheap.
///
/// An exponent vector is a list of nonnegative integers, index `i` carrying
/// the exponent of `x_{i+1}`, with trailing zeros dropped. Pairs come back
/// ordered lexicographically by exponent vector, each vector once, no zero
/// coefficients.
///
/// ```text
/// >>> symfn.schubert_expand([([1, 3, 2], 1)])
/// [((0, 1), 1), ((1,), 1)]
/// ```
///
/// So `S_{132} = x_1 + x_2`, and the `(0, 1)` term is `x_2` — the reading
/// that a 1-based exponent vector would not give.
///
/// # Raises
///
/// Raises `ValueError` unless every term is a permutation.
#[pyfunction]
fn schubert_expand(a: SchubTerms) -> PyResult<Terms> {
    interruptible(move || {
        let a = schub_terms(&a)?;
        Ok(escalate(
            || {
                let x: Schubert<Guarded> = build_schubert(&a)?;
                let e = guarded(|| x.expand())?;
                Some(
                    e.into_iter()
                        .map(|(v, c)| (v.into(), c.to_coeff()))
                        .collect(),
                )
            },
            || {
                let x: Schubert<BigInt> = build_schubert_wide(&a);
                x.expand()
                    .into_iter()
                    .map(|(v, c)| (v.into(), c.to_coeff()))
                    .collect()
            },
        ))
    })
}

/// Write a polynomial in the Schubert basis (the greedy triangular peel).
///
/// The argument is `(exponent vector, coefficient)` pairs in
/// [`schubert_expand`]'s encoding, and the two are inverse on any polynomial
/// that is a combination of Schubert polynomials. Result ordered
/// lexicographically by word.
///
/// ```text
/// >>> symfn.polynomial_to_schubert([([2], 1)])
/// [((3, 1, 2), 1)]
/// ```
///
/// So `x_1² = S_{312}`, which the 1-based reading of the exponent vector
/// would report as `S` of a different permutation.
///
/// # Raises
///
/// Raises `ValueError` if an exponent vector is the code of a permutation
/// past the representation's ceiling.
#[pyfunction]
fn polynomial_to_schubert(terms: Terms) -> PyResult<SchubTerms> {
    interruptible(move || {
        for (e, _) in &terms {
            code_arg(e)?;
        }
        Ok(escalate(
            || {
                let t: Vec<(Vec<u32>, Guarded)> = terms
                    .iter()
                    .map(|(e, c)| Guarded::from_coeff(c).map(|g| (e.to_vec(), g)))
                    .collect::<Option<_>>()?;
                Some(dump_schubert(&guarded(|| Schubert::from_polynomial(&t))?))
            },
            || {
                let t: Vec<(Vec<u32>, BigInt)> = terms
                    .iter()
                    .map(|(e, c)| (e.to_vec(), BigInt::from_coeff_wide(c)))
                    .collect();
                dump_schubert(&Schubert::from_polynomial(&t))
            },
        ))
    })
}

/// The Poincaré pairing on `H*(Fl(n))`.
///
/// `n` is **explicit**. Symmetrica's `scalarproduct_schubert` reads it off
/// however long the stored vectors happen to be, so the same mathematical
/// inputs give different answers depending on prior padding.
///
/// Returns one `int`: the coefficient of `S_{w0(n)}` in the product, which is
/// zero unless the two degrees sum to `ℓ(w0) = n(n−1)/2`.
///
/// ```text
/// >>> symfn.schubert_pairing([([1, 3, 2], 1)], [([3, 1, 2], 1)], 3)
/// 1
/// >>> symfn.schubert_pairing([([2, 1, 3], 1)], [([2, 1, 3], 1)], 3)
/// 0
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments lies in `S_n`. A
/// `w` outside it cannot contribute to the coefficient of `w0(n)`, so the
/// answer would be a `0` that no caller could tell from an honest one — which
/// is the confusion the explicit `n` exists to remove. Also raises if `n` is
/// past the permutation representation's ceiling.
#[pyfunction]
fn schubert_pairing(a: SchubTerms, b: SchubTerms, n: u32) -> PyResult<Coeff> {
    interruptible(move || {
        let (a, b, n) = (schub_terms(&a)?, schub_terms(&b)?, rank_arg(n)?);
        in_flag(&a, n, "a")?;
        in_flag(&b, n, "b")?;
        Ok(escalate(
            || {
                let (x, y): (Schubert<Guarded>, Schubert<Guarded>) =
                    (build_schubert(&a)?, build_schubert(&b)?);
                Some(guarded(|| x.pairing(&y, n))?.to_coeff())
            },
            || {
                let (x, y): (Schubert<BigInt>, Schubert<BigInt>) =
                    (build_schubert_wide(&a), build_schubert_wide(&b));
                x.pairing(&y, n).to_coeff()
            },
        ))
    })
}

/// `∂_{w0(n)}(a·b)`, Symmetrica's `scalarproduct_schubert` — what Sage exposes
/// as `SchubertPolynomial.scalar_product`.
///
/// ⚠️ **It returns a Schubert polynomial, not a scalar.** The two operations
/// meet at one term: the coefficient of `S_id` here is `schubert_pairing`, and
/// everything of higher degree survives it. Reading the name as the pairing is
/// the error `docs/record/schubert.md` records.
///
/// `n` is **explicit**, as it is for `schubert_pairing`; the incumbent reads it
/// off however long its stored vectors happen to be. What a Sage caller reaches
/// is `n` = the longest one-line form among the terms of both arguments, since
/// Sage strips trailing fixed points first. Terms come back ordered
/// lexicographically by word, each word once, no zero coefficients.
///
/// Every surviving term drops `n(n−1)/2` in degree, so the answer is empty
/// whenever that exceeds the product's degree; at `n ≤ 1` the operator is empty
/// and the product comes back unchanged.
///
/// ```text
/// >>> symfn.schubert_scalar_product([([2, 1], 1)], [([2, 1], 1)], 2)
/// [((1, 3, 2), 1)]
/// >>> symfn.schubert_scalar_product([([2, 1], 1)], [([2, 1], 1)], 3)
/// []
/// ```
///
/// `S_{21}·S_{21} = x_1² = S_{312}`, and the single pass at `n = 2` leaves
/// `S_{132}` — the value that separates this from the pairing, which answers 0
/// on the same input.
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a permutation, or
/// if `n` is past the permutation representation's ceiling. Unlike
/// `schubert_pairing` it accepts terms outside `S_n`: they contribute honestly
/// here, rather than to a zero no caller could read.
#[pyfunction]
fn schubert_scalar_product(a: SchubTerms, b: SchubTerms, n: u32) -> PyResult<SchubTerms> {
    interruptible(move || {
        let (a, b, n) = (schub_terms(&a)?, schub_terms(&b)?, rank_arg(n)?);
        Ok(escalate(
            || {
                let (x, y): (Schubert<Guarded>, Schubert<Guarded>) =
                    (build_schubert(&a)?, build_schubert(&b)?);
                Some(dump_schubert(&guarded(|| x.scalar_product(&y, n))?))
            },
            || {
                let (x, y): (Schubert<BigInt>, Schubert<BigInt>) =
                    (build_schubert_wide(&a), build_schubert_wide(&b));
                dump_schubert(&x.scalar_product(&y, n))
            },
        ))
    })
}

/// `S_w(1,…,1)`: the number of pipe dreams, i.e. the size `schubert_expand`
/// would produce. Cheap — it never builds the expansion.
///
/// Returns one `int`, at least 1 for any permutation, since the identity has
/// the single monomial 1.
///
/// ```text
/// >>> symfn.schubert_dimension([1, 4, 3, 2])
/// 5
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless `w` is a permutation.
#[pyfunction]
fn schubert_dimension(w: Vec<u32>) -> PyResult<u128> {
    interruptible(move || Ok(crate::schubert::dimension(&perm_arg(&w)?)))
}

/// A single structure constant `c^w_{uv}`, **without building the product**.
///
/// No other package offers this (`docs/record/schubert.md`). `S_u · S_v` can
/// have a monomial mass of 4.3×10¹⁶, an answer that fits on no machine, while
/// one of its coefficients stays reachable. Positivity searches and
/// rule-hunting want particular constants, not the whole expansion.
///
/// Returns 0 immediately unless `ℓ(w) = ℓ(u)+ℓ(v)` and `u ≤ w`, `v ≤ w` in
/// Bruhat order. Zero is an answer here, not a refusal — a caller sweeping a
/// range of `w` depends on getting it.
///
/// ```text
/// >>> symfn.schubert_coefficient([2, 1, 3], [2, 1, 3], [3, 1, 2])
/// 1
/// >>> symfn.schubert_coefficient([2, 1, 3], [2, 1, 3], [1, 3, 2])
/// 0
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless all three arguments are permutations.
#[pyfunction]
fn schubert_coefficient(u: Vec<u32>, v: Vec<u32>, w: Vec<u32>) -> PyResult<Coeff> {
    interruptible(move || {
        let (pu, pv, pw) = (perm_arg(&u)?, perm_arg(&v)?, perm_arg(&w)?);
        Ok(escalate(
            || {
                guarded(|| crate::schubert::schubert_coeff::<Guarded>(&pu, &pv, &pw))
                    .map(|c| c.to_coeff())
            },
            || crate::schubert::schubert_coeff::<BigInt>(&pu, &pv, &pw).to_coeff(),
        ))
    })
}

/// The **Stanley symmetric function** `F_w` in the Schur basis.
///
/// ⚠️ Symmetrica exports this as `t_SCHUBERT_SCHUR` and Sage inherits the
/// name, but it is not "Schubert → Schur": `F_w = S_w` only in the stable
/// range. `newtrans([2,1,4,3]) = s₂ + s₁₁` while `S_{2143}` is not symmetric
/// at all. The name here says what it computes; the replacement for
/// Symmetrica's `newtrans` is this function, under an honest label.
///
/// Returns `(partition, coefficient)` pairs in the element order. The
/// coefficients are nonnegative, since `F_w` is Schur-positive.
///
/// ```text
/// >>> symfn.schubert_to_stanley_schur([2, 1, 4, 3])
/// [((1, 1), 1), ((2,), 1)]
/// ```
///
/// That is `s_2 + s_{11}`, and `S_{2143}` itself is not symmetric — the value
/// that separates this from a Schubert-to-Schur reading of the name.
///
/// # Raises
///
/// Raises `ValueError` unless `w` is a permutation.
#[pyfunction]
fn schubert_to_stanley_schur(w: Vec<u32>) -> PyResult<Terms> {
    interruptible(move || {
        let p = perm_arg(&w)?;
        Ok(escalate(
            || {
                let f = guarded(|| crate::schubert::stanley::<Guarded>(&p))?;
                Some(dump(&f))
            },
            || dump(&crate::schubert::stanley::<BigInt>(&p)),
        ))
    })
}

/// The product's total monomial mass `S_u(1,…,1)·S_v(1,…,1)`, a count.
///
/// Exposed rather than enforced: it is the only cheap quantity
/// that flags an out-of-family pair — 4.3×10¹⁶ for the one product no engine
/// completes, against 4.4×10¹² for everything else on the ladder — but it is
/// not a runtime predictor, so refusing on it would be guesswork. Callers who
/// want to know before committing can ask.
///
/// Returns one `int`, the product of the two dimensions, saturating rather
/// than wrapping if it exceeds what a 128-bit count holds.
///
/// ```text
/// >>> symfn.schubert_monomial_mass([1, 3, 2], [2, 1, 3])
/// 2
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless both arguments are permutations.
#[pyfunction]
fn schubert_monomial_mass(u: Vec<u32>, v: Vec<u32>) -> PyResult<u128> {
    interruptible(move || {
        let (pu, pv) = (perm_arg(&u)?, perm_arg(&v)?);
        Ok(crate::schubert::dimension(&pu).saturating_mul(crate::schubert::dimension(&pv)))
    })
}

/// Drop every memo cache.
///
/// Exposed for benchmarking rather than for normal use: the caches are
/// referentially transparent, so clearing them cannot change a result, only a
/// timing. A comparison that reuses inputs measures the cache on the second
/// call and not the algorithm — which is exactly the trap
/// `scripts/compare_sage.py` documents on Sage's side.
///
/// Takes no arguments, returns `None`, and raises nothing.
///
/// ```text
/// >>> symfn.clear_caches() is None
/// True
/// ```
#[pyfunction]
fn clear_caches() -> PyResult<()> {
    interruptible(move || {
        crate::memo::clear_caches();
        Ok(())
    })
}

/// A single Littlewood–Richardson coefficient c^λ_{μν}.
///
/// Deliberately [`NaiveLr`] and not [`AutoLr`](crate::strip_lr::AutoLr), which
/// looks backwards — `NaiveLr` is the naive *reference* backend — but is what
/// the measurements say. A targeted backtrack costs roughly the coefficient's
/// own size, while `AutoLr` builds a whole expansion and indexes into it, so
/// the naive search wins whenever the coefficient is small — the
/// overwhelmingly common case — and loses only when it is large, because it
/// then enumerates that many tableaux
/// (`docs/record/littlewood-richardson.md`). Whole *products* are a different
/// question and go through `AutoLr` (see `schur_multiply`).
///
/// **Zero is an answer here, not a refusal.** `c^λ_{μν} = 0` whenever
/// `|λ| ≠ |μ| + |ν|` or λ fails to contain a factor, and that is a theorem
/// rather than a convention standing in for a malformed question — a caller
/// sweeping a range of λ depends on getting it. Contrast
/// [`character_value`] and [`kronecker_coefficient`], where an off-degree
/// argument has no referent at all and raises.
///
/// Returns one nonnegative `int`.
///
/// ```text
/// >>> symfn.lr_coefficient([3, 2, 1], [2, 1], [2, 1])
/// 2
/// >>> symfn.lr_coefficient([4], [2, 1], [2, 1])
/// 0
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless all three arguments are partitions.
#[pyfunction]
fn lr_coefficient(la: Vec<u32>, mu: Vec<u32>, nu: Vec<u32>) -> PyResult<u128> {
    interruptible(move || Ok(NaiveLr.lr_coeff(&part_arg(&la)?, &part_arg(&mu)?, &part_arg(&nu)?)))
}

// --- conversions out of Schur ----------------------------------------------

/// Re-read this module's **own** output as input again.
///
/// The entry points that compose two conversions — [`expand_alphabet`],
/// [`convert_indexed`] — feed one conversion's `Terms` to the next. Those
/// partitions came out of [`Partition::parts`] and so are valid by
/// construction; validating them a second time would be checking this crate
/// rather than the caller, which is what [`terms_arg`] is for.
fn relay(t: &Terms) -> Parsed<'_> {
    t.iter()
        .map(|(p, c)| (Partition::new(p.iter().copied()), c))
        .collect()
}

/// A conversion `Schur -> $basis`, run fixed-width first and re-run exactly if
/// that overflows.
///
/// Generates the inner conversion over already-validated terms alongside the
/// `#[pyfunction]`, so the composing entry points can chain without either
/// re-validating or bypassing validation. The entry point's own doc is passed
/// in, because each pair states a different convention.
macro_rules! out_of_schur {
    ($(#[$doc:meta])* $name:ident, $inner:ident, $basis:ident) => {
        fn $inner(a: &Parsed) -> Terms {
            escalate(
                || {
                    let s: Schur<Guarded> = build(a)?;
                    Some(dump(&guarded(|| $basis::<Guarded>::from_schur(&s))?))
                },
                || {
                    let s: Schur<BigInt> = build_wide(a);
                    dump(&$basis::<BigInt>::from_schur(&s))
                },
            )
        }

        $(#[$doc])*
        #[pyfunction]
        fn $name(a: Terms) -> PyResult<Terms> {
            Ok($inner(&terms_arg(&a)?))
        }
    };
}

out_of_schur!(
    /// `s → h`: rewrite a Schur-basis element in the complete homogeneous
    /// basis.
    ///
    /// Coefficients are integers and may be negative; the inverse is
    /// [`homogeneous_to_schur`].
    ///
    /// ```text
    /// >>> symfn.schur_to_homogeneous([([2, 1], 1)])
    /// [((2, 1), 1), ((3,), -1)]
    /// ```
    ///
    /// # Raises
    ///
    /// Raises `ValueError` unless every term is a partition.
    schur_to_homogeneous,
    s_to_h,
    Homogeneous
);
out_of_schur!(
    /// `s → e`: rewrite a Schur-basis element in the elementary basis.
    ///
    /// Coefficients are integers and may be negative; the inverse is
    /// [`elementary_to_schur`].
    ///
    /// ```text
    /// >>> symfn.schur_to_elementary([([2, 1], 1)])
    /// [((2, 1), 1), ((3,), -1)]
    /// ```
    ///
    /// # Raises
    ///
    /// Raises `ValueError` unless every term is a partition.
    schur_to_elementary,
    s_to_e,
    Elementary
);
out_of_schur!(
    /// `s → m`: rewrite a Schur-basis element in the monomial basis.
    ///
    /// The coefficients are the Kostka numbers `K_{λμ}`, so they are
    /// nonnegative when the argument is. The inverse is
    /// [`monomial_to_schur`].
    ///
    /// ```text
    /// >>> symfn.schur_to_monomial([([2, 1], 1)])
    /// [((1, 1, 1), 2), ((2, 1), 1)]
    /// ```
    ///
    /// # Raises
    ///
    /// Raises `ValueError` unless every term is a partition.
    schur_to_monomial,
    s_to_m,
    Monomial
);
out_of_schur!(
    /// `s → f`: rewrite a Schur-basis element in the forgotten basis.
    ///
    /// The forgotten basis is `ω(m)`, so this is [`schur_to_monomial`] of the
    /// conjugated argument; the inverse is [`forgotten_to_schur`].
    ///
    /// ```text
    /// >>> symfn.schur_to_forgotten([([2, 1], 1)])
    /// [((1, 1, 1), 2), ((2, 1), 1)]
    /// ```
    ///
    /// # Raises
    ///
    /// Raises `ValueError` unless every term is a partition.
    schur_to_forgotten,
    s_to_f,
    Forgotten
);

/// A power-sum element on its way out: rational coefficients as
/// `(numerator, denominator)`.
fn split<C: BoundaryRat>(p: &PowerSum<C>) -> RatTerms {
    p.terms()
        .iter()
        .map(|(part, c)| (part.parts().to_vec().into(), c.split()))
        .collect()
}

/// A conversion **into** the power-sum basis, which is the one family of
/// targets whose coefficients divide.
///
/// h → p and e → p are direct (`crate::convert`); everything else composes
/// through Schur, which is where its own conversion already goes.
macro_rules! into_power {
    ($inner:ident, $basis:ident) => {
        fn $inner(a: &Parsed) -> RatTerms {
            escalate(
                || {
                    let x: $basis<GuardedRat> = build_rat(a)?;
                    Some(split(&guarded(|| {
                        convert::<GuardedRat, _, PowerSum<GuardedRat>>(&x)
                    })?))
                },
                || {
                    let x: $basis<BigRational> = build_rat_wide(a);
                    split(&convert::<BigRational, _, PowerSum<BigRational>>(&x))
                },
            )
        }
    };
}

into_power!(s_to_p, Schur);
into_power!(h_to_p, Homogeneous);
into_power!(e_to_p, Elementary);

/// A conversion into the power-sum basis, from any basis
/// [`convert_indexed`] names.
///
/// Coefficients are rational and come back as `(numerator, denominator)`; every
/// other conversion this module exposes lands in ℤ, which is why this one has
/// an entry point of its own rather than a `dst` on
/// [`convert_terms`](convert_terms).
///
/// `src` is a basis name or one-letter code — `"Schur"` or `"s"`,
/// `"homogeneous"` or `"h"`, `"elementary"` or `"e"`, `"powersum"` or `"p"`,
/// `"monomial"` or `"m"`, `"forgotten"` or `"f"`; the same six spellings every
/// basis argument at this boundary accepts. `"powersum"` is the identity, and
/// is accepted so that a caller dispatching on a basis name does not need a
/// special case for it.
///
/// Returns `(partition, (numerator, denominator))` triples ordered
/// lexicographically by partition, with no zero terms. The fraction is in
/// lowest terms and the denominator is positive.
///
/// ```text
/// >>> symfn.to_power([([2], 1)], "Schur")
/// [((1, 1), (1, 2)), ((2,), (1, 2))]
/// ```
///
/// That is `s_2 = (p_11 + p_2)/2`, and the halves are why this entry point
/// exists apart from [`convert_terms`], which lands in ℤ.
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition and `src` is one of
/// the spellings above.
#[pyfunction]
fn to_power(a: Terms, src: &str) -> PyResult<RatTerms> {
    interruptible(move || {
        let a = terms_arg(&a)?;
        Ok(match Basis::parse(src)? {
            Basis::Schur => s_to_p(&a),
            Basis::Homogeneous => h_to_p(&a),
            Basis::Elementary => e_to_p(&a),
            Basis::PowerSum => a
                .iter()
                .map(|(p, c)| (p.parts().to_vec().into(), ((*c).clone(), Coeff::Small(1))))
                .collect(),
            Basis::Monomial => s_to_p(&relay(&m_to_s(&a))),
            Basis::Forgotten => s_to_p(&relay(&f_to_s(&a))),
        })
    })
}

// --- conversions into Schur -------------------------------------------------

macro_rules! into_schur {
    ($(#[$doc:meta])* $name:ident, $inner:ident, $basis:ident) => {
        fn $inner(a: &Parsed) -> Terms {
            escalate(
                || {
                    let x: $basis<Guarded> = build(a)?;
                    Some(dump(&guarded(|| x.to_schur())?))
                },
                || {
                    let x: $basis<BigInt> = build_wide(a);
                    dump(&x.to_schur())
                },
            )
        }

        $(#[$doc])*
        #[pyfunction]
        fn $name(a: Terms) -> PyResult<Terms> {
            Ok($inner(&terms_arg(&a)?))
        }
    };
}

into_schur!(
    /// `h → s`: rewrite a homogeneous-basis element in the Schur basis.
    ///
    /// The coefficients are the inverse Kostka numbers and may be negative.
    /// The inverse is [`schur_to_homogeneous`].
    ///
    /// ```text
    /// >>> symfn.homogeneous_to_schur([([2, 1], 1)])
    /// [((2, 1), 1), ((3,), 1)]
    /// ```
    ///
    /// # Raises
    ///
    /// Raises `ValueError` unless every term is a partition.
    homogeneous_to_schur,
    h_to_s,
    Homogeneous
);
into_schur!(
    /// `e → s`: rewrite an elementary-basis element in the Schur basis.
    ///
    /// The inverse of [`schur_to_elementary`]. The value below is the
    /// conjugate of [`homogeneous_to_schur`]'s, which is what distinguishes
    /// the two.
    ///
    /// ```text
    /// >>> symfn.elementary_to_schur([([2, 1], 1)])
    /// [((1, 1, 1), 1), ((2, 1), 1)]
    /// ```
    ///
    /// # Raises
    ///
    /// Raises `ValueError` unless every term is a partition.
    elementary_to_schur,
    e_to_s,
    Elementary
);
into_schur!(
    /// `m → s`: rewrite a monomial-basis element in the Schur basis.
    ///
    /// The coefficients are the inverse Kostka numbers and may be negative;
    /// the inverse is [`schur_to_monomial`].
    ///
    /// ```text
    /// >>> symfn.monomial_to_schur([([2, 1], 1)])
    /// [((1, 1, 1), -2), ((2, 1), 1)]
    /// ```
    ///
    /// # Raises
    ///
    /// Raises `ValueError` unless every term is a partition.
    monomial_to_schur,
    m_to_s,
    Monomial
);
into_schur!(
    /// `p → s`: rewrite a power-sum element in the Schur basis.
    ///
    /// The coefficients are the irreducible characters `χ^λ(μ)`, so they are
    /// integers, and this direction never divides — [`to_power`] is the one
    /// that does, which is why it returns fractions and this does not.
    ///
    /// ```text
    /// >>> symfn.power_to_schur([([2, 1], 1)])
    /// [((1, 1, 1), -1), ((3,), 1)]
    /// ```
    ///
    /// # Raises
    ///
    /// Raises `ValueError` unless every term is a partition.
    power_to_schur,
    p_to_s,
    PowerSum
);
into_schur!(
    /// `f → s`: rewrite a forgotten-basis element in the Schur basis.
    ///
    /// The inverse of [`schur_to_forgotten`].
    ///
    /// ```text
    /// >>> symfn.forgotten_to_schur([([2, 1], 1)])
    /// [((2, 1), 1), ((3,), -2)]
    /// ```
    ///
    /// # Raises
    ///
    /// Raises `ValueError` unless every term is a partition.
    forgotten_to_schur,
    f_to_s,
    Forgotten
);

/// A conversion between two multiplicative bases that **skips the Schur hub**.
///
/// The hub is not merely a longer road for these pairs: p_λ with a handful of
/// terms becomes a Schur element with p(n) of them, and the contraction that
/// follows is p(n) determinants. [`crate::convert`] picks the direct rule on
/// its own; naming the pair here is what lets the boundary reach it in one
/// call instead of composing two.
macro_rules! direct_route {
    ($inner:ident, $from:ident, $to:ident) => {
        fn $inner(a: &Parsed) -> Terms {
            escalate(
                || {
                    let x: $from<Guarded> = build(a)?;
                    Some(dump(&guarded(|| convert::<Guarded, _, $to<Guarded>>(&x))?))
                },
                || {
                    let x: $from<BigInt> = build_wide(a);
                    dump(&convert::<BigInt, _, $to<BigInt>>(&x))
                },
            )
        }
    };
}

direct_route!(p_to_h, PowerSum, Homogeneous);
direct_route!(p_to_e, PowerSum, Elementary);
direct_route!(h_to_e, Homogeneous, Elementary);
direct_route!(e_to_h, Elementary, Homogeneous);

/// Plethysm f[g] of two Schur-basis elements.
///
/// The answer comes back in the Schur basis with integer coefficients — a
/// denominator surviving the computation would be a bug, and is reported
/// rather than truncated.
///
/// **Two routes run underneath, chosen from the shapes.** One expands in the
/// power-sum basis and converts back; the other stays in the Schur basis and
/// never converts. Which is cheaper depends on both arguments and the spread
/// is large in both directions, so the choice is made per call
/// (`docs/record/plethysm.md`). What it buys: `s_6[s_6]` in 0.25s against 28s,
/// and `s_3[s_{10,10}]` at degree 60 — 10,198 terms — in 0.75s, where the
/// conversion route does not finish. Nothing about the call changes — same
/// arguments, same answer, same basis — so this is a note about which inputs
/// are cheap, not about the interface.
///
/// Plethysm is not commutative, which the two values below separate.
///
/// ```text
/// >>> symfn.plethysm([([2], 1)], [([1, 1], 1)])
/// [((1, 1, 1, 1), 1), ((2, 2), 1)]
/// >>> symfn.plethysm([([1, 1], 1)], [([2], 1)])
/// [((3, 1), 1)]
/// ```
///
/// Sage writes the same operation as `f(g)` on symmetric functions.
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a partition,
/// and if the computation produces a non-integral coefficient.
#[pyfunction]
fn plethysm(f: Terms, g: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let (f, g) = (terms_arg(&f)?, terms_arg(&g)?);
        escalate(
            || {
                let (x, y): (Schur<GuardedRat>, Schur<GuardedRat>) =
                    (build_rat(&f)?, build_rat(&g)?);
                let r = guarded(|| crate::plethysm::plethysm(&x, &y))?;
                Some(dump_integral(&r, "plethysm"))
            },
            || {
                let (x, y): (Schur<BigRational>, Schur<BigRational>) =
                    (build_rat_wide(&f), build_rat_wide(&g));
                dump_integral(&crate::plethysm::plethysm(&x, &y), "plethysm")
            },
        )
    })
}

// --- classical quantities ---------------------------------------------------

/// Skew a Schur-basis element by `g`, given in `basis` — the adjoint of
/// multiplication by g under the Hall inner product.
///
/// `basis` is a basis name or one-letter code, the same six spellings
/// [`to_power`] lists, and defaults to `"s"`. It selects which rule runs, not
/// merely how `g` is read: `"h"`, `"e"`, and `"p"` take the native Pieri /
/// dual-Pieri / Murnaghan–Nakayama paths and never touch
/// Littlewood–Richardson, while `"s"`, `"m"`, and `"f"` go through it. Passing
/// the same function in a different basis gives the same answer by a different
/// algorithm, which is exactly what the oracle script checks.
///
/// `f` is always Schur-basis, and the result is Schur-basis in the element
/// order.
///
/// ```text
/// >>> symfn.skew_by([([3, 1], 1)], [([1], 1)], "s")
/// [((2, 1), 1), ((3,), 1)]
/// >>> symfn.skew_by([([3, 1], 1)], [([1], 1)], "p")
/// [((2, 1), 1), ((3,), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a partition and
/// `basis` is a known name or code.
#[pyfunction]
#[pyo3(signature = (f, g, basis = "s"))]
fn skew_by(f: Terms, g: Terms, basis: &str) -> PyResult<Terms> {
    interruptible(move || {
        // Resolved once, before either pass runs, so both matches below are
        // exhaustive with no unreachable "unknown basis" arm.
        let b = Basis::parse(basis)?;
        fn fast(f: &Parsed, g: &Parsed, b: Basis) -> Option<Schur<Guarded>> {
            let sf: Schur<Guarded> = build(f)?;
            Some(match b {
                Basis::Schur => SkewBy::skew_by(&sf, &build::<_, Schur<Guarded>>(g)?),
                Basis::Homogeneous => SkewBy::skew_by(&sf, &build::<_, Homogeneous<Guarded>>(g)?),
                Basis::Elementary => SkewBy::skew_by(&sf, &build::<_, Elementary<Guarded>>(g)?),
                Basis::PowerSum => SkewBy::skew_by(&sf, &build::<_, PowerSum<Guarded>>(g)?),
                Basis::Monomial => SkewBy::skew_by(&sf, &build::<_, Monomial<Guarded>>(g)?),
                Basis::Forgotten => SkewBy::skew_by(&sf, &build::<_, Forgotten<Guarded>>(g)?),
            })
        }
        fn wide(f: &Parsed, g: &Parsed, b: Basis) -> Schur<BigInt> {
            let sf: Schur<BigInt> = build_wide(f);
            match b {
                Basis::Schur => SkewBy::skew_by(&sf, &build_wide::<_, Schur<BigInt>>(g)),
                Basis::Homogeneous => {
                    SkewBy::skew_by(&sf, &build_wide::<_, Homogeneous<BigInt>>(g))
                }
                Basis::Elementary => SkewBy::skew_by(&sf, &build_wide::<_, Elementary<BigInt>>(g)),
                Basis::PowerSum => SkewBy::skew_by(&sf, &build_wide::<_, PowerSum<BigInt>>(g)),
                Basis::Monomial => SkewBy::skew_by(&sf, &build_wide::<_, Monomial<BigInt>>(g)),
                Basis::Forgotten => SkewBy::skew_by(&sf, &build_wide::<_, Forgotten<BigInt>>(g)),
            }
        }
        let (f, g) = (terms_arg(&f)?, terms_arg(&g)?);
        Ok(escalate(
            || Some(dump(&guarded(|| fast(&f, &g, b))??)),
            || dump(&wide(&f, &g, b)),
        ))
    })
}

/// Evaluate a Schur-basis element at the alphabet `xs`.
///
/// Integer alphabet only: this is the bridge to concrete values, and a float
/// one would silently make an exact answer approximate.
///
/// The alphabet's length is the number of variables, so a term whose shape has
/// more rows than that contributes `0` — the same vanishing
/// [`principal_specialization`] reports, and an answer rather than a refusal.
///
/// Returns one `int` of any size.
///
/// ```text
/// >>> symfn.evaluate_schur([([2, 1], 1)], [1, 2])
/// 6
/// >>> symfn.evaluate_schur([([2, 1], 1)], [1])
/// 0
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition. A non-integer
/// alphabet entry is a `TypeError` from the extraction itself.
#[pyfunction]
fn evaluate_schur(a: Terms, xs: Vec<Coeff>) -> PyResult<Coeff> {
    interruptible(move || {
        let a = terms_arg(&a)?;
        Ok(escalate(
            || {
                let s: Schur<Guarded> = build(&a)?;
                let alphabet: Option<Vec<Guarded>> =
                    xs.iter().map(<Guarded as Boundary>::from_coeff).collect();
                let alphabet = alphabet?;
                Some(guarded(|| s.eval(&alphabet))?.to_coeff())
            },
            || {
                let s: Schur<BigInt> = build_wide(&a);
                let alphabet: Vec<BigInt> = xs.iter().map(Coeff::to_big).collect();
                Coeff::Big(s.eval(&alphabet))
            },
        ))
    })
}

/// The expansion of an element in `n` variables, as
/// `[(exponent vector, coefficient), ...]` with each vector of length `n`.
///
/// `src` is the basis the input is written in, spelled as [`convert_indexed`]
/// spells it. The result is a polynomial in normal form — no repeated exponent
/// vector, no zero coefficient — so a caller can hand it straight to a
/// polynomial ring's dict constructor without merging.
///
/// This is what Sage's `SymmetricFunction.expand(n)` needs, and it is the
/// displacement of all five `compute_*_with_alphabet`. Note it is *not*
/// [`evaluate_schur`]: there the alphabet is values in a ring, here it is
/// indeterminates, and only the second can answer `expand`.
///
/// Every basis reaches it through `m`, which is the basis whose expansion is
/// definitional; the conversion is where any escalation happens, since laying
/// out the exponents copies coefficients and does no arithmetic.
///
/// The pairs are grouped by the monomial-basis term they come from, and those
/// groups arrive in the element order; **within a group the order of the
/// rearrangements is not part of the contract**. A shape with more than `n`
/// rows contributes nothing.
///
/// ```text
/// >>> symfn.expand_alphabet([([1], 1)], "Schur", 2)
/// [((1, 0), 1), ((0, 1), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition and `src` is a basis
/// name or one-letter code, the same six spellings [`to_power`] lists.
#[pyfunction]
fn expand_alphabet(a: Terms, src: &str, n: usize) -> PyResult<Vec<(Key, Coeff)>> {
    interruptible(move || {
        let a = terms_arg(&a)?;
        let terms = match Basis::parse(src)? {
            Basis::Monomial => return rows_of(&a, n),
            Basis::Schur => s_to_m(&a),
            Basis::Homogeneous => s_to_m(&relay(&h_to_s(&a))),
            Basis::Elementary => s_to_m(&relay(&e_to_s(&a))),
            Basis::PowerSum => s_to_m(&relay(&p_to_s(&a))),
            Basis::Forgotten => s_to_m(&relay(&f_to_s(&a))),
        };
        rows_of(&relay(&terms), n)
    })
}

/// Lay a monomial-basis element out over `n` variables.
///
/// Split out of [`expand_alphabet`] because the `"monomial"` source is already
/// in that basis and must not be routed through a conversion to reach this.
fn rows_of(terms: &Parsed, n: usize) -> PyResult<Vec<(Key, Coeff)>> {
    fn rows<C: Boundary>(terms: &Parsed, n: usize) -> Option<Vec<(Key, Coeff)>> {
        let m: Monomial<C> = build(terms)?;
        Some(expand_rows(&m, n))
    }
    fn rows_wide<C: Wide>(terms: &Parsed, n: usize) -> Vec<(Key, Coeff)> {
        let m: Monomial<C> = build_wide(terms);
        expand_rows(&m, n)
    }
    fn expand_rows<C: Ring + Boundary>(m: &Monomial<C>, n: usize) -> Vec<(Key, Coeff)> {
        m.expand(n)
            .into_iter()
            .map(|(alpha, c)| (alpha.into(), c.to_coeff()))
            .collect()
    }
    Ok(escalate(
        || rows::<Guarded>(terms, n),
        || rows_wide::<BigInt>(terms, n),
    ))
}

/// Multiply two monomial-basis elements.
///
/// The structure constants count pairs of rearrangements summing to a fixed
/// exponent vector, so they are non-negative and no intermediate is wider than
/// the answer — see [`Monomial::mul`](crate::Monomial::mul). Displaces
/// Symmetrica's `mult_monomial_monomial`.
///
/// Both arguments and the result are monomial-basis elements, in the element
/// order.
///
/// ```text
/// >>> symfn.monomial_multiply([([1], 1)], [([1], 1)])
/// [((1, 1), 2), ((2,), 1)]
/// ```
///
/// The `2` is the value that separates this from the Schur product, where
/// `s_1 · s_1` has both coefficients 1.
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a partition.
#[pyfunction]
fn monomial_multiply(a: Terms, b: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let (a, b) = (terms_arg(&a)?, terms_arg(&b)?);
        Ok(escalate(
            || {
                let (x, y): (Monomial<Guarded>, Monomial<Guarded>) = (build(&a)?, build(&b)?);
                Some(dump(&guarded(|| x.mul(&y))?))
            },
            || {
                let (x, y): (Monomial<BigInt>, Monomial<BigInt>) = (build_wide(&a), build_wide(&b));
                dump(&x.mul(&y))
            },
        ))
    })
}

/// The semistandard Young tableaux of shape λ and weight μ, as lists of rows.
///
/// In Symmetrica's `kostka_tab` order, which Sage's `SemistandardTableaux`
/// doctests pin — see
/// [`semistandard_tableaux`](crate::semistandard_tableaux). Use
/// [`kostka_number`] when only the count is wanted: this returns `K_{λμ}`
/// objects and that returns one integer.
///
/// **μ is a composition here**, unlike everywhere else on this surface: entry
/// `i` is how many `i + 1`s the tableau carries, so an interior zero is a value
/// that goes unused rather than a malformed partition. `(2, 0, 1)` and `(2, 1)`
/// have the same count and different tableaux. Sage reaches this constantly —
/// `SemistandardTableaux(λ)` iterates every content vector of `|λ|`, and most
/// are not weakly decreasing — so validating μ as a partition here refuses the
/// majority of the calls the entry point exists to serve.
///
/// Empty is an answer: off-degree there are no such tableaux, which is the
/// same theorem [`kostka_number`] reports as `0`.
///
/// A tableau is a list of rows, each row a list of entries, top row first.
///
/// ```text
/// >>> symfn.semistandard_tableaux([2, 1], [2, 0, 1])
/// [[[1, 1], [3]]]
/// >>> symfn.semistandard_tableaux([2, 1], [2, 1])
/// [[[1, 1], [2]]]
/// ```
///
/// The interior zero is what the two calls separate: `(2, 0, 1)` uses the
/// value 3 and never the value 2, so it is a different set of tableaux from
/// `(2, 1)` and not a malformed spelling of it.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition. μ is a composition and is not
/// held to that.
#[pyfunction]
fn semistandard_tableaux(la: Vec<u32>, mu: Vec<u32>) -> PyResult<Vec<Vec<Vec<u32>>>> {
    interruptible(move || Ok(crate::kostka::semistandard_tableaux(&part_arg(&la)?, &mu)))
}

/// f^λ — the number of standard Young tableaux of shape λ, i.e. the dimension
/// of the irreducible S_{|λ|} representation. `None` past `u128`.
///
/// `None` is a refusal to guess, not a zero: the hook-length product exceeded
/// what the counter holds, and no approximate answer is returned in its place.
///
/// ```text
/// >>> symfn.dimension([2, 1])
/// 2
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn dimension(la: Vec<u32>) -> PyResult<Option<u128>> {
    interruptible(move || Ok(crate::eval::dimension(&part_arg(&la)?)))
}

/// s_λ(1^n), the dimension of the GL_n irreducible. `None` on overflow.
///
/// Zero is an answer: `s_λ` in `n` variables vanishes when `ℓ(λ) > n`, so a λ
/// with too many rows is a legitimate `0` and not a refusal. `None` is the
/// refusal, and means the product exceeded what the counter holds.
///
/// ```text
/// >>> symfn.principal_specialization([2, 1], 3)
/// 8
/// >>> symfn.principal_specialization([2, 1], 1)
/// 0
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn principal_specialization(la: Vec<u32>, n: u32) -> PyResult<Option<u128>> {
    interruptible(move || Ok(crate::eval::principal_specialization(&part_arg(&la)?, n)))
}

/// s_λ(1, q, …, q^{n−1}) as a coefficient list in q, lowest degree first.
///
/// Entry `i` is the coefficient of `q^i`, so the list is dense and includes
/// its zeros — unlike the exponent-keyed rows the parameter families use. The
/// list is empty when the specialization vanishes, which happens exactly when
/// `ℓ(λ) > n`. Setting `q = 1` recovers
/// [`principal_specialization`](principal_specialization).
///
/// ```text
/// >>> symfn.principal_specialization_q([2, 1], 3)
/// [0, 1, 2, 2, 2, 1]
/// ```
///
/// The list starts at `q^0` always; the leading zero here says the lowest
/// exponent occurring in `s_{21}(1, q, q²)` is 1, which is `n(λ)`.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn principal_specialization_q(la: Vec<u32>, n: u32) -> PyResult<Vec<i128>> {
    interruptible(move || Ok(crate::eval::principal_specialization_q(&part_arg(&la)?, n)))
}

/// Kostka number K_{λμ}.
///
/// Zero is an answer: there are no semistandard tableaux of shape λ and weight
/// μ unless `|λ| = |μ|` and λ dominates μ, so both are `0` rather than errors —
/// see [`lr_coefficient`] on which zeros this module refuses instead.
///
/// **μ may be a composition**, and is sorted on the way in. `K_{λμ}` is
/// symmetric in μ — the Bender–Knuth involutions are a bijection between the
/// tableaux of content μ and of any rearrangement of it — so the count is the
/// same and only the sorting is needed. [`semistandard_tableaux`], which
/// returns the tableaux themselves, may **not** do this: they are relabeled by
/// the rearrangement, not preserved.
///
/// Returns one nonnegative `int`.
///
/// ```text
/// >>> symfn.kostka_number([2, 1], [1, 1, 1])
/// 2
/// >>> symfn.kostka_number([2, 1], [1, 2])
/// 1
/// ```
///
/// The second call is the composition case: `(1, 2)` is sorted to `(2, 1)` and
/// counted, rather than refused.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition. μ is sorted, so any list of
/// nonnegative integers is accepted for it.
#[pyfunction]
fn kostka_number(la: Vec<u32>, mu: Vec<u32>) -> PyResult<u128> {
    interruptible(move || {
        let mut mu = mu;
        mu.sort_unstable_by(|a, b| b.cmp(a));
        Ok(crate::kostka::kostka(&part_arg(&la)?, &part_arg(&mu)?))
    })
}

/// Symmetric-group character χ^λ(μ).
///
/// Exact at every size: `try_character` reports overflow rather than wrapping,
/// and the recursion then re-runs in `BigInt`. |χ^λ(μ)| ≤ √(|λ|!), which passes
/// `i128` around |λ| = 58 — reachable, so this is not hypothetical.
///
/// Returns one `int`, which may be negative and has no ceiling.
///
/// ```text
/// >>> symfn.character_value([2, 1], [1, 1, 1])
/// 2
/// >>> symfn.character_value([2, 1], [3])
/// -1
/// ```
///
/// Those two values are the `S_3` standard representation on the identity and
/// on a 3-cycle, which pin λ as the shape and μ as the class rather than the
/// other way round.
///
/// # Raises
///
/// Raises `ValueError` if `|λ| ≠ |μ|`: `μ` must index a conjugacy class of
/// `S_{|λ|}`, so off-degree there is no value to return — unlike
/// [`lr_coefficient`], whose off-degree zero is a theorem. Also raises unless
/// both arguments are partitions.
#[pyfunction]
fn character_value(la: Vec<u32>, mu: Vec<u32>) -> PyResult<Coeff> {
    interruptible(move || {
        let (l, m) = (part_arg(&la)?, part_arg(&mu)?);
        same_degree(&[("la", &l), ("mu", &m)])?;
        Ok(match crate::character::try_character(&l, &m) {
            Some(v) => Coeff::Small(v),
            None => Coeff::Big(crate::character::character_in::<BigInt>(&l, &m)),
        })
    })
}

/// The internal (Kronecker) product of two Schur-basis elements.
///
/// Both arguments and the result are Schur-basis, in the element order. The
/// product is degree-preserving rather than degree-adding: `s_λ ∗ s_μ` is
/// zero unless `|λ| = |μ|`, and then it lives in that same degree. The route
/// runs over ℚ and the answer is integral, so a denominator surviving it is
/// reported rather than truncated.
///
/// ```text
/// >>> symfn.internal_product([([2, 1], 1)], [([2, 1], 1)])
/// [((1, 1, 1), 1), ((2, 1), 1), ((3,), 1)]
/// ```
///
/// The degrees are what separate this from [`schur_multiply`], which would
/// answer in degree 6.
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a partition,
/// and if the power-sum route produces a non-integral coefficient.
#[pyfunction]
fn internal_product(a: Terms, b: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let (a, b) = (terms_arg(&a)?, terms_arg(&b)?);
        escalate(
            || {
                let (x, y): (Schur<GuardedRat>, Schur<GuardedRat>) =
                    (build_rat(&a)?, build_rat(&b)?);
                let r = guarded(|| ops::internal(&x, &y))?;
                Some(dump_integral(&r, "Kronecker"))
            },
            || {
                let (x, y): (Schur<BigRational>, Schur<BigRational>) =
                    (build_rat_wide(&a), build_rat_wide(&b));
                dump_integral(&ops::internal(&x, &y), "Kronecker")
            },
        )
    })
}

/// A single Kronecker coefficient g^ν_{λμ}, computed **without forming the
/// product**.
///
/// The same relationship to [`internal_product`] that [`lr_coefficient`] has to
/// [`schur_multiply`], and it is worth stating because the answer is not the
/// one the Rust-side naming suggests. `internal_product` is `s → p`, a diagonal
/// multiply, and `p → s` back; that last step expands every `p_ρ` into every λ
/// ⊢ n, which is the p(n)² work and the memory wall. This route sums
/// `χ^λ(ρ)χ^μ(ρ)χ^ν(ρ)/z_ρ` instead: three character rows, no symmetric
/// function ever built, O(p(n)) memory.
///
/// Measured over `BigRational` — which is what the wheel always carries, so it
/// is the comparison that applies here — the character sum wins at *every*
/// degree, by a margin that widens with it (`examples/bench_kron_coeff.rs`,
/// `docs/record/kronecker.md`). Over a fixed-width ring the product route wins
/// below n ≈ 12, but no caller reaches this function that way.
///
/// So `internal_product` remains the right call when more than a few ν are
/// wanted, since it produces them all at once; this is the right call for one.
/// Sage times out past n = 32 on either.
///
/// Returns one nonnegative `int`.
///
/// ```text
/// >>> symfn.kronecker_coefficient([2, 1], [2, 1], [2, 1])
/// 1
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless λ, μ and ν are partitions sharing a degree —
/// `g^ν_{λμ}` is an `S_n` multiplicity and has no meaning across degrees. The
/// Rust-side
/// [`kronecker_via_characters`](crate::ops::kronecker_via_characters) instead
/// returns `0` there by convention, so that composing it with
/// [`internal_product`] stays total; this boundary is stricter on purpose
/// (`docs/policies/failure.md`, R11).
#[pyfunction]
fn kronecker_coefficient(la: Vec<u32>, mu: Vec<u32>, nu: Vec<u32>) -> PyResult<Coeff> {
    interruptible(move || {
        let (l, m, n) = (part_arg(&la)?, part_arg(&mu)?, part_arg(&nu)?);
        same_degree(&[("la", &l), ("mu", &m), ("nu", &n)])?;
        Ok(escalate(
            || {
                let v: GuardedRat = guarded(|| ops::kronecker_via_characters(&l, &m, &n))?;
                (v.denom() == 1).then(|| Coeff::Small(v.numer()))
            },
            || {
                let v: BigRational = ops::kronecker_via_characters(&l, &m, &n);
                assert!(
                    v.is_integer(),
                    "Kronecker coefficient is not an integer over BigRational, \
                     which is a bug rather than an overflow"
                );
                Coeff::Big(v.to_integer())
            },
        ))
    })
}

/// A conversion whose output partitions are returned as **indices** rather than
/// as lists: `[(degree, index, coefficient), ...]`, where `index` is into
/// [`partitions`] of that degree.
///
/// The caller almost always has to turn each output partition into an object of
/// its own — a Sage `Partition`, say — and doing that per term dominates. A
/// list of parts must be copied, hashed and looked up before it can be mapped
/// to a cached object; an index is a direct array access.
/// `docs/record/python-and-sage-interop.md` owns the measurement of that
/// lookup on Sage's conversion shim.
///
/// The order is `partitions(degree)`, which is exposed for exactly this reason,
/// so a caller can build its own table once per degree and never build another
/// partition object.
///
/// The triples themselves come in the element order of the output partitions,
/// which is not the `partitions(degree)` order — the index is a lookup key,
/// not a position in this list.
///
/// ```text
/// >>> symfn.convert_indexed([([2], 1)], "Schur", "monomial")
/// [(2, 1, 1), (2, 0, 1)]
/// ```
///
/// Both terms have degree 2, and `1` and `0` index `(1, 1)` and `(2,)` in
/// `partitions(2)`.
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition and both bases are
/// spelled as [`convert_terms`] says.
#[pyfunction]
fn convert_indexed(a: Terms, src: &str, dst: &str) -> PyResult<Vec<(u32, usize, Coeff)>> {
    interruptible(move || {
        Ok(routed(&terms_arg(&a)?, src, dst)?
            .into_iter()
            .map(|(p, c)| {
                let n: u32 = p.iter().sum();
                (n, index_of(n, &p), c)
            })
            .collect())
    })
}

/// The same conversion as [`convert_indexed`], with each output partition
/// spelled out rather than given as an index.
///
/// This is the entry point for a caller that does not hold a table of
/// partitions per degree — a rational input, say, whose coefficients have to be
/// rebuilt term by term anyway, so the index would save nothing.
///
/// `src` and `dst` are each a basis name or one-letter code: `"Schur"` or
/// `"s"`, `"homogeneous"` or `"h"`, `"elementary"` or `"e"`, `"powersum"` or
/// `"p"`, `"monomial"` or `"m"`, `"forgotten"` or `"f"`. The pair is what
/// selects the route: h, e and p reach each other directly, and everything
/// else composes through Schur. Naming the pair in one call is the point —
/// composing two calls in the caller's own language forces the hub and is what
/// made `p → h` cost p(n) determinants (`docs/record/transitions.md`).
///
/// Every pair lands in ℤ; the conversions that divide are [`to_power`]'s,
/// which is why `dst` may not be the power-sum basis. Result in the element
/// order.
///
/// ```text
/// >>> symfn.convert_terms([([2], 1)], "powersum", "homogeneous")
/// [((1, 1), -1), ((2,), 2)]
/// ```
///
/// That is `p_2 = 2·h_2 − h_11`, the direct rule rather than the composition
/// through Schur.
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition and both names are
/// known.
#[pyfunction]
fn convert_terms(a: Terms, src: &str, dst: &str) -> PyResult<Terms> {
    interruptible(move || routed(&terms_arg(&a)?, src, dst))
}

/// `src → dst` over already-validated terms: the direct rule when the pair has
/// one, otherwise out through Schur and back.
fn routed(a: &Parsed, src: &str, dst: &str) -> PyResult<Terms> {
    let src = Basis::parse(src)?;
    let dst = Basis::parse(dst).map_err(|_| bad_dst_basis(dst))?;
    // Rejected before any conversion runs, not after the source has been
    // carried to Schur; the arm below that repeats it is what keeps the match
    // exhaustive without a wildcard.
    if dst == Basis::PowerSum {
        return Err(bad_dst_basis("powersum"));
    }
    match (src, dst) {
        (Basis::PowerSum, Basis::Homogeneous) => return Ok(p_to_h(a)),
        (Basis::PowerSum, Basis::Elementary) => return Ok(p_to_e(a)),
        (Basis::Homogeneous, Basis::Elementary) => return Ok(h_to_e(a)),
        (Basis::Elementary, Basis::Homogeneous) => return Ok(e_to_h(a)),
        _ => {}
    }
    // Validated, so this is the caller's partition in normal form — the shape
    // `index_of` can look up. Handing the raw list through instead is how the
    // identity conversion used to panic on `[2, 1, 0]`: a partition this
    // boundary accepts everywhere else, but not a key in the table.
    let terms = match src {
        Basis::Schur => dump_parsed(a),
        Basis::Monomial => m_to_s(a),
        Basis::Homogeneous => h_to_s(a),
        Basis::Elementary => e_to_s(a),
        Basis::PowerSum => p_to_s(a),
        Basis::Forgotten => f_to_s(a),
    };
    Ok(match dst {
        Basis::Schur => terms,
        Basis::Monomial => s_to_m(&relay(&terms)),
        Basis::Homogeneous => s_to_h(&relay(&terms)),
        Basis::Elementary => s_to_e(&relay(&terms)),
        Basis::Forgotten => s_to_f(&relay(&terms)),
        Basis::PowerSum => return Err(bad_dst_basis("powersum")),
    })
}

/// Validated terms back out as `Terms`, for the conversion that is the
/// identity on the basis but not on the representation.
fn dump_parsed(a: &Parsed) -> Terms {
    a.iter()
        .map(|(p, c)| (p.parts().to_vec().into(), (*c).clone()))
        .collect()
}

/// One of the six classical bases, as named at this boundary.
///
/// Every entry point that takes a basis argument — [`convert_terms`],
/// [`convert_indexed`], [`to_power`], [`expand_alphabet`], [`skew_by`] —
/// resolves it through [`Basis::parse`], so the accepted spellings and the
/// "unknown basis" error exist once. Each accepts the full name or its
/// one-letter code:
///
/// | full name        | code  |
/// |------------------|-------|
/// | `"Schur"`        | `"s"` |
/// | `"homogeneous"`  | `"h"` |
/// | `"elementary"`   | `"e"` |
/// | `"powersum"`     | `"p"` |
/// | `"monomial"`     | `"m"` |
/// | `"forgotten"`    | `"f"` |
///
/// An enum rather than a `&str` threaded through so every downstream match is
/// exhaustive.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Basis {
    Schur,
    Homogeneous,
    Elementary,
    PowerSum,
    Monomial,
    Forgotten,
}

impl Basis {
    fn parse(name: &str) -> PyResult<Self> {
        Ok(match name {
            "Schur" | "s" => Basis::Schur,
            "homogeneous" | "h" => Basis::Homogeneous,
            "elementary" | "e" => Basis::Elementary,
            "powersum" | "p" => Basis::PowerSum,
            "monomial" | "m" => Basis::Monomial,
            "forgotten" | "f" => Basis::Forgotten,
            other => return Err(bad_basis(other)),
        })
    }
}

fn bad_basis(other: &str) -> PyErr {
    PyValueError::new_err(format!(
        "unknown basis {other:?}; expected one of Schur, homogeneous, elementary, powersum, monomial, forgotten \
         or the one-letter codes s, h, e, p, m, f"
    ))
}

/// The `dst` half of [`convert_terms`] and [`convert_indexed`], which reaches
/// five bases rather than six.
///
/// A separate message because the general one listed `powersum` among the
/// accepted names and then rejected it: every pair here lands in ℤ, and the
/// conversions that divide belong to [`to_power`]. Naming the alternative is
/// the difference between a caller fixing the call and a caller concluding the
/// library is wrong about its own basis list.
fn bad_dst_basis(other: &str) -> PyErr {
    PyValueError::new_err(format!(
        "unknown target basis {other:?}; expected one of Schur, homogeneous, elementary, monomial, forgotten \
         or the codes s, h, e, m, f. \
         The power-sum basis is not a target here because that conversion is rational; use to_power(a, src)"
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
///
/// The order is reverse lexicographic — the one-row partition first, the
/// all-ones partition last — and it is **not** the increasing lexicographic
/// order elements come back in. `partitions(0)` is one empty tuple, not an
/// empty list.
///
/// ```text
/// >>> symfn.partitions(4)
/// [(4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]
/// ```
///
/// Raises nothing.
#[pyfunction]
fn partitions(n: u32) -> PyResult<Vec<Key>> {
    interruptible(move || {
        Ok({
            crate::memo::partitions_cached(n)
                .iter()
                .map(|p| p.parts().to_vec().into())
                .collect()
        })
    })
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
///
/// Both indices run over [`partitions`]`(n)`: `i` is the shape, `j` the
/// class.
///
/// ```text
/// >>> symfn.character_table(3)
/// [[1, 1, 1], [-1, 0, 2], [1, -1, 1]]
/// ```
///
/// Row 0 is the trivial character. Row 1 is the standard one, whose value on
/// the identity class is the `2` at the **end** of the row — the orientation
/// that a class-major reading would put at the start.
///
/// Raises nothing.
#[pyfunction]
fn character_table(n: u32) -> PyResult<Vec<Vec<i128>>> {
    interruptible(move || Ok(crate::character::character_table(n)))
}

/// The full Kostka table of degree `n`: `table[i][j]` = K_{λⁱ λʲ}.
///
/// `u128`-backed, with the same caveat as [`character_table`].
///
/// Both indices run over [`partitions`]`(n)`, so the table is upper
/// triangular in that order: `K_{λμ}` vanishes unless λ dominates μ.
///
/// ```text
/// >>> symfn.kostka_table(3)
/// [[1, 1, 1], [0, 1, 2], [0, 0, 1]]
/// ```
///
/// The `2` is `K_{(2,1),(1,1,1)}`, which fixes the orientation: shape first,
/// weight second.
///
/// Raises nothing.
#[pyfunction]
fn kostka_table(n: u32) -> PyResult<Vec<Vec<u128>>> {
    interruptible(move || Ok(crate::kostka::kostka_table(n)))
}

// --- operations -------------------------------------------------------------

/// The ω involution on a Schur-basis element.
///
/// Conjugates indices and copies coefficients, so it cannot overflow; it runs
/// once, over `BigInt`. Argument and result are Schur-basis, in the element
/// order, and ω is its own inverse.
///
/// ```text
/// >>> symfn.omega([([3], 1)])
/// [((1, 1, 1), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition.
#[pyfunction]
fn omega(a: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let s: Schur<BigInt> = build_wide(&terms_arg(&a)?);
        Ok(dump(&s.omega()))
    })
}

/// The Hall inner product of two Schur-basis elements.
///
/// Returns one `int`. The Schur basis is orthonormal for it, so this is the
/// sum of the products of matching coefficients and zero on shapes of
/// different degree.
///
/// ```text
/// >>> symfn.hall_inner_product([([2, 1], 1)], [([2, 1], 1)])
/// 1
/// >>> symfn.hall_inner_product([([2, 1], 1)], [([3], 1)])
/// 0
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term of both arguments is a partition.
#[pyfunction]
fn hall_inner_product(a: Terms, b: Terms) -> PyResult<Coeff> {
    interruptible(move || {
        let (a, b) = (terms_arg(&a)?, terms_arg(&b)?);
        Ok(escalate(
            || {
                let (x, y): (Schur<Guarded>, Schur<Guarded>) = (build(&a)?, build(&b)?);
                Some(guarded(|| ops::hall::<Guarded, _, _>(&x, &y))?.to_coeff())
            },
            || {
                let (x, y): (Schur<BigInt>, Schur<BigInt>) = (build_wide(&a), build_wide(&b));
                Coeff::Big(ops::hall::<BigInt, _, _>(&x, &y))
            },
        ))
    })
}

// --- Hopf structure ---------------------------------------------------------

/// The skew Schur function s_{λ/μ}.
///
/// Zero is an answer: `s_{λ/μ} = 0` unless μ ⊆ λ, by the standard convention
/// that the skew diagram is empty otherwise, and it comes back as the empty
/// list. Result is Schur-basis in the element order, with nonnegative
/// coefficients.
///
/// ```text
/// >>> symfn.skew_schur([2, 1], [1])
/// [((1, 1), 1), ((2,), 1)]
/// >>> symfn.skew_schur([2, 1], [2, 2])
/// []
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless both arguments are partitions.
#[pyfunction]
fn skew_schur(la: Vec<u32>, mu: Vec<u32>) -> PyResult<Terms> {
    interruptible(move || {
        let s: Schur<BigInt> = hopf::skew_schur(&part_arg(&la)?, &part_arg(&mu)?);
        Ok(dump(&s))
    })
}

/// The coproduct Δ, as `[((mu, nu), coefficient), ...]`.
///
/// The support is an ordered **pair** of partitions, and the pairs come back
/// in increasing lexicographic order of that pair. `Δ(s_λ) = Σ s_μ ⊗ s_ν`
/// with the Littlewood–Richardson coefficient `c^λ_{μν}`, so both ends of the
/// range appear: the empty partition pairs with λ itself.
///
/// ```text
/// >>> symfn.coproduct([([2], 1)])
/// [(((), (2,)), 1), (((1,), (1,)), 1), (((2,), ()), 1)]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn coproduct(a: Terms) -> PyResult<Vec<((Key, Key), Coeff)>> {
    interruptible(move || {
        fn split<C: Boundary>(x: &Schur<C>) -> Vec<((Key, Key), Coeff)> {
            hopf::coproduct(x)
                .terms()
                .iter()
                .map(|((m, n), c)| {
                    (
                        (m.parts().to_vec().into(), n.parts().to_vec().into()),
                        c.to_coeff(),
                    )
                })
                .collect()
        }
        let a = terms_arg(&a)?;
        Ok(escalate(
            || {
                let s: Schur<Guarded> = build(&a)?;
                guarded(|| split(&s))
            },
            || {
                let s: Schur<BigInt> = build_wide(&a);
                split(&s)
            },
        ))
    })
}

/// The antipode S on a Schur-basis element.
///
/// Conjugates and negates, so like [`omega`] it cannot overflow. Result
/// Schur-basis, in the element order.
///
/// ```text
/// >>> symfn.antipode([([2], 1)])
/// [((1, 1), 1)]
/// >>> symfn.antipode([([1], 1)])
/// [((1,), -1)]
/// ```
///
/// The sign is `(−1)^{|λ|}`, which the two degrees above separate.
///
/// # Raises
///
/// Raises `ValueError` unless every term is a partition.
#[pyfunction]
fn antipode(a: Terms) -> PyResult<Terms> {
    interruptible(move || {
        let s: Schur<BigInt> = build_wide(&terms_arg(&a)?);
        Ok(dump(&hopf::antipode(&s)))
    })
}

// --- Hall–Littlewood --------------------------------------------------------

/// `Q'_λ(x; t) = Σ_μ K_{μλ}(t) s_μ`, as `[(mu, [(t_exponent, coefficient),
/// ...])]`.
///
/// The coefficients are polynomials, so this cannot reuse [`Terms`]. Sparse in
/// the exponent, which is how [`QtPoly`](crate::QtPoly) already holds them.
///
/// ⚠️ This is `Q'`, the **transformed** Hall–Littlewood function, not `P` and
/// not `Q` — its Schur coefficients are the Kostka–Foulkes polynomials, and at
/// `t = 0` it is `s_λ`. [`hall_littlewood_p`] is the other normalization.
/// Sage's equivalent is `Sym.hall_littlewood().Qp()`.
///
/// Rows come in the element order of μ; each row's `(exponent, coefficient)`
/// pairs are sparse, in increasing exponent, with no zero coefficients.
///
/// ```text
/// >>> symfn.hall_littlewood([1, 1])
/// [((1, 1), [(0, 1)]), ((2,), [(1, 1)])]
/// ```
///
/// So `Q'_{11} = s_11 + t·s_2`, which is the value that separates `Q'` from
/// `P`: [`hall_littlewood_p`] of the same shape is `s_11` alone.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn hall_littlewood(la: Vec<u32>) -> PyResult<Vec<(Key, Vec<(u32, Coeff)>)>> {
    interruptible(move || {
        let l = part_arg(&la)?;
        Ok(escalate(
            || Some(hl_rows(&guarded(|| crate::hall_littlewood::<Guarded>(&l))?)),
            || hl_rows(&crate::hall_littlewood::<BigInt>(&l)),
        ))
    })
}

/// Every `Q'_λ` for `λ ⊢ n`, sharing the recursion's suffixes across the
/// degree.
///
/// Rows come in `partitions(n)` order, not the element order, and each row is
/// one [`hall_littlewood`] answer for that λ.
///
/// ```text
/// >>> symfn.hall_littlewood_table(2)[1]
/// ((1, 1), [((1, 1), [(0, 1)]), ((2,), [(1, 1)])])
/// ```
///
/// Index 1 is `(1, 1)` because [`partitions`] puts the one-row shape first.
///
/// Raises nothing.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn hall_littlewood_table(n: u32) -> PyResult<Vec<(Key, Vec<(Key, Vec<(u32, Coeff)>)>)>> {
    interruptible(move || {
        Ok({
            crate::hall_littlewood_table::<i128>(n)
                .into_iter()
                .map(|(lambda, hl)| (lambda.parts().to_vec().into(), hl_rows(&hl)))
                .collect()
        })
    })
}

/// `K_{λμ}(t)` as `[(t_exponent, coefficient), ...]`.
///
/// Zero unless `|λ| = |μ|` and λ ⊵ μ. The zero polynomial crosses as an empty
/// list, so off-degree arguments return `[]` rather than raising.
///
/// Sparse and in increasing exponent, with no zero coefficients. λ is the
/// **shape** and μ the weight, the orientation [`kostka_table`] uses, and
/// `t = 1` recovers the Kostka number.
///
/// ```text
/// >>> symfn.kostka_foulkes([2, 1], [1, 1, 1])
/// [(1, 1), (2, 1)]
/// >>> symfn.kostka_foulkes([2], [3])
/// []
/// ```
///
/// `K_{(2,1),(1^3)}(t) = t + t²`. Its value at `t = 1` is 2, which is
/// `kostka_number([2, 1], [1, 1, 1])`.
///
/// # Raises
///
/// Raises `ValueError` unless both arguments are partitions.
#[pyfunction]
fn kostka_foulkes(la: Vec<u32>, mu: Vec<u32>) -> PyResult<Vec<(u32, Coeff)>> {
    interruptible(move || {
        Ok(t_poly(&crate::kostka_foulkes::<i128>(
            &part_arg(&la)?,
            &part_arg(&mu)?,
        )))
    })
}

/// Every `K_{λμ}(t)` for a fixed μ, as `[(lambda, [(t_exponent,
/// coefficient)])]`.
///
/// One `Q'_μ` *is* the column, so this costs what a single value costs — see
/// [`crate::kf`]. Rows come in the element order of λ, and the zero entries
/// are omitted rather than listed.
///
/// ```text
/// >>> symfn.kostka_foulkes_column([1, 1])
/// [((1, 1), [(0, 1)]), ((2,), [(1, 1)])]
/// ```
///
/// μ is the **weight** and the λ are the shapes, so this is a column of
/// [`kostka_foulkes_table`] and not a row.
///
/// # Raises
///
/// Raises `ValueError` unless μ is a partition.
#[pyfunction]
fn kostka_foulkes_column(mu: Vec<u32>) -> PyResult<Vec<(Key, Vec<(u32, Coeff)>)>> {
    interruptible(move || {
        Ok(crate::kostka_foulkes_column::<i128>(&part_arg(&mu)?)
            .into_iter()
            .map(|(lambda, k)| (lambda.parts().to_vec().into(), t_poly(&k)))
            .collect())
    })
}

/// `P_λ(x; t)` in the Schur basis — the other Hall–Littlewood normalization.
///
/// Costs the whole degree: the inversion needs every dominance-smaller `P`, so
/// use [`hall_littlewood_p_table`] when more than one shape is wanted.
///
/// Rows in the element order of μ; Sage's equivalent is
/// `Sym.hall_littlewood().P()`. Coefficients may be negative here, which they
/// never are for [`hall_littlewood`].
///
/// ```text
/// >>> symfn.hall_littlewood_p([2])
/// [((1, 1), [(1, -1)]), ((2,), [(0, 1)])]
/// ```
///
/// So `P_2 = s_2 − t·s_11`. At `t = 0` this is `s_2`, as `Q'` also is — the
/// two normalizations separate only at higher order in `t`.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn hall_littlewood_p(la: Vec<u32>) -> PyResult<Vec<(Key, Vec<(u32, Coeff)>)>> {
    interruptible(move || Ok(hl_rows(&crate::hall_littlewood_p::<i128>(&part_arg(&la)?))))
}

/// Every `P_λ` for `λ ⊢ n`, from one inversion of the Kostka–Foulkes matrix.
///
/// Rows in `partitions(n)` order, each one a [`hall_littlewood_p`] answer.
///
/// ```text
/// >>> symfn.hall_littlewood_p_table(2)[1]
/// ((1, 1), [((1, 1), [(0, 1)])])
/// ```
///
/// Raises nothing.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn hall_littlewood_p_table(n: u32) -> PyResult<Vec<(Key, Vec<(Key, Vec<(u32, Coeff)>)>)>> {
    interruptible(move || {
        Ok({
            crate::hall_littlewood_p_table::<i128>(n)
                .into_iter()
                .map(|(lambda, hl)| (lambda.parts().to_vec().into(), hl_rows(&hl)))
                .collect()
        })
    })
}

/// The whole `K_{λμ}(t)` matrix for degree `n`, indexed as `partitions(n)` is.
///
/// Same orientation as [`kostka_table`], of which this is the t-analogue:
/// `table[i][j]` is `K_{λⁱλʲ}(t)`, and `t = 1` recovers that table entry for
/// entry. Asking for the p(n)² values one at a time would recompute each column
/// p(n) times.
///
/// A zero entry is the empty list, so the matrix is dense in its indices and
/// sparse in each entry's exponents.
///
/// ```text
/// >>> symfn.kostka_foulkes_table(2)
/// [[[(0, 1)], [(1, 1)]], [[], [(0, 1)]]]
/// ```
///
/// The lone `[]` is `K_{(1,1),(2)}(t) = 0`, below the diagonal in dominance —
/// the orientation check, since the transpose would put it at `[0][1]`.
///
/// Raises nothing.
#[pyfunction]
fn kostka_foulkes_table(n: u32) -> PyResult<Vec<Vec<Vec<(u32, Coeff)>>>> {
    interruptible(move || {
        Ok({
            crate::kostka_foulkes_table::<i128>(n)
                .into_iter()
                .map(|row| row.iter().map(t_poly).collect())
                .collect()
        })
    })
}

/// A `Schur<QtPoly>` as `[(mu, [(t_exponent, coefficient), ...])]`.
fn hl_rows<C: Ring + ToCoeff>(hl: &Schur<crate::QtPoly<C>>) -> Vec<(Key, Vec<(u32, Coeff)>)> {
    hl.terms()
        .iter()
        .map(|(mu, c)| (mu.parts().to_vec().into(), t_poly(c)))
        .collect()
}

/// A `QtPoly` known not to involve q, as `[(t_exponent, coefficient)]`.
fn t_poly<C: Ring + ToCoeff>(p: &crate::QtPoly<C>) -> Vec<(u32, Coeff)> {
    p.terms()
        .map(|((a, b), v)| {
            debug_assert_eq!(*a, 0, "Hall-Littlewood must not involve q");
            (*b, v.to_coeff())
        })
        .collect()
}

// --- Macdonald --------------------------------------------------------------

/// One Macdonald expansion: per μ, the numerator's `(q_exp, t_exp, coeff)`
/// terms and the denominator's `(q_exp, t_exp, multiplicity)` **factors**.
///
/// The denominator is handed over factored rather than expanded, which is both
/// cheaper and what a caller wants: `prod((1 - q^a*t^b)^m)` builds the element
/// directly in a fraction field, where expanding here and re-factoring there
/// would be work done twice. See [`Frac`](crate::Frac) for why the factored
/// form is the representation and not an optimization.
type MacTerms = Vec<(Key, Vec<(u32, u32, Coeff)>, Vec<(u32, u32, u32)>)>;

fn mac_terms<C: Ring + ToCoeff>(f: &Monomial<crate::Frac<C>>) -> MacTerms {
    f.terms()
        .iter()
        .map(|(mu, c)| {
            let (num, den) = c.parts();
            (
                mu.parts().to_vec().into(),
                num.terms()
                    .map(|((a, b), v)| (*a, *b, v.to_coeff()))
                    .collect(),
                den.map(|(&(a, b), &m)| (a, b, m)).collect(),
            )
        })
        .collect()
}

/// Macdonald `P_λ(x; q, t)` in the monomial basis.
///
/// Escalates: the fixed-width pass reports rather than wrapping, and the call
/// re-runs over `BigInt`, so there is no wall here. There is one underneath —
/// at the extremal one-row shape `λ = (n)`, `i128` gives out at n = 30
/// (`docs/record/failure-and-overflow.md`).
///
/// Each triple is `(mu, numerator terms, denominator factors)`: a numerator
/// term is `(q exponent, t exponent, coefficient)`, a denominator factor is
/// `(q exponent, t exponent, multiplicity)` standing for `(1 − q^a t^b)^m`. An
/// empty factor list means the coefficient is a polynomial. Rows in the
/// element order of μ. Sage's equivalent is `Sym.macdonald().P()`.
///
/// ```text
/// >>> symfn.macdonald_p([1])
/// [((1,), [(0, 0, 1)], [])]
/// ```
///
/// `P_(1) = m_1` with coefficient 1, which is the normalization: `P` is monic
/// in the monomial basis, where `Q` and `J` are not.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn macdonald_p(la: Vec<u32>) -> PyResult<MacTerms> {
    interruptible(move || {
        let l = part_arg(&la)?;
        Ok(escalate(
            || Some(mac_terms(&guarded(|| crate::macdonald_p::<Guarded>(&l))?)),
            || mac_terms(&crate::macdonald_p::<BigInt>(&l)),
        ))
    })
}

/// Macdonald `Q_λ = b_λ · P_λ`. Escalates, as [`macdonald_p`] does; the `i128`
/// wall underneath is n = 26 at λ = (n).
///
/// Same encoding as [`macdonald_p`]; Sage's equivalent is
/// `Sym.macdonald().Q()`.
///
/// ```text
/// >>> symfn.macdonald_q([1])
/// [((1,), [(0, 0, 1), (0, 1, -1)], [(1, 0, 1)])]
/// ```
///
/// So `Q_(1) = (1 − t)/(1 − q) · m_1`, where [`macdonald_p`] gives `m_1` — the
/// value that separates the two normalizations at the smallest shape.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn macdonald_q(la: Vec<u32>) -> PyResult<MacTerms> {
    interruptible(move || {
        let l = part_arg(&la)?;
        Ok(escalate(
            || Some(mac_terms(&guarded(|| crate::macdonald_q::<Guarded>(&l))?)),
            || mac_terms(&crate::macdonald_q::<BigInt>(&l)),
        ))
    })
}

/// Macdonald `J_λ = c_λ · P_λ`, the integral form — every coefficient is a
/// polynomial, so the denominator list comes back empty. Escalates, as
/// [`macdonald_p`] does; the `i128` wall underneath is n = 26 at λ = (n).
///
/// Same encoding as [`macdonald_p`]; Sage's equivalent is
/// `Sym.macdonald().J()`.
///
/// ```text
/// >>> symfn.macdonald_j([1, 1])
/// [((1, 1), [(0, 0, 1), (0, 1, -1), (0, 2, -1), (0, 3, 1)], [])]
/// ```
///
/// The empty third slot on every term is the integral form's signature: `J`
/// clears the denominators `Q` carries.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn macdonald_j(la: Vec<u32>) -> PyResult<MacTerms> {
    interruptible(move || {
        let l = part_arg(&la)?;
        Ok(escalate(
            || Some(mac_terms(&guarded(|| crate::macdonald_j::<Guarded>(&l))?)),
            || mac_terms(&crate::macdonald_j::<BigInt>(&l)),
        ))
    })
}

/// The Schur functions of degree `n` in the Macdonald `J` basis, as
/// `[(lambda, [(mu, numerator, denominator), ...]), ...]` — the **inverse** of
/// the `J → s` transition, in the cell encoding [`MacTerms`] already carries.
///
/// The whole degree, because that is the unit of work: the projection route in
/// [`schur_in_j_table`](crate::schur_in_j_table) reads every row off one
/// `(q,t)`-Kostka table, which costs the same whether one row is wanted or all
/// of them.
///
/// Escalates, as [`macdonald_p`] does. This route divides by `z_ν` where `J`
/// itself does not, so the fixed-width pass here carries denominators; the
/// answers do not, and a numerator that survives with one is a bug rather than
/// a representable result.
///
/// Rows in `partitions(n)` order, one per λ; inside a row the μ come in the
/// element order and the zero entries are omitted. The encoding of each entry
/// is [`macdonald_p`]'s.
///
/// ```text
/// >>> symfn.schur_in_macdonald_j(1)
/// [((1,), [((1,), [(0, 0, 1)], [(0, 1, 1)])])]
/// ```
///
/// So `s_(1) = J_(1)/(1 − t)`, which is the inverse direction: the transition
/// out of `J` would have no denominator here.
///
/// # Raises
///
/// Raises `ValueError` if a numerator coefficient is not an integer.
/// `H_μ[X(1−q)]` is integral in the Schur basis — the `(q,t)`-Kostka entries
/// are, and `s_λ[X(1−q)]` is a `ℤ[q]`-combination of Schur functions — so a
/// fraction reaching here means the power-sum round trip did not cancel.
#[pyfunction]
fn schur_in_macdonald_j(n: u32) -> PyResult<Vec<(Key, MacTerms)>> {
    interruptible(move || {
        fn rows<C: BoundaryRat>(
            table: &[Vec<crate::Frac<C>>],
            n: u32,
        ) -> PyResult<Vec<(Key, MacTerms)>> {
            let parts = crate::partitions_of(n);
            let mut out = Vec::with_capacity(parts.len());
            for (i, lambda) in parts.iter().enumerate() {
                let mut row = MacTerms::new();
                for (j, mu) in parts.iter().enumerate() {
                    if table[i][j].is_zero() {
                        continue;
                    }
                    let (num, den) = table[i][j].parts();
                    let mut terms = Vec::with_capacity(num.len());
                    for (&(a, b), v) in num.terms() {
                        let (numer, denom) = v.split();
                        if !denom.to_big().is_one() {
                            let (x, y) = (numer.to_big(), denom.to_big());
                            return Err(PyValueError::new_err(format!(
                                "non-integral s_{lambda} in J_{mu} coefficient {x}/{y}"
                            )));
                        }
                        terms.push((a, b, numer));
                    }
                    row.push((
                        mu.parts().to_vec().into(),
                        terms,
                        den.map(|(&(a, b), &m)| (a, b, m)).collect(),
                    ));
                }
                out.push((lambda.parts().to_vec().into(), row));
            }
            Ok(out)
        }

        escalate(
            || {
                guarded(|| crate::schur_in_j_table::<GuardedRat>(n))
                    .map(|table| rows::<GuardedRat>(&table, n))
            },
            || rows::<BigRational>(&crate::schur_in_j_table::<BigRational>(n), n),
        )
    })
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
/// the object is. The denominator is handed over **factored**, for the same
/// reason [`MacTerms`] hands over its binomials factored: a caller rebuilding
/// this in ℚ(α) wants `prod(u*a + v)`, and expanding here to re-factor there is
/// work done twice.
///
/// The atoms are *primitive* (`gcd(u, v) = 1`), so the factorization is
/// canonical — unlike the (q,t) family, where `1 − q²` is reducible. See
/// [`AFrac`](crate::afrac::AFrac).
type JackCell = (Vec<Coeff>, Vec<(u32, u32, u32)>, u128);

/// One Jack expansion: per basis index μ, a [`JackCell`].
type JackTerms = Vec<(Key, Vec<Coeff>, Vec<(u32, u32, u32)>, u128)>;

fn jack_cell<C: Boundary>(c: &crate::AFrac<C>) -> JackCell {
    let (num, den, scale) = c.parts();
    (
        num.iter().map(ToCoeff::to_coeff).collect(),
        den.map(|(&(u, v), &m)| (u, v, m)).collect(),
        scale,
    )
}

fn jack_terms<C: Boundary>(f: &Monomial<crate::AFrac<C>>) -> JackTerms {
    f.terms()
        .iter()
        .map(|(mu, c)| {
            let (n, d, s) = jack_cell(c);
            (mu.parts().to_vec().into(), n, d, s)
        })
        .collect()
}

fn jack_terms_p<C: Boundary>(f: &PowerSum<crate::AFrac<C>>) -> JackTerms {
    f.terms()
        .iter()
        .map(|(mu, c)| {
            let (n, d, s) = jack_cell(c);
            (mu.parts().to_vec().into(), n, d, s)
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

/// Jack `P_λ(x; α)` in the monomial basis: monic in `m_λ`,
/// dominance-triangular.
///
/// Computed by the Laplace–Beltrami eigenoperator recursion, which enumerates
/// no tableaux at all. Sage has no whole-degree entry point and walls at
/// n = 12; see [`jack_table`].
///
/// Each row is `(mu, numerator, denominator atoms, scale)`, the [`JackCell`]
/// encoding: the numerator dense in the α-exponent, the denominator factored
/// as `(u, v, multiplicity)` for `(u·α + v)^m`. Rows in the element order of
/// μ. Sage's equivalent is `Sym.jack().P()`.
///
/// ```text
/// >>> symfn.jack_p([2])
/// [((1, 1), [2], [(1, 1, 1)], 1), ((2,), [1], [], 1)]
/// ```
///
/// So `P_(2) = 2/(α + 1)·m_11 + m_2`: monic in `m_λ`, which is what separates
/// `P` from [`jack_q`] and [`jack_j`], both of which scale that leading 1.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn jack_p(la: Vec<u32>) -> PyResult<JackTerms> {
    interruptible(move || {
        let l = part_arg(&la)?;
        Ok(jack_escalate_m(|| crate::jack_p(&l), || crate::jack_p(&l)))
    })
}

/// Jack `Q_λ = (H_λ/H'_λ)·P_λ`, the basis dual to `P` under `⟨·,·⟩_α`.
///
/// Same [`JackCell`] encoding as [`jack_p`], rows in the element order.
///
/// ```text
/// >>> symfn.jack_q([1, 1])
/// [((1, 1), [2], [(1, 0, 1), (1, 1, 1)], 1)]
/// ```
///
/// `Q_{11} = 2/(α(α + 1))·m_11`, where [`jack_p`] of the same shape is `m_11`
/// — the value that separates the two normalizations.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn jack_q(la: Vec<u32>) -> PyResult<JackTerms> {
    interruptible(move || {
        let l = part_arg(&la)?;
        Ok(jack_escalate_m(|| crate::jack_q(&l), || crate::jack_q(&l)))
    })
}

/// Jack `J_λ = H_λ·P_λ`, the integral form.
///
/// Every coefficient is a polynomial in α with non-negative integer
/// coefficients, divisible by `u_μ = ∏ m_i(μ)!` (\[KS\] Thm 1.1) — so the
/// denominator list comes back empty and `scale` comes back 1. None of that is
/// arranged: the coefficients arrive through fraction arithmetic and cancel.
///
/// Same [`JackCell`] encoding as [`jack_p`], rows in the element order.
///
/// ```text
/// >>> symfn.jack_j([2])
/// [((1, 1), [2], [], 1), ((2,), [1, 1], [], 1)]
/// ```
///
/// `J_(2) = 2·m_11 + (1 + α)·m_2`. Every third slot is empty and every fourth
/// is 1, which is the integral form's signature.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn jack_j(la: Vec<u32>) -> PyResult<JackTerms> {
    interruptible(move || {
        let l = part_arg(&la)?;
        Ok(jack_escalate_m(|| crate::jack_j(&l), || crate::jack_j(&l)))
    })
}

/// Every `P_λ` of degree `n` — the unit of work Sage has no entry point for,
/// and the one `docs/record/jack.md` measures the walls in.
///
/// Rows in `partitions(n)` order, each one a [`jack_p`] answer.
///
/// ```text
/// >>> symfn.jack_table(2)[1]
/// ((1, 1), [((1, 1), [1], [], 1)])
/// ```
///
/// Raises nothing.
#[pyfunction]
fn jack_table(n: u32) -> PyResult<Vec<(Key, JackTerms)>> {
    interruptible(move || {
        Ok({
            crate::partitions_of(n)
                .into_iter()
                .map(|l| {
                    let rows = jack_escalate_m(|| crate::jack_p(&l), || crate::jack_p(&l));
                    (l.parts().to_vec().into(), rows)
                })
                .collect()
        })
    })
}

/// `J_λ` in the **power-sum** basis — the Jack character table, and the unit
/// the Goulden–Jackson pipeline consumes.
///
/// Same [`JackCell`] encoding as [`jack_p`], but the support indexes `p_μ`
/// rather than `m_μ`. Rows in the element order, and here the fourth slot is
/// a genuine scale rather than 1.
///
/// ```text
/// >>> symfn.jack_j_powersum([2])
/// [((1, 1), [2], [], 2), ((2,), [0, 2], [], 2)]
/// ```
///
/// So `J_(2) = p_11 + α·p_2`, after dividing each row by its scale of 2.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn jack_j_powersum(la: Vec<u32>) -> PyResult<JackTerms> {
    interruptible(move || {
        let l = part_arg(&la)?;
        Ok(escalate(
            || guarded(|| jack_terms_p(&crate::jack_j_powersum::<Guarded>(&l))),
            || jack_terms_p(&crate::jack_j_powersum::<BigInt>(&l)),
        ))
    })
}

/// `⟨J_λ, J_λ⟩_α = H_λ·H'_λ`, returned **factored** as `[(u, v, mult)]`.
///
/// A product of `2|λ|` linear forms and no pairing at all, where Sage prices
/// the same table like a full expansion (`docs/record/jack.md`).
///
/// Each `(u, v, m)` stands for `(u·α + v)^m`, the atoms primitive and in
/// increasing `(u, v)` order. The answer is the whole product, with no
/// numerator and no scale — the value is a polynomial in α, never a fraction.
///
/// ```text
/// >>> symfn.jack_norm_j([1])
/// [(0, 1, 1), (1, 0, 1)]
/// ```
///
/// So `⟨J_(1), J_(1)⟩_α = 1 · α`, the two factors `H_λ` and `H'_λ` contribute.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn jack_norm_j(la: Vec<u32>) -> PyResult<Vec<(u32, u32, u32)>> {
    interruptible(move || {
        Ok(crate::jack_norm_j(&part_arg(&la)?)
            .into_iter()
            .map(|((u, v), m)| (u, v, m as u32))
            .collect())
    })
}

/// `⟨J_λ J_μ, J_ν⟩_α` — **Stanley's object**, whose membership in `ℕ[α]` is his
/// 1989 conjecture and still open.
///
/// A negative coefficient is a result to report, not a bug: nothing here
/// asserts positivity. `J[3,2,1]²` is out of Sage's range
/// (`docs/record/jack.md`).
///
/// Zero is likewise an answer: the pairing is graded, so `|λ| + |μ| ≠ |ν|`
/// vanishes by orthogonality rather than being a malformed question, and the
/// zero polynomial comes back as an empty numerator.
///
/// Returns one [`JackCell`] — `(numerator, denominator atoms, scale)` — not an
/// element.
///
/// ```text
/// >>> symfn.jack_structure_constant([1], [1], [2])
/// ([0, 0, 2], [], 1)
/// >>> symfn.jack_structure_constant([1], [1], [3])
/// ([], [], 1)
/// ```
///
/// The first is `2α²`, dense in the α-exponent, so the two leading zeros are
/// the absent `α⁰` and `α¹` terms rather than padding.
///
/// # Raises
///
/// Raises `ValueError` unless all three arguments are partitions.
#[pyfunction]
fn jack_structure_constant(la: Vec<u32>, mu: Vec<u32>, nu: Vec<u32>) -> PyResult<JackCell> {
    interruptible(move || {
        let (a, b, c) = (part_arg(&la)?, part_arg(&mu)?, part_arg(&nu)?);
        Ok(jack_escalate(
            || crate::jack_structure_constant(&a, &b, &c),
            || crate::jack_structure_constant(&a, &b, &c),
        ))
    })
}

/// Stanley's **whole table**: every `⟨J_λ J_μ, J_ν⟩_α` with `|λ| = |μ| = k`,
/// as `(lambda, mu, nu, numerator, denominator atoms, scalar)`.
///
/// Zero entries are omitted. Prefer this over looping
/// [`jack_structure_constant`], which recomputes the same p-expansions on
/// every call. It runs to k = 8, degree 16 and 111 804 triples, past anything
/// Sage reaches for even one entry (`docs/record/jack.md`).
///
/// Positivity is Stanley's 1989 conjecture and is **open**. This returns the
/// values and asserts nothing about them.
///
/// The triples run λ, then μ, then ν, each in `partitions` order — λ and μ
/// over `partitions(k)`, ν over `partitions(2k)`. Each row's last three slots
/// are the [`JackCell`] encoding.
///
/// ```text
/// >>> symfn.stanley_table(1)
/// [([1], [1], [2], [0, 0, 2], [], 1), ([1], [1], [1, 1], [0, 0, 2], [], 1)]
/// ```
///
/// Raises nothing. `k = 0` is the single empty triple, not an empty table.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn stanley_table(
    k: u32,
) -> PyResult<
    Vec<(
        Vec<u32>,
        Vec<u32>,
        Vec<u32>,
        Vec<Coeff>,
        Vec<(u32, u32, u32)>,
        u128,
    )>,
> {
    interruptible(move || {
        Ok({
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
        })
    })
}

/// `⟨f, g⟩_α` for two monomial-basis elements whose coefficients are
/// **integer** polynomials in α, given densely: `[(partition, [c0, c1, …])]`.
///
/// The restriction to integral coefficients is the honest boundary: a general
/// `AFrac` input would need the atoms marshalled in too, and every element a
/// caller actually pairs — `J_λ`, integer combinations of them — is already of
/// this shape. Clear denominators on the Python side first if yours is not.
///
/// Returns one [`JackCell`]. The coefficient lists are dense in the
/// α-exponent, so `[1]` is the constant 1 and `[0, 1]` is α.
///
/// ```text
/// >>> symfn.jack_scalar([([1], [1])], [([1], [1])])
/// ([0, 1], [], 1)
/// ```
///
/// `⟨m_1, m_1⟩_α = α`, which is the α-deformed pairing rather than the Hall
/// one, where it would be 1.
///
/// # Raises
///
/// Raises `ValueError` unless every support in both arguments is a partition.
#[pyfunction]
fn jack_scalar(f: Vec<(Vec<u32>, Vec<i128>)>, g: Vec<(Vec<u32>, Vec<i128>)>) -> PyResult<JackCell> {
    interruptible(move || {
        fn shapes(rows: &[(Vec<u32>, Vec<i128>)]) -> PyResult<Vec<Partition>> {
            rows.iter().map(|(mu, _)| part_arg(mu)).collect()
        }
        fn build<C: Ring>(
            shapes: &[Partition],
            rows: &[(Vec<u32>, Vec<i128>)],
        ) -> Monomial<crate::AFrac<C>> {
            let mut out = Monomial::zero();
            for (mu, (_, coeffs)) in shapes.iter().zip(rows) {
                let num: Vec<C> = coeffs.iter().map(|&v| C::from_i128(v)).collect();
                out.add_term(mu.clone(), crate::AFrac::from_coeffs(num));
            }
            out
        }
        let (sf, sg) = (shapes(&f)?, shapes(&g)?);
        Ok(jack_escalate(
            || crate::jack_scalar(&build::<Guarded>(&sf, &f), &build::<Guarded>(&sg, &g)),
            || crate::jack_scalar(&build::<BigInt>(&sf, &f), &build::<BigInt>(&sg, &g)),
        ))
    })
}

/// The zonal polynomial, as exact `(numerator, denominator)` pairs.
///
/// `integral_form = True` returns `J^{(2)}`; `False` returns `P^{(2)}`.
///
/// ⚠️ Sage's `zonal()` is `P^{(2)}` and \[GJ\]'s `Z_λ` is `J^{(2)}`; the two
/// differ by `H_λ(2)`. Measured, not assumed. Both are reachable rather than
/// one under an ambiguous name, because a caller that picks the wrong one
/// still gets plausible-looking output.
///
/// Rows are `(mu, numerator, denominator)` over the monomial basis, in the
/// element order, each fraction in lowest terms. There is no α here: the
/// zonal case is `α = 2`, already substituted.
///
/// ```text
/// >>> symfn.zonal([2], True)
/// [((1, 1), 2, 1), ((2,), 3, 1)]
/// >>> symfn.zonal([2], False)
/// [((1, 1), 2, 3), ((2,), 1, 1)]
/// ```
///
/// `J^{(2)}_(2) = 2·m_11 + 3·m_2` against `P^{(2)}_(2) = 2/3·m_11 + m_2`.
/// Those two values are what a caller who picked the wrong flag would see, so
/// they are the pin.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition.
#[pyfunction]
fn zonal(la: Vec<u32>, integral_form: bool) -> PyResult<Vec<(Key, Coeff, Coeff)>> {
    interruptible(move || {
        let l = part_arg(&la)?;
        let f = if integral_form {
            crate::zonal_j(&l)
        } else {
            crate::zonal_p(&l)
        };
        Ok(f.terms()
            .iter()
            .map(|(mu, c)| {
                (
                    mu.parts().to_vec().into(),
                    Coeff::Small(c.numer()),
                    Coeff::Small(c.denom()),
                )
            })
            .collect())
    })
}

/// The Goulden–Jackson connection tables `c^λ_{μν}(b)` and `h^λ_{μν}(b)` at
/// degree `n`, as `(lambda, mu, nu, [b-coefficients], denominator)`.
///
/// Returns `(c, h)`. Two open conjectures live here — Matchings-Jack on `c`,
/// the b-conjecture on `h` — and no package computes either table
/// (`docs/research-gaps.md`). `ℚ[b]`-polynomiality and `c`'s integrality are
/// theorems and are enforced (a failure raises); **positivity is the open
/// question and is only observed**, so a negative coefficient comes back as
/// data rather than an exception.
///
/// Both tables are lists of `(lambda, mu, nu, numerator, denominator)`, the
/// numerator dense in the `b`-exponent and the denominator one positive
/// integer. Triples run λ, μ, ν in increasing lexicographic order of the
/// triple, and the zero entries are omitted.
///
/// ```text
/// >>> symfn.gj_connection_tables(2)[1][0]
/// ([1, 1], [2], [2], [1], 1)
/// ```
///
/// # Raises
///
/// Raises `ValueError` if `ℚ[b]`-polynomiality or `c`'s integrality fails at
/// this degree — both are theorems, so either would be a defect in this
/// library rather than a caller error.
#[pyfunction]
#[allow(clippy::type_complexity)]
fn gj_connection_tables(
    n: u32,
) -> PyResult<(
    Vec<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<Coeff>, u128)>,
    Vec<(Vec<u32>, Vec<u32>, Vec<u32>, Vec<Coeff>, u128)>,
)> {
    interruptible(move || {
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
    })
}

/// `a^λ_{μν}`, the class-algebra connection coefficient of `S_n`, from
/// characters alone.
///
/// The independent object the `b = 0` slice of [`gj_connection_tables`] is
/// pinned against — no Jack polynomial and no fraction field anywhere in it.
///
/// Returns one nonnegative `int`, and zero is an answer: the product of two
/// class sums need not meet a given class.
///
/// ```text
/// >>> symfn.class_algebra_coefficient([1, 1, 1], [2, 1], [2, 1])
/// 3
/// >>> symfn.class_algebra_coefficient([2, 1], [2, 1], [2, 1])
/// 0
/// ```
///
/// λ is the **target** class and μ, ν the two being multiplied: the three
/// transpositions of `S_3` multiply to the identity in 3 ways and to a
/// transposition in none.
///
/// # Raises
///
/// Raises `ValueError` unless all three are partitions of one `n`: these
/// index conjugacy classes of the same symmetric group, so a mismatch is a
/// malformed question.
#[pyfunction]
fn class_algebra_coefficient(la: Vec<u32>, mu: Vec<u32>, nu: Vec<u32>) -> PyResult<Coeff> {
    interruptible(move || {
        let (l, m, n) = (part_arg(&la)?, part_arg(&mu)?, part_arg(&nu)?);
        same_degree(&[("la", &l), ("mu", &m), ("nu", &n)])?;
        Ok(Coeff::Small(crate::class_algebra_coefficient(&l, &m, &n)))
    })
}

// --- (q,t)-Kostka -----------------------------------------------------------

/// A `QtPoly` over ℤ, as `[(q_exp, t_exp, coeff)]`.
///
/// No integrality assertion, and none is needed: the Bergeron–Haiman recursion
/// never divides by an integer. This path runs over `i128`, so a non-integral
/// value is not representable rather than merely unexpected. `Rat::into_poly`
/// refuses a surviving denominator, and `divide_exact` refuses an inexact
/// division, which is where Macdonald's theorem is enforced.
///
/// `i128` is a ceiling far above this family: `K̃_{λμ}` has non-negative
/// coefficients summing to `f^λ`, and `Σ_λ (f^λ)² = n!`, so nothing here
/// exceeds `√(n!)`. The wall that ceiling produces sits around degree 57.
fn qt_poly<C: Ring + ToCoeff>(p: &crate::QtPoly<C>) -> Vec<(u32, u32, Coeff)> {
    p.terms().map(|(&(a, b), v)| (a, b, v.to_coeff())).collect()
}

/// The (q,t)-Kostka polynomial `K_{λμ}(q,t)`, from `J_μ = Σ_λ K_{λμ} S_λ(x;t)`.
///
/// Computes the whole of `J_μ`; use [`qt_kostka_column`] for more than one λ at
/// a fixed μ, and [`qt_kostka_table`] for a whole degree.
///
/// ⚠️ This is `K_{λμ}(q,t)`, **not** the modified `K̃_{λμ}(q,t)` the modern
/// literature writes; [`macdonald_ht`] carries that one. Terms are
/// `(q exponent, t exponent, coefficient)`, sparse, in increasing exponent
/// pair, with the zero polynomial an empty list.
///
/// ```text
/// >>> symfn.qt_kostka([2], [1, 1])
/// [(0, 1, 1)]
/// ```
///
/// `K_{(2),(11)} = t`, where the modified form has `K̃_{(2),(11)} = 1` —
/// `t^{n(μ)} K_{λμ}(q, 1/t)` with `n(μ) = 1`. That is the smallest pair that
/// separates the two conventions.
///
/// # Raises
///
/// Raises `ValueError` if `|λ| ≠ |μ|`: `K_{λμ}(q,t)` is an entry of one
/// degree's matrix, and off-degree there is no entry rather than a zero one.
/// Also raises unless both arguments are partitions.
#[pyfunction]
fn qt_kostka(la: Vec<u32>, mu: Vec<u32>) -> PyResult<Vec<(u32, u32, Coeff)>> {
    interruptible(move || {
        let (l, m) = (part_arg(&la)?, part_arg(&mu)?);
        same_degree(&[("la", &l), ("mu", &m)])?;
        Ok(qt_poly(&crate::qt_kostka::<i128>(&l, &m)))
    })
}

/// Every `K_{λμ}(q,t)` for a fixed μ — one `J_μ`, which is what a single
/// [`qt_kostka`] costs anyway.
///
/// Unlike a Kostka–Foulkes column this one is **dense**: every λ of the degree
/// appears, since `K_{λμ}` is generally nonzero without λ dominating μ. Rows
/// come in `partitions(n)` order, not the element order.
///
/// ```text
/// >>> symfn.qt_kostka_column([1, 1])
/// [((2,), [(0, 1, 1)]), ((1, 1), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless μ is a partition.
#[pyfunction]
fn qt_kostka_column(mu: Vec<u32>) -> PyResult<Vec<(Key, Vec<(u32, u32, Coeff)>)>> {
    interruptible(move || {
        Ok(crate::qt_kostka_column::<i128>(&part_arg(&mu)?)
            .iter()
            .map(|(lambda, k)| (lambda.parts().to_vec().into(), qt_poly(k)))
            .collect())
    })
}

/// The modified Macdonald polynomial `H̃_μ(x;q,t)` in the **Schur** basis, as
/// `[(lambda, [(q_exp, t_exp, coeff), ...])]`.
///
/// Its coefficients are the modified (q,t)-Kostka polynomials
/// `K̃_{λμ}(q,t) = t^{n(μ)} K_{λμ}(q, 1/t)` — the form the modern literature
/// uses, and where Haiman's positivity reads "non-negative integers" with no
/// normalizing power in the way. Rows in the element order of λ.
///
/// ```text
/// >>> symfn.macdonald_ht([2])
/// [((1, 1), [(1, 0, 1)]), ((2,), [(0, 0, 1)])]
/// >>> symfn.macdonald_ht([1, 1])
/// [((1, 1), [(0, 1, 1)]), ((2,), [(0, 0, 1)])]
/// ```
///
/// `H̃_{(2)} = s_2 + q·s_{11}` and `H̃_{(11)} = s_2 + t·s_{11}`. Conjugating μ
/// swaps `q` and `t`, which pins μ as the index and λ as the Schur shape.
/// [`qt_kostka`] of the same degree gives `K_{(2),(11)} = t` where `K̃` here
/// gives 1, so the two are distinguishable on any of these four values.
///
/// # Raises
///
/// Raises `ValueError` unless μ is a partition.
#[pyfunction]
fn macdonald_ht(mu: Vec<u32>) -> PyResult<Vec<(Key, Vec<(u32, u32, Coeff)>)>> {
    interruptible(move || {
        Ok(crate::macdonald_ht::<i128>(&part_arg(&mu)?)
            .terms()
            .iter()
            .map(|(lambda, k)| (lambda.parts().to_vec().into(), qt_poly(k)))
            .collect())
    })
}

/// The whole `K_{λμ}(q,t)` matrix for degree `n`, indexed as `partitions(n)` is
/// — the same orientation as [`kostka_table`] and [`kostka_foulkes_table`], of
/// which this is the two-variable analogue. `q = 0` recovers the latter.
///
/// `table[i][j]` is `K_{λⁱλʲ}(q,t)`, dense in its indices, each entry sparse
/// in its exponents. Unlike [`kostka_foulkes_table`] it is not triangular.
///
/// ```text
/// >>> symfn.qt_kostka_table(2)
/// [[[(0, 0, 1)], [(0, 1, 1)]], [[(1, 0, 1)], [(0, 0, 1)]]]
/// ```
///
/// The `q` below the diagonal is the entry a Kostka–Foulkes table has as 0,
/// which is the check that this is the two-variable family and not a
/// relabeling of that one.
///
/// Raises nothing.
#[pyfunction]
fn qt_kostka_table(n: u32) -> PyResult<Vec<Vec<Vec<(u32, u32, Coeff)>>>> {
    interruptible(move || {
        Ok({
            crate::qt_kostka_table::<i128>(n)
                .into_iter()
                .map(|row| row.iter().map(qt_poly).collect())
                .collect()
        })
    })
}

// --- the Macdonald operator algebra -----------------------------------------

/// A Schur element with `(q,t)`-polynomial coefficients, as `[(lambda, [(q_exp,
/// t_exp, coeff), ...]), ...]` — the same shape [`macdonald_ht`] already
/// returns, so an `H̃` row can be fed straight back in.
type QtSchur = Vec<(Key, Vec<(u32, u32, Coeff)>)>;

fn qt_schur_out<C: Ring + ToCoeff>(f: &Schur<crate::QtPoly<C>>) -> QtSchur {
    f.terms()
        .iter()
        .map(|(lambda, c)| (lambda.parts().to_vec().into(), qt_poly(c)))
        .collect()
}

/// The same, from the ℚ-bounded operators.
///
/// Every operator here maps `ℤ[q,t]`-Schur combinations to `ℤ[q,t]` ones, so a
/// surviving denominator is a bug and is raised rather than rounded. The
/// general path runs over `Rational` only because `s → p` divides by `z_ρ`; the
/// answer is integral by the time it reaches this boundary.
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
        out.push((lambda.parts().to_vec().into(), row));
    }
    Ok(out)
}

/// Read a Schur element from Python, for the operators that require it to be
/// **homogeneous**.
///
/// Every Macdonald operator acts on one degree at a time — its eigenvalue is a
/// function of the shape's size — so a mixed-degree argument has no meaning
/// rather than an awkward one, and `degree_of` asserts on it. The term-list
/// shape gives a caller no hint of that, and every *other* term-list entry
/// point in this module accepts mixed degrees happily, so the mistake is easy
/// and the answer must be an exception rather than a panic.
fn qt_schur_in(rows: &QtSchur) -> PyResult<Schur<crate::QtPoly<crate::Rational>>> {
    let mut out = Schur::zero();
    let mut degree: Option<(u32, Partition)> = None;
    for (lambda, terms) in rows {
        let lambda = part_arg(lambda)?;
        let n = lambda.size();
        match &degree {
            Some((d, first)) if *d != n => {
                return Err(PyValueError::new_err(format!(
                    "a Macdonald operator requires a homogeneous argument: \
                     {first} has degree {d} but {lambda} has degree {n}"
                )))
            }
            Some(_) => {}
            None => degree = Some((n, lambda.clone())),
        }
        let mut c = crate::QtPoly::zero();
        for (a, b, v) in terms {
            let v = v.as_i128().ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err("coefficient does not fit in i128")
            })?;
            c.add_term(*a, *b, crate::Rational::from_int(v));
        }
        out.add_term(lambda, c);
    }
    Ok(out)
}

/// `∇e_n` in the Schur basis — the shuffle theorem's object.
///
/// The closed form, which never divides by an integer and so runs over ℤ.
/// Prefer this to `nabla(elementary)`: it skips the change of basis entirely.
///
/// Rows are `(lambda, [(q exponent, t exponent, coefficient), ...])` in the
/// element order, each row sparse in its exponent pairs. `n = 0` is the empty
/// partition with coefficient 1.
///
/// ```text
/// >>> symfn.nabla_e(2)
/// [((1, 1), [(0, 1, 1), (1, 0, 1)]), ((2,), [(0, 0, 1)])]
/// ```
///
/// So `∇e_2 = (q + t)·s_11 + s_2`. The `q + t` is symmetric, which is what
/// separates `∇` from an operator carrying the `q`/`t` asymmetry.
///
/// Raises nothing.
#[pyfunction]
fn nabla_e(n: u32) -> PyResult<QtSchur> {
    interruptible(move || Ok(qt_schur_out(&crate::nabla_e::<i128>(n))))
}

/// `Δ'_{e_k} e_n` in the Schur basis — the Delta conjecture's object.
///
/// Same encoding as [`nabla_e`], which this recovers at `k = n − 1`.
///
/// ```text
/// >>> symfn.delta_prime_e(1, 2)
/// [((1, 1), [(0, 1, 1), (1, 0, 1)]), ((2,), [(0, 0, 1)])]
/// ```
///
/// Raises nothing.
#[pyfunction]
fn delta_prime_e(k: u32, n: u32) -> PyResult<QtSchur> {
    interruptible(move || Ok(qt_schur_out(&crate::delta_prime_e::<i128>(k, n))))
}

/// `∇F` for an arbitrary homogeneous `F`, given in the Schur basis.
///
/// The argument uses the same encoding as [`nabla_e`]'s result, so an `H̃` row
/// or a `∇e_n` answer can be fed straight back in. `F` must be homogeneous —
/// every Macdonald operator acts one degree at a time.
///
/// ```text
/// >>> symfn.nabla([([2], [(0, 0, 1)])])
/// [((1, 1), [(1, 1, -1)])]
/// ```
///
/// `∇s_2 = −qt·s_11`. The negative coefficient is real: `∇` is not
/// Schur-positive on a general argument, only on `e_n`.
///
/// # Raises
///
/// Raises `ValueError` if the terms are not all of one degree, if a support
/// is not a partition, or if the answer is not integral.
#[pyfunction]
fn nabla(f: QtSchur) -> PyResult<QtSchur> {
    interruptible(move || qt_schur_out_rat(&crate::nabla(&qt_schur_in(&f)?), "nabla"))
}

/// `∇^r F`, sharing one change of basis across the powers — the object
/// Qiu–Zhang's 2026 theorem is about.
///
/// Same encoding and the same homogeneity requirement as [`nabla`]. `r = 0`
/// is the identity and `r = 1` agrees with [`nabla`].
///
/// ```text
/// >>> symfn.nabla_power([([2], [(0, 0, 1)])], 2)
/// [((1, 1), [(1, 2, -1), (2, 1, -1)]), ((2,), [(1, 1, -1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` on the same three conditions [`nabla`] does.
#[pyfunction]
fn nabla_power(f: QtSchur, r: u32) -> PyResult<QtSchur> {
    interruptible(move || {
        qt_schur_out_rat(&crate::nabla_power(&qt_schur_in(&f)?, r), "nabla_power")
    })
}

/// `Δ_{e_k} F`, with eigenvalue `e_k[B_μ]`.
///
/// Same encoding and homogeneity requirement as [`nabla`], and the degree is
/// preserved.
///
/// ```text
/// >>> symfn.delta_ek(1, [([1, 1], [(0, 0, 1)])])
/// [((1, 1), [(0, 0, 1), (0, 1, 1), (1, 0, 1)]), ((2,), [(0, 0, 1)])]
/// ```
///
/// The constant term is the `1` that [`delta_prime_ek`]'s `B_μ − 1`
/// eigenvalue removes, which is what separates the two operators.
///
/// # Raises
///
/// Raises `ValueError` on the same three conditions [`nabla`] does.
#[pyfunction]
fn delta_ek(k: u32, f: QtSchur) -> PyResult<QtSchur> {
    interruptible(move || {
        let ek = crate::deltaop::elementary(k);
        qt_schur_out_rat(&crate::delta(&ek, &qt_schur_in(&f)?), "delta")
    })
}

/// `Δ'_{e_k} F`, with eigenvalue `e_k[B_μ − 1]`.
///
/// Same encoding and homogeneity requirement as [`nabla`].
///
/// ```text
/// >>> symfn.delta_prime_ek(1, [([1, 1], [(0, 0, 1)])])
/// [((1, 1), [(0, 1, 1), (1, 0, 1)]), ((2,), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` on the same three conditions [`nabla`] does.
#[pyfunction]
fn delta_prime_ek(k: u32, f: QtSchur) -> PyResult<QtSchur> {
    interruptible(move || {
        let ek = crate::deltaop::elementary(k);
        qt_schur_out_rat(&crate::delta_prime(&ek, &qt_schur_in(&f)?), "delta_prime")
    })
}

/// `Θ_{e_k} F`, which raises the degree by `k`.
///
/// Note the cost: Θ expands at degree `n + k`, so it pays for the larger degree
/// and not the input's.
///
/// Same encoding and homogeneity requirement as [`nabla`]; unlike the Δ
/// family the answer sits in degree `n + k`.
///
/// ```text
/// >>> symfn.theta_ek(1, [([1], [(0, 0, 1)])])
/// [((1, 1), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` on the same three conditions [`nabla`] does.
#[pyfunction]
fn theta_ek(k: u32, f: QtSchur) -> PyResult<QtSchur> {
    interruptible(move || {
        let ek = crate::deltaop::elementary(k);
        qt_schur_out_rat(&crate::theta(&ek, &qt_schur_in(&f)?), "theta")
    })
}

/// `ΠF`, with eigenvalue `Π_μ`.
///
/// `Π⁻¹` is deliberately absent: it is genuinely not a polynomial, so it cannot
/// cross this boundary. Only the composite `Θ` can.
///
/// Same encoding and homogeneity requirement as [`nabla`].
///
/// ```text
/// >>> symfn.big_pi([([1], [(0, 0, 1)])])
/// [((1,), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` on the same three conditions [`nabla`] does.
#[pyfunction]
fn big_pi(f: QtSchur) -> PyResult<QtSchur> {
    interruptible(move || qt_schur_out_rat(&crate::big_pi(&qt_schur_in(&f)?), "big_pi"))
}

/// The combinatorial side of the Delta conjecture, in the **monomial** basis,
/// for every `k` at once — entry `k` of the returned list.
///
/// `side` is `"rise"` (a theorem) or `"valley"` (open). One enumeration serves
/// the whole ladder, so asking for one `k` would cost the same.
///
/// The two sides do not cost the same. `"rise"` factors through the per-path
/// LLT polynomials ([`crate::llt`], and `dyck.rs`'s module docs for why), and
/// runs to n = 9. `"valley"` keeps the `(n+1)^{n−1}`-ish labeled enumeration,
/// because `Val` reads the labels: ⚠️ orders of magnitude more, and one degree
/// further is another such step (`docs/record/dyck-paths.md`).
///
/// The list has `n` entries, index `k` holding the side for that `k`, each
/// one a monomial-basis element in the [`nabla_e`] row encoding.
///
/// ```text
/// >>> symfn.delta_conjecture_side(2, "rise")[0]
/// [((1, 1), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless `side` is `"rise"` or `"valley"`.
#[pyfunction]
fn delta_conjecture_side(n: u32, side: &str) -> PyResult<Vec<QtSchur>> {
    interruptible(move || {
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
                    .map(|(mu, c)| (mu.parts().to_vec().into(), qt_poly(c)))
                    .collect()
            })
            .collect())
    })
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
type QtMon = Vec<(Key, Vec<(u32, u32, Coeff)>)>;

/// `i128` for the same reason the (q,t)-Kostka family above uses it: every
/// coefficient here counts tableaux, so it is a non-negative integer bounded by
/// `n!` — 8.7e10 at n = 14, where `i128` holds 1.7e38. `bench_llt` runs the
/// whole ladder at `i64` *and* `i128` and asserts they agree term for term, so
/// the narrower width is checked rather than assumed.
fn qt_mon_out<C: Ring + ToCoeff>(f: &Monomial<crate::QtPoly<C>>) -> QtMon {
    f.terms()
        .iter()
        .map(|(mu, c)| (mu.parts().to_vec().into(), qt_poly(c)))
        .collect()
}

/// A tuple of straight shapes with content offsets, as Python passes it.
fn skew_tuple(shapes: &[Vec<u32>], offsets: Option<Vec<i32>>) -> PyResult<crate::llt::SkewTuple> {
    let ps: Vec<Partition> = parts_arg(shapes)?;
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
    let cells: u32 = ps.iter().map(|p| p.size()).sum();
    cells_arg(cells as usize, "this tuple")?;
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

/// `G̃^(k)_λ(x;q)`, the **cospin** ribbon generating function of \[LLT\] (26),
/// in the monomial basis.
///
/// Empty when λ has no k-ribbon tableaux (nonempty k-core). Sage's
/// `llt(k).cospin(Partition(λ))` is the same object; `docs/record/llt.md` has
/// the comparison.
///
/// Rows are `(weight, [(q exponent, t exponent, coefficient), ...])` in the
/// element order of the weight. The `t` slot is always 0: this family lives
/// in `q` alone.
///
/// ```text
/// >>> symfn.llt_gtilde([2], 2)
/// [((1,), [(0, 0, 1)])]
/// >>> symfn.llt_gtilde([1], 2)
/// []
/// ```
///
/// The empty answer is the nonempty 2-core of `(1)`, a theorem rather than a
/// refusal.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition, `k ≥ 1`, and the shape fits
/// the abacus representation.
#[pyfunction]
#[pyo3(signature = (la, k))]
fn llt_gtilde(la: Vec<u32>, k: u32) -> PyResult<QtMon> {
    interruptible(move || {
        let (l, k) = (part_arg(&la)?, level_arg(k)?);
        abacus_arg(&l, k)?;
        Ok(escalate(
            || {
                Some(qt_mon_out(&guarded(|| {
                    crate::llt::llt_gtilde::<Guarded>(&l, k)
                })?))
            },
            || qt_mon_out(&crate::llt::llt_gtilde::<BigInt>(&l, k)),
        ))
    })
}

/// `H^(k)_μ(x;q) = Σ_R q^{s(R)} x^{w(R)}`, the **spin** family of \[LLT\] (28).
///
/// Takes a partition and a level, never a tuple, and that is a mathematical
/// constraint rather than an API choice: the k-quotient of a shape does not
/// determine `s*`, so there is no honest `H` of a bare tuple. Sage's
/// `llt(k).hspin()[μ]`.
///
/// Same row encoding as [`llt_gtilde`]. Spin and cospin differ by
/// `q^{s*} → q^{−s}`, so the two disagree on any shape with a positive spin
/// range even though they agree below.
///
/// ```text
/// >>> symfn.llt_h([1], 2)
/// [((1,), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless μ is a partition, `k ≥ 1`, and the shape fits
/// the abacus representation.
#[pyfunction]
#[pyo3(signature = (mu, k))]
fn llt_h(mu: Vec<u32>, k: u32) -> PyResult<QtMon> {
    interruptible(move || {
        let (m, k) = (part_arg(&mu)?, level_arg(k)?);
        abacus_arg(&m, k)?;
        Ok(escalate(
            || {
                Some(qt_mon_out(&guarded(|| {
                    crate::llt::llt_h::<Guarded>(&m, k)
                })?))
            },
            || qt_mon_out(&crate::llt::llt_h::<BigInt>(&m, k)),
        ))
    })
}

/// `H̃^(k)_μ = G̃^(k)_{kμ}` (\[LLT\] (27)) — Sage's `llt(k).hcospin()[μ]`.
///
/// Same row encoding as [`llt_gtilde`]. The argument is μ and the shape
/// evaluated is `kμ`, which is the step a caller passing `kμ` directly to
/// [`llt_gtilde`] would duplicate.
///
/// ```text
/// >>> symfn.llt_h_tilde([1], 2)
/// [((1,), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless μ is a partition, `k ≥ 1`, and the shape fits
/// the abacus representation.
#[pyfunction]
#[pyo3(signature = (mu, k))]
fn llt_h_tilde(mu: Vec<u32>, k: u32) -> PyResult<QtMon> {
    interruptible(move || {
        let (m, k) = (part_arg(&mu)?, level_arg(k)?);
        abacus_arg(&m, k)?;
        Ok(escalate(
            || {
                Some(qt_mon_out(&guarded(|| {
                    crate::llt::llt_h_tilde::<Guarded>(&m, k)
                })?))
            },
            || qt_mon_out(&crate::llt::llt_h_tilde::<BigInt>(&m, k)),
        ))
    })
}

/// `Σ_R q^{2s(R)} x^{w(R)}`, the spin-generating grading of \[LT\] (43).
///
/// The rawest of the four normalizations, and the one [`llt_kl_column`] is
/// pinned against. Sage has no entry point for this grading.
///
/// Same row encoding as [`llt_gtilde`]. The exponent is `2s(R)`, so every
/// exponent here is even where the other three normalizations' need not be.
///
/// ```text
/// >>> symfn.llt_g_lt([2], 2)
/// [((1,), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition, `k ≥ 1`, and the shape fits
/// the abacus representation.
#[pyfunction]
#[pyo3(signature = (la, k))]
fn llt_g_lt(la: Vec<u32>, k: u32) -> PyResult<QtMon> {
    interruptible(move || {
        let (l, k) = (part_arg(&la)?, level_arg(k)?);
        abacus_arg(&l, k)?;
        Ok(escalate(
            || {
                Some(qt_mon_out(&guarded(|| {
                    crate::llt::llt_g_lt::<Guarded>(&l, k)
                })?))
            },
            || qt_mon_out(&crate::llt::llt_g_lt::<BigInt>(&l, k)),
        ))
    })
}

/// `H^(k)_μ` for **every** μ ⊢ n — the whole degree, which is the unit
/// `docs/record/llt.md` measures the walls in.
///
/// This is the entry point Sage lacks: there it is `p(n)` separate per-element
/// conversions.
///
/// Rows are `(mu, H)` in the element order of μ, each `H` a [`llt_h`] answer.
///
/// ```text
/// >>> symfn.llt_h_table(1, 2)
/// [((1,), [((1,), [(0, 0, 1)])])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless `k ≥ 1`.
#[pyfunction]
#[pyo3(signature = (n, k))]
fn llt_h_table(n: u32, k: u32) -> PyResult<Vec<(Key, QtMon)>> {
    interruptible(move || {
        Ok(crate::llt::llt_h_table::<i128>(n, level_arg(k)?)
            .iter()
            .map(|(mu, f)| (mu.parts().to_vec().into(), qt_mon_out(f)))
            .collect())
    })
}

/// `G̃^(k)_λ` for **every** λ ⊢ k·n with empty k-core, from a single walk.
///
/// Rows are `(lambda, G)` in the element order of λ, each `G` a
/// [`llt_gtilde`] answer. The shapes with a nonempty k-core are absent rather
/// than present with an empty value, so the list is shorter than
/// `partitions(k·n)`.
///
/// ```text
/// >>> symfn.llt_gtilde_table(1, 2)
/// [((1, 1), [((1,), [(0, 0, 1)])]), ((2,), [((1,), [(0, 0, 1)])])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless `k ≥ 1` and the degree fits the abacus
/// representation.
#[pyfunction]
#[pyo3(signature = (n, k))]
fn llt_gtilde_table(n: u32, k: u32) -> PyResult<Vec<(Key, QtMon)>> {
    interruptible(move || {
        let k = level_arg(k)?;
        abacus_table_arg(n, k)?;
        Ok(crate::llt::llt_gtilde_table::<i128>(n, k)
            .iter()
            .map(|(lambda, f)| (lambda.parts().to_vec().into(), qt_mon_out(f)))
            .collect())
    })
}

/// `G̃^(k)_λ` in the **Schur** basis.
///
/// The same object [`llt_gtilde`] returns in the monomial basis, so the
/// supports here are Schur shapes and not weights. Rows in the element order.
///
/// ```text
/// >>> symfn.llt_schur([2], 2)
/// [((1,), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition, `k ≥ 1`, and the shape fits
/// the abacus representation.
#[pyfunction]
#[pyo3(signature = (la, k))]
fn llt_schur(la: Vec<u32>, k: u32) -> PyResult<QtSchur> {
    interruptible(move || {
        let (l, k) = (part_arg(&la)?, level_arg(k)?);
        abacus_arg(&l, k)?;
        Ok(escalate(
            || {
                Some(qt_schur_out(&guarded(|| {
                    crate::llt::llt_schur::<Guarded>(&l, k)
                })?))
            },
            || qt_schur_out(&crate::llt::llt_schur::<BigInt>(&l, k)),
        ))
    })
}

/// `G_ν(x;q)` for a tuple of shapes, in the monomial basis and the **raw** inv
/// grading.
///
/// `offsets` defaults to all zero. ⚠️ The floor is **not** divided out: `min_T
/// inv(T)` can be positive, and Sage's `llt(k).cospin(tuple)` returns `q^{−min
/// inv} G_ν` instead. Divide by `q^{llt_min_inv(...)}` to compare — exposing
/// the floor is deliberate, since it is real data about ν and hiding it is how
/// the quotient dictionary gets misread (`docs/record/llt.md`).
///
/// `shapes` is a list of straight shapes and `offsets` shifts each component's
/// content, one integer per shape. Rows in the element order of the weight.
///
/// ```text
/// >>> symfn.llt_g([[1], [1]])
/// [((1, 1), [(0, 0, 1), (1, 0, 1)]), ((2,), [(0, 0, 1)])]
/// >>> symfn.llt_g([[1], [1]], [0, 1])
/// [((1, 1), [(0, 0, 2)]), ((2,), [(0, 0, 1)])]
/// ```
///
/// The offsets change the answer, which is what makes them part of the
/// argument rather than a normalization detail.
///
/// # Raises
///
/// Raises `ValueError` unless every shape is a partition, `offsets` has one
/// entry per shape, and the total cell count fits.
#[pyfunction]
#[pyo3(signature = (shapes, offsets=None))]
fn llt_g(shapes: Vec<Vec<u32>>, offsets: Option<Vec<i32>>) -> PyResult<QtMon> {
    interruptible(move || {
        Ok(qt_mon_out(&crate::llt::llt_g::<i128>(&skew_tuple(
            &shapes, offsets,
        )?)))
    })
}

/// `min_T inv(T)` over the semistandard fillings of a tuple — the forced
/// `q`-floor that [`llt_g`] does not divide out.
///
/// Returns one nonnegative `int`. Dividing [`llt_g`]'s answer by `q` to this
/// power is what recovers Sage's `llt(k).cospin(tuple)`.
///
/// ```text
/// >>> symfn.llt_min_inv([[1], [1]])
/// 0
/// ```
///
/// # Raises
///
/// Raises `ValueError` on the same conditions [`llt_g`] does.
#[pyfunction]
#[pyo3(signature = (shapes, offsets=None))]
fn llt_min_inv(shapes: Vec<Vec<u32>>, offsets: Option<Vec<i32>>) -> PyResult<u32> {
    interruptible(move || Ok(crate::llt::llt_min_inv(&skew_tuple(&shapes, offsets)?)))
}

/// The **fundamental quasisymmetric** expansion of `G_ν`, as
/// `[(composition, [(q_exp, t_exp, coeff), ...]), ...]`.
///
/// \[HHL\] (82)'s descent buckets read directly. The crate has no QSym type —
/// the compositions carry their own meaning and nothing here multiplies them.
///
/// The support is a **composition**, so its entries need not decrease, and
/// that is what separates this from the monomial expansion [`llt_g`] returns.
///
/// ```text
/// >>> symfn.llt_fundamental([[1], [1]])
/// [((1, 1), [(1, 0, 1)]), ((2,), [(0, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` on the same conditions [`llt_g`] does.
#[pyfunction]
#[pyo3(signature = (shapes, offsets=None))]
fn llt_fundamental(
    shapes: Vec<Vec<u32>>,
    offsets: Option<Vec<i32>>,
) -> PyResult<Vec<(Key, Vec<(u32, u32, Coeff)>)>> {
    interruptible(move || {
        Ok(
            crate::llt::llt_fundamental::<i128>(&skew_tuple(&shapes, offsets)?)
                .iter()
                .map(|(comp, c)| (comp.clone().into(), qt_poly(c)))
                .collect(),
        )
    })
}

/// The k-core and k-quotient of λ, as `(core, [component, ...])`.
///
/// The abacus primitives the ribbon model rests on. Component **order** (runner
/// 0 first) matters — `G_ν` is not symmetric in its components — and
/// agrees with Sage's `Partition(λ).quotient(k)`, which `scripts/check_llt.py`
/// checks.
///
/// The quotient always has exactly `k` components, some of them empty.
///
/// ```text
/// >>> symfn.k_core_quotient([3, 1], 2)
/// ((), [(2,), ()])
/// ```
///
/// The empty core is what makes `(3, 1)` a shape [`llt_gtilde`] answers on,
/// and the component order is runner 0 first — reversing it gives a different
/// `G_ν`.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition and `k ≥ 1`.
#[pyfunction]
#[pyo3(signature = (la, k))]
fn k_core_quotient(la: Vec<u32>, k: u32) -> PyResult<(Key, Vec<Key>)> {
    interruptible(move || {
        let (l, k) = (part_arg(&la)?, level_arg(k)?);
        Ok((
            l.k_core(k).parts().to_vec().into(),
            l.k_quotient(k)
                .iter()
                .map(|p| p.parts().to_vec().into())
                .collect(),
        ))
    })
}

/// `∇e_n = Σ_D t^{area(D)} G_D(x;q)`, as `[(area_sequence, G_D), ...]`.
///
/// The by-path Schur-positive refinement of the shuffle theorem — `∇e_n`
/// written as a positive sum of positive pieces. No package emits this
/// decomposition (`docs/record/dyck-paths.md`), and it is what makes the rise
/// side of the Delta conjecture cheap (see [`delta_conjecture_side`]).
///
/// ⚠️ `C_n` pieces and `#SYT` work each: n = 10 is 16 796 pieces.
/// Use [`nabla_e`] for the total, which is far cheaper.
///
/// Each entry is `(area sequence, G_D)`, the area sequence a tuple of `n`
/// integers and `G_D` a monomial-basis element in the [`llt_g`] encoding.
/// Summing `t^{area} · G_D` over the list reproduces [`nabla_e`].
///
/// ```text
/// >>> symfn.nabla_e_by_path(2)[1]
/// ((0, 1), [((1, 1), [(0, 1, 1)])])
/// ```
///
/// # Raises
///
/// Raises `ValueError` if the degree is past what the tuple representation
/// holds.
#[pyfunction]
fn nabla_e_by_path(n: u32) -> PyResult<Vec<(Key, QtMon)>> {
    interruptible(move || {
        cells_arg(n as usize, &format!("a degree-{n} path tuple"))?;
        Ok(crate::llt::nabla_e_by_path::<i128>(n)
            .iter()
            .map(|(area, g)| (area.clone().into(), qt_mon_out(g)))
            .collect())
    })
}

/// One **column** of the Schur-expansion table: `c^λ_μ` for every shape μ ⊢ k|λ|,
/// in the \[KMS\] variable `v`, as `[(mu, [(v_exp, 0, coeff), ...]), ...]`.
///
/// These are parabolic affine Kazhdan–Lusztig polynomials (\[LT\] Thm 4.2),
/// computed by exact Fock-space straightening with no Hecke algebra in sight.
/// Sage has no entry point of this shape.
///
/// ⚠️ The variable is `v`, and the ribbon side's grading is recovered at
/// **`q = −v`**. Coefficients are signed for that reason.
///
/// Rows are `(mu, [(v exponent, 0, coefficient), ...])` in `partitions(k|λ|)`
/// order with the zero entries omitted, the second slot always 0 since there
/// is no `t` here.
///
/// ```text
/// >>> symfn.llt_kl_column([2], 2)
/// [((4,), [(0, 0, 1)]), ((3, 1), [(1, 0, -1)]), ((2, 2), [(2, 0, 1)])]
/// ```
///
/// The alternating signs are the `v` convention; under `q = −v` they all
/// become positive, which is the check that this is not already the `q` form.
///
/// # Raises
///
/// Raises `ValueError` unless λ is a partition, `k ≥ 1`, and the shape fits
/// the abacus representation.
#[pyfunction]
#[pyo3(signature = (la, k))]
fn llt_kl_column(la: Vec<u32>, k: u32) -> PyResult<Vec<(Key, Vec<(u32, u32, Coeff)>)>> {
    interruptible(move || {
        let (l, k) = (part_arg(&la)?, level_arg(k)?);
        abacus_arg(&l, k)?;
        Ok(crate::llt::llt_kl_column::<i128>(&l, k)
            .iter()
            .map(|(mu, c)| (mu.parts().to_vec().into(), qt_poly(c)))
            .collect())
    })
}

/// `G_Γ(x;q) = Σ_κ q^{asc(κ)} x^κ` over the colorings of a decorated graph.
///
/// Vertices are `0 … n−1`. `weak` edges are ordered pairs whose *ascents* are
/// counted (`κ(u) < κ(v)` scores a `q`); `strict` edges must run `u < v` and
/// *constrain* (`κ(u) < κ(v)` or the coloring does not count), never scoring.
/// The two sets must be disjoint as unordered pairs, or the statistic silently
/// gains a `q` per strict edge — which is why that is raised rather than
/// tolerated.
///
/// Rows in the [`llt_g`] encoding, the support a weight.
///
/// ```text
/// >>> symfn.llt_graph(2, [(0, 1)], [])
/// [((1, 1), [(0, 0, 1), (1, 0, 1)]), ((2,), [(0, 0, 1)])]
/// ```
///
/// The `(2,)` term is the monochromatic coloring, which a strict edge would
/// remove — the value that separates the two edge kinds.
///
/// # Raises
///
/// Raises `ValueError` if a strict edge runs `u ≥ v`, if an endpoint is not
/// below `n`, or if an edge appears in both sets.
#[pyfunction]
#[pyo3(signature = (n, weak, strict))]
fn llt_graph(n: u32, weak: Vec<(u32, u32)>, strict: Vec<(u32, u32)>) -> PyResult<QtMon> {
    interruptible(move || {
        Ok(qt_mon_out(&crate::llt::llt_graph::<i128>(
            &decorated_graph(n, weak, strict)?,
        )))
    })
}

/// The Shareshian–Wachs chromatic quasisymmetric function `X_Γ(x;q)` of Γ, from
/// its LLT polynomial by the `(q−1)`-plethysm of \[CM\] Prop 3.5.
///
/// Γ should carry the \[CM\] presentation — natural orientation, no strict
/// edges — for the answer to be the chromatic function of the graph rather than
/// of a decorated relative of it. **Isolated vertices are part of Γ and must be
/// counted in `n`**; dropping them is a well-trodden route to a plausible wrong
/// answer.
///
/// Runs over ℚ because the plethysm goes through the power-sum basis, which
/// divides by `z_ρ`; the answer is integral by the time it crosses, and a
/// surviving denominator is raised rather than rounded.
///
/// Rows in the [`llt_g`] encoding.
///
/// ```text
/// >>> symfn.chromatic_from_llt(2, [(0, 1)], [])
/// [((1, 1), [(0, 0, 1), (1, 0, 1)])]
/// ```
///
/// The absent `(2,)` term is the point: a proper coloring of an edge cannot
/// be monochromatic, where [`llt_graph`] of the same Γ keeps it.
///
/// # Raises
///
/// Raises `ValueError` on the same three edge conditions [`llt_graph`] does,
/// and if the plethysm leaves a denominator.
#[pyfunction]
#[pyo3(signature = (n, weak, strict))]
fn chromatic_from_llt(n: u32, weak: Vec<(u32, u32)>, strict: Vec<(u32, u32)>) -> PyResult<QtMon> {
    interruptible(move || {
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
            out.push((mu.parts().to_vec().into(), row));
        }
        Ok(out)
    })
}

/// The \[AS\] **e-expansion** of `Ĝ_Γ(x; q+1)`: `Σ_θ q^{asc(θ)} e_{λ(θ)}` over
/// orientations of the free edges, as `[(partition, poly), ...]`.
///
/// By \[DA\]'s theorem the coefficients are non-negative, so this is
/// **certified positive output** rather than a conjecture to check. Weak edges
/// are read as unordered pairs: \[AS\]'s formula orients them itself.
///
/// ⚠️ `2^{#free edges}` terms.
///
/// Rows are `(partition, poly)` in the element order, the support an
/// **elementary** index rather than a weight.
///
/// ```text
/// >>> symfn.llt_e_expansion(2, [(0, 1)], [])
/// [((1, 1), [(0, 0, 1)]), ((2,), [(1, 0, 1)])]
/// ```
///
/// # Raises
///
/// Raises `ValueError` on the same three edge conditions [`llt_graph`] does,
/// and if the graph has more free edges than the orientation mask holds.
#[pyfunction]
#[pyo3(signature = (n, weak, strict))]
fn llt_e_expansion(
    n: u32,
    weak: Vec<(u32, u32)>,
    strict: Vec<(u32, u32)>,
) -> PyResult<Vec<(Key, Vec<(u32, u32, Coeff)>)>> {
    interruptible(move || {
        let g = decorated_graph(n, weak, strict)?;
        let free = crate::llt::free_edges(&g).len();
        if free >= crate::llt::MAX_FREE_EDGES {
            return Err(PyValueError::new_err(format!(
                "the orientation sum over {free} free edges is 2^{free} terms, \
                 past the {} this mask holds",
                crate::llt::MAX_FREE_EDGES
            )));
        }
        Ok(crate::llt::llt_e_expansion::<i128>(&g)
            .iter()
            .map(|(lambda, c)| (lambda.parts().to_vec().into(), qt_poly(c)))
            .collect())
    })
}

/// `H̃_μ(x;q,t) = Σ_D q^{−a(D)} t^{maj(D)} G_{ν(μ,D)}(x;q)` — the \[HHL\]
/// decomposition, in the monomial basis.
///
/// The fourth route to `H̃` in this crate and the only one positively graded at
/// every intermediate step: Macdonald positivity *is* LLT positivity, and this
/// is where that becomes a computation.
///
/// ⚠️ A reference route, not a fast one — `2^{|μ|−μ₁}` LLT evaluations. Use
/// [`macdonald_ht`] to *compute* `H̃`; use this to check it independently.
///
/// Rows in the [`llt_g`] encoding, and here both exponent slots are used:
/// this is the one LLT-side entry point carrying `t`.
///
/// ```text
/// >>> symfn.htilde_by_llt([2])
/// [((1, 1), [(0, 0, 1), (1, 0, 1)]), ((2,), [(0, 0, 1)])]
/// ```
///
/// The monomial basis is what separates this from [`macdonald_ht`], which
/// answers the same question in the Schur basis.
///
/// # Raises
///
/// Raises `ValueError` unless μ is a partition whose size fits the tuple
/// representation.
#[pyfunction]
fn htilde_by_llt(mu: Vec<u32>) -> PyResult<QtMon> {
    interruptible(move || {
        let m = part_arg(&mu)?;
        cells_arg(m.size() as usize, &format!("the shape {m}"))?;
        Ok(qt_mon_out(&crate::llt::htilde_by_llt::<i128>(&m)))
    })
}

/// A kernel for computing with symmetric functions, exactly.
///
/// Every function here does a **whole-object** operation — multiply two
/// complete elements, convert an entire expansion, return a whole table — and
/// there are no per-monomial accessors, because a cross-language call has
/// overhead and looping one is how a caller loses the speed this library
/// exists for. Where the natural unit of work is larger than one value, the
/// larger unit is the entry point: `kostka_foulkes_column`, the `*_table`
/// family, and `convert_indexed`, which keys by index into `partitions(n)`
/// rather than by lists of parts.
///
/// ## How data crosses
///
/// An **element** is a list of `(support, coefficient)` pairs, keyed by
/// whatever indexes its basis — a partition for a symmetric function, a
/// permutation in one-line notation for a Schubert polynomial:
///
/// ```text
/// >>> symfn.schur_multiply([([2], 1)], [([1], 1)])
/// [((2, 1), 1), ((3,), 1)]
/// ```
///
/// **Every element comes back in one order**: increasing lexicographic by
/// support, with no zero coefficients and no repeated key. So `(2, 1)`
/// precedes `(3,)`, and a support absent from the list has coefficient zero
/// rather than an unknown value. Where a return value is not an element — a
/// coefficient list, a table, an enumeration — the entry point states its own
/// order.
///
/// A **parameter family** crosses as exponent-keyed rows rather than as a
/// polynomial object: `(exponent, coefficient)` for one variable, and
/// `(a, b, coefficient)` for `q^a t^b`. Nothing here returns a type you must
/// import something to unpack, and nothing assumes a coefficient ring on the
/// far side — which is what lets Sage, SymPy and a bare interpreter each
/// rebuild elements in their own ring.
///
/// Partitions are given as weakly decreasing lists of positive integers.
/// Trailing zeros are tolerated, because that is the fixed-width form Sage
/// hands over; anything else malformed raises `ValueError` naming the
/// requirement it violated.
///
/// ## Coefficients have no ceiling
///
/// They cross as Python `int`s of arbitrary size in both directions, and no
/// value returned is ever rounded, truncated or wrapped. Internally a call
/// runs over a fixed-width type that *reports* overflow and, if anything
/// overflowed, again over arbitrary precision — so the width is an
/// implementation detail rather than a wall a caller can hit. Where a
/// computation cannot be exact it raises instead of approximating.
///
/// ## Sage
///
/// This module does not import Sage, depend on it, or know it exists, and it
/// works in any CPython 3.9+. The adapter that makes Sage use it lives on the
/// Sage side of the boundary; see `docs/policies/python.md`, which is this
/// surface's rulebook.
#[pymodule]
fn symfn(m: &Bound<'_, PyModule>) -> PyResult<()> {
    // Sourced from the crate version so the two cannot drift
    // (`docs/policies/python.md`, P10).
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    // At import, so that no entry point can run before the kernel has a way to
    // notice Ctrl-C. Both calls are idempotent under a re-import.
    crate::interrupt::set_checker(check_python_signals);
    install_quiet_cancellation_hook();
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
    m.add_function(wrap_pyfunction!(schur_in_macdonald_j, m)?)?;
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
    m.add_function(wrap_pyfunction!(st_multiply, m)?)?;
    m.add_function(wrap_pyfunction!(reduced_kronecker_product, m)?)?;
    m.add_function(wrap_pyfunction!(reduced_kronecker, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_st, m)?)?;
    m.add_function(wrap_pyfunction!(st_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(ht_multiply, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_ht, m)?)?;
    m.add_function(wrap_pyfunction!(ht_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(lr_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_homogeneous, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_elementary, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_monomial, m)?)?;
    m.add_function(wrap_pyfunction!(schur_to_forgotten, m)?)?;
    m.add_function(wrap_pyfunction!(to_power, m)?)?;
    m.add_function(wrap_pyfunction!(homogeneous_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(elementary_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(monomial_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(forgotten_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(power_to_schur, m)?)?;
    m.add_function(wrap_pyfunction!(plethysm, m)?)?;
    m.add_function(wrap_pyfunction!(kostka_number, m)?)?;
    m.add_function(wrap_pyfunction!(evaluate_schur, m)?)?;
    m.add_function(wrap_pyfunction!(expand_alphabet, m)?)?;
    m.add_function(wrap_pyfunction!(monomial_multiply, m)?)?;
    m.add_function(wrap_pyfunction!(semistandard_tableaux, m)?)?;
    m.add_function(wrap_pyfunction!(dimension, m)?)?;
    m.add_function(wrap_pyfunction!(principal_specialization, m)?)?;
    m.add_function(wrap_pyfunction!(principal_specialization_q, m)?)?;
    m.add_function(wrap_pyfunction!(character_value, m)?)?;
    m.add_function(wrap_pyfunction!(internal_product, m)?)?;
    m.add_function(wrap_pyfunction!(partitions, m)?)?;
    m.add_function(wrap_pyfunction!(kronecker_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(convert_indexed, m)?)?;
    m.add_function(wrap_pyfunction!(convert_terms, m)?)?;
    m.add_function(wrap_pyfunction!(character_table, m)?)?;
    m.add_function(wrap_pyfunction!(kostka_table, m)?)?;
    m.add_function(wrap_pyfunction!(omega, m)?)?;
    m.add_function(wrap_pyfunction!(hall_inner_product, m)?)?;
    m.add_function(wrap_pyfunction!(skew_schur, m)?)?;
    m.add_function(wrap_pyfunction!(skew_by, m)?)?;
    m.add_function(wrap_pyfunction!(coproduct, m)?)?;
    m.add_function(wrap_pyfunction!(antipode, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_multiply, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_multiply_variable, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_divided_difference, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_divided_difference_perm, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_expand, m)?)?;
    m.add_function(wrap_pyfunction!(polynomial_to_schubert, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_pairing, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_scalar_product, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_dimension, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_coefficient, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_monomial_mass, m)?)?;
    m.add_function(wrap_pyfunction!(schubert_to_stanley_schur, m)?)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Below the wall the fast pass answers, so the ladder must not change any
    /// value a caller already had — and must return `Small` coefficients, i.e.
    /// it really did take the fixed-width route.
    #[test]
    fn llt_h_below_the_wall_is_unchanged_and_stays_narrow() {
        let got = llt_h(vec![1; 8], 3).expect("(1^8) at level 3 is a valid abacus argument");
        let want = qt_mon_out(&crate::llt::llt_h::<BigInt>(&Partition::new([1; 8]), 3));
        assert_eq!(got.len(), want.len());
        assert!(
            got.iter()
                .flat_map(|(_, p)| p.iter())
                .all(|(_, _, c)| matches!(c, Coeff::Small(_))),
            "below the wall nothing should have escalated"
        );
    }
}
