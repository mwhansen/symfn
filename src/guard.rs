//! Fixed-width coefficients that **report** overflow instead of wrapping, and
//! the scope that turns a report into a re-run.
//!
//! `impl Ring for i128` uses plain `*`, so a coefficient past the fixed width
//! aborts the computation rather than continuing. `overflow-checks = true`
//! reaches every profile (`docs/policies/failure.md`, R3), so that abort is a
//! panic and never a wrapped answer. A panic is still not an escalation.
//! [`guarded`] watches for `None`, so a panic escapes it and crashes the caller
//! on an input the wide pass answers exactly. The characters module already
//! refuses to do that ([`try_character`] returns
//! `None`, [`character_in`] re-runs the recursion in the coefficient ring), and
//! this module generalizes that pattern to every coefficient:
//!
//! ```text
//!   guarded(|| .. compute over Guarded ..)   ->  Some(answer)  or  None
//! ```
//!
//! `None` means "some operation left the fixed width"; the caller then re-runs
//! over `BigInt` and gets an exact answer. The fast path
//! is unchanged — measured at 0–1% against unchecked arithmetic
//! (`examples/bench_guarded.rs`) — because a `checked_mul` is one predictable
//! branch and the expensive path is never taken.
//!
//! [`try_character`]: crate::character::try_character
//! [`character_in`]: crate::character::character_in
//!
//! ## Why a counter and not a flag
//!
//! Overflow is recorded by **incrementing a global counter**, and [`guarded`]
//! compares its value before and after. The obvious alternative — clear a flag,
//! run, test the flag — is wrong under concurrency in a way that loses answers:
//! if two computations overlap, one can clear the flag *after* the other set
//! it, and the second then reports success on a wrapped result. A monotone
//! counter cannot do that. Its failure direction is the safe one: an unrelated
//! thread's overflow makes this call escalate needlessly, costing time and
//! never correctness.
//!
//! That is also why this is global rather than thread-local. A thread-local
//! counter would miss an overflow on a worker thread — under-reporting, the
//! unsafe direction. Nothing in the parallel Littlewood–Richardson path carries
//! `Ring` coefficients today (it accumulates `u128` counts), but the guarantee
//! should not depend on that staying true.

// Every `u128 → i128` here is guarded by an explicit `n > i128::MAX as u128`
// test immediately above it, which is what those functions are for.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

use core::sync::atomic::{AtomicU64, Ordering};

use crate::coeff::{
    rat_add, rat_div, rat_mul, rat_normalize, Overflow, Plethystic, QAlgebra, Ring,
};

static OVERFLOWS: AtomicU64 = AtomicU64::new(0);

#[cold]
#[inline(never)]
fn note_overflow() {
    OVERFLOWS.fetch_add(1, Ordering::Relaxed);
}

/// Run `f` over guarded coefficients, returning `None` if anything overflowed.
///
/// `None` is not an error: it is the signal to re-run the same computation over
/// an arbitrary-precision ring. A `Some` is a promise that every intermediate
/// stayed inside the fixed width, so the answer is exact.
///
/// # The closure's obligation
///
/// That promise is **the closure's to keep, not this function's**. All this
/// does is compare a counter before and after; an operation that never calls
/// `note_overflow` is invisible to it. So every arithmetic operation reachable
/// from `f` must do one of three things (`docs/policies/failure.md`, R6):
///
/// 1. run through [`Guarded`] / [`GuardedRat`], which report; or
/// 2. be `checked_*` with an explicit refusal or fallback on `None`; or
/// 3. carry a bound proof at the site.
///
/// Native arithmetic that does none of these breaks the ladder in one of two
/// ways, and the second is easy to miss because it is *loud*:
///
/// * without `overflow-checks` it wraps, and this returns `Some(garbage)`;
/// * with `overflow-checks` — which the release profile now carries — it
///   **panics**, and a panic is not something the `escalate` two-pass helper
///   can catch. The caller gets a crash where the wide pass had the answer.
///
/// `powersum_scalar` was the live instance: it formed z_μ in native `u128`, so
/// `jack_scalar` at |μ| ≥ 35 panicked instead of escalating, while |μ| = 34
/// reported correctly because only the *injection* was too wide. It is pinned
/// by `the_z_wall_reports_inside_a_guarded_scope_rather_than_panicking`
/// (`tests/bignum.rs`).
pub fn guarded<T>(f: impl FnOnce() -> T) -> Option<T> {
    let before = OVERFLOWS.load(Ordering::Relaxed);
    let value = f();
    (OVERFLOWS.load(Ordering::Relaxed) == before).then_some(value)
}

/// How many overflows have been recorded process-wide. Diagnostics only.
pub fn overflow_count() -> u64 {
    OVERFLOWS.load(Ordering::Relaxed)
}

/// An `i128` that records overflow rather than wrapping.
///
/// An operation that leaves the width increments the counter rather than
/// storing its result. So `.0` is an answer only inside a [`guarded`] scope
/// that returned `Some`.
#[derive(Clone, Copy, PartialEq, Eq, Debug, Default)]
pub struct Guarded(pub i128);

impl Ring for Guarded {
    #[inline]
    fn zero() -> Self {
        Guarded(0)
    }
    #[inline]
    fn one() -> Self {
        Guarded(1)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
    #[inline]
    fn add_assign(&mut self, other: &Self) {
        match self.0.checked_add(other.0) {
            Some(v) => self.0 = v,
            None => note_overflow(),
        }
    }
    #[inline]
    fn mul(&self, other: &Self) -> Self {
        match self.0.checked_mul(other.0) {
            Some(v) => Guarded(v),
            None => {
                note_overflow();
                Guarded(0)
            }
        }
    }
    /// Reports on `i128::MIN`, which has no negation inside the width.
    /// `checked_mul` can legitimately land on it, and `wrapping_neg` there is a
    /// silently wrong *sign* — `Some(garbage)` out of a scope whose whole
    /// promise is that it never returns one.
    #[inline]
    fn neg(&self) -> Self {
        match self.0.checked_neg() {
            Some(v) => Guarded(v),
            None => {
                note_overflow();
                Guarded(0)
            }
        }
    }
    #[inline]
    fn from_i64(n: i64) -> Self {
        Guarded(n as i128)
    }
    #[inline]
    fn from_u128(n: u128) -> Self {
        if n > i128::MAX as u128 {
            note_overflow();
            return Guarded(0);
        }
        Guarded(n as i128)
    }
    #[inline]
    fn from_i128(n: i128) -> Self {
        Guarded(n)
    }
    /// Exact in ℤ, and it cannot overflow: `|a/b| ≤ |a|` whenever the division
    /// is exact.
    ///
    /// Needed because [`AFrac`](crate::afrac::AFrac) cancels its atoms by
    /// synthetic division through `div_exact`. Declining — the [`Ring`] default
    /// — is *safe* for every other caller, since "not divisible" is an ordinary
    /// outcome there, but for `AFrac` it would silently turn `reduce` into a
    /// no-op and let denominators grow without bound. So the guarded ladder for
    /// Jack needs this and the plain `i128` one does not notice it.
    #[inline]
    fn div_exact(&self, other: &Self) -> Option<Self> {
        // `MIN / -1` is the one exact quotient that leaves the width; it is
        // reported through `neg`, as every other sign flip here is.
        if other.0 == -1 {
            return Some(self.neg());
        }
        crate::coeff::div_exact_i128(self.0, other.0).map(Guarded)
    }
}

/// An exact rational over guarded `i128`s, in lowest terms with `den > 0`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct GuardedRat {
    num: i128,
    den: i128,
}

/// [`GuardedRat`]'s [`Overflow`] policy: an operation that leaves the width
/// increments the counter and yields `0`, and [`GuardedRat::reduced`] then
/// refuses the pair it lands in. The arithmetic itself is `coeff::rat_*`,
/// shared with [`Rational`](crate::coeff::Rational).
struct Reporting;

impl Overflow for Reporting {
    #[inline]
    fn add(a: i128, b: i128) -> i128 {
        ck(a.checked_add(b))
    }
    #[inline]
    fn mul(a: i128, b: i128) -> i128 {
        ck(a.checked_mul(b))
    }
}

#[inline]
fn ck(v: Option<i128>) -> i128 {
    match v {
        Some(x) => x,
        None => {
            note_overflow();
            0
        }
    }
}

impl GuardedRat {
    /// Normalizes to lowest terms with `den > 0`, reporting instead of
    /// returning a value it cannot represent.
    ///
    /// `i128::MIN` is reported rather than stored: normalizing needs `|num|`,
    /// `|den|` and possibly a sign flip, none of which `MIN` has. A report costs
    /// one needless escalation on a measure-zero input; storing it costs a
    /// wrong sign in a scope that promised exactness.
    fn new(num: i128, den: i128) -> Self {
        if den == 0 {
            note_overflow();
            return GuardedRat { num: 0, den: 1 };
        }
        if num == 0 {
            return GuardedRat { num: 0, den: 1 };
        }
        if num == i128::MIN || den == i128::MIN {
            note_overflow();
            return GuardedRat { num: 0, den: 1 };
        }
        let (num, den) = rat_normalize(num, den);
        GuardedRat { num, den }
    }

    /// A pair the shared arithmetic has already put in lowest terms with
    /// `den > 0` — a second gcd would learn nothing — passed through the two
    /// refusals [`GuardedRat::new`] makes. `den ≤ 0` can only come from a
    /// `ck` that reported and returned `0`, and a `MIN` part is what a checked
    /// product can legitimately land on but the representation does not hold.
    /// Every value of this type is built through here or through `new`.
    #[inline]
    fn reduced((num, den): (i128, i128)) -> Self {
        if den <= 0 || num == i128::MIN {
            note_overflow();
            return GuardedRat { num: 0, den: 1 };
        }
        if num == 0 {
            return GuardedRat { num: 0, den: 1 };
        }
        GuardedRat { num, den }
    }
    /// The numerator, in lowest terms. The sign of the rational lives here.
    ///
    /// A value whose construction left the width reads as `0/1`, because the
    /// report is the scope's overflow flag rather than the pair: read it with
    /// [`guarded`], never by inspecting these two.
    pub fn numer(&self) -> i128 {
        self.num
    }
    /// The denominator, in lowest terms and always positive; 1 for an integer,
    /// and 1 for zero.
    pub fn denom(&self) -> i128 {
        self.den
    }
}

impl Ring for GuardedRat {
    fn zero() -> Self {
        GuardedRat { num: 0, den: 1 }
    }
    fn one() -> Self {
        GuardedRat { num: 1, den: 1 }
    }
    fn is_zero(&self) -> bool {
        self.num == 0
    }
    fn add_assign(&mut self, other: &Self) {
        *self = GuardedRat::reduced(rat_add::<Reporting>(
            self.num, self.den, other.num, other.den,
        ));
    }
    fn mul(&self, other: &Self) -> Self {
        GuardedRat::reduced(rat_mul::<Reporting>(
            self.num, self.den, other.num, other.den,
        ))
    }
    /// Reports on a numerator of `i128::MIN` rather than wrapping its sign —
    /// see [`Guarded::neg`]. Normalization already refuses to store one, so
    /// this only fires on a value built past it, and costs one compare.
    fn neg(&self) -> Self {
        match self.num.checked_neg() {
            Some(num) => GuardedRat { num, den: self.den },
            None => {
                note_overflow();
                <Self as Ring>::zero()
            }
        }
    }
    fn from_i64(n: i64) -> Self {
        GuardedRat {
            num: n as i128,
            den: 1,
        }
    }
    fn from_u128(n: u128) -> Self {
        if n > i128::MAX as u128 {
            note_overflow();
            return <Self as Ring>::zero();
        }
        GuardedRat {
            num: n as i128,
            den: 1,
        }
    }
    fn from_i128(n: i128) -> Self {
        // `MIN` cannot be normalized later, so it is refused at the seam rather
        // than stored and reported by whichever operation first needs its sign.
        if n == i128::MIN {
            note_overflow();
            return <Self as Ring>::zero();
        }
        GuardedRat { num: n, den: 1 }
    }
    // Kept, so `convert::integral_sweep` still applies. Its own internal
    // overflow is *not* an overflow of the answer — it bails to the generic
    // path — and it works in raw `i128`, never through this impl, so it cannot
    // trip the counter.
    fn as_ratio(&self) -> Option<(i128, i128)> {
        Some((self.num, self.den))
    }
    fn from_ratio(num: i128, den: i128) -> Option<Self> {
        (den != 0).then(|| GuardedRat::new(num, den))
    }
}

impl QAlgebra for GuardedRat {
    fn div_u128(&self, n: u128) -> Self {
        if n == 0 {
            note_overflow();
            return <Self as Ring>::zero();
        }
        if n > i128::MAX as u128 {
            note_overflow();
            return <Self as Ring>::zero();
        }
        let n = n as i128;
        // `self.num` is never `MIN` (`new` and `from_i128` both refuse one) and
        // `n` is at most `i128::MAX` by the check above.
        GuardedRat::reduced(rat_div::<Reporting>(self.num, self.den, n))
    }
}

impl Plethystic for GuardedRat {
    /// ℚ has no variables to raise.
    fn frobenius(&self, _n: u32) -> Self {
        *self
    }
}

/// The lock every test that constructs a [`Guarded`] must hold.
///
/// The counter is global and monotone, so a test that overflows on purpose
/// makes any *concurrently* running [`guarded`] scope report `None` too. That
/// is the intended fail-safe direction — see the module docs — but it means
/// such tests must not run in parallel with each other, wherever they live:
/// `memo::tests` builds a reporting value to check that an overflowing
/// computation is not cached.
#[cfg(test)]
pub(crate) fn serial() -> std::sync::MutexGuard<'static, ()> {
    static SERIAL: std::sync::Mutex<()> = std::sync::Mutex::new(());
    // Poisoning ignored so one failing test does not cascade into the rest.
    SERIAL.lock().unwrap_or_else(|e| e.into_inner())
}

#[cfg(test)]
mod tests {
    use super::serial;
    use super::*;
    use crate::coeff::Rational;
    use crate::convert::FromSchur;
    use crate::partition::Partition;
    use crate::sym::{PowerSum, Schur, SymFn};

    #[test]
    fn arithmetic_inside_the_range_reports_nothing() {
        let _g = serial();
        let r = guarded(|| {
            let a = Guarded(1_000_000);
            let b = a.mul(&a);
            let mut c = b;
            c.add_assign(&b);
            c
        });
        assert_eq!(r, Some(Guarded(2_000_000_000_000)));
    }

    #[test]
    fn multiplication_past_the_width_is_reported() {
        let _g = serial();
        let r = guarded(|| Guarded(i128::MAX / 2).mul(&Guarded(4)));
        assert_eq!(r, None, "overflowing product must not be returned");
    }

    #[test]
    fn addition_past_the_width_is_reported() {
        let _g = serial();
        let r = guarded(|| {
            let mut a = Guarded(i128::MAX);
            a.add_assign(&Guarded(1));
            a
        });
        assert_eq!(r, None);
    }

    #[test]
    fn u128_injection_past_i128_is_reported() {
        let _g = serial();
        let r = guarded(|| <Guarded as Ring>::from_u128(u128::MAX));
        assert_eq!(r, None, "a structure constant past i128 must escalate");
    }

    /// `i128::MIN` is the one value a `checked_mul` can return that the *rest*
    /// of the arithmetic cannot handle: `-MIN` and `MIN.abs()` both leave the
    /// width. Negating it used to wrap, which put a silently wrong sign inside
    /// a scope whose whole promise is that it never returns one.
    #[test]
    fn negating_the_width_minimum_is_reported_rather_than_wrapped() {
        let _g = serial();
        let reached = guarded(|| Guarded(i128::MIN / 2).mul(&Guarded(2)));
        assert_eq!(
            reached,
            Some(Guarded(i128::MIN)),
            "MIN is a legitimate product, not an overflow"
        );

        let r = guarded(|| Guarded(i128::MIN).neg());
        assert_eq!(r, None, "negating MIN must escalate, not wrap to MIN");
    }

    /// The rational side of the same corner: normalizing needs `|num|`, `|den|`
    /// and possibly a sign flip, and `MIN` has none of them.
    #[test]
    fn rational_width_minimum_is_reported_at_every_seam() {
        let _g = serial();
        assert_eq!(
            guarded(|| <GuardedRat as Ring>::from_i128(i128::MIN)),
            None,
            "injecting MIN must escalate"
        );
        assert_eq!(
            guarded(|| GuardedRat::new(i128::MIN, 3)),
            None,
            "a MIN numerator must escalate"
        );
        assert_eq!(
            guarded(|| GuardedRat::new(3, i128::MIN)),
            None,
            "a MIN denominator must escalate"
        );
        // gcd(MIN, MIN) = 2^127, which is not an i128 either: the magnitudes are
        // taken in u128, so this reports rather than overflowing inside `gcd`.
        assert_eq!(guarded(|| GuardedRat::new(i128::MIN, i128::MIN)), None);
    }

    /// The counter is monotone, so a report inside a nested scope propagates to
    /// the outer one. Anything else would let an inner call swallow the signal.
    #[test]
    fn overflow_inside_a_nested_scope_still_reaches_the_outer_one() {
        let _g = serial();
        let outer = guarded(|| {
            let inner = guarded(|| Guarded(i128::MAX).mul(&Guarded(2)));
            assert_eq!(inner, None);
            7
        });
        assert_eq!(outer, None, "outer scope must not report success");
    }

    #[test]
    fn rationals_stay_exact_and_report_at_the_edge() {
        let _g = serial();
        let ok = guarded(|| {
            let a = GuardedRat::new(1, 3);
            let b = GuardedRat::new(1, 6);
            let mut c = a;
            c.add_assign(&b);
            c
        });
        assert_eq!(ok, Some(GuardedRat::new(1, 2)));

        let bad = guarded(|| GuardedRat::new(i128::MAX, 1).mul(&GuardedRat::new(3, 1)));
        assert_eq!(bad, None);
    }

    /// The guarded shortcuts must give exactly what `Rational`'s give on the
    /// same grid `coeff`'s test uses — every branch of Henrici's addition and
    /// the cross-cancelled product, both signs, zeros, and integers — and
    /// report nothing while doing so.
    #[test]
    fn guarded_shortcut_arithmetic_matches_rational() {
        let _g = serial();
        let nums = [-7i128, -6, -1, 0, 1, 2, 3, 5, 6, 12, 35, 1 << 40];
        let dens = [1i128, 2, 3, 4, 6, 7, 12, 30, 1 << 40, (1 << 40) + 1];
        let got = guarded(|| {
            for &a in &nums {
                for &b in &dens {
                    for &c in &nums {
                        for &d in &dens {
                            let (x, y) = (GuardedRat::new(a, b), GuardedRat::new(c, d));
                            let (rx, ry) = (Rational::new(a, b), Rational::new(c, d));
                            let mut sum = x;
                            sum.add_assign(&y);
                            let mut rsum = rx;
                            rsum.add_assign(&ry);
                            assert_eq!(
                                (sum.num, sum.den),
                                (rsum.numer(), rsum.denom()),
                                "{a}/{b} + {c}/{d}"
                            );
                            let (prod, rprod) = (x.mul(&y), rx.mul(&ry));
                            assert_eq!(
                                (prod.num, prod.den),
                                (rprod.numer(), rprod.denom()),
                                "{a}/{b} * {c}/{d}"
                            );
                            for n in [1u128, 2, 6, 7, 1 << 40] {
                                let (q, rq) = (x.div_u128(n), rx.div_u128(n));
                                assert_eq!(
                                    (q.num, q.den),
                                    (rq.numer(), rq.denom()),
                                    "{a}/{b} / {n}"
                                );
                            }
                        }
                    }
                }
            }
        });
        assert!(got.is_some(), "nothing on this grid leaves i128");
    }

    /// The shortcuts build values without a second normalization, so the two
    /// things `new` refuses must be refused on that path too: a product that
    /// lands exactly on `MIN`, and a denominator a reported overflow left at
    /// zero — neither may come out of the scope as a value.
    #[test]
    fn shortcut_paths_refuse_what_normalization_refuses() {
        let _g = serial();
        let half_min = GuardedRat::new(i128::MIN / 2, 1);
        assert_eq!(
            guarded(|| half_min.mul(&GuardedRat::new(2, 1))),
            None,
            "an integer product landing on MIN is reported"
        );
        assert_eq!(
            guarded(|| half_min.mul(&GuardedRat::new(2, 3))),
            None,
            "a cross-cancelled numerator landing on MIN is reported"
        );
        assert_eq!(
            guarded(|| {
                let mut x = GuardedRat::new(1, i128::MAX / 2);
                x.add_assign(&GuardedRat::new(1, i128::MAX / 2 - 1));
                x
            }),
            None,
            "coprime denominators whose product leaves i128 are reported"
        );
    }

    /// A real conversion over the guarded ring must agree with the unguarded
    /// one wherever it succeeds — the point being that guarding changes only
    /// the reporting, never a value.
    #[test]
    fn guarded_conversions_agree_with_the_plain_ring() {
        let _g = serial();
        for n in 0..=7u32 {
            for lambda in crate::partitions_of(n) {
                let g: Schur<GuardedRat> =
                    Schur::monomial(lambda.clone(), <GuardedRat as Ring>::one());
                let pg: PowerSum<GuardedRat> =
                    guarded(|| PowerSum::from_schur(&g)).expect("degree 7 fits i128");

                let r: Schur<crate::Rational> =
                    Schur::monomial(lambda.clone(), crate::Rational::from_int(1));
                let pr: PowerSum<crate::Rational> = PowerSum::from_schur(&r);

                assert_eq!(pg.terms().len(), pr.terms().len(), "support at {lambda}");
                for (mu, c) in pg.terms() {
                    let want = pr.coeff(mu);
                    assert_eq!(
                        (c.numer(), c.denom()),
                        (want.numer(), want.denom()),
                        "s_{lambda} -> p at {mu}"
                    );
                }
            }
        }
    }

    /// Characters past `i128` must escalate rather than come back wrapped. This
    /// is the case the whole mechanism exists for.
    #[test]
    fn a_character_past_i128_is_reported() {
        let _g = serial();
        let lam = Partition::new([18u32, 16, 14, 12, 10, 8]);
        let n = lam.size();
        let ones = Partition::new(std::iter::repeat_n(1, n as usize));
        let got = guarded(|| crate::character::character_in::<Guarded>(&lam, &ones));
        assert_eq!(got, None, "χ^λ(1^n) for |λ| = {n} exceeds i128");
    }
}
