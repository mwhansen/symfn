//! Fixed-width coefficients that **report** overflow instead of wrapping, and
//! the scope that turns a report into a re-run.
//!
//! `impl Ring for i128` uses plain `*`, so in release a coefficient past the
//! fixed width wraps silently — an answer that is wrong with no signal. The
//! characters module already refuses to do that ([`try_character`] returns
//! `None`, [`character_in`] re-runs the recursion in the coefficient ring), and
//! this module generalises that pattern to every coefficient:
//!
//! ```text
//!   guarded(|| .. compute over Guarded ..)   ->  Some(answer)  or  None
//! ```
//!
//! `None` means "some operation left the fixed width"; the caller then re-runs
//! over [`BigInt`](num_bigint::BigInt) and gets an exact answer. The fast path
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
//! if two computations overlap, one can clear the flag *after* the other set it,
//! and the second then reports success on a wrapped result. A monotone counter
//! cannot do that. Its failure direction is the safe one: an unrelated thread's
//! overflow makes this call escalate needlessly, costing time and never
//! correctness.
//!
//! That is also why this is global rather than thread-local. A thread-local
//! counter would miss an overflow on a worker thread — under-reporting, the
//! unsafe direction. Nothing in the parallel Littlewood–Richardson path carries
//! `Ring` coefficients today (it accumulates `u128` counts), but the guarantee
//! should not depend on that staying true.

use core::sync::atomic::{AtomicU64, Ordering};

use crate::coeff::{Plethystic, QAlgebra, Ring};

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
    #[inline]
    fn neg(&self) -> Self {
        Guarded(self.0.wrapping_neg())
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
}

/// An exact rational over guarded `i128`s, in lowest terms with `den > 0`.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub struct GuardedRat {
    num: i128,
    den: i128,
}

#[inline]
fn gcd(mut a: i128, mut b: i128) -> i128 {
    a = a.abs();
    b = b.abs();
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
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
    fn new(num: i128, den: i128) -> Self {
        if den == 0 {
            note_overflow();
            return GuardedRat { num: 0, den: 1 };
        }
        if num == 0 {
            return GuardedRat { num: 0, den: 1 };
        }
        let g = gcd(num, den);
        let (mut n, mut d) = (num / g, den / g);
        if d < 0 {
            n = -n;
            d = -d;
        }
        GuardedRat { num: n, den: d }
    }
    pub fn numer(&self) -> i128 {
        self.num
    }
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
        let a = ck(self.num.checked_mul(other.den));
        let b = ck(other.num.checked_mul(self.den));
        let num = ck(a.checked_add(b));
        let den = ck(self.den.checked_mul(other.den));
        *self = GuardedRat::new(num, den);
    }
    fn mul(&self, other: &Self) -> Self {
        GuardedRat::new(
            ck(self.num.checked_mul(other.num)),
            ck(self.den.checked_mul(other.den)),
        )
    }
    fn neg(&self) -> Self {
        GuardedRat {
            num: self.num.wrapping_neg(),
            den: self.den,
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
        // Cancel against the numerator before growing the denominator: the
        // divisor is z_μ, which reaches |μ|!.
        let g = gcd(self.num, n).max(1);
        GuardedRat::new(self.num / g, ck(self.den.checked_mul(n / g)))
    }
}

impl Plethystic for GuardedRat {
    /// ℚ has no variables to raise.
    fn frobenius(&self, _n: u32) -> Self {
        *self
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::convert::FromSchur;
    use crate::partition::Partition;
    use crate::sym::{PowerSum, Schur, SymFn};
    use std::sync::Mutex;

    /// The counter is global and monotone, so a test that overflows on purpose
    /// makes any *concurrently* running `guarded` scope report `None` too. That
    /// is the intended fail-safe direction — see the module docs — but it means
    /// these tests must not run in parallel with each other. They are the only
    /// place in the crate that constructs `Guarded`, so serialising them here is
    /// enough; poisoning is ignored so one failure does not cascade.
    static SERIAL: Mutex<()> = Mutex::new(());

    fn serial() -> std::sync::MutexGuard<'static, ()> {
        SERIAL.lock().unwrap_or_else(|e| e.into_inner())
    }

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

    /// A real conversion over the guarded ring must agree with the unguarded
    /// one wherever it succeeds — the point being that guarding changes only
    /// the reporting, never a value.
    #[test]
    fn guarded_conversions_agree_with_the_plain_ring() {
        let _g = serial();
        for n in 1..=7u32 {
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
        let ones = Partition::new(std::iter::repeat(1).take(n as usize));
        let got = guarded(|| crate::character::character_in::<Guarded>(&lam, &ones));
        assert_eq!(got, None, "χ^λ(1^n) for |λ| = {n} exceeds i128");
    }
}
