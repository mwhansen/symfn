//! What does overflow-checked coefficient arithmetic cost?
//!
//! The question is which shape the Python boundary should take. Today it
//! instantiates `i128` / `Rational`, whose `Ring` impls use plain `*` — so past
//! the fixed width they wrap **silently**, which is not something a Symmetrica
//! replacement can ship. Three ways out: keep i128 and refuse loudly, always use
//! bignums, or compute in i128 with checked arithmetic and re-run the whole call
//! in a bignum ring if anything overflowed.
//!
//! The third is the crate's existing pattern (`try_character` → `character_in`)
//! and keeps the fast path fast — but only if "checked" is nearly free. That is
//! the number this measures, and nothing else: `G` and `GRat` below are `i128`
//! and `Rational` with every arithmetic op replaced by its `checked_` form,
//! recording overflow in a flag rather than wrapping.
//!
//! Rounds are interleaved **and the order is rotated**, so each variant spends
//! an equal share of rounds in each position. That is not fussiness: with a
//! fixed order this benchmark reported the checked rationals as 36% *faster*
//! than the unchecked ones, and swapping the two blocks moved the 36% to the
//! other type. Whichever rational pass runs first after the much larger
//! integral passes pays ~60% for arriving with a cold cache, and a fixed order
//! silently charges that to one type. The project has now made a version of
//! this mistake three times; rotating is the fix.
//!
//!   cargo run --release --example bench_guarded

use std::sync::atomic::{AtomicBool, Ordering};
use std::time::Instant;

use symfn::convert::{FromSchur, ToSchur};
use symfn::{
    Elementary, Homogeneous, Monomial, Partition, Plethystic, PowerSum, QAlgebra, Rational, Ring,
    Schur, SymFn,
};

static OVERFLOW: AtomicBool = AtomicBool::new(false);

#[inline(always)]
fn flag() {
    OVERFLOW.store(true, Ordering::Relaxed);
}

// --- guarded i128 ------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Debug)]
struct G(i128);

impl Ring for G {
    #[inline]
    fn zero() -> Self {
        G(0)
    }
    #[inline]
    fn one() -> Self {
        G(1)
    }
    #[inline]
    fn is_zero(&self) -> bool {
        self.0 == 0
    }
    #[inline]
    fn add_assign(&mut self, other: &Self) {
        match self.0.checked_add(other.0) {
            Some(v) => self.0 = v,
            None => flag(),
        }
    }
    #[inline]
    fn mul(&self, other: &Self) -> Self {
        match self.0.checked_mul(other.0) {
            Some(v) => G(v),
            None => {
                flag();
                G(0)
            }
        }
    }
    #[inline]
    fn neg(&self) -> Self {
        G(-self.0)
    }
    #[inline]
    fn from_i64(n: i64) -> Self {
        G(n as i128)
    }
    #[inline]
    fn from_u128(n: u128) -> Self {
        if n > i128::MAX as u128 {
            flag();
            return G(0);
        }
        G(n as i128)
    }
    #[inline]
    fn from_i128(n: i128) -> Self {
        G(n)
    }
}

// --- guarded Rational --------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Debug)]
struct GRat {
    num: i128,
    den: i128,
}

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

impl GRat {
    fn new(num: i128, den: i128) -> Self {
        if den == 0 {
            flag();
            return GRat { num: 0, den: 1 };
        }
        let g = gcd(num, den).max(1);
        let (mut n, mut d) = (num / g, den / g);
        if d < 0 {
            n = -n;
            d = -d;
        }
        GRat { num: n, den: d }
    }
    #[inline]
    fn ck(v: Option<i128>) -> i128 {
        match v {
            Some(x) => x,
            None => {
                flag();
                0
            }
        }
    }
}

impl Ring for GRat {
    fn zero() -> Self {
        GRat { num: 0, den: 1 }
    }
    fn one() -> Self {
        GRat { num: 1, den: 1 }
    }
    fn is_zero(&self) -> bool {
        self.num == 0
    }
    fn add_assign(&mut self, other: &Self) {
        let a = Self::ck(self.num.checked_mul(other.den));
        let b = Self::ck(other.num.checked_mul(self.den));
        let num = Self::ck(a.checked_add(b));
        let den = Self::ck(self.den.checked_mul(other.den));
        *self = GRat::new(num, den);
    }
    fn mul(&self, other: &Self) -> Self {
        GRat::new(
            Self::ck(self.num.checked_mul(other.num)),
            Self::ck(self.den.checked_mul(other.den)),
        )
    }
    fn neg(&self) -> Self {
        GRat {
            num: -self.num,
            den: self.den,
        }
    }
    fn from_i64(n: i64) -> Self {
        GRat {
            num: n as i128,
            den: 1,
        }
    }
    fn from_u128(n: u128) -> Self {
        if n > i128::MAX as u128 {
            flag();
            return GRat::zero();
        }
        GRat {
            num: n as i128,
            den: 1,
        }
    }
    fn from_i128(n: i128) -> Self {
        GRat { num: n, den: 1 }
    }
    // Same opt-in as `Rational`, so the integral p → s sweep is available to
    // both sides and the comparison is like for like.
    fn as_ratio(&self) -> Option<(i128, i128)> {
        Some((self.num, self.den))
    }
    fn from_ratio(num: i128, den: i128) -> Option<Self> {
        Some(GRat::new(num, den))
    }
}

impl QAlgebra for GRat {
    fn div_u128(&self, n: u128) -> Self {
        let n = n as i128;
        let g = gcd(self.num, n).max(1);
        GRat::new(self.num / g, Self::ck(self.den.checked_mul(n / g)))
    }
}

impl Plethystic for GRat {
    fn frobenius(&self, _n: u32) -> Self {
        *self
    }
}

// --- unguarded twin of GRat --------------------------------------------------
//
// Byte-for-byte `GRat` with the `checked_` calls replaced by plain operators.
// Comparing `GRat` against the library's own `Rational` conflates the checking
// with every other difference between two independently written types; this
// isolates it, since the two are compiled in the same crate under the same
// conditions and differ in nothing else.

#[derive(Clone, Copy, PartialEq, Debug)]
struct PRat {
    num: i128,
    den: i128,
}

impl PRat {
    fn new(num: i128, den: i128) -> Self {
        let g = gcd(num, den).max(1);
        let (mut n, mut d) = (num / g, den / g);
        if d < 0 {
            n = -n;
            d = -d;
        }
        PRat { num: n, den: d }
    }
}

impl Ring for PRat {
    fn zero() -> Self {
        PRat { num: 0, den: 1 }
    }
    fn one() -> Self {
        PRat { num: 1, den: 1 }
    }
    fn is_zero(&self) -> bool {
        self.num == 0
    }
    fn add_assign(&mut self, other: &Self) {
        let num = self.num * other.den + other.num * self.den;
        let den = self.den * other.den;
        *self = PRat::new(num, den);
    }
    fn mul(&self, other: &Self) -> Self {
        PRat::new(self.num * other.num, self.den * other.den)
    }
    fn neg(&self) -> Self {
        PRat {
            num: -self.num,
            den: self.den,
        }
    }
    fn from_i64(n: i64) -> Self {
        PRat {
            num: n as i128,
            den: 1,
        }
    }
    fn from_u128(n: u128) -> Self {
        PRat {
            num: n as i128,
            den: 1,
        }
    }
    fn from_i128(n: i128) -> Self {
        PRat { num: n, den: 1 }
    }
    fn as_ratio(&self) -> Option<(i128, i128)> {
        Some((self.num, self.den))
    }
    fn from_ratio(num: i128, den: i128) -> Option<Self> {
        Some(PRat::new(num, den))
    }
}

impl QAlgebra for PRat {
    fn div_u128(&self, n: u128) -> Self {
        let n = n as i128;
        let g = gcd(self.num, n).max(1);
        PRat::new(self.num / g, self.den * (n / g))
    }
}

impl Plethystic for PRat {
    fn frobenius(&self, _n: u32) -> Self {
        *self
    }
}

// --- workloads ---------------------------------------------------------------

fn part(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

fn staircase(n: u32) -> Vec<u32> {
    let mut v: Vec<u32> = Vec::new();
    let mut left = n;
    let mut k = (n as f64).sqrt() as u32 + 2;
    while left > 0 {
        let take = k.min(left);
        v.push(take);
        left -= take;
        if k > 1 {
            k -= 1;
        }
    }
    v.sort_unstable_by(|a, b| b.cmp(a));
    v
}

/// One timed pass over the integral (`Ring`-only) workloads.
///
/// Caches are cleared first, so both sides pay the full cost rather than one of
/// them reading the other's memo table. That is also what a cold call from Sage
/// looks like, and what `scripts/compare_sage.py` does for the same reason.
fn integral_pass<C: Ring>() -> usize {
    symfn::clear_caches();
    let mut n = 0;
    // Schur products: the LR backend returns u128 structure constants which are
    // injected and multiplied, so this is the product path's coefficient cost.
    for parts in [
        &[8u32, 7, 6, 5, 4, 3][..],
        &[10, 8, 6, 4],
        &[7, 6, 5, 4, 3, 2, 1],
    ] {
        let a: Schur<C> = Schur::monomial(part(parts), C::one());
        n += a.mul(&a).terms().len();
    }
    // Conversions that stay over ℤ.
    for deg in [18u32, 20, 22] {
        let lam = part(&staircase(deg));
        let s: Schur<C> = Schur::monomial(lam.clone(), C::one());
        n += Monomial::<C>::from_schur(&s).terms().len();
        n += Homogeneous::<C>::from_schur(&s).terms().len();
        n += Elementary::<C>::from_schur(&s).terms().len();
        let m: Monomial<C> = Monomial::monomial(lam, C::one());
        n += m.to_schur().terms().len();
    }
    n
}

/// One timed pass over the dividing workloads.
fn rational_pass<C: Plethystic>() -> usize {
    symfn::clear_caches();
    let mut n = 0;
    for deg in [14u32, 16, 18] {
        let s: Schur<C> = Schur::monomial(part(&staircase(deg)), C::one());
        let p: PowerSum<C> = PowerSum::from_schur(&s);
        n += p.terms().len();
        n += p.to_schur().terms().len();
    }
    for (f, g) in [
        (&[4u32][..], &[2u32, 1][..]),
        (&[3, 1], &[2, 1]),
        (&[2, 2], &[3]),
    ] {
        let sf: Schur<C> = Schur::monomial(part(f), C::one());
        let sg: Schur<C> = Schur::monomial(part(g), C::one());
        n += symfn::plethysm(&sf, &sg).terms().len();
    }
    for parts in [&[6u32, 5, 4][..], &[7, 5, 3, 1]] {
        let a: Schur<C> = Schur::monomial(part(parts), C::one());
        n += symfn::internal(&a, &a).terms().len();
    }
    n
}

fn main() {
    const ROUNDS: usize = 9;

    // Warm-up: fill the memo caches so neither side pays for building them, and
    // let the CPU settle before the first timed round.
    let _ = integral_pass::<i128>();
    let _ = rational_pass::<Rational>();
    let _ = integral_pass::<G>();
    let _ = rational_pass::<GRat>();

    let _ = rational_pass::<PRat>();
    let (mut plain_i, mut guard_i) = (0.0f64, 0.0f64);
    let (mut plain_r, mut guard_r) = (0.0f64, 0.0f64);
    let mut lib_r = 0.0f64;
    let (mut ci, mut cg, mut cr, mut cgr, mut cpr) = (0, 0, 0, 0, 0);

    for round in 0..ROUNDS {
        // Alternate which of the pair goes first. The position immediately
        // after another pass is measurably penalised -- see the note in main's
        // header -- so a fixed order attributes that penalty to whichever type
        // happens to sit there.
        for k in 0..2 {
            if (round + k) % 2 == 0 {
                let t = Instant::now();
                ci = integral_pass::<i128>();
                plain_i += t.elapsed().as_secs_f64();
            } else {
                let t = Instant::now();
                cg = integral_pass::<G>();
                guard_i += t.elapsed().as_secs_f64();
            }
        }
        for k in 0..3 {
            match (round + k) % 3 {
                0 => {
                    let t = Instant::now();
                    cpr = rational_pass::<PRat>();
                    plain_r += t.elapsed().as_secs_f64();
                }
                1 => {
                    let t = Instant::now();
                    cgr = rational_pass::<GRat>();
                    guard_r += t.elapsed().as_secs_f64();
                }
                _ => {
                    let t = Instant::now();
                    cr = rational_pass::<Rational>();
                    lib_r += t.elapsed().as_secs_f64();
                }
            }
        }
    }

    // Same amount of work on both sides, or the ratio means nothing.
    assert_eq!(ci, cg, "integral passes disagree");
    assert_eq!(cr, cgr, "rational passes disagree");
    assert_eq!(cpr, cgr, "plain-twin pass disagrees");
    assert!(
        !OVERFLOW.load(Ordering::Relaxed),
        "a guarded op overflowed; the comparison would not be like for like"
    );

    println!("{ROUNDS} interleaved rounds, on AC\n");
    println!(
        "{:<34} {:>10} {:>10} {:>9}",
        "workload", "plain", "checked", "cost"
    );
    println!(
        "{:<34} {:>9.4}s {:>9.4}s {:>8.1}%",
        "integral (products, s->m/h/e, m->s)",
        plain_i,
        guard_i,
        (guard_i / plain_i - 1.0) * 100.0
    );
    println!(
        "{:<34} {:>9.4}s {:>9.4}s {:>8.1}%",
        "dividing (s->p, plethysm, kron)",
        plain_r,
        guard_r,
        (guard_r / plain_r - 1.0) * 100.0
    );
    println!(
        "{:<34} {:>9.4}s {:>9.4}s {:>8.1}%",
        "  ^ vs library Rational",
        lib_r,
        guard_r,
        (guard_r / lib_r - 1.0) * 100.0
    );
    let (tp, tg) = (plain_i + plain_r, guard_i + guard_r);
    println!(
        "{:<34} {:>9.4}s {:>9.4}s {:>8.1}%",
        "combined",
        tp,
        tg,
        (tg / tp - 1.0) * 100.0
    );
}
