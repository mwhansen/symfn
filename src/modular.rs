//! Prime-field arithmetic, CRT, rational reconstruction, interpolation.
//!
//! The toolkit for **evaluate–interpolate–lift**: run a computation at numeric
//! points over several small primes, interpolate the answer's coefficients, and
//! reconstruct the exact rationals. It exists because it was measured to be the
//! only surviving way to make `gjmod.rs` fast (`docs/record/jack.md`), but
//! nothing in it knows about Jack polynomials, so it lives here rather than
//! there.
//!
//! ## Where else this applies
//!
//! The crate has several other places that carry a rational function through a
//! long computation and only need the *answer* to be exact:
//!
//! - [`macop`](crate::macop) builds Macdonald operator matrices over ℚ(q,t) and
//!   extracts eigenvectors — the eigenvector is the answer, the matrix is not.
//! - [`deltaop`](crate::deltaop) and [`frac`](crate::frac) divide by factored
//!   atoms exactly as `afrac` does; the same "the intermediate is a fraction
//!   but the answer is a polynomial" shape holds.
//! - [`qtkostka`](crate::qtkostka) is polynomial in q,t by theorem, which is
//!   the same license `gjmod` uses for b.
//!
//! ⚠️ **None of those have been converted, and none has been measured.** The
//! shape matching is not evidence that it would pay: `gjmod` won only after
//! two rounds of sampling corrected three wrong guesses about where its time
//! went, and the win came from the *representation* (31-bit primes, one
//! precomputed reconstruction matrix), not from modularity as such
//! (`docs/record/jack.md`).
//!
//! ## What makes it usable rather than a guess
//!
//! Rational reconstruction does **not** report failure when the true value is
//! bigger than the bound — it returns some spurious small rational instead. So
//! any caller has to answer "was the bound big enough?" independently. The
//! pattern that works, and the one `gjmod` uses, is to hold one prime back from
//! the CRT and require the reconstructed answer to reduce to what that prime
//! computed. See [`Md`] on why the primes are small.

// Every value in this module is a residue mod `p < 2^31`, or a product of two
// such in `u128`, or a Barrett reciprocal `2^62 / p`. Each cast lands inside
// the target by that bound, stated again at the sites that are not obvious.
// This is the module whose ring *is* modular (R4), so wrapping would be
// exactness here — but none of these casts wrap.
#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_sign_loss,
    clippy::cast_possible_wrap
)]

// ------------------------------------------------------------- modular ------

/// Arithmetic in `ℤ/p` for a prime `p < 2^31`, so a product fits `u64`.
///
/// ⚠️ **The prime size is a performance decision, not a taste one.** A 61-bit
/// prime makes a product a `u128` — and `u128 %`
/// is not an instruction on aarch64, it is a call into
/// `compiler_builtins::u128_div`, and sampling put nearly the whole engine
/// inside `__umodti3`. Dropping under `2^31` puts the product in a `u64`,
/// where the remainder is one hardware `udiv`, and Barrett reduction removes
/// even that (`docs/record/jack.md`).
///
/// ⚠️ A 128-bit divide instruction would *not* change this. x86-64's `DIV r64`
/// is 128÷64→64 and Rust cannot emit it for `u128 % u128` — it cannot prove the
/// divisor fits — and at ~30–90 cycles it would still lose to two multiplies.
/// What such an instruction *would* change is the prime size: 61-bit primes
/// would then cost the same per multiply and need one fewer CRT prime, hence
/// one fewer evaluation pass.
///
/// The cost of small primes is range: one 31-bit prime lifts nothing, so
/// answers are CRT'd across several. That is arithmetic needed anyway to get
/// past the `2^30` reconstruction bound a single 61-bit prime gives.
#[derive(Clone, Copy, Debug)]
pub struct Md {
    /// The modulus itself.
    pub p: u64,
    /// `⌊2^62 / p⌋`, the Barrett constant.
    r: u64,
}

impl Md {
    /// Build a modulus.
    ///
    /// # Panics
    ///
    /// Panics unless `p < 2^31`, which is what makes a product of two
    /// residues fit a `u64` with room for the Barrett step. Panics unless `p`
    /// is prime: [`inv`](Md::inv) is Fermat inversion, which for a composite
    /// modulus returns a non-inverse with no signal, and [`crt`] would build a
    /// wrong residue out of it — checked here, where it costs one
    /// Miller–Rabin per modulus, rather than trusted to the type's name.
    pub fn new(p: u64) -> Self {
        assert!(p < (1 << 31), "the modulus must fit 31 bits");
        assert!(
            is_prime(p),
            "the modulus must be prime: Fermat inversion is silently wrong mod {p}"
        );
        Md {
            p,
            r: ((1u128 << 62) / p as u128) as u64,
        }
    }
    /// `x mod p` for `x < 2^62`, with no division.
    #[inline]
    pub fn red(&self, x: u64) -> u64 {
        debug_assert!(x < 1 << 62);
        let q = ((x as u128 * self.r as u128) >> 62) as u64;
        let mut v = x - q * self.p;
        while v >= self.p {
            v -= self.p;
        }
        v
    }
    #[inline]
    pub fn mul(&self, a: u64, b: u64) -> u64 {
        self.red(a * b)
    }
    #[inline]
    pub fn add(&self, a: u64, b: u64) -> u64 {
        let s = a + b;
        if s >= self.p {
            s - self.p
        } else {
            s
        }
    }
    #[inline]
    pub fn sub(&self, a: u64, b: u64) -> u64 {
        if a >= b {
            a - b
        } else {
            a + self.p - b
        }
    }
    pub fn pow(&self, mut a: u64, mut e: u64) -> u64 {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = self.mul(r, a);
            }
            a = self.mul(a, a);
            e >>= 1;
        }
        r
    }
    /// `a⁻¹` by Fermat.
    ///
    /// # Panics
    ///
    /// Panics on zero — for a caller evaluating a rational function that means
    /// the point hit a pole, which is a bug in the choice of points rather than
    /// an arithmetic failure.
    pub fn inv(&self, a: u64) -> u64 {
        assert!(a != 0, "inverse of zero mod {}", self.p);
        self.pow(a, self.p - 2)
    }
    // `wrong_self_convention` reads `from_*` as a constructor whose `self` is
    // the thing converted. Here `self` is the *modulus* and the argument is
    // the value, and the name is deliberately the free functions' name below,
    // so a call site reads the same whether the modulus is a machine word or
    // this struct. Renaming would break that pairing to satisfy a heuristic
    // about the other operand.
    #[allow(clippy::wrong_self_convention)]
    #[inline]
    pub fn from_i128(&self, v: i128) -> u64 {
        v.rem_euclid(self.p as i128) as u64
    }
    #[allow(clippy::wrong_self_convention)]
    #[inline]
    pub fn from_u128(&self, v: u128) -> u64 {
        (v % self.p as u128) as u64
    }
}

// Free-function forms, so a call site reads the same whether the modulus is a
// machine word or a struct.
#[inline]
pub fn mul(a: u64, b: u64, p: Md) -> u64 {
    p.mul(a, b)
}
#[inline]
pub fn add(a: u64, b: u64, p: Md) -> u64 {
    p.add(a, b)
}
#[inline]
pub fn sub(a: u64, b: u64, p: Md) -> u64 {
    p.sub(a, b)
}
#[inline]
pub fn pow(a: u64, e: u64, p: Md) -> u64 {
    p.pow(a, e)
}
#[inline]
pub fn inv(a: u64, p: Md) -> u64 {
    p.inv(a)
}
#[inline]
pub fn from_i128(v: i128, p: Md) -> u64 {
    p.from_i128(v)
}

/// Deterministic Miller–Rabin. The listed bases are a complete witness set for
/// every `u64`, so this decides primality rather than guessing at it — which is
/// why callers can *find* their primes rather than recalling them.
pub fn is_prime(n: u64) -> bool {
    if n < 2 {
        return false;
    }
    for q in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if n.is_multiple_of(q) {
            return n == q;
        }
    }
    let mut d = n - 1;
    let mut s = 0;
    while d.is_multiple_of(2) {
        d /= 2;
        s += 1;
    }
    let mulmod = |a: u64, b: u64| ((a as u128 * b as u128) % n as u128) as u64;
    let powmod = |mut a: u64, mut e: u64| {
        let mut r = 1u64;
        while e > 0 {
            if e & 1 == 1 {
                r = mulmod(r, a);
            }
            a = mulmod(a, a);
            e >>= 1;
        }
        r
    };
    'witness: for a in [2u64, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        let mut x = powmod(a, d);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = mulmod(x, x);
            if x == n - 1 {
                continue 'witness;
            }
        }
        return false;
    }
    true
}

/// The largest prime strictly below `hi`.
pub fn prime_below(hi: u64) -> u64 {
    let mut n = if hi.is_multiple_of(2) { hi - 1 } else { hi - 2 };
    while !is_prime(n) {
        n -= 2;
    }
    n
}

/// A deterministic sequence of distinct 31-bit primes.
///
/// Windows spaced `2^20` apart, each descending to the largest prime below it,
/// so the primes do not depend on how many are asked for — a caller that grows
/// its prime count keeps the ones it already has. `k = 0` lands on `2^31 − 1`.
pub fn nth_prime(k: usize) -> Md {
    Md::new(prime_below((1u64 << 31) - (k as u64) * (1 << 20)))
}

// --------------------------------------------------------- reconstruction ---

/// Recover `num/den ≡ x (mod m)` with `|num|, den ≤ bound`, or `None`.
///
/// The half-extended Euclid on `(m, x)`, stopped when the remainder drops below
/// the bound. Unique when `2·bound² < m`, which is what [`bound_for`] gives —
/// so a successful reconstruction is the *only* small rational congruent to
/// `x`, not merely one of them.
///
/// ⚠️ Uniqueness is not self-detection. If the true value is **larger** than
/// the bound, this can still return some spurious small rational rather than
/// failing. Callers must check the answer against an independently computed
/// prime; see the module docs.
pub fn reconstruct(x: u128, m: u128, bound: u128) -> Option<(i128, u128)> {
    debug_assert!(m < (1u128 << 127), "the modulus must fit a signed i128");
    let (mut r0, mut r1) = (m as i128, x as i128);
    let (mut t0, mut t1) = (0i128, 1i128);
    while r1 > bound as i128 {
        if r1 == 0 {
            return None;
        }
        let q = r0 / r1;
        let r = r0 - q * r1;
        r0 = r1;
        r1 = r;
        let t = t0 - q * t1;
        t0 = t1;
        t1 = t;
    }
    let (mut num, mut den) = (r1, t1);
    if den < 0 {
        num = -num;
        den = -den;
    }
    if den == 0 || den as u128 > bound {
        return None;
    }
    // No congruence check here, deliberately. The extended Euclid maintains
    // `num ≡ x·den (mod m)` by construction, so re-checking it would verify the
    // arithmetic rather than the mathematics — and doing so needs `x·den`, which
    // is 182 bits and does not fit `u128`. The first version accumulated it with
    // a 128-iteration doubling loop.
    //
    // ⚠️ That loop was previously recorded as costing **2.5×** of the engine that
    // uses this. It did not: removing it moved the ladder 7.46s → 6.85s at
    // n = 10, inside the noise. The claim was a guess, and the correction is
    // left in place rather than deleted. The reason to drop the check is the one
    // above — it checks the wrong thing.
    Some((num, den as u128))
}

/// The unique residue mod `∏ pᵢ` agreeing with `aᵢ` mod each `pᵢ`, and that
/// product.
///
/// The primes must be distinct; [`nth_prime`] supplies such a sequence.
///
/// Incremental Garner: carry `x` and `M = ∏` so far, and correct by
/// `M·((aᵢ − x)·M⁻¹ mod pᵢ)`.
///
/// This is what buys the range. A single 31-bit prime lifts nothing at all —
/// the reconstruction bound would be `2^15` — and even a 61-bit one gives only
/// `2^30`. Three 31-bit primes give `2^46`.
///
/// # Panics
///
/// Panics if `residues` and `primes` differ in length — `zip` would silently
/// drop the excess and answer for a *smaller* reconstruction than the caller
/// asked, with a modulus their [`bound_for`] does not expect. Panics if both
/// are empty. Panics through [`inv`] if a prime repeats, because the product
/// so far is then zero mod that prime.
pub fn crt(residues: &[u64], primes: &[Md]) -> (u128, u128) {
    assert!(
        residues.len() == primes.len(),
        "crt needs one residue per prime: {} residues against {} primes",
        residues.len(),
        primes.len()
    );
    assert!(!residues.is_empty(), "crt needs at least one residue");
    let mut x = residues[0] as u128;
    let mut m = primes[0].p as u128;
    for (&a, &q) in residues.iter().zip(primes.iter()).skip(1) {
        let m_mod = q.from_u128(m);
        let diff = q.sub(a % q.p, q.from_u128(x));
        let k = q.mul(diff, q.inv(m_mod));
        x += m * k as u128;
        m *= q.p as u128;
    }
    (x, m)
}

/// The largest reconstruction bound a modulus supports: `⌊√(m/2)⌋`.
pub fn bound_for(modulus: u128) -> u128 {
    isqrt128(modulus / 2)
}

/// Integer square root, by Newton.
pub fn isqrt128(n: u128) -> u128 {
    if n == 0 {
        return 0;
    }
    let mut r = 1u128 << ((128 - n.leading_zeros()).div_ceil(2));
    loop {
        let next = (r + n / r) / 2;
        if next >= r {
            break;
        }
        r = next;
    }
    r
}

pub fn gcd128(mut a: u128, mut b: u128) -> u128 {
    while b != 0 {
        let t = a % b;
        a = b;
        b = t;
    }
    a
}

// --------------------------------------------------------- interpolation ----

/// Lagrange interpolation through `(x_i, y_i)`, returning dense coefficients.
///
/// ⚠️ For many keys sampled at the *same* points, do not call this per key —
/// build [`lagrange_matrix`] once instead. Per key it rebuilds the same basis
/// polynomials and runs a Fermat inversion once per output entry, which
/// dominates any engine calling it that way (`docs/record/jack.md`).
///
/// # Panics
///
/// Panics unless the `xs` are distinct mod `p`, since a repeated point makes
/// the denominator zero and [`inv`] rejects that. Panics if `ys` is shorter
/// than `xs`.
// Dead outside the tests, and deliberately: no engine calls it, because the
// warning above is why [`lagrange_matrix`] exists. It is kept as
// the obviously-correct form that composition is checked against
// (`docs/policies/validation.md`).
#[allow(dead_code)]
pub fn interpolate(xs: &[u64], ys: &[u64], p: Md) -> Vec<u64> {
    let n = xs.len();
    let mut out = vec![0u64; n];
    for i in 0..n {
        // The basis polynomial ∏_{j≠i} (X − x_j) / (x_i − x_j).
        let mut denom = 1u64;
        for j in 0..n {
            if j != i {
                denom = mul(denom, sub(xs[i], xs[j], p), p);
            }
        }
        let scale = mul(ys[i], inv(denom, p), p);
        if scale == 0 {
            continue;
        }
        let basis = basis_poly(xs, i, p);
        for d in 0..n {
            out[d] = add(out[d], mul(basis[d], scale, p), p);
        }
    }
    out
}

/// `∏_{j≠i}(X − x_j)`, dense.
fn basis_poly(xs: &[u64], i: usize, p: Md) -> Vec<u64> {
    let n = xs.len();
    let mut basis = vec![0u64; n];
    basis[0] = 1;
    let mut deg = 0usize;
    for j in 0..n {
        if j == i {
            continue;
        }
        for d in (0..=deg).rev() {
            let v = basis[d];
            basis[d + 1] = add(basis[d + 1], v, p);
            basis[d] = mul(v, sub(0, xs[j], p), p);
        }
        deg += 1;
    }
    basis
}

/// `f(X) ↦ f(shift + X)`, dense — the binomial transform.
// The other half of what [`lagrange_matrix`] folds into one matrix, kept for
// the same reason as [`interpolate`]: the tests compose these two and demand
// the matrix agree.
#[allow(dead_code)]
pub fn shift_by(coeffs: &[u64], shift: u64, p: Md) -> Vec<u64> {
    let n = coeffs.len();
    let mut out = vec![0u64; n];
    for (k, &c) in coeffs.iter().enumerate() {
        if c == 0 {
            continue;
        }
        // Σ_j C(k,j)·shift^{k−j}·X^j
        let mut binom = 1u64;
        for j in 0..=k {
            let w = mul(binom, pow(shift, (k - j) as u64, p), p);
            out[j] = add(out[j], mul(c, w, p), p);
            // C(k, j+1) = C(k, j)·(k−j)/(j+1)
            binom = mul(
                mul(binom, (k - j) as u64 % p.p, p),
                inv((j + 1) as u64 % p.p, p),
                p,
            );
        }
    }
    out
}

/// The composite map: values at `xs` ↦ coefficients of `f(shift + X)`.
///
/// Interpolation and the shift are **both linear, and both identical for every
/// key**, so composing them once turns per-key reconstruction from `N²`
/// multiplies *plus `N` Fermat inversions* into `N²` multiplies and nothing
/// else. `T[j][i]` is the weight of the value at `xs[i]` in output coefficient
/// `j`; apply it with [`apply_matrix`].
///
/// Interpolating per key puts the whole engine inside [`interpolate`] and the
/// shift, with the actual pipeline absent from the profile
/// (`docs/record/jack.md`).
///
/// # Panics
///
/// Panics unless the `xs` are distinct mod `p`, since a repeated point makes
/// the denominator zero and [`inv`] rejects that.
pub fn lagrange_matrix(xs: &[u64], shift: u64, p: Md) -> Vec<Vec<u64>> {
    let n = xs.len();
    // A[d][i]: the coefficient of X^d in the i-th Lagrange basis polynomial.
    let mut a = vec![vec![0u64; n]; n];
    for i in 0..n {
        let mut denom = 1u64;
        for j in 0..n {
            if j != i {
                denom = mul(denom, sub(xs[i], xs[j], p), p);
            }
        }
        let scale = inv(denom, p);
        let basis = basis_poly(xs, i, p);
        for d in 0..n {
            a[d][i] = mul(basis[d], scale, p);
        }
    }
    // Pascal's triangle, so the binomial transform needs no division either.
    let mut binom = vec![vec![0u64; n]; n];
    for d in 0..n {
        binom[d][0] = 1;
        for j in 1..=d {
            binom[d][j] = add(binom[d - 1][j - 1], binom[d - 1][j], p);
        }
    }
    // T[j][i] = Σ_d C(d, j)·shift^{d−j}·A[d][i]
    let mut t = vec![vec![0u64; n]; n];
    for j in 0..n {
        for d in j..n {
            let c = mul(binom[d][j], pow(shift, (d - j) as u64, p), p);
            if c == 0 {
                continue;
            }
            for i in 0..n {
                t[j][i] = add(t[j][i], mul(c, a[d][i], p), p);
            }
        }
    }
    t
}

/// Apply a [`lagrange_matrix`] to one key's values.
pub fn apply_matrix(t: &[Vec<u64>], ys: &[u64], p: Md) -> Vec<u64> {
    t.iter()
        .map(|row| {
            let mut acc = 0u64;
            for (r, y) in row.iter().zip(ys.iter()) {
                if *y != 0 {
                    acc = add(acc, mul(*r, *y, p), p);
                }
            }
            acc
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_primes_are_actually_prime() {
        for k in 0..8 {
            let q = nth_prime(k).p;
            assert!(is_prime(q), "{q} is not prime");
            assert!(q < (1 << 31), "{q} does not fit 31 bits");
        }
        // Distinct, so a CRT across them is meaningful.
        let all: Vec<u64> = (0..8).map(|k| nth_prime(k).p).collect();
        for i in 0..all.len() {
            for j in i + 1..all.len() {
                assert_ne!(all[i], all[j]);
            }
        }
        assert!(!is_prime(3215031751), "a strong pseudoprime to 2,3,5,7");
        assert!(!is_prime(1u64 << 40));
        assert!(is_prime(2));
    }

    /// Barrett reduction must agree with the division it replaces, including at
    /// the top of its declared range.
    #[test]
    fn barrett_agrees_with_division() {
        for k in 0..4 {
            let p = nth_prime(k);
            for x in [0u64, 1, p.p - 1, p.p, p.p + 1, (1 << 62) - 1] {
                assert_eq!(p.red(x), x % p.p, "red({x}) mod {}", p.p);
            }
            // Every product of two residues, at the extremes.
            for &(a, b) in &[(0u64, 0u64), (1, 1), (p.p - 1, p.p - 1), (p.p - 1, 2)] {
                assert_eq!(p.mul(a, b), (a as u128 * b as u128 % p.p as u128) as u64);
            }
        }
    }

    #[test]
    fn rational_reconstruction_round_trips() {
        let primes: Vec<Md> = (0..3).map(nth_prime).collect();
        let m: u128 = primes.iter().map(|q| q.p as u128).product();
        let bound = bound_for(m);
        assert!(bound > (1u128 << 45), "three primes must buy a 2^46 bound");
        for &(num, den) in &[
            (1i128, 1u128),
            (-7, 1),
            (3, 4),
            (123456, 7),
            (-99991, 65537),
            // Past 2^30: exactly the range a single 61-bit prime could not
            // reach, and the coefficients pass 2^30 at n = 11.
            (1234567890123, 1),
            (-12345678901234, 1000003),
        ] {
            assert!(num.unsigned_abs() <= bound && den <= bound);
            let residues: Vec<u64> = primes
                .iter()
                .map(|q| q.mul(q.from_i128(num), q.inv(q.from_u128(den))))
                .collect();
            let (x, mm) = crt(&residues, &primes);
            assert_eq!(mm, m);
            assert_eq!(reconstruct(x, m, bound), Some((num, den)), "{num}/{den}");
        }

        // And the other side of the bound, which is why a caller
        // must hold a prime back. `-98765432109876/1000003` has a numerator of
        // 9.88e13 against a bound of 7.03e13, so it is out of range: the
        // reconstruction does **not** recover it, and does not report failure
        // either. Only an independent prime can tell.
        let (num, den) = (-98765432109876i128, 1000003u128);
        assert!(num.unsigned_abs() > bound);
        let residues: Vec<u64> = primes
            .iter()
            .map(|q| q.mul(q.from_i128(num), q.inv(q.from_u128(den))))
            .collect();
        let (x, _) = crt(&residues, &primes);
        assert_ne!(reconstruct(x, m, bound), Some((num, den)));
    }

    /// A composite modulus makes [`Md::inv`] silently wrong, so construction
    /// is where it must die.
    #[test]
    #[should_panic(expected = "the modulus must be prime")]
    fn a_composite_modulus_is_refused_at_construction() {
        Md::new(15);
    }

    /// A short residue list once truncated silently — `zip` answered for the
    /// primes it could pair and dropped the rest.
    #[test]
    #[should_panic(expected = "one residue per prime")]
    fn crt_refuses_a_residue_list_shorter_than_the_primes() {
        let primes: Vec<Md> = (0..3).map(nth_prime).collect();
        crt(&[1, 2], &primes);
    }

    #[test]
    fn crt_combines_residues() {
        let primes: Vec<Md> = (0..3).map(nth_prime).collect();
        for v in [0u128, 1, 42, 1u128 << 90, 123456789012345678901234567890] {
            let residues: Vec<u64> = primes.iter().map(|q| q.from_u128(v)).collect();
            let (x, m) = crt(&residues, &primes);
            for q in &primes {
                assert_eq!(x % q.p as u128, v % q.p as u128, "mod {}", q.p);
            }
            assert_eq!(x % m, v % m);
        }
    }

    #[test]
    fn interpolation_recovers_a_known_polynomial() {
        let p = nth_prime(0);
        // 5 − 3X + 2X³, sampled at X = 1..6
        let f = |x: u64| -> u64 {
            let x = x as i128;
            p.from_i128(5 - 3 * x + 2 * x * x * x)
        };
        let xs: Vec<u64> = (1..=6).collect();
        let ys: Vec<u64> = xs.iter().map(|&x| f(x)).collect();
        let got = interpolate(&xs, &ys, p);
        let want = [5i128, -3, 0, 2, 0, 0].map(|c| p.from_i128(c));
        assert_eq!(got, want.to_vec());
    }

    /// The matrix must be exactly the composition it replaces, shift and all.
    #[test]
    fn the_matrix_composes_interpolation_with_the_shift() {
        let p = nth_prime(1);
        let xs: Vec<u64> = (1..=7).collect();
        for shift in [0u64, 1, 5] {
            let t = lagrange_matrix(&xs, shift, p);
            for seed in [1i128, 2, 3] {
                let ys: Vec<u64> = xs
                    .iter()
                    .map(|&x| p.from_i128(seed * (x as i128) * (x as i128) + 4 * x as i128 - 1))
                    .collect();
                let direct = shift_by(&interpolate(&xs, &ys, p), shift, p);
                assert_eq!(apply_matrix(&t, &ys, p), direct, "shift {shift}");
            }
        }
    }

    /// `f(1+b)` is the specific shift the \[GJ\] engine needs, and `α = 1` must
    /// land on `b = 0`.
    #[test]
    fn the_shift_by_one_sends_alpha_to_one_plus_b() {
        let p = nth_prime(0);
        // f(α) = α² ↦ (1+b)² = 1 + 2b + b²
        let got = shift_by(&[0, 0, 1], 1, p);
        assert_eq!(got, vec![1, 2, 1]);
        assert_eq!(shift_by(&[7, 0, 0], 1, p), vec![7, 0, 0]);
    }
}
