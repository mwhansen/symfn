//! How big do coefficients actually get, and how much does coefficient
//! arithmetic cost?
//!
//! These decide what the escalated ring needs to be. If values past `i128` are
//! a handful of limbs, then every bignum library runs *schoolbook* on them and
//! GMP's Karatsuba/Toom/FFT never engage — its advantage there is assembly
//! tuning, not algorithms. If coefficient arithmetic is also a minority of
//! runtime, then "track the extra digits correctly" is the requirement and
//! "fast bignum" is not.
//!
//! The answers below are what chose `num-bigint` over GMP for the escalation
//! ring, so this is the example to re-run if that decision is ever revisited.
//!
//!   cargo run --release --features bignum --example coeff_sizes

use std::time::Instant;

use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::Signed;
use symfn::convert::{FromSchur, ToSchur};
use symfn::{Monomial, Partition, PowerSum, Rational, Ring, Schur, SymFn};

fn part(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

/// Widest coefficient in an integer-basis element, in bits and decimal digits.
fn widest(terms: impl Iterator<Item = BigInt>) -> (u32, usize) {
    let mut best = BigInt::from(0);
    for c in terms {
        let a = c.clone().abs();
        if a > best {
            best = a;
        }
    }
    (best.bits() as u32, best.to_string().len())
}

fn report(label: &str, bits: u32, digits: usize) {
    let limbs = bits.div_ceil(64);
    let past = if bits > 127 { "  <-- past i128" } else { "" };
    println!("{label:<44} {bits:>5} bits {digits:>5} digits {limbs:>3} limbs{past}");
}

fn main() {
    println!("=== how wide do coefficients get? (i128 holds 127 bits / 38 digits)\n");

    // Littlewood-Richardson: structure constants of a big square.
    for parts in [
        &[8u32, 7, 6, 5, 4, 3][..],
        &[12, 10, 8, 6],
        &[10, 9, 8, 7, 6, 5],
        &[16, 13, 10, 7],
    ] {
        let a: Schur<BigInt> = Schur::monomial(part(parts), BigInt::from(1));
        let prod = a.mul(&a);
        let (b, d) = widest(prod.terms().values().cloned());
        report(&format!("LR  s_{parts:?}^2"), b, d);
    }

    // Kostka numbers: s -> m at high degree. These are counts of tableaux and
    // are the most plausible source of a genuinely large integer here.
    for deg in [20u32, 26, 32] {
        let lam = staircase(deg);
        let s: Schur<BigInt> = Schur::monomial(part(&lam), BigInt::from(1));
        let m: Monomial<BigInt> = Monomial::from_schur(&s);
        let (b, d) = widest(m.terms().values().cloned());
        report(&format!("Kostka row  s_{lam:?} -> m"), b, d);
    }

    // The single largest Kostka number of a degree: K_{lambda,1^n} = f^lambda.
    for deg in [30u32, 45, 60] {
        let lam = staircase(deg);
        match symfn::dimension(&part(&lam)) {
            Some(v) => report(
                &format!("f^lambda for {lam:?}"),
                128 - (v.leading_zeros()),
                v.to_string().len(),
            ),
            None => println!("f^lambda for {lam:?}: past u128"),
        }
    }

    // Characters: p -> s carries chi^lambda(mu).
    for deg in [24u32, 32, 40] {
        let mu = part(&vec![1u32; deg as usize]);
        let p: PowerSum<BigInt> = PowerSum::monomial(mu, BigInt::from(1));
        let s = p.to_schur();
        let (b, d) = widest(s.terms().values().cloned());
        report(&format!("characters  p_1^{deg} -> s"), b, d);
    }

    // Plethysm: rational intermediates, the worst case in the library.
    for (f, g) in [
        (&[3u32][..], &[2u32, 1][..]),
        (&[4], &[2, 1]),
        (&[2, 2], &[3, 1]),
    ] {
        let sf: Schur<BigRational> = Schur::monomial(part(f), BigRational::from(BigInt::from(1)));
        let sg: Schur<BigRational> = Schur::monomial(part(g), BigRational::from(BigInt::from(1)));
        let r = symfn::plethysm(&sf, &sg);
        let (b, d) = widest(r.terms().values().map(|c| c.numer().clone()));
        report(&format!("plethysm  s_{f:?}[s_{g:?}] numerators"), b, d);
    }

    // The intermediate power-sum form, where the denominators live.
    for deg in [16u32, 20, 24] {
        let lam: Vec<u32> = staircase(deg);
        let s: Schur<BigRational> = Schur::monomial(part(&lam), BigRational::from(BigInt::from(1)));
        let p: PowerSum<BigRational> = PowerSum::from_schur(&s);
        let (b, d) = widest(p.terms().values().map(|c| c.denom().clone()));
        report(&format!("s -> p denominators, degree {deg}"), b, d);
    }

    // --- what does coefficient arithmetic actually cost? ---------------------
    //
    // i128 against BigInt on identical work. This bounds how much of the
    // runtime is coefficient arithmetic at all: BigInt is heap-allocated and
    // several times the cost per operation, so if the whole workload slows by
    // only a little, the arithmetic is a small share of it.

    println!("\n=== what share of runtime is coefficient arithmetic?\n");
    println!("(i128 vs BigInt, identical work, rotated order)\n");

    fn pass<C: Ring>() -> usize {
        symfn::clear_caches();
        let mut n = 0;
        for parts in [&[8u32, 7, 6, 5, 4, 3][..], &[12, 10, 8, 6]] {
            let a: Schur<C> = Schur::monomial(part(parts), C::one());
            n += a.mul(&a).terms().len();
        }
        for deg in [18u32, 20] {
            let s: Schur<C> = Schur::monomial(part(&staircase(deg)), C::one());
            n += Monomial::<C>::from_schur(&s).terms().len();
        }
        n
    }

    let _ = pass::<i128>();
    let _ = pass::<BigInt>();
    let (mut t_i, mut t_b) = (0.0f64, 0.0f64);
    let (mut ci, mut cb) = (0, 0);
    for round in 0..6 {
        for k in 0..2 {
            if (round + k) % 2 == 0 {
                let t = Instant::now();
                ci = pass::<i128>();
                t_i += t.elapsed().as_secs_f64();
            } else {
                let t = Instant::now();
                cb = pass::<BigInt>();
                t_b += t.elapsed().as_secs_f64();
            }
        }
    }
    assert_eq!(ci, cb, "passes disagree");
    println!("  i128          {t_i:.4}s");
    println!("  BigInt        {t_b:.4}s   ({:.2}x)", t_b / t_i);
    println!(
        "\n  => coefficient arithmetic is at most ~{:.0}% of this workload,",
        (1.0 - t_i / t_b) * 100.0
    );
    println!("     since replacing it entirely with heap bignums costs that much.");

    // Rational is the expensive case (a gcd per operation), shown for scale.
    let _ = ratpass::<Rational>();
    let _ = ratpass::<BigRational>();
    let (mut r_i, mut r_b) = (0.0f64, 0.0f64);
    for round in 0..6 {
        for k in 0..2 {
            if (round + k) % 2 == 0 {
                let t = Instant::now();
                let _ = ratpass::<Rational>();
                r_i += t.elapsed().as_secs_f64();
            } else {
                let t = Instant::now();
                let _ = ratpass::<BigRational>();
                r_b += t.elapsed().as_secs_f64();
            }
        }
    }
    println!("\n  i128 Rational {r_i:.4}s");
    println!("  BigRational   {r_b:.4}s   ({:.2}x)", r_b / r_i);
}

fn ratpass<C: symfn::Plethystic>() -> usize {
    symfn::clear_caches();
    let mut n = 0;
    for deg in [14u32, 16] {
        let s: Schur<C> = Schur::monomial(part(&staircase(deg)), C::one());
        let p: PowerSum<C> = PowerSum::from_schur(&s);
        n += p.terms().len();
        n += p.to_schur().terms().len();
    }
    let sf: Schur<C> = Schur::monomial(part(&[3]), C::one());
    let sg: Schur<C> = Schur::monomial(part(&[2, 1]), C::one());
    n += symfn::plethysm(&sf, &sg).terms().len();
    n
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
