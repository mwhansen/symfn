//! Diagnostics for the `phi_slice` bottleneck, before optimising anything.
//!
//! ```text
//!   cargo run --release --example probe_gj -- 10
//! ```
//!
//! Three questions, each of which changes what the fix should be:
//!
//! 1. **Where do the atoms come from?** If `J_θ` in the p basis already carries
//!    linear atoms, the inner loop is stuck with ℚ(α) arithmetic. If they all
//!    come from the single `1/⟨J_θ,J_θ⟩` factor, the triple product is really
//!    integer-polynomial work wearing a fraction's clothes.
//! 2. **How much does an `AFrac` multiply cost against a scalar one?** That
//!    ratio is the ceiling on any evaluate-then-interpolate scheme.
//! 3. **Are there poles?** Every atom is `uα + v` with `u, v ≥ 0`, so `α > 0`
//!    should never hit one — which would make evaluation unconditionally safe
//!    rather than a gamble.

use std::time::Instant;

use symfn::afrac::AFrac;
use symfn::coeff::Rational;
use symfn::sym::SymFn;
use symfn::Ring;

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(9);

    println!("== 1. where the atoms in phi_slice come from ==");
    println!(
        "{:>3} {:>10} {:>14} {:>14} {:>12}",
        "n", "shapes", "J->p atoms", "J->p scales>1", "norm atoms"
    );
    for n in 1..=top {
        let (mut jatoms, mut jscale, mut natoms) = (0usize, 0usize, 0usize);
        for theta in symfn::partitions_of(n) {
            let j = symfn::jack_j_powersum::<i128>(&theta);
            for c in j.terms().values() {
                let (_, den, scale) = c.parts();
                jatoms += den.map(|(_, &m)| m as usize).sum::<usize>();
                jscale += usize::from(scale != 1);
            }
            natoms += symfn::jack_norm_j(&theta)
                .values()
                .map(|&m| m as usize)
                .sum::<usize>();
        }
        println!(
            "{n:>3} {:>10} {jatoms:>14} {jscale:>14} {natoms:>12}",
            symfn::partitions_of(n).len()
        );
    }

    println!();
    println!("== 2. cost of the inner-loop operation, AFrac vs scalar ==");
    // A representative pair: two J-coefficients at degree `top`, and the same
    // values at a numeric alpha.
    let theta = symfn::partitions_of(top).into_iter().next().unwrap();
    let j = symfn::jack_j_powersum::<i128>(&theta);
    let vals: Vec<AFrac<i128>> = j.terms().values().cloned().collect();
    let inverse: std::collections::BTreeMap<(u32, u32), i32> = symfn::jack_norm_j(&theta)
        .into_iter()
        .map(|(k, m)| (k, -m))
        .collect();
    let scaled: Vec<AFrac<i128>> = vals
        .iter()
        .map(|c| {
            let mut v = c.mul_factors(&inverse);
            v.reduce();
            v
        })
        .collect();

    const REPS: usize = 20000;
    let t0 = Instant::now();
    let mut acc = <AFrac<i128> as Ring>::zero();
    for i in 0..REPS {
        let a = &scaled[i % scaled.len()];
        let b = &vals[(i + 1) % vals.len()];
        let c = &vals[(i + 2) % vals.len()];
        acc.add_assign(&a.mul(b).mul(c));
        if i % 64 == 0 {
            acc = <AFrac<i128> as Ring>::zero(); // keep the accumulator realistic
        }
    }
    let afrac_ns = t0.elapsed().as_nanos() as f64 / REPS as f64;

    // The same shape of work at a numeric alpha, which is what an
    // evaluate-then-interpolate scheme would actually run.
    let alpha = Rational::from_int(3);
    let sv: Vec<Rational> = scaled
        .iter()
        .map(|c| c.eval_i128(&alpha).expect("alpha > 0 is never a pole"))
        .collect();
    let vv: Vec<Rational> = vals
        .iter()
        .map(|c| c.eval_i128(&alpha).expect("alpha > 0 is never a pole"))
        .collect();
    let t0 = Instant::now();
    let mut racc = Rational::zero();
    for i in 0..REPS {
        let a = sv[i % sv.len()];
        let b = vv[(i + 1) % vv.len()];
        let c = vv[(i + 2) % vv.len()];
        racc.add_assign(&a.mul(&b).mul(&c));
        if i % 64 == 0 {
            racc = Rational::zero();
        }
    }
    let scalar_ns = t0.elapsed().as_nanos() as f64 / REPS as f64;
    println!("  AFrac<i128>  mul,mul,add : {afrac_ns:>9.1} ns");
    println!("  Rational     mul,mul,add : {scalar_ns:>9.1} ns");
    println!(
        "  ratio                    : {:>9.1}x",
        afrac_ns / scalar_ns
    );

    // Modular arithmetic: the same shape of work with the cheapest possible
    // scalar. This is the number that decides whether evaluate-then-interpolate
    // is worth building, because `Rational` is a SLOW scalar -- it runs a
    // 128-bit gcd on every operation -- and it is not the one such a scheme
    // would have to use.
    const P: u64 = (1 << 61) - 1; // a Mersenne prime, so reduction is cheap
    let mv: Vec<u64> = vv
        .iter()
        .map(|r| r.numer().unsigned_abs() as u64 % P)
        .collect();
    let t0 = Instant::now();
    let mut macc: u64 = 0;
    for i in 0..REPS {
        let a = mv[i % mv.len()] as u128;
        let b = mv[(i + 1) % mv.len()] as u128;
        let c = mv[(i + 2) % mv.len()] as u128;
        let ab = (a * b % P as u128) as u128;
        let abc = (ab * c % P as u128) as u64;
        macc = (macc + abc) % P;
        if i % 64 == 0 {
            macc = 0;
        }
    }
    std::hint::black_box(macc);
    let modp_ns = t0.elapsed().as_nanos() as f64 / REPS as f64;
    println!("  mod p        mul,mul,add : {modp_ns:>9.1} ns");
    println!(
        "  AFrac / mod p            : {:>9.1}x",
        afrac_ns / modp_ns.max(1e-9)
    );

    println!();
    println!("== 2b. where the AFrac time goes, by component ==");
    // The inner loop is `a'.mul(b).mul(c)` then `+=`. `a'` carries every atom
    // (they all come from 1/<J,J>); `b` and `c` carry none. So the question is
    // how much of the cost is the atom bookkeeping rather than the polynomial.
    let plain: Vec<AFrac<i128>> = vals.clone();
    let t0 = Instant::now();
    for i in 0..REPS {
        let b = &plain[i % plain.len()];
        let c = &plain[(i + 1) % plain.len()];
        std::hint::black_box(b.mul(c));
    }
    let no_atoms = t0.elapsed().as_nanos() as f64 / REPS as f64;
    let t0 = Instant::now();
    for i in 0..REPS {
        let a = &scaled[i % scaled.len()];
        let c = &plain[(i + 1) % plain.len()];
        std::hint::black_box(a.mul(c));
    }
    let with_atoms = t0.elapsed().as_nanos() as f64 / REPS as f64;
    println!("  one mul, NO atoms on either side  : {no_atoms:>9.1} ns");
    println!("  one mul, atoms on one side        : {with_atoms:>9.1} ns");
    println!(
        "  the atoms cost                    : {:>9.1}x",
        with_atoms / no_atoms.max(1e-9)
    );

    println!();
    println!("== 2c. the global common denominator D = lcm_theta <J_theta,J_theta> ==");
    // If deg D stays near 2n, every phi_slice entry could be accumulated as a
    // plain integer polynomial over one fixed D -- no atoms, no lcm, no reduce.
    // If it grows like n^2, that restructuring blows up the degrees instead.
    println!(
        "{:>3} {:>8} {:>10} {:>12} {:>10}",
        "n", "deg D", "2n", "distinct", "deg a_rho"
    );
    for n in 1..=top {
        let mut lcm: std::collections::BTreeMap<(u32, u32), u32> = Default::default();
        let mut maxdeg = 0usize;
        for theta in symfn::partitions_of(n) {
            for (k, m) in symfn::jack_norm_j(&theta) {
                let e = lcm.entry(k).or_insert(0);
                *e = (*e).max(m as u32);
            }
            for c in symfn::jack_j_powersum::<i128>(&theta).terms().values() {
                maxdeg = maxdeg.max(c.degree().unwrap_or(0));
            }
        }
        let deg_d: u32 = lcm.iter().map(|(&(u, _), &m)| u.min(1) * m).sum();
        println!(
            "{n:>3} {deg_d:>8} {:>10} {:>12} {maxdeg:>10}",
            2 * n,
            lcm.len()
        );
    }

    println!();
    println!("== 3. poles: every atom is u*alpha + v with u,v >= 0 ==");
    let mut worst = (0u32, 0u32);
    let mut atoms = 0usize;
    for n in 1..=top {
        for theta in symfn::partitions_of(n) {
            for (&(u, v), _) in symfn::jack_norm_j(&theta).iter() {
                atoms += 1;
                worst = worst.max((u, v));
                assert!(u > 0 || v > 0, "the zero form is not an atom");
            }
            let j = symfn::jack_j_powersum::<i128>(&theta);
            for c in j.terms().values() {
                let (_, den, _) = c.parts();
                for (&(u, v), _) in den {
                    assert!(u > 0 || v > 0);
                }
            }
        }
    }
    println!("  {atoms} atoms inspected through degree {top}, largest {worst:?}");
    println!("  none can vanish for alpha > 0, so every positive alpha is a");
    println!("  legal evaluation point -- a proof, not a sampling argument.");
}
