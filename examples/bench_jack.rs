//! The Jack ladder: whole-degree tables, single shapes, norms, and the two
//! engines against each other.
//!
//! ```text
//!   cargo run --release --example bench_jack -- 13
//! ```
//!
//! The Sage side is `docs/record/jack.md`, measured on this machine:
//! whole degree `P → m` costs 1.63 s at n = 8, 5.64 s at 9, 21.29 s at 10,
//! 84.52 s at 11 and dies at 12; the norms table costs >360 s at n = 12. ⚠️
//! Those are mains rows; compare mains to mains, since this machine drifts
//! about 1.8× on battery.
//!
//! Every row runs the two-width ladder — `i128` against `Rational` — and
//! asserts they agree term for term, so a wrapped fixed-width intermediate is a
//! failure and not a silently different number.

use std::time::Instant;

use symfn::afrac::AFrac;
use symfn::coeff::Rational;
use symfn::sym::{Monomial, SymFn};
use symfn::{Partition, Ring};

/// Peak numerator degree, coefficient bits, denominator atoms and scalar bits
/// over a whole expansion — the four columns `scripts/spec_jack_swell.py`
/// predicted before any of this existed.
fn swell(f: &Monomial<AFrac<i128>>) -> (usize, u32, usize, u32) {
    let mut peak = (0usize, 0u32, 0usize, 0u32);
    for c in f.terms().values() {
        let (num, den, scale) = c.parts();
        peak.0 = peak.0.max(num.len().saturating_sub(1));
        peak.1 = peak.1.max(
            num.iter()
                .map(|x| 128 - x.unsigned_abs().leading_zeros())
                .max()
                .unwrap_or(0),
        );
        peak.2 = peak.2.max(den.map(|(_, &m)| m as usize).sum::<usize>());
        peak.3 = peak.3.max(128 - scale.leading_zeros());
    }
    peak
}

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(11);
    // The `Rational` column is the two-width honesty check, and it is the slow
    // one; past this degree only `i128` runs, and the table says so with a `—`.
    let widths: u32 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(top);

    println!("== whole-degree P -> m, the unit Sage has no entry point for ==");
    println!(
        "{:>3} {:>6} {:>11} {:>12} {:>8} {:>7} {:>6} {:>7}",
        "n", "p(n)", "i128(s)", "Rational(s)", "deg", "bits", "atoms", "scale"
    );
    for n in 1..=top {
        symfn::clear_caches();
        let t0 = Instant::now();
        let z: Vec<_> = symfn::jack_table::<i128>(n);
        let integral = t0.elapsed().as_secs_f64();

        let exact = if n <= widths {
            symfn::clear_caches();
            let t0 = Instant::now();
            let q: Vec<_> = symfn::jack_table::<Rational>(n);
            let secs = t0.elapsed().as_secs_f64();
            // The two widths must agree, or `i128` wrapped somewhere.
            for ((la, a), (lb, b)) in z.iter().zip(q.iter()) {
                assert_eq!(la, lb);
                for (mu, ca) in a.terms() {
                    let want = ca.clone().eval_at(&Rational::from_int(7));
                    let got = b.coeff(mu).eval(&Rational::from_int(7));
                    assert_eq!(
                        got, want,
                        "i128 and Rational disagree at degree {n}, {la} / {mu}"
                    );
                }
            }
            format!("{secs:>12.4}")
        } else {
            format!("{:>12}", "—")
        };

        let mut peak = (0usize, 0u32, 0usize, 0u32);
        for (_, a) in z.iter() {
            let p = swell(a);
            peak = (
                peak.0.max(p.0),
                peak.1.max(p.1),
                peak.2.max(p.2),
                peak.3.max(p.3),
            );
        }
        println!(
            "{n:>3} {:>6} {integral:>11.4} {exact} {:>8} {:>7} {:>6} {:>7}",
            z.len(),
            peak.0,
            peak.1,
            peak.2,
            peak.3
        );
    }

    println!();
    println!("== the other three units Sage prices identically (see scripts/bench_jack.py) ==");
    println!(
        "{:>3} {:>12} {:>12} {:>12} {:>12}",
        "n", "P->m(s)", "J->m(s)", "J->p(s)", "norms(s)"
    );
    for n in 1..=top {
        symfn::clear_caches();
        let t0 = Instant::now();
        let _: Vec<_> = symfn::jack_table::<i128>(n);
        let pm = t0.elapsed().as_secs_f64();

        symfn::clear_caches();
        let t0 = Instant::now();
        let _: Vec<_> = symfn::jack_j_table::<i128>(n);
        let jm = t0.elapsed().as_secs_f64();

        symfn::clear_caches();
        let t0 = Instant::now();
        let _: Vec<_> = symfn::jack_powersum_table::<i128>(n);
        let jp = t0.elapsed().as_secs_f64();

        let t0 = Instant::now();
        let mut atoms = 0usize;
        for lambda in symfn::partitions_of(n) {
            atoms += symfn::jack_norm_j(&lambda)
                .values()
                .map(|&m| m as usize)
                .sum::<usize>();
        }
        assert!(atoms > 0);
        let nm = t0.elapsed().as_secs_f64();
        println!("{n:>3} {pm:>12.5} {jm:>12.5} {jp:>12.5} {nm:>12.6}");
    }

    println!();
    println!("== the single shape (n), which IS the degree for Sage ==");
    println!("{:>3} {:>11} {:>9}", "n", "P_(n) (s)", "terms");
    for n in 1..=top {
        symfn::clear_caches();
        let lambda = Partition::new([n]);
        let t0 = Instant::now();
        let p: Monomial<AFrac<i128>> = symfn::jack_p(&lambda);
        println!(
            "{n:>3} {:>11.4} {:>9}",
            t0.elapsed().as_secs_f64(),
            p.terms().len()
        );
    }

    println!();
    println!("== the two engines: eigenoperator vs branching, whole degree ==");
    println!(
        "{:>3} {:>12} {:>12} {:>8}",
        "n", "LB(s)", "branch(s)", "ratio"
    );
    for n in 1..=top.min(10) {
        symfn::clear_caches();
        let parts = symfn::partitions_of(n);
        let t0 = Instant::now();
        let a: Vec<Monomial<AFrac<i128>>> = parts.iter().map(symfn::jack::jack_p_lb).collect();
        let lb = t0.elapsed().as_secs_f64();
        let t0 = Instant::now();
        let b: Vec<Monomial<AFrac<i128>>> =
            parts.iter().map(symfn::jack::jack_p_branching).collect();
        let br = t0.elapsed().as_secs_f64();
        assert_eq!(a, b, "the two engines disagree at degree {n}");
        println!("{n:>3} {lb:>12.4} {br:>12.4} {:>8.2}", br / lb.max(1e-9));
    }

    println!();
    println!("== Stanley's table: every <J_la J_mu, J_nu> with |la| = |mu| = k ==");
    println!(
        "{:>3} {:>10} {:>12} {:>10} {:>9}  positivity",
        "k", "triples", "total(s)", "nonzero", "max deg"
    );
    for k in 1..=top.min(9) {
        symfn::clear_caches();
        let parts = symfn::partitions_of(k);
        let big = symfn::partitions_of(2 * k);
        let triples = (parts.len() * parts.len() * big.len()) as u64;
        let t0 = Instant::now();
        let table: Vec<_> = symfn::stanley_table::<i128>(k);
        let (mut nonzero, mut maxdeg, mut negative) = (0u64, 0usize, 0u64);
        for (_, _, _, g) in &table {
            nonzero += 1;
            match g.clone().into_natural_poly() {
                Some(c) => maxdeg = maxdeg.max(c.len().saturating_sub(1)),
                // OPEN CONJECTURE. A violation is a result to report, not a bug
                // to fix.
                None => negative += 1,
            }
        }
        let verdict = if negative == 0 {
            "all in N[alpha]".to_string()
        } else {
            format!("*** {negative} LEFT N[alpha] -- REPORT ***")
        };
        println!(
            "{k:>3} {triples:>10} {:>12.3} {nonzero:>10} {maxdeg:>9}  {verdict}",
            t0.elapsed().as_secs_f64()
        );
    }

    println!();
    println!("== norms: closed products of 2n linear factors (Sage: >360s at n=12) ==");
    println!("{:>3} {:>12} {:>8}", "n", "table(s)", "shapes");
    for n in 1..=top.max(12) {
        let t0 = Instant::now();
        let parts = symfn::partitions_of(n);
        let mut atoms = 0usize;
        for lambda in &parts {
            atoms += symfn::jack_norm_j(lambda)
                .values()
                .map(|&m| m as usize)
                .sum::<usize>();
        }
        assert!(atoms > 0);
        println!(
            "{n:>3} {:>12.6} {:>8}",
            t0.elapsed().as_secs_f64(),
            parts.len()
        );
    }
}

/// `AFrac<i128>` has no `eval` (that needs a field), so the two-width check
/// goes through ℚ on both sides.
trait EvalAt {
    fn eval_at(self, x: &Rational) -> Option<Rational>;
}

impl EvalAt for AFrac<i128> {
    fn eval_at(self, x: &Rational) -> Option<Rational> {
        let (num, den, scale) = self.parts();
        let mut acc = Rational::zero();
        for c in num.iter().rev() {
            acc = acc.mul(x);
            acc.add_assign(&Rational::from_int(*c));
        }
        let mut d = Rational::from_u128(scale);
        for (&(u, v), &m) in den {
            let mut f = Rational::from_u128(u as u128).mul(x);
            f.add_assign(&Rational::from_u128(v as u128));
            if f.is_zero() {
                return None;
            }
            for _ in 0..m {
                d = d.mul(&f);
            }
        }
        Some(symfn::coeff::Field::div(&acc, &d))
    }
}
