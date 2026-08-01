//! The principal-specialization checksum on Littlewood–Richardson products,
//! with the negative control that makes a PASS mean something.
//!
//! Evaluating `s_μ · s_ν = Σ_λ c^λ_{μν} s_λ` at `x = (1,1,…,1)` turns the
//! identity into a scalar equation, and each `s_λ(1ⁿ)` follows from the
//! hook-content formula `∏_{(i,j)∈λ} (n + j − i) / h(i,j)` — which needs only
//! the diagram and its conjugate, and so shares no code with the LR machinery.
//! That independence is the point: both orientations of one frontier engine
//! fail together, so a second orientation is consistency rather than evidence
//! (`docs/policies/validation.md` V3).
//!
//! This is the check that holds the shapes no external oracle finishes
//! (`docs/record/littlewood-richardson.md`), so it is also a **lossy** one — a
//! weighted sum, not a proof. V7 requires a lossy check to ship with its
//! control, and `the_checksum_detects_every_single_coefficient_error` is it:
//! a checksum that silently always passed would be worse than no check.
//! `examples/verify_specialization.rs` runs the same arithmetic on the shapes
//! that take minutes; this runs on every `cargo test`.

use symfn::{partitions_of, AutoLr, LrBackend, Partition};

/// A Mersenne prime, so a product of two residues fits in `u128`.
const P: u64 = (1 << 61) - 1;

fn mul(a: u64, b: u64) -> u64 {
    ((a as u128 * b as u128) % P as u128) as u64
}

fn pow(mut a: u64, mut e: u64) -> u64 {
    let mut r = 1;
    while e > 0 {
        if e & 1 == 1 {
            r = mul(r, a);
        }
        a = mul(a, a);
        e >>= 1;
    }
    r
}

/// `s_λ(1ⁿ) mod P` by the hook-content formula.
///
/// Zero when `ℓ(λ) > n`, which the formula produces on its own: the cell
/// `(n, 1)` then has content `1 − n − 1`, making the numerator `n + j − i = 0`.
fn schur_at_ones(lambda: &Partition, n: u64) -> u64 {
    let conj = lambda.conjugate();
    let (mut num, mut den) = (1u64, 1u64);
    for (i, &row) in lambda.parts().iter().enumerate() {
        for j in 0..row as usize {
            let hook = (row as usize - j) + (conj.part(j) as usize - i) - 1;
            let content = n as i64 + j as i64 - i as i64;
            if content <= 0 {
                return 0;
            }
            num = mul(num, content as u64 % P);
            den = mul(den, hook as u64 % P);
        }
    }
    mul(num, pow(den, P - 2))
}

/// The right-hand side `Σ_λ c^λ_{μν} s_λ(1ⁿ)`, with term `skip` shifted by
/// `delta` — `None` leaves the expansion alone.
fn rhs_at(prod: &[(Partition, u128)], n: u64, tweak: Option<(usize, i128)>) -> u64 {
    let mut acc = 0u64;
    for (i, (lambda, c)) in prod.iter().enumerate() {
        let c = match tweak {
            Some((k, delta)) if k == i => (*c as i128 + delta) as u128,
            _ => *c,
        };
        acc = (acc + mul((c % P as u128) as u64, schur_at_ones(lambda, n))) % P;
    }
    acc
}

/// The `n` at which a product of two shapes is evaluated. Past `ℓ(μ)+ℓ(ν)` every
/// λ in the support has a nonzero `s_λ(1ⁿ)`, so no term is invisible to the sum.
fn points(mu: &Partition, nu: &Partition) -> Vec<u64> {
    let base = (mu.len() + nu.len()) as u64;
    [1, 3, 8, 21, 55].iter().map(|d| base + d).collect()
}

/// **The checksum holds on every product through `|μ| + |ν| ≤ 8`.**
#[test]
fn principal_specialization_agrees_with_the_lr_expansion() {
    let mut checked = 0usize;
    for a in 0..=8u32 {
        for b in 0..=(8 - a) {
            for mu in partitions_of(a) {
                for nu in partitions_of(b) {
                    let prod = AutoLr.schur_product(&mu, &nu);
                    for n in points(&mu, &nu) {
                        let want = mul(schur_at_ones(&mu, n), schur_at_ones(&nu, n));
                        assert_eq!(rhs_at(&prod, n, None), want, "s{mu}*s{nu} at x = 1^{n}");
                        checked += 1;
                    }
                }
            }
        }
    }
    assert!(checked > 1000, "expected a real sweep, got {checked}");
}

/// **The negative control: every single-coefficient error is caught.**
///
/// Perturbing one coefficient by ±1, for every term in turn, shows that no
/// position in the expansion is a blind spot — which a checksum built from one
/// evaluation point would not give. The assertion is `caught == total`, so a
/// blind spot is a red test rather than a number nobody reads.
#[test]
fn the_checksum_detects_every_single_coefficient_error() {
    let cases = [
        (Partition::new([2u32, 1]), Partition::new([2u32, 1])),
        (Partition::new([3u32, 2, 1]), Partition::new([2u32, 1])),
        (Partition::new([3u32, 2, 1]), Partition::new([3u32, 2, 1])),
    ];
    let (mut caught, mut total) = (0usize, 0usize);
    for (mu, nu) in &cases {
        let prod = AutoLr.schur_product(mu, nu);
        let ns = points(mu, nu);
        let lhs: Vec<u64> = ns
            .iter()
            .map(|&n| mul(schur_at_ones(mu, n), schur_at_ones(nu, n)))
            .collect();
        for k in 0..prod.len() {
            for delta in [1i128, -1] {
                if prod[k].1 as i128 + delta < 0 {
                    continue;
                }
                total += 1;
                let detected = ns
                    .iter()
                    .zip(&lhs)
                    .any(|(&n, &want)| rhs_at(&prod, n, Some((k, delta))) != want);
                assert!(
                    detected,
                    "blind spot: s{mu}*s{nu}, term {k} ({}) shifted by {delta}",
                    prod[k].0
                );
                caught += 1;
            }
        }
    }
    assert_eq!(caught, total, "every perturbation must be detected");
    assert!(total > 100, "the control must be a real sweep, got {total}");
}
