//! Where does a (q,t)-Kostka table spend its time?
//!
//! ```text
//!   cargo run --release --example profile_qtk -- 9
//! ```
//!
//! The route has five separable phases and they are large, so this splits them
//! by the clock rather than by the sampler — a flat profile would name
//! `QtPoly::mul_binomial` and `divide_by_factor` and leave the question of
//! *which phase called them* open, which is the question that decides what to
//! change. `examples/profile_mac.rs` is the sampler harness for when it is the
//! leaf that is in doubt.

use std::time::Instant;

use symfn::coeff::Rational;
use symfn::convert::{FromSchur, ToSchur};
use symfn::sym::{Monomial, PowerSum, Schur, SymFn};
use symfn::{Frac, QtPoly};

type C = Rational;

fn main() {
    let top: u32 = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(8);
    // `-- 9 loop` runs p -> s alone, forever, for the sampler.
    if std::env::args().nth(2).as_deref() == Some("loop") {
        let parts = symfn::partitions_of(top);
        let staged: Vec<PowerSum<Frac<C>>> = parts
            .iter()
            .map(|mu| {
                let j: Monomial<Frac<C>> = symfn::macdonald_j(mu);
                phi_t(&PowerSum::from_schur(&j.to_schur()))
            })
            .collect();
        loop {
            for p in &staged {
                std::hint::black_box(p.to_schur());
            }
        }
    }

    println!(
        "{:>3} {:>6} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8}",
        "n", "p(n)", "J", "m->s", "s->p", "phi_t", "p->s", "reduce", "total"
    );
    for n in 1..=top {
        symfn::clear_caches();
        let parts = symfn::partitions_of(n);
        let (mut tj, mut tms, mut tsp, mut tphi, mut tps, mut tred) =
            (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        let mut terms = 0usize;

        for mu in &parts {
            let t0 = Instant::now();
            let j: Monomial<Frac<C>> = symfn::macdonald_j(mu);
            tj += t0.elapsed().as_secs_f64();

            let t0 = Instant::now();
            let s: Schur<Frac<C>> = j.to_schur();
            tms += t0.elapsed().as_secs_f64();

            let t0 = Instant::now();
            let p: PowerSum<Frac<C>> = PowerSum::from_schur(&s);
            tsp += t0.elapsed().as_secs_f64();

            let t0 = Instant::now();
            let scaled = phi_t(&p);
            tphi += t0.elapsed().as_secs_f64();

            let t0 = Instant::now();
            let back: Schur<Frac<C>> = scaled.to_schur();
            tps += t0.elapsed().as_secs_f64();

            let t0 = Instant::now();
            for c in back.terms().values() {
                let k: QtPoly<C> = c.clone().into_poly().expect("integral");
                terms += k.len();
            }
            tred += t0.elapsed().as_secs_f64();
        }
        let total = tj + tms + tsp + tphi + tps + tred;
        println!(
            "{n:>3} {:>6} {tj:>8.4} {tms:>8.4} {tsp:>8.4} {tphi:>8.4} {tps:>8.4} {tred:>8.4} \
             {total:>8.4}   ({terms} monomials)",
            parts.len()
        );
    }
}

/// The same loop `qtkostka::invert_s_basis` runs, duplicated here because it is
/// private and the point is to time it in isolation.
fn phi_t(p: &PowerSum<Frac<C>>) -> PowerSum<Frac<C>> {
    let mut out = PowerSum::zero();
    for (nu, c) in p.terms() {
        let mut factors = std::collections::BTreeMap::new();
        for &part in nu.parts() {
            *factors.entry((0u32, part)).or_insert(0i32) -= 1;
        }
        let mut v = c.mul_factors(&factors);
        v.reduce();
        out.add_term(nu.clone(), v);
    }
    out
}
