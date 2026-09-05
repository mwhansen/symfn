//! The two routes of s → m against each other, by term count.
//!
//! `Monomial::from_schur` runs one `kostka(λ, μ)` per pair when a degree has
//! few Schur terms and one Pieri trie over every μ ⊢ n when it has many. The
//! per-pair cost grows with the term count and the trie's does not, so the
//! question this harness answers is where they cross, and by how much the
//! trie wins on the inputs the hub routing produces — X → s → m hands this
//! function a Schur element with most of p(n) terms.
//!
//! ```text
//!   cargo build --release --example bench_s2m
//!   ./target/release/examples/bench_s2m 16 20 24
//! ```
//!
//! For each degree, from cold caches every time:
//!
//! * `single_sum` — every s_λ, λ ⊢ n, converted alone, and the times summed.
//!   A single term dispatches per-pair at n ≥ 8, so this is what the per-pair
//!   route costs on full support, whatever the dispatch says.
//! * `terms=k` — `k` shapes thinned evenly across `partitions_of`, converted
//!   as one element, for `k` on both sides of the dispatch threshold. This is
//!   what the dispatch does at that count; which route it took reads back
//!   from the time — a batched row is flat in `k`.
//! * `full` — all p(n) shapes with coefficient 1, i.e. h_{1ⁿ}'s Schur
//!   expansion divided by f^λ, which is the many-term case.
//! * `h1n_to_m`, `e1n_to_m`, `p1n_to_m` — h_{1ⁿ}, e_{1ⁿ}, p_{1ⁿ} → m through
//!   `convert`, the hub route end to end, which is what Sage's `m(h[1]^n)`
//!   reaches through the adapter.
//!
//! Self-relative, no oracle: `tests/sage_oracle.rs` and the batched-versus-
//! per-pair test in `src/convert.rs` say whether the answers agree.

use std::time::Instant;

use symfn::{
    clear_caches, convert, partitions_of, Elementary, FromSchur, Homogeneous, Monomial, Partition,
    PowerSum, Schur, SymFn,
};

fn time<T>(f: impl FnOnce() -> T) -> (f64, T) {
    clear_caches();
    let t = Instant::now();
    let out = f();
    (t.elapsed().as_secs_f64(), out)
}

fn main() {
    let degrees: Vec<u32> = std::env::args()
        .skip(1)
        .map(|a| a.parse().expect("degree"))
        .collect();
    let degrees = if degrees.is_empty() {
        vec![12, 16, 20]
    } else {
        degrees
    };
    println!("n\tcase\tseconds\tterms_in\tterms_out");
    for n in degrees {
        let parts = partitions_of(n);
        let pn = parts.len();

        let mut single_sum = 0.0;
        for lam in &parts {
            let s: Schur<i128> = Schur::monomial(lam.clone(), 1);
            let (dt, _) = time(|| Monomial::<i128>::from_schur(&s));
            single_sum += dt;
        }
        println!("{n}\tsingle_sum\t{single_sum:.6}\t{pn}\t-");

        // Around the dispatch threshold `3n − 30`, on both sides of it.
        let t = (3 * n as usize).saturating_sub(30).max(2);
        for k in [t / 2, t - 1, t, 2 * t, 4 * t] {
            let k = k.clamp(1, pn);
            let mut s: Schur<i128> = Schur::zero();
            for i in 0..k {
                s.add_term(parts[i * pn / k].clone(), 1);
            }
            let (dt, m) = time(|| Monomial::<i128>::from_schur(&s));
            println!("{n}\tterms={k}\t{dt:.6}\t{k}\t{}", m.terms().len());
        }

        let mut full: Schur<i128> = Schur::zero();
        for lam in &parts {
            full.add_term(lam.clone(), 1);
        }
        let (dt, m) = time(|| Monomial::<i128>::from_schur(&full));
        println!("{n}\tfull\t{dt:.6}\t{pn}\t{}", m.terms().len());

        let ones = Partition::new(std::iter::repeat_n(1, n as usize));
        let h: Homogeneous<i128> = Homogeneous::monomial(ones.clone(), 1);
        let (dt, m) = time(|| convert::<i128, _, Monomial<i128>>(&h));
        println!("{n}\th1n_to_m\t{dt:.6}\t1\t{}", m.terms().len());
        let e: Elementary<i128> = Elementary::monomial(ones.clone(), 1);
        let (dt, m) = time(|| convert::<i128, _, Monomial<i128>>(&e));
        println!("{n}\te1n_to_m\t{dt:.6}\t1\t{}", m.terms().len());
        let pw: PowerSum<i128> = PowerSum::monomial(ones, 1);
        let (dt, m) = time(|| convert::<i128, _, Monomial<i128>>(&pw));
        println!("{n}\tp1n_to_m\t{dt:.6}\t1\t{}", m.terms().len());
    }
}
