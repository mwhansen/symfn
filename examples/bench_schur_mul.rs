//! Time the consumer side of a Schur product: what [`Schur::mul`] and
//! [`AutoLr`] add on top of the expansion the LR engine already holds.
//!
//! Every case here runs with the product expansions **warm**, so the engine's
//! own traversal is out of the picture and what remains is copying, term
//! accumulation and dispatch — the costs a caller pays again on every product
//! after the first, and the whole cost of a repeated one.
//!
//!   cargo run --release --example bench_schur_mul [reps] [case ...]
//!
//! Cases (all run when none is named):
//!
//! * `pair` — `s_μ · s_μ` as one product: `AutoLr::schur_product` (the owned
//!   expansion) against `Schur::mul` on one-term elements.
//! * `pair2` — `(s_μ + s_ν) · s_μ`: a second large product landing on a
//!   map the first one filled, with the count of terms new to it.
//! * `sum` — `(Σ_{λ ⊢ n} s_λ)²`: many small products accumulating into one
//!   element, where most terms land on a partition already present.
//! * `power` — `s_μ^k` by repeated multiplication: a growing accumulator
//!   times one term.
//! * `dispatch` — a rectangle, a two-row and a three-row pair, the products
//!   `AutoLr` computes off the general engine: first call, second call, and a
//!   sweep of `lr_coeff` over every term of the product.
//!
//! Times are the minimum over `reps` in-process runs. Compare two builds out
//! of process and interleaved, on AC power (the checklist at the end of
//! `docs/record/memory.md`).

use std::time::Instant;
use symfn::{clear_caches, partitions_of, AutoLr, LrBackend, Partition, Schur, SymFn};

fn p(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

fn one(mu: &Partition) -> Schur<i64> {
    Schur::monomial(mu.clone(), 1)
}

fn min_time<T>(reps: usize, mut f: impl FnMut() -> T) -> (f64, T) {
    let mut best = f64::INFINITY;
    let mut last = f();
    for _ in 1..reps {
        let t = Instant::now();
        last = f();
        best = best.min(t.elapsed().as_secs_f64());
    }
    if reps == 1 {
        best = 0.0;
    }
    (best, last)
}

fn ms(t: f64) -> String {
    format!("{:>9.3}ms", t * 1e3)
}

fn pair(reps: usize) {
    println!("pair: s_μ·s_μ, product warm");
    println!(
        "  {:<20} {:>9} {:>12} {:>12}",
        "μ", "terms", "schur_prod", "Schur::mul"
    );
    for mu in [
        p(&[8, 7, 6, 5, 4, 3]),
        p(&[10, 8, 6, 4]),
        p(&[12, 10, 8, 6]),
        p(&[16, 13, 10, 7]),
    ] {
        clear_caches();
        let warm = AutoLr.schur_product(&mu, &mu);
        let (t_owned, v) = min_time(reps, || AutoLr.schur_product(&mu, &mu));
        assert_eq!(v, warm);
        let (a, b) = (one(&mu), one(&mu));
        let (t_mul, prod) = min_time(reps, || a.mul(&b));
        assert_eq!(prod.terms().len(), warm.len());
        println!(
            "  {:<20} {:>9} {} {}",
            format!("{mu}"),
            warm.len(),
            ms(t_owned),
            ms(t_mul)
        );
    }
}

fn pair2(reps: usize) {
    println!("pair2: (s_μ + s_ν)·s_μ, products warm — the second pair lands on a full map");
    println!(
        "  {:<34} {:>9} {:>9} {:>12}",
        "μ, ν", "terms", "new", "Schur::mul"
    );
    for (mu, nu) in [
        (p(&[10, 8, 6, 4]), p(&[11, 8, 5, 4])),
        (p(&[8, 7, 6, 5, 4, 3]), p(&[10, 7, 6, 5, 3, 2])),
    ] {
        clear_caches();
        let mut left = one(&mu);
        left.add_term(nu.clone(), 1);
        let right = one(&mu);
        let warm = left.mul(&right);
        // How many of the second pair's terms are new to the map: the ones
        // outside the first pair's support.
        let first: std::collections::BTreeSet<Partition> = AutoLr
            .schur_product(&mu, &mu)
            .into_iter()
            .map(|(l, _)| l)
            .collect();
        let new = AutoLr
            .schur_product(&nu, &mu)
            .iter()
            .filter(|(l, _)| !first.contains(l))
            .count();
        let (t, prod) = min_time(reps, || left.mul(&right));
        assert_eq!(prod, warm);
        println!(
            "  {:<34} {:>9} {:>9} {}",
            format!("{mu}, {nu}"),
            warm.terms().len(),
            new,
            ms(t)
        );
    }
}

fn sum(reps: usize) {
    println!("sum: (Σ_{{λ⊢n}} s_λ)², products warm");
    println!(
        "  {:<6} {:>6} {:>9} {:>12}",
        "n", "pairs", "terms", "Schur::mul"
    );
    for n in [8u32, 10, 12] {
        clear_caches();
        let mut s: Schur<i64> = Schur::zero();
        for lam in partitions_of(n) {
            s.add_term(lam, 1);
        }
        let pairs = s.terms().len() * s.terms().len();
        let warm = s.mul(&s);
        let (t, prod) = min_time(reps, || s.mul(&s));
        assert_eq!(prod, warm);
        println!(
            "  {:<6} {:>6} {:>9} {}",
            n,
            pairs,
            warm.terms().len(),
            ms(t)
        );
    }
}

fn power(reps: usize) {
    println!("power: s_μ^k by repeated multiplication, products warm");
    println!(
        "  {:<10} {:>3} {:>9} {:>12}",
        "μ", "k", "terms", "Schur::mul"
    );
    for (mu, k) in [(p(&[2, 1]), 10usize), (p(&[3, 2, 1]), 6), (p(&[4, 2]), 6)] {
        clear_caches();
        let base = one(&mu);
        let run = || {
            let mut acc = base.clone();
            for _ in 1..k {
                acc = acc.mul(&base);
            }
            acc
        };
        let warm = run();
        let (t, prod) = min_time(reps, run);
        assert_eq!(prod, warm);
        println!(
            "  {:<10} {:>3} {:>9} {}",
            format!("{mu}"),
            k,
            warm.terms().len(),
            ms(t)
        );
    }
}

fn dispatch(reps: usize) {
    println!("dispatch: AutoLr products off the general engine; each rep starts cold");
    println!(
        "  {:<26} {:>9} {:>12} {:>12} {:>12} {:>12}",
        "s_μ·s_ν", "terms", "first", "second", "sweep", "sweep again"
    );
    for (mu, nu) in [
        (p(&[6, 6, 6, 6]), p(&[5, 5, 5])),
        (p(&[20, 16, 12]), p(&[20, 16])),
        (p(&[16, 12, 8, 4]), p(&[12, 10, 8])),
    ] {
        let mut best = [f64::INFINITY; 4];
        let mut terms = 0;
        for _ in 0..reps {
            clear_caches();
            let t = Instant::now();
            let first = AutoLr.schur_product(&mu, &nu);
            let t_first = t.elapsed().as_secs_f64();
            let t = Instant::now();
            let second = AutoLr.schur_product(&mu, &nu);
            let t_second = t.elapsed().as_secs_f64();
            assert_eq!(first, second);
            let sweep = || {
                first
                    .iter()
                    .map(|(lam, _)| AutoLr.lr_coeff(lam, &mu, &nu))
                    .sum::<u128>()
            };
            let t = Instant::now();
            let total = sweep();
            let t_sweep = t.elapsed().as_secs_f64();
            let t = Instant::now();
            let again = sweep();
            let t_again = t.elapsed().as_secs_f64();
            assert_eq!(total, first.iter().map(|(_, c)| c).sum::<u128>());
            assert_eq!(total, again);
            terms = first.len();
            for (b, t) in best.iter_mut().zip([t_first, t_second, t_sweep, t_again]) {
                *b = b.min(t);
            }
        }
        println!(
            "  {:<26} {:>9} {} {} {} {}",
            format!("s{mu}·s{nu}"),
            terms,
            ms(best[0]),
            ms(best[1]),
            ms(best[2]),
            ms(best[3])
        );
    }
}

fn main() {
    let mut args = std::env::args().skip(1).peekable();
    let reps: usize = match args.peek() {
        Some(a) if a.parse::<usize>().is_ok() => args.next().unwrap().parse().unwrap(),
        _ => 5,
    };
    let cases: Vec<String> = args.collect();
    let want = |name: &str| cases.is_empty() || cases.iter().any(|c| c == name);
    if want("pair") {
        pair(reps);
    }
    if want("pair2") {
        pair2(reps);
    }
    if want("sum") {
        sum(reps);
    }
    if want("power") {
        power(reps);
    }
    if want("dispatch") {
        dispatch(reps);
    }
}
