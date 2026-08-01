//! The ten basis transitions of [`symfn::convert`], swept across the shape
//! families that separate them, and repeated to a deadline so a sampling
//! profiler can attribute the time.
//!
//! `bench_ops` already times all ten — on **one** shape, `[6,5,4,3,2]`. That is
//! the trap this file exists to avoid: the Jacobi–Trudi directions are
//! exponential in a matrix whose size is ℓ(λ) for s → h and λ₁ for s → e, so a
//! single balanced shape is precisely the input on which the two look alike and
//! both look healthy. A 200x regression in s → e survived every degree of the
//! `compare_sage.py` ladder for exactly that reason
//! ([docs/record/transitions.md](../docs/record/transitions.md), "s → h and
//! s → e"). The families below are chosen to be the extremes of ℓ(λ) against
//! λ₁ — row, column, hook, staircase, rectangle — so no direction can hide.
//!
//! Two modes. The sweep, one tab-separated line per (transition, shape):
//!
//! ```text
//!   cargo build --release --example profile_convert
//!   ./target/release/examples/profile_convert sweep 20
//! ```
//!
//! Columns: tag, transition, shape, seconds, input terms, output terms. Time
//! reads against work, as in `bench_ops`; a transition that is slow because its
//! input has 627 terms is not the same finding as one slow on a single term.
//!
//! And one case held under a profiler, which is what actually attributes the
//! cost to a function:
//!
//! ```text
//!   cargo build --profile profiling --example profile_convert
//!   ./target/profiling/examples/profile_convert loop s2e 24 20 &
//!   sample $! 15 -mayDie -f /tmp/convert.txt
//! ```
//!
//! Caches are cleared every iteration, as in `profile_wide`: without that the
//! second and every later pass is a memo hit and the profile is all hashing.
//! The per-iteration time is printed so it can be checked against the sweep row
//! for the same case — if they disagree, the profile is of something else and
//! the attribution under it means nothing.
//!
//! Both modes are self-relative. There is no oracle here; `compare_sage.py` and
//! `check_backend.py` are what say whether a number is good.

use std::time::{Duration, Instant};

use symfn::{
    clear_caches, convert, Elementary, Forgotten, Homogeneous, Monomial, Partition, PowerSum,
    Rational, Ring, Schur, SymFn,
};

/// The shape families, as (name, builder). Each is a different corner of the
/// ℓ(λ)-against-λ₁ trade the determinant directions turn on: the row and the
/// column are the two extremes, the hook is bad for *both* at once (each about
/// n/2, the family the h ↔ e flip cannot rescue), and the staircase and
/// rectangle are the balanced controls the older ladder measured exclusively.
fn shapes(n: u32) -> Vec<(&'static str, Partition)> {
    let row = Partition::new([n]);
    let column = Partition::new(std::iter::repeat_n(1, n as usize));
    let hook = {
        let arm = n / 2;
        let mut v = vec![arm + n % 2];
        v.extend(std::iter::repeat_n(1, (n - arm - n % 2) as usize));
        Partition::new(v)
    };
    // The largest staircase (k, k-1, …, 1) inside degree n, padded with the
    // remainder as extra 1s so every family is the same degree and the term
    // counts stay comparable.
    let staircase = {
        let mut v = Vec::new();
        let mut left = n;
        let mut k = 1;
        while left >= k {
            v.push(k);
            left -= k;
            k += 1;
        }
        v.reverse();
        v.extend(std::iter::repeat_n(1, left as usize));
        Partition::new(v)
    };
    // Widest rectangle of the degree, falling back to the row when n is prime.
    let rectangle = {
        let rows = (1..=n)
            .filter(|&d| n.is_multiple_of(d))
            .find(|&d| d * d >= n)
            .unwrap_or(1);
        Partition::new(std::iter::repeat_n(n / rows, rows as usize))
    };
    vec![
        ("row", row),
        ("column", column),
        ("hook", hook),
        ("staircase", staircase),
        ("rectangle", rectangle),
    ]
}

/// Time one transition on one shape, printing `tag name shape secs in out`.
///
/// The *input* is built untimed by converting s_λ out into the source basis,
/// so a reverse direction is measured on the element it would really be handed
/// rather than on a one-term artefact. Caches are cleared after that setup and
/// before the timer, so the row is cold.
fn row<A, B>(tag: &str, name: &str, shape: &Partition, input: A)
where
    A: SymFn<Rational> + symfn::ToSchur<Rational>,
    B: SymFn<Rational> + symfn::FromSchur<Rational>,
{
    let n_in = input.terms().len();
    clear_caches();
    let t = Instant::now();
    let out: B = convert::<Rational, A, B>(&input);
    let dt = t.elapsed().as_secs_f64();
    println!(
        "{tag}\t{name}\t{shape}\t{dt:.6}\t{n_in}\t{}",
        out.terms().len()
    );
}

/// Every transition on one shape.
///
/// All ten run over ℚ rather than ℤ, which is not the cheapest choice: only
/// s → p needs division, and the integral directions would be faster in `i128`
/// (`bench_ops` times them there). Uniform ℚ is deliberate here because this
/// harness compares *transitions against each other*, and a table mixing two
/// coefficient rings would attribute the ring's cost to the transition. The
/// integral rows are the ones to read against `bench_ops`, not across this
/// table.
fn sweep_shape(tag: &str, shape: &Partition) {
    let s: Schur<Rational> = Schur::monomial(shape.clone(), Rational::one());

    row::<Schur<Rational>, Monomial<Rational>>(tag, "s2m", shape, s.clone());
    row::<Schur<Rational>, Homogeneous<Rational>>(tag, "s2h", shape, s.clone());
    row::<Schur<Rational>, Elementary<Rational>>(tag, "s2e", shape, s.clone());
    row::<Schur<Rational>, PowerSum<Rational>>(tag, "s2p", shape, s.clone());
    row::<Schur<Rational>, Forgotten<Rational>>(tag, "s2f", shape, s.clone());

    // Reverse directions, each fed the real expansion of s_λ in its basis.
    let m: Monomial<Rational> = convert(&s);
    row::<Monomial<Rational>, Schur<Rational>>(tag, "m2s", shape, m);
    let h: Homogeneous<Rational> = convert(&s);
    row::<Homogeneous<Rational>, Schur<Rational>>(tag, "h2s", shape, h);
    let e: Elementary<Rational> = convert(&s);
    row::<Elementary<Rational>, Schur<Rational>>(tag, "e2s", shape, e);
    let p: PowerSum<Rational> = convert(&s);
    row::<PowerSum<Rational>, Schur<Rational>>(tag, "p2s", shape, p);
    let f: Forgotten<Rational> = convert(&s);
    row::<Forgotten<Rational>, Schur<Rational>>(tag, "f2s", shape, f);
}

/// One transition, repeated until `secs` have passed, for a sampling profiler.
fn profile_loop(which: &str, shape: &Partition, secs: u64) {
    let s: Schur<Rational> = Schur::monomial(shape.clone(), Rational::one());
    // Reverse directions need their input built once, outside the loop: it is
    // setup, not the thing being profiled, and rebuilding it every iteration
    // would put the *forward* transition in the samples too.
    let m: Monomial<Rational> = convert(&s);
    let h: Homogeneous<Rational> = convert(&s);
    let e: Elementary<Rational> = convert(&s);
    let pp: PowerSum<Rational> = convert(&s);
    let f: Forgotten<Rational> = convert(&s);

    eprintln!(
        "profiling {which} on s{shape} for {secs}s (pid {})",
        std::process::id()
    );
    let deadline = Instant::now() + Duration::from_secs(secs);
    let (mut iters, mut terms) = (0u32, 0usize);
    let start = Instant::now();
    while Instant::now() < deadline {
        clear_caches();
        terms = match which {
            "s2m" => convert::<Rational, _, Monomial<Rational>>(&s).terms().len(),
            "s2h" => convert::<Rational, _, Homogeneous<Rational>>(&s)
                .terms()
                .len(),
            "s2e" => convert::<Rational, _, Elementary<Rational>>(&s)
                .terms()
                .len(),
            "s2p" => convert::<Rational, _, PowerSum<Rational>>(&s).terms().len(),
            "s2f" => convert::<Rational, _, Forgotten<Rational>>(&s)
                .terms()
                .len(),
            "m2s" => convert::<Rational, _, Schur<Rational>>(&m).terms().len(),
            "h2s" => convert::<Rational, _, Schur<Rational>>(&h).terms().len(),
            "e2s" => convert::<Rational, _, Schur<Rational>>(&e).terms().len(),
            "p2s" => convert::<Rational, _, Schur<Rational>>(&pp).terms().len(),
            "f2s" => convert::<Rational, _, Schur<Rational>>(&f).terms().len(),
            other => panic!("unknown transition {other}"),
        };
        iters += 1;
    }
    let per = start.elapsed().as_secs_f64() / f64::from(iters);
    eprintln!("{iters} iterations, {terms} terms, {per:.6}s per conversion");
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let mode = args.first().map_or("sweep", String::as_str);

    if mode == "loop" {
        let which = args.get(1).map_or("s2e", String::as_str);
        let n: u32 = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(20);
        let secs: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(20);
        // A bare degree selects the row, the family that broke s → e; an
        // explicit comma list overrides it.
        let shape = args.get(2).filter(|a| a.contains(',')).map_or_else(
            || Partition::new([n]),
            |a| Partition::new(a.split(',').map(|x| x.trim().parse::<u32>().expect("part"))),
        );
        profile_loop(which, &shape, secs);
        return;
    }

    let n: u32 = args.get(1).and_then(|s| s.parse().ok()).unwrap_or(20);
    let tag = args.get(2).map_or("profile", String::as_str);
    for (family, shape) in shapes(n) {
        eprintln!("--- {family} {shape} ---");
        sweep_shape(tag, &shape);
    }
}
