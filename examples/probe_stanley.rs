//! Is there sharing in `stanley`'s transition tree?
//!
//! `stanley` currently walks it as a tree — the very defect §3.1 records for
//! Symmetrica ("nothing is memoized anywhere in the module"). Memoizing is
//! only worth it if distinct nodes are far fewer than visits, so measure that
//! before writing any code. Unlike E3's compression, the per-node value here
//! is a *small* Schur element (1068 terms at S₁₆), so a node-count ratio has
//! a much better chance of converting.

use std::collections::HashMap;
use symfn::permutation::Perm;

fn walk(w: &Perm) -> (u64, usize, u64) {
    let mut visits = 0u64;
    let mut leaves = 0u64;
    let mut seen: HashMap<Perm, u64> = HashMap::new();
    let mut stack = vec![*w];
    while let Some(p) = stack.pop() {
        visits += 1;
        *seen.entry(p).or_insert(0) += 1;
        if p.is_grassmannian().is_some() {
            leaves += 1;
            continue;
        }
        let (_, _, ups) = p.transition().unwrap();
        if ups.is_empty() {
            let m = p.support_len();
            let mut v = vec![1u32];
            for i in 1..=m {
                v.push(p.at(i) + 1);
            }
            stack.push(Perm::new(v).unwrap());
        } else {
            stack.extend(ups);
        }
    }
    (visits, seen.len(), leaves)
}

fn main() {
    let mut seed: u64 = 0xDEAD_BEEF_1234_5678;
    let mut rng = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    println!(
        "{:<16} {:>10} {:>10} {:>10} {:>9}",
        "case", "visits", "distinct", "leaves", "share"
    );
    for n in [8u32, 10, 12, 14, 16, 18] {
        for _ in 0..2 {
            let mut v: Vec<u32> = (1..=n).collect();
            for i in (1..n as usize).rev() {
                let j = (rng() % (i as u64 + 1)) as usize;
                v.swap(i, j);
            }
            let w = Perm::new(v).unwrap();
            let (visits, distinct, leaves) = walk(&w);
            // Rust-side cost, to separate the algorithm from FFI marshalling:
            // the binding hands back one (Vec<u32>, int) pair per Schur term.
            let t = std::time::Instant::now();
            let mut reps = 0;
            let mut terms = 0;
            while t.elapsed().as_secs_f64() < 0.30 {
                let f: symfn::Schur<i128> = symfn::schubert::stanley(&w);
                terms = symfn::SymFn::terms(&f).len();
                reps += 1;
            }
            let per = t.elapsed().as_secs_f64() / reps as f64;
            println!(
                "S_{n} l={:<10} {visits:>10} {distinct:>10} {leaves:>10} {:>8.1}x {terms:>6} terms {:>10.6}s",
                w.length(),
                visits as f64 / distinct as f64,
                per
            );
            if std::env::var("EMIT").is_ok() {
                println!("PERM {:?}", (1..=n).map(|i| w.at(i)).collect::<Vec<_>>());
            }
        }
    }
}
