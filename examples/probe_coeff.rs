//! Does the pruned single-coefficient query buy anything on the cases where
//! the whole product is out of reach?
//!
//! Two measurements:
//!   1. a row that *does* complete (`S_13.2`, 3.24M terms, 5.1s) — query a
//!      coefficient that is genuinely in the support, and compare;
//!   2. `S_13.0`, whose 4.3e16 monomial mass no engine can materialise — ask
//!      for one coefficient anyway.

use std::time::Instant;
use symfn::permutation::Perm;
use symfn::schubert::{schubert_coeff, Schubert};

fn compose(a: &Perm, b: &Perm, n: u32) -> Perm {
    // (a·b)(i) = a(b(i))
    Perm::new((1..=n).map(|i| a.at(b.at(i)))).unwrap()
}

fn main() {
    let s132 = (
        "S_13.2 l=40,41",
        Perm::new([6, 7, 10, 1, 12, 11, 8, 3, 9, 2, 5, 13, 4]).unwrap(),
        Perm::new([2, 12, 9, 11, 4, 10, 6, 3, 1, 8, 7, 5, 13]).unwrap(),
    );
    let s130 = (
        "S_13.0 l=25,36",
        Perm::new([1, 8, 2, 12, 9, 3, 4, 6, 11, 7, 5, 13, 10]).unwrap(),
        Perm::new([1, 12, 7, 3, 2, 11, 8, 6, 13, 10, 5, 9, 4]).unwrap(),
    );

    // --- 1. a computable row: full product vs one coefficient -------------
    let (label, u, v) = &s132;
    let a: Schubert<i128> = Schubert::monomial(*u, 1);
    let b: Schubert<i128> = Schubert::monomial(*v, 1);
    let t = Instant::now();
    let full = a.mul_e2(&b);
    let dt_full = t.elapsed().as_secs_f64();
    // pick a target actually in the support, and one that is not
    let (w_hit, c_hit) = full
        .terms()
        .iter()
        .max_by_key(|(_, c)| **c)
        .map(|(w, c)| (*w, *c))
        .unwrap();
    let t = Instant::now();
    let got: i128 = schubert_coeff(u, v, &w_hit);
    let dt_one = t.elapsed().as_secs_f64();
    assert_eq!(got, c_hit, "pruned query disagreed on {w_hit}");
    println!(
        "{label}: full product {} terms in {dt_full:.3}s; one coefficient \
         (c={c_hit}) in {dt_one:.4}s  -> {:.0}x",
        full.terms().len(),
        dt_full / dt_one.max(1e-9)
    );

    // --- 2. the row nothing can complete ----------------------------------
    let (label, u, v) = &s130;
    let n = 13;
    let w = compose(u, v, n);
    println!(
        "{label}: target w = u·v, l(u)+l(v)={}, l(w)={} -> {}",
        u.length() + v.length(),
        w.length(),
        if w.length() == u.length() + v.length() {
            "lengths add, a valid target"
        } else {
            "lengths do NOT add; c^w = 0 by degree"
        }
    );
    let t = Instant::now();
    let c: i128 = schubert_coeff(u, v, &w);
    println!(
        "{label}: c^w_uv = {c}  in {:.4}s   (zero by degree, no work done)",
        t.elapsed().as_secs_f64()
    );

    // A real target needs l(w) = 61 with u <= w and v <= w. Search for some by
    // random walk from w0 downward -- deterministic LCG, no rand dependency.
    let deg = u.length() + v.length();
    let mut seed: u64 = 0x2545_F491_4F6C_DD1D;
    let mut rng = move || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        seed
    };
    let mut found = 0;
    let mut tried = 0;
    println!("{label}: searching for targets with l(w)={deg}, u<=w, v<=w ...");
    while found < 5 && tried < 400_000 {
        tried += 1;
        let mut one: Vec<u32> = (1..=n).collect();
        for i in (1..n as usize).rev() {
            let j = (rng() % (i as u64 + 1)) as usize;
            one.swap(i, j);
        }
        let w = Perm::new(one).unwrap();
        if w.length() != deg || !u.bruhat_le(&w) || !v.bruhat_le(&w) {
            continue;
        }
        found += 1;
        let t = Instant::now();
        let c: i128 = schubert_coeff(u, v, &w);
        println!("  c^{w}_uv = {c:<6} in {:.4}s", t.elapsed().as_secs_f64());
    }
    if found == 0 {
        println!("  none found in {tried} random draws");
    }
}
