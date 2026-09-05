//! Why is one input pair pathological when a *longer* one is easy?
//!
//! `S_13 ℓ=25,36` defeats both the C `schubmult` (>120s) and E2 (20+ min),
//! while `S_13 ℓ=40,41` — longer, and with 3.4× the output — finishes in 5s.
//! Since both engines fail the same row, the difficulty is a property of the
//! input pair rather than of either implementation, so it should be visible
//! *without* running a product.
//!
//! Everything printed here is computed from permutations alone:
//! `dimension` is `S_w(1,…,1)`, the pipe-dream count, and `transition_tree`
//! walks E2's own recursion without building any element. Both are
//! microseconds even where the product is hopeless.

use symfn::permutation::Perm;
use symfn::schubert::{dimension, transition_tree};

fn stair(k: u32) -> Vec<u32> {
    (1..=k)
        .map(|i| 2 * i)
        .chain((1..=k).map(|i| 2 * i - 1))
        .collect()
}

fn main() {
    let cases: Vec<(&str, Vec<u32>, Vec<u32>, &str)> = vec![
        ("stair6^2", stair(6), stair(6), "0.070s"),
        ("stair7^2", stair(7), stair(7), "1.51s"),
        (
            "S_11.1 l=35,26",
            vec![4, 8, 5, 11, 10, 6, 9, 1, 7, 3, 2],
            vec![7, 5, 3, 2, 9, 6, 11, 10, 1, 4, 8],
            "0.40s",
        ),
        (
            "S_12.0 l=37,24",
            vec![12, 1, 4, 6, 11, 3, 10, 8, 9, 5, 7, 2],
            vec![3, 2, 5, 4, 12, 6, 9, 11, 8, 7, 1, 10],
            "6.58s",
        ),
        (
            "S_12.2 l=43,32",
            vec![7, 12, 10, 3, 9, 5, 2, 8, 4, 11, 6, 1],
            vec![3, 4, 12, 2, 10, 11, 6, 5, 8, 1, 9, 7],
            "0.63s",
        ),
        (
            "S_13.1 l=36,41",
            vec![9, 13, 3, 7, 2, 1, 4, 10, 11, 5, 12, 8, 6],
            vec![9, 5, 8, 10, 2, 12, 11, 4, 3, 6, 1, 13, 7],
            "1.90s",
        ),
        (
            "S_13.2 l=40,41",
            vec![6, 7, 10, 1, 12, 11, 8, 3, 9, 2, 5, 13, 4],
            vec![2, 12, 9, 11, 4, 10, 6, 3, 1, 8, 7, 5, 13],
            "5.10s",
        ),
        (
            "S_13.0 l=25,36 *",
            vec![1, 8, 2, 12, 9, 3, 4, 6, 11, 7, 5, 13, 10],
            vec![1, 12, 7, 3, 2, 11, 8, 6, 13, 10, 5, 9, 4],
            "NEITHER",
        ),
    ];
    println!(
        "{:<17} {:>5} {:>16} {:>16} {:>9} {:>9} {:>9}",
        "case", "deg", "dim(u)", "dim(v)", "tree(u)", "tree(v)", "E2 time"
    );
    for (label, u, v, t) in cases {
        let (pu, pv) = (Perm::new(u).unwrap(), Perm::new(v).unwrap());
        let (nu, _) = transition_tree(&pu);
        let (nv, _) = transition_tree(&pv);
        println!(
            "{:<17} {:>5} {:>16} {:>16} {:>9} {:>9} {:>9}",
            label,
            pu.length() + pv.length(),
            dimension(&pu).map_or(String::from("-"), |d| d.to_string()),
            dimension(&pv).map_or(String::from("-"), |d| d.to_string()),
            nu,
            nv,
            t
        );
    }
}
