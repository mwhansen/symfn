//! Calibration for `three_row::prefer_counting`: where does counting overtake
//! the frontier?
//!
//! Every row times [`SkewLr`] against `three_row_product` directly — not
//! through `AutoLr` — so the crossover is visible on both sides of the current
//! dispatch bound. Outputs are asserted equal before anything is reported.
//!
//! Both sides run `REPS` times with the order alternating each repetition, and
//! the minimum per side is reported: an A/B that always runs the same side
//! second hands it a warm allocator, which put the two-row crossover off by a
//! factor of two once (see `two_row::prefer_counting`'s docs).
use std::io::Write;
use std::time::Instant;
use symfn::three_row::{prefer_counting, three_row_product};
use symfn::{LrBackend, Partition, SkewLr};

macro_rules! p {
    ($($a:tt)*) => {{ println!($($a)*); std::io::stdout().flush().ok(); }};
}

const REPS: usize = 4;

fn main() {
    p!(
        "{:<26} {:>4} {:>7} {:>10} {:>10} {:>7}",
        "mu x nu",
        "n",
        "terms",
        "SkewLr",
        "counting",
        "ratio"
    );
    let cases: Vec<(Vec<u32>, Vec<u32>)> = vec![
        // squares of wide three-row shapes, spanning the dispatch bound
        (vec![8, 6, 4], vec![8, 6, 4]),
        (vec![10, 8, 6], vec![10, 8, 6]),
        (vec![12, 10, 8], vec![12, 10, 8]),
        (vec![14, 12, 10], vec![14, 12, 10]),
        (vec![16, 14, 12], vec![16, 14, 12]),
        (vec![20, 16, 12], vec![20, 16, 12]),
        (vec![22, 18, 14], vec![22, 18, 14]),
        (vec![24, 20, 16], vec![24, 20, 16]),
        // asymmetric: nu supplies the strips, mu has 3-5 rows
        (vec![18, 14, 10], vec![9, 7, 5]),
        (vec![16, 13, 10, 7], vec![8, 6, 4]),
        (vec![14, 12, 10, 8, 6], vec![7, 5, 3]),
        (vec![20, 16, 12], vec![10, 8, 6]),
        // shapes the predicate excludes, as the error bar
        (vec![30, 24, 18], vec![3, 2, 1]),
        (vec![6, 5, 4, 3, 2, 1], vec![6, 5, 4]),
    ];
    for (m, n) in cases {
        let mu = Partition::new(m.iter().copied());
        let nu = Partition::new(n.iter().copied());
        let deg = mu.size() + nu.size();

        let (mut t_dp, mut t_c) = (f64::MAX, f64::MAX);
        let (mut want, mut got) = (Vec::new(), None);
        for rep in 0..REPS {
            for side in 0..2 {
                symfn::clear_caches();
                if (rep + side) % 2 == 0 {
                    let t = Instant::now();
                    want = SkewLr.schur_product(&mu, &nu);
                    t_dp = t_dp.min(t.elapsed().as_secs_f64());
                } else {
                    let t = Instant::now();
                    got = three_row_product(&mu, &nu);
                    t_c = t_c.min(t.elapsed().as_secs_f64());
                }
            }
        }
        assert_eq!(got.as_ref(), Some(&want), "s{mu} · s{nu}");
        p!(
            "s{mu}·s{nu} {:>4} {:>7} {:>9.4}s {:>9.4}s {:>6.2}x{}",
            deg,
            want.len(),
            t_dp,
            t_c,
            t_dp / t_c.max(1e-12),
            if prefer_counting(&mu, &nu) {
                "  <-- DISPATCHED"
            } else {
                ""
            }
        );
    }
}
