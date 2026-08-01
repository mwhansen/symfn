//! Calibration for `two_row::prefer_counting`: where does counting overtake the
//! layer?
//!
//! Every row runs `AutoLr` (which dispatches) against `SkewLr` (which does not),
//! and asserts they agree. Rows marked DISPATCHED are the ones the predicate
//! fires on; the rest route to `SkewLr` on both sides, so their ratios are pure
//! measurement noise and serve as the error bar — roughly ±20% below a
//! millisecond, much tighter above it.
use std::io::Write;
use std::time::Instant;
use symfn::two_row::prefer_counting;
use symfn::{AutoLr, LrBackend, Partition, SkewLr};

macro_rules! p {
    ($($a:tt)*) => {{ println!($($a)*); std::io::stdout().flush().ok(); }};
}

fn main() {
    p!(
        "{:<30} {:>5} {:>8} {:>10} {:>10} {:>8}",
        "mu x nu",
        "n",
        "terms",
        "SkewLr",
        "counting",
        "ratio"
    );
    let cases: Vec<(Vec<u32>, Vec<u32>)> = vec![
        // 3-row mu, proportional nu
        (vec![10, 8, 6], vec![10, 8]),
        (vec![20, 16, 12], vec![20, 16]),
        (vec![28, 22, 17], vec![28, 22]),
        (vec![34, 27, 20], vec![34, 27]),
        (vec![40, 32, 24], vec![40, 32]),
        (vec![50, 40, 30], vec![50, 40]),
        // 2-row mu
        (vec![20, 12], vec![20, 12]),
        (vec![40, 24], vec![40, 24]),
        (vec![70, 42], vec![70, 42]),
        (vec![110, 66], vec![110, 66]),
        // 4-row mu
        (vec![16, 13, 10, 7], vec![16, 13]),
        (vec![24, 20, 16, 12], vec![24, 20]),
        (vec![34, 28, 22, 16], vec![34, 28]),
        // 1-row mu (Pieri-like)
        (vec![60], vec![30, 20]),
        (vec![160], vec![80, 50]),
        // lopsided nu
        (vec![30, 24, 18], vec![40, 2]),
        (vec![30, 24, 18], vec![6, 5]),
        // tall mu: `SkewLr`'s best regime, so counting should lose
        (vec![20, 18, 16, 14, 12, 10, 8, 6], vec![20, 16]),
        (vec![14, 13, 12, 11, 10, 9, 8, 7, 6, 5], vec![14, 11]),
        (vec![12, 11, 10, 9, 8, 7, 6], vec![12, 9]),
        // wide-ish mu with 3 rows near the boundary
        (vec![24, 19, 14], vec![24, 19]),
        (vec![30, 24, 18], vec![30, 24]),
        // nu much larger than mu
        (vec![8, 6, 4], vec![40, 32]),
        (vec![12, 9, 6], vec![60, 48]),
    ];
    for (m, n) in cases {
        let mu = Partition::new(m.iter().copied());
        let nu = Partition::new(n.iter().copied());
        let deg = mu.size() + nu.size();

        symfn::clear_caches();
        let t = Instant::now();
        let want = SkewLr.schur_product(&mu, &nu);
        let t_dp = t.elapsed();

        symfn::clear_caches();
        let t = Instant::now();
        let got = AutoLr.schur_product(&mu, &nu);
        let t_c = t.elapsed();

        assert_eq!(got, want, "s{mu} · s{nu}");
        let r = t_dp.as_secs_f64() / t_c.as_secs_f64().max(1e-12);
        p!(
            "s{mu}·s{nu}{:<4} {:>5} {:>8} {:>9.4?} {:>9.4?} {:>7.2}x{}",
            "",
            deg,
            want.len(),
            t_dp,
            t_c,
            r,
            if prefer_counting(&mu, &nu) {
                "  <-- DISPATCHED"
            } else {
                ""
            }
        );
    }
}
