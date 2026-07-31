//! Where does each `(q,t)` family leave `i128`, and does anyone get there?
//!
//! Every `(q,t)` entry point on the Python boundary instantiates at plain
//! `<i128>` with no escalation ladder. Since `[profile.release]` carries
//! `overflow-checks` (`docs/policies/failure.md`, R3) those walls are loud
//! instead of wrong — but R9 wants them *stated*, and a wall nobody has
//! measured is stated as unmeasured. This measures them:
//!
//! ```text
//!   cargo run --release --example probe_qt_walls
//!   BUDGET=1200 cargo run --release --example probe_qt_walls   # seconds per family
//! ```
//!
//! Each family is walked up by degree until one of two things happens: an
//! arithmetic overflow (the wall this is looking for), or the per-family time
//! budget (the *runtime* wall, which for most of these arrives first and is
//! itself the answer — "the i128 wall is past anything anyone can compute" is a
//! reach statement, and a more useful one than a degree nobody will reach).
//!
//! Overflow is caught rather than fatal, so one family's wall does not end the
//! run. That is the only reason `catch_unwind` appears here; a library path
//! would escalate instead.
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::time::{Duration, Instant};

use symfn::Partition;

enum Stop {
    /// The arithmetic wall: this degree overflowed `i128`.
    Overflow(u32),
    /// The runtime wall: the budget ran out first, at this degree.
    Budget(u32, Duration),
}

/// Walk `f` up the degrees from 1 until it overflows or the budget expires.
fn walk(name: &str, budget: Duration, f: impl Fn(u32)) {
    let start = Instant::now();
    let mut last = 0;
    let mut n = 1;
    loop {
        let t = Instant::now();
        let ok = catch_unwind(AssertUnwindSafe(|| f(n))).is_ok();
        let each = t.elapsed();
        if !ok {
            report(name, Stop::Overflow(n), last);
            return;
        }
        last = n;
        // Printed per degree, not only at the end: these curves get steep, so a
        // run stopped by hand or by the machine still leaves its evidence.
        println!(
            "  {name:<28} n = {n:<3} clean  ({:.1}s)",
            each.as_secs_f64()
        );
        // Stop before starting a degree that cannot finish inside the budget:
        // these curves are steeply superlinear, so the next one costs several
        // times this one.
        if start.elapsed() + each * 4 > budget {
            report(name, Stop::Budget(last, start.elapsed()), last);
            return;
        }
        n += 1;
    }
}

fn report(name: &str, stop: Stop, last: u32) {
    match stop {
        Stop::Overflow(n) => println!("{name:<28} overflows at n = {n} (clean through n = {last})"),
        Stop::Budget(n, t) => {
            println!(
                "{name:<28} no overflow through n = {n}, in {:.1}s",
                t.as_secs_f64()
            )
        }
    }
}

fn main() {
    let budget = Duration::from_secs_f64(
        std::env::var("BUDGET")
            .ok()
            .and_then(|s| s.parse().ok())
            .unwrap_or(600.0),
    );
    println!("per-family budget {:.0}s\n", budget.as_secs_f64());

    // Silence the panic printer: an overflow here is the measurement, not a
    // failure, and the backtrace would bury the table.
    let prev = std::panic::take_hook();
    std::panic::set_hook(Box::new(|_| {}));

    walk("hall_littlewood_table", budget, |n| {
        let _ = symfn::hall_littlewood_table::<i128>(n);
    });
    walk("hall_littlewood_p_table", budget, |n| {
        let _ = symfn::hall_littlewood_p_table::<i128>(n);
    });
    walk("kostka_foulkes_table", budget, |n| {
        let _ = symfn::kostka_foulkes_table::<i128>(n);
    });
    walk("qt_kostka_table", budget, |n| {
        let _ = symfn::qt_kostka_table::<i128>(n);
    });
    walk("macdonald_j (every λ ⊢ n)", budget, |n| {
        for lambda in symfn::partitions_of(n) {
            let _ = symfn::macdonald_j::<i128>(&lambda);
        }
    });
    walk("nabla_e", budget, |n| {
        let _ = symfn::nabla_e::<i128>(n);
    });
    walk("delta_prime_e (k = n/2)", budget, |n| {
        let _ = symfn::delta_prime_e::<i128>(n / 2, n);
    });
    walk("llt_h_table (k = 3)", budget, |n| {
        let _ = symfn::llt::llt_h_table::<i128>(n, 3);
    });
    walk("llt_gtilde_table (k = 3)", budget, |n| {
        let _ = symfn::llt::llt_gtilde_table::<i128>(n, 3);
    });
    walk("llt_schur (staircase, k = 3)", budget, |n| {
        let lambda = Partition::new((1..=n).rev());
        let _ = symfn::llt::llt_schur::<i128>(&lambda, 3);
    });

    std::panic::set_hook(prev);
}
