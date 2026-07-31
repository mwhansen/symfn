//! Where does each `(q,t)` family leave `i128`, and how fast is it heading
//! there?
//!
//! Every `(q,t)` entry point on the Python boundary instantiates at plain
//! `<i128>` with no escalation ladder. Since `[profile.release]` carries
//! `overflow-checks` (`docs/policies/failure.md`, R3) those walls are loud
//! instead of wrong — but R9 wants them *stated*, and a wall nobody has
//! measured is stated as unmeasured.
//!
//! ```text
//!   cargo run --release --example probe_qt_walls
//!   BUDGET=1200 cargo run --release --example probe_qt_walls   # seconds per family
//! ```
//!
//! **The largest coefficient is the measurement, not the panic.** Walking the
//! degree up until something overflows answers "did it overflow by n = N?",
//! which is not a reach statement: it says nothing about whether N+1 is fine or
//! whether the family was one degree from the wall the whole time. So each
//! degree reports the **bit width of its largest coefficient**, and the growth
//! in that width per degree is what extrapolates to the wall — the form
//! `memo.rs` already uses for its own bounds ("25 bits at degree 12, growing
//! about 3 per degree, so `i128` holds past degree 45").
//!
//! Two walls are in play and the run distinguishes them:
//!
//! * the **arithmetic** wall, 127 bits, extrapolated from the measured slope;
//! * the **runtime** wall, where the enumeration stops finishing — the budget
//!   below. For most of these families it arrives first, and that is itself the
//!   reach statement a researcher needs.
//!
//! The bit widths are of the *answers*. Intermediates are covered by the
//! profile flag rather than by inspection: with `overflow-checks` on, an
//! intermediate that left `i128` panics, so a degree reported clean is one
//! where nothing overflowed anywhere — which is why the panic is caught rather
//! than fatal here. A library path would escalate instead of catching.
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::time::{Duration, Instant};

use symfn::{Frac, Monomial, Partition, QtPoly, Schur, SymFn};

/// The widest coefficient in a result, in bits — `0` for an empty or zero one.
///
/// Bits rather than the value itself because that is what extrapolates: these
/// families grow geometrically, so the width is roughly linear in the degree
/// and a slope in bits/degree is a number a reader can carry forward.
trait Widest {
    fn bits(&self) -> u32;
}

impl Widest for i128 {
    fn bits(&self) -> u32 {
        128 - self.unsigned_abs().leading_zeros()
    }
}

impl Widest for QtPoly<i128> {
    fn bits(&self) -> u32 {
        self.terms().map(|(_, c)| c.bits()).max().unwrap_or(0)
    }
}

impl Widest for Frac<i128> {
    /// Numerator only: the denominator of a Macdonald `J` is a *factored*
    /// multiset of binomials `1 − qᵃtᵇ`, never expanded, so it carries no large
    /// integer to measure.
    fn bits(&self) -> u32 {
        self.parts().0.bits()
    }
}

impl<T: Widest> Widest for Vec<T> {
    fn bits(&self) -> u32 {
        self.iter().map(Widest::bits).max().unwrap_or(0)
    }
}

impl<T: Widest> Widest for (Partition, T) {
    fn bits(&self) -> u32 {
        self.1.bits()
    }
}

impl<T: Widest> Widest for (Vec<u32>, T) {
    fn bits(&self) -> u32 {
        self.1.bits()
    }
}

impl<C: Widest + symfn::Ring> Widest for Schur<C> {
    fn bits(&self) -> u32 {
        self.terms().values().map(Widest::bits).max().unwrap_or(0)
    }
}

impl<C: Widest + symfn::Ring> Widest for Monomial<C> {
    fn bits(&self) -> u32 {
        self.terms().values().map(Widest::bits).max().unwrap_or(0)
    }
}

/// One family's curve: `(degree, widest coefficient in bits, seconds)`.
struct Curve {
    name: &'static str,
    rows: Vec<(u32, u32, f64)>,
    overflowed_at: Option<u32>,
}

impl Curve {
    /// Bits gained per degree, over the last third of the walk.
    ///
    /// The tail rather than the whole curve: the first few degrees are
    /// dominated by the shape's own smallness, and it is the asymptotic slope
    /// that extrapolates.
    fn slope(&self) -> Option<f64> {
        let n = self.rows.len();
        if n < 6 {
            return None;
        }
        let tail = &self.rows[n - n / 3..];
        let (a, b) = (tail.first()?, tail.last()?);
        if b.0 == a.0 {
            return None;
        }
        Some(f64::from(b.1 - a.1) / f64::from(b.0 - a.0))
    }

    /// The degree at which the widest coefficient would reach 127 bits, if the
    /// tail slope held. An extrapolation, and labelled as one wherever it is
    /// quoted.
    fn projected_wall(&self) -> Option<u32> {
        let s = self.slope()?;
        if s <= 0.0 {
            return None;
        }
        let last = self.rows.last()?;
        Some(last.0 + ((127.0 - f64::from(last.1)) / s).ceil().max(0.0) as u32)
    }

    fn report(&self) {
        let last = self.rows.last();
        print!("{:<30}", self.name);
        match (self.overflowed_at, last) {
            (Some(n), _) => print!(" overflows at n = {n:<3}"),
            (None, Some(l)) => print!(" clean to n = {:<3} ({} bits)", l.0, l.1),
            (None, None) => print!(" nothing measured"),
        }
        match (self.slope(), self.projected_wall()) {
            (Some(s), Some(w)) => println!("  +{s:.1} bits/degree -> 127 bits near n = {w}"),
            (Some(s), None) => println!("  +{s:.1} bits/degree"),
            _ => println!("  (too few degrees to slope)"),
        }
    }
}

/// The degree to stop at even when nothing else stops the walk.
///
/// Needed, not belt-and-braces: a shape family can be genuinely *cheap* at
/// every degree — `Q'_{(n-1,1)}` is, because the Morris recursion peels one
/// part and the expansion stays tiny — so a time budget alone never trips. The
/// first run of this probe reached n = 85 799 505 that way, having long since
/// stopped measuring anything a caller would ask for.
const MAX_DEGREE: u32 = 200;

/// Walk `f` up the degrees until it overflows, the budget expires, or
/// [`MAX_DEGREE`].
fn walk(name: &'static str, budget: Duration, f: impl Fn(u32) -> u32) -> Curve {
    let start = Instant::now();
    let mut curve = Curve {
        name,
        rows: Vec::new(),
        overflowed_at: None,
    };
    let mut n = 1;
    loop {
        let t = Instant::now();
        let bits = catch_unwind(AssertUnwindSafe(|| f(n)));
        let each = t.elapsed();
        let Ok(bits) = bits else {
            curve.overflowed_at = Some(n);
            return curve;
        };
        // Printed per degree, not only at the end: these curves get steep, so a
        // run stopped by hand or by the machine still leaves its evidence.
        println!(
            "  {name:<30} n = {n:<3} {bits:>3} bits  ({:.1}s)",
            each.as_secs_f64()
        );
        curve.rows.push((n, bits, each.as_secs_f64()));
        // Stop before starting a degree that cannot finish inside the budget:
        // these curves are steeply superlinear, so the next one costs several
        // times this one.
        if start.elapsed() + each * 4 > budget || n >= MAX_DEGREE {
            return curve;
        }
        n += 1;
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

    // Which *unit of work* to walk. The two have different runtime walls and
    // almost the same arithmetic one, so they are different reach statements
    // about the same family: `hall_littlewood_table(n)` is p(n) polynomials,
    // while a researcher asking for one `Q'_λ` reaches much further in degree
    // before the enumeration stops finishing — and can therefore meet the
    // coefficient wall first, which the table never does.
    //
    // `kostka_foulkes` and `qt_kostka_column` deliberately have no single-shape
    // row: both document that one value costs the whole column or the whole
    // table, so the table *is* their unit whatever the caller asks for.
    let unit = std::env::var("UNIT").unwrap_or_else(|_| "table".into());
    let mut curves = Vec::new();

    if unit != "table" {
        // λ = 1ⁿ is where the √(n!) bound on Kostka–Foulkes coefficients is
        // attained (`K_{μ,1ⁿ} = f^μ`), so it is the worst case for growth as
        // well as a natural single-shape ask.
        curves.push(walk("hall_littlewood (λ = 1^n)", budget, |n| {
            symfn::hall_littlewood::<i128>(&Partition::new(std::iter::repeat_n(1, n as usize)))
                .bits()
        }));
        curves.push(walk("hall_littlewood (λ = (n-1,1))", budget, |n| {
            let l = if n < 2 {
                Partition::new([n])
            } else {
                Partition::new([n - 1, 1])
            };
            symfn::hall_littlewood::<i128>(&l).bits()
        }));
        curves.push(walk("macdonald_j (λ = 1^n)", budget, |n| {
            symfn::macdonald_j::<i128>(&Partition::new(std::iter::repeat_n(1, n as usize))).bits()
        }));
        curves.push(walk("llt_h (μ = 1^n, k = 3)", budget, |n| {
            symfn::llt::llt_h::<i128>(&Partition::new(std::iter::repeat_n(1, n as usize)), 3).bits()
        }));
    }

    if unit == "table" || unit == "both" {
        curves.extend([
            walk("hall_littlewood_table", budget, |n| {
                symfn::hall_littlewood_table::<i128>(n).bits()
            }),
            walk("hall_littlewood_p_table", budget, |n| {
                symfn::hall_littlewood_p_table::<i128>(n).bits()
            }),
            walk("kostka_foulkes_table", budget, |n| {
                symfn::kostka_foulkes_table::<i128>(n).bits()
            }),
            walk("qt_kostka_table", budget, |n| {
                symfn::qt_kostka_table::<i128>(n).bits()
            }),
            walk("macdonald_j (every lambda)", budget, |n| {
                symfn::partitions_of(n)
                    .iter()
                    .map(|l| symfn::macdonald_j::<i128>(l).bits())
                    .max()
                    .unwrap_or(0)
            }),
            walk("nabla_e", budget, |n| symfn::nabla_e::<i128>(n).bits()),
            walk("delta_prime_e (k = n/2)", budget, |n| {
                symfn::delta_prime_e::<i128>(n / 2, n).bits()
            }),
            walk("llt_h_table (k = 3)", budget, |n| {
                symfn::llt::llt_h_table::<i128>(n, 3).bits()
            }),
            walk("llt_gtilde_table (k = 3)", budget, |n| {
                symfn::llt::llt_gtilde_table::<i128>(n, 3).bits()
            }),
            // No staircase row here. `λ = (n, n-1, …, 1)` has size `n(n+1)/2`,
            // which carries no 3-ribbon tableaux unless 3 divides it, so two
            // degrees in three return zero instantly. That fits no slope — and
            // it defeats the budget's lookahead as well, since the walk read
            // n = 13 as free (it is vacuous) and then spent over an hour inside
            // n = 14, which is not. A row that alternates empty and enormous
            // measures neither.
            walk("llt_h (μ = (n), k = 3)", budget, |n| {
                symfn::llt::llt_h::<i128>(&Partition::new([n]), 3).bits()
            }),
        ]);
    }

    std::panic::set_hook(prev);

    println!("\n== summary ==");
    for c in &curves {
        c.report();
    }
    println!(
        "\nBit widths are of the answers; intermediates are covered by\n\
         overflow-checks, which panics on one that leaves i128 -- so a degree\n\
         reported clean had no overflow anywhere. The projected wall is an\n\
         extrapolation of the tail slope, not a measured degree."
    );
}
