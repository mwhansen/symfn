//! The catalog of measured workloads, and their memory budgets.
//!
//! **One entry here gives you both** the exploratory report
//! (`cargo run --release --example heapstat -- <name>`) and a regression test
//! (`cargo test --release --test memory`). Adding a subsystem to the memory
//! story is a single `Workload` — that is why the catalog exists
//! rather than each harness carrying its own list.
//!
//! ## Adding one
//!
//! 1. Write the entry with `peak: 0, allocs: 0` and any tolerance.
//! 2. `cargo run --release --example heapstat -- <name>` — it prints the
//!    measured numbers and the `Budget` line to paste back.
//! 3. Round the numbers *up* a little and paste. A budget is a ceiling with
//!    headroom, not a golden value.
//!
//! Size the workload where the work actually lives. A degree that finishes in
//! milliseconds usually measures allocator warm-up rather than the algorithm,
//! and `docs/record/oracles-and-comparisons.md` carries the comparison ladder
//! where sizing at a real input moved the answer by three orders of magnitude.

use super::Budget;
use crate::partition::Partition;

/// A named piece of work, its budget, and how to run it.
pub struct Workload {
    pub name: &'static str,
    /// Runs the work and returns a short note for the report (term counts and
    /// the like), which doubles as a check that the work was not optimized away.
    pub run: fn() -> String,
    pub budget: Budget,
}

fn p(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

/// The `μ²` skew shape the LR benchmarks use: outer = (μ+w) ++ μ, inner = wˡ.
fn square(mu: &[u32]) -> (Partition, Partition) {
    let mu = p(mu);
    let w = mu.part(0);
    let outer: Vec<u32> = mu
        .parts()
        .iter()
        .map(|x| x + w)
        .chain(mu.parts().iter().copied())
        .collect();
    let inner: Vec<u32> = mu.parts().iter().map(|_| w).collect();
    (p(&outer), p(&inner))
}

/// Every workload, in rough order of size.
///
/// Budgets were measured on the numbers recorded in
/// `docs/record/memory.md`. Single-threaded workloads repeat bit-exactly and
/// carry a tight tolerance; `skew-big` fans out across threads, so its shard
/// split — and with it the allocation count — follows thread scheduling.
pub const WORKLOADS: &[Workload] = &[
    Workload {
        name: "skew-mid",
        run: || {
            let (o, i) = square(&[10, 8, 6, 4]);
            format!("{} terms", crate::expand_skew(&o, &i).len())
        },
        budget: Budget {
            name: "skew-mid",
            peak: 4_400_000,
            allocs: 50_000,
            tolerance: 0.10,
        },
    },
    Workload {
        name: "skew-big",
        run: || {
            let (o, i) = square(&[8, 7, 6, 5, 4, 3]);
            format!("{} terms", crate::expand_skew(&o, &i).len())
        },
        budget: Budget {
            name: "skew-big",
            peak: 88_000_000,
            allocs: 340_000,
            tolerance: 0.10,
        },
    },
    Workload {
        // What `expand_skew` costs to *return* an expansion the cache already
        // holds. Guards the `expand_skew_shared` split: if the owned version
        // ever becomes the internal default again, this is where it shows.
        name: "skew-clone",
        run: || {
            let (o, i) = square(&[8, 7, 6, 5, 4, 3]);
            let warm = crate::skew_lr::expand_skew_shared(&o, &i);
            let n = warm.len();
            drop(warm);
            super::reset();
            format!("{} terms, {} returned", n, crate::expand_skew(&o, &i).len())
        },
        budget: Budget {
            name: "skew-clone",
            peak: 16_000_000,
            allocs: 170_000,
            tolerance: 0.10,
        },
    },
    Workload {
        name: "product",
        run: || {
            use crate::LrBackend;
            let b = crate::AutoLr;
            format!(
                "{} terms",
                b.schur_product(&p(&[9, 7, 5, 3, 1]), &p(&[9, 7, 5, 3, 1]))
                    .len()
            )
        },
        budget: Budget {
            name: "product",
            peak: 8_000_000,
            allocs: 58_000,
            tolerance: 0.10,
        },
    },
    Workload {
        // What `Schur::mul` costs on top of an expansion the cache already
        // holds — the consumer side of a product. Guards the shared read in
        // `mul_with`: if the expansion is ever copied whole again, the peak
        // doubles here.
        name: "schur-mul",
        run: || {
            use crate::{LrBackend, SymFn};
            let mu = p(&[8, 7, 6, 5, 4, 3]);
            let n = crate::AutoLr.schur_product(&mu, &mu).len();
            super::reset();
            let s: crate::Schur<i64> = crate::Schur::monomial(mu, 1);
            format!("{} terms, {} in the product", s.mul(&s).terms().len(), n)
        },
        budget: Budget {
            name: "schur-mul",
            peak: 18_000_000,
            allocs: 180_000,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "coproduct",
        run: || {
            use crate::SymFn;
            let s: crate::Schur<i64> = crate::Schur::monomial(p(&[7, 6, 5, 4, 3]), 1);
            format!("{} terms", crate::coproduct(&s).terms().len())
        },
        budget: Budget {
            name: "coproduct",
            peak: 2_500_000,
            allocs: 95_000,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "htilde",
        run: || format!("{} rows", crate::qt_kostka_table::<i128>(10).len()),
        budget: Budget {
            name: "htilde",
            peak: 18_000_000,
            allocs: 760_000,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "s-in-j",
        run: || {
            format!(
                "{} rows",
                crate::schur_in_j_table::<crate::Rational>(9).len()
            )
        },
        budget: Budget {
            name: "s-in-j",
            peak: 10_500_000,
            allocs: 580_000,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "m-in-p",
        run: || {
            let m: crate::Monomial<crate::Frac<crate::Rational>> = crate::sym::SymFn::monomial(
                crate::Partition::new([5, 3, 1]),
                <crate::Frac<crate::Rational> as crate::Ring>::one(),
            );
            format!("{} terms", crate::monomial_to_macdonald_p(&m).len())
        },
        budget: Budget {
            name: "m-in-p",
            peak: 7_350_000,
            allocs: 624_100,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "m-in-jack-p",
        run: || {
            let m: crate::Monomial<crate::AFrac<crate::Rational>> = crate::sym::SymFn::monomial(
                crate::Partition::new([5, 3, 1]),
                <crate::AFrac<crate::Rational> as crate::Ring>::one(),
            );
            format!("{} terms", crate::monomial_to_jack_p(&m).len())
        },
        budget: Budget {
            name: "m-in-jack-p",
            peak: 520_000,
            allocs: 72_000,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "hl",
        run: || format!("{} rows", crate::hall_littlewood_table::<i64>(12).len()),
        budget: Budget {
            name: "hl",
            peak: 2_000_000,
            allocs: 38_000,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "llt",
        run: || format!("{} rows", crate::llt_h_table::<i64>(9, 3).len()),
        budget: Budget {
            name: "llt",
            peak: 400_000,
            allocs: 85_000,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "jack",
        run: || format!("{} rows", crate::jack_table::<i64>(9).len()),
        budget: Budget {
            name: "jack",
            peak: 400_000,
            allocs: 42_000,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "kostka-foulkes",
        run: || format!("{} rows", crate::kostka_foulkes_table::<i64>(12).len()),
        budget: Budget {
            name: "kostka-foulkes",
            peak: 2_000_000,
            allocs: 41_000,
            tolerance: 0.05,
        },
    },
    Workload {
        name: "character",
        run: || format!("{} rows", crate::character::character_table(24).len()),
        budget: Budget {
            name: "character",
            peak: 42_000_000,
            allocs: 13_000,
            tolerance: 0.05,
        },
    },
    Workload {
        // A session under a budget: the sweep `session_sweep` describes, held
        // to SESSION_BUDGET, checking at the end that the caches sit under
        // the budget plus tier 0 — the invariant every insert restores when
        // no other thread holds a table. The peak below is what the sweep
        // costs *while* bounded, so a budget that stopped biting would show
        // here as a peak that grew.
        name: "session",
        run: || {
            const SESSION_BUDGET: usize = 4 << 20;
            crate::set_cache_budget(Some(SESSION_BUDGET));
            session_sweep((16, 12, 8, 7), |_| {});
            let stats = crate::cache_stats();
            let held: usize = stats.iter().map(|r| r.bytes).sum();
            let tier0: usize = stats.iter().filter(|r| r.tier == 0).map(|r| r.bytes).sum();
            crate::set_cache_budget(None);
            assert!(
                held <= SESSION_BUDGET + tier0,
                "{held} bytes held against a budget of {SESSION_BUDGET} plus {tier0} of tier 0"
            );
            format!("{held} bytes held under a {SESSION_BUDGET}-byte budget")
        },
        budget: Budget {
            name: "session",
            peak: 8_400_000,
            allocs: 1_200_000,
            tolerance: 0.10,
        },
    },
    Workload {
        name: "schubert",
        run: || {
            let w = crate::permutation::Perm::new((1..=9u32).rev()).expect("w0 is a permutation");
            let a = crate::schubert::Schubert::<i128>::monomial(w, 1);
            format!("{} terms", a.mul(&a).terms().len())
        },
        budget: Budget {
            name: "schubert",
            peak: 200_000,
            allocs: 2_000,
            tolerance: 0.05,
        },
    },
];

/// A session rather than a benchmark: rising degrees, and at each one the
/// calls a research session makes, one value at a time — every character and
/// Kostka number of the degree, every coefficient of every Schur square, every
/// (q,t)-Kostka number, and each Schur function written in Macdonald `J` —
/// with `after(degree)` called between degrees so a census can read the
/// tables.
///
/// One value at a time is the point. Each of those calls reads a table the
/// crate builds whole — the character memo, a product expansion, the `H̃`
/// table of the degree, the `s → J` matrix — so what a session gets from
/// the caches is the difference between reading that table p(n) times and
/// building it p(n) times, which is what a budget below the working set
/// costs (`docs/record/qt-kostka.md`, `docs/record/macdonald.md`).
///
/// The four ceilings are the last degree at which each family is asked for;
/// the families are the ones whose tables dominate a session, one per tier
/// the budget distinguishes.
pub fn session_sweep(ceilings: (u32, u32, u32, u32), mut after: impl FnMut(u32)) {
    use crate::coeff::Ring;
    use crate::sym::SymFn;
    let (chars, kostka, products, transition) = ceilings;
    let top = chars.max(kostka).max(products).max(transition);
    for n in 1..=top {
        let parts = crate::partitions_of(n);
        if n <= chars {
            for a in &parts {
                for b in &parts {
                    std::hint::black_box(crate::character(a, b));
                }
            }
        }
        if n <= kostka {
            for a in &parts {
                for b in &parts {
                    std::hint::black_box(crate::kostka(a, b));
                }
            }
        }
        if n <= products {
            use crate::LrBackend;
            let doubled = crate::partitions_of(2 * n);
            for a in &parts {
                for c in &doubled {
                    std::hint::black_box(crate::AutoLr.lr_coeff(c, a, a));
                }
            }
        }
        if n <= transition {
            for a in &parts {
                for b in &parts {
                    std::hint::black_box(crate::qt_kostka::<i128>(a, b));
                }
                let s: crate::Schur<crate::QtPoly<crate::Rational>> =
                    crate::Schur::monomial(a.clone(), crate::QtPoly::one());
                std::hint::black_box(crate::schur_to_macdonald_j(&s).len());
            }
        }
        after(n);
    }
}

/// Look a workload up by name.
pub fn find(name: &str) -> Option<&'static Workload> {
    WORKLOADS.iter().find(|w| w.name == name)
}

/// Every workload name, for a usage message.
pub fn names() -> Vec<&'static str> {
    WORKLOADS.iter().map(|w| w.name).collect()
}
