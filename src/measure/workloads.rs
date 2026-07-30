//! The catalogue of measured workloads, and their memory budgets.
//!
//! **One entry here gives you both** the exploratory report
//! (`cargo run --release --example heapstat -- <name>`) and a regression test
//! (`cargo test --release --test memory`). Adding a subsystem to the memory
//! story is a single `Workload` — that is the whole point of the catalogue
//! existing rather than each harness carrying its own list.
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
//! milliseconds usually measures allocator warm-up rather than the algorithm —
//! the same mistake the Sage comparison ladder made before it was resized, where
//! picking a real-sized case moved the answer by three orders of magnitude.

use super::Budget;
use crate::partition::Partition;

/// A named piece of work, its budget, and how to run it.
pub struct Workload {
    pub name: &'static str,
    /// Runs the work and returns a short note for the report (term counts and
    /// the like), which doubles as a check that the work was not optimised away.
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
            allocs: 510_000,
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
            let b = crate::AutoLr::default();
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
        name: "schubert",
        run: || {
            let w = crate::permutation::Perm::new((1..=9u32).rev()).unwrap();
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

/// Look a workload up by name.
pub fn find(name: &str) -> Option<&'static Workload> {
    WORKLOADS.iter().find(|w| w.name == name)
}

/// Every workload name, for a usage message.
pub fn names() -> Vec<&'static str> {
    WORKLOADS.iter().map(|w| w.name).collect()
}
