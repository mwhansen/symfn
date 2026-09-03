//! The bytes the caches report against the bytes the allocator releases when
//! they are cleared.
//!
//! ```text
//!   cargo test --release --test cache_accounting -- --nocapture
//! ```
//!
//! `cache_stats` is the crate's own arithmetic over capacities and `size_of`
//! (`src/memo.rs`, `HeapSize`). What it cannot see is the allocator's rounding
//! of each block to a size class, so the check is a band rather than an
//! equality, and the band is what that rounding measured at
//! (`docs/record/memory.md`). A reading outside it means the arithmetic has
//! drifted from the types it describes.
//!
//! Its own test binary, because the counting allocator is process-wide and
//! `tests/memory.rs` runs its workloads under the same counters.
//!
//! **Release only**, for the same reason as `tests/memory.rs`.

#[global_allocator]
static ALLOC: symfn::measure::Counting = symfn::measure::Counting::new();

use symfn::measure;
use symfn::SymFn;

fn counted() -> usize {
    symfn::cache_stats().iter().map(|r| r.bytes).sum()
}

/// What the caches held, measured as the drop in live bytes when they are
/// cleared — which isolates the tables from anything else a call leaves
/// behind, such as a thread's retained scratch.
fn held_by_caches() -> usize {
    let before = measure::live();
    symfn::clear_caches();
    before.saturating_sub(measure::live())
}

fn check(what: &str, held: usize, counted: usize) -> Result<(), String> {
    let ratio = counted as f64 / held.max(1) as f64;
    println!("{what:<16} held {held:>10}   counted {counted:>10}   ratio {ratio:.3}");
    if (0.95..=1.05).contains(&ratio) {
        Ok(())
    } else {
        Err(format!(
            "{what}: counted {counted} bytes against {held} released by clear_caches, ratio {ratio:.3}"
        ))
    }
}

/// Over the two table shapes that dominate a session: one `Arc` to a large
/// expansion per entry (`skews`), and many small fixed-width entries
/// (`character_masks`). The first half — that every row reads empty after
/// `clear_caches` — holds in both profiles and lives here rather than in
/// `memo.rs`'s unit tests because those share the tables across threads.
#[test]
fn cache_accounting_tracks_the_allocator() {
    let p = |v: &[u32]| symfn::Partition::new(v.iter().copied());
    let _ = symfn::kostka(&p(&[3, 1]), &p(&[2, 1, 1]));
    let _ = symfn::character(&p(&[3, 1]), &p(&[2, 1, 1]));
    let _: symfn::Elementary<i64> =
        symfn::convert(&symfn::Homogeneous::<i64>::monomial(p(&[3, 1]), 1));
    let filled: Vec<&str> = symfn::cache_stats()
        .iter()
        .filter(|r| r.bytes > 0)
        .map(|r| r.name)
        .collect();
    assert!(
        ["kostka", "character_masks", "flip_rows"]
            .iter()
            .all(|n| filled.contains(n)),
        "the calls above fill three named tables; filled: {filled:?}"
    );
    symfn::clear_caches();
    for row in symfn::cache_stats() {
        assert_eq!(
            (row.entries, row.bytes),
            (0, 0),
            "{} after clear_caches",
            row.name
        );
    }

    if cfg!(debug_assertions) {
        eprintln!("skipped: the accounting is calibrated for --release");
        return;
    }
    let mut failures = Vec::new();

    symfn::clear_caches();
    measure::reset();
    let outer = p(&[16, 15, 14, 13, 12, 11, 8, 7, 6, 5, 4, 3]);
    let inner = p(&[8, 8, 8, 8, 8, 8]);
    let terms = symfn::skew_lr::expand_skew_shared(&outer, &inner).len();
    assert!(
        terms > 100_000,
        "the skew expansion is what fills the table: {terms} terms"
    );
    for row in symfn::cache_stats().iter().filter(|r| r.bytes > 0) {
        println!(
            "  {:<22} {:>8} entries {:>10} bytes",
            row.name, row.entries, row.bytes
        );
    }
    let counted_skews = counted();
    if let Err(e) = check("skews", held_by_caches(), counted_skews) {
        failures.push(e);
    }

    measure::reset();
    let parts = symfn::partitions_of(16);
    let mut nonzero = 0usize;
    for lam in &parts {
        for mu in &parts {
            if symfn::character(lam, mu) != 0 {
                nonzero += 1;
            }
        }
    }
    assert!(
        nonzero > 10_000,
        "the character sweep is what fills the table"
    );
    for row in symfn::cache_stats().iter().filter(|r| r.bytes > 0) {
        println!(
            "  {:<22} {:>8} entries {:>10} bytes",
            row.name, row.entries, row.bytes
        );
    }
    let counted_masks = counted();
    if let Err(e) = check("character_masks", held_by_caches(), counted_masks) {
        failures.push(e);
    }

    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
