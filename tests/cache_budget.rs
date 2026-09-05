//! The byte budget on the memo caches: what it clears, in what order, and
//! that clearing changes no result.
//!
//! Its own binary because the budget is process-wide and every other test
//! binary runs its tests on several threads; here the one test sets, uses and
//! lifts the budget alone.

use std::collections::BTreeMap;

use symfn::{
    cache_budget, cache_stats, clear_caches, partitions_of, set_cache_budget, LrBackend, Partition,
    SymFn,
};

fn held() -> usize {
    cache_stats().iter().map(|r| r.bytes).sum()
}

fn row(name: &str) -> (usize, usize) {
    let r = cache_stats()
        .into_iter()
        .find(|r| r.name == name)
        .unwrap_or_else(|| panic!("no cache row named {name}"));
    (r.entries, r.bytes)
}

fn p(v: &[u32]) -> Partition {
    Partition::new(v.iter().copied())
}

/// Every answer a sweep of the cached kernels gives, keyed by what was asked.
fn sweep() -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    let parts: Vec<Partition> = (0..=6).flat_map(partitions_of).collect();
    for a in &parts {
        for b in &parts {
            if a.size() + b.size() <= 8 {
                let prod = symfn::AutoLr.schur_product(a, b);
                out.insert(format!("s{a} * s{b}"), format!("{prod:?}"));
            }
            if a.size() == b.size() {
                out.insert(format!("K {a} {b}"), symfn::kostka(a, b).to_string());
                out.insert(format!("chi {a} {b}"), symfn::character(a, b).to_string());
            }
            if a.contains(b) {
                out.insert(
                    format!("skew {a}/{b}"),
                    format!("{:?}", symfn::expand_skew(a, b)),
                );
            }
        }
        let h: symfn::Homogeneous<i64> = symfn::Homogeneous::monomial(a.clone(), 1);
        let s: symfn::Schur<i64> = symfn::convert(&h);
        out.insert(format!("h{a} in s"), format!("{s}"));
    }
    out
}

#[test]
fn a_budget_clears_whole_tables_in_tier_order_and_changes_no_answer() {
    clear_caches();
    assert_eq!(cache_budget(), None, "the crate starts unbounded");

    // Unbounded first: the reference answers, and the tables they fill —
    // `skews` in tier 2, `kostka` and `character_masks` in tier 1.
    let reference = sweep();
    let unbounded = held();
    let (skew_entries, _) = row("skews");
    let (mask_entries, _) = row("character_masks");
    let (kostka_entries, _) = row("kostka");
    assert!(
        skew_entries > 100 && mask_entries > 100 && kostka_entries > 100,
        "the sweep fills the three tables the test reads: {skew_entries}, {mask_entries}, {kostka_entries}"
    );

    // A budget just under what the sweep holds, then one more skew: the only
    // tier-2 table with bytes is cleared whole, and tier 1 is untouched.
    set_cache_budget(Some(unbounded - 1));
    assert_eq!(cache_budget(), Some(unbounded - 1));
    let _ = symfn::expand_skew(&p(&[5, 4, 3, 2]), &p(&[2, 1]));
    assert!(held() < unbounded, "held {} over the budget", held());
    assert_eq!(row("skews").0, 0, "tier 2 is what a budget clears first");
    assert_eq!(row("character_masks").0, mask_entries, "tier 1 survives");
    assert_eq!(row("kostka").0, kostka_entries, "tier 1 survives");

    // Tier 0 is never cleared, even by a budget of zero, under which every
    // insert clears everything a budget may touch.
    set_cache_budget(Some(0));
    let _ = symfn::kostka::kostka_table(5);
    let _ = symfn::character(&p(&[3, 2, 1]), &p(&[2, 2, 1, 1]));
    assert!(row("partitions").0 > 0, "tier 0 survives a zero budget");
    assert_eq!(
        row("character_masks").0,
        0,
        "a zero budget keeps nothing in tier 1"
    );

    // Under a budget that evicts constantly, every answer is the same, and
    // what is held never exceeds the budget by more than tier 0.
    set_cache_budget(Some(64 * 1024));
    let bounded = sweep();
    let tier0: usize = cache_stats()
        .iter()
        .filter(|r| r.tier == 0)
        .map(|r| r.bytes)
        .sum();
    assert!(
        held() <= 64 * 1024 + tier0,
        "held {} against a budget of 65536 plus {tier0} of tier 0",
        held()
    );
    for (question, answer) in &reference {
        assert_eq!(
            bounded.get(question),
            Some(answer),
            "{question} changed under a budget"
        );
    }

    set_cache_budget(None);
    assert_eq!(cache_budget(), None);
    clear_caches();
}
