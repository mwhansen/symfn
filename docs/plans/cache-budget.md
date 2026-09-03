# Plan: a byte budget for the global caches

Written 2026-09-03, from the code review recorded in
[release-readiness.md](../release-readiness.md), Phase 9, and from Rule 4 of
[memory.md](../record/memory.md), which already names the shape of the
solution. This file is the hand-off: what is built, in what order, and which
decisions are made so they are not re-made. All four stages land before the
0.9 release, because stage 1 adds entry points and the shape of an entry
point is decided before the first tag.

## Why this matters

Every table in `src/memo.rs` grows without eviction — twenty-one of them, plus
two thread-local tables in `src/convert.rs` that `clear_caches` does not
reach. There is no way to ask how much any of them holds, and the one control
is `clear_caches`, which drops everything. No consumer calls it: not the
convenience layer, not Sage's adapter, only the benches. A Sage session that
runs for a day computes a few thousand products and characters at rising
degree and keeps every one of them; the record measured `character_table` at
39.7 MB by degree 24 and past 1 GB by degree 32, and `skew_table` holding a
5M-term expansion. The session's owner sees the process grow and has no name
for the thing growing.

The caches are also what makes the library fast on the workloads the record
measures: a transition matrix kept across a degree sweep is a 17× speedup on
`s → J` and 22× on `m → P`, and the product table is what makes a second
`s_μ · s_ν` free. Dropping them is not the answer; bounding them is.

## Decisions made

1. **Bytes, not entries.** Bytes per entry vary a thousandfold between
   `lr_table` (three partitions and a `u128`) and `skew_table` (an `Arc` to
   an expansion with millions of terms). An entry cap on either is wrong for
   the other.
2. **Whole-table eviction, largest first within a tier.** Per-entry eviction
   needs a recency order, a recency order needs a write on every hit, and the
   hits are what the `RwLock` exists to keep cheap. Clearing one whole table
   is a single write lock and no bookkeeping on the read path. Which table: the
   one holding the most bytes in the lowest evictable tier. Repeat until
   under budget.
3. **Tiers by value per byte**, from the retention figures in
   [memory.md](../record/memory.md):

   | Tier | Tables | Evicted |
   |---|---|---|
   | 0 | `partitions_table`, `lex_parts_table`, and the two `convert.rs` tables once moved here | never |
   | 1 | the three character tables, `kostka_table`, `lr_table`, `bh_pieri_table`, `bh_ell_table` | third |
   | 2 | `product_table`, `skew_table`, `jt_row_table`, `inverse_kostka_row_table`, the four st/ht tables, `reduced_kronecker_table`, `bold_p_table`, `htilde_table` | first |
   | 3 | `transition_store` | second |

   Tier 2 holds large expansions that cost a few times their size to rebuild.
   Tier 3 holds whole-degree matrices whose retention is 13–36% of their
   build cost but which buy the largest speedups; they go after the
   expansions and before the scalar tables. Tier 1 is what makes characters
   and Kostka numbers usable at all above degree 20 — the record's ceilings
   on those functions are the table, not the arithmetic — so it goes last.
   Tier 0 is O(p(n)) and structural.
4. **Enforcement runs on the inserting thread, after its write guard drops.**
   Never while holding a table's lock, so the clearing loop can take each
   table's write lock in turn with no ordering to get wrong. It uses
   `try_write` and skips a table another thread is reading: a character
   recursion holds its mask table's read guard for its whole run, and an
   unrelated insert must not stall behind it. A skipped table is tried again
   at the next enforcement.
5. **Clearing is safe for in-flight callers.** Every large value is behind an
   `Arc`; a caller holding one keeps it and the bytes go when the last clone
   drops. The accounting counts a table's own holdings, so the counter drops
   at clear time and the allocator follows later. That gap is bounded by what
   callers hold, which is bounded by the work in flight.
6. **The estimate is calibrated, not guessed.** A `HeapSize` trait gives each
   key and value type its bytes — headers, capacity, and the heap behind a
   bignum coefficient — and `tests/memory.rs` checks the sum of the counters
   against `measure::live()` before and after filling a table. The tolerance
   is the allocator's per-allocation overhead, measured rather than assumed.
7. **The crate default stays unbounded.** No Rust caller's behavior changes at
   0.9. The wheel sets a budget at module init, because a Sage session's
   lifetime belongs to someone else and nothing in it will ever call
   `clear_caches`. It is overridable from Python by `set_cache_budget` and
   from the environment by `SYMFN_CACHE_BUDGET` (bytes; `0` means unbounded),
   the environment variable being for the Sage user who never sees a symfn
   call. The number comes from stage 3's measurements, not from this file.
8. **`clear_caches` keeps its meaning.** It drops every table including tier
   0, as it does today, because a timing run wants a cold start.
9. **Scoped cache contexts are the right design and are not this plan.** A
   context object passed to every entry point would make two workloads in one
   process independent and would need no global budget at all. It changes
   the signature of nearly every public function. It is a 2.0 question, and
   this plan does not pretend the budget settles it.

## Stage 1 — ~~accounting and introspection~~ — done 2026-09-03

- A `HeapSize` trait in `memo.rs`, implemented for every key and value type
  the tables hold: `Partition`, tuples of them, `u128`/`i128`, `Rat<i128>`,
  `Frac<C>`, `AFrac<C>`, `QtPoly<C>`, `Schur<QtPoly<i128>>`,
  `PowerSum<GuardedRat>`, `Vec<T>`, `BTreeMap<K, V>`, `Arc<T>`, and the
  bignum types under the feature. Inline types report zero heap; the trait
  method is `heap_bytes`, and the table adds `size_of` of the key and value
  itself.
- A per-table `AtomicUsize` beside each `OnceLock`, added on insert and
  zeroed on clear. The `table!` macro grows a name, a tier and the counter.
- `cache_stats() -> Vec<CacheStat>` at the crate root, one row per table:
  name, tier, entries, bytes. Re-exported beside `clear_caches` for the same
  reason it is.
- The two `convert.rs` thread-local tables move into `memo.rs` as tier 0,
  so `clear_caches` and `cache_stats` see them.
- `tests/memory.rs` gains the calibration test of decision 6, over the two
  table shapes that dominate: `skew_table` (large `Arc` values) and
  `character_mask_table` (many small entries).
- Python: `cache_stats` in the contract layer, returning a list of
  `(name, tier, entries, bytes)` tuples in table order, stubbed and
  documented like `clear_caches`. Its home in the policy's table is "a new
  computation": whole-object, plain data.

No eviction in this stage. `cargo test`, `preflight.sh` and
`preflight_python.sh` green; the speed check of stage 4 is not needed here
because the read path is untouched.

Landed as written, with two findings recorded in
[memory.md](../record/memory.md) under Rule 4: the old `clear_caches` kept
every bucket array, and a skew expansion leaves about 8.4 MB live that no
table owns. The calibration ratio was 1.000 on both shapes.

## Stage 2 — ~~the budget and eviction~~ — done 2026-09-03

- `set_cache_budget(Option<usize>)` and `cache_budget() -> Option<usize>` at
  the crate root; the budget is an `AtomicUsize` with `0` for unbounded.
- After `lookup` and every other insert path drops its write guard, it calls
  `enforce()`, which compares the sum of the counters to the budget and
  clears the largest table in the lowest evictable tier until under, with
  the `try_write` rule of decision 4.
- Tests: an insert that crosses a small budget clears exactly the expected
  table and no other; the result of every cached function is equal before
  and after an eviction, over the same sweep `tests/algebra_laws.rs` uses.
- Python: `set_cache_budget` and `cache_budget` in the contract layer; the
  `#[pymodule]` init reads `SYMFN_CACHE_BUDGET` and applies it, else the
  wheel default, which is a named constant in `python.rs` and is
  **unbounded until stage 3 sets it**.

Landed as written. `tests/cache_budget.rs` pins the tier order over the
tables a degree-6 sweep fills (`skews` goes, `kostka` and `character_masks`
stay), that a zero budget leaves tier 0 alone, and that every answer of the
sweep is unchanged under a 64 KB budget. A malformed `SYMFN_CACHE_BUDGET`
fails the import with `ValueError` rather than being ignored.

## Stage 3 — ~~the session workload, the default~~ — done 2026-09-03

- One entry in `src/measure/workloads.rs`: a sweep over rising degree of
  products, skews, characters, Kostka numbers and one transition matrix,
  shaped like a session rather than a benchmark. Under a budget, its live
  bytes must stay within the budget plus the largest single insert.
- Run it unbounded first, on AC power, and record per-table bytes at each
  degree in [memory.md](../record/memory.md) — the census Rule 4 asks for.
- From that census pick the wheel default: the smallest budget at which the
  sweep's speedup over a cold run stays within a stated fraction of the
  unbounded speedup. Record the fraction, the number, and the harness.

Landed, with one change to the recipe: there is no fraction to state,
because the cost has no slope — a budget above the working set costs
nothing measurable and one below it costs multiples at once, since eviction
is by whole table. The default is therefore chosen for headroom over the
working set rather than from a curve: 1 GiB, twenty times what the census
sweep holds, biting only at the degrees where the record's memory walls sit
anyway. The census, the ladder and the argument are in
[memory.md](../record/memory.md) Rule 4; the run was on battery, and the AC
rerun is the one open item of this stage.

## Stage 4 — ~~the speed check~~ — done 2026-09-03

- `bench_lr` and `bench_ops`, interleaved, the tree before stage 1 against
  the tree after stage 2, both unbounded. The miss path gains one atomic add
  and one compare; the hit path gains nothing, and the interleaved runs are
  the proof. If a row moves past the noise floor, the accounting is wrong
  somewhere on the hit path and the stage is not done.
- The numbers land in [memory.md](../record/memory.md) with the power state.

Landed. `bench_lr` is at parity throughout; `bench_ops` is at parity except
the character sweeps, which the record traces to the old `clear_caches`
retaining bucket arrays between cases, and a fresh-process comparison
confirms the parity. On battery; the AC rerun is owed with stage 3's.

## Rejected, with the premise

- **LRU or clock eviction per table.** A write on every read, on the path
  that carries all the reuse.
- **Weak references as the cache.** A value drops the moment its caller
  drops it, so nothing is shared across calls.
- **Entry caps.** Decision 1.
- **A budget derived from physical memory.** Needs a platform query and a
  dependency; the crate depends on nothing. A fixed default with an
  environment override is what a packaging audit can read.
- **One bounded default for crate and wheel.** Rejected for 0.9 only because
  it would change every Rust caller's behavior in a release whose point is
  to change nothing they see. Re-examine at 1.0 with the census in hand.

## Surfaces this touches

`src/memo.rs`, `src/convert.rs` (the two tables), `src/lib.rs` (re-exports),
`src/python.rs` (three entry points, module init), `python/symfn/symfn.pyi`,
`src/measure/workloads.rs`, `tests/memory.rs`,
[memory.md](../record/memory.md) (Rule 4 is amended, the census and the
numbers land in it), [public-api.md](../public-api.md) (the re-exports beside
`clear_caches`), and the Phase 9 cache item in
[release-readiness.md](../release-readiness.md), which points here.
