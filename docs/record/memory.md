# Memory: how we measure it, and what the measurements say

Memory numbers in this project have been inconsistent, and the reason is that
"memory" has meant three different quantities that move independently:

| quantity | what it is | what changes it |
|---|---|---|
| **peak live heap** | bytes the computation holds at its high-water mark | the representation, and how much is alive at once |
| **total allocated** | every byte ever handed out ("churn") | how many temporaries the code builds |
| **peak RSS** | pages the process holds from the OS | the two above, *plus* allocator retention and history |

RSS is the one that is easy to measure and the one that drifts. It depends on
allocation *history* — the order and size mix of every request the allocator has
seen — so two runs of the same binary on the same input can disagree, and a
change that provably reduces data can raise it. Every claim below is therefore
anchored to peak live heap and allocation count first, with RSS quoted only as
the end-to-end check.

## The harness

`src/measure/` is the shared accounting: a `GlobalAlloc` wrapper that counts,
`Stats` (peak / total / allocs), a size-class histogram, and `Budget`. It has no
dependencies — `criterion` and `dhat` are both out on the crate's
dependency-free rule — and costs nothing in a binary that does not install it.

`src/measure/workloads.rs` is the **single catalog**, and it is what makes this
cheap to keep up. One `Workload` entry gives you both halves:

```bash
cargo run --release --example heapstat -- htilde   # explore: report + histogram
cargo test --release --test memory                 # enforce: all budgets, 0.7s
```

Adding a subsystem is one entry: write it with `peak: 0, allocs: 0`, run
`heapstat`, and paste back the `Budget` line it prints. `HEAPSTAT_HIST=1` adds
the size-class histogram, which is what attributes churn to a specific buffer: a
spike in the 1–8 KB classes is polynomial arithmetic, a spike in the 32-byte
class is one `Partition` per output term.

The workloads span the subsystems, plus `skew-clone` and `schur-mul` for the two
call patterns of §Rule 2. `examples/lrheap.rs` predates this and stays because it
divides by the layer counter to report bytes-per-state; it now uses the same
allocator instead of its own copy.

### The counter was unsigned, and `reset()` had already falsified that

Turning on `overflow-checks` in the release profile
([failure-and-overflow.md](failure-and-overflow.md)) failed the budget test
immediately, inside the harness itself. `measure::reset()` zeroes the live-bytes
counter while memory allocated *before* the measurement is still held — the test
harness's own capture buffer is the reliable example — and freeing those blocks
then drives the counter below zero. Unsigned, that underflowed to ~2^64 and the
high-water mark latched it, so a budget could report a peak that was never
allocated.

Live bytes are now `isize`, floored at 0 on the way out, and `tests/memory.rs`
opens with a region that only frees. **The unsigned type was not a detail that
happened to be wrong** — it was the claim "this only goes up", which `reset()`
contradicts by design.

### Why memory can be a test when time cannot

This project has learned the hard way that timings need AC power, rotated
passes, and min-of-three, and are still not assertable. Peak live bytes and
allocation count are different: **a single-threaded workload reports them
bit-exactly across runs.** Three consecutive runs of the degree-10 (q,t)-Kostka
table give 16.5 MB / 786.9 MB / 729 802 allocations, identical every time, while
its RSS moves. Workloads that fan out across threads move about 0.1% on
allocation count and 1% on peak, because the shard split follows scheduling —
hence `Budget::tolerance`, tight (5%) for the single-threaded workloads and
looser (10%) for the LR ones.

Budgets are **ceilings with headroom, not golden values**, so ordinary churn
never touches them. A failure prints the measured numbers and the `Budget` line
to paste back, because the right answer is often "this legitimately costs more"
and the point is to make that a deliberate, visible edit:

```text
memory budget exceeded for `htilde`
  peak   16.5 MB  (budget 11.4 MB, +5% = 12.0 MB)
  allocs 729801  (budget 700000, +5% = 735000)
  If the increase is intended, update the budget to:
    Budget { name: "htilde", peak: 17325012, allocs: 729801, tolerance: 0.05 }
```

Two constraints the design has to respect. The counters are **process-wide**,
not thread-local, because the LR expansion spawns workers and a thread-local
counter would silently report the merge as free — the same class of error as the
`PEAK_LIVE_STATES` counter that sampled after the merge and made sharding look
memory-neutral. So measurements must not overlap: every budget runs inside one
sequential `#[test]`, since the harness runs test *functions* in parallel. And
the budgets are calibrated for `--release`; the test asserts nothing under
`debug_assertions`.

Baseline, as measured:

| workload | peak | total | allocs | churn |
|---|---|---|---|---|
| `skew-mid` | 4.1 MB | 16.0 MB | 48 339 | 3.9x |
| `skew-big` | 82.5 MB | 486.2 MB | 499 499 | 5.9x |
| `product` | 7.3 MB | 25.7 MB | 55 013 | 3.5x |
| `coproduct` | 2.2 MB | 5.8 MB | 90 750 | 2.6x |
| `htilde` (deg 10) | 16.5 MB | 786.9 MB | 729 802 | **47.6x** |
| `s-in-j` (deg 9) | 9.7 MB | 943.8 MB | 564 971 | **97.7x** |
| `m-in-p` (deg 9) | 7.0 MB | 5515.9 MB | 624 097 | **787.0x** |
| `hl` (deg 12) | 1.7 MB | 4.4 MB | 35 654 | 2.6x |
| `llt` (9, 3) | 0.2 MB | 7.7 MB | 80 565 | **43.3x** |
| `jack` (deg 9) | 0.2 MB | 1.8 MB | 40 149 | 11.7x |
| `kostka-foulkes` (deg 12) | 1.8 MB | 4.9 MB | 38 439 | 2.7x |
| `character` (deg 24) | 38.7 MB | 499.4 MB | 11 221 | 12.9x |

`s-in-j` is the `s → J` transition matrix of a degree, and it was the one
workload here whose result was **retained**: `memo::schur_in_j_cached` holds it
after the call. That is measured separately, since the table above reports the
computation and not what survives it — a second, warm call allocates one copy
and nothing else, at 0.5 MB for degree 8, 1.4 MB for degree 9 and 3.7 MB for
degree 10, against cold peaks of 3.8, 9.4 and 25.0 MB. So the cache keeps
13–15% of what building it costs, which is the trade `docs/record/qt-kostka.md`
records against a 17× speedup for a caller expanding one shape at a time.

`m-in-p` is `m_{(5,3,1)}` written in the Macdonald `P` basis, and it is
retained the same way, by `memo::mac_p_inverse_cached`. Measured with
`measure::live()` after a cold call — which is the retention quantity, since
`peak` is a high-water mark and cannot tell a table that was built and kept
from one built and dropped: 0.33 MB at degree 7, 0.99 MB at 8 and 2.64 MB at
9, against cold peaks of 1.12, 3.02 and 7.35 MB. So this cache keeps 29–36% of
what building it costs, against the 22× speedup
[macdonald.md](macdonald.md) records for a caller sweeping a degree one shape
at a time. Both figures are higher than `s-in-j`'s 13–15%, because the table
here is a solve over ℚ(q,t) rather than a matrix read off one projection.

Its churn — **787×, the highest here** — is the `Frac` arithmetic of the
back-substitution, which builds and drops a numerator polynomial at every step
of every entry. Rule 1 below says that is a CPU cost and not a memory one
while the sizes stay uniform and the buffers are freed promptly, which at 7.0
MB peak against 5.5 GB allocated they evidently are. It has not been profiled;
[macdonald.md](macdonald.md) records it as the first place to look if this
direction is worth another pass.

## Rule 1: churn costs memory only when sizes are diverse or buffers retained

`htilde` allocates 787 MB to hold 16.5 MB. That looks like the obvious target,
and it was measured and rejected.

`QtPoly::add_shifted` merged into a **fresh output `Vec` per call**, so folding
`m` terms into an accumulator allocated `m` times, each pass copying the whole
accumulator and dropping the previous buffer — quadratic bytes for a linear
result. Replacing it with a standard in-place backward merge (open the new slots
at the end, fill from the back, close the gap left by cancellations with one
`drain`) is a genuine algorithmic improvement and did exactly what it promised:

| degree-12 `bench_htilde` | churn | peak live | peak RSS | time |
|---|---|---|---|---|
| as shipped | 787 MB | 16.5 MB | 177.6 MB | 1.040s |
| in-place, room for `n` | 397 MB | 20.9 MB | **182.1 MB** | 1.079s |
| in-place, exact new-key count | **366 MB** | 16.6 MB | 166.8 MB | **1.204s** |

**Halving the churn moved RSS by 6% and cost 16% of the run time.** The reverted
buffers were all the same size and were freed immediately, which is the best case
for a size-class allocator: it hands the same block straight back, so the churn
never becomes residency. The variant that reserved room for every incoming term
rather than counting the new ones first was *worse than baseline on RSS*, because
the slack it left behind is retained by every stored polynomial.

Both variants are reverted. Churn matters when the sizes are **diverse** (the
allocator cannot recycle a block into a differently-shaped request) or when the
buffers are **retained**. Uniform, promptly-freed churn is a CPU cost, not a
memory cost, and should be judged as one.

This is the same shape of result as the layer-pooling experiment in
[littlewood-richardson.md](littlewood-richardson.md#-tried-the-obvious-fix-it-made-things-worse):
the allocator is usually doing a better job than a hand-rolled scheme, and
"reuse the allocations" is not a memory argument on its own.

## Rule 2: peak is the constraint, and duplication is where it goes

The changes that pay are the ones that stop the same data existing twice.

**Shipped: `expand_skew_shared`.** The memo holds an `Arc<Vec<(Partition, u128)>>`
and `expand_skew` returned a **deep clone of it**. Measured on `[8,7,6,5,4,3]²`
(`heapstat skew-clone`), handing back an expansion the cache already held cost
**164 041 allocations and 14.1 MB per call** — one allocation per term, because
every `Partition` owns a `Vec` — while an identical 14.1 MB sat in the cache. At
the 5.3M-term scale of `[24,20,16,12]²` that is hundreds of MB of pure
duplication, and it was the standing item 5 ("output residency") on the LR list.

`expand_skew_shared` returns the `Arc`. Every in-crate caller only iterates, so
`hopf` (three sites), `llt`, and `SkewLr::lr_coeff` take it; `StripLr` gained the
same split as `product_shared`. `expand_skew` stays for callers that want an
owned vector. `SkewLr::lr_coeff` was the worst case: reading **one** coefficient
deep-copied the entire expansion.

The general rule: **a memoized value returned by clone is a design error.**
Return the `Arc` and let callers copy only if they must.

**Shipped 2026-08-18: `LrBackend::schur_product_shared`, and a copy-free
`Schur::mul`.** The rule above had one violator left, and it was the main
entry point: `Schur::mul_with` took the trait's owned `schur_product` — the
deep clone — and then `SymFn::add_term` copied every key again to serve its
rare cancel-and-remove path, so a product's terms existed three times over
during the loop. `mul_with` now reads the shared expansion, builds the first
pair's terms into the map in one pass through an exactly-sized vector, and
accumulates later pairs by reference; `add_term` goes through `Entry` and never
copies its key. Measured on `[8,7,6,5,4,3]²` with the product warm (`heapstat
schur-mul`, the workload added for it): **peak 24.1 → 17.0 MB, 355 418 →
178 959 allocations**, and 3.5–4.2x faster on the same products
([littlewood-richardson.md](littlewood-richardson.md), "the consumer side").
The one-pass build's vector is itself a transient duplicate — 32 bytes per
term beside the ~70 the map holds, and 3 MB more than that when it was left to
grow by doubling — accepted for a 6x on this stage over sorted insertion; the
threshold that would give the largest shapes their peak back is in that
file's open tail.

## Rule 3: know which allocations are structural

`skew-big` makes 499 505 allocations, of which **471 907 are in the 32-byte size
class** — that is `Partition`'s `Vec<u32>`, one heap allocation per partition,
and partitions are the unit of output everywhere in the crate. They are only
20.4 MB of bytes, but allocation *count* is what drives fragmentation and the
RSS-over-peak gap, and this is 94% of the count.

An inline representation — parts stored in the struct up to some length, spilling
to the heap beyond it, exactly what `skew_lr::Key` already does for layer
states — would remove essentially all of them. It is invasive (`Partition` is the
crate's most-used type) and the byte-level trade is close to neutral, so it is a
**prototype-behind-the-harness** job, not a refactor to start on faith. `Key` is
the precedent that it works: inlining layer states removed a malloc per state
and took the allocator from a measured 38% of wall time to 5%.

## Rule 4: the caches are unbounded, and that is a policy, not an oversight

Every table in `memo.rs` grows without eviction; `clear_caches()` is all-or-
nothing. That suits interactive research at modest degrees, and the (q,t)-Kostka
measurements show the sharing is close to free — degree-12 peak RSS moves
156 MB → 157 MB for a 1.5x speedup. It is a genuine hazard for a long-running
Sage session, where nothing ever calls `clear_caches`, and for the dense tables
(`character_table` is p(n)² `i128`: 39.7 MB at degree 24, past 1 GB at degree 32,
which is the real ceiling on that function, well before the `i128` one).
`kostka_table` has the same shape and the same ceiling: 1.1 GB at degree 32 and
22 GB at degree 40, against a precision wall at n ≈ 58 where the table would be
8 TB. Both rustdocs state the memory wall and not the precision one, because it
is the memory wall a caller meets.

If this needs to change, the useful step is per-table byte accounting and a
budget, not an LRU on every table — the tables have very different value per
byte, and `partitions_cached` should never be evicted while `skew_table` is
holding a 5M-term expansion.

## Checklist for the next memory change

0. **If you added a subsystem, add a workload.** One entry in
   `src/measure/workloads.rs`, sized where the work actually lives. Everything
   below then happens for free, and keeps happening.
1. **Measure with `heapstat`, not `/usr/bin/time -l`.** Quote peak live and
   allocation count. RSS is the last check, not the first.
2. **Say which of the three quantities you are moving**, and check the other two
   did not move the wrong way — the reverted merge above raised time 16% and the
   reserve-`n` variant raised RSS while lowering churn.
3. **Churn needs a second argument**: diverse sizes, or retention. Uniform
   promptly-freed churn is a CPU question.
4. **Look for duplication first.** A clone of a cached value, a buffer live at
   the same time as the thing it will become, output held alongside the layer
   that produced it.
5. **Interleave and rotate benchmark passes, on AC power.** Position and power
   state have both produced double-digit phantom differences in this project.
6. **Record negative results here.** Two "obvious" memory fixes have now been
   implemented and reverted; each cost a day that the note would have saved.
