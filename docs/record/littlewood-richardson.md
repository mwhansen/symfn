# Littlewood–Richardson: backends, counting, parallelism, memory

The LR engine is the deepest single piece of work in the crate: four
backends, two per-output counting routes, an orientation dispatch, a
byte-packed layer, and a parallel traversal. It is also where the
project's measurement discipline was learned, so most of the warnings in
this file are about benchmarking rather than about mathematics.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

## Two findings worth remembering

**Sage memoizes.** A benchmark that repeats an identical computation measures
Sage's cache, not its algorithm — under that condition Sage originally beat
symfn (0.89x). Adding our own memoization turned that into ~8x. Always benchmark
with distinct, cold inputs (`scripts/bench_vs_sage.py` does).

**The DP backend was not a universal win — and the fix was to ask a different
question.** `StripLr` has better asymptotics than `NaiveLr` (DP states merge)
but higher constants (hashing and allocating per state), so it *lost* on small
inputs and the default had to dispatch on |μ|+|ν| against an empirical
threshold.

What removed the trade-off was noticing that both backends answered the wrong
question. Both compute "what is c^λ_{μν}?" and then repeat the entire
computation for every candidate λ. `SkewLr` expands the shape itself, letting
each filling report its own content — so one pass produces every nonzero
coefficient, and the cost stops scaling with p(n).

**But one traversal is necessary, not sufficient.** The shape behind
s[8,7,6,5,4,3]² has 2.1×10⁸ LR tableaux across 164 037 output terms; touching
them individually is hopeless no matter how tight the inner loop. The real win
is that the state carried between rows is only *(content so far, the row above
clipped to the columns the next row overlaps)*. Everything undecided depends on
exactly those, so partial fillings that agree on them merge into one weighted
state.

A later pass took a further ~2.2x on top of that, uniform across every
non-trivial shape, by filling each row **by runs** rather than cell by cell — a
weakly increasing row is a sequence of runs, and both the column-strictness and
ballot constraints reduce to O(1) per run — and by packing the layer key into
a single buffer, so a transition that *merges* (the common case, and the whole
point of a layer) allocates nothing. `examples/bench_shapes.rs` is the
interleaved A/B harness for measuring that kind of change.

Measured (`cargo run --release --example bench_lr`). **Re-measured 2026-07-30**
— min of 3 consecutive runs, bracketed by an external control (an unmodified
`schubmult` binary: 0.32s before, 0.34s after, against 0.355s recorded the
previous day) so that machine drift could be distinguished from code change.
The three runs agreed to within 0.8% on the dominant row.

| product | NaiveLr | StripLr | SkewLr | vs best |
|---|---|---|---|---|
| s[5,4,3,2,1]² | 0.0097s | 0.0091s | 0.0012s | **7.6x** |
| s[6,5,4,3,2]² | 0.0786s | 0.0619s | 0.0036s | **17.2x** |
| s[6,5,4,3,2,1]² | 0.4154s | 0.1406s | 0.0080s | **17.6x** |
| s[7,6,5,4,3]² | 0.6023s | 0.4497s | 0.0121s | **37.2x** |
| s[8,7,6,5,4,3]² | 151.65s | 14.23s | 0.1287s | **110.6x** |

⚠️ **What changed since the previous version of this table, and what it does
and does not license concluding.** Every cell got faster: `NaiveLr` 327.69s →
151.65s and `StripLr` 30.14s → 14.23s, both **≈2.1× uniformly**, while
`SkewLr` moved 1.2× at the smallest shape and 6.5× at the largest. The
external control rules out the machine having drifted *between yesterday and
today* — it does **not** rule out this machine being faster than whatever ran
the original table, whose date is not recorded. So:

* the **uniform ≈2.1× across all three backends** is not attributable and
  should be read as ambient (machine, toolchain, allocator);
* the **differential** — `SkewLr` improving progressively more with size,
  which no uniform factor can produce — is real code improvement, and it is
  what moved `vs best` on the largest shape from 36.2× to **110.6×**.

The lesson is procedural: a benchmark table with no date and no control is not
a baseline, and this one silently understated the library by up to 6× for an
unknown period. Tables here should carry both, as this one now does.

The win is larger still on skew expansions, where the old path ran a full
backtrack per candidate content. ⚠️ **These rows are now below this harness's
measurement floor and their speedups have been withdrawn.** `SkewLr` reads
0.0000–0.0001s on all three, and the *same* case varied 2.2×–6.1× across the
three runs — so the previously published 124×, 3040× and 257× are not
reproducible figures, they are noise divided by noise. The direction (`SkewLr`
much faster) is not in doubt; the magnitudes need a harness that loops each
case rather than timing it once.

| skew shape | per-ν | SkewLr | speedup |
|---|---|---|---|
| s[8,7,6,5,4]/[3,2,1] | 0.0001s | <0.0001s | below floor |
| s[9,9,8,8,7]/[2,1] | 0.0012s | <0.0001s | below floor |
| s[10,9,8,7,6,5]/[4,3,2,1] | 0.0009s | 0.0001s | below floor |

`SkewLr` wins at every size measured — checked explicitly at small sizes, where
the layer map's hashing could have dominated; it does not. So `AutoLr` no
longer dispatches. It stays a distinct type as the one place to reintroduce
dispatch if a future backend wins only in some regime.

**Where that leaves us against lrcalc: ahead everywhere except moderate wide
shapes.** Run `scripts/compare_lrcalc.py`, which drives both as CLI processes on
identical inputs, verifies the outputs agree, and times best-of-N. Only cases
well clear of the ~6ms process-startup floor say anything about either
algorithm:

Measured on **AC power** (see the variance warning below — an earlier run of
this same table on battery was uniformly ~2x pessimistic and distorted
`[16,13,10,7]²` in particular):

| product | lrcalc | symfn | |
|---|---|---|---|
| s[8,7,6,5,4,3]² | 4.19s | 0.378s | **us 11.1x** |
| wide [16,13,10,7]² | 7.15s | 1.26s | **us 5.7x** |
| s[9,8,7,6,5]² | 0.421s | 0.169s | us 2.5x |
| s[6,5,4,3,2,1]² | 0.033s | 0.015s | us 2.2x |
| s[8,7,6,5,4]² | 0.141s | 0.069s | us 2.0x |
| rectangle [14⁷]² | 0.113s | 0.069s | us 1.7x |
| s[7,6,5,4,3]² | 0.037s | 0.023s | us 1.6x |
| rectangle [12⁶]² | 0.018s | 0.012s | us 1.4x |
| asym [18,14,10]×[9,7,5] | 0.007s | 0.008s | lrcalc 1.2x |
| wide [12,10,8]² | 0.006s | 0.008s | lrcalc 1.3x |
| asym [14,12,10,8,6]×[7,5,3] | 0.015s | 0.019s | lrcalc 1.3x |
| asym [16,13,10,7]×[8,6,4] | 0.012s | 0.015s | lrcalc 1.3x |
| wide [14,12,10]² | 0.010s | 0.013s | **lrcalc 1.3x** |
| wide [20,16,12]² | 0.090s | 0.143s | **lrcalc 1.6x** |
| wide [24,20,16,12]² | >200s | **58.8s** | we finish, lrcalc doesn't |

**The asymmetric losses are not a new regime.** Every one of them has ℓ(ν) = 3
— `[8,6,4]`, `[7,5,3]`, `[9,7,5]` — and the chain's depth is ℓ(ν), so they are
the *same* few-row deficit as the wide three-row band, seen in a second family.
`[20,16,12]×[10,5]` has ℓ(ν) = 2 and wins 1.10x. This matters for
prioritisation: it is one defect with two symptoms, not two. It also means the
three-row counting prototype would in principle cover both — though its
crossover sits near 64k terms and these cases have ~12k, so it would not help
*these* without further work.

**We buy speed with memory, and that had been invisible.** `RSS=1` on the
comparison reports peak resident set per side (a separate invocation, so the
`/usr/bin/time -l` wrapper never contaminates timings):

| case | lrcalc | symfn | |
|---|---|---|---|
| rectangle `[14^7]²` | 15.6 MB | 21.5 MB | 1.4x |
| wide `[16,13,10,7]²` | 46.6 MB | **229.3 MB** | **4.9x** |
| skew `[13,12..2]/[5,4,3,2,1]` | 7.6 MB | 21.7 MB | 2.9x |

So the 5–16x time wins come with 1.4–4.9x the memory. The layer trades
residency for speed by construction, and on `[24,20,16,12]²` that reaches
2.36 GB — which is why lrcalc's inability to finish that case is not purely a
speed result.

⚠️ **Except where a section marks otherwise, the timings in this file were
taken on battery**, so absolute times are not comparable across runs and the
per-case ratios are the durable quantity —
see "Power state" in [README.md](README.md), whose sweep-throttling finding
came from `[16,13,10,7]²` here. The memory numbers above stand; they are not
timing-sensitive.

**Skew expansions and single coefficients had never been measured** — every
case for both sat inside the ~4ms process-startup floor, so those rows timed
`exec` rather than either implementation. Sizing them is not like sizing a
product: a skew expansion's term count tracks the diagram's **row count**, not
its width (`[24,20,16,12]/[8,6,4,2]` is 52 cells and yields 168 terms, while
`s[8,7,6,5,4,3]²` — 66 cells over 12 rows — yields 164 037). A first attempt
used wide four- and five-row shapes and stayed at the floor.

Once sized by rows, skew is the regime we are *strongest* in, and it had been
invisible for the whole project:

| skew case | terms | lrcalc | symfn | |
|---|---|---|---|---|
| `[12,11..3]/[4,3,2,1]` | 4 527 | 0.011s | 0.009s | 1.26x |
| `[14,13..7]/[6,5,4,3,2,1]` | 3 828 | 0.023s | 0.010s | 2.30x |
| `[13,12..2]/[5,4,3,2,1]` | 42 325 | 1.07s | 0.069s | **15.41x** |

**And the coefficient rows found a real defect.** `lr_coeff` always expanded
λ/μ, which has |ν| cells. On `c^λ_{μν}` with λ = `[13,12..2]`, μ = `[5,4,3,2,1]`,
|ν| = 75, that built all 42 325 terms of a 75-cell expansion to read a single
coefficient: **91ms against lrcalc's 4.4ms, a 21x loss**. Since
c^λ_{μν} = c^λ_{νμ}, peeling off the *larger* factor instead leaves |μ| = 15
cells. Fixed; the same query is now 3.9ms, i.e. 1.34x ahead of lrcalc — a 23x
improvement that no product benchmark could ever have surfaced.

Asymmetric products were also newly measured (the sweep was almost entirely
`s_μ²`, which cannot see a cost depending on which factor supplies the strips):
0.80–1.00x, so mild losses, in a regime previously invisible.

The sweep drives `AutoLr`, the backend a caller actually gets, not `SkewLr`
directly — `examples/lr_cli.rs` named `SkewLr` until the rectangle path landed,
which would have made that path invisible to every row here.

The ratios moved against the previous run of this table (13.6x → 10.5x,
8.2x → 5.2x) entirely because of lrcalc's own drift: our times are unchanged to
three digits (2.8056s vs 2.81s on `[16,13,10,7]²`) while lrcalc went 23.0s →
14.5s on the same binary and input. Treat both as the same result.

The `[24,20,16,12]²` row needs `TIMEOUT≥200`; at the default 120s (or the 90s
used in one sweep) *both* sides time out and the row says nothing.

⚠️ The rectangle rows in `CASES` are useless as written for the rectangle path:
`[5^5]²` (252 terms) and `[7^7]²` (3432) both sit inside the ~5ms process
startup floor, so they measure `exec`. `[12^6]²` and `[14^7]²` were added for
this reason but have not yet been run out-of-process; in-process the closed
form is 38–47x (see `src/rect.rs`).

Noise is ±30% run to run; treat anything inside ±20% as a tie. lrcalc's own
timing on an unchanged binary drifted 6.3s → 8.6s → 11.2s across this project's
sweeps, so single-digit-percent differences mean nothing. The startup floor also
moves with machine load — it was ~14ms in the sweep above, so every row under
~0.02s there is measuring `exec`.

**The remaining weakness is few-row factors, and it is structural rather than a
tuning problem.** We lose 1.2–1.6x whenever the factors have two or three rows.
It does *not* extend to wide shapes generally — `[16,13,10,7]²` is wide and we
win it — so width is not the variable.

Profiled with macOS `sample` (`examples/profile_wide.rs`, `--profile profiling`),
`[20,16,12]²` spends 54% of its time in `fill_runs` and 37% committing states,
while `[16,13,10,7]²` — which we win 5.2x — spends only 16% in `fill_runs` and
32% in the outer merge. Instrumenting the layer explains the difference, and
it is not the one previously recorded here:

| factor rows | shapes | LR tableaux ÷ states produced |
|---|---|---|
| 2 | `[24,12]²`, `[20,10]²` | **1.0x** |
| 3 | `[12,10,8]²`, `[20,16,12]²` | **1.0–1.1x** |
| 4 | `[16,13,10,7]²` | 40.8x |
| 5 | `[9,8,7,6,5]²`, `[12,10,8,6,4]²` | 11–164x |
| 6–7 | `[8,7,6,5,4,3]²`, `[7,6,5,4,3,2,1]²` | 57–80x |

The layer's entire advantage is that one state stands for many tableaux. At
two or three rows it stands for **one** — `[20,16,12]²` produces 2 614 952
states against 2 802 764 LR tableaux. So the DP walks exactly what a
tableau-at-a-time enumerator walks and then pays key assembly, encoding and
hashing on top of it. That is the whole deficit, and no amount of tuning the
layer removes it; the layer is the cost.

⚠️ This corrects an earlier claim here that few rows mean *little merging*.
Merging within the layer is in fact strongest exactly where we lose:
`[20,16,12]²` collapses 2.6M productions into 104 876 live states (24.9x), the
highest measured, while `[16,13,10,7]²` manages only 2.2x. Those are different
quantities — productions-per-live-state says the layer stays small,
tableaux-per-production says whether the enumeration was avoided — and only the
second is a saving. The earlier note conflated them.

**A layer-free enumerator was prototyped and does not help — measured, and
worth not repeating.** The obvious consequence of 1.0x compression is that
few-row shapes should skip the layer: walk the same run decompositions
depth-first, keep the previous row on the stack instead of in a hashed key, and
hash only completed tableaux keyed on content. That was built and verified
against `AutoLr` on every shape below. It is **parity at best**:

| shape | layer | direct | |
|---|---|---|---|
| `[12,10,8]²` | 8.35ms | 7.99ms | 1.05x |
| `[14,12,10]²` | 13.7ms | 14.9ms | 0.92x |
| `[20,16,12]²` | 266ms | 283ms | 0.94x |
| `[9,8,7,6,5]²` | 293ms | 1.95s | 0.15x |
| `[8,7,6,5,4,3]²` | 657ms | 22.9s | 0.03x |

(A first version measured 0.75x on `[20,16,12]²` purely because it used std's
`HashMap`; SipHash against the layer's own fast hasher measures the hasher,
not the algorithm. The table is after fixing that.)

The row-fill counts confirm the compression reading rather than contradicting
it — direct does 2 846 571 fills against the layer's 2 614 952 on
`[20,16,12]²` (+9%), and 243 534 430 against 3 678 951 on `[8,7,6,5,4,3]²`
(66x worse, exactly the compression the layer is buying there).

The reasoning error was treating the profile's 37% "state commit" as removable.
Removing the layer **relocates** that cost rather than deleting it: every
completed tableau still has to be hashed to bin it by λ, and there are 2.8M of
them either way. The layer hashes 2.6M longer keys spread across rows;
direct hashes 2.8M shorter keys at the last row.

So the real statement is: **both approaches touch all 2.8M LR tableaux while the
answer has only 64 335 terms.** Beating lrcalc here needs an algorithm that does
not enumerate tableaux at all — not a cheaper enumeration. Two candidate
directions, neither validated:

**Counting per output instead of enumerating chains — validated on a model
problem.** `examples/model_count_vs_chains.rs` computes `s_μ·h_a·h_b` both ways:
by building every chain μ ⊂ λ¹ ⊂ λ² (what the layer does, minus the lattice
condition), and by iterating over candidate λ² and counting the λ¹ directly.
Interlacing pins λ¹ᵢ to `[max(μᵢ, λ²ᵢ₊₁), min(μᵢ₋₁, λ²ᵢ)]` *independently*, with
Σλ¹ fixed, so the coefficient is a lattice-point count in a box on a hyperplane
— a bounded-composition count, closed form by inclusion–exclusion. Both agree
with the library:

| problem | terms | chains | candidates | work | time |
|---|---|---|---|---|---|
| `s[6,4,2]·h₆·h₄` | 173 | 869 | 174 | 5.0x | 1.6x |
| `s[20,16,12]·h₂₀·h₁₆` | 11 714 | 523 966 | 11 725 | 44.7x | 6.7x |
| `s[30,24,18]·h₃₀·h₂₄` | 52 725 | 6 203 378 | 52 771 | 117.6x | 17.5x |

Candidate enumeration is near waste-free (11 725 tested for 11 714 terms), and
the margin grows with size. So the *strategy* is sound where the fibre is a box.

**The lattice condition does not break it — measured, two-row ν.**
`examples/lr2_count_vs_layer.rs` runs the smallest LR case that carries the
real difficulty. The lattice condition constrains *prefix sums* of λ¹ rather
than individual λ¹ᵢ, so the fibre stops being a box and the closed form above
does not apply. But the admissible range for λ¹ⱼ given the running prefix Lⱼ₋₁
stays **contiguous** — `λ¹ⱼ ≥ Λⱼ + Mⱼ₋₁ − 2Lⱼ₋₁` — so a DP over rows keyed on
that prefix sum works, with each step a range-add on a difference array rather
than an enumeration. Verified against the library on every case:

| product | terms | AutoLr | counting | |
|---|---|---|---|---|
| `s[6,4,2]·s[6,4]` | 139 | 139µs | 60µs | 2.30x |
| `s[20,16,12]·s[20,16]` | 7 909 | 6.9ms | 8.2ms | 0.84x |
| `s[30,24,18]·s[30,24]` | 35 557 | 44.6ms | 51.8ms | 0.86x |
| `s[40,32,24]·s[40,32]` | 105 817 | 240ms | 214ms | 1.12x |
| `s[60,48,36]·s[60,48]` | 504 157 | **5.55s** | **1.51s** | **3.68x** |

This is an asymptotic crossover, not noise. Counting costs O(terms × rows ×
span) and its per-term time grows roughly linearly (0.43 → 1.04 → 1.46 → 2.02 →
2.99 µs); the layer costs O(tableaux) and its per-term time grows far faster
(1.0 → 0.87 → 1.25 → 2.27 → **11.0** µs). Crossover is near
`[40,32,24]·[40,32]` and the gap widens after it.

⚠️ **Measure against `AutoLr`, not against chain enumeration.** The same file
also implements the naive chain enumerator, which counting beats by 5–244x —
a meaningless number, since the layer DP exists precisely to beat chain
enumeration. An earlier version of this note quoted a "44.7x less work" figure
that compared incommensurable units: it counted *candidates tested* while hiding
a ~340-operation DP inside each candidate.

**Three-row ν: the 2-D state is tractable** — now landed as `src/three_row.rs`
and dispatched from `AutoLr`. Three strips need the prefix sums of both λ¹
and λ². The saving move is to key the state on *cells added so far by each
strip* rather than on absolute prefix sums: those are bounded by ν₁ and ν₂,
not by |μ|+ν₁. Chaining the interlacings also puts λ¹ in known bounds,
`λ¹ᵢ ∈ [max(μᵢ, λᵢ₊₂), min(μᵢ₋₁, λᵢ)]`, leaving state (λ¹ⱼ, aⱼ, bⱼ).

I predicted this would cost thousands of operations per term and be hopeless.
Measured, it is 26–2192, because the reachable state space is far smaller than
its bounding box. Self-contained (candidates enumerated, no oracle), verified
against the library:

| product | terms | layer | counting | |
|---|---|---|---|---|
| `[12,10,8]²` | 6 579 | 7.6ms | 7.6ms | 1.00x |
| `[14,12,10]²` | 12 068 | 12.5ms | 14.8ms | 0.85x |
| `[20,16,12]²` | 64 335 | 241ms | 194ms | 1.24x |
| `[24,20,16]²` | 145 505 | 1.15s | 518ms | **2.22x** |
| `[30,24,18]²` | 419 032 | 14.0s | 3.56s | **3.94x** |

A first pass using a `HashMap` keyed on the state tuple measured 3.3x *slower*
while reporting the same operation counts — the algorithm was fine and the data
structure was wrong. Generation-stamped dense tables fixed it. Worth remembering
before concluding an approach has failed.

Landed with `n = |μ|+|ν| ≥ 90`. Verified by **interleaved** A/B of the two
`lr_cli` binaries, alternating builds, min of 5 each, both repetitions agreeing
to three digits:

| case | before | after | |
|---|---|---|---|
| `[20,16,12]²` | 0.1489 / 0.1486 | 0.1184 / 0.1190 | **1.26x** |
| `[22,18,14]²` | 0.2659 / 0.2673 | 0.1906 / 0.1886 | **1.40x** |
| `[12,10,8]²` | 0.0168 / 0.0168 | 0.0170 / 0.0167 | unchanged (excluded) |

Against lrcalc that moves `[20,16,12]²` — the worst case in the sweep — from
0.59x to 0.74x. Still a loss, but the deficit is roughly halved.

⚠️ **Two calibration traps, both of which produced wrong thresholds first.**
An in-process A/B put the crossover at n ≈ 60, but it ran `SkewLr` first and
counting second every time, so counting inherited a warm allocator; measured
out-of-process, n = 60 was a *regression* (`[12,10,8]²` 0.79x → 0.53x). And an
apparent regression on that same case turned out to be measurement context —
the interleaved A/B above shows it unchanged. Neither is visible without
alternating builds in one process-per-run harness.

A second, independent idea, not yet tested: because ν has 3 rows, entries come
from `{1,2,3}` and every *column* is one of 7 subsets, so a DP keyed on
(content, previous column) would have a tiny state space. The obstacle is again
mathematical — the ballot condition is defined on the row reading word, and
whether it survives a column-wise reformulation is unresolved.

**Resolved 2026-07-31, in both directions at once: the mathematics says yes,
the measurement says it does not pay.** The ballot condition *does* survive —
see "The column question is settled" below for the argument, the exhaustive
check, and the prototype — but the resulting DP merges almost nothing, for the
same reason the layer compresses 1.0x here, so the tiny-state-space hope
was wrong. The section below records why that is structural.

The conjugate dispatch does not help either (it deliberately does not fire
here), because conjugating trades few rows for few columns and the state still
pins the tableau.

## 2026-07-31: the wide-band deficit closed — constants, a stale dispatch, and two negative results

⚠️ **Every number in this section was measured on battery (64% → 50%,
discharging)**, so its absolute times are not comparable to the AC tables
above. All conclusions rest on interleaved, same-condition, out-of-process
ratios (min of 5, one cold process per run). The lrcalc binary is conda's
`lrcalc 2.1` from the `sage-dev`
environment; the measured process-startup floor today was ~2.5ms per side,
not the ~6ms of earlier sweeps — the floor moves with machine load.

**The n ≥ 90 dispatch had gone stale: it preferred the slower backend on every
shape it admitted.** A new order-alternating in-process harness
(`examples/calibrate_three_row.rs`) measured the dispatched band at 0.76–0.89x
— counting *behind* the layer — and an out-of-process A/B of two builds
(dispatch on vs forced off) confirmed the direction: turning the dispatch off
was 1.02–1.07x faster on `[20,16,12]²`, `[22,18,14]²`, `[24,20,16]²`. The
calibration was sound when taken; what moved is recorded above — the ambient
2026-07-30 re-measurement found `SkewLr` improving most at exactly these sizes.
A dispatch bound is a measurement with a shelf life, and nothing in the tree
re-checks it; this is the second table in this file to be silently invalidated
by drift, after the undated-baseline lesson above.

**Negative result: checkpointing the fibre DP along the candidate DFS, with
difference-array transitions, loses 5–10x — measured, instrumented, reverted.**
The premises looked airtight: row j of the fibre DP reads nothing of λ beyond
λ_{j+2}, so DP layers can be checkpointed per DFS depth and shared across
every candidate extending the prefix, and the admissible (λ¹ⱼ, λ²ⱼ) transitions
form contiguous b-intervals, so a difference array can replace per-cell adds.
Both were implemented and verified (the module's exhaustive oracle sweep stayed
green); on `[20,16,12]²` the combination measured 1.42s against 159ms for the
per-candidate DP it replaced, degrading with size. Instrumentation counted, for
64 335 terms: 1 680 388 step calls (the DFS visits ~26 prefix-attempts per
completing candidate, and the checkpointed version pays a DP step on every one,
where the leaf-only design pays nothing for a prefix that dies on the size
constraint), 485M difference-array cells scanned at drain against 89M states
recovered (integration scans the whole b-run per touched group), and 117M
range-adds whose windows average *under one cell* — there was never a wide
interval to collapse; the recorded 26–2192 ops/term should have said so in
advance. What this licenses: the dense per-cell generation-stamped table is
well matched to this workload, and sharing schemes must first prune the
non-completing DFS bush before they can pay.

**The column question is settled: the ballot condition survives column order —
proved, verified, and measured not to matter.** The argument is three plactic
facts chained: the column word (columns left to right, bottom to top) is Knuth
equivalent to the row word; a word is ballot iff its rectification is the
superstandard tableau of its content; rectification is a plactic invariant. So
scanning columns **right to left, top to bottom** tests LR-ness exactly. The
μ-bounded form survives too: prepending the superstandard word of μ turns
"μ-floored ballot" into plain ballot, and concatenation respects Knuth moves,
so the μ-bounded fillings of ν that lrcalc enumerates admit the same column
scan with counters initialized at μ. Verified exhaustively (2 808 semistandard
fillings across 10 skew shapes × 2–4 letters: row-ballot ⟺ column-ballot with
no exception) and end-to-end (a prototype column transfer-matrix DP over ν's
diagram — state (partial content κ, previous column pattern τ), coefficients
read off final states as λ = μ + κ with no candidate enumeration — reproduced
`lr_cli mult` exactly on seven products including asymmetric and lopsided
factors).

The engineering answer is no. On the shapes that matter the DAG barely merges:

| product | terms | states | edges | paths (= tableaux) | paths/edge |
|---|---|---|---|---|---|
| `[8,6,4]²` | 1 185 | 5 249 | 8 305 | 7 488 | 0.9 |
| `[12,10,8]²` | 6 579 | 37 880 | 58 361 | 58 102 | 1.0 |
| `[14,12,10]²` | 12 068 | 75 657 | 112 743 | 114 081 | 1.0 |
| `[14,12,10,8,6]·[7,5,3]` | 12 279 | 70 130 | 165 655 | 291 161 | 1.8 |

One edge per tableau is enumeration wearing a hash map. This is the same 1.0x
compression the row layer measures at ℓ(ν) ≤ 3, now seen from the transposed
sweep direction, and the shared cause is now plain: **any DP that produces all
outputs in one traversal must carry partial content in its state — the output
is binned by content — and at three rows the partial content pins the filling
almost uniquely, so state-merging cannot beat enumeration no matter which way
the diagram is scanned.** Only per-output counting escapes, because fixing λ
turns content from state into constraint. That closes both "one big traversal"
directions (rows: the layer-free prototype above; columns: this one) and
leaves the fibre count as the only lane that scales past enumeration here.

**What won instead: three constants in the fibre count and one in the CLI.**
The per-candidate DP was kept exactly as designed and made ~2x cheaper:

1. **The state decode was two integer divisions.** The dense table's `touched`
   list held flat cell indices, and unpacking one cost two divisions by
   *runtime* strides on every state visit. States now pack `(λ¹, a, b)` into
   one u32 (12/10/10 bits); decode is three shift-masks. `three_row_product`
   declines shapes wider than the packing (ν₁ ≥ 1024, or first-row candidates
   ≥ 4096) and the caller falls back to `SkewLr`, which owns that regime
   regardless.
2. **The inner loop re-derived its bounds per cell.** All five clamps on the
   admissible b-interval (strip, `b ≤ ν₂`, ballot `b ≤ aⱼ₋₁`, `c ≥ 0`,
   `c ≤ bⱼ₋₁`) are monotone in λ²ⱼ, so the window is computed once per
   (state, λ¹ⱼ) and walked without checks; the `a ≤ ν₁` cut folds into the λ¹ⱼ
   loop bound. (The *difference-array* version of this same observation is the
   negative result above — the window is real, it is just too narrow to encode.)
3. **Row constants were re-read per state.** μⱼ, λⱼ, λⱼ₊₁ and the running
   Λⱼ−Mⱼ are hoisted out of the state loop; the intermediate `live` vector is
   gone in favor of iterating the touched list directly.
4. **`lr_cli` spent ~360ns a line printing.** `format!` per term plus a
   `to_string` per part is ~8 allocations a line, the same order as the whole
   DP on mid-sized products, and every comparison row pays the print path.
   One reused buffer, zero per-term allocations. (lrcalc prints through
   stdio's buffer; this only removes a handicap, it does not add an edge.)

In-process, order-alternating, min of 4 (`examples/calibrate_three_row.rs`):
counting beats `SkewLr` on every three-row-ν case measured, 1.11–1.72x —
the dispatched band that read 0.76–0.89x before the change. Out-of-process,
counting-forced-on vs forced-off, min of 5: 0.97–1.36x (the n = 36 square is
the one tie). End-to-end against the morning's HEAD binary, same five-rep
interleaved discipline, identical outputs:

| case | before | after | |
|---|---|---|---|
| `[12,10,8]²` | 8.0ms | 6.1ms | **1.31x** |
| `[14,12,10]²` | 13.6ms | 9.5ms | **1.43x** |
| `[20,16,12]²` | 113.8ms | 80.0ms | **1.42x** |
| `[22,18,14]²` | 186.1ms | 128.4ms | **1.45x** |
| `[16,13,10,7]·[8,6,4]` | 14.3ms | 10.2ms | **1.40x** |

`prefer_counting` was recalibrated from the out-of-process numbers: the
crossover drops n ≥ 90 → **n ≥ 48** (n = 36 stays out as a tie), and the
balance clause widens 3|ν| ≥ |μ| → 4|ν| ≥ |μ| to admit the measured five-row
win `[14,12,10,8,6]·[7,5,3]` (|μ|/|ν| = 3.3 at 1.36x) while still excluding
the ratio-12 tie `[30,24,18]·[3,2,1]`. Six-row μ measured 1.24x in-process but
has no out-of-process number, so `rows ≤ 5` stands until it does.

**Where that leaves us against lrcalc: ahead everywhere clear of the startup
floor.** Full `scripts/compare_lrcalc.py` sweep, 38 cases, every output
verified equal, `[24,20,16,12]²` excluded on battery (lrcalc exceeds the
timeout; the AC result above stands):

| former loss | was | now |
|---|---|---|
| wide `[12,10,8]²` | 0.80x | **1.06x** |
| wide `[14,12,10]²` | 0.76x | **1.09x** |
| wide `[20,16,12]²` | 0.71x | **1.11x** |
| asym `[18,14,10]·[9,7,5]` | 0.88x | **1.10x** |
| asym `[16,13,10,7]·[8,6,4]` | 0.80x | **1.15x** |
| asym `[14,12,10,8,6]·[7,5,3]` | 0.77x | **1.33x** |

("was" is this session's pre-change baseline under the same battery
conditions, matching the AC table above to within its noise.) Every remaining
sub-1.0 row in the sweep — three coef rows at 0.77–0.96x, two small skews at
0.88–0.93x, tall `[2⁸]²` at 0.99x — sits at 2.6–3.7ms total against a ~2.5ms
exec floor, the regime the sizing notes above already classify as measuring
`exec` rather than either algorithm. No above-floor case loses.

The counting path is single-threaded, so every ratio in the two tables above
compares algorithms; the paragraph below about parallel layer rows concerns
only the staircase and four-row-factor rows of the full sweep, whose large
layers can cross the row-parallel threshold.

**AC re-validation, same day: the bound and the conclusion hold, and the
battery caveat above is discharged.** This is the session that produced the
block-sequential-vs-interleaved artifact recorded under "Power state" in
[README.md](README.md); the rule it established is that a sweep row
disagreeing with a dedicated interleaved A/B loses. Confirmed on AC, all
out-of-process, min of 5–7:

* Crossover: counting over `SkewLr` 1.03–1.42x across the whole
  three-row family; the n = 36 square reads 0.98x in-process and 1.09x
  out-of-process — still the tie zone, still excluded — and the lopsided
  control is 0.96x, correctly excluded. The n ≥ 48 bound stands unchanged.
* Against lrcalc with the shipped binary, interleaved min of 7:
  `[12,10,8]²` 1.06x, `[14,12,10]²` 1.08x, `[16,13,10,7]·[8,6,4]` 1.16x,
  `[14,12,10,8,6]·[7,5,3]` 1.43x, `[18,14,10]·[9,7,5]` 1.02x; the larger
  squares measured 1.15–1.49x in the crossover probe. The full 38-case sweep
  agrees on every output.
* `[24,20,16,12]²` was not re-run: four-row factors do not dispatch to
  counting, so nothing in this change touches its path, and the standing AC
  result above stands.

In-process, the same AC session puts counting at 1.03–2.06x over `SkewLr`
on every dispatched row (`examples/calibrate_three_row.rs`), and the six-row-μ
case reads 1.32x in-process for the third time — the `rows ≤ 5` widening still
waits on an out-of-process number.

**Both implementations are single-threaded, and that is what makes this table
mean something.** lrcalc runs at ~99% of one core; symfn uses no threads at all.
So these ratios compare *algorithms*, not core counts. If symfn is ever
parallelised, this comparison must keep reporting a single-threaded number —
a multi-threaded wall-clock figure set against a single-threaded lrcalc would
conflate an algorithmic win with a hardware one, and note that a
parallel build pinned to one thread is not the same as a sequential build
(per-thread structures and merge machinery cost something even at N=1). lrcalc
being single-threaded is a property of its implementation, not of the problem;
its enumeration is at least as parallelisable as our layer, so threads are a
real engineering win for users but not a durable claim of algorithmic
superiority.

**The memory wall, and how it fell.** `[24,20,16,12]²` was once killed at
27m27s wall with 2.0+ GB resident and still climbing, **38% of it system
time** — allocation and page faults, not combinatorics. Peak memory is
(states) × (bytes per state), and the second factor was soft. Three changes
(same commit series, measured by interleaved A/B with
`examples/bench_shapes.rs`, which now also reports peak live layer states):

1. **Byte-packed inline keys.** Every element of a layer key is bounded by
   the shape's cell count, so keys serialize at one byte per element for
   anything practically computable and live inline in a 32-byte enum — no
   heap allocation per state at all. (Narrowing is where overflow bugs live:
   the width bound is proven in `elem_width` and tested across the 255/256
   and 65535/65536 boundaries against Pieri.)
2. **u64 map values with a checked u128 fallback.** Multiplicities are
   tableau counts with no provable narrow bound, so every merge is a
   `checked_add` and the expansion transparently reruns wider on saturation
   (exercised in tests with a u8 accumulator).
3. **Single-table layers.** Between rows the layer is drained into an
   exactly-sized `Vec`; the hash table — whose power-of-two bucket array can
   run 2–4× the payload — only exists on the side being merged into.

A trap worth remembering: shorter keys meant fewer hash-mixer rounds, which
exposed weak low-bit dispersion in `MixHasher` (multiply and a small left
rotation only move entropy upward) — a 2.7× probe-clustering slowdown on
`[20,16,12]²` until a murmur-style avalanche finalizer fixed it.

Results (peak RSS via `/usr/bin/time -l`, min-of-3 interleaved times):

| case | RSS before | RSS after | time before | time after |
|---|---|---|---|---|
| [16,13,10,7]² | 222.7 MB | 134.0 MB | 6.97s | 6.98s |
| [8,7,6,5,4,3]² | 114.0 MB | 78.8 MB | 1.04s | 0.93s |
| [20,16,12]² | 20.2 MB | 19.0 MB | 0.259s | 0.259s |
| **[24,20,16,12]²** | **killed at 27m, 2.0+ GB, climbing** | **completes: 1072s, 2.06 GB peak** (148s with the orientation dispatch below) | | |

(The table's [16,13,10,7]² time is the packed layer alone, same
orientation; the dispatch below then takes it to 2.4s.)

`[24,20,16,12]²` = 5 313 471 terms, peak 23.0M live states. Independently
re-measured through the ordinary `lr_cli mult` path: **125.6s wall, 2.36 GB
peak RSS, and 5% system time** — down from 38% when the case was
allocation-bound, which is the diagnosis confirming itself rather than just
the symptom improving. Neither lrcalc (>200s timeout, still running at 27m in
earlier sweeps) nor the old representation finishes it on this machine. Much
of the remaining 2 GB is the 5.3M-term *output* (two copies: the memoized
`Arc` plus the caller's clone), not the layer.

**This result now has an independent oracle** (`examples/verify_specialization.rs`).
lrcalc cannot finish the case, and for a long time its 5.3M terms were checked
only against our own conjugate orientation — a real consistency check, but not
an independent one, since a bug in the shared layer code reproduces itself in
both orientations.

Principal specialization closes that. Evaluating `s_μ·s_ν = Σ c^λ s_λ` at
`x = (1,…,1)` makes a scalar identity whose weights come from the hook-content
formula, which shares no code with the LR machinery, and whose left-hand side
never touches a coefficient. **5 313 471 terms, five independent n, all OK**
(peak RSS 2.69 GB).

It is a weighted checksum, not a proof, so it ships with a negative control:
perturbing each coefficient of a correct expansion by ±1 in turn is detected
412/412 times. A checksum that silently always passed would be worse than no
check, so that number is the one that makes the PASS meaningful.

Both halves are in `tests/lr_specialization.rs`, not only in the example: the
checksum over every product through `|μ| + |ν| ≤ 8` at five n, and the control
over three shapes, asserting `caught == total` so a blind spot is a red test
rather than a printed number. `examples/verify_specialization.rs` keeps the
shapes that take minutes. V7 wants the control to be part of the check rather
than a one-time experiment, and an example nothing runs is the experiment.

**Transposition, re-measured on the right axes — and now dispatched.**
Since c^λ_{μν} = c^{λ'}_{μ'ν'} the walk can run on the transposed diagram.
An earlier note said wide shapes prefer the original orientation — true at
moderate size, but it inverts exactly where it matters. Peak *states* are
nearly orientation-independent (±25% both ways on every case measured; 23.0M
direct vs 21.3M conjugate on the big one — the layer is the same
information either way). Time is not: the per-row run fill enumerates
fillings combinatorially in row width, and the conjugate bounds row width by
the original row count. Measured conjugate speedups: `[16,13,10,7]²` 2.1×,
`[18,15,12,9]²` 4.0×, `[10,9,8,7,6,5]²` 2.1×, `[11,10,9,8,7,6]²` 2.7×,
`[24,20,16,12]²` **7.9×** (135s vs 1072s, same 5 313 471 terms — a
cross-orientation agreement check as well). Three-row wide shapes still
prefer the direct orientation ~2× ([12,10,8], [20,16,12]); staircases tie at
moderate size (a staircase's conjugate is itself) and swing conjugate when
large. Memory does not decide the orientation; the fill cost per row does.

`expand_skew` now dispatches (`prefer_conjugate`: ≥ 8 rows, wider than
tall, ≥ 60 cells — empirical thresholds; every case the rule fires on
measured ≥ 2× or a tie, its known losses are rectangles like `[12⁶]²` at
~1.4× on a 45 ms case). With dispatch, `[24,20,16,12]²` runs **148s /
2.15 GB peak RSS** end to end through the default path, and terms are
conjugated back one at a time so no vector of unconjugated partitions is
materialized.

**`lr_coeff` answers sweeps from a cached product.** A caller sweeping many
λ against one (μ,ν) — the natural way to read coefficients off a product —
used to pay one fresh λ/μ traversal per λ. `lr_coeff` now peeks the skew
cache under the juxtaposed (μ,ν) shape first and answers by binary search:
sweeping all p(36) = 17 977 λ against μ = ν = [6,5,4,3] with the product
warm went 0.481s → 0.004s. One-shot queries still take the λ/μ route
(the smaller expansion) and nothing is computed speculatively.

## Parallel LR (`src/skew_lr.rs`)

The layer traversal is now multi-threaded. A row is a barrier — row r+1 needs
row r complete — so this is bulk-synchronous, and the only question is how to
split the states *within* a row. Two things mattered more than the threading
itself, and neither was obvious up front:

**Heterogeneous cores.** This machine is 4 performance + 6 efficiency cores, the
latter roughly a third the throughput. An even split leaves the row barrier
waiting on whichever chunk landed on the slowest core. Work is therefore claimed
from a shared counter in 2,048-state chunks, so a fast core takes three while a
slow one takes one.

**The merge was the Amdahl ceiling.** Combining each worker's table into one was
measured at **35–50% of wall time** on the large shapes — a hard 2x limit
regardless of core count, and exactly why the first version topped out at 1.73x.
The layer is now *sharded*: each key is routed to a shard by a cheap hash of
its tail bytes, so every copy of a key lands in the same shard whoever produced
it, and shard j can be combined independently of shard k. The merge became
parallel and the ceiling went with it.

| product | serial | parallel | |
|---|---|---|---|
| `[8,7,6,5]²` | 0.0043s | 0.0044s | 0.98x *(below threshold, untouched)* |
| `[10,8,6,4]²` | 0.0205s | 0.0182s | 1.13x |
| `[12,10,8,6]²` | 0.1063s | 0.0559s | 1.90x |
| `[8,7,6,5,4,3]²` | 0.2722s | 0.1319s | 2.06x |
| `[16,13,10,7]²` | 1.0399s | 0.3633s | **2.86x** |

Rows below 24,576 states stay on one thread, so small shapes are bit-for-bit the
old code path. Verified by checksumming every coefficient of every shape against
the serial expansion, not just the term counts — and the lrcalc oracle still
passes.

⚠️ **These numbers are AC-only, and that is not a formality.** This A/B is the
2.02x-vs-1.73x case behind "Power state" in [README.md](README.md): parallel
results measured on battery are not comparable to anything.

## Memory: two thirds of RSS is allocator retention, not data

`examples/lrheap.rs` wraps the global allocator to count live bytes, which
separates what the traversal actually holds from what the process has not given
back. The two differ by a lot, and the difference grows with the shape:

| shape | | live heap | peak RSS | retention |
|---|---|---|---|---|
| `[8,7,6,5,4,3]²` | serial | 84.5 MB | 157.8 MB | 1.9x |
| | parallel | 86.4 MB | 168.7 MB | 2.0x |
| `[16,13,10,7]²` | serial | 122.6 MB | 364.4 MB | **3.0x** |
| | parallel | 123.2 MB | 364.1 MB | **3.0x** |

Live data is ~66 bytes per layer state, which is about right for a 40-byte
`(Key, u64)` entry plus hash-table slack — the representation is not the
problem. **RSS is 3x that because every row allocates a fresh layer and frees
the previous one**, and after 32 rows of multi-megabyte churn the allocator is
holding the difference. That reframes the standing "we use 1.4–4.9x lrcalc's
memory" line: on live data the gap is far smaller, and most of what was being
compared is retention.

### ⚠️ Tried the obvious fix; it made things worse

Carrying the layer tables across rows and `clear()`ing them — keeping
capacity instead of reallocating — was implemented and **reverted**. It did
exactly what it was supposed to and still lost:

| `[16,13,10,7]²` | live heap | peak RSS | retention |
|---|---|---|---|
| as shipped | 123.2 MB | 332.2 MB | 2.7x |
| pooled | 417.6 MB | **492.2 MB** | 1.2x |
| pooled (tables only) | 399.2 MB | 547.0 MB | — |

Retention fell from 2.8x to 1.2x as predicted. RSS still rose 45%, because
**pooling pins each of ~110 buffers at its own high-water mark**, and the sum of
per-buffer peaks is much larger than the peak of the sum. The allocator was
doing the better job: a block freed by one row can be handed to a differently-
shaped request in the next, which a dedicated pool by construction cannot do.
Isolating the two halves showed the tables, not the layer vectors, were
responsible (399 MB with only the tables pooled).

So the retention is real but it is *not* free to reclaim, and the naive reading
— "reuse the allocations" — is wrong here. Worth knowing before anyone tries it
again.

What is still untried, and is a genuine reduction rather than a reshuffle: the
per-shard entry buffers and the output vector are both live at once during the
merge, holding the layer twice. Having the shards write into disjoint ranges
of a single output vector would remove one full copy — roughly 58 MB on this
shape — and helps the current code, pooled or not.

Two incidental findings:

* **Parallelism is memory-neutral in bytes** (86.4 vs 84.5 MB, 123.2 vs
  122.6 MB) even though it holds **2.25x more live entries** before the merge,
  because a hundred small shard tables carry less absolute slack than one giant
  power-of-two table.
* That 2.25x was invisible until `PEAK_LIVE_STATES` was corrected. It sampled
  `cur.len() + next.len()` *after* the merge, so it could not see pre-merge
  duplication at all, and reported the sharded path as free. Sharding silently
  invalidated the metric's documented meaning — the counter now samples the true
  pre-merge total.

## Next, in priority order

1. **Parallelism.** Deliberately deferred until after the memory work
   (per-thread layers multiply residency); now that bytes-per-state is
   ~4× smaller, a row-parallel merge is the next big lever.
2. ~~**Few-row factors below the counting crossover**~~ **Done 2026-07-31**,
   by exactly the route this item named: a cheaper fibre count (packed state,
   window-form inner loop) lowered the crossover to n ≥ 48, and the whole
   three-row band plus the ℓ(ν) = 3 asymmetric family now measures ahead of
   lrcalc — 1.06–1.33x where it was 0.71–0.88x. See "the wide-band deficit
   closed" above; the AC re-validation the first version of this item asked
   for ran the same day and confirmed both the bound and the sweep (1.02–1.49x
   against lrcalc, interleaved). What remains is one widening question:
   `rows ≤ 5` → 6 has measured 1.24–1.32x in-process three times but still
   has no out-of-process number.
3. **Extend counting to four-row factors.** The state gains one dimension per
   strip, so ℓ(ν) = 4 needs (λ¹ⱼ, aⱼ, bⱼ, cⱼ). Whether that stays affordable
   is unknown — the three-row case cost 26–2192 ops/term against a predicted
   "thousands, hopeless", so the bounding-box estimate is not trustworthy here
   and it should be measured rather than reasoned about. `[24,20,16,12]²`, our
   largest case, has four-row factors. Two lessons from 2026-07-31 apply: the
   packed-state decode trick is worth ~2x before any algorithm work, and
   one-traversal alternatives are now ruled out in both scan directions, so
   the fibre count is the only lane.
4. **Shape preprocessing** — factoring a skew diagram into connected
   components and expanding each separately, since the expansion of a
   disconnected shape is the product of its pieces.
5. ~~**Output residency.**~~ **Done.** On `[24,20,16,12]²` a growing share of
   peak RSS was the 5.3M-term *output* (the memoized `Arc<Vec>` plus the
   caller's clone), not the layer. `expand_skew_shared` returns the `Arc`;
   every in-crate caller only iterates, so none of them copy any more. Measured
   on `[8,7,6,5,4,3]²` (`heapstat skew-clone`), the clone alone was **164 041
   allocations and 14.1 MB per call** on top of the identical 14.1 MB in the
   cache. `SkewLr::lr_coeff` was the worst case — an entire expansion copied to
   read one coefficient. See [memory.md](memory.md).
