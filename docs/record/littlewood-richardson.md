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

This table carried neither a date nor a control, and silently understated the
library by up to 6× for an unknown period. It now carries both.

A **second copy of this table** — the same five shapes, the older numbers —
sat in `AutoLr`'s rustdoc, where it had neither date nor control and no reason
to be re-run when this one was. It is gone; this is the only copy. What the
rustdoc keeps is the conclusion that outlives the numbers: `SkewLr` wins at
every size measured, including every pair with |μ|+|ν| ≤ 12, so there is no
crossover left to dispatch on.

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
prioritization: it is one defect with two symptoms, not two. It also means the
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
this reason but have not yet been run out-of-process. In-process, against
`SkewLr` on the same terms:

```text
  s(4⁴)·s(4⁴)      70 terms   0.099ms →  0.010ms   10x
  s(8⁵)·s(8⁵)    1287 terms   2.06ms  →  0.116ms   18x
  s(10⁸)·s(6⁴)    210 terms   0.83ms  →  0.016ms   52x
  s(12⁶)·s(12⁶) 18564 terms  40.9ms   →  1.09ms    38x
  s(14⁷)·s(14⁷)116280 terms 356ms     →  7.52ms    47x
```

A *single* coefficient gains far more, because the predicate is O(ℓ(λ)) and
replaces a whole search outright: one `c^λ_{μν}` with μ = ν = (12⁶) goes from
467 ms to 1.1 µs. This table lived in `AutoLr`'s rustdoc, and this file pointed
at `src/rect.rs` for it — which never held it.

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
structure was wrong. Generation-stamped dense tables fixed it.

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
A dispatch bound is a measurement, and nothing in the tree re-checks it; this
is the second table in this file to be silently invalidated by drift, after the
undated baseline above.

**Negative result: checkpointing the fibre DP along the candidate DFS, with
difference-array transitions, loses 5–10x — measured, instrumented, reverted.**
The premises looked sound: row j of the fibre DP reads nothing of λ beyond
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
non-completing DFS branches before they can pay.

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

One edge per tableau means the DP enumerates the tableaux and pays hashing on
top of that. This is the same 1.0x compression the row layer measures at
ℓ(ν) ≤ 3, now seen from the transposed sweep direction, and the shared cause is
now plain: **any DP that produces all outputs in one traversal must carry
partial content in its state — the output is binned by content — and at three
rows the partial content pins the filling almost uniquely, so state-merging
cannot beat enumeration no matter which way the diagram is scanned.** Only
per-output counting escapes, because fixing λ turns content from state into
constraint. That closes both "one big traversal"
directions (rows: the layer-free prototype above; columns: this one) and leaves
the fibre count as the only approach that scales past enumeration here.

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
parallelized, this comparison must keep reporting a single-threaded number —
a multi-threaded wall-clock figure set against a single-threaded lrcalc would
conflate an algorithmic win with a hardware one, and note that a
parallel build pinned to one thread is not the same as a sequential build
(per-thread structures and merge machinery cost something even at N=1). lrcalc
being single-threaded is a property of its implementation, not of the problem;
its enumeration is at least as parallelizable as our layer, so threads are a
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
than a one-time experiment, and an example that nothing runs is not part of the
check.

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

### One coefficient: the naive backend is the right one

`lr_coefficient` at the Python boundary is deliberately `NaiveLr` rather than
`AutoLr`, which reads backwards until the shape of the cost is clear: a
targeted backtrack costs roughly the coefficient's own size, while `AutoLr`
builds a whole expansion and indexes into it.

```text
  c^[16,14,12,10,8,6]_{[8,7,6,5,4,3],[8,7,6,5,4,3]} = 1        5µs vs  771µs
  c^[24,20,16,12]_{[12,10,8,6],[12,10,8,6]}         = 1        6µs vs  145µs
  c^[13,12..2]_{[5,4,3,2,1],[12,11..3]}         = 14080     2446µs vs  347µs
```

So the naive search wins whenever the coefficient is small — the overwhelmingly
common case — and loses only when it is large, since it then enumerates that
many tableaux. This table was in `python.rs`'s docstring, which is the Sage
user's `help()` output; the docstring keeps the rule and this file keeps the
measurement.

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

## 2026-08-18: the allocator costs 8%, and every way of paying less costs more

Sample profiling put a number on what the traversal spends inside the
allocator, which until now was only known by its symptoms. `sample` at 1 ms on
`[22,18,14,10]²` (2 841 490 terms, 14.9M peak states, `examples/lrmem.rs`,
**on battery** — the shares below are within one profile, so power state does
not enter, but none of the wall times in this section are comparable to the AC
figures above):

| where the busy samples are | share |
|---|---|
| symfn itself | 83.4% |
| `libsystem_platform` (memmove, memset, memcmp) | 9.4% |
| `libsystem_malloc` | **8.2%** |
| `libsystem_kernel` | 0.3% |

So 8.2% is the ceiling on any allocation work at all, down from the 38% that
motivated the inline `Key` (see [`skew_lr::Key`](../../src/skew_lr.rs)). Of
that 8.2%, 37% enters from `fill_runs` — the innermost loop, which the
reference says allocates nothing.

**It does allocate, on exactly the shapes that matter.** Counting inline
against spilled keys at the insert site:

| shape | longest key | states spilling to `Key::Heap` |
|---|---|---|
| `[10,8,6,4]²` | 14 B | 0 |
| `[8,7,6,5,4,3]²` | 17 B | 0 |
| `[16,13,10,7]²` | 33 B | 3.4% |
| `[20,16,12,8]²` | 41 B | **50.8%** (18.9M of 37.2M) |

`INLINE = 30` holds for every shape small enough to measure quickly and fails
on the ones that take minutes, so the property was never false where anyone
looked. The cause is orientation: `prefer_conjugate` fires on the large shapes,
and the conjugate of `[20,16,12,8]²` has 40 rows, so the content prefix of a
key is 40 elements before the clipped row is appended.

### ⚠️ Widening the key to remove those 18.9M allocations made it 52% slower

Interleaved, four passes each, `[20,16,12,8]²`, same binary except `INLINE`:

| `INLINE` | `size_of::<Key>()` | map entry | spill rate | time |
|---|---|---|---|---|
| 22 | 24 B | 32 B | 99.8% | 5.67s |
| **30** | **32 B** | **40 B** | **50.8%** | **2.88s** |
| 38 | 40 B | 48 B | 1.1% | 4.43s |
| 46 | 48 B | 56 B | 0% | 4.38s |

Removing *every* allocation (`INLINE = 46`) costs half again as much time as
leaving half of them in. The layer is the working set, the traversal is bound
by how much of it fits in cache, and 8 more bytes per entry outweighs a malloc
and a free per new state. The shipped value sits near the optimum by accident,
not by design — it was chosen to make `Key` 32 bytes.

The reading to take from this: **entry size, not allocation, is the binding
constraint**, and the two are traded against each other. It also explains the
pooling failure above without appealing to allocator internals.

### ⚠️ Pre-sizing the layer tables from the previous row's growth: no effect

Layers grow 1.5–8x per row (`SKEW_TRACE=1`), and each shard table is sized
from the *input* layer, so it doubles two or three times while filling and
rehashes each time. Predicting the next row's size from the last row's ratio
(clamped at 4x) was implemented and **reverted**: 2.89s vs 2.86s and 1094 MB
vs 1124 MB peak RSS, both inside the run-to-run spread, and the allocator's
share of the profile moved 8.2% → 8.1%. Table growth is amortized and the
re-inserts are on hot lines; there was nothing there to win.

### What did land

One allocation per output term, on conjugated shapes only. The output map
built a `Partition` from the decoded content and immediately dropped it for
its conjugate; [`partition::conjugate_parts`](../../src/partition.rs)
transposes the parts in the buffer instead, so one `Vec` is allocated per term
where two were. On `skew-big` (`examples/heapstat.rs`) that is **499 552 →
335 516 allocations**, a third of the total, at 164 037 terms — and the
`skew-big` budget in [workloads.rs](../../src/measure/workloads.rs) was lowered
to 340 000 so the old number cannot come back unnoticed. Wall time is unchanged
(2.86s vs 2.89s on `[20,16,12,8]²`, inside the spread), which the 8.2% ceiling
predicts.

`examples/lrmem.rs` now takes the shape on the command line. The shapes that
expose any of this — `[20,16,12,8]²` and up — are larger than the four presets
it used to carry, and `lrheap` keeps the presets and the allocator wrapper.

## 2026-08-18, later: the layer key is two lattice-path bitmaps, 1.14–2.10x

The `INLINE` sweep above said the layer is bound by entry size and that the
byte key could not be made narrower because the keys themselves were 33–49
bytes on the shapes that matter. Both halves of a state are monotone
sequences of bounded integers, and such a sequence is a set of distinct bit
positions: `content` is a partition with at most `rows` parts each at most
`width` (the 1's form a horizontal strip), so part `i` of `k` sits at bit
`content[i] + (k − 1 − i)`; the clipped row above is weakly increasing with
values at most `k + 1`, so cell `j` sits at bit `above[j] + j`. Every
position is below `rows + width`, so a state is two `u64`s whenever
`rows + width ≤ 64` — the conjugate walk of `[20,16,12,8]²` (40 × 8) packs
its 41-byte key into 16 bytes — two 128-bit words up to 128, and the byte
`Key` only past that. A layer entry is 24 bytes where it was 40, nothing
spills, and hashing is two word rounds. [`skew_lr::PackedKey`](../../src/skew_lr.rs)
is the type; `LayerKey` is the trait the byte `Key` now also implements, and
`expand_oriented` picks the representation once per shape from `rows + width`.
The parallel merge also moved from a probe-then-insert to `entry`, one hash
per key instead of two on every key new to the accumulator.

Verified: every existing `skew_lr` test runs through the bitmap key (all of
them are small enough), and new tests pin the encoding as a bijection over a
6 × 6 box and 126 rows, force all three representations onto the same merging
shapes and compare them with each other and with `NaiveLr`, and place a
content bit at position 63 and 127 exactly (`[62,3]/∅`, `[126,3]/∅`) with the
shapes one past each boundary dispatched to the next representation. Both
oracle suites pass. `lr_cli` output is byte-identical to the previous binary
on `[16,13,10,7]·[8,6,4,2]` and `[20,16,12,8]²`.

Measured **on battery (31% → 24%)**, so nothing here is comparable to the AC
tables; every row is an interleaved out-of-process A/B of the HEAD binary
against this one (`examples/bench_shapes.rs`, min of 3, min of 2 on the
largest), which is the comparison that survives the power state:

| shape | before | after | | terms | peak states |
|---|---|---|---|---|---|
| `[12,10,8]²` | 3.30 ms | 2.70 ms | 1.22x | 6 579 | 9 400 |
| `[20,16,12]²` | 61.3 ms | 46.6 ms | 1.32x | 64 335 | 215 131 |
| `[16,13,10,7]²` | 356 ms | 296 ms | 1.20x | 390 075 | 1.96M |
| `[8,7,6,5,4,3]²` | 127 ms | 104 ms | 1.22x | 164 037 | 1.20M |
| `[9,8,7,6,5]²` | 62.6 ms | 54.5 ms | 1.15x | 105 533 | 628 134 |
| `[8,7,6,5,4]²` | 31.6 ms | 27.6 ms | 1.14x | 45 791 | 213 500 |
| `[6,5,4,3,2,1]²` | 8.2 ms | 7.1 ms | 1.15x | 10 873 | 26 489 |
| `[7⁵]²`, `[3¹²]²` | 0.4 ms | 0.4 ms | tie | | |
| `[20,16,12,8]²` | 2.83 s | 1.60 s | **1.77x** | 1 393 833 | 6.90M |
| `[22,18,14,10]²` | 8.30 s | 3.96 s | **2.10x** | 2 841 490 | 14.97M |
| `[32,26,20]²` (128-bit band, `rows + width` 70) | 2.36 s | 1.91 s | 1.23x | 568 289 | 5.26M |
| `[70,66]²` (byte band, 144) | 12.9 ms | 11.7 ms | 1.10x | 20 234 | 20 569 |
| `[33,31]²` (128-bit band) | 0.7 ms | 0.7 ms | floor | 2 576 | 2 672 |

Peak live states are unchanged on every row, as they must be — the states are
the same, only their bytes moved — so the whole gain is bytes per state, and
it grows with the layer exactly as the `INLINE` sweep predicted: the two
shapes where half or more of the states spilled to `Key::Heap` gain most.
Peak RSS on `[20,16,12,8]²` (`/usr/bin/time -l`, one run each): 1113 MB →
1013 MB. That is the ~110 MB the entry narrowing accounts for (6.9M states ×
16 bytes); the rest of the resident set is the output and allocator
retention, per the "two thirds of RSS is retention" section above.

Not re-run: `[24,20,16,12]²`, on battery. Its conjugate walk is 48 × 8, in
the `u64` band with the same key shape as `[22,18,14,10]²`, so the direction
is expected to carry, but that is an expectation and this file does not
record expectations as results. The `prefer_conjugate` thresholds were
calibrated when key length depended on orientation (a direct walk of a wide
shape carried its 40-cell row at a byte a cell); with the bitmap it no longer
does, so the rule may now be conservative and should be re-measured before it
is next relied on.

## 2026-08-18, later still: product orientation — direct for asymmetric pairs, 1.14–1.92x

A product `s_a·s_b` is one skew expansion of the juxtaposed shape, and four
walks compute it: either factor can be the enumerated block (the other is the
ballot offset, its own block having one canonical filling), on the diagram or
on its transpose — and transposing the juxtaposition swaps the roles, so the
four walks enumerate `b`, `a`, `a'` and `b'`. Until today the enumerated
factor was chosen by *lexicographic* order (`mu >= nu` put the lex-larger on
top), which is unrelated to cost, and the transpose by `prefer_conjugate` on
the juxtaposed shape, whose thresholds predate the bitmap key. The lrcalc
sweep is almost entirely squares, where the first choice is moot, which is
how this stayed invisible.

`examples/calibrate_orientation.rs` times all four walks of a pair
in-process (order rotated, min of `reps`, every walk's expansion checked equal
to the library's), and reports each walk's **productions** — row fillings
committed to the layer, merged or not, now counted by
`skew_lr::take_productions` — and peak states. Two grids, 33 pairs, on
battery (14–12%), so ratios only:

| product | walk `b` | walk `a` | walk `a'` | walk `b'` | library took | best |
|---|---|---|---|---|---|---|
| `[16,13,10,7]·[8,6,4,2]` | 9.9 ms | **9.8** | 21.0 | 15.6 | `a'` (2.15x) | direct |
| `[20,16,12,8]·[8,6,4,2]` | 12.8 | **11.8** | 31.2 | 19.6 | `a'` (2.65x) | direct |
| `[18,15,12,9]·[9,7,5,3]` | **24.4** | 32.8 | 39.6 | 28.9 | `a'` (1.62x) | direct |
| `[24,20,16,12]·[10,8,6,4]` | **61.1** | 63.3 | 92.2 | 74.2 | `a'` (1.51x) | direct |
| `[20,16,12,8]·[10,8,6,4]` | 58.3 | **51.5** | 76.5 | 67.3 | `a'` (1.49x) | direct |
| `[16,13,10,7]·[5,4,3,2,1]` | 2.19 | **1.86** | 5.34 | 2.80 | `a'` (2.87x) | direct |
| `[20,16,12,8]·[6,5,4,3,2,1]` | **12.2** | 12.4 | 32.0 | 16.9 | `a'` (2.62x) | direct |
| `[16,13,10,7]·[10,8,6,4]` | 42.6 | **36.4** | 51.5 | 44.8 | `a'` (1.42x) | direct |
| `[12,10,8,6]·[10,8,6,4]` (64 cells) | **22.2** | 23.1 | 29.1 | 26.7 | `a'` (1.31x) | direct |
| `[12,10,8,6]·[9,8,7,6,5]` (ratio 0.97) | **36.5** | 44.4 | 46.7 | 44.0 | `a'` (1.28x) | direct |
| `[16,13,10,7]·[10,9,8,7,6,5]` | **324** | 572 | 390 | 324 | `a'` (1.20x) | direct |
| `[16,13,10,7]·[12,10,8,6]` | 149 | 137 | 104 | **103** | `a'` (1.01x) | transpose |
| `[10,9,8,7,6,5]·[9,8,7,6,5,4]` | 954 | 874 | **614** | 642 | `a'` | transpose |
| `[10,9,8,7,6,5]·[8,7,6,5,4,3]` | 364 | 313 | **272** | 287 | `a'` | transpose |
| `[16,13,10,7]²` | 655 | 650 | 291 | **290** | `a'` | transpose |
| `[10,9,8,7,6,5]²` | 2063 | 2065 | **1267** | 1272 | `a'` | transpose |
| `[12,10,8,6]²` (72 cells) | 51.6 | 52.7 | **45.7** | 46.0 | `a'` | transpose |
| `[8,7,6,5,4,3]²` (66) | 108 | 106 | 104 | **103** | `a'` | tie |
| `[9,8,7,6,5]²` (70) | 54.1 | 54.2 | **52.2** | 52.3 | `a'` | tie |
| `[12,9,6,3]²` (60) | **21.7** | 21.9 | 30.1 | 30.1 | `a'` (1.39x) | direct |
| `[11,9,7,5]·[10,8,6,4]` (60) | **16.5** | 17.7 | 24.4 | 22.2 | `a'` (1.49x) | direct |
| `[9,8,7,6,5]·[9,7,5,3,1]` (60) | **19.6** | 24.7 | 28.5 | 28.1 | `a'` (1.45x) | direct |

(`[8,7,6,5,4]²`, `[10,8,6,4]²`, `[12⁶]²`, and nine further asymmetric pairs
in the harness's default list read the same way and are omitted for space.)

**What the productions column says.** Time is productions times a
per-production cost, and both move with the walk. Transposing pays only
where the direct fill is loose enough to leave compression on the table:
`[16,13,10,7]²` commits 40.2M productions directly (103 per term) and 8.9M
transposed, at roughly twice the cost each because the content is longer
(up to `a₁ + b₁` parts against `ℓ(a) + ℓ(b)`), net 2.25x. An asymmetric pair
is already constrained by its small offset — 8–30 productions per term
directly, `[16,13,10,7]·[8,6,4,2]` at 200k for 24k terms — so the transpose
has little to compress and pays its per-production premium for nothing:
1.5–2.9x slower. That is the regime the lex rule plus `prefer_conjugate` was
putting every asymmetric pair into. Between the two, near-squares of equal
row count and comparable size still want the transpose (1.06–1.55x), and the
60-cell squares that used to fire the shape rule now prefer the direct walk
by 1.3–1.5x — the bitmap key removed the direct orientation's key-length
penalty, so the crossover moved up.

**The rule that landed** — `skew_lr::product_walk`, replacing both
`mu >= nu` and the shape rule for products; general skew expansions keep
`prefer_conjugate` unchanged. The larger factor by cells (then lex, so the
pair is a total order and both argument orders share one cache entry) goes
on top and the smaller is enumerated; the walk transposes iff the factors
have the same number of rows, the smaller is at least 0.7 of the larger by
cells, the product has at least 72 cells, and the juxtaposed shape has at
least eight rows and is wider than tall. Fitted to the 33 pairs above; the
tests pin 24 of them. What the fit leaves on the table is second-order: the
choice between enumerating `a` and `b` in the direct regime is within 1.21x
of the best everywhere for the smaller factor and 1.34x for the larger, with
no clean predictor in the data (row width and offset flatness both matter and
pull against each other), so the smaller factor is enumerated as the plainer
rule.

Out-of-process, `lr_cli mult`, interleaved against the previous binary, min
of 3, outputs identical on every row:

| product | before | after | |
|---|---|---|---|
| `[16,13,10,7]·[8,6,4,2]` | 28.6 ms | 16.5 ms | **1.73x** |
| `[20,16,12,8]·[8,6,4,2]` | 40.0 | 20.8 | **1.92x** |
| `[20,16,12,8]·[6,5,4,3,2,1]` | 42.1 | 21.9 | **1.92x** |
| `[16,13,10,7]·[5,4,3,2,1]` | 8.7 | 6.0 | 1.45x |
| `[18,15,12,9]·[9,7,5,3]` | 52.4 | 37.1 | 1.41x |
| `[24,20,16,12]·[10,8,6,4]` | 121 | 89 | 1.36x |
| `[14,12,10,8,6]·[7,5,3,1]` | 22.0 | 16.4 | 1.34x |
| `[16,13,10,7]·[8,7,6,5,4,3]` | 148 | 111 | 1.33x |
| `[11,9,7,5]·[10,8,6,4]` | 34.0 | 26.0 | 1.31x |
| `[12,9,6,3]²` | 41.2 | 33.1 | 1.25x |
| `[20,16,12,8]·[10,8,6,4]` | 104 | 84 | 1.24x |
| `[12,10,8,6]·[10,8,6,4]` | 39.8 | 33.2 | 1.20x |
| `[12,10,8,6]·[9,8,7,6,5]` | 66.7 | 56.8 | 1.18x |
| `[16,13,10,7]·[10,9,8,7,6,5]` | 523 | 459 | 1.14x |
| `[12,11,10,9,8,7]·[6,5,4,3]` | 22.6 | 19.8 | 1.14x |
| `[16,13,10,7]·[12,10,8,6]`, `[10,9,8,7,6,5]·[8,7,6,5,4,3]` | | | 0.99–1.00x (same walk) |
| `[16,13,10,7]²`, `[12,10,8,6]²`, `[8,7,6,5,4,3]²`, `[10,9,8,7,6,5]²` | | | 0.98–1.01x (same walk) |

The two `prefer_counting` routes and the rectangle path sit in front of
`SkewLr` in `AutoLr`, so nothing here touches a three-row or rectangular
product.

## 2026-08-18, last: skew-shape transposition follows the outer shape's steps; the direct-regime tiebreak is closed

Two items the product-orientation section left open, both measured on
battery (11%), out of process through `lr_cli skew` with `SKEW_ORIENT`
forced, interleaved, min of 3, outputs identical.

**`prefer_conjugate` for genuine skew shapes.** Sixteen shapes shaped like
the rule's remaining callers — coproduct-style λ/μ with λ large, and
`lr_coeff`-style λ over the larger factor — all with `rows ≥ 8`,
`width > rows`, and 60–99 cells, so all firing the old rule:

| shape | cells | direct | transposed | direct/transposed |
|---|---|---|---|---|
| `[16,15,…,9]/[8,7,…,1]` | 64 | 16.4 ms | 8.2 ms | **2.00** |
| `[14,13,…,7]/[6,5,…,1]` | 63 | 8.4 | 5.9 | 1.41 |
| `[16,15,…,7]/[6,5,…,1]` | 94 | 46.4 | 36.9 | 1.26 |
| `[14,13,…,5]/[6,5,…,1]` | 74 | 38.2 | 31.6 | 1.21 |
| `[12,11,…,3]/[5,4,3,2,1]` | 60 | 13.1 | 11.8 | 1.11 |
| `[12,11,…,3]/[4,3,2,1]`, `/[3,2,1]`, `[13,12,…,2]/[5,4,3,2,1]`, `[15,…,6]/[3,2,1]`, `[12,12,11,11,10,10,9,9]/[6,6,5,5]`, `[10⁸]/[4⁴]` | 60–99 | | | 0.95–1.04 |
| `[18,16,…,4]/[8,6,4,2]` | 68 | 25.5 | 28.7 | 0.89 |
| `[16,14,…,2]/[6,4,2]` | 60 | 5.1 | 7.0 | 0.73 |
| `[24,21,18,15,12,9,6,3,1]/[16,13,10,7]` | 63 | 103 | 153 | 0.67 |
| `[30,26,22,17,13,9,5,3,1]/[20,16,12,8]` | 70 | 217 | 373 | **0.58** |
| `[20,17,14,11,8,5,2,2]/[8,5,2]` | 64 | 8.9 | 15.5 | **0.57** |

Cell count does not separate the two groups, and neither does the average
coefficient (the two `lr_coeff`-style shapes have the largest, 10⁵–10⁶
tableaux per term, and want the direct walk). What does is the outer shape's
descent: every λ that steps down by one cell per row wants the transpose,
every λ that steps by two or more wants the direct walk. The productions
counter says why (`take_productions`, in-process, min of 2):

| shape | direct prod. | ns each | transposed prod. | ns each |
|---|---|---|---|---|
| `[16,15,…,9]/[8,…,1]` | 332 556 | 41.5 | **105 846** | 51.4 |
| `[14,13,…,7]/[6,…,1]` | 143 641 | 51.9 | **66 204** | 58.7 |
| `[16,15,…,7]/[6,…,1]` | 1 574 295 | 23.9 | **826 833** | 35.1 |
| `[16,14,…,2]/[6,4,2]` | **56 088** | 43.3 | 59 469 | 76.5 |
| `[20,17,14,11,8,5,2,2]/[8,5,2]` | **127 067** | 42.5 | 148 917 | 83.7 |
| `[24,21,…,1]/[16,13,10,7]` | **3 927 441** | 21.8 | 6 610 534 | 20.8 |
| `[30,26,22,…,1]/[20,16,12,8]` | **11 773 231** | 16.0 | 15 402 288 | 21.9 |
| `[18,16,…,4]/[8,6,4,2]` | 968 888 | 20.8 | **474 924** | 50.6 |

On the step-one staircases the direct rows overlap almost entirely, the
direct layer barely merges, and the transpose commits 1.3–3.1x fewer
fillings at ~1.3x the cost each: a win. On the steep shapes the transpose
commits as many or more — the compression the rule assumed is not there —
and pays 1.3–2x per filling for its longer content: a loss, and it grows
with the shape. (`[18,16,…,4]/[8,6,4,2]` is the one shape where the
transpose does compress, 2x, and still loses on the 2.4x cost each.)

**Landed:** `prefer_conjugate` now also requires every step of `outer` to
be at most one. That keeps the transpose exactly on the shapes it was
calibrated on and measured to win, and turns it off where it was measured
to lose, up to 1.75x. The remaining thresholds are unchanged. A shape with
mixed steps gets the direct walk, which is the untransposed default and not
a measured loss anywhere; a finer rule wants a productions model this data
does not supply.

**The direct-regime tiebreak, closed without a rule.** Across the 25
asymmetric pairs in the calibration harness, enumerating the smaller factor
is within 1.21x of the best walk everywhere and the larger within 1.34x, in
different places. A perfect predictor over "always the smaller" would gain
21% on one pair (`[10,9,8,7,6,5]·[7,6,5,4]`), 9–18% on seven, and nothing
on the rest — about 4% on average — and the two components pull against
each other: enumerating the larger factor commits fewer fillings (its
offset is flatter, so the fill is more constrained; 1.6–4x fewer) but pays
1.1–4x more per filling (its rows are wider, so the run fill visits more
partial runs per completed row). `[18,15,12,9]·[9,7,5,3]` and
`[20,16,12,8]·[10,8,6,4]` are the same shape family and fall on opposite
sides. Not worth a fitted rule; the smaller factor stays enumerated.

## 2026-08-18, mixed-step outer shapes: the transpose follows the inner shape's depth, up to 6x

The step-one rule above left mixed-step outer shapes on the direct walk
unmeasured. Measured (AC, charging at 11–12%; `lr_cli skew` with
`SKEW_ORIENT` forced, out of process, interleaved, min of 3, outputs
identical), 41 further shapes in five sweeps: one big step among ones, a few
ones among big steps, alternating, zero steps among big ones, product terms
as λ, wide bands, and `lr_coeff`-style shapes with the inner shape four rows
short of the outer. Neither the step pattern nor the cell count nor the
average coefficient predicts the winner; the depth of the inner shape does —
**gap = ℓ(outer) − ℓ(inner)**, the number of full rows under the inner
shape's last row, which the direct walk meets last and the transposed walk
meets first:

| gap | transposed wins | direct wins | ties |
|---|---|---|---|
| 0–2 | 6.06, 3.56, 2.33, 2.22, 1.64, 1.47, 1.38, 1.12, 1.08 (`[26,23,20,17,14,11,8,5,2]/[14,12,10,8,6,4,2]` is the 6.06x, 1810 → 299 ms) | 0.77 (`[24,20,16,12,8,4,3,2]/[9,7,5,3,1,1]`) | 0.98 |
| 3 | 1.85, 1.08 | 0.95 | 1.03 |
| 4 | 3.95 (`[26,23,20,17,14,11,8,5,2]/[12,10,8,6,4]`, 2316 → 586 ms), 1.26, 1.25, 1.24, 1.21, 1.10 | 0.89, 0.86, 0.84, 0.83, 0.78 | 1.04, 1.01 |
| ≥ 5 | 1.11 (a step-one staircase) | 0.57, 0.58, 0.63, 0.67, 0.70, 0.73, 0.73, 0.76, 0.81, 0.84, 0.87, 0.89 | six |

Below the size floor the picture is a coin flip with a bad tail
(`[24,20,16,12,8,4,2]/[12,10,8,6,4,2]`, seven rows, transposes 2.2x
*slower*), so `rows ≥ 8`, `width > rows`, `cells ≥ 60` stand. Two a-priori
models were tried against the 52 shapes and both fail: the overlap between
consecutive rows (what the state carries; the transposed walk's is smaller
almost everywhere, yet it loses on many) and a per-row state estimate
`Σ_r p_{≤r+1}(cells so far) · C(overlap_r + r, overlap_r)` (its value-range
factor calls everything direct; the ballot condition constrains far more
than it knows). The rule is therefore empirical: **transpose iff the size
floor holds and either every step of `outer` is at most one, or
`gap ≤ 4`.** Gap 4 is a genuine coin flip (six wins against five losses); it
transposes because the wins sit on the slow cases (3.95x on 2.3 s) and the
losses on 60–130 ms ones, and the calibration criterion, stated by the
project's owner this session, is that a large win on a slow case outweighs
a small loss on a fast one. Gap 5 loses on the slow cases too
(`[30,26,22,17,13,9,5,3,1]/[20,16,12,8]` 0.58 at 217 ms) and stays direct.

Out of process against HEAD (the step-one rule), same discipline:

| shape | HEAD | band rule | |
|---|---|---|---|
| `[26,23,20,17,14,11,8,5,2]/[14,12,10,8,6,4,2]` | 1751 ms | 295 ms | **5.95x** |
| `[26,23,20,17,14,11,8,5,2]/[12,10,8,6,4]` | 2444 | 573 | **4.27x** |
| `[20,18,16,14,12,10,8,6]/[10,8,6,4,2,1]` | 242 | 72 | **3.38x** |
| `[20,15,14,13,12,11,10,9]/[8,7,6,5,4,3,2,1]` | 36.9 | 16.4 | 2.24x |
| `[16,15,14,12,11,10,9,8]/[7,6,5,4,3,2,1]` | 28.8 | 13.5 | 2.13x |
| `[16,14,12,10,9,8,7,6]/[6,5,4,3,2,1]` | 21.8 | 14.0 | 1.56x |
| `[30,26,22,17,13,9,5,3,1]/[20,16,12,8,4,2,1]` | 322 | 211 | 1.52x |
| `[16,15,14,13,9,8,7,6]/[6,5,4,3,2,1]` | 19.8 | 14.5 | 1.37x |
| `[20,18,16,14,12,10,8,6]/[10,8,6,4]` | 77.6 | 65.6 | 1.18x |
| `[18,16,14,12,10,8,6,4]/[8,6,4,2]` | 27.5 | 31.9 | 0.86x |
| `[28,24,20,16,13,9,5,3]/[20,16,12,8]` | 97.9 | 125.6 | 0.78x |
| `[24,20,16,12,8,4,3,2]/[9,7,5,3,1,1]` | 33.1 | 43.1 | 0.77x |
| step-one staircases; steep shapes with gap ≥ 5 | | | 0.99–1.00x (same walk) |

The product rule (`product_walk`) is untouched: a juxtaposed product shape
has gap = ℓ(smaller factor), and its top block's single canonical filling
makes it a different problem — `[20,16,12,8]·[8,6,4,2]` has gap 4 and wants
the direct walk by 2.65x there.

## 2026-08-18, the consumer side: one cache entry per product on every route, and `Schur::mul` without the copies

Three costs sat between the engine's memoized expansion and a caller.
`AutoLr::schur_product` memoized only its `SkewLr` branch: a rectangle,
two-row or three-row product was recomputed on every call, and
`AutoLr::lr_coeff`'s peek at the product cache — the route that answers a
sweep of coefficients from an expansion already built — looked under a key
those routes never wrote, so every coefficient off a counting-route product
ran its own λ/μ traversal. `Schur::mul_with` went through the trait's owned
`schur_product`, a deep clone of the memoized vector — the one caller
`expand_skew_shared` left behind, and the main entry point, since
`schur_multiply` at the Python boundary is `Schur::mul`. And
`SymFn::add_term` copied every key it was handed, `map.entry(p.clone())`,
to serve the rare cancel-and-remove path. Per output term: two allocations
and two frees around one map descent, with a third live copy of the product
for the length of the loop.

**Built.** `LrBackend::schur_product_shared`, provided as
`Arc::new(self.schur_product(..))` and overridden by every memoizing backend
to hand out its cached vector. `skew_lr::memoized_product(mu, nu, shortcut)`:
the one place a product enters the skew table, under the shape
`product_walk` chooses; `AutoLr` passes its closed form and counting routes
as the shortcut and `SkewLr` passes none, so every route stores under the
entry the peek reads. `Schur::mul_with` reads each pair's expansion in
place; the first pair's terms — sorted by λ, distinct — build the `BTreeMap`
in one pass through an exactly-sized vector (`BTreeMap::from_iter` takes a
vector's buffer as it is; fed an iterator of unknown length it grows one by
doubling, which cost 3 MB of peak on the 164k-term case below), and later
pairs accumulate by reference, copying a partition only when its term is
new. `add_term` goes through `Entry` and never copies its key.

**Measured** (AC, charging at 42%; `examples/bench_schur_mul`, before =
`68769e4` built with the same harness, out of process and interleaved, min
of 5 in-process reps, three rounds for the consumer cases and two for
dispatch; every product warm, so only the consumer side is timed):

| `Schur::mul`, products warm | terms | HEAD | now | |
|---|---|---|---|---|
| `s_μ·s_μ`, μ = `[10,8,6,4]` | 23 973 | 2.95 ms | 0.82 ms | 3.6x |
| `[12,10,8,6]` | 79 241 | 10.4 | 2.62 | 4.0x |
| `[8,7,6,5,4,3]` | 164 037 | 27.5 | 6.5 | 4.2x |
| `[16,13,10,7]` | 390 075 | 58.7 | 16.6 | 3.5x |
| `(s_μ + s_ν)·s_μ`, `[10,8,6,4]`, `[11,8,5,4]` — the second pair puts 2 730 of its 26 703 terms on partitions the first did not | 26 703 | 5.67 | 2.60 | 2.2x |
| `[8,7,6,5,4,3]`, `[10,7,6,5,3,2]` — 63 970 of 228 007 new | 228 007 | 59.9 | 41.8 | 1.43x |
| `(Σ_{λ⊢8} s_λ)²`, 484 pairs | 231 | 1.13 | 0.34 | 3.4x |
| n = 10, 1 764 pairs | 627 | 10.1 | 3.17 | 3.2x |
| n = 12, 5 929 pairs | 1 575 | 80.7 | 26.3 | 3.1x |
| `s_{21}^{10}` | 5 410 | 19.0 | 10.1 | 1.9x |
| `s_{321}^{6}` | 15 388 | 166.6 | 80.0 | 2.1x |
| `s_{42}^{6}` | 12 050 | 102.7 | 49.6 | 2.1x |

The owned `AutoLr::schur_product` — the deep clone by itself — is 4.2 ms on
the 164k-term product on either side, so of `Schur::mul`'s 27.5 ms, 23 ms
was the map and its key copies; the copy-free path is 6.5 ms all in.

Dispatch, each rep from an empty cache (`first` = the product, `second` =
the same call again, `sweep` = `AutoLr::lr_coeff` over every term of the
product, `sweep again` = once more; min of 3, two rounds):

| | terms | first | second | sweep | sweep again |
|---|---|---|---|---|---|
| `[6,6,6,6]·[5,5,5]` (rectangle), HEAD | 56 | 3 µs | 3 µs | 1 µs | 1 µs |
| now | | 3 µs | 1 µs | 1 µs | 1 µs |
| `[20,16,12]·[20,16]` (two-row route), HEAD | 7 909 | 3.62 ms | 3.55 ms | 1 456 ms | 11.2 ms |
| now | | 3.63 | 0.094 | 1.34 | 1.34 |
| `[16,12,8,4]·[12,10,8]` (three-row route), HEAD | 41 105 | 38.3 | 37.8 | 18 764 | 471 |
| now | | 38.2 | 0.53 | 8.46 | 8.34 |

The first call is unchanged — storing costs nothing measurable — and a
repeat is the copy (38x, 71x). The cold sweep is the item that mattered: on
HEAD it cost 400x and 490x the product it was reading, one λ/μ expansion per
coefficient; now 0.4x and 0.2x.

Memory (`heapstat schur-mul`, new: `Schur::mul` on `[8,7,6,5,4,3]²` with
the product warm, bit-exact across runs): peak 24.1 → 17.0 MB, allocations
355 418 → 178 959, total 30.7 → 22.0 MB. What remains is the result itself,
164k map entries each owning a partition, plus the build's temporary vector
at 32 bytes per term over the tree's ~70.

**Two choices decided by measurement**, three builds interleaved over three
rounds, every case agreeing to the millisecond between rounds:

- The first pair's terms could enter the map by sorted insertion instead of
  in one pass. That was **40.8 ms** on the 164k-term pair against 6.5 —
  slower than HEAD's 27.5, because a term new to the map costs the
  by-reference path two descents (`get_mut`, then `insert`), and on the
  first pair every term is new. The one-pass build has no descent at all.
- Later pairs could go through `add_term` with a copied key (one descent,
  one allocation and one free per term) instead of by reference (one
  descent, plus a second and the copy only when the term is new). By
  reference was **1.85x** faster on `(Σ_{λ⊢12} s_λ)²` (26.3 against 48.8 ms)
  and 1.5x on `s_{321}^6` (80 against 122), where nearly every term lands on
  a partition already present. The second-pair case with 28% of its terms
  new is where the two descents show, and it still gains 1.43x over HEAD.

**Not done.** The one-pass build's temporary vector is a transient peak
contributor at about 40% of the tree it builds; sorted insertion would
remove it at 6x on this stage — about 1% of the expansion's time on the
multi-million-term shapes where peak binds, and the whole of a warm repeat's
time everywhere else. A size threshold would serve both regimes; unmeasured
at that scale, and in the open tail. The Python boundary's own copies stand:
`dump` builds a `Vec<(Key, Coeff)>` from the map and PyO3 builds the list of
tuples from that, so `schur_multiply`'s peak is cache + map + terms + Python
objects, and the clone removed here was never the binding one there.

## 2026-08-18, a Pieri route in `AutoLr`: measured and declined

The proposal: `s_λ·s_{(k)}` and `s_λ·s_{(1^k)}` with λ not a rectangle go
to the layer (a rectangle times a row or column is Okada's form), and a
horizontal- or vertical-strip enumeration generates the answer with no LR
machinery under it — h → s and e → s in `convert.rs` already run that way,
off `AutoLr`. The question was what the layer costs on such a shape.

⚠️ **Measured on battery (84%)**, ratios only. Ad-hoc probe, not kept: a
scratch crate against the tree at `575c217`; per rep, `clear_caches()`, then
`SkewLr::schur_product`, then `AutoLr::schur_product`, then the direct route
— the strip enumeration of `convert.rs`'s `horizontal_strips` (on the
conjugate for a vertical strip), sorted into the product's form; min of 5;
all three equal on every case.

| product | terms | `SkewLr` | direct | ratio |
|---|---|---|---|---|
| `[8,7,6,5,4,3]·s_10` | 128 | 0.029 ms | 0.022 ms | 1.3x |
| `[16,13,10,7]·s_20` | 512 | 0.114 | 0.088 | 1.3x |
| `[20,16,12,8,4]·s_30` | 3 125 | 0.500 | 0.359 | 1.4x |
| `[30,20,10]·s_30` | 1 331 | 0.155 | 0.073 | 2.1x |
| `[80,50]·s_160` | 1 581 | 0.256 | 0.376 | 0.7x |
| `[10,9,…,1]·s_30` | 1 024 | 0.206 | 0.228 | 0.9x |
| `[15,14,…,1]·s_15` | 32 768 | 9.20 | 4.65 | 2.0x |
| `[10,9,…,1]·e_10` | 1 024 | 0.499 | 0.234 | 2.1x |
| `[12,11,…,1]·e_6` | 2 510 | 0.524 | 0.548 | 1.0x |

`s_1` and `e_1` read 3–8x, on 3–7 µs. Everything else is 0.7–2.1x, and the
layer's time is proportional to the output. Nothing is there to remove: the
partial fillings the layer carries on a one-row or one-column skew shape are
the LR fillings of a Pieri product of a prefix of λ, which is
multiplicity-free, so their number is a strip count, never a tableau count,
and the direct route wins only its constant factor. `two_row.rs`'s
`rows ≥ 3` clause records the same fact from the other side.

The consumer pattern where a route would show most — `Schur::mul` of a
many-term element by `s_k` or `e_k`, one cold `AutoLr` product per pair,
which is what a Sage `X * s[k]` crosses as — against one direct strip step
over the terms into a `BTreeMap`:

| X · s_k | pairs | `Schur::mul` | direct step | ratio |
|---|---|---|---|---|
| `Σ_{λ⊢12} s_λ · s_3` | 77 | 0.62 ms | 0.18 ms | 3.4x |
| `Σ_{λ⊢20} s_λ · s_5` | 627 | 4.62 | 2.40 | 1.9x |
| `Σ_{λ⊢20} s_λ · e_5` | 627 | 7.06 | 3.20 | 2.2x |

About 7 µs a pair against 4, and a route would keep the store
(`memoized_product`'s `Arc` and table insert), so it would land nearer
1.5–2x, on milliseconds. Not written: it would cost a strip enumerator
shared out of `convert.rs` or a third copy, a predicate, tests, and one more
shortcut under `memoized_product`, for at most 2x on products under 10 ms
and nothing on the products where time is spent.

## 2026-08-18, four-row factors: the layer is not enumeration there, and the counting bands have moved

Item 3 of the tail — extend the fibre count to four-row factors, for
`[24,20,16,12]²` — rested on two premises: that the case costs 148 s, and
that at four rows the layer enumerates tableaux the way it does at two and
three, so that only per-output counting escapes. Both were checked today.

⚠️ **Measured on battery (78–84%)**, one process per shape, in-process
timing (`examples/bench_shapes`, which now also reports productions and the
coefficient sum, the number of LR tableaux; `/usr/bin/time -l` for CPU and
RSS). The default path throughout — `SkewLr` through `product_walk`, which
transposes the four-row squares (eight-row juxtaposed shapes) and walks the
three-row ones directly, parallel fill on — so wall and CPU differ by the
parallel speedup:

| square | terms | peak states | productions | tableaux ÷ productions | wall | CPU | peak RSS |
|---|---|---|---|---|---|---|---|
| `[20,16,12]²` (three-row) | 64 335 | 215 042 | 2 614 952 | **1.07x** | 0.072 s | 0.17 s | 24 MB |
| `[30,24,18]²` (three-row) | 419 032 | 3.87M | 56 974 473 | **1.09x** | 1.36 s | 10.5 s | 485 MB |
| `[16,13,10,7]²` | 390 075 | 1.96M | 8 898 576 | **41x** | 0.33 s | 1.63 s | 242 MB |
| `[20,16,12,8]²` | 1 393 833 | 6.90M | 40 977 079 | **131x** | 1.77 s | 9.7 s | 1.01 GB |
| `[22,18,14,10]²` | 2 841 490 | 14.95M | 92 294 835 | **189x** | 4.26 s | 24.8 s | 1.36 GB |
| `[24,20,16,12]²` | 5 313 471 | 29.69M | 185 983 383 | **243x** | 10.5 s | 57.5 s | 1.96 GB |

**The 148 s is 8.5–10.5 s wall** (three runs today; 45–58 s CPU) — the
bitmap key did to this case what it measured on `[22,18,14,10]²`, and the
standing number predates it. So the ceiling on anything aimed at this case
is about ten seconds of wall.

**At four rows the layer is not enumeration.** Tableaux over productions is
1.07–1.09x at three rows — one production per tableau, which is why the
fibre count won there — and 41x, 131x, 189x, 243x on the four four-row
squares, growing with size. A four-row fibre count would replace
O(productions), 23–35 per term, at 4–11 µs of CPU per term, not O(tableaux).
Its state is also five components, not the four item 3 listed:
(λ¹ⱼ, λ²ⱼ, aⱼ, bⱼ, cⱼ). With three strips λ³ = λ is the candidate itself, so
`three_row.rs` applies λ³ⱼ₊₁ ≤ λ²ⱼ a row early and drops λ²ⱼ from the
state; with four, λ³ⱼ₊₁ ≤ λ²ⱼ has to be checked when λ³ⱼ₊₁ is chosen, so λ²ⱼ
is carried. On `[24,20,16,12]²` the per-row box is (λ¹−μⱼ ≤ 24) × (λ²−λ¹ ≤ 20)
× (a ≤ 24) × (b ≤ 20) × (c ≤ 16) ≈ 4.7M cells against ~13k for the three-row
`[24,20,16]²`, and today's three-row route runs at 1.2–2.7 µs per term
(`examples/calibrate_three_row`, below). To tie it would have to come in
under ~9 µs of CPU per term on this case with a state two dimensions
larger, and under ~1.6 µs of wall unless it is parallelized over candidates
as the layer is over rows. The caveat that the three-row box estimate was
wrong by an order of magnitude stands, and so does the conclusion of the
2026-07-31 section that no one-traversal method escapes enumeration at few
rows — but at four rows the layer does not need to escape it. **Dropped as a
build target.**

**The counting bands have moved, the second time a dispatch bound has gone
stale.** The same session ran both calibration harnesses of record,
in-process, on battery, order-alternating in the three-row one:

| `calibrate_three_row`, dispatched rows | 2026-07-31 | today |
|---|---|---|
| `[10,8,6]²`, `[12,10,8]²`, `[14,12,10]²`, `[16,14,12]²` | 1.31x, 1.43x end to end for the second and third; the in-process sweep read 1.03–2.06x on every dispatched row | 0.86, 0.97, 1.04, 1.09x |
| `[20,16,12]²`, `[22,18,14]²`, `[24,20,16]²` | 1.42x, 1.45x end to end for the first two | **0.60, 0.56, 0.69x** |
| `[18,14,10]·[9,7,5]`, `[20,16,12]·[10,8,6]` | in the 1.03–2.06x sweep | 0.86, 0.87x |
| `[16,13,10,7]·[8,6,4]`, `[14,12,10,8,6]·[7,5,3]` | 1.40x end to end, 1.36x in the crossover probe | 1.06, 1.22x |

(The 2026-07-31 end-to-end figures are out-of-process `lr_cli` A/Bs of the
layer against counting; the in-process sweep is this same harness on AC.)

`calibrate_two_row` (not order-alternating) reads its dispatched rows at
0.69–2.73x: `[28,22,17]·[28,22]` 0.71x, `[16,13,10,7]·[16,13]` 0.76x,
`[24,19,14]·[24,19]` 0.69x, `[20,16,12]·[20,16]` 1.21x, ties near n = 116–126,
and 2.0–2.7x at the top (`[50,40,30]·[50,40]`, `[34,28,22,16]·[34,28]`).

One change has touched the layer on these shapes since 2026-07-31, and none
has touched the counting routes: the bitmap key, 1.22–1.32x on exactly the
three-row squares (the product-walk and transposition rules leave a
three-row square on the direct walk it always took, and the allocator
study's one landing was on conjugated shapes only). That does not account
for the whole swing — `[20,16,12]²` read 1.42x end to end then and 0.60x
in-process now — and the rest is unattributed until the AC run: battery favors the parallel side of a
parallel-against-serial comparison (the 2.02x-vs-1.73x finding under "Power
state" in [README.md](README.md)), and the counting routes are
single-threaded — on CPU they are ahead of the layer by 1.6x on `[20,16,12]²`
(0.077 s against 0.12 s), 2.2x on `[22,18,14]²`, 3.2x on `[24,20,16]²` and
9x on `[30,24,18]²`, where counting is 2.7 µs a term against the layer's
25 µs of CPU and 3.2 µs of wall. What the tail item asks for is the
protocol the bounds were set by: on AC, out of process, `lr_cli` builds with
counting forced on and off, interleaved, min of 5 — before `prefer_counting`
moves in either direction. The likelier fix than narrowing the bands is to
parallelize the fibre count over candidates, which the CPU column says would
put it well ahead on wall; that is a build, and it waits on the same AC
number.

## 2026-08-18, the counting routes go parallel over candidates: 1.1–10x over the layer, and both bands widen

The AC number arrived the same day, and it went the way the CPU column said.
Protocol throughout: **AC (charging)**, out of process, one `lr_cli` binary
carrying two environment switches — counting forced off, and counting forced
on wherever a route applies — arms alternating per case, min of 5, output
digests compared (the script and binary were session scratch; the switches
were a temporary edit to both `prefer_counting`s, not shipped).

**Before the change, the stale bands confirmed on AC.** Layer against
counting as dispatched, every case single-threaded on the counting side:

| dispatched today | terms | layer | count | |
|---|---|---|---|---|
| `[10,8,6]²`, `[12,10,8]²`, `[14,12,10]²`, `[16,14,12]²` | 3k–20k | 3.3–15.9 ms | 3.5–15.6 ms | 0.95–1.02x |
| `[20,16,12]²` | 64 335 | 58.2 ms | 90.0 ms | **0.65x** |
| `[22,18,14]²` | 99 208 | 101.1 | 146.4 | **0.69x** |
| `[24,20,16]²` | 145 505 | 170.9 | 222.8 | **0.77x** |
| `[30,24,18]²` | 419 032 | 1 448 | 1 558 | 0.93x |
| `[18,14,10]·[9,7,5]`, `[16,13,10,7]·[8,6,4]`, `[14,12,10,8,6]·[7,5,3]`, `[20,16,12]·[10,8,6]` | 5k–12k | 5.2–11.9 | 5.8–10.2 | 0.87–1.17x |
| two-row: `[20,16,12]·[20,16]`, `[28,22,17]·[28,22]`, `[16,13,10,7]·[16,13]`, `[24,19,14]·[24,19]` | 8k–27k | 5.5–18.6 | 6.7–22.1 | **0.80–0.84x** |
| two-row: `[34,27,20]·[34,27]`, `[40,32,24]·[40,32]`, `[24,20,16,12]·[24,20]`, `[30,24,18]·[30,24]` | 36k–106k | 35–134 | 34–109 | 1.03–1.23x |
| two-row: `[50,40,30]·[50,40]`, `[34,28,22,16]·[34,28]` | 249k, 384k | 641, 1 205 | 306, 440 | 2.09x, 2.74x |

So the three-row band was a net loss out of process on AC — nothing above
1.17x, the three mid squares at 0.65–0.77x — and the two-row band lost its
lower third. The battery in-process numbers of the section above had the
direction right and the magnitude 1.2–1.5x too pessimistic.

**Built: `candidates.rs`.** The two routes' candidate walks were one function
written twice with a different depth (`strips` = 2 or 3: at most that many
new rows, λⱼ ≤ μ_{j−strips}); it is now one `pub(crate)` walk, and
`count_all` drives it: the candidates are split by their first two rows into
work items (hundreds to thousands for a dispatched product), workers claim
items from an atomic counter — lexicographic order, so the deepest subtrees
go first — and each worker owns a `Fibre` (the route's per-product scratch,
now per worker: two dense tables for three rows, two vectors for two) and an
output vector; the vectors are concatenated and sorted, so the output does not
depend on the thread count. A worker's panic is re-raised on the caller with
`resume_unwind`, so `two_row_product`'s documented `i128` panic reaches the
caller as itself rather than as a "worker panicked" message. One worker per
32 items and never more than `available_parallelism`; a product with a few
dozen items stays on the calling thread. The counting routes still do not
poll for interrupts (they never did; `interrupt.rs` says why a worker may
not, and the calling thread now spends its time in `join`).

**After, same protocol, counting forced on:** every dispatched case wins, and
the row clauses of both predicates turn out to be excluding the largest
wins:

| | terms | layer | count | |
|---|---|---|---|---|
| `[10,8,6]²` (n = 48) | 3 114 | 3.1 ms | 2.7 ms | 1.12x |
| `[12,10,8]²`, `[14,12,10]²`, `[16,14,12]²` | 7k–20k | 5.3–15.1 | 3.7–7.2 | 1.44–2.09x |
| `[20,16,12]²`, `[22,18,14]²`, `[24,20,16]²` | 64k–146k | 55–142 | 27–61 | 1.99–2.31x |
| `[30,24,18]²` | 419 032 | 1 206 | 327 | **3.69x** |
| three-row asymmetric, 3–5-row μ (four cases) | 5k–12k | 4.5–10.4 | 3.7–7.6 | 1.23–1.48x |
| **three-row, six-row μ**: `[12,11,10,9,8,7]·[6,5,4]`, `[14,12,10,8,6,4]·[9,7,5]` | 10k, 95k | 8.2, 115 | 6.6, 44 | 1.24x, 2.62x |
| **seven-row μ**: `[12,…,6]·[6,5,4]`, `[16,14,…,4]·[10,8,6]` | 22k, 520k | 20, 2 586 | 13, 260 | 1.58x, **9.96x** |
| **eight-, ten-, twelve-row μ**: `[10,…,3]·[8,6,4]`, `[10,…,1]·[8,6,4]`, `[12,…,1]·[9,7,5]` | 67k, 159k, 1.53M | 50, 162, 9 507 | 35, 88, 1 238 | 1.44x, 1.84x, **7.68x** |
| **ten-row μ, ratio 6.3**: `[14,…,5]·[6,5,4]` | 215k | 324 | 151 | 2.15x |
| two-row band as dispatched, ten cases | 8k–384k | 5.4–1 039 | 4.1–129 | 1.31–8.04x (`[50,40,30]·[50,40]` 6.44x, `[34,28,22,16]·[34,28]` 8.04x) |
| **two-row μ**: `[40,24]²`, `[70,42]²`, `[110,66]²`, `[60,30]·[50,40]` | 10k–200k | 5.9–273 | 4.0–71 | 1.47x, 2.38x, 3.87x, 2.08x |
| **seven- to sixteen-row μ, two-row ν**: `[12,…,6]·[12,9]`, `[12,…,6]·[16,12]`, `[16,14,…,2]·[16,12]`, `[14,…,5]·[14,11]`, `[10,…,1]·[12,10]`, `[12,…,1]·[12,10]`, `[15,…,1]·[10,8]`, `[16,…,1]·[16,12]`, `[20,18,…,6]·[20,16]` | 20k–20.7M | 11.6 ms – 66.6 s | 10.6 ms – 25.1 s | 1.10, 1.18, 3.53, 1.63, 1.18, 1.62, 2.83, 2.66, **7.42x** |
| controls that stay out: `[8,6,4]²` (n = 36), `[6,5,4,3,2,1]·[6,5,4]` (n = 36), `[10,8,6]·[10,8]` (n = 42), `[30,24,18]·[3,2,1]` (ratio 12) | | | | 1.02, 1.07, 0.97, 0.98x — all at the 2–3 ms floor |
| `[30,24,18]·[6,5]` (ratio 6.5, now inside the two-row bound) | 504 | 2.9 | 2.9 | 1.00x, at the floor |
| controls that stay out, and lose: `[160]·[80,50]` (one-row μ), `[30,24,18]·[40,2]` (lopsided ν), `[8,6,4]·[40,32]` (μ small against ν) | | | | 0.79x, 0.64x, 0.89x |

**The predicates that landed.** Two-row: `rows ≥ 2` (was `3..=6`) and
`8·|ν| ≥ |μ|` (was 3); the lopsided-ν, μ-small and n ≥ 75 clauses stand on
today's losses and floor ties. Three-row: `rows ≥ 3` with no upper bound (was
`3..=5`) and `8·|ν| ≥ |μ|` (was 4); balanced-ν, μ-small and n ≥ 48 stand.
Every case above is admitted or excluded as measured; the ratio-8 bound sits
between the ratio-6.7 win (`[15,…,1]·[10,8]`, 2.83x on 3.4M terms) and the
ratio-12 floor tie, with nothing measured between — a small product at ratio
6–8 may tie at the floor, which is the trade the calibration criterion
accepts. Both pinning tests carry the new cases with their ratios.

In-process, the two harnesses of record on AC after the change: three-row
2.4–4.0x on every dispatched row (`[8,6,4]²` at n = 36 reads 1.25x here and
1.02x out of process, so it stays out); two-row 2.4–12x.

Memory and CPU, `/usr/bin/time -l`, one cold process each: `[30,24,18]²`
counting 0.33 s wall, 2.47 s CPU, 149 MB against the layer's 1.31 s, 10.5 s,
482 MB; `[34,28,22,16]·[34,28]` counting 0.13 s, 0.69 s CPU, 114 MB against
1.08 s, 1.06 s, 103 MB — that layer walk runs nearly serial (its rows never
reach the parallel threshold) and the count's parallel overhead shows in its
CPU. Per-worker scratch is the three-row route's two tables, sized
`(μ₁+|ν|+1)(ν₁+1)(ν₂+1)` cells at 20 bytes: 3.2 MB a worker on `[30,24,18]²`,
so tens of MB across ten workers on the largest dispatched products. λ¹ never
exceeds μ₁+ν₁, so `max_l1` could shrink the table about |ν|/ν₁-fold;
unmeasured, and only worth it if a dispatched product ever runs the count into
memory, which none measured does.

## Next, in priority order

1. ~~**Parallelism.**~~ **Done** — the row-parallel fill with a sharded merge
   is the "Parallel LR" section above, at 2.86x on `[16,13,10,7]²`; this item
   predates it and was left standing by mistake.
2. ~~**Few-row factors below the counting crossover**~~ **Done 2026-07-31**,
   by exactly the route this item named: a cheaper fibre count (packed state,
   window-form inner loop) lowered the crossover to n ≥ 48, and the whole
   three-row band plus the ℓ(ν) = 3 asymmetric family now measures ahead of
   lrcalc — 1.06–1.33x where it was 0.71–0.88x. See "the wide-band deficit
   closed" above; the AC re-validation the first version of this item asked
   for ran the same day and confirmed both the bound and the sweep (1.02–1.49x
   against lrcalc, interleaved). The widening question it left — `rows ≤ 5`
   → 6, 1.24–1.32x in-process three times with no out-of-process number —
   closed 2026-08-18 with the parallel count: 1.24x and 2.62x out of process
   on AC, and the upper bound is gone altogether ("the counting routes go
   parallel" above).
3. ~~**Extend counting to four-row factors.**~~ **Dropped 2026-08-18**
   ("four-row factors" above). Both premises failed on measurement:
   `[24,20,16,12]²` is 8.5–10.5 s of wall today, not 148 s, and at four rows
   the layer compresses 41–243x tableaux per production, so it is not the
   enumeration a fibre count escapes — the count would replace 23–35
   productions per term at 4–11 µs of CPU with a five-component state
   (λ¹ⱼ, λ²ⱼ, aⱼ, bⱼ, cⱼ), not the four this item listed. What the item got
   right stands: the estimate is untrustworthy in both directions, and only a
   prototype would settle it; the bar it would have to clear is recorded
   there.
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
6. **`AutoLr`'s product-cache peek, versus a per-λ DP.** Replacing the peek
   would be a regression for the many-λ sweep pattern and a win for one-shot
   queries; which dominates has not been measured. Lifted here from
   [two_row.rs](../../src/two_row.rs), where it sat as rustdoc — the reference
   describes the present, so an unmeasured trade-off belongs in this tail.
   Since 2026-08-18 the peek finds products from every route ("the consumer
   side" above), so the sweep pattern it protects now holds off the counting
   routes too — a cold sweep of a two-row product's terms was 400x the
   product's own cost before that. The one-shot question is unchanged: a
   `two_row_coeff` after the peek would cost the sweep nothing.
7. ~~**Keys out of line, in a per-shard arena.**~~ **Superseded 2026-08-18**
   by the bitmap key (the section above): the density this item wanted — a
   16-byte entry of `(u32 offset, u32 len, C)` plus the bytes in an arena —
   is now a 24-byte entry with the whole 16-byte state inline and no
   indirection, and moving that state into an arena would add an offset
   without removing any bytes. What remains open on this axis is the
   accumulator: 8 of the 24 bytes are the `u64` multiplicity, and a `u32`
   first pass with the existing widen-and-rerun fallback would make an entry
   20 bytes — but the rerun costs a whole traversal, and it would fire on
   exactly the largest shapes, whose partial-filling multiplicities are
   tableau counts far past 2³². Unmeasured, and not obviously a win.
8. ~~**`prefer_conjugate` for skew expansions is calibrated against the byte
   key.**~~ **Done, same day, in two steps** — first conditioned on the
   outer shape descending by at most one cell per row ("skew-shape
   transposition follows the outer shape's steps"), then, once the mixed-step
   shapes were measured, on the inner shape reaching to within four rows of
   the bottom ("mixed-step outer shapes"): 1.2–6x on bands, steep deep
   shapes 1.1–1.75x faster than the original rule. Open within it: the
   gap-4 boundary is a coin flip decided by the calibration criterion, and
   bands under 60 cells (four of five measured wanted the transpose, by
   1.1–1.8x, at 5–16 ms) sit below the floor unmeasured further.
9. ~~**The direct regime's second-order choice.**~~ **Closed, same day,
   without a rule** — measured worth ~4% on average and 21% at most against
   "always the smaller factor", with the two cost components pulling
   opposite ways (fewer fillings against costlier ones); recorded above.
10. **`Schur::mul`'s one-pass build at the multi-million-term scale.** The
    first pair's terms go through a temporary vector, 32 bytes per term for
    an `i64` coefficient, alive alongside the map it builds ("the consumer
    side" above). At `[24,20,16,12]²` that is a few hundred MB next to the
    ~380 MB cached expansion and the ~500 MB map; sorted insertion has no
    such vector and measured 6x slower on this stage, which at that scale
    is about 1% of the expansion. A term-count threshold between the two
    would cost the small, repeated products nothing and give the huge ones
    the peak back — unmeasured there, since one run is 148 s and 2 GB.
11. ~~**A Pieri route in `AutoLr`.**~~ **Measured and declined 2026-08-18**
    ("a Pieri route" above): the layer is within 0.7–2.1x of a direct strip
    enumeration on every one-row or one-column product tried, its time
    proportional to the output, so a route would buy a constant factor on
    sub-10 ms products. Battery numbers; a re-measurement on AC would need a
    harness, since the probe was not kept.
12. ~~**Re-calibrate the two- and three-row counting bands, on AC.**~~
    **Done 2026-08-18, the same day** ("the counting routes go parallel"
    above): the AC out-of-process run confirmed the loss (three-row band
    0.65–1.17x, nothing above 1.17x; two-row lower third 0.80–0.84x), the
    fibre count was parallelized over candidates (`candidates.rs`), and the
    bands were re-fitted around that — wider, not narrower: two-row `rows ≥ 2`
    and `8·|ν| ≥ |μ|`, three-row `rows ≥ 3` with no upper bound and
    `8·|ν| ≥ |μ|`, every dispatched case 1.1–10x over the layer. Left open
    inside it: the ratio-6.7 to 12 gap on both size clauses (a floor tie at
    12, a 2.8x win at 6.7, nothing between), and two-row n between 42 (a
    floor tie) and 75 (the bound), unmeasured since 2026-07-27.
