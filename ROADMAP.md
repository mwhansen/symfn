# symfn roadmap — base functionality to par

**Target ("at par"):** the five classical bases {m, e, h, p, s}, multiplication in
each, all conversions between them, the ω involution, the Hall inner product, the
partition conjugate — **plus the Hopf structure**: skew Schur functions,
comultiplication, counit, and antipode. All tested against Sage as oracle.

## Current state

**Phases 0–4 complete; Phase 5 partially complete.** 45 tests green
(34 unit + 6 algebra-law + 4 oracle + 1 doctest), no external dependencies.

All five bases exist as real types with multiplication; every ordered pair of
bases converts; ω, the Hall inner product, and the full Hopf structure (skew
Schur, coproduct, counit, antipode) are implemented and law-tested.

Remaining: the three Phase 5 items that need network access to crates.io
(`rug`, `pyo3`, `proptest`) — see below.

## Design decisions (settled)

1. **Introduce ℚ** via a `Ring` / `Field` trait split. Integral bases (s, h, e, m)
   stay exact over ℤ; only the power-sum basis, ω on p, and the inner product pull
   in ℚ. (Later: `rug::Rational` under the `gmp` feature.)
2. **Hub-and-spoke conversions through Schur** for integral bases (compose to
   cover all 20 ordered pairs) rather than a full pairwise matrix.
3. **Tensor type** `SymTensor` (map `(Partition, Partition) → C` in the Schur
   basis) as the codomain of the coproduct.

---

## Phase 0 — coefficient & partition groundwork
*Small; unblocks everything.*
- [x] Split `Coeff`: `Ring` (add/mul, today's trait) and `Field: Ring` (+ `inv`).
- [x] Pure-Rust `Rational` implementing `Field` (i128 numer/denom, normalized).
- [x] `Partition::conjugate()` (transpose) and `Partition::z()` = ∏ iᵐⁱ·mᵢ!.
- **Done when:** rationals arithmetic tested; `conjugate` involutive; `z` matches
  known values (z of (2,1,1)=4, etc.).

## Phase 1 — the multiplicative bases
*Easy; mostly free.*
- [x] Add `Elementary<C>` (e) and `Homogeneous<C>` (h) basis types.
- [x] Multiplication for e, h = multiset union of parts (same rule as p).
- [x] m-multiplication (routed through Schur, in Phase 2).
- **Done when:** eλ·eμ, hλ·hμ tested; all five basis types exist.

## Phase 2 — the transition engine  *(the core)*
Per-basis-element conversion to/from the Schur hub, extended linearly; compose
for all pairs. A uniform `convert::<Target>()` on top.
- [x] `e ↔ h`: recursion Σ(−1)ⁱ eᵢ h₍ₙ₋ᵢ₎ = 0 (integral).
- [x] `s ↔ h`: Jacobi–Trudi `sλ = det(h₍λᵢ−i+j₎)`; inverse via Pieri (integral).
- [x] `s ↔ e`: dual Jacobi–Trudi (integral).
- [x] `s ↔ m`: Kostka `sλ = Σ Kλμ mμ`; inverse Kostka (integral).
- [x] `s ↔ p`: Murnaghan–Nakayama characters; `sλ = Σ zμ⁻¹ χλ(μ) pμ` (needs ℚ).
- [x] `m` multiplication via convert→Schur→multiply→convert-back.
- **Done when:** every ordered pair converts; round-trips are identity on a sweep
  of small partitions; transition values match Sage.

## Phase 3 — standard operations
*Small once Phase 2 lands.*
- [x] ω involution per basis: `ω(pλ)=(−1)^{|λ|−ℓ(λ)}pλ`, `ω(hλ)=eλ`, `ω(sλ)=s_{λ'}`.
- [x] Hall inner product `⟨·,·⟩` via p-expansion + `z` (check vs Schur orthonormality).
- **Done when:** ω²=id; ⟨sλ,sμ⟩=δλμ; ⟨pλ,pμ⟩=zλ δ across small degrees.

## Phase 4 — Hopf structure
*Leverages the existing LR machinery.*
- [x] `SymTensor<C>` type (element of Sym ⊗ Sym in the Schur basis) with linear ops.
- [x] Skew Schur `s_{λ/μ} = Σν c^λ_{μν} sν` (direct from `LrBackend`).
- [x] Coproduct Δ: on multiplicative bases `Δ(hₙ)=Σ hᵢ⊗h₍ₙ₋ᵢ₎`, `Δ(eₙ)` dual,
      `Δ(pₙ)=pₙ⊗1+1⊗pₙ`; on Schur `Δ(sλ)=Σ c^λ_{μν} sμ⊗sν`.
- [x] Counit ε (project to degree-0 coefficient) and antipode `S(sλ)=(−1)^{|λ|}s_{λ'}`.
- **Done when:** coassociativity on small inputs; `Δ` agrees across bases;
      antipode satisfies the Hopf axiom `m∘(S⊗id)∘Δ = ε·1` on small degrees.

## Phase 5 — validation & interop
*Makes it trustworthy and actually usable from Sage.*
- [x] `Ring::from_u128` injection, so structure constants (LR, Kostka, z_λ) are
      never silently truncated by a cast — the seam a bignum type needs.
- [x] Cross-cutting algebraic-law suite (`tests/algebra_laws.rs`): conversions are
      ring homomorphisms, round-trips are identity, ω is an involutive algebra
      map, Hall pairings ⟨s,s⟩/⟨h,m⟩/⟨p,p⟩ correct, Δ is an algebra map.
- [x] Sage oracle: `scripts/gen_sage_oracle.sage` generates a 743-line fixture
      (Kostka, characters, Schur products, all four conversions, skew Schur);
      `tests/sage_oracle.rs` checks symfn against it. Fixture committed, so
      `cargo test` needs no Sage.
- [x] `gmp` feature: `impl Ring for rug::Integer`, `Ring`+`Field` for
      `rug::Rational`. `from_u128` exact; verified on coefficients beyond i64.
- [x] `python` feature: PyO3 module (abi3, so one artifact serves CPython 3.9+),
      coarse-grained whole-object API. Verified importing into Sage (Py 3.14).
- **Done:** `import symfn` works inside Sage; 49 Schur products cross-checked
      in-process with zero mismatches; benchmarks show 2–20x over the pure-Python
      path (see below).

## Phase 6 — performance (memoization → plethysm → optimized LR)

- [x] **Memoization** (`src/memo.rs`). Thread-safe caches of *pure* function
      results: partition lists, characters, Kostka numbers, LR coefficients,
      inverse-Kostka data, and whole Schur products. Referentially transparent,
      so unlike Symmetrica's mutable globals nothing is observable in results.
      `clear_caches()` reclaims the memory.
- [x] **Plethysm** (`src/plethysm.rs`). Computed through the power-sum basis,
      where `p_n[g]` is just part-scaling and plethysm is multiplicative in f.
      Cross-checked against Sage on 27 values.
- [x] **Optimized LR backend** (`src/strip_lr.rs`). `StripLr`: a row-level DP
      over horizontal strips with state merging, replacing cell-by-cell
      enumeration. Verified against `NaiveLr` on every product with |μ|+|ν| ≤ 7.
- [x] **Whole-shape LR backend** (`src/skew_lr.rs`). `SkewLr`: expand a skew
      shape in one traversal, advancing a *merged frontier* of partial fillings
      rather than enumerating tableaux individually. Products go through the
      disconnected skew shape whose skew Schur function is s_μ·s_ν. Now the
      default; verified against both other backends exhaustively and against
      `lrcalc` on large shapes. Written clean-room — see `NOTICE.md` and
      `docs/cleanroom-spec-skew-lr.md`.

### Benchmark summary (vs Sage, same machine)

| workload | speedup |
|---|---|
| 144 small Schur products | **15.2x** |
| repeated large products | **8–10x** |
| single large products, cold | **1.2–4x** |
| 225 Kostka numbers | **62x** |
| plethysm (12 cases, cold) | **13.9x** |

### Known ceiling: the Python boundary is `i128`, Python integers are not

`src/python.rs` crosses coefficients as `i128`. That is wide enough for every
structure constant this library computes in practice, and it replaced an `i64`
boundary that silently truncated plethysm numerators — but it is still a
fixed width, while Python integers are arbitrary precision and the `gmp`
feature is exact.

So `--features gmp` and `--features python` do not compose: a Sage caller
cannot get bignum coefficients even though the core supports them. Closing
that needs a decision about how arbitrary-precision integers cross PyO3
(`num-bigint` via pyo3's feature, a decimal-string representation, or making
the module generic and exposing two variants), which is an API question rather
than an implementation one. Until then, the ceiling on that boundary is
symfn's, not Python's.

### Two findings worth remembering

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
ballot constraints reduce to O(1) per run — and by packing the frontier key into
a single buffer, so a transition that *merges* (the common case, and the whole
point of a frontier) allocates nothing. `examples/bench_shapes.rs` is the
interleaved A/B harness for measuring that kind of change.

Measured (`cargo run --release --example bench_lr`):

| product | NaiveLr | StripLr | SkewLr | vs best |
|---|---|---|---|---|
| s[5,4,3,2,1]² | 0.0113s | 0.0141s | 0.0014s | **8.1x** |
| s[6,5,4,3,2]² | 0.1416s | 0.1220s | 0.0071s | **17.2x** |
| s[6,5,4,3,2,1]² | 0.9398s | 0.2991s | 0.0160s | **18.7x** |
| s[7,6,5,4,3]² | 1.3488s | 0.9215s | 0.0321s | **28.7x** |
| s[8,7,6,5,4,3]² | 327.69s | 30.14s | 0.83s | **36.2x** |

The win is larger still on skew expansions, where the old path ran a full
backtrack per candidate content:

| skew shape | per-ν | SkewLr | speedup |
|---|---|---|---|
| s[8,7,6,5,4]/[3,2,1] | 0.0036s | <0.0001s | **124x** |
| s[9,9,8,8,7]/[2,1] | 0.0545s | <0.0001s | **3040x** |
| s[10,9,8,7,6,5]/[4,3,2,1] | 0.0442s | 0.0002s | **257x** |

`SkewLr` wins at every size measured — checked explicitly at small sizes, where
the frontier map's hashing could have dominated; it does not. So `AutoLr` no
longer dispatches. It stays a distinct type as the one place to reintroduce
dispatch if a future backend wins only in some regime.

**Where that leaves us against lrcalc: ahead everywhere except moderate wide
shapes.** Run `scripts/compare_lrcalc.py`, which drives both as CLI processes on
identical inputs, verifies the outputs agree, and times best-of-N. Only cases
well clear of the ~6ms process-startup floor say anything about either
algorithm:

| product | lrcalc | symfn | |
|---|---|---|---|
| s[8,7,6,5,4,3]² | 8.64s | 0.82s | **us 10.5x** |
| wide [16,13,10,7]² | 14.55s | 2.81s | **us 5.2x** |
| s[9,8,7,6,5]² | 0.848s | 0.360s | us 2.4x |
| s[6,5,4,3,2,1]² | 0.068s | 0.030s | us 2.3x |
| s[8,7,6,5,4]² | 0.284s | 0.136s | us 2.1x |
| s[7,6,5,4,3]² | 0.075s | 0.048s | us 1.6x |
| wide [12,10,8]² | 0.014s | 0.017s | lrcalc 1.2x |
| wide [14,12,10]² | 0.021s | 0.028s | **lrcalc 1.3x** |
| wide [20,16,12]² | 0.183s | 0.293s | **lrcalc 1.6x** |
| wide [24,20,16,12]² | >200s | **118s** | we finish, lrcalc doesn't (see below) |

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
32% in the outer merge. Instrumenting the frontier explains the difference, and
it is not the one previously recorded here:

| factor rows | shapes | LR tableaux ÷ states produced |
|---|---|---|
| 2 | `[24,12]²`, `[20,10]²` | **1.0x** |
| 3 | `[12,10,8]²`, `[20,16,12]²` | **1.0–1.1x** |
| 4 | `[16,13,10,7]²` | 40.8x |
| 5 | `[9,8,7,6,5]²`, `[12,10,8,6,4]²` | 11–164x |
| 6–7 | `[8,7,6,5,4,3]²`, `[7,6,5,4,3,2,1]²` | 57–80x |

The frontier's entire advantage is that one state stands for many tableaux. At
two or three rows it stands for **one** — `[20,16,12]²` produces 2 614 952
states against 2 802 764 LR tableaux. So the DP walks exactly what a
tableau-at-a-time enumerator walks and then pays key assembly, encoding and
hashing on top of it. That is the whole deficit, and no amount of tuning the
frontier removes it; the frontier is the cost.

⚠️ This corrects an earlier claim here that few rows mean *little merging*.
Merging within the frontier is in fact strongest exactly where we lose:
`[20,16,12]²` collapses 2.6M productions into 104 876 live states (24.9x), the
highest measured, while `[16,13,10,7]²` manages only 2.2x. Those are different
quantities — productions-per-live-state says the frontier stays small,
tableaux-per-production says whether the enumeration was avoided — and only the
second is a saving. The earlier note conflated them.

**A frontier-free enumerator was prototyped and does not help — measured, and
worth not repeating.** The obvious consequence of 1.0x compression is that
few-row shapes should skip the frontier: walk the same run decompositions
depth-first, keep the previous row on the stack instead of in a hashed key, and
hash only completed tableaux keyed on content. That was built and verified
against `AutoLr` on every shape below. It is **parity at best**:

| shape | frontier | direct | |
|---|---|---|---|
| `[12,10,8]²` | 8.35ms | 7.99ms | 1.05x |
| `[14,12,10]²` | 13.7ms | 14.9ms | 0.92x |
| `[20,16,12]²` | 266ms | 283ms | 0.94x |
| `[9,8,7,6,5]²` | 293ms | 1.95s | 0.15x |
| `[8,7,6,5,4,3]²` | 657ms | 22.9s | 0.03x |

(A first version measured 0.75x on `[20,16,12]²` purely because it used std's
`HashMap`; SipHash against the frontier's own fast hasher measures the hasher,
not the algorithm. The table is after fixing that.)

The row-fill counts confirm the compression reading rather than contradicting
it — direct does 2 846 571 fills against the frontier's 2 614 952 on
`[20,16,12]²` (+9%), and 243 534 430 against 3 678 951 on `[8,7,6,5,4,3]²`
(66x worse, exactly the compression the frontier is buying there).

The reasoning error was treating the profile's 37% "state commit" as removable.
Removing the frontier **relocates** that cost rather than deleting it: every
completed tableau still has to be hashed to bin it by λ, and there are 2.8M of
them either way. The frontier hashes 2.6M longer keys spread across rows;
direct hashes 2.8M shorter keys at the last row.

So the real statement is: **both approaches touch all 2.8M LR tableaux while the
answer has only 64 335 terms.** Beating lrcalc here needs an algorithm that does
not enumerate tableaux at all — not a cheaper enumeration. Two candidate
directions, neither validated:

* Because ν has 3 rows, entries come from `{1,2,3}` and every *column* is one of
  7 subsets, so a DP keyed on (content, previous column) would have a tiny state
  space. The obstacle is mathematical, not incidental: the ballot condition is
  defined on the row reading word, and whether it survives a column-wise
  reformulation is an open question that should be settled before any code.
* The last row's fillings are already in bijection with their content
  increments (a weakly increasing row is determined by its per-value counts), so
  the final row could in principle be replaced by a convolution over frontier
  states rather than an enumeration.

The conjugate dispatch does not help either (it deliberately does not fire
here), because conjugating trades few rows for few columns and the state still
pins the tableau.

**Both implementations are single-threaded, and that is what makes this table
mean something.** lrcalc runs at ~99% of one core; symfn uses no threads at all.
So these ratios compare *algorithms*, not core counts. If symfn is ever
parallelised, this comparison must keep reporting a single-threaded number —
a multi-threaded wall-clock figure set against a single-threaded lrcalc would
conflate an algorithmic win with a hardware one, and note that a
parallel build pinned to one thread is not the same as a sequential build
(per-thread structures and merge machinery cost something even at N=1). lrcalc
being single-threaded is a property of its implementation, not of the problem;
its enumeration is at least as parallelisable as our frontier, so threads are a
real engineering win for users but not a durable claim of algorithmic
superiority.

**The memory wall, and how it fell.** `[24,20,16,12]²` was once killed at
27m27s wall with 2.0+ GB resident and still climbing, **38% of it system
time** — allocation and page faults, not combinatorics. Peak memory is
(states) × (bytes per state), and the second factor was soft. Three changes
(same commit series, measured by interleaved A/B with
`examples/bench_shapes.rs`, which now also reports peak live frontier states):

1. **Byte-packed inline keys.** Every element of a frontier key is bounded by
   the shape's cell count, so keys serialize at one byte per element for
   anything practically computable and live inline in a 32-byte enum — no
   heap allocation per state at all. (Narrowing is where overflow bugs live:
   the width bound is proven in `elem_width` and tested across the 255/256
   and 65535/65536 boundaries against Pieri.)
2. **u64 map values with a checked u128 fallback.** Multiplicities are
   tableau counts with no provable narrow bound, so every merge is a
   `checked_add` and the expansion transparently reruns wider on saturation
   (exercised in tests with a u8 accumulator).
3. **Single-table frontiers.** Between rows the frontier is drained into an
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

(The table's [16,13,10,7]² time is the packed frontier alone, same
orientation; the dispatch below then takes it to 2.4s.)

`[24,20,16,12]²` = 5 313 471 terms, peak 23.0M live states. Independently
re-measured through the ordinary `lr_cli mult` path: **125.6s wall, 2.36 GB
peak RSS, and 5% system time** — down from 38% when the case was
allocation-bound, which is the diagnosis confirming itself rather than just
the symptom improving. Neither lrcalc (>200s timeout, still running at 27m in
earlier sweeps) nor the old representation finishes it on this machine. Much
of the remaining 2 GB is the 5.3M-term *output* (two copies: the memoized
`Arc` plus the caller's clone), not the frontier.

⚠️ **This result has no external oracle.** lrcalc cannot finish the case, so
those 5.3M terms are checked only against ourselves — the conjugate and direct
orientations independently produce identical output, which is a real
consistency check but not an independent one. Every *other* case in the
comparison table is verified term-by-term against lrcalc. Treat the largest
result as unconfirmed until something else can compute it.

**Transposition, re-measured on the right axes — and now dispatched.**
Since c^λ_{μν} = c^{λ'}_{μ'ν'} the walk can run on the transposed diagram.
An earlier note said wide shapes prefer the original orientation — true at
moderate size, but it inverts exactly where it matters. Peak *states* are
nearly orientation-independent (±25% both ways on every case measured; 23.0M
direct vs 21.3M conjugate on the big one — the frontier is the same
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

**Next**, in priority order:

1. **Parallelism.** Deliberately deferred until after the memory work
   (per-thread frontiers multiply residency); now that bytes-per-state is
   ~4× smaller, a row-parallel merge is the next big lever.
2. **Few-row factors** — the last regime where lrcalc beats us, 1.2–1.6x
   ([20,16,12]², [14,12,10]², [12,10,8]², and 2-row shapes). Measured cause:
   at ≤3 rows the frontier compresses 1.0x, so the DP does a naive
   enumerator's work plus hashing. **Skipping the frontier does not fix it** —
   prototyped, parity at best, see above; the hashing relocates rather than
   disappears. Any real win has to stop enumerating all 2.8M tableaux to
   produce 64K terms, so this is now a research question (column-wise DP, or
   a convolution for the final row) rather than an engineering one. Lowest
   confidence of anything on this list; do not schedule it as a known fix.
3. **Shape preprocessing** — factoring a skew diagram into connected
   components and expanding each separately, since the expansion of a
   disconnected shape is the product of its pieces.
4. **Output residency.** On `[24,20,16,12]²` a growing share of peak RSS is
   the 5.3M-term *output* (the memoized `Arc<Vec>` plus the caller's clone),
   not the frontier. An `Arc`-returning variant of `expand_skew` would halve
   that.

---

### Dependency order
`Phase 0` → `Phase 1` → `Phase 2` (needs 0,1) → `Phase 3` (needs 2) →
`Phase 4` (needs LR + conjugate + ω) → `Phase 5` (needs it all).

### Notes
- Keep `NaiveLr` as the in-house LR oracle: it is the most obviously-correct
  backend, and the faster ones are held to exhaustive agreement with it.
  (Skew Schur and the Schur coproduct now go through `SkewLr` instead.)
- Every phase ships with tests; Sage supplies ground truth wherever a value is
  non-trivial (Kostka tables, characters, transition matrices).
