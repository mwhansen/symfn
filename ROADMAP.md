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

So the 5–16x time wins come with 1.4–4.9x the memory. The frontier trades
residency for speed by construction, and on `[24,20,16,12]²` that reaches
2.36 GB — which is why lrcalc's inability to finish that case is not purely a
speed result.

⚠️ **Every timing in this file was taken on a laptop running on battery, and
that is an uncontrolled confound.** The same binary on `[16,13,10,7]²` measured
1.37s in a short run and 2.91s inside a 25-minute sweep — 2.1x. I first
attributed that to being memory-bound. That was wrong: re-measured on AC power
the case runs 1.20–1.58s, matching the *short* battery run, and the long sweep
had by then drained the battery to 14%, where macOS throttles aggressively. The
2.91s is a power/thermal artifact.

Two consequences, both practical:

* **Absolute times here are not comparable across runs**, and a long sweep
  throttles progressively, so its later rows are systematically pessimistic.
  Per-case *ratios* are the durable quantity, since the harness runs both sides
  adjacently under the same conditions.
* **Re-baseline before believing any change under ~2x.** Interleave builds,
  keep the machine on AC, and never compare a number from a long sweep against
  one from a short run.

The memory numbers above stand — they are not timing-sensitive — but the
inference that memory explains the variance does not.

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

**Counting per output instead of enumerating chains — validated on a model
problem.** `examples/model_count_vs_chains.rs` computes `s_μ·h_a·h_b` both ways:
by building every chain μ ⊂ λ¹ ⊂ λ² (what the frontier does, minus the lattice
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
`examples/lr2_count_vs_frontier.rs` runs the smallest LR case that carries the
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
2.99 µs); the frontier costs O(tableaux) and its per-term time grows far faster
(1.0 → 0.87 → 1.25 → 2.27 → **11.0** µs). Crossover is near
`[40,32,24]·[40,32]` and the gap widens after it.

⚠️ **Measure against `AutoLr`, not against chain enumeration.** The same file
also implements the naive chain enumerator, which counting beats by 5–244x —
a meaningless number, since the frontier DP exists precisely to beat chain
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

| product | terms | frontier | counting | |
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

**This result now has an independent oracle** (`examples/verify_specialization.rs`).
lrcalc cannot finish the case, and for a long time its 5.3M terms were checked
only against our own conjugate orientation — a real consistency check, but not
an independent one, since a bug in the shared frontier code reproduces itself in
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

### Non-LR baselines (`examples/bench_ops.rs`)

Every non-LR operation was unmeasured; `scripts/bench_vs_sage.py` covers only
Schur products. The first run of the new harness found the hot spot immediately:

| operation | time | work |
|---|---|---|
| `character_table_n16` | 0.036s | 53 361 values |
| `coproduct_[7,6,5,4,3]` | 0.059s | 8 336 terms |
| `convert_m_to_s` (deg 12) | 0.254s | inverts the 77×77 Kostka matrix |
| `kostka_all_pairs_n12` | 0.264s | 5 929 values |
| `kostka_warm_50x` | 0.00027s | 3 850 cached lookups (66 ns each) |

**Kostka is the bottleneck, and it is exponential.** Cost of one K_{λμ}:

| λ | degree | per value |
|---|---|---|
| `[3,2,1]` | 6 | 1.1 µs |
| `[4,3,2,1]` | 10 | 12 µs |
| `[5,4,3,2]` | 14 | 527 µs |
| `[5,4,3,2,1]` | 15 | 1.9 ms |
| `[6,5,4,3,2]` | 20 | **638 ms** |

At degree 20 that makes `convert_s_to_m` take **400 seconds** for a single
Schur function, since it needs K_{λμ} for all p(20) = 627 partitions μ. Compare
characters at 0.6 µs per value — Kostka is ~500x slower per number, and it is
the s ↔ m transition, so it also gates `convert_m_to_s` and anything routed
through the monomial basis.

**Fixed.** `kostka_uncached` enumerated SSYT one cell at a time. K_{λμ} instead
counts chains ∅ ⊆ λ¹ ⊆ … ⊆ λ with λⁱ/λⁱ⁻¹ a horizontal strip of size μᵢ, and
chains through the same intermediate shape *merge* — `strip_lr`'s DP without the
lattice condition. Pruning every intermediate to λ bounds the state space by the
partitions inside λ instead of by the tableaux of shape λ.

Measured by interleaved A/B of two `bench_ops` binaries, min of 2 rounds:

| operation | before | after | |
|---|---|---|---|
| `kostka_row_[5,4,3,2,1]` | 0.3246s | 0.0006s | **538x** |
| `kostka_row_[5,4,3,2]` | 0.0704s | 0.0003s | **216x** |
| `convert_s_to_m` | 0.0058s | 0.00013s | **46x** |
| `convert_m_to_s` | 0.2506s | 0.0075s | **34x** |
| `kostka_all_pairs_n12` | 0.2504s | 0.0080s | **31x** |

The speedup grows with degree, as exponential → polynomial should. On the case
that prompted this, `convert_s_to_m` for `s[6,5,4,3,2]` went from **400 s to
11.7 ms** (~34 000x), and degrees that were simply unreachable now run: degree
30 in 154 ms (2 317 terms), degree 40 in 2.70 s (16 306 terms).

Nothing else moved — characters, plethysm and the Hopf operations are unchanged
to within noise, which is the expected result since none routes through Kostka,
and is worth stating because a rewrite that quietly perturbed them would be a
regression hiding behind a headline.

Correctness: the Sage fixture covers only small degrees, so it cannot exercise
the new implementation where it now operates. `K_{λ,1ⁿ}` counts standard
tableaux, so the hook-length formula `n!/∏h(i,j)` checks it independently at
degrees up to 25.

**The coproduct had the same defect, one function away from its own fix.**
`skew_schur`'s doc comment records that it used to sweep all p(n) partitions
running a full LR backtrack per candidate, "the same answer for orders of
magnitude more work". `coproduct` was still doing exactly that: sweeping every
(μ, ν) pair of the right total degree and calling `NaiveLr::lr_coeff` on each.
Computing Δ(s_λ) = Σ_{μ⊆λ} s_μ ⊗ s_{λ/μ} instead gets every ν from one
traversal per μ — **0.057s → 0.0039s, ~15x**, identical 8 336 terms.

Worth noting as a pattern rather than a one-off: the fix already existed in the
same file, and the benchmark is what made the second instance visible. Nothing
about reading the code had surfaced it, in a codebase this well-commented,
because the comment explaining the mistake sat on the function that no longer
made it.

**After the Kostka rewrite the profile changed shape.** Rescaled to degree 20
(`examples/bench_ops.rs`), the leaders are now:

| operation | time | work |
|---|---|---|
| `convert_m_to_s` | **2.03s** | **1 term** |
| `kostka_all_pairs_n20` | 1.97s | 393 129 values |
| `convert_p_to_s` | 0.195s | 1 term |
| `hall_s_p_degree17` | 0.173s | 1 pairing |
| `coproduct_[8,7,6,5,4]` | 0.0092s | 19 758 terms |

`convert_m_to_s` builds the entire p(20)×p(20) Kostka matrix and inverts it to
use **one row**, and the matrix construction is ~97% of that (2.03s against
1.97s for the same number of Kostka values). The inversion itself is nearly
free by comparison.

A dominance early-out (`K_{λμ} ≠ 0` iff λ ⊵ μ) was added and is worth only
**1.1x** — I expected much more, on the reasoning that dominance is far sparser
than the lex triangle the matrix is built in. The chain DP's capacity pruning
was already rejecting those pairs quickly, so the test mostly replaces a fast
zero with a faster one.

**Fixed by computing only the row that is needed**, via the forward recurrence
`w[j] = 1`, `w[jj] = −Σ_{m=j}^{jj−1} w[m]·K[m][jj]`, cached per μ rather than per
degree. Measured **2.02s → 0.527s, 3.8x** — better than the ~2x predicted from
the triangle-versus-square argument alone, because skipping `m` where
`w[m] = 0` avoids the Kostka call entirely rather than computing a value and
multiplying it by zero. Inverse-Kostka rows are sparse enough for that to be
the larger half of the win.

**p → s without computing characters.** `p_μ = Σ_λ χ^λ(μ) s_λ` was implemented
by asking for χ^λ(μ) once per λ — p(n) independent recursions per term. But
Murnaghan–Nakayama is itself a multiplication rule, `p_k·s_λ = Σ (−1)^{ht} s_{λ∪ξ}`
over k-rim-hooks ξ, so multiplying successively by each part of μ builds the whole
expansion in ℓ(μ) passes and never evaluates a character.

| operation | characters | iterated MN | |
|---|---|---|---|
| `convert_p_to_s` | 0.203s | 0.0997s | **2.0x** |
| `hall_s_p_degree17` | 0.179s | 0.0991s | **1.8x** |
| `plethysm_[3,2][[2,1]]` | 0.00464s | 0.00188s | **2.5x** |
| `plethysm_[4][[3]]` | 0.00266s | 0.00109s | **2.4x** |

⚠️ **The first implementation of this measured 0.8x — a regression — and the
algorithm was not what changed.** Keying the frontier on `Vec<i64>` β-numbers
meant hashing a 160-byte key per rim hook, and the character path it was
competing against has a memo cache that shares subproblems across every term.
β values here are below 2n, so for n ≤ 32 the whole set is a `u64` bitmask:
adding a rim hook becomes two shifts and a popcount, and the key is one word.
Same algorithm, 2.5x swing. Degrees past 32 fall back to characters, which stay
exact in `C` for bignum rings.

That is the second time this session a data structure inverted an algorithmic
conclusion — the other was `HashMap` versus generation-stamped dense tables in
`three_row`, also worth 3.3x. Both times the operation counts were unchanged and
only the constants moved. **A structural idea that benchmarks badly on its first
implementation has not been tested yet.**

### Against Sage (`scripts/compare_sage.py`)

Sage is now a second comparison harness alongside lrcalc, covering the
operations lrcalc cannot. It differs in one respect and that respect is the
whole design: lrcalc runs **out of process**, one invocation per case, so both
sides pay startup and neither reuses a warm cache. Sage's interpreter costs
seconds, so both sides run **in one process** and startup is excluded by
construction — which forfeits the cache isolation a fresh process gave free, and
that has to be bought back explicitly:

* Sage memoizes, so every case uses *distinct inputs of comparable size*, each
  computed once. Repeating one input measures Sage's cache from the second call.
* symfn memoizes too, so `clear_caches()` (newly exposed through PyO3) runs
  before each timed call.
* Sage's first touch of a basis is lazy, so an untimed warm-up runs first or the
  first case absorbs setup belonging to all of them.

Every case is verified, not merely timed. **All cases agree with Sage at every
degree measured**, which independently cross-checks the rewrites above.

#### ⚠️ Half these rows are not benchmarks against Python

Sage's conversions between the five classical bases **are Symmetrica**:
`sage.combinat.sf.classical.init` populates `conversion_functions` with
`t_<FROM>_<TO>_symmetrica`, confirmed at runtime. Everything else in the table
is Sage's own Python. The harness now prints a `via` column (`C` / `py`) so a
ratio is never read without its baseline, because the two mean very different
things: **3.7x on a `C` row is a stronger result than 9x on a `py` row.**

| case | via | deg 8 | deg 12 | deg 16 | deg 20 |
|---|---|---|---|---|---|
| p → s | C | 6.3x | 6.8x | 7.5x | 7.6x |
| s → m | C | 3.2x | 3.7x | 4.3x | 5.6x |
| **m → s** | C | 5.7x | 6.4x | 5.5x | **3.7x** *(was 0.004x)* |
| **s → e** | C | 6.6x | 8.1x | 4.9x | **3.7x** *(was 0.4x)* |
| s → h | C | 4.8x | 4.6x | 4.2x | 3.8x |
| **s → p** | C | 3.1x | 2.0x | 1.6x | **1.4–5.1x** |
| Kostka | py | 1.3x | 8.6x | 7.2x | 6.0x |
| plethysm | py | 9.0x | 8.4x | 8.5x | 9.0x |
| coproduct | py | 2.6x | 2.7x | 3.1x | 2.4x |
| skew | py | 7.1x | 6.4x | 4.5x | 4.5x |
| Hall | py | 2.5x | 1.8x | 1.9x | 1.7x |

Both former scaling failures are fixed (see the commit "Fix the two scaling
failures the Sage ladder found"); the ladder is what exposed them, since each
was *faster than Sage at degree 8* and only lost as degree climbed.

#### The plethysm row is measuring the wrong baseline

Symmetrica also ships C implementations of **plethysm** and Schur products that
Sage does not use, so the `py` rows above compare against the weaker of two
available baselines. Called directly through
`sage.libs.symmetrica.all.plethysm` (min of 3 fresh processes, all outputs
agreeing):

| case | Symmetrica C | symfn | ratio |
|---|---|---|---|
| s_3[s_{21}] | 0.000105s | 0.000133s | 0.79x |
| s_4[s_{21}] | 0.000297s | 0.000462s | 0.64x |
| s_5[s_{21}] | 0.000829s | 0.001985s | 0.42x |
| s_3[s_{31}] | 0.000200s | 0.000657s | 0.30x |
| s_4[s_{22}] | 0.000571s | 0.004050s | **0.14x** |

So plethysm is **slower than Symmetrica, and the gap widens with output size** —
another scaling problem, hidden because the visible baseline was an interpreter.
Caveat: Symmetrica only supports a *single-row outer* here ("for the moment only
for outer S_n"), which is the easy case and may use a specialised path, so this
is not like-for-like on generality — but on the inputs where both run, we lose.

#### s → p never got the fix that p → s did

`PowerSum::from_schur` calls `character_in(λ, μ)` once per μ — **p(n)
independent Murnaghan–Nakayama recursions per λ**. That is precisely the pattern
removed from `p → s`, which replaced p(n) character queries with one iterated-MN
sweep and now runs at 7.6x. `s → p` is the mirror image and is our weakest
conversion row (1.4–1.6x at degree 16–20).

The fix is less mechanical than `p → s`'s was: `p_expand` produces a *column* of
the character table (one μ, all λ) in one sweep, whereas `s → p` needs a *row*
(one λ, all μ). Options: build the table column-by-column and cache it per
degree — a loss for a single cold conversion but a large win for the memoized
many-query pattern that is the real Sage usage — or find a direct row sweep.

⚠️ The plethysm row is capped at degree 10 regardless of the ladder setting
(its input is the *outer* partition and the result reaches degree 30), so that
row does not vary across the columns above — its four entries are the same
measurement repeated.

**Next**, in priority order:

1. **Parallelism.** Deliberately deferred until after the memory work
   (per-thread frontiers multiply residency); now that bytes-per-state is
   ~4× smaller, a row-parallel merge is the next big lever.
2. **Few-row factors below the counting crossover** — the remaining regime
   where lrcalc beats us, now 0.74–0.85x rather than 0.59–0.79x. Cause is
   understood: at ℓ(ν) ≤ 3 the frontier compresses 1.0x, so it does a naive
   enumerator's work plus hashing. Per-output counting fixes that *above* a
   crossover (`src/two_row.rs`, `src/three_row.rs`, both dispatched), but
   below it the fibre DP's own cost dominates and the frontier still wins.
   Closing the rest needs either a cheaper fibre count or a lower crossover;
   note that removing the frontier outright was prototyped and is parity at
   best, so that door is shut.
3. **Extend counting to four-row factors.** The state gains one dimension per
   strip, so ℓ(ν) = 4 needs (λ¹ⱼ, aⱼ, bⱼ, cⱼ). Whether that stays affordable
   is unknown — the three-row case cost 26–2192 ops/term against a predicted
   "thousands, hopeless", so the bounding-box estimate is not trustworthy here
   and it should be measured rather than reasoned about. `[24,20,16,12]²`, our
   largest case, has four-row factors.
4. **Shape preprocessing** — factoring a skew diagram into connected
   components and expanding each separately, since the expansion of a
   disconnected shape is the product of its pieces.
5. **Output residency.** On `[24,20,16,12]²` a growing share of peak RSS is
   the 5.3M-term *output* (the memoized `Arc<Vec>` plus the caller's clone),
   not the frontier. An `Arc`-returning variant of `expand_skew` would halve
   that.

---

## Beyond the core (deferred, but intended)

The target above is the symmetric-function core. Symmetrica — the library this
one succeeds in spirit — also covers, and symfn does not:

- modular and projective representation theory of the symmetric group
- Schubert polynomials, commutative and non-commutative
- Hecke algebras of type A
- finite group operations
- ordinary representation theory of the classical groups

That is the eventual scope, not the current one. Nothing here is scheduled and
none of it should be read as implied by the phases above; the core comes first
and is where all the depth work (LR backends, Kostka, conversions) lives.

Recorded now for one practical reason: **Symmetrica is public domain**, so when
this work does begin it is a legitimate source of algorithms, not merely a
reference point — see `NOTICE.md`. It would also serve as a third test oracle
alongside Sage and lrcalc, with no licensing friction and coverage of
operations neither of those makes convenient.

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
