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
| s → h | C | 4.8x | 4.6x | 4.2x | 4.0x |
| **s → p** | C | 3.1x | 3.5x | 2.1x | **2.0x** *(was 1.2–1.6x)* |
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

| case | Symmetrica C | symfn (before) | | symfn (now) | |
|---|---|---|---|---|---|
| s_3[s_{21}] | 0.000098s | 0.000133s | 0.79x | 0.000055s | **1.78x** |
| s_4[s_{21}] | 0.000282s | 0.000462s | 0.64x | 0.000111s | **2.54x** |
| s_5[s_{21}] | 0.000824s | 0.001985s | 0.42x | 0.000368s | **2.24x** |
| s_3[s_{31}] | 0.000197s | 0.000657s | 0.30x | 0.000120s | **1.64x** |
| s_4[s_{22}] | 0.000572s | 0.004050s | 0.14x | 0.000551s | **1.04x** |

**Ahead on every case**, from behind on every case. Against Sage, plethysm went
9x -> ~40x.

Plethysm was **slower than Symmetrica and the gap widened with output size** —
a scaling problem hidden because the visible baseline was an interpreter. Three
fixes, all in `p → s`, which profiling showed was 98–99% of plethysm's runtime
(the module doc's claim that cost was split between "the two conversions at the
ends" was wrong: s → p is ~0.5%):

1. **The DP now accumulates in i128, not the caller's ring.** Every value in it
   is a character, but plethysm runs over ℚ, so each rim hook was doing a
   rational add — a gcd — plus a temporary from negating. Safe because the β
   mask already caps l at 32 and √(32!) ≈ 1.6·10¹⁸.
2. **Terms are batched by degree and share their sweep.** Expanding each p_μ
   separately from ∅ discards all the Murnaghan–Nakayama work that p_μ and p_ν
   share whenever they share parts; the batch sorts by part sequence and
   continues one frontier per common prefix.
3. **Accumulation is keyed on the β-mask, not the partition.** Every leaf
   touches the whole frontier, so a partition key allocated, sorted and hashed a
   fresh `Vec` once per (μ, mask) pair — ~9,000 allocations to produce 63 terms.

A fourth, worth its two lines: **`p_step` pre-sizes its output map.** The
frontier grows monotonically through a sweep, so a default-capacity map rehashed
several times per step. Worth 1.22x, measured over 8 interleaved rounds.

Net ~3–4x, and it moved plethysm from 9x to ~25x against Sage. **We are now
ahead on three of five cases, at parity on a fourth, and behind only on
`s_4[s_{22}]` (0.50x)** — the scaling issue is reduced, not eliminated.

#### ⚠️ Retracted: "the rational leaf arithmetic is not the bottleneck"

This section previously recorded a *negative* result — that killing the rational
arithmetic at each (μ, mask) leaf had only a 24% ceiling, and so was not worth a
`Ring` hook, an lcm with overflow guards, and a fallback path. **That conclusion
was wrong, and it was wrong because the measurement was a single un-interleaved
run.**

A sampling profile (2,496 samples) said otherwise:

| symbol | self | share |
|---|---|---|
| `p_step` (frontier DP) | 852 | 34% |
| `u128_div_rem` (gcd) | 733 | 29% |
| `Rational::add_assign` | 366 | 15% |
| `Rational::mul` | 146 | 6% |
| `__modti3` / `__divti3` | 139 | 6% |

Rational arithmetic was **~55%**, not 24%. Re-running the same ceiling
experiment over 6 interleaved rounds gave 0.000875s → 0.000434s: a **2.02x**
ceiling. The original number came from one run of each side, on a machine that
has repeatedly been shown to drift 2x as it warms — the exact failure this
document had already warned about two paragraphs earlier, committed anyway.

The lesson is not "always build it". It is that a *cheap* experiment used to
**cancel** work needs the same rigour as one used to justify it, and it did not
get it. A profiler would have settled it in one shot for less effort than the
experiment cost.

#### The fix that followed: a common-denominator integer sweep

`p → s` computes Σ_μ c_μ χ^λ(μ) — a sum of (coefficient × integer) terms. Over
ℚ that is a rational multiply and a rational add per leaf, each normalising by a
gcd. Putting every c_μ over one denominator D makes the whole accumulation
integer, with a single conversion back per output term.

D is the lcm of the denominators, and measurement said that is the right shape:
across these plethysms **the lcm equalled the largest denominator every time**,
never exceeding ~5·10⁵. Overflow is still checked at every step rather than
argued away, since the guarantee only covers the cases measured, and any failure
falls back to the generic path having written nothing.

Exposed as two provided `Ring` methods (`as_ratio` / `from_ratio`) defaulting to
`None`, so rings that cannot answer — `i128`, and the `gmp` types — simply keep
the existing path.

**3.0x on p → s**, which beat the "ceiling" above because that experiment still
built a `Map<u64, C>` and still called `Rational::new` on denominator-1 values.

### Parallel LR (`src/skew_lr.rs`)

The frontier traversal is now multi-threaded. A row is a barrier — row r+1 needs
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
The frontier is now *sharded*: each key is routed to a shard by a cheap hash of
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

⚠️ **These numbers are AC-only, and that is not a formality.** An earlier run of
the same A/B on battery reported 2.02x where AC says 1.73x for the identical
binaries — the power state changes the *ratio*, not just the absolute times,
because throttling hits ten busy cores differently from one. Parallel results
measured on battery are not comparable to anything.

### The internal (Kronecker) product

The third product on symmetric functions, after the ordinary one and plethysm
(Macdonald I.7). Under the Frobenius characteristic it is the tensor product of
S_n representations: `s_λ * s_μ = Σ_ν g^ν_{λμ} s_ν` with
`g^ν_{λμ} = ⟨χ^λ χ^μ, χ^ν⟩`. The Kronecker coefficients are famously harder than
Littlewood–Richardson and have **no known positive combinatorial rule**.

Computing them is nonetheless nearly free, because the internal product is
*diagonal in the power-sum basis*: `p_λ * p_μ = δ_{λμ} z_λ p_λ`. So it is s → p
on both sides, a coefficientwise multiply weighted by z_λ, and p → s back —
nothing enumerates anything. That a hard object falls out of a diagonal basis is
the payoff for keeping the power-sum route fast.

**434/434 products agree with Sage's `itensor` across degrees 1–7**, and the
implementation is checked in-crate against the character-theoretic definition
`g^ν_{λμ} = Σ_ρ χ^λ(ρ)χ^μ(ρ)χ^ν(ρ)/z_ρ` for every (λ, μ, ν) triple through
degree 6. That second check is the one that matters: the implementation rests
entirely on the p-basis identity plus s ↔ p, so a wrong identity or a wrong
conversion would produce a self-consistent but false answer. Symmetries
(g invariant under permuting the three indices, and under conjugating any two)
and the dimension count Σ_ν g·dim ν = dim λ · dim μ are asserted too.

| degree | products | Sage | symfn | ratio |
|---|---|---|---|---|
| 8 | 36 | 0.1057s | 0.0019s | 55.5x |
| 10 | 36 | 0.3854s | 0.0047s | 82.2x |
| 12 | 36 | 1.1893s | 0.0130s | 91.2x |

⚠️ That is against **Sage's Python**, not C. Symmetrica exposes no internal
product at all, so unlike the conversions there is no compiled baseline to check
against — and a `py` ratio is exactly the kind that read 9x for plethysm while
the truth against C was 0.14x. Treat 55–91x as "not obviously slow", not as a
result.

### The forgotten basis: Macdonald's sixth

`f_λ = ω(m_λ)` (Macdonald I.2). This is the last of the six classical bases and
the only one with no independent combinatorial description — it is *defined* as
the image of the monomial basis under ω, which is where the name comes from.
Adding it makes `convert` total over the standard set.

The implementation is two lines each way, because ω is an involution:

```text
  to_schur:    Σ c_λ f_λ = ω(Σ c_λ m_λ)         →  ω(monomial_to_schur(c))
  from_schur:  x = Σ c_λ f_λ ⟺ ω(x) = Σ c_λ m_λ →  monomial coeffs of ω(x)
```

and ω on the Schur basis is conjugation of every index. So both directions
inherit Muir's rule and the Kostka machinery — and their asymptotics — for one
transpose per term.

Being nearly free makes the *testing* the real work, since a test phrased in
terms of ω would only restate the implementation. Three checks that don't:

* **Endpoints, hand-derived.** `f_(n) = (−1)^{n−1} p_n` and `f_(1^n) = h_n`,
  from `m_(n) = p_n` and `m_(1^n) = e_n` pushed through ω. Neither mentions the
  forgotten basis; both are facts about the other bases.
* **Duality**, the structural characterisation: `⟨f_λ, e_μ⟩ = δ_{λμ}` for every
  pair through degree 7, reached through `hall` and the e → s expansion — code
  the conversion never touches. Wiring `Forgotten` to the wrong involution, or
  to conjugation of the *index* rather than of the Schur expansion, would fail
  this while a round trip still closed.
* **Round trips** out through m, e, and h and back, all partitions to degree 7.

Both directions agree with Sage at degrees 8 and 12.

⚠️ The measured ratios are 146x/989x (s → f) and 379x/888x (f → s), and they
should not be quoted. Symmetrica has **no forgotten basis** — verified, not
assumed: `sage.combinat.sf.classical.conversion_functions` holds exactly 20
entries, the ordered pairs among {Schur, elementary, homogeneous, monomial,
powersum}, and forgotten appears in none. So Sage falls back to a generic
Python basis-change through its own machinery, and this is a `py` row against an
unoptimised path. It says the basis is not a bottleneck; it says nothing more.

### The Python boundary's integer ceiling — decided: compute-and-escalate

Two corrections to what this file previously implied. **`gmp` and `python` do
compose** — `maturin build --features "gmp,python"` produces a wheel; that note
was stale. But enabling `gmp` changes *nothing* about the Python API, because
`python.rs` hardcodes `i128` in 20 places and every entry point builds
`Schur<i128>` or the i128-backed `Rational`. The GMP wheel is behaviourally
identical to the plain one.

The forcing issue is soundness rather than capability. `impl Ring for i128` uses
plain `*`, so in release **coefficient arithmetic wraps silently**. Characters
are already guarded (`character()` panics, `try_character` returns `None`) and
`integral_sweep` is fully checked with a clean bail-out — but nothing protects
generic coefficient arithmetic in products, plethysm, or conversions. At the
degrees a Sage user reaches that returns wrong answers with no signal, which is
not something a Symmetrica replacement can ship.

Three options: (A) keep i128 and refuse loudly, (B) always bignum, (C) compute
in i128 with checked arithmetic and re-run the whole call in `rug::Integer` if
anything overflowed. C is the crate's existing pattern — `try_character` →
`character_in` — and keeps the fast path fast, but only if "checked" is nearly
free. `examples/bench_guarded.rs` measures exactly that: `i128` and `Rational`
with every arithmetic op replaced by its `checked_` form, over Schur products,
s→m/h/e, m→s, s→p, plethysm and the internal product.

**Checked arithmetic costs 0–1%, indistinguishable from run-to-run noise.**

| workload | plain | checked | cost |
|---|---|---|---|
| integral (products, s→m/h/e, m→s) | 2.428s | 2.422s | −0.2% |
| dividing (s→p, plethysm, Kronecker) | 0.095s | 0.095s | −0.3% |

So **C is decided**. Escalation is the rare path, the common path is unchanged,
and the ceiling stops being visible to callers.

⚠️ **The methodology mattered more than the result, and this is the third time.**
With a fixed pass order the benchmark first reported checked rationals as **36%
faster** than unchecked. Swapping the two blocks moved the 36% to the other
type: whichever rational pass runs immediately after the (25x larger) integral
passes pays ~60% for arriving with a cold cache, and a fixed order silently
charges that to one type. Interleaving alone was not enough — the fix is
**rotating** the order so each variant spends an equal share of rounds in each
position. Two hypotheses were tested and discarded on the way: missing
`#[inline]` on `Rational`'s `Ring` impl (adding 25 of them changed nothing,
since generic instantiation already inlines) and a difference between `Rational`
and the hand-written twin (a byte-identical unchecked twin in the same crate
showed the same anomaly, which is what localised it to position).

#### What the escalated ring needs to be: digits, not speed

`examples/coeff_sizes.rs` measures the two things that decide it.

**Coefficients are small.** Across LR squares, Kostka rows, characters, plethysm
numerators and s → p denominators, nothing in these workloads exceeds **2
limbs**, and most are one:

| workload | widest coefficient |
|---|---|
| LR `s_[10,9,8,7,6,5]²` | 18 bits |
| Kostka row, degree 32 | 55 bits |
| characters `p_1^40 → s` | 76 bits |
| f^λ at degree 60 | 118 bits |
| s → p denominators, degree 24 | 43 bits |
| plethysm numerators | 2–5 bits |

Even f^λ for a degree-60 staircase — 36 digits — still fits `i128`. Past the
ceiling we are therefore in the 2–5 limb range, where **every** bignum library
runs schoolbook: Karatsuba engages around 10–30 limbs, Toom near 100, FFT in the
thousands. GMP's advantage at this size is assembly tuning, not algorithms.

**Coefficient arithmetic is a small share of runtime.** Replacing `i128` with
heap-allocated `rug::Integer` *entirely* costs **1.04x** on products and s → m,
so at most ~4% of that workload is coefficient arithmetic; the combinatorial
traversal dominates. The rational path is the exception at **3.39x**, since
`rug::Rational` pays a gcd per operation — but on a 100x smaller absolute base,
and its hot case (`integral_sweep`) already runs in raw `i128` under a common
denominator.

So the requirement is **exactness, not speed**, and a pure-Rust bignum is
enough. That also settles the wheel: `rug` and `gmp-mpfr-sys` are **LGPL-3.0+**,
and statically linking them into a distributed wheel would attach LGPL terms to
a crate that is deliberately MIT OR Apache-2.0. `num-bigint`/`num-rational` are
permissive, need no C toolchain or `m4` to build, and PyO3 0.29 has
`num-bigint`/`num-rational` features that convert them to Python `int`/`Fraction`
natively — which removes the decimal-string encoding step as well. `gmp` stays
an optional feature for Rust callers who want it.

Implementing C, in order: a guarded coefficient type with an overflow flag; the
`python.rs` entry points made generic over it and re-run in `BigInt`/`BigRational`
on the flag; and the boundary handing those to PyO3 directly, so callers always
see a plain arbitrary-precision Python `int` with no mixed-type list to branch on.

### Coefficient rings that are not fields (ℚ[t], ℚ[q,t])

Prompted by the question of what it would take to serve a Sage user working over
`QQ['t']`. Most of the answer is "nothing": **Sage already factors the
coefficient ring out of the conversion path.** From `sage/combinat/sf/classical.py`:

```python
if R == QQ and P.base_ring() == QQ:
    return self._from_dict(t(m)._monomial_coefficients, coerce=True)   # one bulk call
f = lambda part: self._from_dict(t({part: ZZ.one()})._monomial_coefficients)
return self._apply_module_endomorphism(x, f)                            # integer rows, R-recombine
```

Over any ring that is not ℚ, Symmetrica is asked only for the **integer
transition row of a single basis element**, and Sage does the arithmetic in R
itself. Structure constants are ring-independent, so `python.rs` already exposes
what is needed. (The non-QQ path is one FFI call per partition in the support,
which a bulk entry point would collapse — but measured at degree 12, 77 terms,
`QQ['t']` took 0.053s against `QQ`'s 0.066s, so this is structural tidiness and
not a bottleneck.)

What *did* block ℚ[t] was our own bound. **Every division in this library is by
z_μ, a positive integer** — `s → p` carries z_μ⁻¹ and the internal product and
plethysm inherit it through the power-sum basis. Nothing divides by a general
ring element. But those paths were bounded on `Field`, which demands the ability
to invert *t*, so ℚ[t] and ℚ[q,t] were excluded from operations that are
perfectly well defined over them — and those are exactly the coefficient rings
Hall–Littlewood and Macdonald need.

`Field` is replaced in those bounds by `QAlgebra` — a ring containing ℚ — with
the single method `div_u128`. `Field` is kept (it is a real thing to name, and
`Rational` is one) but nothing in the library requires it. The implication
"characteristic-0 field ⟹ ℚ-algebra" is deliberately *not* a blanket impl: that
would occupy the impl for every downstream type, and a polynomial ring could
then never implement it.

**Plethysm needed a second, less obvious fix, and it was silently wrong.**
Everything else in the library is linear with integer structure constants, so a
coefficient is only ever multiplied and added. `p_n` is different: it substitutes
into the alphabet, and a coefficient ring's variables belong to that alphabet.
Sage:

```text
sage: R.<t> = QQ[];  p[2](t*p[1])  ->  t^2*p[2]        p[2](t*p[1], exclude=[t]) -> t*p[2]
```

`scale_parts` carried coefficients through unchanged — correct over ℚ, where
there is nothing to raise, which is why the omission was invisible for as long
as ℚ was the only coefficient ring in use. It is now `Plethystic::frobenius`, a
separate trait above `QAlgebra` so that `s → p` (which divides but never
substitutes) does not demand it. Which variables get raised is a genuine
convention choice, so it is the implementor's to state rather than a default.

`tests/qalgebra.rs` is the proof: a ℚ[t] implementing `Ring + QAlgebra +
Plethystic` and deliberately **not** `Field`, exercising `s → p`, the round
trip, the Kronecker product, and plethysm. If any of those paths regressed to a
`Field` bound the file would stop compiling, which is a stronger assertion than
its `assert_eq!`s. Plethysm values are Sage's, including the degree-3 cases
(`s[3](t*s[1]) = t^3*s[3]`) that a doubling bug could not fake. Over ℚ nothing
changed: the full Sage ladder still agrees.

Remaining for a real Sage backend, in order: arbitrary-precision integers across
the FFI boundary (the known `i128` ceiling, plus `--features gmp` and
`--features python` still not composing), then the bulk expansion entry point.

### Skewing by an arbitrary symmetric function

`g^⊥`, the adjoint of multiplication by g under the Hall inner product:
⟨g^⊥ f, h⟩ = ⟨f, g·h⟩. `skew_schur` was only the case g = s_μ, f = s_λ — which
is why the classical notation for it is a quotient, s_μ^⊥ s_λ = s_{λ/μ}.

`SkewBy<C, G>` is generic over **the basis g is written in**, and that is the
design rather than a convenience. g ↦ g^⊥ is linear, so any g could be expanded
into Schur and handed to Littlewood–Richardson — but three bases have adjoints
with direct rules, each a Pieri or Murnaghan–Nakayama rule read backwards:

| g in basis | g^⊥ on s_λ | machinery |
|---|---|---|
| h | remove a horizontal strip | Pieri |
| e | remove a vertical strip | dual Pieri |
| p | remove a rim hook, signed by height | M–N |
| s, m, f | Σ_μ d_μ s_{λ/μ} | LR |

`e` is implemented as ω ∘ h^⊥ ∘ ω rather than as a second strip enumerator: ω
is an isometry with ω(h_r) = e_r, so ⟨e_r^⊥ s_λ, s_ν⟩ = ⟨h_r^⊥ s_{λ'}, s_{ν'}⟩.
An identity, and it costs one transpose per term against a duplicate enumerator
with its own separate bugs. `p` reuses `character::border_strips` — the same
β-number bit tricks that make the character table fast.

**The p case is the one that most repays the native path, and not only for
speed: it needs no division.** Expanding p_μ into Schur requires 1/z_μ, so the
LR route forces ℚ. Rim-hook removal stays in ℤ, so `Schur<i64>` can be skewed
by a power sum — something the generic route could not have offered at all.

Native path vs the same g expanded into Schur, interleaved in one process:

| g | shapes | ratio |
|---|---|---|
| p | λ of 4–7 rows, degree 42–52 | **43–268x** |
| h | same | 2.2–3.3x |
| e | *tall* λ (10–12 rows) | 3.0–7.2x |

⚠️ On battery. Same-process single-threaded ratios are far less power-sensitive
than the parallel LR numbers were, but they are not re-measured on AC.

The `e` row needed the shape family changed, and the first attempt is worth
recording as a measurement error rather than a result. On the *wide* shapes used
for h and p it read 0.8–1.6x — because a vertical strip needs a distinct row per
cell, so on a 4-row λ the answer is nearly empty and both routes finish in under
a millisecond. The ratio was measuring harness noise on a trivial answer, not
the algorithms. Tall shapes make the operation non-trivial and it behaves like h.

Verified against Sage (`scripts/check_skew.py`): **32448 checks through degree
8, zero mismatches**. Every (λ, μ, basis) is put to symfn twice — once in its
own basis, once expanded into Schur — so agreeing with Sage says the answer is
right and agreeing with each other says the fast path is a shortcut and not a
different operation. In-crate, the defining adjointness ⟨g^⊥ f, h⟩ = ⟨f, g·h⟩ is
checked for every triple through degree 6; that test mentions no algorithm at
all, and its two sides share no code.

### Evaluation at an alphabet, and the principal specializations

`src/eval.rs`. Everything else in the crate computes *with* symmetric functions
as formal objects; this is the bridge back to concrete numbers. Two different
things live there and the distinction is the design:

**Evaluation at an arbitrary alphabet**, generic over `Ring`, one algorithm per
basis rather than "convert, then evaluate". p, e, h are products of one-row
generators and cost a linear DP each. Schur uses the **branching rule**: a
tableau is a chain ∅ = ν⁰ ⊆ … ⊆ νⁿ = λ of horizontal strips, so sweeping
variable by variable with a frontier of *shapes* collapses every tableau sharing
a prefix into one number. Cost is the shapes inside λ, not the tableaux — the
same chain DP as `kostka.rs`, carrying ring elements instead of counts.

**The bialternant is deliberately absent.** s_λ = a_{λ+δ}/a_δ is the textbook
formula and would be an O(n³) determinant, but it needs *division* — so it is
not generic over `Ring` — and it is 0/0 whenever two x_i coincide. The branching
rule is slower on generic input and always right. The oracle script exercises
exactly that hole: one of its alphabets is `[2,2,2,0,5]`.

**Closed forms** for the two special alphabets, which enumerate nothing:

| | formula |
|---|---|
| `dimension(λ)` = f^λ | \|λ\|! / ∏ h(u) |
| `principal_specialization(λ,n)` = s_λ(1ⁿ) | ∏ (n + c(u)) / h(u) |
| `principal_specialization_q(λ,n)` | q^{n(λ)} ∏ (1−q^{n+c(u)}) / (1−q^{h(u)}) |

The q-analogue returns a coefficient vector. Neither product divides the other
cell-by-cell, so the quotient is taken once at the end; both have constant term
1, which makes it a truncated power-series inversion — no leading-coefficient
case analysis, and exact in ℤ because the quotient is known in advance to be a
polynomial.

`dimension` and `principal_specialization` interleave their divisions with their
multiplications rather than forming the factorial first, and that is not a
micro-optimisation: for the staircase λ = (10,9,…,1), f^λ has **35 digits** and
fits `u128`, while 55! has **74**. Forming the numerator first would overflow by
thirty-five orders of magnitude on an answer that is comfortably representable.
Both return `Option`, `None` on genuine overflow.

Verified against Sage (`scripts/check_eval.py`): **1440 checks through degree 9,
zero mismatches** — evaluation against `expand(n)` substituted, both
specializations against `principal_specialization`, and f^λ against Sage's own
tableaux count. In-crate, the five bases are checked to agree with each other at
a shared alphabet through degree 6, which is the check that catches a wrong
recurrence in any one of them: they share no code, so they can only agree by all
being right.

No timings are quoted. This is new surface rather than a faster route to
something we already had, and Symmetrica has no equivalent entry point.

### Memory: two thirds of RSS is allocator retention, not data

`examples/lrheap.rs` wraps the global allocator to count live bytes, which
separates what the traversal actually holds from what the process has not given
back. The two differ by a lot, and the difference grows with the shape:

| shape | | live heap | peak RSS | retention |
|---|---|---|---|---|
| `[8,7,6,5,4,3]²` | serial | 84.5 MB | 157.8 MB | 1.9x |
| | parallel | 86.4 MB | 168.7 MB | 2.0x |
| `[16,13,10,7]²` | serial | 122.6 MB | 364.4 MB | **3.0x** |
| | parallel | 123.2 MB | 364.1 MB | **3.0x** |

Live data is ~66 bytes per frontier state, which is about right for a 40-byte
`(Key, u64)` entry plus hash-table slack — the representation is not the
problem. **RSS is 3x that because every row allocates a fresh frontier and frees
the previous one**, and after 32 rows of multi-megabyte churn the allocator is
holding the difference. That reframes the standing "we use 1.4–4.9x lrcalc's
memory" line: on live data the gap is far smaller, and most of what was being
compared is retention.

#### ⚠️ Tried the obvious fix; it made things worse

Carrying the frontier tables across rows and `clear()`ing them — keeping
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
Isolating the two halves showed the tables, not the frontier vectors, were
responsible (399 MB with only the tables pooled).

So the retention is real but it is *not* free to reclaim, and the naive reading
— "reuse the allocations" — is wrong here. Worth knowing before anyone tries it
again.

What is still untried, and is a genuine reduction rather than a reshuffle: the
per-shard entry buffers and the output vector are both live at once during the
merge, holding the frontier twice. Having the shards write into disjoint ranges
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

### Against Symmetrica directly (`scripts/compare_symmetrica.py`)

`compare_sage.py` can only see half of what it measures: Sage's five classical-
basis conversions dispatch into Symmetrica's C, but its *other* operations are
pure Python — even where Symmetrica ships a C implementation Sage never calls.
Plethysm read as **9x faster** than Sage while being **0.14x** against
`sage.libs.symmetrica.all.plethysm` on the same inputs. That defect sat inside a
green benchmark.

This harness drives Symmetrica directly through the bindings Sage installs but
mostly leaves unused, and verifies every case rather than only timing it — which
also makes it a second independent oracle for LR alongside lrcalc.

| case | deg 10 | deg 12 | deg 14 |
|---|---|---|---|
| LR product s_μ·s_ν | 3.5x | 3.8x | **4.7x** |
| skew s_{λ/μ} | 3.4x | 3.6x | 4.1x |
| Kostka K_{λμ} (one value) | 3.0x | 4.8x | 5.1x |
| character χ^λ(μ) (one value) | 7.9x | 0.6x | 0.5x |
| Hall ⟨s_λ, s_λ⟩ | 8.0x | 8.2x | 9.4x |
| **character table** | — | 3.6x *(was 0.6x)* | **2.9x** *(was 0.7x)* |
| **Kostka table** | — | 3.6x *(was 0.5x)* | **2.5x** *(was 0.39x)* |
| plethysm (single-row outer) | 1.5x | 1.7x | 2.1x |

At degree 16: character table **3.8x**, Kostka table **2.8x**, LR product 5.3x.

**LR is comfortably ahead of Symmetrica** — and the small-degree number badly
understates it. The ladder above tops out where a whole product is ~0.1 ms,
which says nothing about the regime this library targets. At real sizes:

| product | degree | Symmetrica | symfn | ratio |
|---|---|---|---|---|
| `[4,3,2,1]²` | 20 | 0.0026s | 0.0008s | 3.2x |
| `[6,5,4,3]²` | 36 | 0.259s | 0.0028s | **91x** |
| `[8,7,6,5]²` | 52 | 7.50s | 0.0170s | **442x** |
| `[10,8,6,4]²` | 56 | 530s | 0.0589s | **9002x** |

Every coefficient agrees, which also gives the LR engine a second independent
oracle alongside lrcalc — the first it has had. Symmetrica's
`outerproduct_schur` degrades violently over this range (0.0026s → 530s while
symfn goes 0.0008s → 0.059s), so `run_big_lr` stops as soon as it passes the
budget. Note the direction of the lesson: at toy sizes this reads 3.2x, and
sizing the benchmark where the work actually lives changed the answer by three
orders of magnitude.

**Two new deficits, both on whole tables.** Note the shape of it: our *per-value*
Kostka and character are 3–5x faster, but the *whole table* is 2–2.5x slower and
the Kostka gap widens with degree (1.4x → 0.51x → 0.39x). That is the signature
of Symmetrica computing a table **as a table**, sharing work across entries,
while we answer p(n)² independent memoized queries. It is exactly the "produce
the whole answer in one sweep rather than query it entry by entry" pattern this
library has already applied to Kostka rows, the coproduct, p → s and m → s — and
has not applied to either table.

**Both are fixed, by the same observation: a table is p(n) *sweeps*, not p(n)²
numbers.**

* `kostka` bounds its chain DP by λ and reads one entry out of the final
  frontier, discarding everything else that frontier holds. Drop the bound and
  the frontier at the end of μ's chain **is** the whole column — every λ with its
  K_{λμ} — for barely more than the single-value cost.
* `p_expand` already computed p_μ = Σ_λ χ^λ(μ) s_λ in one Murnaghan–Nakayama
  sweep, which *is* a column of the character table. The machinery was there; it
  had simply never been pointed at the table.

Columns then share with each other. K_{λμ} depends on μ only as a multiset, so
its parts can be consumed in any order, and characters likewise — taking them
descending lets partitions with a common prefix share the whole initial segment
of their chain. One traversal covers every μ. That is the same trie as
`convert::p_expand_shared`, which both now use.

Roughly **6–7x** on each table, turning both from losses into 2.5–3.8x wins.
Verified entry-by-entry against the per-pair functions through degree 12 (Kostka)
and 11 (characters) — a prefix-grouping slip would misattribute a column, and the
β-mask row index would drop rows rather than corrupt them, neither of which a
spot check catches.

⚠️ Single-value character rows above are unreliable and should not be read as a
deficit on their own: at 20–50 µs they are dominated by harness overhead, and
Symmetrica caches internally across calls while symfn calls `clear_caches()`
before each timed one. The *table* rows are the trustworthy comparison, which is
why they exist — a big symmetric unit of work with no caching asymmetry.

#### The comparison above is on Symmetrica's home turf

Symmetrica's plethysm **refuses a multi-row outer partition** — it reports "for
the moment only for outer S_n" and computes nothing. Every case in the table is
therefore a single-row outer, the easy case, and possibly a specialised path.

symfn has no such restriction. `s_{21}[s_{21}]` (17 terms, 0.0011s),
`s_{22}[s_2]`, `s_{32}[s_{11}]`, `s_{21}[s_{31}]` (39 terms) and
`s_{31}[s_{22}]` (104 terms) all compute and all agree with Sage — 48/48 cases
across both shapes of input. So the honest summary is that symfn is faster than
Symmetrica on the inputs Symmetrica accepts, and is the only one of the two that
handles the rest.

#### s → p: fixed, but not the way expected

`PowerSum::from_schur` calls `character_in(λ, μ)` once per μ, which looked like
the pattern removed from `p → s` — p(n) independent queries where one sweep
would do. Two measurements redirected the work:

1. **The obvious fix is a loss.** `p_expand` produces a *column* of the character
   table (one μ, all λ); `s → p` needs a *row*. Running the batched column sweep
   over every μ of the degree and reading off one row costs 0.0398s at degree 20
   against 0.0117s for the three shapes it would serve — **3.4x worse**, with
   break-even only around ten shapes.
2. **Characters were 84–100% of `s → p`**, so the transition arithmetic was never
   worth touching.

The real cost was representation, not algorithm. `border_strips` runs at *every
node* of the recursion and allocated a `Vec` for β, a heap **`HashSet`** for
membership, then per strip another `Vec` plus a sort — thousands of heap
allocations per character, to do arithmetic that fits in registers. Holding the
β-set in a u64 makes membership a bit test and the height a masked
`count_ones`. Worth **1.5x** at every degree (1.51x, 1.57x, 1.51x, over 6
interleaved rounds), taking `s → p` from 1.2–2.0x to **2.0–3.5x**.

Same β-number mathematics as before, and the general form is retained both as
the fallback past |λ| = 32 and as the reference the masked path is checked
against. That test earned its place immediately: asserting on strip *order*
failed, because the masked form walks β upward and the general one downward —
invisible to callers, which only sum, but noticed rather than assumed harmless.

Still open: `character_cached` clones both partitions into its key and takes a
global `RwLock` at every node, and `character_uncached` allocates a fresh
`Partition` for the μ-suffix at every node. Interning shapes to integer ids
would remove both.

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
