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

| product | lrcalc | SkewLr | |
|---|---|---|---|
| s[8,7,6,5,4,3]² | 8.55s | 1.31s | **us 6.5x** |
| wide [16,13,10,7]² | 14.46s | 5.99s | **us 2.4x** |
| s[6,5,4,3,2,1]² | 0.073s | 0.034s | us 2.2x |
| s[9,8,7,6,5]² | 1.098s | 0.653s | us 1.7x |
| s[7,6,5,4,3]² | 0.079s | 0.052s | us 1.5x |
| s[8,7,6,5,4]² | 0.349s | 0.272s | us 1.3x |
| rectangle [7^7]² | 0.023s | 0.024s | tie |
| wide [14,12,10]² | 0.029s | 0.036s | lrcalc 1.3x |
| wide [12,10,8]² | 0.023s | 0.030s | lrcalc 1.3x |
| wide [20,16,12]² | 0.236s | 0.358s | **lrcalc 1.5x** |
| wide [24,20,16,12]² | >60s | >60s | neither finishes |

Noise is ±30% run to run; treat anything inside ±20% as a tie. lrcalc's own
timing on an unchanged binary drifted 6.3s → 8.6s → 11.5s across this project's
sweeps, so single-digit-percent differences mean nothing.

**The remaining weakness is moderate-size wide shapes** — few rows, large parts,
in the 0.02–0.4s band — where we lose ~1.3–1.5x. It does *not* extend to large
wide shapes: `[16,13,10,7]²` is wide and we win it 2.4x. So this is a
constant-factor problem at moderate size, not a scaling one. Few rows means
little merging, so the frontier's hashing and allocation overhead is paid
without collecting its benefit, against lrcalc's very tight per-tableau loop.

**Two findings worth keeping.**

*Transposition does not fix wide shapes.* Since c^λ_{μν} = c^{λ'}_{μ'ν'}, the
walk can run on the transposed diagram, and "few rows become few columns" looks
like the obvious fix. Measured, it does the opposite: wide shapes prefer the
*original* orientation ([20,16,12] 0.65x, [12,10,8] 0.57x, [20,10] 0.15x), and
only staircases mildly prefer the conjugate — i.e. it helps only where we
already win handily. What actually fixed wide shapes was collapsing the row fill
into runs.

*The biggest shapes are memory-bound, not compute-bound.* `[24,20,16,12]²` was
killed at 27m27s wall with 2.0+ GB resident and still climbing, **38% of it
system time** — allocation and page faults, not combinatorics. Per-transition
optimization cannot reach it; the state set itself has to shrink. One candidate
was ruled out: pruning states whose `above` row exceeds the next row's value cap
never fires, because any value `w` in `above` was actually placed, so
content[w−1] ≥ 1 and hence w ≤ content.len() < cap, always.

**Next**, in priority order:

1. **Shrink the state set.** This is the real barrier and it needs a different
   decomposition, not tuning. Everything large is bounded by it.
2. **Moderate wide shapes** — the last regime where lrcalc beats us, now a
   ~1.3–1.5x constant factor rather than 3.4x.
3. **Shape preprocessing** — factoring a skew diagram into connected components
   and expanding each separately, since the expansion of a disconnected shape is
   the product of its pieces. `SkewLr` already exploits that fact in one
   direction (to *build* a product); using it in reverse, to decompose, should
   show up most on skew inputs with gaps.

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
