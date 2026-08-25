# Vertical-strip and ribbon LLT polynomials

Three engines over `QtPoly` — SYT enumeration, beta-set ribbon strips, and
Fock-space straightening — sharing nothing but the coefficient ring and
cross-checking each other. Sage's pure-Python implementation is the only
other one anywhere and walls at whole-degree n = 10–11; Symmetrica has no LLT
at all, a capability gap rather than a backend swap. Speed against Sage runs
from ~700× (whole-degree tables) to ~19 100× (k = 4); two
deliverables — the `∇e_n` by-path Schur-positive refinement and
Kazhdan–Lusztig polynomial columns by exact straightening — have no Sage
entry point at any speed.

---

`src/llt.rs` — the ribbon model ([LLT] §6) and the tuple model ([HHL] Def
3.2), and the dictionary between them, pinned by 1517 comparisons against
Sage at a degree-5 dump (`examples/llt_dump.rs` → `scripts/check_llt.py`, 0
failures) and 30 in-crate tests (plus 4 for the abacus primitives in
`partition.rs`), all green. Python bindings are 16 functions in `python.rs`,
checked by 29 boundary tests in `scripts/check_bindings.py`, 0 failures.

## Where Sage stops, and where the crate already is

This machine, 2026-07-30, against `llt_h_table` (`i64`) — the authoritative
ratios; the battery walls below are provenance, not arithmetic.

Whole-degree spin tables (every μ ⊢ n, `HSp_k[μ] → s`):

```text
  k   n      Sage      ours     ratio     Sage worst shape μ=(n)
  2   7      0.33s   0.0005s      660×     0.31s  (94% of degree)
  2   8      1.08s   0.0012s      900×     1.03s  (95%)
  2   9      4.84s   0.0031s     1560×     4.73s  (98%)
  2  10     24.16s   0.0078s     3100×    23.86s  (99%)
  2  11     >120s    0.0190s    >6300×        —
  3   7      1.11s   0.0009s     1230×     1.08s  (98%)
  3   8      6.77s   0.0028s     2420×     6.70s  (99%)
  3   9     47.34s   0.0080s     5920×    47.17s  (100%)
  3  10     >120s    0.0222s    >5400×        —
  4   7      8.29s   0.0017s     4880×     8.26s  (100%)
  4   8    101.45s   0.0053s    19100×   101.36s  (100%)
  4   9     >120s    0.0162s    >7400×        —
```

Tuple LLTs (`llt(k).cospin(tuple)` against `llt_g`):

```text
  ((2,2),(2,1),(2))     n=9     0.58s   0.0001s   ~5800×
  ((2,2),(2,2),(2))     n=10    1.45s   0.0002s   ~7250×
  ((3,2),(2,2),(1))     n=10    1.67s   0.0002s   ~8350×
  ((2,2),(2,2),(2,1))   n=11   10.36s   0.0010s  ~10400×
  ((3,2),(2,2),(2,1))   n=12   66.99s   0.0049s  ~13700×
```

**Where Sage is competitive, stated plainly.** `e[n].nabla()` is a real
implementation and the LLT route is not the way to beat it:

```text
  n    Sage e[n].nabla()   deltaop::nabla_e   llt::nabla_e_by_path
  9         4.64s               0.212s             0.956s
 10        14.64s               0.751s            19.278s
 11        41.10s               1.801s               —
```

`deltaop` is the route for the total (~20× over Sage — `macdonald-operators.md`'s
number, not this module's). At n = 10 the by-path route is *slower than
Sage's total* because it is not computing the total: it emits all 16 796
per-path Schur-positive pieces, an object Sage has no entry point for at any
speed. Quoting 19.3s against 14.6s as a loss, or against nothing as a win,
would both be wrong — it is a different deliverable at comparable cost.

**The one number that is not a speed claim.** Sage's single shape
`HSp3[[n]] → s` cost 0.25 / 1.06 / 6.76 / 47.96s at n = 6…9 and timed out at
10, while `llt_h(&(n), 3)` cost 10–70 microseconds across that whole range —
a ratio of 10⁵–10⁶ that says nothing about arithmetic throughput: the
one-row shape is 94–100% of a Sage degree and under 1% of the pruned walk's,
because chains to a single row pass only through single rows. **10–70
microseconds is timer resolution here, so the ratio records that the shape
left the walk's critical path, not how fast the arithmetic is.**

## Provenance: the same walls, on battery

SageMath 10.9, this machine, 2026-07-29 — every row on battery (~96%,
discharging), superseded by the re-run above and kept for provenance, not
arithmetic. ⚠️ The re-run measured ~2× faster across the board, and is the
third of the three drift measurements under "Power state" in
[README.md](README.md).

```text
  k \ n      7        8        9        10       11
   2       0.65s    2.22s   10.48s    51.16s    >120s
   3       2.11s   13.01s   90.01s    >120s        —
   4      18.72s   >120s       —         —        —
```

Worst shape μ=(n): 94% of the degree at k=2, n=7, rising to 98–100% by n=9–10
across k — the same cost-concentration finding the table above
restates in seconds instead of shares.

Tuple LLTs, the same five tuples as the table above, roughly 2× slower:
1.14s / 2.80s / 3.14s / 19.69s / >120s. Single coefficients: `(6,6,4,2)` wt
`1⁹` k=2 0.15s; `(8,8,6,4,2)` wt `1¹⁴` k=2 >120s; `(12,9,6,3)` wt `1¹⁰` k=3
0.86s.

## Specializations and degenerations, verified

```text
  k = 1:       G̃^(1)_λ = s_λ                    (1-ribbons = cells, spin 0; definitional)
  q = 1:       G_ν(x; 1) = ∏_i s_{ν^(i)}        (verified vs Sage skew Schurs)
  k ≥ ℓ-bound: H^(k)_μ = Q'_μ(x; q)              ([LLT] Thm 6.6; verified H^(4)_{3211} = Sage HL Qp)
  ω-duality:   G_{ν'}(x; q) = q^{A(ν)} ω G_ν(x; 1/q), A(ν) = #attacking pairs,
               ν' = (transpose each component, reverse order)  (verified, 4 tuples)
  coproduct:   G_ν(X+Y) = Σ G_inner(X) G_{outer/inner}(Y), the [BHMPS]
               inv-offset bookkeeping                      (verified on ((2),(1)))
  Macdonald:   H̃_μ(x; q, t) = Σ_D q^{−a(D)} t^{maj(D)} G_{ν(μ,D)}(x; q)
               ([HHL]; verified against Sage Ht for all μ ⊢ n ≤ 5)
  shuffle:     ∇e_n = Σ_D t^{area(D)} G_D(x; q)   (verified n ≤ 5 vs Sage nabla)
  chromatic:   X_Γ(x; q) = (q−1)^{−n} G_Γ[x(q−1); q]   ([CM] Prop 3.5 / [AP] Lemma 47;
               verified vs Sage chromatic QSF, all unit-interval graphs n ≤ 5)
  e-expansion: Ĝ_P(x; q+1) = Σ_θ q^{asc(θ)} e_{λ(θ)}, blocks by highest-reachable
               vertex                                   ([AS]; verified on the fixtures)
```

Nine identities, each checked against Sage or a cited theorem rather than
assumed; the shuffle and Macdonald rows are what tie this module to
`deltaop.rs` and `qtkostka.rs` respectively.

## Recorded dead ends

Four, each a wrong first guess before the verified law replaced it.

- **[HHL]'s remark-shaped `q^e G(1/q)` as the k-quotient dictionary.** [HHL]'s
  `G̃_λ` is spin-flavored, [LLT]'s `G̃` is cospin — coding the remark directly
  fails. The direct, min-inv-floored equality between the tuple and ribbon
  sides is the actual law. Cost one debugging round.
- **Recovering spin from a bare tuple.** Impossible, by fixture: λ = (1,1,1,1)
  at k=2 has one ribbon tableau with `s* = 1`, but its 2-quotient tuple has
  `max inv = min inv = 0` — no statistic of the tuple recovers `s*`. The
  spin-graded `H` is exposed on the partition-plus-level side only.
- **The cell-below-path graph as the dinv dictionary.** It is a real graph
  but the *other* one: `a = (0,0)` has no cells below the path yet its rows
  form a primary dinv pair, so this graph and the dinv-faithful decorated
  graph disagree. Both are kept, for different jobs.
- **[LLT] §7's printed straightening rule.** Contains two misprints; [KMS]
  (43)/(45) is normative, and is what the Fock route implements.

Three of the four are literature misprints or near-miss symbols rather than
implementation bugs: the formula was transcribed correctly and was still
wrong, so what caught each one was a comparison against something outside the
paper, not a failing test.

## Targets, guessed and measured

Stated before the code existed, so the measurement could contradict them —
house precedent by then: the `st` basis guessed 50× and got 3400×, the
Δ-operators guessed 100× and got ~21×, Jack guessed 200× and got 8340×.
Measured 2026-07-30 (`cargo run --release --example bench_llt --
14`, single-threaded, `i64`):

| target | guessed | ours | Sage (mains) | verdict |
|---|---|---|---|---|
| `H^(3)` table, n=9 | ≥100× (≤0.9s) | 0.008s | 47.34s | **5900×** |
| `H^(2)` table, n=10 | — | 0.014s | 24.16s | **1700×** |
| `H^(4)` table, n=8 | — | 0.0053s | 101.45s | **19 100×** |
| k=2,3 tables through n=14 | "minutes each" | 0.25 / 1.08s | dies at n=11/10 | met, ~250× |
| tuple `((3,2),(2,2),(2,1))`, n=12 | <1s | 0.006s | 66.99s | **13 700×** |
| `nabla_e_by_path(10)`, 16 796 pieces | "minutes" | 19.3s + 3.7s | no entry point | met |
| KL columns, all λ⊢6, k=2 | seconds *per column* | 0.041s for all 11 | no entry point | met |

The two "no entry point" rows are the honest half of the table: not
speedups, but objects Sage cannot produce at any speed.

## What the profiler found

`examples/profile_llt.rs`, one attributable workload per route, `sample` on
macOS. The first working version measured 0.019 / 0.053 / 0.115s at n=10 for
k=2/3/4 and 80.4s for `nabla_e_by_path(10)`; the numbers above are *after*
this pass:

| route | first profile | fixed | gain |
|---|---|---|---|
| R1 (SYT walk) | 85% in the walk, 6% hashing | bit-mask poset (`pred_mask`, `attack_mask`), direct-indexed descent table | **4.1×** |
| R2 (abacus strips) | ~35% `malloc`/`free`; `strip_rec` 70% of the rest | scratch buffers, all-weights-in-one-walk, O(changed bits) containment, runner-mask block scan | **2.2–2.5×** |
| R3 (Fock straightening) | 63% `malloc`/`free`, `straighten` itself 6% | in-place wedge mutation, generated offset ladder, one-pass coefficient helpers, probe-before-insert | **3.0×** |

Three things generalize past this module:

- **In two of the three routes, the first profile pointed at the allocator
  rather than at the mathematics.** R3 spent two-thirds of its time there
  because the straightening recursion cloned its wedge per branch; the fix
  mutates two positions and restores them rather than cloning. R2 was ~35%
  `malloc`/`free`. The same shape of finding recurs at `jack.rs` and
  `macop.rs`.
- **The naive form of each hot loop stayed on as a test oracle** where one
  existed (`strip_blocks_agree_with_naive_subsets`,
  `abacus_containment_is_partition_containment`,
  `the_all_weights_walk_agrees_with_the_fixed_weight_one`) — the fast form
  carries the reason in a comment, the slow form keeps checking it, and the
  last of those tests exists specifically because the optimization moved
  production code off the path the first test covered.
- **R1 and R2 are now compute-bound; R3 is still ~38% allocator**, from the
  `QtPoly` temporary each straightening branch creates. Left there: the next
  step would be threading a coefficient stack through the recursion, and the
  route is not on the critical path of any sweep.

**The 6% that came back on the fallback.** The direct-indexed descent table
above removes R1's per-leaf hash only for tuples that fit
`FLAT_TABLE_BUDGET`; past it `MapSink` reinstates one hash per standard
filling, and it was still SipHash on a bare `u64`. Switching it to
`fasthash::Map` is **1.51×**, min-of-3:

| tuple | n | A(ν) | SipHash | `MixHasher` |
|---|---|---|---|---|
| `((2,2),(2,2),(2,2),(2,2))` | 16 | 60 | 19.40s | 12.81s |
| `((3,3),(3,3),(3,3))` | 18 | 54 | 44.89s | 29.81s |

Larger than the 1.12× the same one-line change bought
[the character sweep](transitions.md), and for the same reason the layer
maps that key on `Vec<u32>` got almost nothing: what the change is worth
tracks whether the key is *already a word*. `MapSink`'s is; a `Partition`
key's is not, and the
`to_vec()` behind it costs more than the hasher either way. Measured
non-results, same harness pattern: `kostka_uncached`'s layer 1.08×,
`strip_lr`'s state map and `convert::jt_terms`' accumulator both nil.

A trap the A/B turned up on the way: the two sinks **disagree on row
length** and always have — `FlatSink` pads every row to `A(ν) + 1`,
`MapSink` grows a row only to the largest `inv` it saw. Both consumers skip
zero counts, so this was invisible, but the fallback had no test at all
(every tuple that reaches it costs tens of seconds).
`both_syt_bucket_sinks_agree_up_to_trailing_zeros` now drives the budget to
zero to run small tuples down the map path, and pins the trailing-zero
latitude rather than hiding it.

**The algorithm not taken.** R1's walk costs `#SYT(ν)` per path, which is now
essentially all of it. A subset DP over (assigned set, last cell) would cost
`2^n · n²` per content instead — better once `#SYT` passes
`p(n) · 2^n · n²`, around n=12–13 for these tuples, so *worse* at the n ≤ 10
the shuffle refinement wants. Recorded rather than built.

Peak RSS for the whole bench is 220 MB, dominated by `nabla_e_by_path(10)`
holding all 16 796 pieces at once rather than by any ribbon walk. The
per-phase split is instrumented for the by-path table only; the rest is
recorded as owed.

## The rise ladder, factored through the LLT engine

`dyck::ladder`'s rise side is a function of the area sequence alone — no
labels — so the whole labeled-path walk factors into one LLT evaluation per
area sequence plus a knapsack (`dyck::rise_ladder_via_llt`, dispatched from
`Side::Rise`). The labeled walk survives as `ladder_at_content`: it is the
oracle, and it stays cheaper for one coarse content, where μ=(n) is a single
labeling.

Measured (`examples/probe_llt_ladder.rs`), mains:

```text
  n           5      6      7       8       9
  llt (s)   .0006  .0037  .0174   .0987   1.310
  walk (s)  .0013  .0171  .1350   2.873  73.699
  speedup    2.1×   4.7×   7.7×   29.1×   56.3×
```

~2× per degree, so n=10 extrapolates to ~110× and the walk's ~25 minutes to
~13s. The same table on battery gave 26.1×/56.6× at n=8,9 with the absolutes
roughly doubled — an independent confirmation of the mains drift, since the
*ratio* is power-independent even though the seconds are not. It also yields
a per-path Schur-positive refinement of `Δ'_{e_k}e_n` for every k, where the
by-path decomposition alone had only reached the k=n−1 top.

⚠️ **A number first reported in this session was wrong, and the mistake is
worth keeping.** The first probe measured 140× at n=8 by looping
`dyck::side(n, k, Rise)` over k — but `side` computes the *whole* ladder and
discards all but one slot, so looping it charged the labeled walk n times
over. The honest figure is 29×. **`side`'s own rustdoc already said that
asking for one slot costs the same as asking for all of them, and the A/B was
built by looping it anyway.**

The valley side does not factor — `Val` reads the labels, and its weights
are per-labeling — so it keeps the full enumeration and stays the open,
expensive half.

**A fourth H̃ engine was tried and is not competitive.** `htilde_by_llt`
beat `bh`'s route by ~1.3× through n=7, tied at n=8, and lost 3× at n=9,
growing — the `2^{|μ|−μ₁}` descent subsets it competes against eventually
win, and there is no pruning left to add. It keeps its value as an
independent, positively-graded cross-check on the H̃ assembly, and nothing
more.

## Two open questions, answered

**[LLT] Conjecture 6.4** — is `H^(k+1)_μ − H^(k)_μ` Schur-positive? — was
open on every reading available when this module was specified. Swept
2026-07-30 (`bench_llt`'s Conj 6.4 section): k=1…4, every μ⊢n for n ≤ 14, 56
(k,n) rows, **120 943 nonzero Schur coefficients of the difference, none
negative**, the whole sweep costing 34.9s. The conjecture is still
open — a sweep is evidence, not a proof — and the bench prints
`*** COUNTEREXAMPLE -- REPORT ***` rather than failing an assertion, the same
report-not-assert rule `macdonald-operators.md` records. At 34.9s for n ≤ 14
the live question is what range would be *informative*, not what range is
affordable.

**Does R3 (Fock straightening) beat R2 (abacus strips) for Schur output at
scale?** Partly answered, and the answer is no in the measured range. All
λ⊢6 at k=2 costs R3 0.041s of straightening; the whole-degree R2 table it
would have to beat (n=6, k=2) costs 0.0007s, and R2's own range is n=14 at
0.25s. R3 also grows faster — about 12× per degree at k=3: 0.004 → 0.051 →
0.65s for λ⊢4,5,6 — because the wedge straightening branches where R2's
weight-trie shares. R2 is therefore the engine for whole-degree sweeps; R3
is kept for what it produces rather than what it costs — its columns
are Kazhdan–Lusztig polynomials, an output R2 cannot produce at all.

## Offline oracle fixture

The three ribbon dictionaries, committed and checked with no Sage: 36 `H^(k)`
and 36 `H̃^(k)` expansions for k = 1..3 through degree 4, and 120 `G̃^(k)`
expansions wherever k divides |λ|.

Two conventions are pinned. All three dictionaries are carried rather than
derived from one another, because four normalizations of `G̃` circulate and
agree on the easy cases. And the fixture writes the grading exponent into
symfn's **q** slot, where Sage names the same parameter `t` — the translation a
generator written without thinking about it gets backwards. `k = 1` is in the
sweep because `H^(1)` is the Schur function, the cheapest place a spin/cospin
swap shows.

## Next

- **Is the fundamental expansion really unshipped elsewhere?** `llt_fundamental`
  and its Python binding both claimed nobody ships it, and neither pointed at a
  survey; grepping `docs/record/` and `research-gaps.md` for any backing found
  none. Both claims were deleted in the 2026-08-07 audit rather than left
  unbacked. If a survey confirms it, the claim is worth making again — properly
  this time, with a dated incumbent comparison behind it.
- **Nonzero-core quotients.** The k-quotient dictionary is verified for
  empty-core λ with offset 0. What offset vector makes the tuple model match
  `G̃` for a general k-core, and is the min-inv floor still the only
  normalization needed? The literature's stated conventions were already wrong
  against the fixtures once here; measure this rather than trust it.
- **[AP] Conj 25, [AS] Problem 6.20, and Shareshian–Wachs e-positivity** —
  the unicellular open problems the chromatic bridge opens onto. The sweep
  infrastructure belongs to the `chromatic-corpus` branch, not here; this
  module's bridge is verified and waiting.
- **The [BHMPS] Catalanimal route** — `∇` of a general LLT by
  raising-operator combinatorics, no Macdonald-basis pass. The v2 route past
  what `deltaop.rs` plus a basis conversion can do; not specced further.
- **Does R3 win at a large k with a small λ?** Open — that is where R2's
  abacus has many runners and few beads per runner, the shape least
  favorable to the weight-trie sharing that wins everywhere measured so far.
- **Positivity bookkeeping for [GH]-only cases.** For skew (not straight)
  tuples, Schur-positivity rests on an unpublished 2007 preprint. The engine
  should flag a negative coefficient on a skew-tuple expansion as a
  first-order finding — it would be either our bug or new mathematics.

## Negative result: both LLT routes are already enumeration-bound

Checked in the same sweep as [jack.md](jack.md) — looking for hot maps keyed on
a freshly allocated `Partition`, the defect worth 2.0–4.0x in
[transitions.md](transitions.md). `llt.rs` has two such maps, so it looked like
a candidate.

It is not, in either route. Sampling `profile_llt r1 9`: `SkewTuple::walk_rec`
**81.3%** and `for_each_area::rec` 11.5% — 93% in the combinatorial walk, with
the allocator under 3%. Sampling `profile_llt r2 8`: `strip_any_rec` **52.8%**,
`collect_blocks` 7.2%, `QtPoly::add_shifted` 7.0%, allocator ~14%.

Both profiles are the shape a healthy engine has: the enumeration dominates and
the data structures are noise. Nothing here resembles the 53.8%-allocator,
11.3%-mathematics profile the Pieri layer had before its rewrite. Recorded so
the two `HashMap<Partition, _>` sites are not "fixed" on the strength of
grepping for them.

## The skew tuples reach the surface (2026-08-25)

Stage 4 of
[convenience-surface-review.md](../plans/convenience-surface-review.md), the
optional item. The kernel has always computed on the \[HHL\] Def 3.2 object —
`SkewTuple::from_skews` predates this change — but the boundary took straight
shapes only, so `G_ν` could not be asked for the object its mathematics is
defined on. `llt_g`, `llt_min_inv` and `llt_fundamental` now accept each
component as either a plain shape or an `(outer, inner)` pair (permissive-in,
P1; a pair is two sequences, which a list of parts never reads as, so no tag
is needed), and `llt.G` and `llt.min_inv` take the same spellings. An inner
not contained in its outer is refused at the entry point with a `ValueError`
naming both shapes, rather than reaching `from_skews`'s panic (R11).

Pinned by ten `lltgskew` fixture rows against Sage's `cospin` on skew
partitions (`llt_skew_tuples_match_sage` in `tests/sage_oracle.rs`,
regenerated with `SAGE_DISABLE_SYMFN=1`): Sage divides out the floor
`q^{min inv}`, so the comparison multiplies it back through `llt_min_inv` on
the same tuples, and two of the rows are content translations of each other —
`(2,1)/(1) ∪ (2)/(1)` against `(2,1)/(1) ∪ (1)` — which the generator emits
separately and Sage values identically, the invariance the offset model
predicts. The convenience sweep holds the spellings to each other: empty
inner is the straight shape, pairs and plain shapes mix, translation moves
nothing (`check_llt_skew_tuples`).
