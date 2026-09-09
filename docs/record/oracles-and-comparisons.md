# Oracles and comparison harnesses

Sage, lrcalc and Symmetrica are all used as oracles and as baselines, and
the three differ in what they can prove. This file holds the harnesses
themselves and the findings that came out of building them — including
two whole-table deficits that only a like-for-like C comparison could
show.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

## Against Sage (`scripts/compare_sage.py`)

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

### ⚠️ Half these rows are not benchmarks against Python

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

## Against Symmetrica directly (`scripts/compare_symmetrica.py`)

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
`outerproduct_schur` degrades steeply over this range (0.0026s → 530s while
symfn goes 0.0008s → 0.059s), so `run_big_lr` stops as soon as it passes the
budget. At toy sizes this same comparison reads 3.2x, so sizing the benchmark
where the work actually lives changed the answer by three orders of magnitude.

**Two new deficits, both on whole tables.** Our *per-value*
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
  layer, discarding everything else that layer holds. Drop the bound and
  the layer at the end of μ's chain **is** the whole column — every λ with its
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

### The plethysm comparison covers only the inputs Symmetrica accepts

Symmetrica's plethysm **refuses a multi-row outer partition** — it reports "for
the moment only for outer S_n" and computes nothing. Every case in the table is
therefore a single-row outer, the easy case, and possibly a specialized path.

symfn has no such restriction. `s_{21}[s_{21}]` (17 terms, 0.0011s),
`s_{22}[s_2]`, `s_{32}[s_{11}]`, `s_{21}[s_{31}]` (39 terms) and
`s_{31}[s_{22}]` (104 terms) all compute and all agree with Sage — 48/48 cases
across both shapes of input. So the honest summary is that symfn is faster than
Symmetrica on the inputs Symmetrica accepts, and is the only one of the two that
handles the rest.

### s → p: fixed, but not the way expected

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
against. That test caught something immediately: asserting on strip *order*
failed, because the masked form walks β upward and the general one downward —
invisible to callers, which only sum, but noticed rather than assumed harmless.

What was still open here — `character_cached` cloning both partitions into its
key and taking a global `RwLock` at every node, `character_uncached`
allocating a fresh `Partition` for the μ-suffix at every node — is closed in
[transitions.md](transitions.md), "The character recursion on β-masks": the
recursion now runs on β-masks end to end with a word-keyed memo, 4-7x on a
character and 1.6-4.6x on `s → p`.

⚠️ The plethysm row is capped at degree 10 regardless of the ladder setting
(its input is the *outer* partition and the result reaches degree 30), so that
row does not vary across the columns above — its four entries are the same
measurement repeated.

---

## The Sage job in CI (2026-09-04)

Every live harness in this file, and every `check_*.py` that imports Sage,
ran only when someone had Sage and remembered — the state
[../policies/validation.md](../policies/validation.md) V1 describes. From
2026-09-04 they also run unattended in `.github/workflows/sage.yml`, a
workflow separate from `ci.yml` so that nothing in it can gate a pull request.
It runs on pushes to main, on tags, on a weekly schedule, and by hand.

**What it installs.** Sage and lrcalc from conda-forge, through
`mamba-org/setup-micromamba`, unpinned. conda-forge's Sage was 10.9 on the day
this was written (published 2026-05-29; the fixtures were generated by
10.10.beta4). Unpinned on purpose: the weekly schedule exists to meet the next
Sage, and a value Sage changes between releases arrives as a fixture that no
longer regenerates. The conda-forge `lrcalc` package ships the `lrcalc`
executable `gen_lrcalc_oracle.py` drives, checked on the local sage-dev
environment, whose binary comes from that package. The scripts' usage lines
say `sage -python`; the workflow says `python`, because in a conda
environment that is Sage's interpreter, and because the development checkout
of Sage here has no `-python` at all — the first timing run below failed
fourteen times on that before it ran once.

**What it runs, in order.** `scripts/sage_guard.py` as a program, which
prints which Sage answered. Then both fixture generators, with `diff` against
the committed file and a nonzero exit on any difference: a Sage that changed
an answer, a fixture edited by hand, and a generator edited without
regenerating all land there. Then the seven dump-then-check harnesses at the
sizes their usage lines name, then the six that import the module into
Sage's interpreter, after `pip install .` into that interpreter. Every step
after the first fixture carries `if: !cancelled()`, so one red family does not
hide the others.

**What it costs.** Timed once on this machine, on battery, with the local
backend disabled and the dumps prebuilt, by a shell loop around each command
(`time_sage_checks.sh`, not committed). Sage startup is about 3 s of every
figure. Wall times are what a runner schedule needs, not a benchmark, so
they are rounded:

| step | local wall |
|---|---|
| `gen_sage_oracle.sage` | 7 s |
| `gen_lrcalc_oracle.py` | under 1 s |
| the six dumps together | about 1 s |
| `check_jack`, `check_hl`, `check_llt`, `check_macdonald`, `check_kf`, `check_eval`, `check_skew`, `check_qt_kostka`, `check_schubert_bindings` | 3–5 s each |
| `check_hl_p` | 6 s |
| `check_bindings 7` | 18 s |
| `check_st 6 4` | 44 s |
| `check_deltaop 7` | 52 s |

About two and a half minutes of scripts. The Sage install dominates the job
and is cached across runs by the action.

Both regenerated fixtures were byte-identical to the committed files, so the
workflow uses a plain `diff` rather than a sorted one; the generators'
output order is stable across the two Sage builds that have produced it.

**What stays out, and why.**

- `check_backend.py` drives the adapter through Sage's own dispatch and asserts
  which side answered. The adapter lives on the Sage branch
  ([../sage-backend.md](../sage-backend.md)) and in no Sage conda-forge ships,
  so on the runner it has nothing to drive. It is tested where it lives, by
  Sage's doctests on that branch. This is why the Phase 5b item that names
  the CI job stays open in
  [../release-readiness.md](../release-readiness.md).
- The `bench_*` and `compare_*` scripts are timings, and a shared runner's
  clock is not a measurement.
- `verify_deltaop_formulas.py` and the `spec_*.py` trio verify formulas on
  Sage alone, before any Rust existed; `check_deltaop.py` holds the same
  identities against the implementation.

**Open.** Two things this job has not yet shown.

- It has not run. Written and validated as YAML, like every workflow in this
  tree was before its first push; the first run will say whether conda-forge's
  Sage 10.9 regenerates a fixture that 10.10.beta4 produced.
- `SAGE_DISABLE_SYMFN=1` is set for the whole job against a Sage that does
  not yet need it. The day conda-forge's Sage carries the backend, the
  variable makes Sage answer out of Symmetrica and its own Python, and the
  first step's line changes from "stock, no symfn backend module" to "symfn
  backend present and disabled". That is the intended state and the job keeps
  working. What the variable cannot survive is a Sage that drops its own
  implementations rather than keeping them behind the switch; then Sage
  stops being an oracle for those families, and the job's Sage must be pinned
  to the last release that still has them. The fixtures stay valid evidence
  either way, because they record what an independent implementation
  answered.

---

## Laws over inputs nobody chose (`tests/random_laws.rs`, 2026-09-05)

`tests/algebra_laws.rs` sweeps every partition of degree at most 5 and stops
its product laws at degree 3, so a defect needing a shape outside that list
is never asked about. `tests/random_laws.rs` asks the same laws — the three
Schur round trips, the power-sum round trip over ℚ, `to_schur` as a ring
homomorphism on `h` and `e`, ω as an involutive algebra map, the three Hall
pairings, Δ as an algebra map — about generated partitions instead: round
trips to degree 10, pairings to degree 10, products with factors to degree 7
each, the coproduct to degree 4 each, a thousand cases per law by default.

No `proptest`. The default build has no dependencies and `cargo test` on the
published tarball runs offline ([../release-readiness.md](../release-readiness.md),
Phase 4), so the generator is a xorshift64* seeded by a constant mixed with
the law's name, the pattern `examples/bench_schubert_wall.rs` already used
for permutations. What that gives up is shrinking; every assertion names the
input, and `SYMFN_LAWS_SEED` and `SYMFN_LAWS_CASES` rerun any seed at any
size. The partitions are drawn part by part from the largest down, which is
not uniform over partitions of `n` and is not meant to be: it reaches long
thin shapes and short wide ones in one run.

Both sides of each comparison share `i128`, the shape
[../policies/failure.md](../policies/failure.md) R10 names, and the escape
clause is the a-priori bound `widest_value_stays_far_below_the_width`
measures on every run: the largest coefficient any law produces is held
under `i64::MAX`, leaving 64 bits between the values and the width.

First run: eight seeds at a thousand cases per law, no failure, a quarter of
a second per seed in the debug build. The suite is in `cargo test` and so in
every CI lane.

## symfn against its own history (`examples/bench_suite.rs`, 2026-09-05)

The `bench_*` examples each measure one subsystem against an external
baseline, as deep as that comparison needs; none of them says whether a
commit made the library slower than it was. `bench_suite` is the shallow and
wide instrument for that: every workload in `symfn::measure::workloads` —
already spanning the subsystems, already sized where the work lives for the
memory budgets — timed with the caches cleared before each repetition and
the minimum of three kept, printed as tab-separated rows.
`scripts/bench_compare.py` diffs two such files and marks the ratios outside
±20%, the band the LR record draws for one run on this machine; rows under
5 ms are printed and never flagged, since at that scale the timer and the
scheduler are what is compared, and `schubert` in the catalog is such a row
by design (one term, sized for its memory budget).

Not a test, and not made one: wall time is not assertable
([memory.md](memory.md), "Why memory can be a test when time cannot"), which
is why the catalog carries memory budgets and this carries none. The
committed run is [bench_suite.tsv](bench_suite.tsv), whose header names the
machine, the date, the power state and the commit; a comparison against it
means something on that machine and nowhere else. Two runs on the day it
was committed, on battery with low power mode off, agreed within 13% on
every row above the floor and within 6% on all but one.

## Coverage: what the test suite never reaches (2026-09-05)

`cargo llvm-cov --features bignum` over the test suite, in the debug
profile, read once for the question Phase 7 asked — which paths no check
reaches — and run on every push by the `coverage` job in `ci.yml`, which
uploads the report and is never red: a percentage is not a claim this tree
makes, and a threshold would either pass today and check nothing or fail on
the next workload added. Lines executed, 95%; the number is here for the
date, not as a bar. Two caveats shape what the report can say. Doctests are
not instrumented, so a function whose only check is its doctest reads as
uncovered although `cargo test` runs it. And `tests/memory.rs` asserts
nothing under `debug_assertions`, so `measure/` reads as unexecuted in a
debug run and is not.

What the report found, after those two are discounted:

- **Two root re-exports no test calls.** `two_row_coeff`, the single
  coefficient of the two-row strategy, has no caller in the suite, no
  doctest and no example; the product form `two_row_product` is checked, the
  coefficient form is not. `modified_qt_kostka` likewise.
- **Two functions checked only by an example.** `qt_kostka_table_via_operator`
  agrees with the other two routes wherever `examples/bench_qtk_routes.rs`
  asserts it, and nowhere in `cargo test`; this file's "a benchmark asserts
  all three agree" is true and is not a test. `gj::double_coset_coefficient`
  is reached only from `examples/gj_tables.rs`.
- **Size-gated paths the suite never crosses.** `convert.rs`'s prefix-group
  walk `expand_shared` and the `JT_LIMIT` branch of the Schur sweep, taken
  for shapes with more than 14 rows and columns, run only on inputs the
  oracle scripts and the LR examples reach; every fixture is below the gate.
- **The `bignum` feature's `BigRational` impls** of `Integral` and
  `Plethystic`, and `Guarded`'s `gcd`, are unexercised although the feature
  was on: `tests/bignum.rs` drives `BigInt`, not the rational.
- **Branches with one side never taken:** `deltaop::Ratio` equality for two
  different denominators (the lcm branch), `gjmod`'s degree-0 table and its
  `g > 1` reduction, and the `Display` impls for `AFrac`, `QtPoly` and
  `Perm`, which no test formats.

None is a defect found; each is a place where the claim "every value a
public family can produce is covered by a check that does not share its
mathematics" ([../policies/validation.md](../policies/validation.md)) rests
on a doctest, an example, or nothing. The first two bullets are the ones to
close, and they are cheap: an agreement test between `two_row_coeff` and
`okada_coeff` or `SkewLr::lr_coeff` on two-row shapes, and the three-route
agreement moved from the benchmark into a test at one small degree. Open
until then.

---

## Spot rows above the sweeps (2026-09-05)

Every sweep in `scripts/gen_sage_oracle.sage` stopped at degree 5 or 6, so
above that the families were checked only by agreement between this crate's
own engines, which [../policies/validation.md](../policies/validation.md) V5
accepts and its own "does not share its mathematics" clause does not, for
the range the rustdoc advertises. The generator now ends with a block of
spot rows — single shapes and pairs, not degree sweeps — at degrees 8 to 15,
each family asked about a long shape, a wide one and a balanced one, and
about a zero where it has zeros:

| family | rows | degrees |
|---|---|---|
| Kostka numbers, characters, Kostka–Foulkes | 11 pairs each, two of them zeros | 10, 12, 15 |
| Hall–Littlewood `Q'` and `P` | 7 shapes each | 10, 12 |
| `(q,t)`-Kostka | 6 pairs | 8, 9, 10 |
| `H̃` in the Schur basis | 6 shapes | 8, 10, 12 |
| Macdonald `P` | 5 shapes | 8, 9, 10 |
| Macdonald `J` | 5 shapes | 8, 10, 12 |

The tests read the fixture row by row and assume no degree is complete, so
the rows needed no parser change; every one passed on the first run, and
the regenerated fixture differs from the committed one by exactly those 69
added lines.

**What the floors cost, and why they stop where they do.** Measured on Sage
10.10.beta4 with `SAGE_DISABLE_SYMFN=1`, on battery. Kostka numbers,
characters and Kostka–Foulkes are instant at degree 15; Hall–Littlewood is
under 2 s a shape at degree 12; `H̃` is under 0.3 s a shape at 12. Macdonald
`P` and `J` cost by the degree rather than the shape: Sage builds the
degree's transition matrix on the first shape asked for and answers the rest
from it — 2 s for degree 8, 7 s for 9, 18 s for 10 and 185 s for 12. The
generator went from 7 s to 3 min 50 s, and the degree-12 `J` rows are that
difference. Kept, because they are the only Sage answers above degree 10 the
fixture holds for the family; the Sage workflow regenerates the file weekly
and pays it there.

On the Rust side the debug-build suite went from 0.15 s to 11 s, in two
tests that run in parallel: the Macdonald expansions and `H̃`, each about
10 s, which is symfn computing the degree-12 rows unoptimized. `cargo test`
still runs in seconds.

The after-the-tag item this closes named degree 10 to 15 for characters,
Kostka numbers, Hall–Littlewood and Macdonald; the `(q,t)`-Kostka rows and
the zeros are additions, the first because the table had no external row
above degree 5 and the second because a family that is right on its support
and wrong about the support passes a nonzero-only fixture.

## The incumbent's walls, before any of this was built (2026-07-28)

The survey that the capability claims in `src/lib.rs` and the family module
docs rest on. It asked what a symmetric-function package cannot do at all,
and measured SageMath 10.9 on this machine to find out — single runs,
per-item `SIGALRM` timeouts of 90–120 s, power state not recorded. It was
compiled before the crate had Jack, the Macdonald operators, LLT or the
Orellana–Zabrocki bases, so every row below is the incumbent alone.

| operation | input | result |
|---|---|---|
| `internal_product` (Kronecker) | `s[10,7,5]·s[9,7,6]`, n=22 | 2.92s → 721 terms |
| | `s[14,10,8]·s[13,10,9]`, n=32 | **>90s timeout** |
| | `s[12,9,7,4]·s[11,10,7,4]`, n=32 | **>90s timeout** |
| | `s[20,15,10]·s[18,15,12]`, n=45 | **>90s timeout** |
| single coefficient `g(λ,μ,ν)` | 3-row, n=45 | **>90s timeout** |
| plethysm | `h₅[h₅]`, deg 25 | 2.13s → 245 |
| | `h₄[h₇]`, deg 28 | 7.66s → 173 |
| | `h₅[h₆]`, deg 30 | 12.03s → 492 |
| | `h₆[h₆]`, deg 36 | 110.87s → 2002 |
| | `h₄[h₁₀]`, deg 40 | **>120s timeout** |
| chromatic symmetric function | path `P₁₆` | 0.34s |
| | path `P₂₀` | 3.48s |
| | path `P₂₄` | 58.48s |
| | random tree, 26 vertices | **>120s timeout** |
| Jack `P → s` | n=9 | 5.59s |
| | n=15 | **>90s timeout** |
| Macdonald `H̃ → s` | n=15 | 2.10s |
| | n=18 | 10.55–15.41s |
| | n=21 | 66.82s |
| ∇ (nabla) | `∇e₆` | 0.21s |
| | `∇e₈` | 2.13s |
| stable Kronecker via `st` basis | `st[3,1]·st[3,1]` | 0.31s |
| | `st[6,4]·st[6,4]` | **>90s timeout** |
| | `st[8,5]·st[7,4]` | **>90s timeout** |
| e-/Schur-positivity certification | — | **no API exists** |

⚠️ Single runs on a laptop with coarse timeouts. These are order-of-magnitude
walls, not benchmarks — the caution this file's first table carries applies
here too.

**The `st` row is the sharpest.** The Orellana–Zabrocki
irreducible-character basis is the modern tool for reduced Kronecker
coefficients — its outer-product structure constants are the stable Kronecker
coefficients — and Sage's implementation dies on two two-row partitions of 10.
What came of it is in [kronecker.md](kronecker.md), "The Orellana–Zabrocki
character bases".

**No single-coefficient path existed anywhere.** Every package computes the
whole product to read one number, which is the `lr_coeff` defect
[littlewood-richardson.md](littlewood-richardson.md) records ("peeling off the
*larger* factor instead", 23x, a defect no product benchmark surfaces). It was
unexercised for Kronecker and plethysm then;
[kronecker.md](kronecker.md), "A single coefficient, without the product",
is what it became.

**What does exist, to be fair.** Sage has `internal_product`,
`reduced_kronecker_product`, `plethysm`, Macdonald (`P Q J H Ht S`), Jack, LLT,
`nabla`, `theta_qt`, the Orellana–Zabrocki `st` basis,
`chromatic_symmetric_function` and `chromatic_quasisymmetric_function` — broad
coverage, with speed as the gap and only the modern operators absent outright.
`lrcalc` is fast and narrow: LR and quantum LR. `barvikron` implements the
Christandl–Doran–Walter lattice-point algorithm for Kronecker coefficients of
bounded height, in polynomial time, as an unmaintained Python prototype;
Baldoni–Vergne–Walter distribute Maple code for the bounded-length case.
Stembridge's `SF` computes Kronecker naively and takes very general
user-defined bases.

Four of these rows were re-measured later, sharper and with the methodology
fixed, in the files that own them: Jack in [jack.md](jack.md) — where Sage's
whole-degree tables die at n = 12, not 15, and where the note is recorded that
these walls cannot be reproduced in a single Sage process at all, because it
memoizes the transition matrices and `SIGALRM` corrupts them mid-build — ∇ and
`H̃` in [macdonald-operators.md](macdonald-operators.md), Kronecker in
[kronecker.md](kronecker.md), and plethysm in
[plethysm.md](plethysm.md).
