# Plethysm

Plethysm has **two routes**. The general one goes through the power-sum basis,
where `p_n[g]` is part scaling and plethysm is multiplicative in `f`; a
one-row inner argument takes a Schur-basis recursion instead and never
converts at all. Most of what follows is the first route being made faster and
then having its wall found; the last sections are the second route, which
arrived because an outside tester asked why `s[6](s[6])` finishes nowhere.

The early sections are also about which baseline was being measured against —
the `py` rows the first one opens by criticizing are the Sage ladder in
[oracles-and-comparisons.md](oracles-and-comparisons.md).

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

## The plethysm row is measuring the wrong baseline

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

**Ahead on every case**, having been behind on every case. Against Sage,
plethysm went 9x -> ~40x.

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
   continues one layer per common prefix.
3. **Accumulation is keyed on the β-mask, not the partition.** Every leaf
   touches the whole layer, so a partition key allocated, sorted and hashed a
   fresh `Vec` once per (μ, mask) pair — ~9,000 allocations to produce 63 terms.

A fourth, worth its two lines: **`p_step` pre-sizes its output map.** The
layer grows monotonically through a sweep, so a default-capacity map rehashed
several times per step. Worth 1.22x, measured over 8 interleaved rounds.

Net ~3–4x, and it moved plethysm from 9x to ~25x against Sage. **We are now
ahead on three of five cases, at parity on a fourth, and behind only on
`s_4[s_{22}]` (0.50x)** — the scaling issue is reduced, not eliminated.

## ⚠️ Retracted: "the rational leaf arithmetic is not the bottleneck"

This section previously recorded a *negative* result — that killing the rational
arithmetic at each (μ, mask) leaf had only a 24% ceiling, and so was not worth a
`Ring` hook, an lcm with overflow guards, and a fallback path. **That conclusion
was wrong, and it was wrong because the measurement was a single un-interleaved
run.**

A sampling profile (2,496 samples) said otherwise:

| symbol | self | share |
|---|---|---|
| `p_step` (layer DP) | 852 | 34% |
| `u128_div_rem` (gcd) | 733 | 29% |
| `Rational::add_assign` | 366 | 15% |
| `Rational::mul` | 146 | 6% |
| `__modti3` / `__divti3` | 139 | 6% |

Rational arithmetic was **~55%**, not 24%. Re-running the same ceiling
experiment over 6 interleaved rounds gave 0.000875s → 0.000434s: a **2.02x**
ceiling. The original number came from one run of each side, on a machine that
has repeatedly been shown to drift 2x as it warms — the exact failure this
document had already warned about two paragraphs earlier, committed anyway.

The experiment was cheap and its result was used to **cancel** work, and it got
one un-interleaved run — less scrutiny than the same number would have needed
to justify building something. A profiler would have settled it in one shot,
for less effort than the experiment itself cost.

## The fix that followed: a common-denominator integer sweep

`p → s` computes Σ_μ c_μ χ^λ(μ) — a sum of (coefficient × integer) terms. Over
ℚ that is a rational multiply and a rational add per leaf, each normalizing by a
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

## The degree-32 cliff, and what was behind it

The first outside tester asked why `s[6](s[6])` does not finish anywhere,
while `s[5](s[5])` takes under a second. The answer was not the size of the
answer. It was a wall at degree 32 that nothing in the tree named as one.

β values run from 0 to under 2l, so a `u64` β-mask holds a degree-l sweep only
for l ≤ 32. Past that `to_schur` fell back to `character_in` **per (λ, μ)
pair** — p(n) character recursions in the coefficient ring for every term of
the p-element, against one shared sweep below the wall. At degree 36 that is
17,977 partitions against a p-element with up to 17,977 terms, in rational
arithmetic. Measured on this laptop:

| | degree | time |
|---|---|---|
| `s5[s5]` | 25 | 0.174s |
| `s4[s7]` | 28 | 0.919s |
| `s5[s6]` | 30 | 1.769s |
| `s4[s8]` | 32 | 5.717s |
| `s3[s11]` | 33 | **did not finish in 9 minutes** |
| `s6[s6]` | 36 | **did not finish in 25 minutes**, 2.3 GB resident |

The growth up to 32 is a factor of ~2 per two degrees. The step from 32 to 33
is not a step on that curve at all; it is a different algorithm.

**The sweep is now generic over the mask width**, `u64` through degree 32 and
`u128` above it, and the same laptop:

| | degree | before | after |
|---|---|---|---|
| `s3[s11]` | 33 | did not finish | **88.1s** |
| `s6[s6]` | 36 | did not finish | **28.1s** |

Degree 33 costing more than degree 36 is not a mistake. The sweep's cost is
driven by how many p_μ the batch carries and how much prefix they share, not
by the degree alone: `s3[s11]` expands an inner s_11 whose p-expansion has 56
terms, and the products of those fill degree 33 far more densely than
`s6[s6]` fills 36.

**The ceiling is the accumulator, not the mask.** A `u128` mask would hold
l ≤ 64, but the layer accumulates in `i128` and every value in it is a
character, so |χ^λ(μ)| ≤ √(l!), with one slot taking at most l contributions
before it settles: l·√(l!) is 6.2·10³⁷ at l = 55 and 1.5·10³⁹ at l = 56,
against an `i128` ceiling of 1.7·10³⁸. Hence `WIDE_MASK_LIMIT = 55`, which
covers `s7[s7]` at 49. Past it the character fallback still stands.

**The width is chosen per degree because widening is not free.** Interleaved
A/B, 4 rounds, min per (build, case), battery with low power mode off, forcing
every degree through `u128` against the shipped routing:

| case | `u64` | `u128` | |
|---|---|---|---|
| `convert_p_to_s` | 0.0124s | 0.0185s | **1.49x** |
| `hall_s_p_degree17` | 0.0123s | 0.0181s | **1.48x** |
| everything not on the p → s sweep | | | 0.98–1.01x |

So an unconditional widening would have cost about half again on every p → s
below the wall to fix the degrees above it. Against HEAD the shipped routing
measures 0.97–1.02x across the same suite: the `u64` path monomorphizes to
what it compiled to before the `Beta` trait existed.

**A first attempt at this measurement was wrong and said 12x.** It forced the
wide path by setting `MASK_LIMIT = 0`, which also moves the threshold
`character.rs` reads for *its* fallback — so `character_beta_sweep_n28` was
timing the per-entry character recursion, not a `u128` sweep. The constant is
shared; a build flag that looks local is not. The corrected experiment changes
only the dispatch in `to_schur`.

Correctness at these degrees has no oracle — Sage cannot compute them, which
is the whole reason the tester asked. Two independent checks stand in.
`the_two_mask_widths_agree_where_both_apply` sweeps every partition of every
degree up to 14 both ways, which is where a transcription slip in the wider
`Beta` impl would show. `the_wide_mask_expands_a_power_sum_to_its_hooks`
checks p_n against the closed form Σ_{r<n} (−1)^r s_{(n−r,1^r)} at degrees 33,
40 and 55 — exact, independent, and cheap because a one-part μ is a single
Murnaghan–Nakayama step. Above those, `s_2[g] + s_{1,1}[g] = g²` at g = s_18
reproduced exactly at degree 36 against a Littlewood–Richardson product, which
shares none of the sweep (737s, so it is an experiment and not a test).

## The one-row route: no conversion at all

The widening moved the wall without changing the shape of the cost, which was
still one p → s at degree n·m. For a **one-row inner argument** that
conversion turns out to be avoidable entirely, and avoiding it is worth two
orders of magnitude. This is now what `plethysm` runs whenever g is `s_m`;
every other g takes the power-sum route unchanged.

The recursion is Newton's, and it stays in the Schur basis:

```text
  n·h_n[g] = Σ_{k=1..n} p_k[g] · h_{n−k}[g]
```

Every product is an ordinary Littlewood–Richardson product — the crate's
fastest primitive — and the outer argument is carried onto the resulting
ladder by its **h-expansion**, since `h_μ[g] = ∏_i h_{μ_i}[g]`. The h-expansion
is what replaces the Jacobi–Trudi determinant the same identity suggests: a
determinant of symmetric functions is factorial in the number of rows, and
`Homogeneous::from_schur` already exists.

What makes the recursion usable is that `p_k[s_m]` has a closed combinatorial
form with no conversion behind it. When every part of the cycle type is
divisible by k, χ^λ vanishes unless λ has empty k-core and otherwise factors
through the k-quotient; for h_m the surviving quotients are exactly the
k-tuples of one-row partitions summing to m:

```text
  p_k[h_m] = Σ ±s_λ   over (a_0, …, a_{k−1}) with Σ a_i = m
```

`C(m+k−1, k−1)` terms — 462 at k = m = 6, against p(36) = 17,977. λ is read
off an abacus with **one bead per runner**, β_i = k·a_i + i, which is enough
because these λ never have more than k parts and keeps the sign O(k²) rather
than O((km)²).

### The sign went wrong three times, the same way each time

Every failure was a **correct support with flipped signs** — the set of λ was
right on the first attempt and stayed right — which is exactly what a
convention error looks like and nothing like what a wrong algorithm looks
like. In order:

1. the bead count was chosen per composition, so terms were compared across
   different abacuses;
2. with the count fixed, the whole sum still carried a constant depending on
   it;
3. the one-bead form reintroduced the same constant, `(−1)^{k(k−1)/2}`, which
   is why it was wrong for k = 2, 3, 6 and right for k = 4, 5 — a pattern that
   reads as a deep bug and is a missing normalization.

The fix each time was to measure the sign against the empty configuration,
which must give λ = ∅ with sign +1. The rule is now bead-count independent,
and `adams_one_row`'s doctest asserts values at k = 2, where the constant is
−1, rather than a k where a missing normalization would pass.

This is the fourth sign slip in this subsystem's history and the third in one
sitting. The lesson the tree already states holds exactly:
a plausible answer is the failure mode, so the check has to be a value that
distinguishes the conventions, never a shape or a count.

### What it is checked against

There is no oracle at these degrees — Sage cannot compute them, which is why
the question was asked. `the_two_routes_are_one_operation` runs both routes
over four inner rows and eight outer shapes and demands term-for-term
agreement; the two share no machinery, one being LR products of an abacus rule
and the other z_μ, characters and a β-mask sweep. `s_2[s_3]` and `s_3[s_2]`
are pinned to Sage-checked values because they have equal degree and different
answers, so a route that transposed its arguments would pass either alone.

The Sage fixture caught the one defect: `s_∅[s_1]` came back **zero**. A
constant outer argument has no largest part, and reading that absent maximum
as "no work to do" dropped the answer instead of returning 1. It is a
degenerate input rather than a mathematical error, which is the kind the
committed fixtures exist for.

### Measured

Release, **AC power**, caches cleared per case, min of 3 for the fast cases
and a single run for the slow ones:

| case | degree | p-route | ladder | | terms |
|---|---|---|---|---|---|
| `s2[s2]` | 4 | 0.000s | 0.000s | 0.7x | 2 |
| `s3[s3]` | 9 | 0.000s | 0.000s | 0.8x | 5 |
| `s4[s4]` | 16 | 0.001s | 0.000s | 3.3x | 28 |
| `s5[s5]` | 25 | 0.162s | 0.007s | **21.9x** | 245 |
| `s4[s8]` | 32 | 5.718s | 0.004s | **1,452x** | 254 |
| `s3[s11]` | 33 | 98.923s | 0.001s | **159,828x** | 72 |
| `s6[s6]` | 36 | 26.932s | 0.271s | **99.4x** | 2002 |
| `s7[s7]` | 49 | not run | 17.9s | | 15293 |

Both routes agree on the term count in every row. The sub-1.0x rows are
sub-millisecond and are noise, but they are the honest shape of the trade: the
ladder builds every rung up to the outer degree, so on an outer argument small
enough that the conversion was never the cost, it does slightly more work.
Nothing reachable makes that matter.

The speedups do not order by degree, because neither route's cost does.
`s3[s11]` is the extreme case in both directions: degree 33 with only 72
terms in the answer, where the p-route pays for a degree-33 conversion over a
dense p-element and the ladder pays for eleven rungs of an inner s_11 whose
outer is a 3.

A first version of the ladder was **2.1x slower than the prototype it came
from**: it rebuilt `adams_one_row(k, m)` inside the rung loop, so each p_k was
recomputed once per remaining rung. Hoisting them out is the whole difference,
and it is the sort of thing a prototype gets right by accident — the Python
version cached them in a dict without anyone deciding to.

### The profile says the new code is not where the time goes

`examples/profile_plethysm.rs`, `ladder 7 7`, 15s of `sample` against the
`profiling` build, leaves bucketed by subsystem (12,552 leaf samples):

| bucket | share |
|---|---|
| Littlewood–Richardson (`skew_lr` fill and merge) | 34.9% |
| Littlewood–Richardson (`three_row` strategy) | 23.2% |
| allocator and `memmove` | 15.8% |
| unresolved or deduplicated symbols | 25.8% |
| **the plethysm code added here** | **0.2%** |

24 samples out of 12,552 for the Adams rule, the composition enumeration, the
division by n and the h-expansion combined. That is the answer to "is there
unnecessary overhead": the route spends its time in Littlewood–Richardson,
which is the primitive it was designed to spend it in, and the combinatorics
that made the route possible cost nothing measurable. `interrupt::poll` does
not appear in the profile at all.

The next lever, if one is wanted, is the 15.8% in the allocator rather than
anything in this section — the same place
[littlewood-richardson.md](littlewood-richardson.md) already found 38% of wall
time on large shapes before `Key` was packed inline.

### Still open

- **The ladder is rebuilt per call.** Within a call it is shared across every
  term of the outer argument, which is where the speedup comes from, but a
  second plethysm over the same g starts again. In the prototype, caching it
  made every later outer shape over the same inner free — `s_5[s_6]` and
  `s_4[s_6]` in 0.000s, and degree-36 shapes like `s_{3,2,1}[s_6]` in 0.085s
  where the general route needs ~30s. The cache wants the `bold_guarded`
  shape: a narrow tier that stores only when the overflow counter agrees.
- **The inner argument must be one row.** A general inner needs the general
  Adams operation, whose k-quotient form is a k-tuple of arbitrary partitions
  rather than rows, and that is **not** verified here.
- Coefficient growth along the ladder is untested at large degree; `s_7[s_7]`
  is the largest case run.
