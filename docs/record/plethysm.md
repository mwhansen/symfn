# Plethysm

Plethysm is computed through the power-sum basis, where `p_n[g]` is part
scaling and plethysm is multiplicative in `f`. What it records is almost
entirely about which baseline was being measured against —
the `py` rows it opens by criticizing are the Sage ladder in
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

## Open: computing `s_n[s_m]` without the degree-`nm` conversion

The widening moves the wall; it does not remove the shape of the cost, which
is still one p → s at degree n·m. The route that would remove it is the
Newton recursion, which stays in the Schur basis throughout:

```text
  n·h_n[g] = Σ_{k=1..n} p_k[g] · h_{n-k}[g]
```

with s_6 = h_6, and the products ordinary Littlewood–Richardson — the crate's
fastest primitive. What it needs is p_k[s_m] in the Schur basis without a
general conversion, and that appears to exist. When every part of the cycle
type is divisible by k, χ^λ vanishes unless λ has empty k-core and otherwise
factors through the k-quotient, which collapses the Adams operation to

```text
  p_k[h_m] = Σ ±s_λ  over λ with empty k-core whose k-quotient is a
             k-tuple of one-row partitions summing to m
```

— C(m+k−1, k−1) terms, so 462 at k = m = 6, against p(36) = 17,977.

**Checked numerically, not assumed**, against the existing route (p_k in the
Schur basis is the hook sum Σ_r (−1)^r s_{(k−r,1^r)}, and plethysm is linear
in its outer argument): all nine (k, m) with k ∈ {2,3,4}, m ∈ {2,3,4} agree
exactly. The support was right on the first attempt and the **signs were not**,
twice, which is the normalization trap this tree keeps meeting:

1. the bead count was chosen per composition, so terms were compared across
   different abacuses — a wrong sign on a right support;
2. with the count fixed, the whole sum still carried a constant depending on
   it, until the sign was measured against the empty configuration, which must
   give λ = ∅ with sign +1.

Both showed up as a *correct set of λ with some signs flipped*, which is
exactly what a convention error looks like and nothing like what a wrong
algorithm looks like. The rule is now bead-count independent, which is
asserted in the experiment rather than argued.

Not built. What it needs before it is: the inverse k-quotient map (build λ
from an empty core and a k-tuple of rows) does not exist in the tree —
`k_core_quotient` goes the other way — and the recursion's cost is then
dominated by LR products of degree-36 Schur elements, which is a different
profile from anything measured here and could be worse. The experiment is
`scripts/`-shaped work, not a kernel change, until those two are answered.
