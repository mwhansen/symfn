# Plethysm

Plethysm is computed through the power-sum basis, where `p_n[g]` is part
scaling and plethysm is multiplicative in `f`. Its performance story is
almost entirely a story about which baseline was being measured against —
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

The lesson is not "always build it". It is that a *cheap* experiment used to
**cancel** work needs the same rigour as one used to justify it, and it did not
get it. A profiler would have settled it in one shot for less effort than the
experiment cost.

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
