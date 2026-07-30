# Macdonald polynomials

`P_λ(x; q, t)` by the branching formula, plus `Q`, `J`, and Hall–Littlewood
`P` by inversion. Symmetrica has no Macdonald polynomials at all, so Sage
is the only external oracle. The interesting engineering is the
coefficient field: a fraction field over ℚ(q,t) that never needs a
bivariate gcd.

Split out of [docs/record/README.md](../../docs/record/README.md), which carries the phase plan
and a summary of this file.

---

`src/macdonald.rs` — `P_λ(x; q, t)` in the monomial basis, by the branching
formula (Macdonald VI (6.24), (7.13')). Symmetrica has **no Macdonald
polynomials at all**, so unlike Hall–Littlewood there is no C implementation to
port from or measure against; Sage is the only external oracle.

Verified against Sage on **every λ through degree 10** — 138 partitions, 1719
coefficients — compared as elements of the fraction field rather than as
strings, which matters because the representations legitimately differ (see
below). In-crate: monic and dominance-triangular, `q = t` gives the Schur
function, `t = 1` gives `m_λ`, and `P_(2)` is checked against a hand computation
that pins both halves of ψ's row/column condition.

**Speed: ~94× Sage** — 0.189s against 17.85s for every shape through degree 9
(`scripts/bench_macdonald.py`). Per degree the ratio runs 81-139x, 81x at the
top degree where Sage takes 12.1s and symfn 0.15s. It was 16x when the layer
first worked; see *Six times over, from three wrong assumptions* below.

## ℚ(q,t) without a gcd

The blocker was the coefficient ring. A general fraction field needs a gcd in
ℤ[q,t] — content and primitive parts over ℤ[q][t], with the coefficient swell
that implies — which would have been the bulk of the work and none of the point.

It is also unnecessary. Every denominator Macdonald produces is a product of
**binomials `1 − qᵃtᵇ`**, and that class is closed under both the product and
the **lcm**, which is all that multiplication and addition need. So `Frac` keeps
the denominator *factored*, as a multiset of exponent pairs, and never expands
it. Cancellation is trial division by a binomial — a short exact loop, using the
fact that multiplying by `qᵃtᵇ` strictly increases the lexicographic key
`QtPoly` already sorts on.

One consequence had to be handled rather than assumed away: **the factored form
is not canonical**, because `1 − qᵃtᵇ` need not be irreducible. `(1+q)/(1−q²)`
and `1/(1−q)` are the same element, both fully reduced against whole binomial
factors, and structurally different. `PartialEq` therefore cross-multiplies. A
derived `PartialEq` would have silently called equal things unequal, and the
Sage comparison would have failed for a reason that had nothing to do with the
mathematics.

## What the profile said this time

Naively — one `Frac::mul` per ψ factor — degree 9 took 4.17s. Three changes,
each from a sampling profile, took it to **1.11s (3.8×)** with byte-identical
output at every step:

| | samples before | what it was |
| --- | --- | --- |
| `memmove` | 2654 | `QtPoly::mul` accumulating with `add_term` |
| `divide_by_factor` | 825 (after the above) | `reduce` running on every `add_assign` |

1. **ψ is a product of ratios of binomials, so count them instead of
   multiplying.** `Frac::from_factors` sums signed exponents first, so factors
   appearing on both sides cancel before anything is expanded. This also
   produced *better-reduced* answers — one coefficient went from 6 denominator
   factors to 5 — which is why the output changed and had to be re-checked
   against Sage rather than assumed equivalent.
2. **`QtPoly::mul` collects, sorts, and combines in one pass.** The double loop
   visits keys in no useful order, so every `add_term` shifted the tail: 2654
   samples in `memmove` against 247 in the multiplication. This is the
   counter-case to the Hall–Littlewood measurement — the `Vec` is still right,
   but *only* with a bulk insertion pattern. Hall–Littlewood's numbers are
   unchanged and its output byte-identical, so the two workloads now agree.
3. **`Frac::add_assign` no longer reduces.** Trial division is the expensive
   operation, and a running sum reduced after every addition pays it once per
   term for a cancellation that can only be decided once the sum is complete.
   Callers accumulate and call `reduce` once; correctness does not depend on it,
   since `is_zero` reads the numerator and equality cross-multiplies.

## Q, J, and Hall–Littlewood P

`Q_λ = b_λ · P_λ` and `J_λ = c_λ · P_λ` are single scalar multiples, with

```text
  b_λ = ∏_{s∈λ} (1 − q^a t^{l+1}) / (1 − q^{a+1} t^l)      c_λ = ∏_{s∈λ} (1 − q^a t^{l+1})
```

so `c_λ` is the **numerator** of `b_λ`, not its denominator (`c'_λ`). Writing the
denominator instead gives a `J_(2)` with leading coefficient `(1−q²)(1−q)` where
Sage has `(1−t)(1−qt)`. That was a real bug, and it was caught in-crate by
asserting the property `J` exists for — that its monomial coefficients are
polynomials — rather than by the Sage comparison.

Hall–Littlewood `P` came from the other direction entirely. The same
Kostka–Foulkes matrix works both ways:

```text
  Q'_λ = Σ_μ K_{μλ}(t) s_μ            s_μ = Σ_λ K_{μλ}(t) P_λ
```

so `P` is what comes out of **inverting** it, with no new enumeration. `K` is
unitriangular in dominance order, so the inverse stays in ℤ[t] and the solve
never divides; visiting the partitions lex-ascending is a linear extension of
dominance, so every `P` the sum needs is already known.

**All four checked against Sage**: `P`, `Q`, `J` through degree 8 (471
coefficients each), Hall–Littlewood `P` through degree 12 (4688 coefficients).

## The cross-check this was all for

`q = 0` turns Macdonald `P` into Hall–Littlewood `P`, and the two are computed
by nothing in common — a branching formula over rational functions on one side,
inversion of the Morris recursion's transition matrix on the other. They agree
for every λ through degree 6.

Compared at several numeric values of t rather than symbolically: `Frac::eval`
substitutes numbers, and holding t formal while setting q = 0 would need a
separate exact division in ℤ[t]. Sage covers the symbolic case.

One assertion had to be weakened, and it is worth recording which: `P_λ(x;1) = m_λ`
holds coefficient by coefficient, but the *supports* differ. A `P` coefficient
can be a nonzero polynomial that vanishes at t = 1, because unlike `Q'` its
coefficients are not sign-definite. Comparing `terms().len()` — which is right
for `Q'` — fails there on a correct answer.

## Six times over, from three wrong assumptions

A second profiling pass over the finished layer. Degree 10, one process, every
shape: **6.45s to 0.80s**; every shape through degree 9: **1.11s to 0.189s
(5.9x)**. All five dumps stayed **byte-identical**, and the Sage checks were
re-run rather than assumed.

Each of the three findings contradicted something the previous pass had left in
place, which is the reason to profile a *finished* thing and not only a new one.

**1. Every multiplication was by a binomial — 39% of the profile was sorting
two already-sorted runs.**

The profile put `quicksort` + `small_sort_general` at 1850 of 4738 samples,
inside `QtPoly::mul`. Rather than infer the operand shape, it was counted:
**100% of the 100k `mul` calls at degree 9 had a two-term operand**, averaging
129 terms on the other side. `Frac` multiplies by `1 − qᵃtᵇ` and by nothing
else, because `from_factors`, `lift` and `denominator` build products of
binomials — so the previous pass's "collect, sort, combine" was quicksorting a
concatenation of two sorted sequences, every time.

`QtPoly::mul_binomial` merges them instead, in one pass over `self` read at two
offsets. **1.8x.** The general `mul` stays for the general case; it is simply
not the case that occurs here.

**2. The reduction inside `from_factors` cost 4x and bought nothing.**

`reduce` trial-divides the numerator by every denominator factor, and
instrumenting it found **72% of those divisions fail**. It ran once per
*tableau*. Removing the call took degree 10 from **3.42s to 0.84s** — and the
output was byte-identical, because the one `reduce` at the end of each
coefficient already reaches the same form.

This is the same mistake the previous pass fixed in `add_assign` and did not
finish: reduction is a per-coefficient operation, and it had been left on a
per-term path. Worth stating as a rule — *in a non-canonical representation,
normalise where the result is consumed, not where it is built.*

**3. Divisibility by `1 − qᵃtᵇ` is a statement about chains.**

Matching coefficients in `N = Q·(1 − qᵃtᵇ)` gives `Q[k] = N[k] + Q[k − δ]`, so
`Q` along a chain `k, k+δ, k+2δ, …` is a running sum of `N`, the chains are
independent, and **the division is exact iff every chain sums to zero**. The
failing 72% then cost what the successes cost, and the `BTreeMap` remainder —
which popped the least key and inserted a larger one per step — is gone.

Worth **1.07x** by the time (2) had removed most of the calls. It was worth much
more before that, and the honest ordering is that (2) superseded it.

Two smaller ones: `mul_binomial` no longer clones before merging (1.05x), and
`Q`/`J` apply their scalar as *factors* rather than expanding `b_λ`/`c_λ` and
running the general product against it (1.22x on those two).

**Two things measured and not done**, recorded because the measurement is the
result:

- Memoising `Frac::from_factors` on the ψ exponent multiset. **61% of the
  multisets are distinct**, and the key is 12.8 factors wide — the hash would
  cost more than the 39% it could save.
- A fast path in `Frac::add_assign` skipping the lift when the accumulator's
  denominator already is the lcm. Measured **flat**, so it was checked whether
  it fires rather than kept on the reasoning that it should.

| | samples | share |
| --- | --- | --- |
| `quicksort` + `small_sort` (was 39%) | 0 | gone |
| `QtPoly::mul_binomial` | 1398 | 46% |
| `divide_by_factor` | 701 | 23% |

What is left is arithmetic that is actually being asked for.

## At the Python boundary

Everything above is now reachable from Python, which is what makes it usable
from Sage: `macdonald_p` / `macdonald_q` / `macdonald_j`, `hall_littlewood_p`
and `hall_littlewood_p_table`, and `kostka_foulkes_table`.

**The denominator crosses the boundary factored**, as `(q_exp, t_exp,
multiplicity)` triples alongside the numerator's terms. That is not a detail of
the encoding — it is the whole design surfacing. A caller writes

```python
d = prod((1 - q**a * t**b) ** m for a, b, m in den)
```

which is the form a fraction field wants anyway; expanding here would mean
factoring again on the other side, and `Frac` exists precisely so nothing has to
factor a bivariate polynomial.

`i128` is not the ceiling it might look like. The widest Macdonald numerator
coefficient through degree 10 is 31594374 — **25 bits against 127**, growing
about 3.5 bits per degree (`examples/mac_coeff_sizes.rs`, which runs each degree
in both widths and compares, since a wrapped `i128` is otherwise silent). The
enumeration becomes impractical long before the width does, so the escalation
path the classical bases carry is not needed here.

`scripts/check_bindings.py` tests **the boundary rather than the mathematics**,
which the dumps already cover. It calls the bindings the way Sage would and
rebuilds the answers as Sage objects, because the failures available here are
different in kind: a denominator marshalled unfactored, a `(q, t)` exponent pair
swapped, a table returned transposed. The last of those is the reason the
Kostka–Foulkes check also asserts the matrix is **not symmetric** — an
orientation test on a symmetric matrix proves nothing, and would have passed
while the table was wrong.
