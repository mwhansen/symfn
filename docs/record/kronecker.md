# The Kronecker product, ordinary and reduced

Two related subsystems. The ordinary internal product falls out of the
power-sum basis being diagonal for it; the reduced (stable) Kronecker
coefficients come from the Orellana–Zabrocki character bases, where the
structure constants of the `st` basis *are* the stable coefficients.

Split out of [docs/record/README.md](../../docs/record/README.md), which carries the phase plan
and a summary of this file.

---

## The internal (Kronecker) product

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

## The Orellana–Zabrocki character bases, and reduced Kronecker coefficients

`docs/research-gaps.md` §2.2 asked for the `st` basis as a first-class ring, on
the grounds that its outer-product structure constants *are* the stable Kronecker
coefficients and that Sage dies on two two-row partitions of 10.
`docs/record/st-basis-spec.md` is the specification; this is what came out.

The source is Orellana–Zabrocki, [arXiv:1605.06672](https://arxiv.org/abs/1605.06672)
v5, read on 2026-07-29. Sage was used as a black-box oracle and never read, the
same clean-room posture as the LR work.

### Theorem 14 is a map, not a formula

The paper recommends its Eq (23) for computing `s̃_λ`, and read as a formula it is
one: a sum over p(n) partitions of `𝐩_γ/z_γ` weighted by characters. Read as what
it says — that **Γ with `Γ(s_λ) = s̃_λ` is the linear map with `Γ(p_γ) = 𝐩_γ`** —
it is much more, because `𝐩_γ = Π_i 𝐩_{i^{m_i(γ)}}` and each factor is a
univariate polynomial in `P_i = (1/i)Σ_{d|i} μ(i/d) p_d` of degree `m_i` with
leading coefficient `i^{m_i}`.

So the change of basis between monomials in the `P_i` and the `𝐩_γ` is a **tensor
product of univariate triangular matrices**. Γ⁻¹ is then a back substitution
against an (m+1)×(m+1) matrix per part size — not a p(n)×p(n) inversion, and not
a peel over one row at a time. Both directions become whole-element operations:

```text
    s̃_λ · s̃_μ = Γ⁻¹( Γ(s_λ) · Γ(s_μ) ),  read in the Schur basis
```

and the multiplication in the middle is in the **power-sum basis**, where a
product is a multiset union. **There is no Littlewood–Richardson coefficient
anywhere in a reduced Kronecker calculation.**

That is the second time on this project that the win came from *not* using the
most optimized thing available. The spec (§3.3) routed the product through
`Schur::mul` — three backends, `AutoLr`, memoized whole expansions — at ~59²
cached LR products for the target case. It would have been correct and it would
have been slower. The correction is recorded in the spec rather than deleted,
because reaching for the good hammer is the natural mistake.

### Measured, against the wall it was built for

Single runs, this machine, 2026-07-29; ours release-built with caches cleared per
row, Sage 10.9 timed in one process with a 90s `SIGALRM`. ⚠️ Order-of-magnitude
comparisons, not benchmarks.

```text
  case                        terms       ours       sage      ratio
  st[2,1] · st[2,1]              26     0.0003     0.0221        74x
  st[3,1] · st[3,1]              53     0.0006     0.0577        96x
  st[3,2] · st[3,2]             101     0.0019     1.4418       759x
  st[4,2] · st[4,2]             186     0.0030     8.6438      2881x
  st[4,3] · st[4,3]             308     0.0089    30.2569      3400x
  st[5,3] · st[5,3]             525     0.0241        >90         —
  st[6,4] · st[6,4]            1282     0.1447        >90         —
  st[8,5] · st[7,4]            2845     1.0683        >90         —
  st[3,2,1] · st[3,2,1]         211     0.0010     0.1467       147x
  st[4,3,2] · st[4,3,2]        1047     0.0470          —         —
```

`st[4,3]·st[4,3]` at **3400×** is the honest headline, because it is the largest
case where Sage still answers. The three `>90` rows are the point of the exercise
and are deliberately *not* given a ratio: a number divided by a timeout is not a
measurement.

⚠️ **Sage's timings are not stable, in both directions.** Re-measured inside
`scripts/check_st.py`, where a session accumulates state, `st[4,2]²` took 2.42s
instead of 8.64s (3.6× faster) while `st[3,2,1]²` took 10.99s instead of 0.147s
(75× *slower*). The table above uses the one-process sweep for both sides. The
lesson is the one `bench_ops` already carries — interleave and take the min — and
it applies to the oracle as much as to us.

### Where it stops: `z_γ`, not the answer

The engine is fixed-width by default and the wall is at **total degree 24**:
`st[8,5]·st[7,4]` completes, `st[8,5]·st[8,5]` does not.

The overflow is entirely in the **intermediate** rationals. Routing through power
sums divides by `z_γ`, which passes 10²⁶ around degree 26 — and the product
multiplies two of those together, so the denominators square. The *answers* are
nowhere near: the largest coefficient in `st[8,5]·st[7,4]` is 32835, **16 bits**,
growing about 1 bit per unit of `|λ|+|μ|`. Under `bignum` the same call escalates
and answers, and the coefficients stay small there too:

```text
  case                 deg     terms       sec    max coeff
  st[8,5] · st[7,4]     24      2845      1.08        32835   (fixed width)
  st[8,5] · st[8,5]     26      4150     18.18        82994   (escalated)
  st[9,5] · st[9,5]     28      5902     48.00       171959   (escalated)
  st[9,6] · st[9,6]     30      8193    114.02       340166   (escalated)
  st[10,6] · st[10,6]   32     11354    278.76       679234   (escalated)
```

The 1.08s → 18.18s step is not the mathematics getting harder; it is the whole
computation re-running over `BigRational` because one intermediate did not fit.
An obvious future fix is a common-denominator integer formulation — the seam
`Ring::as_ratio` was added for exactly this shape of problem on the `p → s` path —
which would keep the fast path fast well past degree 26.

**One bug worth recording, because it made the escalation dead code.** The first
version detected a bad intermediate by asserting the final coefficient was an
integer. That assert fired *inside* the closure passed to `guarded`, so the panic
escaped before `guarded` could report `None`, and the `bignum` build failed
exactly where the plain build did — the escalation existed and was unreachable.
Detection has to be a returned `Option`, not a panic. Worse, the same code had a
real silent-wrong-answer hole: `Rational` wraps in release, and a wrapped
intermediate can perfectly well land on denominator 1 and be accepted. Running
over `GuardedRat` inside `guarded` closes it.

`memo`'s `𝐩_γ` cache needed a matching change and is the only table here that is
peek-and-store rather than `lookup`: a value computed during an overflowing call
must not be cached, or a later reader sees a clean counter and accepts garbage.

### What checks it

Six layers, because at these sizes there is nothing left to ask:

1. **Published values.** OZ Eq (20) and Eq (21) are printed in the paper and are
   unit tests — a check against the literature rather than against ourselves.
2. **Thm 1(3)**, `s̃_{1^r} = Σ(−1)^i e_{r−i}`, which pins the *embedding* rather
   than any structure constant.
3. **Round trips** on every partition to size 8 (`st`) and 6 (`ht`).
4. **Top degree is Littlewood–Richardson**: `ḡ^ν_{λμ} = c^ν_{λμ}` when
   `|ν| = |λ|+|μ|`, which ties the newest engine to the oldest one with no step
   in common.
5. **Stability against `ops::kronecker`** — the reduced coefficient must be the
   ordinary one once n is large. The only test touching the existing Kronecker
   code.
6. **Sage, on 208 comparisons, 0 mismatches** (`scripts/check_st.py` +
   `examples/stdump`): both transitions to degree 6, all 100 product pairs to
   degree 4 a side, and `st[3,2]²`, `st[4,2]²`, `st[3,2,1]²`, `st[4,3]²`
   individually. `st[5,3]²` is reported as skipped rather than silently dropped.

The escalated path has no oracle at all — nothing on this machine computes total
degree 26 — so it is held to (4) and to non-negativity, in an `#[ignore]`d test
(`cargo test --release --features bignum -- --ignored`); it is ~17s in release
and ~220s in debug, and nothing escalates below degree 26 to make it cheaper.

### The second route is weaker than intended

The `h̃` basis carries an independent product rule, derived here rather than read:
tensoring two Young permutation modules gives orbits indexed by non-negative
integer matrices with row sums `(n−|λ|, λ)` and column sums `(n−|μ|, μ)`, so

```text
    h̃_λ · h̃_μ = Σ_A h̃_{entries of A, minus the (0,0) corner}
```

with only the corner depending on n. Checked against Sage on 9 pairs before it
was implemented, exact on all of them.

It was meant to be the cross-check that survives past Sage's wall, and it does
not. The matrix count explodes on long partitions — λ = μ = (1¹⁰) is 11¹⁰ — and
`s̃_λ → h̃` *always* produces `(1^k)` terms, so the explosion is reached from any
input of interest. `reduced_kronecker_via_ht` refuses (budgeted) rather than
hanging, and the two routes are held to agreement only at degree 3 a side. Making
an independent route reach the frontier is the open problem, and it matters more
than usual: past `st[4,3]·st[4,3]` there is no third-party package left to ask.

### Next

- A common-denominator integer formulation, to push the fixed-width wall past
  degree 26 and delete most of the escalation cost.
- Profile which of the two Γ's dominates. Neither has been measured; the whole
  cost model above is "the Γ's are what cost anything", which is an inference
  from the multiplication being free, not an observation.
- An independent route that reaches the sizes the engine reaches (see above).
- Python bindings, following the whole-object rule: `reduced_kronecker_product`
  and the two transitions, not per-coefficient calls.
