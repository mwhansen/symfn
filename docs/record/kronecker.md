# The Kronecker product, ordinary and reduced

Two related subsystems. The ordinary internal product falls out of the
power-sum basis being diagonal for it; the reduced (stable) Kronecker
coefficients come from the Orellana–Zabrocki character bases, where the
structure constants of the `st` basis *are* the stable coefficients.

Split out of [the record index](README.md), which carries the phase plan
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
nothing enumerates anything. So a hard object costs exactly what s ↔ p costs,
and no more.

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

## A single coefficient, without the product

`docs/research-gaps.md` §2.1 records that no package has a single-coefficient
Kronecker query — every one of them computes the whole product to read one
number, the same defect already found and fixed for Littlewood–Richardson. Here
that defect was explicit: `kronecker` documented itself as costing the same as
`internal`, because the power-sum route produces every ν at once.

The fix was already in the test suite. The orthogonality formula
`g^ν_{λμ} = Σ_ρ χ^λ(ρ)χ^μ(ρ)χ^ν(ρ)/z_ρ` was already there as the *oracle* for
`internal` — it exists because the product route rests entirely on p-basis
diagonality plus s ↔ p, so a wrong identity would be self-consistently wrong.
Read the other way it is an algorithm. `ops::kronecker_via_characters` is that
reading: three character rows and a weighted dot product, no symmetric function
ever built.

What it drops is `p → s`, which expands every `p_ρ` into every λ ⊢ n — the
p(n)×p(n) work and the 1.1 GB. Cost becomes 3·p(n) Murnaghan–Nakayama
evaluations, heavily shared because `try_character` memoizes and the recursion
re-enters itself; memory becomes O(p(n)).

| n | g | product route | character sum | ratio |
|---|---|---|---|---|
| 8 | 4 | 266µs | 145µs | 1.84x |
| 12 | 28 | 1.57ms | 208µs | 7.55x |
| 16 | 28 | 22.0ms | 630µs | 34.9x |
| 20 | 28 | 133ms | 1.39ms | 95.4x |
| 24 | 28 | 990ms | 2.97ms | 334x |
| 28 | 28 | 6.94s | 26.7ms | 260x |
| 32 | 28 | 43.6s | 75.3ms | 579x |
| 36 | 28 | (budget) | 186ms | — |
| 40 | 28 | (budget) | 463ms | — |
| 44 | 28 | (budget) | 1.32s | — |

`examples/bench_kron_coeff.rs`, release with `bignum`, λ = (n−5,3,2),
μ = (n−6,4,2), ν = (n−4,3,1). **Both routes run over `BigRational`**, so every
row is exact and the comparison is like-for-like. `(budget)` is the bench
declining to spend minutes on a p(n)² route, not a limit of it.

**The crossover is ring-dependent, and an earlier version of this section got
that wrong.** Measured first over fixed-width `Rational` — where both routes are
exact only to n ≈ 24 — the product route won below n ≈ 12 (0.34x at n = 8). Over
`BigRational` it never wins: bignum arithmetic costs the product route far more,
because it does asymptotically more of it. Both statements are true of their own
ring, and neither generalizes. The product route is still the right one whenever
more than a few ν are wanted, since it produces them all at once.

The 7x step between n = 24 and n = 28 in the character-sum column is the
escalation switching on: past it every call pays a discarded fixed-width pass
before the bignum one.

Note the coefficient itself is 28 from n = 12 up — Murnaghan stability, visible
in the table for free.

Two things this is not. It is **not asymptotic** — p(n) grows like exp(c√n), so
this is subexponential, and computing Kronecker coefficients is #P-hard either
way. And it is **not** `docs/research-gaps.md` §2.1, which asks for the
polynomial-time bounded-row algorithms (Christandl–Doran–Walter lattice-point
counting; Panova, arXiv:2502.20253, still unimplemented anywhere). This routine
is the intended *oracle* for those if they get built: it is exact, it shares no
code with a lattice-point method, and it reaches sizes where the product route
cannot answer at all.

### At the Python boundary

Exposed as `symfn.kronecker_coefficient(lambda, mu, nu)`, standing to
`internal_product` exactly as `lr_coefficient` stands to `schur_multiply` — the
same defect, fixed the same way, and worth the symmetry in the API for that
reason.

No feature gymnastics were needed: `python` already implies `bignum`, so the
gating below is invisible from Python. The return is `Coeff`, not `i128`, so a
coefficient past the fixed width comes back as a Python `int` rather than
hitting the Rust signature's ceiling.

Against Sage's own `itensor` on the same coefficient, same machine, through the
wheel:

| n | symfn | Sage `itensor` | ratio |
|---|---|---|---|
| 12 | 0.0002s | 0.008s | 40x |
| 16 | 0.0004s | 0.124s | 310x |
| 20 | 0.0011s | 0.929s | 845x |
| 24 | 0.0028s | 6.570s | 2350x |
| 28 | 0.0272s | 43.818s | 1610x |
| 32 | 0.0729s | **>120s timeout** | — |

λ = (n−5,3,2), μ = (n−6,4,2), ν = (n−4,3,1); `SIGALRM` at 120s. This is the
comparison a user actually faces, and unlike the in-crate table it is not
like-for-like on purpose: Sage has no single-coefficient path to offer, which is
the point `docs/research-gaps.md` §2.1 was making.

⚠️ Sage's `itensor` is Python, not C, so the caveat above this file's first
table applies here too — treat these as "the wall is in a different place", not
as a compiled-baseline result.

`scripts/check_bindings.py` checks the binding separately from the library, on
four **deliberately asymmetric** triples: g is symmetric in its three indices,
so three same-shaped arguments marshalled in the wrong order give a plausible
number rather than an error, and symmetric inputs would hide it. Plus the
trivial-character identity at n = 40, where Sage cannot answer at all.

### The z_λ ceiling, and removing it

`Partition::z` returns `u128`, and z_{1^n} = n!. **34! ≈ 2.95e38 is the last one
that fits** (the ceiling is 3.40e38); 35! ≈ 1.03e40 is not. Past that it wrapped
in release and panicked in debug, undocumented. An earlier version of the test
here asserted that 34! wrapped — it does not, and the pin now records the
boundary instead of an inequality around it.

That ceiling was capping far more than this routine. `PowerSum::from_schur` —
the `s → p` half of *every* conversion — divided by `mu.z()`, and `ops::internal`
multiplied by `C::from_u128(lambda.z())`. Both were therefore silently wrong
above degree 34 for any coefficient ring, bignum included.

Two escapes, because the two directions want different things:

- **Multiplying by z_λ** — `Partition::z_in::<C>()`, accumulating in the
  coefficient ring. The same seam as `character_in`: exact for a bignum ring,
  and for a fixed-width one the limit is the caller's choice of ring rather than
  the method's.
- **Dividing by z_λ** — `Partition::div_by_z`, which divides by each part and
  each multiplicity separately and so never forms z_λ at all. Every divisor is
  ≤ n.

A bignum `z()` alone would not have sufficed, and this is the reason that
decided it: `QAlgebra::div_u128` takes a `u128` *by design*, because the trait's
whole point is that the library never divides by a ring element — that is what
keeps ℚ[t] and ℚ[q,t] eligible as coefficient rings. Widening it to accept a
bignum divisor would have bought degree 35 at the cost of the trait. The
division schedule buys it for nothing.

The run-length scan the three routes share is `for_each_part_multiplicity`,
taking a closure rather than returning a `Vec`: `div_by_z` runs once per term of
every `s → p`, so an allocation there would be a real cost paid for tidiness.
`tests/memory.rs` holds the allocation counts that would have caught it.

The running sum is the real ceiling, and it is far lower than the n ≈ 58
character ceiling: measured, plain `Rational` returns confident nonsense from
**n ≈ 26**, because the partial sums are rationals whose denominators divide
lcm(z_ρ) even though the answer is a small integer. Same shape as the `st`-basis
wall recorded below — intermediates, not answers. So `kronecker_coeff`
runs over `GuardedRat` and escalates to `BigRational`.

### Why `kronecker_coeff` requires `bignum` rather than panicking without it

The first version existed in every build and panicked when it could not
escalate, matching `character_basis::escalating`. That is right *there* — the
`st` basis is useful over its whole fixed-width range and the wall is at total
degree 24. It is wrong here. A single-coefficient query is wanted precisely at
the degrees where the whole product does not fit, and the fixed-width path stops
being trustworthy at n ≈ 26 — so a non-bignum build could serve almost none of
the function's reason for existing. An absent function states that; one that
panics on most of its inputs does not.

The default build is unaffected in the way that matters: still zero
dependencies, and it keeps `kronecker`, which over a fixed-width ring is the
faster route below the crossover anyway. `kronecker_via_characters` also stays
ungated, for a caller bringing its own exact ring.

## The Orellana–Zabrocki character bases, and reduced Kronecker coefficients

`docs/research-gaps.md` §2.2 asked for the `st` basis as a first-class ring, on
the grounds that its outer-product structure constants *are* the stable Kronecker
coefficients and that Sage dies on two two-row partitions of 10. This is what
came out.

The source is Orellana–Zabrocki, [arXiv:1605.06672](https://arxiv.org/abs/1605.06672)
v5, read on 2026-07-29. Sage was used as a black-box oracle and never read, the
same clean-room rule as the LR work.

### Prior art, and what is actually ours

| who | what they have | where it stops |
|---|---|---|
| **Sage** | `st`, `ht`, `o`, `sp` bases; transitions; product | The product: `st[4,3]²` takes 30s, `st[5,3]²` exceeds 90s (below). |
| **Stembridge `SF`** | user-defined bases, so this is expressible | No reduced-Kronecker support out of the box; Maple-speed. |
| **`barvikron`** | Christandl–Doran–Walter lattice points, polynomial time for *bounded height* | Kronecker, not reduced Kronecker; Python prototype, unmaintained. |
| **Baldoni–Vergne–Walter** | vector partition functions, bounded length | Maple, distributed as research code. |
| **lrcalc** | LR only | No Kronecker coefficient of any kind. |

The mathematics is entirely prior art and mostly Orellana–Zabrocki's own: Eq (23)
is their recommended computational route, and the stable Kostka transition is
their Eq (7)–(8). Nothing above reaches past where Sage already stops; what
follows is the product, past that wall.

### Transitions are cheap, the product is not

Measured before any Rust existed, to find the actual wall: SageMath 10.9, this
machine, 2026-07-29, single runs, per-item `SIGALRM` of 90s. ⚠️
Order-of-magnitude walls, not benchmarks.

```text
  st[λ] → s                          s[λ] → st
  λ            terms      sec        λ            terms      sec
  [2,1]            4    0.0040       [2,1]            5    0.0049
  [4,2]           16    0.0077       [4,2]           18    0.0296
  [3,3]           15    0.0274       [5,3]           33    0.1355
  [5,3]           33    0.0260       [6,4]           60    0.5738
  [6,4]           59    0.0693
  [4,3,2]         50    0.0385
```

Both directions cost hundredths of a second through degree 6, in Sage itself —
**the transition was never the problem**. An earlier sketch of this work assumed
otherwise and planned to attack the transition; measuring first redirected the
whole design toward the product, which is where the wall in fact sits (below).

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
most optimized thing available. The design draft routed the product through
`Schur::mul` — three backends, `AutoLr`, memoized whole expansions — at ~59²
cached LR products for the target case. It would have been correct and it would
have been slower. The correction is recorded here rather than dropped, because
routing a product through the most optimized product code is the natural first
choice.

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
hanging, and the two routes are held to agreement only at degree 3 a side.
Getting an independent route out to that same computational wall is the open
problem, and it matters more than usual: past `st[4,3]·st[4,3]` there is no
third-party package left to ask.

### A recorded dead end: Littlewood's triple-LR formula

An earlier draft of the design was going to build the product on Littlewood's
triple-LR formula, `ḡ^ν_{λμ} = Σ_{α,β,γ} c^λ_{αβ} c^μ_{αγ} c^ν_{βγ}`, on the
strength of it being the formula everyone quotes for reduced Kronecker
coefficients. **It is wrong in that generality**, and the check that refuted it
is one line: the formula forces `|ν| ≡ |λ|+|μ| (mod 2)`, while `s̃_1 · s̃_1`
contains `s̃_1`. Tested against Sage on 24 pairs, it disagreed on all 24 — a
comparison the draft had not run before building the design on the formula.

### The engine reaches Python, and Sage

Five entry points, following the whole-object rule: `st_multiply` (two whole
elements), `reduced_kronecker_product` (one column), `reduced_kronecker` (one
coefficient, read out of that column rather than by a route of its own),
`schur_to_st` and `st_to_schur`.

Checked against Sage, which has the `st` basis and is therefore an oracle here
rather than only a comparison: the product agrees on every pair of shapes
through degree 3, and both transitions on every shape through degree 4.

Timed in one process against Sage's own `st` product, caches cleared per row,
⚠️ **on battery** — so these are the order-of-magnitude ratios this file's
earlier table already deals in, not benchmarks:

```text
  case                    terms      symfn        sage      ratio
  st[2,1] · st[2,1]          26     0.0001      0.0508       386x
  st[3,2] · st[3,2]         101     0.0017      0.5225       307x
  st[4,2] · st[4,2]         186     0.0013      1.9479      1453x
  st[4,3] · st[4,3]         308     0.0046      6.2729      1355x
```

`st[4,3]²` is the largest case Sage answers at all, and the three rows past it
in the table above are why the bindings exist. These are the *kernel* numbers;
the end-to-end figure through Sage's own dispatch is smaller and is the honest
one to quote, on the same Amdahl argument as everywhere else in
[python-and-sage-interop.md](python-and-sage-interop.md).

### Next

- A common-denominator integer formulation, to push the fixed-width wall past
  degree 26 and delete most of the escalation cost.
- Profile which of the two Γ's dominates. Neither has been measured; the whole
  cost model above is "the Γ's are what cost anything", which is an inference
  from the multiplication being free, not an observation.
- An independent route that reaches the sizes the engine reaches (see above).
- **A single-coefficient query, still unbuilt.** OZ Lemma 20 gives the
  coefficient of `s̃_λ` in `f`, for `r > 2·deg(f)`, as
  `Σ_{μ⊢r} (1/z_μ) s̃_λ[Ξ_μ] f[Ξ_μ]` — both evaluations are integers
  (`s̃_λ[Ξ_μ] = χ^{(r−|λ|,λ)}(μ)` by Theorem 1(1); `h_n[Ξ_μ]` counts weak
  compositions by Prop 24, so anything expanded in `h` evaluates integrally).
  Not built: the sum runs over p(r) partitions with `r > 2(|λ|+|μ|)`, which is
  p(41) ≈ 4.5·10⁴ for the target case — plausibly worse than computing the
  whole product — and now that the whole column is fast enough (above), a
  per-coefficient route has little room to win.

## Offline oracle fixture — the ordinary product

All 505 `(λ, μ, ν)` triples through degree 5, against Sage's `itensor` —
**302 of them zero**.

The zeros are why the sweep runs over triples rather than over the nonzero
support. This route rests entirely on p-basis diagonality plus s ↔ p, so a
wrong identity would be self-consistently wrong, and getting the support wrong
is exactly the error an in-tree cross-route check reproduces on both sides.

## Offline oracle fixture — the reduced product

Sage's `st()` basis is the irreducible-character basis of Orellana–Zabrocki, so
by their Theorem 7 an ordinary product there has the reduced Kronecker
coefficients as its structure constants — a direct oracle for
`reduced_kronecker`, and one this file previously had no equivalent of. 49
products, every pair of shapes through degree 3, checked on every `cargo test`.

The in-tree evidence — the published expansions, the LR top degree, stability
against the ordinary Kronecker product — is real, but none of it is an
independent implementation of ḡ. The gap it left is the one the sweep is shaped
around: the s~ expansion is **inhomogeneous**, terms of every degree up to
|λ|+|μ| appear, and a route that dropped the lower-degree tail would still look
like a plausible product. So the test sweeps the zeros as well, over every ν the
product could reach, and the negative control confirms that a truncated tail
fails.

## Profiling the `h̃` route (`examples/bench_ht_product.rs`)

`reduced_kronecker_via_ht` calls `ht_product_terms` once per pair of `h̃` rows,
and that function's leaf enumerates matrices under `HT_PRODUCT_BUDGET` — up to
4·10⁶ of them — so the per-leaf cost is the route's cost. Nothing had a harness;
`bench_htilde` measures the modified-Macdonald table, a different thing.

Sampled at degree 7: leaf `block` 13.7%, `horizontal_strips::rec` 12.6%,
allocator **24%**, `ht_to_st_row` 8.1%, `memmove` 4.8%, SipHash 3.6%.

The leaf allocated **twice** per matrix — `entries.clone()`, then
`Partition::new` collecting a second vector out of its zero-filter — to produce
a key the map has to own once. It now builds in a reused buffer, sorts and
strips zeros in place, and allocates once.

| case | before | after | |
|---|---|---|---|
| `ht_kronecker_n5` | 0.1023s | 0.0992s | 1.03x |
| `ht_kronecker_n6` | 1.0853s | 1.0631s | 1.02x |
| `ht_kronecker_n7` | 9.1686s | 8.7475s | **1.05x** |

Interleaved A/B, min of 4 rounds, binaries verified distinct.

⚠️ **1.02–1.05x, against a 24% allocator share.** Halving the *count* of
allocations does not halve allocator time: the surviving allocation is the same
size, and much of that 24% belongs to `ht_to_st_row` and the `BTreeMap`, not to
the leaf. The change is kept because it is strictly less work and consistent
across all three degrees, not because it is a win worth repeating the analysis
for. What is left is dominated by the enumeration itself.

### The real defect was a memoization gap, not allocation: 3.0-3.8x

Looking past the leaf at the rest of that profile found something the grep for
allocation patterns could not have. `ht_to_st_row(μ)` sweeps every partition of
every size up to |μ| computing `ht_to_st_coeff`, and `st_to_ht_row(λ)` runs a
back substitution that calls it once per pivot. **Neither was cached**, while
three of the module's sibling row functions — `st_to_schur_cached`,
`schur_to_st_cached`, `reduced_kronecker_cached` — already were. The
reduced-Kronecker route recomputes both once per pair it is asked for, so over
the p(n)² pairs of a degree each row is rebuilt p(n) times.

They now go through `memo::ht_to_st_cached` / `st_to_ht_cached`, the same shape
as the three that were already there, and both are released by `clear_caches`.

| case | before | after | |
|---|---|---|---|
| `ht_kronecker_n5` | 0.1017s | 0.0269s | **3.79x** |
| `ht_kronecker_n6` | 1.0985s | 0.2988s | **3.68x** |
| `ht_kronecker_n7` | 8.9606s | 3.0181s | **2.97x** |

Interleaved A/B, min of 3 rounds, binaries verified distinct, taken
against the single-allocation leaf above rather than against the original — so
these two results compose rather than overlap.

⚠️ **The allocation fix was 1.05x and this is 3.0-3.8x, and the profile pointed
at the allocation.** A sampling profile attributes time to where it is *spent*,
which is inside the recomputation; it cannot show that the recomputation should
not have happened at all. The 24% allocator share was real and was still the
wrong thing to fix first. The question the profile could not raise is the one
that mattered: why `ht_to_st_row` was being rebuilt p(n) times per degree when
three sibling row functions already cached theirs.

Verification is unchanged and independent: `two_product_routes_agree`,
`product_agrees_with_the_schur_route`, `agrees_with_the_ordinary_kronecker_once_
stable` and `ht_and_schur_round_trip` all exercise this route against ones that
do not share its mathematics.

## `s → s̃` is now the whole cost of the character bases in Sage

Sage's `sage.combinat.sf.character` reached these bases by **peeling** — one
leading term removed and expanded per step — and that peel is now intercepted:
both directions of `s ↔ s̃` and `h ↔ h̃` are single whole-element conversions
through `schur_to_st` / `st_to_schur` / `schur_to_ht` / `ht_to_schur`. 6134
small conversions for one degree-16 element became one
([transitions.md](transitions.md)).

Which puts the entire remaining cost inside `schur_to_st_row`, and it is the
`p(n)²` shape: for each ν it builds `PowerSum::from_schur(s_ν)` — p(n)
characters — applies Γ⁻¹, and converts back. One degree-16 element with 199
Schur terms is 1.23s, of which essentially all of it is that.

The fix is the one [Theorem 14 is a map, not a formula](#theorem-14-is-a-map-not-a-formula)
already names from the other side: `r_{νμ}` has a direct description, and the
route through the power sums is a convenience the whole-element caller no longer
needs. Repeat conversions at a degree are already free — 0.006s against
Symmetrica's 0.983s, since the rows memoize — so this is the cold call only, and
it is the last of it.
