# (q,t)-Kostka polynomials

Three independent routes to the same table — the branching formula, the
Lapointe–Lascoux–Morse eigenvector, and the Bergeron–Haiman Pieri
recursion Sage actually uses. All three are kept, and a benchmark asserts
they agree at every degree it times.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

`J_μ = Σ_λ K_{λμ}(q,t) S_λ(x;t)` with `S_λ(x;t) = s_λ[X(1−t)]`, so reading
`K` off means inverting that basis change: `φ_t : p_n ↦ p_n/(1−t^n)`, applied to
`J_μ`'s power-sum expansion. Expand `J_μ` in power sums, divide the `p_ν`
coefficient by `∏_i (1−t^{ν_i})`, expand back into Schur.

## The plan said `Plethystic`, and the plan was wrong

The line this section replaces read: *"That needs plethystic substitution over
`Frac`, i.e. a `Plethystic` impl — straightforward, since `p_n` raises q and t
exponents in both numerator and denominator factors."*

It needs no such thing, and writing one and using it would have produced wrong
answers. `φ_t` is **ℚ(q,t)-linear**: it acts on the alphabet and holds the
coefficients fixed. A plethysm does the opposite — `Plethystic::frobenius`
exists precisely because `p_n` substitutes into the coefficient ring too, which
is the one thing that must *not* happen to the q and t already sitting inside
`J_μ`. The two get conflated because both are written `f[X/(1−t)]`. The
`(1−t^n)` in the divisor really is the Frobenius image of `1−t`, but it comes
from `S_λ`'s definition, where the coefficients are rational constants and there
is nothing else to raise.

So the whole step is one loop over the p-basis calling `Frac::mul_factors` with
`(0, ν_i) ↦ −1`. The divisor is never expanded: it is already in the encoding
`Frac` speaks.

## What the tests caught

Two of the ones written first were wrong, both by assuming the Kostka–Foulkes
shape carries over:

- **`K_{μμ} = 1`.** It is not. `K_{(21),(21)} = 1 + qt`.
- **`K_{λμ}(1,1)` is the Kostka number.** It is not — `K_{(11),(2)} = q`, which
  is 1 at `q = 1` while the Kostka number is 0. The table is **dense**;
  triangularity is a `q = 0` phenomenon.

The replacement is much better than what it replaced: `K_{λμ}(1,1) = f^λ`, the
number of standard tableaux of shape λ, **independent of μ** — these are graded
multiplicities in the Garsia–Haiman module, which is a graded regular
representation. Independence of μ is a sharp check, and it reads the whole
polynomial where the `q = 0` test reads one slice. `f^λ > 0` also makes it the
place the density claim is checked.

The other two: `K_{λμ}(0,t)` is Kostka–Foulkes, which crosses to the Morris
recursion in ℤ[t] sharing no code below `Partition`; and
`K_{λμ}(q,t) = K_{λ'μ'}(t,q)`, which is the only one of the four that says
anything about how the two variables are stored relative to each other.
Integrality is checked too — the route divides by `z_ν` and by `1−t^n`, so
arriving back in ℤ[q,t] is Macdonald's theorem and not an invariant this code
maintains. `Frac::into_poly` refuses a surviving denominator; the test refuses a
surviving `1/2`.

Every (λ,μ) pair through degree 7 agrees with Sage's `qt_kostka` — 225 pairs at
degree 7 alone.

## Speed: 1.0× Sage, then 2.8×

Whole table per degree, each in its own process so both sides are cold:

```text
  n   values     before      after       sage    ratio
  6      121     0.0138     0.0061     0.1530    25.1x
  7      225     0.0641     0.0251     0.2556    10.2x
  8      484     0.3283     0.1146     0.5988     5.2x
  9      900     1.4445     0.5113     1.4364     2.8x
```

Sage's first four degrees are ~0.1s of fixed setup, so the trend only means
anything from n=6. The first measurement landed at parity at degree 9, which is
not where this crate usually lands; 2.8× is where two changes put it.

The benchmark had to be rebuilt to see any of this. Sage caches the transition
matrices behind `qt_kostka`, so the obvious per-pair timing loop measures
`dict.__getitem__` for every pair after the first. The cache is not even confined
to one degree: degrees 5, 6, 7 in one process time 0.106s, 0.056s, 0.144s, where
cold they are 0.107s, 0.150s, 0.250s — degree 6 comes out *faster* than degree 5
because degree 5 paid for machinery both share. Hence one process per degree.

## Reduce once, before a value is used many times

The phase split at degree 9 (`examples/profile_qtk.rs`, which times the five
phases rather than sampling, because a flat profile names `divide_by_factor` and
leaves *which phase called it* open):

```text
        J     m->s     s->p    phi_t     p->s   reduce    total
   0.3846   0.0169   0.1112   0.0396   0.0907   0.0001   0.6431
```

`p → s` started at 0.88s of a 1.39s total. Sampling it gave `divide_by_factor`
at 38%, reached from `Frac`'s `Ring::mul` — and the multiplier there is a
symmetric-group **character**, a constant, which over ℚ is a unit and so cannot
make a binomial newly divide anything. Every trial division that ran was doomed
before it started.

So the reduce came out of `Frac::mul`, and the run got **1.7× slower**. The
reasoning about that product was right and the conclusion was wrong. The
coefficient arriving there had never been reduced by anything else —
`mul_factors` does not reduce, `div_u128` does not — so this was where a
denominator first got cut down, and `PowerSum::to_schur` then uses each
coefficient p(n) times, lifting the running sum to an lcm every time. Reducing
was cheap; skipping it was quadratic.

The fix was to move it, not remove it: reduce in `invert_s_basis`, once per `p_ν`
coefficient, before `to_schur` sees it. `p → s` went 0.88s → 0.10s and the total
1.39s → 0.66s. `Frac::mul` keeps its reduce, now with a comment recording why the
apparent inconsistency with `add_assign` and `from_factors` is not one: the
policy is not *never reduce eagerly*, it is **reduce once, before a value is used
many times**.

## A `Rational` that notices it is holding an integer

With `p → s` fixed, the sampler's other two entries were `u128_div_rem` at 29%
and `Rational::add_assign` at 21% — a 128-bit Euclidean gcd on every operation.
(These figures are for the branching route over `Frac<Rational>`, which was the
route at the time and is now the check route; the shipped Bergeron–Haiman route
below runs over `i128` and has no rational arithmetic in it. The shared fix
that finally reached what remained here — the gcd and quotients narrowed to 64
bits, the redundant renormalizations dropped — is in
[coefficient-arithmetic.md](coefficient-arithmetic.md), 1.44–1.46x on this
route.)

Almost none of that arithmetic ever leaves ℤ. A Macdonald `J` over ℚ(q,t) is
integral throughout; the fractions appear only at `s → p`, where `z_ν⁻¹` enters.
`gcd(n, 1) == 1` needs no Euclid to discover, so `add_assign`, `mul` and `new`
now check for `den == 1` first. Degree 9 went 0.66s → 0.48s.

This is deliberately **not** advertised as a general win. `bench_ops` is flat on
it — `hall_s_p_degree17`, `omega_on_h` and all three plethysms move by less than
noise — because those genuinely work in ℚ, where `z_ν⁻¹` is in every coefficient
and the guard never fires. It pays where a ℚ-algebra is being used to hold
integers, which is what `Frac<Rational>` does.

The Macdonald dumps are byte-identical across both changes.

## The curve is real, and tuning will not fix it

Extending the benchmark two degrees settles what four points could only suggest:

```text
  n   values      symfn       sage    ratio
   8      484     0.1138     0.6064     5.3x
   9      900     0.4874     1.4271     2.9x
  10     1764     2.2332     3.8607     1.7x
  11     3136    10.5855     9.4227     0.9x
```

symfn grows **4.3×, 4.6×, 4.7×** per degree; Sage grows a steady **~2.5×**. The
crossover is at n=11. The phase split says where: `macdonald_j` is 67% at degree
11 (6.93s of 10.31s) and is itself growing at 5.5× per degree, while every other
phase grows ~3.5×. So the branching enumeration is both the majority and the
curve.

## Four candidate fixes, all measured, all dead

The plan was to hoist the ψ cache across a table — each `macdonald_p` builds its
own from empty, and different λ share sub-shapes. Sampling `macdonald_j` at
degree 11 ruled it out before a line was written:

```text
  QtPoly::mul_binomial     5212      54%
  divide_by_factor         1462      15%
  memmove                   765       8%
  QtPoly::add_shifted       617       6%
  charge::build             100       1%
  charge::strips             35
  DefaultHasher::write       14
```

`build` + `strips` + hashing is **1.6%**. A perfect cache hoist wins 1.6%.

The sample pointed at `mul_binomial` instead, which looked like `from_factors`
building ψ per tableau. Replacing the accumulation with a `black_box` says
otherwise: `from_factors` alone is **0.13s** of degree 10's **1.23s**. The
binomial multiplications are in `Frac::lift`, inside `add_assign` — 23.7k lifts
at ~9.5 binomials each, every one of them re-expanding a term's numerator up to
the accumulator's denominator.

Three ways to attack that, and instrumentation killed each:

- **Skip the lift when the denominators already match.** They match on **0.8%**
  of calls (212 of 27,502). `lift` clones the numerator even then, which is the
  8% of `memmove` — but 0.8% of the calls is not worth a branch.
- **Bucket tableaux by denominator and sum within a bucket first.** At degree 10
  there are 24,537 tableaux and **16,693 distinct denominators** — 1.5 per
  bucket. There is nothing to group.
- **Reduce the accumulator periodically** rather than once at the end, so its
  denominator stops growing and later lifts are cheaper. Measured across periods:

  ```text
    period      1      4     16     64    256    never
    J (s)   4.926  2.192  1.406  1.157  1.206   1.239
  ```

  Best case **6.6%**, for a tuning constant on a code path that needs replacing.
  Declined. (Period 1 being 4× worse reproduces the earlier `from_factors`
  finding from the other direction.)

The conclusion the measurements force: the lifting **is** the algorithm, not an
inefficiency in it. 24,537 tableaux each carry a rational function with an
essentially unique denominator, and summing them over a common denominator is
what the branching formula asks for. The constant is already close to the floor;
the exponent is the problem.

## Next: Lapointe–Lascoux–Morse

L. Lapointe, A. Lascoux, J. Morse, *Determinantal expressions for Macdonald
polynomials*, International Mathematics Research Notices **1998** no. 18,
957–978 (arXiv:math/9808050) — cited below as **[LLM]**, by their numbering.
Their Theorem 3.1 gives `J_λ` as a determinant over partitions `μ ⊵ λ` whose
entries are explicit Laurent polynomials — no tableau enumeration anywhere.

Two things to be clear about before building on it:

- **Evaluating it as a determinant would be a mistake.** [LLM] 3.9 is Cramer's
  rule for an eigenvector: their 3.7 says the Macdonald operator `M₁` is
  triangular on
  `S_μ[X(t−1)/(q−1)]` with distinct eigenvalues `[|λ|] = Σ q^{λᵢ}t^{n−i}`, so
  back-substitution computes the same thing in O(p(n)²) per shape without the
  intermediate degree swell a `k×k` determinant over ℤ[q,t] carries.
- **It breaks the premise `Frac` rests on.** Back-substitution divides by
  `[|λ|] − [|μ|] = Σᵢ (q^{λᵢ} − q^{μᵢ}) t^{n−i}`, which is *not* a product of
  `1 − qᵃtᵇ`. Every denominator in this library is, which is exactly why no
  bivariate gcd is needed. The way out is to clear denominators through
  `v_λ = ∏_{μ>λ}([|λ|] − [|μ|])`, work over `QtPoly`, and divide exactly at the
  end — which needs `QtPoly::divide_exact` by an arbitrary polynomial.
  `divide_by_factor` is already the binomial special case, and exact division
  with a known-existing quotient is leading-term elimination, not gcd.

## Step one: exact division, without a field

`QtPoly::divide_exact` is in, and it is division and not gcd: the quotient is
assumed to exist and the routine only finds it. Lex order on the exponent pair is
a monomial order — the same well-ordering `mul_binomial` already rests on — so
`lt(qd) = lt(q)·lt(d)`, the remainder's leading term must be divisible at every
step, and the quotient monomials come out in descending order and are reversed
once rather than sorted.

Coefficients needed a seam: eliminating a leading term divides one. `Ring` now
has `div_exact`, defaulting to `None` in the shape `as_ratio` established. It
asks for far less than `Field::inv` — not an inverse, only the answer where one
exists — which is exactly the distinction that lets ℤ have it. Declining is
conservative and never wrong: a caller already has to handle `None`, since *not
divisible* is an ordinary outcome, so a ring that abstains loses divisions it
could have done and cannot produce a false quotient.

**The route needs no field, because every coefficient of `[|λ|] − [|μ|]` is
±1.** Two monomials could only collide if `λ_i = μ_i` at the
same `i`, and then they cancel to nothing instead of accumulating. Lex order is
multiplicative, so `v_λ = ∏([|λ|] − [|μ|])` is lex-monic up to sign and the
elimination never divides by anything but a unit. `QtPoly<i64>` is enough;
ℚ is not required anywhere in the plan. A test asserts this directly through
degree 6 and round-trips `v_λ` against a spread of polynomials — if it ever
fails, `div_exact` over ℤ starts declining and the route needs ℚ after all.

The cross-check that makes the new routine trustworthy is against
`divide_by_factor`, which shares no code and no idea with it — one walks
arithmetic progressions of exponents, the other descends a monomial order. Both
the divisible and the non-divisible case, since a divider that silently returned
a truncated quotient would pass a round trip.

The authors' own basis is `S_μ[X(t−1)/(q−1)]`, not the `S_λ[X(1−t)]` the
(q,t)-Kostka live in; [LLM] Theorem 3.3 reaches the monomial basis but with
entries that are scalar products, and the paper says plainly *"We skip the
problem of computing efficiently all the scalar products in the matrix."* Their
Corollary 3.2's ordinary-Schur form is the useful one, and `φ_t` still has to run
afterwards — but `φ_t` is 3% of the profile, so that is fine.

## Step two: the operator, and an indexing that is wrong in silence

`macop::operator_matrix` builds `M₁` on `{S_μ[X^{tq}]}` from [LLM] 3.6, and
`macop::eigenvector` solves for `J_λ` by back substitution over ℤ[q,t] using
`divide_exact`. Both are verified: the matrix reproduces [LLM] Theorem 3.7
(dominance-triangular, `[|μ|]` on the diagonal) and the eigenvector agrees with
`macdonald_j` — two routes sharing no code, one enumerating tableaux and
multiplying out ψ, the other solving a linear system built from permutation
signs.

Getting there needed one correction that no structural test could have caught.
`[|α|] = Σ_i q^{α_i}t^{n−i}` reads α positionally, and the natural reading — lay
α out by **row** of the Jacobi–Trudi determinant — is wrong. [LLM] 3.5 defines
`M₁` through the formal operators of their 2.4, which add the alphabet `X^t` to
one **column** and sum over which; expanding the determinant, row `j`'s factor
picks up `q^{α_j}` exactly when `σ(j)` is that column, so the power of `t`
travels with the column. Their α is the rearrangement `σ(μ+ρ)−ρ`, which is that
indexing. (Their §1 also numbers rows bottom-to-top, which sends you looking in
the wrong place.)

The row indexing yields distinct eigenvalues, a dominance-triangular matrix with
the right diagonal, and it **satisfies `M₁b = [|λ|]b`** — that equation only says
the solve agrees with the matrix it was handed, never that the matrix is `M₁`.
All three of those tests passed on it. What caught it was a hand computation at
`λ = (1,1)`: the coefficient ratio must be `−(q−t)/(1−qt)` and the row indexing
forces 1 no matter which power of `t` is paired with which position. That is the
argument for the comparison test being against a real Macdonald polynomial and
not against the theory's own internal consistency.

## Step three: measured, and 29× the wrong way

`J_λ` for every shape of a degree, both routes (`examples/bench_llm.rs`):

```text
  n   p(n)    branching  eigenvector
  9     30      0.2380s      7.0724s
 10     42      1.2386s     43.2393s
 11     56      7.3173s    211.5943s
```

Not like for like — the branching side lands in the monomial basis and the
eigenvector side in `S_μ[X^{tq}]`, and the crossing is not written — so this
flatters the new route, and it still loses by 29×.

**The split says the idea is fine and the implementation is not.** At degree 10
the operator matrix takes **0.0038s** and holds 5,630 terms in total; the solve
takes **42.14s**. The matrix is 0.01% of the run. That matrix is the whole of
LLM's contribution — `J_λ` reached with no tableau enumeration — and against the
branching formula's 1.24s for the same degree it is not close.

What destroys it is a decision made three paragraphs after warning against it.
The note above says a determinant over ℤ[q,t] swells to degree ~kn and that back
substitution avoids it — then the back substitution was written on
`b_κ = a_κ · v` with `v = ∏_{κ ▷ λ}([|κ|] − [|λ|])`, clearing every denominator
at once, which is the same swell by another route:

```text
  n   p(n)   |v| terms   max |b_κ|   total b terms
   8     22       8 667       9 699         690 283
   9     30      20 266      22 251       2 746 303
  10     42      48 419      52 097      11 637 890
```

`v` alone is 48,419 terms at degree 10 — an order of magnitude more than the
entire operator matrix — and every one of the p(n)³ multiplications in the solve
runs against polynomials that size, to produce coefficients with a handful of
terms. `divide_exact` is what made the global clearing *possible*, and that
capability is what invited the mistake.

## Step four: two fixes, 36× then 3.3×

**Factored denominators.** `Coeff` in `macop.rs` is [`Frac`](../../src/frac.rs)'s design
over a different family: the divisors are the p(n) known polynomials
`gap_k = [|κ_k|] − [|λ|]`, enumerated before the solve starts, so a denominator
is a multiset of indices into `gap` and combining two needs no gcd. Reduced after
every row, by the argument this session already learned twice — each `b_κ` is
used p(n) times by later rows, so cutting it down once is p(n) lifts avoided.

Degree 11 went **211.6s → 5.9s**, and `v` at degree 10 went 48,419 terms → 2,002.

**Merging instead of sorting.** With the swell gone, sampling put **64%** of the
run in `quicksort` + `small_sort`, inside `QtPoly`'s general `mul`. A uniform
shift is monotone for the lexicographic key — the fact `mul_binomial` already
rests on — so `q^a t^b · other` is *already sorted* and a product is a merge of
`self.len()` sorted runs, not a sort of their concatenation.

That collect-and-sort was itself a fix, for accumulating with `add_term` (which
put 2654 samples in `memmove` against 247 in the multiplication), and it was the
right answer for the workload of the time: one operand was almost always a
binomial, and `mul_binomial` now handles that case without either. Here both
operands are general. Degree 12 went **19.1s → 5.8s**, and the Macdonald dumps
are byte-identical across it.

## Where it stands

`J_λ` for every shape of a degree, both routes over ℤ:

```text
  n   p(n)    branching  eigenvector    ratio
   9     30      0.1802s      0.1091s     1.7x
  10     42      0.8705s      0.3816s     2.3x
  11     56      4.5279s      1.1891s     3.8x
  12     77     27.0584s      3.6336s     7.4x
```

**The curves are what matter**: the eigenvector route grows a steady
**3.1×** per degree, the branching formula **5.5×** and rising. That is the
scaling problem this was started to fix, fixed.

Still not like for like — the branching side lands in the monomial basis and the
eigenvector side in `S_μ[X^{tq}]`, and the crossing is not written. So the 7.4×
is an upper bound, not a result.

Both sides run over `i128` rather than `Rational`, worth ~1.6× to each: the route
never leaves ℤ, and `Rational` was carrying two `i128` fields and a branch to
represent denominators that are always 1. `examples/llm_coeff_sizes.rs` runs it
in two widths and compares, since `impl_ring_for_int` wraps silently — ~6 bits
per degree, 50 bits at degree 12, so `i64` holds to about degree 14 and `i128`
far past where the enumeration is feasible.

## Step five: the crossing, and the route loses anyway

`qt_kostka_table_via_operator` completes the chain, and it agrees with the
branching route term for term through degree 6 — which inherits the Sage check,
since the branching route is held to `qt_kostka` every pair through degree 7.

**The two plethysms compose into one.** Getting from `S_κ[X^{tq}]` to `K` looked
like two substitutions — cross to an ordinary alphabet, then apply `φ_t`:

```text
  p_k ↦ p_k (t^k−1)/(q^k−1)   then   p_k ↦ p_k/(1−t^k)   =   p_k / (1 − q^k)
```

The `t` half cancels outright. Doing them separately would build and then destroy
every `(1 − t^k)` in the expansion, so this is one pass through the power sums,
and `1 − q^k` is a binomial [`Frac`] already holds.

One ordering bug worth recording: `c_{μ'}(t,q)` has to be applied **before**
leaving `Frac`, because the denominators `Ψ` leaves behind are cancelled by it
and by nothing else. At μ = (1) the expansion is `s_1/(1−q)` and `c' = 1−q`, so
asking for a polynomial first fails on the smallest case there is.

**And the whole table is 0.5×.**

```text
  n   p(n)    branching     operator    ratio
   8     22      0.1124s      0.3129s     0.4x
   9     30      0.4724s      1.1596s     0.4x
  10     42      2.2433s      4.8675s     0.5x
```

The `J` half is 2× ahead and the crossing is **87.5%** of the operator route at
degree 10 (4.15s of 4.75s). Both routes run the *same* crossing code; the
operator route just feeds it far more:

```text
  n    branching feeds        operator feeds
   8   max  203, 24 376      max 1 030,  67 505
   9   max  290, 63 392      max 1 618, 182 056
  10   max  404, 166 804     max 2 348, 489 534
```

3× the terms, 5.8× the largest. The cause is the same `v` as before: `solve`
returns `b_κ = a_κ · v`, and `|b_κ|` tracks `|v|` rather than the size of the
answer. Reducing inside the solve stopped `v` from growing quadratically, but it
is still there in the output, and the crossing pays for it.

## Bergeron–Haiman, and the algorithm Sage actually uses

A report on Sage's internals settled what all of the above had been guessing at:
**Sage has not used Lapointe–Lascoux–Morse for the Kostka matrix since 2015.**
That route survives in Sage only for the `J`/`P`/`Q` bases. The engine behind
`qt_kostka` is the Pieri recursion of

> F. Bergeron, M. Haiman, *Tableaux formulas for Macdonald polynomials*, Int. J.
> Algebra Comput. **23** (2013), 833–852 — cited as **[BH]**.

So the LLM work above was aimed at the algorithm Sage abandoned, and the honest
reading of "0.5× Sage" was never a verdict on LLM against Sage's method. It was
LLM against a *third* thing.

### The recursion

`H̃_μ = Σ_ν L_{μν} m_ν` with `L_{μν} = ⟨H̃_μ, h_ν⟩`. Peeling the smallest part `r`
off ν and expanding `h_r^⊥ H̃_μ = Σ_γ c⁽ʳ⁾_{μγ} H̃_γ` gives

```text
  L_{μν} = Σ_{γ ⊆ μ, |γ| = |ν̂|} c⁽ʳ⁾_{μγ} L_{γν̂},        L_{μ,(n)} = 1.
```

One box has a closed form — ratios of `(q,t)`-hook weights over the cells of ν in
the row and column the box vacated, everything else cancelling. More boxes go
through [BH] Proposition 5 and the bi-exponent generator
`B_{μ/ν} = Σ_{(i,j) ∈ μ/ν} t^i q^j`.

This is *not* a per-shape enumeration: it is a recursion over pairs of partitions
ordered by containment, and every value is shared by every μ and ν whose
recursion reaches it. That is the whole difference.

### A second factored fraction field

The hook weights are `q^a − t^b`, which [`Frac`](../../src/frac.rs) cannot hold — it is
closed under `1 − qᵃtᵇ`, and `t² − q³` is not of that shape for any exponent
pair. `bh::Rat` is the same design over the family that does close.

**No field is needed anywhere.** `q^a − t^b` is lex-monic up to sign, so
`divide_exact` never divides a coefficient by anything but a unit; and `m → s` is
the integral inverse Kostka transition. So this is the only one of the three
routes bounded on `Ring` rather than `QAlgebra` — it never divides by an integer,
where the other two carry `z_ν⁻¹` through ℚ.

### The bug at degree 11

`B_{μ/ν}` stopped dividing the Pieri sum at μ = (5,2,2,2), ν = (4,2,1). The cause
is the non-canonicity `Frac` already documents: the atoms are not irreducible —
`q⁴ − t²` is `(q² − t)(q² + t)` — so reducing can cancel a *proper factor* of an
atom against the numerator and leave `B` no longer dividing it. Dividing by `B`
**before** reducing fixes it. Everything below degree 11 is clean either way,
which is how a defect like this survives to ship.

### Where it lands

Whole table per degree, one fresh process each so both sides are cold:

```text
  n   values      symfn       sage    ratio
   9      900     0.0637     1.5303    24.0x
  10     1764     0.2183     4.1120    18.8x
  11     3136     0.6035     9.7494    16.2x
  12     5929     1.7814    24.6699    13.8x
```

From parity to **13.8×**, and the measured Sage times match the report's
independently (24.7s against 22.4s at degree 12). Growth is ~2.9× per degree
against Sage's ~2.45×, so the gap narrows slowly — Sage's curve is genuinely
slightly better and the lead is a constant factor.

Part of that constant is free: the report notes 80% of Sage's time is the
`t → 1/t` substitution applied in the fraction field, once per (μ,ν) pair. This
crate never does it — `K̃` comes out of the recursion directly and `K` is the
`t`-reversal of an integer polynomial, which `macdonald_ht` already did in the
other direction.

`qt_kostka_table` is now this route. The branching and operator versions are kept
as `qt_kostka_table_via_branching` and `qt_kostka_table_via_operator`, and a
benchmark asserts all three agree at every degree it times — three algorithms
sharing nothing above `Partition`.

### Keeping the slow routes, and making the fast one reach everything

Three implementations of the same table is a lot to carry, and the crate already
has a policy for it — stated for [`NaiveLr`](../../src/lr.rs) ("the most
obviously-correct of the three, and the faster ones are held to exhaustive
agreement with it") and again for `kostka_foulkes_by_charge` ("it is here because
it shares no code with the recursion, which makes agreement between the two
evidence rather than tautology"). Both slow (q,t)-Kostka routes are justified
by that standard, and the operator one twice over: it shares no
*mathematics* with either alternative, being an eigenvector problem where the
others sum over tableaux or recurse on containment. `macop::operator_matrix` is
also the only place the Macdonald operator `M₁` exists as an explicit matrix.

What was **not** defensible was the split: `qt_kostka_table` took the fast route
while `qt_kostka`, `qt_kostka_column` and `macdonald_ht` still took the slow one.
That is a trap for callers, not a design, and the measurement says there was
never a trade to make:

```text
  n    one column (branching)    whole table (BH)    crossover
   8            0.0099s               0.0266s        2.7 columns
  10            0.0438s               0.1805s        4.1 columns
  12            1.2952s               1.7050s        1.3 columns
```

The whole table costs 1.3 branching columns at degree 12 and the crossover is
falling, so by degree 13 the entire table is cheaper than a single column the
other way. Every entry point now reads out of the recursion.

`macdonald_ht` got simpler rather than faster-and-more-complicated: Bergeron–
Haiman produces `K̃` **natively**, so `H̃` is now the direct read and `K` is the
one paying for a reflection — the opposite of the arrangement from when `K` came
first.

Deleted: `examples/bench_llm.rs`, which timed `J_λ` in the monomial basis against
`J_λ` in `S_μ[X^{tq}]`. It carried a paragraph explaining that it was not
like-for-like, which is not a fix — a benchmark whose caveat is "this number is
not the comparison you want" can only mislead. `bench_qtk_routes.rs` is the
honest one, and it asserts the three routes agree at every degree it times.

### One cache answers both open questions

The two things left were "what does instantiating at `Rational` cost" and "where
does Sage's slightly better curve come from". They turned out to be the same
question, and the answer was neither.

**`i128` against `Rational` is only 1.2×** (1.29× at degree 6, 1.22× at degree
12) — far less than the 1.6× the same swap was worth on the eigenvector route,
because the `den == 1` fast path added earlier already absorbs most of it when
the arithmetic never leaves ℤ.

**The real cost was p(n) recomputations.** The recursion's unit of work is the
degree, but `qt_kostka` and `qt_kostka_column` are asked for one value or one
column, and each rebuilt the whole thing:

```text
  n    every column one at a time    whole table    wasted
   7            0.0615s                0.0038s       16.3x
   8            0.3519s                0.0158s       22.3x
   9            1.5453s                0.0521s       29.7x
```

That is the mistake `kostka_table` documents, arrived at from the other
direction — and `check_qt_kostka.py`, which asks per pair, was paying it.
`memo::htilde_cached` takes it to **1.7×**, the remainder being the per-call
conversion out of the cache.

It also makes the type question moot. The cache holds `i128`, so the recursion
now runs at `i128` whatever the caller asks for, and `C` only decides the
conversion on the way out. The 1.2× is collected without anyone choosing a ring.

`i128` is safe here by a **bound**, not a measurement: `K̃_{λμ}` has non-negative
coefficients (Haiman) summing to `K̃_{λμ}(1,1) = f^λ`, and `Σ_λ (f^λ)² = n!`, so
no coefficient exceeds `√(n!)` — past `i128` only around degree 57. Measured,
they are 9 bits at degree 12. A test compares the cached path against the
uncached generic one in three rings including `i64`, which is narrower than the
cache, so a value that failed to fit would show up as a disagreement rather than
a silent wrap.

Against Sage, end to end through the bindings, one fresh process per degree:

```text
  n   values      symfn       sage    ratio
   9      900     0.0567     1.5026    26.5x
  10     1764     0.1560     3.8073    24.4x
  11     3136     0.4455     9.3057    20.9x
  12     5929     1.4196    24.4714    17.2x
```

### At the Python boundary

Four entry points: `qt_kostka`, `qt_kostka_column`, `qt_kostka_table` and
`macdonald_ht`, all reaching the recursion. `modified_qt_kostka` is deliberately
**not** bound — the module's rule is whole-object operations only, and a single
`K̃` is one filter away from `macdonald_ht`, which returns the whole column.
(`schur_in_macdonald_j` is the fifth, added later; it is the one that keeps the
`QAlgebra` bound the four below shed, because it round-trips through the power
sums and so cannot run over `i128`. It escalates instead.)

The bounds were wrong until this point and the Python layer was paying for it.
Every one of these was declared `C: QAlgebra` — inherited from the branching
route, which divides by `z_ν` — while the recursion that now backs them divides
by nothing. Relaxing them to `Ring` lets the bindings run over `i128` directly
instead of computing at `Rational` and converting back, and takes degree 12 from
17.2× Sage to **18.4×**.

That also retires an assertion. The old marshalling refused a non-integral
coefficient, on the grounds that landing in ℤ[q,t] was Macdonald's theorem rather
than something the code arranged. Over `i128` a non-integral value is not
representable rather than merely unexpected, and the theorem is enforced where it
belongs: `Rat::into_poly` refuses a surviving denominator and `divide_exact`
refuses an inexact division.

```text
  n   values      symfn       sage    ratio
   9      900     0.0547     1.4695    26.9x
  10     1764     0.1538     3.8005    24.7x
  11     3136     0.4419     9.9721    22.6x
  12     5929     1.3912    25.5895    18.4x
```

`check_qt_kostka.py` now runs to degree 8 — 484 pairs plus 918 `H̃` coefficients
against Sage's own `Ht` — in 4.4s. It stopped at 7 before because the per-pair
cost made 8 impractical.

### Cross-degree cache sharing

`c⁽ʳ⁾_{μν}` and `L_{μν}` are indexed by pairs of partitions of **every** size
below `n`, so a degree-12 run rebuilds most of what a degree-11 run already knew.
Moving both caches out of a per-call `Recursion` and into `memo` shares them:

```text
  degrees 1..=12    clearing between    shared
  time                      2.1203s     1.4649s      1.45x
  peak RSS                    211MB       229MB      +9%
```

The memory is nearly free because the degree-12 call was building that cache
*inside itself* either way — a single degree 12 alone already peaks at 200MB, so
the sharing adds 18MB on top of a cost that was already being paid. All it does
is stop the smaller degrees rebuilding it.

This required making the recursion engine `i128`-only internally, with the
generic conversion moved to the edges, since a `static` cannot be generic — the
same shape `htilde_cached` already had. Safe by measurement rather than by the
`√(n!)` bound, because these are intermediate rational functions and not the
coefficients Haiman's theorem constrains: **25 bits at degree 12**, growing about
3 per degree, so `i128` holds past degree 45.

**It does not explain Sage's curve, and the earlier note claiming it might was
wrong.** Sharing helps a *walk up the degrees*; it does nothing for a single cold
degree, which is what `bench_qt_kostka.py` measures with one fresh process each.
The per-degree growth is unchanged at ~2.9× against Sage's ~2.45×, and the
headline stays where it was:

```text
  n   values      symfn       sage    ratio
   9      900     0.0557     1.4132    25.4x
  10     1764     0.1578     3.8057    24.1x
  11     3136     0.4840     9.6639    20.0x
  12     5929     1.4092    24.5137    17.4x
```

### Next

- The per-degree curve. ~2.9× against Sage's ~2.45×, and where Sage's better
  exponent comes from is still unexplained — cross-degree cache sharing was the
  one candidate tested, and the section above rules it out.
- ~~**`schur_in_j_table` is not memoized**, where `bh::htilde_table` is~~ —
  **closed 2026-08-21**. The obstacle stated here was that the table is
  `Frac<C>` over a generic `C` where `memo::htilde_cached` stores one concrete
  ring; the answer was to key the cache by the ring rather than narrow the
  value, which turned out to be *required* and not merely convenient, because
  `schur_in_macdonald_j` escalates. The resident cost was the other open
  question and is now measured. See "The element-wise form" below.

## The inverse of `J → s` is a projection, not a solve

`schur_in_j_table` expands every Schur function of a degree in the Macdonald `J`
basis — the inverse of the matrix `macdonald_j` fills — and it never inverts
anything. `{J_μ}` is orthogonal for the `(q,t)` scalar product, so each
coefficient is a projection:

```text
  s_λ = Σ_μ ⟨s_λ, J_μ⟩_{q,t} / (c_μ c'_μ) · J_μ
```

and the numerator is a Schur coefficient of something this module already
builds. Against the Hall product `⟨f, g⟩_{q,t} = ⟨f, g[X(1−q)/(1−t)]⟩`, and
`H_μ = J_μ[X/(1−t)] = Σ_λ K_{λμ} s_λ` is the `(q,t)`-Kostka column, so

```text
  ⟨s_λ, J_μ⟩_{q,t} = [s_λ] H_μ[X(1−q)]
```

— one power-sum round trip per μ on top of one [`qt_kostka_table`], against the
`O(p(n)³)` triangular solve over ℚ(q,t) that the matrix inverse costs. The
`X ↦ X(1−q)` step is the same alphabet scaling `invert_s_basis` performs in the
other direction, and it is **not** the plethysm of that name for the same reason
(see the module docs); nothing in `J`'s coefficients gets raised.

Every numerator comes out in ℤ[q,t]: `H_μ[X(1−q)]` divides by nothing, and only
the hook products put anything under the line. The `p(n)` round trips do divide
by `z_ν`, so the route runs over ℚ(q,t) and the boundary refuses a surviving
denominator rather than rounding it.

**Two failed guesses, before the scalar product.** The first was that the table
would be the `(q,t)`-Kostka matrix read transposed over the hook products —
the shape that made Hall–Littlewood `Q'` work, where `s → P` *is* `K` and
`Q' → s` is `K` again. It is not: that would need `⟨s_λ, S_ν⟩_{q,t} = δ`, and
`S_ν = s_ν[X(1−t)]` is dual to `s_ν[X/(1−q)]`, not to `s_ν`. The error is
invisible on the diagonal and wrong everywhere else — degree 2 already separates
them. The second was that a `q ↔ t` swap would repair it; it does not, and the
gap is a whole extra matrix rather than a substitution.

**The check that matters is the matrix product.** `the_schur_table_inverts_the_j_expansion`
multiplies the table against `macdonald_j`'s own Schur expansion through degree
5 and asks for the identity. The two sides share `J` and nothing of how the
inverse is reached, so a wrong hook product, a wrong plethysm or a transposed
index all surface as an off-diagonal entry that fails to cancel.
`the_diagonal_is_one_over_c_lambda` pins `c` against `c'` separately, because
conjugating λ and swapping the variables exchanges them and most of this
module's checks are symmetric under exactly that.

### The element-wise form

Added 2026-08-21, the third of the inverse expansions
([parametric-basis-inverses.md](../plans/parametric-basis-inverses.md); the
others are [hall-littlewood.md](hall-littlewood.md), "The inverse direction",
and [macdonald-operators.md](macdonald-operators.md), "The expansion on its
own"). `schur_in_j_table` holds a whole degree, which is the right unit for
the projection but the wrong shape for a caller holding one element:
`schur_to_macdonald_j` groups the argument by degree, calls the table once per
degree, and returns a `BTreeMap<Partition, Frac<C>>`. No new mathematics at
all — this is the plan's cheapest item, and it is the whole of it.

The encoding at the boundary is the existing `MacTerms`, not the pair form
`s → H̃` needed: these denominators **are** products of `1 − qᵃtᵇ`, since they
are the hook products `c_μ c'_μ`, so they cross factored as `macdonald_p`'s do
and the convenience type is `QtFrac`. `macdonald.to_J` tags the result `McdJ`.
Escalation follows `schur_in_macdonald_j`: guarded rational, then
`BigRational`. Factoring `mac_cell` out of that function's inner loop is what
let both share the integrality refusal.

**Pinned by** the existing `the_schur_table_inverts_the_j_expansion`, which is
the round trip and belongs to the table, plus three tests on the element-wise
form: `s2_and_s11_in_j_are_the_hand_values` (`s_11 = J_11/((1−t)(1−t²))`,
`s_2 = J_2/((1−t)(1−qt)) + (t−q)/((1−t)(1−t²)(1−qt))·J_11`, confirmed against
Sage 10.9 with `SAGE_DISABLE_SYMFN=1`); `the_j_expansion_is_the_table_read_by_rows`
through degree 5, which is the only check that can catch the wrapper reading
the table transposed — a transposed read is still triangular and still
plausible; and a mixed-degree linearity test with the zero element. On the
Python side `check_convenience.py` holds `to_J(s_λ)` against the corresponding
row of `schur_in_macdonald_j`, through a different entry point, for every λ
through degree 5.

**Measured: 45–70× Sage, after a memo that was worth 17×.**
`scripts/bench_inverse.py`, on **AC power**, 2026-08-21; same harness and same
conditions as the `s → H̃` numbers in
[macdonald-operators.md](macdonald-operators.md), including
`SAGE_DISABLE_SYMFN=1` in the Sage arm and one process per degree. Sage
dispatches this to its own Python.

```text
  n  p(n)       sage   per-call      ratio   whole-degree      ratio
  4     5     0.0350     0.0005      70.0x         0.0004      87.5x
  5     7     0.0797     0.0014      56.9x         0.0012      66.4x
  6    11     0.2509     0.0061      41.1x         0.0050      50.2x
  7    15     0.7283     0.0143      50.9x         0.0125      58.3x
  8    22     2.3542     0.0522      45.1x         0.0513      45.9x
```

"per-call" is `schur_to_macdonald_j` once per λ; "whole-degree" is
`schur_in_j_table(n)` once. The first measurement of this had them 17.4× apart
at degree 8 — close to p(8) = 22 — with the per-call ratio at **2.6×** and
falling with every degree, because `schur_in_j_table` was rebuilding the whole
table on every call. `bh::htilde_table` had been memoized all along, which is
the entire reason `s → H̃`'s curve widened away from Sage while this one
narrowed toward it; the two entry points are the same shape and only one of
them happened to sit on a cached callee.

`memo::schur_in_j_cached` closes it, and the two columns are now the same
number. It is not shaped like the other tables here, in two ways that both
matter:

* **Keyed by the coefficient ring** (`(TypeId, degree)`), where
  `htilde_cached` stores one concrete `i128` table and converts at the edges.
  These values are `Frac`s over ℚ(q,t) rather than integers, so that trick has
  nothing to narrow to — and `schur_in_macdonald_j` *escalates*, so a
  ring-blind cache would hand the `BigRational` pass the guarded pass's values
  and make the escalation nominal. That is the failure the module docs above
  describe for `htilde_table`, arrived at from the other side.
* **The store is conditional.** Over `GuardedRat` a value computed while the
  overflow counter moved is garbage, and caching it is worse than recomputing
  it: a later reader inside a clean `guarded` window sees an untouched counter
  and accepts it. The counter is read on both sides of the computation and the
  entry stored only if it did not move — `bold_p_peek`'s discipline, done
  inside the computing function rather than by the caller.
  `memo::tests::a_reported_overflow_is_not_cached` pins both halves, since a
  cache that never stored anything would pass the negative one.

**What it costs in memory**, `examples/heapstat.rs` and the `s-in-j` workload,
same machine and power state:

```text
   n     retained     cold peak
   8       0.5 MB       3.8 MB
   9       1.4 MB       9.4 MB
  10       3.7 MB      25.0 MB
```

"retained" is one copy of the table — measured as the allocation total of a
second, warm call, which does nothing but clone what the cache holds. It is
13–15% of what computing the table transiently costs, so the entry is cheap
against the run that produced it. Degree 10 is the largest measured;
`clear_caches` drops it.

**No new fixture was needed.** The offline oracle below already carries 53
`s → J` coefficients through degree 5, and the row test ties the element-wise
form to the table those cover. A live comparison against Sage over the same
range was run once while building this and agreed on all 53 — which is what
the fixture asserts, so it added nothing and is not kept as a script.

## Offline oracle fixture

`check_qt_kostka.py` is the wider, live check; the offline half is 89
`(q,t)`-Kostka pairs through degree 5, 30 `H̃_μ` Schur expansions through
degree 6, and 53 `s → J` coefficients through degree 5, committed and checked
on every `cargo test`.

The `s → J` rows are compared by **cross-multiplying** and not at a generic
`(q,t)`, which is how every other Macdonald row in that fixture is checked.
The denominators here are the hook products — ten binomials at degree 5 — and
`2^a·3^b` over that degree overflows the `i128` inside `Rational` before any
pole is reached, so the point has to go rather than be moved. `a·d = c·b` in
ℤ[q,t] is exact and needs neither.

⚠️ **The generator now refuses to run with the backend enabled.** Sage's
Macdonald, Hall–Littlewood, Jack and character bases reach this library through
it, so a fixture taken without `SAGE_DISABLE_SYMFN` would be symfn quoting
itself — committed, and passing forever. Regenerating under the guard produced
a diff containing *only* the new `s → J` rows, so what was already there is
verifiably Sage's own.

The matrix has no zero entries, so the risk here is orientation rather than
support: the two indices enter asymmetrically and a transposed table is
otherwise plausible. `H̃` carries its own pin — `H̃_{(2)} = s_2 + q·s_{11}`
against `H̃_{(11)} = s_2 + t·s_{11}` is the smallest pair separating it from
`H`, from `J`, and from a `q ↔ t` transpose, all of which agree on `H̃_{(1)}`.
