# Jack polynomials, and the Goulden–Jackson tables

Three engines over a canonical factored fraction field in ℚ(α), and the
Matchings-Jack / b-conjecture tables that no other package computes.
Symmetrica has no Jack at all, so this is a capability gap rather than a
backend swap.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

`src/jack.rs` and `src/afrac.rs` — `research-gaps.md` §2.6. Symmetrica has **no
Jack at all** (zonal only, `zo.c`), so unlike Schubert this is a capability gap
rather than a backend swap.

Three engines, sharing `Partition` and `AFrac` and nothing else:

| fn | route | role |
|---|---|---|
| `jack_p_lb` | the [MOPS] Laplace–Beltrami eigenoperator recursion | the engine |
| `jack_p_branching` | chains of horizontal strips, ψ^α | the cross-check |
| `jack_j_tableaux` | Knop–Sahi's tableau formula | the reference, and the *positive* route |

`research-gaps.md` asked for "Knop–Sahi, the Lassalle recurrences, or the
Laplace–Beltrami eigenoperator, **rather than Gram–Schmidt**". The eigenoperator
route wins, and the reason is arithmetic rather than combinatorics: it
enumerates nothing at all, and its denominator at every step is a *single* atom.

## Sources

Read directly rather than recalled, and cited by tag throughout this file and
in `src/jack.rs` / `src/gj.rs`:

| tag | source | used for |
|---|---|---|
| **[KS]** | Knop, Sahi, *A recursion and a combinatorial formula for Jack polynomials*, [arXiv:q-alg/9610016](https://arxiv.org/abs/q-alg/9610016) (Invent. Math. 128, 1997) | Thm 1.1 (positivity and `u_μ` divisibility); Thm 4.6 (nonsymmetric recursion); Thm 5.1 (tableau formula) |
| **[MOPS]** | Dumitriu, Edelman, Shuman, *MOPS: Multivariate orthogonal polynomials (symbolically)*, [arXiv:math-ph/0409066](https://arxiv.org/abs/math-ph/0409066) (J. Symb. Comput. 2007) | the Laplace–Beltrami eigenvalue `ρ^α_κ`; the moving-box recursion; Lemma 2.16 |
| **[GJ]** | Goulden, Jackson, *Connection coefficients, matchings, maps and combinatorial conjectures for Jack symmetric functions*, TAMS 348 (1996) 873–892 | series (1) `Φ`, (2) `Ψ`; coefficient definitions (4), (5); both conjectures verbatim |
| **[BD]** | Ben Dali, *Integrality in the Matching-Jack conjecture and the Farahat–Higman algebra*, [arXiv:2203.14879](https://arxiv.org/abs/2203.14879) | integrality half of Matchings-Jack |
| **[DF]** | Dołęga, Féray, *Gaussian…/Cumulants of Jack symmetric functions and the b-conjecture*, [arXiv:1601.01501](https://arxiv.org/abs/1601.01501) | ℚ[b]-polynomiality of both [GJ] coefficient families |

## The coefficient field is what the design turns on

Every scalar in the Jack calculus — hooks, ψ-ratios, eigenvalue differences,
norms — is a ratio of **integer-linear forms `uα + v`**. Normalize the atoms to
primitive and three things become true that are false for `Frac`'s `1 − qᵃtᵇ`:
distinct atoms are irreducible and pairwise coprime (so the factored form is
canonical, where `Frac`'s `PartialEq` must cross-multiply); taking the larger
exponent of each atom gives the **exact** lcm rather than a common multiple, so
addition grows the denominator no more than it has to; and a failed
cancellation is refuted on its first step by Gauss's lemma.

⚠️ Skip the primitive part and the answers leave the ring. `E(κ)−E(λ)` is
genuinely non-primitive — κ = (2,2), λ = (1,1,1,1) gives `2α+4` — and dividing
ℤ[α] by a non-primitive form leaves ℤ[α].

`AFrac<C>` is a **ℚ-algebra for any `C`**, including `i128`, which is not one:
dividing by an integer multiplies the scalar denominator and needs nothing from
`C`. That is what lets the engine run over `AFrac<i128>` and still be handed to
`s → p`, which asks for `QAlgebra` because it divides by `z_ν`.

## Measured: 3000–8000×, against a target of 200×

⚠️ **Battery to battery, one session, isolated processes.** The
pre-implementation measurements had `P → m` at n = 11 as 84.5 s on mains, and
this machine drifts ~1.8× — confirmed here rather than assumed: the same cell
re-measured on battery is 155.1 s, a ratio of 1.84. This is one of the three
measurements behind "Power state" in [README.md](README.md).

```text
  whole degree, n = 11 (p(11) = 56 shapes)
  unit          Sage (s)   symfn ℤ (s)        ratio
  P → m          155.119       0.01857        8353×
  J → m          147.428       0.02161        6822×
  J → p          149.242       0.13546        1102×
  norms          148.937      0.000091     1636670×

  the P → m ladder
  n   p(n)     Sage (s)   symfn ℤ (s)     ratio
   9     30       11.104        0.0035    3173×
  10     42       44.297        0.0086    5151×
  11     56      155.119        0.0186    8340×
  12     77        >180         0.0416       —
```

**Sage prices all four units identically** — 147–155 s at n = 11 — which is the
sharper version of the spec's finding that the single shape λ = (n) costs as
much as the whole degree. It is not the *unit* that costs: the `P → m`
transition is, and `J`, the p-expansion and even the closed-form norm all route
through it. Our four differ by four orders of magnitude, because they are
actually different computations.

`J → p` is our slowest because it goes `m → s → p` through the generic
`convert` hub — Murnaghan–Nakayama and Kostka, not Jack work at all.

⚠️ **"So it is the one worth attacking next" is what this paragraph said, and
it was wrong.** Sampling its two consumers says otherwise. In the [GJ] pipeline
`J → p` does not reach the top fourteen frames at all — 96.5% is `phi_slice`,
i.e. the `p(n)³` triple product, split across `reduce` (49%), `Ring::mul` (37%)
and `add_assign` (36%). In Stanley's table it is **94.6%** — but of *repeated*
calls, not slow ones: 27951 conversions for 99 distinct values, because
`jack_structure_constant` recomputes all three p-expansions per triple.
Hoisting them into `stanley_table` took the degree-12 table from **46.6 s to
1.44 s (32×)** without touching `J → p` at all, and put degree 16 in reach:
111804 triples in 74 s.

Deliberately **not** solved by memoizing `jack_j_powersum`. The Python boundary
runs over `Guarded` precisely so an overflowing intermediate is detected, and a
cache filled at `i128` and handed out to other widths would hide exactly
that — the `memo::bold_p` hazard. Hoisting the loop has no correctness
question in it.

The real "next factor" is therefore the `p(n)³` triple product in `gj.rs`, and
`J → p` stays as it is: 1102× is the least impressive ratio in the table and
also the one that costs nothing.

The ratio *grows*, because Sage costs ~3.5× per degree and this costs ~2.0×.
The 200× target was beaten by more than an order of magnitude, and — unlike the
Δ-operator spec, which guessed 100× and got 21× — the guess was low for a
reason worth recording: it priced the arithmetic correctly and the *incumbent*
wrongly.

The wall shift, run to exhaustion:

```text
  n     16      18      20      22      24      26
  s   0.73    2.68   15.17   49.55  161.56  520.66
```

n = 16 is the degree Stembridge's SF ships as **precomputed archives**; it is
0.73 s live here. Sage cannot do n = 12 at all. Coefficients are 89 bits at
n = 26 against `i128`'s 127, growing ~4.4 bits/degree, so the fixed-width wall
is around n = 34 and the two-width ladder (`i128` against `Rational`, term for
term) is what would catch it.

Norms are not benchmarked past making the point: `⟨J_λ,J_λ⟩ = H_λH'_λ` is a
closed product of `2|λ|` linear factors, so the whole n = 11 table is 91 µs
against Sage's 148.9 s. Sage prices a product of 22 linear factors like a full
expansion.

Scored against the targets set before any of this was built:

| requirement | target | measured | verdict |
|---|---|---|---|
| headline `J → m` at n = 11 | ≥ 200× (≤ 0.45 s) | 0.0216 s, **6822×** | **beaten by 34×** |
| whole tables through n = 16 live | "in minutes" | **0.73 s**; n = 26 in 520 s | beaten |
| single shape `P_(20)` | "single-digit seconds", flagged a wild guess | **0.0071 s** | beaten by ~10³ |
| norms table at n = 12 | "microseconds" | 134 µs | met |
| Stanley table, `\|λ\|=\|μ\|=6` | "minutes" | **1.44 s**, 9317 triples | beaten |
| [GJ] `c` and `h` complete for n ≤ 10 | — | **26.3 s** at n = 10 | met |

## What the sampling said

`sample`, per the `Cargo.toml` workflow, on `profile_jack` at n = 18. **58% of
the profile was inside `reduce_at`, and roughly half of *that* was `malloc` and
`free` rather than arithmetic.** `divide_by_linear` built the quotient buffer
before knowing whether the division succeeded, and `reduce` trial-divides by
every denominator atom, so most of those allocations were thrown away.

The arithmetic argument above predicted a failed cancellation would cost "one
dot product". It cost one dot product **and two heap allocations**, and the
allocations dominated. Splitting off an allocation-free predicate
(`divides_by_linear`, a single running scalar — the recurrence never needs the
whole quotient array) and rewriting in place was worth **1.6×**: n = 18 went
4.33 s → 2.72 s.

Three follow-up changes — in-place `lift`, an allocation-free
`content_reduce`, and cached partitions with precomputed eigenvalue statistics
— were worth **~2%, i.e. nothing measurable**, and are kept only because they
strictly allocate less. That is `QtPoly::mul_diff` again: sound reasoning about
a real cost that turns out not to be on the critical path. Worse, the
`content_reduce` rewrite silently **grew** the scalar denominators (reading the
loop bound off the live scale makes the final step try the uncancelled scalar
as one lump, so `202 = 2·101` against a numerator of content 101 keeps its 101
forever). The regression was invisible in the timing column and visible in the
scale-width column, which is why `bench_jack` prints that column.

## The Goulden–Jackson pipeline

`src/gj.rs` computes `c^λ_{μν}(b)` and `h^λ_{μν}(b)` — the Matchings-Jack and
b-conjecture coefficients. **No package computes either table.**
ℚ[b]-polynomiality is a theorem (Dołęga–Féray) and `c`'s integrality is a
theorem (Ben Dali), so both are enforced; **positivity is open for both and is
only observed**, with any negative coefficient reported as a finding rather
than debugged away — the same rule the valley-Delta conjecture gets.

The transcription is pinned at **both** known specializations, each against an
object with an independent definition and neither touching a Jack polynomial:

- **`b = 0` is the class algebra of `S_n`.** `class_algebra_coefficient`
  computes `a^λ_{μν}` from characters alone. Exhaustive through n = 6.
- **`b = 1` is the double coset algebra of `(S_2n, H_n)`.**
  `double_coset_coefficient` counts matchings: `b^λ_{μν}` is the number of `δ`
  with `type(δ₀,δ) = μ` and `type(δ,δ₁) = ν`, for fixed `δ₀,δ₁` of relative type
  λ. ⚠️ The **normalization was measured, not read from [GJ]** — the ratio came
  back exactly 1 on all 285 live triples through n = 5, with no factor of `z_λ`
  or `2^ℓ`, and the test then requires it through n = 6, which is 484 triples the
  constant was not fitted on. Computing this side from *zonal* polynomials would
  have re-used Jack at α = 2 and checked nothing; counting matchings is what
  makes it independent. Cost is `(2n−1)!!` per λ — 10,395 at n = 6, 135,135 at
  n = 7 — so it is a pin, not an engine.

Two engines now compute the tables (see below), and where both are affordable
they are required to agree entry for entry. The modular one carries the ladder:

```text
  n   p(n)  modular(s)  exact(s)   c terms   h terms  b deg  bits
  10     42       2.130         —     54108     28752      9    26
  12     77      20.020         —    358478    193603     11    35
  14    135     186.370         —   2045553   1121377     13    44
```

Every coefficient computed lies in ℕ[b].

**⚠️ But a bulk sign check is close to worthless here, and the ladder now says
so.** Positivity and integrality are already theorems, the degree bound is
characterized (Promyslov), and Ben Dali's marginal sums `Σ_{ℓ(ν)=m} c^λ_{μν}`
are *already known* b-positive with a matchings interpretation — so a
counterexample has to hide inside a marginal sum, with the other terms of that
sum cancelling it. The bar for "verified through n = N" as a remark worth
making is n ≥ 25.

So `gj_tables` reports **how much of its output is not already a theorem**, via
`matchings_jack_coverage`:

```text
  n       open  proved [GJ]  (n)-variant    open %
   4         15           19           31     23.1%
   6        494           55          221     64.2%
   8       5879          137         1026     83.5%
```

`proved [GJ]` is λ = [1ⁿ] or [2,1^{n−2}], where Goulden and Jackson built the
statistic and proved the conjecture outright. `(n)-variant` is one of the three
partitions equal to `(n)`: Kanunnikov–Vassilieva proved μ = ν = (n), and with
Promyslov any one of the three — but for a **variation involving labeled
matchings**, so it is weaker than the column beside it and is counted
separately. Everything else is open; the smallest such triple is
λ = μ = ν = (2,2) at n = 4, and nothing below n = 4 is open at all.

## Closing out §5

Every item on the spec's correctness list is met except two, both recorded
rather than quietly dropped: the Knop–Sahi tableau route stops at n = 5 rather
than 6 (it is exponential, and its job is to be a *different* algorithm rather
than a wide one), and [KS] Thm 1.1 is checked to n = 8 rather than 10.

Three things were nearly missed and are worth naming:

- **The `ω_α` duality** — `ω_α P_λ^{(α)} = Q_{λ'}^{(1/α)}` — is the only law
  relating `P` to `Q`, to conjugation, and to the parameter inversion at once,
  so it is the only independent check on `jack_q`'s normalization; `⟨P,Q⟩ = δ`
  would survive a compensating error in both. It needed `AFrac::invert_alpha`,
  which turns out to stay inside the family exactly: `(uα+v) ↦ (vα+u)/α`, still
  primitive, except that the `α` atom becomes the constant 1 and leaves. So the
  substitution is symbolic and the law never evaluates anything. The test ships
  with its own negative control, because the α-twist is silent when wrong.
- **Offline oracle fixtures.** `check_jack.py` is wider but only runs when
  someone has Sage and remembers; `gen_sage_oracle.sage` now emits 132 Jack
  expansions to degree 7 (776 coefficients) and `cargo test` checks them with
  no Sage installed. Sage's `numerator()` over ℚ(α) returns a polynomial with
  *rational* coefficients, so the generator clears them — a `1/2` token in an
  integer fixture is a parse error, not a wrong answer, but only because the
  parser was strict.
- **`AFrac` over bignum was never instantiated by a test.** It compiled,
  because the Python boundary uses it, and `coeff.rs` warns exactly about this:
  a bound is only checked where it is instantiated. `BigInt` implements
  `div_exact` as exact division and `BigRational` as a field — two different
  meanings, both correct here, neither previously exercised.

## The inverse direction: `m → P`, `m → Q`, `m → J`

Built 2026-08-21, item 4 of `docs/plans/parametric-basis-inverses.md`, after
the Macdonald pair. `monomial_to_jack_p`, `_q` and `_j` in `src/jack.rs`, the
pyfunctions of the same names, `jack.to_P` / `to_Q` / `to_J`, tagged `JackP`,
`JackQ` and `JackJ` as Sage prints them.

Structurally it is the Macdonald pair with `AFrac` in place of `Frac`: `P` is
monic and dominance-unitriangular in the monomial basis, so `monomial_in_p_table`
back-substitutes through `jack_table(n)` and the entry point applies one
degree's table to that degree's terms. `Q` and `J` are the same solve rescaled
by a *factored* product of linear forms — `jack_norm_p(λ) = H'_λ/H_λ` for `Q`,
the reciprocal of `hook_lower(λ)` for `J` — so neither costs a second solve.
The source basis is the monomial one for all three, because that is what all
three are expanded in, so a forward answer feeds straight back.

`AFrac::div_int` is new and exists for the boundary: `AFrac::parts` hands out
`num / (scale · ∏ atoms)` and nothing put a `scale` *back*, which an inbound
coefficient needs. It is also what lets `_jack_rows` in the convenience layer
skip the least-common-denominator round trip `_mac_rows` performs — every Jack
row already carries an integer denominator of its own, so a rational numerator
folds into that row's `scale` and crosses unchanged.

### Memoizing the solve, from the start

The Macdonald item's advice was to memoize from the beginning rather than
measure the loss first, and it was right. `memo::jack_p_inverse_cached` shares
`transition_cached` with the Macdonald and `s → J` tables; the key carries the
table's own Rust type, which carries shape and coefficient ring together, so
the three cannot collide.

What it is worth, cold per degree, `--release`, AC power (harness: a sweep of
every λ of the degree, against the same sweep with `clear_caches()` between
shapes):

| n | p(n) | memoized | uncached | ratio |
|---|---|---|---|---|
| 6 | 11 | 0.0007s | 0.0066s | 9.3× |
| 7 | 15 | 0.0020s | 0.0187s | 9.6× |
| 8 | 22 | 0.0034s | 0.0513s | 15.0× |
| 9 | 30 | 0.0056s | 0.1592s | 28.4× |
| 10 | 42 | 0.0154s | 0.6397s | 41.5× |

The ratio approaches p(n) for the reason it did in the Macdonald case: without
the cache the sweep rebuilds the whole-degree table once per shape.

**Where the cached build's time goes is *not* what Macdonald's was.** There the
solve was 98% of a cold call and the forward table almost free; here the two
are about even — `jack_table(8)` is 0.0028s of a 0.0054s cold call, and
`jack_table(10)` 0.0134s of 0.0241s. The eigenoperator recursion is expensive
relative to a back-substitution over linear forms, where the branching-based
Macdonald table is cheap relative to a back-substitution over binomial
fractions. Neither half is negligible, so there is no single place to optimize.

### Against Sage

`scripts/bench_inverse.py 8`, 2026-08-21, **AC power**, one process per degree
*and per arm*, `SAGE_DISABLE_SYMFN=1` in the Sage environment. Sage's arm is
`Sym.jack(t=a).P()(m(λ))` and its siblings, over the fraction field of
`QQ[a]`; the workload is every λ of the degree.

| n | p(n) | m→P sage | symfn | ratio | m→Q sage | symfn | ratio | m→J sage | symfn | ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 7 | 0.0642s | 0.0002s | 321× | 0.0696s | 0.0002s | 348× | 0.0673s | 0.0002s | 337× |
| 6 | 11 | 0.1377s | 0.0004s | 344× | 0.1542s | 0.0005s | 308× | 0.1488s | 0.0004s | 372× |
| 7 | 15 | 0.3812s | 0.0007s | 545× | 0.3958s | 0.0009s | 440× | 0.3800s | 0.0007s | 543× |
| 8 | 22 | 1.3637s | 0.0022s | 620× | 1.4511s | 0.0027s | 537× | 1.4310s | 0.0023s | 622× |

The Jack arms alone continue past the degree the Macdonald table stops at
(`s → H̃` costs 19 s at degree 8 and the run stops being cheap):

| n | p(n) | m→P sage | symfn | ratio | m→Q sage | symfn | ratio | m→J sage | symfn | ratio |
|---|---|---|---|---|---|---|---|---|---|---|
| 9 | 30 | 4.7828s | 0.0058s | 824× | 4.8662s | 0.0071s | 685× | 4.7499s | 0.0060s | 792× |
| 10 | 42 | 17.4891s | 0.0161s | 1086× | 18.0653s | 0.0188s | 961× | 17.7939s | 0.0167s | 1065× |

The margin is an order of magnitude above the Macdonald directions' 27–51×,
and unlike `m → P` over `ℚ(q,t)` it *grows* with degree. One parameter and
linear denominators is the whole of the reason, and the memory numbers below
say the same thing in a second currency.

### Two harness defects the Jack arms exposed

Both were in `scripts/bench_inverse.py` and both **change numbers already
recorded**, so they are stated here rather than fixed silently.

1. **Sage shares a family's transition matrix between its normalizations, and
   the arms shared a process.** The script's docstring claimed one process per
   degree *and per arm*; the code ran every arm in one process per degree. The
   Jack arms made it visible: `m → Q` and `m → J` read 0.08s and 0.05s at
   degree 8 against `m → P`'s 1.36s, and in isolation all three cost ~1.4s.
   The later arms were reading the matrix the first one built. `one(n)` now
   takes an arm index and the driver spawns one process per arm.

2. **Sage builds a family's coercion machinery on first use, and the first arm
   paid for it.** One untimed degree-1 conversion now runs before the timed
   region. It touches no transition matrix of the degree being measured.

**Correction to `docs/record/macdonald.md`, "The inverse direction".** With
both fixed, Sage's Macdonald `m → P` at degree 8 is 4.5055s rather than the
2.7233s recorded there, and the ratio is **51.4×** rather than 30.9×; `m → Q`
is 4.5957s and **44.1×** rather than 2.8571s and 27.4×. `s → H̃` (19.14s,
122×) and `s → J` (2.31s, 44.7×) are unchanged within noise, being the arms
that ran first. The claim in that section that the `m → P` ratio *falls* with
degree where `s → H̃`'s rises does not survive: on isolated arms it runs
82×, 101×, 79×, 66×, 69×, 51× from degree 3 to 8 — still falling, but from a
much higher start, and the degree-8 figure is now above `s → J`'s. The
sentence there that Sage's `m → P` overtakes its own `s → J` between degrees 7
and 8 remains true and is in fact stronger: 4.51s against 2.31s.

### Memory

`m-in-jack-p` is the same shape as the Macdonald `m-in-p` workload —
`monomial_to_jack_p` of `m_{(5,3,1)}`, 22 terms out — so the two are directly
comparable (`cargo run --release --example heapstat`):

| workload | peak | total | allocs | churn |
|---|---|---|---|---|
| `m-in-jack-p` | 0.4 MB | 5.7 MB | 68 590 | 12.9× |
| `m-in-p` | 7.0 MB | 5515.9 MB | 624 097 | 787.0× |

Same answer shape, **970× less total allocation**. `m-in-p` has the highest
churn in the catalog and this has one of the lowest; the difference is entirely
`AFrac` against `Frac` — a dense `Vec<C>` in one variable with primitive linear
atoms, against a bivariate term map over a binomial multiset.

Retention, read with `measure::live()` immediately after a cold call (the
quantity `peak` cannot report, since it is a high-water mark):

| n | cold peak | retained | share |
|---|---|---|---|
| 7 | 0.11 MB | 0.06 MB | 55% |
| 8 | 0.25 MB | 0.14 MB | 56% |
| 9 | 0.46 MB | 0.26 MB | 57% |

A higher share than the Macdonald table's 29–36%, on a tenth of the absolute
size: less scratch is thrown away because there is less fraction arithmetic to
throw away.

### What holds it up

The plan asked for the orthogonality route as a **second engine sharing no
mathematics**, and it is cheap here, so it is in `cargo test` rather than only
in the record. `⟨P_λ, Q_μ⟩_α = δ_λμ` makes the `P`-coefficient of `f` equal
`⟨f, P_λ⟩_α / ⟨P_λ, P_λ⟩_α`; that route goes `m → s → p` and pairs diagonally,
where the solve never leaves the monomial basis and never pairs anything, and
the norm is a closed product of `2|λ|` linear forms rather than a computation.
`orthogonality_gives_the_same_coefficients_as_the_solve` checks every `m_μ`
against every λ through degree 5.

Also committed: the round trip for all three normalizations through degree 7;
`at_alpha_one_the_p_expansion_is_the_schur_expansion`, which is `P_λ(x;1) = s_λ`
read backwards and compares against the ordinary `m → s` transition through
degree 6; the hand values at free α; and linearity across three degrees.

⚠️ **The `α = 1` check is blind to the `α → 1/α` twist**, which fixes it, and
so is every `Q`-versus-`P` comparison there: `m_11` is `P_11` outright and
`[α(α+1)/2] Q_11`, and both are 1 at α = 1. The free-α hand values are the
only thing separating either pair, which is why they are doctests as well as
tests. Sage prints `-(2/(a+1))*JackP[1, 1] + JackP[2]` for `m_2`, confirmed
directly in the sage-dev environment.

**Confirmed against Sage** by extending the committed `scripts/check_jack.py`
rather than by a one-off dump — which is the gap the Macdonald item's record
flagged, closed here. `examples/jack_dump.rs` emits three new kinds `mp`, `mq`,
`mj`, and the check compares them against `Sym.jack().P()(m(λ))` and its
siblings by value in the fraction field. 351 coefficients — every λ through
degree 6, in all three normalizations — 0 mismatches, alongside the 1212 the
script already compared.

And it is fixtured as well as scripted, so `validation.md`'s "committed
fixtures, not scripts someone must remember to run" is met outright.
`gen_sage_oracle.sage` emits `jminp`, `jminq` and `jminj` — `m_λ` in all three
normalizations for every λ through degree 7: 135 rows, 702 coefficients — and
`monomial_in_jack_matches_sage` in `tests/sage_oracle.rs` reads them with no
Sage installed, at three generic α. ⚠️ Not at α = 1, which separates neither
the normalizations nor the twist.

## Next

- **Neither half of the `m → P` solve is negligible**, unlike Macdonald's,
  where the back-substitution was 98% of it. `jack_table(n)` is about half a
  cold call at degrees 8 and 10. If this direction is ever worth optimizing,
  it needs both the eigenoperator recursion and the solve, and the profile
  above ("Jack is coefficient-bound") says `AFrac::reduce_at` is where the
  first half's time goes.
- **Push the [GJ] tables past n = 10** — the deliverable, and what the engine
  exists for. The cost is `phi_slice`: `Σ_θ` of a rank-1 tensor over `p(n)³`
  entries, so `p(n)⁴` coefficient operations, measured growing ~4.6×/degree
  (faster than `p(n)⁴`'s 3.8×, because the coefficients grow too).

  `examples/probe_gj.rs` measures the four things that decide the fix, and
  **three of them ruled out the design they were testing**:

  | probe | result | verdict |
  |---|---|---|
  | where the atoms come from | `J → p` carries **zero** atoms at every degree; all of them come from one `1/⟨J_θ,J_θ⟩` per θ | true, and useless on its own |
  | do the atoms dominate a multiply? | atom-carrying × plain is **0.6×** the cost of plain × plain | **no** — the numerator polynomial dominates, and atom-carrying values have *smaller* numerators |
  | one global denominator `D = lcm_θ⟨J_θ,J_θ⟩`? | `deg D` = 73 at n = 10 against `2n` = 20, growing ~n² | **no** — accumulating over one fixed `D` blows the degrees up |
  | evaluate at numeric α, interpolate? | ℚ: **5.4×** per op, needing ~n+2 points | **no** — a net *loss* |

  So "hoist the atoms out of the inner loop" — the obvious move given the first
  row — is worth nothing, and the fourth row is worth less than nothing.

  **The one that survives is the same idea with a cheaper scalar.** `Rational`
  is a *slow* scalar: it runs a 128-bit gcd per operation. Modular arithmetic
  does not:

  ```text
    mul, mul, add        AFrac<i128>   3512 ns
                         Rational       648 ns      5.4×
                         mod p (2⁶¹−1)  9.6 ns    366×
  ```

  366× per operation against ~n+2 points is **~26× net at n = 10, and it grows**
  — `AFrac` operations get more expensive with degree while a modular one stays
  flat, and the point count only grows linearly.

  Two things make this safe rather than a guess. **There are no poles for
  α > 0**: every atom is `uα + v` with `u, v ≥ 0` and not both zero, so any
  positive α is a legal evaluation point — a proof, not a sampling argument.
  And the object to interpolate is the *final* `c` and `h`, which are
  polynomials in `b` of measured degree ≤ n−1, not the intermediate `Φ`, which
  is not a polynomial at all. So the whole pipeline — including the `log`
  recurrence — runs at a numeric α and only the answers are reconstructed.

  **This was built** — `src/gjmod.rs`, with the reusable half extracted to
  `src/modular.rs` — and it took the ladder from n = 10 to n = 14. Three of the
  planned details turned out wrong, all corrected by sampling:

  | planned | measured | what shipped |
  |---|---|---|
  | a 61-bit prime, `u128` remainder | **87.6% of the engine inside `__umodti3`** — `u128 %` is a call into `compiler_builtins` on aarch64, not an instruction | 31-bit primes so the product fits a `u64`, with Barrett reduction. **3.3×** |
  | interpolate each output entry | **78% in `interpolate`, 10.5% in the shift**, with the pipeline absent from the profile | one Lagrange×shift matrix per prime, not per key. **8–15×** |
  | one prime plenty, second is the check | coefficients grow ~3 bits/degree and hit **44 at n = 14** against the `2^46` three primes give | three primes CRT'd, a fourth held back, and the count **grows on demand** rather than reporting our own range limit as non-polynomiality |

  ⚠️ A 128-bit divide instruction would not have changed the first row. x86-64's
  `DIV r64` is 128÷64→64 and Rust cannot emit it for `u128 % u128` — it cannot
  prove the divisor fits — and at ~30–90 cycles it would still lose to two
  multiplies. What it *would* change is the prime size: 61-bit primes would then
  cost the same per multiply and need one fewer CRT prime, hence one fewer
  evaluation pass.

  The scalar is where the whole design comes from — `Rational` runs a 128-bit
  gcd per operation and a prime field does not — and what the prime size then
  cost, on the `gj_tables` ladder:

```text
  mul, mul, add     AFrac<i128>   3512 ns
                    Rational       648 ns      5.4×
                    mod p          9.6 ns    366×

  gj_tables ladder, seconds     61-bit primes    31-bit + Barrett
  n = 10                             7.46              1.99      3.7×
  n = 11                            19.5               5.72      3.4×
  n = 12                            59.9              17.96      3.3×
```

  The 61-bit column is not a naive baseline — it is the same engine after the
  `reconstruction_matrix` fix in the row above. Both tables were in
  `gjmod.rs`'s module doc; the rustdoc keeps the reason (`u128 %` is a function
  call on aarch64 and `u64 %` is not) and this file keeps the numbers.

  ⚠️ One recorded number was simply false and is corrected in place: a doc
  comment claimed the removed 128-iteration `mulmod` cost "2.5×" of the engine.
  It cost 7.46 → 6.85 s at n = 10, inside the noise. That was the third wrong
  guess of the session, and the reason to drop the check is that it verifies the
  arithmetic rather than the mathematics.

  It is a second engine, not a tweak, and the exact one stays as its
  cross-check — the `qtkostka.rs` "three routes" standard. It also loses one
  free law (that the denominators collapse, which is how [DF] is currently
  enforced), so three checks have to do that job instead: the `b = 0` class
  algebra, the `b = 1` double coset algebra, and the held-back prime. Rational
  reconstruction returns a *spurious small rational* rather than failing when a
  value exceeds the bound, so a large bound is not by itself evidence it was
  large enough.
- **⚠️ The next rung on the [GJ] ladder is not n = 15.** Pushing degree was the
  plan and it is now the wrong plan: see the coverage table above. What the
  literature has not settled is the **statistic** `wt_λ` itself — Matchings-Jack
  asserts one function of λ and a matching δ simultaneously produces the right
  polynomial for *every* pair (μ,ν), and that is a constraint-satisfaction
  problem rather than a sign check. The search space grows fast: matchings on
  2n points number `(2n−1)!!`, so 2.0×10⁶ at n = 8 and 6.5×10⁸ at n = 10 —
  exhaustive is comfortable through 9, painful at 10.

  The result worth having is not yes/no but the **shape of the solution
  space**. Is the statistic unique for a given λ? If it is rigid, the definition
  can be read off the data and then proved. If there is enormous slack, hunting
  for "the" canonical statistic is the wrong framing. Likewise: for which λ does
  La Croix's θ work verbatim, where does it need patching, and is the patch
  systematic?

  Everything needed is already here — `double_coset_table` enumerates and types
  matchings, and both engines produce the target polynomials on the open triples.
  ⚠️ Not started, and it competes directly with LLT for the next slot. The
  framing above comes from a secondary summary rather than from the community
  itself, so it is a hypothesis about what would be worth reading, not a
  reported view.
- **Stanley's table now reaches degree 16** (111804 triples, 74 s, all in
  ℕ[α]); Sage cannot do the single degree-12 product `J[3,2,1]²`. Degree 18 is
  the next rung.
- **`jack_p_branching` for a single coefficient.** E1 fills the whole row
  whatever you asked for; E2 computes one μ. The `lr_coeff` lesson says the
  peeling order matters. Candidate, not plan.
- **The "Lassalle recurrences" named above** are Lassalle–Schlosser Pieri
  inversion — an explicit expansion by inverting Pieri. Left unexplored rather
  than rejected: E1 is already enumeration-free with unit-cost denominators, so
  no advantage was identified, but no numbers argue against it either.
- Shifted / interpolation Jack (Knop–Sahi's other family, with its own open
  positivity conjecture on structure constants) is the natural v2, and the
  reason `AFrac` is its own module rather than buried in `jack.rs`.
- Nonsymmetric `E_η` via [KS] Thm 4.6, which is what Cherednik-operator methods
  are built on.
- `AFrac` is now the **fourth** factored fraction field after `Frac`,
  `bh::Rat` and `deltaop::Ratio`, and the only canonical one. The
  `FactoredFrac<A>` refactor the Macdonald spec argued for now has a fourth
  witness and its cleanest instantiation.

## Negative result: Jack is coefficient-bound, not container-bound

Checked while sweeping the tree for the defect that gave `m → s` 2.35x and the
Pieri directions 2.0–4.0x (`docs/record/transitions.md`) — a hot map keyed on a
freshly allocated `Partition`. `jack.rs` has one, `HashMap<Partition, AFrac<C>>`
in the coefficient loop, so it looked like the same shape.

It is not. Sampling `profile_jack 12`: `AFrac::reduce_at` **15.5%**,
`AFrac::lift` 5.4%, `AFrac as Ring` 5.3%, `compiler_builtins` integer division
5.0%, allocator ~18%, and `jack_p_lb` — the algorithm — 4.7%. The cost is
rational-function arithmetic in the coefficient ring, and the gcd reduction
inside it above all. The partition keying does not appear.

So the transitions fix does not port here, and the place to look, if there is
one, is `AFrac` — not the container. Recorded so the pattern match is not made
a second time from the code alone.

## The forward direction, for a whole element (2026-08-21)

`jack_p_to_monomial`, `jack_q_to_monomial` and `jack_j_to_monomial` in
`src/jack.rs` expand a `P`, `Q` or `J` element — the coefficient map the
inverse expansions return — back into the monomial basis. One `jack_p` per
shape *present*, grouped so nothing rebuilds a degree it does not need, where
`jack_table` is the whole-degree route. They exist because the Python surface
now names a shape in its own basis and expands on request
(`docs/record/python-and-sage-interop.md`).

`every_monomial_comes_back_as_itself` closes the composite the other way
round from `every_jack_polynomial_comes_back_as_itself`: `m_μ → P → m_μ` at
every shape through degree 7, in all three normalizations. A table inverted
correctly in one direction only passes the older test and fails this one.

## Jack has no plethysm, and the obstruction is `AFrac`, not the engine (2026-08-24)

The other nine operations `Sym` has reached `Param` over all four coefficient
rings (`docs/plans/element-model.md`). Plethysm reached three of them —
`plethysm_qt`, `plethysm_macdonald`, `plethysm_ht`, over new `Plethystic` impls
for `Frac` and `Ratio` — and stopped at ℚ(α).

`Plethystic::frobenius` is the nth plethystic Frobenius, the map that raises
every variable of the coefficient ring. Over `ℚ[q,t]` and its two fraction
types that is `q^a t^b ↦ q^{an} t^{bn}`, which carries `1 − qᵃtᵇ` to
`1 − q^{an}t^{bn}` and `qᵃ − tᵇ` to `q^{an} − t^{bn}`: both denominator
families are closed, so the impls are a key remapping. Over ℚ(α) it is
α ↦ α^n, and that is **not** closed on `AFrac`, whose denominator is an integer
times a product of primitive *linear* forms `uα + v`. At `n = 2` a factor
`α + 1` becomes `α² + 1`, irreducible over ℚ.

The values are real rather than an artifact of the encoding. Asked of Sage with
`SAGE_DISABLE_SYMFN=1`, over `SymmetricFunctions(FractionField(QQ['alpha']))`:

    p[2].plethysm((1/(alpha+1))*p[1])   : (1/(alpha^2+1))*p[2]
    JackP[2].plethysm(JackP[2])         : (alpha^2+1) divides three of the
                                          five coefficients' denominators

So Jack plethysm needs a general ℚ(α) — a univariate rational function ring
with a polynomial gcd — which is a new coefficient ring, not an impl on an
existing one. `AFrac` cannot be widened to it without giving up the factored
form the Jack engine's speed rests on: `AFrac::reduce_at` is already 15.5% of
`profile_jack 12` (the negative result above), and a gcd over ℚ[α] is the more
expensive operation.

`Param.plethysm` therefore refuses ℚ(α) by name and points at `.at()`, which
specializes α and hands back a `Sym` where the integer route applies. That is
the one remaining gap in the ten, and it is recorded in the plan's deferred
section rather than left as a silent absence.
