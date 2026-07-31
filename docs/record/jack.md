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

## The coefficient field is the whole story

Every scalar in the Jack calculus — hooks, ψ-ratios, eigenvalue differences,
norms — is a ratio of **integer-linear forms `uα + v`**. Normalise the atoms to
primitive and three things become true that are false for `Frac`'s `1 − qᵃtᵇ`:
distinct atoms are irreducible and pairwise coprime (so the factored form is
canonical, where `Frac`'s `PartialEq` must cross-multiply); atom-wise max
multiplicity is the **exact** lcm rather than a common multiple; and a failed
cancellation is refuted on its first step by Gauss's lemma.

⚠️ Primitivity is load-bearing, not cosmetic. `E(κ)−E(λ)` is genuinely
non-primitive — κ = (2,2), λ = (1,1,1,1) gives `2α+4` — and dividing ℤ[α] by a
non-primitive form leaves ℤ[α].

`AFrac<C>` is a **ℚ-algebra for any `C`**, including `i128`, which is not one:
dividing by an integer multiplies the scalar denominator and needs nothing from
`C`. That is what lets the engine run over `AFrac<i128>` and still be handed to
`s → p`, which asks for `QAlgebra` because it divides by `z_ν`.

## Measured: 3000–8000×, against a target of 200×

⚠️ **Battery to battery, one session, isolated processes.** The
pre-implementation measurements had `P → m` at n = 11 as 84.5 s on mains, and
this machine drifts ~1.8× — confirmed here rather than assumed: the same cell
re-measured on battery is 155.1 s, a ratio of 1.84.

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
sharper version of the spec's finding that the single shape λ = (n) is the whole
degree. It is not the *unit* that costs: the `P → m` transition is, and `J`, the
p-expansion and even the closed-form norm all route through it. Our four differ
by four orders of magnitude, because they are actually different computations.

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
cache filled at `i128` and handed out to other widths would launder exactly
that away — the `memo::bold_p` hazard. Hoisting the loop has no correctness
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
forever). It was caught only because `bench_jack` prints the scale width as a
column. Print the shape of the data, not just the time.

## The Goulden–Jackson pipeline

`src/gj.rs` computes `c^λ_{μν}(b)` and `h^λ_{μν}(b)` — the Matchings-Jack and
b-conjecture coefficients. **No package computes either table.**
ℚ[b]-polynomiality is a theorem (Dołęga–Féray) and `c`'s integrality is a
theorem (Ben Dali), so both are enforced; **positivity is open for both and is
only observed**, with any negative coefficient reported as a finding rather
than debugged away — the valley-Delta posture.

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
counterexample has to hide inside a marginal sum with its siblings cancelling
it. The bar for "verified through n = N" as a remark worth making is n ≥ 25.

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
Promyslov any one of the three — but for a **variation involving labelled
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

## Next

- **Push the [GJ] tables past n = 10** — the deliverable, and what the engine
  exists for. The cost is `phi_slice`: `Σ_θ` of a rank-1 tensor over `p(n)³`
  entries, so `p(n)⁴` coefficient operations, measured growing ~4.6×/degree
  (faster than `p(n)⁴`'s 3.8×, because the coefficients grow too).

  `examples/probe_gj.rs` measures the four things that decide the fix, and
  **three of them killed the design they were testing**:

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

  Two things make this safe rather than a gamble. **There are no poles for
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

  ⚠️ One recorded number was simply false and is corrected in place: a doc
  comment claimed the removed 128-iteration `mulmod` cost "2.5×" of the engine.
  It cost 7.46 → 6.85 s at n = 10, inside the noise. That was the third wrong
  guess of the session, and the reason to drop the check is that it verifies the
  arithmetic rather than the mathematics.

  It is a second engine, not a tweak, and the exact one stays as its
  cross-check — the `qtkostka.rs` "three routes" standard. It also loses one
  free law (that the denominators collapse, which is how [DF] is currently
  enforced), so three checks become load-bearing: the `b = 0` class algebra, the
  `b = 1` double coset algebra, and the held-back prime. Rational reconstruction
  returns a *spurious small rational* rather than failing when a value exceeds
  the bound, so a large bound is not by itself evidence it was large enough.
- **⚠️ The next rung on the [GJ] ladder is not n = 15.** Pushing degree was the
  plan and it is now the wrong plan: see the coverage table above. What the
  literature has not fenced in is the **statistic** `wt_λ` itself — Matchings-Jack
  asserts one function of λ and a matching δ simultaneously produces the right
  polynomial for *every* pair (μ,ν), and that is a constraint-satisfaction
  problem rather than a sign check. It dies fast: matchings on 2n points number
  `(2n−1)!!`, so 2.0×10⁶ at n = 8 and 6.5×10⁸ at n = 10 — exhaustive is
  comfortable through 9, painful at 10.

  The payload is not yes/no but the **shape of the solution space**. Is the
  statistic unique for a given λ? If it is rigid, the definition can be read off
  the data and then proved. If there is enormous slack, hunting for "the"
  canonical statistic is the wrong framing. Likewise: for which λ does La Croix's
  θ work verbatim, where does it need patching, and is the patch systematic?

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
- Nonsymmetric `E_η` via [KS] Thm 4.6 — the door to Cherednik-operator methods.
- `AFrac` is now the **fourth** factored fraction field after `Frac`,
  `bh::Rat` and `deltaop::Ratio`, and the only canonical one. The
  `FactoredFrac<A>` refactor the Macdonald spec argued for now has a fourth
  witness and its cleanest instantiation.
