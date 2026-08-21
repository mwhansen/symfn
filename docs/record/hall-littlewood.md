# Past the classical core: Hall–Littlewood and Kostka–Foulkes

The first subsystem past the classical target: a `(q,t)` coefficient ring,
the charge statistic as an independent reference, the Morris recursion
read out of Symmetrica's C, and the Kostka–Foulkes polynomials that fall
out of the resulting transition.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

## The `(q,t)` coefficient ring — done

`src/qt.rs`. `QtPoly<C>` is a sparse bivariate polynomial in q and t over any
`Ring`, so one type covers ℤ[q,t] (`QtPoly<i64>`, where Kostka–Foulkes
coefficients live), ℚ[q,t] (`QtPoly<Rational>`) and an unbounded variant
(`QtPoly<BigInt>`). Sparse because Hall–Littlewood and Macdonald expansions are:
a Kostka–Foulkes polynomial has few terms relative to its degree.

It exists only because the dividing paths were re-bounded on `QAlgebra` rather
than `Field` — ℚ[q,t] is not a field, and z_μ⁻¹ is all they ever need.
`s → p → s` round-trips over it for every partition through degree 6.

**`Plethystic` is bounded on `QAlgebra`, not on the coefficient ring.**
`QtPoly<i64>` is a perfectly good ring for *holding* Hall–Littlewood
coefficients and cannot do plethysm, because plethysm routes through the
power-sum basis and carries z_μ⁻¹. The bound rejects `QtPoly<i64>` at compile
time instead of failing at runtime.

The Frobenius raises both variables, q^a t^b ↦ q^{an} t^{bn}. Checked as a ring
homomorphism, as the identity at n = 1, and against Sage on the case a single
monomial cannot distinguish — a *sum* in the coefficient:

```text
sage: s[2]((q+t)*s[1])   ->   q*t*s[1,1] + (q^2+q*t+t^2)*s[2]
```

A frobenius scaling the whole polynomial by q^n t^n, rather than raising each
variable, passes every monomial test and fails that one.

## Charge and Kostka–Foulkes — the reference

`src/charge.rs`. `K_{λμ}(t) = Σ_{T ∈ SSYT(λ,μ)} t^{charge(T)}`, by enumerating
tableaux. **918 polynomials through degree 8 agree with Sage exactly.**

This is deliberately the slow one — the role `NaiveLr` plays for
Littlewood–Richardson. It is ~4x per degree (0.0013s for the whole degree-8
table, 1.08s at degree 13), so it is a usable oracle to about degree 14 and
nothing more.

The value of this route is that it shares no code with the fast one.
Hall–Littlewood will come from a recursion over skewing and straightening, so
agreement between them is evidence rather than tautology. Reading Symmetrica is
what established that: its `hall_littlewood` does not use charge at all, so the
two routes are genuinely disjoint. Had charge been on the fast path — as the
earlier plan in this file assumed — this check would have been worth much less.

Two conventions were pinned by hand-computation before any code ran, since the
literature differs by reversal: the reading word is **bottom row to top row,
left to right**, and `index(i) = index(i−1) + 1` when i lies to the *right* of
i−1. For λ = (2,1), μ = (1,1,1) that gives charges 2 and 1, i.e. `t² + t`, which
is Sage's answer. Peeling a general word into standard subwords **wraps** at the
left end; omitting the wrap gives increasing rather than standard subwords and
is the easy mistake.

Beyond Sage, `K_{λμ}(1)` is checked against the ordinary Kostka number for every
pair through degree 7, and the support against dominance — both tying it to
machinery already tested by other means.

## What Symmetrica does for Hall–Littlewood (read before building)

`Symmetrica_2.0/sr.c` and `rest.c`, read on the same footing as `tms.c`/`muir.c`
earlier — it is public domain, see `NOTICE.md`. Three findings, and the first
changes the plan.

**1. The default `hall_littlewood` is a recursion on λ, not a charge sum.**
Morris (1963). Peel the last part, recurse on λ⁻, then rebuild:

```text
  base  ℓ(λ) = 1            ->  a single term, coefficient t⁰
  step  HL(λ) = Σ over ν in HL(λ⁻) of
           ν extended by a new row λ_last                     (the i = 0 term)
         + Σ_{i=1..|λ⁻|} t^i · (h_i^⊥ ν) extended by λ_last + i
```

then straighten. The inner operation is `part_part_skewschur(ν, (i))` — skewing
by a **one-row partition**, which is exactly `h_i^⊥`. We already have that as a
native path (Pieri run backwards, horizontal-strip removal, no Littlewood–
Richardson), measured at 2.2–3.3x the LR route.

**2. `reorder_hall_littlewood` is β-number straightening**, the same rules we
already implement twice: a negative part kills the term, two adjacent entries
differing by exactly −1 kill it, and a descending pair is fixed by negating the
coefficient and swapping with a ±1 adjustment.

**3. `charge_word` is there, but Hall–Littlewood does not use it**, and
Symmetrica has no Kostka–Foulkes entry point at all. The statistic itself:

```text
  standard word (content all 1s):
      index(1) = 0;  index(i) = index(i−1) + 1 if i lies right of i−1, else index(i−1)
      charge = Σ index
  general word:
      peel standard subwords by a cyclic right-to-left scan (find 1, then 2, …,
      wrapping at the start), remove, recurse; charge = Σ over subwords
```

⚠️ **This corrects the plan given earlier in this file.** I wrote that charge was
"the one genuinely new combinatorial primitive" and put Kostka–Foulkes first,
with Hall–Littlewood built on top. Symmetrica's default path does the opposite
and uses no charge: HL comes from a recursion over machinery we already have and
have already made fast, and K_{λμ}(t) then falls out of the s ↔ P transition.
Charge is better used as an *independent check* on the result than as the way to
compute it — which is the more valuable role anyway, since the two routes would
share no code.

Revised order: **Hall–Littlewood by the Morris recursion over `QtPoly`**, using
`SkewBy<Homogeneous>` and the existing straightening; then Kostka–Foulkes from
the transition; then charge as an independent check; then Macdonald, which
needs a fraction field ℚ(q,t) over `QtPoly` and degenerates to HL at q = 0.
The t = 0 and t = 1 specializations remain the first tests, and
`QtPoly::eval` exists for them.

## Hall–Littlewood: built, and where the time actually went

`src/hl.rs` returns `Q'_λ = Σ_μ K_{μλ}(t) s_μ`. The specialization that pins
*which* Hall–Littlewood this is turned out to be t = 1, not t = 0: both `P` and
`Q'` give `s_λ` at t = 0, while `Q'_λ(x;1) = h_λ` and `P_λ(x;1) = m_λ`. Only
the t = 1 test distinguishes them, and it is the one worth writing first.

Three oracles, of decreasing independence:

| oracle | scope | what it can catch |
| --- | --- | --- |
| `charge::kostka_foulkes` | degrees 1–8, in-crate | shares no code — real evidence |
| Sage `hall_littlewood().Qp()` | 507 shapes, 17977 coefficients, ≤ deg 14 | a convention both of ours got wrong |
| Symmetrica's C `hall_littlewood` | every λ ≤ deg 15 | a port error — same algorithm, so weakest |

**Against Symmetrica end to end: 2.2–3.0×, growing with degree**, both sides
charged for building the Sage object (`scripts/bench_hl.py`).

### The predicted optimization was the wrong one

The plan said the win would be **sharing the recursion's suffixes across a
degree**, the pattern that took `kostka_table` from 0.39× to 2.5× and the
character table from 0.65× to 3.8×. It was implemented, it works, and it is
worth **1.1–1.2×** — not nothing, but nowhere near the earlier sweeps. The
reason is structural: for a single λ the recursion already calls itself only
once per part, so the top level dominates and there is little below it to share.

Sampling the release binary (`sample`, 5s at degree 21) found the real cost, and
it was not combinatorial at all:

| | before | after |
| --- | --- | --- |
| `Schur::add_term` | 660 | 105 |
| malloc/free | ~750 | ~250 |
| `QtPoly::add_assign` / `mul` | 97 | 464 → (see below) |
| **`remove_horizontal` / `removals`** | **44** | 71 |

44 samples out of ~1800 in the actual strip enumeration. Everything else was
temporaries. Two changes, each pointed at directly by a profile:

1. **Stop materializing `h_i^⊥ prev` for each i.** Written the way the recursion
   reads, each i built a whole `Schur<QtPoly>` map, walked it once and dropped
   it. Removing a horizontal strip of *any* size from ν is a single interlacing
   walk with `i = |ν| − |μ|` falling out at the leaf, so one pass replaces
   `|λ⁻| + 1` passes and nothing is built in between.
2. **Stop allocating a map to multiply by `t^i`.** `shift_t` is a rename of the
   exponents; adding the shifted terms straight into the destination slot
   removed the 464-sample `add_assign`.

Together: **2.15×** on the Rust path at degree 17 (0.1639s → 0.0764s), with
byte-identical output to the dump Sage had already verified.

The predicted optimization paid 1.15× and the unpredicted one paid 2.15×; the
profile named the second, and the plan had not mentioned it. The Jacobi–Trudi
row order is the same failure from the other direction — there a 200×
regression hid because the benchmark shapes were too uniform to expose it.

### The `QtPoly` representation, and a premise that was wrong twice

With those fixed, `QtPoly::add_term` rose to the top of the profile (219
samples). The obvious move — a `BTreeMap` of a handful of terms should be a
sorted `Vec` — was implemented and measured **7% slower** (0.0755s → 0.0812s at
degree 17, three runs each way, confirmed against a stashed build).

Measuring the premise instead of assuming it explained both halves. Over the
148 448 coefficients of degree 18:

```text
  terms per coefficient   mean 16.3   median 12   p90 35   p99 69   max 99
  density over the support                                          0.999
```

Two facts, each contradicting something I had assumed:

* **Not a handful.** At 16 terms with a tail to 99, per-term binary-search-plus-
  memmove insertion is worse than a B-tree rebalance. That is the 7%.
* **Essentially gapless.** Kostka–Foulkes polynomials have no interior holes, so
  the terms arrive as a *sorted run*, and shifting by t^i preserves that order.

So the representation was fine and the insertion pattern was not. `add_shifted`
merges the two sorted sequences in one pass instead of inserting term by term:
**1.40×** on top (0.0755s → 0.0538s), and the `Vec` now wins clearly.

Totals for Hall–Littlewood, all with byte-identical output to the dump Sage
verified:

| | degree 17, Rust | vs Symmetrica, degree 15, end to end |
| --- | --- | --- |
| as first written | 0.1639s | 2.58× |
| no intermediate `h_i^⊥` map, no `shift_t` temporary | 0.0764s | 2.77× |
| `Vec` + merging accumulation | **0.0538s** | **3.36×** |

**3.05× on the Rust path**, from three changes, none of which was the one the
plan predicted. The prediction — sharing the recursion across a degree — is real
but worth 1.06–1.2×, and it is now the *smallest* of the four effects measured.

## Kostka–Foulkes, from the transition

`src/kf.rs`. `Q'_μ = Σ_λ K_{λμ}(t) s_λ`, so the polynomials *are* the
coefficients Hall–Littlewood already produces and the module is the entry point
that says so. Symmetrica has none — `hall_littlewood` is the only way to reach
these from it, and the transition has to be read off by hand.

Checked against Sage's `kfpoly` on **every (λ, μ) pair through degree 11** —
3136 pairs at degree 11 alone — including the zero pairs, since a transition
that is right on its support and wrong about where the support *is* would pass a
nonzero-only comparison. Also against `kostka_foulkes_by_charge` (degrees 1–8,
no shared code), against `kostka_table` at t = 1, and against the identity at
t = 0.

| degree | pairs | symfn | Sage `kfpoly` | |
| --- | --- | --- | --- | --- |
| 9 | 900 | 0.0084s | 0.2604s | 31× |
| 10 | 1764 | 0.0220s | 0.8199s | 37× |
| 11 | 3136 | 0.0554s | 2.2102s | **40×** |

Those ask per pair, the way Sage is asked, so both sides answer the same
question. But **one `Q'_μ` is an entire column**, so per-pair is the wrong unit:
`kostka_foulkes_column` produces all 3136 values of degree 11 in **0.0025s**, a
further 22× on our own per-pair number and 884× on Sage's. This is the same
effect as the `kostka_table` result — the cost is in answering p(n)² independent
queries, not in the mathematics — and it is why `kf` exposes the column and the
table, not only the single value.

Still unexploited: the 0.999 density means a coefficient could be a dense `Vec<C>`
with a base offset, making accumulation O(1) index arithmetic. That is a
bigger change to `QtPoly` and would need to stay honest about the bivariate
case, where nothing guarantees density in q. Worth revisiting when Macdonald
gives a second workload to measure against — both premises above were formed on
the Hall–Littlewood workload alone, and both were wrong.

**Macdonald has since provided that second workload, and it was worth waiting
for.** Its numerators hold hundreds of terms where Hall–Littlewood's hold a
dozen, and the two agree that the sorted `Vec` is right — but only once
`QtPoly::mul` stopped accumulating with `add_term`. See below.

## Offline oracle fixture

`check_hl.py`, `check_hl_p.py` and `check_kf.py` are wider but only run when
someone has Sage and remembers. `gen_sage_oracle.sage` emits 30 `Q'` and 30 `P`
expansions in the Schur basis through degree 6, plus all 210 Kostka–Foulkes
pairs through degree 6 — **92 of them zero**. A transition right on its support
and wrong about where the support *is* passes a nonzero-only comparison; the 92
zeros are what catch it. `cargo test` checks all of it with no Sage installed.

`Q'`, `P` and `K(t)` are carried separately rather than derived from one
another. They differ by `b_λ(t)` and by the transition matrix, all three answer
to "the Hall–Littlewood polynomial", and deriving two from one would make a
convention swap self-consistent instead of visible.

## The inverse direction: `s → P` and `s → Q'`

Until 2026-08-21 every parametric family ran one way — `P_λ`, `Q'_λ`, the
Macdonald and Jack polynomials, all returned *expanded* in a classical basis —
and nothing took an element back into the family's own basis. The research
questions these families serve mostly run the other way: a symmetric function
defined some other way, expanded in `P` or `Q'` or `H̃` to ask whether its
coefficients are positive. A Sage user had the inverse the whole time, because
Sage's triangularity machinery solves for it from symfn's forward expansion;
the direct-Python user did not. The Hall–Littlewood pair is built first
because it needs no new linear algebra, and it fixes the design the other
families will follow ([parametric-basis-inverses.md](../plans/parametric-basis-inverses.md)).

**No solve is needed.** Both transitions are already in the crate as data and
only had to be re-indexed:

* `s_μ = Σ_λ K_{μλ}(t) P_λ`, so the `s → P` coefficients are the
  Kostka–Foulkes matrix read by rows — `kostka_foulkes_table(n)`, the table
  `hall_littlewood_p_table` inverts to get `P → s`.
* `⟨P_λ, Q'_μ⟩ = δ_{λμ}` under the Hall inner product, so the coefficient of
  `Q'_λ` in `s_ν` is the coefficient of `s_ν` in `P_λ`: the `s → Q'` matrix is
  the transpose of `hall_littlewood_p_table(n)`.

Both stay in `ℤ[t]`. `schur_to_hall_littlewood_p` and
`schur_to_hall_littlewood_qp` in `src/hl.rs` take a `Schur<QtPoly<C>>` of any
mixture of degrees, apply each degree's table to its own terms, and return a
`BTreeMap<Partition, QtPoly<C>>` — a plain map, because the crate has no
`P`-basis type and a `Schur` holding `P`-coefficients would be the basis
confusion the types exist to prevent. Cost is one `hall_littlewood_table(n)`
per degree present (plus the inversion for `Q'`), which is what a single
`hall_littlewood_p(λ)` already costs.

**The boundary** takes and returns the `[(lambda, [(t_exponent,
coefficient)])]` encoding `hall_littlewood` emits, so a `P` or `Q'` answer
feeds straight back in; it escalates from guarded `i128` to `BigInt` like the
forward entry points. **The convenience layer** adds `hl.to_P` and `hl.to_Qp`,
which accept a `Sym` in the Schur basis, a `Param` in `t`, or the rows, scale
rational coefficients through the boundary the way `Sym.to` does, and return a
`Param` tagged `HLP` or `HLQp`. The tags are the names Sage prints, so
`hl.to_P(s([2]))` reads `t*HLP[1,1] + HLP[2]` here and `t*HLP[1, 1] + HLP[2]`
there. `Param.at` refuses the parametric tags — a `Sym` carries only the six
classical codes, and a parametric basis has no meaning at a fixed `t` — and
`ParamBasis` in `_types.py` is the `Literal` that makes the new codes a type
error before a `ValueError` (python.md, the typing section).

**Pinned by** three tests in `src/hl.rs`: the inverse of each forward
expansion is the unit element for every λ through degree 8
(`the_inverse_expansions_undo_the_forward_ones`); the hand values
`s_2 = P_2 + t·P_11` and `s_11 = Q'_11 − t·Q'_2`
(`s2_in_p_and_s11_in_qp_are_the_hand_values`), which are the smallest values
that tell the two directions apart — the `t` lands on the smaller shape with a
plus sign in one and on the larger with a minus in the other; and
`ℤ[t]`-linearity across a mixed-degree argument. Both hand values were checked
against Sage 10.9 (`SAGE_DISABLE_SYMFN=1`, the sage-dev environment):
`HLP(s[2])`, `HLQp(s[1,1])`, and further `HLP(s[2,1]) = (t^2+t)*HLP[1,1,1] +
HLP[2,1]`, `HLQp(s[2,1]) = HLQp[2,1] - t*HLQp[3]`,
`HLP(s[3]) = t^3*HLP[1,1,1] + t*HLP[2,1] + HLP[3]`. On the Python side
`scripts/check_convenience.py` holds `to_P(P(λ))` and `to_Qp(Qp(λ))` to the
unit, the wrappers to the contract rows, and the coefficient of `P_λ` in `s_μ`
to `kostka_foulkes(μ, λ)` through a different entry point, for every λ through
degree 6. Not timed: each call costs what the forward table costs, which the
sections above already measure.

### Open tail

* **The other families.** `s → H̃` landed 2026-08-21
  ([macdonald-operators.md](macdonald-operators.md), "The expansion on its
  own"). Macdonald (`m → P`, `m → Q`, `s → J`),
  Jack (`m → P`, `m → Q`, `m → J`) and the choice of tags, encodings and
  oracles for each are laid out in
  [parametric-basis-inverses.md](../plans/parametric-basis-inverses.md).
  The `s → J` direction already exists as a whole-degree table
  (`schur_in_macdonald_j`; [qt-kostka.md](qt-kostka.md), "The inverse of
  `J → s` is a projection, not a solve") and only wants the element-wise
  form and the wrapper.
* **A fixture for the inverse direction.** The Sage values above were checked
  by hand in one session and are pinned as Rust constants, not read from
  `tests/fixtures/`. The round trip against the fixtured forward expansions
  covers the mathematics; a committed `HLP(s_μ)` / `HLQp(s_μ)` table through
  degree 6 would cover the orientation independently, and belongs in
  `gen_sage_oracle.sage` when that file is next regenerated.

