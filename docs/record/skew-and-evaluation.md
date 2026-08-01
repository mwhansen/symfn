# Skewing by an arbitrary symmetric function, and evaluation

Two pieces of surface that are new rather than faster: the adjoint of
multiplication under the Hall inner product, generic over the basis the
skewing function is written in; and the bridge from formal symmetric
functions back to concrete numbers.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

## Skewing by an arbitrary symmetric function

`g^⊥`, the adjoint of multiplication by g under the Hall inner product:
⟨g^⊥ f, h⟩ = ⟨f, g·h⟩. `skew_schur` was only the case g = s_μ, f = s_λ — which
is why the classical notation for it is a quotient, s_μ^⊥ s_λ = s_{λ/μ}.

`SkewBy<C, G>` is generic over **the basis g is written in**, and that is the
design rather than a convenience. g ↦ g^⊥ is linear, so any g could be expanded
into Schur and handed to Littlewood–Richardson — but three bases have adjoints
with direct rules, each a Pieri or Murnaghan–Nakayama rule read backwards:

| g in basis | g^⊥ on s_λ | machinery |
|---|---|---|
| h | remove a horizontal strip | Pieri |
| e | remove a vertical strip | dual Pieri |
| p | remove a rim hook, signed by height | M–N |
| s, m, f | Σ_μ d_μ s_{λ/μ} | LR |

`e` is implemented as ω ∘ h^⊥ ∘ ω rather than as a second strip enumerator: ω
is an isometry with ω(h_r) = e_r, so ⟨e_r^⊥ s_λ, s_ν⟩ = ⟨h_r^⊥ s_{λ'}, s_{ν'}⟩.
An identity, and it costs one transpose per term against a duplicate enumerator
with its own separate bugs. `p` reuses `character::border_strips` — the same
β-number bit tricks that make the character table fast.

**The p case is the one that most repays the native path, and not only for
speed: it needs no division.** Expanding p_μ into Schur requires 1/z_μ, so the
LR route forces ℚ. Rim-hook removal stays in ℤ, so `Schur<i64>` can be skewed
by a power sum — something the generic route could not have offered at all.

Native path vs the same g expanded into Schur, interleaved in one process:

| g | shapes | ratio |
|---|---|---|
| p | λ of 4–7 rows, degree 42–52 | **43–268x** |
| h | same | 2.2–3.3x |
| e | *tall* λ (10–12 rows) | 3.0–7.2x |

⚠️ On battery. Same-process single-threaded ratios are far less power-sensitive
than the parallel LR numbers were, but they are not re-measured on AC.

The `e` row needed the shape family changed, and the first attempt is worth
recording as a measurement error rather than a result. On the *wide* shapes used
for h and p it read 0.8–1.6x — because a vertical strip needs a distinct row per
cell, so on a 4-row λ the answer is nearly empty and both routes finish in under
a millisecond. The ratio was measuring harness noise on a trivial answer, not
the algorithms. Tall shapes make the operation non-trivial and it behaves like h.

Verified against Sage (`scripts/check_skew.py`): **32448 checks through degree
8, zero mismatches**. Every (λ, μ, basis) is put to symfn twice — once in its
own basis, once expanded into Schur — so agreeing with Sage says the answer is
right and agreeing with each other says the fast path is a shortcut and not a
different operation. In-crate, the defining adjointness ⟨g^⊥ f, h⟩ = ⟨f, g·h⟩ is
checked for every triple through degree 6; that test mentions no algorithm at
all, and its two sides share no code.

## Evaluation at an alphabet, and the principal specializations

`src/eval.rs`. Everything else in the crate computes *with* symmetric functions
as formal objects; this is the bridge back to concrete numbers. Two different
things live there and the distinction is the design:

**Evaluation at an arbitrary alphabet**, generic over `Ring`, one algorithm per
basis rather than "convert, then evaluate". p, e, h are products of one-row
generators and cost a linear DP each. Schur uses the **branching rule**: a
tableau is a chain ∅ = ν⁰ ⊆ … ⊆ νⁿ = λ of horizontal strips, so sweeping
variable by variable with a frontier of *shapes* collapses every tableau sharing
a prefix into one number. Cost is the shapes inside λ, not the tableaux — the
same chain DP as `kostka.rs`, carrying ring elements instead of counts.

**The bialternant is deliberately absent.** s_λ = a_{λ+δ}/a_δ is the textbook
formula and would be an O(n³) determinant, but it needs *division* — so it is
not generic over `Ring` — and it is 0/0 whenever two x_i coincide. The branching
rule is slower on generic input and always right. The oracle script exercises
exactly that hole: one of its alphabets is `[2,2,2,0,5]`.

**Closed forms** for the two special alphabets, which enumerate nothing:

| | formula |
|---|---|
| `dimension(λ)` = f^λ | \|λ\|! / ∏ h(u) |
| `principal_specialization(λ,n)` = s_λ(1ⁿ) | ∏ (n + c(u)) / h(u) |
| `principal_specialization_q(λ,n)` | q^{n(λ)} ∏ (1−q^{n+c(u)}) / (1−q^{h(u)}) |

The q-analogue returns a coefficient vector. Neither product divides the other
cell-by-cell, so the quotient is taken once at the end; both have constant term
1, which makes it a truncated power-series inversion — no leading-coefficient
case analysis, and exact in ℤ because the quotient is known in advance to be a
polynomial.

`dimension` and `principal_specialization` interleave their divisions with their
multiplications rather than forming the factorial first, and that is not a
micro-optimisation: for the staircase λ = (10,9,…,1), f^λ has **35 digits** and
fits `u128`, while 55! has **74**. Forming the numerator first would overflow by
thirty-five orders of magnitude on an answer that is comfortably representable.
Both return `Option`, `None` on genuine overflow.

Verified against Sage (`scripts/check_eval.py`): **1440 checks through degree 9,
zero mismatches** — evaluation against `expand(n)` substituted, both
specializations against `principal_specialization`, and f^λ against Sage's own
tableaux count. In-crate, the five bases are checked to agree with each other at
a shared alphabet through degree 6, which is the check that catches a wrong
recurrence in any one of them: they share no code, so they can only agree by all
being right.

No timings are quoted. This is new surface rather than a faster route to
something we already had, and Symmetrica has no equivalent entry point.

## Offline oracle fixture

Both halves of this file's subject now carry oracle evidence `cargo test`
re-establishes on its own: 30 coproducts, 30 antipodes and 30 counits through
degree 6, and 210 principal specializations in each of the plain and graded
forms, over alphabets from 0 to 6, plus 30 dimensions.

Three things the sweep is shaped by:

- **Laws were the only evidence here, and laws are convention-blind.** The
  Hopf axioms and the LR-route agreement hold under any self-consistent
  normalization, so what the oracle adds is the values.
- **The antipode is the pin.** `S(s_λ) = (−1)^{|λ|} s_{λ'}` couples a sign to a
  conjugation, and both halves are silent on a self-conjugate shape of even
  size — so the sweep runs over whole degrees, and `S(s_{(2,1)}) = −s_{(2,1)}`
  is the case where the shape is fixed by conjugation and only the sign is left
  to be wrong.
- **The coproduct's tensor orientation is not pinned, and cannot be.** Sym is
  cocommutative: all 30 coproducts in the fixture are invariant under swapping
  the two sides, so no value distinguishes this convention from its transpose.
  That is a property of the object, not a hole in the sweep — recorded so the
  sweep is not later mistaken for evidence it cannot carry.

`f^λ` is checked against `StandardTableaux(λ).cardinality()`, a direct count,
rather than Sage's hook-length formula: symfn computes `f^λ` by hooks, and an
oracle using the same formula would check the arithmetic and nothing else. The
alphabet size sweeps from 0, where `s_∅(1⁰) = 1` and every other shape gives 0
— a convention over a well-posed question, so it is swept rather than skipped
(V6), and 30-plus of the cases are that vanishing.
