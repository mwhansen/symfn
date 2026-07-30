# Research gaps in symmetric-function software

*Survey compiled 2026-07-28. Question: what capabilities related to symmetric
functions would be compelling for research and do not currently exist in any
package (Sage, Magma, Symmetrica, Stembridge's SF, lrcalc)?*

The walls below were **measured** on this machine against SageMath 10.9, not
recalled. The literature section was gathered by search on the same date.

Companion document: `sage-acceleration.md` profiles *why* those walls are where
they are, and asks the separate question of how to speed up capabilities that
already exist. This document is about capabilities that do not exist at all.

---

## 1. Measured walls in the incumbent

SageMath 10.9, single runs, per-item `SIGALRM` timeouts of 90–120s.

| operation | input | result |
|---|---|---|
| `internal_product` (Kronecker) | `s[10,7,5]·s[9,7,6]`, n=22 | 2.92s → 721 terms |
| | `s[14,10,8]·s[13,10,9]`, n=32 | **>90s timeout** |
| | `s[12,9,7,4]·s[11,10,7,4]`, n=32 | **>90s timeout** |
| | `s[20,15,10]·s[18,15,12]`, n=45 | **>90s timeout** |
| single coefficient `g(λ,μ,ν)` | 3-row, n=45 | **>90s timeout** |
| plethysm | `h₅[h₅]`, deg 25 | 2.13s → 245 |
| | `h₄[h₇]`, deg 28 | 7.66s → 173 |
| | `h₅[h₆]`, deg 30 | 12.03s → 492 |
| | `h₆[h₆]`, deg 36 | 110.87s → 2002 |
| | `h₄[h₁₀]`, deg 40 | **>120s timeout** |
| chromatic symmetric function | path `P₁₆` | 0.34s |
| | path `P₂₀` | 3.48s |
| | path `P₂₄` | 58.48s |
| | random tree, 26 vertices | **>120s timeout** |
| Jack `P → s` | n=9 | 5.59s |
| | n=15 | **>90s timeout** |
| Macdonald `H̃ → s` | n=15 | 2.10s |
| | n=18 | 10.55–15.41s |
| | n=21 | 66.82s |
| ∇ (nabla) | `∇e₆` | 0.21s |
| | `∇e₈` | 2.13s |
| stable Kronecker via `st` basis | `st[3,1]·st[3,1]` | 0.31s |
| | `st[6,4]·st[6,4]` | **>90s timeout** |
| | `st[8,5]·st[7,4]` | **>90s timeout** |
| e-/Schur-positivity certification | — | **no API exists** |

⚠️ Single runs on a laptop with coarse timeouts. These are order-of-magnitude
walls, not benchmarks — the same caution `ROADMAP.md` applies to its own tables
applies here.

**The headline is the `st` row.** The Orellana–Zabrocki irreducible-character
basis is *the* modern tool for reduced Kronecker coefficients — its outer-product
structure constants literally are the stable Kronecker coefficients — and Sage's
implementation dies on two 2-row partitions of 10.

**The second headline is that no single-coefficient path exists anywhere.**
Every package computes the entire product to read one number. This is precisely
the `lr_coeff` defect already found and fixed for Littlewood–Richardson
(ROADMAP: "peeling off the *larger* factor instead", a 23x improvement no
product benchmark could surface), still unexercised for Kronecker and plethysm.

### What *does* exist, to be fair

- Sage has `internal_product`, `reduced_kronecker_product`, `plethysm`,
  Macdonald (`P Q J H Ht S`), Jack, LLT, `nabla`, `theta_qt`, the
  Orellana–Zabrocki `st` basis, `chromatic_symmetric_function` and
  `chromatic_quasisymmetric_function`. Coverage is broad; **speed is the gap**,
  and in a few places the modern operators are simply absent.
- `lrcalc` — LR and quantum LR, fast, narrow.
- `barvikron` — Christandl–Doran–Walter lattice-point algorithm for Kronecker
  coefficients of *bounded height*, polynomial time. Python prototype, niche,
  effectively unmaintained.
- Baldoni–Vergne–Walter distribute Maple code for bounded-length Kronecker.
- Stembridge's `SF` — Kronecker naively; very general user-defined bases.

---

## 2. The gaps worth building, ranked

### 2.1 A Kronecker engine with a single-coefficient query

Three algorithmic routes are under-exploited, and the third is unimplemented:

1. **Christandl–Doran–Walter lattice-point counting** for bounded height —
   polynomial time, only `barvikron` implements it.
2. **Vector partition functions** (Baldoni–Vergne–Walter) — Maple only.
3. **Panova 2025** polynomial-time classical algorithms for bounded parameters
   ([arXiv:2502.20253](https://arxiv.org/abs/2502.20253)) — as far as could be
   determined, **nobody has implemented these**.

Architecturally this is the best fit for existing machinery. Bounded-height
Kronecker is a lattice-point count in a box on a hyperplane — the exact regime
where `three_row.rs` ("counting per output instead of enumerating chains") won,
including the same trap: the first attempt there was 3.3x *slower* because of a
`HashMap` keyed on the state tuple, fixed by generation-stamped dense tables.

Downstream demand: Saxl / tensor-square conjecture, geometric complexity theory,
Kronecker positivity (NP-hard to decide, so instances matter).

### 2.2 Reduced / stable Kronecker as a first-class ring

Murnaghan stability, the rate of stabilization, restriction coefficients. Making
`st`-basis multiplication fast is the whole task — it is a *product* in a
nonhomogeneous basis rather than a character-table blowup, so it sidesteps the
p(n)×p(n) memory ceiling documented in `ROADMAP.md` (1.1 GB at n=32).

Note: computing reduced Kronecker coefficients is #P-hard and deciding their
positivity is NP-hard, so the target is a fast engine for reachable instances,
not a general algorithm.

### 2.3 Chromatic symmetric functions at scale, with positivity certification

**Highest leverage item in this document.**

- Stanley's tree-isomorphism conjecture is verified to 29 vertices only by
  bespoke C code. Sage times out on a random tree at **26**.
- The strongest argument for building this landed in 2026: a paper killed
  Stanley's claw-free Schur-positivity conjecture *and* Monical's saturated-
  Newton-polytope conjecture by **finding counterexamples**
  ([arXiv:2607.21508](https://arxiv.org/abs/2607.21508)). Search over this space
  is where the results are.
- No package ships `is_e_positive` / `is_schur_positive` with a witness, let
  alone a counterexample-search driver.

Algorithms to draw on: deletion–contraction, the p-expansion over connected
partitions of the edge set, and the modular relations of Orellana–Scott and
Gebhard–Sagan.

### 2.4 The Macdonald operator algebra: Δ_f, Δ'_f, Θ_f

**Cheapest item on the list.** Sage has `nabla` and the classical plethystic
`theta_qt` only. The Delta and Theta operators that the entire diagonal-harmonics
community works with are absent from every package and live in personal Maple
files. The valley Delta conjecture is still open.

*Both versions of the Delta conjecture are now checked to n = 9 as whole
symmetric functions (`src/dyck.rs`), the valley one being open. The wall is no
longer the operator but the `(n+1)^{n−1}` path enumeration.*

*Built: `src/deltaop.rs`, 2026-07-29. ∇, ∇^r, Δ_f, Δ'_f, Π, Π⁻¹ and Θ_f, all
against `spec-macdonald-operators.md`. `∇e_13` takes 16s where Sage takes 5m38s
(~21×, mains-to-mains); Δ, Δ' and Θ exist nowhere else to compare against. The
constraint on the valley Delta conjecture is now the labelled-Dyck-path
enumeration, not the operator — see `ROADMAP.md`.*

*Specified in `spec-macdonald-operators.md` (2026-07-29), where the claim above
is re-measured and holds: `theta_qt` and `scalar_qt` are the near-misses and are
both genuinely different operators, and there is no Δ, Δ', Θ or Π anywhere. Two
findings from that document belong here. First, **Sage's ∇ is slow at the change
of basis, not at the Macdonald polynomials** — `Ht(e[10])` is 20.2s of the 29.2s
`∇e_10` costs, while one H̃_μ reaches the Schur basis in 0.056s; the crate's
whole degree-10 table takes 0.294s. Second, the expansion into `{H̃_μ}` needs
**no matrix inversion**: `H̃` is orthogonal for the star scalar product, so the
coefficient is `⟨F,H̃_μ⟩_*/w_μ`.*

This is algebra over ℚ(q,t), which `QtPoly` plus a fraction field nearly supports
already, and it needs no new enumeration engine. The one real requirement is
**honest plethystic substitution at formal alphabets** — `X(1-q)/(1-t)`,
`1/(1-t)`, and friends — which is the single most error-prone corner of Sage's
symmetric-function library.

Relevant existing note in `ROADMAP.md`: the `(q,t)` Frobenius already raises both
variables (`q^a t^b ↦ q^{an} t^{bn}`) and was checked against the case a single
monomial cannot distinguish. That is the right foundation.

### 2.5 Witnesses, not just numbers

Explicit LR tableaux, charge and cocharge, crystal isomorphisms, jeu-de-taquin
paths, RSK. Every package returns integers; researchers testing conjectural
bijections need the objects. `NaiveLr` already enumerates the tableaux — the gap
is that nothing surfaces them.

Note: *charge* is already on the critical path for Kostka–Foulkes (ROADMAP:
"the one genuinely new combinatorial primitive"), so this partly comes for free.

### 2.6 Jack polynomials at scale

Sage times out at n=15. The Goulden–Jackson matchings-Jack conjecture and the
b-conjecture are open and computationally starved. Routes: Knop–Sahi, the
Lassalle recurrences, or the Laplace–Beltrami eigenoperator, rather than
Gram–Schmidt.

*Specified in `spec-jack.md` (2026-07-29). Walls re-measured there and
sharper than the row above: Sage's whole-degree tables (P→m, J→m, J→p) all
die at **n = 12**, the single shape λ=(n) carrying essentially the entire
degree's cost, and the Stanley-conjecture product `J[3,2,1]·J[3,2,1]` is
already >120s. Three shared-nothing routes were verified against Sage before
entering the spec (`scripts/spec_jack_verify.py`, "all formulas verified"):
the [MOPS] Laplace–Beltrami moving-box recursion — the engine, honoring the
"rather than Gram–Schmidt" warning — the branching formula via the per-atom
`q = t^α` limit of the Macdonald ψ, and the Knop–Sahi tableau formula. The
decisive measurement (`scripts/spec_jack_swell.py`): the whole calculus
lives over factored integer-linear atoms `uα+v` — canonical, exact lcm,
integer root-test cancellation — and the recursion shows **no denominator
swell at all** (peak numerator degree 6 at n = 12). Status of the targets:
Lassalle's conjecture is now a theorem (Ben Dali–Dołęga 2305.07966);
matchings-Jack has polynomiality (Dołęga–Féray) and integrality (Ben Dali
2203.14879) settled, positivity open; the b-conjecture likewise open.*

*Built 2026-07-29: `src/afrac.rs`, `src/jack.rs`, `src/gj.rs`. All three routes
implemented and agreeing (E1 ≡ E2 to n = 8, E3 ≡ both to n = 5), 1565 values
checked against Sage with no mismatches. Measured battery-to-battery in one
session with isolated Sage processes: whole-degree `P → m` at n = 11 is
**0.0186 s against Sage's 155.1 s (8340×)**, and Sage prices `P → m`, `J → m`,
`J → p` and the norms table identically — it is the `P → m` transition that
costs, not the unit. Whole tables run to **n = 26** (520 s); n = 16, the degree
Stembridge ships as precomputed archives, is 0.73 s. Stanley's full
degree-12 table (9317 triples, including the `J[3,2,1]²` case Sage cannot do)
is 46.6 s and lies entirely in ℕ[α]. The **Goulden–Jackson `c` and `h` tables
are complete through n = 14** (2026-07-30, 186 s, 2045553 + 1121377
coefficients, all in ℕ[b]) — computed nowhere else — by a second engine
(`src/gjmod.rs`) that runs the whole pipeline at numeric α over prime fields and
reconstructs the answers, sharing no arithmetic with the exact ℚ(α) one and
required to agree with it wherever both are affordable. Its reusable half is
`src/modular.rs`. The transcription is pinned at **both** known specializations,
neither touching a Jack polynomial: b = 0 against the `S_n` class algebra from
characters, b = 1 against the double coset algebra of `(S_2n, H_n)` by counting
matchings.

⚠️ **The degree is the wrong dial, and this is the correction that matters most
here.** Positivity and integrality are theorems, the degree bound is
characterized (Promyslov), and Ben Dali's marginal sums are already known
b-positive with a matchings interpretation — so a counterexample must hide inside
a marginal sum with its siblings cancelling it, and the bar for a bulk sign check
as a remark worth making is n ≥ 25. Most of the 2.0M coefficients also fall in
already-proved cases; `matchings_jack_coverage` now separates them, and at n = 8
only 83.5% of live triples are open at all. What is *not* fenced in is the
**statistic `wt_λ`** — one function of λ and a matching that must produce the
right polynomial for every (μ,ν) simultaneously — where the payload is the
rigidity of the solution space rather than a yes/no, and where n ≤ 9 is enough
because `(2n−1)!!` is 2.0×10⁶ at n = 8. Not started. ⚠️ This framing comes from
a secondary summary, not from the community, so it is a hypothesis about what
would be worth reading. ⚠️ Two methodological corrections are recorded in
`spec-jack.md` §6.1: the walls above cannot be reproduced in a single Sage
process (it memoizes the transition matrices, and `SIGALRM` corrupts them
mid-build), and the spec's "a failed cancellation costs one dot product" was
wrong — it cost two heap allocations, worth 1.6× once removed.*

### 2.7 Cylindric / affine and quantum LR, k-Schur, Catalan functions

Underserved computationally; see the Blasiak–Haiman–Morse–Pun–Seelinger results
below, several of which are *raising-operator formulas* and therefore directly
implementable.

---

## 3. Recent results worth building on

| result | reference | why it matters here |
|---|---|---|
| **Hikita's proof of Stanley–Stembridge** | [arXiv:2410.12758](https://arxiv.org/abs/2410.12758) (Oct 2024, rev. Dec 2025) | Proves e-positivity for (3+1)-free graphs via a *probabilistic* interpretation of the e-coefficients of the chromatic quasisymmetric function of unit interval graphs. Those probabilities are a brand-new computable object with no implementation anywhere. |
| **Claimed proof of Saxl's conjecture** | [arXiv:2512.15035](https://arxiv.org/pdf/2512.15035) (Dec 2025) | Staircase-minimality theorem + Ikenmeyer + Bessenrodt–Bowman–Sutton lifting. ⚠️ Preprint, treat as unverified. Makes Kronecker positivity testing newly interesting either way. |
| **Panova, classical vs. quantum multiplicities** | [arXiv:2502.20253](https://arxiv.org/abs/2502.20253) (2025) | Polynomial-time *classical* algorithms for Kronecker and plethysm in many bounded-parameter cases, refuting claimed quantum speedups. Unimplemented algorithms sitting in a paper. |
| **BHMPS: LLT in the Schiffmann algebra** | Crelle 811 (2024) 93–133 | Explicit raising-operator formula for ∇ applied to any LLT polynomial. ⚠️ The LLT side of this now exists (`src/llt.rs`, `docs/spec-llt.md`); the Catalanimal route for `∇` of a general LLT is deferred there as spec §3.7, and this row is the open half. |
| **BHMPS: Demazure crystals and Schur positivity of Catalan functions** | Invent. Math. 236 (2024) 483–547 | |
| **BHMPS: raising-operator formula for Macdonald polynomials** | Forum Math. Sigma (2025) | Plausibly beats Sage's Macdonald path; directly implementable. |
| **Nonsymmetric shuffle theorem** | [arXiv:2509.24040](https://arxiv.org/pdf/2509.24040) (Sep 2025) | |
| **Extended / compositional Delta theorems** | D'Adderio–Mellit; Blasiak–Haiman–Morse–Pun–Seelinger | Rise version is a theorem by two independent routes; **valley version still open**. Neither extended nor compositional implies the other. |
| **Claw-free CSFs are not Schur positive** | [arXiv:2607.21508](https://arxiv.org/abs/2607.21508) (2026) | Counterexamples to Stanley (claw-free Schur positivity) and Monical (SNP). Found by search — the existence proof for §2.3. |
| **Two stability theorems on plethysms** | [arXiv:2505.06104](https://arxiv.org/pdf/2505.06104) | |
| **Geometric / generating-function approach to plethysm** | [arXiv:2511.02649](https://arxiv.org/pdf/2511.02649) | Gutiérrez, Orellana, Saliola, Schilling, Zabrocki. |
| **Schur positivity of ∇m_μ** | [arXiv:2607.00940](https://arxiv.org/pdf/2607.00940) (2026) | |
| **Orellana–Zabrocki character basis** | [arXiv:1605.06672](https://arxiv.org/abs/1605.06672); Hopf structure in *Alg. Comb.* | The basis whose outer-product structure constants are the stable Kronecker coefficients. Underpins §2.2. |

Foulkes' conjecture, for calibration: known for a ≤ 4, verified for a = 5
(Cheung–Ikenmeyer–Mkrtchyan, symmetrizing tableaux). Open since 1949. The Sage
plethysm wall measured above sits essentially exactly at the frontier.

---

## 4. Recommendation

1. **§2.4 — Δ/Θ operators plus plethystic calculus over ℚ(q,t).** Weeks, not
   months, given `QtPoly`. Hands an active community a tool that exists nowhere.
   *Specified 2026-07-29; every formula verified against Sage by
   `scripts/verify_deltaop_formulas.py`. The one engineering risk it identifies
   is denominator swell in `Σ_μ (num_μ / w_μ) H̃_μ`, with `macop.rs` as the
   standing precedent for that going badly.*
2. **§2.3 — CSF at scale with positivity search.** Where a fast engine most
   plausibly produces a *result* rather than a convenience.
3. **§2.1 — Kronecker single-coefficient engine.** The deepest and best fit for
   the frontier/counting machinery, but the largest build.

All three are outside the "Beyond the core (deferred, but intended)" list in
`ROADMAP.md` — that list is Symmetrica's remaining scope (modular/projective
representations, Schubert, Hecke algebras). This document is about scope that
*no* package covers, which is a different and more interesting target.
