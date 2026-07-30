# Specification: LLT polynomials (`llt`)

> **Status, 2026-07-30: implemented, profiled, and tuned.** `src/llt.rs`
> (+ `k_core`/`k_quotient`/β-number helpers in `partition.rs`), 30 in-crate
> tests plus 4 for the abacus primitives, all green. The Sage conformance pair
> is `examples/llt_dump.rs` → `scripts/check_llt.py` (**1517 comparisons, 0
> failures** at a degree-5 dump). Benchmarks are `examples/bench_llt.rs`,
> profiling workloads `examples/profile_llt.rs`; §6 carries measured numbers
> instead of guesses and §6.1 the profiling pass that made them **2.2–4.1×**
> better than the first working version.
>
> Implementing the spec falsified two of its own claims — §1.3(a)'s floor
> witness and §1.3(g)/§5.12's location of the multi-term KL entries — both
> corrected in place and marked ⚠️. §7's two answerable open questions moved:
> Conj 6.4 is swept to n = 14 (§7.2) and the R2/R3 crossover is measured
> (§7.4). **Python bindings built** (§4, 16 functions; `check_bindings.py` has 29
> boundary checks, 0 failures). Not built: the [BHMPS] route (§3.7, deferred by
> design).

*Written 2026-07-29. Implements the `docs/record/dyck-paths.md` "Next" item under the Dyck
ladder work — "Vertical-strip LLT polynomials as a first-class object, which
subsumes the decomposition above and is a `docs/research-gaps.md` item in its own
right" — and closes the two obstructions recorded there (§3.6). LLT appears
in `docs/research-gaps.md` twice: the Sage-coverage list (line 66) and the
[BHMPS] row (line 235, "Explicit raising-operator formula for ∇ applied to
any LLT polynomial").*

Sources, all fetched and read today, 2026-07-29 — the PDFs, not
recollections:

| tag | source | used for |
|---|---|---|
| **[LLT]** | Lascoux, Leclerc, Thibon, *Ribbon tableaux, Hall–Littlewood functions, quantum affine algebras and unipotent varieties*, [arXiv:q-alg/9512031](https://arxiv.org/abs/q-alg/9512031) (J. Math. Phys. 38, 1997) | spin (15), cospin (24)–(25), G̃ (26), H̃ (27), H (28), φ-image (29); Thm 6.6 (k ≥ ℓ-bound ⇒ Hall–Littlewood); Ex 6.8(i)/(ii) fixtures; Conj 6.4 (open); §7's straightening rule — **misprinted, see §1.3(g)** |
| **[LT]** | Leclerc, Thibon, *Littlewood–Richardson coefficients and Kazhdan–Lusztig polynomials*, [arXiv:math/9809122](https://arxiv.org/abs/math/9809122) | the spin-generating G (their (43)); Ex 4.1 fixture; Thm 4.2 (c^λ_μ(−q⁻¹) = parabolic affine KL, Varagnolo–Vasserot); Prop 5.11 (49)–(51) three-term straightening; Lemma 6.5 (beta-strips) |
| **[KMS]** | Kashiwara, Miwa, Stern, *Decomposition of q-deformed Fock spaces*, [arXiv:q-alg/9508006](https://arxiv.org/abs/q-alg/9508006) (Selecta Math. 1, 1996) | the **normative** straightening rules (43) and (45) — the version the verify script implements |
| **[HHL]** | Haglund, Haiman, Loehr, *A combinatorial formula for Macdonald polynomials*, [arXiv:math/0409538](https://arxiv.org/abs/math/0409538) (JAMS 18, 2005) | Def 3.2 (tuple model, attacking inversions); the standardization identity (their (82)); the Macdonald decomposition H̃_μ = Σ_D q^{−a} t^{maj} G_ν(μ,D) |
| **[BHMPS]** | Blasiak, Haiman, Morse, Pun, Seelinger, *LLT polynomials in the Schiffmann algebra*, [arXiv:2112.07063](https://arxiv.org/abs/2112.07063) (Crelle 811, 2024) | Def 2.2.1 (their tuple conventions); the coproduct law; eq (4): ω∇^m G_ν as a Catalanimal — the v2 route (§3.7) |
| **[CGKM]** | Corteel, Gitlin, Keating, Meza, *A vertex model for LLT polynomials*, [arXiv:2012.02376](https://arxiv.org/abs/2012.02376) | coinversion-LLT conventions (content j−i), Cauchy identity — recorded as convention datapoints, not an engine route |
| **[GH]** | Grojnowski, Haiman, *Affine Hecke algebras and positivity of LLT and Macdonald polynomials* (preprint, 2007 — **never published**) | positivity of G_ν for arbitrary skew tuples; cited as the status of the theorem, nothing implemented from it |
| **[CM]** | Carlsson, Mellit, *A proof of the shuffle conjecture*, [arXiv:1508.06239](https://arxiv.org/abs/1508.06239) (JAMS 31, 2018) | Prop 3.5 (chromatic bridge); the Dyck path algebra the [DA] recursion lives in |
| **[DA]** | D'Adderio, *e-positivity of vertical strip LLT polynomials*, [arXiv:1906.02633](https://arxiv.org/abs/1906.02633) (JCTA 172, 2020) | the theorem G_ν[X; q+1] is e-positive; Remark 2.2 (the tuple **reversal**, §1.3(e)); Ex 5.6 fixture; the path-algebra recursion (5.1)–(5.4) |
| **[AS]** | Alexandersson, Sulzgruber, *A combinatorial expansion of vertical-strip LLT polynomials in the basis of elementary symmetric functions*, [arXiv:2004.09198](https://arxiv.org/abs/2004.09198) (Adv. Math. 400, 2022) | the orientation/e-expansion formula (hrv blocks); Schröder-path picture (diagonal steps = strict edges); Ex 6.1 fixture; Problem 6.20 (open) |
| **[AP]** | Alexandersson, Panova, *LLT polynomials, chromatic quasisymmetric functions and graphs with cycles*, [arXiv:1705.10353](https://arxiv.org/abs/1705.10353) (Discrete Math. 341, 2018) | Def 16 (coloring model); Lemma 47 (chromatic bridge, their form); Conj 25 (open); the area-word convention clash with [CM], §1.3(f) |

Every formula below was **verified numerically before being written down**,
by `scripts/spec_llt_verify.py`, committed alongside this document — 26
PASS lines, ~19 s, ends `all formulas verified`. The walls script is
`scripts/spec_llt_walls.py` (round 2 of 2; round 1's exploratory sweep
lived in the session scratchpad and its surviving numbers are folded into
§2). **Sage is an oracle, never a source** — SageMath 10.9 (2026-05-04), this machine, today; its
`sage/combinat/ribbon_tableau.py` and LLT code were not read and must not
be, the same clean-room posture as `docs/record/jack-spec.md` and
`docs/record/macdonald-operators-spec.md`. Everything Sage-shaped below (basis
dictionaries, oracle warts) was established by black-box probing.

---

## 1. What the thing is

One family, two combinatorial presentations, and a dictionary between them
that is nowhere written down in one place correctly — pinning that
dictionary operationally is a large fraction of this spec's value (§1.3).

### 1.1 The ribbon model ([LLT] §6)

Fix a level `k ≥ 1`. A **k-ribbon** is a connected skew shape of k cells
containing no 2×2 square; its **spin** is `s = (h−1)/2 ∈ ½ℕ` where `h` is
its height ([LLT] (15)). A k-ribbon tableau of shape λ is a chain peeling
λ by **horizontal k-ribbon strips**; `s(R)` = total spin, `s*(λ)` = max
spin over tableaux of shape λ, **cospin** `s̃(R) = s*(λ) − s(R) ∈ ℤ`
([LLT] (24)–(25); integrality verified on every shape the script touches).

```text
  G̃^(k)_λ(x; q) = Σ_R q^{s̃(R)} x^{w(R)}      cospin generating   [LLT] (26)
  H̃^(k)_μ       = G̃^(k)_{kμ}                                     [LLT] (27)
  H^(k)_μ        = Σ_R q^{s(R)} x^{w(R)} = q^{s*} H̃^(k)_μ(x; 1/q) [LLT] (28)
  G_LT,λ(x; q)   = Σ_R q^{2s(R)} x^{w(R)} = q^{2s*} G̃_λ(x; q⁻²)   [LT] (43)
```

`kμ = (kμ_1, kμ_2, …)`. G̃ and H are symmetric functions (nontrivial —
[LLT] Thm 6.1); the script verifies symmetry implicitly by expanding in
the m-basis against Sage's Schur arithmetic throughout.

Existence of k-ribbon tableaux of shape λ ⟺ λ has empty k-core. Beta-set
mechanics: peeling a horizontal k-ribbon strip of weight m from λ ⟺ adding
k to m distinct beta numbers of the **conjugate** partition, with the spin
read off the crossing count ([LT] Lemma 6.5; implemented as `strips_down`
in the verify script and cross-checked against fixtures, §5).

### 1.2 The tuple model ([HHL] Def 3.2)

Input: a tuple `ν = (ν^(1), …, ν^(r))` of skew shapes, each with an
integer **content offset**. Cells have content `c = col − row + offset`
(classical, per component). The **reading order** sorts all cells by
(content, then component index, then row). Cells `u ∈ ν^(i)`, `v ∈ ν^(j)`
**attack** iff (a) same content and `i < j`, or (b) `c(v) = c(u) + 1` and
`j < i` — i.e. between distinct components only. For a semistandard
filling `T`, `inv(T)` = attacking pairs `(u, v)` with `u` earlier in
reading order and `T(u) > T(v)`:

```text
  G_ν(x; q) = Σ_{T ∈ SSYT(ν)} q^{inv(T)} x^T
```

This model is **not** spin-normalized: `min_T inv(T)` can be strictly
positive (§1.3(a)), and the tuple loses the absolute spin grading
entirely (§1.3(c)). `q = 1` gives `∏_i s_{ν^(i)}` (verified, §5.E).
Schur-positivity of every `G_ν` for arbitrary skew tuples is [GH] — a
2007 preprint that was never published; for tuples of partitions it is
[LT]/Varagnolo–Vasserot KL-positivity. Treat "positive for all skew
tuples" as morally settled but cite-carefully.

**Standardization** ([HHL] (82)): `G_ν = Σ_{S ∈ SYT(ν)} q^{inv(S)}
Q_{n,D(S)}` where `Q_{n,D}` is the fundamental quasisymmetric function
and `D(S)` the descent set read along **content reading order** — ties in
content broken by the fixed reading order. `[x^μ] Q_{n,D} = 1` iff
`D ⊆ {partial sums of μ}`. Verified against direct SSYT enumeration on
every tuple in the script (§5.A). This identity is what dissolves the
first `docs/record/dyck-paths.md` obstruction (§3.6).

### 1.3 The convention minefield, measured

Every entry below cost a wrong first guess either in this session or in
the literature itself. Each is pinned by a PASS line.

- **(a) The min-inv floor.** The k-quotient dictionary is
  `G_{quot_k(λ), offsets 0} = q^{min inv} · G̃^(k)_λ` — *direct*, cospin to
  cospin, after dividing the tuple side by `q^{min inv}`. The floor is
  **forced**: λ = (2,2,2), k = 2 has 2-quotient `((1),(1,1))` and
  `min inv = 1`, and no offset choice removes it. Verified k = 2, 3 over all
  empty-core λ with |λ| ∈ {k, 2k, 3k}.
  ⚠️ **Corrected 2026-07-30** (`src/llt.rs::the_min_inv_floor_is_forced`):
  this entry previously named λ = (2,2) as the witness, which is wrong — its
  2-quotient is `((1),(1))`, whose `G` is `m₂ + (1+q)m₁₁`, floor **0**. (The
  convention gate in §5.1 is that same tuple and depends on the floor being
  zero there, so the two claims were never consistent.) The floored shapes at
  k = 2 through |λ| = 8 are (2,2,2), (4,2,2), (3,3,2), (2,2,2,2) and
  (2,2,2,1,1). The law itself is unaffected and still verified.
- **(b) The tilde collision.** [HHL]'s symbol "G̃_λ" is the
  **spin-flavored** function; [LLT]'s G̃ is **cospin**. Coding [HHL]'s
  remark-shaped `q^e G(1/q)` as the quotient dictionary fails; (a) is the
  truth. The script's section D records the correction in place.
- **(c) Tuples lose absolute spin.** λ = (1,1,1,1), k = 2: one ribbon
  tableau, `s* = 1`, but the 2-quotient tuple is two single cells with
  `max inv = min inv = 0`. No statistic of the tuple recovers `s*`; the
  spin-graded H is *shape data*. Any API that promises H from a bare
  tuple is lying — H is exposed on the partition-plus-level side only.
- **(d) Sage's dictionaries** (black-box, all verified):
  `llt(k, t=q).hspin()[μ] = Σ_R t^{s(R)}` on shape kμ (= [LLT] (28) H);
  `.hcospin()[μ] = H̃^(k)_μ`; `.cospin(Partition λ) = G̃^(k)_λ`;
  `.cospin(tuple) = q^{−min inv} G_ν` — Sage floors the tuple model.
  Sage's LLT **parameter is named t**, and the level-k family constructor
  needs `t` coercible into the base ring (`QQ['t']` works; a fraction
  field in a variable named q needs an explicit `t=q`).
- **(e) The shuffle-world reversal** ([DA] Remark 2.2). The
  vertical-strip tuple of a Dyck path (area sequence `a`, maximal
  +1-runs as columns, run starting at row r placed at offset `−a_r`,
  i.e. the cell of row i has content `−a_i`) must list components in
  **reverse row order** for [HHL] attacking inversions to coincide with
  [HRW] dinv pairs **pointwise** (verified per-path, per-filling-count,
  n ≤ 6). Row order is off by exactly the primary/secondary orientation.
- **(f) Two graphs per path, and two area words.** The cell-below-path
  rule (edge {j, i}, j < i, iff j ≥ i − a_i — the [CM]/[AS] picture)
  builds a *different* graph than the dinv-faithful one: a = (0,0) has
  no cells below the path but its rows form a primary dinv pair. The
  dinv-faithful decorated graph is: weak edges = dinv pairs (primary
  `a_x = a_y` oriented (x,y), secondary `a_x = a_y + 1` oriented (y,x)),
  strict edges = consecutive same-run rows. Additionally [AP] and [CM]
  write area words in opposite orders (P₃ is (1,1,0) in [AP], (0,1,1)
  in [CM]). Both graphs are real objects — §3.5 uses each where it is
  the right one.
- **(g) The straightening rules as printed in [LLT] §7 contain two
  misprints.** [KMS] (43)/(45) is normative: `u_l ∧ u_m = −u_m ∧ u_l`
  when `l ≡ m (mod k)`; otherwise, with `i = (m−l) mod k`,
  `u_l ∧ u_m = −q u_m∧u_l + (q²−1)(u_{m−i}∧u_{l+i} − q u_{m−k}∧u_{l+k}
  + q² u_{m−k−i}∧u_{l+k+i} + ⋯)`, offsets `i, k, k+i, 2k, 2k+i, …`,
  truncated at `2a < m − l`. [LT] Prop 5.11 (49)–(51) is the three-term
  equivalent. `q_KMS = q_LT^{−1}`; the operational dictionary the script
  pinned end-to-end (§5.M): Fock coefficients in the KMS variable v
  equal ribbon `c^λ_μ(q)` (the [LT] (43) grading) at **q = −v**. Single
  monomials cannot distinguish q = −v from its mirrors, so the sweep must
  reach a **multi-term** KL entry.
  ⚠️ **Corrected 2026-07-30**
  (`src/llt.rs::fock_straightening_is_the_ribbon_column_at_q_eq_minus_v`):
  those entries are **not at k = 2**. Every Schur coefficient of `G_LT` at
  k = 2 is a single monomial through |μ| = 10, so `spec_llt_verify.py`'s
  k = 2-only section M did *not* pin the variable dictionary — its range was
  consistent with the mirrors after all. The first multi-term entries live at
  **k = 3**: `s_{21}` in `G_LT,(3,3,3)` is `q² + q⁴`, and [LT] Ex 4.1's
  `q³ + q⁵` on `s_{211}` is the k = 3 shape (3,3,3,2,1) — the fixture was in
  the spec the whole time, at the wrong level. The Rust test sweeps k = 2
  (λ ⊢ ≤ 4) **and** k = 3 (λ ⊢ ≤ 3) and asserts a multi-term entry was
  actually seen, so the range cannot silently shrink below the one that pins
  the dictionary.

### 1.4 Specializations and degenerations, all verified

```text
  k = 1:      G̃^(1)_λ = s_λ                    (1-ribbons = cells, spin 0; definitional)
  q = 1:      G_ν(x; 1) = ∏_i s_{ν^(i)}        (verified vs Sage skew Schurs)
  k ≥ ℓ-bound: H^(k)_μ = Q'_μ(x; q)            ([LLT] Thm 6.6; verified H^(4)_{3211} = Sage HL Qp)
  ω-duality:  G_{ν'}(x; q) = q^{A(ν)} ω G_ν(x; 1/q),  A(ν) = #attacking pairs,
              ν' = (transpose each component, reverse order)   (verified, 4 tuples)
  coproduct:  G_ν(X+Y) = Σ G_{inner}(X) G_{outer/inner}(Y) with the [BHMPS]
              inv-offset bookkeeping                       (verified on ((2),(1)))
  Macdonald:  H̃_μ(x; q, t) = Σ_D q^{−a(D)} t^{maj(D)} G_{ν(μ,D)}(x; q)
              ([HHL]; verified against Sage Ht for ALL μ ⊢ n ≤ 5)
  shuffle:    ∇e_n = Σ_D t^{area(D)} G_D(x; q)  (verified n ≤ 5 vs deltaop-style Sage nabla)
  chromatic:  X_Γ(x; q) = (q−1)^{−n} G_Γ[x(q−1); q]   ([CM] Prop 3.5 / [AP] Lemma 47;
              verified vs Sage chromatic QSF, all unit-interval graphs n ≤ 5)
  e-expansion: Ĝ_P(x; q+1) = Σ_θ q^{asc(θ)} e_{λ(θ)} over orientations, blocks by
              highest-reachable-vertex ([AS]; verified on the fixtures)
```

Fixtures hit: [LLT] Ex 6.8(i) (G̃_{(3,3,3,2,1)}, k=3, in m **and** s),
Ex 6.8(ii) (H^(2)_{3211}), [LT] Ex 4.1 (q⁷s₃₁ + q⁵s₂₂ + (q⁵+q³)s₂₁₁ +
qs₁₁₁₁ in the (43) grading), [AS] Ex 6.1 (G_nndee = q²s₁₁₁ + qs₂₁),
[DA] Ex 5.6 (G_ndndee = qe₁e₃ + q(q−1)e₄).

### 1.5 Why anyone wants this

1. **LLT polynomials are the atoms of Macdonald theory.** [HHL] writes
   every H̃_μ as a positive sum of level-μ₁ LLTs; Macdonald positivity
   *is* LLT positivity. A fast H̃ engine that emits the LLT
   decomposition refines `qtkostka.rs`'s three routes with a fourth that
   is positively graded at every intermediate step.
2. **The ∇e_n refinement `dyck.rs` already wants.** `∇e_n = Σ_D t^{area}
   G_D` with every `G_D` Schur-positive: the by-path Schur refinement of
   the shuffle theorem. `docs/record/dyck-paths.md` records this as "the win left on
   the table" and names two obstructions; §3.6 shows both are resolved
   by the tuple model's standardization. No package emits this
   decomposition.
3. **Chromatic symmetric functions.** The unicellular bridge (§1.4) puts
   Shareshian–Wachs e-positivity (open) and [AP] Conj 25 (open) and [AS]
   Problem 6.20 (open) within sweep range of a fast engine — three named
   open problems that are explicitly computationally starved; [DA]'s
   theorem (G[X; q+1] is e-positive) plus [AS]'s formula make the
   e-expansion *certified output* rather than a conjecture check.
4. **Parabolic affine Kazhdan–Lusztig polynomials for free.** The Fock
   route (§3.4) computes columns of the `c^λ_μ` table — actual KL
   polynomials ([LT] Thm 4.2) — by exact straightening, no Hecke algebra
   in sight.
5. **The k-interpolation.** `H^(k)` walks from Schur (k=1) to
   Hall–Littlewood (k large, [LLT] Thm 6.6) through genuinely new
   symmetric functions; [LLT] Conj 6.4 (Schur-positivity of
   H^(k+1) − H^(k), still open per today's reading) has, to our
   knowledge, never been swept beyond small degrees.

### 1.6 Prior art, and what is actually ours

| who | what they have | where it stops |
|---|---|---|
| **Sage 10.9** (`llt.py` + `ribbon_tableau.py`, pure Python) | the only maintained implementation anywhere: `llt(k)` family with `hspin`/`hcospin` bases, `cospin`/`spin` on shapes and tuples, `cospin_polynomial` single coefficients | whole-degree walls at **n = 11/10/8** for k = 2/3/4 (§2.2); tuples die at n = 12; the worst shape *is* the degree (94–100%); a global cache makes fresh-parent measurement misleading (§2.1); parameter/coercion warts (§1.3(d)) |
| **Symmetrica** (the C library this repo wraps) | **no LLT at all** — it predates the subject's computational demand; nothing in the 2026-07-29 doc inventory | a capability gap, not a backend swap (same situation as Jack) |
| **[CGKM]** and the vertex-model line | transfer-matrix formulation, Cauchy identity | paper + research worksheets; no packaged engine; ⚠️ unmeasured here |
| **[BHMPS]** | ∇ of any LLT as a Catalanimal | theory; no implementation known to us |
| **everyone** | — | **nobody** ships: tuple-LLT as a documented API at scale, vertical-strip/unicellular LLTs as first-class objects, e-expansions via [AS], F-expansions via standardization, ∇e_n by-path Schur refinements, or KL columns via Fock straightening |

Honest accounting, per the standing house questions:

1. **The mathematics is entirely prior art** — every formula above has a
   1995–2022 citation. What no one has: an engine where the **dictionary
   web of §1.3–1.4 is pinned by machine-checked fixtures**, at speeds
   that reach the open-problem frontier. The convention soup is not a
   nuisance footnote; it is *the* reason no two implementations
   interoperate and half the remarks in the literature are off by a
   twist. 26 PASS lines say ours compose.
2. **The superpower is scale plus the web.** Sage's unit of work
   (whole-degree H^(k) tables) dies at n = 10–11 for k = 2; the single
   worst shape is the whole cost. The Rust engine targets the same
   ~100–8000× class of wins the Jack and Δ-operator builds measured, on
   an object whose incumbents are pure Python.
3. **The machinery fit is near-total** (§3.9): `QtPoly` is the
   coefficient ring as-is; `dyck.rs` already enumerates vertical-strip
   fillings with exactly the verified statistics (`d_row` *is* the
   attacking-inversion count); `charge.rs`/`hl.rs` are the beta-set
   strip template; `deltaop.rs` provides the ∇ cross-checks; `skew_lr`
   the q=1 gate. The genuinely new pieces are small: k-cores/quotients
   on `Partition`, a `SkewTuple` type, and the Fock straightening.
4. **The record's debt is paid** (§3.6): the recorded dinv-tie-breaking
   obstruction is dissolved by [HHL] (82) — verified, not argued.

---

## 2. Measured walls

SageMath 10.9, this machine, 2026-07-29, single runs, per-item `SIGALRM`
120 s. ⚠️ Order-of-magnitude walls, not benchmarks.

**Power provenance: every row in this section ran on battery (~96%,
discharging).** The Jack spec measured battery-at-100% rows within noise
of mains on this machine, and the Macdonald spec measured **1.8× drift at
lower charge** — so these tables are internally consistent but a
mains-vs-battery correction of up to ~2× must be assumed before comparing
against §6's eventual Rust numbers. §6 comparisons must be re-run
mains-to-mains (the standing rule; recorded here so nobody "discovers" a
2× speedup that is a power cable).

### 2.1 What Sage has, and the oracle warts

By introspection and black-box probing only: `Sym.llt(k, t=…)` with
bases `.hspin()`, `.hcospin()`; methods `.cospin(shape-or-tuple)`,
`.spin(...)`; `sage.combinat.ribbon_tableau.cospin_polynomial(shape,
weight, k)` for single coefficients. Warts, all hit today:

- **A global cache defeats fresh-parent timing.** With a *fresh*
  `SymmetricFunctions` parent per item, `HSp3[[6]]` costs 0.40 s but
  `HSp3[[7]]`, `[[8]]`, `[[9]]` then cost 0.00–0.01 s — some layer
  (ribbon-tableau/Kostka machinery) caches globally across parents. Yet
  `HSp3[[10]]` still times out: the cache helps its own reruns, not the
  frontier. Walls below are first-touch numbers; reruns in one process
  can look arbitrarily better and must not be quoted.
- **The t-parameter coercion trap** (§1.3(d)): `llt(3)` over
  `Frac(QQ[q])` fails (`t` not coercible) — construct over `QQ['t']` or
  pass `t=q` explicitly.
- **`cospin` rejects plain int lists** for shapes
  (`'float' object has no attribute 'floor'`) — wrap in `Partition`.
- **No whole-degree entry point**; per-element conversions, so the
  tables below are the per-degree sums a user actually pays.

### 2.2 Where Sage stops

Whole-degree spin tables (every μ ⊢ n, `HSp_k[μ] → s`), fresh parent per
(k, n):

```text
  k \ n      7        8        9        10       11
   2       0.65 s   2.22 s  10.48 s   51.16 s   >120 s
   3       2.11 s  13.01 s  90.01 s   >120 s      —
   4      18.72 s   >120 s     —         —        —

  worst shape μ = (n), share of whole degree:
   k=2: 94% (n=7) → 99% (n=10);  k=3: 98–100%;  k=4: 100%
```

**The finding that shapes the design: the one-row shape is the entire
degree**, exactly the Jack pattern. Touching `H^(k)_{(n)}` costs the whole
degree; growth is ~4.3–6.9× per degree (k = 2 steps: 3.4×, 4.7×, 4.9×;
k = 3: 6.2×, 6.9×).

Tuple LLTs (`llt(3).cospin(tuple)`, 3-tuples, fresh parents):

```text
  ((2,2),(2,1),(2))    n=9    1.14 s
  ((2,2),(2,2),(2))    n=10   2.80 s
  ((3,2),(2,2),(1))    n=10   3.14 s
  ((2,2),(2,2),(2,1))  n=11  19.69 s
  ((3,2),(2,2),(2,1))  n=12  >120 s
```

Single coefficients (`cospin_polynomial`, standard weight):

```text
  k=2  (6,6,4,2)      1⁹    0.15 s     k=3  (9,6,3)      1⁶   0.01 s
  k=2  (6,6,6,4,2)    1¹²   6.67 s     k=3  (9,9,6,3)    1⁹   0.30 s
  k=2  (8,8,6,4,2)    1¹⁴  >120 s      k=3  (12,9,6,3)   1¹⁰  0.86 s
```

So: the whole *family* — spin tables, tuples, single coefficients — walls
at n ≈ 10–14 in pure Python, and the objects §1.5 wants (∇e_n
decompositions at n = 10+, Conj 6.4 sweeps, chromatic corpora) sit on the
far side.

### 2.3 The same walls on mains, and the mains-to-mains comparison

Re-run 2026-07-30 with `pmset -g ps` reporting `AC Power`, same script, same
machine (`scripts/spec_llt_walls.py`). **Mains is ~2× faster than §2.2's
battery rows across the board** — the drift figure the provenance notes had been
assuming, now measured on this exact workload — and two items that timed out on
battery complete here, so the bounds become numbers. These are the rows every
§6 ratio is computed against; §2.2 is kept for provenance, not for arithmetic.

Whole-degree spin tables (every μ ⊢ n, `HSp_k[μ] → s`), against `llt_h_table`
(`i64`, mains, same session):

```text
  k   n      Sage      ours     ratio     Sage worst shape μ=(n)
  2   7      0.33s   0.0005s      660x     0.31s  (94% of degree)
  2   8      1.08s   0.0012s      900x     1.03s  (95%)
  2   9      4.84s   0.0031s     1560x     4.73s  (98%)
  2  10     24.16s   0.0078s     3100x    23.86s  (99%)
  2  11     >120s    0.0190s    >6300x        —
  3   7      1.11s   0.0009s     1230x     1.08s  (98%)
  3   8      6.77s   0.0028s     2420x     6.70s  (99%)
  3   9     47.34s   0.0080s     5920x    47.17s  (100%)
  3  10     >120s    0.0222s    >5400x        —
  4   7      8.29s   0.0017s     4880x     8.26s  (100%)
  4   8    101.45s   0.0053s    19100x   101.36s  (100%)
  4   9     >120s    0.0162s    >7400x        —
```

Tuple LLTs (`llt(k).cospin(tuple)` against `llt_g`):

```text
  ((2,2),(2,1),(2))     n=9     0.58s   0.0001s   ~5800x
  ((2,2),(2,2),(2))     n=10    1.45s   0.0002s   ~7250x
  ((3,2),(2,2),(1))     n=10    1.67s   0.0002s   ~8350x
  ((2,2),(2,2),(2,1))   n=11   10.36s   0.0010s  ~10400x
  ((3,2),(2,2),(2,1))   n=12   66.99s   0.0049s  ~13700x
```

Single coefficients (`cospin_polynomial`). Note the units differ in our favour
and the table says so: Sage returns **one** coefficient, `llt_gtilde` returns
the whole m-expansion and the requested coefficient is one of its 114 terms.

```text
  (8,8,6,4,2) wt 1^14 k=2   Sage 80.43s (1 coeff)   ours 0.0081s (114 coeffs)
```

**Where Sage is competitive, stated plainly.** `e[n].nabla()` is a real
implementation and the LLT route is *not* the way to beat it:

```text
  n    Sage e[n].nabla()   deltaop::nabla_e   llt::nabla_e_by_path
  9         4.64s               0.212s             0.956s
 10        14.64s               0.751s            19.278s
 11        41.10s               1.801s               —
```

`deltaop` is the route for the **total** (~20× over Sage, and that is
`docs/record/macdonald-operators-spec.md`'s number, not this module's). At n = 10 the
by-path route is *slower than Sage's total* — because it is not computing the
total: it emits all **16 796** per-path Schur-positive pieces, which is an
object Sage has no entry point for at any speed. Quoting 19.3 s against 14.6 s
as a loss, or against nothing as a win, would both be wrong; the honest
statement is that it is a different deliverable at comparable cost.

**The one number that is not a speed claim.** Sage's single shape `HSp3[[n]] → s`
costs 0.25 / 1.06 / 6.76 / 47.96 s at n = 6…9 and times out at 10, while
`llt_h(&(n), 3)` costs 10–70 **microseconds** across that whole range — a ratio
of 10⁵–10⁶ that says nothing about arithmetic throughput. It is the
cost-concentration finding of §2.2 turned inside out: the one-row shape is
94–100% of a degree for Sage and under 1% for the pruned walk, because chains to
a single row pass only through single rows. Our absolute times there are at
timer resolution, so the ratio should be read as "this shape stopped being the
bottleneck", not as a benchmark.

---

## 3. The algorithms

Three engine routes for the general object, one specialization family,
and one deliberately deferred route. All were prototyped in
`spec_llt_verify.py` and agree with each other and the oracle on
everything §5 lists. The pattern is `qtkostka.rs`'s: shared-nothing
routes in one module, cross-checking each other.

### 3.1 The coefficient ring: `QtPoly`, unchanged

LLT coefficients live in ℕ[q] (once normalized); the ∇-decomposition and
the Macdonald assembly want ℕ[q, t]. That is exactly `qt.rs`'s
`QtPoly<C>` (qt.rs:52) — sparse, sorted, already the crate's workhorse.
No new arithmetic, no fraction field anywhere in the core (contrast
Jack's `AFrac`): every algorithm below is subtraction-free in the main
path except Fock straightening, which is signed but exact over ℤ[v].
The min-inv floor (§1.3(a)) is bookkeeping on exponents, not arithmetic.

### 3.2 R1 — SYT + descent buckets: the reference engine and the F-expansion

[HHL] (82). Enumerate standard fillings of the tuple ν (linear
extensions, backtracking over the content reading order), bucket
`q^{inv(S)}` by descent set `D(S)`; then

```text
  G_ν = Σ_D ( Σ_{S: D(S)=D} q^{inv(S)} ) · Q_{n,D}
  [x^μ] G_ν = Σ_{D ⊆ partialsums(μ)} bucket(D)
```

Cost `#SYT(ν)` — for the n ≤ 12 tuple range this is the honest reference,
and the buckets *are* the fundamental-quasisymmetric expansion, an output
nobody ships. Monomial extraction from buckets is a subset-sum filter,
trivially cheap. This is the engine the verify script trusts everywhere
(`syt_buckets` ≡ direct SSYT enumeration: §5.A).

### 3.3 R2 — beta-set ribbon strips: the whole-degree H^(k) engine

The [LT] Lemma 6.5 mechanics, which is `charge.rs`/`hl.rs`'s shape: work
on beta numbers of the conjugate; a horizontal k-ribbon strip of weight m
= choose m beta numbers, add k to each (all results distinct), spin from
the crossing count. `H^(k)_μ` and whole-degree tables enumerate chains of
such moves down the weight word — a Morris-recursion-shaped walk
(hl.rs:57's `hall_littlewood` is the k-large degeneration, [LLT] Thm 6.6,
so the two modules must agree on a shared range: §5.C). Memoize on
(intermediate beta-set, remaining weight) exactly as `hl.rs`/`charge.rs`
do; the per-degree table shares all sub-chains across μ, which is where
the μ = (n)-is-the-degree cost concentration (§2.2) gets amortized away
instead of paid per shape.

**As built (2026-07-30), two refinements the prototype did not need.**

1. *Decompose the strip by runner.* "Choose m beta numbers and add k to
   each" reads as a `C(rows, m)` subset enumeration, and at the shapes
   §2.2 says are the whole cost that is fatal: λ = (kn) has
   `rows = kn`, so weight-m strips would be `C(kn, m)` candidates —
   `C(28,14) ≈ 4·10⁷` at the k = 2, n = 14 target. But `β ↦ β + k` never
   leaves its residue class, so the choice **factors over the k runners**,
   and within a runner a set of beads can move up one slot each and stay
   distinct exactly when it is a **prefix of a maximal block** of occupied
   slots (move a bead with an occupied slot above it and the two collide; a
   maximal block's top always has room). A strip is therefore one prefix
   length per block, and for λ = (kn) that is *one* candidate where the
   subset reading has 4·10⁷. Crossings come out as `k−1` popcounts —
   an unmoved bead is passed iff it sits at `B+i` for some `0 < i < k`.
   The naive subset form is kept as a test oracle
   (`strip_blocks_agree_with_naive_subsets`), since the two share nothing.
2. *Walk up, not down, and let the weight trie be the memo.* Building from
   ∅ by strips of weakly **decreasing** weight makes the recursion tree the
   partition trie, so every weight *prefix* is walked once and shared by all
   the partitions extending it — the sharing this section wanted, without a
   memo table. Two entry points fall out of the same walk: pruned to
   subshapes of one λ (exact — a chain to λ never leaves `⊆ λ`) it is
   `llt_gtilde`; unpruned it is `llt_gtilde_table`, **every** `G̃^(k)_λ` of a
   degree from one pass, which no incumbent exposes. `the_two_walks_agree`
   holds them together.

The abacus is a `u128` bitmask (occupancy, the `+k` move, and the crossing
count are then bit ops), which caps `ℓ(λ) + λ₁ + k < 128` — comfortably past
the n = 14 targets at k ≤ 4, and an asserted limit rather than a silent one.

### 3.4 R3 — Fock straightening: the Schur/KL engine

[KMS] (43)/(45) verbatim (§1.3(g)). States are strictly decreasing
integer wedges; `V_m` (the h_m-shaped boson) adds k to a multiset of m
positions; straightening normal-orders with the two rules; then

```text
  S_λ |ρ⟩ = Σ_ν κ_{λν} V_ν |ρ⟩   (κ: s_λ in the h-basis, i.e. one
                                   determinant row per λ)
  ⟨μ + ρ| S_λ |ρ⟩ = c^λ_μ(q) at q = −v     (verified deg ≤ 4, k = 2,
                                            multi-term entries included)
```

One run with fixed λ yields **one column of the Schur-expansion table
for every μ at once** — the transpose of what tableau enumeration
produces — and the entries are parabolic affine KL polynomials ([LT]
Thm 4.2). This is the route for "give me `G̃_λ` in the Schur basis
without enumerating fillings", and its wedge arithmetic (sorted integer
vectors, signed ℤ[v] coefficients) is `charge.rs`-adjacent code, not a
new algebra layer.

### 3.5 Vertical strips and graphs: the specials

Four presentations of the same objects, each earning its keep (the
graph dictionary subtleties are §1.3(e)/(f)):

- **dinv model** (the shuffle-world native): `dyck.rs` *already
  enumerates this* — `for_each_labelling` (dyck.rs:117, strictness at
  :130) with `d_row` (dyck.rs:148) is precisely `Σ q^{dinv} x^label`
  per area sequence. `G_D` as a first-class object is a regrouping of
  code that exists and is measured.
- **Tuple form**: runs-in-reverse-order (§1.3(e)) hands `G_D` to R1 for
  F-expansions and standardization.
- **Coloring form** (the [AP]/[DA] shape): vertices, oriented weak
  edges, strict edges; equals the dinv model on the dinv-faithful graph
  (§1.3(f), verified n ≤ 6). This is the form the **chromatic bridge**
  consumes (p-basis twist `p_r ↦ (q^r − 1)-scaled`, verified against
  Sage's chromatic QSF on all unit-interval graphs n ≤ 5) — the
  `chromatic-corpus` branch's client.
- **[AS] orientation formula** for the *e-expansion*: sum over
  orientations of free edges, `q^{asc}`, block partition by
  highest-reachable-vertex along strict+ascending edges; then
  `G(x; q+1) = Σ_θ q^{asc(θ)} e_{λ(θ)}` — 2^{#free edges} terms,
  certified e-positive output by [DA]'s theorem. Fine to n ≈ 20 for
  path-shaped graphs (sparse); the [DA] path-algebra recursion
  ((5.1)–(5.4), `d_− φ^m d_+` = multiplication by e_{m+1}) is the
  recorded alternative if orientation counting walls.

### 3.6 The ∇e_n by-path decomposition — the record's debt, paid

`docs/record/dyck-paths.md` (the Dyck-ladder postmortem) recorded two obstructions to
computing `Σ_labelings q^{dinv}` without enumeration: (i) dinv is not
invariant under either tie-breaking convention when standardizing
labelled paths, and (ii) `Val`'s tie clause is destroyed by
standardization. Resolution, verified:

- (i) dissolves **in the tuple model**: [HHL] (82)'s standardization is
  on tuple fillings with the content-reading-order tie-break, and it is
  an *identity* (§5.A), not a convention choice. The path's labellings
  map to tuple fillings by the reverse-run dictionary (§1.3(e)),
  pointwise in q (§5.I). So `G_D` in F- or m-basis comes from `#SYT`
  standard objects, not `#labelings` — the exponential-to-polynomial
  drop per path the postmortem wanted.
- (ii) is **out of scope correctly**: `Val` is not an LLT statistic;
  the valley side of the Delta conjecture stays with `dyck.rs` and its
  honest enumeration. Nothing here claims otherwise.
- The assembled statement `∇e_n = Σ_D t^{area} G_D` is verified n ≤ 5
  against Sage's nabla, and `deltaop.rs`'s own `nabla_e` (deltaop.rs:979
  ff.) gives the in-crate cross-check at every degree the engine
  reaches.

### 3.7 The [BHMPS] route — recorded, deferred

Their eq (4) writes `ω∇^m G_ν` as an explicit Catalanimal — ∇ of *any*
LLT by raising-operator combinatorics, no Macdonald basis pass. This is
the v2 route for pushing `∇(LLT)` past what `deltaop.rs` + basis
conversion can do, and the natural follow-on once the `st`-basis work
and this module coexist. Not specced further; the row exists so the
research-gaps table can point here.

### 3.8 Recorded dead ends

- **[HHL]'s `q^e G(1/q)` remark as the quotient dictionary** — the
  tilde collision (§1.3(b)). Cost one debugging round; the direct
  min-inv-floored equality is the law.
- **Recovering spin from tuples** (§1.3(c)) — impossible, by fixture.
- **The cell-below-path graph as the dinv dictionary** (§1.3(f)) — it
  is the *other* graph; a = (0,0) separates them. Both are kept, for
  different jobs.
- **[LLT] §7's printed straightening rule** — two misprints; [KMS] is
  normative (§1.3(g)).
- **Building the k-quotient with nonzero offsets** to chase [HHL]'s
  stated content conventions — offsets 0 with the min-inv floor is the
  verified normal form for empty-core shapes; general-core offsets are
  an open question (§7), not a v1 blocker.

### 3.9 Crate fit

Verified against the source today, with line numbers:

| need | existing | file |
|---|---|---|
| coefficient ring ℕ[q,t], sparse | `QtPoly<C>` | qt.rs:52 |
| vertical-strip enumeration + dinv | `for_each_area` / `for_each_labelling` / `d_row` | dyck.rs:95/117/148 |
| the ∇e_n consumers to refine | `ladder`, `side` | dyck.rs:215/299 |
| beta-set strip walking (R2 template) | `charge::build` / `strips` | charge.rs:154/177 |
| the k-large degeneration oracle | `hall_littlewood`, `_p`, tables | hl.rs:57/124/68 |
| ∇ cross-check for §3.6 | `nabla`, `nabla_e`, `delta` | deltaop.rs:979/1002 |
| q=1 gate (products of skew Schurs) | `expand_skew` | skew_lr.rs:56 |
| arm/leg/content helpers | `arm`, `leg` | macdonald.rs:214/219 |
| m-basis container, Schur container | `Monomial`, `Schur` (`basis!`) | sym.rs |
| H̃ cross-check for the [HHL] assembly | the `qtkostka` routes | qtkostka.rs |
| partitions, dominance, z_λ | `partition.rs` | partition.rs |
| per-degree memo tables | `table!` pattern | memo.rs:42 |
| Python boundary + guards | `guarded`/`escalate`, `#[pymodule]` list | guard.rs / python.rs |
| Sage conformance harness | `check_backend.py` pattern | scripts/check_backend.py |
| committed oracle fixtures | `gen_sage_oracle.sage` → `tests/fixtures/` | tests/sage_oracle.rs:18 |

Genuinely new, all small: **k-core/k-quotient on `Partition`**
(`partition.rs` has neither — grep-verified today), a **`SkewTuple`**
type (shapes + offsets + attacking-pair iteration), **wedge
straightening** (R3), and a **composition-indexed container for the
F-expansion** (the crate has no QSym anywhere — grep-verified; a
`Vec<(Vec<u32>, QtPoly)>` with documented reading-word convention is
enough, no Hopf structure needed). `lib.rs` slot: `pub mod llt;` between
`kostka` and `lr` (lib.rs:70/71); no name collisions (`llt` appears
nowhere in `src/`).

---

## 4. API

```rust
// llt.rs — coefficients in QtPoly<C>; q is the LLT variable, t reserved
// for area/maj gradings in the assembly functions.

/// A tuple of skew shapes with content offsets — the [HHL] object.
pub struct SkewTuple { /* (outer, inner, offset) per component */ }

impl SkewTuple {
    pub fn from_partitions(shapes: &[Partition], offsets: &[i32]) -> Self;
    /// The k-quotient of λ (empty k-core required), offsets 0 — the
    /// §1.3(a) normal form.
    pub fn quotient(lambda: &Partition, k: u32) -> Option<Self>;
    /// Vertical strips of a Dyck path, reverse-run order (§1.3(e)).
    pub fn from_area(area: &[u32]) -> Self;
    pub fn conjugate(&self) -> Self;           // for the ω-duality law
    pub fn attacking_pairs(&self) -> usize;
}

/// G_ν, raw inv grading (min-inv floor NOT divided out; callers get
/// `min_inv` to normalize — the honest exposure of §1.3(a)).
pub fn llt_g<C: Ring>(nu: &SkewTuple) -> Monomial<QtPoly<C>>;
pub fn llt_min_inv(nu: &SkewTuple) -> u32;
/// Cospin ribbon functions: G̃^(k)_λ, H̃^(k)_μ, spin H^(k)_μ (partition
/// side only — spin is shape data, §1.3(c)).
pub fn llt_gtilde<C: Ring>(lambda: &Partition, k: u32) -> Monomial<QtPoly<C>>;
pub fn llt_h_tilde<C: Ring>(mu: &Partition, k: u32) -> Monomial<QtPoly<C>>;
pub fn llt_h<C: Ring>(mu: &Partition, k: u32) -> Monomial<QtPoly<C>>;
/// Whole degree — the unit §2.2's walls are measured in.
pub fn llt_h_table<C: Ring>(n: u32, k: u32) -> Vec<(Partition, Monomial<QtPoly<C>>)>;
/// Schur expansion via R3 (one λ-column per straightening run); ribbon
/// side for partition shapes, R1+LR for general tuples.
pub fn llt_schur<C: Ring>(lambda: &Partition, k: u32) -> Schur<QtPoly<C>>;
/// Fundamental quasisymmetric expansion via R1's buckets.
pub fn llt_fundamental<C: Ring>(nu: &SkewTuple) -> Vec<(Vec<u32>, QtPoly<C>)>;
/// The by-path shuffle refinement: (area sequence, G_D) pairs with
/// Σ_D t^{area} G_D = ∇e_n — checked against deltaop::nabla_e in tests.
pub fn nabla_e_by_path<C: QAlgebra>(n: u32) -> Vec<(Vec<u32>, Monomial<QtPoly<C>>)>;
/// Unicellular / vertical-strip via decorated graphs (weak+strict,
/// oriented), the coloring model; consumed by the chromatic bridge.
pub struct DecoratedGraph { /* n, weak (oriented), strict */ }
pub fn llt_graph<C: Ring>(g: &DecoratedGraph) -> Monomial<QtPoly<C>>;
/// X_Γ(x;q) from G_Γ by the (q−1)-twist — the chromatic-corpus client.
pub fn chromatic_from_llt<C: QAlgebra>(g: &DecoratedGraph) -> Monomial<QtPoly<C>>;
/// [AS] e-expansion of G(x; q+1) — certified positive ([DA]).
pub fn llt_e_expansion<C: Ring>(g: &DecoratedGraph) -> Vec<(Partition, QtPoly<C>)>;
```

Python bindings mirror the Macdonald/Jack pattern (`python.rs`), through
`guarded`/`escalate`; boundary payloads are the existing
`(partition, [(qexp, texp, coeff)])` term shape. Not in v1: ∇ of general
LLTs ([BHMPS], §3.7), nonzero-core quotients (§7), compositional/Val
refinements (explicitly out of scope, §3.6), Schröder-word parsing
(graphs are the API; words are a constructor the tests use).

**As built, 2026-07-30.** Every function above exists with these signatures,
plus four the implementation wanted: `llt_g_lt` ([LT] (43)'s grading, the one
the Fock route is pinned against — the spec's §1.1 lists four normalizations
and the API exposed three), `llt_max_inv` (beside `llt_min_inv`, same
argument), `llt_gtilde_table` (§3.3's unpruned walk), and
`SkewTuple::from_cells` / `from_skews` (the [HHL] Macdonald components are
arbitrary cell sets that walk left out of the first column, not skew shapes,
so the general constructor is load-bearing rather than a convenience).
`llt_kl_column` is the name `llt_schur`'s doc-comment promised for the column
half. `htilde_by_llt` is the [HHL] assembly of §5.8, which the API section
forgot to list.

**The Python bindings, built 2026-07-30.** Sixteen functions in `python.rs`,
following the (q,t)-Kostka family's `i128` boundary rather than the
`guarded`/`escalate` ladder — every coefficient here counts tableaux, so it is a
non-negative integer bounded by `n!` (8.7e10 at n = 14 against `i128`'s 1.7e38),
and `bench_llt` runs the whole ladder at both `i64` and `i128` asserting term-for-
term agreement, so the narrower width is checked rather than assumed.

Payloads reuse the existing `(partition, [(qexp, texp, coeff)])` shape, under a
distinct `QtMon` alias: structurally identical to `QtSchur` and deliberately
named apart, because these partitions index **weights** and feeding one to an
operator expecting Schur shapes would type-check in Python and be wrong. The `t`
slot is zero for the LLT families proper and carries the area/maj grading only in
`nabla_e_by_path` and `htilde_by_llt` — which the boundary tests assert rather
than assume, since a stray `t` exponent would otherwise be silently dropped.

Three things get *validated* at the boundary instead of trusted, because each
would otherwise return a plausible wrong answer: strict edges must run `u < v`,
weak and strict edge sets must be disjoint as unordered pairs (or the ascent
statistic gains a `q` per strict edge), and edges must land inside `0..n`. All
three raise `ValueError`, and `check_bindings.py` asserts they do.

§5.13's boundary half is now covered: 29 checks appended to
`scripts/check_bindings.py` — hspin/hcospin/cospin against Sage in **both** m and
s, the floored tuple dictionary with `llt_g` and `llt_min_inv` crossing
separately (and an assertion that the sweep contains a nonzero floor, or a broken
pair would look fine), table index sets and per-shape agreement, the k-core and
k-quotient **component order** against Sage's own, `nabla_e_by_path` summing to
Sage's `nabla` with `C_n` pieces, path graphs against Sage's chromatic QSF, the
isolated-vertex fixture, and the three rejection paths. **0 failures.**

## 5. Correctness requirements

Layered as the crate does it. Items 1–13 **already pass in the Python
prototype** (`spec_llt_verify.py`, 26 PASS lines); the requirement is
that the Rust port reproduces them and extends the ranges.

> **Done 2026-07-30.** All thirteen layers are covered by the 28 tests in
> `src/llt.rs` (plus 4 in `partition.rs` for the abacus primitives), and item 4
> — the Sage oracle dictionaries — became `examples/llt_dump.rs` →
> `scripts/check_llt.py` rather than committed fixtures, following the
> `check_hl.py`/`check_jack.py` pattern that already exists: dump from a release
> binary, compare in Sage, never read Sage's source. That run is **1517
> comparisons, 0 failures**, and it checks two things nothing in-crate can — the
> quotient's *component order* against Sage's own `Partition.quotient`, and the
> floored tuple dictionary with the floor crossing the boundary separately, so a
> tuple whose floor happens to be zero cannot make the check pass by accident.
> Item 13's Python-boundary half is **not** built (§4).
>
> Two layers came back with corrections rather than confirmations: item 2's
> floor witness (§1.3(a) ⚠️) and item 12's range (§1.3(g) ⚠️). Both were spec
> errors, not port errors.

1. **The convention gate, first.** `G_{((1),(1))} = m₂ + (1+q)m₁₁` — one
   attacking pair, and the line that dies if reading order, attack rule,
   or inv orientation drifts.
2. **Model equivalences** (mine vs mine, shared-nothing): R1 ≡ direct
   SSYT (8 tuples); ribbon G̃ ≡ floored tuple G on k-quotients (51
   empty-core shapes, k = 2,3 — the §1.3(a) law); dinv model ≡
   reverse-run tuple ≡ coloring-on-dinv-graph, pointwise per path,
   n ≤ 6; offsets are load-bearing (a shifted-offset tuple must differ).
3. **Published fixtures**: [LLT] Ex 6.8(i) in m and s; Ex 6.8(ii); [LT]
   Ex 4.1 (the (43) grading); [AS] Ex 6.1; [DA] Ex 5.6. Rust must carry
   these as unit tests with the values inlined.
4. **Sage oracle dictionaries** (§1.3(d)): hspin/hcospin sweeps k = 2, 3
   for |μ| ≤ 4; `cospin(λ)` on the k = 3 shape list; `cospin(tuple)` =
   `q^{−min inv} G_ν` on the 9-tuple list — as committed fixtures per
   the `gen_sage_oracle.sage` pattern, regenerated only with the wart
   list (§2.1) in view.
5. **Internal gradings**: H = q^{s*} G̃(1/q) on even-spin shapes; cospin
   integrality everywhere; `s* ` consistency between `ribbon_g`,
   `ribbon_g_spin`, `ribbon_g_cospin` ports.
6. **Degenerations**: q = 1 products of skew Schurs (vs `skew_lr`);
   [LLT] Thm 6.6 vs `hl.rs`'s own Q′ (in-crate, not Sage) on a shared
   range; k = 1 ≡ Schur.
7. **Symmetries**: ω-duality with the `q^{A(ν)}` twist; the [BHMPS]
   coproduct on split alphabets.
8. **The [HHL] assembly**: Σ_D q^{−a} t^{maj} G_{ν(μ,D)} = H̃_μ for all
   μ ⊢ n ≤ 5 against `qtkostka.rs` (in-crate) — this is the fourth
   H̃-route and must agree with the existing three.
9. **The shuffle refinement**: Σ_D t^{area} G_D = `deltaop::nabla_e(n)`
   for n ≤ 8 in Rust (n ≤ 5 in the prototype); per-path Schur
   positivity **recorded, never repaired** — a negative coefficient is
   a finding.
10. **The chromatic bridge**: p-twist against a committed Sage
    chromatic-QSF fixture set (all unit-interval graphs n ≤ 5), with
    isolated vertices in the graph fixture (the empty-edge-list trap
    cost this session a debugging round: a two-vertex edgeless graph
    must not collapse to the empty graph).
11. **The [AS] e-expansion**: equals `llt_graph` after q ↦ q+1 twist on
    every graph n ≤ 7; [DA] positivity of the q+1 coefficients asserted
    (theorem, so a violation is our bug — the *reverse* of item 9's
    posture, deliberately).
12. **R3 Fock**: KMS rules only; c^λ_μ at q = −v vs the ribbon side,
    k = 2 deg ≤ 4 in the prototype, k = 2, 3 deg ≤ 5 in Rust; the
    multi-term KL entries must be present in the range or the range is
    too small to pin the dictionary (§1.3(g)).
13. **Fixed-width honesty + bindings**: i64/i128/bignum ladders agree
    (`macop.rs` two-width pattern); `check_bindings.py`-style slot tests
    at the Python boundary.

## 6. Performance requirements — and the measurement

Stated before the code existed, so the measurement could embarrass them
(house precedent: `st` guessed 50×, got 3400×; Δ-operators guessed 100×,
got ~21×; Jack guessed 200×, got 8340×). **Measured 2026-07-30 on mains**
(`cargo run --release --example bench_llt -- 14`, single-threaded, `i64`
column), and **both sides of every ratio are mains** — §2.3's Sage re-run, not
§2.2's battery rows. The standing mains-to-mains rule is satisfied here rather
than caveated; §2.2 is kept for provenance, not for arithmetic.

| target | guessed | ours | Sage (mains) | verdict |
|---|---|---|---|---|
| `H^(3)` table, n = 9 | ≥ 100× (≤ 0.9 s) | **0.008 s** | 47.34 s | **5900×** |
| `H^(2)` table, n = 10 | — | **0.014 s** | 24.16 s | **1700×** |
| `H^(4)` table, n = 8 | — | **0.0053 s** | 101.45 s | **19 100×** |
| k = 2, 3 tables through n = 14 | "minutes each" | **0.25 / 1.08 s** | dies at n = 11 / 10 | met, by ~250× |
| tuple `((3,2),(2,2),(2,1))`, n = 12 | < 1 s | **0.006 s** | 66.99 s | **13 700×** |
| `nabla_e_by_path(10)`, 16 796 pieces | "minutes" | **19.3 s + 3.7 s** | *no entry point* | met — see §2.3 |
| KL columns, all λ ⊢ 6, k = 2 | seconds *per column* | **0.041 s for all 11** | *no entry point* | met |

The two "no entry point" rows are the honest half: they are not speedups, they
are objects Sage cannot produce at any speed. And §2.3 records the one place
Sage is competitive — `e[n].nabla()`, where `deltaop` rather than this module is
the right answer, and where the by-path route costs *more* than Sage's total
because it is not computing the total.

Whole-degree tables, mains, seconds:

```text
  k \ n      9        10       11       12       13       14    growth/degree
   2      0.008    0.014    0.022    0.045    0.107    0.249     1.6–2.5×
   3      0.008    0.022    0.058    0.157    0.409    1.085     2.6–2.9×
   4      0.016    0.049    0.142    0.412    1.140    3.111     2.7–3.1×
```

### 6.1 The profiling pass

The table above is *after* a sampling pass (`examples/profile_llt.rs`, one
attributable workload per route; `sample` on macOS). The first working version
measured 0.019 / 0.053 / 0.115 s at n = 10 and **80.4 s** for
`nabla_e_by_path(10)`. What the profile said, per route, and what came of it:

| route | first profile | fixed | gain |
|---|---|---|---|
| **R1** SYT walk | 85% in the walk, 6% hashing | bit-mask poset (`pred_mask`, `attack_mask`), direct-indexed descent table | **4.1×** |
| **R2** abacus strips | ~35% `malloc`/`free`; `strip_rec` 70% of the rest | scratch buffers, all-weights-in-one-walk, O(changed bits) containment, runner-mask block scan | **2.2–2.5×** |
| **R3** Fock straightening | **63% `malloc`/`free`**, `straighten` itself 6% | in-place wedge mutation, generated offset ladder, one-pass coefficient helpers, probe-before-insert | **3.0×** |

Three things worth keeping from that:

- **Every route's first profile pointed somewhere different**, and in two of
  three cases not at the mathematics. R3 was spending two thirds of its time in
  the allocator because the straightening recursion cloned its wedge per branch;
  the fix is that it mutates two positions and restores them. This is the
  [`crate::jack`]/[`crate::macop`] finding again — the standing note to
  re-sample the Δ-operator paths for allocator churn was right about the shape
  of the problem.
- **The naive form of each hot loop was the readable one**, so the fast form
  carries the reason in a comment and the slow form survives as a test oracle
  where one existed (`strip_blocks_agree_with_naive_subsets`,
  `abacus_containment_is_partition_containment`,
  `the_all_weights_walk_agrees_with_the_fixed_weight_one`). The last of those
  exists because the optimization moved production off the path the first test
  covered — worth checking for after any such change.
- **R1 and R2 are now compute-bound; R3 is still ~38% allocator**, from the
  `QtPoly` temporary each straightening branch creates. Left there: the next step
  would be threading a coefficient stack through the recursion, and the route is
  not on the critical path of any sweep (§7.4).

Where the remaining R1 time goes, and the algorithm *not* taken: the walk is
`#SYT(ν)` per path and that is now essentially all of it. A subset DP over
(assigned set, last cell) would be `2^n · n²` per content instead, i.e. better
once `#SYT` passes `p(n) · 2^n · n²` — around n = 12–13 for these tuples, so
**worse** at the n ≤ 10 the shuffle refinement wants. Recorded rather than
built.

Notes the numbers force:

- **The growth factor moved, not just the constant.** Sage's measured
  4.3–6.9× per degree is 1.6–3.1× here. That is the weight-trie sharing
  (§3.3): the marginal degree adds chains, not whole re-walks. "The first
  degree over an hour" is not reached in range — extrapolating k = 4 at
  2.7× puts it around n = 21.
- **`μ = (n)` stopped being the degree.** In Sage the one-row shape is
  94–100% of the whole-degree cost (§2.2); here it is `0.0002 s` at
  n = 10, k = 2 — under 1% — because the pruned walk restricts to
  subshapes of a single row, of which there are `n`. The design finding
  that shaped §3.3 dissolved the cost concentration it was aimed at.
- **The i64/i128 columns are within noise of each other** (0.626 vs 0.635
  at k = 2, n = 14), and are asserted equal term-for-term in the bench
  itself, so the fixed-width column is honest rather than merely fast.
- **`nabla_e_by_path` is enumeration-bound, as predicted.** 19.4 s to
  build all 16 796 pieces, 3.7 s to convert them all to the Schur basis and
  sum: the `#SYT` walk is 84% of it, and the m → s conversion the postmortem
  worried about is 16%. Verified exactly against `deltaop::nabla_e(10)`,
  and every one of the 16 796 pieces is Schur-positive.
- **Peak RSS is 220 MB for the whole bench**, dominated by
  `nabla_e_by_path(10)` holding all 16 796 pieces at once rather than by any
  ribbon walk. The per-phase split (enumeration vs bucketing vs conversion) is
  instrumented for the by-path table only; the rest is recorded as owed.

## 7. Open questions

1. **Nonzero-core quotients.** The §1.3(a) dictionary is verified for
   empty-core λ with offsets 0. What offset vector makes the tuple model
   match `G̃` for general k-core — and is the min-inv floor still the
   only normalization needed? (The literature states conventions that
   did not survive contact with the fixtures; measure, don't trust.)
2. **[LLT] Conj 6.4** — is `H^(k+1)_μ − H^(k)_μ` Schur-positive?
   Open per today's reading. A k-sweep at n ≤ 12 is a first-week
   deliverable of the R2 engine and would be the first serious
   computational evidence either way that we know of.
   **Swept 2026-07-30** (`bench_llt`, the Conj 6.4 section): `k = 1 … 4`,
   every `μ ⊢ n` for `n ≤ 14` — 56 (k, n) rows, **120 943 nonzero Schur
   coefficients of the difference, none negative**. The whole sweep costs 34.9 s
   on mains, which is
   the point — this was out of reach and is now a rounding error, so the
   honest next question is what range would be *informative* rather than what
   range is affordable. Still open, and still a conjecture: a sweep is
   evidence, and the bench prints `*** COUNTEREXAMPLE -- REPORT ***` rather
   than failing an assertion, so a violation stays a result.
3. **[AP] Conj 25 / [AS] Problem 6.20 / Shareshian–Wachs** — the
   unicellular e-positivity constellation. The `chromatic-corpus`
   branch is the consumer; the bridge (§3.5) is verified; the sweep
   infrastructure is that branch's spec, not this one's.
4. **Does R3 beat R2 for Schur output at scale?** R3 emits columns
   (fixed λ, all μ), R2+conversion emits rows. The crossover matters
   for the Conj 6.4 sweep, which wants rows.
   **Partly answered 2026-07-30, and the answer is no in the measured
   range.** All λ ⊢ 6 at k = 2 costs R3 0.041 s of straightening; the
   whole-degree R2 table it would have to beat (n = 6, k = 2) costs
   0.0007 s, and R2's own reach is n = 14 at 0.25 s. R3 also grows faster
   (~12× per degree at k = 3: 0.004 → 0.051 → 0.65 s for λ ⊢ 4, 5, 6),
   because the wedge straightening branches where the trie shares. So R2
   is the engine for the sweeps and R3 earns its place for what it *is*
   rather than what it costs: the columns are Kazhdan–Lusztig polynomials,
   an output R2 cannot produce at all. Still open: whether R3 wins at a
   large `k` with a small λ, where R2's abacus has many runners and few
   beads per runner.
5. **The [BHMPS] Catalanimal route** (§3.7) — v2; whether `∇G_ν` for
   vertical strips beats `deltaop`-side computation is the question the
   research-gaps row exists to answer.

7. **Where else in the crate is this engine an engine?** Measured
   2026-07-30 (`examples/probe_llt_ladder.rs`), two candidates, one win:

   - **`dyck::ladder`'s rise side factors through per-path LLTs — verified,
     and it is the big one.** [HRW]'s rise statistic selects
     `e_{n−1−k}` of the weights `t^{−a_i}` over `i ∈ Rise(D)`, and both the
     rise set and those weights are functions of the **area sequence
     alone** (`dyck.rs:255`) — no labels. So the whole `z`-extraction
     factors out of the labelling sum:

     ```text
       Rise_{n,k} = Σ_D [ Σ_{S ⊆ Rise(D), |S| = n−1−k} t^{area(D) − Σ_{i∈S} a_i} ] · G_D(x; q)
     ```

     which replaces one labelled-path walk per content with `C_n` LLT
     evaluations plus a knapsack.

     **Implemented 2026-07-30** in `dyck::ladder`, which now dispatches
     `Side::Rise` to `rise_ladder_via_llt`. The labelled walk stays as
     `ladder_at_content` — it is the oracle
     (`the_rise_ladder_via_llt_agrees_with_the_labelled_walk`, every k and
     every content, n ≤ 6) and it is genuinely cheaper for one coarse
     content, `μ = (n)` being a single labelling. Measured
     (`examples/probe_llt_ladder.rs`), **mains**:

     ```text
       n           5     6     7      8      9
       llt (s)  .0006 .0037 .0174 .0987  1.310
       walk (s) .0013 .0171 .1350 2.873 73.699
       speedup   2.1x  4.7x  7.7x  29.1x  56.3x
     ```

     ~2× per degree, so n = 10 extrapolates to ~110× — the walk goes from
     ~25 min to ~13 s. The same table on battery gave 26.1× / 56.6× at
     n = 8, 9 with the absolutes roughly doubled, which independently
     confirms §2's ~1.8× drift figure: the A/B *ratios* are
     power-independent, only the seconds are not.

     It also yields something the old ladder could not: a
     **per-path Schur-positive refinement of `Δ'_{e_k} e_n` for every k**,
     where §3.6 only claimed the `k = n−1` top.

     ⚠️ **Correction to a number first reported in this session.** The first
     probe measured 140× at n = 8 by looping `dyck::side(n, k, Rise)` over k —
     and `side` computes the *whole* ladder and discards all but one slot, so
     that charged the labelled walk n times over. The honest figure is 29×.
     Worth recording as the failure mode: a per-`k` API whose docs already say
     "asking for a single k costs the same as asking for all of them" will
     silently inflate any A/B built out of it.

     ⚠️ The **valley** side does *not* factor — `Val` reads the labels, and its
     weights are `q^{d_i+1}` per labelling — so it keeps the enumeration, and it
     is the open one. This makes the rise half of the comparison nearly free
     rather than moving the conjecture.
   - **The [HHL] route to `H̃_μ` is not a competitive fourth engine.**
     `htilde_by_llt` beats `bh`'s route by ~1.3× through n = 7, ties at
     n = 8 and loses **3×** at n = 9, growing — the `2^{|μ|−μ₁}` descent
     subsets win eventually and there is no pruning to add. It keeps its
     value as §5.8's independent, positively-graded cross-check, and that is
     all it should be asked to be.
6. **Positivity bookkeeping for [GH]-only cases.** For skew (not
   straight) tuples, Schur-positivity rests on an unpublished preprint.
   The engine should *flag* skew-tuple Schur expansions with a negative
   coefficient as findings of the first order (they would be either our
   bug or mathematics).
