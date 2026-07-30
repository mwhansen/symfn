# Specification: Schubert polynomials (`schubert`)

*Written 2026-07-29. Implements the first item scheduled from docs/record/README.md's
"Beyond the core (deferred, but intended)" list. This is the one Symmetrica
subsystem Sage still actively routes through — retiring it completes the
displacement the README opens with.*

Sources, pinned:

- **Knutson, *Schubert polynomials, pipe dreams, equivariant classes, and a
  co-transition formula*, [arXiv:1909.13777v2](https://arxiv.org/abs/1909.13777)**
  — the v2 PDF read today, 2026-07-29, not recalled. Proposition and section
  numbers below are that paper's.
- **Symmetrica's C** (`Symmetrica_2.0/sb.c`, `mss.c`, `perm.c`), read today at
  function level. Symmetrica is public domain, so it is a legitimate *source of
  algorithms*, not merely an oracle — the same standing `NOTICE.md` records for
  the Morris recursion taken from `sr.c`. Function:line references below.
- Primary literature for provenance (not fetched; the formulas were instead
  **verified numerically today**, see below): Lascoux–Schützenberger 1982
  (definition), Monk 1959, Lascoux–Schützenberger 1985 (transition), Macdonald,
  *Notes on Schubert polynomials* (1991), Billey–Jockusch–Stanley 1993.

**Every formula in §3 was checked against Sage before it entered this
document.** SageMath 10.9 (2026-05-04), this machine, today: the transition
formula (exhaustive on S₅ — 119/119 — plus 20 random elements of S₈), the
signed Monk rule (exhaustive on S₄ × i ∈ 1..4), the pipe-dream example C₁₄₂₃
(Knutson p.6), the dominant base case, stability under S_n ↪ S_{n+1}, the
lex-triangularity `x^{code(w)} = lex-min monomial of S_w` (40 random S₇), and
Macdonald's reduced-word evaluation `S_w(1,…,1) = (1/ℓ!) Σ_{red words} Π aᵢ`
(4 cases). The harnesses are committed as `scripts/spec_schubert_verify.py`
(formula checks), `scripts/spec_schubert_walls1.py` / `_walls2.py` (Sage
walls), and `scripts/spec_schubert_schubmult.py` (the C incumbent); the
formula checks become committed tests in §5.

**Revision 2026-07-30.** The two items this document left open as
prerequisites are now closed, both before any Rust exists:

- **The Grassmannian flag-truncation identity holds** (was the one ⚠️
  unverified formula; §7 Q2). Exhaustive over same-descent Grassmannian pairs
  for k ≤ 4 with |λ|,|μ| ≤ 5 — 760 pairs, 0 failures — plus nine larger cases
  through k = 5, with a different-descent negative control that fails as it
  should. `scripts/spec_schubert_grass.py`. §3.5 is rewritten accordingly, and
  the measurement turned up a *cost* consequence the identity's truth hides
  (the flag discards 98.7% of the LR expansion at k = 7), which changes what
  leaf dispatch needs from `AutoLr`.
- **The E3 compression ratio is measured** — §3.5 called this "the first
  implementation step, because it predicts whether E3 beats E2 before either
  is tuned." It is 687× on `stair7`, 26 338× on `stair8`, and 125 599× on the
  S₁₃ element where the C incumbent times out.
  `scripts/spec_schubert_peel.py`. §3.5, §6 and §7 Q1 are rewritten
  accordingly. These are **exact combinatorial counts, not timings**, so
  unlike §2 they carry no power-state caveat.

**Clean-room boundaries.** Three Schubert implementations exist and their
licenses differ, so the posture differs per source and must be kept straight:

| code | license | may we read it? | role |
|---|---|---|---|
| Symmetrica `sb.c` | public domain | **yes** (read today) | algorithm source + oracle |
| lrcalc's `schubmult` (Buch, C; ships inside Sage's tree) | GPL | **no** | out-of-process oracle, measured |
| `schubmult` 4.1.0 (Samuel, Python/SymPy, PyPI 2026-07-14) | GPL-3 | **no** | prior art; oracle if ever needed |
| Sage's `schubert_polynomial.py` wrapper | GPL | **no** (runtime introspection only) | oracle through the public API |

---

## 1. What the thing is

Schubert polynomials `S_w`, indexed by permutations `w ∈ S_∞` (eventually-
identity permutations of ℕ₊), form a ℤ-basis of the *whole* polynomial ring
ℤ[x₁, x₂, …] — not of Sym. Two equivalent definitions:

**Divided differences** (Lascoux–Schützenberger; Knutson §2). With
`∂ᵢ f := (f − sᵢf)/(xᵢ − xᵢ₊₁)`,

```text
    S_{w₀⁽ⁿ⁾} = x₁ⁿ⁻¹ x₂ⁿ⁻² ⋯ x_{n−1}          (w₀⁽ⁿ⁾ = n, n−1, …, 1)
    ∂ᵢ S_w = S_{w sᵢ}   if w(i) > w(i+1),   else 0
```

**Dual basis** (Knutson Prop 2.1): the `{S_w}` are the ℤ-basis of ℤ[x̲] dual
to the nil Hecke algebra basis `{∂_w}` under a perfect pairing. `deg S_w =
ℓ(w)`, and `S_w` is stable under `S_n ↪ S_{n+1}` (verified), which is what
lets one index by `S_∞` and is the invariant our `Perm` type will carry by
construction.

**Why anyone wants this.** The structure constants

```text
    S_u · S_v = Σ_w  c^w_{uv} S_w
```

are the intersection numbers of Schubert varieties in full flag manifolds:
non-negative for geometric reasons, with **no combinatorial rule known in
general** — one of the standing open problems of algebraic combinatorics. When
u, v are Grassmannian (single descent), `c^w_{uv}` specializes to
Littlewood–Richardson. Computing these numbers *is* the use case; positivity
searches and rule-hunting are computationally starved in exactly the way
`docs/research-gaps.md` documents for Kronecker coefficients.

Facts verified today and kept as tests (§5):

- **Dominant** w (weakly decreasing Lehmer code): `S_w = x^{code(w)}`, one
  monomial.
- **Grassmannian** w (single descent at k): `S_w = s_λ(x₁,…,x_k)`, the flagged
  Schur polynomial, λ = code sorted increasing → partition.
- `S_{1423} = x₁² + x₁x₂ + x₂²` (Knutson p.6, pipe-dream example).
- **Monk, signed form**: `xᵢ·S_w = Σ_{j>i} S_{w t_{ij}} − Σ_{j<i} S_{w t_{ji}}`,
  both sums over Bruhat covers only. (Equivariant version: Knutson Prop 2.3.)
- **Transition** (Lascoux–Schützenberger 1985; the choice below matches
  Symmetrica `mss.c`): let r = w's **last** descent, s = the **largest** s > r
  with w(s) < w(r), v = w·t_{rs} (so ℓ(v) = ℓ(w)−1). Then

  ```text
      S_w = x_r · S_v + Σ_{q < r, ℓ(v t_{qr}) = ℓ(v)+1} S_{v t_{qr}}
  ```

- **Lex-triangularity**: `x^{code(w)}` is the lex-*minimal* monomial of `S_w`.
  This is what makes the greedy from-polynomial peel (§3.4) correct, and it is
  the invariant Symmetrica's `t_POLYNOM_SCHUBERT` silently relies on.
- **Macdonald's reduced-word identity**: `S_w(1,…,1) = (1/ℓ(w)!) Σ Π aᵢ` over
  reduced words (a₁,…,a_ℓ) of w. Shares no machinery with any product engine —
  this is the independent checksum, playing the role hook-content plays in
  `examples/verify_specialization.rs`.

### 1.1 Prior art, and what is actually ours

| who | what they have | where it stops |
|---|---|---|
| **Symmetrica** (= what Sage ships) | product, expand, from-polynomial, ∂, ∂_{w₀} pairing, Stanley/Schur via transition | Product walls at S₁₀ (§2): 6.9s on `stair5²`, timeout on `stair6²` and every random S₁₁ pair tried. No memoization; cost = #pipe dreams. Aborts at interpreter exit with "permutation memory not freed?" after heavy use (observed today, `mem_counter_perm = 215164`). |
| **lrcalc's `schubmult`** (Buch, C — installed inside Sage's own tree at `local/bin/schubmult`) | single-Schubert products, CLI | The actual frontier. `stair6²` 0.35s, random S₁₂ ~2–21s, S₁₃ up to 106s / 3.2M terms, and one S₁₃ case (ℓ(u)=25, ℓ(v)=36) **>120s**. |
| **`schubmult` 4.1.0** (Samuel, PyPI, GPL-3, released 2026-07-14 — actively maintained) | single, double, quantum, quantum-double, parabolic; per its site, "transitioning to the elementary symmetric monomial basis (with a twist) and then using the appropriate Pieri formula" (Sottile) | Python/SymPy; breadth is its point, not headroom on single products. Not measured here (not installed). |
| **Sage native** | nothing — `SchubertPolynomialRing` is a pure Symmetrica wrapper (`combinat/schubert_polynomial.py` is the only Sage library module using the Schubert bindings) | its engine is the Symmetrica column |

Honest accounting: **the mathematics is entirely prior art** — transition,
Monk, the peel, the pipe-dream expansion are all in the literature and all
already implemented in the three columns above. What is ours is narrower:

1. **Retiring Symmetrica from Sage is worth 400–1000× to Sage users by
   itself.** Sage's Schubert ring dies on `stair6²` while the C sitting in
   Sage's *own tree* does it in 0.35s. Nothing about that gap is deep; it is
   the same story as the classical-basis conversions, and the same fix.
2. **A frontier formulation of the Monk chain** (§3.5, E3). Symmetrica's
   engine touches one pipe dream at a time, twice over (expand, then re-walk
   per monomial). The crate's repeated experience — `SkewLr`, `kostka`,
   `coproduct` — is that the win is *merging partial states that agree on the
   undecided part*, not enumerating. **It transfers** (measured 2026-07-30,
   §3.5): the peel DAG has 2 740 states where the tree has 31 877 880 leaves
   on the very pair the C incumbent cannot finish. Still a node count and not
   a runtime — but it is no longer a hypothesis about whether the structure
   is there.
3. **Leaf dispatch into the LR engine.** Transition trees bottom out in
   Grassmannian permutations; products of those are flag-truncated LR
   expansions (identity **verified 2026-07-30**, §3.5), and the fastest LR
   machinery anywhere is already in this crate (`AutoLr`, `strip_lr.rs:264`).
   None of the incumbents can have this, since none of them has the LR
   frontier. Contingent on the row-bounded product of §7 Q6: the flag discards
   98.7% of the expansion at k = 7, so dispatching to an unbounded LR product
   hands back most of the advantage.
4. **A single-coefficient query** `c^w_{uv}` (§3.7). No package has one — the
   same "whole product to read one number" defect `docs/research-gaps.md` §2.1
   records for Kronecker. **Built and measured 2026-07-30**, and it is the
   item on this list that turned out to matter most. E2 with Bruhat pruning:
   **578×** faster than the whole product on `S_13.2` (0.0091s against
   5.10s), and — the part no amount of engine tuning could have delivered —
   it answers `S_13.0`, whose full expansion has a monomial mass of 4.3×10¹⁶
   and **cannot be materialised on any machine**, in 0.02–0.16s.

   That reframes the whole §6 exercise. Beating the C `schubmult` by 1.9–30×
   on whole products is a real result but an incremental one; being able to
   ask a question that no existing package can express, on inputs where the
   whole product provably does not fit in memory, is the difference in kind.
   Per `lr_coeff`'s lesson the query recurses on whichever factor has the
   smaller transition tree (`transition_tree` builds no elements, so asking
   costs microseconds). ⚠️ That side choice is *directionally* right — it
   picks `S_13.0`'s 1 477-node factor over its 4 909-node one — but measured
   only as 0.414s → 0.277s across five targets with one row getting worse;
   the absolute times are small and noisy, so it is not the clean 21× the LR
   case was.
5. **The only maintained permissive implementation.** Symmetrica is public
   domain but dead (and leaks); both maintained engines are GPL. A
   MIT/Apache-2.0 wheel matters to anyone embedding this outside Sage.

## 2. Measured walls

SageMath 10.9, this machine, 2026-07-29. Single runs, per-case alarm of
120s. ⚠️ **On battery power (68%, discharging)** — per docs/record/README.md's standing
warning these are order-of-magnitude walls, not benchmarks; ratios between the
two tables are safe (adjacent runs, same conditions), absolute times are not.

**Products through Sage** (`SchubertPolynomialRing(ZZ)`, i.e. Symmetrica):

```text
  case                                   sec      terms   max|coeff|
  stair4² (S_8)                        0.078         16       4
  stair5² (S_10)                       6.85          59      16
  stair6² (S_12)                       >120           —       —
  rand S_10  ℓ=21,22                   47.2        1070       5
  rand S_10  ℓ=19,18                   0.0004        12       1
  rand S_10  ℓ=28,20                   0.0067        25       1
  rand S_11  ℓ=32,27 / 35,26 / 31,21   >120 (all three)
  rand S_12  (three pairs)             >120 (all three)
```

`stairk` is the Grassmannian permutation [2,4,…,2k,1,3,…,2k−1] whose Schubert
polynomial is s_{(k,k−1,…,1)}(x₁..x_k) — chosen because its stable shadow is
the staircase LR family this crate measures everything on. For calibration:
the *unflagged* s₍₅,₄,₃,₂,₁₎² costs `AutoLr` **1.4ms** (docs/record/README.md) against
Symmetrica's 6.85s for the flagged version, which has *fewer* terms.

**Same cases, lrcalc's C `schubmult`** (out of process, min-of-3 below 0.5s).
⚠️ Every row carries ~3 ms of process startup and rows ≤ ~0.03s are mostly
that; see §6. **Re-measured 2026-07-30** and reproducible: 10 of 11 rows
within 5% of the values below (the outlier, `S_12 ℓ=27,21`, went 1.61 → 2.25s,
inside the record's stated ±30% noise (`docs/record/littlewood-richardson.md`)). That reproducibility is what licenses
using this table as a baseline at all — and it independently confirms the
machine was in the same state for both sessions, so the §2 Sage rows stand
too:

```text
  case                                   sec      terms
  stair4² (S_8)                        0.005         16
  stair5² (S_10)                       0.016         59
  stair6² (S_12)                       0.355        247
  stair7² (S_14)                      20.6         1111
  rand S_10  ℓ=21,22                   0.066       1070
  rand S_11  ℓ=32,27                   0.222      30143
  rand S_11  ℓ=35,26                   0.748     118822
  rand S_11  ℓ=31,21                   3.86      185284
  rand S_12  ℓ=37,24                  21.0       230901
  rand S_12  ℓ=27,21                   1.61        8280
  rand S_12  ℓ=43,32                  19.0       242507
  rand S_13  ℓ=25,36                  >120            —
  rand S_13  ℓ=36,41                  20.9       947254
  rand S_13  ℓ=40,41                 105.8      3241903
```

Things to notice:

- **Sage ships the slow one.** 6.85s vs 0.016s on `stair5²` (430×), timeout vs
  0.355s on `stair6²` — and the fast binary is *in Sage's own install*. The
  displacement win requires no algorithmic novelty at all.
- **The C frontier is S₁₃–S₁₄, and its cost does not track output size.**
  `S_13 ℓ=25,36` times out while `S_13 ℓ=40,41` finishes with 3.2M terms.
  Whatever it recurses on, some shapes are pathological for it — the same
  "cost tracks something other than the answer" signature the `st` spec found
  in Sage, and a hint that a different recursion order has room.
- **Coefficients are small.** max|c^w_{uv}| ≤ 16 across everything measured,
  ≤ 9 on the random S₁₀–S₁₂ cases. Structure constants at reachable sizes are
  nowhere near `i64`; the `guard.rs` escalation still wraps everything on
  principle (the record (`docs/record/python-and-sage-interop.md`) records why fixed-width without a guard is not
  shippable).
- **Expansion has its own wall** (Symmetrica `expand`): 11.9s / 84 084
  monomials for one random S₁₂ element of ℓ=33. #monomials = S_w(1,…,1) grows
  super-exponentially; any engine that expands before multiplying inherits it.
- **The exit-time abort.** After the product sweep, Sage's process printed
  Symmetrica's `ERROR: permutation memory not freed?` banner with a quarter
  of a million live objects — in a library embedded in Sage for twenty years.
  Manual-memory C with mutable globals is the failure mode; it is the same
  argument `lib.rs` already makes against `OP`.

Two methodology traps hit today, recorded so they are not hit again:

- **Rectangle Grassmannians are dominant.** The first ladder used
  w = [k+1,…,2k,1,…,k] (S_w = s_{(kᵏ)}(x₁..x_k)) as a "hard LR case" — but
  s_{(kᵏ)} in exactly k variables is the single monomial (x₁⋯x_k)ᵏ, so every
  product collapsed to one term and the rows measured nothing. A Schubert
  benchmark must check the *flagged* object is non-trivial, not import
  intuition from the stable limit.
- **The peel memo key must carry the level, not just `(perm, stufe)`.** Added
  2026-07-30 from the implementation. `level = n − len(p) + 1` pins the level
  from the permutation length — but only *within one top-level call*, because
  `n` differs between permutations. A memo shared across `S_{132}` (n = 3) and
  `S_{1423}` (n = 4) returns a length-3 subtree computed in `x₁, x₂` to a
  caller expecting `x₂, x₃`. The bug is invisible to the obvious test: the
  single-permutation `expand → from_polynomial` round trip passed exhaustively
  through S₆, because each call built its own memo. Only `from_polynomial` on
  a *sum of permutations of different sizes* shares one memo across differing
  `n`, and it is the only test that failed. Any engine that memoizes the peel
  inherits this, E3 included — and it is the concrete form of the general
  warning in §3.5 that shift bookkeeping goes wrong in ways tests notice late.
- **Variable-index conventions are a minefield.** Sage's
  `multiply_variable(i)` is 0-based; Symmetrica's `mult_schubert_variable` is
  0-based; Symmetrica's `divdiff_schubert` is **1-based** (`sb.c`: the module
  is internally inconsistent). The first transition verification failed
  100% of cases on exactly this. §4 fixes one convention (1-based, matching
  the mathematics) and tests it at the boundary.

## 3. The algorithms

### 3.1 What the incumbent actually does (read today, public domain)

Condensed from `sb.c` / `mss.c` / `perm.c`; adoptable pieces marked ✓, defects
marked ✗.

- **Expansion** (`m_perm_schubert_monom_summe`, `sb.c:123`; worker
  `algorithmus2`, `sb.c:323`): iterated *inverse Monk at position 1* — peel
  x₁'s exponent, recurse on covers `w ⋖ w·t_{1,i}` via a running-minimum scan
  (✓ the scan is the right cover enumeration; reused in E1/E3). One leaf per
  pipe dream, so cost = S_w(1,…,1). ✗ No sharing between the leaves.
- **Product** (`mult_schubert_schubert`, `sb.c:965`): expand one factor to
  monomials, then for each monomial run |α| single-variable Monk passes over
  the whole accumulating sum (`mult_schubert_variable`, `sb.c:1003` — ✓ its
  two-scan signed cover enumeration is exactly the verified signed Monk rule).
  ✗ Cost = (#pipe dreams of one factor) × (Monk chains re-run per monomial),
  the two compounding blowups behind the S₁₀ wall. ✗ The "expand the smaller
  factor" heuristic compares only the *head term's* vector length.
- **From-polynomial** (`t_POLYNOM_SCHUBERT`, `sb.c:232`): greedy peel on the
  lex-smallest monomial, pad to a valid Lehmer code, subtract `c·S_w`, repeat
  (✓ triangularity verified today; adopted in §3.4). ✗ Re-expands S_w from
  scratch every iteration, no caching.
- **Divided differences**: on the Schubert basis `∂ᵢS_w = S_{wsᵢ}` or drop
  (`divdiff_schubert`, `sb.c:1334`) — ✓ trivially cheap, adopted. ✗ 1-based
  index disagreeing with the 0-based product routine; ✗ latent out-of-range
  read when a stored permutation is shorter than the letter applied.
- **Pairing** (`scalarproduct_schubert`, `sb.c:1840`): Poincaré pairing
  computed as full product followed by n(n−1)/2 divided-difference passes,
  where **n is read off the stored vector lengths** — the answer depends on
  how padded the inputs happen to be. ✗ Semantics a caller cannot predict;
  §4 takes n as an explicit argument and reads one coefficient instead.
- **Stanley/Schur** (`newtrans`, `mss.c:44` — the one well-engineered routine):
  the Lascoux–Schützenberger transition, iterative with an explicit stack,
  Grassmannian base case emitting a single s_λ, and the `F_w = F_{1×w}` shift
  when the left sum is empty (✓ the exact mechanics of §3.8, adopted). ✗ Static
  `char[1000]` state, not reentrant. ⚠️ Its Sage name `t_SCHUBERT_SCHUR` and
  `schur.doc`'s phrasing oversell it: it computes the **Stanley symmetric
  function** F_w, which equals S_w only in the stable range. `newtrans([2,1,4,3])
  = s₂ + s₁₁` (verified) while S₂₁₄₃ is not even symmetric. §4's name says
  what it is.
- ✗ **Nothing is memoized anywhere in the module**, every routine mutates its
  inputs in place (padding, exponent growth), and the identity permutation is
  represented as `[1,2]`. Stability is patched in the comparator (mixed-length
  lists compare fixed-point-padded) rather than by normal forms.

### 3.2 `Perm`: the index type

A `Perm` is a permutation of ℕ₊ fixing all but finitely many points, stored in
one-line notation **with trailing fixed points stripped at construction** —
stability as an invariant, exactly the move `Partition` makes for
weakly-decreasing (`partition.rs:14–52`), and the opposite of Symmetrica's
compare-time patch. Invariants: `w[i] ≥ 1`, values distinct, last entry not
fixed. The empty `Perm` is the identity.

Needed operations, all standard and all testable against Sage: `length()`
(inversions), `code()` (Lehmer) and `from_code()` (inverse bijection —
round-trip test), `descents()`, the two cover scans (right-multiplication by
t_{ij} increasing length by 1, enumerated by the running-max/min idiom
verified in Symmetrica's two routines), `is_dominant()`, `is_grassmannian()`,
`inverse()`. No reduced words are needed on any hot path (`∂_w` composes
one-letter steps; Macdonald's checksum enumerates reduced words only in tests,
via Sage).

### 3.3 The element type

```rust
pub struct Schubert<C: Ring> { terms: BTreeMap<Perm, C> }   // zeros never stored
```

mirroring `sym.rs:24`'s storage decision (deterministic iteration order — the
peel in §3.4 *requires* an ordered map). It is deliberately **not** a `SymFn`:
that trait is `Partition`-indexed (`sym.rs:26–32`) and everything generic over
it assumes Sym. `Forgotten` is the precedent for a type that opts out of a
trait rather than half-satisfying it (`sym.rs:184` — there it was
`SymAlgebra`); here the whole trait is declined and the handful of shared
helpers (`add_term`, `scale`, `format`) are reimplemented in the module. If a
third basis-indexed-by-something-else ever arrives, *then* extract the common
shape; not before.

### 3.4 Verified primitives

- **`mul_variable(i)`** — the signed Monk rule, 1-based, whole-element: one
  pass over the terms, two cover scans per term, merged into a fresh map.
  Verified form (S₄ × i exhaustive, today):
  `xᵢ·S_w = Σ_{j>i, cover} S_{w t_{ij}} − Σ_{j<i, cover} S_{w t_{ji}}`.
- **`divided_difference(i)`** — on the basis: swap-or-drop per term
  (`∂ᵢ S_w = S_{wsᵢ}` if w(i) > w(i+1), else 0). A second, independent
  implementation on *expanded polynomials* — `(f − sᵢf)/(xᵢ−xᵢ₊₁)` with exact
  division — exists only in tests, so the two can disagree loudly (§5.3).
- **`expand()`** — `algorithmus2`'s peel, shared-prefix: the recursion is a
  tree whose nodes are (remaining permutation, current variable index); nodes
  merge when equal. Even without merging it is the reference expansion; with
  merging it is E3's skeleton. Output is exponent-vector terms
  `Vec<(Vec<u32>, C)>`.
- **`from_polynomial(terms)`** — greedy peel on the lex-minimal exponent
  vector, justified by the verified triangularity; each subtracted `S_w`
  expansion comes from the (memoized) `expand`. This is `t_POLYNOM_SCHUBERT`
  minus the recompute-everything defect.
- **`transition(w) -> (r, v, Vec<Perm>)`** — the verified last-descent
  transition step; shared by E2, `stanley`, and tests.

### 3.5 The product — engine candidates

The crate's own history warns twice over here: *the obvious fix is often a
loss* (batched character sweep, `docs/record/README.md`), and *a structural idea that
benchmarks badly on its first data structure has not been tested* (three_row's
HashMap → dense-table 3.3×). So the spec commits to the reference route and
the cross-checks, states two engine candidates with the reasoning, and defers
the choice to measurement — targets in §6 either way.

**E1 — `NaiveSchubert`** (the in-house oracle, role of `NaiveLr`): Symmetrica's
route — expand the factor with the *smaller pipe-dream count* (estimated by
`dimension()`, which is `algorithmus3`'s substitution and costs one peel
sweep — fixing the head-length heuristic), then iterated `mul_variable`. Kept
forever, exhaustively agreed against whatever wins.

**E2 — memoized transition recursion.** Recurse on the product:

```text
    S_u · S_v:
      if ℓ(v) = 0            → S_u
      if v dominant          → x^{code(v)} · S_u        (iterated mul_variable)
      else (r, v', {v''}) := transition(v)
           S_u·S_v = x_r · (S_u·S_v')  +  Σ (S_u·S_v'')
```

memoized on the (u, ·) pair like `memo::product_cached` (`memo.rs:281`).
Every v'' has ℓ(v'') = ℓ(v) but the recursion terminates — verified
computationally (the transition identity held exhaustively on S₅ and the
recursion is exactly `newtrans`'s, whose termination Symmetrica has relied on
for thirty years; Macdonald's notes prove the tree finite). ⚠️ The recursion
tree on v is shared across every caller with the same v-side subproblems only
if the memo key is the *pair* — whether the hit rate justifies the table is a
measurement, not an argument.

**E3 — the frontier bet, now measured.** E1's two costs are (i) touching every
pipe dream of the expanded factor and (ii) re-running Monk chains per
monomial. But `algorithmus2`'s recursion tree shares prefixes massively, and
the Monk chain for a monomial x^α is the *same operator sequence* for every
monomial sharing a prefix of α. So: walk the peel DAG of the smaller factor
**once**, merging at equal states. That is precisely the `SkewLr` move —
"partial fillings that agree on the undecided part merge into one weighted
state" — transplanted from tableaux to Monk chains.

Read off `algorithmus2` (sb.c:323), a state is `(perm, alphabetindex, stufe)`
and the transitions are: `len(perm) = 2` → emit a monomial [LEAF];
`perm[0] = len(perm)` → `x_alphabetindex^stufe ·` recurse on `perm[1:]` at
`alphabetindex+1` [DESCEND]; otherwise sum over the Bruhat covers
`w ↦ w·t_{1,i}` at `stufe−1` [BRANCH]. Two facts make the merge cheap, and
both are why the state is effectively just the permutation:

- **`len(perm)` determines `alphabetindex`**, since DESCEND is the only
  transition that changes either, and it changes both by one.
- **`stufe` enters only as an overall factor.** It is *read* only at DESCEND,
  so the subtree value satisfies
  `V(p, ℓvl, stufe) = x_ℓvl^{stufe − s₀} · V(p, ℓvl, s₀)`. A node reusing a
  cached perm therefore pays a shift by a power of **one** variable — a Monk
  pass, not a re-walk.

So E3 evaluates the DAG bottom-up on states keyed by `perm` alone, each state
holding `V(state)·S_u` in the Schubert basis. Measured on the §2 cases
(`scripts/spec_schubert_peel.py`; leaf counts validated against Sage's
`expand()` exhaustively on S₂–S₆, 872 permutations, 0 mismatches):

```text
  case                     ℓ     pipe dreams      states  edges  peak live  compress
  stair4  (S_8)           10              64          88    113        12      0.7×
  stair5  (S_10)          15           1 024         283    382        28      3.6×
  stair6  (S_12)          21          32 768         923  1 286        69     35.5×
  stair7  (S_14)          28       2 097 152       3 052  4 340       182       687×
  stair8  (S_16)          36     268 435 456      10 192 14 693       499    26 338×
  rand S_11.0u            32         142 818         521    744        27       274×
  rand S_12.1u            27      15 468 012       5 092  7 706       279     3 038×
  rand S_13.0u            25      31 877 880       2 740  4 094       133    11 634×
  rand S_13.0v            36   1 355 962 209      10 796 17 006       546   125 599×
  w0(S_12), dominant      66               1          11     10         2       0.1×
```

`compress = pipe dreams ÷ states`, the direct analogue of the LR record's (`docs/record/littlewood-richardson.md`)
`LR tableaux ÷ states produced`. What the table says:

- **The leaf count is exponential and the state count is not.** On the
  staircase family pipe dreams are exactly `2^C(k,2)` while states grow ~3.3×
  per rung — the gap is the whole engine bet, and it is not marginal.
- **The compression is largest exactly where the incumbents die.** §2's
  `rand S_13 ℓ=25,36` is the pair the C `schubmult` does **not** finish in
  120s. Peeling its smaller factor is 31 877 880 pipe dreams against **2 740
  states**. That is a mechanism for §2's "cost does not track output size"
  observation, not just a correlation with it.
- **Memory is not the trade this time.** `peak live` is the maximum number of
  simultaneously-live state values under DFS post-order with refcounting (each
  live state holds a whole Schubert element). It never exceeds 546, including
  on the 1.4-billion-pipe-dream row. The LR work bought speed with memory and
  said so; on this evidence E3 does not have to, which makes the §6 RSS curve
  a confirmation rather than a disclosure. ⚠️ Element *sizes* are not measured
  here — only how many are live.
- **Compression below 1× exists and is bounded.** Dominant and short
  permutations have more states than pipe dreams. But their state counts are
  ≤ 68 in absolute terms, so the loss is a small constant, not a regime; E1
  stays the oracle and the dominant case already short-circuits (E2's second
  line).

⚠️ Still a node count, not a runtime. Per-state work is an element-sized Monk
operation, and E1's per-leaf work is also element-sized, so the ratio is
*indicative* of the speedup and not equal to it. The §6 ladder settles it.

**It was settled, 2026-07-30 — and the answer is split.** E3 is built
(`src/schubert.rs`, keyed on `(perm, level, stufe)`) and measured
(`examples/bench_schubert.rs`, `examples/profile_schubert.rs`):

```text
  case        compression   E1        E3        E3 vs E1   schubmult C   vs C
  stair4²           0.7×    0.0047s   0.0049s      1.0×       0.005s      1.1×
  stair5²           3.6×    0.1767s   0.0674s      2.6×       0.016s      0.24×
  stair6²          35.5×   27.7062s   2.4751s     11.2×       0.355s      0.14×
  stair7²           687×    (~8h est) 141.54s        —        20.6s       0.15×
  S_11.1 ℓ=35,26     50×          —    38.05s        —        0.748s      0.02×
```

- **The engine bet paid off against E1**, at roughly a third of the nominal
  compression (35.5× → 11.2×). The shortfall has two named causes: the
  `(perm, level, stufe)` key merges less than the `perm`-only key the table
  measures, and per-state work is element-sized where E1's per-leaf work is a
  monomial chain. So the ratio was a *real* signal and an *inflated* one —
  worth recording as the calibration for the next time a state-count ratio is
  used to predict a speedup.
- **The displacement bar is met.** `stair6²` and `stair7²` complete; Sage's
  engine completes neither. Where both run, `stair5²` is 6.85s → 0.067s, 102×.
- **The frontier bar is missed**, and in two different ways, which is the
  useful part:
  - On the staircase family it is a flat ~7× — and a sampling profile of
    `stair7²` (15 660 samples) puts **~47% of the time in `BTreeMap` and the
    allocator** (malloc/free 25%, `BTreeMap::insert`/`remove` 15%,
    memmove/bzero 7%), because `Schubert` stores `BTreeMap<Perm, C>` and every
    `Perm` is a heap `Vec`. Every Bruhat cover allocates; every key comparison
    chases a pointer. This is `three_row`'s lesson exactly — a sound
    structural idea benchmarked on the wrong data structure — and it is
    plainly worth 3–5×, which lands the staircase rows at the bar.
  - On `S_11.1` it is 51×, and that is **not** a constant-factor story. With
    118 822 output terms, E3's cost is (Monk passes) × (element size), and the
    element is the size of the answer for most of the DAG. schubmult produces
    those 118 822 terms in 0.748s; E3 performs on the order of 10⁹ term
    updates to do it. No data-structure work closes that.

**First round of engine-independent tuning, same day.** Two changes, one that
paid and one that did not — both recorded, because the one that did not is the
more useful entry:

```text
                                     stair5²   stair6²   stair7²   vs schubmult C
  E3 as first written                 0.0674s   2.4751s   141.54s      0.15×
  + inline Perm ([u8; 32], Copy)      0.0345s   1.3862s    81.12s      0.25×
  + in-place add + Rc-keyed memo      0.0565s   1.3540s        —       0.26×
  schubmult C                         0.016s    0.355s     20.6s        —
```

- **Inline `Perm` is worth ~1.8×** and confirms the profile's reading: with
  `[u8; 32]` + length, `Copy` and zero-filled so `Ord` is unchanged, the
  allocator falls from 25% to ~8% of samples and `BTreeMap` insert/remove from
  15% to ~4%. Below the 3–5× predicted, which is itself the lesson: the
  profile attributed 47% to map-and-allocator, but removing the *allocations*
  does not remove the *inserts*.
- **In-place `add_assign` plus an `Rc`-keyed memo bought ~2%, i.e. nothing.**
  The hypothesis was that `acc = acc.add(&part)` copying the accumulator per
  cover, and the memo cloning a whole element per hit, were the residue after
  the profile's `clone_subtree` line. They were not: with those gone the time
  stayed put, which locates the remaining cost in the **number of term
  insertions**, not the cost of each. A recorded dead end in the sense
  the record uses — the change is kept because it is strictly better code, but it
  bought nothing and no one should re-derive the idea expecting a win.

What that leaves, in order: the accumulate-then-sort representation
(`HashMap` + `fasthash::MixHasher` during accumulation, ordered only at the
boundary — `BTreeMap`'s ordering is needed by `from_polynomial`, not by the
engine), then merging on `perm` alone (1.3–1.8× by the §3.5 table). Together
those plausibly reach the staircase bar; neither touches the `S_11.1` regime.

### E2, built the same day — and it is the engine

Implemented as `Schubert::mul_e2` and verified against E1 exhaustively on
S₅ × S₅ in both argument orders (E2 recurses on the *second* factor, so the
orders walk different trees), plus staircases. The transition recursion
terminated everywhere; the depth cap in `mul_e2_depth` exists so that a
counterexample would be a clean panic rather than a hang.

```text
  case              terms      E3 nodes  E2 nodes  E2 passes   E3 time    E2 time   schubmult C   E2 vs C
  stair5²              59           641       191        333    0.0345s   0.0066s     0.016s        2.4×
  stair6²             247         2 559       633        931    1.3540s   0.0665s     0.355s        5.3×
  stair7²           1 111        10 052     2 085      2 667   81.1190s   1.4821s     20.6s        13.9×
  S_10.0            1 070           307       152        305          —   0.0125s     0.066s        5.3×
  S_11.0           30 143           181        91        204          —   0.0848s     0.222s        2.6×
  S_11.1          118 822           571     1 102      5 477   38.0513s   0.3923s     0.748s        1.9×
  S_11.2          185 284           265     1 156      2 064          —   0.7605s     3.86s         5.1×
  S_12.0          230 901         7 067     1 862      3 154          —   6.4134s    21.0s          3.3×
  S_12.1            8 280         1 787       556      1 117          —   0.6541s     1.61s         2.5×
  S_12.2          242 507           267     1 126      1 718          —   0.6221s    19.0s         30.5×
  S_13.1          947 254           734       351      1 368          —   1.8852s    20.9s         11.1×
  S_13.2        3 241 903         2 288       447        943          —   5.0792s   105.8s         20.8×
```

**The §6 frontier bar has two clauses, and E2 meets one of them.**

- ✅ *Within 2× of the C `schubmult` on every row it finishes.* E2 is **ahead
  on every such row**, by 1.9× to 30.5×. (`stair4²` is excluded — at 0.005s it
  is the startup floor, per §6.) The escalated `stair7²` clause of ≥1× is met
  at 13.9×.
- ❌ *Complete `S_13 ℓ=25,36`, which schubmult does not.* **E2 does not
  complete it either.** It runs 428s and dies at **6.56 GB** peak RSS.

The cause is not the mathematics, it is `E2`'s memo: `HashMap<Perm,
Rc<Schubert<C>>>` with **no eviction**, so every node's whole product is
retained for the entire run. On `S_13.2` that is 447 nodes against a 3.2M-term
answer; nodes holding even a fraction of that are hundreds of megabytes each.

This is precisely the failure §6 predicted — *"a rising RSS curve means the
implementation is retaining the whole DAG rather than refcounting it, and that
is a bug, not a trade"* — written about E3 and then shipped in E2 anyway.

**Fixed the same day, and the fix is free.** `count_transition_parents` walks
the transition tree using only `Perm` operations (no products), counting how
many parents each node has; `E2::release` drops a child's cached value the
moment its last parent has consumed it. Results:

- **No time regression anywhere.** Every ladder row is within noise of the
  pre-eviction numbers (`stair7²` 1.482s → 1.515s, `S_13.2` 5.079s → 5.099s).
  Eviction costs nothing because the pre-pass is `Perm`-only and the drops
  replace work the allocator would have done at the end regardless.
- **`S_13.2`, the largest row that finishes** — 3 241 903 output terms in
  5.27s — peaks at **1.52 GB**.
- **`S_13.0 ℓ=25,36` no longer dies, but still does not finish.** RSS now
  oscillates between 1.2 GB and 2.4 GB for 20+ minutes instead of climbing
  monotonically to 6.56 GB and aborting at 428s. So the memory bug is real and
  gone; what is left on that row is genuine time. **§6's second clause remains
  unmet**, now for a different and more honest reason.

### What the pathological row actually is

`ℓ=25,36` defeats both engines while `ℓ=40,41` — longer, 3.4× the output —
finishes in 5s. Since *both* implementations fail the same row and no other,
the difficulty has to be a property of the input pair, and therefore visible
without running a product. `examples/probe_schubert.rs` looks, using only
permutation-level quantities (`dimension` is `S_w(1,…,1)`;
`schubert::transition_tree` walks E2's own recursion building no elements, so
it costs microseconds even where the product is hopeless):

```text
  case             deg           dim(u)           dim(v)  tree(u)  tree(v)   E2 time
  stair6²           42            32 768           32 768      633      633    0.070s
  stair7²           56         2 097 152        2 097 152    2 085    2 085     1.51s
  S_11.1 ℓ=35,26    61            13 305           11 067      376    1 102     0.40s
  S_12.0 ℓ=37,24    61         3 352 856        1 215 900    1 384    1 862     6.58s
  S_12.2 ℓ=43,32    75             4 228        1 396 206      157    1 126     0.63s
  S_13.1 ℓ=36,41    77           496 776            9 310      856      351     1.90s
  S_13.2 ℓ=40,41    81           585 004          832 723      663      447     5.10s
  S_13.0 ℓ=25,36 *  61        31 877 880    1 355 962 209    1 477    4 909   NEITHER
```

Three conclusions, all measured rather than argued:

- **It is not the node count.** `S_13.0`'s transition tree is 1 477 nodes,
  *smaller* than `stair7²`'s 2 085 (which takes 1.51s) and comparable to
  `S_12.0`'s 1 384 (6.58s). E2 visits few nodes there; each is enormous. The
  difficulty is element size, confirming from the opposite direction the same
  thing the E2-vs-E3 result showed.
- **It is not the degree.** `deg = 61` is shared with `S_11.1` (0.40s) and
  `S_12.0` (6.58s), and is *lower* than `S_13.2`'s 81 (5.10s). Degree
  predicts nothing, which is why `ℓ` was a misleading label on these rows all
  along.
- **The all-ones specialization singles it out by 10⁴.** `S_u(1,…,1) ·
  S_v(1,…,1) = Σ_w c^w_{uv} S_w(1,…,1)` is an exact identity with every term
  non-negative, so `dim(u)·dim(v)` is the product's total monomial mass. It is
  **4.32×10¹⁶** for `S_13.0` against a ladder maximum of 4.40×10¹² everywhere
  else — the only quantity in the table that separates the row at all.

So `S_13.0` is not an algorithmic pathology: **its answer is astronomically
large**, and both engines fail it for the same reason. That also explains the
§2 observation it was filed under — "the C frontier's cost does not track
output size" — the cost tracks monomial mass, which `ℓ` and the Schubert-basis
term count both fail to reflect.

⚠️ Two limits on this, stated because the quantity is tempting to over-use.
`dim(u)·dim(v)` is **not** a predictor of runtime in general: `stair7²`
(4.40×10¹², 1.51s) and `S_12.0` (4.08×10¹², 6.58s) have the same mass and 4×
different times. And it does **not** predict the Schubert-basis term count —
mass ÷ terms ranges from 1.2×10³ to 4.0×10⁹ across the completed rows. What it
does is flag a row as out of family, and it costs microseconds, which makes it
the right thing for `mul` to consult before committing to a case. Whether to
*refuse* such a product, or attempt it and let the caller run out of memory,
is an API question (§7 Q9) rather than an engine one.

**Why it wins is not the node count, and that is the finding.** On `S_11.1`
E2 uses *more* nodes than E3 — 1 102 against 571 — and is still 97× faster.
Both engines cost (nodes) × (size of the running element); what differs is the
second factor. E1 and E3 expand a factor into monomials and push the other
through Monk chains, so the running element inflates toward the size of the
answer early and stays there. The transition recursion never expands: each
node holds a product that is only as large as that subproblem needs. The
element-size term, not the node count, was the whole gap to `schubmult` — and
the compression measurement, which only ever counted nodes, was structurally
incapable of seeing it. That is the sharper version of §3.9's rule: a cost
model that omits a factor will rank engines confidently and wrongly.

A sampling profile of E2 on `stair7²` (~13 400 samples) is healthier than
E3's ever was — the time is in the mathematics, not the plumbing:
`mul_variable` 50%, `add_assign` 11%, the cover scans 10%, `BTreeMap`
insert/remove/drop 8%, allocator 8%, memmove/memset 8%. Two specific
follow-ups, neither needed to hold the bar: the cover scans still allocate a
`Vec<(u32, Perm)>` per term per pass, and `BTreeMap::remove` alone is 504
samples — that is signed Monk's cancellation churning entries in and out,
which an accumulate-then-sort representation would absorb.

**Consequence for §7 Q1: E2 is back, and for a reason rather than for
symmetry.** The whole E1/E3 family expands one factor into monomials and
pushes the other through Monk chains, so its cost is tied to the size of the
running element. Buch's `schubmult` does not expand at all. A memoized
transition recursion is the shape that avoids the element-size term, and it is
now the *next engine to build*, not the fallback. E3 stays: it is 11× over E1,
it wins where output is small relative to the DAG, and a hybrid dispatching on
predicted output size is the likely end state.

**Leaf dispatch, both engines — verified, with a cost caveat.** Products where
both factors are Grassmannian with the same descent k reduce to LR with a
k-row flag truncation: `AutoLr::schur_product` (`strip_lr.rs:281`) followed by
dropping output partitions with more than k rows and re-indexing as
Grassmannian `Perm`s. **Verified 2026-07-30** (`spec_schubert_grass.py`):
exhaustive over same-descent pairs for k ≤ 4, |λ|,|μ| ≤ 5 (760 pairs, 0
failures), nine larger cases through k = 5, and the shape ↔ permutation
convention `code(w) = (λ_k,…,λ_1)`, `w(i) = λ_{k+1−i} + i` pinned against
`s_λ(x₁..x_k)` over 128 shapes. Independent corroboration: the surviving term
counts are 2, 5, 16, 59, 247, 1111 for k = 2..7, and 16/59/247/1111 are
exactly §2's measured `stair4²`/`stair5²`/`stair6²`/`stair7²` product sizes.
A different-descent pair produces non-Grassmannian output, so "same descent"
is load-bearing rather than decorative.

The caveat the identity's truth conceals is what it costs to use:

```text
  k       LR terms   kept (≤k rows)   dropped   waste
  2              7                2         5   71.4%
  3             34                5        29   85.3%
  4            206               16       190   92.2%
  5          1 433               59     1 374   95.9%
  6         10 873              247    10 626   97.7%
  7         87 452            1 111    86 341   98.7%
```

`AutoLr::schur_product` has no row bound, so dispatching a `stair7` leaf to it
computes 87 452 terms to keep 1 111 — 79× redundant work, which would eat the
whole point of dispatching to the fast engine. So leaf dispatch is contingent
on a **row-bounded product**, `schur_product_bounded(μ, ν, max_rows)`:
`SkewLr`'s recursion adds horizontal strips row by row (`strip_lr.rs:47`), so
a partial shape already exceeding `max_rows` can be cut where it is generated
rather than filtered at the end. That is a new, small piece of LR work this
spec now depends on (§7 Q6) — it is not a Schubert change. Until it exists,
leaves stay on E1/E2 and nothing else moves.

**Not chosen: cotransition.** Knutson's co-transition formula (1909.13777,
the Lemma on p.2) recurses *upward* — `(xᵢ − y_{π(i)}) P_π = Σ P_σ` over
covers — and computing S_π from it means dividing a sum of Schubert
polynomials by a linear form. `Ring::div_exact` (`coeff.rs:116`) is the seam
it would need. The 2026 specialization paper (arXiv:2603.20104) compares
descent/transition/cotransition recursions *for specializations*; for whole
products nothing here needs division, so the division-shaped route is recorded
and skipped, per the standing rule that being able to do a thing is not a
reason to.

### 3.6 Coefficients, overflow, and the boundary

Structure constants are non-negative integers, small at every measured size
(≤16). The engine runs generic over `Ring` exactly like `Schur::mul_with`
(`sym.rs:244`), coefficients injected via `C::from_u128`. The Python entry
points wrap in `guarded` + `escalate` (`guard.rs:57`, `python.rs:222`)
unchanged — intermediates in E2/E3 are signed (Monk's left scan subtracts), so
`Guarded`'s checked arithmetic is doing real work even while the answers stay
tiny. Permutations cross the FFI as one-line `Vec<u32>` with values 1..n,
normalized on entry exactly as `part()` re-normalizes partitions
(`python.rs:111`).

### 3.7 The single-coefficient query — built, and it is the answer for the
### out-of-family cases

*Rewritten 2026-07-30. This section was "candidate, not plan"; the
characterisation of `S_13.0` promoted it, because for a product whose answer
cannot exist in memory a single coefficient is the only question left that has
an answer.*

`schubert_coeff(u, v, w)` runs E2 with **Bruhat pruning**: at every node,
terms not `≤ w` are discarded. Sound because the signed Monk rule moves
strictly *up* the Bruhat order — both of its sums run over covers — so nothing
below the cut can reach `w`, directly or by cancelling against something that
does. Two necessary conditions are checked first and answer most queries with
no work at all: `ℓ(w) = ℓ(u) + ℓ(v)`, and `u ≤ w`, `v ≤ w`.

All three premises are *verified*, not assumed
(`product_support_lies_above_both_factors`,
`monk_covers_move_up_the_bruhat_order`,
`bruhat_is_a_partial_order_with_known_bounds`), and the query itself agrees
with reading the coefficient off the full product exhaustively on S₄ × S₄ with
targets ranging over S₆ — **including every zero answer**, since a pruning bug
shows up as a spurious zero rather than a wrong number.

```text
  S_13.2 ℓ=40,41   full product, 3 241 903 terms   5.08s
                   one coefficient (c = 130)       0.0088s      578×
  S_13.0 ℓ=25,36   full product                    IMPOSSIBLE (mass 4.3×10¹⁶)
                   one coefficient                 0.028–0.190s
```

**The row no engine can complete is queryable in a tenth of a second.** That
is the practical answer to §3.5's out-of-family case, and it is what the use
case in §1 — positivity searches, rule-hunting — actually needs: nobody wants
2.9×10¹¹ terms, they want particular ones.

**Nonzero constants exhibited on the impossible pair** (`examples/find_nonzero.rs`).
Random `w` are essentially never in the support, so sampling the group only
ever returned 0. Walking the support works: **propose** by applying the Monk
chain of `x^{code(v)}` to `S_u` with a beam of 4 000, which yields candidates
automatically of the right degree (each pass raises `ℓ` by one and
`Σ code(v) = ℓ(v)`) and automatically `≥ u`; **dispose** by re-checking each
one with the exact pruned query. The beam only decides where to look — its own
coefficients are wrong, since truncation discards terms that would have
cancelled.

```text
  S_13.0 ℓ=25,36   3 545 candidates in 0.06s, 2 238 with v ≤ w
    c^[8,16,12,10,4,6,5,7,13,2,3,11,1,9,14,15] = 5
    c^[8,16,12,10,4,6,5,7,14,2,3,9,1,11,13,15] = 5
    c^[8,16,12,13,4,6,5,7,10,2,3,9,1,11,14,15] = 5
    c^[8,16,13,10,4,6,5,7,12,2,3,9,1,11,14,15] = 3
    c^[9,16,12,10,4,6,5,7,13,2,3,8,1,11,14,15] = 18
```

**Refereed on the computable neighbour**, which is the only way to trust
answers on a pair that by construction has no referee: the same pipeline run
on `S_13.2` and checked against its full 3 241 903-term product agrees on
**400 candidates with 0 mismatches, 214 of them nonzero**. The zeros matter as
much as the nonzeros there — a pruning bug shows up as a spurious zero, not a
wrong number.

So the claim is now the strong one: **structure constants of a product that
cannot exist in memory, computed and checked.** `c^w_{uv} = 18` above is a
number no other package can produce by any route.

### 3.7.1 Superseded sketch (kept for the reasoning)

`c^w_{uv}` without the whole product. Three shapes, none validated:
read one coefficient off the (memoized) E2/E3 product — the baseline;
`c^w_{uv} = c^w_{vu}`, so peel the factor with the smaller transition tree
(the `lr_coeff` lesson — peeling the wrong side was a 21× loss there);
or extract via the pairing, `c^w_{uv} = [S_{w₀}] (S_u·S_v·S_{w₀ w w₀...})` —
degree bookkeeping makes this a real formula only in the right range, and it
computes a *bigger* product, so it is probably the dead end it looks like.
Deliverable is the honest measurement, after the engine exists. No package
has any such query; even a 2× win over whole-product would be new.

### 3.8 `stanley(w)`: the bridge into Sym

The transition recursion with the x_r-term dropped and the `1×w` shift when
the left sum is empty computes the Stanley symmetric function F_w expanded in
Schur functions — `newtrans`'s semantics (verified today:
`newtrans([2,1,4,3]) = s₂ + s₁₁`), landing in the existing `Schur<C>`. This is
the one function whose output lives in Sym, it is what Sage's
`t_SCHUBERT_SCHUR` misleadingly calls Schubert→Schur, and it gives the
module a second connection to the rest of the crate (vexillary w → a single
s_λ — a test; and F_w's Schur positivity — a law test). Iterative with an
explicit stack, per the incumbent, minus the static buffers.

### 3.9 A recorded dead end

The first benchmark family (§2's rectangle trap) *was* the dead end this time:
half a day of "hard cases" that were dominant monomials in disguise, caught
only because 1-term outputs are conspicuous. The rule it re-teaches is the
`st` spec's: compute the thing in Sage and *look at it* before building
anything on it. It is why every formula in this document carries a
verified-today marker or an explicit ⚠️ unverified flag. **As of 2026-07-30
the count of unverified flags is zero** — the last one (Grassmannian flag
truncation) went through its Sage session and survived, and the flags that
remain in §3.5 and §6 mark *predictions about cost*, which only running code
can settle, not unchecked mathematics.

---

## 4. API

Following the house shape — one module, whole-object operations, conversions
where they mathematically exist.

```rust
// permutation.rs (new; the index type is not Schubert-specific)
pub struct Perm(/* one-line, trailing fixed points stripped */);
impl Perm {
    pub fn new(one_line: impl IntoIterator<Item = u32>) -> Result<Perm, PermError>;
    pub fn from_code(code: &[u32]) -> Perm;          // Lehmer inverse
    pub fn code(&self) -> Vec<u32>;
    pub fn length(&self) -> u32;                     // ℓ(w) = |code|
    pub fn is_dominant(&self) -> bool;
    pub fn is_grassmannian(&self) -> Option<u32>;    // Some(descent)
    pub fn inverse(&self) -> Perm;
}

// schubert.rs (new)
pub struct Schubert<C: Ring> { /* BTreeMap<Perm, C> */ }

impl<C: Ring> Schubert<C> {
    pub fn monomial(w: Perm, c: C) -> Self;
    pub fn grassmannian(lambda: &Partition, descent: u32) -> Self;   // S_w = s_λ(x_1..x_k)
    pub fn mul(&self, other: &Self) -> Self;                          // the engine
    pub fn mul_variable(&self, i: u32) -> Self;                       // x_i·f, 1-based, signed Monk
    pub fn divided_difference(&self, i: u32) -> Self;                 // ∂_i, 1-based
    pub fn divided_difference_perm(&self, w: &Perm) -> Self;          // ∂_w composed
    pub fn expand(&self) -> Vec<(Vec<u32>, C)>;                       // monomials, exponent vectors
    pub fn from_polynomial(terms: &[(Vec<u32>, C)]) -> Self;          // greedy peel (§3.4)
    /// Poincaré pairing ⟨f,g⟩ on H*(Fl(n)) — n EXPLICIT, unlike the incumbent.
    pub fn pairing(&self, other: &Self, n: u32) -> C;
}

pub fn schubert_coeff<C: Ring>(u: &Perm, v: &Perm, w: &Perm) -> C;    // §3.7
pub fn stanley<C: Ring>(w: &Perm) -> Schur<C>;                        // §3.8
pub fn dimension(w: &Perm) -> Option<u128>;            // S_w(1,…,1), #pipe dreams
pub fn principal_specialization_q(w: &Perm) -> Vec<i128>;             // S_w(1,q,q²,…)
```

**Python bindings — built 2026-07-30**, ten entry points, all through
`guarded`/`escalate`, permutations as 1-based one-line lists normalized on
entry exactly as `part()` normalizes partitions:
`schubert_multiply`, `schubert_multiply_variable`,
`schubert_divided_difference`, `schubert_divided_difference_perm`,
`schubert_expand`, `polynomial_to_schubert`, `schubert_pairing(a, b, n)`,
`schubert_dimension`, `schubert_coefficient`, `schubert_monomial_mass`.

Checked by `scripts/check_schubert_bindings.py` against Sage: **415 checks, 0
failures** — products, the 1-based/0-based `multiply_variable` boundary,
divided differences, `expand` and its round trip, `dimension`, every
single-coefficient value over S₅ for three pairs, and stability under padded
input. The harness is separate from the engine tests on purpose (§5.7): every
one of those already passed, so a failure here can only be marshalling. It
earned its keep immediately — four apparent failures turned out to be the
*harness* looking up a padded key in a dict Sage keys unpadded, which is
exactly the class of fault it exists to localise.

End to end from Python: `S_13.2` is 3 241 903 terms in 6.0s against the C
`schubmult`'s 105.8s; `schubert_monomial_mass` reports 4.323×10¹⁶ for
`S_13.0` and `schubert_coefficient` answers on that same pair in 0.038s.

**`stanley()` and `schubert_to_stanley_schur` are built too**, closing §4 —
the binding check is now **427 checks, 0 failures**, including agreement with
Symmetrica's own `newtrans` on 11 permutations. Tests: the §3.8 hand value
(`F_{2143} = s₂ + s₁₁`), Grassmannian `w` giving exactly one Schur term of the
shape `grassmannian_perm` was built from, Schur-positivity, correct degree,
and invariance under the `1×w` shift — a law test for the one move in the
recursion that is not a transition step.

**Against `newtrans`: parity, then 1.9× from removing allocations.** As first
written it was 2.7× at S₁₄ and **0.7× at S₁₆ — slower** — which is the
expected outcome, since §3.1 calls `newtrans` the one well-engineered routine
in the module. Two measurements fixed that:

- **Memoization would have bought nothing, and was not written.** The obvious
  fix was the defect §3.1 records for Symmetrica — `stanley` walked the
  transition tree as a *tree*. But `examples/probe_stanley.rs` measured the
  sharing first: **1.0×–2.1×, mostly 1.0–1.3×**. The Stanley transition tree
  really is a tree. Recorded as a dead end that cost one measurement instead
  of an implementation — the discipline E3 taught the same day.
- **Marshalling was not the cost either**: 0.00241s in Rust against 0.00287s
  through the binding on the same input, so the FFI is 16%.
- What was left is ~180 ns/node, and per node `transition()` heap-allocated
  **three** times: `descents()` built a whole `Vec` for its last element, and
  `covers_left` built a `Vec<(u32, Perm)>` that callers re-collected into a
  `Vec<Perm>`. `Perm::last_descent` and `Perm::for_each_cover_left` remove all
  three, with the ups pushed straight onto the work stack. **1.9× uniformly**
  (S₁₈ 0.00241s → 0.00133s).

Now ahead of `newtrans` everywhere measured: **1.3× at S₁₆** (was 0.7×), 4.4×
at S₁₄, 2.5× at S₁₂. ⚠️ The first row's 205× is Sage warmup, not a result —
its very next call costs 0.0003s.

Original binding plan (whole-object rule, `python.rs` module list at :882),
mirroring the seven `sage.libs.symmetrica` Schubert entry points plus the
from-polynomial pair: `schubert_multiply`, `schubert_multiply_variable`,
`schubert_expand`, `polynomial_to_schubert`, `schubert_divided_difference`
(int and permutation forms), `schubert_pairing(a, b, n)`,
`schubert_to_stanley_schur`, `schubert_dimension`. All escalate through
`guard`; permutations as 1-based one-line lists.

**Where it lives**: a module in symfn, not a sibling crate. Reasons: it reuses
`Ring`/`guard`/`memo`/`python` wholesale and `Schur`/`AutoLr` at the leaves;
the wheel must keep serving Sage as one artifact for the backend-substitution
story; and the crate boundary would cut exactly through `stanley` and the leaf
dispatch, the two places the domains genuinely touch. The alternative (a
`schubert` crate depending on symfn) is real and revisitable the day someone
wants Schubert without Sym — nothing in the module layout below forecloses it.
`lib.rs` gains `pub mod permutation; pub mod schubert;` in the alphabetical
block (lib.rs:52–82) with re-exports alongside (lib.rs:84–117).

### Crate fit

| need | existing | file |
|---|---|---|
| index newtype w/ invariant at construction | the `Partition` pattern | `partition.rs:14–52` |
| term storage, zero-elision, ordered iteration | `BTreeMap` per `SymFn` | `sym.rs:24` |
| structure constants via `from_u128` injection | `Schur::mul_with` shape | `sym.rs:244` |
| LR at Grassmannian leaves | `AutoLr::schur_product` — **needs a `max_rows` variant**, §7 Q6 | `strip_lr.rs:281`, strip recursion `:47` |
| memoized whole products | `table!` + `product_cached` shape | `memo.rs:42, :281` |
| overflow guard + bignum re-run at the boundary | `guarded` / `escalate` | `guard.rs:57`, `python.rs:222` |
| coefficients as Python ints, never strings | `Coeff` enum | `python.rs:65` |
| Stanley output | `Schur<C>` + `convert` hub | `sym.rs:153`, `convert.rs:44` |
| S_w(1..1), S_w(1,q,..) closed forms | `eval.rs` precedent (`dimension`, `principal_specialization_q`) | `eval.rs:359, :440` |

Deliberately **not** reused: `SymFn`/`SymAlgebra` (Partition-indexed; §3.3),
`LrBackend` as a trait bound (the engine is not an LR backend; it *calls* one
at leaves), `QtPoly` (single Schuberts need no second alphabet — double
Schuberts would, and that is the v2 seam), and `convert.rs`'s hub (Schubert is
not a basis of Sym; the only bridges are `grassmannian` and `stanley`, both
explicit functions rather than `ToSchur` impls, so the type system never
suggests a conversion that does not exist).

## 5. Correctness requirements

Layered as the crate does, each layer catching what the previous cannot:

1. **Hand values.** `S_{132} = x₁+x₂`, `S_{321} = x₁²x₂`, `S_{1423} =
   x₁²+x₁x₂+x₂²` (Knutson p.6), the dominant and Grassmannian families on
   small cases, identity/unit laws.
2. **Structural laws.** `code ↔ Perm` round trip; stability (constructing from
   padded one-line input yields the same element); `deg = ℓ`; `∂ᵢ∂ᵢ = 0`; the
   braid relations on ∂; `∂ᵢ S_w = S_{wsᵢ}`-or-0 — with the basis-side and
   expanded-polynomial-side implementations of ∂ compared against each other
   (they share nothing).
3. **Defining recursion.** `S_{w₀⁽ⁿ⁾} = x^δ` and downward ∂-recursion
   reproduces `expand()` for all of S₅.
4. **Engines agree.** E1 ≡ E2 (≡ E3 if built) on all products in S₄ and S₅,
   random S₇/S₈ pairs; Monk and transition each re-derived from the other's
   output on the same sweep. Grassmannian leaf dispatch against E1 on every
   same-descent pair through S₈ — the truncation statement has now passed its
   Sage check (§3.5), so this is a regression test for the implementation
   rather than a test of the mathematics.
   Plus the two checks that came out of the 2026-07-30 measurements:
   `expand()`'s leaf count equals `S_w(1,…,1)` exhaustively on S₂–S₆ (the
   validation `spec_schubert_peel.py` already runs, 872 permutations), and
   **E3's merged DAG agrees with the unmerged tree** on small cases — the
   merge is the one place E3 can be subtly wrong while staying plausible, and
   `peel_stats_nomemo` exists precisely to catch it (`[1,4,2,3]`: 15 tree
   nodes, 10 merged states, 3 leaves).
5. **External oracles, committed fixtures.** Sage/Symmetrica for products and
   expansions below its wall (≤ S₉); **lrcalc's `schubmult` for large ones** —
   generated by a script driving the binary out of process, exactly the
   `tests/lrcalc_oracle.rs` pattern including the licensing rationale in
   `NOTICE.md` (GPL program invoked as an external oracle, output committed,
   no code read). `stanley` against `newtrans` on a sweep including
   non-vexillary and non-stable-range cases.
6. **The scale checksum.** Macdonald's reduced-word identity plus
   `S_w(1,…,1) = Σ coeffs of expand` at sizes where no third-party oracle
   finishes, with the ±1-perturbation negative control — a checksum that
   silently always passed would be worse than none
   (`examples/verify_specialization.rs` sets the standard: 412/412 detected).
7. **Bindings separately** (`scripts/check_bindings.py` pattern): a correct
   answer in the wrong slot is a different failure from a wrong answer.

## 6. Performance requirements

Stated before the code exists, so the measurement can embarrass them:

- **Displacement bar**: complete `stair6²` and every seed-1 S₁₁/S₁₂ row of §2
  (Sage's current engine completes none). Through-Sage end-to-end comparison
  against the Symmetrica backend on the S₉-and-below range where both run,
  like-for-like per the backend-shim methodology already in `docs/record/README.md`.
- **Frontier bar**: the C `schubmult` rows of §2. Target: within **2×** of it
  on every row it finishes, and **complete `S_13 ℓ=25,36`**, which it does
  not (in 120s). The compression measurement has now been done, and it
  triggers the escalation clause this bullet carried as a guess: `stair7`'s
  ratio is **687×**, far past the "~10×" trip point, so the target on that row
  is raised to **≥1× schubmult C** — beat 20.6s, do not merely approach it.
  On `S_13 ℓ=25,36` (11 634× on the smaller factor) mere completion is now the
  *weak* reading of the bar; if E3 is built and that row is not comfortable,
  the compression did not convert and that is the finding to report.
- **Curves, not points**: the §2 ladder re-run per engine, plus peak RSS
  (`RSS=1` harness pattern) — the LR work bought speed with memory and said so;
  this must too. Predicted here: E3's live frontier is ≤ 546 Schubert elements
  on every §2 case (§3.5), so the RSS curve should come out *flat* in the
  state count and track element size instead. A rising RSS curve means the
  implementation is retaining the whole DAG rather than refcounting it, and
  that is a bug, not a trade.
- **Single-threaded numbers only**, per the standing rule: schubmult C is
  single-threaded, and a parallel symfn number set against it would conflate
  an algorithmic claim with a hardware one.
- **Subtract the out-of-process floor, and know it flatters us.** `schubmult`
  is driven as a subprocess, so its timings bundle work symfn does not do.
  Measured 2026-07-30 (`scripts/spec_schubert_floor.py`, min of 5):

  ```text
    startup floor (trivial input, 1 term out)        ~3.0–3.3 ms
    printing: pipe vs /dev/null, 30 143 terms              +1.0%
    printing: pipe vs /dev/null, 118 822 terms             −0.6%
    printing: pipe vs /dev/null, 185 284 terms             −0.4%
  ```

  Two consequences, opposite in sign:

  - **Output formatting is not a confound.** This was the one worth worrying
    about, because it would scale with term count — the very axis the engine
    comparison rides on, and a "win" over `printf` would be worthless. At
    185 284 terms it is inside the noise. The comparison axis is clean.
  - **Startup is ~3 ms of every row, and it inflates the incumbent.** So the
    §2 schubmult rows at or below ~0.03s — `stair4²` (0.005s), `stair5²`
    (0.016s), and both fast random S₁₀ rows — are substantially measuring
    `exec` and are **not valid targets**; beating them proves nothing about
    the engine. Use rows ≥ ~0.1s, or subtract the floor explicitly. This cuts
    against symfn, which is exactly why it has to be stated: an in-process
    Rust call compared against a process launch wins the small rows for free.
- **Report where it stops.** ✅ Done 2026-07-30, and it is **S₁₅**
  (`examples/bench_schubert_wall.rs`). S₁₄ is routine — 2.6s, 52.9s and 83.4s
  for 0.37M, 7.1M and 6.7M terms — while the first S₁₅ pair takes **421.9s**
  for 12.4M terms. Recorded in `docs/record/README.md` alongside the other subsystems.
  Coefficient growth on the same run settles §7 Q4: **130** at S₁₃, **591** at
  S₁₄, **863** at S₁₅, against the "≤16" this document originally carried from
  Symmetrica's easy rows.

## 7. Open questions

1. **E2 vs E3 vs hybrid.** *Answered 2026-07-30, then re-opened the other way
   by the same day's measurements.* Sequence, kept because the reversal is the
   lesson: the compression measurement pointed hard at E3 (687× on `stair7`),
   so E3 was built first and **E2 was demoted to "must justify itself"**. E3
   then came in 11.2× over E1 — the bet was real — but 7× *behind* the C
   `schubmult` on staircases and 51× behind on a large-output random pair. The
   node-count ratio converted at about a third, and the residue is not a
   constant: E1 and E3 are both `expand × Monk`, so both carry an element-size
   factor that `schubmult` does not pay.

   So: **E2 is the next engine**, on the merits. Open sub-questions —
   (a) does a memoized transition recursion actually avoid the element-size
   term, or does it just relocate it? (b) what is the pair-key hit rate
   (Q4)? (c) is the end state a hybrid dispatching on predicted output size,
   with E3 taking the small-output/large-DAG cases it already wins? Measure
   (a) on `S_11.1 ℓ=35,26`, the row where E3 is worst, before building
   anything on top of it.
2. **The hot representation.** `Schubert` is `BTreeMap<Perm, C>` with
   `Perm(Vec<u32>)`, so every Bruhat cover heap-allocates and every key
   comparison chases a pointer; a `stair7²` profile puts ~47% of the time in
   the map and the allocator (§3.5). An inline fixed-size `Perm` (`[u8; 32]`
   plus a length, `Copy`, n ≤ 32 — the ladder's worst case is ~25) removes the
   allocation without changing any signature except `one_line`. Worth 3–5× on
   the profile's evidence, and it is engine-independent: E2 will want it too,
   which is a reason to do it **before** E2 rather than after. ⚠️ `three_row`'s
   warning applies in both directions — measure it, do not assume it.
2. ~~**Does the Grassmannian flag-truncation identity hold as stated?**~~
   **Closed 2026-07-30: it holds**, exhaustively for k ≤ 4 plus larger spot
   cases, with the indexing convention pinned and a negative control (§3.5).
   The follow-on is Q6, not this.
4. **Coefficient growth.** ✅ **Answered.** 130 at S₁₃, 591 at S₁₄, 863 at
   S₁₅ — growing steadily, still four orders from `i64`, and `guard` covers the
   rest. ⚠️ **The ≤16 figure was wrong** — it came from §2's
   Sage/Symmetrica rows, which only ever reached the products Symmetrica could
   finish. Our own engine finds **c = 130** in `S_13.2 ℓ=40,41`. Still nowhere
   near `i64`, but the trend is now measured on the frontier rather than on
   the incumbent's easy cases, and it should be re-checked at S₁₄–S₁₅. One more rung (S₁₄–S₁₅ via our own
   engine once it exists) before concluding anything; `guard` covers the
   meantime.
5. **Memo policy for pair-keyed products.** `product_cached` works on
   (Partition, Partition); permutation pairs are a much bigger key space with
   much lower hit rates. Measure the hit rate before paying the table.
6. **Double Schubert polynomials** (two alphabets; Symmetrica has them as a
   copy-paste second recursion, `sb.c:1451`; Samuel's package centers them).
   The v2 item, and the real reason `schubmult`-the-package exists — parity
   there means a two-alphabet coefficient story (`QtPoly`'s shape, different
   ring). Out of scope until single Schuberts hold the §6 bars.
7. **A row-bounded LR product** (new, 2026-07-30, and the only piece of this
   spec that asks for work *outside* the Schubert module).
   `schur_product_bounded(μ, ν, max_rows)` on `SkewLr`, pruning partial shapes
   that exceed `max_rows` inside the strip recursion (`strip_lr.rs:47`) rather
   than filtering the finished expansion. Without it the verified leaf
   dispatch does 79× redundant work at `stair7` (§3.5). Questions: does the
   pruning actually cut the recursion, or does `SkewLr`'s state merging mean
   the excess rows cost little to begin with? Does `rect`/`two_row`/
   `three_row` dispatch survive a row bound, or does bounding force everything
   down the general path and lose more than it saves? Measure before building
   the dispatch — this is a small LR-side experiment, and if it fails the
   answer is simply that leaves stay on E1/E3.
9. **Refuse, or attempt and die?** (new, 2026-07-30.) `dim(u)·dim(v)` is the
   product's total monomial mass, costs microseconds, and flags out-of-family
   inputs like `S_13 ℓ=25,36` (4.32×10¹⁶ against a ladder maximum of
   4.40×10¹²) whose answers cannot be materialised. So `mul` *can* know in
   advance that a case is hopeless. Whether it should refuse, warn, or
   attempt it anyway is an API decision and not obviously "refuse": the
   threshold is soft (mass is not a runtime predictor, §3.5), a caller with
   a big machine may legitimately want the attempt, and silently declining a
   product a competitor also cannot do is a worse failure than running out of
   memory loudly. Leaning toward: compute it, expose it
   (`schubert::monomial_mass`), refuse nothing.
10. **Type B/C/D** (`bar.c` shows type B living entirely in expanded-polynomial
   land, sharing nothing) and **Grothendieck/K-theory** (Knutson §7: same
   recursions, isobaric operators, `1−exp` base case): recorded, out of scope,
   and the module layout should not have to move for either.
