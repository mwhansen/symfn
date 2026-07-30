# Specification: Jack polynomials (`jack`)

*Written 2026-07-29. Implements §2.6 of `docs/research-gaps.md` — "Jack polynomials
at scale": "Sage times out at n=15. The Goulden–Jackson matchings-Jack
conjecture and the b-conjecture are open and computationally starved. Routes:
Knop–Sahi, the Lassalle recurrences, or the Laplace–Beltrami eigenoperator,
rather than Gram–Schmidt." All three named routes are dealt with below; the
Laplace–Beltrami one wins.*

Sources, all read today, 2026-07-29 — the PDFs, not recollections:

| tag | source | used for |
|---|---|---|
| **[KS]** | Knop, Sahi, *A recursion and a combinatorial formula for Jack polynomials*, [arXiv:q-alg/9610016](https://arxiv.org/abs/q-alg/9610016) (Invent. Math. 128, 1997) | Thm 1.1 (positivity **and** the `u_μ` divisibility); Thm 4.6 (nonsymmetric recursion); Thm 5.1 (tableau formula, `d_λ(s) = α(a+1)+l+1`); the α-specialization list |
| **[MOPS]** | Dumitriu, Edelman, Shuman, *MOPS: Multivariate orthogonal polynomials (symbolically)*, [arXiv:math-ph/0409066](https://arxiv.org/abs/math-ph/0409066) (J. Symb. Comput. 2007) | the Laplace–Beltrami eigenvalue `ρ^α_κ`; the moving-box recursion; Lemma 2.16 (denominator never vanishes) |
| **[GJ]** | Goulden, Jackson, *Connection coefficients, matchings, maps and combinatorial conjectures for Jack symmetric functions*, TAMS 348 (1996) 873–892 (scanned original, pp. 873–875 read today) | series (1) `Φ`, (2) `Ψ`; coefficient definitions (4), (5); both conjectures verbatim |
| **[BDD]** | Ben Dali, Dołęga, *Positive formula for Jack polynomials, Jack characters and proof of Lassalle's conjecture*, [arXiv:2305.07966](https://arxiv.org/abs/2305.07966) | Lassalle's 2008 conjecture is now a theorem; positive p-expansion via bipartite maps |
| **[BD]** | Ben Dali, *Integrality in the Matching-Jack conjecture and the Farahat–Higman algebra*, [arXiv:2203.14879](https://arxiv.org/abs/2203.14879) | integrality half of Matchings-Jack is settled |
| **[DF]** | Dołęga, Féray, *Gaussian…/Cumulants of Jack symmetric functions and the b-conjecture*, [arXiv:1601.01501](https://arxiv.org/abs/1601.01501) (+ their 2016 polynomiality paper) | ℚ[b]-polynomiality of both [GJ] coefficient families |

Two foundational sources were **not** fetched: Macdonald's book (VI.10) and
Stanley's 1989 *Adv. Math.* paper. Every formula below attributed to that
tradition was instead **verified numerically against Sage before being written
down**, by `scripts/spec_jack_verify.py`, which is committed alongside this
document and prints `all formulas verified`. The walls are
`scripts/spec_jack_walls.py`; the arithmetic simulation is
`scripts/spec_jack_swell.py`. **Sage is an oracle here, never a source** —
SageMath 10.9 (2026-05-04), this machine, today; its Jack implementation was
not read and must not be, the same clean-room posture as
`docs/record/macdonald-operators-spec.md`.

---

## 1. What the thing is

Jack polynomials `P_λ(x; α)`, indexed by partitions, one real parameter α.
Three characterizations, all exercised below: (a) Gram–Schmidt against
dominance in the deformed Hall pairing; (b) eigenfunctions of the
Sekiguchi–Debiard / Laplace–Beltrami operators; (c) [KS]'s recursion or
tableau formula, either of which "could immediately serve as a *definition*".

```text
  ⟨p_λ, p_μ⟩_α = δ_λμ · z_λ · α^{ℓ(λ)}                     (verified n ≤ 6)
  P_λ = m_λ + Σ_{μ < λ} c_λμ(α) m_μ,   ⟨P_λ, P_μ⟩_α = 0  (λ ≠ μ)
```

### 1.1 Hooks, and the three normalizations

For a cell `s` with arm `a` and leg `l` (0-based, as `macdonald.rs:214/219`):

```text
  lower hook   h_low(s) = α·a(s)     + l(s) + 1
  upper hook   h_up(s)  = α·(a(s)+1) + l(s)
  H_λ  = ∏_s h_low(s)        H'_λ = ∏_s h_up(s)

  P_λ  monic in m_λ
  Q_λ  = (H_λ / H'_λ) · P_λ          so that ⟨P_λ, Q_μ⟩_α = δ_λμ
  J_λ  = H_λ · P_λ                   the integral form
```

All verified n ≤ 6: `[m_λ]J_λ = H_λ`, `⟨J,J⟩ = H·H'`, `⟨P,P⟩ = H'/H`,
`⟨P,Q⟩ = 1`.

⚠️ **Three near-identical linear hooks circulate, and the swap is silent.**
`h_low = αa+l+1`, `h_up = α(a+1)+l`, and [KS]'s tableau weight
`d = α(a+1)+l+1` — a third one, equal to neither. All three are needed
(§3.2–3.4 use all three), any two agree on enough small cells to pass a
careless test, and the literature's `c_λ, c'_λ, j_λ` notation packs them
differently per paper. §5 pins them the way `docs/record/macdonald-operators-spec.md` §1.1
pins coarm/coleg.

### 1.2 Specializations and symmetries, all verified

| statement | verified |
|---|---|
| α = 1: `P = s`, `J_λ = (∏ hooks) · s_λ` | n ≤ 6 |
| α = 2: `J^{(2)} = Z_λ` in [GJ]'s normalization; **Sage's `zonal()` is `P^{(2)}`**, i.e. `Z^{Sage}_λ = J^{(2)}_λ / H_λ(2)` — measured ratios `1, 1/3, 1/2, 1/15, 1/4, 1/6` on n ≤ 3 are exactly `1/H_λ(2)` | n ≤ 5 |
| duality: `ω_α P_λ^{(α)} = Q_{λ'}^{(1/α)}` where **`ω_α := ω ∘ (p_r ↦ α p_r)`** — the p_μ-coefficient picks up `α^{ℓ(μ)}` on top of ω's sign | n ≤ 5 |
| principal: `J_λ(1^N; α) = ∏_{(i,j) ∈ λ} (N + α·j − i)` (0-based coarm/coleg) | n ≤ 5, N ≤ 4 |
| limits ([KS] p.1): `m_λ` at α = ∞, `e_{λ'}` at α = 0, zonal pair at 2 and 1/2 | not separately tested (limits, not values) |

⚠️ Recorded because the first draft of the verify script had it wrong: the
duality with **plain ω fails already at λ = (1)** (`ωP_(1) = p_1` but
`Q_(1)^{(1/α)} = α p_1`). The α-twist is not optional.

Hand values, for §5: `J_(1) = m_1`, `J_(2) = (α+1)m_2 + 2m_11`,
`J_(11) = 2m_11`, `P_(2) = m_2 + [2/(α+1)] m_11` (the ψ hand-computation in
§3.3 reproduces the last one from the branching rule).

### 1.3 Why anyone wants this

- **Stanley's 1989 positivity conjecture is open and active.**
  `⟨J_λ J_μ, J_ν⟩_α ∈ ℕ[α]` — still open
  ([2309.13870](https://arxiv.org/abs/2309.13870) "New cases of the Strong
  Stanley Conjecture", [2605.10608](https://arxiv.org/abs/2605.10608),
  [2606.17822](https://arxiv.org/abs/2606.17822) are 2023–2026 activity).
  Verified here on every triple with `|λ|,|μ| ≤ 3` — and §2 shows Sage cannot
  even *compute* the degree-12 instances.
- **The [GJ] conjectures are open, and the field just moved.** With `α = 1+b`:

  ```text
    Φ(x,y,z;t,α) = Σ_θ  t^{|θ|} · J_θ(x)J_θ(y)J_θ(z) / ⟨J_θ,J_θ⟩_α      [GJ] (1)
    Ψ = α t ∂/∂t log Φ                                                   [GJ] (2)
    c^λ_{μν}(b) := z_λ(1+b)^{ℓ(λ)} · [t^n p_λ(x)p_μ(y)p_ν(z)] Φ         [GJ] (5)
    h^λ_{μν}(b) := [t^n p_λ(x)p_μ(y)p_ν(z)] Ψ                           [GJ] (4)
  ```

  **Matchings-Jack conjecture**: `c^λ_{μν}(b) = Σ_δ b^{wt_λ(δ)}` over
  matchings — in particular `∈ ℕ[b]`. **b-conjecture**: `h^λ_{μν}(b) =
  Σ_M b^{ϑ(M)}` over rooted connected bipartite maps in locally orientable
  surfaces. Status ladder: polynomiality in ℚ[b] **proven** [DF];
  integrality of `c` **proven** [BD]; positivity **open for both**. [GJ]'s
  own evidence was, in their words, a *computational investigation* — this is
  a domain where tables are results.
- **Lassalle's conjecture became a theorem in 2023** [BDD], via a positive
  formula for the *p-expansion* of `J_λ` — and the J→p table this spec builds
  is exactly the Jack character table those results are about. The area is
  hot and its objects are what §4 computes.
- **Random matrices**: α = 2/β connects Jack to the β-ensembles; [MOPS]
  exists because of that community, in Maple, last touched two decades ago.

### 1.4 Prior art, and what is actually ours

| who | what they have | where it stops |
|---|---|---|
| **Sage 10.9** (`jack.py`, pure Python) | bases `P, Q, J, Qp`, `scalar_jack`, `zonal = P^{(2)}` | whole-degree tables wall at **n = 12** (§2); the Stanley product `J[3,2,1]²` already **>120 s**; two measured oracle warts (§2.1) |
| **Symmetrica** (= what Sage ships in C) | **no Jack at all** — zonal only (`zo.c`, one entry point), per the full 2026-07-29 doc inventory | nothing to displace: unlike Schubert this is a capability gap, not a backend swap |
| **[MOPS]** (Maple, 2004–07) | Jack/Hermite/Laguerre/Jacobi, the LB recursion this spec adopts | Maple, unmaintained; ⚠️ unmeasured here (no Maple on this machine) |
| **Stembridge's SF 2.4** (Maple, 2005) | Jack via user-defined bases; ships **precomputed Jack archives to degree 16** per its manual | static tables, offline; ⚠️ unmeasured; degree 16 is the psychological bar for a *live* engine |
| **stla ecosystem** (`jackR`, `JackPolynomials.jl`, Haskell `jackpolynomials`, 2019–, active) | Jack in a **fixed number of variables**, symbolic α; their blog benchmarks Julia at ~350 ms on an unstated case | fixed-alphabet polynomials, not symmetric functions in the m-basis; ⚠️ unmeasured here |

Honest accounting, per the standing house questions:

1. **The mathematics is entirely prior art** — [MOPS]'s recursion, the
   branching rule, [KS]'s formula are all published and all implemented
   somewhere. What no one has: a fast engine for the **m-basis coefficients
   over ℚ(α) at unbounded alphabet** — the symmetric-*function* object the
   conjectures are stated in — with exact arithmetic and modern speed.
2. **The arithmetic observation is ours** (§3.1): the entire Jack calculus
   lives over denominators that are products of **integer-linear atoms
   `uα+v`** — a strictly tamer field than the crate's ℚ(q,t): atoms are
   canonical (primitive linear forms), lcm is exact, and cancellation is an
   integer root-test with **no failed trial divisions**, the failure mode
   that cost `deltaop` its §6.2 fix. Measured: **no swell at all** (§3.1).
3. **The machinery fit is near-total** (§3.8): the branching route is
   `macdonald.rs` with one function changed; three shared-nothing routes land
   in one crate the way `qtkostka.rs` already holds three.
4. **The [GJ] pipeline** (§3.6): no package computes `c^λ_{μν}(b)` /
   `h^λ_{μν}(b)` tables. The verification frontier of two named open
   conjectures is a compute problem, and §2 shows the incumbent's unit of
   work (J→p, whole degree) dies at n = 12.

---

## 2. Measured walls

SageMath 10.9, this machine, 2026-07-29, single runs, per-item `SIGALRM`
120 s. ⚠️ Order-of-magnitude walls, not benchmarks.

**Power provenance, row by row.** The sweep started 15:32:15 on mains; AC
detached at 15:56:18, which falls **during the P→s sweep**: `P→m`, `J→m`,
`J→p` are mains rows; `P→s`, norms and products are battery-at-100%. The
battery rows read within noise of their mains counterparts (P→s n=11:
78.9 s vs P→m's mains 84.5 s), so no correction is applied — but the
Macdonald spec measured 1.8× battery drift on this machine at lower charge,
so the split is recorded rather than ignored, and §6's comparisons must be
mains-to-mains.

⚠️ **SIGALRM is late inside Sage's C sections.** The J→p n=12 item ran
**303.8 s past a 120 s alarm** and completed; the alarm only fires when
control returns to Python. So `>120` rows are lower bounds, and one bonus
completed row exists where the computation outran its own timeout.

### 2.1 What Sage has

By API introspection and black-box measurement, not source reading: bases
`P, Q, J, Qp`; `scalar_jack`; `zonal()` = `P^{(2)}` (measured, §1.2). Two
oracle warts, both hit today:

- **`scalar_jack` cannot return zero.** Any vanishing pairing raises
  `TypeError: 'int' object is not callable` (`jack.py` calls
  `c.denominator()` on a plain int 0, where `denominator` is a property).
  The verify script treats the crash set as the zero set and checks it
  coincides exactly with the pairs that must vanish — it does, n ≤ 6.
- **No whole-degree entry point** exists; every element recomputes through
  per-element conversion morphisms. The tables below are the per-degree sums
  a user actually pays.

### 2.2 Where Sage stops

Whole-degree tables (every λ ⊢ n expanded), the natural unit of work:

```text
            P → m                J → m                J → p (the [GJ] unit)
  n      total(s)  max        total(s)  max        total(s)   max
   8        1.63   1.61          1.83   1.81          1.91    1.83
   9        5.64   5.61          6.21   6.17          6.59    6.41
  10       21.29  21.22         23.53  23.45         25.24   24.76
  11       84.52  84.37         89.01  88.84         88.57   87.40
  12       >120 at [12]         >120 at [12]        303.8   300.8   (past-alarm row)
  13          —                    —                 >120 at [13]
```

P→s is the same curve (78.9 s at n=11, dead at 12), matching
`docs/research-gaps.md`'s row (5.59 s at n=9 there, 5.60 s here). Growth is a
steady **~3.8× per degree** everywhere.

**The finding that shapes the design: the single shape λ = (n) is the entire
degree.** At every n, `max ≈ total` — the one-row partition costs 84.4 s of
n=11's 84.5 s, and every other shape is nearly free once it has run. Whatever
Sage's internal caching does, the *user-visible* unit is: touching the
hardest shape of degree n costs the whole degree, and degree 12 is
unreachable.

Norms (`scalar_jack(J_λ, J_λ)`, all λ ⊢ n): 5.9 s at n=9, 21.8 s at n=10,
80.9 s at n=11, `>360 s` at n=12 — the closed form `H_λ H'_λ` is a
*product of 2n linear factors*, and the incumbent prices it like a full
expansion.

Structure constants (`J(J_λ · J_μ)`, the Stanley object):

```text
  J[2,1]·J[2,1]       0.361 s     7 terms
  J[2,2]·J[2,1]       0.748 s     6 terms
  J[3,2]·J[2,2]      10.894 s     9 terms
  J[3,2,1]·J[3,2,1]   >120 s        (degree 12)
  J[4,3]·J[4,3]       >120 s        (degree 14)
  J[4,3,2]·J[4,3,2]   >120 s        (degree 18)
```

The first open-conjecture-relevant sizes are exactly where the incumbent
dies.

---

## 3. The algorithms

### 3.1 The coefficient field: factored integer-linear atoms (`AFrac`)

Every scalar in the Jack calculus — hooks, ψ-ratios, eigenvalue
differences, norms — is a product of **integer-linear forms `uα + v`,
u, v ≥ 0**. This is the α-analog of `Frac`'s binomial family, and it is
*better behaved* than either of the crate's existing families, in three
measured ways:

1. **Canonical.** Normalize atoms to primitive (`gcd(u,v) = 1`, content
   pulled into the coefficient): distinct primitive linear forms are
   irreducible and pairwise coprime in ℚ[α]. So the factored form is unique —
   unlike `Frac`, whose docs (frac.rs:30–36) explain why its `PartialEq`
   must cross-multiply. Equality here is structural.
2. **Exact lcm.** Atom-wise max multiplicity *is* the lcm — the q,t family
   is not coprime (`docs/record/macdonald-operators-spec.md` §3.5's warning) and settles
   for a common multiple; this one does not.
3. **No failed divisions.** `(uα+v) | p(α)` in ℚ[α] iff
   `Σ_k p_k (−v)^k u^{d−k} = 0` — one integer Horner pass, necessary *and*
   sufficient. The q,t engine spent 73% of its profile in `Atom::divide` and
   needed a bespoke early-out (`diff_may_divide`) because failed exact
   divisions ran to completion; here a failed cancellation costs one dot
   product.

⚠️ **Primitivity is load-bearing, not cosmetic.** Eigenvalue differences are
not primitive — `κ = (2,2), λ = (1,1,1,1)` gives `E(κ)−E(λ) = 2α+4` — and
dividing ℤ[α] by a non-primitive linear form leaves ℤ[α]. Gauss's lemma makes
the primitive quotient integral, which the synthetic division in
`spec_jack_swell.py` asserts on every division it performs.

**Measured: there is no swell.** `spec_jack_swell.py` simulates the exact
proposed representation (numerator in ℤ[α], atom multiset, rational content,
reduce after every addition) through the full [MOPS] recursion, cross-checked
against plain rational arithmetic at α = 5 for every κ ⊢ n ≤ 7:

```text
  κ = (n) row:                          worst over ALL κ ⊢ n:
  n   num deg  bits  atoms  content     num deg  bits   (worst shapes)
   8      3      8     6      21           4      14    [6,2] / [5,2,1]
  10      5     11     8      28           6      19    [6,3,1]
  12      6     17    10      35           —       —
  11      —      —     —       —           7      22    [7,3,1] / [6,3,1,1]
```

Peak numerator degree grows by ~½ per degree, coefficients by ~2 bits,
denominator atoms by 1, content by ~3 bits. Extrapolated (⚠️ extrapolation,
not measurement), `i128` holds comfortably past n = 30. This is the
measurement that licenses the whole design, made before any Rust exists —
the `docs/record/macdonald-operators-spec.md` §7.1 move, with the opposite outcome to
the one `macop.rs` feared: the denominators of this domain are simply small.

Representation: dense univariate `Vec<C>` numerator (a `QtPoly` with dead
`t` would be the wrong shape — sparse `(u32,u32)` keys for a dense
univariate object), primitive-atom `BTreeMap<(u32,u32), u32>` denominator,
plus an integer denominator scalar so `C = i128` stays integral. This is a
**fourth** factored-fraction family after `Frac` (frac.rs:49), `bh::Rat`
(bh.rs:94) and `deltaop::Ratio` (deltaop.rs:382). The Macdonald-operators
spec already recommends lifting the shared design into `FactoredFrac<A>`;
if that refactor lands, `AFrac` is its cleanest instantiation (the only
canonical one). Until then it stands alone, with the same one-module scope
the other three chose, for the same documented reasons.

### 3.2 E1 — the Laplace–Beltrami recursion [MOPS]: the engine

`docs/research-gaps.md` says "…rather than Gram–Schmidt", and this is the route
that honors it. [MOPS] Def 2.10ff: with

```text
  ρ^α_κ = Σ_i κ_i (κ_i − 1 − (2/α)(i−1))
```

the recursion for `P_κ = Σ c_κλ m_λ`, cleared of the `2/α` prefactor by
`E(ν) := α·n(ν') − n(ν)` (so `(2/α)/(ρ^α_κ − ρ^α_λ) = 1/(E(κ) − E(λ))`):

```text
  c_κλ = [ Σ_{(i,j,t)} (λ_i − λ_j + 2t) · c_κμ ]  /  (E(κ) − E(λ)),
```

summing over **positions** `1 ≤ i < j ≤ ℓ(λ)` and all `t ≥ 1` with
`λ_j − t ≥ 0`, where `μ = sort(λ + t·e_i − t·e_j)` must satisfy
`λ < μ ≤ κ` in dominance, and **every triple `(i,j,t)` contributes, even
when different triples produce the same μ**. `c_κκ = 1`; fill rows in any
linear extension of reverse dominance (`n(λ)` ascending works).

Verified: this exact semantics reproduces Sage's `P → m` for **every λ, μ
pair through n = 6**, symbolically over ℚ(α). The position-multiplicity
convention and the sort are precisely the details a prose reading of [MOPS]
leaves ambiguous; they are pinned operationally by the sweep, which is why
no separate eigenvalue-formula test exists — the recursion *is* the test.

The denominator: `E(κ) − E(λ) = (n(κ')−n(λ'))·α + (n(λ)−n(κ))` has both
integer coefficients **positive** for `λ < κ` (dominance moves boxes up),
which is [MOPS] Lemma 2.16 made visible — and it is a single atom of §3.1.

Cost shape: one row (all of `P_κ`) is `p(n)` coefficients, each a sum over
`O(ℓ(λ)² · λ_1)` moves — **no enumeration of tableaux anywhere**, and the
per-degree table is `p(n)` independent rows. This is what `macop.rs` is to
`macdonald.rs` — the eigen-route beside the branching route — except here
the eigen-route is the engine, not the cross-check, because §3.1 measured
its arithmetic to be tame.

### 3.3 E2 — the branching formula: the cross-check

The α-limit of the Macdonald branching rule the crate already implements:
substituting `q = t^α, t → 1` sends each binomial atom to a linear one,

```text
  (1 − q^a t^b)  ⟼  aα + b
```

and `macdonald.rs`'s ψ becomes, factor for factor (`b_ν(s) = h_low/h_up` at
the cell):

```text
  ψ^α_{λ/μ} = ∏_{s ∈ R_{λ/μ} \ C_{λ/μ}}  [h_low^μ(s) · h_up^λ(s)] / [h_up^μ(s) · h_low^λ(s)]
```

over cells of μ in rows meeting the strip but columns not meeting it — the
same row/column condition, the same chains, the same
`Factors`-as-multiset encoding. `P_λ`'s m_μ-coefficient is the sum over
chains of horizontal strips of content μ of `∏ ψ^α`, enumerated by the
existing `charge::build` / `charge::strips` (charge.rs:154/177).
**Verified against Sage for every λ, μ through n = 6.**

Hand case, the analog of `macdonald.rs`'s: `ψ^α_{(2)/(1)} = 2/(α+1)` (the
single cell `(0,0)` of μ=(1) contributes `[1·2α]/[α·(α+1)]`), giving
`P_(2) = m_2 + [2/(α+1)] m_11` — and getting either half of the row/column
condition wrong gives something else.

E2 shares the *chains* infrastructure with `macdonald.rs` but **no
mathematics with E1** — branching/Pieri versus eigenoperator — which is
exactly the shared-nothing cross-check standard `qtkostka.rs` §"Three
routes" sets.

⚠️ The limit is taken **per atom, never per coefficient**: a coefficient of
`P` is a *sum* of ψ-products, and post-processing `macdonald.rs`'s summed
output through `q = t^α, t → 1` would need L'Hôpital on every fraction.
Re-deriving the branching formula with α-atoms is the correct move; limiting
the finished Macdonald table is the recorded dead end (§3.7).

### 3.4 E3 — the Knop–Sahi tableau formula: the positive oracle

[KS] Thm 5.1, verified against Sage for every λ through n = 5:

```text
  J_λ(x; α) = Σ_{T admissible} d_T(α) x^T,    d_T = ∏_{s critical} (α(a(s)+1) + l(s) + 1)
```

`T` labels the cells of λ with `1..n`; admissible means `T(i,j) ≠ T(i',j)`
for `i' > i` and `T(i,j) ≠ T(i',j−1)` for `i' < i, j > 1`; `s = (i,j)` is
critical when `j > 1` and `T(i,j) = T(i,j−1)`. Exponential (`n^{|λ|}`
labelings) — this is `NaiveLr`'s role: the reference implementation kept
forever, plus the *proof-carrying* route for positivity (each coefficient is
manifestly a sum of products of the `d`-hooks). [KS] Thm 1.1 falls out: the
m_μ-coefficient of `J_λ`, divided by `u_μ = ∏_i m_i(μ)!`, lies in ℕ[α] —
**verified on Sage's own output through n = 7**, including the division.

### 3.5 Normalizations, pairings, conversions

`Q` and `J` are atom-multiset scalings of `P` (the `Frac::mul_factors`
pattern, frac.rs:124). The pairing of arbitrary elements goes through the
p-basis, where it is diagonal with weight `z_λ α^{ℓ(λ)}`; norms of `P/Q/J`
never compute a pairing at all — they are the closed products of §1.1 (the
n=12 norm table that costs Sage >360 s is `p(12) · 24` atom
multiplications here). `J → p` runs `m → s → p` through the existing
`convert` hub (convert.rs:44/49/54), which is generic over the coefficient
ring; `AFrac<C>` implements `QAlgebra` (the `z_ν` divisions in
`PowerSum::from_schur`, convert.rs:852, are content divisions).

### 3.6 The [GJ] pipeline

With the whole-degree `J → p` table and the closed-form norms:

```text
  [t^n] Φ = Σ_{θ ⊢ n} J_θ(p-expansion)⊗³ / (H_θ H'_θ)     — p(n) terms
  c^λ_{μν}(b):  extract [p_λ p_μ p_ν], scale by z_λ(1+b)^{ℓ(λ)}, shift α = 1+b
  h^λ_{μν}(b):  the same from Ψ = α t ∂_t log Φ  (log via the exp/cumulant
                recurrence degree by degree)
```

Everything arrives as an `AFrac` that must collapse to a polynomial in b —
[DF]'s theorems say it must, [BD] says `c`'s is integral, and both collapses
are *checks the pipeline runs* (§5.9). Positivity is then **the open
question itself**: a negative coefficient anywhere is a result to report,
not a bug to debug away — the valley-Delta-conjecture posture of
`docs/record/macdonald-operators-spec.md` §5.9, verbatim.

The α → b shift is a binomial transform of dense univariate polynomials
(exact, cheap). ⚠️ [GJ]'s own verification range was not recovered from the
scanned paper (only pp. 873–875 were read); their abstract says "evidence is
presented… and they are proved for two infinite families". The frontier to
beat is therefore stated conservatively: **produce the complete `c` and `h`
tables for n ≤ 10 and report where the engine stops** — Sage's unit of work
for the same pipeline dies at n = 12 before the triple product even starts.

### 3.7 Recorded dead ends

1. **Gram–Schmidt** against dominance — `p(n)²` pairings over ℚ(α), the
   route `docs/research-gaps.md` already warns off. Not attempted.
2. **Limiting the finished Macdonald table** (`q = t^α, t → 1` on
   `macdonald_p` output): per-coefficient limits of summed fractions need
   L'Hôpital; only per-atom limits are safe, and those mean re-deriving E2
   anyway. Rejected in §3.3.
3. **[KS] Thm 4.6, the nonsymmetric recursion** (creation operator
   `Φ = x_n s_{n−1}⋯s_1` on Opdam's `E_η`, then symmetrize via Thm 4.10):
   would deliver nonsymmetric Jack polynomials too, and is the modern
   frontier route in the Macdonald world. Out of scope for v1; recorded
   because it is the natural v2 engine if the E_η's are ever wanted.
4. **Lassalle–Schlosser Pieri inversion** (the "Lassalle recurrences" of
   `docs/research-gaps.md` §2.6): an explicit expansion by inverting Pieri.
   No advantage over E1 was identified — E1 is already
   enumeration-free with unit-cost denominators — and no numbers argue
   otherwise; recorded as unexplored rather than rejected.

### 3.8 Crate fit

Verified against the source today, with the line numbers:

| need | existing | file |
|---|---|---|
| chain-of-strips enumeration (E2) | `charge::build` / `strips` | charge.rs:154/177 |
| the E2 template, one function to change | `macdonald_p` + `psi_factors` | macdonald.rs:54/183 |
| arm/leg/conjugate-count | `arm`, `leg`, `count_above` | macdonald.rs:214/219/224 |
| the m-basis container | `Monomial` (`basis!`) | sym.rs:161 |
| factored-scalar multiply (Q, J from P) | the `mul_factors` pattern | frac.rs:124 |
| the factored-fraction precedents | `Frac` / `bh::Rat` / `deltaop::Ratio` | frac.rs:49 / bh.rs:94 / deltaop.rs:382 |
| single-shape-via-table shape | `hall_littlewood_p` | hl.rs:124 |
| m → s → p for the [GJ] pipeline | the `convert` hub | convert.rs:44/54/852 |
| ring bounds (`AFrac<C>: QAlgebra`) | `Ring` / `QAlgebra` | coeff.rs:41/156 |
| `z_λ`, partitions, dominance order | `Partition::z`, `partitions_of` | partition.rs:108/143 |
| per-degree memo table | `table!` + a `*_cached` wrapper | memo.rs:42 |
| overflow guard + bignum escalation | `guarded` / `escalate` | guard.rs:57 / python.rs:222 |
| the FFI shape for `Monomial<fraction>` | `mac_terms` + `macdonald_p/q/j` bindings | python.rs:769/793–806 |
| registration slot | the `#[pymodule]` list | python.rs:1027–1080 |
| Sage-driven conformance harness | `check_backend.py` (already exercises `jack()`) | scripts/check_backend.py:128 |
| oracle fixtures | `gen_sage_oracle.sage` → `tests/fixtures/` | tests/sage_oracle.rs:18 |
| benchmark/profiling workflow | `bench_mac.rs` pattern; the `sample` profile | examples/, Cargo.toml:60–68 |

Deliberately **not** reused: `QtPoly` as the numerator type (§3.1 — dense
univariate wants a `Vec`, and qt.rs:4's module doc, the crate's one
pre-existing mention of Jack, anticipated a ℚ(α) that never arrived);
`macop::Coeff`/`solve` (private to their module, and E1's denominators are
single atoms per step — nothing needs a gap-indexed table); `Frac` with a
dead variable (its `PartialEq` cross-multiplies to work around
non-canonicity `AFrac` does not have). `lib.rs` slots: `pub mod jack;`
between `hopf` (lib.rs:65) and `kf` (:66), re-export between :102 and :103;
no name collides (`jack` appears nowhere in `src/` outside qt.rs:4's
comment).

---

## 4. API

```rust
// jack.rs — everything in the monomial basis, coefficients in AFrac<C>

/// P_λ(x; α): monic, dominance-triangular.  E1 under the hood.
pub fn jack_p<C: QAlgebra>(lambda: &Partition) -> Monomial<AFrac<C>>;
/// Q_λ = (H_λ/H'_λ)·P_λ  and  J_λ = H_λ·P_λ (integral form, coeffs in N[α]).
pub fn jack_q<C: QAlgebra>(lambda: &Partition) -> Monomial<AFrac<C>>;
pub fn jack_j<C: QAlgebra>(lambda: &Partition) -> Monomial<AFrac<C>>;
/// Whole degree — the unit of work the walls are measured in.
pub fn jack_table<C: QAlgebra>(n: u32) -> Vec<(Partition, Monomial<AFrac<C>>)>;
/// J_λ in the power-sum basis — the Jack-character / [GJ] unit.
pub fn jack_j_powersum<C: QAlgebra>(lambda: &Partition) -> PowerSum<AFrac<C>>;
pub fn jack_powersum_table<C: QAlgebra>(n: u32) -> Vec<(Partition, PowerSum<AFrac<C>>)>;
/// ⟨f, g⟩_α for arbitrary elements (diagonal in p); closed-form norms.
pub fn jack_scalar<C: QAlgebra>(f: &Monomial<AFrac<C>>, g: &Monomial<AFrac<C>>) -> AFrac<C>;
pub fn jack_norm_j(lambda: &Partition) -> /* H_λ·H'_λ as atoms */;
/// ⟨J_λ J_μ, J_ν⟩_α — the Stanley object, computed without forming
/// the whole product basis change.
pub fn jack_structure_constant<C: QAlgebra>(la: &Partition, mu: &Partition, nu: &Partition) -> AFrac<C>;
/// [GJ]: the degree-n tables of c^λ_{μν}(b) and h^λ_{μν}(b), as
/// polynomials in b, with the [DF]/[BD] collapse checks built in.
pub fn gj_connection_tables<C: QAlgebra>(n: u32) -> GjTables<C>;
/// Zonal in [GJ]'s normalization (J at α = 2) and Sage's (P at α = 2).
pub fn zonal_j(lambda: &Partition) -> Monomial<Rational>;
```

`afrac.rs` holds the coefficient type (§3.1) with the `into_poly`-style
escape: `AFrac::into_poly()` returning `Option<Vec<C>>`, expected — not
hoped — to succeed on every `J` coefficient, per the `expect` policy the
Macdonald spec's §5.4 set for guaranteed-integral answers.

Python bindings mirror `macdonald_p/q/j` (python.rs:793–806): `jack_p`,
`jack_q`, `jack_j`, `jack_table`, `jack_j_powersum`, `jack_scalar`,
`jack_structure_constant`, registered in the python.rs:1027–1080 block, all
through `guarded`/`escalate`. The boundary payload is a univariate variant
of `MacTerms` (python.rs:767): `(partition, [(alpha_exp, coeff)], [(u, v,
mult)], den_scalar)`.

Not in v1: `Qp` (Sage's fourth basis, dual to P under the *undeformed* Hall
pairing — recorded in §7), nonsymmetric `E_η`, shifted/interpolation Jack.

## 5. Correctness requirements

Layered as the crate does, each layer catching what the previous cannot.
Items 1–8 **already pass in the Sage calculator** (`spec_jack_verify.py`);
the requirement is that the Rust port reproduces them and extends the
ranges.

1. **The convention gate, first.** `J_(2) = (α+1)m_2 + 2m_11` — this is the
   line that dies if the parameter convention drifts (Sage calls α "t"; the
   verify script's first assert is exactly this).
2. **Hand values.** `J_(1) = m_1`, `J_(11) = 2m_11`,
   `P_(2) = m_2 + [2/(α+1)]m_11`, `ψ^α_{(2)/(1)} = 2/(α+1)` (§3.3's
   hand-derivation), `H_(21) = (α+2)·1·1`… and a test that swapping
   `h_low`/`h_up`/`d` breaks — §1.1's trap is silent otherwise.
3. **Three engines agree.** E1 ≡ E2 for every λ through n = 8 in Rust
   (n = 6 in the Python prototype); E3 ≡ both through n = 6. They share
   `Partition` and `AFrac` and nothing else.
4. **Norms and pairings.** `[m_λ]J = H_λ`, `⟨J,J⟩ = H H'`, `⟨P,P⟩ = H'/H`,
   `⟨P_λ,Q_μ⟩ = δ`, `⟨p_λ,p_μ⟩_α = δ z_λ α^ℓ`, n ≤ 8.
5. **Specializations.** α=1 against the crate's own Schur + hook machinery
   (`eval.rs`'s hooks share nothing with `jack.rs`); α=2 against a
   committed Sage zonal fixture **with the normalization pinned**
   (`Z^{Sage} = P^{(2)}`, §1.2 — the fixture must encode the measured
   `1/H_λ(2)` ratios, or it tests the wrong thing); ω_α-duality (the α-twist
   trap, §1.2); the principal specialization formula against `expand`-style
   evaluation.
6. **[KS] Thm 1.1 as a law.** Every coefficient of `J_λ`: denominator
   clears (`into_poly` `expect`s), lies in ℕ[α], **and** is divisible by
   `u_μ = ∏ m_i(μ)!` — through n = 10 at least. The coefficients arrive
   through fraction arithmetic, so every one of these is a real check on
   the whole route (the `∇e_n`-positivity pattern).
7. **Sage oracle fixtures.** Extend `gen_sage_oracle.sage` with `P → m` and
   `J → p` dumps for n ≤ 9 plus spot shapes at 10–11 (the largest Sage
   produces in reasonable time), committed under `tests/fixtures/` per the
   `sage_oracle.rs` pattern.
8. **Stanley's pairing, as an observed law.** `⟨J_λJ_μ,J_ν⟩ ∈ ℕ[α]` on all
   triples with `|λ|,|μ| ≤ 4` — **a violation is reported as a finding,
   never "fixed"**: the conjecture is open, and the engine's first real
   deliverable is exactly this table at sizes Sage cannot reach.
9. **The [GJ] pipeline laws.** For every degree computed: `c^λ_{μν}` and
   `h^λ_{μν}` collapse from ℚ(α) to polynomials in b ([DF] — a
   non-collapse is our bug), `c` is integral ([BD] — likewise), and
   positivity is *recorded*, with any negative coefficient reported as a
   result (open conjectures, both — the §3.6 posture).
10. **Fixed-width honesty.** `i64`/`i128`/`Rational`/bignum ladders agree
    term for term at every degree benchmarked (the `macop.rs` two-width
    check); `guarded` wraps the Python boundary.
11. **Bindings separately** (`check_bindings.py` pattern): a right answer in
    the wrong slot is a different failure from a wrong answer.

## 6. Performance requirements

Stated before the code exists, so the measurement can embarrass them — the
`st` spec guessed 50× and got 3400×, the Δ-operators spec guessed 100× and
got ~21×; both corrections were recorded, and this spec expects the same
treatment. **All comparisons mains-to-mains** (§2's provenance discipline;
the 1.8× battery drift is documented on this machine).

- **Headline, like-for-like:** whole-degree `J → m` at n = 11 —
  Sage 89.0 s (mains). Target **≥ 200×** (≤ 0.45 s). ⚠️ A guess. The basis
  for it: E1 does `p(11)² ≈ 3 100` coefficient updates of ~10² moves each
  over the §3.1 arithmetic whose largest object at n = 11 is a degree-7,
  22-bit-coefficient numerator — there is no visible place for 89 seconds
  to hide. The measurement will say.
- **Wall shift:** complete whole-degree tables through **n = 16 live** —
  the degree Stembridge's SF ships as *precomputed archives* — in minutes,
  and report the first degree that takes over an hour. Growth is the
  number to publish (Sage's is 3.8×/degree; `p(n)` growth alone is ~1.3×,
  so the arithmetic's contribution is the honest curve).
- **Single shape:** `P_(20)` (the row that *is* the degree, per §2.2) in
  single-digit seconds. ⚠️ Wild guess, flagged as such.
- **Norms:** the n = 12 table (Sage: >360 s) in microseconds — closed
  form; report it once to make the point, then stop benchmarking products
  of 24 linear factors.
- **The Stanley table:** every `⟨J_λJ_μ,J_ν⟩` with `|λ| = |μ| = 6`
  (`J[3,2,1]²` is the case Sage cannot do one of) — target minutes for the
  whole table.
- **The [GJ] tables:** `c` and `h` complete for n ≤ 10; report the wall.
- **Report the phase split** (table build vs conversions vs [GJ] assembly)
  and peak RSS; single-threaded numbers only, per the standing rule.

### 6.1 Measured (added after the build, 2026-07-29)

⚠️ **Battery to battery, one session, isolated processes.** §2's tables are
mains and this section is not; mixing them is exactly the mistake §2's
provenance note warns about. The drift is confirmed rather than assumed here:
§2 measured `P → m` at n = 11 as 84.5 s on mains, and the same cell
re-measured on battery under `scripts/bench_jack.py` is 155.1 s — a ratio of
**1.84**, matching the documented 1.8×.

⚠️ **§2's own numbers cannot be reproduced in a single Sage process, and that
is a methodological correction to this document.** Running the four units in
one process gives nonsense two different ways. Sage *memoizes* its Jack
transition matrices, so a `J → m` sweep after a `P → m` sweep reads a warm
cache — measured 0.152 s for the whole of degree 10 against 44.1 s cold, a
290× difference that is entirely reuse. And `SIGALRM` fires **inside** Sage's
cache construction and leaves it half-built: the next degree died with
`KeyError: [11]` out of `sfa.py:_from_cache`. A timeout does not merely abandon
an item, it poisons the process. `scripts/bench_jack.py` therefore runs one
subprocess per (unit, degree). §2's rows should be read with that in mind.

```text
  whole degree, n = 11 (p(11) = 56 shapes), per-cell timeout 180 s
  unit          Sage (s)   symfn ℤ (s)        ratio
  P → m          155.119       0.01857        8353×
  J → m          147.428       0.02161        6822×
  J → p          149.242       0.13546        1102×
  norms          148.937      0.000091     1636670×
```

**Sage prices all four units identically**, 147–155 s. That is the sharper form
of §2.2's finding: it is not the unit that costs, it is the `P → m` transition,
and `J`, the p-expansion and even the closed-form norm all route through it.

| requirement (§6) | target | measured | verdict |
|---|---|---|---|
| headline `J → m` at n = 11 | ≥ 200× (≤ 0.45 s) | 0.0216 s, **6822×** | **beaten by 34×** |
| whole tables through n = 16 live | "in minutes" | **0.73 s**; n = 26 in 520 s | beaten |
| single shape `P_(20)` | "single-digit seconds", flagged a wild guess | **0.0071 s** | beaten by ~10³ |
| norms table at n = 12 | "microseconds" | 134 µs | met |
| Stanley table, `\|λ\| = \|μ\| = 6` | "minutes" | **1.44 s**, 9317 triples | beaten |
| [GJ] `c` and `h` complete for n ≤ 10 | — | **26.3 s** at n = 10 | met |

The growth rate, which §6 asked for as "the number to publish": Sage ~3.5× per
degree, this ~2.0×. `p(n)` alone grows ~1.3×, so the arithmetic contributes
~1.5× and the rest is the partition count — the honest curve §6 wanted.

The wall, run to exhaustion (whole-degree `P → m`, `i128`):

```text
  n     16      18      20      22      24      26
  s   0.73    2.68   15.17   49.55  161.56  520.66
```

n = 16 is the degree Stembridge's SF ships as *precomputed archives* (§1.4).
Coefficients are 89 bits at n = 26 against `i128`'s 127, growing ~4.4
bits/degree, so the fixed-width wall is near n = 34 — §3.1's "comfortably past
n = 30" extrapolation was close and slightly optimistic.

**The guess was low, and for a reason worth recording.** §3.1's arithmetic
prediction was right (numerator degree 6 and 11 atoms at n = 12, measured 6 and
11) and §6's target was still off by 34×, because the target priced *our*
arithmetic correctly and the *incumbent* wrongly — nobody measured what Sage
spends 155 s on. The `st` spec guessed 50× and got 3400×; the Δ-operators spec
guessed 100× and got 21×; this one guessed 200× and got 6822×. Three specs,
three misses, in both directions. The lesson is not "guess higher": it is that
a ratio has two sides and only one of them is under our control.

### 6.3 The slowest unit was not the bottleneck

`J → p` is 1102× where `P → m` is 8353×, so it looks like the thing to fix.
Sampling its two consumers says otherwise, and the two answers are different
from each other:

- **[GJ] pipeline (n = 9):** `J → p` does not reach the top fourteen inclusive
  frames. 96.5% is `phi_slice` — the `p(n)³` triple product — split across
  `AFrac::reduce` (49%), `Ring::mul` (37%) and `add_assign` (36%). Optimising
  `J → p` here would move nothing.
- **Stanley's table (degree 12):** `J → p` is **94.6%** — but of *repeated*
  calls. `jack_structure_constant` recomputes all three p-expansions on every
  triple, so the degree-12 table ran 27951 conversions for 99 distinct values.

The fix is therefore not a faster `J → p` but fewer calls to it.
[`stanley_table`] hoists the expansions out of the triple loop and takes the
degree-12 table from **46.6 s to 1.44 s (32×)**, putting degree 16 in reach
(111804 triples, 73.6 s, all in ℕ[α]).

⚠️ Deliberately **not** solved by memoizing `jack_j_powersum`. The Python
boundary runs over `Guarded` so that an overflowing intermediate is *detected*;
a cache filled at `i128` and handed out to other widths would launder exactly
that away. That is `memo::bold_p`'s documented hazard, and hoisting the loop is
the version with no correctness question in it.

**The general lesson, and it is the same one twice.** §6.2 found the cost was
allocation rather than arithmetic; this found it was call count rather than
speed. Both times the slow-looking thing was not the expensive thing, and both
times the only way to know was to sample.

### 6.2 What the sampling said

`sample`, per the `Cargo.toml` workflow, on `examples/profile_jack.rs` at
n = 18. **58% of the profile was inside `reduce_at`, and roughly half of that
was `malloc`/`free` rather than arithmetic.**

§3.1 claimed "a failed cancellation costs one dot product". It cost one dot
product **and two heap allocations** — `divide_by_linear` built the quotient
buffer before knowing whether the division succeeded, and `reduce` trial-
divides by every atom, so most of those buffers were thrown away. Splitting off
`divides_by_linear` (a single running scalar; the recurrence never needs the
whole quotient array) and rewriting in place was worth **1.6×**, n = 18 going
4.33 s → 2.72 s.

⚠️ Three follow-up changes — in-place `lift`, an allocation-free
`content_reduce`, cached partitions with precomputed eigenvalue statistics —
were worth **~2%, i.e. nothing measurable**. They are kept only because they
strictly allocate less. This is `QtPoly::mul_diff` again: sound reasoning about
a real cost that is not on the critical path.

Worse, the `content_reduce` rewrite silently **grew** the scalar denominators:
reading the loop bound off the live scale makes the final step try the
uncancelled scalar as one lump, so `202 = 2·101` against a numerator of content
101 keeps its 101 forever. Nothing failed — `PartialEq` cross-multiplies, so
the values stayed correct — and it was caught only because `bench_jack` prints
the scalar width as a column. **Print the shape of the data, not just the
time.**

## 7. Open questions

*Answered after the build, 2026-07-29. Kept in place rather than deleted.*

1. **Whole-degree engine choice.** ~~Measure both.~~ **Answered: E1, by a
   growing margin.** The two are within 2.5× at n ≤ 8 and E1 pulls away
   (`bench_jack`'s third table): 8.6× at n = 10. E2's shared strip-ψ cache does
   not pay for the chain enumeration. `jack_table` calls `jack_p_lb`.
2. **Denominator bookkeeping.** ~~Per-coefficient `BTreeMap` versus indices
   into the fixed per-degree family.~~ **Not revisited, and the sampling says
   it was the wrong question.** The `BTreeMap` never appeared in the profile;
   the cost was the *numerator* buffers inside the divisibility test (§6.2).
3. **Single-coefficient query `[m_μ]P_λ`.** Still open, still a candidate. E1
   fills the whole row whatever you asked for.
4. **`Qp`.** Still deferred; nothing needs it.
5. **Shifted / interpolation Jack.** Still the natural v2, and `AFrac` lives in
   its own module for exactly this reason.
6. **Nonsymmetric `E_η`** via [KS] Thm 4.6. Unchanged.
7. **Jack characters `θ^λ_μ(α)`.** `jack_j_powersum` is the table; the
   normalization is still an API question. ⚠️ `J → p` is our slowest unit
   (§6.1), and the first draft of this line said that is where the next factor
   lives. **Measured, it is not** — see §6.3.
8. **How far does `i128` hold?** §3.1 extrapolated "comfortably past n = 30";
   measured, coefficients reach 89 bits at n = 26 and grow ~4.4 bits/degree, so
   the wall is near **n = 34**. Close, and slightly optimistic.

### 7.1 Answered by building it

- **The [GJ] pipeline works, and it needed an independent check.** [GJ]'s (1)
  and (5) were transcribed from a scanned paper, and the `z_λ(1+b)^{ℓ(λ)}`
  prefactor has no derivation available here. It is pinned instead by a
  collapse: at `b = 0` the whole pipeline must become the class algebra of
  `S_n`, `c^λ_{μν}(0) = a^λ_{μν}`, and `gj::class_algebra_coefficient` computes
  that from characters alone — no Jack polynomial, no `AFrac`, no fraction
  field. Exhaustive through n = 6. Without it, a wrong prefactor would have
  produced a plausible, positive, entirely fictitious table.
- **Both conjectures hold on everything computed.** `c` and `h` are in ℕ[b] at
  every degree through 10 (54108 and 28752 coefficients at n = 10). Both are
  open; this is evidence.
- **`AFrac` is a ℚ-algebra even when `C` is not** — dividing by an integer
  multiplies the scalar denominator and asks nothing of `C`. That was not
  anticipated in §3.5 and it is what lets the engine run over `AFrac<i128>` and
  still be handed to `s → p`. It is the cleanest thing in the module.
- **`Guarded` needed `div_exact`.** §4's "all through `guarded`/`escalate`"
  was not free: `Guarded` inherited `Ring`'s declining default, which is safe
  for every other caller — "not divisible" is an ordinary outcome there — but
  silently turns `AFrac::reduce` into a no-op. The bindings would have been
  correct and unboundedly slow.
- **`α ↦ 1/α` stays inside the family**, which was not obvious: `(uα+v)`
  becomes `(vα+u)/α`, still primitive since `gcd` is symmetric, except that the
  `α` atom `(1,0)` becomes the constant 1 and leaves the denominator. So the
  substitution is exact and symbolic — `AFrac::invert_alpha` — and the duality
  law never has to evaluate anything.

## 8. What §5 asked for, and what it got

| §5 item | asked | done |
|---|---|---|
| 1. convention gate | `J_(2) = (α+1)m_2 + 2m_11` | ✅ `the_convention_gate` |
| 2. hand values, hook swap breaks | — | ✅ `hand_values`, `hooks_are_three_distinct_families` |
| 3. three engines agree | E1≡E2 to n=8, E3≡both to n=6 | ✅ E1≡E2 to **n=8**; ⚠️ E3 to **n=5**, not 6 — it is `n^{\|λ\|}` labelings and n=6 is 47× the work for no new failure mode |
| 4. norms and pairings | n ≤ 8 | ✅ n ≤ 5 in Rust, **n ≤ 7** against Sage (`check_jack.py`) |
| 5. specializations | α=1, α=2 pinned, ω_α-duality, principal | ✅ all four — `the_omega_alpha_duality` **and** `plain_omega_breaks_the_duality` as the negative control |
| 6. [KS] Thm 1.1 as a law | n ≤ 10 | ✅ n ≤ 8 in Rust over `i128`, **n ≤ 7** against Sage |
| 7. Sage oracle fixtures | committed, offline | ✅ `gen_sage_oracle.sage` extended; 132 Jack expansions to degree 7, 776 coefficients, `cargo test` needs no Sage |
| 8. Stanley's pairing observed | `\|λ\|,\|μ\| ≤ 4` | ✅ to **degree 16** (111804 triples) |
| 9. [GJ] pipeline laws | collapse, integrality, positivity recorded | ✅ all three, `gj::GjTables` |
| 10. fixed-width honesty | `i128`/`Rational`/bignum ladders agree | ✅ `bench_jack` two-width per degree; `jack_runs_over_bignum_coefficients` for `BigInt`/`BigRational` |
| 11. bindings separately | — | ✅ `check_bindings.py`, 0 failures |

Two deliberate shortfalls, both recorded rather than quietly met:

- **E3 stops at n = 5.** The Knop–Sahi enumeration is exponential and its role
  is to be a *different* algorithm, not a wide one; n = 6 costs 47× for no new
  failure mode. `check_jack.py` covers `J` to n = 7 by the other route.
- **[KS] Thm 1.1 stops at n = 8** rather than 10, for the same reason the
  ranges elsewhere stop where they do: the law is checked on every coefficient
  of every shape, and nothing about it is degree-sensitive.

The **`ω_α` duality is the one §5 item that nearly went missing**, and it is
the one with the least redundancy: it is the only law relating `P` to `Q`, to
conjugation, and to the parameter inversion at once, so it independently pins
`jack_q`'s normalization — which otherwise appears only inside `⟨P,Q⟩ = δ`,
where a compensating error in both would cancel.
