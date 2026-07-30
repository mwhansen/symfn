# Specification: the Orellana–Zabrocki character bases (`st`, `ht`)

*Written 2026-07-29. Implements §2.2 of `docs/research-gaps.md` — "reduced / stable
Kronecker as a first-class ring".*

The source is Orellana–Zabrocki, *Symmetric group characters as symmetric
functions*, [arXiv:1605.06672v5](https://arxiv.org/abs/1605.06672) (published in
*Adv. Math.*). Equation and theorem numbers below are that paper's, from the v5
PDF read on 2026-07-29 — not recalled.

**Sage is an oracle here, never a source.** Every number in §2 was measured on
this machine against SageMath 10.9 (release 2026-05-04) today; every formula in
§3 was checked against Sage's answers before being written down. Sage's
implementation of these bases was not read, and must not be. The same clean-room
posture as `docs/cleanroom-spec-skew-lr.md`, for the same reason.

---

## 1. What the basis is

`Sym` has the six classical bases, all homogeneous. The two bases here are
**inhomogeneous**: `s̃_λ` is a sum of terms of every degree from 0 to |λ|.

**The irreducible character basis `s̃_λ` (`st`).** Theorem 1(1): for a fixed
partition λ, `s̃_λ` is the unique symmetric function such that for all n ≥ |λ|+λ₁
and all γ ⊢ n,

```text
    s̃_λ(x₁,…,x_n) = χ^{(n−|λ|,λ)}(γ)
```

where the xᵢ are the *eigenvalues of a permutation matrix* of cycle type γ. The
paper writes this evaluation `f[Ξ_γ]`. So Schur functions are the characters of
GLₙ and these are the characters of Sₙ, in the same ring — that is the whole
idea, and the reason λ indexes the shape `(n−|λ|,λ)` with the long first row left
implicit.

**The induced trivial character basis `h̃_λ` (`ht`).** Definition 4 / Eq (6):
`h̃_λ[Ξ_μ] = ⟨h_{|μ|−|λ|} h_λ, p_μ⟩`, the character of the trivial representation
induced from a Young subgroup — the permutation module `M^{(n−|λ|,λ)}`. It is to
`s̃` what `h` is to `s`, and it exists because the transition between the two is
the tool that makes everything else computable.

**Why anyone wants this (Theorem 7, Eq 16).** The product is

```text
    s̃_λ · s̃_μ = Σ_{|ν| ≤ |λ|+|μ|} ḡ^ν_{λμ} · s̃_ν
```

where `ḡ` are the **reduced (stable) Kronecker coefficients** — the eventual
value of `g((n−|λ|,λ), (n−|μ|,μ), (n−|ν|,ν))` once n is large. This is an
ordinary product of symmetric functions that happens to compute a Kronecker-type
invariant, which is why it sidesteps the p(n)×p(n) character-table blowup
`docs/record/kronecker.md` records (1.1 GB at n=32): nothing here is ever indexed by a
partition of n.

Three facts worth keeping as tests, all from the paper:

- Theorem 1(3): `s_{1^r} = s̃_{1^r} + s̃_{1^{r−1}}`, equivalently
  `s̃_{1^r} = Σ_{i=0}^{r} (−1)^i e_{r−i}`.
- Theorem 1(2): `s_λ = Σ_μ r_{λμ} s̃_μ`, where `r_{λμ}` is the multiplicity of
  `S^{(n−|μ|,μ)}` in the restriction of the GLₙ-module `W^λ` to Sₙ. So the
  transition matrix *is* the restriction problem, and its entries are
  non-negative.
- Eq (20): `h_{21} = s̃_3 + s̃_{21} + 4s̃_2 + 3s̃_{11} + 7s̃_1 + 4s̃_∅`.

### 1.1 Prior art, and what is actually ours

| who | what they have | where it stops |
|---|---|---|
| **Sage** | `st`, `ht`, `o`, `sp` bases; transitions; product | The product. Measured below: `st[4,3]²` takes 30s, `st[5,3]²` exceeds 90s. |
| **Stembridge `SF`** | user-defined bases, so this is expressible | No reduced-Kronecker support out of the box; Maple-speed. |
| **`barvikron`** | Christandl–Doran–Walter lattice points, polynomial time for *bounded height* | Kronecker, not reduced Kronecker; Python prototype, unmaintained. |
| **Baldoni–Vergne–Walter** | vector partition functions, bounded length | Maple, distributed as research code. |
| **lrcalc** | LR only | Does not know what a Kronecker coefficient is. |

Honest accounting: **the mathematics is entirely prior art and mostly the source
paper's.** Eq (23) is OZ's recommended computational route and we are taking it.
The stable Kostka transition is their Eq (7)–(8). Nothing in §3.1–§3.2 is a new
idea.

What is ours is narrower and worth stating exactly:

1. **The product is the wall, and the product is what we are good at.** Sage's
   transitions are already fast (§2) — `st[6,4] → s` costs 0.07s. It is
   `st_λ · st_μ` that dies. Our product is a Schur-basis product, and the Schur
   product is the single most optimized thing in this crate: three LR backends,
   `AutoLr` selected per shape, `memo::product_cached` memoizing whole
   expansions. The transition-then-multiply route reduces the open problem to
   ~10³–10⁴ *small* LR products, all cached, all integer.
2. **A second product route that shares no code with the first** (§3.4), derived
   here rather than read: the `ht`-basis product is a sum over non-negative
   integer matrices. This is the crate's established pattern — three LR
   backends, three (q,t)-Kostka algorithms — and it is what makes a claim about
   reduced Kronecker coefficients believable when no third-party package can
   check the sizes we are targeting.
3. **A single-coefficient query.** `docs/research-gaps.md` §2.1: no package computes
   one Kronecker-type coefficient without computing the whole product. Lemma 20
   of the paper is an inner-product formula for exactly that (§3.6). ⚠️ It is
   *not* obviously faster than the whole product; it is a candidate, not a claim.

---

## 2. Measured walls

SageMath 10.9, this machine, 2026-07-29. Single runs, per-item `SIGALRM` of 90s.
⚠️ Order-of-magnitude walls, not benchmarks — the same caution `docs/record/README.md`
applies to its own tables.

**Transitions are cheap.** This is the finding that redirected the design; an
earlier sketch of this document assumed the transition was the bottleneck and
planned to attack it.

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

**The product is not.**

```text
  product                       terms       sec
  st[2,1] · st[2,1]                26    0.0221
  st[3,1] · st[3,1]                53    0.0577
  st[3,2] · st[3,2]               101    1.4418
  st[4,2] · st[4,2]               186    8.6438
  st[4,3] · st[4,3]               308   30.2569
  st[5,3] · st[5,3]                 —      >90
  st[6,4] · st[6,4]                 —      >90
  st[8,5] · st[7,4]                 —      >90
  st[2,1,1] · st[2,1,1]            54    0.0108
  st[3,2,1] · st[3,2,1]           211    0.1467
  ht[6,4] · ht[6,4]               120   18.3357
```

Two things to notice. The wall is at **degree 10–12 in |λ|+|μ|**, which is small.
And the cost is *not* a function of the output size: `st[3,2,1]²` returns 211
terms in 0.15s while `st[4,2]²` returns 186 terms in 8.6s — 58× slower for fewer
terms, at the same total degree. Whatever Sage's product does, its cost tracks
something other than the answer. We do not know what, and per the clean-room rule
we will not look; we only need to not reproduce it.

---

## 3. The algorithms

### 3.1 `st → p → s`, by Eq (23)

Theorem 14: for λ ⊢ n,

```text
    s̃_λ = Σ_{γ ⊢ n} χ^λ(γ) · 𝐩_γ / z_γ                                    (23)

    𝐩_{i^r} = Σ_{k=0}^{r} (−1)^{r−k} i^k C(r,k) · ( P_i )_k                (24)
    P_i     = (1/i) Σ_{d|i} μ(i/d) p_d          (μ = number-theoretic Möbius)
    ( x )_k = x(x−1)···(x−k+1)                  (falling factorial, in Sym)
    𝐩_γ     = Π_i 𝐩_{i^{m_i(γ)}}
```

The paper recommends this one explicitly ("The reader interested in computing
`s̃_λ` from the formulae in this paper is recommended to use Equation (23)").

**Verified**: implemented against Sage as a calculator and compared to Sage's own
`st`, exact agreement on all 12 of ∅, (1), (2), (11), (3), (21), (111), (31),
(22), (211), (42), (321). Note the falling factorial is taken in the *ring*, so
`(P_i)_k` is a product of symmetric functions and `x − j` means `x − j·1`.

Everything after that is existing machinery: `p → s` is iterated
Murnaghan–Nakayama, already in `convert`.

Eq (23) divides by `z_γ` and by `i`, so the intermediate ring must be a
[`QAlgebra`]. The composite `st_λ → s` is integral (Theorem 1(2): the inverse
transition has non-negative integer entries, and both are unitriangular), so the
public entry point should hand back integers and the ℚ must not leak into the
signature the way it did on the Macdonald branching route.

**⚠️ Superseded during implementation, and worth keeping.** Eq (23) is correct
and was verified, but reading it as a *formula per λ* is the wrong reading.
Theorem 14's actual content is that **Γ is a linear map with `Γ(p_γ) = 𝐩_γ`**,
and `𝐩_γ = Π_i 𝐩_{i^{m_i(γ)}}` where each factor is a univariate polynomial in
`P_i` of degree `m_i` with leading coefficient `i^{m_i}`. So the change of basis
between monomials in the `P_i` and the `𝐩_γ` is a **tensor product of univariate
triangular matrices**, invertible one variable at a time. What follows in §3.2
and §3.3 was rewritten around that; `src/character_basis.rs` implements the
rewritten version.

### 3.2 `s → st`, by peeling

`s̃_λ = s_λ + (terms of strictly smaller degree)` — Theorem 1(2) gives
`r_{λλ} = 1` and `r_{λμ} = 0` for |μ| > |λ|, so the transition is unitriangular
with respect to size. Hence no matrix inversion is needed:

```text
    R := the Schur element to convert
    while R ≠ 0:
        pick a term c·s_λ of R of maximal degree
        emit c·s̃_λ
        R := R − c·(s̃_λ expanded in s)
```

This needs exactly one primitive, `st_λ → s` from §3.1, and it terminates because
each step kills a top-degree term without creating one.

**⚠️ Not what was built, and the peel is not needed at all.** Γ⁻¹ is available
directly, by the tensor-product structure noted above: expand into monomials in
the `P_i` (using `p_i = Σ_{d|i} d·P_d`, which has no denominators), invert
per-variable into the `𝐩` basis, and reinterpret the resulting indices as `p_γ`.
That is a back substitution against an (m+1)×(m+1) triangular matrix per part
size, not a peel over p(≤20) partitions — and it converts a whole element in one
pass instead of one row at a time, which is what makes §3.3 collapse. The peel
above would have worked and would have been slower; it is recorded rather than
deleted because it is the obvious design and the reason it loses is not obvious.

Sage examples to pin either version:

```text
    s̃_1  = s_1 − s_∅
    s̃_2  = s_2 − 2s_1
    s̃_11 = s_11 − s_1 + s_∅
    s̃_21 = s_21 − 2s_2 − 2s_11 + 3s_1
    s_21 = s̃_21 + 2s̃_2 + 2s̃_11 + 3s̃_1 + s̃_∅
```

### 3.3 The product — the primary route

`s̃_λ · s̃_μ` is the **ordinary** product of symmetric functions (Theorem 7 is a
statement about the outer product), so:

```text
    st_λ ⊗ st_μ  =  peel_to_st( to_schur(st_λ) · to_schur(st_μ) )
```

Cost, for the target case |λ| = |μ| = 10:

- `to_schur` on each side: ≤ P(≤10) = 139 terms — measured at 59 for (6,4).
- the Schur product: ≤ 59×59 ≈ 3.5·10³ small LR products, every one of them a
  product of partitions of size ≤ 10, every one memoized by
  `memo::product_cached`. Output lives in degrees ≤ 20, so ≤ P(≤20) = 2714
  distinct terms.
- the peel: ≤ 2714 iterations, each subtracting one cached `st_λ → s` row.

That is an integer computation of roughly 10⁵–10⁶ term-updates against a Sage
timeout of >90s. ⚠️ Predicted, not measured — this document is written before
the code exists, and the prediction is the thing to check first.

**⚠️ Wrong, and replaced. There is no Littlewood–Richardson coefficient in a
reduced Kronecker calculation.** Because Γ is linear and available in both
directions (§3.1, §3.2), the product can be taken *before* returning to the Schur
basis:

```text
    s̃_λ · s̃_μ = Γ⁻¹( Γ(s_λ) · Γ(s_μ) ),  read in the Schur basis
```

and the multiplication in the middle happens in the **power-sum basis**, where a
product is a multiset union of indices — free. The route above instead pays
~3.5·10³ LR products to arrive at the same answer. Both are correct; the
power-sum one is what `src/character_basis.rs` implements, and the measured
result is in `docs/record/kronecker.md`. The mistake is instructive: the draft reached for the
crate's best-optimized primitive (`Schur::mul`, three backends, memoized) when
the right move was to not need it. Optimized machinery is not an argument for
routing through it.

### 3.4 `ht`, and the independent second route

The `ht` product has a rule of its own, and it is the cross-check.

**Derivation** (mine, from the module description in Eq (6); nothing was read to
obtain it). `h̃_λ` is the character of the permutation module `M^{(n−|λ|,λ)}`. A
tensor product of two permutation modules is the permutation module on the
product of the two coset spaces; its orbits are the double cosets, indexed by
non-negative integer matrices A with row sums `(n−|λ|, λ₁, λ₂, …)` and column
sums `(n−|μ|, μ₁, μ₂, …)`; and the stabiliser of an orbit is the Young subgroup
on the entries of A. Therefore

```text
    h̃_λ · h̃_μ = Σ_A h̃_{sorted multiset of entries of A, minus the (0,0) corner}
```

Only `A₀₀` grows with n, so dropping it leaves a finite, n-free sum. Concretely,
the free data is the block `A_{ij}` for i,j ≥ 1 with row sums ≤ λᵢ and column
sums ≤ μⱼ; the border entries are the slack.

**Verified**: exact agreement with Sage on 9/9 pairs — (1,1), (2,1), (11,1),
(2,2), (21,1), (21,2), (21,21), (31,22), (22,11).

The transition to `st` is the **stable Kostka matrix**, Eq (7):

```text
    h̃_μ = Σ_{|λ| ≤ |μ|} K_{(n−|λ|,λ)(n−|μ|,μ)} · s̃_λ                       (7)
    s̃_λ = Σ_{|μ| ≤ |λ|} K⁻¹_{(n−|λ|,λ)(n−|μ|,μ)} · h̃_μ    (n ≥ 2|λ|)        (8)
```

and the paper gives an **n-free** form of the same coefficient, which is the one
to implement: the coefficient of `s̃_λ` in `h̃_μ` is `Σ_γ K_{γμ}` over partitions
γ ⊢ |μ| such that γ/λ is a horizontal strip of size |μ|−|λ|. Ordinary Kostka
numbers of degree |μ| only — no partitions of 2|λ| anywhere. Checked by hand
against Sage's `h̃_{21} = s̃_∅ + 2s̃_1 + s̃_11 + 2s̃_2 + s̃_21 + s̃_3`: for λ = (1),
γ ranges over (3) and (21) — both have γ/(1) a horizontal strip of size 2 — so
the coefficient is K₍₃₎,₍₂₁₎ + K₍₂₁₎,₍₂₁₎ = 1 + 1 = 2. ✓ For λ = (11) only
γ = (21) contributes, since K₍₁₁₁₎,₍₂₁₎ = 0, giving 1. ✓

Inverting Eq (7) is again unitriangular by size — the |μ|=|λ| block is the
ordinary Kostka matrix — so the same peeling structure as §3.2 applies.

⚠️ **This route does not scale as written and is a checker, not the engine.**
The matrix count explodes on long partitions: for λ = μ = (1^10) the free block
is 10×10 with row sums 1, giving 11^10 ≈ 2.6·10^10 matrices. It is cheap exactly
where the partitions are short — λ = μ = (6,4) is a 2×2 free block, at most
7·7·5·5 = 1225 matrices — which happens to be the shape of the open cases. Use
it for cross-checking two-row and three-row inputs and refuse the rest.

### 3.5 Crate fit

Every component maps to something that already exists. Verified against the
source today, not assumed:

| need | existing | file |
|---|---|---|
| basis types with linear structure | the `basis!` macro + `SymFn` | `sym.rs:132` |
| a hub to convert through | `ToSchur` / `FromSchur` / `convert` | `convert.rs:41` |
| `χ^λ(γ)` for Eq (23) | `character_in::<C>`, `character_table` | `character.rs:72` |
| `z_γ` | `Partition::z()` | `partition.rs:108` |
| `p → s` | `PowerSum: ToSchur`, iterated MN | `convert.rs` |
| the product | `Schur::mul` → `AutoLr` | `sym.rs:229`, `strip_lr.rs` |
| memoized whole products | `memo::product_cached` | `memo.rs:195` |
| Kostka numbers for Eq (7) | `kostka`, `memo::kostka_cached` | `kostka.rs:39` |
| the unitriangular solve pattern | `inverse_kostka_row` | `convert.rs:1063` |
| ring bounds (divide by `z_γ`, `i`) | `QAlgebra`, not `Field` | `coeff.rs` |
| caching `st_λ → s` rows | a new `memo` table, same `table!` shape | `memo.rs:55` |
| bignum + Python boundary | the `guarded` compute-and-escalate path | `guard.rs`, `python.rs` |

Two things it does **not** reuse, and both are deliberate:

- `ops::internal` / `ops::kronecker` (`ops.rs:97`) compute the *ordinary*
  Kronecker product through the power-sum basis at a fixed degree n. They are
  useful as a **stability cross-check** at small n — `ḡ^ν_{λμ}` must equal
  `g((n−|λ|,λ),(n−|μ|,μ),(n−|ν|,ν))` once n is large enough — and useless as an
  engine, because the n we would need is where the 1.1 GB character table lives.
- `Forgotten` is the precedent for a basis that is deliberately not a
  `SymAlgebra` (`sym.rs:174`). `St` and `Ht` are the opposite case: both *are*
  closed under multiplication, and `St::mul` is the entire point.

### 3.6 The single-coefficient query

Lemma 20: for r > 2·deg(f), the coefficient of `s̃_λ` in f is
`Σ_{μ⊢r} (1/z_μ) s̃_λ[Ξ_μ] f[Ξ_μ]`. Both evaluations are integers —
`s̃_λ[Ξ_μ] = χ^{(r−|λ|,λ)}(μ)` by Theorem 1(1), and `h_n[Ξ_μ]` is a count of weak
compositions by Prop 24, so any f expanded in h evaluates integrally.

⚠️ This is recorded as a **candidate**, not a plan. The sum runs over p(r)
partitions with r > 2(|λ|+|μ|), which is p(41) ≈ 4.5·10⁴ for the target case —
plausibly worse than computing the whole product, which §3.3 predicts costs about
the same. Implement §3.3 first, measure, and only then decide whether this is the
`lr_coeff` story ("peel off the larger factor", 23×) repeating for Kronecker or a
dead end. `eval.rs` already evaluates at a finite alphabet, but not at `Ξ`.

### 3.7 A recorded dead end

An earlier draft of this document was going to build the product on Littlewood's
triple-LR formula, `ḡ^ν_{λμ} = Σ_{α,β,γ} c^λ_{αβ} c^μ_{αγ} c^ν_{βγ}`, on the
strength of it being the formula everyone quotes for reduced Kronecker
coefficients. **It is wrong in that generality**, and the check that killed it is
one line: the formula forces |ν| ≡ |λ|+|μ| (mod 2), while `s̃_1 · s̃_1` contains
`s̃_1`. Tested against Sage on 24 pairs, it disagreed on all 24. Recorded because
it is the obvious thing to reach for, and because it demonstrates the rule: the
formula went into a Sage comparison *before* it went into a document.

---

## 4. API

Following the existing shape — one type per basis, conversions through the hub,
an inherent `mul`.

```rust
// sym.rs
basis!(St, "st");        // irreducible character basis  s̃_λ
basis!(Ht, "ht");        // induced trivial character basis  h̃_λ

impl<C: Ring> St<C> {
    /// s̃_λ · s̃_μ, whose structure constants are the reduced Kronecker
    /// coefficients (OZ Thm 7).
    pub fn mul(&self, other: &Self) -> Self;
}
impl<C: Ring> Ht<C> { pub fn mul(&self, other: &Self) -> Self; }

impl<C: Ring> SymAlgebra<C> for St<C> { ... }
impl<C: Ring> SymAlgebra<C> for Ht<C> { ... }

// a new module, character_basis.rs
pub fn reduced_kronecker<C: Ring>(lam: &Partition, mu: &Partition) -> St<C>;
pub fn reduced_kronecker_coeff<C: Ring>(lam: &Partition, mu: &Partition, nu: &Partition) -> C;

/// The second, independent route (§3.4). Cross-check only; refuses inputs
/// whose matrix count would explode.
pub fn reduced_kronecker_via_ht<C: Ring>(lam: &Partition, mu: &Partition) -> Option<St<C>>;
```

Conversions land as `impl ToSchur for St` / `FromSchur for St`, and likewise for
`Ht`, so `convert::<_, St<_>, Schur<_>>` and every other ordered pair against the
existing six bases come for free through the hub.

`degree()` already returns the max over terms (`sym.rs:58`), which is the right
answer for an inhomogeneous element. Nothing in `SymFn` assumes homogeneity —
worth stating because it is the one place these bases could have needed surgery
and do not.

## 5. Correctness requirements

Layered the way the rest of the crate is:

1. **Unit, hand-checkable.** `s̃_{1^r} = Σ(−1)^i e_{r−i}` (Thm 1(3));
   `s̃_1 = s_1 − 1`; Eq (20) `h_{21} = s̃_3 + s̃_{21} + 4s̃_2 + 3s̃_{11} + 7s̃_1 + 4s̃_∅`;
   Eq (21) `h_{1^4}` in the `s̃` basis, which the paper prints in full.
2. **Round trips.** `st → s → st` and `ht → st → ht` are the identity on every
   partition up to size 10.
3. **Algebraic laws.** Commutativity and associativity of `St::mul`;
   `s̃_∅ = 1` is the unit; `ḡ^ν_{λμ} = c^ν_{λμ}` whenever |ν| = |λ|+|μ|
   (stated in the proof of Thm 1(3), and a direct link to the existing LR
   engine).
4. **Two routes agree.** §3.3 against §3.4 on every pair of partitions of size
   ≤ 6, plus the two-row cases up to size 10 where the matrix rule stays cheap.
5. **Stability against the ordinary Kronecker.** `ops::kronecker` at n large
   enough must reproduce `ḡ`. This is the only check that ties the new code to
   the *existing* Kronecker implementation, so it is the one that would catch a
   consistent misreading of the paper. Small n only.
6. **Sage fixtures.** `scripts/gen_sage_oracle.sage` extended with `st` and `ht`
   products and transitions in the range Sage can still answer (degree ≤ 8 for
   products, per §2), committed under `tests/fixtures/`.
7. **Bindings checked separately** (`scripts/check_bindings.py`), per the
   existing rule that a right answer in the wrong slot is a different failure.

## 6. Performance requirements

- `st[5,3] · st[5,3]` and `st[6,4] · st[6,4]` must **complete**. Sage does not.
  That alone is the deliverable; a ratio against a timeout is not a measurement.
- `st[4,3] · st[4,3]`, where Sage takes 30.26s, is the honest like-for-like
  headline. Target ≥ 50×. ⚠️ Guess, stated in advance so the measurement can
  embarrass it.

  *Measured: 3400× (0.0089s). The guess was low by a factor of 68, which is what
  a guess made before knowing the product would leave the Schur basis is worth.
  Both targets met; see `docs/record/kronecker.md`.*
- Degree ladder in |λ|+|μ|, both routes, reported as a growth *curve* and not a
  point — the ladder is what caught three bad conversions before (`convert.rs`
  header), and a point measurement would have caught none of them.
- Report where it stops. Every previous engine here has a wall; find this one's
  and put it in `docs/record/kronecker.md` rather than quoting only the cases that work.

## 7. Open questions

Answered by the implementation, recorded here rather than deleted:

1. ~~**Does the peel dominate?**~~ Moot: there is no peel, and no LR (§3.2,
   §3.3). What dominates instead has not been profiled — the two Γ's are the
   only candidates left, and which one costs more is unmeasured.
2. **Is Eq (23) the right primitive?** Yes as a *map*, no as a *formula*. See the
   correction at the end of §3.1.
3. ~~**Coefficient growth.**~~ Measured, not assumed: the largest coefficient in
   `st[8,5]·st[7,4]` is 32835 — **16 bits**, growing about 1 bit per unit of
   `|λ|+|μ|`. `i128` is not close to being the constraint on the answers.
   ⚠️ This measures the *output*. The intermediate rationals are a different
   quantity, and they are where it breaks: the fixed-width wall is at total
   degree 24, entirely because of `z_γ`. `guard.rs` is now wired in — the engine
   runs over `GuardedRat` and re-runs over `BigRational` under the `bignum`
   feature, reaching degree 32 so far. `docs/record/kronecker.md` has the numbers and the two
   bugs found doing it.

Still open:

4. **`ht` is exposed at all?** It ships public, because it earns its place as the
   independent second route. But `reduced_kronecker_via_ht` refuses long
   partitions (matrix explosion, §3.4), and `s̃_λ → h̃` always involves `(1^k)`
   terms — so the cross-check is limited to small degree in practice, which is
   less coverage than intended. Whether a cheaper independent route exists is
   the open question, and it matters because past `st[4,3]·st[4,3]` there is no
   third-party package left to ask.
5. **The single-coefficient path** (§3.6) is still unbuilt, and now less
   attractive: the whole column is fast enough that a per-coefficient route has
   little room to win. Revisit only if a search driver wants one coefficient of a
   pair it will never ask about again.
