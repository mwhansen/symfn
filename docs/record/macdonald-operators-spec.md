# Specification: the Macdonald operator algebra (∇, Δ_f, Δ'_f, Π, Θ_f)

*Written 2026-07-29. Implements §2.4 of `docs/research-gaps.md` — "the Macdonald
operator algebra: Δ_f, Δ'_f, Θ_f", the item that document ranks first.*

Sources, all read from the arXiv PDFs on 2026-07-29, not recalled:

| tag | paper | used for |
|---|---|---|
| **[DM]** | D'Adderio, Mellit, *A proof of the compositional Delta conjecture*, [arXiv:2011.11467](https://arxiv.org/abs/2011.11467) (Invent. Math.) | (5)–(8) `M, B_μ, T_μ, Π_μ`; (10) `f*`; (11) ∇; (12) Δ_f, Δ'_f; (21) **Π**; (22) Θ_f |
| **[DIV]** | D'Adderio, Iraci, Vanden Wyngaerd, *Theta operators, refined Delta conjectures, and coinvariants*, [arXiv:1906.02623](https://arxiv.org/abs/1906.02623) | (9) `w_μ`; (27) **Π**; (28) Θ_f |
| **[IR]** | Iraci, Romero, *Delta and Theta operator expansions*, [arXiv:2203.10342](https://arxiv.org/abs/2203.10342) | the star scalar product `⟨F,G⟩_* = ⟨F,(ωG)[MX]⟩` |
| **[HRW]** | Haglund, Remmel, Wilson, *The Delta Conjecture*, [arXiv:1509.07058](https://arxiv.org/abs/1509.07058) | the rise and valley versions; `Val(P)`, `d_i(P)`, `dinv`, `area` |
| **[QZ]** | Qiu, Zhang, *The Schur positivity of ∇m_μ*, [arXiv:2607.00940](https://arxiv.org/abs/2607.00940) | ∇ without a sign; the 2026 result |

**Sage is an oracle here, never a source.** Every number in §2 was measured on
this machine against SageMath 10.9 (release 2026-05-04) today. Every formula in
§1 and §3 was checked numerically against Sage *before* being written down, by
`scripts/verify_deltaop_formulas.py`, which is committed alongside this document
and prints `all formulas verified`. Sage's implementation was not read and must
not be — the same clean-room posture as `docs/cleanroom-spec-skew-lr.md` and
`docs/record/st-basis-spec.md`.

---

## 1. What the operators are

All of them are **diagonal in the modified Macdonald basis** `{H̃_μ}`, which the
crate already computes (`bh::htilde_table`, `qtkostka::macdonald_ht`). That is
the whole reason this is the cheapest item in `docs/research-gaps.md`: no new
enumeration engine is required, only a change of basis and a scalar per μ.

### 1.1 Cell statistics

For a cell `c = (i,j)` of μ with `i` the 0-based row and `j` the 0-based column:

```text
  coarm  a'(c) = j            arm  a(c) = μ_i − j − 1
  coleg  l'(c) = i            leg  l(c) = μ'_j − i − 1
```

and then, [DM] (5)–(8) and [DIV] (9):

```text
  M    = (1 − q)(1 − t)
  B_μ  = Σ_{c ∈ μ}       q^{a'(c)} t^{l'(c)}
  T_μ  = ∏_{c ∈ μ}       q^{a'(c)} t^{l'(c)}   = q^{n(μ')} t^{n(μ)}
  Π_μ  = ∏_{c ∈ μ, c ≠ (0,0)} (1 − q^{a'(c)} t^{l'(c)})
  w_μ  = ∏_{c ∈ μ} (q^{a(c)} − t^{l(c)+1}) (t^{l(c)} − q^{a(c)+1})
```

Note which statistics each uses: `B, T, Π` are **co**arm/**co**leg, `w` is
arm/leg. They are easy to swap and the swap is not loud — `T_μ` built from arms
and legs is still a monomial, still distinct across μ, and still gives a
triangular-looking answer. §5 pins it.

### 1.2 The operators

```text
  ∇ H̃_μ      = T_μ H̃_μ                                        [DM] (11), [QZ]
  Δ_f H̃_μ    = f[B_μ(q,t)] H̃_μ                                [DM] (12)
  Δ'_f H̃_μ   = f[B_μ(q,t) − 1] H̃_μ                            [DM] (12)
  Π H̃_μ      = Π_μ H̃_μ           (μ ≠ ∅;  Π 1 = 1)            [DM] (21)
  Θ_f F      = Π f* Π^{-1} F     (F homogeneous, f* = f[X/M])  [DM] (22)
```

`∇ = Δ_{e_n}` on degree `n`, because `e_n[B_μ] = ∏_c q^{a'}t^{l'} = T_μ`.

**⚠️ A sign, recorded rather than assumed.** The HTML rendering of [DM] (11)
reads `∇H̃_μ = (−1)^{|μ|} T_μ H̃_μ`. Measured against Sage's own `nabla` on every
μ up to n = 5, **there is no sign**: `∇H̃_μ = T_μ H̃_μ`, which is also how [QZ]
states it. The two spellings differ by a global `(−1)^n` on degree n, so they
agree on `∇e_2` and disagree in sign on `∇e_3` — invisible to exactly the first
case anyone checks. We implement the unsigned one.

**Δ' needs no virtual alphabet.** `B_μ − 1` looks like a formal alphabet with a
negative letter, and `Frac` cannot hold one. It does not have to: the cell
`(0,0)` contributes `q^0t^0 = 1` to `B_μ`, so `B_μ − 1` is just `B_μ` with that
cell deleted. `Δ'_f` is `Δ_f` over the other `|μ|−1` cells, and the same cell set
`Π_μ` runs over. Verified for every f in the p-basis, n ≤ 5.

**Θ_f raises degree.** `f ∈ Λ^{(k)}` and `F ∈ Λ^{(n)}` give `Θ_f F ∈ Λ^{(n+k)}`,
because the middle step is an honest **product** of symmetric functions. So Θ is
the only one of these that is not a scalar-per-μ: it is *diagonal, multiply,
diagonal again at a higher degree*.

### 1.3 Why anyone wants this

- **The shuffle theorem** (Carlsson–Mellit) is a statement about `∇e_n`.
- **The Delta conjecture** [HRW] is about `Δ'_{e_k} e_n`. The **rise version is a
  theorem** (D'Adderio–Mellit; Blasiak–Haiman–Morse–Pun–Seelinger, two
  independent routes). The **valley version is open**. Both were verified
  numerically here for n ≤ 5, so the machinery to search past that is exactly
  what is missing.
- **Θ operators** are how the modern refinements are stated, and the identity
  `Θ_{e_k} ∇ e_{n−k} = Δ'_{e_{n−k−1}} e_n` ties all three families together. It
  is verified here for 2 ≤ n ≤ 5 and every k, and it is the single best
  self-consistency test the design has (§5).
- **[QZ] 2026** proves signed Schur positivity of `∇^r m_μ`, settling a 1999
  BGHT conjecture — a live area where the objects are cheap to compute and the
  next conjecture is one search away.

---

## 2. Measured walls

SageMath 10.9, this machine, 2026-07-29. Single runs, per-item `SIGALRM`.
⚠️ Order-of-magnitude walls, not benchmarks — the same caution `docs/record/README.md`
applies to its own tables.

⚠️⚠️ **The tables in §2.2 and §2.3 were taken on battery and are not the
comparison numbers.** Mains power is worth about 1.8× on this machine — Sage's
`∇e_12` is 234.6s on battery and 126.9s plugged in — so every ratio computed
across those two tables is inflated by roughly that factor. §6 carries the
like-for-like ladder, both sides re-measured on mains after the implementation
landed, and that is the one to quote. These are kept because they are what the
design was decided against.

### 2.1 What Sage has

Measured by API introspection and by numerically identifying each method — not
by reading source:

| method | what it actually is | verified |
|---|---|---|
| `nabla()` | ∇, on any basis | `∇H̃_μ = T_μ H̃_μ`, every μ, n ≤ 5 |
| `theta_qt(q,t)` | the **classical** plethysm `f[X(1−q)/(1−t)]` | equals the p-basis substitution, exactly |
| `theta(a)` | `p_ρ ↦ a^{ℓ(ρ)} p_ρ` (the Jack-flavoured one) | exactly |
| `scalar_qt` | **Macdonald's** `⟨,⟩_{q,t}`, weight `∏(1−q^r)/(1−t^r)` | exactly |
| `omega_qt`, `plethysm`, `macdonald().{P,Q,J,H,Ht,S}` | as documented | — |

and what is **absent**: `Δ_f`, `Δ'_f`, `Θ_f`, `Π`, and the star scalar product
`⟨,⟩_*`. There is no `delta`, `Delta`, `delta_prime`, `Theta` or `big_pi`
anywhere on a symmetric function, and nothing in `sage.all` global namespace
either. `theta_qt` and `scalar_qt` are the near-misses — both are genuinely
different operators, and both are the ones a quick look would mistake for these.

So `docs/research-gaps.md` §2.4's claim holds as measured: **∇ exists, everything
else does not.**

### 2.2 Where ∇ stops

```text
  s(nabla(e[n]))                  Ht(e[n])              s(Ht[[n]])
  n     terms      sec            (the H~ expansion)    (one H~, to Schur)
   6       11     0.287                  0.196              0.005
   7       15     0.854                  0.544              0.009
   8       22     3.014                  1.972              0.018
   9       30     9.158                  6.064              0.032
  10       42    29.165                 20.241              0.056
  11       56    81.506                 56.992              0.111
  12       77   234.631                183.545                  —
```

The ratio is a steady ~2.8× per degree and `Ht(e[n])` is 70–78% of it at every
size measured.

```text
  s(nabla(s[2,1]))     0.011      s(nabla(s[5,4]))       6.229
  s(nabla(s[3,2]))     0.050      s(nabla(s[4,3,2]))     6.214
  s(nabla(s[4,3]))     0.583
```

**The finding that shapes the design is the third column.** Sage's ∇ is not
slow because modified Macdonald polynomials are expensive — one of them lands in
the Schur basis in 0.056s at degree 10. It is slow because **expanding an input
into the `H̃` basis** costs 20.2s of the 29.2s, and 57.0s of the 81.5s.

An earlier sketch of this document assumed the Macdonald polynomials were the
bottleneck and planned to attack them. They are not; the change of basis is.
This is the same redirect `docs/record/st-basis-spec.md` §2 records ("transitions are cheap,
the product is not"), arriving at the opposite conclusion for the opposite
reason, which is why it had to be measured rather than carried over.

### 2.3 Where the crate already is

`cargo run --release --example bench_htilde`, same machine, same day. This is
the **whole** table `{H̃_μ : μ ⊢ n}` in the Schur basis — every polynomial the
operators need for one degree, not one of them:

```text
  n   p(n)    htilde(s)      terms
  8     22       0.0279        484
  9     30       0.0864        900
 10     42       0.2938       1764
 11     56       0.8524       3136
 12     77       2.6458       5929
 13    101       7.0548      10201
```

At degree 10 the crate produces all 42 modified Macdonald polynomials in
**0.294s** while Sage spends **20.2s** expanding a single element into that same
basis. The head start is roughly **70×**, and it exists before one line of
operator code is written. That is the honest starting position, and §6 is about
not squandering it.

**Re-measured 2026-07-30, and 1.3× of it is a change made that day.** The table
now reads 0.0210 / 0.0421 / 0.1119 / 0.3245 / 1.0100 / 2.7360 s at n = 8…13 — a
2.6× improvement on the rows above. ⚠️ **Only part of that is attributable**: the
same-session before/after of the one change described below was **1.27–1.34×**
(n = 12: 1.344 → 1.032 s). The rest predates it — other work since the original
table, and possibly a power-state difference, neither of which was controlled
for. The measured 1.3× is the claim; the 2.6× is just where the number stands.

**Where the 1.3× came from.** A sampling profile of
`htilde_table(11)` put **35% of the run inside `QtPoly::divide_exact`**, reached
from `bh::Rat::reduce`'s trial-division loop. `bh` divides by the `qᵃ − tᵇ`
family, and this module had *already* learned twice over that the generic
`divide_exact` is the wrong tool for it — the `diff_may_divide` filter (most
trial divisions fail, and `divide_exact` is expensive about failing, because the
lex-leading monomial `qᵃ` means its early exit never fires on the `t` exponent)
and the `divide_by_diff` chain flow (which stays in a sorted `Vec` rather than a
B-tree with a rebalance per elimination). Neither had reached `bh`. Both now live
in `frac`, next to `divide_by_factor` which is the same job for the other
binomial family, with three callers instead of one.

The rest of the table below is unchanged by this: `∇e_n` spends only ~13% of
itself in `htilde_table`, and the remaining 65% is `Atom::divide`, which is
already through this treatment. ⚠️ Its sort is 15.6% of a `∇e_11` profile and
resisted the obvious fix — `frac::divide_by_diff` records the measurement, so
nobody need repeat it.

---

## 3. The algorithms

### 3.1 The expansion into `H̃`, by the star scalar product

Everything else is a scalar. This is the operation, and it is where a naive
design dies.

The obvious route is to invert the `K̃` matrix — a `p(n)×p(n)` linear solve over
ℚ(q,t), 77×77 at degree 12. **Do not.** The modified Macdonald polynomials are
orthogonal with respect to the star scalar product [IR], which makes the
expansion diagonal and needs no solve at all:

```text
  ⟨p_ρ, p_σ⟩_*  =  δ_{ρσ} · z_ρ · (−1)^{|ρ|−ℓ(ρ)} · ∏_i (1 − q^{ρ_i})(1 − t^{ρ_i})

  ⟨H̃_μ, H̃_ν⟩_* =  δ_{μν} · w_μ

  F = Σ_{μ ⊢ n}  ( ⟨F, H̃_μ⟩_* / w_μ ) · H̃_μ
```

**Verified**: the power-sum form agrees with [IR]'s `⟨F,(ωG)[MX]⟩` on every
Schur pair through n = 4; `⟨H̃_μ,H̃_μ⟩_* = w_μ` and the off-diagonal vanishes for
every pair through n = 5; the round trip `F → coefficients → F` is exact for
`s_λ`, `e_n`, `h_n` and `m_{1^n}` through n = 5.

**The unit of work is the degree, and it is a Gram matrix.** `⟨·,·⟩_*` is a
bilinear form; in the Schur basis its matrix is

```text
  G_{λκ} = ⟨s_λ, s_κ⟩_* = Σ_{ρ ⊢ n} z_ρ^{-1} χ^λ_ρ χ^κ_ρ (−1)^{|ρ|−ℓ(ρ)}
                                    ∏_i (1 − q^{ρ_i})(1 − t^{ρ_i})
```

computed **once per degree** from the memoized character table, and then reused
for every μ and every input. **Measured, not assumed**: `G` is symmetric, its
entries lie in **ℤ[q,t]** — not merely ℚ[q,t] — and every entry is divisible by
`M`, for all n ≤ 7. So the pairing never leaves the polynomial ring, and the
only division in the whole route is the one by `w_μ`.

This is structurally `macop::operator_matrix`: a `p(n)×p(n)` matrix over `QtPoly`
that depends on nothing but the degree, built once and used p(n) times.

### 3.2 The closed form for `e_n`, and when to prefer it

The single most-wanted input has its expansion in closed form, with no pairing
at all. Both of these were verified for n ≤ 6:

```text
  e_n        = Σ_{μ ⊢ n}  M B_μ Π_μ / w_μ · H̃_μ
  e_n[X/M]   = Σ_{μ ⊢ n}          1 / w_μ · H̃_μ
```

So `∇e_n`, `Δ'_{e_k} e_n` and the whole Delta-conjecture family skip §3.1
entirely: `p(n)` scalars, each a product of known factors over `w_μ`. The
general pairing is needed for `∇s_λ`, `∇m_μ` ([QZ]'s object), and for the inner
step of Θ.

### 3.3 The eigenvalues

`f[B_μ]` is `f` evaluated at the multiset of cell monomials. Through the power
sums it is one line and needs no plethysm machinery:

```text
  p_k[B_μ]     = Σ_{c ∈ μ} q^{k a'(c)} t^{k l'(c)}   =  B_μ(q^k, t^k)
  p_k[B_μ − 1] = the same sum over c ≠ (0,0)
  f[·]         = Σ_ρ c_ρ ∏_i p_{ρ_i}[·]   for  f = Σ_ρ c_ρ p_ρ
```

Every eigenvalue is therefore a **polynomial** in ℤ[q,t] when `f` is integral —
`T_μ` is a monomial, `Δ_{e_k}`'s eigenvalue is `e_k` of `|μ|` monomials. Only
`Π_μ^{±1}` and `1/w_μ` are fractions.

⚠️ **The coefficient convention, stated because it is the trap.** `f[B_μ]` and
`f[X/M]` are plethysms *at a formal alphabet*, and q, t are letters of that
alphabet. For `f` with coefficients in ℚ — which is every `f` the literature
applies these to — the distinction is vacuous and `Plethystic::frobenius` is the
identity. For `f` with coefficients in ℚ(q,t) it is not, and the two readings
give different operators. **We define `Δ_f` and `Θ_f` for `f` with constant
coefficients only**, and the API enforces it rather than silently picking a
reading. `qtkostka.rs`'s module docs ("it is not a plethysm, and that matters")
record the crate getting this exact distinction wrong once already, in the other
direction.

### 3.4 Θ, and the one place a product appears

```text
  Θ_f F  =  Π ( f[X/M] · ( Π^{-1} F ) )
```

Three steps: expand `F` (degree n) in `H̃`, divide each coefficient by `Π_μ`, and
sum back into the Schur basis; multiply by `f[X/M]`, an ordinary Schur product
against a symmetric function with ℚ(q,t) coefficients; expand the result (degree
n+k) in `H̃` again and multiply each coefficient by `Π_μ`. `f[X/M]` is the
diagonal substitution `p_k ↦ p_k/((1−q^k)(1−t^k))` — the identical shape to
`qtkostka::invert_s_basis`, down to the factors being `Frac` atoms.

The degree-0 case is a genuine special case, not an accident: [DM] (22) sets
`Θ_f F = 0` when `deg F = 0` and `deg f ≥ 1`, and `f·F` when both are 0. It must
be written explicitly; the general formula divides by `Π_∅` and does not reach
it.

**Verified**: `Θ_{e_k} ∇ e_{n−k} = Δ'_{e_{n−k−1}} e_n` for every 2 ≤ n ≤ 5 and
1 ≤ k < n. This is the test worth having — it exercises Θ's product step, both
Π directions, the Δ' cell-deletion, ∇'s monomial, and two different degrees of
the §3.1 expansion, and it is a published theorem rather than a self-consistency
tautology.

### 3.5 Arithmetic: two factored fraction fields, both already in the crate

Only two families of denominator arise, and each is exactly the closed class one
of the crate's existing types was built around:

| denominator | atom family | existing type |
|---|---|---|
| `M`, `Π_μ`, `f[X/M]`, `⟨,⟩_*` weights | `1 − qᵃtᵇ` | **`Frac`** (`frac.rs:49`) |
| `w_μ` | `qᵃ − tᵇ` | **`bh::Rat`** (`bh.rs:94`) |

`frac.rs`'s module docs explain why the first family works — closed under
products *and lcms*, so addition never needs a gcd; `bh.rs`'s explain the second
in the same words for the other family. This spec is the first thing in the
crate that needs **both at once**.

The recommendation is to **generalise, not to add a third type**: lift the shared
design into a `FactoredFrac<A>` parameterised by the atom family (a trait giving
`atom(a,b) -> QtPoly` and reusing `QtPoly::divide_exact` for reduction), and let
the operator module instantiate it over the union `{1 − qᵃtᵇ} ∪ {qᵃ − tᵇ}`.
`Frac` and `bh::Rat` become instantiations, which is a refactor with two existing
test suites already pointed at it.

⚠️ **The union family is not coprime, and that is fine.** `q² − t²` factors as
`(q−t)(q+t)`, and `q^a − 1` is `−(1 − q^a)`, an atom of the *other* family. So
"lcm = max multiplicity per atom" produces a common multiple that is not always
the least one. Correctness does not depend on minimality — `lift`/`add_assign`
need a common multiple only — and `reduce`'s trial division cuts the excess back
down. `bh::Rat` already lives with exactly this and says so.

### 3.6 Recorded dead ends

Two, both worth stating in advance because both are the obvious design.

1. **Inverting `K̃`.** A `p(n)×p(n)` solve over ℚ(q,t) per degree, when the star
   product makes the expansion diagonal. Not attempted; §3.1 is strictly better
   and needs no fraction-free elimination.
2. **Clearing all denominators at once.** `Σ_μ (num_μ / w_μ) H̃_μ` over a common
   denominator `∏_μ w_μ` is the natural spelling and it will not survive. The
   precedent is exact and in this crate: `macop::solve` first cleared every
   denominator at once and put **48,419 terms** in the common denominator at
   degree 10, against 5,630 terms in the entire operator matrix (`docs/record/README.md`).
   The fix there is the design to copy — `macop::Coeff` (`macop.rs:229`) holds a
   denominator as a *multiset of indices into a fixed family* and calls `reduce`
   after every step, so factors that cancel never get multiplied out. §3.5's
   `FactoredFrac` is the same object one level of generality up.

### 3.7 Crate fit

Every component maps to something that exists. Verified against the source
today, not assumed:

| need | existing | file |
|---|---|---|
| `H̃_μ` in the Schur basis, whole degree | `bh::htilde_table` | `bh.rs:408` |
| `H̃_μ` for one shape | `qtkostka::macdonald_ht` | `qtkostka.rs:190` |
| ℤ[q,t] / ℚ[q,t] coefficients | `QtPoly` | `qt.rs` |
| exact division, no gcd | `QtPoly::divide_exact` | `qt.rs:319` |
| denominators of shape `1 − qᵃtᵇ` | `Frac` | `frac.rs:49` |
| denominators of shape `qᵃ − tᵇ` | `bh::Rat` | `bh.rs:94` |
| factored denominator over a fixed family, reduced eagerly | `macop::Coeff` | `macop.rs:229` |
| `χ^λ_ρ` for the Gram matrix | `character_in`, memoized table | `character.rs` |
| `z_ρ` | `Partition::z()` | `partition.rs:108` |
| `s ↔ p` | `PowerSum: ToSchur` / `FromSchur` | `convert.rs:545` |
| diagonal p-basis substitution (`f[X/M]`) | the shape of `qtkostka::invert_s_basis` | `qtkostka.rs:378` |
| the Schur product for Θ's middle step | `Schur::mul` → `AutoLr` | `sym.rs:262` |
| a per-degree memo table | the `table!` macro | `memo.rs:42` |
| ring bounds (`z_ρ⁻¹`) | `QAlgebra`, not `Field` | `coeff.rs` |
| bignum + Python boundary | `guarded` compute-and-escalate | `guard.rs`, `python.rs` |

Two things it deliberately does **not** reuse:

- **`macop.rs`.** The name collides and the subject does not. That module is
  Lapointe–Lascoux–Morse's *first Macdonald operator* `M₁`, an eigenvector route
  to `J_λ`. It shares `QtPoly` and the factored-denominator idea with this work
  and nothing else. The new module is therefore `deltaop.rs`, not an extension
  of `macop.rs`.
- **`sym::Ht`.** Also a name collision, also unrelated: `Ht` is
  Orellana–Zabrocki's induced-trivial *character* basis (`sym.rs:209`). The
  modified Macdonald functions do not need a basis type of their own here —
  everything enters and leaves in the Schur basis, and `{H̃_μ}` appears only as
  an internal per-degree table. ⚠️ If one is ever wanted publicly it must not be
  called `Ht`.

---

## 4. API

A new module `deltaop.rs`. Everything takes and returns Schur elements, because
that is what the rest of the crate and the Python layer speak.

```rust
// deltaop.rs

/// ∇F. [DM] (11).
pub fn nabla<C: QAlgebra>(f: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>>;

/// ∇^r F, sharing the per-degree table across the powers ([QZ]'s object).
pub fn nabla_power<C: QAlgebra>(f: &Schur<QtPoly<C>>, r: u32) -> Schur<QtPoly<C>>;

/// Δ_f F and Δ'_f F. [DM] (12). `f` carries **integer** coefficients — see
/// §3.3 on why the signature refuses ℚ(q,t) there.
pub fn delta<C: QAlgebra>(f: &Schur<i128>, x: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>>;
pub fn delta_prime<C: QAlgebra>(f: &Schur<i128>, x: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>>;

/// Π F and Π⁻¹F. [DM] (21).
pub fn big_pi<C: QAlgebra>(x: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>>;
pub fn big_pi_inverse<C: QAlgebra>(x: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>>;

/// Θ_f F. [DM] (22); zero on degree 0 unless f is a scalar.
pub fn theta<C: QAlgebra>(f: &Schur<i128>, x: &Schur<QtPoly<C>>) -> Schur<QtPoly<C>>;

/// The shortcut of §3.2: e_n's H̃-coefficients in closed form, no pairing.
pub fn nabla_e<C: QAlgebra>(n: u32) -> Schur<QtPoly<C>>;
pub fn delta_prime_e<C: QAlgebra>(k: u32, n: u32) -> Schur<QtPoly<C>>;

/// The per-degree shared work, exposed because it is the unit of work.
pub fn star_gram<C: Ring>(n: u32) -> Vec<Vec<QtPoly<C>>>;
pub fn htilde_coefficients<C: QAlgebra>(f: &Schur<QtPoly<C>>) -> Vec<Frac...>;
```

Inputs must be homogeneous; a mixed-degree argument is a caller error and should
say so rather than pick a degree. (`SymFn::degree` returns the max over terms,
which is right for `St` and wrong to lean on here.)

## 5. Correctness requirements

Layered the way the rest of the crate is. Every item below **already passes in
the Sage calculator** (`scripts/verify_deltaop_formulas.py`); the requirement is
that the Rust port reproduces them.

1. **Unit, hand-checkable.** `∇e_2 = s_2 + (q+t)s_{11}`; `H̃_{(2)} = s_2 + q s_{11}`
   and `H̃_{(11)} = s_2 + t s_{11}` (already tested in `qtkostka.rs`);
   `Δ'_{e_0} = id` on `e_n`; `Π_{(1)} = 1`.
2. **The statistics.** `T_μ = q^{n(μ')}t^{n(μ)}`, and a test that swapping
   arm/leg for coarm/coleg in `B, T, Π` breaks — §1.1 is the failure mode and it
   is silent otherwise.
3. **Orthogonality.** `⟨H̃_μ,H̃_ν⟩_* = δ_{μν} w_μ` for every pair, n ≤ 8. This
   pins `w_μ`, the star weights and `H̃` against each other simultaneously.
4. **The Gram matrix is integral.** Every entry of `star_gram(n)` lies in ℤ[q,t]
   and is divisible by `M`, n ≤ 8. Measured true for n ≤ 7 already; a ℚ leaking
   in is a bug, not a fallback, and should be an `expect` (the `into_poly`
   policy).
5. **Round trip.** `F → H̃-coefficients → F` is the identity on `s_λ`, `e_n`,
   `h_n`, `m_μ`, every λ and μ, n ≤ 8.
6. **∇ against Sage.** `nabla` on `e_n`, `h_n`, `s_λ`, `m_μ` versus Sage's own,
   through the largest degree Sage will answer (n = 11 for `e_n`, per §2.2).
   This is the only external check available — Δ, Δ' and Θ have **no oracle
   anywhere**, which is precisely why they are worth building and why item 7
   carries the weight.
7. **The Θ/Δ'/∇ identity.** `Θ_{e_k} ∇ e_{n−k} = Δ'_{e_{n−k−1}} e_n` for every
   k < n, n as far as it runs. A published theorem relating all three families;
   see §3.4 on why it is the sharpest test available.
8. **Δ_{e_n} = ∇** on degree n, and `Δ'_{e_{n−1}} e_n = ∇e_n`. Cheap, and they
   catch an off-by-one in the cell deletion.
9. **The Delta conjecture, both versions.** `⟨Δ'_{e_k} e_n, h_1^n⟩` against a
   direct enumeration of labelled Dyck paths [HRW]: the **rise** version (a
   theorem — a mismatch is our bug) and the **valley** version (open — a mismatch
   is a *result*, and must be reported as such rather than debugged away).
   Verified here n ≤ 5; the Rust side should reach further, and that reach is the
   deliverable.
10. **Positivity.** `∇e_n` and `Δ'_{e_k} e_n` have coefficients in ℕ[q,t];
    `(−1)^{|μ|−ℓ(μ)}∇^r m_μ` likewise ([QZ]). Not something the code arranges —
    it arrives over ℚ(q,t) and every denominator must cancel — so it is a real
    check on the whole route.
11. **`dim DH_n = (n+1)^{n−1}`** via `⟨∇e_n, h_1^n⟩` at q = t = 1. One line, and
    it is the check a wrong `T_μ` cannot survive.
12. **Fixed-width exactness**, by running `i64` and `i128` and comparing — the
    `macop.rs` pattern. Coefficient growth here is unmeasured (§7).
13. **Bindings checked separately** (`scripts/check_bindings.py`), per the
    existing rule that a right answer in the wrong slot is a different failure.

## 6. Performance: the requirement, and what it actually did

The requirements as written, before the code existed:

- **`∇e_13` must complete, and then whatever is past it.** Since Sage's growth
  (2.8× per degree) and the crate's `H̃` table growth (3.0×) are within noise of
  each other, the win is a **constant factor and a shifted wall**, not a change
  of asymptotics — worth saying plainly, because a spec that promised otherwise
  would be wrong.
- **`∇e_10` is the honest like-for-like headline. Target ≥ 100×.** ⚠️ A guess,
  stated in advance so the measurement can embarrass it — the `st` spec guessed
  50× and got 3400×.
- **Report the phase split.** §3.6 predicts the final denominator cancellation
  dominates. That prediction is the first thing to check and the most likely
  thing to be wrong.
- **Degree ladder as a growth curve, not a point**; **report where it stops.**

### 6.1 Measured

Both sides re-measured on **mains power**, 2026-07-29, after the implementation
landed: `cargo run --release --example bench_deltaop` against
`s(e[n].nabla())` in SageMath 10.9.

```text
  n   p(n)     Sage (s)   symfn ℤ (s)   symfn ℚ (s)   ratio (ℤ)
  8     22        1.681        0.0632        0.0769        26.6×
  9     30        4.722        0.2000        0.2428        23.6×
 10     42       20.246        0.8197        1.0524        24.7×
 11     56       55.083        2.0076        2.5329        27.4×
 12     77      126.853        8.0598       10.4485        15.7×
 13    101      337.766       16.0083       20.7937        21.1×
```

**The 100× target was missed; the real figure is ~20–27×.** Recorded as the
guess it was. The `st` spec's guess was low by 68× and this one is high by ~4×,
and the difference is instructive: there the design changed underneath the guess
(the product left the Schur basis entirely), while here the design is exactly
what §3 describes and the constant factor is simply what exact division in
ℚ(q,t) costs.

`∇e_13` completes in 16s where Sage takes 5m38s, so the requirement that it
"complete, and then whatever is past it" is met — the wall moves by about two
degrees at a fixed time budget.

### 6.2 The phase split, and four fixes

§3.6 predicted the denominator cancellation would dominate. **It does** — that
much was right. Instrumented at degree 12, the split is `H̃` table 2.6s,
coefficients 0.02s, eigenvalues 0.0005s, and everything else in the
accumulation; inside that, `reduce` was 39.6s against `add_mul`'s 6.2s.

Sampling (`sample`, the workflow `Cargo.toml` documents) then drove four
changes, in order, taking `∇e_12` from 44.4s to 8.1s:

| change | why | effect |
|---|---|---|
| `diff_may_divide`, a one-pass necessary condition | `divide_exact` cannot detect a failed `qᵃ − tᵇ` division early — the leading term is `qᵃ`, so a doomed division runs the whole elimination | 44.4 → 14.8s |
| `divide_by_diff`, a specialised chain walk | the general routine keeps its remainder in a `BTreeMap`; 83% of the profile was in its node rebalancing | 1.35× |
| `QtPoly<i128>` for the closed form | `Rational` where every value is an integer | 1.30× |
| multiply **before** lifting in `add_mul`, and `sum_tree` | a running sum pays a full trial-division sweep at final size for all `p(n)` steps; a balanced tree does most of its work near the leaves | 31.5 → 12.3s |

Two things that did **not** work, recorded rather than deleted:

- **`QtPoly::mul_diff`.** Written to stop `Ring::mul` quicksorting a merge of two
  sorted runs — sound reasoning, and worth 0.13s out of 48s, i.e. nothing. The
  sort really was 30% of the profile, but of a *different* product; reordering
  `add_mul` is what removed it.
- **Breaking the chain walk in `divide_by_diff` as soon as the running sum hits
  zero.** Correct (the flow onward is zero, so the terms ahead are a fresh
  chain) and **slower**, 12.3s → 14.6s: each restart resets the search window to
  the full term list, and losing that shrink costs more than the walking it
  saves.

After all four, 73% of the profile is in `Atom::divide`, and it is now genuine
*successful* divisions — the exact-division work the design is made of, not
overhead around it. Further gains need a different algorithm (evaluation and
interpolation over a modular grid is the obvious candidate), not more tuning.

⚠️ **A methodological note worth keeping.** Three performance claims in this
document and in the source were written *before* the measurement meant to
support them, and all three were wrong (2.9× → 1.30×, 2.4× → 1.00×, 8.5s →
12.3s). They were corrected in place. The rule the house style already states —
measured, not recalled — has to cover numbers about one's own code just as much
as numbers about Sage.

## 7. Open questions

Answered by the implementation, recorded here rather than deleted:

1. ~~**Does the denominator swell?**~~ **No** — measured before a line of Rust
   was written, by simulating whole-atom-only cancellation in Sage (the weaker
   thing a factored denominator can do, not the gcd Sage would use). Reducing
   after every term holds the peak numerator to 1393 terms and the denominator
   to 15 atoms at degree 9, and the denominator cancels to nothing every time.
   That measurement is what licensed the design; the risk `macop.rs` warned
   about is real but the `macop::Coeff` policy defuses it.
2. ~~**Is the closed form faster than the pairing?**~~ Barely: 8.06s against
   10.76s at degree 12, and most of that gap is the [`Ring`] bound letting the
   closed form run over `i128`. The pairing is **not** the cost of the general
   path — both spend essentially all their time in the accumulation they share
   (§6.2), which is why skipping it buys so little.
3. ~~**Coefficient growth.**~~ `i128` holds through degree 13 at least;
   `bench_deltaop` asserts the `i128` and `Rational` ladders agree term for term
   at every degree it runs, which is the `macop.rs` two-width check in another
   spelling. Where it *stops* is still unmeasured, and `guard.rs` is the escape
   hatch when it does.
4. **Θ at higher degree** is still the expensive case, and unimproved:
   `Θ_{e_k} F` for `F` of degree n builds a degree-(n+k) table, so
   `Θ_{e_1} ∇e_9` at degree 10 costs what degree 10 costs. Whether the composite
   identities let one avoid ever forming it is not known here.
5. **Should `{H̃_μ}` become a public basis type?** Still argued against in §3.7.
   Nothing in the implementation needed it; `Ratio` and the Schur basis were
   enough, and the name `Ht` is still taken.

Still open, and now the point:

6. ~~**The valley Delta conjecture.**~~ Partly answered: `src/dyck.rs` now builds
   both combinatorial sides, and **both agree with the operator for every k at
   every n ≤ 9**, as whole symmetric functions rather than at `h_1ⁿ` only. The
   rise version is a theorem, so that half is a check on us; the valley version
   is open, so that half is evidence.
   The constraint moved again, exactly one step: the enumeration is
   `(n+1)^{n−1}` paths, ~1×10⁸ at n = 9 (12 minutes for the full ladder) and 2.4×10⁹ at n = 10, so **the
   combinatorial side is now the wall and a recursive decomposition of the
   generating function is what would move it** — not a faster operator and not a
   faster loop. `docs/record/macdonald-operators.md` has the ladder.
   Worth recording: Sage's `ParkingFunctions` turned out to be a genuine external
   oracle for the `k = n−1` slice (its `dinv`/`area` are [HRW]'s, confirmed
   before use), which is the only independent check this enumeration has.
7. **A different algorithm, if 20× is not enough.** 73% of the profile is now
   genuine exact division. Evaluation at many `(q,t)` points with modular
   arithmetic and interpolation would replace the fraction field entirely; the
   degrees are known in advance (`n(n−1)/2`), so the grid size is known. That is
   a rewrite of §3, not a tuning pass, and nothing here says it would win.