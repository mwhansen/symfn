# The Macdonald operator algebra: ∇, Δ_f, Δ'_f, Π, Θ_f

The operators the diagonal-harmonics community works with. Sage has
`nabla` and nothing else, so most of this has no oracle anywhere and is
tied to the one oracled operator by published identities.

Split out of [ROADMAP.md](../../ROADMAP.md), which carries the phase plan
and a summary of this file.

---

`src/deltaop.rs` — the operators the diagonal-harmonics community works with,
specified in `docs/spec-macdonald-operators.md` and implemented against it.
**Sage has `nabla` and nothing else**: no Δ, no Δ', no Θ, no Π, and no star
scalar product. That was measured, not assumed — `theta_qt` and `scalar_qt` are
the near-misses and were both identified numerically as different operators (the
classical plethysm `f[X(1−q)/(1−t)]`, and Macdonald's `⟨,⟩_{q,t}`).

Every formula went into a Sage comparison *before* it went into the document:
`scripts/verify_deltaop_formulas.py` implements all five operators on top of
Sage's `Ht` basis from the papers and checks them, and it caught two things a
recollection would have shipped. The HTML of D'Adderio–Mellit (11) renders ∇'s
eigenvalue with a `(−1)^{|μ|}`; against Sage there is none, and the two differ by
a global sign on odd degrees — invisible on `∇e_2`, which is the first case
anyone checks. And a remembered closed form for `⟨Δ'_{e_k}e_n⟩` at q=t=1 was
simply wrong, replaced by a direct labelled-Dyck-path enumeration that now
verifies both the rise version (a theorem) and the valley version (open) in full
`(q,t)`.

## No linear solve, because H̃ is orthogonal

The obvious route expands the input in `{H̃_μ}` by inverting the `K̃` matrix —
77×77 over ℚ(q,t) at degree 12. It is not needed: `H̃` is orthogonal for the star
scalar product, so the coefficient is `⟨F,H̃_μ⟩_*/w_μ` and the expansion is
diagonal. Better, pairing against the *Schur* basis rather than shape by shape
cancels the `z_ρ` — `⟨F,s_κ⟩_* = Σ_ρ F_ρ χ^κ_ρ ε_ρ W_ρ` has no division in it —
so `H̃_μ` is never converted to the power sums at all. The Schur-basis Gram
matrix of `⟨,⟩_*` was measured integral, symmetric, and divisible by `M` for
n ≤ 7.

## A third factored fraction field, and why the atoms must normalise

Two denominator families arise: `1 − qᵃtᵇ` (from `M`, `Π_μ`, `f[X/M]`, the star
weights) and `qᵃ − tᵇ` (from `w_μ`). `Frac` closes over the first and `bh::Rat`
over the second; this is the first thing in the crate needing both at once, so
`Ratio` holds a denominator as a multiset over the union.

The families **overlap**, and normalising the overlap is load-bearing rather
than tidy: `q^a − 1` is `−(1 − q^a)` and `q⁰ − t^b` is `1 − t^b`, and `w_μ`
produces both while the star weights produce them independently. Left as
distinct atoms they never cancel against each other.

The reduction policy is the **opposite** of `Frac`'s — reduce after every term,
not once at the end — because the `w_μ` are largely coprime and an unreduced
running sum grows to their lcm. That is `macop::Coeff`'s policy for
`macop::Coeff`'s reason, and it was justified by simulation before any Rust
existed: with whole-atom-only cancellation (weaker than the gcd Sage would use),
the peak numerator is 1393 terms and the denominator 15 atoms at degree 9, and
the denominator cancels to nothing every time.

## Measured: ~20–27×, against a target of 100×

Both sides on mains power, after the implementation landed. ⚠️ Mains is worth
about 1.8× on this machine and the spec's earlier tables are on battery; mixing
them inflates every ratio.

```text
  n   p(n)     Sage (s)   symfn ℤ (s)   ratio
  8     22        1.681        0.0632    26.6×
 10     42       20.246        0.8197    24.7×
 11     56       55.083        2.0076    27.4×
 12     77      126.853        8.0598    15.7×
 13    101      337.766       16.0083    21.1×
```

`∇e_13` completes in 16s where Sage takes 5m38s. **The 100× target was missed**
and is recorded as the guess it was: the `st` spec guessed low by 68× because
its design changed underneath it, and this one guessed high by ~4× because the
design did *not* change and the constant factor is simply what exact division in
ℚ(q,t) costs.

## What the sampling said, and what it cost to ignore it

`sample`, per the workflow in `Cargo.toml`. Four changes took `∇e_12` from 44.4s
to 8.1s:

- **A one-pass necessary condition for `qᵃ − tᵇ` divisibility** (44.4 → 14.8s).
  `divide_exact` cannot fail early here: the leading term of `qᵃ − tᵇ` is `qᵃ`,
  so the "not a multiple" exit never fires on the `t` exponent and a doomed
  division runs the whole elimination. Substituting `q ↦ s^{b'}, t ↦ s^{a'}`
  kills the factor, so a nonzero image proves non-divisibility in one pass.
- **`divide_by_diff`, a specialised chain walk** (1.35×). Sampling put 83% of the
  profile in `divide_exact`'s `BTreeMap` rebalancing — the identical finding
  `frac.rs` records for Macdonald `P`, in the other atom family. Eliminating
  against `qᵃ − tᵇ` is a *flow*: `Q` is a running sum of `N` along chains under
  `σ(x,y) = (x−a, y+b)`.
- **`QtPoly<i128>` for the closed form** (1.30×). `nabla_e`/`delta_prime_e` never
  divide by an integer, so they take a `Ring` bound, not `QAlgebra` — the move
  `qt_kostka_table_via_bh` already makes. Smaller than expected because
  `Rational::add_assign` already short-circuits on integers.
- **Multiply before lifting, and sum as a balanced tree** (31.5 → 12.3s). A
  running sum reaches final size after a few terms and then pays a full
  trial-division sweep at that size for all `p(n)` steps; a tree does most of its
  work near the leaves.

Two failures worth keeping. `QtPoly::mul_diff` was written to stop `Ring::mul`
quicksorting a merge of two sorted runs — sound reasoning, worth 0.13s out of
48s. The sort really was 30% of the profile, but of a *different* product; the
big operand was there because `add_mul` lifted before multiplying. And breaking
the chain walk as soon as the running sum hits zero is correct and **slower**
(12.3 → 14.6s): each restart resets the search window to the full term list.

⚠️ Three performance claims here and in the source were written before the
measurement meant to support them, and all three were wrong (2.9× → 1.30×,
2.4× → 1.00×, 8.5s → 12.3s). Corrected in place. "Measured, not recalled" has to
cover numbers about one's own code too.

## What checks it

`nabla` and `nabla²` are held to Sage on `e_n` and **every** Schur function
through degree 7 — 1126 coefficients across 144 expansions
(`scripts/check_deltaop.py`). Δ, Δ' and Θ have **no oracle anywhere**, so they
are tied to the oracled one by published identities instead:
`Θ_{e_k} ∇ e_{n−k} = Δ'_{e_{n−k−1}} e_n` and `Δ'_{e_{n−1}} e_n = ∇e_n`, both
exact for every k through degree 7. The rise version of the Delta conjecture is
checked against a direct labelled-Dyck-path enumeration, and `dim DH_n =
(n+1)^{n−1}` falls out of `⟨∇e_n, h_1^n⟩` at q = t = 1. `divide_by_diff` is held
to `QtPoly::divide_exact` on multiples and non-multiples alike, the same way
`frac.rs` holds `divide_by_factor`.

## Next

- **The valley Delta conjecture is the point, and the operator is no longer the
  constraint.** `Δ'_{e_k}e_n` is 0.1s at degree 8; the labelled-Dyck-path
  enumeration is what walls out, around n = 9. A search driver wants that
  enumeration written properly, next to `research-gaps.md` §2.3's positivity
  certification.
- Θ still costs a degree-(n+k) table; whether the composite identities avoid
  ever forming it is unknown.
- Generalise `Frac`, `bh::Rat` and `Ratio` into one `FactoredFrac<A>`. Three
  copies of one design is two too many, and the spec argued for it before
  `Ratio` existed.
- ~~Python bindings~~ — done: `nabla`, `nabla_e`, `nabla_power`, `delta_ek`,
  `delta_prime_ek`, `delta_prime_e`, `theta_ek`, `big_pi` and
  `delta_conjecture_side`, all whole-object, all checked against Sage by
  `scripts/check_bindings.py`. `Π⁻¹` is deliberately not exposed: it genuinely is
  not a polynomial, so it cannot cross that boundary; only the composite `Θ` can.
- If 20× is not enough, the next step is a different algorithm — evaluation at
  many `(q,t)` points with modular arithmetic and interpolation, with the degree
  bound `n(n−1)/2` fixing the grid — not more tuning. 73% of the profile is now
  genuine exact division.
