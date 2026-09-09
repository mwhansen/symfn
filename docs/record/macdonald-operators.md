# The Macdonald operator algebra: ∇, Δ_f, Δ'_f, Π, Θ_f

The operators the diagonal-harmonics community works with. Sage has
`nabla` and nothing else, so most of this has no oracle anywhere and is
tied to the one oracled operator by published identities.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

`src/deltaop.rs` — the operators the diagonal-harmonics community works with,
with every formula checked against Sage before it was coded. **Sage has
`nabla` and nothing else**: no Δ, no Δ', no Θ, no Π, and no star
scalar product. That was measured, not assumed — `theta_qt` and `scalar_qt` are
the near-misses and were both identified numerically as different operators (the
classical plethysm `f[X(1−q)/(1−t)]`, and Macdonald's `⟨,⟩_{q,t}`).

Every formula went into a Sage comparison *before* it went into the document:
`scripts/verify_deltaop_formulas.py` implements all five operators on top of
Sage's `Ht` basis from the papers and checks them, and it caught two things a
recollection would have shipped. The HTML of D'Adderio–Mellit (11) renders ∇'s
eigenvalue with a `(−1)^{|μ|}`; against Sage there is none, and the two differ by
a global sign on odd degrees — invisible on `∇e_2`, which is the first case
anyone checks. A direct labeled-Dyck-path enumeration now
verifies both the rise version (a theorem) and the valley version (open) in full
`(q,t)`.

## Before implementation: Sage's wall, on battery

⚠️⚠️ **Both tables below ran on battery and are not the comparison numbers.**
Sage's `∇e_12` at 234.6s here against 126.9s from the mains ladder further down
is one of the three measurements behind the ~1.8× drift recorded under "Power
state" in [README.md](README.md); a ratio taken across these two tables alone
is inflated by roughly that factor. They are kept because they are what the
design was decided against.

SageMath 10.9, this machine, 2026-07-29, battery:

```text
  n   terms   s(∇e_n)   Ht(e_n)   s(H̃_(n))   power
   6     11     0.287     0.196      0.005   battery
   7     15     0.854     0.544      0.009   battery
   8     22     3.014     1.972      0.018   battery
   9     30     9.158     6.064      0.032   battery
  10     42    29.165    20.241      0.056   battery
  11     56    81.506    56.992      0.111   battery
  12     77   234.631   183.545          —   battery
```

`s(∇e_n)` is the total cost of `∇e_n` in the Schur basis; `Ht(e_n)` is the
cost of expanding `e_n` into the `H̃` basis alone; `s(H̃_(n))` is one modified
Macdonald polynomial converted to Schur. Sage's `∇` was not slow because
modified Macdonald polynomials are themselves expensive — one `H̃_μ` reached
the Schur basis in 0.056s at degree 10. It was slow because *expanding an
input into the `H̃` basis* cost 20.2s of the 29.2s that `∇e_10` cost (57.0s of
81.5s at degree 11). **The expensive step here was the change of basis, not
the object being computed** — an earlier sketch of the design assumed the
Macdonald polynomials themselves were the bottleneck and planned to attack
them; they were not.

The crate's own whole-degree `H̃` table was already most of the way there
before any operator existed (`cargo run --release --example bench_htilde`,
same machine, same day, battery):

```text
  n   p(n)   htilde(s)   terms   power
   8     22    0.0279       484  battery
   9     30    0.0864       900  battery
  10     42    0.2938      1764  battery
  11     56    0.8524      3136  battery
  12     77    2.6458      5929  battery
  13    101    7.0548     10201  battery
```

At degree 10 the crate produced all 42 modified Macdonald polynomials in
0.294s against Sage's 20.2s for the change of basis alone — a head start of
roughly 70×, before one line of operator code existed. ⚠️ That 70× is
battery-to-battery and not the number to quote; it is the starting position the
operators were built on top of. The ladder below is the one to cite.

## No linear solve, because H̃ is orthogonal

The obvious route expands the input in `{H̃_μ}` by inverting the `K̃` matrix —
77×77 over ℚ(q,t) at degree 12. It is not needed: `H̃` is orthogonal for the star
scalar product, so the coefficient is `⟨F,H̃_μ⟩_*/w_μ` and the expansion is
diagonal. Better, pairing against the *Schur* basis rather than shape by shape
cancels the `z_ρ` — `⟨F,s_κ⟩_* = Σ_ρ F_ρ χ^κ_ρ ε_ρ W_ρ` has no division in it —
so `H̃_μ` is never converted to the power sums at all. The Schur-basis Gram
matrix of `⟨,⟩_*` was measured integral, symmetric, and divisible by `M` for
n ≤ 7.

## A third factored fraction field, and why the atoms must normalize

Two denominator families arise: `1 − qᵃtᵇ` (from `M`, `Π_μ`, `f[X/M]`, the star
weights) and `qᵃ − tᵇ` (from `w_μ`). `Frac` closes over the first and `bh::Rat`
over the second; this is the first thing in the crate needing both at once, so
`Ratio` holds a denominator as a multiset over the union.

The families **overlap**, and normalizing the overlap is necessary, not
tidiness: `q^a − 1` is `−(1 − q^a)` and `q⁰ − t^b` is `1 − t^b`, and `w_μ`
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

Both sides measured after the implementation landed. ⚠️ The spec's earlier
tables are on battery; mixing them with these inflates every ratio.

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
- **`divide_by_diff`, a specialized chain walk** (1.35×). Sampling put 83% of the
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
measurement meant to support them, and all three were wrong: 2.9× against a
measured 1.30×, 2.4× against 1.00×, and 8.5s against 12.3s. Each was a
prediction about code in this tree, written in the same session that wrote the
code, which is where the estimate felt safest to state without checking.
Corrected in place.

## What checks it

`nabla` and `nabla²` are held to Sage on `e_n` and **every** Schur function
through degree 7 — 1126 coefficients across 144 expansions
(`scripts/check_deltaop.py`). Δ, Δ' and Θ have **no oracle anywhere**, so they
are tied to the oracled one by published identities instead:
`Θ_{e_k} ∇ e_{n−k} = Δ'_{e_{n−k−1}} e_n` and `Δ'_{e_{n−1}} e_n = ∇e_n`, both
exact for every k through degree 7. The rise version of the Delta conjecture is
checked against a direct labeled-Dyck-path enumeration, and `dim DH_n =
(n+1)^{n−1}` falls out of `⟨∇e_n, h_1^n⟩` at q = t = 1. `divide_by_diff` is held
to `QtPoly::divide_exact` on multiples and non-multiples alike, the same way
`frac.rs` holds `divide_by_factor`.

## Theorems get asserts; open conjectures get reports

One rule governs every check in this crate on an object tied to an open
problem: **a theorem gets an assert — a mismatch is our bug — and an open
conjecture gets collected and reported, never asserted or "fixed."** Silently
repairing a mismatch on the open side would convert a potential
counterexample into a passing test, and the crate would never know it
happened.

`src/dyck.rs`'s `Side` enum is the executable statement of it — one type,
two doc comments, each stating the rule at its point of use: `Rise` is
"[HRW]'s rise version — a theorem, so a mismatch is our bug," `Valley` is
"[HRW]'s valley version — open, so a mismatch is a result and must be
reported as one rather than debugged away." This crate's own Delta-conjecture
check applies exactly that split: the rise side is an assertion, the valley
side — still open — is collected and reported instead.

The same asymmetry recurs in at least seven other sites (`src/gj.rs`,
`src/jack.rs`, `src/python.rs`, `examples/bench_jack.rs`,
`examples/bench_llt.rs`, `scripts/README.md`) — this is the first record
entry to state it as a general rule rather than repeat it at each call site.

## The expansion on its own: `s → H̃`

Added 2026-08-21, the second of the inverse expansions; the first is
[hall-littlewood.md](hall-littlewood.md), "The inverse direction". Every
operator above already writes its argument in `H̃` — that is what
`coefficients` does, `c_μ = ⟨f,H̃_μ⟩_* / w_μ` — and nothing exposed it. "Is
this `H̃`-positive" is the question the modified basis is mostly asked, and it
is a question about coefficients in a basis one has to convert *into*.

`schur_to_macdonald_ht` in `src/deltaop.rs` is that half factored out. It
takes a `Schur<QtPoly<C>>` of any mixture of degrees, groups by degree, calls
`coefficients` once per degree, and returns a `BTreeMap<Partition, Ratio<C>>`.
No new mathematics and no new arithmetic: the section above is why there is no
`K̃` inversion here, and `Ratio` is the type the operators already carry the
answer in. Cost is one `bh::htilde_table` per degree present, which one `nabla`
call already pays — and that table is memoized, so p(n) calls at one degree
build it once.

**The coefficients are genuinely not polynomials.** `K̃` is unitriangular in
neither direction and its inverse divides by `w_μ`, whose atoms `qᵃ − tᵇ` do
not cancel — `s_2 = q/(q−t)·H̃_11 − t/(q−t)·H̃_2` at the smallest shape. That
decided the boundary encoding. `MacTerms`, which `macdonald_p` uses, hands the
denominator over *factored* as `1 − qᵃtᵇ` powers, and these denominators are
not products of those; rather than grow that encoding a second atom family —
which every existing consumer of the first would then have to read — the
denominator crosses **expanded**, as an ordinary `QtPoly` term list. So
`schur_to_macdonald_ht` returns `(mu, numerator, denominator)` rows with both
halves polynomials, the denominator never empty and `[(0, 0, 1)]` when the
coefficient is a polynomial. That is the pair form, for the reason above. The
convenience type is `QtRatio` in `python/symfn/_param.py`, and
`macdonald.to_Htilde` returns a `Param` tagged `McdHt`, the name Sage prints.

The numerators arrive integral and the boundary raises rather than rounding if
one does not (`failure.md`, P8). That is not an invariant this code maintains:
it follows from the `z_ρ` cancellation the module docs derive, and it was
checked over every λ through degree 7 before being relied on.

**Pinned by** three tests in `src/deltaop.rs`:
`every_htilde_comes_back_as_itself` runs each `H̃_μ` back through the
expansion for every μ through degree 8 and demands `H̃_μ` alone with
coefficient 1 — the whole `K̃` matrix, and the proof that this is an inverse;
`s2_and_s11_in_htilde_are_the_hand_values` pins the orientation, which the
round trip cannot see, against Sage 10.9 (`SAGE_DISABLE_SYMFN=1`, the sage-dev
environment): `Ht(s[2]) = q/(q−t)·H̃_11 − t/(q−t)·H̃_2` and
`Ht(s[1,1]) = −1/(q−t)·H̃_11 + 1/(q−t)·H̃_2`. All three degree-3 expansions
were compared against Sage in the same session and agree; they are not in the
test. The `q ↔ t` swap gives a different and perfectly plausible answer, which
is the trap
[qt-kostka.md](qt-kostka.md) records as an indexing that was wrong in silence.
The third takes a mixed-degree argument and the zero element.

On the Python side `scripts/check_convenience.py` holds `to_Htilde(Htilde(λ))`
to the unit, the wrapper to the contract rows, and — the check that shares no
code with the solve — specializes the coefficients at `q = 3, t = 2/7` and
recombines them with the `H̃_μ` at the same point, which must rebuild `s_λ`.
All three for every λ through degree 5.

### Measured: 59–132× Sage, and widening

`scripts/bench_inverse.py`, on **AC power**, 2026-08-21. The unit is the
degree — every λ ⊢ n expanded — one process per degree and per arm, because
both sides cache their transition across degrees. The Sage arm ran with
`SAGE_DISABLE_SYMFN=1`, without which it reaches this library and reports
about 1.0×; the script refuses rather than trusting that. Sage dispatches this
to **its own Python**: Symmetrica has no Macdonald bases, so there is no C
arm to lose to here.

```text
  n  p(n)       sage      symfn    ratio
  4     5     0.0837     0.0014     59.8x
  5     7     0.2020     0.0034     59.4x
  6    11     0.9732     0.0126     77.2x
  7    15     3.8273     0.0416     92.0x
  8    22    19.4286     0.1469    132.3x
```

Degrees 1–3 are omitted: Sage's ~0.04s of fixed setup is most of them, and the
trend only means anything above it. The ratio **widens** with the degree —
Sage costs 5.1× more from n=7 to n=8 where this costs 3.5× more — which is the
star-orthogonality expansion against `K̃` matrix machinery.

### A silent truncation the new method inherited

`_schur_rows` in `python/symfn/_families.py` read a `Sym`'s coefficients with
`int(c)`, and `int()` on a `Fraction` **truncates**: a Schur-basis `Sym` with
coefficient `1/2` crossed as 0, so `macdonald.nabla(Sym("s", {(2,): 1/2}))`
returned the zero element rather than an error. Six methods were affected —
`nabla`, `nabla_power`, `delta_ek`, `delta_prime_ek`, `theta_ek`, `big_pi` —
and `to_Htilde` would have been the seventh. These rows are `ℤ[q,t]`, so the
fix is to refuse: the helper now raises `ValueError` naming the shape and the
coefficient. Scaling up by the least common denominator, as `_t_schur_rows`
does for the Hall–Littlewood inverses, would also be exact and is what a
caller may want; it is not done here because each of the seven would need its
own way of dividing the scale back out of a different return type.

## Offline oracle fixture

`∇e_n` in the Schur basis for n = 0..6, committed and checked on every
`cargo test`. ∇ is the one operator here with an external oracle, and Δ_f, Δ'_f
and Θ_f are all tied back to it — so ∇ is the piece that most needed evidence
`cargo test` can re-establish on its own.

## Next

- ~~**A fixture for `s → H̃`.**~~ Done 2026-08-21. `gen_sage_oracle.sage` emits
  `sinht` — `s_λ` in `H̃` for every λ through degree 6, 190 coefficients — and
  `schur_in_macdonald_ht_matches_sage` reads it. ⚠️ These denominators are
  products of `qᵃ − tᵇ`, so the fixture hands them over **expanded** and the
  comparison cross-multiplies; a generic point does not work here, because
  `q = t` is a pole of every coefficient. The hand-checked Rust constants stay
  as the convention pin.
- **The valley Delta conjecture is the point, and the operator is no longer the
  constraint.** `Δ'_{e_k}e_n` is 0.1s at degree 8; the labeled-Dyck-path
  enumeration is what walls out, around n = 9. A search driver wants that
  enumeration written properly.
- ~~**`nabla_e`'s fixed-width guard has no test.**~~ Its rustdoc cited
  `nabla_e_is_exact_in_fixed_width` as the pin; the 2026-08-07 audit found that
  name exists nowhere in the tree, so the citation was removed and the guard
  was uncovered. **Closed 2026-08-20**: the test now exists in
  `tests/bignum.rs` — `i128` held to `BigInt` through degree 8, a different
  width rather than shared mathematics — and the citation is back. The dead
  citation remains the more interesting half: a doc naming a test nobody wrote
  reads exactly like a doc naming a test that passes, and no gate compares the
  two.
- Θ still costs a degree-(n+k) table; whether the composite identities avoid
  ever forming it is unknown.
- Generalize `Frac`, `bh::Rat` and `Ratio` into one `FactoredFrac<A>`. The
  three carry the same design over different atom families, so a fix to one
  has to be found and repeated in the other two; the spec argued for merging
  them before `Ratio` existed.
- ~~Python bindings~~ — done: `nabla`, `nabla_e`, `nabla_power`, `delta_ek`,
  `delta_prime_ek`, `delta_prime_e`, `theta_ek`, `big_pi` and
  `delta_conjecture_side`, all whole-object, all checked against Sage by
  `scripts/check_bindings.py`. `Π⁻¹` is deliberately not exposed: it genuinely is
  not a polynomial, so it cannot cross that boundary; only the composite `Θ` can.
- If 20× is not enough, the next step is a different algorithm — evaluation at
  many `(q,t)` points with modular arithmetic and interpolation, with the degree
  bound `n(n−1)/2` fixing the grid — not more tuning. 73% of the profile is now
  genuine exact division.
