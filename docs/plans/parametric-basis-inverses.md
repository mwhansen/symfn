# Plan: expanding into the parametric bases

Written 2026-08-21, after the Hall–Littlewood pair shipped. This file is the
hand-off for the remaining families. The Hall–Littlewood section of
[hall-littlewood.md](../record/hall-littlewood.md) ("The inverse direction")
records what was built and why; this file says what is left, in what order,
and which decisions are already made so they are not re-made.

## Why this matters

Every parametric family — Hall–Littlewood, Macdonald, Jack, LLT — returned its
polynomials *expanded in a classical basis* and nothing went the other way. A
researcher's question is usually the other way: take a symmetric function
defined by some other construction and write it in `P`, `Q'`, `H̃` or `J` to
see whether the coefficients are positive, or what they count. Positivity
conjectures are statements about coefficients in a basis one has to *convert
into*. Sage users have always had this, because Sage's triangularity machinery
solves for the inverse from symfn's forward expansion; the direct-Python user
had a lookup table rather than a tool. Closing that gap is the point.

## Decisions already made (by the Hall–Littlewood pair)

These were settled when `schur_to_hall_littlewood_p` / `_qp` and `hl.to_P` /
`hl.to_Qp` were built. Follow them unless a family forces a change, and if it
does, change them everywhere in the same commit.

1. **The contract entry point is element-wise, not a table.** It takes a whole
   element in the *forward family's* encoding and returns the same encoding
   in the target basis, so a forward answer feeds straight back in. Mixed
   degrees are accepted and handled degree by degree; the zero element gives
   the empty list. (`schur_in_macdonald_j` is the older whole-degree-table
   shape; it stays, and gets an element-wise sibling.)
2. **The source basis is the one the forward family returns.** Hall–Littlewood
   returns Schur, so the inverse takes Schur. Macdonald `P`/`Q`/`J` return
   monomial, so their inverses take monomial; `H̃` returns Schur, so `s → H̃`.
   Jack returns monomial. A caller holding the element in another classical
   basis converts first with `Sym.to` (or the coefficient-slice trick below).
3. **The crate function returns a plain `BTreeMap<Partition, Coefficient>`.**
   The crate has no `P`-basis type, and a `Schur<…>` holding `P`-coefficients
   is the basis confusion the types exist to prevent. Do not add basis types
   for the parametric families unless a ring structure on them is being
   built too — it is not, and should not be (Sage is that ring).
4. **The convenience layer tags the result with Sage's printed name.**
   `HLP`, `HLQp` today; `McdP`, `McdQ`, `McdJ`, `McdHt` and `JackP`, `JackQ`,
   `JackJ` to come — exactly as Sage spells them, so a reader sees the same
   print form in both systems. The codes live in `ParamBasis`
   (`python/symfn/_types.py`) and `PARAM_BASES` / `check_param_basis`
   (`_bases.py`); `Basis` and `Sym` never learn them.
5. **`Param.at` refuses a parametric tag.** A `Sym` carries only the six
   classical bases, and a parametric basis has no meaning once its parameter
   is fixed. Do not relax this for the "`P_λ(x;0) = s_λ`" special case; it
   would make `at` basis-dependent in a way no other method is.
6. **Methods live on the family object as `to_<Basis>`**: `hl.to_P`,
   `hl.to_Qp`; so `macdonald.to_P`, `macdonald.to_Htilde`, `jack.to_P`. Not
   `Param.to("HLP")` — `to` on `Sym`/`Param` converts between classical bases
   and should stay that.
7. **Escalation and exactness as everywhere else.** The boundary runs guarded
   `i128` then `BigInt` (or `GuardedRat` then `BigRational` where the family
   divides); a surviving denominator where the mathematics promises none is a
   raised error, never rounding (`failure.md`, P8).
8. **Validation is the round trip plus a hand-checked orientation pin plus
   linearity**, in `cargo test`, with the orientation values confirmed
   against Sage once and written into the test. The round trip
   (`to_X(X(λ)) == X_λ` for every λ through degree 8) proves the inverse; it
   cannot catch a transposed matrix on its own, which is what the hand values
   are for. See `s2_in_p_and_s11_in_qp_are_the_hand_values` in `src/hl.rs`.

## What is left, in order

### 1. ~~Macdonald `s → H̃`~~ — done 2026-08-21

`schur_to_macdonald_ht` (`src/deltaop.rs`), `symfn.schur_to_macdonald_ht`,
`macdonald.to_Htilde`, tagged `McdHt`. It was where the plan said it was: the
operators' own `coefficients`, factored out. What the plan left open and this
change settled, for the families still to come:

* **The encoding is the pair form.** The denominators are products of
  `qᵃ − tᵇ`, not of `1 − qᵃtᵇ`, so `MacTerms`' factored denominator does not
  reach; the denominator crosses expanded, as an ordinary `QtPoly` term list,
  and is `[(0, 0, 1)]` rather than empty when the coefficient is a
  polynomial. `HtElement` in `symfn.pyi` names it; `QtRatio` in `_param.py`
  is the convenience type. A family whose denominators *are* `1 − qᵃtᵇ`
  products keeps `MacdonaldElement` and `QtFrac`.
* **The coefficient ring is `Rational`, and there is no escalation.** This
  entry point sits in the operator family, which runs over fixed-width
  `Rational` throughout and refuses rather than wrapping; adding an escalation
  here alone would make one operator behave unlike its neighbors.
* **A round trip through the forward table plus a specialize-and-recombine
  check is enough Python-side evidence.** `check_convenience.py` sets `q` and
  `t` to unrelated values, evaluates the coefficients, recombines with the
  `H̃_μ` at the same point, and demands `s_λ` back.

Recorded in [macdonald-operators.md](../record/macdonald-operators.md), "The
expansion on its own: `s → H̃`", together with a silent `int(Fraction)`
truncation in `_schur_rows` that the work uncovered and closed.

### 2. Macdonald `s → J` (element-wise form of what exists)

`schur_in_macdonald_j(n)` already expands every `s_λ` of a degree in `J`,
by projection rather than a solve (`qt-kostka.md`, "The inverse of `J → s`
is a projection, not a solve"). Add `schur_to_macdonald_j(f)` that groups
`f` by degree, calls the table once per degree, and combines rows; output in
the existing `MacdonaldElement` encoding (denominators *are* hook products
here). Tag `McdJ`; Sage `Sym.macdonald().J()(f)`. The existing test
`the_schur_table_inverts_the_j_expansion` is the round trip; add the
mixed-degree linearity test and the `s_2`, `s_11` pins.

### 3. Macdonald `m → P` and `m → Q`

`P_λ` is monic and dominance-unitriangular in the monomial basis, so `m → P`
is a back-substitution over `ℚ(q,t)` with the `macdonald_p` table as input —
the same shape as `hall_littlewood_p_table`'s loop in `src/hl.rs`, but over
`Frac<C>` instead of `QtPoly<C>`, so every step divides and reduces. Measure
before optimizing: `src/frac.rs` reduction is where Macdonald's time goes
(`docs/record/macdonald.md`). `m → Q` is `m → P` followed by division by
`b_λ`, which is a product of atoms — do it as one function with a flag on the
crate side, two entry points at the boundary. Tags `McdP`, `McdQ`. The
orientation trap is `P` versus `Q` (only `P` is monic) and `q ↔ t`; pin with
`m_11` and `m_2` in both.

Open question to settle here: whether `s → P` should also be offered. It is
`s → m` (integer, `convert_terms`) followed by `m → P`, and the convenience
layer can compose the two *only if* the `Sym` argument is converted by
`Sym.to("m")` before crossing — which it can, since the input is a `Sym`. A
`Param` in `s` over `(q,t)` would need the slice trick below. Decide when the
`m → P` entry exists; do not build a second solve.

### 4. Jack `m → P`, `m → Q`, `m → J`

Identical structure to Macdonald over `AFrac<C>` (the `α`-rational type) with
`jack_table(n)` as the triangular input. Tags `JackP`, `JackQ`, `JackJ`; Sage
`Sym.jack().P()(f)` etc. Orientation trap is `α → 1/α`; pin with `m_11` at
`α = 1` (`P_λ(x;1) = s_λ`) *and* one value with `α` free, since the `α = 1`
check is blind to the twist. `jack_scalar` and orthogonality give a second
route to the same coefficients (`⟨f, P_λ⟩_α / ⟨P_λ, P_λ⟩_α`) that shares no
code with the solve — that is the validation.md row "second engine sharing no
mathematics", and it is cheap here, so use it in the test rather than
only the round trip.

### 5. Not planned: LLT

LLT polynomials `G̃^{(k)}_λ` are not a basis of `Λ` in general (they are
indexed by tuples / `k`-quotients and are linearly dependent across `k`), so
"expand in the LLT basis" is not a well-posed request. The record should say
so in `docs/record/llt.md`'s open tail if anyone asks; nothing to build.

## The coefficient-slice trick (a convenience-layer conversion for `Param`)

A `Param` in a classical basis with coefficients in `ℤ[t]` or `ℤ[q,t]` can be
converted between classical bases without any new Rust: every conversion is
`ℤ`-linear on coefficients, so `f = Σ_k t^k f_k` with integer `f_k`, and
`convert_terms` on each slice gives the answer. That is a loop of contract
calls, which `python.md` P4 permits (it computes nothing the contract layer
does not; compare `Sym.__mul__`). It would let `macdonald.to_Htilde` accept a
`Param` in `m`, and `hl.to_P` a `Param` in `h`. Worth adding as `Param.to`
restricted to classical targets once a second family wants it; not before.

## Surfaces each addition touches (the Hall–Littlewood change as template)

* `src/<family>.rs`: the function, its rustdoc with a Rust doctest pinning the
  orientation, three tests (round trip through degree 8, hand values,
  mixed-degree linearity and the zero element). `src/lib.rs` export.
* `src/python.rs`: the `#[pyfunction]` with ` ```text ` example, escalation,
  registration in the `#[pymodule]` list; `python/symfn/symfn.pyi` stub.
* `python/symfn/_types.py` (`ParamBasis`), `_bases.py` (`PARAM_BASES`),
  `_families.py` (the `to_X` method and, if needed, a reader), doctests in
  each.
* `scripts/check_convenience.py` (round trip, wrapper-equals-contract, tag),
  `check_python_marshalling.py` (shape), `check_python_boundary.py` (a bad
  partition inside the element).
* `docsite/conventions.md` (the family's section gets the inverse example);
  the record file for the family (what was built, the Sage confirmation with
  its environment, what is pinned); this file (strike the item).
* `scripts/preflight.sh` and `scripts/preflight_python.sh` both green.

## Oracle and fixture policy for the inverse direction

The round trip against the *fixtured* forward expansion covers the
mathematics; the hand values cover orientation. When `gen_sage_oracle.sage`
is next regenerated (it needs `SAGE_DISABLE_SYMFN=1` and the sage-dev
environment — `docs/record/oracles-and-comparisons.md`), add a block of
`HLP(s_μ)`, `HLQp(s_μ)`, `Ht(s_μ)`, `J(s_μ)` through degree 6 so the
orientation is fixtured too, and read it from `tests/sage_oracle.rs`. Until
then the Sage values live in the test names and the record, which is the
state `validation.md` calls "committed fixtures, not scripts someone must
remember to run" only half-met; say so in the record when each family lands.
