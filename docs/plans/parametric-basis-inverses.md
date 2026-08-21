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

### 2. ~~Macdonald `s → J`~~ — done 2026-08-21

`schur_to_macdonald_j` (`src/qtkostka.rs`), `symfn.schur_to_macdonald_j`,
`macdonald.to_J`, tagged `McdJ`. Exactly the item as written: group by degree,
call `schur_in_j_table` once per degree, combine. The denominators here *are*
products of `1 − qᵃtᵇ` — they are the hook products `c_μ c'_μ` — so this keeps
the `MacdonaldElement` encoding and `QtFrac`, where `s → H̃` needed the pair
form. `the_schur_table_inverts_the_j_expansion` is the round trip, and the new
`the_j_expansion_is_the_table_read_by_rows` is what catches a transposed read,
which the round trip cannot because it belongs to the table rather than to the
wrapper.

Recorded in [qt-kostka.md](../record/qt-kostka.md), "The element-wise form".

### 3. ~~Macdonald `m → P` and `m → Q`~~ — done 2026-08-21

`monomial_to_macdonald_p` and `monomial_to_macdonald_q` (`src/macdonald.rs`),
the pyfunctions of the same names, `macdonald.to_P` / `to_Q`, tagged `McdP`
and `McdQ`. The item as written: a back-substitution through the new
`macdonald_p_table`, with `m → Q` the same solve divided by `b_λ` applied
factored. What it settled, for Jack:

* **The measure-before-optimizing instruction paid, and not where the item
  said.** The item pointed at `Frac` reduction; the measurement pointed at the
  *number of solves*. The solve is 98% of an `m → P` call and its unit is the
  degree, so a sweep rebuilt it p(n) times — memoized, and degree 8 went
  1.947s to 0.088s. Do the same for Jack from the start.
* **The shared ψ cache in `macdonald_p_table` saves 3%**, measured, not the
  large factor the strip-sharing argument suggests. The table is the
  whole-degree unit, not a faster route to one.
* **`MacdonaldElement` crosses inbound unchanged.** The argument is the
  forward encoding read the other way (`MacdonaldElementArg`), so a `P`, `Q`
  or `J` value feeds straight back. The one thing the parse step must do
  beyond the partition check is reject a `(0, 0)` denominator factor, which
  is the zero binomial and which `Frac::mul_factors` asserts on.

**The open question is settled: there is no `s → P`.** An element in another
classical basis is refused rather than converted, on `BasisError`'s grounds —
a silent conversion picks a basis the caller did not choose and hides its
cost. The caller writes `macdonald.to_P(f.to("m"))`, which is one call, one
visible conversion, and no second solve. The coefficient-slice trick below is
still what a `Param` in another basis would need, and is still not built.

Recorded in [macdonald.md](../record/macdonald.md), "The inverse direction:
`m → P` and `m → Q`".

### 4. ~~Jack `m → P`, `m → Q`, `m → J`~~ — done 2026-08-21

`monomial_to_jack_p`, `_q` and `_j` (`src/jack.rs`), the pyfunctions of the
same names, `jack.to_P` / `to_Q` / `to_J`, tagged `JackP`, `JackQ`, `JackJ`.
The item as written: the Macdonald structure over `AFrac<C>` with
`jack_table(n)` as the triangular input, memoized from the start on item 3's
advice (`memo::jack_p_inverse_cached`, worth 41.5× at degree 10). The
orthogonality route is in `cargo test` as the item asked, not only in the
record. What it settled:

* **The `α = 1` pin really is blind, and so is the `Q`-versus-`P` comparison
  there.** `m_11` is `P_11` outright and `[α(α+1)/2] Q_11`, and both are 1 at
  α = 1. The free-α hand values are the only thing separating either pair.
* **`AFrac::div_int` was missing and the boundary needed it.** `parts` hands
  out `num / (scale · ∏ atoms)` and nothing put a `scale` back. It also lets
  the convenience layer skip `_mac_rows`' least-common-denominator round trip
  — every Jack row carries its own integer denominator already.
* **The Sage confirmation is in `scripts/check_jack.py`**, not a one-off dump:
  three new dump kinds and 351 coefficients through degree 6 in all three
  normalizations, 0 mismatches. That is the gap item 3's record flagged.
* **Adding the arms found two defects in `scripts/bench_inverse.py`** that
  changed numbers already recorded — Sage shares a family's transition matrix
  between its normalizations, and the arms shared a process. The correction to
  item 3's Sage table is in `docs/record/jack.md` and flagged in
  `docs/record/macdonald.md`.

Recorded in [jack.md](../record/jack.md), "The inverse direction: `m → P`,
`m → Q`, `m → J`".

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

## ~~Oracle and fixture policy for the inverse direction~~ — settled 2026-08-21

The round trip against the *fixtured* forward expansion covers the
mathematics and the hand values cover orientation, but neither is an outside
check on the inverse itself: the round trip is blind to any error the forward
direction shares. That is what the fixture is for, and it is now built.
`gen_sage_oracle.sage` was regenerated (`SAGE_DISABLE_SYMFN=1`, the sage-dev
environment — `docs/record/oracles-and-comparisons.md`) with **every** inverse
this plan produced, not the four the item listed:

| tag | what | degrees | coefficients | compared |
|---|---|---|---|---|
| `sinhlp`, `sinhlqp` | `s_λ` in HL `P`, `Q'` | 6 | 225 | exactly, in ℤ[t] |
| `sinht` | `s_λ` in `H̃` | 6 | 190 | cross-multiplied, denominator expanded |
| `sinj` | `s_λ` in Macdonald `J` | 5 | — | cross-multiplied (already present) |
| `minp`, `minq` | `m_λ` in Macdonald `P`, `Q` | 5 | 108 | cross-multiplied |
| `jminp`, `jminq`, `jminj` | `m_λ` in Jack `P`, `Q`, `J` | 7 | 702 | three generic α |

263 new fixture lines; the existing ones are unchanged, so the regeneration is
a pure addition. Four new tests in `tests/sage_oracle.rs` read them, and each
was negative-controlled by swapping the normalization or conjugating the
argument. **Every normalization of every family is present**, because the pair
is what separates them: `m_11` is `P_11` outright in both families, and its
`Q` coefficient is what a dropped or doubled `b_λ` would change.

Two encodings do not reach: `H̃`'s denominators are products of `qᵃ − tᵇ` and
cross the wire expanded, and the Macdonald ones are hook products where a
generic point overflows the `i128` under `Rational` before a pole is reached —
so those three compare by cross-multiplication, which needs no point and no
division. Jack evaluates at three generic α, ⚠️ **never at α = 1**, which
separates neither the normalizations nor the `α → 1/α` twist.

`validation.md`'s "committed fixtures, not scripts someone must remember to
run" is now met outright for the inverse direction in all three families.
