# Plan: one kind of element, in fifteen bases

Written 2026-08-24. This file is the plan for the convenience layer's element
model: what a `Sym` and a `Param` are, and why they should be one thing. It was
Phase 8 of [release-readiness.md](../release-readiness.md) until 2026-08-24,
and moved here because it is a design change with its own dependency order
rather than a release gate — one item in it precedes 0.1.0 and the rest is
additive.

[parametric-basis-inverses.md](parametric-basis-inverses.md) is the direct
predecessor: it introduced the nine parametric tags, and its closing section
records that a tag you can only convert *into* is not yet a basis. This plan is
the rest of that — giving the nine the operations the six classical codes have.

**Not blocking a release. One decision inside it is, and it is the first
item.**

## The model a user brings

A user of this layer holds three independent facts about a value: the **base
ring** its coefficients live in, the **basis** it is written in, and the
element itself. The base ring is chosen once and is a property of the whole
computation; the basis is a way of writing the element down. Neither one
narrows what the value can do, because it is an element of a ring either way.

The `Sym`/`Param` split cuts across all three. `Param` is not a base ring and
it is not a basis: it means "the coefficients are not `int` or `Fraction`",
which is a fact about representation. Nothing in the model a user brings has a
slot for that, and every complaint collected below follows from the one
mismatch.

Sage is where most users of this layer form the model, and it has *more*
element classes than symfn — one per basis — and nobody notices, because they
are interchangeable in every observable way (Sage's own code, run with
`SAGE_DISABLE_SYMFN=1` on 2026-08-22):

    type of P[2]      : MacdonaldPolynomials_p_with_category.element_class
    type of q*m[2]    : SymmetricFunctionAlgebra_monomial_with_category.element_class
    P[2]*P[1]         : -((q^3*t^2-q*t^2-q^2+1)/(-q^3*t^2+q^2*t+q*t-1))*McdP[2,1] + McdP[3]
    s(P[2])           : -((q-t)/(q*t-1))*s[1,1] + s[2]
    P[2].omega()      : ((q^2*t^2-q^2-t^2+1)/(q^2*t^2-2*q*t+1))*McdP[1,1] - ((q-t)/(q*t-1))*McdP[2]
    P[2].degree()     : 2
    (q*m[2]).degree() : 2

**The target, in the user's vocabulary: the nine tags are bases like `s` and
`m`, the parameters are the base ring, and neither one makes a different kind
of object.**

**A parametric basis is not narrower than a classical one.** `McdP` is a basis
of the ring: it has products, it converts to every other basis, and ω is
defined on it — the values above are Sage's. So a refusal in this layer that
cites a basis being parametric is a statement about symfn's coverage, never
about the mathematics, and the message has to say which it is.

`q * m([2])` is the sharpest case and the one to reason from. Multiplying by a
scalar is the most innocuous thing that can be done to a ring element, and here
it deletes `degree`, `is_homogeneous`, `omega`, `antipode`, `**`, `to` and
every product — from a value in the monomial basis with no restriction on it.
The product refusal is worse than absent: it cites a parametric basis's
structure constants for an element that is not in one, because the check is
`isinstance(other, Param)` with no basis test.

## The decision

- [ ] **`q * m([2])` returns whatever `m([2])` returns, and the end state is
      one class.** In the user's model there is no second kind of element, so
      the return type of a scalar multiplication cannot be the place a
      distinction appears. `Sym` carries `Poly`, `QtPoly`, `QtFrac`, `QtRatio`
      and `AlphaFrac` alongside `int` and `Fraction`; the nine tags join the
      six codes as bases it can be written in; `Param` becomes an alias for it.
      `_lift` in [python/symfn/_sym.py](../../python/symfn/_sym.py) is the site.

      **What precedes 0.1.0 is the commitment, not the merge.** An earlier
      draft of this phase held that the return type could not change in 0.1.1
      without breaking a caller. That is weaker than it looked. If the merged
      class is named `Sym` and `Param` is kept as an alias of it, then
      `isinstance(x, Param)` stays true of everything it is true of today, and
      the only observable change is that it becomes true of elements it used to
      be false of — which is exactly the belief this phase is retiring. So the
      real ordering constraint is that the two interfaces converge — same
      accessors, same operations, same error shapes — after which the merge is
      invisible to anything but `type()`. The commitment is written into
      [docs/policies/python.md](../policies/python.md) (P10) as of 2026-08-24,
      which is the part that precedes 0.1.0; the merge lands once the items
      below have removed the reasons to tell the two apart.

      The cost, whenever it lands, is that `Coefficient` in
      [python/symfn/_types.py](../../python/symfn/_types.py) becomes a public
      union of seven types, in `symfn.pyi`, under `mypy --strict`, and read by
      the adapter.

## The work, ordered by what a user hits first

- [x] **Accessors that touch no coefficient arithmetic**: `degree`,
      `is_homogeneous`. Pure Python over `terms`, no contract call, no entry
      point. First because they are the entire gap for anyone who only ever
      scales a classical element. Done 2026-08-24, same bodies as `Sym`'s;
      `macdonald.P([2]).degree()` is 2 without expanding, which is the point —
      the degree is a fact about the partitions alone.
- [x] **Six-way `to`**, done 2026-08-24, for all nine tags and every classical
      basis but `p`.

      `Param.to` used to refuse everything but a parametric basis's own pivot,
      so the only route from a family to the Schur basis ran through
      specializing the parameter first. It now takes any of `s`, `h`, `e`, `m`,
      `f`, and a parametric basis reaches them by expanding into its pivot and
      converting from there. `jack.P([2]).to("s")` is
      `s[2] + ((1 − α)/(1 + α))·s[1,1]`, the value this item predicted, and
      `macdonald.P([2]).to("s")` matches Sage's `s(P[2])` term for term.

      **It is four boundary entry points, not the two this item first said.**
      A `Param` carries five coefficient classes over four Rust rings, each
      with its own encoding: `QtPoly` for Hall-Littlewood, LLT and any scaled
      classical element (`convert_qt_terms`), `Frac` for Macdonald `P`, `Q`,
      `J` (`convert_macdonald_terms`), `AFrac` for Jack
      (`convert_jack_terms`), and `Ratio` for `H̃` (`convert_ht_terms`). All
      four are the same routing over a different ring, because
      `crate::convert_named` is generic in `C: Ring` — `convert` resolves its
      route from the *types* of its two ends, and a caller holding a basis code
      has no types to offer, so the same routing is restated with the pair
      resolved at runtime. The three rational ones reuse the row builders their
      family's inverse expansion already had.

      **`p` is not among the five, and it is not simply a coverage gap.**
      Conversions into the power-sum basis divide by z_μ and so need a ring
      containing ℚ. `QtPoly<i128>` is a `Ring` and nothing more, and
      `Frac<Guarded>` is a `QAlgebra` only when its own coefficients are — so
      of the four, only `AFrac` qualifies at every width. Opening `p` for Jack
      alone would make the surface uneven; it is left closed, and the integer
      path states the same restriction and routes `p` through `to_power`.

      This is also what makes `to` between two *parametric* bases reachable,
      since both sides expand into a classical one.
- [x] **`omega` and `antipode`**, done 2026-08-24 for the same coefficient
      classes `to` reaches. Both are defined on the Schur basis — ω conjugates
      the index, the antipode conjugates and signs — so the boundary pair
      `omega_qt_terms` and `antipode_qt_terms` needs no basis argument and no
      escalation: conjugation is a bijection on the partitions of a degree, so
      no two terms meet and the coefficient ring is never added in. The one
      place a sign is not total is `i128::MIN`, where `Coeff::negated` widens.

      An element in a family's own basis comes back in it, which is three legs
      rather than `Sym`'s one: expand into the pivot, convert to Schur, act,
      and travel back. `hl.Qp([2]).omega()` is `HLQp[1,1] − t·HLQp[2]`, and
      that last leg is an inverse expansion the classical route never runs, so
      `check_parametric_hopf` asserts the basis as well as the value. It works
      for `HLP`, `HLQp` and `McdHt`.
- [ ] **ω and the antipode over the rational-function rings.** Both act in the
      Schur basis and the entry point that does so reads the polynomial
      encoding, so `macdonald.P([2]).omega()` refuses where `to` no longer
      does. Two ways out, and the choice is not made: three more entry-point
      pairs on the model of the converters, or the observation that ω sends
      `h_μ` to `e_μ`, so it is a relabeling of the basis tag with `to` on
      either side and needs no new boundary at all. The second is free and puts
      a mathematical identity in the convenience layer, which P4 in
      [docs/policies/python.md](../policies/python.md) is about not doing.
- [ ] **Products, in every basis.** This is the item that makes the nine tags
      bases rather than labels. A product in a parametric basis is reachable
      and always was: expand to the pivot the basis expands in, multiply,
      invert back — which is what
      [python/symfn/_sym.py](../../python/symfn/_sym.py) already does for `h`,
      `e`, `p` and `f` through Schur. Both directions exist for all nine tags.
      A product in a classical basis with parametric coefficients is the same
      route without the inverse, so the multiply is the only new capability and
      one entry point per coefficient kind covers both.

      **Hall-Littlewood is the case to build first.** Its pivot is Schur, so
      the multiply is the memoized Littlewood-Richardson path, and its
      coefficients are polynomials in `t` rather than rational functions —
      no denominators and no reduction. Measured 2026-08-22 (release, AC
      power, caches cleared per case, Apple M4, scratch example not kept):

      | `P[n]²`, degree | HL invert | Jack invert |
      |---|---|---|
      | 8 | 0.30 ms | 7.0 ms |
      | 12 | 7.5 ms | 119 ms |
      | 16 | 82 ms | overflow |
      | 20 | 676 ms | overflow |
      | 24 | 7.6 s | overflow |

      The multiply is never the cost — 2.4 ms at degree 24 against 7.6 s for
      the inverse, which is a `p(n)` triangular solve over every partition of
      the product degree and does not care that the answer has 13 terms.
- [ ] **Exponentiation follows from products** and needs nothing else:
      `Sym.__pow__` is repeated squaring over the same multiply.
- [ ] **The single-coefficient route, for Jack.**
      [src/jack.rs](../../src/jack.rs)'s `jack_structure_constant` computes
      `⟨J_λ J_μ, J_ν⟩_α` in the power-sum basis, where the product is a
      multiset union and the pairing is diagonal, so no basis change of the
      product is formed. It wins for **one** coefficient and loses for a whole
      product: 9.4 ms against 119 ms at degree 12 for a single constant, but
      317 ms against 123 ms for the full expansion, because pairing against
      every `J_ν` of the degree needs all 77 of their p-expansions while the
      triangular solve returns all 77 coefficients at once. `stanley_table` is
      the amortized form and its own doc measures the redundancy at 282× for
      k = 6. So: p-route for a single coefficient or a whole-degree batch,
      inverse for a one-off product. Nothing analogous is exposed for
      Macdonald — [src/hl.rs](../../src/hl.rs) and
      [src/macdonald.rs](../../src/macdonald.rs) have no scalar product, norm or
      power-sum route — though `⟨,⟩_t` and `⟨,⟩_{q,t}` are diagonal in `p` for
      the same reason.
- [ ] **The overflow is in the product, not the expansion**, and it needs a
      row of the mechanism table in
      [docs/policies/failure.md](../policies/failure.md). `P[n].to("m")` is clean
      to at least `P[20]` (627 terms); `P[8]²` at degree 16 panics with
      `attempt to multiply with overflow` at [src/coeff.rs](../../src/coeff.rs).
      Coefficients in ℚ(α) grow much faster than the integer coefficients a
      classical product produces, so this route leaves `i128` at a degree
      where `Schur<i64>` is comfortable. The boundary row — two-pass
      escalation over `BigInt`/`BigRational` — is the one that applies, so the
      Python entry point escalates rather than refusing; a Rust caller keeps
      the loud panic.
- [ ] **Pin the normalization before shipping any of the products.** Sweeping
      every `P_μ · P_ν` with `|μ| = |ν| ≤ 5` gives 1871 Hall-Littlewood
      coefficients, **331 of them negative** — `P[2,1]² → P[3,1,1,1]` is
      `1 + t − t³ − t⁴`, whose constant term is the Littlewood-Richardson
      coefficient `c^{3111}_{21,21} = 1` and whose value at `t = 1` is 0,
      matching `m[2,1]²` having no `m[3,1,1,1]` term. So these are in ℤ[t],
      while the classical Hall polynomials counting subgroups of abelian
      p-groups are in ℕ[t]; the two differ by a normalization twist that is
      not restated here from memory. A doctest pinning the convention must use
      a coefficient where the twist shows — a negative one — and not
      `P[1]² = P[2] + (1 + t)·P[1,1]`, which every convention in circulation
      agrees on. This is the trap CLAUDE.md names as the house failure mode.
- [ ] **Evidence**, per the table in
      [docs/policies/validation.md](../policies/validation.md): products in these
      families are a row Sage covers, so an offline fixture sweep against it
      plus the specialization pins — `t = 0` to Littlewood-Richardson, `t = 1`
      to monomial, `α = 1` to Schur — plus the convention pin above. The
      specializations are the check that shares no mathematics with the route.

## The refusals to keep, and they are the only two

Everything else this layer refuses today is coverage, and its message should
say so rather than naming the basis's parameters as the reason.

- [ ] **Different bases do not add.** `s([1]) + m([1])` raises `BasisError`
      here; Sage silently coerces, returning `McdP[1,1] + McdP[2]` for
      `P[2] + m[1,1]` in the run above. The refusal is the better answer and
      stays, for the nine tags exactly as for the six codes.
- [ ] **Different base rings do not combine, and the error has to say that.**
      `alpha * m([2]) + q * m([2])` raises `TypeError: unsupported operand
      type(s) for +: 'Poly' and 'QtPoly'` — a coefficient-class accident
      leaking through a domain question, for two elements in the same basis.
      It needs the same treatment `BasisError` gets: one error type, naming
      ℚ(α) and ℚ(q,t). `Param.__add__` compares `_params` already but only
      reaches that check when the bases differ, which is why the mismatch
      above escapes it.
- [ ] Audit every remaining refusal against these two. The product message is
      the known offender; the `to` message — "already classical; substitute
      with at() and convert the result" — is the other, and it stops being
      true the moment six-way `to` lands.

## Deferred, with no work planned

`scalar`, `skew_by`, `coproduct`, `expand`, `evaluate`,
`principal_specialization`, `dimension`, `plethysm`. Each is a pair of entry
points on the same pattern; plethysm over ℚ(q,t) is the only one that is real
work rather than wiring. These are the reason the merge does not by itself
deliver the target: a merged `Sym` would carry all eight as methods that raise
for parametric coefficients, which relocates the refusal rather than removing
it. That is an argument for doing the work above *before* the merge, not for
keeping two classes. Revisit once `to` is six-way, which is expected to absorb
most of what the gap feels like.

**Done when:** an element with parameters answers every question an element
without them answers, in all fifteen bases, and the only refusals left are the
two above.

There is no Pieri rule in the tree for `P_λ · P_(r)`, which has a closed form
for Jack (Stanley 1989) and would beat both routes for a one-row multiplier
without touching anything of size `p(n)`. `jack_p_branching` is a different
object — branching for the m-expansion, not a product rule. Left open here
rather than in a record tail because it is a capability, not a dead end.
