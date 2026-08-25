# Plan: fixing what the convenience-surface review found

Written 2026-08-25, from an interface review of `python/symfn/` conducted
against the built module — every defect below was reproduced, not read off
the source. The review's findings are ordered here into four stages. Stages
1–3 change observable behavior or documented semantics of the surface that
[python.md](../policies/python.md) freezes hardest, so they precede `0.1.0`;
stage 4 is capability work that can follow a release. Each item names the
defect, the fix, and the gate that holds it.

## Why this matters

The layer's stated purpose is that misusing it raises rather than computing
garbage (P7 in [python.md](../policies/python.md)). Several findings are
exactly the failure mode that rule exists against, one layer up from where it
was aimed: an element that prints `0` but does not equal `0`, a constructor
that silently drops a term, and a pairing whose correct answer reads as a bug
to anyone who knows the Jack family is orthogonal — under a different pairing
than the one the method computes.

## The recorded decisions, re-evaluated

Two earlier decisions sat in this plan's way. Both were re-examined on
2026-08-25 rather than followed blind, because the premises of each had
moved; the first is overturned and the second is given an expiry. Two others
were checked and stand unchanged.

1. **Overturned: `to` converts between the six classical codes only**
   ([parametric-basis-inverses.md](parametric-basis-inverses.md), decision
   6). The decision predates the `Sym`/`Param` merge, and its load-bearing
   premise was the two classes: `Sym.to("HLP")` would have returned a
   `Param`, an argument-dependent return type across two types. The merge
   removed that — `to` returns `Sym` unconditionally — and two later
   commitments contradict the decision outright: P10's "the nine tags are
   bases, on the same footing as the six codes", and
   [element-model.md](element-model.md)'s "`McdP` is a basis of the ring: it
   converts to every other basis". The sharpest symptom: `x.to(x.basis)`
   raises `unknown basis 'McdP'` for a tag the same layer's own `basis`
   property just returned — the identity conversion fails. The machinery is
   already built and running: `_back_to` over `_INVERSE` is every
   parametric operation's return leg, so parametric targets are dispatch,
   not new mathematics. Decision 5 of the same file set the precedent — "at
   refuses a parametric tag" was already evolved into expand-first `at` —
   and the file itself says to change a decision everywhere in the same
   commit when one forces a change. Stage 3 carries the work; the family
   `to_*` methods stay, as the convention-pinned, row-accepting entries.
2. **Kept until stage 4, then opened: the `p` conversion for parametric
   coefficients** ([element-model.md](element-model.md), "six-way `to`").
   The recorded grounds were that only the Jack ring is a ℚ-algebra at every
   width, and that opening `p` for one family makes the surface uneven. The
   first is softer than recorded: the s → p transition's coefficients are
   `χ^λ(μ)/z_μ` — parameter-free rationals — so converting a parametric
   element into `p` multiplies its coefficients by rational constants, the
   clear-scale/restore pattern this boundary already uses everywhere; no new
   ring type, only entry points with a scale in the encoding. The evenness
   ground has an expiry: the stage 4 pairings must build exactly this route
   inside the crate for all four rings, and once they land, *closed* is the
   uneven state. So the surface stays closed through stage 3 — with its
   message corrected to the true wall — and opens for all four rings in the
   same stage 4 change that lands the pairings.
3. **Kept: no `__mul__` dispatch change for Jack.** The single-coefficient
   route loses on a whole product (317 ms against 123 ms at degree 12;
   [element-model.md](element-model.md), "The single-coefficient route").
   What remains of that open item is subsumed by stage 4's pairing work —
   see the correction there.
4. **Kept: equality computes nothing.** `__eq__` stays structural: no
   conversion, no ring arithmetic beyond what the coefficient classes
   already answer. Stage 2 narrows what "structural" means for constants; it
   does not add a computed comparison.

## Stage 1 — documentation and messages, no behavior changes

- [x] **Purge the `Param` merge artifacts.** Done 2026-08-25. The b417d8a
      rename left
      sentences that now read as nonsense on the public surface: "a `Sym` or
      a `Sym` in the monomial basis"
      ([_families.py](../../python/symfn/_families.py) near lines 268, 323,
      610, 786), "returns a `Sym` rather than a `Sym`" (`jack.zonal`, near
      line 673), "three legs where `Sym` needs one" (`_hopf`, near line
      1203), the internal row-builder docstrings (near lines 2424, 2471,
      2523, 2577), and "the codes a `Param` accepts"
      ([_bases.py](../../python/symfn/_bases.py), line 36). The distinction
      those sentences carried needs different words, not a rename: "a
      parameter-free element" against "an element with parameters set".
      Gate: a grep for `` a `Sym` or a `Sym` `` and for `` `Param` `` outside
      the `ParamBasis`/`ParamCoefficient` type names comes back empty, and
      `scripts/check_convenience_docs.py` stays green.
- [x] **`Sym.scalar` states which pairing it is.** Done 2026-08-25. It is
      the classical Hall
      pairing, `⟨p_λ, p_μ⟩ = z_λ δ_{λμ}`, and the parametric families are
      not orthogonal under it: `jack.P([2]).scalar(jack.P([1, 1]))` is
      `(1 − α)/(1 + α)`, which is correct and which a reader expecting
      `⟨,⟩_α` will read as a defect — the deformed pairing gives 0 there.
      The docstring says so and carries that value as the doctest, so the
      distinction sits exactly where the confusion arises. No roadmap
      language; when stage 4's siblings exist, they are named here in the
      same change that adds them.
- [x] **The `p` refusal states the true wall.** Done 2026-08-25. The old
      message said the
      coefficient rings "are not all closed" under dividing by z_μ — false:
      z_μ is an integer and every ring here is closed under dividing by one.
      The wall is the boundary encoding (integer numerators, no scale slot
      on the `p` route). The message states the encoding fact in a sentence,
      with no tree path, per P11. It is interim: decision 2 opens the
      conversion in stage 4, and the message lives only until then.

      Two review findings that were items here dissolved when decision 1 was
      overturned, and this records why so they are not re-added. The
      "unknown basis 'McdP'" message and the unfollowable `BasisError` hint
      ("convert one with .to()" for `McdP + McdQ`, which `to` could not
      honor) are both symptoms of `to` refusing parametric targets; once
      stage 3 makes `to` fifteen-way, the first message can no longer arise
      and the hint becomes true for every pair.
- [x] **The LLT class docstring matches the signature.** Done 2026-08-25. It
      said "`G` takes
      a tuple of skew shapes"; `llt_g` takes straight partitions plus
      offsets ([python.rs](../../src/python.rs), `fn llt_g`). The sentence
      is corrected to what the surface accepts; exposing the kernel's
      skew-tuple model is the optional stage 4 item.

Gate for the stage: `scripts/preflight_python.sh`.

## Stage 2 — behavior defects, before 0.1.0

- [x] **Constants and zero compare as the values they are.** Done
      2026-08-25. Fixing the hash sweep surfaced one defect beyond the
      review's list: `AlphaFrac.__eq__` compared `num[:1]` against a number,
      so `3 + 5α == 3` was `True`; the higher coefficients are now required
      to vanish, and the sweep pins the counterexample. Before the fix,
      `q*m([2]) - q*m([2]) == 0`, `hl.P([1])**0 == 1` and `jack.P([]) == 1`
      were all `False` while printing `0` and `1` — an assertion written
      against them was silently always false, the failure mode this layer
      exists to prevent. The semantics after the fix:

      - an element with no terms equals `0` and equals every other empty
        element, whatever the basis and parameters — the reading `_same_ring`
        already takes ("the zero of whichever ring the other one names");
      - an element supported on the empty partition alone compares to a
        scalar, and to another such element, by its one coefficient — the
        empty partition indexes 1 in all fifteen bases, so a constant is the
        same element wherever it is written, and this is also what makes
        equality transitive again (`s([]) * 3 == 3` and `h([]) * 3 == 3`
        were both `True` while `s([]) * 3 == h([]) * 3` was `False`);
      - everything else is unchanged: differing basis or parameters is
        `False`, never an error.

      Hash follows equality, which it did not even for the cases that
      worked: `s([]) * 3 == 3` was `True` with differing hashes, and the
      coefficient classes have the same violation (`QtPoly` of a constant
      equals the `int` and hashes apart). The empty element hashes as `0`; a
      constant element hashes as its coefficient; each coefficient class
      hashes a constant value as the number it equals.

      The `Sym` class docstring is corrected in the same change: it claimed
      comparison raises `BasisError`, and comparison has never raised.

      Gate: the three values above as doctests, and a `check_convenience.py`
      sweep asserting `hash(x) == hash(y)` wherever `x == y` across the
      constant and zero cases of all five coefficient classes.
- [x] **The constructor validates what it accepts.** Done 2026-08-25. Two
      holes in
      `Sym.__init__` ([_sym.py](../../python/symfn/_sym.py)):

      - a repeated shape whose coefficient carries a parameter silently
        last-won — `Sym("m", [((2,), q), ((2,), t)], ("q", "t"))` was
        `t*m[2]` and the `q` term was gone. Numeric coefficients accumulate;
        parameter-carrying ones raise `ValueError` naming the shape, since
        the alternative to accumulating is refusing, not dropping;
      - the coefficient kinds were not checked against `parameters`, so
        `Sym("m", {(2,): q})` constructed and died later as `'<' not
        supported between instances of 'QtPoly' and 'int'` from inside
        `repr`. After the fix: empty `parameters` requires `int`/`Fraction`
        coefficients (the `TypeError` names the class found and the
        `parameters` argument); non-empty requires one coefficient class
        across all terms, with a `Poly`'s variable equal to the single
        parameter; a mix of classes raises rather than `_kind` reading the
        first term and the rest surfacing wherever they surface.

      Gate: a malformed-construction battery in `check_convenience.py`, in
      the style of `scripts/check_python_boundary.py`.
- [x] **`evaluate` refuses a non-integer alphabet in a typed way.** Done
      2026-08-25.
      `s([2]).evaluate([Fraction(1, 2), 1])` leaked PyO3's `'Fraction' object
      cannot be interpreted as an integer` — the exact leak class delta 2 of
      [python.md](../policies/python.md) closed at the contract layer.
      Both routes check the alphabet and raise naming the requirement.
      Rational alphabets as a capability are recorded in stage 4's closing
      note, not silently implied here.

## Stage 3 — ergonomics

- [x] **`to` takes all fifteen codes.** Done 2026-08-25. The work decision
      1's overturn
      opens: a parametric target dispatches through `_INVERSE`/`_back_to` —
      convert to the basis the family's inverse expansion reads, then run
      it — which is the return leg every parametric operation already
      takes, exposed as a conversion. `x.to(x.basis)` is the identity for
      every tag; `s([2]).to("HLP")` equals `hl.to_P(s([2]))`; two
      parametric tags convert through the pivot both expansions share. `p`
      from parametric coefficients stays refused until stage 4 (decision
      2), with the stage 1 message. Ring mismatches keep raising through
      the row builders (`to("HLP")` from `ℚ(q,t)` coefficients refuses, as
      `hl.to_P` does today). The `to` docstring names the family methods as
      the home of each convention rather than restating the pins (P11
      delegation). The dated correction at decision 6 of
      [parametric-basis-inverses.md](parametric-basis-inverses.md) is
      already stamped, 2026-08-25. Gate:
      `check_convenience.py` asserts `x.to(tag) == family.to_*(x.to(pivot))`
      over every tag and the identity `x.to(x.basis) == x` over all
      fifteen.
- [x] **The parameter symbols reach their families.** Done 2026-08-25. An HL
      element printed its coefficients as `t` while `t * hl.P([2])` raised,
      because the exported `t` is the two-variable class and the fix was the
      undiscoverable `t_hl`, which the message did not name. LLT was worse:
      elements print `q` and no exported symbol produced their parameter at
      all. The fix, in four parts:

      - export `q_llt`, a one-variable `Poly` in `q`, beside `t_hl`;
      - `_scale`'s polynomial branch accepts a `QtPoly` supported on a
        single variable when it matches the element's `Poly` variable,
        demoting it — `t * hl.P([2])` and `q * llt.H([1, 1], 2)` then just
        work and stay over their own one-variable ring, which also
        dissolves the `t_hl * x + t * y` mismatch;
      - the Hall-Littlewood row builders (`_t_schur_rows`) demote the same
        way, so `hl.to_P(t * s([2]))` stops refusing the element the
        exported `t` builds;
      - where demotion is impossible the message names `t_hl` or `q_llt`
        and both variables involved.

      `_lift` on a classical element is unchanged: `t * s([2])` stays over
      `ℚ(q,t)`, because which one-variable ring is meant is not decidable
      there. Gate: `check_convenience.py` asserts `t * hl.P(la) ==
      t_hl * hl.P(la)` and the LLT twin across shapes.
- [x] **Division by a rational scalar.** Done 2026-08-25. `s([2]) / 2` was a
      bare `TypeError`
      although `Fraction` coefficients are native. `Sym.__truediv__` takes
      `int` and `Fraction` on both routes (`ZeroDivisionError` on zero), and
      `Poly` and `QtPoly` get the same for coefficient work. `Schub` stays
      over ℤ; its refusal already says why. Gate: `x / 2 ==
      Fraction(1, 2) * x` in the convenience sweep.
- [x] **Extracted coefficients get arithmetic that keeps canonical form.**
      Done 2026-08-25.
      `QtFrac`, `QtRatio` and `AlphaFrac` had no `+`, `-`, `*`, so a
      coefficient pulled out of an element was read-only. They gain the four
      by routing through the element entry points at the empty partition —
      `macdonald_element_add`/`_scale`, `jack_element_add`/`_scale`,
      `macdonald_ht_element_add`/`_scale` — the same reduction `_add` uses,
      so structural `==` keeps working; nothing is computed in Python (P4).
      In the same change `Sym.__mul__` lets a parameter-carrying coefficient
      scale a *classical* element by lifting the element into the
      coefficient's ring through `_lift_to` — today `c * m([1])` for a
      `QtFrac` `c` says "cannot multiply QtFrac with a symmetric function"
      while `c * macdonald.P([1])` works. Gate: build one value along two
      routes and assert `==`; spot-check `(c1 + c2).at(...) == c1.at(...) +
      c2.at(...)`.
- [x] **`principal_specialization` honors its own signature on the classical
      route.** Done 2026-08-25. The parameter is typed `AnyCoefficient` and
      the parametric route takes a `Poly` `q` (there is a doctest), but a
      classical element died in `exact()`:
      `s([2]).principal_specialization(3, q=t_hl)` raised "coefficient must
      be an int or a Fraction". The accumulation no longer forces `exact`
      and returns what the coefficient arithmetic produces.
      Gate: the `Poly("q")` route agrees with
      `principal_specialization_q(n)` coefficient by coefficient.

## Stage 4 — capabilities, not release-blocking

- [x] **The deformed pairings, named as Sage spells them:**
      `Sym.scalar_jack`, `Sym.scalar_t`, `Sym.scalar_qt`. Done 2026-08-25;
      the record entries are "The pairing the bases are orthogonal under"
      in [hall-littlewood.md](../record/hall-littlewood.md), "`⟨·,·⟩_{q,t}`,
      and the power-sum route it forced open" in
      [macdonald.md](../record/macdonald.md), and "`jack_scalar` reaches the
      whole ring" in [jack.md](../record/jack.md).

      **A correction to the review's framing.** `stanley_table` already
      crossed the boundary — it is a `#[pyfunction]` with a doctest
      ([python.rs](../../src/python.rs), `fn stanley_table`) — so "expose
      the batch form" was already done. What remained of
      [element-model.md](element-model.md)'s open item was exactly this
      item's Jack half: `jack_scalar` took `(partition, dense numerator)`
      rows, so it reached integral α-coefficients only — `J` but not `P` or
      `Q`. Landing `scalar_jack` closed that open item, and the same change
      says so there.

      - **Jack first**: `jack_scalar` widened to the full
        `(numerator, atoms, scale, tail)` row encoding the other Jack entry
        points share, and wrapped as `Sym.scalar_jack`. Pinned: `⟨P_λ, P_μ⟩_α
        = 0` for `λ ≠ μ` — the orthogonality that motivated the item, and a
        value the classical pairing gets wrong — the `P`/`Q` duality, and
        the norm against Sage's `scalar_jack` (the `scalarj` fixture rows).
        Two defects surfaced beyond the review's list and were fixed in the
        same change: `_alpha_scalar` dropped a coefficient's `tail`, so
        scaling by a plethysm-produced coefficient multiplied by a different
        value; and `_scalar` accepted a mixed pair in one order only —
        `s([2]).scalar(jack.P([2]))` raised while the reverse worked. The
        shared front leg now lifts whichever side is parameter-free, for
        `scalar` and all three siblings.
      - **Hall-Littlewood and Macdonald crate work**, as planned: both
        pairings are diagonal in `p`, and the crate computes them there —
        `scalar_t` in [hl.rs](../../src/hl.rs), `scalar_qt` in
        [macdonald.rs](../../src/macdonald.rs), each mirroring
        `jack::powersum_scalar` over `Frac`, plus `scalar_qt_ratio` in
        [deltaop.rs](../../src/deltaop.rs) for the `H̃` ring's own
        coefficient field, pinned against the `Frac` form cross-type. No new
        coefficient type: the deformation weights are products of the
        binomials `Frac` already holds factored. Pinned: `⟨P_λ, Q_μ⟩ =
        δ_{λμ}` in both families, `⟨J_λ, J_λ⟩_{q,t} = c_λ·c'_λ`, the
        convention-separating hand values, and Sage's `scalar_t`/`scalar_qt`
        over 39 Schur pairs each (`deformed_pairings_match_sage`, fixture
        rows regenerated with `SAGE_DISABLE_SYMFN=1`).
      - **The same change opened `p` at the surface, for all four rings** —
        decision 2's expiry. `convert_named_to_power` is the sixth
        destination `convert_named` could not carry (`QAlgebra`, because
        only this one divides); `to_power_jack` returns the `JackTerms`
        rows whose own scale slot takes the z_μ division, and `to_power_qt`,
        `to_power_macdonald` and `to_power_ht` return their encodings with
        each numerator coefficient split `(numerator, denominator)`, the
        way the classical `to_power` returns rational coefficients.
        `Sym.to("p")` reaches every ring and the stage 1 interim message is
        gone. Pinned: `to("p")` then `at` equals `at` then the classical
        `to("p")`, and the round trip back into all six family bases,
        across shapes (`check_to_power_parametric`).
      - `Sym.scalar`'s docstring now names the siblings — the sentence the
        stage 1 item deliberately left for this change.
- [ ] **Partial specialization of `at`.** `at` demands every parameter at
      once, so `macdonald.P([2]).at(q=0)` (the Hall-Littlewood
      degeneration) and `qt_kostka(la, mu).at(t=1)` (the one-variable
      polynomial) are unreachable — the standard specializations of these
      families. `at` accepts a subset of the declared parameters; each
      coefficient class substitutes the given ones and keeps the rest.
      The classes are closed under this: a `QtFrac` atom `(a, b)` at
      `t = 1` is `(a, 0)`; a `QtRatio` kind-1 atom `q^a − t^b` at `t = 1`
      is `−(1 − q^a)` with the sign absorbed into the numerator, and
      `a = b` gives the vanishing denominator it should. The element keeps
      its declared `parameters` — narrowing the ring would need coefficient
      classes that do not exist — and an element in a parametric basis
      expands first, as full `at` already does. Pin:
      `macdonald.P(la).at(q=0)` then `.at(t=k)` equals `hl.P(la).at(t=k)`
      across shapes through degree 5 — the q = 0 degeneration is a theorem
      and fails under the `q ↔ t` twist. Symbolic values (`at(q=t)`) are a
      separate decision deferred until this lands; kind-1 atoms under
      `q → t` need `t^a − t^b` factored, and whether that earns its keep is
      not decided here.
- [ ] **Optional: the LLT skew tuples.** The kernel models them
      ([llt.rs](../../src/llt.rs), the \[HHL\] Def 3.2 object) and the
      surface does not. Widen `llt_g` and `llt_min_inv` to accept
      `(outer, inner)` pairs beside plain partitions — permissive-in, P1 —
      and `llt.G` then computes what its mathematics is defined on. If not
      taken, the stage 1 docstring fix already stops the overpromise.
- [ ] **Recorded open, no work planned: rational alphabets in `evaluate`.**
      Exact by per-degree scaling — `s_λ(x/k) = k^{−|λ|} s_λ(x)` on each
      homogeneous component — but nothing asks for it yet, and the stage 2
      typed refusal states the integer requirement.

## What this plan does not change

- Cross-basis equality of nonconstant elements stays `False` — not an
  error, not a conversion.
- The family `to_*` methods stay, as the convention-pinned entries that
  also accept raw contract rows; fifteen-way `to` is a second spelling, not
  a replacement.
- The two refusals [element-model.md](element-model.md) kept — different
  bases do not add, different base rings do not combine — stay exactly as
  they are.

Every stage ends with `scripts/preflight.sh` and
`scripts/preflight_python.sh` green; the stage 4 crate items also carry
their Sage-oracle fixtures, regenerated with `SAGE_DISABLE_SYMFN=1` in the
oracle's environment. Stages 1–3 precede the Phase 2 freeze in
[release-readiness.md](../release-readiness.md); when the first of them
lands, that file's Phase 2 gains the pointer here.
