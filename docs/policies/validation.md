# Validation policy — the evidence a new capability must carry

What this file governs: how much independent evidence a new module, public
family, or engine must have before it ships, and what counts as independent.
What it does not govern: how the evidence is written up — doctest form, the
circulating conventions, fixture provenance and the pre-ship checklist are
[style.md](../style.md); failure mechanics are [failure.md](failure.md); the
Python surface and its layers are [python.md](python.md), whose Sage-free
boundary suite is where this file's edge pins land on that side (V6); and
measurement methodology and the numbers themselves belong to the record
([oracles-and-comparisons.md](../record/oracles-and-comparisons.md) holds the
harness discipline).

Like the other policies, this is a rulebook, not a description of current
practice; where the two disagree, the gap is listed in
[What this changes](#what-this-changes). Unlike the others, most of it has
been in force for months as folklore: the record cites "the crate's standing
policy" ([qt-kostka.md](../record/qt-kostka.md)), and five source comments
restate the core rule locally ([lr.rs](../../src/lr.rs),
[charge.rs](../../src/charge.rs), [kf.rs](../../src/kf.rs),
[jack.rs](../../src/jack.rs), [gjmod.rs](../../src/gjmod.rs) — the last citing
another module as "the `qtkostka.rs` standard"). A policy that exists only as
scattered quotes has no home to extend; this file is the home.

## The invariant

**Every value a public family can produce is covered by at least one check
that does not share its mathematics — over the range the documentation
advertises.**

The two clauses block the two ways a plausible wrong value ships. Without the
first, a family validates itself: the Kronecker product "rests entirely on
p-basis diagonality plus s ↔ p, so a wrong identity would be self-consistently
wrong" ([kronecker.md](../record/kronecker.md)) — which is why it is
separately checked against the character-theoretic definition. Without the
second, correctness established at toy sizes gets asserted at real ones: the
5 313 471-term LR product was for a long time "checked only against our own
conjugate orientation — a real consistency check, but not an independent one"
([littlewood-richardson.md](../record/littlewood-richardson.md)).

The invariant is a floor. How much evidence a family needs above it is set by
one gradient: **the less the outside world can check, the more of the evidence
must come from inside this tree.** A family Sage computes on the whole intended
range ships on a fixture sweep, a convention pin, and its laws; a family
nothing computes — the Goulden–Jackson tables — ships with two engines sharing
no mathematics, specialization pins against objects with independent
definitions, and its theorems enforced as assertions
([jack.md](../record/jack.md)). Rigor is allocated by checkability and by
what a wrong answer would cost, not by how large the feature feels.

The reason is the failure policy's reason: who the caller is. The
researcher publishes under their own name, and the failure that destroys them
is the plausible wrong value ([style.md](../style.md), "The readers"). Every
check this crate performs is one the researcher no longer has to.

## What counts as evidence

Five classes, in descending strength. A family carries the strongest class
that exists for it, plus enough of the lower ones to cover what the top class
cannot see.

1. **An external oracle.** Sage, lrcalc, Symmetrica, schubmult — an
   implementation sharing no authorship with this one. Two forms, not
   interchangeable: an **offline fixture**, committed with its generating
   script, is durable evidence `cargo test` re-establishes on every run; a
   **live harness** (`scripts/check_*.py`, the compare scripts) is wider but
   "only runs when someone has Sage and remembers"
   ([jack.md](../record/jack.md)). V4 governs the split. The exact side of
   any comparison runs wide, per [failure.md](failure.md) R10.
2. **An independent in-tree route.** A second engine sharing no code with the
   first — and preferably no mathematics: charge against the Morris recursion
   ([charge.rs](../../src/charge.rs)), the character sum against the p-basis
   product ([kronecker.md](../record/kronecker.md)), three (q,t)-Kostka
   routes "sharing nothing above `Partition`"
   ([qt-kostka.md](../record/qt-kostka.md)). Agreement between such routes is
   evidence; between routes sharing a core it is only consistency (V3).
3. **Identities and specializations.** A published identity tying the new
   object to an oracled one — `q = 0` collapses Macdonald P to
   Hall–Littlewood P, "two routes with nothing in common"
   ([macdonald.md](../record/macdonald.md)); Δ_f, Δ′_f and Θ_f are tied to
   ∇, the one operator with an oracle
   ([macdonald-operators.md](../record/macdonald-operators.md)) — or a
   specialization equal to an object with its own definition: the GJ tables
   at `b = 0` (the class algebra, from characters alone) and `b = 1`
   (matchings, counted directly — "computing this side from zonal
   polynomials would have re-used Jack at α = 2 and checked nothing",
   [jack.md](../record/jack.md)). Each is chosen for what it uniquely pins:
   `ω_α` duality guards `jack_q` because it is "the only law relating P to Q,
   to conjugation, and to the parameter inversion at once", where
   `⟨P,Q⟩ = δ` "would survive a compensating error in both".
4. **Algebraic laws.** Ring homomorphisms, round-trip identities, ω an
   involution, the Hopf axioms
   ([algebra_laws.rs](../../tests/algebra_laws.rs)). Powerful, cheap, and
   convention-blind: every circulating normalization of a family satisfies
   its laws, so laws never discharge the convention question.
5. **Pins.** Hand-checkable values, including the mandatory convention pin
   ([style.md](../style.md), "Examples are convention pins"). The floor: the
   only class that distinguishes conventions, and never sufficient alone.

## The rules

### V1 — A published formula is a claim until a machine has checked it

The literature carries misprints and the implementer carries false memories,
and both produce code that faithfully computes the wrong thing. A formula
taken from a paper is verified numerically against an oracle **before** it is
implemented: `scripts/verify_deltaop_formulas.py` checked all five operators
on top of Sage's `Ht` basis before `deltaop.rs` existed and "caught two
things a recollection would have shipped" — a rendered `(−1)^{|μ|}` in a
published eigenvalue, invisible on `∇e_2`, "the first case anyone checks",
and a wrong remembered closed form
([macdonald-operators.md](../record/macdonald-operators.md)). When sources
disagree, the normative one is named with the reason — [KMS] over [LLT] §7,
whose printing of the straightening rules carries two misprints
([llt.rs](../../src/llt.rs)). Witnesses get the same treatment: the LLT
min-inv floor's witness "was written down as λ = (2,2) twice before anyone
checked", and the true smallest is (2,2,2) ([llt.md](../record/llt.md)).

### V2 — Know the incumbent before building against it

What computes the object today, under which convention, out to what wall —
measured, not assumed. Reading Symmetrica's C before Hall–Littlewood changed
the plan and freed charge to be the oracle
([hall-littlewood.md](../record/hall-littlewood.md)); a report on Sage's
internals showed the (q,t)-Kostka comparison had been "against a third thing
entirely" — Sage has not used LLM for that table since 2015
([qt-kostka.md](../record/qt-kostka.md)); `theta_qt` and `scalar_qt` were
"identified numerically as different operators" before Δ shipped
([macdonald-operators.md](../record/macdonald-operators.md)). A capability
claim then cites the measured survey, per [style.md](../style.md) ("Say what
it opens"). Reading source is bounded by license, and the boundary lives in
[NOTICE.md](../../NOTICE.md): public domain may be read and ported, with the
consultation recorded; GPL is used as a black-box oracle only, and where
access happened anyway, the clean room is the remedy
([cleanroom-spec-skew-lr.md](../cleanroom-spec-skew-lr.md)).

### V3 — Independence is the absence of a shared failure

The test of a proposed check: name one bug both sides could contain. If one
exists, the two sides are one check. Both orientations of one `SkewLr` traversal
fail together — "a bug in the shared layer code reproduces itself in both
orientations" ([littlewood-richardson.md](../record/littlewood-richardson.md)).
An oracle route that reuses the code under test checks nothing (zonal
polynomials would have re-used Jack, V-class 3 above). **Sage is such a route
by default now**: with the backend installed its conversions, Hall–Littlewood,
Jack, Macdonald and character bases all dispatch here, so a harness or
generator that calls Sage as the oracle has to run under `SAGE_DISABLE_SYMFN`
— and *refuse* rather than trust the invoker to have set it, as
[gen_sage_oracle.sage](../../scripts/gen_sage_oracle.sage) and
[check_qt_kostka.py](../../scripts/check_qt_kostka.py) do. A run in the other
state reported 9547 comparisons and 0 mismatches and meant nothing
([python-and-sage-interop.md](../record/python-and-sage-interop.md)). Two
sides on one fixed width "wrap identically and agree on the same wrong answer"
([failure.md](failure.md) R10). A law that a consistent normalization error
satisfies pins nothing (`⟨P,Q⟩ = δ`, above). Each added check is chosen for
what the existing ones cannot see, and its doc says so — the house form is
already in the tree: "agreement between them is evidence rather than
tautology" ([charge.rs](../../src/charge.rs), [kf.rs](../../src/kf.rs),
[jack.rs](../../src/jack.rs), [gjmod.rs](../../src/gjmod.rs)).

### V4 — Durable oracle evidence is a committed fixture

The offline form: generating script committed, fixture committed, `cargo
test` needs nothing installed
([gen_sage_oracle.sage](../../scripts/gen_sage_oracle.sage) →
[sage_oracle.rs](../../tests/sage_oracle.rs);
[gen_lrcalc_oracle.py](../../scripts/gen_lrcalc_oracle.py) →
[lrcalc_oracle.rs](../../tests/lrcalc_oracle.rs), with the license terms of
recording a GPL program's output in [NOTICE.md](../../NOTICE.md)). Live
harnesses serve width, in-process comparison, and driving Sage itself
(`sage_backend.py`'s 8647 computations) — but a family whose only oracle
evidence is a live script has no oracle evidence on most days; that is the
state the Jack fixture was added to correct ([jack.md](../record/jack.md),
"Offline oracle fixtures"). Fixture parsers are strict, so generator drift is
a parse error rather than an absorbed wrong value — "a `1/2` token in an
integer fixture is a parse error, not a wrong answer, but only because the
parser was strict" (same entry). Fixture cases are cold and distinct wherever
memoization on either side would turn a repeat into a cache hit
([oracles-and-comparisons.md](../record/oracles-and-comparisons.md)).

### V5 — Evidence reaches where the claims reach

[failure.md](failure.md) R9 requires the wall stated; this rule requires the
territory inside it defended. The lrcalc fixture exists because Sage
"already validates symfn on small degrees" and lrcalc "earns a separate
fixture by reaching much larger shapes"
([lrcalc_oracle.rs](../../tests/lrcalc_oracle.rs)). Past every oracle's wall,
an independent check that scales: the principal-specialization checksum holds
the 5 313 471-term product no external tool finishes
([littlewood-richardson.md](../record/littlewood-richardson.md)); the modular
second engine took the GJ ladder from n = 10 to n = 14
([jack.md](../record/jack.md)). Where no check reaches, the claim shrinks to
what the evidence covers — the honest summary of the plethysm comparison is
"faster than Symmetrica on the inputs Symmetrica accepts"
([oracles-and-comparisons.md](../record/oracles-and-comparisons.md)).

### V6 — Evidence starts where the domain starts

V5 defends the top of the advertised range; this rule defends the bottom.
[style.md](../style.md) classifies the degenerate inputs — ∅, degree 0,
`k = 1`, equal shapes in a skew pair — as convention choices, not corner
cases, and its own distinction applies to them: a stated convention with no
pin is a promise. So the edges carry pins like any convention —
`degenerate_shapes` in [skew_lr.rs](../../src/skew_lr.rs),
`empty_shape_gives_one_tableau_and_mismatched_sizes_give_none` in
[kostka.rs](../../src/kostka.rs) (one function pinning the convention and the
refusal), "everything contains ∅" in
[partition.rs](../../src/partition.rs) — and two obligations go past the pin:

- **Sweeps include their floor.** Every evidence class quantifies over a
  range, and the degenerate cases sit at its bottom — where recursion base
  cases and convention forks live, so an agreement sweep from `n = 1` checks
  neither engine's base case. The law suite is the model: `(0..=5)` in
  [algebra_laws.rs](../../tests/algebra_laws.rs) puts ∅ and degree 0 inside
  every law rather than beside them. A floor above zero is a decision, stated
  at the site — not a default inherited from an example.
- **Every edge is classified.** A degenerate value that is a **theorem** is
  an answer, and stays pinned so a later validation pass cannot quietly eat
  it — `check_theorem_zeros_still_answer` in
  [check_python_boundary.py](../../scripts/check_python_boundary.py) pins
  six. A **convention** over a well-posed question is pinned and swept (the
  kostka test above). A question the family leaves undefined **refuses** —
  the five conventional zeros that became errors — and the refusal is
  itself pinned ([failure.md](failure.md) R2 and R11; the `None` half of
  the same kostka test).

### V7 — A verifier must be shown able to fail

A lossy check — a checksum, an invariant, a law — ships with its negative
control, committed beside it. The LR specialization checksum detects 412/412
single-coefficient perturbations, and the record says why that number exists:
"a checksum that silently always passed would be worse than no check, so that
number is the one that makes the PASS meaningful"
([littlewood-richardson.md](../record/littlewood-richardson.md)). The `ω_α`
duality test "ships with its own negative control, because the α-twist is
silent when wrong" ([jack.md](../record/jack.md)). The control is part of the
check, not a one-time experiment.

### V8 — The slow engine stays, and the fast ones answer to it

The standing policy, stated once, here. Every family keeps its most
obviously-correct implementation as the in-house oracle — `NaiveLr`
([lr.rs](../../src/lr.rs)), charge ([charge.rs](../../src/charge.rs)),
Schubert E1 "kept forever as the oracle"
([schubert.md](../record/schubert.md)), the general border-strip form behind
the masked one
([oracles-and-comparisons.md](../record/oracles-and-comparisons.md)) — and
every faster engine is held to exhaustive agreement with it on everything it
can finish, as a `cargo test`, before any benchmark ("verified against
`NaiveLr` on every product with |μ|+|ν| ≤ 7",
[the record index](../record/README.md)). An engine earns its keep by
independence and range, not speed: the operator (q,t)-Kostka route stays
"twice over, since it shares no mathematics with either alternative"
([qt-kostka.md](../record/qt-kostka.md)); an engine that duplicates a kept
one's mathematics and extends no wall is a maintenance cost, not evidence.
A benchmark that times several routes also asserts their agreement at every
degree it times ([bench_qtk_routes.rs](../../examples/bench_qtk_routes.rs)).

### V9 — Theorems are enforced, conjectures observed

Every property the output is known to satisfy is classified and treated by
class. A **theorem** is an assertion: ℚ[b]-polynomiality (Dołęga–Féray) and
integrality (Ben Dali) of the GJ coefficients "are enforced". An **open
question** is observed and reported, never corrected: "positivity is open for
both and is only observed, with any negative coefficient reported as a
finding rather than debugged away" ([jack.md](../record/jack.md)) — the same
rule as the valley side of the Delta conjecture, whose driver states which
disagreement is a bug and which is a discovery
([dyck-paths.md](../record/dyck-paths.md); [style.md](../style.md),
"Research drivers"). A **convention** is a pin (class 5). The
misclassification to fear is enforcing a conjecture: it converts the one
output that would matter into a crash — or worse, into a "fix".

### V10 — Coverage is a number with a bound

"Verified" is a quantified claim: 434/434 Kronecker products, 8647
backend-driven computations, 32 448 skewing checks, 48/48 plethysm cases
across both input shapes, exhaustive through `|μ| + |ν| ≤ 7` — the record's
existing forms. The bound is part of the claim, and so is what lies outside
it: `matchings_jack_coverage` reports how much of its own output is not
already a theorem — 83.5% open at n = 8 — and the same file states the bar a
wider sweep must clear before it is worth reporting ("n ≥ 25",
[jack.md](../record/jack.md)). A sweep with an unstated bound reads as
exhaustive; a coverage nobody computed reads as 100%.

### V11 — A timed case is a verified case, against the baseline that matters

A benchmark verifies every case it times, in the same run. The
comparison harnesses do exactly that ("Every case is verified,
not merely timed",
[oracles-and-comparisons.md](../record/oracles-and-comparisons.md)) — which
is how the Sage ladder cross-checks every rewrite for free. The baseline is
like-for-like and labeled (the `via C/py` column; the plethysm 9x that was
0.14x against the C that Sage never calls), and a benchmark that cannot be
like-for-like is deleted rather than caveated — "a benchmark whose caveat is
'this number is not the comparison you want' can only mislead"
([qt-kostka.md](../record/qt-kostka.md)). The numbers, methodology, and
caveats land in the record, which owns them ([style.md](../style.md), "Range
and performance in rustdoc").

## Choosing the evidence

### Two questions before any code

**Who else computes this object?** The answer sets the top of the ladder:
an oracle over the whole intended range, an oracle that walls early, a
near-miss (a similarly-named object under another convention — establish the
difference numerically, V2), or nothing. Everything below the top class is
then chosen to cover what it cannot see.

**What would the plausible wrong answer look like here?** Wrong-by-a-twist —
a convention family — demands distinguishing pins and the dictionary
([style.md](../style.md), "Conventions get their own section").
Self-consistently wrong — a route resting on one identity — demands a
definitional check ([kronecker.md](../record/kronecker.md)). Right at small
degree, wrong at range — sharing, overflow, and DP bugs with onset — demands
V5. Faithfully wrong — a transcribed misprint — demands V1. The battery is
chosen against the failure the family is actually exposed to.

### The table

| the new thing | minimum evidence before it ships | in-tree model |
|---|---|---|
| family an oracle covers on the whole intended range | offline fixture sweep + convention pin + laws | the classical layer in [sage_oracle.rs](../../tests/sage_oracle.rs) |
| family whose oracle walls below the intended range | fixture in the overlap + an independent scaling check past the wall, with negative control | LR: [lrcalc_oracle.rs](../../tests/lrcalc_oracle.rs) + `verify_specialization` |
| family with only a near-miss oracle | the difference identified numerically + the dictionary + a distinguishing pin | `theta_qt`/`scalar_qt`; the four circulating G̃'s in [llt.rs](../../src/llt.rs) |
| family with no oracle but identities to one | identity suite chosen for what each pins + the anchor object fixtured | Δ, Δ′, Θ tied to ∇; `q = 0` Macdonald → HL |
| family nothing else computes | second engine sharing no mathematics + specialization pins vs independent definitions + theorems enforced | the GJ tables ([jack.md](../record/jack.md)) |
| new engine for an existing family | exhaustive `cargo test` agreement with the in-house oracle + external spot checks at scale | `SkewLr` vs `NaiveLr` + lrcalc |
| conjecture-checking driver | the theorem side as control + the interpretation contract | [delta_conjecture.rs](../../examples/delta_conjecture.rs) |

### The distinctions that get miscalled

- **Consistency vs independence.** Two runs of one engine — another
  orientation, another degree, another thread count — share every bug. They
  are worth having and worth nothing alone (V3).
- **Wide vs different.** A redundant engine's job is to be a *different*
  algorithm, not a far-reaching one: Knop–Sahi tableaux stop at n = 5 and
  that is accepted, "it is exponential, and its job is to be a different
  algorithm rather than a wide one"; the `b = 1` matchings count is
  `(2n−1)!!` per λ, "a pin, not an engine"
  ([jack.md](../record/jack.md)). Range is V5's job.
- **A law vs a pin.** Round-trips, homomorphisms, and involutions hold under
  every self-consistent convention. Laws catch broken arithmetic; only a
  pinned value catches the wrong normalization shipped fluently (class 4
  vs 5).
- **The oracle vs a namesake.** An incumbent function with the right name may
  compute a different object (`theta_qt`), a different engine than believed
  (Sage's post-2015 (q,t)-Kostka), or a different normalization (four
  circulating G̃'s). Until the difference is measured, agreement and
  disagreement are both uninterpretable (V2).
- **More degrees vs another kind.** Once theorems and existing checks cover a
  region, widening a sweep of it buys little — the GJ file sets its own bar
  at n ≥ 25 for the next sweep to be worth reporting. A new *kind* of check
  widens what is covered; a wider run of an old kind covers nothing new (V10).

### Defaults when unsure

New family → extend
[gen_sage_oracle.sage](../../scripts/gen_sage_oracle.sage) first. Sage cannot
compute it → hunt the identity to something oracled. Nothing computes it →
build the obviously-correct slow engine first and keep it (V8). New engine →
the agreement test lands before the first benchmark (V8). New sweep or
fixture block → floor at 0, and a floor above 0 states its reason at the
site (V6). New formula from a paper → script it against an oracle before
implementing it (V1). Unsure whether two checks are independent → name the
bug both could share (V3). A family that fits no row of the table is a
policy gap: extend this file in the same change, rather than improvising
silently.

## Where the policy surfaces

- the committed fixtures and their generating scripts
  (`tests/fixtures/`, [scripts/README.md](../../scripts/README.md));
- the oracle, law, and engine-agreement tests `cargo test` runs with nothing
  installed ([sage_oracle.rs](../../tests/sage_oracle.rs),
  [lrcalc_oracle.rs](../../tests/lrcalc_oracle.rs),
  [algebra_laws.rs](../../tests/algebra_laws.rs), the per-module agreement
  tests);
- the live harnesses, one role statement each in
  [scripts/README.md](../../scripts/README.md);
- the Sage-free boundary suite
  ([check_python_boundary.py](../../scripts/check_python_boundary.py)),
  where the surface's theorem zeros and refusals are pinned — the surface
  itself is [python.md](python.md)'s;
- convention-pin doctests, per [style.md](../style.md);
- the record: coverage counts, negative-control numbers, and the measured
  surveys behind capability claims.

## What this changes

Current practice built this policy — the rules above quote the tree rather
than aspire for it. The three deltas it opened with are closed; each is kept
here with what closing it cost, because the gate is what a later drift is
caught by:

1. **The five local statements of the standing policy point here — done.**
   [lr.rs](../../src/lr.rs), [charge.rs](../../src/charge.rs),
   [kf.rs](../../src/kf.rs), [jack.rs](../../src/jack.rs) and
   [gjmod.rs](../../src/gjmod.rs) each keep their one-line local conclusion
   and carry a `docs/policies/validation.md` pointer, so the next engine's
   author finds the rule rather than the folklore; `gjmod.rs`'s citation of
   "the `qtkostka.rs` standard" is now a citation of V3. Gate: the five files
   above, re-grepped at edit time.

2. **Every sweep floor is a decision — done.** Of 148 `for n in 1..=` sites
   in `src/` and `tests/`, 137 lowered to 0 with the suite still green: the
   ∅ case was simply never being swept. The 11 that did not are the ones
   whose object has no degree-0 case, and each now says so in one line at the
   site — the Adams operations (`ψ⁰` is not multiplicative), `f_(n)`'s sign
   `(−1)^{n−1}`, the Delta/shuffle ladder indexed by `k = n−1`, and the [GJ]
   tables, which `gj_connection_tables` returns empty at n = 0. That last
   edge is classified rather than merely skipped, per V6's second bullet:
   `the_tables_are_empty_at_zero` pins it. The generator sweeps from 0
   throughout. Gate: grep `for n in 1..=` in `src/` and `tests/`, and
   `range(1,` in the generator — every surviving hit carries its reason, and
   the one in [skew_lr.rs](../../src/skew_lr.rs) is a run length rather than
   a degree, which it says.

3. **The offline fixture reaches every family — done.**
   [sage_oracle.rs](../../tests/sage_oracle.rs) now covers Hall–Littlewood
   `Q'` and `P`, Kostka–Foulkes, the (q,t)-Kostka table, `H̃`, `∇e_n`,
   Macdonald `P`/`Q`/`J`, the Kronecker product, the three LLT ribbon
   dictionaries and the Schubert structure constants, alongside the classical
   layer and Jack, all checked with nothing installed. Case counts are in each
   family's record file. Cases were chosen to distinguish
   conventions, so the fixture doubles as the pin: `H̃_{(2)}` against
   `H̃_{(11)}` separates `H̃` from a `q ↔ t` transpose, `Q'` against `P`
   separates the two Hall–Littlewoods, and the LLT block writes the grading
   into symfn's `q` slot where Sage names it `t`. Every family's fixture was
   perturbation-tested before it landed — a wrong value in each must fail the
   suite, per V7. Sage's own walls keep the degrees small, which is the
   design: the fixture is the durable floor, the live harnesses keep the
   width, and the in-tree second routes carry the range (V5).

A pass over the table afterwards found three more gaps, all closed:

- **A lossy check whose control nothing ran.** The LR principal-specialization
  checksum is the only independent check past lrcalc's wall — this file's own
  model for the second row of the table — but it and its 412/412 negative
  control lived only in `examples/verify_specialization.rs`. V7 wants the
  control to be part of the check, and an example nothing runs is exactly the
  one-time experiment V7 excludes. Both halves are now
  [lr_specialization.rs](../../tests/lr_specialization.rs).
- **The Hopf structure had laws but no values.** `coproduct`, `antipode` and
  `counit` are a first-row family — Sage computes all three — resting on the
  Hopf axioms and the LR route, and laws are convention-blind. Fixtured, with
  the antipode as the pin. One limit is recorded rather than left unstated:
  Sym is cocommutative, so no value can distinguish the coproduct's tensor
  orientation from its transpose.
- **The specializations likewise.** `principal_specialization`, its graded
  form and `dimension` had strong in-tree checks and a live script but no
  offline oracle; `f^λ` is now held to a direct standard-tableau count rather
  than to Sage's hook formula, which would have re-used symfn's own.

V1 and V2 describe practice from the Hall–Littlewood work onward and bind
prospectively; nothing here is retroactive beyond items 1–3.
