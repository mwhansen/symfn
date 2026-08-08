# Python policy — how this library is consumed

What this file governs: the Python surface — which entry points exist, which
layer of the boundary each belongs to, what shape their data has, and what
the surface owes its two audiences: a bare CPython caller, and Sage through
the adapter. What it does not govern: how a computation may fail
([failure.md](failure.md), whose invariant applies at this boundary in
full), prose style ([style.md](../style.md), which owns docstrings as it
owns every other surface — P11 below states how its rules read at this
one), or measurements — those belong to the record,
and [python-and-sage-interop.md](../record/python-and-sage-interop.md) is
the subsystem's file. Numbers appear below only where a rule leans on one,
and each names its owner.

Like the other rulebooks, this is policy, not a description of current
practice. Where the two disagree, the gap is listed in
[What this changes](#what-this-changes) at the bottom, in execution order.

## The invariant

**The Python surface is the contract this library is consumed through:
plain, ring-free data, crossed whole-object, serving bare CPython and Sage
alike — and Sage is one consumer on the far side of the boundary, never a
dependency.**

It is stated this hard because of who leans on it. The Rust caller is rare
and served for free ([style.md](../style.md), "The readers"); every other
consumer reaches this library through Python, and the heaviest intended one
is Sage itself. Once the adapter — and eventually Sage
([release-readiness.md](../release-readiness.md), Phase 5c) — pins this
surface, a break here is a Sage bug. So the asymmetry is deliberate: the
crate's internals stay free to churn, and the Python surface freezes
hardest of anything in the tree.

The boundary has exactly **three layers**:

1. the **contract layer** — the `#[pyfunction]`s of
   [python.rs](../../src/python.rs): coarse entry points over plain data;
   documented, stubbed, versioned, and what both other layers are built on;
2. the **convenience layer** — pure Python inside the wheel: ergonomics
   defined entirely by contract calls, computing nothing (`Sym`, the basis
   factories, the family namespaces and `Schub`, in `python/symfn/`);
3. the **adapter** — Sage-side, outside the wheel
   ([sage_backend.py](../../scripts/sage_backend.py) and the compiled loop
   [symfn_cy.pyx](../../scripts/symfn_cy.pyx)): everything Sage-shaped
   lives here and nowhere else.

There is no fourth layer: no Sage import inside the wheel, no convenience
class with an algorithm in it, no entry point whose caller nobody can name.

## The rules

### P1 — Elements cross as plain, ring-free data

An element is a list of `(support, coefficient)` pairs keyed by its
combinatorial support — partition tuples for symmetric functions,
permutations for Schuberts. A parameter family crosses as exponent-keyed
rows: `(exponent, coefficient)` for one variable (`t_poly` in
[python.rs](../../src/python.rs)), `(a, b, coefficient)` for `q^a t^b`.

**A support crosses out as a `tuple` and in as any sequence.** Outbound the
type is part of the contract, not an accident of PyO3's default: every
consumer uses a support as a key — Sage's adapter interns partitions on the
tuple of parts and builds `{Partition: coefficient}` dicts, and a bare
CPython caller reaches for a `dict` or a `set` just as fast — and a list has
to be copied into a tuple before either can hold it, once per output term,
in the loop this library's marshalling budget is spent in
([python-and-sage-interop.md](../record/python-and-sage-interop.md)).
Inbound, a caller may hand back what it received or pass the list it already
has. The boundary is strict about what it promises and permissive about what
it accepts; `Key` in [python.rs](../../src/python.rs) is the one type that
holds both halves.
**Plain data** means a caller reads a result with nothing but the standard
library, and every value is exact: `int`s of any size, and a documented
integer encoding wherever a denominator exists. Nothing Sage-shaped,
nothing that requires importing a class to unpack.

Ring-freedom is the deeper half, and it is what serves both audiences at
once: Sage factors the coefficient ring out of its conversion path and asks
only for integer transition rows, recombining over `ℚ[t]` on its own side
([python-and-sage-interop.md](../record/python-and-sage-interop.md)). A
surface that ships ring-independent structure data lets every consumer —
Sage, SymPy, a bare interpreter — rebuild elements in its own ring, and is
the only thing a Sage-free wheel can promise without lying. A new
coefficient kind gets its documented plain-data encoding before any
function ships it.

### P2 — Calls cross whole objects; tables and indices are first-class

Every entry point does a whole-object operation — multiply two complete
elements, convert an entire expansion — never a per-monomial call; the
module doc of [python.rs](../../src/python.rs) opens with this rule, and it
is not relaxed here. Where the natural unit of work is larger than one
value, the larger unit is the entry point: `kostka_foulkes_column` exists
because the Morris recursion produces a column at a time
([hall-littlewood.md](../record/hall-littlewood.md) owns the measurement),
and the `*_table` family exists because consumers ask per degree. The hot
path's currency is indices into `partitions(n)` rather than lists of parts
(`convert_indexed`), because marshalling, not mathematics, is where
end-to-end time went: the first backend measurement was 91% Python glue
([python-and-sage-interop.md](../record/python-and-sage-interop.md)). A
fine-grained accessor is never the fix for a slow caller; a bulk entry
point is.

**A pair of one-way entry points lets the caller build the accessor P2
forbids.** `schur_to_homogeneous` and `power_to_schur` invite the caller to
compose them, and a caller that has composed them has already chosen the route:
`p → h` becomes `p → s → h`, and no direct rule in the kernel can be
reached. So a conversion names its *pair* — `convert_terms(a, src, dst)`,
`convert_indexed`, `to_power(a, src)` — and the routing decision stays on
this side of the boundary, where it can be measured
([transitions.md](../record/transitions.md)).

### P3 — Sage is a consumer, never a dependency

The wheel builds and imports with no Sage anywhere: not in
`[build-system]`, not at import time, not as an optional extra.
[python.rs](../../src/python.rs) names Sage nowhere today
([release-readiness.md](../release-readiness.md), Phase 5) and keeps it
that way deliberately. Everything Sage-shaped — its `Partition` and
`Integer` types, its ZZ-vs-QQ element contract, its orderings, its
exception surface — lives in the adapter. The enforcement is CI's `wheel`
job, which installs the wheel in a bare interpreter, exercises both layers and
fails if any `sage` module reaches `sys.modules` — because the failure mode
this rule exists against is a convenience import creeping in silently.

### P4 — The convenience layer computes nothing

Pure Python, defined entirely in terms of contract calls: it may batch,
memoize, type-tag, pretty-print, and overload operators; it may not
implement mathematics. A convenience that computes is a second engine — the
trap [failure.md](failure.md) names for fallback paths ("the wide pass is
the same code"), one language out: an untested second copy, exercised only by
the callers least equipped to notice a twist. The adapter never imports this
layer; it exists for humans.

**"Computes nothing" is a claim about equality, so it is checked as one.**
`scripts/check_convenience.py` runs each method against the contract sequence
it claims to be, over every partition to degree 6 and every basis. A method
with no such sequence to compare against does not belong in this layer.

Two things the rule does *not* forbid, because both were needed and both stay
inside it. **Routing is allowed**: `h_λ · h_μ` has no product entry point, so
the multiplication converts to Schur, multiplies and converts back — the
decision is made on this side, where P2 says it belongs, and the check above
asserts the route changes no value. **Arithmetic on a returned coefficient is
allowed**: the parameter families cross as exponent-keyed rows (P1), and
substituting numbers into those rows — `Param.at`, and the `at` on each
coefficient type — is arithmetic on plain data, not symmetric-function
mathematics. It can produce no coefficient the contract layer did not. It also
earns its place: `P_λ(x; q, q) = s_λ` and `Q'_λ(x; 1) = h_λ` are how a
convention is checked from Python at all, and without evaluation the families
arrive as rows nothing on this side can test.

The line is that neither may reach a value the contract layer cannot. A
convenience that answers where the contract layer has no answer is a second
engine whatever it is called.

### P5 — Totality lives in the wheel; fidelity lives in the adapter

The wheel answers everything it correctly can, in the literature's
conventions; being more total than an incumbent is capability, the point of
the library. The adapter narrows that to the consumer's exact contract: its
values, its element types, its list orders, its exception surface. The
record's Schubert chapter is the grounding: symfn correctly answers 187
inputs where Symmetrica raises `ValueError`, Sage's doctests assert those
errors, and a faithful drop-in must reproduce them
([python-and-sage-interop.md](../record/python-and-sage-interop.md)) — so
the refusals belong to the adapter, and the answers stay in the wheel.
An adapter entry is designed against what its caller does with the result,
which the entry point's name does not tell you: only driving the consumer —
[check_backend.py](../../scripts/check_backend.py), not a dump — revealed
Sage's two calling conventions.

### P6 — Every list states its order, or states that it has none

A list-returning entry point documents its output order, zero-freeness, and
deduplication, or says explicitly that the order is unspecified. The rule
has teeth because order becomes contract the moment a consumer prints it:
`semistandard_tableaux` must return increasing lex in the row-major reading
word because Sage's doctests print the list
([python-and-sage-interop.md](../record/python-and-sage-interop.md)). Leave
the order undocumented and it is still a promise, made without anyone
deciding to make it, and found when changing it breaks a consumer.

### P7 — Conventions are pinned where Python executes them

The conventions in circulation cross the boundary intact, and Python is where
it bites hardest: no type distinguishes a basis, so a wrong convention is a
plausible wrong answer, not an error. Every family's Python-facing doc
states its convention and names the Sage equivalent or its absence
([style.md](../style.md), "Conventions get their own section"), and the pin
is executable on the Python side: a value in the Sage-free boundary suite
(and, once the convenience layer exists, its doctests) chosen to
distinguish the shipped normalization from its rivals. The convenience
layer also carries basis identity — misusing a basis raises; it never
computes garbage — which is the crate's "basis confusion is a compile
error" ([README.md](../../README.md), the design table) translated into a
language with no compiler.

### P8 — Exact or loud crosses the boundary intact

[failure.md](failure.md)'s invariant applies unchanged; this surface adds
its readings. A refusal is a typed Python exception naming the violated
requirement — `perm_arg` in [python.rs](../../src/python.rs) is the model —
and a `PanicException` reaching Python from input a caller can type is a
bug by definition ([failure.md](failure.md), R2), never an interface. No
ceiling is visible from Python: coefficients cross as arbitrary-size `int`
in both directions, and element-valued entry points run fixed-width and
re-run over `BigInt` on overflow — at 0–1% on the fast path
(`examples/bench_guarded.rs`), allocation-free in the common case (the
`Coeff` enum; [python-and-sage-interop.md](../record/python-and-sage-interop.md)
owns both measurements). The shipping bar that follows: a family is
supported only once crossing its wall escalates or raises a typed
exception. The `(q,t)` and `α` families meet neither today; delta 4 below
states the bar, and [failure.md](failure.md) owns the execution.

### P9 — One wheel, one behavior

There is exactly one behavior a user can install: the `python` feature
always carries the escalation ring ([Cargo.toml](../../Cargo.toml)), the
build is `abi3` so one artifact serves CPython 3.9+, and no build variant,
platform, or support tier may change an answer. The tier policy
([release-readiness.md](../release-readiness.md), Phase 5) decides where a
prebuilt wheel exists — availability, never semantics. The wheel has no
Python-level dependencies, matching the crate's zero-dependency default.

### P10 — The supported surface is a deliberate list

Membership is decided, not accumulated: **108 entry points** exist, and each is
either **supported** — stubbed in `symfn.pyi`, documented to
[style.md](../style.md)'s checklist, held stable — or **harness-only** —
underscore-prefixed, absent from the stubs, free to change, kept for
`scripts/check_*.py`.

⚠️ **Count it from the `#[pymodule]` block, never by grepping
`#[pyfunction]`.** This file said 91 and
[release-readiness.md](../release-readiness.md) has said 87 and 92; all three
were attribute greps, and the attribute undercounts twice over — two of them
sit inside `out_of_schur!` and `into_schur!` and expand to nine conversion
entry points between them, and one apparent match is the string
`#[pyfunction]` inside a doc comment. `dir(symfn)` and the registration block
agree at 108.

The sort ran at 98, and at 108 after the Cython branch merged; **the
harness-only set came out empty either way**, which is a
result rather than a deferral: every entry point does a whole-object operation
with a named audience, because P2 has been enforced since the file was
written, so the fine-grained probes this category exists to absorb never
accumulated. The two that looked like probes are not — `schubert_monomial_mass`
documents a caller who wants the out-of-family flag before committing, and
`clear_caches` is what any consumer timing a run needs. If the set is still
empty at the next addition, that is P2 working, not the sort being skipped.

**A convenience name may never take a contract name.** The supported names sit
flat at `symfn.*`, and `symfn/__init__.py` re-exports the contract layer and
then the convenience layer, so a collision does not raise — it replaces. The
Hall-Littlewood namespace was `hall_littlewood` until this was noticed, which
removed the entry point of that name from the surface entirely; it is `hl`, and
`scripts/check_convenience.py` asserts every contract name still resolves to
the contract object. When the readable name is taken, the convenience one
yields, because only one of the two is frozen.

Low-level is not a third category:
the indexed and bulk entry points are supported *and* documented as
low-level, because they are precisely what the adapter — and, at Phase 5c,
Sage — pins. The supported names live flat at `symfn.*` and survive any
package layout change; `__version__` is sourced from the crate version so
the two cannot drift.

### P11 — Docstrings are rustdoc in translation

[style.md](../style.md)'s rustdoc rules govern the Python docstrings — for
the contract layer literally, since the `///` on a `#[pyfunction]` is the
`__doc__` PyO3 ships. This rule states the translation: the text is written
for the Python reader, who meets it as `help(...)` and a stub tooltip,
never as docs.rs — while remaining valid rustdoc.

- **The checklist crosses intact**: a first sentence that stands alone; the
  contract — what comes back and in what state (P6's order statement
  included), argument requirements, degenerate inputs, cost in shape terms;
  the convention stated, with the Sage equivalent named or its absence
  stated; backticked Unicode math; no roadmaps, no history, no seconds, no
  `×`-ratios.
- **`# Panics` becomes `Raises`.** The contract documents typed exceptions
  (P8); a panic is never documented as interface here, because it is not
  one.
- **Examples are Python doctests** — `>>> symfn...` blocks with
  hand-checkable, convention-distinguishing values — because this reader
  pastes Python, not Rust. In `///` they sit in ` ```text ` fences (the
  house form, already in use in [python.rs](../../src/python.rs)) so
  cargo's doctest runner does not compile them as Rust; the runner that
  executes them is `scripts/check_python_docs.py`; the convenience layer's
  are ordinary doctests, run by `scripts/check_convenience_docs.py`. An
  example no runner executes is not a pin: nothing fails when it stops being
  true.
- **Pointers are backticked repo paths**, which read identically in all
  three renderings — docs.rs, `help()`, the stubs — where an intra-doc
  link resolves only in the first.
- **The module docstring owns the model**, as in rustdoc: the
  `#[pymodule]`'s doc carries the data representation, the escalation
  contract, and the pointer to this file; `symfn/__init__.py`'s docstring is
  the package's front page and names both layers.
  Each family's entry point states its own convention and may delegate
  depth, never the convention itself, to the Rust module doc it names. In
  the convenience layer the division repeats one level down: a class
  docstring owns the type's invariant and representation; its methods own
  their contracts. No framework markup anywhere — plain sections and
  standard doctest format, readable as text with no Sphinx to render it.

## Choosing a home

Two questions before any addition:

**Who calls it?** A human at a REPL wants the convenience layer — over an
existing contract call, which is added first if missing. The adapter's
per-term loop wants an indexed or bulk contract entry. A check script alone
wants a harness-only function. If no caller can be named, it is not added.

**Does it cross once?** If callers would loop it, the loop is the entry
point: expose the column, the table, or the whole-object form instead.

| the new thing | where it lands | model |
|---|---|---|
| a new computation | contract layer: whole-object, plain data, escalating | `schur_multiply`, `macdonald_p` |
| ergonomics over an existing computation | convenience layer, computing nothing | `Sym.omega`, `Sym.to`, `macdonald.P` |
| anything Sage-shaped: types, orders, exceptions | the adapter | the ZZ/QQ element rule in [sage_backend.py](../../scripts/sage_backend.py) |
| a hot-loop marshalling win | an indexed/bulk contract entry; optionally the compiled shim, adapter-side | `convert_indexed`; [symfn_cy.pyx](../../scripts/symfn_cy.pyx) |
| a probe only a check script calls | harness-only: underscore-prefixed, no stub | the set came out empty; P10 records why |
| a new coefficient kind | a documented plain-data encoding, before any function ships it | `t_poly` rows; the `(a, b, coefficient)` triples |

### The distinctions that get miscalled

- **Low-level is not unsupported.** The indexed entry points are among the
  most heavily used functions in the file — hiding them would unpin the
  exact surface Sage is meant to pin. They are supported, and documented as
  low-level.
- **Convenience is not a second engine.** The moment a convenience method
  computes something the contract layer does not, it stops being
  convenience; the computation moves down and the method becomes a call.
- **Total is not faithful.** The wheel answering where the incumbent errors
  is capability; the adapter doing the same is a broken drop-in. Which one
  is being written decides which behavior is the bug (P5).
- **A documented wall is not a visible ceiling.** Rustdoc stating a
  fixed-width wall makes the *crate* honest ([failure.md](failure.md), R9);
  the *wheel* additionally may not list the family as supported until
  crossing the wall is a typed exception or an escalation (P8). The two
  bars are different, and the second one is this file's.

### Defaults when unsure

A new entry point is contract-layer: whole-object, plain data, escalating,
order stated, stubbed, convention pinned. When the caller is not yet known,
add the contract call and wait for the caller before adding convenience —
only the caller reveals what it must return. Anything that would import
Sage goes to the adapter, whatever else it is.

## Where the policy surfaces

- the module doc of [python.rs](../../src/python.rs) opens with the
  coarse-grained rule and points here;
- `python/symfn/symfn.pyi` is the supported list made machine-readable, held
  to the module by `scripts/check_python_stubs.py`;
- `python/symfn/` is the convenience layer, held to the contract layer by
  `scripts/check_convenience.py` and to P11 by
  `scripts/check_convenience_docs.py`;
- `docsite/` renders both layers from the objects themselves, so `help()` and
  the website cannot drift, and `scripts/check_docs_complete.py` fails when a
  supported name reaches no page;
- `scripts/preflight_python.sh` runs all of the above, and CI runs it;
- each supported function's docstring carries contract, convention, and the
  Sage equivalent per [style.md](../style.md)'s checklist, read as P11
  translates it;
- [release-readiness.md](../release-readiness.md) Phases 2 and 5 keep the
  execution checklists and name this file as their bar.

## What this changes

Current practice already embodies most of this policy — the plain-data
encodings, the whole-object rule, the escalation ladder, and the Sage-free
`python.rs` all exist and are kept as-is. The deltas, in execution order;
each names its gate:

1. **Sort the entry points (P10, P6, P11).** *Done.* This item said "the 91"
   until P10 above established that
   every count taken by grepping `#[pyfunction]` was wrong. All **108** are
   supported and none is harness-only (P10 records
   why), `symfn.pyi` exists and is held to the module by
   `scripts/check_python_stubs.py` — the gate this item named — the
   `#[pymodule]` carries the docstring it lacked, and `__version__` comes from
   `CARGO_PKG_VERSION` so the two cannot drift. No name needed an underscore,
   so `scripts/check_*.py` needed no edit.

   **The sort found a defect first.** 27 of the 98 named
   their first argument `lambda`, a Python keyword, so those calls could not
   use keyword arguments at all and no stub could be written for them — `def
   jack_p(lambda: list[int])` is as much a `SyntaxError` as the call was. All
   27 are now `la`, matching the μ→mu, ν→nu transliteration already in the
   file. Writing the stubs is what forced it into the open, which is the
   argument for the stub file being a gate rather than a courtesy.

   **The docstring half is now done too.** Every one of the 108 carries a
   `# Raises` section naming the exception and the requirement, and an
   executed example: 132 of them, run by `scripts/check_python_docs.py`
   against the built module, which fails both on a wrong value and on an
   entry point that has none. P6's order statement is stated once in the
   `#[pymodule]` doc — increasing lexicographic by support, zero-free,
   deduplicated — and repeated at an entry point only where it differs, which
   three do. The first run of the runner found the module docstring's own
   example printing lists where the boundary returns tuples;
   [python-and-sage-interop.md](../record/python-and-sage-interop.md) has
   what else the sweep turned up.
2. **Typed exceptions at the boundary (P8).** *Done.* Every precondition a
   Python caller can violate is checked at the boundary and raised as a typed
   exception naming the requirement, pinned by
   `scripts/check_python_boundary.py` — 116 malformed calls over 85
   `#[pyfunction]`s, plus a completeness check that fails when a new one
   appears uncovered.

   **A correction to what this item said.** Partition arguments did not meet a
   `PanicException`: `part()` called `Partition::new`, which *normalizes* —
   drops zeros and sorts — so `[1, 3]` reached the mathematics as `[3, 1]` and
   the caller got a well-formed answer to a question they had not asked, across
   roughly 50 entry points. The defect was real and worse than a crash, and
   auditing by "does it panic?" would have missed it entirely, because it
   arrived through the convenience constructor. Five clusters *did* panic —
   non-permutation Schubert term lists, indices past `MAX_SUPPORT`, `k = 0`,
   inhomogeneous Macdonald operator arguments, and the LLT capacity walls —
   and all now
   raise. [failure.md](failure.md) R11 carries the mechanism half, including
   the distinction between a zero that is a theorem and a zero that was a
   convention over an undefined question;
   [python-and-sage-interop.md](../record/python-and-sage-interop.md) has the
   measurements and the dead ends.
3. **The Sage-free gates (P3, P7, P11).** *All but the round-trip half is
   done, and all of it is in CI.* `scripts/check_python_docs.py` executes
   every example on the supported surface against the built module;
   `scripts/check_convenience_docs.py` does the same one layer up;
   `scripts/check_convenience.py` carries the convention pins as
   *degenerations* — `P_λ(x;q,q) = s_λ`, `P_λ(x;α=1) = s_λ`, `Q'_λ(x;0) =
   s_λ`, `Q'_λ(x;1) = h_λ`, `K_{λμ}(1) = K_{λμ}` — each of which fails under
   the twist its family's rivals use, which is a stronger pin than a committed
   value because it is a theorem rather than a transcript. The
   bare-interpreter assertion is CI's `wheel` job: it installs the wheel,
   computes, and fails if any `sage` module is on `sys.modules`.
   `scripts/preflight_python.sh` runs the set.

   **What remains** is the round-trip half — values computed in Rust and
   asserted from Python. Everything above compares the two Python layers to
   each other, which catches a convenience defect and would not catch a kernel
   one; the fixtures under `tests/fixtures/` are the Rust side's answer and
   nothing on this side reads them yet.
4. **The parameter families reach the bar (P8).** The `(q,t)` and `α`
   entry points stop being able to wrap in release; execution is owned by
   [failure.md](failure.md), "What this changes". This file adds the
   supported-list consequence: those families are listed as supported only
   once their walls raise typed exceptions or escalate — a
   `PanicException` from the profile backstop is interim honesty, not an
   interface.
5. **The convenience layer (P4, P7).** *Done.* Pure Python in the wheel:
   `Sym` with a basis tag, the factories `s`, `h`, `e`, `p`, `m`, `f`, the
   family namespaces `macdonald`, `jack`, `hl`, `llt`, and `Schub` over
   permutations. Coefficients are `int` and `Fraction`; the `repr` is
   readable and says so where it is not `eval`-able. It is deliberately not a
   coercion framework and not a `SymmetricFunctions` re-creation — mixing
   bases raises and names `.to()`, which is P7 executed rather than
   described; users who want the environment have Sage. The gate is
   `scripts/check_convenience.py`, 2177 checks, and
   `scripts/check_convenience_docs.py` for the doctests.

   **The layer needed a packaging change to exist**, which is why it landed
   with Phase 5's `pyproject.toml`: a wheel that is only a compiled module has
   nowhere to put Python. The layout is `python-source`, the compiled half is
   `symfn.symfn`, and the supported names stay flat at `symfn.*`.

   **Two rules came out of building it**, both now above rather than here: P4
   says what routing and coefficient evaluation are allowed to be, because the
   layer needed both and neither is a second engine; P10 says a convenience
   name may never take a contract name, because one did.
6. **The adapter becomes installable (P3, P5).** `install()` as a real
   entry point instead of import-time patching, no `sys.path.insert`,
   locating `symfn` as an ordinary installed package. Owned by Phase 5b of
   [release-readiness.md](../release-readiness.md); named here because P3
   and P5 are its bar.
