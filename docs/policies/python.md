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
   defined entirely by contract calls, computing nothing (pending — the
   layer is a delta below);
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

### P3 — Sage is a consumer, never a dependency

The wheel builds and imports with no Sage anywhere: not in
`[build-system]`, not at import time, not as an optional extra.
[python.rs](../../src/python.rs) names Sage nowhere today
([release-readiness.md](../release-readiness.md), Phase 5) and keeps it
that way deliberately. Everything Sage-shaped — its `Partition` and
`Integer` types, its ZZ-vs-QQ element contract, its orderings, its
exception surface — lives in the adapter. The enforcement is a CI job that
imports the wheel in a bare interpreter and exercises the supported surface
(delta 3 below), because the failure mode this rule exists against is a
convenience import creeping in silently.

### P4 — The convenience layer computes nothing

Pure Python, defined entirely in terms of contract calls: it may batch,
memoize, type-tag, pretty-print, and overload operators; it may not
implement mathematics. A convenience that computes is a second engine — the
trap [failure.md](failure.md) names for fallback paths ("the wide pass is
the same code"), one language out: an untested twin exercised only by the
callers least equipped to notice a twist. The adapter never imports this
layer; it exists for humans.

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
Adapter entries are designed against the caller, not the name: "an entry
point's name tells you what it computes; only its caller tells you what it
must return" (same file), and only the consumer driving —
[check_backend.py](../../scripts/check_backend.py), not a dump — reveals
contracts like Sage's two calling conventions.

### P6 — Every list states its order, or states that it has none

A list-returning entry point documents its output order, zero-freeness, and
deduplication, or says explicitly that the order is unspecified. The rule
has teeth because order becomes contract the moment a consumer prints it:
`semistandard_tableaux` must return increasing lex in the row-major reading
word because Sage's doctests print the list
([python-and-sage-interop.md](../record/python-and-sage-interop.md)). An
undocumented order is a promise made by accident and discovered by a
breakage.

### P7 — Conventions are pinned where Python executes them

The convention minefield crosses the boundary intact, and Python is where
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

Membership is decided, not accumulated: 91 `#[pyfunction]`s exist at last
count (re-grep at sort time), and each is either **supported** — stubbed in
`symfn.pyi`, documented to [style.md](../style.md)'s checklist, held stable
— or **harness-only** — underscore-prefixed, absent from the stubs, free to
change, kept for `scripts/check_*.py`. Low-level is not a third category:
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
  executes them is the Sage-free boundary suite (delta 3 below), and
  `doctest` itself for the convenience layer. An example no runner executes
  is prose wearing a doctest's costume.
- **Pointers are backticked repo paths**, which read identically in all
  three renderings — docs.rs, `help()`, the stubs — where an intra-doc
  link resolves only in the first.
- **The module docstring owns the model**, as in rustdoc: the
  `#[pymodule]`'s doc — empty today; delta 1 below — carries the data
  representation, the escalation contract, and the pointer to this file.
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
| ergonomics over an existing computation | convenience layer, computing nothing | pending (delta 5 below) |
| anything Sage-shaped: types, orders, exceptions | the adapter | the ZZ/QQ element rule in [sage_backend.py](../../scripts/sage_backend.py) |
| a hot-loop marshalling win | an indexed/bulk contract entry; optionally the compiled shim, adapter-side | `convert_indexed`; [symfn_cy.pyx](../../scripts/symfn_cy.pyx) |
| a probe only a check script calls | harness-only: underscore-prefixed, no stub | delta 1 below establishes the set |
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
- `symfn.pyi` is the supported list made machine-readable (delta 1 below);
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

1. **Sort the 91 (P10, P6, P11).** Every `#[pyfunction]` becomes supported
   or harness-only; harness-only names gain the underscore and leave the
   stubs; each supported function's doc is brought to
   [style.md](../style.md)'s checklist as P11 reads it — `Raises`, Python
   doctest examples, P6's order statement; the module docstring lands on
   the `#[pymodule]`, which carries none today; and `__version__` lands,
   sourced from the crate version. Absorbs the Phase 2 Python item in
   [release-readiness.md](../release-readiness.md). Gate: `symfn.pyi`
   exists and lists exactly the supported set, and `scripts/check_*.py`
   still pass using the renamed probes.
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
   auditing by "does it panic?" would have missed it entirely; the convenience
   constructor was the leak. Five clusters *did* panic — non-permutation
   Schubert term lists, indices past `MAX_SUPPORT`, `k = 0`, inhomogeneous
   Macdonald operator arguments, and the LLT capacity walls — and all now
   raise. [failure.md](failure.md) R11 carries the mechanism half, including
   the distinction between a zero that is a theorem and a zero that was a
   convention over an undefined question;
   [python-and-sage-interop.md](../record/python-and-sage-interop.md) has the
   measurements and the dead ends.
3. **The Sage-free gates (P3, P7, P11).** A boundary suite that runs on a
   stock runner — [check_bindings.py](../../scripts/check_bindings.py)'s
   job with committed values in place of the Sage oracle — carrying the
   convention pins and executing the supported surface's docstring
   examples, plus the bare-interpreter import assertion, both in CI.
   Absorbs two Phase 5 items. Gate: CI is red if the wheel imports Sage, a
   pin moves, or a docstring example fails.
4. **The parameter families reach the bar (P8).** The `(q,t)` and `α`
   entry points stop being able to wrap in release; execution is owned by
   [failure.md](failure.md), "What this changes". This file adds the
   supported-list consequence: those families are listed as supported only
   once their walls raise typed exceptions or escalate — a
   `PanicException` from the profile backstop is interim honesty, not an
   interface.
5. **The convenience layer (P4, P7).** Pure Python in the wheel:
   basis-tagged elements, operator overloading, exact `int`/`Fraction`
   coefficients, a readable repr — and deliberately not a coercion
   framework or a `SymmetricFunctions` re-creation; users who want the
   environment have Sage. Gate: every method is a composition of contract
   calls, equality-tested against them, with doctests that run Sage-free.
6. **The adapter becomes installable (P3, P5).** `install()` as a real
   entry point instead of import-time patching, no `sys.path.insert`,
   locating `symfn` as an ordinary installed package. Owned by Phase 5b of
   [release-readiness.md](../release-readiness.md); named here because P3
   and P5 are its bar.
