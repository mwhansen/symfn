# CLAUDE.md — the session bootstrap

This file is loaded at the start of every session. Most of this library was
written by coding agents and is maintained the same way: every session starts
with no memory of the last one, so the tree is the interface between sessions
([docs/style.md](docs/style.md), "The readers"). The rules below are the ones
that bite before you would think to look them up; everything else this file
does is route.

symfn is a zero-dependency Rust kernel for symmetric functions — the six
classical bases through Hall–Littlewood, Macdonald, LLT, Jack and Schubert
polynomials — that also builds as a PyO3 wheel and can stand in for
Symmetrica underneath Sage. The [README](README.md) states what the library
opens, and [docs/layout.md](docs/layout.md) lists every module in the tree;
[docs/sage-backend.md](docs/sage-backend.md) is the one account of the Sage
side — what is covered, where the adapter lives, how to turn it off — and
nothing else in the tree restates it.

## The five rulebooks

1. **[docs/style.md](docs/style.md)** governs every prose surface: rustdoc,
   `//` comments, test names, commit messages, docs/. Its closing checklist
   is the bar a new `pub` item must meet. The genre rule matters most:
   rustdoc describes the present — no roadmaps, no history, no seconds, no
   ×-Sage ratios; narrative and measured numbers live in the record.
2. **[docs/policies/failure.md](docs/policies/failure.md)** governs every way
   a computation may fail. The invariant: **every value that leaves this
   library is exact, or the call fails loudly — in every build profile.**
   Overflow escalates or refuses; it never wraps. New code picks its
   mechanism from that file's table; a site that fits no row is a policy gap
   — extend the policy in the same change, not silently.
3. **[docs/policies/python.md](docs/policies/python.md)** governs the Python
   surface, the contract nearly every consumer builds on: plain ring-free
   data crossed whole-object, three layers (contract, convenience, adapter),
   and Sage a consumer on the far side of the boundary, never a dependency.
   Before adding or changing a `#[pyfunction]` — or anything in
   `python/symfn/` — pick its home from that file's table. Two traps that
   only bite in this layer: a convenience name may never take a contract
   name (P10 — one did, and it deleted the entry point), and a wrapper must
   carry the basis its entry point actually returns, which is not always the
   one the family is usually written in (the LLT entry points return
   monomial).
4. **[docs/policies/validation.md](docs/policies/validation.md)** governs the
   evidence a capability must carry. The invariant: **every value a public
   family can produce is covered by at least one check that does not share
   its mathematics, over the advertised range** — and the less the outside
   world can check, the more the tree carries itself. New families pick
   their battery from that file's table; the slow obviously-correct engine
   is kept as the in-house oracle; oracle agreement lives in committed
   fixtures, not in scripts someone must remember to run.
5. **[docs/record/](docs/record/)** is the long-term memory: one file per
   subsystem recording what was built, what was measured, and what failed.
   **Before working in a subsystem, read its record file** — dead ends are
   recorded with their premises exactly so they are not re-explored at full
   price. Every measurement lands there with its harness named; negative
   results are first-class. Future work goes in a record file's open tail or
   in [docs/todo-1.0.md](docs/todo-1.0.md), which holds what is left after
   0.9.0 — nowhere else.

## Commands

The Rust side needs no dependencies, no network, no Sage — the oracle tests
read committed fixtures under `tests/fixtures/`:

    cargo test                        # default suite; runs in seconds
    cargo test --features bignum      # + arbitrary-precision coefficients
    scripts/preflight.sh              # the gate: fmt check + both suites
    cargo fmt --all                   # the pre-commit hook checks, never fixes
    cargo doc --no-deps               # render the reference

The Python surface has its own gate, because everything in it needs
`cargo build --features python` and the last two steps need ruff and Sphinx —
which is why it is not inside `preflight.sh`. Run it when anything under
`src/python.rs`, `python/symfn/` or `docsite/` changes — or the README's
Python example, which runs in its page-examples step:

    scripts/preflight_python.sh       # stubs, typed exceptions, both layers'
                                      # doctests, the convenience layer against
                                      # the contract layer, ruff, mypy --strict,
                                      # docs coverage

What does need externals:

- **Sage** — nearly every script in `scripts/` imports it. The verification
  pattern is dump-then-check:

      cargo run --release --example jack_dump -- 7 5 > /tmp/jack.txt
      sage -python scripts/check_jack.py /tmp/jack.txt

- **Regenerating a fixture** needs the oracle it came from
  (`scripts/gen_sage_oracle.sage`, `scripts/gen_lrcalc_oracle.py`) — and the
  Sage one needs `SAGE_DISABLE_SYMFN=1`, or Sage answers out of this library
  and the fixture is symfn quoting itself. It refuses rather than trusting you.
- **The wheel** needs maturin; the build lines are in the README.

Once per clone:

    git config core.hooksPath .githooks
    git config blame.ignoreRevsFile .git-blame-ignore-revs

## Conventions that bite

- **Normalization traps are the house failure mode.** Every `K_{λμ}` variant,
  LLT `G`, and `H̃` in circulation differs by a twist that yields plausible
  wrong answers, not errors — and an agent fluent in all of them will supply
  the wrong one confidently. Before touching a family, read the module doc's
  convention section ("The conventions in circulation" in [src/llt.rs](src/llt.rs)
  is the model); pin any convention you add with a doctest whose value
  distinguishes it from its rivals, not one every convention agrees on.
- **No TODO, FIXME, or commented-out code.** The tree has zero; keep it at
  zero. Future work goes to the record's open tail, rejected code to git.
- **Test names are propositions** (`to_schur_is_a_ring_homomorphism`), and
  every assertion inside a sweep names its counterexample input.
- **A `//` comment must state** an invariant, the mathematical
  fact that licenses the step, a measured reason (naming its harness), or a
  trap. Never narrate mechanics, never address the reviewer.
- **Pointers are greppable file paths** — `docs/record/llt.md` — never
  section numbers, which nothing checks and which have already drifted here
  once ([docs/style.md](docs/style.md), "Specs, and how they end").
- **Bracketed math goes inside backticks in rustdoc** — a bare `ℚ[q,t]`
  parses as a markdown link and becomes a warning. Bare `Σ_ν`, `μ ⊢ n`, `λ'`
  are safe and stay unmarked, and `//` comments need nothing at all. Citation
  keys resolve from a module-doc `## References` block and are escaped
  (`\[HHL\]`) in item docs, where definitions do not reach. `cargo doc
  --no-deps --all-features` is silent; keep it that way. ASCII in identifiers
  (`lambda`, never `λ`). American English. Prose wraps at 80 columns.

## Commit messages

Plain language, for someone reading `git log` cold — a report, not a record
entry; the record is where the narrative goes. The title says what changed,
with the number when there is one, and a failed experiment gets a title of
the same form. The body says what was wrong before, what changed (naming the
functions), how it was measured and what the numbers are, what was tried and
rejected, what is left open and where it is recorded, and what the tests
pin. No narrator's voice, no coined words, no metaphor. A correction to an
earlier claim gets its own paragraph, never a silent fix. The rules and the
model commits are in [docs/style.md](docs/style.md), "Commit messages".
Before committing, run `scripts/preflight.sh`; the pre-commit hook re-checks
only formatting.

## Measurement discipline

- **On AC power, and record the power state.** Battery drifts ~1.8× and
  changes ratios, not just times (scripts/README.md;
  [docs/record/jack.md](docs/record/jack.md) straddled an AC detach mid-run).
- **Sage memoizes** — time distinct, cold inputs. **symfn memoizes too** —
  `clear_caches()` between timed runs ([src/memo.rs](src/memo.rs)).
- A ×-Sage figure states what Sage dispatched to — its own Python or
  Symmetrica's C — because the two differ by an order of magnitude
  ([oracles-and-comparisons.md](docs/record/oracles-and-comparisons.md)).
- **The Sage A/B needs `SAGE_DISABLE_SYMFN=1` in the control arm's
  environment.** Sage's own branch now defaults to symfn, so a bare
  `classical.init()` control compares symfn to symfn — it passes, it proves
  nothing, and a run of ratios all near 1.0x is the symptom
  ([python-and-sage-interop.md](docs/record/python-and-sage-interop.md)).
- "Memory" is three quantities that move independently — peak live bytes,
  total allocated, RSS. Use `symfn::measure` and `examples/heapstat.rs`;
  only the first two are assertable
  ([memory.md](docs/record/memory.md)).
- The number, its harness, and its caveats land in the record file that owns
  the subsystem — never in rustdoc.

## Definition of done

- `scripts/preflight.sh` is green.
- A new `pub` item meets the checklist at the bottom of
  [docs/style.md](docs/style.md) — summary sentence, contract, degenerate
  inputs, `# Panics`, a convention-pinning doctest, resolving citations.
- A new failure path fits a row of the mechanism table in
  [docs/policies/failure.md](docs/policies/failure.md).
- A new or changed `#[pyfunction]` sits in a row of the home table in
  [docs/policies/python.md](docs/policies/python.md).
- A new family or engine carries the evidence its row of the table in
  [docs/policies/validation.md](docs/policies/validation.md) requires; a
  timed case is a verified case.
- Measurements are in the record, with harness and power state.
- If the change closes or opens an item in
  [docs/todo-1.0.md](docs/todo-1.0.md) or a record file's open tail, the same
  change says so there.
