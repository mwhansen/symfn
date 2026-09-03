# Release readiness — from "works on this machine" to a package

[The record](record/) tracks what the library *computes*, and how it came to.
This file is the forward-looking half: what stands between the current tree and
something a stranger can depend on — continuous integration, a curated API
surface, a stated contract for failure,
and two publishable artifacts — a crate and a **Sage-free** wheel — with Sage
interoperability layered on top of the wheel rather than baked into it.

The mathematics is not the gap. The suite passes across every feature
combination, the oracles are committed, and the licensing is clean and audited.
What is missing is the operational layer — and one fact framed the whole
document when it was written:

> **Nothing in this repository has ever been built or tested anywhere except one
> macOS arm64 machine running rustc 1.96.** There is no `.github/`. Every claim
> in the README is, today, a claim about one laptop.

**That is no longer true, and what replaced it is worth keeping.** The tree is
at `github.com/mwhansen/symfn`, `.github/workflows/` holds a CI and a release
workflow, and the first run of the first one found 43 gate failures across
clippy, Sphinx and mypy — every one of them a tool version this laptop did not
have, and not one of them a regression
([record/python-and-sage-interop.md](record/python-and-sage-interop.md)). The
premise was right about the size of the exposure.

The `v0.1.0` tag this document once called "already tagged" is a local tag on a
commit 171 behind, from before any of that existed. It was never pushed and it
does not describe a shippable tree; the versions that reach anyone start at
`v0.1.0-rc.1`.

That is Phase 0, and almost everything else is easier once it exists.

---

## Phase 0 — CI, so the claims become checkable
*Blocks everything. Nothing below can be verified without it.*

- [x] `.github/workflows/ci.yml`: `cargo test` across {default, `bignum`} ×
      {ubuntu, macos, windows}, plus a **release-profile lane** — the tests
      pinning `overflow-checks` only carry information there
      ([policies/failure.md](policies/failure.md), R3) — and a `python` build
      job on Linux and macOS. `python` is built rather than tested: an
      `extension-module` test binary has no interpreter to resolve `_PyExc_*`
      against, which is why that boundary's pin lives in
      `scripts/check_schubert_bindings.py`. **Never executed** — this
      repository has no remote, so the workflow is written and unverified until
      one exists.
- [x] `cargo fmt --all --check` as a CI gate. `.githooks/pre-commit` already
      does this, but it is opt-in per clone (`git config core.hooksPath`), so it
      is a convenience, not an enforcement.
- [x] `rustup component add clippy`, then triage. ⚠️ The backlog was **180, not
      the 155 recorded here** — 155 was one `--all-features` run, and `#[cfg]`
      changes which code exists, so the union over {default, `bignum`, `python`,
      all} is the real number. 124 were fixed (83 by `--fix`, the rest by hand),
      2 lints were exempted crate-wide with the reason in `Cargo.toml`
      (`needless_range_loop`, `type_complexity` — both wrong for this domain),
      and the 56 `cast_*` in `tests/`/`examples/` stay visible-but-not-fatal.
      The advisory job is now a gate at `-D warnings` over all four feature
      sets. See [record/failure-and-overflow.md](record/failure-and-overflow.md)
      for the census and what each exemption buys.
- [x] Declare `rust-version` in `Cargo.toml` and add an MSRV job pinned to it.
      **1.87**, chosen rather than discovered: it is the floor that buys
      `u32::is_multiple_of` (1.87) and `iter::repeat_n` (1.82), which is what 37
      of the warnings above wanted. Verified by running both suites on 1.87.0,
      not inferred.
- [x] `cargo doc --no-deps --all-features` gated with `-D warnings`.
      `scripts/build_docs.sh` sets `RUSTDOCFLAGS=-D warnings` and the `docs`
      job in CI runs it, which is also what bundles the rendered pair. Until
      that job existed, CLAUDE.md's standing claim that this command is silent
      was checked by nothing.
- [x] Fix the 5 dead-code warnings `cargo package` surfaces. All four functions
      turned out to be **exercised by tests and dead only outside them**, and
      each documents something the live code no longer says, so they are kept
      with `#[allow(dead_code)]` and that reason; the unread `room` field was
      dead only because its reader was. Four more warnings went with them (two
      unused imports, an unused `mut`, and a non-local `impl` inside a test
      body — which implemented `Acc for u8` across the whole test module while
      looking local). The tree is now warning-free across
      {default, `bignum`, `python`} × all targets, which is what lets CI run
      `-D warnings`.
- [ ] A separate, non-blocking job for the Sage-dependent checks. **Nearly
      every script in `scripts/` imports Sage**, so they cannot run on a normal
      runner; put them behind a container image or a nightly schedule and let
      the fast suite gate PRs.

**Done when:** a push runs the full non-Sage suite on three platforms and three
feature sets, and a red build blocks merge.

---

## Phase 1 — documentation that renders

**`cargo doc --no-deps --all-features` is silent.** It emitted 186 warnings
when this phase was written (the 173 recorded here had grown), and on docs.rs
those were broken links scattered across the whole public API. The breakdown
as executed, which differed from the estimate above in two ways worth keeping:

| count | warning | resolution |
|---|---|---|
| 105 | citation key | link definitions in the module doc; `\[KEY\]` escapes in item docs — see [style.md](style.md), "Citations" |
| 54 | math read as a link | backticked: `` `ℚ[q,t]` ``, `` `ℤ[α]` ``, `` `f[g]` `` |
| 17 | public docs link to a private item | delinked to a plain code span; none of the 17 was a target worth making public |
| 14 | `X` is both a function and a module | `mod@` prefix on `crate::{kostka,character,charge,convert,plethysm}` |
| 7 | redundant explicit link target | dropped the explicit path |
| 8 | genuinely broken link | five needed a real path (`crate::lr::NaiveLr`, `crate::coeff::{QAlgebra,Field,Plethystic::frobenius}`, `Self::from_beta_numbers`); three point at feature-gated or optional-dependency items that cannot resolve in a default build and are now code spans |

Two corrections to the estimate. **Citations were the largest class, not
math** — 105 against 54, where this table had them lumped together at 142.
And **link definitions fix the module doc only**: rustdoc scopes them to the
doc comment they sit in, so the "real links are the better answer" plan works
for `//!` and leaves every `///` mention warning. Item docs escape instead;
[style.md](style.md), "Citations", records the split and why.

- [x] Math backticked — the rule narrowed to *bracketed* expressions only,
      since no bare glyph ever warned ([style.md](style.md), "Mathematical
      notation").
- [x] Citation keys given real targets, with `## References` entries added for
      the keys that lacked them. `[GJ]`, `[BH]` and `[GH]` are journal-only or
      unpublished and are escaped everywhere; there is nothing to link to.
- [x] The five function/module collisions disambiguated.
- [x] The public→private links resolved.

**Remaining for this phase:** the CI gate from Phase 0 that switches
`-D warnings` on, so this cannot regress.

**Seven measurements live only in rustdoc**, found by the 2026-08-07 audit
while moving figures to the record. The rule assumed the record already held
whatever rustdoc mentions, so the only moves available were "delete" or "keep",
and deleting would have destroyed the sole copy. Each needs re-measuring with a
named harness and a record entry, after which the rustdoc sentence goes:

| site | figure | record file that should own it |
|---|---|---|
| [rect.rs](../src/rect.rs) | 37.0 ms → 3.8 ms on `s(12⁶)²` | littlewood-richardson.md, which has a *different* run of the same workload |
| [skew_lr.rs](../src/skew_lr.rs) | "tens to hundreds of MB" transient | memory.md, whose "hundreds of MB" is the unrelated deep-clone figure |
| [jack.rs](../src/jack.rs) | 8.6× at n = 10, E1 vs branching | jack.md |
| [qtkostka.rs](../src/qtkostka.rs) | 8.8× at degree 9, BH vs branching | qt-kostka.md |
| [fasthash.rs](../src/fasthash.rs) | ~1.3× on `s[8,7,6,5,4,3]²` | llt.md, which has only the later MixHasher sweep |
| [frac.rs](../src/frac.rs) | 782 of 3300 profile samples | macdonald.md |
| [schubert.rs](../src/schubert.rs) | 1.3–1.6× for the perm-only peel memo | schubert.md, which discusses the merge but records no speedup |

In all seven the rustdoc holds the only copy, so applying the rule "the record
owns measurements" by deleting the figure would delete the measurement. Moving
each one needs the record file checked first, not the rustdoc edited first.

### Prose a reader outside this project can follow

Raised 2026-08-08, by the one reader who has read the whole tree: the code is
distributable and the writing is not. Figurative and maxim-shaped prose —
banned as of [style.md](style.md), "No aphorisms, no metaphors" — is dense
enough in places that a sentence's meaning depends on already knowing the
answer. The survey, by marker grep over the tree and `git log`:

| surface | figures | note |
|---|---|---|
| `src/` rustdoc | 4 sites, 0 maxims | ⚠️ **wrong — 16, see below** |
| `README.md` | 6 | the shop-window framing itself is one |
| `docs/` | ~49 marker hits, 298 bolded leads | [style.md](style.md) is the densest single file, at 23 |
| commit bodies | 45 of 219 | across all 16 branches, not just `main` |

⚠️ **The `src/` row was wrong, and the way it was wrong is worth recording.**
The survey was a marker grep run before the ban existed, so its markers were
the figures someone had already noticed. Re-scanned against the vocabulary the
sweep itself collected — `posture`, `gamble`, `launder`, `papering over`, `in
disguise`, `wearing`, `the payoff`, `the whole reason` — `src/` has **16**,
six of them `posture` alone. The conclusion drawn from the 4, that the shipped
surface was nearly clean and could go last, is what put this item at the
bottom of the list.

The surface a downstream consumer reads on docs.rs is `src/` plus the README,
so this blocks a contributor or a packager reading the repository well before
it blocks distribution.

The log is being rewritten rather than left alone: all 219 commit messages
across all 16 branches, by `git filter-repo`, decided 2026-08-08. That trades
the capture stage's contemporaneity for legibility, and the terms of the
trade are stated where the rewrite lands. Measurements are carried over
verbatim; only the prose around them is rewritten.

- [x] `src/` and `README.md` — **done**, 17 sites rather than the 10 the
      survey predicted. `cargo doc --no-deps --all-features` is still silent.
- [x] `docs/style.md` — **done**, and it was using four of the figures its own
      ban section names as examples: "filing the number off", "the shop
      window", "a verdict with no premise is a permanent wall", and "the record
      is the working agent's only long-term memory", the last two of which the
      ban section quotes *with their replacements already written*. It also
      used the construction its Voice section bans in the same breath, "the
      whole point of X", twice. Four of its claims were stale, including the
      spelling bullet still saying "converge on touch, don't sweep" after the
      sweep and its gate had landed.
- [x] The rest of `docs/` — **done**, all 12 record files plus the three
      policies, the digest, the four audits and the clean-room spec. Two
      findings from doing it in one pass rather than converging on touch.
      **Figurative density tracks how much of a file is about process rather
      than measurement**: `llt.md` (four bolded maxims) and `kronecker.md`
      (three) argue about how to work, while `macdonald.md` and
      `skew-and-evaluation.md` are almost entirely measurement and needed two
      edits and none. And **the same maxim kept turning up in three or four
      files at once** — "an entry point's name tells you what it computes;
      only its caller tells you what it must return" was in an audit, a
      record and a policy; "a counter that is reset while its subject is still
      live is signed" was in three; "a cost model with a factor missing will
      rank engines confidently and wrongly" was in two and is the sentence
      style.md quotes as its specimen aphorism. Each is now stated once, as
      the finding, in the file that owns it.
- [x] A marker-word grep in `scripts/preflight.sh` — **done**,
      `scripts/check_figures.py`, and it was cheap because the sweep produced
      its own vocabulary: the words that recurred, with the count each had
      when it was found. It fires on 48 lines at the pre-sweep revision and
      all 48 are lines the sweep rewrote, which is the measurement behind
      gating on it rather than reporting. It also found five sites the reading
      pass had missed, in `tests/`, `examples/` and two files already
      declared done. It covers roughly a third of what reading found and its
      docstring says so; the four bolded maxims in `llt.md` are invisible to
      it. The same argument was applied to spelling first, since there the
      grep is exact: `scripts/check_spelling.py` now scans `docs/`,
      `README.md` and `CLAUDE.md` as well as `src/`, and found 51 British
      forms in 19 files that the `.rs`-only version could not see, plus four
      stem gaps that had been invisible in `src/` since it was written.

---

## Phase 2 — decide what the API *is*

`src/lib.rs` declares **39 `pub mod`s and exactly one private one**. Publishing
that freezes every one of them under semver — including `gjmod`, `rect`,
`three_row`, `two_row`, `memo`, `modular`, which read as implementation strategy
rather than interface. The cost of getting this wrong is paid forever; the cost
of getting it right is one afternoon before the first crates.io release.

⚠️ **40, not 39** — the count above omitted the `#[cfg(feature = "python")]`
one. The sort landed at **29 API / 9 `#[doc(hidden)]` / 2 `pub(crate)`**, and
the tiers with their membership test are recorded as a comment above the module
list in [lib.rs](../src/lib.rs) so a new module is sorted the same way.

- [x] Sort the 39 into **API** (documented, semver-stable) and
      **implementation** (`pub(crate)`). The test for API membership: would a
      caller who only wants symmetric functions ever name it? Only `memo` and
      `modular` reached `pub(crate)`: nothing outside `src/` names either, and
      everything else demoted is named by a test or an example, so it landed in
      the hidden tier instead. Nothing outside `src/` needed an edit — every
      test and example already reached the demoted modules through paths that
      survive.
- [x] For anything that must stay public but is not a stable promise — the
      alternative LR backends, the cross-check engines that exist to disagree
      with each other — mark it `#[doc(hidden)]` or gate it behind an
      `unstable-internals` feature, and say so in the README. **`#[doc(hidden)]`
      chosen**; the feature would have cost a fifth CI feature set and
      `required-features` on the tests that reach these modules, for no wall the
      README does not already state. The nine: `bh`, `gjmod`, `macop` (the
      cross-check engines), `rect`, `two_row`, `three_row`, `strip_lr` (the
      product strategies), `measure`, `python`. The strategy modules keep their
      root re-exports **documented** — `okada_product`, `two_row_product`,
      `three_row_product` and `AutoLr` are mathematical results and the default
      LR backend, and hiding the module path is the whole demotion.
      `gjmod::engines_agree` and the `macop` re-exports are hidden with their
      modules, because a predicate that runs two engines and compares them is a
      test helper rather than a computation; `qtkostka`'s three `_via_*` routes
      stay API on the same distinction, since they return the table.
- [x] Add `#![deny(missing_docs)]` once the surface is small enough to hold that
      line. It found **53**, of which the sort had already retired 17 — the lint
      skips `#[doc(hidden)]` items as well as private ones, which is the
      concrete return on sorting first. The 36 that remained were accessors,
      six struct fields, two `PartitionError` variants, and unwritten
      trait-method signatures on `SymFn`, `Ring`, `ToSchur`/`FromSchur`,
      `SkewBy` and `Dual`; not one basis conversion was among them.
- [x] Write the semver policy into the README: what 0.x means here, what will
      break, and that the coefficient-ring traits (`Ring`, `QAlgebra`,
      `Plethystic`) are the ones consumers build on. Under "The public API, and
      what a version number promises", which also names the two breaks that do
      not look like breaks: a method added to `Ring`/`SymFn`/`LrBackend`/
      `SkewBy` breaks external implementors while breaking no caller, and the
      Python surface freezes harder than the crate rather than in step with it.

      Correction, 2026-08-25: the strict 0.x promise is withdrawn until the
      first release — with the crate at 0.1.0-rc there is nothing for a
      minor-is-breaking rule to govern yet, and committing to one now would
      promise on a number line no consumer holds. The enumeration work above
      stands; its home moved to [public-api.md](public-api.md) in the README
      restructure, whose pre-release section now says the policy is decided at
      the first release and keeps the two not-look-like-breaks as facts. Deciding
      the policy is now an open item in Phase 4, beside the CHANGELOG item,
      not a standing README section.
- [x] Same exercise for Python: ⚠️ **108 entry points** — 98 when the sort
      ran, and 108 once the Cython branch merged — not the 87 this file
      first counted, the 91 it then said, or the 92 a re-grep gives. Every one
      of those came from grepping `#[pyfunction]`, which undercounts twice
      over: two attributes sit inside `out_of_schur!` and `into_schur!` and
      expand to nine conversion entry points between them, and one apparent
      match is the string `#[pyfunction]` inside a doc comment. Count from the
      `#[pymodule]` registration block; `dir(symfn)` agrees with it.

      **All of them are supported and none is harness-only** — a result, not a
      deferral. Every entry point does a whole-object operation with a named
      audience, because [policies/python.md](policies/python.md) P2 has been
      enforced since the file was written, so the fine-grained probes the
      category exists to absorb never accumulated. The indexed and bulk entries
      the adapter needs are supported *and* documented as low-level, which is
      the carve-out P10 already made. `symfn.pyi` is the list made
      machine-readable and `scripts/check_python_stubs.py` holds the two
      together; `__version__` comes from the crate version, and the
      `#[pymodule]` carries the docstring it had none of.

      **The sort found a defect first**: 27 of the 98 named their first
      argument `lambda`, a Python keyword, so none of them could be called with
      keyword arguments and none could be stubbed. All 27 are now `la`. See
      [record/python-and-sage-interop.md](record/python-and-sage-interop.md).

**Done when:** the public module list is a deliberate list, and every item on it
has module-level docs. **Done** for the crate: 29 API modules, 9
`#[doc(hidden)]`, 2 `pub(crate)`, `#![deny(missing_docs)]` holding the line, and
the policy in the README. **Done on the Python side too**, membership and
docstrings both: all 108 carry a `# Raises` section and an example, 132 of
them executed by `scripts/check_python_docs.py`, and P6's order statement
sits in the `#[pymodule]` doc with the three departures from it named at the
entry points that depart. [policies/python.md](policies/python.md) delta 1 is
closed; the sweep is recorded in
[record/python-and-sage-interop.md](record/python-and-sage-interop.md).

A 2026-08-25 interface review of the convenience layer reopened the Python
half in one narrow sense: the entry-point list stands, but some behaviors
and documented semantics beneath it are defects — an element that prints
`0` without equaling `0`, a constructor that silently drops a term, an
identity conversion that raises. The fixes are staged in
[plans/convenience-surface-review.md](plans/convenience-surface-review.md);
its stages 1–3 precede this phase's freeze, and all four stages are done as
of 2026-08-25 — the deformed pairings, the power-sum opening, partial `at`
and the LLT skew tuples all landed the same day, so what the plan leaves is
one recorded open (rational alphabets in `evaluate`) with no work planned.

---

## Phase 3 — a stated contract for failure

The count was **138 `panic!` / `unwrap()` / `expect()` sites in `src/`**; at
audit time, outside tests, it was 93. That is not automatically wrong — for a
library whose inputs are partitions, "this partition is not a partition" is a
programming error, and panicking is the right answer. What was wrong is that a
caller reading the docs could not tell which inputs panic, which return
`Result`, and what happens on overflow.

- [x] Write the policy down — [policies/failure.md](policies/failure.md):
      contract violations panic and say so, reachable states refuse loudly,
      overflow escalates or refuses and never wraps. Its caller-facing half is
      now the crate front page ([lib.rs](../src/lib.rs), "Exactness").
- [x] Audit them against that rule. The live defect was the Python boundary
      panicking on a malformed permutation; fixed structurally, so the
      escalation path has no `unwrap` to make. The `# Panics` sweep that
      followed took every remaining `pub fn`: 46 sections added, from four in
      the whole crate, and no bare `.unwrap()` left in `src/`. It found a
      second defect — `Perm::at` returned `w(0) = 0` rather than panicking.
      See [record/failure-and-overflow.md](record/failure-and-overflow.md).
- [x] Document the overflow story properly: what a caller gets from each
      coefficient type, on the front page rather than only in the README, with
      each fixed-width family's own wall stated where that family lives.
- [ ] Confirm the two `unsafe` blocks are justified with `// SAFETY:` comments,
      or add `#![forbid(unsafe_code)]` to the modules that do not need them.
      (`skew_lr.rs`'s has one; `measure/`'s `GlobalAlloc` impl forwards to
      `System` and has not been reviewed under this heading.)

**Done when:** every public function that can panic says so, and the overflow
contract is on the type, not only in the README. **Both halves are now done** —
the front page carries the contract, and the `# Panics` sweep
([style.md](style.md), delta 2) covered every `pub fn` in `src/` outside
`python.rs`. What is left in this phase is the `unsafe` review above.

---

## Phase 4 — the crate as a publishable artifact

`cargo package` used to warn *manifest has no documentation, homepage or
repository*. It no longer warns about anything.

- [x] Fill in `repository`, `homepage`, `documentation`, `readme`, and
      `rust-version`. `documentation` is `docs.rs/symfn`, which is deliberately
      not the Read the Docs URL in `pyproject.toml` — two surfaces, two
      references.
- [x] `exclude` added, and the question it hinged on decided explicitly:
      **`cargo test` on a published crate is meant to work**, because that is
      what distro packagers do with the tarball, and a suite that cannot run is
      worse than one never shipped. So `tests/` and its 440 KB of fixtures
      stay — excluding them would not slim the crate but break it, since both
      oracle files are `include_str!`d into the test binaries.

      Out go `scripts/`, `docs/`, `docsite/`, `python/` and the dot-directories:
      measured, nothing in `src/` or `tests/` reads any of them at build or
      test time, and the `docs/` pointers throughout the rustdoc are prose
      rather than `include_str!`. `python/` leaves the *crate* only — maturin
      builds the wheel and sdist from the working tree, which was confirmed by
      rebuilding the sdist and running `scripts/check_sdist_offline.sh`
      against it.

      **234 files and 3.5 MiB before, 124 and 2.2 MiB after** — 600 KiB
      compressed against 1.0 MiB. Verified by unpacking the packaged crate and
      running `cargo test` inside it: both oracle suites ran there.
- [x] `cargo publish --dry-run` is clean, and
      `[package.metadata.docs.rs] features = ["bignum"]` is set. docs.rs builds
      `--all-features` by default, which would build PyO3's `extension-module`
      cdylib with no interpreter to resolve `_Py*` against — the same link
      failure this tree hit twice in one week from two other causes. It would
      have failed *after* publishing, on a machine nothing local reproduces.
- [x] `cargo-deny` is a CI gate, configured by `deny.toml`. Measured at the
      revision it landed: 19 dependencies across all features, every one
      permissive, so the allow-list is the exact set rather than a category —
      `MIT`, `Apache-2.0`, `Apache-2.0 WITH LLVM-exception` (target-lexicon)
      and `Unicode-3.0` (unicode-ident, in a conjunction, so it genuinely
      applies).

      The gate was negative-tested rather than assumed. Removing `MIT` from the
      allow-list still passes, because every dependency is `MIT OR Apache-2.0`
      and the disjunction is satisfied by the other half — so that test proves
      nothing. Removing `Unicode-3.0` fails, because it arrives in a
      conjunction. That is the one that shows the gate bites.

      It also denies any source but crates.io, which would defeat the offline
      sdist as well as the licensing claim.
- [ ] `CHANGELOG.md`, starting from `v0.1.0-rc.1` — the first tag that is
      pushed, built and downloadable. Not from the local `v0.1.0`, which this
      item used to name: it sits 171 commits back, predates CI, and nothing was
      ever built from it.
- [ ] Decide the versioning policy at the first release: what a minor and a
      patch may change, and whether the Python surface's harder freeze is a
      stated rule. The material for the decision is in
      [public-api.md](public-api.md) — the tiers, and the two breaks that do
      not look like breaks. Deferred from the Phase 2 semver item, whose
      2026-08-25 correction says why.

**Done when:** `cargo publish --dry-run` is clean and the docs.rs build is
verified.

---

## Phase 5 — the wheel, and it does not know Sage exists

**The invariant: the `symfn` wheel has no Sage dependency, at build time or at
import time.** It is a standalone symmetric-function library for any Python.
Sage is one consumer of it, and the adapter that makes that work lives on the
Sage side of the boundary (Phase 5b) — never inside the package. That
invariant, and the three-layer boundary that enforces it, are standing policy
now — [policies/python.md](policies/python.md) — and several items below are
that file's deltas; this checklist remains the execution plan.

This is already true and worth keeping true deliberately: `src/python.rs`
contains **zero** references to Sage. All the coupling is in
`scripts/sage_backend.py` and `scripts/symfn_cy.pyx` (deleted since; the
adapter lives in Sage now), the second of which `cimport`s Sage's
`Integer` type and therefore cannot build without Sage present. The work here is
to make that separation enforced and packaged rather than incidental.

- [x] `pyproject.toml` with `[build-system] requires = ["maturin>=1.5,<2.0"]`
      and a `[project]` table: description, README, license, classifiers,
      `requires-python = ">=3.9"` (matching the `abi3-py39` build).
      **`[project.urls]` and `Cargo.toml`'s `repository` both name
      `github.com/mwhansen/symfn`.** They were left empty until there was a
      published home, on the grounds that two URLs written separately disagree;
      they were filled in together, in one commit, when the repository existed.

      **It also carries `[tool.maturin]`, which is what made the convenience
      layer possible**: `python-source = "python"` and
      `module-name = "symfn.symfn"` turn the wheel into a mixed project, so
      there is somewhere to put Python at all. A wheel that is only a compiled
      module has no such place, which is why `docs/policies/python.md` delta 5
      waited on this item rather than the other way round.

      And `[tool.ruff]`, which lints the convenience layer — `E`, `F`, `I`,
      `B`, `C4`, `UP`, `D`, with the pydocstyle rules that disagree with
      `docs/style.md` turned off and each exclusion carrying its reason.
      **No Sage in `dependencies`, and no optional extra that pulls it in** —
      Sage is not pip-installable in the normal case, and an extra advertising
      it would be a lie.
- [x] `__version__` in the Python package, sourced from the crate version so
      the two cannot drift. `#[pymodule]` adds it from `CARGO_PKG_VERSION` and
      `symfn/__init__.py` re-exports it; `pyproject.toml` takes it from Cargo
      through maturin's `dynamic = ["version"]`, so all three are one value.
- [x] Wheel matrix via `maturin-action`, in
      [.github/workflows/release.yml](../.github/workflows/release.yml). It is
      the full `rpds_py` platform set — macOS x86_64 and arm64; manylinux
      x86_64, aarch64, armv7l, ppc64le, s390x, i686; musllinux x86_64, aarch64,
      i686; Windows win32, amd64, arm64 — as
      [docs/sage-packaging-audit.md](sage-packaging-audit.md) argued it should
      be, and `abi3-py39` turns it into **14 artifacts** rather than that
      package's 55.

      **The QEMU legs this item expected are not there, and are not needed.**
      `abi3` builds without an interpreter for the target, so maturin
      cross-compiles armv7l, ppc64le, s390x and i686 inside its manylinux
      containers at the cost of an ordinary compile. What emulation would have
      bought is *testing*, not building; the native legs run an
      import-and-compute step and the cross legs do not, which
      [support-tiers.md](support-tiers.md) records as accepted exposure.
- [x] **A support-tier policy, written down**, at
      [docs/support-tiers.md](support-tiers.md): the fourteen Tier 1 platforms
      by wheel tag, what Tier 2 requires of a builder, and what moving a
      platform between them costs.
- [x] **An sdist that builds offline.** `scripts/build_sdist.sh` vendors the
      `python` feature's four dependencies and writes the
      `.cargo/config.toml` that redirects crates-io at them, both of them
      untracked and present only inside the artifact; `pyproject.toml`'s
      `[tool.maturin] include` is what carries them in. 4.4 MB, against 1.4 MB
      for the wheel.

      `scripts/check_sdist_offline.sh` is the assertion, and it runs in the
      `sdist` job of the ordinary CI workflow rather than only at release:
      the way this breaks is a new dependency landing in `Cargo.toml`, which is
      an ordinary commit. Offline is asserted rather than simulated —
      `CARGO_NET_OFFLINE=true` and pip's `--no-index` make a missed dependency
      a hard error naming the crate, where cutting the runner's network would
      have tested the runner.
- [x] `symfn.pyi` type stubs. The API is coarse-grained and takes
      list-of-`(partition, coefficient)` pairs; stubs are the difference between
      that being discoverable and being guesswork. Landed with Phase 2's sort,
      since the stub file is what makes the supported list machine-readable:
      all 108 entry points, named type aliases where a Rust alias already drew
      the distinction a structural type loses, and
      `scripts/check_python_stubs.py` failing on any drift in name, arity or
      parameter name. It now sits at `python/symfn/symfn.pyi`, beside the
      compiled module it describes, which is where a type checker looks for it
      — at the repository root it typed nothing once the layout became mixed.
      `py.typed` ships beside it, so the convenience layer's inline annotations
      are read too.
- [x] A Python test suite that runs **without Sage** — round-trip the marshalling
      layer against values computed in Rust. `check_bindings.py` already tests
      the boundary rather than the library, which is the right idea; it just
      needs a Sage-free sibling that CI can run on a stock runner.
      `scripts/check_python_boundary.py` is the first such sibling and covers
      the *failure* half: 128 malformed calls over 94 pyfunctions, asserting
      only typed exceptions come back. `scripts/check_convenience.py` is the
      second, and holds the convenience layer to the contract layer over 2177
      checks. The round-trip half is `scripts/check_python_marshalling.py`
      (2026-08-21), and it is a marshalling test rather than an oracle: a
      value handed in comes back out intact, at the widths and shapes
      `docs/policies/python.md` P1 promises — every export's return shape,
      coefficients at and past both `i128` edges, the permissive inbound
      spellings. Catching a kernel defect is not this surface's job —
      `docs/policies/validation.md` owns that, and the Rust suites and
      `tests/fixtures/` discharge it. Its first run found two boundary
      defects at `i128::MIN`, one of them a silently wrong value
      ([record/python-and-sage-interop.md](record/python-and-sage-interop.md)).
- [x] A CI assertion that the invariant holds: import `symfn` in a bare
      interpreter with no Sage on the path and exercise the public API. That is
      the test that stops a convenience import from creeping in later. The
      `wheel` job does it against an *installed* wheel rather than the source
      tree, because the failure it guards is the packaging one — a pure-Python
      module left out of the wheel imports fine from the tree and not at all
      from a `pip install` — and it fails if any `sage` module reaches
      `sys.modules`.
- [x] Decide how coefficients cross the boundary for non-Sage callers. Ints and
      `Fraction` come through natively via PyO3's `num-bigint`/`num-rational`
      conversions; the `q,t` and `α` types cross as the exponent-keyed rows
      `docs/policies/python.md` P1 specifies — factored denominators for
      Macdonald, primitive atoms for Jack — and the convenience layer's `Poly`,
      `QtPoly`, `QtFrac` and `AlphaFrac` wrap exactly those rows with a `repr`
      and an `at`. No Sage ring is assumed anywhere, and the specializations
      `at` makes expressible are what check the conventions from Python.
- [x] A publish workflow for PyPI and crates.io, with trusted publishing on
      both sides — PyPI's OIDC exchange and `rust-lang/crates-io-auth-action`,
      so neither registry needs a long-lived token stored in the repository.

      **It is not publish-on-*tag*, which is what this item asked for.** A `v*`
      tag builds all fifteen artifacts and attaches them to a GitHub Release,
      and stops there; publishing to either registry takes a deliberate
      `workflow_dispatch` naming the registry, behind a GitHub environment
      whose required reviewer holds even when the dispatch is wrong. The reason
      is that the first testers install from the Release page before the name
      goes to a registry at all, and a registry publish is the one step here
      that cannot be undone — a version yanked from PyPI or crates.io can never
      be reused, so a wrong 0.1.0 costs the number permanently. Making the
      not-yet state structural is cheaper than remembering it.
- [x] **The README's Install section, rewritten for the registry names.** It
      opened "Not yet on crates.io or PyPI" and gave a `cargo add --git` line
      and a `pip install` of one wheel URL from the `v0.1.0-rc.2` Release page;
      all three are gone, replaced by `cargo add symfn` and `pip install
      symfn`. The pinned URL was the part that failed quietly rather than
      loudly: it named one asset of one tag, so it kept working — and kept
      installing a release candidate — for as long as that Release exists,
      which is indefinitely. The same change adds the crates.io, PyPI and
      docs.rs badges to the header, and drops the same not-yet-published note
      from `docsite/install.md`.

      **These lines resolve only once the registry dispatch runs.** The
      rewrite went in ahead of it rather than in the same change, so the
      window in which the README names a package neither registry carries is
      the window between this commit and the first dispatch; closing it is
      what the dispatch is for.
- [x] **The rendered reference**, at `docsite/`, published by Read the Docs.
      Sphinx with MyST, building the wheel first so both layers are documented
      from the objects themselves and `help()` cannot drift from the website.
      `scripts/check_docs_complete.py` compares `symfn.__all__` against the
      inventory Sphinx writes and fails when a supported name reaches no page:
      Sphinx warns about references that do not resolve and says nothing about
      an entry point nobody wrote a directive for, which is the hole that opens
      as the surface grows. All 108 contract entry points and 79 convenience
      names are covered.

      **Both references also ship as a downloadable bundle**,
      `symfn-docs-<version>.tar.gz` — 5.1 MB, self-contained, no server.
      `scripts/build_docs.sh` renders the Sphinx site and `cargo doc` and puts
      a landing page over the pair, because Read the Docs publishes the Python
      half and nothing publishes the Rust half; a tester who downloaded a wheel
      would otherwise have no reference to read beside it. It is attached to
      every Release and uploaded from every CI run.

      That script is also the first thing holding `cargo doc --no-deps
      --all-features` to being silent — a claim CLAUDE.md has made since before
      there was CI and nothing checked. It denies rustdoc warnings.

**Done when:** `pip install symfn` works on Linux, macOS and Windows without a
Rust toolchain, and the package imports and computes with no Sage anywhere.

---

## Phase 5b — the Sage extension module, on the far side of the boundary

*What is true today — the coverage claim, where the adapter lives, how to turn
it off — is stated once in [sage-backend.md](sage-backend.md). This phase and
the next are the plan around it.*

Sage integration stays possible, but as a **separate artifact** that depends on
`symfn` rather than the other way round. When this phase was written it existed
as two files in `scripts/` that were built to run experiments, not to be
installed by anyone: `sage_backend.py` (monkey-patches
`sage.combinat.sf.classical.conversion_functions` in a live session) and
`symfn_cy.pyx` (the compiled per-term loop, which needs Sage's headers to
build). Both are deleted; neither path resolves any more.

**That is no longer where the adapter lives.** It is a branch of Sage —
`mwhansen/sage`, branch `symfn` — carrying `src/sage/libs/symfn/`
(`backend.py`, `extras.py`, `terms.pyx`), `src/sage/features/symfn.py`,
`build/pkgs/symfn/`, and the call-site changes across `combinat/sf/`,
`partition.py` and `schubert_polynomial.py`. It builds and its doctests pass;
what has not happened is the *upstream* part, which is Phase 5c.

- [x] **Decided: the third option**, the adapter inside the Sage codebase. The
      deciding argument is the one this list already gave — `terms.pyx`
      `cimport`s Sage's `Integer` and so must be compiled against a specific
      Sage build, which is a build-per-version problem maintained externally
      and a non-problem maintained inside. The first two options are recorded
      below as the paths not taken:
      - **An in-tree module Sage imports** — the adapter lives here under, say,
        `sage/`, is not part of the wheel, and is installed by pointing a Sage
        session at it. Lower ceremony; the natural next step from where the
        scripts already are, and the right *interim* answer.
      - **A separate distribution** (`symfn-sage`) that depends on `symfn` and
        is built inside a Sage environment, shipping the Cython shim compiled
        against Sage's `Integer`. Cleaner boundary; needs its own build and
        release path, and the wheel matrix problem is *harder*, because a wheel
        carrying compiled Sage-linked code must match a Sage build.
      - **Upstream, in the Sage codebase** — see below.
- [x] Resolved, though not the way this item describes. The adapter was not
      lifted and relocated — it was **rewritten inside Sage**, as
      `sage/libs/symfn/backend.py` at 1001 lines against the script's 404. It
      imports `symfn` as an ordinary installed package, and the entry point is
      not `install()` but `sage.combinat.sf.classical.init()`, which populates
      the conversion table conditionally at import instead of monkey-patching a
      live session.

      **And the copy that was left behind is gone.**
      `scripts/sage_backend.py`, `scripts/symfn_cy.pyx` and
      `scripts/setup_cy.py` were a second, older implementation of the same
      adapter that no preflight ran — the drift this tree records elsewhere.
      `scripts/check_backend.py` and `scripts/bench_backend.py` now drive
      Sage's own adapter instead, switching arms through `SAGE_DISABLE_SYMFN`
      as they already did and **asserting** which backend answered rather than
      installing one. That assertion is new and is the stronger arrangement: an
      unverified control is the failure CLAUDE.md records as "a run of ratios
      all near 1.0x".

      The repointed harness reproduces the recorded figure exactly —
      **8647 computations at degree 8, 0 mismatches** — against
      `sage/libs/symfn/` rather than the deleted script, which is what says the
      deletion cost no coverage.
- [x] Documented in two places on the Sage side: `build/pkgs/symfn/SPKG.rst`
      states both, and `sage/libs/symfn/__init__.py` names which module serves
      which. The two modes are **backend replacement** (drop into Sage's conversion table and
      accelerate everything Sage already does) and **direct use** (call symfn
      for the things Sage has no equivalent for — reduced Kronecker via the `st`
      basis, the Macdonald operator algebra, the (q,t)-Kostka tables).

### Becoming a complete drop-in

The adapter as this section was written displaced **one** of Symmetrica's six
consumer sites —
the conversion table in `combinat/sf/classical.py`. That is why it can claim
4678 comparisons and still not be a replacement. The other five call sites reach
Symmetrica directly, and
[docs/record/python-and-sage-interop.md](record/python-and-sage-interop.md)
enumerates them. **These are ordinary library and binding tasks — they need no
upstream involvement and can be done at any point**, which is why they belong
here rather than in Phase 5c:

- [x] **`compute_*_with_alphabet`** (`combinat/sf/sfa.py:5656`). ⚠️ The audit
      called this a binding gap around `eval()`; it was not. `sfa._expand` needs
      a **polynomial in `n` indeterminates**, which it hands to `resPR(...)`,
      while `eval()` returns a single value from a `Ring` alphabet. Closed by a
      new operation — `Monomial::expand` in `eval.rs`, returning exponent
      vectors — bound as `expand_alphabet`, with the other five bases reaching
      it through `m`.
- [x] **`mult_monomial_monomial`** (`combinat/sf/monomial.py:129`). ⚠️ The
      pre-existing `Monomial::mul` was in `convert.rs`, not `sym.rs`, and went
      `m → s → LR → m` — which costs the whole degree, because `m → s` inverts
      the Kostka matrix. Replaced by the direct overlay rule in `sym.rs`; the
      Schur route stays as the reference oracle, the role `NaiveLr` plays for
      LR.
- [x] **`kostka_tab`** (`combinat/tableau.py:7016,7036`) — the audit's one
      correctly-identified gap, and the only one. `kostka::semistandard_tableaux`
      enumerates the chains the counting DP merges. Its **order is contract**:
      Sage prints the list in its doctests, so it is increasing lexicographic in
      the row-major reading word, checked against Symmetrica over all 1818 (λ, μ)
      pairs through degree 9.
- [x] **Wire `hall_littlewood`** (`combinat/sf/hall_littlewood.py:28`). This one
      site binds the function at import rather than reaching it through the
      module, so the adapter rebinds the name in `sf/hall_littlewood`'s own
      namespace; the other five all resolve through
      `sage.libs.symmetrica.all` at call time.
- [x] Extend `scripts/check_backend.py` to cover the newly-intercepted sites.
      **8647 computations at degree 8, 0 mismatches**, up from 4678.

Covered and needing no work: the 20 conversions, `kostka_number`.

- [x] **`scalarproduct_schubert`** (`combinat/schubert_polynomial.py:350`) — ⚠️
      **a gap the audit missed**, and then a decision it got wrong. It returns a
      Schubert *polynomial*; `schubert_pairing`, which the audit paired it with,
      returns an integer. The map is `∂_{w₀⁽ⁿ⁾}(S_u·S_v)`, with the pairing as
      its coefficient at the identity, and it ships as
      `Schubert::scalar_product` / `symfn.schubert_scalar_product` with `n`
      explicit.

      **The reversal, since the decision is on the record.** "Won't do" rested
      on nothing in sagelib calling `SchubertPolynomial.scalar_product`, which
      is still true. It is not the deciding fact: the method is public API on a
      user-facing class, not a low-level export, so third-party code can depend
      on it without its author knowing Symmetrica was underneath. See
      [docs/record/schubert.md](record/schubert.md).
- [x] **Wire `schubert_polynomial.py` in the adapter.** Done, on the Sage side
      at `mwhansen/sage` branch `symfn`. Five sites dispatch to symfn whenever
      it is present; the adapter supplies `scalar_product`'s rank as the
      longest one-line form across both supports, which is well defined because
      Sage strips trailing fixed points before the call.

      **The exception-fidelity requirement this item carried was dropped, not
      met.** It read: Sage's doctests assert `ValueError`s symfn answers
      instead, so a drop-in has to reproduce them. Those refusals turned out to
      live behind `divided_difference`'s explicit `algorithm='symmetrica'`, so
      they were never on the default path. symfn is a third `algorithm` value
      rather than a rebinding of that one, and it is now the default when
      installed — a caller who wants Symmetrica's refusals still asks for them
      by name. 1404 comparisons against `algorithm='sage'` agree, 136 of them
      in the range Symmetrica refuses.

      `\delta_i` at `i <= 0` is the case that cut the other way: there is no
      value, so Sage's own message is raised for both `'sage'` and `'symfn'`
      rather than letting a backend phrase it.

      Symmetrica still blocks on an interactive prompt at teardown after
      `divdiff_perm_schubert`, so the comparison stays a standalone probe.
- [x] **Superseded by the shape decision, and worth stating rather than
      quietly dropping.** This item existed because a shim maintained *outside*
      Sage might not build on a given Sage install, so a pure-Python fallback
      was the safety net. Inside Sage there is nothing to fall back from:
      `terms.pyx` is registered in `src/sage/libs/meson.build` and compiles
      whenever Sage does, exactly like every other `.pyx` in the tree.
      `backend.py` imports `build_terms` unconditionally, and that is now
      correct rather than a missing guard.

      The measurement that motivated it stands and is why the shim exists at
      all: the pure-Python per-term loop ran ~185 ns/term, which came to 0.76×
      the entire Rust computation it wrapped.
- [ ] The Sage-dependent CI job from Phase 0 is what tests all of this, and it
      is the only place Sage ever appears in the build graph.

**Done when:** a Sage user installs `symfn` from PyPI, installs or points at the
adapter, and gets both modes — with the adapter's absence costing the wheel
nothing.

**Where that stands:** the adapter half is built and verified; the *installs
from PyPI* half is not, and it is the same blocker the whole release story has.
`build/pkgs/symfn/requirements.txt` asks for `symfn >=1.0.0rc1` and
`SPKG.rst` points at `pypi.org/project/symfn/`, which does not exist yet — so
today the only route is a wheel downloaded from a GitHub Release on a private
repository. Nothing about the Sage side moves until symfn is published.

---

## Phase 5c — upstreaming the adapter into Sage
*The intended end state, and deliberately not the first move. The present-tense
statement is [sage-backend.md](sage-backend.md); this is the route to it.*

The adapter's natural long-term home is **inside the Sage codebase**, for one
concrete reason: the per-term loop `cimport`s Sage's `Integer` and therefore must
be compiled against a specific Sage build. Maintained externally, that means a
build per Sage version and ABI — the version-chasing problem in its purest
form. Compiled as part of Sage, it is just another Cython file, and the problem
does not exist. The pure-Python adapter is far less version-sensitive; the
compiled per-term loop is the piece with no good home outside.

This is a well-trodden arrangement, not a special case. **`lrcalc` is the exact
precedent**: an independently released library with its own versioning, wrapped
by thin code that lives in Sage (`sage/libs/lrcalc/`) and is built and tested
with it. Symmetrica is the same shape.

### The target is displacement, not coexistence

**Decided: symfn replaces Symmetrica in Sage rather than sitting beside it.**
That is a strictly larger project than adding an optional package, and it
changes three things about what has to be true first.

**1. Coverage becomes total rather than selective.** Being faster on the paths
we chose stops being the bar. Every `sage.libs.symmetrica` entry point that
anything in Sage calls needs an answer: covered by symfn, reimplemented, or
deprecated. The symmetric-function core of that is largely done — the 20
conversions, Kostka, characters, plethysm, Hall–Littlewood, Jack, Macdonald —
and `record/schubert.md` already designed the bindings to mirror the seven
Symmetrica Schubert entry points, so `sage/combinat/schubert_polynomial.py` is
anticipated too.

The obvious worry was what Symmetrica does *outside* symmetric functions —
[docs/record/README.md](record/README.md)'s "Beyond the core (deferred, but intended)" list:
modular and projective representation theory of the symmetric group, Hecke
algebras of type A, finite group operations, classical groups. **The audit
retired that worry: Sage calls none of it.** Those entry points are exported but
unreached, which makes them a deprecation question rather than an implementation
one.

- [x] **Do the audit first.** Done, and its census is now in
      [docs/record/python-and-sage-interop.md](record/python-and-sage-interop.md).
      **Sage reaches 36 of the 66 exported entry points, from six files**, and
      the whole representation-theory half of Symmetrica — the part that would
      have been a research programme — **is never called by Sage at all**. Its
      per-entry-point gap list was then corrected by implementing it: **all 36
      are computed by symfn today** and 29 are intercepted.

**The coverage requirement is closed.** Every entry point Sage reaches is
computed here, `scalarproduct_schubert` included — the last one, and the one
this plan and the audit had both agreed to leave behind. What remains is
binding work inherited from Phase 5b: five of the six consumer sites are
intercepted, and the sixth, `schubert_polynomial.py`, is now unblocked and
listed under
[Becoming a complete drop-in](#becoming-a-complete-drop-in).

⚠️ **The accurate claim is narrower than "symfn replaces Symmetrica".** Sage
reaches 36 of the 66 entry points Symmetrica exports; the other 30 are
reachable only by a user who writes `from sage.libs.symmetrica.all import ...`,
and symfn does not cover them. So: **symfn covers every Symmetrica entry point
Sage calls** — which is what the displacement needs, and no more than that.

What is genuinely new at this phase is a policy question, not code:

- [ ] Decide, upstream, what happens to the **30 entry points Sage never
      reaches**
      ([the census](record/python-and-sage-interop.md) lists them): reimplement
      them, deprecate them through Sage's normal cycle, or keep Symmetrica
      installable as an optional package for whoever imports it directly. That
      decision belongs to Sage and is not settled here.

      ⚠️ **This item used to settle it, and what let it is gone.** It read
      "propose demoting Symmetrica from standard to optional", answering 31
      entry points at once — the 30 plus `scalarproduct_schubert` — and it was
      built around that last one existing: the displacement pitch needed an
      asterisk while a method a user could already be calling had no
      replacement, and demotion was the cheapest way to carry it. With the gap
      closed the asterisk goes, and with it the reason to bundle the 30
      unreached exports into the same answer. They are an ordinary deprecation
      question now, on their own merits.

**2. Standard package, not optional — and the goal is that Sage installs a
wheel, never a compiler.** A replacement for a standard package must itself be
standard, so the platform question is unavoidable. The intended answer is
Phase 5's matrix: `abi3-py39` means one wheel per *platform*, and Sage on
Windows is WSL, so the set Sage actually needs is small — manylinux x86_64 and
aarch64, musllinux, macOS x86_64 and arm64. An end user installing Sage should
never need cargo.

**This is settled, and it is settled in our favor** —
[docs/sage-packaging-audit.md](sage-packaging-audit.md) has the evidence:

- **131 of Sage's 272 standard packages are already distributed as prebuilt
  wheels**, with multi-platform wheel support documented in
  `build/sage_bootstrap/package.py` and `*.whl)` install branches in
  `build/bin/sage-spkg`.
- **`rpds_py` is the precedent and it is exact**: `type: standard`, a Rust
  extension, **built with maturin**, shipped as 55 platform wheels plus an
  sdist, with *no* `spkg-install` because nothing is compiled.
- **`clarabel`** ships `cp39-abi3` wheels — symfn's exact build configuration —
  which is why abi3 collapses `rpds_py`'s 55 artifacts to 5.
- **Sage has no `rust` or `cargo` package and no `spkg-install` anywhere invokes
  cargo.** Sage does not build Rust from source; it consumes Rust wheels. A
  symfn spkg introduces no new policy.

Two qualifications that survive the audit:

- **Downstream packagers still build from source.** Debian, conda-forge, Gentoo
  and nix build everything from source on principle. Rust is a *packager*
  problem, not a Sage problem and not an end-user problem — which is what
  Phase 5's offline `cargo vendor` sdist is for.
- **Platform reach is the real risk, and it replaces the toolchain concern.**
  Symmetrica is C and compiles anywhere; a wheel only reaches platforms someone
  built for. `rpds_py` sets the bar for a *standard* package and covers armv7l,
  ppc64le, s390x, i686 and Windows arm64 — beyond what Phase 5 currently lists.
  Match that before proposing: "fewer platforms than the package you are
  displacing" is a concrete, fair objection.

**The staging below already de-risks this**, which is the main reason to keep
the order. Symmetrica remains the fallback through the "flip the default"
landing, so any platform with no symfn wheel and no cargo degrades to exactly
what Sage does today. The platform question only becomes forcing at the final
step, when Symmetrica is removed — by which point there is real deployment data
to answer it with.

**3. Demotion, not removal — which is what retires the deprecation cycle.**
Sage deprecates rather than deletes, and anything that *loses* functionality
goes through the standard period. Demoting Symmetrica from `type: standard` to
`type: optional` loses none: the 30 entry points Sage never reaches stay
reachable for anyone who installs it. That turns what would have been a
multi-release deprecation of public API into a packaging change.

⚠️ **It no longer carries the coverage bar as well.** This paragraph used to end
"and it is the reason the coverage bar in **1.** above can be met without
reimplementing `scalarproduct_schubert`". That operation is implemented, so the
bar is met directly; demotion now answers only the 30 exports no sagelib code
path reaches, and whether it is the right answer for them is an upstream call.

One consequence that disappears with it: `combinat/schubert_polynomial.py` was
going to be the only sagelib file still importing Symmetrica, needing a
`sage.features` gate and `# optional - symmetrica` doctests for
`scalar_product` alone. With all seven of its calls available, wiring the file
takes it off Symmetrica outright and no gate is needed.

### Staging it inside Sage

The de-risking move is to **not** make displacement a single PR. Three landings,
each independently useful and revertible:

- [ ] **Land as optional.** Every artifact this step names now exists on
      `mwhansen/sage`, branch `symfn`, and none of it has been proposed
      upstream — which is what "land" means, so the box stays open. Built:
      `build/pkgs/symfn/` (`type: optional`), `sage/features/symfn.py` with a
      version floor, 11 doctests tagged `# optional - symfn`,
      `classical.init()` populating the conversion table conditionally at
      import, and `terms.pyx` in Sage's meson build. Verified with the wheel
      installed: `s(h[3,2,1])` and `p(s[2,1])` agree with the Symmetrica answers
      under `SAGE_DISABLE_SYMFN=1`, and the doctests of `combinat/sf/` and
      `combinat/schubert_polynomial.py` pass under `--optional=sage,symfn`.

      **One deviation from this staging plan, and a reviewer will find it.**
      `init()` defaults to `is_available()`, so installing the optional package
      switches the backend immediately — which is the *next* step's behavior
      arriving inside this one. The staging argument was that each landing be
      independently revertible, and "installs but stays off by default" is a
      genuinely different review than "installs and takes over". Either the
      default becomes opt-in for the first PR, or the two steps merge and the
      plan says so. Deciding that is the next real piece of work in this phase.
- [ ] **Flip the default.** symfn becomes the backend when present; Symmetrica
      stays as the fallback. This is the release where the performance claim is
      tested by actual users on actual hardware, and the one that generates the
      evidence for the third step.
- [ ] **Promote symfn to standard and demote Symmetrica to optional.** Only
      after the flip has survived a release in the wild. Not a removal and not a
      deprecation — the 30 entry points Sage never reaches stay reachable
      through the optional package, which is what makes this step a packaging
      change rather than an API break.

**The burden inverts rather than vanishes, and this is the thing to plan
around.** Once Sage depends on symfn, symfn's *Python API* becomes the interface
that must not break — Sage pins a version range, and a change here becomes a
Sage bug. Two consequences:

- It raises the stakes on **Phase 2**. The bulk/indexed entry points the shim
  relies on — partition *indices* rather than lists of parts, which is the whole
  reason it beats Sage's own Symmetrica wrapper — become a contract Sage pins.
  They need to be a deliberate, documented, stable subset of the 108 entry
  points, not whatever happened to be exported. Phase 2 made that list — all of
  them, with `symfn.pyi` and `scripts/check_python_stubs.py` holding it — and
  the docstring half of [policies/python.md](policies/python.md) delta 1 has
  since closed as well, so each of the 108 states its contract and pins its
  convention with an example that runs. What this raises the stakes on now is
  keeping those gates running, which needs the CI Phase 0 still owes.
- Independent release cadence is gone. An adapter fix ships when Sage ships, and
  upstream review is measured in months.

Which is why the order is: **publish the wheel, keep the adapter external, prove
it against two or three consecutive Sage releases, then propose.** Proposing to
*replace a standard package* is the hardest version of this ask, and it is
carried by evidence rather than argument: a track record across releases, the
coverage audit showing nothing is lost, and the 4678 comparisons in
`scripts/check_backend.py` — which are already, precisely, symfn answering
Symmetrica's questions with Sage asking them.

Open questions to resolve before writing any of it, in descending order of risk:

- ~~**What does Sage actually call?**~~ **Answered** —
  [the census](record/python-and-sage-interop.md). 36 of 66
  entry points, six files; all six now intercepted, and all 36 computed by
  symfn — `scalarproduct_schubert`, which this plan had left behind, included.
- **What is Symmetrica's current standing?** How much appetite there is upstream
  for demoting an unmaintained C dependency determines whether this is a welcome
  contribution or an uphill one. Needs checking rather than assuming — it is the
  difference between a receptive review and a dead PR. Note the ask is now
  *demotion to optional*, not retirement, which is the version most likely to
  find agreement.
- ~~**Will Sage take a binary wheel for a standard package?**~~ **Answered:
  yes.** See [docs/sage-packaging-audit.md](sage-packaging-audit.md). A
  maturin-built Rust package (`rpds_py`) is already standard and wheel-only, and
  Sage builds no Rust from source at all.
- **Does the wheel matrix reach every platform Sage supports?** The replacement
  for the toolchain question, and now the open packaging risk. Measured against
  `rpds_py`'s platform set, not against Phase 5's default list.
- Licensing is fine in the direction needed: GPL Sage may depend on an
  MIT/Apache wheel, and the clean-room rule in [NOTICE.md](../NOTICE.md)
  protects against contamination the other way. Note the asymmetry that makes
  displacement easier than it looks: **Symmetrica is public domain**, so if the
  audit turns up something Sage calls and symfn lacks, its algorithms are a
  legitimate source rather than merely a reference — the record already records
  this, and displacement is the scenario it was recorded for.

---

## Phase 6 — the files a contributor needs

- [x] `CONTRIBUTING.md`. Done 2026-08-25, built around the item's point:
      **which checks need what** leads the file — `cargo test` needs nothing
      (the oracle tests read committed fixtures), `preflight.sh` needs
      python3 and no Sage, `preflight_python.sh` needs the python build plus
      ruff/mypy/Sphinx, and regenerating a fixture needs the oracle that
      produced it. Then the once-per-clone git config and a where-things-are
      list (layout, public-api, the record, CLAUDE.md as the rulebook
      router, the release plan). The README's Contributing section shrinks
      to the issue-tracker line and a pointer; the file is scanned by the
      spelling, figures and link gates alongside README.md and CLAUDE.md.
- [ ] `SECURITY.md`, `CODE_OF_CONDUCT.md`, issue and PR templates.
- [x] Restructure the README. Done 2026-08-25, in two rounds. The item's
      original complaint — status and benchmarks before anything runnable —
      had already dissolved: Install and Usage sit directly under the intro.
      What remained heavy was the back half, ~250 lines of contributor and
      governance material inlined where a link would do. The Layout tree
      moved to [layout.md](layout.md) (it had already drifted once, missing
      `_bases.py`), the three API tiers and the 0.x break list moved to
      [public-api.md](public-api.md) with a summary paragraph kept in place,
      and the Validation bullets compressed to one paragraph pointing at
      [policies/validation.md](policies/validation.md). CLAUDE.md's routing
      and the README's Contributing section point at the new files; the
      README went from 407 lines to 252.

**Done when:** someone who has never seen the repo can clone it, run the right
tests, and know which ones they cannot run.

---

## Phase 7 — durability
*Not blocking a release; what keeps it good afterwards.*

- [ ] **Property tests.** Validation today is oracle-and-law based over inputs
      that we or Sage chose. The algebraic laws in `tests/algebra_laws.rs` are
      already written as universally-quantified statements — ω is an involution,
      conversions are ring homomorphisms, Δ is an algebra map — so putting
      `proptest`-generated partitions behind them is nearly free and covers the
      space nobody thought to enumerate. This is the same reasoning that made
      `check_backend.py` (Sage choosing the inputs) find the 200× regression the
      degree ladder never generated.
- [ ] **Benchmarks in a harness.** The speedup figures are the crate's headline
      claim and there is no committed criterion suite to reproduce them or to
      catch a regression. The `scripts/bench_*.py` files measure against Sage;
      what is missing is symfn-against-its-own-history.
- [ ] Coverage reporting, if only to find the paths the oracles never reach.

---

## Phase 8 — the element model
*Closed 2026-08-25.*

The convenience layer's `Sym`/`Param` split — what an element with parameters
is, and why the nine parametric tags should be bases like `s` and `m` — moved
to [plans/element-model.md](plans/element-model.md) on 2026-08-24 and **closed
on 2026-08-25**. There is one element class, `Sym`, carrying all seven
coefficient types and all fifteen bases and picking each operation's route from
`parameters`; `Param` no longer exists.
As of 2026-08-25 that commitment holds for all ten operations in all fifteen
bases, with no refusals left. The last two closed differently: the principal
specialization's gap by adding the operation Sage's message names rather than
by widening a ring, and Jack plethysm by giving `AFrac` a general denominator
factor beside its linear atoms — see [record/jack.md](record/jack.md) for the
measurement that rejected replacing them.

---

## Phase 9 — the code review of 2026-09-03
*What a reader of `src/` alone, with the docs closed, would want changed. The
items are sorted by whether they get more expensive after the first tag.*

The review read the tree without the record or the policies and reported what
the code itself shows. Most of what it found is already covered above and is
not repeated; what follows is the remainder, with the reason each item sits on
its side of the tag. The measurements it took are in
[record/littlewood-richardson.md](record/littlewood-richardson.md) only where
they add to what was there; the rest were spot checks that agreed with the
record.

### Before the first tag

Each of these changes a signature, a name, or a promise. Before the tag they
are free; after it each one is a breaking change or a permanent commitment.

- [ ] **Decide the number line with the review's finding in view.** The crate
      is at 1.0.0-rc.1 and the versioning decision is the open Phase 4 item.
      The finding: the tree has no external caller yet, the root re-exports
      roughly 150 names, and the API tier holds about 400 `pub fn`s. A 1.0 tag
      freezes every one of those under semver on the strength of six weeks of
      in-house use. The alternative is a 0.x first release that gathers callers
      and cuts 1.0 once the surface has held still for a few months. Either
      way the decision is made explicitly, in Phase 4, and
      [public-api.md](public-api.md) says which was chosen and why.
- [ ] **Prune the root re-exports to the entry points a consumer names.**
      Phase 2 sorted the *modules*; it left every module's contents re-exported
      flat at the crate root. The membership test is the same one Phase 2
      used — would a caller who only wants symmetric functions ever name
      it? — applied to the `pub use` list in [lib.rs](../src/lib.rs).
      Whatever stays
      at the root is a promise; whatever moves behind its module path is still
      reachable and can be promoted later without a break. Removing a
      re-export after the tag is a break.
- [ ] **A non-panicking twin for every entry point that panics on overflow.**
      Only `Partition::try_new` and `try_character` exist. `character`, the
      `from_u128`/`from_i128` conversions on the fixed-width rings, and
      `Partition::z` all panic when a value leaves the type, and a Rust caller
      has no way to ask first. The failure policy
      ([policies/failure.md](policies/failure.md)) is satisfied — the panic is
      loud — but the shape of the API is decided here: adding `try_` twins
      later is additive, while changing a return type to `Result` or `Option`
      is not. Pick which entry points get a twin and which change shape, and
      do the shape changes now.
- [ ] **`Partition::z` says the wrong thing about release builds.** Its
      `# Range` section reads "wraps in release and panics in debug". With
      `[profile.release] overflow-checks = true` (Phase 3, R3) it panics in
      both. The sentence contradicts the front page's exactness contract; fix
      it in the same change as the twin above, since the twin is the escape
      the section should point at.
- [ ] **`std::ops` on the basis types.** `Schur`, `PowerSum` and the rest
      have `add`, `sub`, `mul`, `neg`, `scale` as inherent methods and no
      `Add`/`Sub`/`Mul`/`Neg` impls, so `a * b` does not compile and every
      example reads `a.mul(&b)`. Adding the impls is additive, but the
      examples, doctests and README are what callers copy, and the style they
      copy is set by the first release. Implement for references at minimum;
      decide whether by-value impls are wanted at the same time, since adding
      them later changes inference for existing callers.
- [ ] **Bound the caches.** [record/memory.md](record/memory.md) Rule 4
      records that every table in `memo.rs` grows without eviction and that a
      long-running Sage session is the hazard, and names the shape of the
      fix: per-table byte accounting and a budget, not an LRU. The plan is
      [plans/cache-budget.md](plans/cache-budget.md), in four stages —
      accounting and `cache_stats`, the budget and eviction, a long-session
      workload that picks the wheel's default, and the speed check — all
      before 0.9, because the first stage adds entry points.
- [ ] **Read `SKEW_TRACE` once.** `expand_layer` in
      [skew_lr.rs](../src/skew_lr.rs) calls `std::env::var_os` on every
      invocation, which is a syscall and a lock on the hot path of every
      Schur product. The facility is used by the record and stays; read the
      variable into a `OnceLock` at first use. Measure before and after with
      `bench_lr`, because the cost is per call and the calls are short.
- [ ] **Run the whole gate on a machine that is not this one, from the
      tarball.** CI runs the crate suites on three platforms; the Python gate
      and the from-tarball `cargo test` (Phase 4) have run only here. Do both
      on a clean Linux checkout before the tag, because the record's own
      history says the first remote run of any gate finds what the laptop
      cannot
      ([record/python-and-sage-interop.md](record/python-and-sage-interop.md)).
- [ ] **Cut the tag from a clean tree.** The working tree carries an ignored
      `symfn_cy.cpython-314-darwin.so` at the root and `build/`, `dist/`,
      `pybuild/` directories from earlier wheel and sdist runs. None reaches
      the crate (Phase 4's `exclude`) but maturin builds from the working
      tree, and a stale `.so` beside the source is the kind of thing an sdist
      picks up. Delete them, rebuild both artifacts, and diff the file lists
      against the ones Phase 4 and Phase 5 recorded.

### After the first tag

Each of these is internal, additive, or needs users to be worth doing. None
changes a signature.

- [ ] **Property tests and a benchmark harness** — already Phase 7's first
      two items, which the review confirmed: the algebraic laws in
      `tests/algebra_laws.rs` sweep degree ≤ 5 by enumeration, and the
      benches are examples with no committed baseline. Nothing to add beyond
      the confirmation.
- [ ] **Oracle rows above degree 6 for the families other than LR.** The
      lrcalc fixture reaches degree 42 for Schur products and skews; the Sage
      fixture stops at degree 6 for everything it covers. Above that, the
      families are checked by agreement between in-house engines, which
      [policies/validation.md](policies/validation.md) accepts but which does
      not satisfy its own "does not share its mathematics" clause for the
      range the rustdoc advertises. A handful of rows at degree 10–15 for
      characters, Kostka numbers, Hall–Littlewood and Macdonald, generated
      with `SAGE_DISABLE_SYMFN=1`, closes that. After the tag because it
      changes no interface and because the generator's floors are a record
      matter ([policies/validation.md](policies/validation.md), V4).
- [ ] **Collapse the per-ring quadruplication in `python.rs`.** The file is
      9,400 lines and 340 functions. `omega`, `antipode`, `skew_by` and
      `multiply` each exist four times, once per parametric ring, with the
      same body modulo the parse and dump helpers; the `out_of_schur!` and
      `into_schur!` macros show the pattern that would absorb them. The
      Python contract does not move — every entry point keeps its name and
      its plain-data shape ([policies/python.md](policies/python.md)) — so
      this is a refactor behind a frozen surface, and the doctest and stub
      gates in `preflight_python.sh` are the proof it changed nothing.
- [ ] **A bridge to `num-traits`, or an implementation of `Ring` for its
      types.** `Ring` is a twelve-method trait of this crate's own, so a Rust
      caller with an existing coefficient type must implement it by hand.
      A blanket impl over `num_traits::{Zero, One, Signed}` behind the
      `bignum` feature (which already depends on `num-traits`) would let most
      types in for free. Additive; after the tag because the trait's method
      list should be frozen first, and a method added to `Ring` breaks
      external implementors ([public-api.md](public-api.md)).
- [ ] **Allocation in the core types.** `Partition` is a `Vec<u32>` used as
      the key of every `BTreeMap`, and `src/` has about 570 `.clone()`
      calls. A small-vector key, or interning, is a measured job with
      `heapstat` attached and the memory record's checklist followed; it is
      the kind of change Rule 3 in [record/memory.md](record/memory.md) says
      has been reverted twice when done on instinct.
- [ ] **Publish the ×-Sage figures for the families where Sage dispatches to
      its own Python**, in the release notes rather than the rustdoc
      ([style.md](style.md), genre rule). Those are the cases where the
      speedup is an order of magnitude and where a Sage user decides whether
      to install the wheel.
- [ ] **Read the first month of issues before touching the surface again.**
      The reports from real Sage sessions will say which of the entry points
      anyone calls, which is the list a later 1.0 (if the first release is
      0.x) or a 2.0 (if it is not) should be built around.

**Done when:** every item in "Before the first tag" is checked or has a
recorded reason not to be, and the tag is cut from a tree that CI has built
from the tarball.

---

### Dependency order

`Phase 0` (CI) → `Phase 1` (docs render) and `Phase 2` (API surface) in
parallel → `Phase 3` (failure contract, needs 2's surface decided) →
`Phase 4` (crate) and `Phase 5` (wheel) → `Phase 5b` (Sage adapter, needs 5's
package to depend on) → `Phase 6` → `Phase 9`'s "Before the first tag" list,
which is the last thing before the tag because its items are the ones that
stop being free once it exists. `Phase 5c` (upstreaming) trails 5b by several
Sage releases, on purpose. `Phase 7` and `Phase 9`'s "After the first tag" list
are continuous.

The shortest path to something publishable is 0 → 1 → 2 → 4. Phase 5 is
independent of the crate release and can be pulled forward if Python users come
first. Phase 5b is the only phase Sage appears in at all, and it is deliberately
last of the packaging work: the wheel must be shippable and useful with the
adapter never written.
