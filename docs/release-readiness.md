# Release readiness — from "works on this machine" to a package

[The record](record/) tracks what the library *computes*, and how it came to.
This file is the forward-looking half: what stands between the current tree and
something a stranger can depend on — continuous integration, a curated API
surface, a stated contract for failure,
and two publishable artifacts — a crate and a **Sage-free** wheel — with Sage
interoperability layered on top of the wheel rather than baked into it.

The mathematics is not the gap. The suite passes across every feature
combination, the oracles are committed, the licensing is clean and audited, and
`v0.1.0` is already tagged. What is missing is the operational layer — and one
fact frames the whole document:

> **Nothing in this repository has ever been built or tested anywhere except one
> macOS arm64 machine running rustc 1.96.** There is no `.github/`. Every claim
> in the README is, today, a claim about one laptop.

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
- [ ] `cargo doc --no-deps --all-features` gated with `-D warnings` — *after*
      Phase 1 clears the existing 173.
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
- [ ] A separate, non-blocking job for the Sage-dependent checks. **38 of the 40
      scripts in `scripts/` import Sage**, so they cannot run on a normal
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

`cargo package` warns: *manifest has no documentation, homepage or repository*.
`repository` is literally `""`.

- [ ] Fill in `repository`, `homepage`, `documentation`, `readme`, and
      `rust-version`.
- [ ] Add `exclude` so the published tarball is the library. Today it carries
      `scripts/` (40 files, most of them Sage harnesses), all of `docs/`, and
      both oracle fixtures — `lrcalc_oracle.txt`, and `sage_oracle.txt`, which
      covers the `(q,t)` layer, Kronecker, LLT and Schubert as well as the
      classical one. They are the largest thing in the tarball. The fixtures
      should stay if `cargo test` on a published crate is meant to work —
      decide that explicitly rather than by default.
- [ ] `cargo publish --dry-run`, and verify the docs.rs build with the right
      feature set (`all-features` will try to build PyO3; configure
      `[package.metadata.docs.rs]` with `features = ["bignum"]` instead).
- [ ] Add `cargo-deny` to CI. The "every dependency is permissive, so the wheel
      carries no copyleft obligation" claim in `NOTICE.md` is what the whole
      licensing story depends on, and nothing currently stops a future
      dependency from quietly breaking it.
- [ ] `CHANGELOG.md`, starting from the already-tagged `v0.1.0`.

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
`scripts/sage_backend.py` and `scripts/symfn_cy.pyx`, which `cimport`s Sage's
`Integer` type and therefore cannot build without Sage present. The work here is
to make that separation enforced and packaged rather than incidental.

- [x] `pyproject.toml` with `[build-system] requires = ["maturin>=1.5,<2.0"]`
      and a `[project]` table: description, README, license, classifiers,
      `requires-python = ">=3.9"` (matching the `abi3-py39` build).
      **`[project.urls]` is deliberately absent** — `Cargo.toml`'s `repository`
      is empty too and the tree has no published home yet; both get filled in
      together so they cannot disagree.

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
- [ ] Wheel matrix via `cibuildwheel` or `maturin-action`. The abi3 build means
      one wheel per *platform* covers 3.9+, so the matrix is platform-only and
      cheap — and GitHub's arm64 runners make the aarch64 legs native rather
      than emulated. **Size it against `rpds_py`**, the standard Sage package
      symfn would be joining (see
      [docs/sage-packaging-audit.md](sage-packaging-audit.md)): macOS x86_64 and
      arm64; manylinux x86_64, aarch64, armv7l, ppc64le, s390x, i686;
      musllinux x86_64, aarch64, i686; Windows win32, amd64, arm64. abi3 turns
      that into ~13 artifacts rather than the 55 `rpds_py` needs. The exotic
      Linux arches need QEMU legs; decide which are Tier 1 and which are Tier 2
      rather than dropping them silently.
- [ ] **A support-tier policy, written down.** *Tier 1* — a prebuilt wheel
      exists, `pip install symfn` needs no toolchain. *Tier 2* — no wheel;
      builds from the sdist, needs cargo. This is the sentence Phase 5c's
      review will turn on, so it should exist before then and be honest about
      which platforms are which.
- [ ] **An sdist that builds offline.** `cargo vendor` the `python` feature's
      dependencies (PyO3, num-bigint, num-rational, num-traits) into the sdist.
      The default build already has zero dependencies and builds offline by
      design; the wheel build does not, and offline source builds are exactly
      the configuration distro packagers use. Test it in CI with the network
      off.
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
- [ ] A Python test suite that runs **without Sage** — round-trip the marshalling
      layer against values computed in Rust. `check_bindings.py` already tests
      the boundary rather than the library, which is the right idea; it just
      needs a Sage-free sibling that CI can run on a stock runner.
      `scripts/check_python_boundary.py` is the first such sibling and covers
      the *failure* half: 128 malformed calls over 94 pyfunctions, asserting
      only typed exceptions come back. `scripts/check_convenience.py` is the
      second, and holds the convenience layer to the contract layer over 2177
      checks. **The round-trip half is still open**, and it is a marshalling
      test rather than an oracle: a value handed in comes back out intact, at
      the widths and shapes `docs/policies/python.md` P1 promises. Catching a
      kernel defect is not this surface's job — `docs/policies/validation.md`
      owns that, and the Rust suites and `tests/fixtures/` discharge it.
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
- [ ] A publish-on-tag workflow for PyPI, with trusted publishing.
- [x] **The rendered reference**, at `docsite/`, published by Read the Docs.
      Sphinx with MyST, building the wheel first so both layers are documented
      from the objects themselves and `help()` cannot drift from the website.
      `scripts/check_docs_complete.py` compares `symfn.__all__` against the
      inventory Sphinx writes and fails when a supported name reaches no page:
      Sphinx warns about references that do not resolve and says nothing about
      an entry point nobody wrote a directive for, which is the hole that opens
      as the surface grows. All 108 contract entry points and 79 convenience
      names are covered.

**Done when:** `pip install symfn` works on Linux, macOS and Windows without a
Rust toolchain, and the package imports and computes with no Sage anywhere.

---

## Phase 5b — the Sage extension module, on the far side of the boundary

Sage integration stays possible, but as a **separate artifact** that depends on
`symfn` rather than the other way round. Today it exists as two files in
`scripts/` that were written to run experiments, not to be installed by anyone:
`sage_backend.py` (monkey-patches
`sage.combinat.sf.classical.conversion_functions` in a live session) and
`symfn_cy.pyx` (the compiled per-term loop, which needs Sage's headers to
build).

That is the right architecture already — it just needs to become a thing with a
name, rather than a script that assumes `sys.path.insert(0, "pybuild")`.

- [ ] Decide the shape. Three options, and the third is the intended
      destination:
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
- [ ] Either way, lift `sage_backend.py` out of `scripts/`: remove the hardcoded
      `sys.path.insert(0, "pybuild")`, give it a real entry point
      (`symfn_sage.install()` rather than import-time patching), and let it
      locate `symfn` as an ordinary installed package.
- [ ] Document the two modes it supports, because they are genuinely different
      products: **backend replacement** (drop into Sage's conversion table and
      accelerate everything Sage already does) and **direct use** (call symfn
      for the things Sage has no equivalent for — reduced Kronecker via the `st`
      basis, the Macdonald operator algebra, the (q,t)-Kostka tables).

### Becoming a complete drop-in

Today's `sage_backend.py` displaces **one** of Symmetrica's six consumer sites —
the conversion table in `combinat/sf/classical.py`. That is why it can claim
4678 comparisons and still not be a replacement. The other five call sites reach
Symmetrica directly, and
[docs/symmetrica-coverage-audit.md](symmetrica-coverage-audit.md) enumerates
them. **These are ordinary library and binding tasks — they need no upstream
involvement and can be done at any point**, which is why they belong here rather
than in Phase 5c:

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
      **a gap the audit missed**, and **decided: won't do.** It returns a
      Schubert *polynomial*; `schubert_pairing`, which the audit paired it with,
      returns an integer. Nothing in sagelib calls
      `SchubertPolynomial.scalar_product` — only its own definition and its own
      doctests — so it is public API with no internal dependents, and Symmetrica
      staying available as an optional package covers it. The other six of the
      seven Schubert calls are verified exact over 1679 comparisons and can be
      wired when the file is taken on; see the audit for the two requirements
      that adds (Sage's doctests assert `ValueError`s symfn answers instead, and
      Symmetrica blocks on an interactive prompt at teardown after
      `divdiff_perm_schubert`).
- [ ] Keep the Cython shim optional. It is a measured 185 ns/term win on the
      per-term loop, but it needs Sage's headers and a working Cython; the pure
      Python path must stay a working fallback.
- [ ] The Sage-dependent CI job from Phase 0 is what tests all of this, and it
      is the only place Sage ever appears in the build graph.

**Done when:** a Sage user installs `symfn` from PyPI, installs or points at the
adapter, and gets both modes — with the adapter's absence costing the wheel
nothing.

---

## Phase 5c — upstreaming the adapter into Sage
*The intended end state, and deliberately not the first move.*

The adapter's natural long-term home is **inside the Sage codebase**, for one
concrete reason: `symfn_cy.pyx` `cimport`s Sage's `Integer` and therefore must
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

- [x] **Do the audit first.** Done —
      [docs/symmetrica-coverage-audit.md](symmetrica-coverage-audit.md).
      **Sage reaches 36 of the 66 exported entry points, from six files**, and
      the whole representation-theory half of Symmetrica — the part that would
      have been a research programme — **is never called by Sage at all**. Its
      per-entry-point gap list was then corrected by implementing it: 35 of the
      36 are covered today and 29 are intercepted, with
      `scalarproduct_schubert` the one operation symfn does not have.

The coverage requirement that remains is small, and **it is inherited from
Phase 5b rather than new here** — the tasks under
[Becoming a complete drop-in](#becoming-a-complete-drop-in) are what close it,
and they need no upstream involvement. Five of the six consumer sites are
intercepted; the sixth, `schubert_polynomial.py`, waits on
`scalarproduct_schubert`. By the time this phase opens that should be done too,
and the adapter proven across releases.

What is genuinely new at this phase is a policy question, not code:

- [ ] Propose **demoting Symmetrica from standard to optional**, rather than
      removing it. That is the answer to the **31 entry points that stay
      behind** — the 30 unreached ones plus `scalarproduct_schubert` — in one
      move: they keep working for anyone who installs the optional package, so
      none of them needs reimplementing *or* deprecating. **Decided here**; it
      is a materially softer ask than removal, and it is what lets Phase 5b
      close with a gap it deliberately did not fill.

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
`type: optional` loses none: the 31 entry points symfn will not displace stay
reachable for anyone who installs it. That turns what would have been a
multi-release deprecation of public API into a packaging change, and it is the
reason the coverage bar in **1.** above can be met without reimplementing
`scalarproduct_schubert`.

One consequence to plan for: `combinat/schubert_polynomial.py` becomes the only
sagelib file still importing Symmetrica, so it needs a `sage.features` gate and
`# optional - symmetrica` doctests — the same mechanism the *first* landing uses
for symfn, pointed the other way. Wiring its other six calls (verified exact,
see the audit) narrows that gate from the whole file to `scalar_product` alone,
which is the argument for doing it before this step rather than never.

### Staging it inside Sage

The de-risking move is to **not** make displacement a single PR. Three landings,
each independently useful and revertible:

- [ ] **Land as optional.** `build/pkgs/symfn/`, a `sage.features` gate so Sage
      builds and runs fine without it, doctests tagged `# optional - symfn`, and
      the conversion table in `sage/combinat/sf/classical.py` populated
      conditionally at import rather than monkey-patched by `sage_backend.py`.
      The Cython shim moves into Sage's build here.
- [ ] **Flip the default.** symfn becomes the backend when present; Symmetrica
      stays as the fallback. This is the release where the performance claim is
      tested by actual users on actual hardware, and the one that generates the
      evidence for the third step.
- [ ] **Promote symfn to standard and demote Symmetrica to optional.** Only
      after the flip has survived a release in the wild. Not a removal and not a
      deprecation — the 31 entry points symfn does not displace stay reachable
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
  [docs/symmetrica-coverage-audit.md](symmetrica-coverage-audit.md). 36 of 66
  entry points, six files; five files intercepted, and
  `scalarproduct_schubert` deliberately left to the optional package.
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

- [ ] `CONTRIBUTING.md`. The single most valuable thing in it: **which checks
      need what.** `cargo test` needs nothing. `tests/lrcalc_oracle.rs` and
      `tests/sage_oracle.rs` run against committed fixtures and also need
      nothing — but *regenerating* them needs lrcalc or Sage. 38 of the 40
      scripts need Sage. `CLAUDE.md` now states that map for agent sessions;
      this file is where it reaches human contributors, for whom it is the
      first question.
- [ ] `SECURITY.md`, `CODE_OF_CONDUCT.md`, issue and PR templates.
- [ ] Restructure the README. It currently spends ~60 lines on status and
      benchmarks before anything a reader can run. Lead with `cargo add symfn`
      and a five-line example; move the achievement narrative below the fold or
      into `record/README.md`, which is where that story already lives in full.

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

### Dependency order

`Phase 0` (CI) → `Phase 1` (docs render) and `Phase 2` (API surface) in
parallel → `Phase 3` (failure contract, needs 2's surface decided) →
`Phase 4` (crate) and `Phase 5` (wheel) → `Phase 5b` (Sage adapter, needs 5's
package to depend on) → `Phase 6`. `Phase 5c` (upstreaming) trails 5b by
several Sage releases, on purpose. `Phase 7` is continuous.

The shortest path to something publishable is 0 → 1 → 2 → 4. Phase 5 is
independent of the crate release and can be pulled forward if Python users come
first. Phase 5b is the only phase Sage appears in at all, and it is deliberately
last of the packaging work: the wheel must be shippable and useful with the
adapter never written.
