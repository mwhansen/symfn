# Todo for 1.0

What is left after 0.9.0. The plan that tracked the whole path from "works on
this machine" to a package closed everything it set out to — CI, the API sort,
the failure contract, the crate and wheel artifacts, the Sage adapter, the
2026-09-03 code review — and what it measured on the way is in
[record/](record/). What remains is here.

Nothing below changes a published signature. The Python surface freezes hardest
of anything in the tree ([policies/python.md](policies/python.md)), and it
freezes harder still once Sage depends on it.

## The Sage landing

The adapter is `sage/libs/symfn/` on `mwhansen/sage`, branch `symfn`. What it
covers and how to turn it off is the README's "Under Sage" section; the
measurements are in
[record/python-and-sage-interop.md](record/python-and-sage-interop.md). Nothing
below has been proposed upstream, which is what "land" means.

- [ ] **A CI job that tests the adapter.** `.github/workflows/sage.yml` exists
      and Sage appears nowhere else in the build graph, but the Sage it
      installs is conda-forge's, which has no adapter — so it runs the oracle
      scripts and not `scripts/check_backend.py`. The adapter is tested where
      it lives, by Sage's own doctests on the branch. This closes when a Sage
      that carries the adapter is installable on a runner.
- [ ] **Decide, upstream, what happens to the 30 entry points Sage never
      reaches.** symfn covers every Symmetrica entry point Sage itself calls —
      36 of the 66 Symmetrica exports, per the census in
      [record/python-and-sage-interop.md](record/python-and-sage-interop.md).
      The other 30 are reachable only by a user who writes `from
      sage.libs.symmetrica.all import ...`. Reimplement them, deprecate them
      through Sage's normal cycle, or keep Symmetrica installable as an
      optional package for whoever imports it directly. That decision belongs
      to Sage. It is an ordinary deprecation question on its own merits: an
      earlier version of this item bundled it with proposing Symmetrica's
      demotion, which was only worth doing while `scalarproduct_schubert` had
      no replacement, and it now has one.
- [ ] **Land as optional.** Built and verified with the wheel installed:
      `build/pkgs/symfn/` (`type: optional`), `sage/features/symfn.py` with a
      version floor, 11 doctests tagged `# optional - symfn`,
      `classical.init()` populating the conversion table conditionally at
      import, and `terms.pyx` in Sage's meson build. `s(h[3,2,1])` and
      `p(s[2,1])` agree with the Symmetrica answers under
      `SAGE_DISABLE_SYMFN=1`, and the doctests of `combinat/sf/` and
      `combinat/schubert_polynomial.py` pass under `--optional=sage,symfn`.

      ⚠️ **One deviation a reviewer will find.** `init()` defaults to
      `is_available()`, so installing the optional package switches the backend
      immediately — the next step's behavior arriving inside this one. The
      staging argument is that each landing be independently revertible, and
      "installs but stays off by default" is a genuinely different review from
      "installs and takes over". Either the default becomes opt-in for the
      first PR, or the two steps merge and this file says so. Deciding that is
      the next real piece of work here.
- [ ] **Flip the default.** symfn becomes the backend when present; Symmetrica
      stays as the fallback. This is the release where the performance claim is
      tested by actual users on actual hardware, and the one that generates the
      evidence for the step after it.
- [ ] **Promote symfn to standard and demote Symmetrica to optional.** Only
      after the flip has survived a release in the wild. Not a removal and not
      a deprecation — the 30 entry points Sage never reaches stay reachable
      through the optional package, which is what makes this a packaging change
      rather than an API break.

**The interim is safe.** Until the last step, Symmetrica remains the fallback,
so a platform with no symfn wheel and no cargo gets exactly what Sage gives it
today. Platform reach is the real risk in the end — Symmetrica is C and
compiles anywhere, a wheel reaches only the platforms someone built for, and
"fewer platforms than the package you are displacing" is a concrete review
objection — which is why the wheels cover exactly `rpds_py`'s platform set
(the README's Install section lists them). It only becomes forcing at the
promotion.

**The burden inverts rather than vanishes, and it is the thing to plan
around.** Once Sage depends on symfn, symfn's *Python API* is the interface
that must not break: Sage pins a version range, and a change here becomes a
Sage bug. The bulk and indexed entry points the shim relies on — partition
indices rather than lists of parts, which is why it beats Sage's own
Symmetrica wrapper — are part of that contract. Independent release cadence
goes with it: an adapter fix ships when Sage ships.

## Release infrastructure

What 0.9.0 left unfinished outside the tree. The account of the release is in
[record/python-and-sage-interop.md](record/python-and-sage-interop.md).

- [ ] **Create the Read the Docs project.** `symfn.readthedocs.io` returns 404,
      and the published PyPI metadata's Documentation link and the release
      notes both point there. The project slug must be exactly `symfn`, pull
      request builds should be enabled, and `stable` should follow `v0.9.0`.
      Creating it needs the maintainer's GitHub sign-in. This closes when
      `symfn.readthedocs.io/en/stable/` renders 0.9.0.
- [ ] **Register the crates.io trusted publisher, and revoke the 0.9.0 token.**
      0.9.0 was published by hand with an API token, because crates.io has no
      pending publishers. On the crate's settings: owner `mwhansen`, repository
      `symfn`, workflow `release.yml`, environment `crates-io`, then
      trusted-publishing-only mode. Until then the `crates` job in
      `.github/workflows/release.yml` cannot publish 0.9.1. This closes when
      the crate's settings list the publisher and no API token exists.

## Internal work

Additive or invisible from outside, and none of it changes a signature.

- [ ] **Collapse the per-ring quadruplication in `python.rs`.** The file is
      9,400 lines and 340 functions. `omega`, `antipode`, `skew_by` and
      `multiply` each exist four times, once per parametric ring, with the same
      body modulo the parse and dump helpers; the `out_of_schur!` and
      `into_schur!` macros show the pattern that would absorb them. The Python
      contract does not move — every entry point keeps its name and its
      plain-data shape ([policies/python.md](policies/python.md)) — so this is
      a refactor behind a frozen surface, and the doctest and stub gates in
      `preflight_python.sh` are the proof it changed nothing.
- [ ] **A bridge to `num-traits`, or an implementation of `Ring` for its
      types.** `Ring` is a twelve-method trait of this crate's own, so a Rust
      caller with an existing coefficient type must implement it by hand. A
      blanket impl over `num_traits::{Zero, One, Signed}` behind the `bignum`
      feature (which already depends on `num-traits`) would let most types in
      for free. It waits on the trait's method list being frozen, because a
      method added to `Ring` breaks external implementors
      ([public-api.md](public-api.md)).
- [ ] **Allocation in the core types.** `Partition` is a `Vec<u32>` used as the
      key of every `BTreeMap`, and `src/` has about 570 `.clone()` calls. A
      small-vector key, or interning, is a measured job with `heapstat`
      attached and the memory record's checklist followed; it is the kind of
      change Rule 3 in [record/memory.md](record/memory.md) says has been
      reverted twice when done on instinct.
- [ ] **Read the first month of issues before touching the surface again.**
      Reports from real Sage sessions will say which entry points anyone
      actually calls, and that is the list 1.0 should be built around.
