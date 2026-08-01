# Can a standard Sage package be a prebuilt wheel? — an audit

*Audit of `sagemath/sage` at `09472ff` (10.10.beta7, 2026-07-26), 441 packages
in `build/pkgs/`, cross-checked against the SageMath 10.9 installation on this
machine.*

[docs/release-readiness.md](release-readiness.md) Phase 5c called this "the
single highest-value fact in this phase": if symfn is to displace Symmetrica it
must become a **standard** package, and the question was whether Sage would
accept one satisfied by prebuilt binary wheels — or insist it be buildable from
source in Sage's own tree, which would put a Rust toolchain in Sage's build
path.

**Answer: yes, and there is a direct precedent that is a Rust extension module
built with maturin.** The concern was misplaced.

---

## 1. Sage's package types

| type | count |
|---|---|
| standard | 272 |
| optional | 157 |
| experimental | 8 |
| base | 1 |

## 2. Standard packages are routinely prebuilt wheels

**131 of the 272 standard packages** (48%) declare a `.whl` as their upstream
artifact in `checksums.ini` rather than a source tarball — `sympy`, `networkx`,
`jupyterlab`, `sphinx`, `pip`, `setuptools`, and so on. Sage downloads the wheel
and installs it; there is no build.

This is a first-class mechanism, not a workaround. `build/sage_bootstrap/package.py`
documents multiple platform-tagged artifacts per package in its own docstring:

```
Supports multiple tarballs with format:
tarball=package-VERSION-cp311-cp311-manylinux_2_17_x86_64.whl
sha256=abc123...
tarball=package-VERSION-cp311-cp311-macosx_11_0_arm64.whl
```

and `build/bin/sage-spkg` has explicit `*.whl)` branches for install.

Of those 131, only **one** is platform-specific rather than
`py3-none-any` — and it is the case that matters.

## 3. `rpds_py` — the precedent

`build/pkgs/rpds_py/` is **`type: standard`**, and it is:

- a **Rust** extension module (bindings to the `rpds` crate — its `SPKG.rst`
  says so outright);
- **built with maturin** — verified locally, the installed
  `rpds_py-0.30.0.dist-info/WHEEL` in Sage 10.9 reads
  `Generator: maturin (1.10.2)`, `Root-Is-Purelib: false`;
- shipped as **55 prebuilt platform wheels** plus one sdist, covering
  macOS x86_64/arm64, manylinux x86_64/aarch64/armv7l/ppc64le/s390x/i686,
  musllinux x86_64/aarch64/i686, and Windows win32/amd64/arm64, for CPython
  3.12, 3.13, 3.14 and 3.14t;
- installed in the Sage on this machine as a compiled
  `rpds/rpds.cpython-314-darwin.so`, shipping `.pyi` stubs, importable from
  `sage -python`.

Its entire `build/pkgs/rpds_py/` directory is five files — `SPKG.rst`,
`checksums.ini`, `dependencies`, `package-version.txt`, `type`. **There is no
`spkg-install`**, because nothing is built.

That is architecturally the same object symfn would be.

## 4. `clarabel` — the closer packaging match

`build/pkgs/clarabel/` is `type: optional`, also Rust, and ships
**`cp39-abi3`** wheels — the *exact* build configuration symfn already uses
(`pyo3` with `abi3-py39`). Five artifacts cover every platform and every Python
version.

The contrast with `rpds_py` is the useful part: **abi3 collapses 55 artifacts
into 5**, because one wheel per platform serves all Python versions instead of
one per (platform × Python). symfn's existing `abi3-py39` choice is already the
right one, and it makes the artifact matrix roughly a tenth the size of the
precedent that is a standard package.

## 5. Sage never builds Rust from source

The finding that removes the objection outright:

- **There is no `rust` or `cargo` package in `build/pkgs/`.** Sage ships no Rust
  toolchain and has no way to acquire one.
- **No `spkg-install` in any of the 441 packages invokes `cargo`.** The strings
  "Rust"/"cargo" appear only in the *prose descriptions* of `rpds_py/SPKG.rst`
  and `clarabel/SPKG.rst`.

So Sage's answer to "how do we handle a Rust dependency" is already settled, and
it is: **consume prebuilt wheels, never compile it**. A symfn spkg would not be
introducing a new burden or a new policy — it would be the third package in a
pattern Sage already relies on for a *standard* dependency.

## 6. What this corrects in the plan

- ❌ *"It puts a Rust toolchain in Sage's build path"* — **wrong for Sage.** The
  toolchain question does not arise upstream at all. It remains true only for
  downstream distro packagers (Debian, conda-forge, Gentoo, nix), who build from
  source on principle. That is a real but much narrower constraint, and it is
  what the offline `cargo vendor` sdist in Phase 5 is for.
- ❌ *"Expect the `cryptography` objection"* — **substantially weaker than
  assumed.** That argument was about a package that forced Rust into build
  pipelines. Sage has already resolved this by not building Rust at all, and has
  accepted a maturin-built Rust package as **standard**.
- ✅ The Phase 5 wheel matrix is the right shape, but should be **widened**.
  `rpds_py` is the bar for a standard package, and it covers armv7l, ppc64le,
  s390x, i686 and Windows arm64 — well beyond the five or six platforms Phase 5
  currently lists. `cibuildwheel` reaches most of these; the exotic Linux arches
  need QEMU legs. Worth matching before proposing, because "fewer platforms than
  the package you are displacing" is a concrete review objection, and the
  displacement target — Symmetrica — is C and therefore builds anywhere.

**That last point is the one genuine risk this audit surfaces**, and it is a
narrower and more tractable one than the toolchain concern it replaces: not
*"will Sage accept a wheel"* (it will), but *"does the wheel reach every platform
Sage supports"*. Symmetrica compiles from source on any platform with a C
compiler; a wheel only reaches platforms someone built for. The staged rollout
in Phase 5c already handles the interim — Symmetrica remains the fallback until
the final landing — so this only becomes forcing at removal.

## 7. Symmetrica's position, for reference

- `type: standard`, C, **public domain**, version 3.1.0, sourced from
  `gitlab.com/sagemath/symmetrica` (Sage maintains its own modernized fork —
  upstream's author died in 2013).
- It has an `spkg-configure.m4`, so Sage can use a system copy if one is
  present and suitably patched.
- **Exactly one package depends on it: `sagelib`.** The dependency surface is
  Sage's own library and nothing else.
- Its `SPKG.rst` scope statement is worth reading against
  [docs/record/README.md](record/README.md)'s "Beyond the core" list, because
  they are nearly the same list: ordinary representation theory of the
  symmetric group and related groups, ordinary representation theory of the
  classical groups, modular and projective representation theory of the
  symmetric group, combinatorics of tableaux, symmetric functions and
  polynomials, commutative and non-commutative Schubert polynomials,
  operations of finite groups, and Hecke algebras of type Aₙ.

That is the *library's* scope, not the part Sage calls. The Phase 5c coverage
audit — enumerating `sage/libs/symmetrica/` entry points and, separately, the
call sites in Sage that reach them — is still the open question, and it is now
the highest-value unknown in the phase. Nothing here measures it.

---

### Reproducing

```
git clone --filter=blob:none --sparse --depth 1 https://github.com/sagemath/sage
cd sage && git sparse-checkout set build
for d in build/pkgs/*/; do
  grep -q '^tarball=.*\.whl' "$d/checksums.ini" 2>/dev/null &&
    echo "$(cat $d/type) $(basename $d) $(grep -c '^tarball=' $d/checksums.ini)"
done | sort
```
