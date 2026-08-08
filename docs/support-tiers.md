# Support tiers

**Tier 1 — a prebuilt wheel exists.** `pip install symfn` needs no compiler, no
Rust toolchain, and no network beyond the index. The wheel is built by
[.github/workflows/release.yml](../.github/workflows/release.yml) and every
release carries one for each platform below.

**Tier 2 — no wheel; the source distribution builds it.** `pip install symfn`
compiles the crate, so it needs a Rust toolchain at the version floor
`Cargo.toml` states. The sdist carries its dependencies vendored, so the build
needs no network — `scripts/build_sdist.sh` produces it and
`scripts/check_sdist_offline.sh` asserts the property, on every push rather than
at release time.

Nothing is unsupported. The line between the tiers is whether someone else
already compiled it for you.

## Tier 1

| Platform | Wheel tag | Runner |
| --- | --- | --- |
| Linux x86\_64 (glibc) | `manylinux_2_17_x86_64` | native |
| Linux aarch64 (glibc) | `manylinux_2_17_aarch64` | cross-compiled |
| Linux i686 (glibc) | `manylinux_2_17_i686` | cross-compiled |
| Linux armv7l (glibc) | `manylinux_2_17_armv7l` | cross-compiled |
| Linux ppc64le (glibc) | `manylinux_2_17_ppc64le` | cross-compiled |
| Linux s390x (glibc) | `manylinux_2_17_s390x` | cross-compiled |
| Linux x86\_64 (musl) | `musllinux_1_2_x86_64` | cross-compiled |
| Linux aarch64 (musl) | `musllinux_1_2_aarch64` | cross-compiled |
| Linux i686 (musl) | `musllinux_1_2_i686` | cross-compiled |
| macOS x86\_64 | `macosx_10_12_x86_64` | native |
| macOS arm64 | `macosx_11_0_arm64` | native |
| Windows amd64 | `win_amd64` | native |
| Windows win32 | `win32` | cross-compiled |
| Windows arm64 | `win_arm64` | cross-compiled |

Fourteen artifacts, one per platform rather than one per platform × Python
version, because the extension is built `abi3-py39`: a single wheel serves every
CPython from 3.9 up, and the floor is the one `pyproject.toml` advertises.

That list is not chosen for convenience. It is exactly the platform set of
`rpds_py`, the Rust extension module built with maturin that Sage already ships
as a **standard** package ([sage-packaging-audit.md](sage-packaging-audit.md)).
Sage ships no Rust toolchain and no `spkg-install` in its 441 packages invokes
cargo, so its answer to a Rust dependency is to consume prebuilt wheels — which
makes the platform set the whole of the compatibility question, and "fewer
platforms than the package you are displacing" a concrete review objection.
Symmetrica, the C library symfn would displace, compiles anywhere a C compiler
runs.

**What "cross-compiled" costs.** A cross-compiled wheel is built but not
imported on the platform it targets, because the runner cannot execute it.
`abi3` is what makes cross-compiling possible at all — the build needs no
interpreter for the target — and it is also what makes the gap narrow, since
the ABI a cross build links against is the stable one rather than a version's
internals. The native legs run an import-and-compute step in the release
workflow. A cross leg that builds and does not load would reach a user first,
and that is a known and accepted exposure rather than an oversight.

## Tier 2

Everything else: FreeBSD, OpenBSD, illumos, Linux on riscv64 or mips, macOS
back past the wheel's minimum, and any platform Rust supports that this list
does not name.

The requirement is a Rust toolchain at or above `Cargo.toml`'s `rust-version`,
plus a C linker. Nothing else — the crate's default build has no dependencies,
and the wheel build's four (PyO3, num-bigint, num-rational, num-traits) are
vendored into the sdist.

This is also the tier the source-building distributions live in by choice.
Debian, conda-forge, Gentoo and nix build from source on principle whatever
wheel exists, which is why the offline property is a gate rather than a
courtesy.

## Moving a platform between tiers

Adding one is a matrix entry in the release workflow and a row here. Removing
one is a breaking change for whoever was installing without a toolchain, so it
belongs in a release note with the reason, not in a quiet matrix edit.

Raising `rust-version` moves nothing between tiers but narrows Tier 2, since it
is the toolchain floor a source build has to clear —
[release-readiness.md](release-readiness.md) records that as a breaking change
for distro packagers.
