# symfn

[![CI](https://github.com/mwhansen/symfn/actions/workflows/ci.yml/badge.svg)](https://github.com/mwhansen/symfn/actions/workflows/ci.yml)
[![Sage oracle](https://github.com/mwhansen/symfn/actions/workflows/sage.yml/badge.svg)](https://github.com/mwhansen/symfn/actions/workflows/sage.yml)
[![docs](https://img.shields.io/badge/docs-symfn.readthedocs.io-blue)](https://symfn.readthedocs.io)
[![crates.io](https://img.shields.io/crates/v/symfn.svg)](https://crates.io/crates/symfn)
[![PyPI](https://img.shields.io/pypi/v/symfn.svg)](https://pypi.org/project/symfn/)
[![docs.rs](https://img.shields.io/docsrs/symfn)](https://docs.rs/symfn)
[![license](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](#license)

A Rust and Python library for computing with **symmetric functions**
and other related objects: the six classical bases and every
transition between them, the Hall–Littlewood, Macdonald, LLT and Jack
families above them, and Schubert polynomials.

Each basis is a distinct type over a coefficient ring the caller chooses, so a
basis mix-up is a compile error and `ℚ[q,t]` is as ordinary a coefficient ring
as `ℤ`. Every value that leaves the library is exact, or the call fails
loudly. The default build has no dependencies.

## Install

```sh
# Rust
cargo add symfn

# Python
pip install symfn
```

Rust 1.87 or later; CPython 3.9 or later. Any platform without a
prebuilt wheel builds from the source distribution and needs a Rust
toolchain — [Building from source](#building-from-source) has the
local build instructions.

### Platforms

**Tier 1 — a prebuilt wheel exists.** `pip install symfn` needs no compiler, no
Rust toolchain and no network beyond the index. Every release carries one wheel
for each platform below, built by
[.github/workflows/release.yml](.github/workflows/release.yml).

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
CPython from 3.9 up. Which platforms the set has to reach is decided by what
shipping under Sage requires — [Under Sage](#under-sage) states that
requirement.

⚠️ **A cross-compiled wheel is built but never imported on the platform it
targets**, because the runner cannot execute it. `abi3` is what makes cross
builds possible at all — no interpreter for the target is needed — and also
what keeps the gap narrow, since the ABI linked against is the stable one
rather than a version's internals. The native legs run an import-and-compute
step in the release workflow; the cross legs do not, so one could reach a user
first. That is accepted exposure, not an oversight.

**Tier 2 — no wheel; the source distribution builds it.** Everything else:
FreeBSD, OpenBSD, illumos, Linux on riscv64 or mips, macOS back past the
wheel's minimum, and any platform Rust supports that the table does not name.
The requirement is a Rust toolchain at or above `Cargo.toml`'s `rust-version`
plus a C linker, and nothing more — the sdist carries the `python` feature's
four dependencies vendored, so the build needs no network.
`scripts/check_sdist_offline.sh` asserts that on every push. It is also the
tier the source-building distributions live in by choice: Debian, conda-forge,
Gentoo and nix build from source whatever wheel exists.

Nothing is unsupported; the line between the tiers is whether someone else
already compiled it for you. Adding a platform is a matrix entry in the release
workflow and a row above. Removing one is a breaking change for whoever was
installing without a toolchain, so it belongs in a release note with the
reason, not in a quiet matrix edit. Raising `rust-version` moves nothing
between tiers but narrows Tier 2, since it is the floor a source build clears.

## Usage

```rust
use symfn::{Partition, Schur, SymFn};

// s_2 · s_1 = s_3 + s_{21}
let s2: Schur<i64> = Schur::monomial(Partition::new([2]), 1);
let s1: Schur<i64> = Schur::monomial(Partition::new([1]), 1);
let prod = &s2 * &s1;
assert_eq!(prod.coeff(&Partition::new([3])), 1);
assert_eq!(prod.coeff(&Partition::new([2, 1])), 1);
```

The Python wheel has two supported layers. The convenience layer is for code a
person reads; the contract layer is plain data, for marshalling in bulk or
building another library on top.

```pycon
>>> import symfn
>>> from symfn import s, h
>>> s([2, 1]) * s([1])
s[2,1,1] + s[2,2] + s[3,1]
>>> h([2]).to("s")
s[2]
>>> symfn.schur_multiply([([2, 1], 1)], [([1], 1)])
[((2, 1, 1), 1), ((2, 2), 1), ((3, 1), 1)]
```

Full reference: [docs.rs](https://docs.rs/symfn) for the crate,
[symfn.readthedocs.io](https://symfn.readthedocs.io) for the Python surface.

## What it computes

- **The six classical bases** {m, e, h, p, s, f}, conversions between every
  ordered pair, and three products: Littlewood–Richardson, plethysm, and the
  internal (Kronecker) product.
- **The full Hopf structure** — coproduct, antipode, skewing by an arbitrary
  symmetric function — plus evaluation, the principal specializations and
  their q-analogue, symmetric-group characters and Kostka numbers, as single
  values and as whole tables.
- **Hall–Littlewood, Macdonald and Jack** in their `P`/`Q`/`J`
  normalizations, the Kostka–Foulkes and (q,t)-Kostka polynomials, the
  modified basis `H̃`, LLT polynomials in the ribbon and tuple models, and
  the delta-operator tower ∇, Δ'_f, Θ.
- **Every family backwards too**: an element rewritten *into* `P`, `Q`,
  `Q'`, `J` or `H̃` — the direction a positivity question asks in — with
  `to` converting among all fifteen basis codes, classical and parametric
  alike.
- **The deformed Hall pairings** `⟨,⟩_t`, `⟨,⟩_{q,t}` and `⟨,⟩_α` beside
  the classical one; each family is orthogonal under its own.
- **Schubert polynomials**, including single structure constants `c^w_{uv}`
  and `g^ν_{λμ}` for products too large to materialize.
- **Reduced Kronecker coefficients** as an outer product in the
  Orellana–Zabrocki bases, and the Matchings–Jack and b-conjecture
  coefficients, which no other package computes.

Coefficients are generic over the ring, and arbitrary precision is automatic
at the Python boundary: a call runs in fixed width and re-runs exactly if
anything overflows.

## Under Sage

Sage can compute its symmetric functions and Schubert polynomials through symfn
instead of Symmetrica. This section is the whole statement of that
relationship; the measurements and the reasoning behind each binding are in
[docs/record/python-and-sage-interop.md](docs/record/python-and-sage-interop.md),
and what is left to do is [docs/todo-1.0.md](docs/todo-1.0.md).

**symfn covers every Symmetrica entry point Sage calls.** That is narrower than
"symfn replaces Symmetrica", and the difference is not rhetorical. Sage reaches
roughly half of what `sage/libs/symmetrica/all.py` exports, from the six files
below, and symfn computes all of it. The rest are public API a user reaches by
writing `from sage.libs.symmetrica.all import ...`, which no sagelib code path
touches — mostly the representation-theory half, the decomposition matrices and
wreath-product and Gupta tables. symfn does not cover those, and retiring them
is a deprecation question for Sage rather than an implementation question here.

| Sage file | what reaches symfn |
| --- | --- |
| `combinat/sf/classical.py` | the basis conversions |
| `combinat/sf/sfa.py` | the `compute_*_with_alphabet` family |
| `combinat/sf/monomial.py` | the monomial product |
| `combinat/sf/hall_littlewood.py` | the `P` and `Q'` caches |
| `combinat/tableau.py` | Kostka numbers, semistandard tableaux |
| `combinat/schubert_polynomial.py` | product, `multiply_variable`, `expand`, polynomial to Schubert, the scalar product, divided differences |

**The adapter is not in this repository.** The wheel contains no Sage code,
imports no Sage module, and gains nothing from Sage being present — a test
asserts it. The adapter lives on the Sage side as `sage/libs/symfn/`:
`backend.py`, `extras.py`, and the compiled per-term loop `terms.pyx`. It is on
`mwhansen/sage` branch `symfn`, which has not been proposed upstream. So
installing symfn into a stock Sage changes nothing by itself; what changes
Sage's behavior is that branch, where `sage.combinat.sf.classical.init()` fills
the conversion table at import if `sage.features.symfn` finds the package, and
where the five other consumer files ask `is_available()` before dispatching.
`build/pkgs/symfn/SPKG.rst` states the two modes on the Sage side: **backend
replacement**, which accelerates what Sage already does, and **direct use**,
which covers what Sage has no equivalent for at all — reduced Kronecker
coefficients through the `st` basis, the Macdonald operator algebra, the
(q,t)-Kostka tables.

### Turning it off

`SAGE_DISABLE_SYMFN=1` makes Sage answer out of Symmetrica and its own Python.
It has to be in the process *environment* before Sage starts: the conversion
table is filled when `sage.combinat.sf.classical` is imported, so setting it
inside a running session is too late.

⚠️ **Everything here that uses Sage as an oracle, or as a benchmark's control
arm, depends on that — and every one of them fails quietly rather than loudly
if it is missing.** An A/B whose control arm does not set it compares symfn to
symfn: it passes, it proves nothing, and the symptom is a run of ratios near
1.0×. `scripts/gen_sage_oracle.sage` refuses to run without it, because Sage
would otherwise quote this library back into its own fixture. Because the
failure is silent, the rule is not left to the invoker to remember:
`scripts/sage_guard.py` states it and its `require_own_sage` refuses to run
without the variable, and `scripts/check_sage_guards.py` — a static scan inside
`scripts/preflight.sh` — fails when a Sage-importing script under `scripts/`
neither calls the guard nor names itself, with a reason, as one that measures
the backend on purpose.

## Validation

Every public family is checked with the strongest class of evidence
available to it, and with at least one check that does not share its
mathematics ([docs/policies/validation.md](docs/policies/validation.md) is
the policy). The five classes, strongest first:

- **External oracles**: independent implementations written by others,
  used two ways: committed fixtures that `cargo test` re-checks on every
  run, and live harnesses where the oracle picks the inputs.
- **Independent in-tree routes**: a second engine that shares no code, and
  preferably no mathematics, with the first, so agreement between the two
  is evidence rather than consistency.
- **Identities and specializations**: published identities tying a new
  family to one that already has an oracle, each chosen to catch an error
  the other checks would miss.
- **Algebraic laws**: ring homomorphisms, round trips, the Hopf axioms.
  Cheap and broad, but every convention satisfies the same laws, so a law
  can never say which convention this is.
- **Convention pins**: hand-checkable values chosen to differ between this
  library's normalization and other conventions appearing in the
  literature.

## Performance

Every figure is in [the record](docs/record/) with its harness and power
state, including the ones that went the wrong way.

| Computation | Ratio | Measured against |
|---|---|---|
| Classical-basis conversions | 2–8× | Symmetrica (C, via Sage) |
| Littlewood–Richardson products | 1–30×, ahead on every above-floor case | lrcalc (C) |
| Schubert products | 2–30× | schubmult (C) |
| Hall–Littlewood `Q'`, end to end | 2.5–3.4× | Symmetrica `hall_littlewood` (C, via Sage) |
| Kostka–Foulkes `K_{λμ}(t)`, per pair | 30–40× | Sage `kfpoly` (its own Python) |
| Kostka–Foulkes, a whole column at once | 880× | Sage `kfpoly` (its own Python) |
| Macdonald `P/Q/J`, through degree 9 | ~94× | Sage (its own Python) |
| (q,t)-Kostka matrix, degree 12 | ~18× | Sage (its own Python) |
| `∇e_n`, degrees 8 to 13 | 16–27× | Sage (its own Python) |
| Jack `m → P`, `m → Q`, `m → J`, degree 10 | 960–1090× | Sage (its own Python) |
| LLT, whole-degree tables to `k = 4` | 700–19 100× | Sage (its own Python) |
| A Kronecker coefficient, degrees 12 to 24 | 40–2350× | Sage `itensor` (its own Python) |
| `st[4,3] · st[4,3]` (reduced Kronecker) | 3400× | Sage (its own Python) |
| Plethysm | 25–40× | Sage (its own Python) |
| `sage.combinat.sf` tests | 1.84× | Symmetrica, like-for-like |

## Building from source

The Rust side needs no dependencies, no network, and no external oracle — the
oracle tests read committed fixtures under `tests/fixtures/`.

```sh
cargo test                      # core suite
CARGO_TARGET_DIR=target/bignum cargo test --features bignum
scripts/preflight.sh            # the commit gate: fmt check + both suites
cargo doc --open                # the reference
```

Each feature set wants its own target directory. Changing the set
invalidates everything built under a shared one, and the rebuild is minutes
where the same suite warm in its own directory is under one; the gate scripts
do this for themselves.

The Python extension module needs maturin:

```sh
CARGO_TARGET_DIR=target/python maturin build --release --features python
mkdir -p pybuild && unzip -q -o target/python/wheels/*.whl -d pybuild
PYTHONPATH=pybuild python -c "import symfn; print(symfn.schur_multiply([([2],1)],[([1],1)]))"
```

### Features

| feature | what it adds |
|---|---|
| *(default)* | the whole library over `i64`/`i128` and an exact `Rational`; no dependencies |
| `bignum` | `BigInt` / `BigRational` coefficients (`num-bigint`, pure Rust) — exact beyond `i128` |
| `python` | PyO3 extension module (abi3, CPython 3.9+); the wheel `pyproject.toml` builds adds the pure-Python convenience layer on top |

## Contributing

Bugs and questions go to
[the issue tracker](https://github.com/mwhansen/symfn/issues).
[CONTRIBUTING.md](CONTRIBUTING.md) has the contributor setup: which checks
need what, the commit gates, and where the rulebooks and the record live.
[SECURITY.md](SECURITY.md) is the private route for a vulnerability, and
[CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md) applies to every discussion of the
project.

## License

**MIT OR Apache-2.0** — see [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).

symfn contains no third-party code. It is validated against two GPL programs,
named in [NOTICE.md](NOTICE.md), by invoking them as external oracles to
generate committed test fixtures. The LR engine was written **clean-room** —
specification and implementation by separate parties, the specification
carrying no implementation technique and the implementer having no access to
`lrcalc`. Every dependency, optional ones included, is permissively licensed,
so the published wheel carries no copyleft obligation. [NOTICE.md](NOTICE.md)
has the details.
