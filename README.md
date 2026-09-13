# symfn

[![CI](https://github.com/mwhansen/symfn/actions/workflows/ci.yml/badge.svg)](https://github.com/mwhansen/symfn/actions/workflows/ci.yml)
[![Sage oracle](https://github.com/mwhansen/symfn/actions/workflows/sage.yml/badge.svg)](https://github.com/mwhansen/symfn/actions/workflows/sage.yml)
[![docs](https://img.shields.io/badge/docs-symfn.readthedocs.io-blue)](https://symfn.readthedocs.io)
[![crates.io](https://img.shields.io/crates/v/symfn.svg)](https://crates.io/crates/symfn)
[![PyPI](https://img.shields.io/pypi/v/symfn.svg)](https://pypi.org/project/symfn/)
[![docs.rs](https://img.shields.io/docsrs/symfn)](https://docs.rs/symfn)
[![license](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](#license)

A Rust and Python library for computing with **symmetric functions**
and other related objects: the six classical bases, the
Hall–Littlewood, Macdonald, LLT and Jack families above them, and
Schubert polynomials.

Each basis is a distinct type over a coefficient ring the caller
chooses. Every value that leaves the library is exact, or the call
fails loudly.

## Install

```sh
# Rust
cargo add symfn

# Python
pip install symfn
```

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

The Python wheel has two supported layers. The convenience layer is best for
interactive use as it provides things like human-readable representations, etc.

```pycon
>>> from symfn import s, h
>>> s([2, 1]) * s([1])
s[2,1,1] + s[2,2] + s[3,1]
>>> h([2]).to("s")
s[2]
```

The contract layer returns results in plain Python types, most useful for building
a library on top of the `symfn` kernel.

```pycon
>>> import symfn
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
  their q-analog, symmetric-group characters and Kostka numbers, as single
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
  and `g^ν_{λμ}` for products too large to compute.
- **Reduced Kronecker coefficients** as an outer product in the
  Orellana–Zabrocki bases, and the Matchings–Jack and b-conjecture
  coefficients.
- ... and more.

Coefficients are generic over the ring, and arbitrary precision is automatic
at the Python boundary: a call that runs in fixed width in Rust will re-run
with a wider coefficient ring if anything overflows.

## Validation

Every public family of symmetric functions is checked with the
strongest class of evidence available to it, and with at least one
check that does not share its mathematics
([docs/policies/validation.md](docs/policies/validation.md) is the
policy). The five classes, strongest first:

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

## Platforms

Rust 1.87 or later; CPython 3.9 or later. Any platform without a
prebuilt wheel builds from the source distribution and needs a Rust
toolchain.

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

A single wheel serves every CPython from 3.9 up.

**Tier 2 — no wheel; the source distribution builds it.** Everything
else: FreeBSD, OpenBSD, illumos, Linux on riscv64 or mips, macOS back
past the wheel's minimum, and any platform Rust supports.  The
requirement is a Rust toolchain at or above `Cargo.toml`'s
`rust-version` plus a C linker.


### Building from source

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

### Rust Features

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
