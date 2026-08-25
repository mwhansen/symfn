# symfn

[![CI](https://github.com/mwhansen/symfn/actions/workflows/ci.yml/badge.svg)](https://github.com/mwhansen/symfn/actions/workflows/ci.yml)
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

## Usage

```rust
use symfn::{Partition, Schur, SymFn};

// s_2 · s_1 = s_3 + s_{21}
let s2: Schur<i64> = Schur::monomial(Partition::new([2]), 1);
let s1: Schur<i64> = Schur::monomial(Partition::new([1]), 1);
let prod = s2.mul(&s1);
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
| (q,t)-Kostka matrix, degree 12 | ~18× | Sage |
| `st[4,3] · st[4,3]` (reduced Kronecker) | 3400× | Sage (its own Python) |
| The whole substitution, through Sage | 1.84× | Symmetrica, like-for-like |

The last row is the one to read twice. End to end through Sage the
substitution is **1.84× like-for-like**; the 4.37× also quotable there includes
a Sage-`Partition` cache that Symmetrica's wrapper does not have and could
equally adopt. The gap between the in-crate ratios and the end-to-end one is
object marshalling, which both backends pay.

A ratio against a timeout is not a measurement, so cases Sage does not finish
are reported as times rather than ratios.

## Features

| feature | what it adds |
|---|---|
| *(default)* | the whole library over `i64`/`i128` and an exact `Rational`; no dependencies |
| `bignum` | `BigInt` / `BigRational` coefficients (`num-bigint`, pure Rust) — exact beyond `i128` |
| `python` | PyO3 extension module (abi3, CPython 3.9+); the wheel `pyproject.toml` builds adds the pure-Python convenience layer on top |

## The public API, and what a version number promises

While the major version is 0, **the minor number is the breaking one**: 0.1 →
0.2 may change or remove anything public, and a patch release will not. The
public module list is decided rather than accumulated, in three tiers — the
API, the `#[doc(hidden)]` cross-check engines and product strategies that are
promised nothing, and the private core — and what consumers build on is the
coefficient-ring layer: `Ring`, with the `QAlgebra` and `Plethystic`
refinements above it. The Python surface is frozen harder than the crate,
because it is the contract nearly every consumer reaches this library
through. The tier lists, and the two breaks that do not look like breaks, are
in [docs/public-api.md](docs/public-api.md).

## Validation

The check that matters most: `scripts/check_backend.py` has an incumbent
computer algebra system drive symfn as its own backend and compares the
answers on every input the incumbent's dispatch reaches, so the inputs are
its choice rather than ours — which is what found a 200x regression on shape
families the degree ladder never generated. Under it sit committed oracle
fixtures from two independent external programs (each regenerable by script,
so an auditor can check rather than trust), exhaustive agreement between the
three Littlewood–Richardson backends, and law suites whose two sides share no
code: conversions are ring homomorphisms, ω is an involutive algebra map, the
Hopf axioms hold, and skewing matches its defining adjunction. The Python
layer is checked against a live oracle separately, because a correct answer
marshalled into the wrong slot is a different failure from a wrong answer.
[docs/policies/validation.md](docs/policies/validation.md) states what
evidence a new family owes before it ships.

## Building from source

The Rust side needs no dependencies, no network, and no external oracle — the
oracle tests read committed fixtures under `tests/fixtures/`.

```sh
cargo test                      # core suite
cargo test --features bignum    # + arbitrary-precision coefficients
scripts/preflight.sh            # the commit gate: fmt check + both suites
cargo doc --open                # the reference
```

The Python extension module needs maturin:

```sh
maturin build --release --features python
mkdir -p pybuild && unzip -q -o target/wheels/*.whl -d pybuild
PYTHONPATH=pybuild python -c "import symfn; print(symfn.schur_multiply([([2],1)],[([1],1)]))"
```

## Contributing

Bugs and questions go to
[the issue tracker](https://github.com/mwhansen/symfn/issues).

Read [CLAUDE.md](CLAUDE.md) first: it routes to the five rulebooks that govern
prose, failure handling, the Python surface, validation, and the record.
[docs/layout.md](docs/layout.md) maps every module in the tree. Before
working in a subsystem, read its file in
[docs/record/](docs/record/) — dead ends are recorded with their premises
exactly so they are not re-explored at full price.

Once per clone:

```sh
git config core.hooksPath .githooks
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

`scripts/preflight.sh` is the gate for a Rust change;
`scripts/preflight_python.sh` is the gate for anything under `src/python.rs`,
`python/symfn/` or `docsite/` — and for an edit to this file's Python
example, which runs there. The release and packaging scripts —
`build_sdist.sh`, `check_sdist_offline.sh`, `build_docs.sh` — are documented in
[scripts/README.md](scripts/README.md).
[docs/release-readiness.md](docs/release-readiness.md) is the release plan,
and [docs/plans/](docs/plans/) holds the design plans that are not release
gates.

## License

**MIT OR Apache-2.0** — see [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).

symfn contains no third-party code. It is validated against two GPL programs,
named in [NOTICE.md](NOTICE.md), by invoking them as external oracles to
generate committed test fixtures. The LR engine was written **clean-room** —
specification and
implementation by separate parties, the implementer having no access to
`lrcalc` — with the spec committed at
[docs/cleanroom-spec-skew-lr.md](docs/cleanroom-spec-skew-lr.md) as the audit
trail. Every dependency, optional ones included, is permissively licensed, so
the published wheel carries no copyleft obligation. [NOTICE.md](NOTICE.md) has
the details.
