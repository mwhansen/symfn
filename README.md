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

```
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

```python
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

**The classical core.** All six bases {m, e, h, p, s, f} as distinct types,
with conversions between every ordered pair. Three products: ordinary
(Littlewood–Richardson), plethysm, and the internal/Kronecker product. Full
Hopf structure — coproduct, counit, antipode, and skewing by an *arbitrary*
symmetric function. Evaluation at a finite alphabet, the principal
specializations and their q-analogue, symmetric-group characters and Kostka
numbers, as single values and as whole tables.

**The parametrized families.** Hall–Littlewood `Q'_λ(x;t)` by the Morris
recursion and `P_λ` by inverting the Kostka–Foulkes matrix; Kostka–Foulkes
`K_{λμ}(t)` from that transition. Macdonald `P_λ`, `Q_λ` and `J_λ` by the
branching formula, over a `ℚ(q,t)` that avoids bivariate gcd by keeping
denominators factored. The (q,t)-Kostka polynomials `K_{λμ}(q,t)` and the
modified form `H̃_μ = Σ_λ K̃_{λμ}(q,t) s_λ`, computed by the Bergeron–Haiman
Pieri recursion with the branching formula and a Lapointe–Lascoux–Morse
eigenvector solve kept as independent cross-checks — three algorithms sharing
nothing above `Partition`. Jack `P/Q/J_λ(x;α)` with Laplace–Beltrami as the
engine. LLT polynomials in both the ribbon and tuple models, the tuple
entries straight or skew. The deformed Hall pairings `⟨,⟩_t`, `⟨,⟩_{q,t}` and
`⟨,⟩_α` beside the classical one — each family is orthogonal under its own.

Every one of those runs **backwards** too — an element rewritten *into* `P`,
`Q`, `Q'`, `J` or `H̃` rather than expanded out of one, which is the direction
a positivity question asks in, and the direction that had no entry point until
the transitions were inverted. Each is a back-substitution through the forward
expansion of its own degree, memoized, so a sweep over a degree costs one
solve rather than p(n). The forward direction takes a whole element too, so
each basis is a place an element can be written rather than a table it is read
out of, and the Python surface names shapes in it: `jack.P([2])` is `JackP[2]`,
and `.to` converts among all fifteen basis codes, classical and parametric
alike, the power-sum expansion included. LLT is the exception and cannot be otherwise: its
polynomials are linearly dependent across the level `k`, so they are not a
basis of Λ.

**Reduced Kronecker coefficients** as an outer product, in the
Orellana–Zabrocki bases `s̃_λ` and `h̃_λ`. The calculation never leaves the
power-sum basis, so no Littlewood–Richardson coefficient enters it at all.
`st[6,4] · st[6,4]` takes 0.14s; Sage does not finish it.

**Single structure constants for products that cannot be materialized.**
`schubert::schubert_coeff` answers `c^w_{uv}` by Bruhat pruning for pairs whose
product no machine holds, and `ops::kronecker_coeff` answers one `g^ν_{λμ}`
where the whole internal product does not fit. Also the Matchings–Jack and
b-conjecture coefficients, which no other package computes.

Coefficients are generic over a `Ring`; the paths that divide ask only for a
`QAlgebra` (a ring containing ℚ), so ℚ[t] and ℚ[q,t] qualify even though
neither is a field. Arbitrary precision is automatic at the Python boundary: a
call runs in fixed width and re-runs exactly if anything overflows.

## Performance

Every figure is in [the record](docs/record/) with its harness and power
state, including the ones that went the wrong way.

| Computation | Ratio | Measured against |
|---|---|---|
| Classical-basis conversions | 2–10× | Symmetrica (C) |
| Hall–Littlewood `Q'`, end to end | 2.5–3.4× | Symmetrica `hall_littlewood` (C) |
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

## Design in one screen

| Decision | What we did | Why |
|---|---|---|
| No untyped object | One **type per basis** (`Schur`, `PowerSum`, `Monomial`) behind the `SymFn` trait | Basis confusion becomes a compile error, not a runtime bug (vs. Symmetrica's `OP`) |
| Coefficients | Generic over `Ring`; dividing paths bounded on `QAlgebra`, not `Field` | Every division is by z_μ, an *integer* — so ℚ[t] and ℚ[q,t] qualify, which is what Macdonald/Hall–Littlewood need |
| Littlewood–Richardson | Behind the `LrBackend` trait; three native backends, `SkewLr` the default (no external C lib) | The trait paid off: each new backend swapped in with no caller changes and is cross-checked against the previous ones |
| Partitions | `Partition` newtype, invariant enforced at construction | Weakly-decreasing/positive guaranteed, not merely assumed |
| Correctness | Known-value + algebraic-law tests, committed oracle fixtures, and Sage driving symfn as its own backend | An oracle only tests inputs you thought of; letting Sage pick them found a 200x regression the benchmark could not see |

## Features

| feature | what it adds |
|---|---|
| *(default)* | the whole library over `i64`/`i128` and an exact `Rational`; no dependencies |
| `bignum` | `BigInt` / `BigRational` coefficients (`num-bigint`, pure Rust) — exact beyond `i128` |
| `python` | PyO3 extension module (abi3, CPython 3.9+); the wheel `pyproject.toml` builds adds the pure-Python convenience layer on top |

## The public API, and what a version number promises

The [Layout](#layout) maps the tree; this maps the *interface*, which is
smaller. The public module list is decided rather than accumulated, and the
test that sorts it is whether a caller who only wants symmetric functions would
ever name the module. Three tiers, checkable with `cargo doc --no-deps`:

- **API.** Every family and every conversion, plus the types and coefficient
  rings they are written over: `sym`, `convert`, `ops`, `hopf`, `plethysm`,
  `eval`, `kostka`, `character`, `character_basis`, `charge`, `kf`, `hl`,
  `jack`, `macdonald`, `qtkostka`, `llt`, `schubert`, `lr`, `skew_lr`,
  `deltaop`, `dyck`, `gj`, `partition`, `permutation`, `coeff`, `guard`, `qt`,
  `frac`, `afrac`. Every item in them carries a doc comment —
  `#![deny(missing_docs)]` is what keeps that true rather than a habit.
- **Reachable, promised nothing** — all `#[doc(hidden)]`. `bh`, `gjmod` and
  `macop` are cross-check engines that exist to disagree with a primary route;
  `rect`, `two_row`, `three_row` and `strip_lr` are product strategies whose
  *results* are API even though their paths are not (`okada_coeff`,
  `okada_product`, `two_row_coeff`, `two_row_product`, `three_row_product`,
  `AutoLr` and `StripLr` are re-exported at the crate root and documented
  there); `measure` is heap accounting and `python` is the PyO3 bridge, and
  neither is symmetric functions. Naming one of these compiles. It is not a
  promise, and it is where a break lands without a version bump.
- **Private.** `memo`, `modular`, `fasthash`, `candidates`. `clear_caches` is
  re-exported at the crate root, because timing a run means clearing the
  caches between them.

**What consumers build on is the coefficient-ring layer**: `Ring`, and the
`QAlgebra` and `Plethystic` refinements above it. Generic code bounded on those
three is what survives a basis or backend being rewritten underneath it — the
bound is deliberately weaker than `Field` so that ℚ[t] and ℚ[q,t] qualify,
which is what Macdonald and Hall–Littlewood need. `LrBackend` is the same shape
one level down: three native backends implement it, and swapping one for
another changed no caller.

### 0.x

While the major version is 0, **the minor number is the breaking one**: 0.1 →
0.2 may change or remove anything in the API tier, and a patch release will
not. The hidden tier is outside that guarantee entirely — it can move in a
patch, which is why it is a separate tier rather than a naming convention.

Two things break API-tier callers that do not look like breaks, so they are
worth naming:

- **A method added to `Ring`, `SymFn`, `LrBackend` or `SkewBy`** breaks any
  code implementing the trait outside this crate, while breaking no caller.
  The coefficient-ring traits are the seams this library changes behavior
  through ([docs/policies/failure.md](docs/policies/failure.md)), so they are
  the ones most likely to gain a method.
- **The Python surface is frozen harder than the crate**, not in step with it.
  It is the contract nearly every consumer reaches this library through, and
  [docs/policies/python.md](docs/policies/python.md) is its rulebook; a
  crate-internal change is free, and the same change at that boundary is not.

## Validation

- **Sage as the driver, not the oracle** — `scripts/check_backend.py` runs Sage
  once on each backend and compares the two answers on every input Sage's own
  dispatch reaches, covering Hall–Littlewood, Jack and Macdonald as well as the
  classical bases, `expand`, the monomial product and semistandard tableaux.
  This is the check that matters most: every other test uses inputs *we* chose,
  so it can only find bugs we thought of. Letting Sage pick them found a 200x
  regression on shape families the degree ladder never generated.
- **Sage oracle** — `tests/sage_oracle.rs` checks values Sage computed
  independently (Kostka numbers as tableau counts, characters, Schur products,
  all four conversions out of Schur, skew Schur). The fixture is committed and
  `scripts/gen_sage_oracle.sage` regenerates it, so an auditor can check rather
  than trust.
- **lrcalc oracle** — `tests/lrcalc_oracle.rs` checks products and skew
  expansions against `lrcalc`, on shapes far larger than Sage can finish. That
  size is the point: the single-traversal backend produces a whole *set* of
  terms at once, and a bug there shows up as a missing or extra term rather
  than a wrong single coefficient. Fixture is committed.
- **Backend agreement** — all three LR backends are checked against each other
  exhaustively over every product with |μ|+|ν| ≤ 7, and `expand_skew` against
  the coefficient-at-a-time path over every skew shape up to degree 8.
- **Algebraic laws** — `tests/algebra_laws.rs`: conversions are ring
  homomorphisms, round-trips are the identity, ω is an involutive algebra map,
  Hall pairings ⟨s,s⟩/⟨h,m⟩/⟨p,p⟩ are correct, Δ is an algebra map.
- **Hopf axioms** — the antipode satisfies `m∘(S⊗id)∘Δ = ε·1`, and skewing is
  checked against its *definition*, ⟨g⊥f, h⟩ = ⟨f, g·h⟩, for every triple
  through degree 6 — a test naming no algorithm, whose two sides share no code.
- **The Python layer separately** — `scripts/check_bindings.py` checks the
  bindings against Sage rather than against a dump, because a correct answer
  marshalled into the wrong slot is a different failure from a wrong answer,
  and only one of the two shows up in a dump.

[docs/policies/validation.md](docs/policies/validation.md) states what evidence
a new family owes before it ships.

## Building from source

The Rust side needs no dependencies, no network, and no Sage — the oracle tests
read committed fixtures under `tests/fixtures/`.

```
cargo test                      # core suite
cargo test --features bignum    # + arbitrary-precision coefficients
scripts/preflight.sh            # the commit gate: fmt check + both suites
cargo doc --open                # the reference
```

The Python extension module needs maturin:

```
maturin build --release --features python
mkdir -p pybuild && unzip -q -o target/wheels/*.whl -d pybuild
PYTHONPATH=pybuild python -c "import symfn; print(symfn.schur_multiply([([2],1)],[([1],1)]))"
```

## Contributing

Bugs and questions go to
[the issue tracker](https://github.com/mwhansen/symfn/issues).

Read [CLAUDE.md](CLAUDE.md) first: it routes to the five rulebooks that govern
prose, failure handling, the Python surface, validation, and the record.
Before working in a subsystem, read its file in
[docs/record/](docs/record/) — dead ends are recorded with their premises
exactly so they are not re-explored at full price.

Once per clone:

```
git config core.hooksPath .githooks
git config blame.ignoreRevsFile .git-blame-ignore-revs
```

`scripts/preflight.sh` is the gate for a Rust change;
`scripts/preflight_python.sh` is the gate for anything under `src/python.rs`,
`python/symfn/` or `docsite/`. The release and packaging scripts —
`build_sdist.sh`, `check_sdist_offline.sh`, `build_docs.sh` — are documented in
[scripts/README.md](scripts/README.md).
[docs/release-readiness.md](docs/release-readiness.md) is the release plan,
and [docs/plans/](docs/plans/) holds the design plans that are not release
gates.

## Layout

```
src/
  coeff.rs      Ring / Field traits, i64+i128 rings, exact Rational field
  partition.rs  Partition newtype, conjugate, z(λ), partition generator
  lr.rs         LrBackend trait + native NaiveLr (incremental pruning)
  strip_lr.rs   StripLr row-strip DP; AutoLr, the backend the library uses
  skew_lr.rs    SkewLr — whole-shape expansion, merged row layer (default)
  two_row.rs    s_μ·s_ν with a two-row factor, counting fibres per output
  three_row.rs  the three-row analogue of the same counting route
  candidates.rs the candidate walk both counting routes share, and its workers
  rect.rs       Okada's closed form for a product of two rectangles
  kostka.rs     Kostka numbers K_{λμ}: SSYT counting, and enumeration
  qt.rs         ℤ[q,t] / ℚ[q,t] coefficients — sparse, sorted, merge-accumulated
  hl.rs         Hall–Littlewood Q'_λ(x;t) by the Morris recursion, and P
  kf.rs         Kostka–Foulkes K_{λμ}(t) from the Hall–Littlewood transition
  frac.rs       ℚ(q,t) with denominators kept factored — no bivariate gcd
  macdonald.rs  Macdonald P/Q/J_λ(x;q,t) by the branching formula
  bh.rs         modified Macdonald H̃_μ by the Bergeron–Haiman Pieri recursion
  qtkostka.rs   the (q,t)-Kostka polynomials K_{λμ}(q,t); three routes kept
  macop.rs      the Macdonald operator M₁ as a matrix on modified Schurs
  deltaop.rs    the operator algebra: ∇, Δ_f, Δ'_f, Π and Θ_f
  dyck.rs       labeled Dyck paths; the Delta conjecture's combinatorial side
  llt.rs        LLT polynomials — ribbon and tuple models, three engines
  jack.rs       Jack P/Q/J_λ(x;α); Laplace–Beltrami is the engine
  afrac.rs      ℚ(α) as factored integer-linear atoms, kept canonical
  gj.rs         the Goulden–Jackson connection-coefficient pipeline c^λ_{μν}(b)
  gjmod.rs      the same tables by modular evaluation; engines_agree
  modular.rs    prime fields, CRT, rational reconstruction, interpolation
  permutation.rs  Perm, permutations moving finitely many points (Schubert)
  schubert.rs   Schubert polynomials S_w; coefficients of products too large
                to materialize
  charge.rs     the charge statistic; K_{λμ}(t) by tableau enumeration (reference)
  character.rs  χ^λ(μ) via Murnaghan–Nakayama (β-number rim hooks)
  character_basis.rs  the OZ bases s̃/h̃; reduced Kronecker via the power-sum route
  sym.rs        SymFn / SymAlgebra traits; all six bases; multiplication
  convert.rs    ToSchur / FromSchur hub; Jacobi–Trudi; Muir's rule; h↔e flip
  ops.rs        ω involution, Hall inner product, internal (Kronecker) product
  hopf.rs       SymTensor, skew Schur, SkewBy, coproduct, counit, antipode
  plethysm.rs   f[g] through the power-sum basis
  eval.rs       evaluation at an alphabet; principal specializations; dim λ
  guard.rs      overflow-reporting coefficients + the escalation scope
  measure/      heap accounting shared by benchmarks, budget tests, heapstat
  fasthash.rs   the DP layers' hasher; memo.rs  the caches
  python.rs     the PyO3 bridge; lib.rs  crate docs and re-exports
python/symfn/  the wheel's pure-Python half — the convenience layer
  __init__.py   the package: the contract layer re-exported flat, then this
  _bases.py     the basis codes; denominator clearing, so rational
                coefficients cross the integer contract layer exactly
  _sym.py       Sym, basis-tagged; the factories s, h, e, p, m, f; skew
  _param.py     Poly, QtPoly, QtFrac, QtRatio, AlphaFrac: the coefficient
                types that carry a parameter
  _families.py  the namespaces macdonald, jack, hl, llt
  _schubert.py  Schub over permutations, and the factory X
  _types.py     the type vocabulary the layer is annotated in
  symfn.pyi     the contract surface as a list, held to the module by
                scripts/check_python_stubs.py
  py.typed      so a checker reads both halves
docsite/   the rendered reference (Sphinx + MyST), published by Read the Docs
tests/
  oracle.rs        known Schur expansions + commutativity/associativity/degree
  algebra_laws.rs  ring-hom conversions, ω algebra map, Hall pairings, Δ algebra map
  sage_oracle.rs   Sage-computed values, from a committed fixture
  lrcalc_oracle.rs products and skew expansions vs lrcalc, past Sage's sizes
  qalgebra.rs      the library over ℚ[t] — a ring that is deliberately not a Field
  bignum.rs        exactness past i128
  memory.rs        peak-bytes and allocation budgets over the measure workloads
  fixtures/        the committed oracle outputs both *_oracle suites read
docs/
  style.md             the prose rulebook, for every documentation surface
  sage-backend.md      standing in for Symmetrica under Sage, in one file
  policies/            failure, the Python surface, validation
  record/              the memory — one file per subsystem; README.md indexes
examples/  research drivers
  *_dump.rs        emitters whose output scripts/check_*.py hold to Sage
  bench_*, profile_*, probe_*  per-subsystem instruments of the record
  delta_conjecture.rs, find_nonzero.rs  conjecture checks
scripts/   nearly all need Sage; scripts/README.md documents the main ones
  preflight.sh      the local gate: fmt check + both test suites, no Sage
  check_backend.py  A/B the two backends through Sage itself; the adapter it
                    drives is not here — docs/sage-backend.md says where
  check_bindings.py the Python layer itself against Sage, not a dump
  preflight_python.sh  the Python gate: stubs, typed exceptions, both layers'
                    docstring examples, the convenience layer against the
                    contract layer, ruff, mypy --strict, docs completeness
  check_python_boundary.py, check_python_stubs.py, check_python_docs.py,
  check_convenience.py, check_convenience_docs.py, check_docs_complete.py
                    the ones needing no Sage; preflight_python.sh runs them
  check_*.py        one Sage oracle per subsystem (hl, kf, macdonald, jack,
                    llt, qt_kostka, deltaop, eval, skew, st, …)
  bench_*.py        the Sage side of each ladder, same work on both sides
  spec_*.py         pre-implementation verification and wall measurement
  gen_*, compare_*  fixture generators; direct lrcalc/Symmetrica comparisons
```

## License

**MIT OR Apache-2.0** — see [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).

symfn contains no third-party code. It is validated against two GPL programs
(Sage and `lrcalc`) by invoking them as external oracles to generate committed
test fixtures. The LR engine was written **clean-room** — specification and
implementation by separate parties, the implementer having no access to
`lrcalc` — with the spec committed at
[docs/cleanroom-spec-skew-lr.md](docs/cleanroom-spec-skew-lr.md) as the audit
trail. Every dependency, optional ones included, is permissively licensed, so
the published wheel carries no copyleft obligation. [NOTICE.md](NOTICE.md) has
the details.
