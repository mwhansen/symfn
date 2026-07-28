# symfn

A modern, fast, testable kernel for computing with **symmetric functions** — a
clean-room Rust successor *in spirit* to Symmetrica (it shares no code with the
old C library), designed to interoperate with Sage without paying Python
overhead on the hot path.

## Status

**v0.1 — complete over the classical core, and validated by replacing
Symmetrica inside Sage.**

All six classical bases {m, e, h, p, s, f} as distinct types, with conversions
between every ordered pair. Three products: ordinary (Littlewood–Richardson),
plethysm, and the internal/Kronecker product. Full Hopf structure — coproduct,
counit, antipode, and skewing by an *arbitrary* symmetric function. Evaluation
at a finite alphabet, the principal specializations and their q-analogue,
symmetric-group characters and Kostka numbers (single values and whole tables).

Beyond the classical core: **Hall–Littlewood** `Q'_λ(x;t)` by the Morris
recursion, and **Kostka–Foulkes** `K_{λμ}(t)` from that transition — 2.5–3.4×
Symmetrica's `hall_littlewood` end to end, and 30–40× Sage's `kfpoly` per pair
(880× when a whole column is asked for at once, which is the recursion's natural
unit of work). **Macdonald** `P_λ(x;q,t)` by the branching formula, ~16× Sage,
over a ℚ(q,t) that avoids bivariate gcd by keeping denominators factored.

Coefficients are generic over a `Ring`; the paths that divide ask only for a
`QAlgebra` (a ring containing ℚ), so ℚ[t] and ℚ[q,t] work even though neither is
a field. Arbitrary precision is automatic at the Python boundary: a call runs in
fixed width and re-runs exactly if anything overflows.

Validation is layered — 147 unit and integration tests, algebraic-law suites,
committed fixtures from Sage and `lrcalc`, and **4678 computations driven by
Sage itself** with symfn substituted for Symmetrica as its conversion backend
(`scripts/check_backend.py`), covering Hall–Littlewood, Jack and Macdonald as
well as the classical bases.

Performance, with the caveats that matter: the classical-basis conversions run
**2–10x** Symmetrica's C. End to end *through Sage* the same substitution is
**1.84x** like-for-like, or 4.37x with a Sage-`Partition` cache that Symmetrica's
wrapper does not have and could equally adopt. The gap between those is object
marshalling, which both backends pay. `ROADMAP.md` carries the full numbers,
including the ones that went the wrong way.

The default build has **zero dependencies**.

```
cargo test                      # core suite (no dependencies needed)
cargo test --features bignum    # + arbitrary-precision coefficients
cargo doc --open                # design docs

# Build the Python/Sage extension module (needs maturin):
maturin build --release --features python
mkdir -p pybuild && unzip -q -o target/wheels/*.whl -d pybuild
PYTHONPATH=pybuild sage -python -c "import symfn; print(symfn.schur_multiply([([2,1],1)],[([2,1],1)]))"

# Optional: compile the Sage shim's per-term loop
sage -python scripts/setup_cy.py build_ext --inplace
```

## Design in one screen

| Decision | What we did | Why |
|---|---|---|
| No untyped object | One **type per basis** (`Schur`, `PowerSum`, `Monomial`) behind the `SymFn` trait | Basis confusion becomes a compile error, not a runtime bug (vs. Symmetrica's `OP`) |
| Coefficients | Generic over `Ring`; dividing paths bounded on `QAlgebra`, not `Field` | Every division is by z_μ, an *integer* — so ℚ[t] and ℚ[q,t] qualify, which is what Macdonald/Hall–Littlewood need |
| Littlewood–Richardson | Behind the `LrBackend` trait; three native backends, `SkewLr` the default (no external C lib) | The trait paid off: each new backend swapped in with no caller changes and is cross-checked against the previous ones |
| Partitions | `Partition` newtype, invariant enforced at construction | Weakly-decreasing/positive guaranteed, not merely assumed |
| Correctness | Known-value + algebraic-law tests, committed oracle fixtures, and Sage driving symfn as its own backend | An oracle only tests inputs you thought of; letting Sage pick them found a 200x regression the benchmark could not see |

## Layout

```
src/
  coeff.rs      Ring / Field traits, i64+i128 rings, exact Rational field
  partition.rs  Partition newtype, conjugate, z(λ), partition generator
  lr.rs         LrBackend trait + native NaiveLr (incremental pruning)
  strip_lr.rs   StripLr row-strip DP; AutoLr, the backend the library uses
  skew_lr.rs    SkewLr — whole-shape expansion, merged frontier (default)
  kostka.rs     Kostka numbers K_{λμ} (SSYT counting)
  qt.rs         ℤ[q,t] / ℚ[q,t] coefficients — sparse, sorted, merge-accumulated
  hl.rs         Hall–Littlewood Q'_λ(x;t) by the Morris recursion
  kf.rs         Kostka–Foulkes K_{λμ}(t) from the Hall–Littlewood transition
  frac.rs       ℚ(q,t) with denominators kept factored — no bivariate gcd
  macdonald.rs  Macdonald P_λ(x;q,t) by the branching formula
  charge.rs     the charge statistic; K_{λμ}(t) by tableau enumeration (reference)
  character.rs  χ^λ(μ) via Murnaghan–Nakayama (β-number rim hooks)
  sym.rs        SymFn / SymAlgebra traits; all six bases; multiplication
  convert.rs    ToSchur / FromSchur hub; Jacobi–Trudi; Muir's rule; h↔e flip
  ops.rs        ω involution, Hall inner product, internal (Kronecker) product
  hopf.rs       SymTensor, skew Schur, SkewBy, coproduct, counit, antipode
  plethysm.rs   f[g] through the power-sum basis
  eval.rs       evaluation at an alphabet; principal specializations; dim λ
  guard.rs      overflow-reporting coefficients + the escalation scope
  memo.rs       the caches; python.rs  the PyO3 bridge
  lib.rs        crate docs, re-exports, roadmap
tests/
  oracle.rs        known Schur expansions + commutativity/associativity/degree
  algebra_laws.rs  ring-hom conversions, ω algebra map, Hall pairings, Δ algebra map
  qalgebra.rs      the library over ℚ[t] — a ring that is deliberately not a Field
  bignum.rs        exactness past i128
scripts/
  sage_backend.py   symfn as Sage's conversion backend, replacing Symmetrica
  check_backend.py  A/B the two backends through Sage itself
  symfn_cy.pyx      the shim's per-term loop, compiled
  check_hl.py       Q'_λ against Sage; bench_hl.py A/Bs Symmetrica's own C
  check_kf.py       K_{λμ}(t) against Sage's kfpoly, every pair including zeros
  check_macdonald.py  P_λ(x;q,t) against Sage, compared in the fraction field
```

## Features

| feature | what it adds |
|---|---|
| *(default)* | the whole library over `i64`/`i128` and an exact `Rational`; no dependencies |
| `bignum` | `BigInt` / `BigRational` coefficients (`num-bigint`, pure Rust) — exact beyond `i128` |
| `python` | PyO3 extension module (abi3, CPython 3.9+) with a coarse-grained API for Sage |

## Validation

- **Sage oracle** — `tests/sage_oracle.rs` checks 743 values Sage computed
  independently (Kostka numbers as tableau counts, characters, Schur products,
  all four conversions out of Schur, skew Schur). Fixture is committed.
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
- **Sage as the driver, not the oracle** — `scripts/check_backend.py` installs
  symfn in place of Symmetrica in Sage's own dispatch table and compares 4678
  computations against the C library it displaces, including Hall–Littlewood,
  Jack and Macdonald. This is the check that matters most: every other test
  uses inputs *we* chose, so it can only find bugs we thought of. Letting Sage
  pick them found a 200x regression on shape families the degree ladder never
  generated.

See [ROADMAP.md](ROADMAP.md) for benchmarks and what's next.

## License

**MIT OR Apache-2.0** — see [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).

symfn contains no third-party code. It is validated against two GPL programs
(Sage and `lrcalc`) by invoking them as external oracles to generate committed
test fixtures. The LR engine was written **clean-room** — specification and
implementation by separate parties, the implementer having no access to `lrcalc`
— with the spec committed at
[docs/cleanroom-spec-skew-lr.md](docs/cleanroom-spec-skew-lr.md) as the audit
trail. Every dependency, optional ones included, is permissively licensed, so the
published wheel carries no copyleft obligation. [NOTICE.md](NOTICE.md) has the
details.
