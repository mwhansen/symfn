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
unit of work). **Macdonald** `P_λ`, `Q_λ` and `J_λ` by the branching formula, **~94× Sage**, over a
ℚ(q,t) that avoids bivariate gcd by keeping denominators factored. Both
Hall–Littlewood bases: `Q'` from the recursion, `P` by inverting the
Kostka–Foulkes matrix — and `q = 0` turns Macdonald `P` into Hall–Littlewood
`P`, checked between two computations that share no code. On top of those, the
**(q,t)-Kostka polynomials** `K_{λμ}(q,t)` from `J_μ = Σ_λ K_{λμ} S_λ(x;t)`,
which `q = 0` sends back to Kostka–Foulkes and `q = t = 1` to the number of
standard tableaux — and the modified form `H̃_μ = Σ_λ K̃_{λμ}(q,t) s_λ`, which is
where the literature states Haiman's positivity. The Kostka matrix comes from
the Bergeron–Haiman Pieri recursion, **~18× Sage** at degree 12, with the
branching formula and a Lapointe–Lascoux–Morse eigenvector solve kept as
independent cross-checks — three algorithms sharing nothing above `Partition`.

Newest, and the first thing here aimed at a gap rather than at parity: the
**Orellana–Zabrocki character bases** `s̃_λ` and `h̃_λ`, whose outer-product
structure constants are the **reduced (stable) Kronecker coefficients**. The
product never leaves the power-sum basis — a consequence of reading their
Theorem 14 as a statement about a linear map, which makes both transitions
invertible one part-size at a time — so no Littlewood–Richardson coefficient is
computed in a reduced Kronecker calculation at all. `st[4,3] · st[4,3]` is
**3400× Sage**, and `st[6,4] · st[6,4]`, the case `docs/research-gaps.md`
recorded Sage timing out on, takes 0.14s. The fixed-width wall is at total
degree 24 and is `z_γ` rather than the answers, which are 16 bits; `bignum`
escalation carries it to 32. Validated on 208 Sage comparisons with no
mismatches, plus two expansions printed in the paper.

Coefficients are generic over a `Ring`; the paths that divide ask only for a
`QAlgebra` (a ring containing ℚ), so ℚ[t] and ℚ[q,t] work even though neither is
a field. Arbitrary precision is automatic at the Python boundary: a call runs in
fixed width and re-runs exactly if anything overflows.

Validation is layered — 202 unit and integration tests, algebraic-law suites,
committed fixtures from Sage and `lrcalc`, and **8647 computations driven by
Sage itself** with symfn substituted for Symmetrica at five of the six places
Sage calls it (`scripts/check_backend.py`), covering Hall–Littlewood, Jack and
Macdonald as well as the classical bases, `expand`, the monomial product and
semistandard tableaux. The sixth site, Schubert polynomials, is still on
Symmetrica — one operation is missing, and
[docs/symmetrica-coverage-audit.md](docs/symmetrica-coverage-audit.md) says
which. The Python layer is checked separately from the
library (`scripts/check_bindings.py`): a correct answer marshalled into the
wrong slot is a different failure from a wrong answer, and only one of the two
shows up in a dump.

Performance, with the caveats that matter: the classical-basis conversions run
**2–10x** Symmetrica's C. End to end *through Sage* the same substitution is
**1.84x** like-for-like, or 4.37x with a Sage-`Partition` cache that Symmetrica's
wrapper does not have and could equally adopt. The gap between those is object
marshalling, which both backends pay. [The record](docs/record/) carries the
full numbers, including the ones that went the wrong way.

The default build has **zero dependencies**.

```
cargo test                      # core suite (no dependencies needed)
cargo test --features bignum    # + arbitrary-precision coefficients
scripts/preflight.sh            # the commit gate: fmt check + both suites
cargo doc --open                # design docs

# Once per clone: the versioned pre-commit hook (rustfmt check), and the
# formatting-only commits that git blame should look straight through.
git config core.hooksPath .githooks
git config blame.ignoreRevsFile .git-blame-ignore-revs

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
  skew_lr.rs    SkewLr — whole-shape expansion, merged row layer (default)
  two_row.rs    s_μ·s_ν with a two-row factor, counting fibres per output
  three_row.rs  the three-row analogue of the same counting route
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
  dyck.rs       labelled Dyck paths; the Delta conjecture's combinatorial side
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
tests/
  oracle.rs        known Schur expansions + commutativity/associativity/degree
  algebra_laws.rs  ring-hom conversions, ω algebra map, Hall pairings, Δ algebra map
  sage_oracle.rs   743 Sage-computed values, from a committed fixture
  lrcalc_oracle.rs products and skew expansions vs lrcalc, past Sage's sizes
  qalgebra.rs      the library over ℚ[t] — a ring that is deliberately not a Field
  bignum.rs        exactness past i128
  memory.rs        peak-bytes and allocation budgets over the measure workloads
  fixtures/        the committed oracle outputs both *_oracle suites read
docs/
  style.md             the prose rulebook, for every documentation surface
  policies/failure.md  how the library is allowed to fail
  record/              the memory — one file per subsystem; README.md indexes
examples/  research drivers — instruments, not demos (docs/style.md)
  *_dump.rs        emitters whose output scripts/check_*.py hold to Sage
  bench_*, profile_*, probe_*  per-subsystem instruments of the record
  delta_conjecture.rs, find_nonzero.rs  conjecture checks that state which
                   disagreement is a bug and which is a discovery
scripts/   nearly all need Sage; scripts/README.md documents the main ones
  preflight.sh      the local gate: fmt check + both test suites, no Sage
  sage_backend.py   symfn as Sage's conversion backend, replacing Symmetrica
  check_backend.py  A/B the two backends through Sage itself
  check_bindings.py the Python layer itself against Sage, not a dump
  symfn_cy.pyx      the shim's per-term loop, compiled (setup_cy.py builds it)
  check_*.py        one Sage oracle per subsystem (hl, kf, macdonald, jack,
                    llt, qt_kostka, deltaop, eval, skew, st, …)
  bench_*.py        the Sage side of each ladder, same work on both sides
  spec_*.py         pre-implementation verification and wall measurement
  gen_*, compare_*  fixture generators; direct lrcalc/Symmetrica comparisons
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
  symfn in place of Symmetrica at five of Sage's six call sites and compares 8647
  computations against the C library it displaces, including Hall–Littlewood,
  Jack and Macdonald. This is the check that matters most: every other test
  uses inputs *we* chose, so it can only find bugs we thought of. Letting Sage
  pick them found a 200x regression on shape families the degree ladder never
  generated.

See [the record](docs/record/) for what was built and measured, subsystem by
subsystem, and [docs/release-readiness.md](docs/release-readiness.md) for what
stands between this tree and a package.

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
