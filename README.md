# symfn

A modern, fast, testable kernel for computing with **symmetric functions** — a
clean-room Rust successor *in spirit* to Symmetrica (it shares no code with the
old C library), designed to interoperate with Sage without paying Python
overhead on the hot path.

## Status

**Complete and Sage-validated.** All five bases {m, e, h, p, s} as real types
with multiplication; conversions between every ordered pair; the ω involution,
Hall inner product, skew Schur functions, coproduct, counit, and antipode.
Plus **plethysm** and a memoized, size-dispatching Littlewood–Richardson engine.
Cross-validated against Sage on 770 independently-computed values. Optional GMP
coefficients and a PyO3 bridge importable into Sage (typically **8–60x** faster
than the pure-Python path). The default build has **zero dependencies**;
59 tests green.

```
cargo test                      # core suite (no dependencies needed)
cargo test --features gmp       # + GMP-backed bignum coefficients
cargo doc --open                # design docs

# Build the Python/Sage extension module:
cargo build --release --features python
mkdir -p pybuild && cp target/release/libsymfn.dylib pybuild/symfn.so   # .so on Linux
PYTHONPATH=pybuild sage -python -c "import symfn; print(symfn.schur_multiply([([2,1],1)],[([2,1],1)]))"
```

## Design in one screen

| Decision | What we did | Why |
|---|---|---|
| No untyped object | One **type per basis** (`Schur`, `PowerSum`, `Monomial`) behind the `SymFn` trait | Basis confusion becomes a compile error, not a runtime bug (vs. Symmetrica's `OP`) |
| Coefficients | Generic over `Coeff` (`i64` now) | `gmp` feature swaps in `rug::Integer`; later carries `(q,t)` for Macdonald |
| Littlewood–Richardson | Behind the `LrBackend` trait; three native backends, `SkewLr` the default (no external C lib) | The trait paid off: each new backend swapped in with no caller changes and is cross-checked against the previous ones |
| Partitions | `Partition` newtype, invariant enforced at construction | Weakly-decreasing/positive guaranteed, not merely assumed |
| Correctness | Known-value + algebraic-law tests; Sage as the eventual oracle | Born tested against an independent implementation |

## Layout

```
src/
  coeff.rs      Ring / Field traits, i64+i128 rings, exact Rational field
  partition.rs  Partition newtype, conjugate, z(λ), partition generator
  lr.rs         LrBackend trait + native NaiveLr (incremental pruning)
  strip_lr.rs   StripLr row-strip DP; AutoLr, the backend the library uses
  skew_lr.rs    SkewLr — one traversal per shape, binned by content (default)
  kostka.rs     Kostka numbers K_{λμ} (SSYT counting)
  character.rs  χ^λ(μ) via Murnaghan–Nakayama (β-number rim hooks)
  sym.rs        SymFn / SymAlgebra traits; all five bases; multiplication
  convert.rs    ToSchur / FromSchur hub; Jacobi–Trudi; inverse Kostka
  ops.rs        ω involution, Hall inner product
  hopf.rs       SymTensor, skew Schur, coproduct, counit, antipode
  lib.rs        crate docs, re-exports, roadmap
tests/
  oracle.rs        known Schur expansions + commutativity/associativity/degree
  algebra_laws.rs  ring-hom conversions, ω algebra map, Hall pairings, Δ algebra map
```

## Features

| feature | what it adds |
|---|---|
| *(default)* | the whole library over `i64`/`i128` and an exact `Rational`; no dependencies |
| `gmp` | `rug::Integer` / `rug::Rational` coefficients — GMP Karatsuba/Toom/FFT, exact beyond `i64` |
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
- **Hopf axioms** — the antipode satisfies `m∘(S⊗id)∘Δ = ε·1`.

See [ROADMAP.md](ROADMAP.md) for benchmarks and what's next.

## License

**MIT OR Apache-2.0** — see [LICENSE-MIT](LICENSE-MIT) and
[LICENSE-APACHE](LICENSE-APACHE).

symfn contains no third-party code. It is validated against two GPL programs
(Sage and `lrcalc`) by invoking them as external oracles to generate committed
test fixtures, and its LR engine uses published algorithmic ideas that `lrcalc`
also uses — but no code from either. The optional `gmp` feature links LGPL
libraries, which does not affect symfn's own terms. [NOTICE.md](NOTICE.md) has
the details.
