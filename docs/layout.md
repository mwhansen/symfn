# Layout

Every module in the tree, in one map. A module added, removed or renamed
updates this file in the same change.

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
  cache_accounting.rs  what cache_stats reports against what clear_caches releases
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
