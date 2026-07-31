# The Python boundary, and running as Sage's backend

What it takes to be a drop-in replacement for Symmetrica underneath Sage:
the PyO3 module, the conversion-table shim, the arbitrary-precision
escalation path, and the coefficient-ring bounds that let ℚ[t] and ℚ[q,t]
through.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

## Known ceiling: the Python boundary is `i128`, Python integers are not

`src/python.rs` crosses coefficients as `i128`. That is wide enough for every
structure constant this library computes in practice, and it replaced an `i64`
boundary that silently truncated plethysm numerators — but it is still a
fixed width, while Python integers are arbitrary precision and the `gmp`
feature is exact.

So `--features gmp` and `--features python` do not compose: a Sage caller
cannot get bignum coefficients even though the core supports them. Closing
that needs a decision about how arbitrary-precision integers cross PyO3
(`num-bigint` via pyo3's feature, a decimal-string representation, or making
the module generic and exposing two variants), which is an API question rather
than an implementation one. Until then, the ceiling on that boundary is
symfn's, not Python's.

## Running as Sage's backend, in place of Symmetrica

`scripts/sage_backend.py` fills `sage.combinat.sf.classical.conversion_functions`
with symfn shims; `scripts/check_backend.py` A/Bs Sage against itself with only
the backend changed. **4678 computations agree at degree 8.**

> Superseded by [From one consumer site to five](#from-one-consumer-site-to-five)
> below: the adapter now displaces five of Sage's six consumer files, not just
> the conversion table, and the A/B is **8647 computations at degree 8**. The
> account below still describes the conversion table accurately.

This is a different kind of test from everything before it. Every earlier script
used Sage as an *oracle* on inputs we chose, so it could only find bugs we
thought to look for. Here Sage drives, and the comparison is against the C
library the shim displaces — same inputs, same code. Coverage is deliberately
indirect as well as direct: the 20 table entries head-on, but also the
operations that merely *reach* a conversion on the way to something else —
products in a non-Schur basis, `scalar`, `expand`, plethysm, `itensor`,
`skew_by`, coproduct, antipode, and the Hall–Littlewood, Jack and Macdonald
bases, which are defined by transitions from the classical ones.

**It found something within minutes that the conversion table alone cannot
show: Sage has two calling conventions.** `classical.py` passes a
`{Partition: coeff}` dict. But `sf.py`'s `SymmetricaConversionOnBasis` — the
wrapper that builds conversion *morphisms*, and therefore what **every non-QQ
base ring goes through** — passes a `CombinatorialFreeModule` element and calls
`dict()` on the result. A backend that handles only dicts passes every direct
table check and then fails the instant a caller touches Macdonald, Jack, HL, or
any ring other than QQ. Nothing in the table's shape hints at it; only letting
Sage drive surfaces it.

Both coefficient regimes are exercised for the same reason, since they are
distinct paths: over QQ the whole element crosses in one call, over ℚ[t] Sage
calls once per partition and recombines. The two halves run in **separate
processes** — Sage memoises conversion morphisms hard enough that swapping the
backend in-process risks comparing a cached answer with a fresh one and calling
it agreement.

## End to end: 1.84x like-for-like, 4.37x with a partition cache

⚠️ **Read both numbers.** The headline 4.37x includes a Sage-`Partition` cache
that Symmetrica's wrapper does not have and could equally well adopt — nothing
about it is specific to this backend. With it disabled
(`SYMFN_NO_PARTITION_CACHE=1`, kept for exactly this measurement) the same
benchmark gives **1.84x**, and that is the like-for-like figure: both sides then
construct a `Partition` per output term, as `symmetrica.pxi` does.

So the cache is worth **2.4x of the 4.37x** — more than half the end-to-end win
is a marshalling trick, not the Rust core. On the marshalling alone it is worth
7-9x (9.16x at degree 10, 7.07x at degree 18). Its cost is a table per degree
built on first touch: 0.19 ms at degree 10, 1.21 ms at 18, 12.7 ms at 30 for
p(30) = 5604 objects. Payback is roughly one conversion at degree 18, so it is
clearly right to keep — but the honest claim against Symmetrica is 1.84x plus
"and it should cache its partitions too".

The gap between 1.84x and the 5-10x the conversions themselves show is the
answer to why: at these sizes a conversion is microseconds of arithmetic wrapped
in milliseconds of object marshalling, and the marshalling is common to both.



`scripts/bench_backend.py` times Sage-level operations with the backend as the
only difference, in alternating separate processes. **Total 3.72x**, every row a
win (1.08x to 11.7x), across both coefficient regimes and degrees 10/14/18.

The first run said **1.91x**, with most rows at ~1.0x and several *below* it —
which is the number worth keeping in mind, because the Rust core is 5-10x faster
and that did not show up. **A faster core does not make Sage faster on its own.**
For anything but the heaviest conversions the shim's Python dominated:

| | time |
|---|---|
| the symfn call | 0.033s |
| + rebuilding Sage `Partition` / `ZZ` objects | 0.309s |
| the full shim entry | 0.389s |

**91% Python glue**, and the Partition rebuild alone was 9.3x the computation.
`_Partitions(list)` validates and interns on every call, once per output term.

The fix is a dict. Every conversion of degree n draws from the same p(n)
partitions, so caching them by part tuple turns construction into a lookup:
0.279s → 0.018s for that step, and 1.91x → 3.72x overall. `element_class`,
which skips validation, only reached 0.154s — so the cost is *construction*,
not checking, and avoiding it entirely is what matters.

A second pass added a **fast path for integral input to an integral basis** —
nearly every call. The general path walks the output four times (build with
`QQ(c)/den`, drop zeros, scan denominators to choose ZZ or QQ, build the dict);
the fast path stays in Python ints and walks it once. **4.35x** total.

The lesson generalises past this shim: at these sizes a classical-basis
conversion is microseconds of arithmetic wrapped in milliseconds of object
marshalling, and optimising the former without the latter is invisible. It is
the same coarse-grained argument `python.rs` opens with, one layer further out.

### The Cython interface, built

`scripts/symfn_cy.pyx` compiles the per-term loop; `scripts/setup_cy.py` builds
it. The import is optional — the pure-Python fallback is the same computation —
so a wheel without it still works.

Two changes made it worth doing:

* **Indices instead of partitions.** `symfn.convert_indexed` returns each output
  partition as its *position* in `symfn.partitions(degree)` rather than as a
  list of parts, so the shim reads a C array instead of building a tuple and
  hashing it. Sage's own wrapper cannot do this — `symmetrica.pxi` gets lists of
  parts back from C and calls `Partition(res)` on each, paying object
  construction per term.
* **`smallInteger` instead of `Integer(...)`.** Sage's internal constructor for
  values fitting a C long, reached through `PyLong_AsLongAndOverflow` with the
  generic parse kept for the rare escalated coefficient. This was the single
  biggest step: **110 → 71 ns/term**.

The loop went **220 → 71 ns/term (3.12x)** and the whole shim, on 40 shapes of
degree 14:

| | time | glue |
|---|---|---|
| original | 0.389s | 91% |
| + Partition cache | 0.134s | — |
| + integral fast path | 0.077s | 58% |
| + indices & compiled loop | **0.053s** | **40%** |

⚠️ **The end-to-end total barely moved — 4.35x to 4.36x — and that is the honest
headline.** `bench_backend.py` is weighted towards heavy conversions where the
Rust computation dominates and there is little glue left to remove. The Cython
win lands on the *light* rows, which is exactly where it should: at degree 10,
`s → h` went 1.48x → 3.40x, `s → e` 2.12x → 3.47x, `m → s` 2.25x → 3.68x.

Two things this establishes about the ceiling. Below about degree 10 the limit
is **Sage's own dispatch**, not either backend: a sweep of many small degree-8
conversions is 1.28x whatever we do. And at the top the limit is now symfn
itself — the Rust call is 60% of the shim, so further glue work has little left
to win. A C ABI (level 2) would attack the remaining 40%, and on this evidence
is worth perhaps another 1.5x on light workloads and nothing on heavy ones.

## The Python boundary's integer ceiling — decided: compute-and-escalate

Two corrections to what this file previously implied. **`gmp` and `python` do
compose** — `maturin build --features "gmp,python"` produces a wheel; that note
was stale. But enabling `gmp` changes *nothing* about the Python API, because
`python.rs` hardcodes `i128` in 20 places and every entry point builds
`Schur<i128>` or the i128-backed `Rational`. The GMP wheel is behaviourally
identical to the plain one.

The forcing issue is soundness rather than capability. `impl Ring for i128` uses
plain `*`, so in release **coefficient arithmetic wraps silently**. Characters
are already guarded (`character()` panics, `try_character` returns `None`) and
`integral_sweep` is fully checked with a clean bail-out — but nothing protects
generic coefficient arithmetic in products, plethysm, or conversions. At the
degrees a Sage user reaches that returns wrong answers with no signal, which is
not something a Symmetrica replacement can ship.

Three options: (A) keep i128 and refuse loudly, (B) always bignum, (C) compute
in i128 with checked arithmetic and re-run the whole call in `rug::Integer` if
anything overflowed. C is the crate's existing pattern — `try_character` →
`character_in` — and keeps the fast path fast, but only if "checked" is nearly
free. `examples/bench_guarded.rs` measures exactly that: `i128` and `Rational`
with every arithmetic op replaced by its `checked_` form, over Schur products,
s→m/h/e, m→s, s→p, plethysm and the internal product.

**Checked arithmetic costs 0–1%, indistinguishable from run-to-run noise.**

| workload | plain | checked | cost |
|---|---|---|---|
| integral (products, s→m/h/e, m→s) | 2.428s | 2.422s | −0.2% |
| dividing (s→p, plethysm, Kronecker) | 0.095s | 0.095s | −0.3% |

So **C is decided**. Escalation is the rare path, the common path is unchanged,
and the ceiling stops being visible to callers.

⚠️ **The methodology mattered more than the result, and this is the third time.**
With a fixed pass order the benchmark first reported checked rationals as **36%
faster** than unchecked. Swapping the two blocks moved the 36% to the other
type: whichever rational pass runs immediately after the (25x larger) integral
passes pays ~60% for arriving with a cold cache, and a fixed order silently
charges that to one type. Interleaving alone was not enough — the fix is
**rotating** the order so each variant spends an equal share of rounds in each
position. Two hypotheses were tested and discarded on the way: missing
`#[inline]` on `Rational`'s `Ring` impl (adding 25 of them changed nothing,
since generic instantiation already inlines) and a difference between `Rational`
and the hand-written twin (a byte-identical unchecked twin in the same crate
showed the same anomaly, which is what localised it to position).

### What the escalated ring needs to be: digits, not speed

`examples/coeff_sizes.rs` measures the two things that decide it.

**Coefficients are small.** Across LR squares, Kostka rows, characters, plethysm
numerators and s → p denominators, nothing in these workloads exceeds **2
limbs**, and most are one:

| workload | widest coefficient |
|---|---|
| LR `s_[10,9,8,7,6,5]²` | 18 bits |
| Kostka row, degree 32 | 55 bits |
| characters `p_1^40 → s` | 76 bits |
| f^λ at degree 60 | 118 bits |
| s → p denominators, degree 24 | 43 bits |
| plethysm numerators | 2–5 bits |

Even f^λ for a degree-60 staircase — 36 digits — still fits `i128`. Past the
ceiling we are therefore in the 2–5 limb range, where **every** bignum library
runs schoolbook: Karatsuba engages around 10–30 limbs, Toom near 100, FFT in the
thousands. GMP's advantage at this size is assembly tuning, not algorithms.

**Coefficient arithmetic is a small share of runtime.** Replacing `i128` with
heap-allocated `rug::Integer` *entirely* costs **1.04x** on products and s → m,
so at most ~4% of that workload is coefficient arithmetic; the combinatorial
traversal dominates. The rational path is the exception at **3.39x**, since
`rug::Rational` pays a gcd per operation — but on a 100x smaller absolute base,
and its hot case (`integral_sweep`) already runs in raw `i128` under a common
denominator.

So the requirement is **exactness, not speed**, and a pure-Rust bignum is
enough. That also settles the wheel: `rug` and `gmp-mpfr-sys` are **LGPL-3.0+**,
and statically linking them into a distributed wheel would attach LGPL terms to
a crate that is deliberately MIT OR Apache-2.0. `num-bigint`/`num-rational` are
permissive, need no C toolchain or `m4` to build, and PyO3 0.29 has
`num-bigint`/`num-rational` features that convert them to Python `int`/`Fraction`
natively — which removes the decimal-string encoding step as well. `gmp` stays
an optional feature for Rust callers who want it.

### Implemented

`src/guard.rs` holds `Guarded` / `GuardedRat` — `i128` arithmetic that *reports*
overflow instead of wrapping — and `guarded(|| …) -> Option<T>`. Every
element-valued entry point in `python.rs` now runs fixed-width first and re-runs
over `BigInt` / `BigRational` if anything overflowed. `character_value` escalates
through `try_character` → `character_in`.

**Overflow is recorded by a monotone global counter, not a flag**, and that is a
correctness requirement rather than a style choice. Clear-run-test loses answers
under concurrency: two overlapping computations, and one can clear the flag
*after* the other set it, so the second reports success on a wrapped result. A
counter compared before/after cannot do that, and its failure direction is the
safe one — an unrelated overflow forces a needless escalation, costing time and
never correctness. Global rather than thread-local for the same asymmetry: a
thread-local would *miss* a worker thread's overflow, which is the unsafe
direction. The guard tests must therefore be serialised against each other,
since one test's deliberate overflow is visible to another's scope; that showed
up immediately as two failures that passed in isolation.

**Coefficients cross as a `Coeff` enum, not as `BigInt`.** Python sees a plain
`int` either way — the enum never escapes Rust, so there is no mixed-type list —
but `BigInt` is heap-allocated and almost every coefficient is small. Measured
against the old `i128` boundary:

| workload | all-`BigInt` | `Coeff` enum |
|---|---|---|
| term-heavy pass (171k terms) | 1.076x | **1.005x** |
| `coproduct` (marshalling-dominated) | 1.224x | **0.969x** |

So the ceiling is gone for free. Verified against Sage on cases that previously
came back **wrapped, with no signal**: a product with 10³⁰ coefficients, χ^λ(1⁷⁸)
= 1.789…e48, s → p with a 10⁴⁰ input, and a Hall product of 10⁵⁰. The full
ladder, 660 evaluation checks and 6660 skew checks still agree.

Still fixed-width, and now documented as such: `character_table` (`i128`) and
`kostka_table` (`u128`) are built on fixed-width accumulators inside the sweep,
so widening their signatures alone would not help. `character_value` escalates
per entry and is the exact route past |λ| ≈ 58.

## Coefficient rings that are not fields (ℚ[t], ℚ[q,t])

Prompted by the question of what it would take to serve a Sage user working over
`QQ['t']`. Most of the answer is "nothing": **Sage already factors the
coefficient ring out of the conversion path.** From `sage/combinat/sf/classical.py`:

```python
if R == QQ and P.base_ring() == QQ:
    return self._from_dict(t(m)._monomial_coefficients, coerce=True)   # one bulk call
f = lambda part: self._from_dict(t({part: ZZ.one()})._monomial_coefficients)
return self._apply_module_endomorphism(x, f)                            # integer rows, R-recombine
```

Over any ring that is not ℚ, Symmetrica is asked only for the **integer
transition row of a single basis element**, and Sage does the arithmetic in R
itself. Structure constants are ring-independent, so `python.rs` already exposes
what is needed. (The non-QQ path is one FFI call per partition in the support,
which a bulk entry point would collapse — but measured at degree 12, 77 terms,
`QQ['t']` took 0.053s against `QQ`'s 0.066s, so this is structural tidiness and
not a bottleneck.)

What *did* block ℚ[t] was our own bound. **Every division in this library is by
z_μ, a positive integer** — `s → p` carries z_μ⁻¹ and the internal product and
plethysm inherit it through the power-sum basis. Nothing divides by a general
ring element. But those paths were bounded on `Field`, which demands the ability
to invert *t*, so ℚ[t] and ℚ[q,t] were excluded from operations that are
perfectly well defined over them — and those are exactly the coefficient rings
Hall–Littlewood and Macdonald need.

`Field` is replaced in those bounds by `QAlgebra` — a ring containing ℚ — with
the single method `div_u128`. `Field` is kept (it is a real thing to name, and
`Rational` is one) but nothing in the library requires it. The implication
"characteristic-0 field ⟹ ℚ-algebra" is deliberately *not* a blanket impl: that
would occupy the impl for every downstream type, and a polynomial ring could
then never implement it.

**Plethysm needed a second, less obvious fix, and it was silently wrong.**
Everything else in the library is linear with integer structure constants, so a
coefficient is only ever multiplied and added. `p_n` is different: it substitutes
into the alphabet, and a coefficient ring's variables belong to that alphabet.
Sage:

```text
sage: R.<t> = QQ[];  p[2](t*p[1])  ->  t^2*p[2]        p[2](t*p[1], exclude=[t]) -> t*p[2]
```

`scale_parts` carried coefficients through unchanged — correct over ℚ, where
there is nothing to raise, which is why the omission was invisible for as long
as ℚ was the only coefficient ring in use. It is now `Plethystic::frobenius`, a
separate trait above `QAlgebra` so that `s → p` (which divides but never
substitutes) does not demand it. Which variables get raised is a genuine
convention choice, so it is the implementor's to state rather than a default.

`tests/qalgebra.rs` is the proof: a ℚ[t] implementing `Ring + QAlgebra +
Plethystic` and deliberately **not** `Field`, exercising `s → p`, the round
trip, the Kronecker product, and plethysm. If any of those paths regressed to a
`Field` bound the file would stop compiling, which is a stronger assertion than
its `assert_eq!`s. Plethysm values are Sage's, including the degree-3 cases
(`s[3](t*s[1]) = t^3*s[3]`) that a doubling bug could not fake. Over ℚ nothing
changed: the full Sage ladder still agrees.

Remaining for a real Sage backend, in order: arbitrary-precision integers across
the FFI boundary (the known `i128` ceiling, plus `--features gmp` and
`--features python` still not composing), then the bulk expansion entry point.

## From one consumer site to five

`docs/symmetrica-coverage-audit.md` found that Sage reaches Symmetrica from six
files and that the conversion table — all `sage_backend.py` displaced — was one
of them. This chapter is what happened when its remaining task list was
implemented: the four tasks landed, but **two of the three gap
classifications were wrong**, and a fourth gap the audit did not see turned up
in the file it had marked complete.

**The portable lesson: an entry point's name tells you what it computes; only
its caller tells you what it must return.** Both misclassifications came from
matching a Symmetrica name to a symfn name and stopping there.

| audit said | actually |
|---|---|
| `compute_*_with_alphabet`: bind `eval()`, "a binding gap, not a mathematics gap" | the reverse — `eval()` cannot produce it at all |
| `mult_monomial_monomial`: "`Monomial::mul` exists in `sym.rs`" | it was in `convert.rs`, routed through Schur |
| `kostka_tab`: the one real gap | correct — and its *order* is contract too |
| the Schubert seven "map one-to-one" | six do; `scalarproduct_schubert` is a different operation |

**`eval` and `expand` answer different questions, and only one of them is
`expand(n)`.** `sfa._expand` does `resPR(e(part, n, alphabet))` — it wants a
polynomial in `n` indeterminates. `Schur::eval` evaluates at an alphabet of
`Ring` elements and returns one value; the indeterminates are not ring elements
the coefficients live in, so no binding bridges the two. The fix was a new
operation, `Monomial::expand`, emitting `(exponent vector, coefficient)` pairs.
Routing the other five bases through `m` costs one conversion and makes `m` the
only basis that needs an expansion rule — which is right, because `m`'s
expansion *is* its definition.

That operation then paid for itself twice: the distinct-rearrangement multiset
walk it needs is also what the monomial product needs, so `Monomial::mul` became
a direct rule — enumerate `α + β` slotwise, keep the sum weakly decreasing, and
the multiplicity with which λ arrives *is* the structure constant. The old
`m → s → LR → m` route costs the whole degree, because `m → s` inverts the
Kostka matrix; it is kept as the reference oracle, the role `NaiveLr` plays for
LR, and `monomial_product_agrees_with_the_schur_route` is the agreement test.

**`kostka_tab`'s order is part of the interface, and no amount of reading the
audit would have said so.** `SemistandardTableaux(λ, μ).list()` returns the
backend's list verbatim and Sage's doctests print it. The order is increasing
lexicographic in the row-major reading word — which the natural
chain-of-horizontal-strips walk does *not* produce, because that walk groups by
value and the reading word orders by position. They first disagree at three-row
shapes. Generating by chains and sorting was chosen over generating in reading
order: it keeps the existing pruned walk, costs a log factor on an enumeration
already paying `K_{λμ}·|λ|`, and turns the contract into a stated property
rather than an artifact of a traversal. 1818 (λ, μ) pairs through degree 9 agree
with Symmetrica exactly, order included.

**Two patching mechanisms, because the six sites bind differently.** Five reach
a function through the `sage.libs.symmetrica.all` *module object* at call time —
by `getattr`, by attribute access, or by `lazy_import` of the module — so
rebinding attributes on that module reaches all five at once.
`sf/hall_littlewood.py` does `from ... import hall_littlewood_symmetrica as
hall_littlewood` at import, so the name has to be rebound in that module's own
namespace. A shim that only did the first would silently miss Hall–Littlewood,
which is exactly the failure mode the earlier `_items` discovery had.

`check_backend.py` grew sections for the five: **8647 computations at degree 8,
0 mismatches**, up from 4678 covering the conversion table alone.

### Why `schubert_polynomial.py` is still on Symmetrica

The audit marked all seven Schubert entry points covered. Six are: 1679
comparisons over `S₁`–`S₄` found **zero disagreements on any input Symmetrica
answers**. The seventh is not covered at all —
`scalarproduct_schubert` returns a Schubert *polynomial* and
`schubert_pairing` returns an integer:

```text
  X([2,1]).scalar_product(X([2,1]))  =  X[1,3,2]     (Symmetrica)
  symfn.schubert_pairing([([2,1],1)], [([2,1],1)], 2)  =  0
```

Two more facts, recorded so the next session does not rediscover them:

- **symfn is more total than Symmetrica here, which is a problem, not a
  feature.** 187 of the 1679 comparisons are cases where Symmetrica or its Sage
  wrapper raises `ValueError` — `∂_i` past the permutation's length, `∂_w` on
  the identity — and symfn returns the correct value. Sage's doctests assert
  those exceptions, so a faithful drop-in must reproduce them.
- **Symmetrica leaks permutations and then blocks on a prompt.** After a run of
  `divdiff_perm_schubert` calls it prints `ERROR: permutation memory not
  freed?: mem_counter_perm = 99` at teardown and drops into an interactive
  menu (`enter a to abort with core dump, g to go, …`). On a non-tty that hangs
  forever — the probe has to close stdin. This cost most of an hour before the
  cause was visible, and it is the strongest robustness argument for
  displacement the project has found so far.

Displacing six of the seven would leave Symmetrica loaded for the seventh, so
the adapter leaves the whole file alone for now.

**Decided: `scalarproduct_schubert` will not be reimplemented.** Symmetrica
stays available to Sage as an *optional* package rather than being removed, so
the operation keeps working for anyone who installs it. The fact that makes this
cheap is that **nothing in sagelib calls
`SchubertPolynomial.scalar_product`** — its only references are its own
definition and its own doctests, so it is public API with no internal
dependents, in exactly the position of the 30 unreached entry points.

The decision generalises past this one function, which is why it is recorded
here rather than as a footnote: **displacement does not have to mean removal.**
Demoting Symmetrica from `type: standard` to `type: optional` answers all 31
entry points symfn will not cover in a single packaging change, and retires the
multi-release deprecation cycle that removing public API would have required. It
converts the hardest part of the upstream ask into the easiest.

It also inverts the "wire six of seven buys nothing" conclusion above. While
Symmetrica is standard, that is true. Once it is optional, wiring the six is
what keeps Schubert polynomials working for users who do not install it, with
only `scalar_product` behind the feature gate — so it becomes worth doing before
the demotion lands, not never. The exception-fidelity requirement is still the
price of admission.

## The boundary raises where it panicked: 5 clusters, 30 entry points

The premise this started from was that `part()` calls `Partition::new`, "which
asserts", so a non-partition reached Sage as a `PanicException`. **`Partition::new`
does not assert.** It normalizes — filters zeros, sorts weakly decreasing — so
`symfn.schur_multiply([([1,3], 1)], …)` returned, cheerfully, the product for
`[3,1]`. The bug was real and worse than the one described: not a crash but a
well-formed answer to a question the caller had not asked, across roughly 50
entry points. Recorded because the correction is the interesting part — the
convenience constructor was the leak, and "does it panic?" was the wrong
question to audit by.

Probed by building the cdylib and driving it from CPython, five clusters
produced a `PanicException` from correctly-typed input. All five are now
`ValueError`, and `docs/policies/failure.md` gained R11 for the rule they share.

| cluster | entry points | trigger | was |
|---|---|---|---|
| non-permutation term list | 6 Schubert | `schubert_multiply([([1,1],1)], …)` | `unwrap()` on `None` |
| index past `MAX_SUPPORT` | 4 Schubert | `schubert_pairing(…, 33)`, `…_variable(…, 32)` | `w0 is a permutation: TooLarge` |
| `k = 0` | 9 LLT + `k_core_quotient` | `llt_gtilde([2,1], 0)` | `a ribbon level needs k ≥ 1` |
| inhomogeneous argument | 6 Macdonald operators | mixed degrees to `nabla` | `degree_of`'s assert |
| capacity walls | `llt_g`, the abacus five, `llt_e_expansion`, `nabla_e_by_path`, `htilde_by_llt` | `llt_h([130], 1)` | `past this abacus's 128 bits` |

Three findings worth keeping:

- **The Schubert unwrap was reachable precisely because coefficients are
  small.** `build_schubert` returned `None` for two unrelated reasons — a
  coefficient too wide for the fixed-width pass, and a word that is not a
  permutation — and the escalation path unwrapped it. The fast pass declined
  the malformed word by returning `None`, so escalation ran and the unwrap
  fired; a *wide* coefficient would have taken the same path legitimately. The
  fix is structural rather than a better message: a `Wide` trait marks the rings
  that cannot decline, so `build_wide` has no `Option` to unwrap, and
  `schub_terms` validates words before either pass. That removed every `unwrap`
  in `python.rs`. This duplicates the fix on `claude/codebase-failure-policy-2854d6`
  (commit 3512be8) and deliberately keeps its trait name and shape, so the two
  branches converge rather than conflict.

- **`convert_indexed` crashed on a *valid* input.** The `Schur → Schur` identity
  passed the caller's raw list to `index_of`, whose table is keyed by normal
  forms, so `[2,1,0]` — padding, the one tolerance this module documents,
  because Sage hands over fixed-width lists — was a missing key. Validating
  restored it: the test asserts padded input still answers, since a validation
  pass that costs the padding tolerance would break the Sage adapter.

- **Capacity walls belong to the module that knows the bound.** The abacus reach
  is `ℓ(λ) + λ₁ + k` for the pruned walk and `2kn + k` for the unpruned one, and
  the note at `assert_abacus_fits` explains why they cannot share. Restating
  either at the boundary would have been a second copy to drift, so `llt.rs`
  exposes `abacus_reach`, `abacus_reach_table`, `MAX_CELLS` and `free_edges`,
  and the asserts and the boundary check the same expression.

The pin is `scripts/check_python_boundary.py` — 111 malformed calls over 85
pyfunctions, asserting each raises a typed exception, plus a completeness check
that fails when a new `#[pyfunction]` appears with neither a case nor an entry
in its `TOTAL` list. It is Sage-free and imports the `cargo build` artifact
directly, so it needs no maturin step; it is not part of `scripts/preflight.sh`,
which builds only the default features. A Rust unit test is not available here:
an `extension-module` binary has no interpreter, so a test that merely
constructs a `PyErr` aborts at load.

### The zeros: which are theorems and which were conventions

The first pass left every precondition violation that returned a plausible `0`
alone, listed as open. Working through them, they are **two kinds**, and the
distinction is the whole answer:

- **A zero that is a theorem.** `c^λ_{μν} = 0` off-degree because the product is
  homogeneous; `K_{λμ} = 0` unless λ dominates μ because there are no such
  tableaux; `s_λ` in `n` variables vanishes when `ℓ(λ) > n`. These are values,
  and a caller sweeping a range depends on getting them. `kostka.rs` already
  documented its own ("returns 0 unless |λ| = |μ|"). **Unchanged**, and each now
  says outright in rustdoc that zero is an answer — plus
  `check_theorem_zeros_still_answer` in the boundary script, so a later pass of
  "validate more" cannot quietly eat them.

- **A zero that was a convention over an undefined question.** `χ^λ(μ)` needs μ
  to index a class of `S_{|λ|}`; `g^ν_{λμ}` and `a^λ_{μν}` need all three in one
  `S_n`; `K_{λμ}(q,t)` is an entry of one degree's matrix. Off-degree there is
  no value, and `ops.rs` said as much in its comment — "unequal degrees pair to
  zero" is labelled a convention, not a theorem. **These five now raise**
  (`character_value`, `kronecker_coefficient`, `class_algebra_coefficient`,
  `qt_kostka`, `schubert_pairing`).

`schubert_pairing` is the sharpest of the five and did not look like a degree
question at all. Its `n` is explicit precisely because Symmetrica's
`scalarproduct_schubert` infers it from padding, so the same inputs give
different answers — but a term outside `S_n` simply could not contribute to the
coefficient of `w0(n)`, so it came back `0`, and the caller could not tell that
from an honest zero. The explicit `n` removed the trap one level down and left
it one level up.

The split is deliberate about **who** is calling. The Rust-side
`kronecker_via_characters` keeps returning `0` off-degree, because
`internal_product` reaches the same zero by having no shared λ and the two
routes must agree — the Schubert `mul`/binding divergence in
[schubert.md](schubert.md) is what disagreeing costs. Only the boundary is
stricter, which is the licence R11 already grants: totality is what a composing
Rust caller needs, and a typed refusal is what a foreign caller needs, and they
are not the same requirement.

Boundary script: 116 malformed calls, plus the theorem-zero pins and the
padding pins.
