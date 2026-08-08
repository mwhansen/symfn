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
processes** — Sage memoizes conversion morphisms hard enough that swapping the
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

At these sizes a classical-basis conversion is microseconds of arithmetic
wrapped in milliseconds of object marshalling, so optimizing the former without
the latter is invisible. It is the same coarse-grained argument `python.rs`
opens with, one layer further out.

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

### Supports cross out as tuples: 16 ns/term off the keying step

Every output partition, permutation, composition and exponent vector used to
arrive as a Python `list`, which is PyO3's default for `Vec<u32>` and was never
a decision. It is now a `tuple` — one newtype, `Key` in
[python.rs](../../src/python.rs), carrying the outbound conversion and an
inbound one that still accepts any sequence, so no caller's existing input
broke.

The argument is that a support is a *key*: every consumer puts it in a `dict`
or a `set` on arrival, and a list has to be copied into a tuple before either
can hold it. On the s → m expansion of `s_[9,5,3,1]` (300 output terms),
building the result dict costs **46.5 ns/term from tuple keys against 62.6 from
list keys — 16.1 ns/term, 1.35x on that step**. Measured with `timeit` over
2000 repetitions, Python 3.12, **on battery** (so treat the absolute numbers as
soft and the ratio as the durable quantity — [README.md](README.md) documents
the 1.8x drift).

For scale: the compiled per-term loop is 71 ns/term and the pure-Python one 185
ns/term, so this is a fifth of the compiled loop's whole budget on any path that
keys by the support. The *indexed* path (`convert_indexed`) does not key at all
— it reads a position out of a table — so it neither gains nor loses here; the
win lands on the general path and on every non-Sage consumer, which has no
index table to read from.

`scripts/check_bindings.py` (0 failures) and `scripts/check_python_boundary.py`
(116 malformed calls, every one a typed exception) both pass unchanged across
the switch, which is what establishes that only the container type moved.

## Supply both directions; never let Sage invert

A technique, established on Hall–Littlewood `P` and expected to apply wherever
this backend meets a non-classical basis.

**The shape of the problem.** Sage stores a change of basis as two dictionaries
and computes one from the other with `_invert_morphism`: it fills the direction
it knows one matrix entry at a time, then recovers the other by triangular
back-substitution **over the fraction field** — `ℚ(t)` for Hall–Littlewood,
`ℚ(q,t)` for Macdonald, `ℚ(α)` for Jack. Both halves are expensive and the
second is usually the larger. On HL `P` at degree 15: 14.4 s to fill, ~6.9 s to
solve, 21.0 s total.

**The move.** Where symfn has both directions, hand Sage both and let
`_invert_morphism` go unused. For HL `P` the two are the same Kostka–Foulkes
matrix read two ways — `s_μ = Σ_λ K_{μλ}(t) P_λ` against
`Q'_λ = Σ_μ K_{μλ}(t) s_μ` — so `s → P` *is* `K`, and `P → s` is the inverse
symfn already takes by back-substitution in **`ℤ[t]`**, which never divides
because `K` is unitriangular in dominance order. Sage's solve does the same
algebra over the fraction field and pays for it.

**It is worth more than making either direction faster.** Supplying only
`P → s` and still letting Sage invert took degree 15 from 21.0 s to 2.32 s;
dropping the inversion as well took it to 0.498 s. The inversion was 95% of what
the first version left.

**Then look at the marshalling, because it becomes the cost.** Building each
entry in `ℤ[t]` from its coefficient dictionary and coercing once — rather than
summing `c·t^e` inside the fraction field, which builds a rational function per
monomial — was the difference between 2.32 s and 0.498 s. That runs once per
nonzero entry of a `p(n) × p(n)` table, so at these sizes it is not a
micro-optimization. This is the Amdahl argument from
[the shim's own history](#the-cython-interface-built) arriving one layer out:
remove the algorithmic cost and the boundary is what is left.

**Check it exactly, not by sampling.** Both cache dictionaries are
bit-identical to the ones the old route produces, for every degree through 8,
compared in separate processes. A change of basis is exactly the kind of object
where a plausible wrong answer survives spot checks, and Sage's own output is
available as the oracle.

### Which half is the cost is not guessable — measure the split first

Four changes of basis, four different answers to "where does the time go":

| basis | the expensive half | what fixed it | result |
|---|---|---|---|
| Hall–Littlewood `P` | the **inverse** (6.9s of 21.0s, and the fill was the rest) | both directions supplied | **42x** |
| Hall–Littlewood `Q'` | the **inverse** | both directions, the second one *transposed* | **15x** |
| Jack `P` | the **forward** fill — Gram–Schmidt, 17.1s against 0.14s to invert | forward only; Sage's inverse untouched | **74x** |
| Macdonald `J` | the **inverse**, and symfn had no `s → J` | forward only, then both once `s → J` existed | 1.4x, then **4.9x** |

The pattern: supplying one direction is worth 1.2–1.5x when the inverse is the
cost, and everything when the fill is. Supplying **both** is what turns 1.3x
into 15x. Before touching one of these, time the fill and the solve separately —
the four rows above would each have been mispredicted.

### Macdonald `J`, the other direction: 1.9x to 4.9x

`s → J` exists now ([qt-kostka.md](qt-kostka.md), "The inverse of `J → s` is a
projection, not a solve"), so `_s_cache` fills both caches and returns without
calling `_invert_morphism` at all — the third and last of the four to reach that
shape.

```text
  scripts/bench_macdonald_cache.py 11 2, ⚠️ on AC power
  n   cells  symmetrica       symfn    ratio
  7     116     0.5037s     0.0604s    8.33x
  8     238     1.7324s     0.2209s    7.84x
  9     430     5.4143s     0.7397s    7.32x
 10     818    17.9231s     2.8850s    6.21x
 11    1426    54.8776s    11.2200s    4.89x
```

The forward-only state measured in the same session on the same power gives
30.2s at degree 11 and 1.91x, so dropping the solve is 2.7x of the symfn arm.
The 1.4x in the table above and the 86s in the earlier standing list were taken
⚠️ on battery and are not comparable term for term; 1.91x is the re-measurement.

**The ratio falls with degree, and what is left is the other half.** At degree 11
the 11.2s splits 5.1s for `s → J` and 6.1s for `J → s` — so the direction just
added is now the cheaper one, and the forward fill through `macdonald_j` in the
monomial basis is what a further pass would have to take. Marshalling is not the
story: of the 5.1s, 4.8s is inside symfn and 0.3s is the 1426 fraction-field
cells being built in Python.

**Printed forms had to be matched, not just values.** ℚ(q,t) does not
canonicalize the sign of a fraction, and each `1 − qᵃtᵇ` factor contributes a
`−1`, so a denominator with an odd number of factors comes back negated relative
to what Sage computes for itself — the same element, printing differently, and
nine `sf` doctests failed on it. `_mac_cell` now normalizes to a positive leading
coefficient, **except** on the diagonal, where Sage forms `1/c_λ` by inverting a
polynomial directly and never reduces it. Both halves are needed: the rule was
found by dumping all 233 cells through degree 7 in each arm and diffing the
printed strings, and with it the whole `sf` suite passes with the backend
installed and without it.

### Three traps in `_invert_morphism`, all found the hard way

- **It recomputes the known direction** unless *both* caches already hold the
  degree. Pre-filling one and then calling it does nothing but add work; the
  first Macdonald attempt measured *slower* than no change at all. Hand the
  table over as the `to_other_function` it calls, or bypass the method entirely.
- **Its triangular branch is `O(p(n)³)` like its dense one.** The flag is a
  constant factor, not an asymptotic one — Macdonald `J → s` *is* triangular and
  setting the flag measured 10.4s against 10.6s. Not worth the diff.
- **Comparing printed forms will lie to you.** Switching Macdonald to the
  triangular branch appeared to change the answer; it had not. The fraction
  field normalizes the sign of numerator and denominator together, so the same
  element prints two ways. Compare values.

### The control arm was symfn, and every ratio read 1.0x

⚠️ **The A/B harness stopped being an A/B, silently.** `check_backend.py` and
`bench_backend.py` build their control arm by calling `classical.init()` — and
`classical.init()` on the Sage branch now *defaults to symfn* whenever the wheel
is importable. Both arms were symfn. The benchmark passed, the checker reported
0 mismatches, and neither meant anything.

What gave it away was not a failure but a **shape**: every one of 60 ratios came
back between 0.88x and 1.25x, including `s → m`, which is a 15x row. A harness
that agrees with itself to within noise on a case known to differ is reporting
that it compared nothing.

The fix is `SAGE_DISABLE_SYMFN`, read by `sage.features.symfn.Symfn` — so it
reaches the conversion table, the Hall–Littlewood, Jack and Macdonald caches,
the character bases, `expand`, the monomial product, `SemistandardTableaux` and
the Schubert polynomials alike, all of which choose through the one feature. It
has to be set **before the process starts**, because `classical` fills its table
at import and `Feature.is_present` caches; both harnesses set it in the child's
environment. Putting it in the feature rather than in `is_available` is what
also makes the doctest framework skip `# optional - symfn` tests, instead of
running them against the backend they are not testing.

Two corrections to what this file previously recorded. The **"9547 comparisons,
0 mismatches"** line from the session before this one cannot be trusted: it was
produced by the harness in this state. It is re-established here at **8647
comparisons of degree 8 and 13838 of degree 9, 0 mismatches**, with a control
arm that is verifiably Symmetrica. The per-route numbers in
[transitions.md](transitions.md) are *not* affected — their two arms differ by
up to 20x, which a self-comparison cannot produce.

### The standing list, ranked by what was measured

Everything below was timed on this machine with the symfn backend installed and
`SAGE_DISABLE_SYMFN` marking the control — so these are the walls that *remain*.
⚠️ On battery except the Macdonald rows, which are on AC and say so.

1. **`s → s̃`, the character-basis floor.** With the peel intercepted, the whole
   of `h → ht` at degree 16 is one `schur_to_ht`, and inside it
   `schur_to_st_row(ν)` runs a full `s → p` and back per ν: `p(n)²` character
   work, once per Schur term. Orellana–Zabrocki give `r_{νμ}` directly.
   Everything above it is now free — the *second* conversion at a degree costs
   0.006s where Symmetrica's peel costs 0.983s.
2. **The h → p and e → p generator table is rebuilt per call.** The only two of
   the six direct routes not ahead of Symmetrica (0.78–0.98x), and the reason is
   the one thing `thp.c` does that this does not: cache the table. It is
   ring-dependent, so the `htilde_cached` rule — cache at a concrete ring and
   convert — is the shape of the answer.
3. **Macdonald `J → s`, now the larger half.** With `s → J` supplied, the
   degree-11 fill is 6.1s forward against 5.1s back (⚠️ AC). The forward table
   builds `J_μ` in the **monomial** basis one shape at a time and converts each
   row; a Schur-native route, or the `S` basis and its creation operators from
   the symfn side, is what would move it.
4. **Jack `Q` and `J`, and Macdonald `P`/`Q`/`H`/`H̃`.** All are defined off the
   two bases that now have fast caches, so they may already be fixed — unmeasured.
5. **The `ht` matrix rule's decline.** `h̃_λ · h̃_μ` refuses on long partitions
   and falls back to the Schur route. A second route for the long case would
   close the last slow corner of that basis.
6. **`itensor` at large degree.** 0.60s at n = 21 and growing; symfn has a
   single-coefficient Kronecker query that Sage has no equivalent for, but the
   whole-product path is already respectable and this is the weakest row here.
7. **The `check_*.py` scripts that use Sage as an oracle do not all disable the
   backend.** `gen_sage_oracle.sage` and `check_qt_kostka.py` now refuse to run
   without `SAGE_DISABLE_SYMFN`; `check_macdonald.py`, `check_jack.py`,
   `check_hl.py` and the rest are exposed to the same self-comparison and have
   not been audited. Nothing says a run of theirs was honest except the date it
   was taken. The committed fixture is clear — regenerating it under the guard
   changed nothing but the rows being added — so this is about future runs.

Off the list because they were measured and are fine: `m → s`, `p → s`,
`e → s`, `scalar`, `omega`, `LLT`, and Macdonald `H̃`; the whole family
`p → h`, `p → e`, `h → e`, `e → h`, which went from *losing* to Symmetrica to
1.1–1.5x ahead ([transitions.md](transitions.md)); and — new — the Macdonald
`s → J` direction, which was item 3 and is now 4.9x at degree 11.

### One call per pair, not two

The routing gap had a second home, in the adapter. `_convert` composed
`_TO_SCHUR[src]` and `_FROM_SCHUR[dst]` in Python, which forces the hub no
matter what the kernel can do — a direct `p → h` in Rust is unreachable if the
caller has already asked for `p → s`. So the per-pair entry points that invited
that composition are gone from the adapter's path, replaced by two that name the
pair:

| entry point | covers |
|---|---|
| `convert_terms(a, src, dst)` | any pair with an integral target |
| `convert_indexed(a, src, dst)` | the same, output partitions as indices |
| `to_power(a, src)` | any source into the power sums, coefficients rational |

`to_power` replaces `schur_to_power`, which could only say one thing. Nothing
released depends on the old name.

## The Python boundary's integer ceiling — decided: compute-and-escalate

Two corrections to what this file previously implied. **`gmp` and `python` do
compose** — `maturin build --features "gmp,python"` produces a wheel; that note
was stale. But enabling `gmp` changes *nothing* about the Python API, because
`python.rs` hardcodes `i128` in 20 places and every entry point builds
`Schur<i128>` or the i128-backed `Rational`. The GMP wheel is behaviorally
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
showed the same anomaly, which is what localized it to position).

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
direction. The guard tests must therefore be serialized against each other,
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

**Both misclassifications came from matching a Symmetrica name to a symfn name
and stopping there**, without checking what the caller does with the result.

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

It is recorded here rather than as a footnote because it applies to all 31
entry points symfn will not cover: **displacement does not have to mean
removal.**
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
entry points. Recorded because the correction is the interesting part — it
arrived through the convenience constructor, so "does it panic?" was the wrong
question to audit by.

This closes delta 2 of [../policies/python.md](../policies/python.md), whose
own statement of the gap said partition arguments met a `PanicException` —
the correction is recorded there too, since the item is that file's.

Probed by building the cdylib and driving it from CPython, five clusters
produced a `PanicException` from correctly-typed input. All five are now
`ValueError`. P8 already owned the rule that they must; `failure.md` gained
R11 for the mechanism half P8 defers — which wrong answer is being guarded
against, and how to tell it from a right one.

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
  (commit 2186e3f) and deliberately keeps its trait name and shape, so the two
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
  zero" is labeled a convention, not a theorem. **These five now raise**
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
stricter, which is the license R11 already grants: totality is what a composing
Rust caller needs, and a typed refusal is what a foreign caller needs, and they
are not the same requirement.

Boundary script: 116 malformed calls, plus the theorem-zero pins and the
padding pins.

## Sorting the surface: 98, not 91, and 27 of them uncallable by keyword

`docs/release-readiness.md` Phase 2's Python half, and the membership half of
[../policies/python.md](../policies/python.md) delta 1. Two things came out of
it that the plan did not anticipate, and one that it did.

**The count was wrong every time it was written, in the same way.** This file's
policy said 91, the checklist has said 87 and 92, and all three came from
grepping `#[pyfunction]`. That undercounts twice over. Two of the attributes
are inside `out_of_schur!` and `into_schur!`, which expand to nine conversion
entry points between them — so the ten conversions out of and into Schur cost
three attributes, not ten. And one apparent match is the string `#[pyfunction]`
inside a doc comment, which is why a naive `grep -c` says 92 where a parser
says 91. The authoritative list is the `#[pymodule]` registration block;
`dir(symfn)` agrees with it at **98**. Anything counting this surface counts it
there.

**27 of the 98 named their first argument `lambda`.** That is a Python keyword,
so `symfn.jack_p(lambda=[2, 1])` was a `SyntaxError` rather than a call, and
every one of those functions was positional-only in practice while advertising
a parameter name. All 27 are now `la` — the same transliteration rule the
surface already used for μ→mu and ν→nu, and already reached for once in
`class_algebra_coefficient(la, mu, nu)`. No caller changed, and none could
have: a keyword call did not parse, so every script reaches them positionally.
The argument names inside `same_degree`'s message were renamed with them, since
that message quotes the name the caller typed and is wrong the moment the
parameter moves without it.

The defect had been there since the boundary was written and no check could
have found it, because every existing check *calls* the functions and calling
them positionally works. What found it was writing `symfn.pyi`: a stub cannot
express `def jack_p(lambda: list[int])` any more than the call could. That is
the argument for the stub file being a gate rather than a courtesy — it is the
only artifact that has to name every parameter, so it is the only one that can
notice a name is unusable.

**The harness-only category came out empty**, and that is the part worth
recording as a result rather than as a gap. P10 anticipated a set of probes
that exist only for `scripts/check_*.py`, to be underscore-prefixed and left
out of the stubs. There are none. Every entry point does a whole-object
operation with a caller who can be named, which is P2 having been enforced
since the file was written — the fine-grained accessors the category exists to
absorb never accumulated in the first place. The two candidates that looked
like probes are not: `schubert_monomial_mass` documents a caller who wants the
out-of-family flag before committing to a product no engine finishes, and
`clear_caches` is what anyone timing a run needs, cache-clearing being the trap
`scripts/compare_sage.py` documents on Sage's side. The check scripts therefore
needed no edit, which is the gate that item named.

**Gates.** `scripts/check_python_stubs.py` fails on three kinds of drift —
exported without a stub, stubbed without being exported, and a parameter name
that differs between the two. The third is the quiet one: PyO3 exports every
argument as keyword-callable, so a parameter name is contract, and a stub
saying `mu` where the module says `nu` type-checks a call that fails. All three
were tested by breaking the file deliberately before it was committed passing.
It needs no Sage and no maturin, taking the artifact `cargo build --features
python` leaves behind, the way the boundary script does.

**The other half of delta 1 — the per-function docstring sweep — is recorded
below, under "The docstring sweep".** When this paragraph was written it was
untouched, and the stub summaries were the existing first sentences.

Reading all 98 of those sentences while generating the stubs did turn up one
outright false one, which is a preview of what the sweep is for:
`schubert_monomial_mass` described its return as **"in microseconds"**, when it
returns `schubert::dimension(u) · schubert::dimension(v)` — a count of
monomials with no time in it anywhere. Corrected. Nothing checked it, and
nothing could have: the summary sentence is the one part of a docstring no test
executes, which is why it is the part a doctest example is worth adding to.

### The merge with the Cython branch: 108, and eight more renames

Recorded because the numbers above are now historical. The Cython/Sage-backend
branch was developed in parallel and merged after this sort; it adds eleven
entry points — the `st` and `ht` character bases, the six direct routes among
h, e and p, `convert_terms`, `schur_in_macdonald_j` — and replaces
`schur_to_power` with the generalized `to_power(a, src)`. **The surface is
108.** Nothing was public, so the replaced name was taken rather than kept
alongside.

Two things the sort's own artifacts caught, which is the argument for having
built them:

- **Eight more entry points named `lambda`**, six of them the ones this branch
  had renamed and git had resolved back to the other side, plus
  `reduced_kronecker` and `reduced_kronecker_product`, which are new. Merging
  the two branches textually left the surface *half* renamed — some functions
  keyword-callable and some not — and one function with a signature from one
  side and a body from the other, which is the only reason it failed to
  compile rather than shipping. All eight are `la` now.
- **`Key` made the stub types wrong in a way no name check would see.** The
  Cython branch changed outbound partitions from lists to Python tuples, so
  every regenerated return type moved from `list[int]` to `tuple[int, ...]`
  while every *argument* stayed permissive — `Key` extracts from any sequence.
  The stub file now carries that asymmetry explicitly, with `*Arg` aliases for
  the inbound halves, because it is a real property of the boundary and not an
  oversight: strict about what it promises, permissive about what it accepts.
  Verified by feeding a result straight back in.

`scripts/check_python_stubs.py` compares names, arity and parameter names, not
types — so it caught the first of these and not the second. The second was
found by regenerating rather than by checking, which is the limit of that gate
and worth knowing before trusting it further than it goes.

## The docstring sweep: 108 entry points, 132 executed examples

Delta 1's other half, 2026-08-08. Every `#[pyfunction]` now carries a
`# Raises` section naming the exception and the requirement it names, and at
least one example that runs.

**The runner came first, and it earned itself on its first run.**
`scripts/check_python_docs.py` loads the cdylib `cargo build --features
python` leaves behind — the same trick `check_python_stubs.py` and
`check_python_boundary.py` use, so it needs neither Sage nor maturin — and
executes every ` ```text ` fence containing `>>>` as a doctest with `symfn`
bound to the loaded module. The fence is what makes the arrangement work in
both directions: cargo's doctest runner skips a text fence rather than trying
to compile it as Rust, and running one fence at a time keeps the prose
between fences from being parsed as expected output.

Before a single entry point had been touched, it failed on the module
docstring's own example, which showed

    >>> symfn.schur_multiply([([2], 1)], [([1], 1)])
    [([2, 1], 1), ([3], 1)]

— lists where the boundary returns tuples. That is the exact claim P1 makes
about the outbound type, printed wrong on the first surface a Python caller
reads, and nothing had executed it since it was written. It is the same
failure mode as `schubert_monomial_mass`'s "in microseconds" above, and the
same argument for the same fix.

**Order got one home.** P6 asks every list-returning entry point to state its
order. Repeating it 108 times would have been 108 chances to drift, so the
rule is stated once in the `#[pymodule]` doc — increasing lexicographic by
support, no zero coefficients, no repeated key — and an entry point states
its own order only where it differs. Three do, and finding them was most of
the value:

- `expand_alphabet` groups its pairs by the monomial term they come from,
  and **the order within a group is not contract**. It looks sorted and is
  not: `s_2` in two variables comes back `(1,1), (2,0), (0,2)`.
- `partitions`, and therefore every table indexed by it, is *reverse*
  lexicographic — the opposite of the element order. A caller who assumed
  one order held everywhere would read `character_table` transposed.
- `stanley_table` runs λ, μ, ν each in `partitions` order, not
  lexicographic; `llt_kl_column` and `llt_gtilde_table` omit their zero
  entries, so neither list is as long as `partitions(n)`.

**What the runner cannot check is the part that matters most.** An example
that every rival convention also satisfies pins nothing. That judgment is not
mechanizable and is stated in the script's own docstring rather than left
implied. The examples written to be distinguishing, and what each separates:

| entry point | the value | what it rules out |
| --- | --- | --- |
| `st_to_schur` | `s̃_(2) = s_2 − 2·s_1` | a homogeneous reading of `s̃` |
| `schubert_expand` | `S_{132} = x_1 + x_2` at `(0,1)` | 1-based exponent vectors |
| `schubert_multiply_variable` | `i = 1` means `x_1` | Symmetrica's 0-based index |
| `character_table` | the `2` at the row's end | a class-major table |
| `semistandard_tableaux` | `(2,0,1)` ≠ `(2,1)` | μ read as a partition |
| `qt_kostka` vs `macdonald_ht` | `K_{(2),(11)} = t`, `K̃ = 1` | the modified form |
| `hall_littlewood` vs `_p` | `Q'_{11} = s_11 + t·s_2` | the other normalization |
| `jack_p`/`q`/`j` | empty denominator and scale 1 | the other two forms |
| `zonal(la, True/False)` | `2m_11 + 3m_2` vs `⅔m_11 + m_2` | Sage's `P^{(2)}` for `[GJ]`'s `Z_λ` |
| `llt_kl_column` | alternating signs | the `q` grading, which is positive |
| `chromatic_from_llt` | no monochromatic term | `llt_graph`, which keeps it |
| `k_core_quotient` | runner 0 first | the reversed component order |

**Two mechanical notes.** The two conversion macros now take the entry
point's doc comment as a macro argument (`$(#[$doc:meta])*`), because nine
generated functions sharing one doc would state nine different conventions
badly or none at all. And one example was shortened for width rather than
content: `nabla_power`'s argument moved from `s_11` to `s_2` so the expected
output fits 80 columns.

**What is still not done.** The runner is not in CI, because there is no CI
([release-readiness.md](../release-readiness.md), Phase 0), and
`scripts/preflight.sh` cannot run it — preflight builds default features and
this needs `--features python`. It sits with `check_python_stubs.py` and
`check_python_boundary.py`, the two other gates that need the cdylib and are
run by hand after touching `src/python.rs`.

## The convenience layer, and what building it found

`policies/python.md` delta 5 is closed: the wheel is a mixed layout with a
pure-Python layer over the contract calls. `Sym` carries a basis tag and refuses
to combine two elements that disagree; the parameter families are namespaces
returning coefficient objects that print readably and specialize; `Schub` is
the permutation-keyed sibling. Nothing in the layer computes — every method is
a contract call with bookkeeping, checked as such below.

The layout change is Phase 5's, not a side effect: `pyproject.toml` now exists
with the maturin backend and `python-source = "python"`, so the compiled module
lands at `symfn.symfn` and the supported names stay flat at `symfn.*`. The
stubs moved from the repository root into the package beside the module they
describe, where a type checker looks for them; at the root they typed nothing
once the layout became mixed.

### Three defects the gates found, in the order they appeared

**A convenience namespace silently deleted a contract entry point.** The
Hall-Littlewood namespace was called `hall_littlewood`, and `symfn.__init__`
imports it after `from .symfn import *`. `symfn.hall_littlewood` was therefore
the namespace and not the entry point of that name, which is a break in a
surface that is supposed to freeze hardest of anything in the tree (P10). The
namespace is now `hl`, and `check_convenience.py` asserts every contract name
still resolves to the contract object — the check exists because of this, and
it is the cheapest of the three to have missed.

**The LLT wrappers were tagged Schur; the entry points return monomial.**
`llt_g`, `llt_gtilde` and `llt_h` all return in the monomial basis, and only
`llt_g` says so in its first line — `llt_schur` naming itself is what makes the
others' basis inferable rather than stated. A wrong basis tag survives every
value check, because the values are right; it fails only a check that asserts
the tag. This is the failure mode `CLAUDE.md` names as the house one, met in a
new place: the convenience layer can now mislabel a correct answer, which the
contract layer could not, because the contract layer had nothing to label it
with. `check_convenience.py` asserts the tag of every family.

**The same wrappers also claimed two parameters where the family has one**,
found the same way and worth recording beside it. LLT rows share the
`(q_exponent, t_exponent, coefficient)` encoding with the Macdonald operators
and use only the first slot, so wrapping them as a `(q, t)` element produced
correct values under a signature demanding a `t` that appears nowhere:
`llt.G(shapes).at(q=2)` raised, and `parameters` said `('q', 't')`. The
projection now builds a one-variable `Poly` and raises if a row ever carries a
nonzero `t` exponent, so a family that grew a second parameter would fail
loudly rather than have it silently dropped. Both defects are the same shape —
a label on a right answer — and neither a value check nor a doctest of a
`repr` can see either.

**`convert_terms` listed a target basis it rejects.** The `dst` arm reaches
five bases; the shared error message named six, `powersum` among them, so a
caller who asked for it was told it was expected and refused in the same
sentence. The rustdoc had it right the whole time — "the conversions that
divide are `to_power`'s, which is why they are not reachable here" — which is
the argument for the message being written against the match arm rather than
against the family. It now names the five and points at `to_power`.

### What the layer is checked by

Four gates, all Sage-free, all in CI, run together by
`scripts/preflight_python.sh`:

| gate | what it holds | size |
| --- | --- | --- |
| `check_convenience.py` | every method equals its contract composition; the families hit their classical limits; no name shadows a contract name; mixing bases raises | 2177 checks |
| `check_convenience_docs.py` | every docstring example runs, and every public item has one | 74 examples over 78 items |
| `check_docs_complete.py` | every supported name reaches a rendered page | 187 names |
| the `wheel` CI job | `pip install symfn` imports and computes with no Sage on the path | — |

The first is the one that makes P4 checkable rather than aspirational.
"Convenience computes nothing" is a claim about equality, so the check runs each
method against the contract sequence it claims to be, over every partition to
degree 6 and every basis, rather than at one shape.

**The degenerations are what check the conventions from Python.** `P_λ(x;q,q) =
s_λ`, `P_λ(x;α=1) = s_λ`, `Q'_λ(x;0) = s_λ`, `Q'_λ(x;1) = h_λ` and `K_{λμ}(1) =
K_{λμ}` are theorems that each fail under a `q ↔ t`, `α → 1/α` or `t → 1/t`
twist, and `at` is what makes them expressible as a Python equality. That is the
argument for the coefficient types carrying an evaluation map at all: without
it the parameter families cross the boundary as rows nothing on this side can
check.

### The rendered documentation

The site is Sphinx with MyST, under `docsite/`, published by Read the Docs,
which builds the wheel first — the reference is generated from the objects, so
`help()` and the website cannot drift.

**The docstrings are Markdown and Sphinx's autodoc feeds RST to docutils**,
which is the one real obstacle and the reason for `markdown_docstrings` in
`docsite/conf.py`. Four constructions in the house form need translating, and
three of them were found by the build failing rather than by reading:

- a ` ```text ` fence becomes a literal block, `pycon` when it holds a doctest;
- `# Raises` becomes a rubric, not a section, so it does not enter the page's
  heading tree;
- a single-backtick span becomes a double-backtick literal — and needs RST's
  escaped space after it when a word runs on, which ``` `int`s ``` in the module
  doc does;
- a bare `|` is escaped, because `|λ|` for the size of a partition is a
  substitution reference to RST and five entry points write it.

The heading rule had a bug worth recording because it looked like a docstring
problem rather than a regex one: the pattern ended `\s*$`, and `\s` matches
newlines, so it ate the blank line after every heading and docutils reported
"explicit markup ends without a blank line" once per entry point — 60-odd
warnings that all had one cause. The class is `[ \t]*$`.

The build runs with warnings as errors, and `check_docs_complete.py` compares
`symfn.__all__` against the `objects.inv` Sphinx writes. Sphinx does not check
this itself: it warns about a reference that does not resolve and says nothing
about an entry point no page documents, which is the hole that actually opens
when the surface grows. All 108 contract entry points and 79 convenience names
are on a page.

### Annotating the layer, and what a type checker found

The convenience layer shipped unannotated for one commit while the wheel
carried `py.typed` and the `Typing :: Typed` classifier. That combination is
not a small overclaim: `py.typed` tells a checker to read the package's inline
types, and finding none it infers `Any` for every `Sym`, `Param` and `Schub`
call, so the classifier promised the opposite of what a user got. The layer is
now annotated throughout, with `ruff`'s `ANN` rules holding completeness and
`mypy --strict` holding correctness — both in `preflight_python.sh` and CI.

**`symfn.pyi` could not be read by a checker at the floor it advertises.** It
imported `TypeAlias` from `typing`, which arrived in 3.10, while
`requires-python` is `>=3.9` to match the `abi3-py39` build. A checker running
as 3.9 stopped at the import. The annotations were never needed — a bare
assignment is a type alias — so they are gone, and this is the second defect
found by pointing a tool at the stub rather than reading it
(`check_python_stubs.py` found the 27 `lambda` parameters the same way).

**Two modeling decisions worth keeping.** `Mapping` is invariant in its key, so
a permissive `Mapping[PartitionArg, Coefficient]` is *rejected* for a
`dict[Partition, Coefficient]` rather than accepting more — the opposite of the
intent. The mapping branches of `TermsArg` therefore name `Partition` and `int`,
which is also just true: a mapping key has to be hashable, so a support arriving
as a key is already a tuple or a bare integer, and a list can only reach the
constructor through the pairs form. And `Basis` is a `Literal`, which is P7 as
far into the type system as Python reaches; the one narrowing from `str` sits in
`check_basis`, where the validation is.

**One contract went from documented to enforced.** `Schub.__init__` said it
raises `TypeError` unless every coefficient is an integer, and did not check —
a `Fraction` was accepted and stored. mypy found it as a type error on the
assignment. It now raises, which is what the docstring already claimed.

The pass also closed the asymmetry between the two factories: `_partition` took
a bare `int` and `_permutation` did not, so `_SchubFactory.__getitem__`
compensated by wrapping non-tuples itself while `_Factory.__getitem__` did not.
Both normalizers take an `int` now and both factories read identically. No
behavior changed — `X[1]` was already the identity — but the compensation was
the kind that stops being applied the next time someone adds a call site.

### Filling out the wrappers, and one the types could not hold

Sixteen more entry points reached the convenience layer: the Delta and
Macdonald operator algebra (`nabla_power`, `delta_prime_e`, `delta_ek`,
`delta_prime_ek`, `theta_ek`, `big_pi`), Jack's `structure_constant`, LLT's
Schur-basis `schur`, `Htilde` and `min_inv`, `Sym.principal_specialization_q`,
and the Schubert half that was thinnest — `multiply_variable`,
`divided_difference_perm`, `pairing`, `dimension`, plus the module-level
`stanley_schur` and `from_polynomial`. 38 of the 108 entry points had a wrapper
before; 54 do now, and the rest stay reachable flat at `symfn.*`, which is the
documented answer rather than a gap.

**`jack_norm_j` was wrapped wrong and then withdrawn**, which is the finding
worth recording. It returns `⟨J_λ, J_λ⟩_α` as a factored list of atoms — a
*numerator*, handed over factored for the same reason everything else here is —
and `AlphaFrac` is a fraction type whose atoms are its *denominator*. Wrapping
it as `AlphaFrac([1], atoms, 1)` produced `1/(α(α+1)·2α)` where the answer is
`α(α+1)(2α+1)`: the reciprocal, printed confidently. The repr is what showed
it, since the value was never evaluated in a check.

The type gap is real and the entry point stays flat until it is closed. A
factored *numerator* has no home in the coefficient types, and the three ways
out are each a decision rather than a wrapper: give `AlphaFrac` a factored
numerator slot, expand the product into the dense numerator (polynomial
arithmetic, further from P4's carve-out than substitution is), or add a
product type for one function. None of them belongs in a pass whose subject was
wrapping existing calls.

**The cross-checks the new wrappers brought** are the ones whose identity is
easy to get wrong rather than whose values are: `∇^1 = ∇`, `Δ'_{e_k} e_n`
against `delta_prime_e(k, n)`, `Θ` raising the degree where `Δ` preserves it,
the graded specialization summing to the plain one at `q = 1`, and `expand` and
`from_polynomial` being inverse over all of `S_4`. `check_convenience.py` is at
2906 checks, from 2183.
## The release pipeline, and why a tag does not publish

The repository got a home — `github.com/mwhansen/symfn` — which unblocked the
two URL fields that had been held empty on purpose, and the four packaging
items Phase 5 still owed.

**The wheel matrix is 14 artifacts**, the full `rpds_py` platform set that
`docs/sage-packaging-audit.md` argued for, built by `maturin-action` in
`.github/workflows/release.yml`. The audit had expected the exotic Linux
architectures to need QEMU legs. **They do not, and the reason is worth
recording because it changes the cost of the matrix by an order of magnitude:**
an `abi3` build needs no interpreter for the target, so maturin cross-compiles
armv7l, ppc64le, s390x and i686 inside its manylinux containers at the price of
an ordinary compile. Emulation would have bought *testing*, not building. What
that leaves is the honest gap — a cross-compiled wheel is built and never
imported on the platform it targets — and `docs/support-tiers.md` states it as
accepted exposure rather than leaving it to be discovered. The native legs
(macOS x86\_64 and arm64, Windows x64) do run an import-and-compute step, and
the exposure is narrow precisely because `abi3` links the stable ABI rather
than a version's internals.

**The sdist vendors its dependencies**, at 4.4 MB against 1.4 MB for the wheel.
`scripts/build_sdist.sh` writes `vendor/` and `.cargo/config.toml`, calls
maturin, and removes both; neither is tracked, and `pyproject.toml`'s
`[tool.maturin] include` is what carries them into the artifact when they exist
and matches nothing when they do not. `cargo vendor --locked` is what ties the
vendored set to the committed `Cargo.lock` rather than to whatever resolves on
the day.

**The trap, paid for once.** `.cargo/config.toml` was already tracked and
carries the `-undefined dynamic_lookup` link arguments a bare
`cargo build --features python` needs on macOS — maturin supplies its own, so
nothing in the wheel or sdist path notices they are gone. The first draft of
`build_sdist.sh` wrote that file outright to add the redirect, which deleted
the flags, and the next `scripts/preflight_python.sh` failed at the link step
with a page of undefined `_Py*` symbols that named neither the file nor the
script. The script appends and restores from a backup now, with the restore on
an EXIT trap so an interrupted run also puts it back, and `.gitignore` carries
`vendor/` and the backup as the second net. The generalization: a build script
that writes a configuration file should assume the file is someone else's.

The measurement that mattered: **offline is asserted, not simulated.**
`scripts/check_sdist_offline.sh` sets `CARGO_NET_OFFLINE=true` and pip's
`--no-index --no-build-isolation`, so a dependency the vendoring missed is a
hard error naming the crate. Cutting the runner's network would have tested the
runner. Verified locally on macOS arm64: the sdist built
`symfn-0.1.0-cp39-abi3-macosx_11_0_arm64.whl` and computed, with both fetchers
refusing. It runs in the ordinary CI workflow rather than only at release,
because the way it breaks is a new dependency landing in `Cargo.toml` — an
ordinary commit, not a release-day event.

**A tag builds; it does not publish.** Phase 5 asked for publish-on-tag and
this is deliberately not that. A `v*` tag builds all fifteen artifacts and
attaches them to a GitHub Release; reaching PyPI or crates.io takes a
`workflow_dispatch` that names the registry, behind a GitHub environment whose
required reviewer holds even when the dispatch input is wrong. Two reasons, and
the second is the one that decides it: testers install from the Release page
before the name exists on any registry, and a registry publish is the only
irreversible step in this pipeline — a yanked version number can never be
reused, so a wrong 0.1.0 costs the number permanently. Making the not-yet state
structural is cheaper than remembering it, and the gate stays useful after the
first publish.

Both registries authenticate by OIDC — PyPI's trusted publishing and
`rust-lang/crates-io-auth-action`, which exchanges the run's identity for a
short-lived token — so no long-lived secret is stored in the repository for
either half. `cargo publish` runs `--dry-run` first and `--locked` in both
passes, so the crate resolves to the same lockfile the wheels were built from.
Measured: the crate packages to 230 files, 1.0 MiB compressed, well inside
crates.io's limit.

### What is still open

- The round-trip half of the Sage-free suite (Phase 5) is still not written:
  a value handed in comes back out intact, at the widths and shapes P1
  promises.

  **A correction to what this entry said above when it was written.** It
  listed the round-trip half as "values computed in Rust and asserted from
  Python" and called the existing gates weak for comparing Python to Python,
  "so a kernel defect would pass them". Catching a kernel defect is not what
  these tests are for. `docs/policies/validation.md` owns that evidence and
  the Rust suites and committed fixtures carry it; a second oracle reached
  through PyO3 would be a slower, narrower copy with the boundary in the way
  of every failure it reported. These gates owe the boundary and the layer
  above it, and Python-against-Python is the right instrument for that.
- **`docs/` is not published, and that is now the decision rather than an
  open question.** The rulebooks and this record are the tree's internal memory
  — written for the next session, dense with dead ends and measurements — and
  the site is an outside reader's manual. They are different documents for
  different readers, and publishing the first as if it were the second would
  mislead about which parts are contract. Revisit only if an outside reader
  asks for something the site cannot say without them.
- The convenience layer wraps most of the families' constructors, not all 108
  entry points. The rest stay reachable flat at `symfn.*`, which is the
  documented answer rather than a gap. The `*_table` and `*_column` bulk family
  is the largest group without one, and `jack_scalar`, `chromatic_from_llt` and
  `llt_e_expansion` are the individual ones most likely to want it next.
- `jack_norm_j` has no wrapper because the coefficient types cannot hold a
  factored numerator; the section above states the three ways out and why none
  of them is a wrapper.
- **Nothing in the release pipeline has run yet.** It is written and its
  scripts are verified locally, but no tag has been pushed, so the fourteen
  cross-compiled legs are untested against GitHub's runners. The three settings
  it depends on are outside the tree and have to be made once on the
  repository: a PyPI trusted publisher for `symfn` naming
  `release.yml`, the same on crates.io, and the `pypi` / `testpypi` /
  `crates-io` environments with required reviewers. Until those exist the
  publish jobs fail at authentication, which is the correct failure.
