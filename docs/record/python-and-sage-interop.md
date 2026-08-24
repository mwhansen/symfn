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

> The adapter this section and the next several describe lived at
> `scripts/sage_backend.py`, with its per-term loop at `scripts/symfn_cy.pyx`.
> Both were deleted once the adapter moved into Sage itself as
> `sage/libs/symfn/`; those paths do not resolve, and everything below is the
> history of how the arrangement was arrived at. `scripts/check_backend.py`
> still exists and now drives Sage's copy.

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
is marshalling, not the Rust core. On the marshalling alone it is worth
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

> **Both of those sentences describe an arrangement that no longer exists.**
> The loop is `sage/libs/symfn/terms.pyx` now, registered in Sage's
> `src/sage/libs/meson.build`, so it compiles whenever Sage does and
> `backend.py` imports it unconditionally. There is nothing to fall back from,
> and the optional-import guard the paragraph describes was removed with the
> files it lived in. The measurement below is why the loop exists at all and
> still stands.

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
- **Printed forms disagree on equal values.** Switching Macdonald to the
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

## The coverage census, and the 30 entry points nobody calls

A scan of 3053 sagelib source files against the 66 entry points
`sage/libs/symmetrica/all.py` exports, at `sagemath/sage` commit `09472ff`
(10.10.beta7, 2026-07-26), found that **Sage reaches 36 of the 66, from six
files**: `combinat/sf/classical.py` (the 20 basis conversions), `sf/sfa.py`
(the five `compute_*_with_alphabet`), `sf/hall_littlewood.py`,
`sf/monomial.py`, `combinat/tableau.py`, and `combinat/schubert_polynomial.py`.
All 36 are now computed by symfn.

This is the census behind the claim in
[../sage-backend.md](../sage-backend.md), and the claim is narrower than it
looks: **covering what Sage calls is not covering what Symmetrica exports.**
The other 30 are public API a user can reach with
`from sage.libs.symmetrica.all import ...`, and no sagelib code path touches
them:

```
bdg  chartafel  charvalue  compute_schur_with_alphabet_det  dimension_schur
dimension_symmetrization  gupta_nm  gupta_tafel  kostka_tafel  kranztafel
mult_schur_schur  ndg  newtrans  odd_to_strict  odg  outerproduct_schur
part_part_skewschur  plethysm  q_core  random_partition  scalarproduct_schur
schur_schur_plet  sdg  specht_dg  start  strict_to_odd_part
t_POLYNOM_ELMSYM  t_POLYNOM_MONOMIAL  t_POLYNOM_POWER  t_POLYNOM_SCHUR
```

(`start` is the library initializer, called by `all.py` itself.)

symfn covers much of that list incidentally — `plethysm`, `mult_schur_schur`,
`outerproduct_schur`, `part_part_skewschur`, `dimension_schur`,
`charvalue`/`chartafel` and `newtrans` all have direct equivalents. What is
genuinely uncovered is the representation-theory group: `bdg`, `sdg`, `odg`,
`ndg`, `specht_dg`, `dimension_symmetrization`, `kranztafel` (wreath products),
`gupta_nm`/`gupta_tafel`, `q_core`, `strict_to_odd_part`/`odd_to_strict`. That
is [README.md](README.md)'s "Beyond the core" territory, and it is exactly the
part Sage never calls.

Retiring those is a **deprecation question, not an implementation question**,
and it belongs upstream: either reimplement them or take them through Sage's
own deprecation cycle. It is the reason "symfn covers every Symmetrica entry
point Sage calls" is the sentence to use, and "symfn replaces Symmetrica" is
not.

## From one consumer site to five

The census above also found that the conversion table — all `sage_backend.py`
displaced — was one of the six files. This chapter is what happened when the
remaining task list was implemented: the four tasks landed, but **two of the
three gap classifications were wrong**, and a fourth gap that scan did not see
turned up in the file it had marked complete.

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
the demotion lands, not never. The exception-fidelity requirement still
applies.

## The boundary raises where it panicked: 5 clusters, 30 entry points

The premise this started from was that `part()` calls `Partition::new`, "which
asserts", so a non-partition reached Sage as a `PanicException`. **`Partition::new`
does not assert.** It normalizes — filters zeros, sorts weakly decreasing — so
`symfn.schur_multiply([([1,3], 1)], …)` silently returned the product for
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

**The runner came first, and it caught a defect on its first run.**
`scripts/check_python_docs.py` loads the cdylib `cargo build --features
python` leaves behind — the same arrangement `check_python_stubs.py` and
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

Five gates, all Sage-free, all in CI, run together by
`scripts/preflight_python.sh`:

| gate | what it holds | size |
| --- | --- | --- |
| `check_convenience.py` | every method equals its contract composition; the families hit their classical limits; no name shadows a contract name; mixing bases raises | 2177 checks |
| `check_convenience_docs.py` | every docstring example runs, and every public item has one | 271 examples over 94 items |
| `check_docsite_docs.py` | every example on the docsite's narrative pages runs, each page one interpreter session | 79 examples over 4 pages |
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
`α(α+1)(2α+1)`: the reciprocal, printed as if correct. The repr is what showed
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

**Both references travel with the artifacts.** `scripts/build_docs.sh` renders
the Sphinx site and `cargo doc --no-deps --all-features`, stages the extension
module the first of those needs, and bundles the pair under a landing page as
`symfn-docs-<version>.tar.gz` — 5.1 MB, 14 MB of Python reference and 12 MB of
Rust reference unpacked, opening from a local filesystem with no server. The
reason it exists is the same stretch of time this whole pipeline is sized for:
Read the Docs publishes the Python half, nothing publishes the Rust half, and
a tester with a downloaded wheel has no reference beside it. It is attached to
every Release and uploaded from every CI run, so a reviewer can also read what
a branch did to the documentation.

Two things fell out of building it. The script denies rustdoc warnings, which
makes CLAUDE.md's standing claim that `cargo doc --no-deps --all-features` is
silent into one that can be false — nothing had checked it before. And the
distributable artifacts had to be renamed to a common `dist-*` prefix: the
publish job downloads with `merge-multiple`, and an unfiltered download would
have handed PyPI `symfn-docs-0.1.0.tar.gz` as a second source distribution
sitting beside the real one. The docs bundle is named `docs` so it cannot match
the pattern; the GitHub Release takes everything.

Both registries authenticate by OIDC — PyPI's trusted publishing and
`rust-lang/crates-io-auth-action`, which exchanges the run's identity for a
short-lived token — so no long-lived secret is stored in the repository for
either half. `cargo publish` runs `--dry-run` first and `--locked` in both
passes, so the crate resolves to the same lockfile the wheels were built from.
Measured: the crate packages to 230 files, 1.0 MiB compressed, well inside
crates.io's limit.

## The first CI run, and what one laptop had been hiding

The tree was pushed to `github.com/mwhansen/symfn` and CI ran for the first
time on 240 commits. Every Rust suite passed on all three platforms, both MSRV
legs passed, the release-profile legs passed, the sdist built offline, and the
wheel imported with no Sage. **Three jobs failed, and all three failed for the
same reason: the gates had only ever met one machine's tool versions.** Nothing
had regressed.

- **clippy, all four feature sets.** 40 errors under clippy 0.1.97 that 0.1.96
  did not have — 27 `doc_overindented_list_items`, 10 `useless_conversion` on
  an identity `Vec<u32>`, one `manual_is_multiple_of`, two
  `wrong_self_convention`. The lint set grew; the code did not change. The
  first three were mechanical and applied from clippy's own suggestion spans
  rather than by hand. The fourth was not: `Md::from_i128` and `from_u128` take
  `&self` as the *modulus* and the value as the argument, and the name is
  deliberately the free functions' name so a call site reads the same either
  way. Clippy reads `from_*` as a constructor and mis-identifies which operand
  is which, so those two carry an `allow` with the reason rather than a rename
  that would break a documented pairing.
- **The Sphinx build, under docutils 0.22.4 against the laptop's 0.21.2.** Two
  warnings, one cause: `docsite/conf.py`'s `CODE_SPAN` excluded newlines, so a
  backticked span that landed across the 80-column wrap kept its single
  backticks while every span around it became double. `stanley_table` and
  `expand_alphabet` were the two that happened to wrap. The rule admits a
  newline now, but still not a blank line — a span running to the next
  paragraph would pair with whatever backtick it found there. Both now render
  as proper inline literals, checked in the built HTML rather than inferred
  from a silent build.
- **mypy, and this one has no version-independent answer.** mypy 2.0 refuses
  `python_version = "3.9"` outright, so the config was a hard error before any
  code was checked; and 2.x narrows `check_basis`'s membership test well enough
  that the `cast` to `Basis` is a `redundant-cast` error, where 1.x needs that
  cast or reports a return-type error. The two cannot both pass. The gate
  targets 2.0, `scripts/preflight_python.sh` checks the major version and skips
  an older checker with a message about itself rather than a failure about the
  code, and CI installs `mypy>=2.0` so the skip cannot silently disable it.
  Dropping `python_version` costs nothing real: ruff's `target-version =
  "py39"`, the `wheel` job on a 3.9 interpreter, and `requires-python` are what
  hold the floor.

The second run left one failure, on macOS alone, and it was not a tool version:
**a `RUSTFLAGS` in the environment replaces `.cargo/config.toml`'s
`target.*.rustflags` rather than adding to them.** The workflow sets
`RUSTFLAGS: -D warnings` globally, which discarded the `-undefined
dynamic_lookup` pair that file carries for macOS, and `cargo build --features
python` failed at the link step with the same page of undefined `_Py*` symbols
the `build_sdist.sh` trap had produced a few commits earlier — a second cause
with an identical symptom, and an error naming neither the workflow line nor
the config file. Reproduced locally with `RUSTFLAGS="-D warnings" cargo build
--features python` before it was fixed, rather than inferred from the log.
**A correction to the fix this paragraph first recorded.** It said every job
that builds that feature empties `RUSTFLAGS` for itself. That does not work, and
the next CI run failed identically: cargo reads the variable as
present-and-empty and still overrides the config file, so `RUSTFLAGS: ""`
changes nothing. Measured directly —

    RUSTFLAGS="" cargo build --features python   # fails to link
    unset RUSTFLAGS; cargo build --features python   # builds

— and YAML has no way to unset an environment variable a workflow already set.
The workflow-level `RUSTFLAGS` is gone; the three pure-Rust lanes (`test`,
`msrv`, `release`) declare `-D warnings` for themselves, and the jobs that
build the `python` feature leave the variable absent. Warnings under that
feature are the `clippy` job's business. The two `RUSTFLAGS: ""` step overrides
in `casts` and `clippy` went with it — they existed only to neutralize the
global, and a mechanism whose reason has expired is worse than none.

The third run was green, and it surfaced one more thing by warning rather than
failing: every `actions/*` step was on the node20 runtime GitHub is retiring.
The majors were bumped to node24 ones — checkout 7, setup-python 7,
upload-artifact 7, download-artifact 8, action-gh-release 3 — after reading
`runs.using` in each action's own `action.yml` rather than assuming the newest
major had moved; `PyO3/maturin-action@v1` and
`rust-lang/crates-io-auth-action@v1` were already node24 inside v1, and
`pypa/gh-action-pypi-publish` is composite and has no JS runtime.

Reading the manifests caught a second, unrelated defect in the release
workflow, which has never run: its macOS and Windows legs asked
`setup-python` for 3.9, and `actions/python-versions` **has no darwin/arm64
build below 3.10**, so the aarch64 leg would have failed at a step nobody had
exercised. Both legs take 3.11 now. Nothing is lost — `abi3` means the build
does not depend on which interpreter is present, and the floor is proven by
ci.yml's `wheel` job, which installs and computes on a real 3.9.

**The first tag, `v0.1.0-rc.1`, and the one leg that never ran.** Sixteen of
the seventeen jobs succeeded on the first attempt: the version guard, all six
manylinux legs, all three musllinux, all three Windows, macOS aarch64, the
offline sdist and the docs bundle. `macos x86_64` sat in `queued` for
twenty-five minutes and never started — `macos-13` is retired and no longer
appears in `actions/runner-images`, so no runner was ever going to claim it,
and the symptom is a queue rather than an error. `macos-14` is deprecated on
the same page. Both legs moved to `macos-15-intel` and `macos-15`, which keeps
them native and therefore keeps the import-and-compute step meaningful.

Because the run never reached `github-release`, no Release was created and no
artifact left the repository, so the tag was re-pointed rather than the
candidate number burned. That is the narrow case where re-pointing is right:
the rule about a version never being reusable is about a version somebody
*received*, and nobody received this one.

The generalization worth keeping: **a green gate is green for the toolchain
that ran it.** Phase 0 said CI existed because every portability claim was a
claim about one laptop. That turned out to be true of the lint and type claims
too, and it took one push to find out.

## The adapter has a home, and two defects only a second user would find

The Sage side is a branch — `mwhansen/sage`, branch `symfn`, 14 commits over
10.10.beta4 — carrying `src/sage/libs/symfn/`, `src/sage/features/symfn.py`,
`build/pkgs/symfn/` and the call-site changes. Phase 5b's open question about
shape is answered by it: the adapter lives *inside* Sage, for the reason that
phase already gave, which is that `terms.pyx` `cimport`s Sage's `Integer` and
so must compile against a specific Sage build.

**A wheel is not enough to run it, and that is the answer to the obvious
question.** The branch has to be built: `terms.pyx` is registered in
`src/sage/libs/meson.build`, so Sage compiles it, and only then does installing
the wheel into that Sage's Python do anything.

Two defects surfaced from asking what a *second* person would have to do, and
neither is visible to the person who wrote it:

- **`build/pkgs/symfn/requirements.txt` asked for `symfn >=0.1.0`, which no
  release candidate satisfies.** `0.1.0rc1` sorts below `0.1.0`, so the
  specifier rejects it even with prereleases enabled — verified against
  `packaging`, not assumed. Upstream's first artifacts are release candidates,
  so the floor is `0.1.0rc1`.
- **`Symfn._is_present()` tested importability, and importability is not the
  question.** Measured against a real leftover install: 5 of the 33 entry
  points `backend.py` calls were missing — `convert_terms`, `ht_multiply`,
  `reduced_kronecker_product`, `schur_in_macdonald_j`, `to_power` — while the
  module imported perfectly well. The feature reported present, so a conversion
  routed into symfn and raised `AttributeError` from inside the basis machinery
  instead of falling back to Symmetrica. It now checks a version floor.

  The subtlety that decides the implementation: the version is read from
  `symfn.__version__`, **not** from the distribution metadata, because the two
  disagree in exactly this case. The stale install kept a `dist-info` claiming
  `0.1.0` while the module it installed predated the attribute entirely, so a
  metadata check would have waved it through.

Both fixed on the branch, with the sf and Schubert doctests passing under
`--optional=sage,symfn` and the Symmetrica answers reproduced under
`SAGE_DISABLE_SYMFN=1`.

**What the reconciliation also turned up, and did not fix.** `classical.init()`
defaults to `is_available()`, so installing the optional package switches the
backend immediately — which is Phase 5c's *second* landing arriving inside its
first. The staging plan's whole argument was that each landing be independently
reviewable and revertible, and "installs but stays off" is a different review
from "installs and takes over". The second finding is settled rather than recorded: `scripts/sage_backend.py`,
`scripts/symfn_cy.pyx` and `scripts/setup_cy.py` were a second, older copy of
the adapter that no preflight ran, and they are deleted.
`scripts/check_backend.py` and `scripts/bench_backend.py` drive Sage's own
adapter now, switching arms through `SAGE_DISABLE_SYMFN` as they already did
and **asserting** which backend answered instead of installing one — an
unverified control being the failure CLAUDE.md records as "a run of ratios all
near 1.0x". Measured after the change: **8647 computations at degree 8, 0
mismatches**, the recorded figure reproduced exactly against `sage/libs/symfn/`
rather than the deleted script. That is what says the deletion cost no
coverage.

## Ctrl-C did nothing, and being fast is what exposed it

The first outside tester reported it: under Sage with symfn installed, a
computation that runs long cannot be interrupted, where the same computation on
stock Sage breaks out of Ctrl-C normally. The report is exactly right and the
cause is structural rather than incidental.

All 101 `#[pyfunction]`s ran as plain Rust bodies holding the GIL for their
whole duration, and `check_signals` appeared nowhere in the tree. CPython
records SIGINT and runs the handler at its next bytecode boundary, so the
signal sat until the call returned. Measured before the fix, with SIGINT
delivered 1.0s into a 5.8s call:

| | signal delivered | `KeyboardInterrupt` raised |
|---|---|---|
| before | 1.0s | **4.83s** — when the call ended |
| after | 1.0s | **1.01s** |

**Stock Sage is interruptible for a reason that goes away when symfn is
installed**, which is why this arrived with the first user rather than earlier.
Sage's plethysm assembles its answer in the p basis in Python and finishes with
one coercion to s, so on stock Sage the work is spread over Python bytecode and
every loop iteration is an interrupt point. With symfn under it, that same
coercion is *one* call that can run for minutes. The speedup did not create the
defect; it merged thousands of small interruptible steps into one large
uninterruptible one, which is the general shape of the hazard and not a fact
about plethysm.

The mechanism is [interrupt.rs](../../src/interrupt.rs): the embedder installs
a `fn() -> bool`, the kernel polls it on loops whose trip count grows with the
input, and a poll that sees a cancellation unwinds with a payload no other
panic uses. `docs/policies/failure.md` carries the reasoning for the panic —
`Ring::mul` returns `Self`, so generic code has no channel to thread an
`Option` cancellation through, and the alternative was every long loop and
every caller of one changing signature.

**Cost, interleaved A/B, 4 rounds, min per (build, case), battery with low
power mode off** (`examples/bench_ops.rs`, `examples/bench_shapes.rs`):

| case | before | after | |
|---|---|---|---|
| `kostka_all_pairs_n20` | 2.0729s | 2.0903s | 1.01x |
| `character_beta_sweep_n28` | 1.3984s | 1.4016s | 1.00x |
| `convert_p_to_s` | 0.0127s | 0.0125s | 0.99x |
| `kostka_row_[8,7,6,5,4]` | 0.1433s | 0.1427s | 1.00x |
| LR `[14,12,10,8]²` | 0.1588s | 0.1566s | 0.99x |
| LR `[12,10,8,6]²` | 0.0598s | 0.0571s | 0.95x |

Worst case 1.01x and nothing outside noise, which is what the shape of `poll`
predicts: with no checker installed it is one relaxed load and a predictable
branch, and with one installed the checker itself runs once per 64 polls. The
sub-1.00x rows are noise, not a speedup. LR term counts were identical on both
sides, which is what says the chunked sequential fill changed only where the
poll sits.

**Two things were already right, and neither was written for this.** `memo`
runs `compute` outside its guard and clears lock poison rather than
propagating it, so an unwind leaves no half-built entry and does not disable a
table; `escalate` tests for `None`, so a panic passes through it instead of
being read as overflow and silently restarting the whole computation over
`BigInt`. Every store in the crate writes an already-finished value —
`bold_guarded` stores only once the overflow counter agrees, the monotone level
tables push finished levels — so a cancellation cannot leave a poisoned answer
behind. `tests/interrupt.rs` cancels at forty different depths and demands the
answers back, because that failure would be silent: a wrong result on the call
*after* the Ctrl-C.

What a cancellation does perturb is diagnostics. `PEAK_LIVE_STATES` and the
`measure` counters accumulate across an abandoned run, so a memory figure taken
right after a Ctrl-C includes work that never finished. Nothing depends on them
for an answer; the measurement discipline does.

### Docstrings that point where a Python reader can go

P11 used to say "pointers are backticked repo paths", written for docs.rs
rendering, and the Python surface followed it: `python/symfn/*.py` cited
`docs/policies/python.md` by rule number in seven places, `_families.py`
sent the reader to `src/llt.rs`, and 24 `#[pyfunction]` docstrings in
`src/python.rs` cited `docs/record/*` for a measurement or `scripts/*` for a
check. None of those is reachable from the wheel, `help()`, a stub tooltip,
or the two published sites, which are all that reader has.

P11 now says the opposite: house material (`docs/policies/*`,
`docs/record/*`, `scripts/*`) is never cited from a Python-facing docstring —
the docstring states the fact and stops — and depth on a family points at the
crate's rendered reference by module path (`symfn::llt` at
https://docs.rs/symfn) or at the docsite's Conventions page. Where a pointer
had been carrying a measured number as its evidence (`s_6[s_6]` in 0.25 s
against 28 s, Sage exceeding 90 s), the number went with it, which is the
rustdoc rule anyway; the shape-terms claim stayed. The private `///` docs in
`src/python.rs` and the `#` comments in the package are the maintainer's and
keep their tree paths.

`scripts/check_python_pointers.py` holds this over every string in
`python/symfn/*.py` and `symfn.pyi` and over the `///` block of every
`#[pyfunction]` and the `#[pymodule]` in `src/python.rs`, and
`scripts/preflight_python.sh` runs it as the "pointers" step. It reads
source and needs nothing built.

### Signatures on the site name their types

Sphinx documents the compiled module from `symfn.pyi`, and the stub's aliases
(`QtElement`, `JackElement`, …) were expanded in every rendered signature —
`delta_conjecture_side` came out as `list[list[tuple[tuple[int, ...],
list[tuple[int, int, int]]]]]`, and the stub spelled several shapes inline
rather than by alias in the first place. Three things changed. The stub now
postpones its annotations (`from __future__ import annotations`), which is
what lets autodoc consult `autodoc_type_aliases` at all; `docsite/conf.py`
builds that table from the stub's own assignments, so an alias added there is
documented by name with no second edit; and every signature in the stub is
written in a vocabulary of about twenty-five aliases, each with a `#:` line
saying what its tuples mean — `Partition`/`PartitionArg`, `Element`,
`TCoefficient`, `QtCoefficient`, `TElement`, `QtElement`, `AlphaAtoms`,
`JackCell` (the Rust boundary's own name), `MacdonaldElement`, `Edges`,
`BTable`, and the `*Arg` halves. The same structural type gets a different
name where it means something different: `AlphaAtoms`, `QtCoefficient` and
`IndexedElement` are all `list[tuple[int, int, int]]`.

Two Sphinx limits needed hooks in `conf.py`. An alias nested inside another
type — `list[QtElement]` — reaches the signature as
`list[TypeAliasForwardRef('...')]` while a bare one is resolved, so
`alias_names` recovers the name on `autodoc-process-signature`. And a type in
a signature is cross-referenced as `py:class`, whose lookup is restricted to
classes and exceptions, so an alias documented as module data is found by
name and rejected by kind; `link_alias` on `missing-reference` looks it up
again by name alone. `python_use_unqualified_type_names` shows `QtElement`
rather than `symfn.symfn.QtElement`.

One rename fell out: the coefficient aliases were going to be `TPoly` and
`QtPoly`, and `from .symfn import *` in the package's `__init__` would then
have shadowed the convenience class `symfn.QtPoly` for mypy. They are
`TCoefficient` and `QtCoefficient`; nothing in the stub's alias namespace may
share a name with anything the package exports.

## The marshalling suite: 198 checks, and two defects at i128::MIN on its first run (2026-08-21)

`scripts/check_python_marshalling.py` is the round-trip half the tail below
had open, written to the correction recorded there: a boundary test, not an
oracle. Three checks, all Sage-free, all against the built extension module,
run as the "marshalling" step of `scripts/preflight_python.sh` and so by CI's
`python` job:

- **shapes** — every exported callable (109 at writing) runs once on a small
  valid input and its return is validated against the stub alias `symfn.pyi`
  declares, `type() is` strict: tuples where tuples are promised, `int`
  coefficients with `bool` excluded, partitions weakly decreasing with no
  zeros, permutations with trailing fixed points dropped, denominators
  positive and in lowest terms, α-atoms primitive and increasing, no
  `(1 - q^0 t^0)` denominator factor. The table must name every export — the
  same completeness device as `check_python_boundary.py` — so a new entry
  point cannot ship with its encoding unvalidated.
- **widths** — coefficients at 1, the `i64` edges, both `i128` edges, and
  past them (`2^127`, `10^40`, `-2^200`) survive identity-shaped calls
  unchanged: the identity conversion, a product with `s_∅`, ω on `s_1`, a
  Schubert product with the identity permutation, `to_power` on `s_1`,
  `∇` on `s_1`. Identity-shaped so that any change in the value is the
  marshalling's, which is the module doc's "no ceiling" claim exercised on
  both sides of the escalation.
- **permissive inbound** — the `*Arg` halves: list against tuple, padded
  against normalized, and the strict outbound form handed straight back in,
  every spelling required to agree.

The first run reported four failures. Two were the suite's own cases
overstepping the contract — a `(q,t)` triple spelled as a list where the stub
promises a tuple, and an integrality assumption `s_2 = (p_11 + p_2)/2`
falsifies — and were fixed in the suite. Two were real, both at `i128::MIN` —
the one `i128` value with no negation in the width, and a value no test in
the tree had ever pushed through the boundary:

- **`to_power` returned `[]`** — a silently wrong value, R1's forbidden
  fourth outcome, on `to_power([([1], -2**127)], "s")`. The mechanism is the
  escalation seam: `GuardedRat::from_i128(i128::MIN)` reports-and-zeroes,
  which is correct *inside* a `guarded` window, but `build_rat` loads
  coefficients before the window opens, so the report was already in the
  counter when `guarded` read its baseline and the fast pass returned `Some`
  with the term zeroed — the wide pass, which had the right answer, never
  ran. `plethysm` and `internal_product` load through the same seam. The fix
  is one line: `BoundaryRat::from_coeff` for `GuardedRat` now declines
  `i128::MIN`, which is what routes `escalate` to the `BigRational` pass.
  The general lesson is the guarded-window protocol's edge: a report is only
  visible if it fires between the baseline read and the check, so a
  reporting *load* must instead decline.
- **the ∇ family panicked** — `PanicException` across the boundary, the
  outcome P8 exists to forbid, on the same coefficient: `qt_schur_in` stored
  `Rational::from_int(i128::MIN)` in the *panicking* ring, and the first sign
  flip inside `nabla` hit `Rational::neg`'s refusal. The extraction now
  refuses `i128::MIN` with the same typed `ValueError` as a value past the
  width — the `(q,t)` inbound window is `i128::MIN < v ≤ i128::MAX` — and
  the ∇-family docstrings name the width refusal in `# Raises`.

Everything else held on the first run: all 109 return shapes, every other
width in both directions, every inbound spelling. Both defects are pinned by
the suite's widths section, which runs `i128::MIN` through every path above.

## The doctest gate ran 80 of 271 examples; 16 unrun were wrong (2026-08-21)

`check_convenience_docs.py` collected examples with `doctest.DocTestFinder`
over the package and its private modules, on the premise that a layer that is
Python all the way down needs no custom extraction. The premise was false.
Every class in the layer says `__module__ = "symfn"` so `help()` reads well,
and that lie fails the finder's ownership tests three ways: the module scan
drops the class (its `__module__` is not the module being scanned), the
package scan accepts the class but drops its methods (a method's `__globals__`
are the private module's), and a class reachable only as the type of an
instance — the basis factories, the four family namespaces — is never reached
at all, because the finder does not recurse into instances. The gate printed
"80 examples pass" while 191 more sat in docstrings it never read. The
completeness half of the same script counted those items as documented, which
is what let the two halves disagree silently.

Sixteen of the unrun examples were wrong, and the cluster is exactly what
`docs/policies/validation.md` predicts for unexecuted convention pins —
plausible rival values, not typos: `qt_kostka([2], [1, 1])` documented as `q`
where the Garsia–Haiman orientation gives `t`; `hl.Qp([1, 1])` documented as
`t·s_11 + s_2`, the cocharge shape, where charge gives `s_11 + t·s_2`; all
four LLT entry points documented "in the Schur basis" with Schur-basis values,
where they return monomial (the P11 trap, again); `llt.Gtilde([2, 1], 2)`
documented with a nonzero value where a shape with a nonempty 2-core has no
ribbon tiling and the sum is empty; `e_1²` documented as `e_11 + e_2`, a false
identity; and `hl.Qp([2, 1]).at(t=1)` documented as `h_111`'s expansion rather
than `h_21`'s. Every corrected value was re-derived by hand or checked against
`docsite/conventions.md`, whose examples were written against a live
interpreter and were almost all correct.

The fix keeps the stock finder for modules and walks the supported classes
explicitly with a finder whose ownership test is waived — safe because every
member of such a class is defined beside it. The gate now runs 271 examples
over the same 94-item surface the completeness half counts. The same sweep
added `check_docsite_docs.py`: the docsite's narrative pages repeat docstring
values so a reader never leaves the page, and repetition is safe only while
both copies execute; each page runs as one interpreter session, the way it
reads. Its first run found two wrong outputs on `docsite/quickstart.md` —
`(q + t)` written in the literature's term order where the `QtPoly` repr
prints `(t + q)` — and nothing else.

The general lesson matches the marshalling suite's, one layer up: a gate that
counts a surface and a gate that executes it must walk the *same* enumeration,
or the counted-but-unexecuted gap rots in the dark precisely because the gate
is green.

### What is still open

- Cancellation latency inside a parallel Littlewood–Richardson row is one row,
  not one poll. A worker thread must not poll — the scoped join reads any
  worker panic as a bug, and the checker takes the GIL, which a worker cannot
  assume it may do — so a parallel section is cancelled at the boundary that
  dispatched it. On a shape whose single row runs for minutes that is the
  wait. Fixing it means a cancellation channel the workers can *read* rather
  than raise on, and it has not been needed yet.

- ~~The round-trip half of the Sage-free suite (Phase 5) is still not
  written: a value handed in comes back out intact, at the widths and shapes
  P1 promises.~~ **Written 2026-08-21** — `scripts/check_python_marshalling.py`;
  see "The marshalling suite" above, including the two `i128::MIN` defects
  its first run found.

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

## The oracle scripts were not refusing a live backend (2026-08-21)

Found while weighing whether to dispatch the new inverse expansions from the
Sage branch. The answer to that question is elsewhere; this is what looking
for it turned up.

**The defect.** `SAGE_DISABLE_SYMFN` is what keeps a Sage comparison from
being symfn quoting itself, and `docs/sage-backend.md` has said so since the
backend landed. Two files enforced it — `gen_sage_oracle.sage`, which refuses,
and `check_qt_kostka.py`, which refuses. Thirty others imported Sage and did
not. Meanwhile `sage.libs.symfn.is_available()` returns `True` in the
development environment, running the adapter out of `~/projects/sage` on
`combinat/symfn-backend`, where Hall–Littlewood `P` and `Q'`, Macdonald `J`,
Jack's forward direction, the character-basis product, plethysm and the
classical conversions all dispatch.

So `check_jack.py`, `check_macdonald.py`, `check_hl.py` and `check_hl_p.py`
were, when run without the variable, holding those families to themselves. The
scripts pass either way, which is the whole difficulty: there is no symptom on
the oracle side. On the benchmark side there is one — a run of ratios near
1.0× — and `bench_hl.py`, `bench_jack.py`, `bench_macdonald.py`,
`bench_qt_kostka.py`, `bench_vs_sage.py`, `compare_sage.py` and
`compare_symmetrica.py` had no guard either.

**Nothing wrong was published.** The committed fixture comes from
`gen_sage_oracle.sage`, which has always refused without the variable, and
every fixture-backed claim in the record files rests on that file rather than
on the check scripts. The check scripts are the *wider* half of each family's
oracle, run by hand; what they were failing to add was independence, not
correctness. Runs made with the variable set — which is the documented
invocation in every one of their docstrings — were honest all along.

**The fix, and why it is not thirty copies of ten lines.**
`scripts/sage_guard.py` holds the rule and one function, `require_own_sage`,
which exits unless Sage will answer out of its own code and passes straight
through on a stock Sage with no backend installed. It cannot *set* the
variable: `Feature.is_present` caches, so the variable has to be in the
environment before Sage starts. Thirty-two scripts call it, one line each.

The part that matters more is `scripts/check_sage_guards.py`, which fails when
a script under `scripts/` imports Sage and neither calls the guard nor names
itself in an allowlist with a reason. Five scripts are on that list, each
measuring what the backend does and needing it on for at least one arm. The
gate parses with `ast`, imports nothing, needs no Sage, and runs inside
`scripts/preflight.sh` — so the hole cannot reopen the way it opened, which
was one script at a time over several months, each one individually
reasonable.

**Checked both ways.** Stripping the guard from `check_jack.py` fails the gate
with the script named; with the backend live, `python scripts/check_jack.py`
now exits 1 saying which comparison would have been vacuous, and
`SAGE_DISABLE_SYMFN=1 python scripts/check_jack.py` still compares its 1212
values and reports 0 failures.

**What this says about the dispatch question.** Sending the inverse
expansions through the backend would widen the set of families a guardless
oracle script cannot see past — the Macdonald and Jack inverses are exactly
what `tests/fixtures/sage_oracle.txt` gained the same day. The gate had to
come first, and now has.

## Jack `P` supplies both directions to Sage (2026-08-21)

The Sage-side change is `e65f661f59e` on `mwhansen/sage` branch
`combinat/symfn-backend`; this is the symfn half of the record, because the
entry point it consumes is `monomial_to_jack_p` and the reason it pays is the
memoization in `src/jack.rs`.

Sage's `JackPolynomials_p._m_cache` filled `P → m` from `jack_p_table` and
handed it to `_invert_morphism`, which recovers `m → P` by a triangular solve
over ℚ(t). The comment there read: "here the inversion is already cheap and
Gram-Schmidt is nearly all of it."

**That was true when written and had already stopped being true.** It was the
same commit that made it false: once symfn took over the fill, the solve
became the larger half. On this branch, with the fill from symfn —

| n | `_m_cache` | of which `_invert_morphism` |
|---|---|---|
| 8 | 0.048s | 0.029s (59%) |
| 11 | 0.491s | 0.336s (72%) |

— and the share grows with the degree. A statement about which half dominates
is only true relative to the other half, so replacing one half invalidates it,
and this one was invalidated by the change that shipped alongside it. Nothing
re-read the comment for a month.

`jack_p_caches(n, ring)` in the adapter now returns both directions and
`_invert_morphism` is not called at all when symfn is present. The inverse
side calls `monomial_to_jack_p` once per shape, p(n) times; each hits the
memoized whole-degree table, so the p(n) calls cost **one** solve. That is
what `mac_p_inverse_cached`'s sibling was for, and this is its first outside
consumer.

Measured, AC power, one process per point:

| n | pure Sage | symfn fill only | both directions | vs fill only |
|---|---|---|---|---|
| 8 | 1.354s | 0.048s | 0.031s | 1.5× |
| 10 | 17.770s | 0.217s | 0.139s | 1.6× |
| 12 | 286.226s | 1.165s | 0.570s | 2.0× |
| 13 | — | 2.568s | 1.137s | 2.3× |
| 14 | — | 6.087s | 2.184s | 2.8× |

Against Sage's own route degree 12 is 502×. This is the fourth of the four
transitions to reach the both-directions shape, after Hall–Littlewood `P` and
`Q'` and Macdonald `J`.

**No normalization pass was needed**, which the Macdonald `J` change had
warned to expect: there the sign of a ℚ(q,t) fraction had to be matched cell
by cell and nine doctests failed on it first. ℚ(t) is univariate and
canonicalizes its fractions, so dividing in the ring lands on the
representative `_invert_morphism` produces. Checked rather than assumed —
both cache dictionaries are bit-identical to the old route's, 942 cells over
every degree through 8, dumped from separate processes and compared as
strings.

The equivalence check needed the guard work above to be meaningful in the
other direction too: the "old route" arm is `SAGE_DISABLE_SYMFN=1`, and
without that variable it would have been the new route compared to itself.

## Macdonald `m → P` and `m → Q` against the `J` route (2026-08-21)

Measured first, then adopted — the measurement and the decision are both
below, in that order. Macdonald `P` and `Q` did not dispatch directly: Sage
registers `P` as a diagonal coercion from `J` through `c2`, so `P(m(λ))` runs
`m → s → J → P` on `J`'s both-directions cache. The question was whether
`monomial_to_macdonald_p` beats that detour.

**It does, by a lot.** Every λ of the degree, from `m` into `P` or `Q`, both
arms ending on a list of basis elements, AC power, one process per point:

| n | `P` now | `P` direct | | `Q` now | `Q` direct | |
|---|---|---|---|---|---|---|
| 6 | 0.222s | 0.007s | 37× | 0.235s | 0.009s | 26× |
| 7 | 0.735s | 0.021s | 35× | 0.768s | 0.029s | 27× |
| 8 | 2.943s | 0.089s | 33× | 3.035s | 0.113s | 27× |
| 9 | 9.529s | 0.376s | 25× | 9.678s | 0.445s | 22× |

"Direct" is `monomial_to_macdonald_p`/`_q` per shape plus `_mac_cell`
marshalling into ℚ(q,t) plus building the elements, so the marshalling the
figure has to survive is in it.

**The win is not a cache-fill accounting artifact.** Building `J`'s `_s_cache`
for the degree first, untimed, barely moves the Sage arm — 9.403s to 8.782s at
degree 9 — so what costs is the per-element `m → s → J → P` arithmetic, not
the table it rides on.

### What it costs: the representative does not match

Values are exact — 234 cells through degree 6, 0 rows differing. Printed form
is another matter, and this is the trap `macdonald_s_to_j_table` documented,
met from the other side:

| `_mac_cell` | cells printing differently from Sage |
|---|---|
| `normalize=True` | 35 of 234 |
| `normalize=False` | 96 of 234 |

Neither setting matches, so this is not a flag to flip. Sage's representative
here is a byproduct of the composition — the `s → J` entry times an
unreduced `1/c2(λ)` — and it follows no rule that could be reproduced: at
degree 4 every denominator has positive leading coefficient, at degrees 3 and
5 some do not, and the split does not follow `|λ|`, `ℓ(μ)`, or the parity of
the factor count (which is what `normalize=False` would give).

So dispatching this changes what Sage *prints* for about 15% of Macdonald
`P` and `Q` coefficients, depending on whether symfn is installed.
That is the one thing the backend has held to throughout
(`_mac_cell`'s docstring: matching the representative "is what keeps
installing the backend from rewriting printed output"). The doctest blast
radius is small — 13 lines in `src/sage/combinat/sf/` print such a
coefficient as a fraction, 2 with a negative leading denominator — but the
divergence is not confined to doctests.

### Adopted the same day, at 4.5–8.6×

The judgment call was made: Sage takes the changed representative.
`836dda1a406` on the Sage branch adds `macdonald_pq_caches` and a mixin
carrying `_m_cache`, `_m_to_self` and `_self_to_m` for both bases. **The
integrated figure is much smaller than the 25–37× above**, and the difference
is worth stating because it is where the rest of the work now is:

| n | `P` before | after | | `Q` before | after | |
|---|---|---|---|---|---|---|
| 6 | 0.225s | 0.035s | 6.4× | 0.242s | 0.054s | 4.5× |
| 8 | 2.932s | 0.347s | 8.4× | 3.017s | 0.508s | 5.9× |
| 9 | 9.580s | 1.119s | 8.6× | 9.790s | 1.549s | 6.3× |

The table-building prototype stopped at the table; the integrated path has to
go through Sage's `_from_cache`, and that is now the whole cost.

**`_from_cache` substitutes `q` and `t` into every cell it reads**, even when
they are the ring's own generators and the substitution is the identity.
Against the same map applied without it: 0.064s → 0.001s at degree 7, 0.197s →
0.004s at 8, 0.519s → 0.006s at 9 — **49× to 86×**. That is where the missing
factor went. It is shared by every parametric basis (Jack, Hall–Littlewood,
Macdonald, `orthotriang`), so a fast path there is a larger change than this
one and is not made here. It is the biggest single number left on the Sage
side.

**Registered as a conversion, not a coercion.** `P` already reaches everything
through `J`, and a second *coercion* changes which composite Sage's discovery
picks for unrelated pairs: `P → H̃` started routing through `m`, which
`sage/structure/parent.pyx` has a doctest printing. Nothing got slower either
way — `P → H̃` was 5.369s against 5.429s at degree 7 — but an explicit
`P(m[2])` is a conversion, and a conversion does not join that graph. Worth
remembering for the next basis: `register_coercion` has effects beyond the
pair it names.

**One doctest printed the old representative.** Rather than pin the new one,
which would make that file's output depend on whether symfn is installed, it
now asserts the coefficient — so it checks the mathematics and passes in both
arms. `src/sage/combinat/sf/`, `sage/structure/parent.pyx` and
`non_symmetric_macdonald_polynomials.py` are green with the backend and with
`SAGE_DISABLE_SYMFN=1`.

## The families name their own basis (2026-08-21)

`jack.P([2])` returned `2/(alpha + 1)*m[1,1] + m[2]`. It now returns
`JackP[2]`, and `.to("m")` returns the old value. Same for the other eight
constructors: `macdonald.P`, `.Q`, `.J`, `.Htilde`, `jack.Q`, `.J`, `hl.Qp`
and `.P`.

**What was wrong before.** The inverse expansions of 2026-08-21 introduced
nine parametric bases — `McdP`, `McdQ`, `McdJ`, `McdHt`, `JackP`, `JackQ`,
`JackJ`, `HLP`, `HLQp` — as the tags their results carry. That left the
surface asymmetric: `jack.to_P(m([2]))` printed `JackP` terms, while
`jack.P([2])` printed an expansion, so the basis a family is *about* was
reachable only by going out into the monomial basis and back. A shape was the
one thing you could not ask for by name.

**Why it needed kernel work.** `Param` is an inert value object — no
arithmetic, and none on `AlphaFrac`, `QtFrac` or `QtRatio` either. So
returning `JackP[2]` alone would have been contentless, and `.to("m")` on a
multi-term element means adding and multiplying rational functions, which
`docs/policies/python.md` puts in the contract layer, not this one. Every
forward entry point took a single partition. Nine new ones take an element:

    jack_p_to_monomial          macdonald_p_to_monomial
    jack_q_to_monomial          macdonald_q_to_monomial
    jack_j_to_monomial          macdonald_j_to_monomial
    hall_littlewood_p_to_schur  hall_littlewood_qp_to_schur
    macdonald_ht_to_schur

Eight of them take the encoding their inverse sibling returns, so the two
directions share a row format and neither converts anything. They cost one
forward polynomial per shape *present*, not per shape of the degree — which is
what makes them right for an element with few terms and the `_table` functions
right for a whole degree.

**What the change bought in evidence.** `m_μ → basis → m_μ` is now a closed
composite at every shape, and it is the direction the existing round trips did
not cover: `every_monomial_comes_back_as_itself` in `jack.rs` and
`macdonald.rs`, `every_schur_function_comes_back_as_itself` in `hl.rs`. A
table inverted correctly in one direction only would pass the old tests and
fail these.

**`at` expands rather than refusing.** It used to raise on a parametric basis
because no `Sym` can carry one. It now goes through `to` first, so
`macdonald.P([2]).at(q=5, t=5)` still gives `m[1,1] + m[2]` — every `at`
example in the tree survived unchanged, which is the check that the two
directions compose.

**`McdHt` is the one basis whose expansion is partial.** Its coefficients
cross the boundary as a numerator over an *expanded* denominator, while the
crate divides by a factored multiset of `q^a − t^b` atoms — so a denominator
that is not 1 cannot be handed back, and `_ht_rows` refuses rather than
dropping it. `Htilde(mu).to("s")` works; the general output of `to_Htilde`
does not. Closing that means giving `HtElement` a factored denominator, which
is a new documented encoding (`docs/policies/python.md`, the last row of the
home table) and a change to a contract type, so it was not made here.

**LLT has no basis of its own** and its constructors still return the monomial
basis directly, because there is no expansion *into* an LLT basis to make the
tag mean anything. The quickstart says so rather than leaving the reader to
notice.

**The Sage adapter is untouched.** It imports contract entry points only —
`symfn.jack_p`, `symfn.macdonald_p`, and the rest — and never the `jack`,
`macdonald` or `hl` namespaces, so none of this reaches
`sage/libs/symfn/backend.py`.

### The bases were places you could name but not compute in (2026-08-21)

`jack.P([2])` returning `JackP[2]` left no way to write `q·H̃_{21}`: `Param`
had no arithmetic and there were no parameter values, so a scalar meant
building `QtRatio([(1, 0, 1)], [(0, 0, 1)])` by hand. `symfn.q`, `symfn.t`,
`symfn.t_hl` and `symfn.alpha` are those values now, and `Param` has `+`, `-`,
unary `-` and scalar `*`.

**Where the arithmetic had to live was decided by one fact**: the Python
coefficient classes compare **structurally**, not by cross-multiplying —
`QtFrac.__eq__` is `num == num and den == den`. So an unreduced sum is the
right number in a representation nothing else produces, and `==` against a
value built another way would be false. `AlphaFrac` does not even normalize its
integer scale on construction: `AlphaFrac([0, 2], (), 2)` prints `2*alpha/2`.

That ruled out adding in Python. Four contract entry points do it instead —
`macdonald_element_add`, `macdonald_element_scale`, `jack_element_add`,
`jack_element_scale` — each reducing through the crate's own `Frac::reduce` and
`AFrac::reduce`. They are basis-blind, because addition in a basis does not
depend on which basis it is; the tag stays in the convenience layer. The check
that this was the right call is in `check_convenience.py`:
`macdonald.to_P(m([2])) + macdonald.to_P(m([1,1]))` equals
`macdonald.to_P(m([2]) + m([1,1]))`, which needs a common denominator on one
side and not the other.

`Poly` and `QtPoly` add and multiply in Python, because a polynomial sum is
already canonical, and they needed the ring operations anyway so that
`1 - q*t` can be written as a scalar.

**`McdHt` is again the partial one.** Multiplying is fine — the denominator is
untouched — but two coefficients at one shape with different denominators
cannot be put over a common one, since the crate divides by factored
`q^a - t^b` atoms and the encoding hands them over multiplied out. It refuses
rather than answering over `b*d`, which would be the right value in a
representation nothing else produces. That is the third thing this encoding has
blocked, after `Param.to` and scaling by a fraction.

`adding_and_scaling_commute_with_expanding` in `src/macdonald.rs` and
`src/jack.rs` ties the four to a route that never touches them: `Monomial` adds
and scales through the ordinary `Ring` operations, so agreement is not the two
sharing an implementation.

### `H̃` stopped being the exception (2026-08-21)

Three separate pieces of work hit the same wall — `Param.to`, `Param.__add__`,
`Param.__mul__` by a fraction — and each time the answer was that `HtElement`
hands its denominator over multiplied out while the crate divides by factored
atoms. The encoding now carries `(kind, a, b, multiplicity)` atoms and all
three work; `docs/record/qt-kostka.md` has the account and the reason the kind
tag cannot be dropped.

What that closed, in the convenience layer: `macdonald.Htilde(mu)` and
`macdonald.to_Htilde(f)` both expand with `to("s")`, `McdHt` elements add and
scale like the other eight bases, and `macdonald.to_Htilde(s(la)).to("s")`
returns `s_λ` exactly rather than only agreeing at a point. That round trip is
checked at every shape in `check_convenience.py`.

**A `Sym` scaled by a parameter is now a `Param` in the same basis.**
`q * m([2])` used to raise — a `Sym` carries `int` and `Fraction` coefficients
and nothing else — which left no way to hand a scaled classical element to an
inverse expansion. It lifts instead, keeping the basis, so
`macdonald.to_P(q * m([2]))` is the ordinary way to write that. `Sym.__pow__`
calls a private `_times` rather than `*`, because a parameter has no place in
the middle of a composed product and the narrower return is what keeps that
loop typed.

## A parameter is a base ring, not a kind of object (2026-08-24)

`Param` held two unrelated things — an element in one of the nine parametric
bases, and an element in a classical basis whose coefficients happen to carry
`q`, `t` or α — and the split was read as a mathematical distinction. It is
not. A user of this layer holds three independent facts about a value: the
base ring its coefficients live in, the basis it is written in, and the element
itself. `Param` is neither of the first two: it means "the coefficients are not
`int` or `Fraction`", which is a fact about representation.

**The check that settled it was Sage's own model**, run with
`SAGE_DISABLE_SYMFN=1` so it was not this library answering. Sage has *more*
element classes than symfn — one per basis — and nobody notices, because they
are interchangeable in every observable way: `P[2]*P[1]` multiplies, `s(P[2])`
converts, `P[2].omega()` is defined, and `(q*m[2]).degree()` is 2. So a
parametric basis is not narrower than a classical one, and the earlier claim
here that tier A had "one legal expansion, no product" was wrong about the
mathematics rather than about symfn's coverage.

**What was decided**, and written into
[docs/policies/python.md](../policies/python.md) under P10: the nine tags are
bases on the same footing as the six codes, the parameters are the base ring,
`q * m([2])` returns whatever `m([2])` returns, and `Sym`/`Param` being two
classes is an artifact of where the coefficient arithmetic lives rather than a
distinction a caller may rely on. Two refusals stay — different bases do not
add, different base rings do not combine.

**The release gate turned out to be smaller than it looked.** The plan had held
that the return type of `q * m([2])` could not change after 0.1.0 without
breaking a caller, which made the whole merge a release blocker. If the merged
class is named `Sym` and `Param` is kept as an alias of it, `isinstance(x,
Param)` stays true of everything it is true of today; the only observable
change is that it becomes true of elements it used to be false of. So what
precedes 0.1.0 is the commitment, and the merge follows the work.

**First items landed.** `Param.degree` and `Param.is_homogeneous`, same bodies
as `Sym`'s — the degree is a fact about the partitions alone, so
`macdonald.P([2]).degree()` is 2 with no expansion and no coefficient
arithmetic. That is the entire gap for anyone who only ever scales a classical
element.

**One defect found while checking, not yet fixed.** `alpha * m([2]) + q *
m([2])` raises `TypeError: unsupported operand type(s) for +: 'Poly' and
'QtPoly'` — a base ring mismatch, in a single basis, leaking as a
coefficient-class accident. `Param.__add__` compares `_params`, but only
reaches that check when the bases differ, so this case escapes it. It needs
the treatment `BasisError` gets: one error type, naming ℚ(α) and ℚ(q,t).

The rest — six-way `to`, ω, products in all fifteen bases, and the merge — is
[docs/plans/element-model.md](../plans/element-model.md).

## Six-way `to`, for polynomial coefficients (2026-08-24)

`Param.to` reached one basis: the classical pivot its family expands in. It now
reaches five — `s`, `h`, `e`, `m`, `f` — for the coefficient classes that are
polynomials, which is Hall-Littlewood, LLT, `H̃` where the denominators cancel,
and any classical element scaled by a parameter.

    >>> hl.Qp([1, 1]).to("h")
    h[1,1] + (-1 + t)*h[2]
    >>> (q * m([2])).to("s")
    -q*s[1,1] + q*s[2]

**What was needed in the crate.** `convert` picks its route — direct rule or
Schur hub — from the *types* of its two ends, and a caller holding a basis code
has no types to offer. `convert_named` in [convert.rs](../../src/convert.rs) is
the same routing with the pair resolved at runtime from `SymFn::SYMBOL`, twelve
arms rather than thirty-six because the destination half is factored out. It is
generic over `C: Ring`, so nothing about it is specific to `q` and `t`.

**`p` is not one of the five.** Conversions into the power-sum basis divide by
z_μ and so want a `QAlgebra`; `QtPoly<i128>`, which carries every `t`- and
`(q,t)`-polynomial coefficient at this boundary, is a `Ring` and nothing more.
The first draft of `convert_named` was bound on `QAlgebra` and would not
compile against `QtPoly<i64>`, which is how this was found. The integer path
states the same restriction and sends `p` through `to_power`.

**One boundary entry point, `convert_qt_terms`**, taking `nabla`'s
`[(lambda, [(q_exp, t_exp, coeff), ...])]` rows with `src` and `dst` names, and
escalating over `BigInt` on the same pattern as the Hall-Littlewood pair. The
Python side packs a one-variable `Poly` into the `t` slot and restores the
variable name on the way back, which is legitimate because the entry point
never asks what the exponents count — P1 in
[python.md](../policies/python.md) is exactly that.

**The count in the plan was wrong, and is corrected there: four entry points,
not two.** A `Param` carries five coefficient classes over four Rust rings —
`QtPoly`, `Frac`, `AFrac`, `Ratio` — and each has its own boundary encoding.
`macdonald.P([2]).to("s")` and `jack.P([2]).to("s")` still refuse, and the
message now says the mathematics is a basis change like any other and the
converter is what is missing.

**The evidence shares no route with the thing it checks.** `Param.to` carries
polynomial coefficients through `convert_qt_terms`; `Sym.to` clears
denominators and calls the integer conversions. So
`f.to(b).at(t=v) == f.at(t=v).to(b)` compares two implementations rather than
one with itself, and `check_parametric_conversions` in
`scripts/check_convenience.py` sweeps it over both Hall-Littlewood
normalizations, every shape to degree 4, all five destinations, and four values
of `t` including 0 and 1 where the family degenerates. The suite went from 4424
to 5140 checks. In the crate, `convert_named_scales_with_the_coefficient_ring`
converts an `i64` element and the same element scaled by `q²t` over every
ordered pair and requires the answers to differ by exactly that factor — `i64`
addition against `QtPoly` addition, so a transposed arm in either dispatch
table fails at the pair that names it.

**The quickstart documented the old refusal** and its doctest is what caught
the behavior change. It now shows the six-way conversion and keeps a refusal
example, pointed at the rational-function families where one still applies.

## ω and the antipode, on the same three families (2026-08-24)

`Param.omega` and `Param.antipode` exist, over the coefficient classes six-way
`to` reaches, and they return in the basis they were handed — including a
parametric one.

    >>> hl.Qp([2]).omega()
    HLQp[1,1] - t*HLQp[2]
    >>> (q * m([2, 1])).antipode()
    q*m[2,1] + 2*q*m[3]

**Neither needed a basis argument or an escalation.** Both are defined on the
Schur basis, so the boundary pair `omega_qt_terms` and `antipode_qt_terms` acts
there and the two changes of basis around it are `convert_qt_terms`. Neither
adds in the coefficient ring: conjugation is a bijection on the partitions of a
degree, so no two terms can meet, and the antipode only copies a sign. The one
input where copying a sign is not total is `i128::MIN`, whose negation leaves
the narrow arm — `Coeff::negated` widens to `BigInt` there, and the marshalling
suite's width round trips cover it.

**Three legs rather than one.** `Sym` applies ω by converting to Schur and
back. A `Param` in a family's own basis expands into its pivot first and
travels back through the inverse expansion, so `hl.Qp([2]).omega()` is an
`HLQp` element rather than a Schur one. That last leg is code the classical
route never runs, which is why `check_parametric_hopf` in
`scripts/check_convenience.py` asserts the basis as well as the value. It works
for `HLP`, `HLQp` and `McdHt`; the six Macdonald and Jack tags wait on the same
three converters six-way `to` waits on.

**`_carry` is the refactor that made it cheap.** Packing an element's
coefficients into exponent rows, clearing denominators, calling once, and
rebuilding in the class it went in as is now one function in
`python/symfn/_families.py`; `_convert` and `_hopf` differ only in the call
they pass it. The next three converters inherit it.

**What the evidence is.** `f.omega().at(t=v) == f.at(t=v).omega()` over both
Hall-Littlewood normalizations and `H̃`, every shape to degree 4, three values
each, plus ω being its own inverse in the parametric basis. The two sides share
no route — `Sym.omega` runs the integer entry points over integer coefficients.
The suite went from 5140 to 5461 checks.

The antipode's doctest uses a shape of odd degree on purpose: at even degree it
equals ω, and a value the two agree on pins neither.

## The other three converters, so every family reaches every basis (2026-08-24)

Six-way `to` shipped for polynomial coefficients only, which left Macdonald and
Jack reaching one basis each. The three rational-function converters close it,
and `Param.to` is now total over the nine tags and the five classical
destinations.

    >>> jack.P([2]).to("s")
    (1 - alpha)/(alpha + 1)*s[1,1] + s[2]
    >>> macdonald.P([2]).to("s")
    (-t + q)/(1 - q*t)*s[1,1] + s[2]

**Both values were predicted before the code existed, by different sources.**
The Jack one is what `docs/plans/element-model.md` wrote down from
`jack_p(&[2]).to_schur()` in the crate; the Macdonald one is Sage's `s(P[2])`
from the run recorded above, term for term. Neither is this library checking
itself.

**Three entry points, and they cost almost nothing.**
`convert_macdonald_terms`, `convert_jack_terms` and `convert_ht_terms` live in
[python.rs](../../src/python.rs).
Each parses with the row builder its family's inverse expansion already had,
hands the term map to `convert_named`, and emits with that family's writer. The
routing is identical across all four converters because `convert_named` is
generic in `C: Ring`; what differs is only the encoding on the wire. Two
helpers were factored out while adding them: `convert_pair`, which parses the
basis names and rejects the power-sum destination once, and `routed_ring`,
which is the call plus the R2 panic for the state `convert_pair` has already
excluded.

**`p` stays closed, and the earlier note here about it was too optimistic.**
That conversion divides by z_μ and needs a ring containing ℚ. Of the four
rings, `AFrac<C>` is a `QAlgebra` for any `C: Ring`, but `Frac<C>` is one only
when `C` is, and `QtPoly<i128>` is not. So opening `p` would reach Jack and
nothing else, which is a worse surface than a uniform refusal.

**The `H̃` converter is the one the convenience layer barely reaches.**
`_expand` collapses an `H̃` expansion to `QtPoly` whenever the atoms cancel,
which is the usual case because `K̃_{λμ}` is a polynomial — a sweep over
`macdonald.to_Htilde(c * s(λ))` for four scalars and six shapes produced no
element with a surviving denominator. So the `QtRatio` leg is exercised by
constructing one directly, `q/(q − t)·s_2`, which converts to
`q/(q − t)·(m_2 + m_11)`. That is in `check_convenience.py` rather than left to
a case that may not arise.

**What the evidence is.** `f.to(b).at(...) == f.at(...).to(b)` over all nine
tags, every shape to degree 4, all five destinations, and two or four parameter
values each; the two sides share no route, since `Sym.to` clears denominators
and calls the integer conversions. Plus the classical limits reached through
the *new* route rather than through the pivot — `P_λ(x; 1) = s_λ` for Jack and
`P_λ(x; q, q) = s_λ` for Macdonald, both of which fail under the `α → 1/α` and
`q ↔ t` twists. The suite went from 5461 to 6861 checks.

**ω and the antipode did not come along.** Both act in the Schur basis and the
entry point that does so reads the polynomial encoding, so
`macdonald.P([2]).omega()` still refuses where `to` no longer does. Recorded as
an open item with two candidate routes, one of which — ω sends `h_μ` to `e_μ`,
so it is a relabeling with `to` on either side — needs no new boundary but puts
a mathematical identity in the convenience layer.

The quickstart documented the refusal that just went away, and its doctest
caught it for the second time in two changes. It now shows both values above.

## The Hall-Littlewood product, and the normalization it pins (2026-08-24)

`Param.__mul__` refused two elements with "a parametric basis has structure
constants this does not compute". It computes them now, for every coefficient
class the polynomial encoding carries: Hall-Littlewood, LLT, and any classical
element scaled by a parameter.

    >>> hl.P([1]) * hl.P([1])
    (1 + t)*HLP[1,1] + HLP[2]

**One entry point, and it is the ordinary Schur product.** `schur_multiply_qt`
is `schur_multiply` over the `(q,t)`-polynomial encoding, reaching the same
Littlewood-Richardson backend, because the structure constants are integers and
carry no parameter — `Schur<C>::mul` is generic in `C: Ring`, so nothing in the
crate changed. The route around it is the one the plan predicted: expand to the
pivot, convert to Schur, multiply, and return the same way. `Param.__pow__` is
repeated squaring over it.

**The normalization, and the numbers reproduced.** Sweeping every `P_μ · P_ν`
with `|μ| = |ν| ≤ 5` gives **1871 coefficients, 331 of them negative** —
exactly what `docs/plans/element-model.md` recorded from a scratch experiment
that was not kept. The implementation and that experiment reached the same
normalization independently, which is the strongest thing available here short
of Sage.

`P[2,1]² → P[3,1,1,1]` is `1 + t − t³ − t⁴`, and that is the doctest on
`Param.__mul__`. `P[1]² = P[2] + (1 + t)·P[1,1]` is deliberately not the pin:
every convention in circulation gives it, so it distinguishes nothing. Two
readings of the negative value confirm the convention, and both were checked
against something else in this tree rather than asserted: its constant term is
`c^{3111}_{21,21} = 1` against `symfn.lr_coefficient`, since `P_λ(x; 0) = s_λ`,
and its value at `t = 1` is 0 against `m([2,1])**2` having no `m[3,1,1,1]`
term, since `P_λ(x; 1) = m_λ`. So these constants sit in ℤ[t], while the
classical Hall polynomials counting subgroups of abelian p-groups sit in ℕ[t];
the two differ by a twist, and this is the one that is being shipped.

**What the checks are, and what they are not.** `check_hall_littlewood_products`
sweeps `t = 0` to the Littlewood-Richardson product and `t = 1` to the monomial
product over every `P_μ · P_ν` with `|μ| = |ν| ≤ 4`, generic `t` against
multiplying the two specialized expansions, `Q'` at `t = 0`, and repeated
squaring against repeated multiplication. It also asserts that the sweep
*reached* a negative coefficient, since a pin that never sees one pins nothing.
The suite went from 6861 to 7066 checks.

These are all this library checking itself. The offline fixture sweep against
Sage that `docs/policies/validation.md` asks for on a family Sage covers is
still owed, and is the open item in the plan.

**No new failure-policy row was needed.** `schur_multiply_qt` escalates to
`BigInt` on the same pattern as everything else at this boundary, and the
Hall-Littlewood coefficients are in ℤ[t] with no denominators. The overflow the
plan records is in the *Jack* product, whose coefficients are in ℚ(α) and grow
much faster; that is still open and still wants the boundary row.

Two refusal messages went away with this: `Param.to`'s "already classical" and
the product's citation of structure constants. Both had been describing this
library's coverage in the language of mathematics, which is what P10 in
[python.md](../policies/python.md) now forbids.

## Sage confirms the Hall-Littlewood normalization (2026-08-24)

The products shipped with specialization pins and a convention doctest, all of
which were this library checking itself. Sage now backs them.

    P[2,1]^2 coefficient of [3,1,1,1]: -t^4 - t^3 + t + 1

That is Sage's own `hall_littlewood().P()`, run with `SAGE_DISABLE_SYMFN=1`, and
it is symfn's `1 + t - t^3 - t^4` written in the other order. `P[1]^2` and
`Qp[1]^2` agree too.

**78 products are committed as fixtures** — `P_μ · P_ν` and `Q'_μ · Q'_ν` for
every pair with `|μ| = |ν| ≤ 4` — as `hlpmul` and `hlqpmul` records in
`tests/fixtures/sage_oracle.txt`. `hall_littlewood_products_match_sage` in
`tests/sage_oracle.rs` reads them through the crate: expand both operands with
`hall_littlewood_p`, multiply with the Littlewood-Richardson backend, and
back-substitute with `schur_to_hall_littlewood_p`.
`check_hall_littlewood_products_against_sage` in `scripts/check_convenience.py`
reads the same records through the Python side, which adds the packing, the
denominator clearing and the rebuild that the crate route never touches. Both
pass over all 78.

Sage reaches these by coercing both operands into the Schur basis and inverting
the transition matrix. symfn expands through its own forward polynomials and
back-substitutes. So the two share the definition of `P` and `Q'` and nothing
about how the product is obtained.

**The regeneration changed nothing else.** Diffing the new fixture against the
committed one with the two new tags filtered out is empty, so the 4179 existing
lines are byte-identical and the 78 new ones are the whole change.

**The fixture test discriminates the two normalizations, and that was checked
rather than assumed.** Swapping `P` for `Q'` in the test fails at the very
first pair, `P_1 · P_1`, with `(1+t)·P_11 + P_2` against
`P_11 + (1−t)·P_2`. Both sweeps also assert they saw a negative coefficient,
because a pin that only meets `P_1² = P_2 + (1 + t)·P_11` — which every
convention in circulation gives — pins nothing.

What is still owed is the same evidence for the other families' products, which
do not exist yet: `α = 1` to Schur for Jack, and the Macdonald pairs.

## Every family multiplies, and a second exception (2026-08-24)

Hall-Littlewood multiplied; the six Macdonald and Jack tags did not, because
the multiply read the polynomial encoding. Three more entry points close it —
`schur_multiply_macdonald`, `schur_multiply_jack`, `schur_multiply_ht` — and
all nine tags now have a product.

    >>> jack.P([1]) * jack.P([1])
    2*alpha/(alpha + 1)*JackP[1,1] + JackP[2]
    >>> macdonald.P([1]) * macdonald.P([1])
    (1 + t - q - q*t)/(1 - q*t)*McdP[1,1] + McdP[2]

**Sage agrees on both, and on the first shape where the answer is not obvious.**
`McdP[2]·P[1]` is `(1 − qt² − q² + q³t²)/((1 − qt)(1 − q²t))` here and
`−(q³t² − qt² − q² + 1)/(−q³t² + q²t + qt − 1)` in Sage, which is the same
after clearing the signs and expanding the factored denominator; `JackP[2]·P[1]`
is `(4α + 2α²)/((α + 1)(2α + 1))` here and `(a² + 2a)/(a² + 3/2·a + 1/2)` in
Sage, the same after scaling by 2.

**`_back_to` is what had to grow, and six-way `to` is what paid for it.** It
knew only the tags whose pivot is Schur. It now converts into whichever pivot
`EXPANDS_IN` names and calls that family's inverse expansion, so all nine tags
are re-enterable — which also means ω and the antipode reach the Macdonald and
Jack tags the moment their Schur-basis half is written.

**28 products are committed as fixtures**, `macpmul` and `jackpmul` for every
pair with `|μ| = |ν| ≤ 3`, read by `parametric_products_match_sage` in
`tests/sage_oracle.rs`. Compared by evaluation at three generic points, not by
representation: symfn keeps denominators factored and Sage expands them, so
`(1+q)(1−t)/(1−qt)` has two correct normal forms. The route under test —
expand, convert to Schur, multiply, come back, re-enter `P` — exercises the
forward expansion and the inverse at once, so a wrong inverse cannot be
absorbed by a matching wrong forward one, which the `macp`/`jackp` records
alone cannot rule out. The regeneration left the other 4257 lines byte-
identical.

`check_parametric_products_degenerate` adds the half that needs no fixture:
`P_λ(x; 1) = s_λ` for Jack and `P_λ(x; q, q) = s_λ` for Macdonald, so a product
of two of them specializes to the Schur product computed over integers. Both
fail under the `α → 1/α` and `q ↔ t` twists. The suite went from 7145 to 7204.

**The Jack overflow this tree recorded is real, and it escalates.**
`docs/plans/element-model.md` had `P[8]²` at degree 16 panicking with "attempt
to multiply with overflow". At the Python boundary it does not: the boundary
row of the mechanism table applies, `schur_multiply_jack` escalates over
`BigInt`, and the answer arrives. Measured 2026-08-24 (debug build, AC power,
Apple M4, caches not cleared between cases): `jack.P([4])²` 0.02 s,
`jack.P([6])²` 0.82 s, `jack.P([8])²` 94 s. Slow and correct is what the policy
asks for, so no new row was needed — the plan item is closed rather than
actioned.

**`symfn.BaseRingError` is new, and it fixes a defect recorded above.**
`alpha * m([2]) + q * m([2])` used to raise `TypeError: unsupported operand
type(s) for +: 'Poly' and 'QtPoly'` — the coefficient classes' own failure
leaking through what is a question about the elements, for two operands in the
*same* basis. It now says "cannot combine an element over Q(alpha) with one
over Q(q, t)". It sits beside `BasisError` and subclasses `TypeError` for the
same reason, and the two are separate because only the basis mismatch is fixed
by `.to()`.

⚠️ **`Param`'s cross-basis refusal changed exception type**, from `ValueError`
to `BasisError`, in the same change. That is what `Sym` has always raised for
the same question, so the two classes now agree — one of the interface
differences the merge in `docs/plans/element-model.md` was waiting on. A caller
catching `ValueError` around `Param` arithmetic is affected; `BasisError`
subclasses `TypeError`, not `ValueError`.

## ω and the antipode over the rational-function rings, 2026-08-24

Six entry points — `omega_macdonald_terms`, `antipode_macdonald_terms`,
`omega_jack_terms`, `antipode_jack_terms`, `omega_ht_terms`,
`antipode_ht_terms` — finish the pair for the four coefficient rings, so all
nine parametric tags answer both. They share one generic `hopf_of`, which
conjugates the index and signs when the operation is the antipode and the
degree is odd. Conjugation is a bijection on the partitions of a degree, so no
two terms meet and the coefficient ring is never added in; the only arithmetic
is `Ring::neg`, and the two rings that can decline it escalate on the boundary
row rather than refusing.

**The route not taken was the free one.** ω sends `h_μ` to `e_μ`, so it could
have been a relabeling of the basis tag with a conversion on either side and no
new boundary at all. It was rejected on P4: that is a mathematical identity,
and putting it in the convenience layer is what P4 exists to prevent. The cost
of the route taken is six names on the contract surface.

**Sage agrees, checked in `ℚ(q,t)` rather than by string.**
`macdonald.P([2]).omega()` is
`(1 − t² − q² + q²t²)/(1−qt)²·McdP[1,1] + (q − t)/(1−qt)·McdP[2]`; Sage writes
the same value with the signs the other way. `macdonald.J([2,1]).omega()` has
three coefficients that Sage writes over expanded denominators, and asking Sage
whether each pair is equal as a fraction gives `[True, True, True]`.
`jack.P([2,1]).antipode()` agrees with Sage after scaling by 2, which is where
Sage puts a `1/2` in the denominator and this library does not.

`jack.P([2]).omega()` is `4α/(α+1)²·JackP[1,1] + (1−α)/(α+1)·JackP[2]`, and it
is the doctest on `Param.omega` because it pins which ω is meant: the plain
involution carries α, and the α-deformed one sends `P_λ^{(α)}` to
`Q_{λ'}^{(1/α)}` and would invert the parameter.

**This found a defect in the products committed earlier the same day.**
`_back_to` converted into `EXPANDS_IN[tag]` before calling the family's inverse
expansion. For `McdJ` those are two different bases: `J` expands in the
monomial basis, but `schur_to_macdonald_j` is triangular the other way and
reads the Schur basis. So `macdonald.J([1])**2` raised "to_J needs a
Schur-basis element, not m". `_INVERSE` now carries the basis each inverse
reads beside the function it calls, and `EXPANDS_IN` is no longer consulted on
the return leg. `macdonald.J([1])**2` is
`(1−q)/(1−qt)·McdJ[1,1] + (1−t)/(1−qt)·McdJ[2]`, and `macdonald.J([2])·J([1])`
is `(1−q²)/(1−q²t)·McdJ[2,1] + (1−t)/(1−q²t)·McdJ[3]`; both are Sage's values.

`_demote` is the one piece of encoding that leg needed. It rewrites
coefficients that are fractions with an empty factored denominator over the
polynomial class their numerators already are, and returns the element
untouched otherwise. `J` is the integral form, so a Schur-basis element on its
way back into it has polynomial coefficients — but the route there passes
through the monomial basis over `ℚ(q,t)` and comes out in that ring's class,
which `to_J` does not read. `_ht_element` already narrowed the same way when no
atoms survived, so this is that rule applied to the second ring rather than a
new one.

**Evidence is the specialization, which shares no entry point with the route.**
`check_parametric_hopf` now covers `McdP`, `McdQ`, `McdJ`, `JackP`, `JackQ` and
`JackJ`: ω is an involution on each, both operations come back in the basis
they were handed, and setting the parameter first and acting over ℚ gives the
same answer as acting first and setting it after. Two more pin the
degenerations against the integer path — `jack.P(λ).omega()` at α = 1 and
`macdonald.P(λ).omega()` at q = t both equal `s(λ).omega()`, which runs the
integer entry points end to end. The suite went from 7204 to 7562 checks.

## Three interface differences closed, 2026-08-24

Found while checking what the merge in `docs/plans/element-model.md` still
waits on.

**`Param` refused a scalar in `+` and `-`, and had no `__radd__`.** So
`0 + q*m([2])` raised, and `sum` over a list of parametric elements raised on
its first term, because `sum` starts from `0`. A scalar now adds as the
constant it names times the unit: `_constant` is `_scale` applied to
`_unit_like`, so a fraction ring reduces the product at the boundary rather
than in Python, which is the same reason `_add` sends the fraction kinds
through the contract layer. `2 + jack.P([1])` is `2 + JackP[1]`, the shape
`Sym` has always given. The scalar set is the one `*` already took, so a scalar
that can multiply an element can add to it.

**`_unit_like` covered only the two polynomial classes**, so
`jack.P([1])**0` and `macdonald.Htilde([1])**0` raised "the unit is not written
for AlphaFrac". It covers all five now. An element with no terms records its
parameters but no coefficient class, and gets the unit over the polynomial ring
in those — the smallest of the five containing both 1 and the parameters.

**`Param.__init__` did not sort its terms where `Sym` does.** The contract
layer returns rows sorted, so this showed only for the two classes this layer
adds itself: `q*m([3]) + q*m([2,1])` and `q*m([2,1]) + q*m([3])` printed
differently for the same element. It sorts now, on the same key.

One refusal came out of this rather than a fix. `H̃`'s boundary encoding takes
integer numerators and puts its denominator in factored `q^a − t^b` atoms, so a
rational numerator has no slot. `(1/2)·McdHt[1] + McdHt[1]` leaked
`TypeError: 'Fraction' object cannot be interpreted as an integer` from PyO3 —
reachable before this change, since scaling by `1/2` succeeds and produces a
value addition cannot take back. `_ht_rows` states it now, naming the shape and
the coefficient. `check_parametric_scalars` asserts the refusal rather than
skipping the case. The suite went from 7562 to 7904 checks.

**What the merge still waits on is now exactly the deferred list.** The
operations `Sym` has and `Param` does not are `scalar`, `skew_by`,
`coproduct`, `expand`, `evaluate`, `principal_specialization`,
`principal_specialization_q`, `dimension`, `internal_product` and `plethysm` —
the ten `docs/plans/element-model.md` defers. Everything else the two classes
answer agrees in shape and in the exceptions it raises.

## skew_by over the four coefficient rings, 2026-08-24

The first of the ten operations `Sym` had and `Param` did not. Four entry
points — `skew_by_qt`, `skew_by_macdonald`, `skew_by_jack`, `skew_by_ht` —
over one generic `skew_ring`, on the pattern the converters and the Hopf pair
already use. The engine needed nothing: `SkewBy<C, G>` is implemented for all
six spellings of `G` at `C: Ring`, because Pieri, dual Pieri,
Murnaghan–Nakayama and Littlewood–Richardson all have integer structure
constants.

`macdonald.P([2,1]).skew_by(s([1]))` is
`(1 − t² − q²t + q²t³)/((1−qt)(1−qt²))·McdP[1,1] + McdP[2]`, which is Sage's
value after clearing signs. Hall–Littlewood, Jack, `H̃` and a scaled monomial
element all match Sage too.

**`g` keeps its own basis, and `Sym.skew_by` was changed to agree.** It used
`_same`, the coercion `+` and `*` use, so `s([2,1]).skew_by(h([1]))` raised
`BasisError` — which defeats the point of the basis argument, since that basis
selects which rule runs and not merely how `g` is read. Sage accepts any basis
here. Skewing is not a combination of two elements of one ring but an operator
built from `g` and applied to the element, so the refusal that is right for `+`
is wrong for this. `check_basis_identity` no longer lists `skew_by` among the
operations that must raise, and checks instead that the six spellings of one
`g` give one answer.

**Two checks, and the second found the bug.** The specialization —
set the parameter, then skew over ℚ through `Sym.skew_by` — shares no entry
point with the parametric route. The basis sweep runs all six rules on the same
`g`. The sweep caught a dropped argument: the Jack and `H̃` branches of `_skew`
called their entry point without the basis code, so `g` was read as a
Schur-basis element whatever it was written in. The specialization check would
not have caught it, since it only ever passed `s`.

⚠️ **The sweep has to compare by subtracting, not by `==`.** Two rules can
reach the same value over different denominators — `(1−t+q−qt)/(1−qt)` from the
Schur path and its multiple by `(1+qt)/(1+qt)` from the `e` and `p` paths — and
the fraction coefficient classes compare structurally. The difference goes
through the contract layer, which reduces. This is the same trap `_add`
already documents, met from the other side: there it forced the addition
through the boundary, here it forces the comparison through it.

The suite went from 7908 to 8993 checks.

## The Hall inner product over the four rings, 2026-08-24

`hall_inner_product_qt`, `hall_inner_product_macdonald`,
`hall_inner_product_jack` and `hall_inner_product_ht`, second of the ten. The
engine needed nothing again: `ops::hall<C, A, B>` is `C: Ring` already, because
the Schur basis is orthonormal for the pairing and the value is the sum of the
products of matching coefficients — bilinear over whatever ring they live in.

These are the first entry points that return **one coefficient** rather than
element rows, so the two cells that had no name got one: `MacCell` and
`HtCell`, a row of `macdonald_p` and of `macdonald_ht` without its partition.
`JackCell` already existed.

Four values against Sage, all exact: `⟨McdP[2], McdP[1,1]⟩` is
`(q − t)/(1 − qt)`, `⟨McdP[2], s[2]⟩` is 1, `⟨JackP[2,1], JackP[2,1]⟩` is
`(8 − 4α + 5α²)/(α+2)²`, and `⟨q·m[2], q·m[2]⟩` is `2q²`.

**`Sym.scalar` was relaxed to take the argument in any basis**, the second
`_same` that had no business being there. `h([2]).scalar(m([2]))` raised
`BasisError` and is now 1 — the `h`/`m` duality, a value rather than a
mismatch. The pairing is defined on the ring, so two spellings of one argument
give one number. Sage agrees, and `s[2].scalar(m[1,1])` is 0 there and here.
`check_basis_identity` no longer lists `scalar` either; it checks the six
spellings agree.

**The strongest check is orthonormality read backwards.** `⟨f, s_μ⟩` is the
coefficient of `s_μ` in `f`, so pairing against every shape of the degree and
comparing with `f.to("s").coefficient(mu)` puts the pairing against a change of
basis, which shares no entry point with it. Over five families and every shape
to size 4 that is most of the new checks; the suite went from 8993 to 9727.

## The coproduct over the four rings, 2026-08-24

`coproduct_qt`, `coproduct_macdonald`, `coproduct_jack` and `coproduct_ht`,
third of the ten, over one generic `coproduct_ring`. `hopf::coproduct<C: Ring>`
was already generic, because `Δ(s_λ) = Σ c^λ_{μν} s_μ ⊗ s_ν` has
Littlewood–Richardson coefficients and those carry no parameter.

`s(P[2]).coproduct()` in Sage is `−((q−t)/(qt−1))·s∅ ⊗ s11 + s∅ ⊗ s2 +
((qt−q+t−1)/(qt−1))·s1 ⊗ s1 + …`, and every coefficient matches after clearing
signs. The Jack pair `2/(α+1)` at `(1),(1)` and `(1−α)/(α+1)` at `∅,(11)` is
Sage's too, and the second is the doctest: the `α → 1/α` mirror gives its
negative.

**Both factors come back in the Schur basis, and Sage's do not.** Sage writes
`P[2].coproduct()` in `McdP ⊗ McdP`. Returning it that way needs the inverse
expansion applied to both factors of a tensor, which is not an operation here
— and `Sym.coproduct` has always returned Schur pairs whatever basis it was
handed, so matching `Sym` is what keeps the two classes converging. The value
is the same; only the basis it is written in differs.

**The check that shares nothing with the coproduct is its defining identity.**
`⟨Δf, g ⊗ h⟩ = ⟨f, gh⟩`, and the Schur basis of each factor is orthonormal, so
the coefficient at `(μ, ν)` is `⟨f, s_μ · s_ν⟩` — read through the product and
the Hall pairing, three entry points, none of them the coproduct's. The
specialization check is the second reading, and it has to drop zeros: a
coefficient can be a nonzero rational function that vanishes at the point it is
specialized to, and the integer route never builds a term for it. The suite
went from 9727 to 10148.

## expand and evaluate over the four rings, 2026-08-24

Fourth and fifth of the ten, eight entry points over two generic helpers,
`expand_ring` and `evaluate_ring`. Neither needed anything from the crate:
`Monomial::expand` is `C: Ring` and only copies coefficients, and `Schur::eval`
is `C: Ring` with the alphabet in the same ring, so an integer alphabet injects
and the parameters ride through.

**The entry points take one basis, not six.** `expand_alphabet` carries a
`src` argument and routes every basis to `m` inside; the parametric ones take
the monomial basis and nothing else, because the caller already has
`convert_qt_terms` and its three siblings to get there. Restating the routing
four more times would have been four more copies of a conversion the boundary
already exposes.

Three values against Sage, all exact after clearing signs.
`macdonald.P([2]).expand(2)` puts `(1 − t + q − qt)/(1 − qt)` on `x0 x1`,
Sage's `(qt − q + t − 1)/(qt − 1)`. `hl.P([2,1]).evaluate([1,1,1])` is
`8 − t − t²`, Sage's `−t² − t + 8`. `jack.P([2,1]).evaluate([1,1,1])` is
`(18 + 6α)/(α + 2)`, Sage's `(6a + 18)/(a + 2)`. Both evaluations degenerate to
`s_21(1,1,1) = 8`, the first at `t = 0` and the second at α = 1.

**The two check each other at the all-ones alphabet.** `f(1,…,1)` is the sum of
the coefficients of the expansion, and the two sides run different engines —
the expansion lays out the monomial basis, the evaluation runs the Schur one.
The sum is taken over the specialized values, so the addition is ℚ's rather
than the coefficient classes'.

`_ring_rows` and `_cell_coeff` came out of this: the pack-call-unpack half that
every one of these operations repeats, with the entry point chosen by the
coefficient class. The suite went from 10148 to 10486.

## dimension and the principal specialization over the four rings, 2026-08-24

Sixth and seventh of the ten. Both answer `Σ_λ c_λ w(λ)` for a weight the shape
alone decides — `f^λ` for one, `s_λ(1^n)` for the other — so the eight entry
points share `combine_ring` and differ only in which `w` they pass.

**The `u128` wall is read before any coefficient arithmetic runs.**
`shape_weights` collects every weight first and raises `OverflowError` naming
the shape, so escalation is about the ring and the wall is about the shape, and
the two cannot be confused. That is why the weights are keyed by partition
rather than positional: the parsed rows and the built term map do not iterate
in the same order.

`hl.P([2,1]).dimension()` is `2 − t − t²`, which reads off the expansion
directly: `HLP[2,1] = s[2,1] − (t + t²)·s[1,1,1]`, and `f^{21} = 2`,
`f^{111} = 1`. `jack.P([2,1]).dimension()` is `6/(α + 2)`, which is 2 at α = 1.

**The principal specialization and the evaluation now check each other.**
`f.principal_specialization(n)` and `f.evaluate([1]*n)` are the same number by
different routes — one weighs each shape by `s_λ(1^n)`, the other lays out an
alphabet and runs the Schur evaluation — and the check compares them as
coefficients, not after specializing, so it is the parametric values that have
to agree. With the expansion's sum from the previous change that is three
routes to one number. The suite went from 10486 to 10714.

## The principal specialization in q, and the wall it runs into, 2026-08-24

Eighth of the ten, and the first that is **not** available over every ring.
`principal_specialization_q_qt` is one entry point rather than four, because
the operation introduces a variable and the coefficient classes here carry at
most two: `Poly` one, `QtPoly` two, and the three fraction classes none to
spare. So it is written where the base ring is a single variable other than
`q` — Hall–Littlewood and LLT, over `ℚ[t]` — and refused elsewhere.

`hl.P([2,1]).principal_specialization_q(3)` is
`q + 2q² + 2q³ − q³t − q³t² + 2q⁴ + q⁵`. Sage writes the same value as
`q^5 + 2q^4 + (−t² − t + 2)q³ + 2q² + q`, over `ℚ(t)` with the default `q`.

**This is Sage's own wall, reported the same way.** Over `ℚ(q,t)` Sage says
"the variable q is in the base ring, pass it explicitly" and takes any ring
element as `q` — `P[2].principal_specialization(3, q=t)` answers. Here the two
refusals are separated, because they are two different facts: `q` already being
a parameter, and a coefficient class having no free variable at all. The Jack
case is the second, and calling it the first would have been wrong — `ℚ(α)`
does not carry `q`.

The boundary guards it too: the entry point refuses rows whose `q` exponent is
nonzero, naming the shape and the exponent pair. That is the same fact stated
where it can be checked rather than trusted.

**What makes the refusal actionable is not written yet.** Sage's advice is to
pass the variable, which for a `ℚ(q,t)` element means substituting an existing
ring element. That is `evaluate` at an alphabet drawn from the base ring, and
`evaluate` here takes integers. Recorded rather than done.

Two checks on the family that does work: at `q = 1` it is
`principal_specialization(3)`, and at `t = 0` it is the classical q-analogue,
since `P_λ(x; 0) = s_λ`. Both are ways of confirming the introduced `q` and the
`t` already there stayed apart. The other four families are checked to refuse.
The suite went from 10714 to 10781.

## The internal product over the four rings, 2026-08-24

Ninth of the ten, and the first that needed the ring widened rather than
carried. `ops::internal<C: QAlgebra>` routes through the power-sum basis and
divides by z_μ, and of the four boundary rings only `AFrac<Guarded>` and
`Ratio<Rational>` are `QAlgebra` — `QtPoly<Guarded>` and `Frac<Guarded>` are
rings without ℚ in them.

**The fix is to compute over ℚ and answer in ℤ, which is what Sage's base ring
does implicitly.** `internal_product_qt` builds over `QtPoly<GuardedRat>` and
`internal_product_macdonald` over `Frac<GuardedRat>`, both escalating to the
`BigRational` width. The answer is a ℤ-bilinear combination of the arguments,
so the denominators cancel; `qt_poly_integral` and `mac_coeff_integral` raise
rather than round if one does not, on the model of `dump_integral`. Jack and
`H̃` needed no widening at all — `AFrac<C>` is a `QAlgebra` for every `C`,
because α is an indeterminate and dividing by z_μ never asks for its inverse.

`build_mac_rat` is the one new builder. Everything else reuses `build_qt` and
`build_qt_wide`, which were already over the rational widths.

Four values against Sage, all exact after clearing signs.
`HLP[2,1] ∗ HLP[2,1]` agrees in all three coefficients, the largest being
`1 + t − t² − 3t³ − 2t⁴ + t⁵ + 2t⁶ + t⁷`. `McdP[2] ∗ McdP[1,1]`,
`McdHt[2] ∗ McdHt[2]` and `JackP[2,1] ∗ JackP[2,1]` likewise; the Jack one
matches after scaling, where Sage writes halves in the denominator.

⚠️ **The Jack contract-layer value is unreduced, and that is by design.**
`internal_product_jack` on `JackP[2,1]` returns `([3], [], 3)` for each
coefficient — three thirds, not one. `AFrac` normalizes its atoms and not its
integer content, because cancelling the content needs a gcd inside `C` that
`Ring` does not offer; `src/afrac.rs`'s module doc records that. The
convenience layer's answer is reduced, because the inverse expansion on the way
back normalizes it, but the entry point's doctest shows the raw form.

**Unlike `scalar` and `skew_by`, this one keeps the same-basis refusal.** It
combines two elements of the ring rather than pairing them or building an
operator, so it is in the family `+` and `*` belong to.

Three checks: `h_n` is the Kronecker identity in degree n — a fact about the
operation and not about any coefficient, so it holds over every ring —
symmetry in the two arguments, and the specialization. The suite went from
10781 to 11921.

## Plethysm, over three of the four coefficient rings (2026-08-24)

The tenth and last of the operations `Sym` had and `Param` did not.
`plethysm_qt`, `plethysm_macdonald` and `plethysm_ht` in `src/python.rs`, on
the same per-ring pattern as the nine before, over new `Plethystic` impls for
`Frac` (`src/frac.rs`) and `Ratio` (`src/deltaop.rs`). `Param.plethysm` and
`_plethysm` in `python/symfn/_families.py` are the convenience half.

**The parameters are part of the alphabet, so `p_n` raises them.** That was
already `QtPoly::frobenius`'s convention and it is Sage's default; the two new
impls extend it to the denominators, where `1 − qᵃtᵇ ↦ 1 − q^{an}t^{bn}` and
`qᵃ − tᵇ ↦ q^{an} − t^{bn}` keep both families closed. Sage's `exclude=`,
which holds a variable constant instead, has no counterpart here.

`plethysm_qt` and `plethysm_macdonald` run over `ℚ[q,t]` and `ℚ(q,t)` and
answer in the integral ring, for the reason the internal product does above.

**There is no `plethysm_jack`.** Over ℚ(α) the Frobenius is α ↦ α^n, which
takes a denominator `α + 1` to `α² + 1` and so leaves the
product-of-linear-forms class `AFrac` holds. `Param.plethysm` refuses that ring
by name and points at `.at()`. The full account, with the Sage values showing
the obstruction is mathematical rather than an encoding artifact, is in
[jack.md](jack.md).

Values against Sage, all exact: `HLP[2][t·HLP[1]] = t²·HLP[2]` — `t²` and not
`t` is the value that separates the raising convention from its rival —
`HLP[2][HLP[1,1]]`, `McdP[2][McdP[1,1]]`, `McdHt[2][McdHt[2]]` (all three
coefficients, after clearing signs from Sage's expanded denominators), and
`s_2[s_1/(1−qt)] = (s_2 + qt·s_11)/((1−qt)(1−q²t²))`, where the `1 − q²t²` is
`p_2`'s raised copy and a Frobenius that left the denominator alone would give
`(1−qt)²`.

**Specializing does not commute with plethysm**, which is what makes the
convention checkable at all: `t·s_1` composed into `p_2` gives `t²p_2`, while
setting `t = 3` first gives `3p_2` and not `9p_2`. So
`check_parametric_plethysm` crosses to the integer route only with an inner
argument whose coefficients carry no parameter, and pins the raising with a law
instead — `f[t^k·g] = t^{kd}·f[g]` for `f` homogeneous of degree `d` — plus
linearity and multiplicativity in the outer argument, both of which run through
machinery plethysm does not share. The suite went from 11921 to 11997.

`_plethysm` restores the outer argument's cleared denominator and refuses the
inner one's, which is what `Sym.plethysm` already does: plethysm is linear in
`f` and not in `g`.

## The principal specialization at a base-ring alphabet (2026-08-24)

`principal_specialization_q` introduces a fresh `q` and so needs a free
variable in the coefficient ring. Hall-Littlewood and LLT have one; `ℚ(q,t)`
and `ℚ(α)` do not, and the entry point refused them. That refusal was accurate
and not actionable — it pointed at "evaluate at an alphabet you name yourself",
and `evaluate` took integers only.

Sage's message says what to do instead: *pass it explicitly*. What that means
is that the alphabet is drawn from **the base ring**, not from a new variable —
`P[2].principal_specialization(3, q=q)` substitutes the ring's own `q`. That is
a different operation from the one this tree had, and it is available in every
ring, including ℚ(α), which has no free variable at all.

`principal_specialization_at_{qt,macdonald,jack,ht}` in `src/python.rs` are
that operation, over one generic `ps_at_ring`. `s_λ(1,q,…,q^{n−1})` is a
polynomial in `q` with non-negative integer coefficients, which
`crate::eval::principal_specialization_q` already returns, so substituting a
ring element for `q` is ring arithmetic and the bound stays at `Ring` — no
widening, no escalation past the usual pair. The powers of the alphabet are
shared across shapes.

The alphabet argument is one coefficient in that ring's encoding, and it is
parsed and built through the same path a row of the element takes: a one-term
element at the empty partition. A malformed cell therefore raises where a
malformed row would, and `one_coefficient` reads the single value back out.

`Param.principal_specialization` and `Sym.principal_specialization` both took
`n` alone and now take `n, q=None`, `q = None` meaning the value at `1^n` they
already answered. Keeping the two signatures identical is one fewer difference
for the merge. `Sym`'s route needs no entry point of its own: its base ring is
ℚ, so substituting is Python arithmetic over the same q-analogue.

⚠️ **A non-integral alphabet is refused, and the encodings force it.** `Frac`'s
denominator is a product of binomials `1 − qᵃtᵇ` and `Ratio`'s a product of
atoms, so neither holds `1/2`; the `Poly`/`QtPoly` path clears denominators by
scaling, and the alphabet enters at every power from 0 to the degree rather
than linearly, so a cleared scale cannot be restored. `_coeff_cell` raises and
says so. `AFrac` carries an integer scale, so Jack does take `q = 1/2`.

Values against Sage, exact: `McdP[2]` at `q = q` in 3 variables is
`(1 + q − 2qt + 3q² − 2q²t + 2q³ − 3q³t + 2q⁴ − q⁴t − q⁵t)/(1 − q*t)`, which is
Sage's after clearing signs; `HLP[2,1]` at `q = t`; `McdHt[2]` at `q = q`; and
`JackP[2]` at `q = 2`, which Sage refuses as a plain integer and this layer
lifts into the ring.

`check_parametric_principal_at` crosses three ways: `q = 1` must give
`principal_specialization(n)`, `q = c` must give `evaluate([1, c, …, c^{n−1}])`
— a route that lays the alphabet out and expands in the monomial basis, sharing
nothing with the q-analogue — and specializing the parameters afterwards must
agree with specializing first, since `at` is a ring homomorphism and the
alphabet is a ring element like any other. The suite went from 11997 to 12189.
