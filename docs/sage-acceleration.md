# Accelerating the Sage ecosystem: where the time actually goes

*Profiling study, 2026-07-28, against SageMath 10.9 on this machine.*

**Question:** if we wanted to accelerate as much of the Sage combinatorics
ecosystem as possible, what would we build? Could an efficient kernel for a few
key operations have broad impact?

**Answer:** yes, but the kernel that matters is not the obvious one. Sage's
hottest symmetric-function path is *already* C, and over half its cost is
converting that C output into Python objects. The broadly-shared bottleneck is
the term store and the key objects, not the mathematics.

Companion to `research-gaps.md`, which covers capabilities that do not exist in
any package. This document is about making capabilities that *do* exist run
faster.

---

## 1. The mathematics is often already C

Profile of 20 evaluations of `s[5,4,3,2] * s[6,4,2]` (`cProfile`, sorted by
`tottime`):

```
 ncalls  tottime  cumtime  function
     20    0.028    0.038  sage/libs/lrcalc/lrcalc.py:195(_lrcalc_dict_to_sage)
     20    0.002    0.041  sage/libs/lrcalc/lrcalc.py:272(mult)        <- the C call
   9100    0.004    0.005  sage/combinat/combinat.py:1528(__init__)
   9100    0.002    0.008  sage/combinat/partition.py:519(__init__)
   9102    0.001    0.002  sage/combinat/partition.py:551(__hash__)
     20    0.001    0.053  sage/combinat/free_module.py:1054(linear_combination)
```

Sage's `schur.product_on_basis` calls **lrcalc**. The Littlewood–Richardson
computation itself is compiled C and is not the bottleneck.
`_lrcalc_dict_to_sage` — which walks the C result dict and builds a Sage
`Partition` for every key — is **0.028s of 0.054s total, ~52%**.

**Consequence:** a faster LR algorithm attacks the remaining ~48%, and roughly
half of *that* is also object layer. This is the single most important fact in
this document.

### The same pattern in the pure-Python paths

Conversions and plethysm show the identical shape, with the coercion entry point
`sf.py:1698(__call__)` replacing lrcalc as the top line:

| workload | top cost | share of cumtime |
|---|---|---|
| `s → p`, `s[7,5,3,1]` | `sf.py:1698(__call__)` | 0.014 / 0.021 (67%) |
| `s → m`, `s[7,5,3,1]` | `sf.py:1698(__call__)` | 0.050 / 0.090 (56%) |
| `h[3][h[3]]` plethysm | `sf.py:1698(__call__)` | 0.019 tottime, 620 calls |

Beneath those, in every case: thousands of `Partition._element_constructor_`
calls, thousands of `Partition.__contains__` validations, and tens of thousands
of the weakly-decreasing-check generator expression at `partition.py:6528`.

The `s → m` profile additionally shows Sage's **category framework running at
computation time** — `dynamic_class.py:338(dynamic_class_internal)` at 0.013s,
plus `category.py:_all_super_categories`, `category.join`, and
`is_subcategory`. That is pure framework overhead on a numeric path.

---

## 2. One substrate underlies the whole ecosystem

Checked by walking the MRO of each parent:

| algebra | backed by `CombinatorialFreeModule` |
|---|---|
| `Sym.s` | ✅ |
| `QSym.F` | ✅ |
| `NCSF.R` | ✅ |
| `NCSym.m` | ✅ |
| `FQSym.F` | ✅ |
| `WQSym.M` | ✅ |
| `FSym.G` | ✅ |
| `SchubertPolynomialRing` | ✅ |
| `WeylCharacterRing("A3")` | ✅ |
| `SymmetricGroup(4).algebra(QQ)` | ✅ |

11 of 11. This is one shared bottleneck, not ten independent ones — which is
precisely what makes a kernel investment leverage the whole ecosystem rather
than one corner of it.

---

## 3. The tax is key-object construction

| operation | µs/op | vs. baseline |
|---|---|---|
| `tuple([5,4,3,2,1])` — baseline | 0.03 | 1x |
| `Partition([5,4,3,2,1])` | **1.71** | **57x** |
| `_Partitions([5,4,3,2,1])` — internal "fast path" | 1.47 | 49x (saves only 14%) |
| `list(Partition)` | 0.15 | 5x |
| `Composition([2,1,2])` | 1.65 | 55x |
| `Permutation([3,1,2])` | 2.04 | 68x |
| `s._from_dict` (2 terms) | 0.43 | — |
| `s[5,4,3,2]` — build a *one-term* element | **4.71** | — |

Every key of every element of every algebra in §2 costs ~1.7–2 µs to construct,
against 0.03 µs for the tuple carrying identical information. Constructing a
single-term Schur element costs 4.71 µs before any arithmetic happens.

Note that Sage's own internal fast path, `_Partitions`, recovers only 14% — the
cost is in `Element.__init__`, parent bookkeeping, and validation, not in the
public constructor's dispatch.

---

## 4. The rest of the ecosystem is worse than Sym

Per-product cost on deliberately tiny inputs:

| product | µs/op |
|---|---|
| **QSym `F[2,1]·F[1,2]`** | **10,229** |
| QSym `M[2,1]·M[1,2]` | 268 |
| WQSym `M` | 247 |
| FQSym `F` | 151 |
| NCSym `m` | 46 |
| S₅ group algebra | 23 |
| Weyl character ring A₃ | 20 |
| Schubert `X[2,1,3]·X[1,3,2]` | 18 |
| NCSF `R[2,1]·R[1,2]` | 8 |

The QSym `F` outlier was re-measured with **60 distinct pairs** to rule out a
memoization artifact; it got *worse* (10,229 µs distinct vs 3,175 µs for one
repeated pair — so the 4,711 µs first measurement was itself flattered by
caching).

Profile of 40 such products:

```
 ncalls  tottime  function
   6334    0.074  integer_lists/lists.py:225(_element_iter)
  82679    0.059  composition.py:187(__init__)
  54063    0.053  composition.py:1766(_element_constructor_)
  82679    0.042  combinat.py:1528(__init__)
  35086    0.038  composition.py:360(__add__)
 220977    0.027  composition.py:1780(<genexpr>)          <- validation
  95516    0.022  combinat.py:1359(__hash__)
  10074    0.021  shuffle.py:605(__iter__)
 296407    0.020  {built-in method builtins.isinstance}
```

**~2,000 `Composition` objects constructed per product of two degree-3
compositions.** Not one entry in the top nine is mathematics.

---

## 5. The Amdahl trap

This project has already measured this effect without naming it. From
`docs/record/README.md`:

> End to end *through Sage* the same substitution is **1.84x** like-for-like, or
> 4.37x with a Sage-`Partition` cache that Symmetrica's wrapper does not have.
> The gap between those is object marshalling, which both backends pay.

That is §3 of this document seen from the other side. The general statement:

> **The faster the kernel gets, the larger the marshalling fraction becomes.**

A 30x algorithmic win behind a boundary that costs 52% delivers under 2x
end-to-end. The object layer is therefore not an afterthought to be tidied up
later — it is the **prerequisite** that determines the ceiling on everything
else.

It also means benchmark hygiene matters here in the specific way
`docs/record/README.md` already documents: an in-process A/B that measures the
kernel alone will overstate the user-visible win, and the honest number is
always end-to-end through Sage on cold, distinct inputs.

---

## 6. What to build, in order

### Tier 0 — the term store *(prerequisite; multiplies every other tier)*

Cross the Python boundary **once per expansion, never per term**.

- Carry results as compact arrays — small-integer parts plus a coefficient
  array — not as lists of Sage objects.
- Materialize Sage `Partition` / `Composition` / `Permutation` objects **lazily**,
  only for terms the caller actually touches.
- **Intern** key objects on a tuple key. At 1.71 µs construction against ~0.1 µs
  for a dict hit, this is a >10x saving on any workload that revisits keys — and
  every conversion, product, and plethysm revisits keys constantly.

This is the only tier that helps QSym, Schubert, and Weyl character rings for
free, because it operates on the substrate identified in §2. It can also be done
entirely on our side of the boundary without asking Sage to change anything.

### Tier 1 — operations with no C behind them at all

Plethysm, Kronecker, `s ↔ m`, `s ↔ p`, Jack, LLT, the `st` (Orellana–Zabrocki)
basis. These are pure Python end to end, so headroom is 10–1000x rather than
2–3x. They are also exactly the operations that timed out in the survey recorded
in `research-gaps.md`:

- Kronecker `n=32`: >90s
- Jack `n=15`: >90s
- `st[6,4]·st[6,4]`: >90s
- plethysm `h₄[h₁₀]`: >120s

### Tier 2 — the other combinatorial Hopf algebras

QSym, WQSym, FQSym, NCSym. Structurally these are **the same problem already
solved here**: a graded connected Hopf algebra with combinatorial keys,
products, coproducts, and basis conversions. The `SymFn` / `SymAlgebra` trait
design generalizes to compositions and set partitions with the
basis-confusion-is-a-compile-error discipline intact.

The 10.2 ms QSym `F` product is the largest single headroom figure in this
investigation.

### Tier 3 — Schur and skew LR

Already lrcalc (§1). Headroom is real — `SkewLr` beats lrcalc on most shapes
measured in `docs/record/python-and-sage-interop.md` — but bounded, and halved
again by the marshalling tax until Tier 0 lands.

### Flagged, but different in kind

- **Nonsymmetric Macdonald in affine types.** Enormous headroom, but requires
  the whole root-system stack and is Ram–Yip-shaped, not `SkewLr`-shaped.
  High cost.
- **Crystals, root systems, Weyl character rings.** `CombinatorialFreeModule`-backed
  (so Tier 0 helps), but the hot loop is weight arithmetic and Weyl-orbit
  enumeration — a genuinely different kernel from anything in this codebase.
- **Schubert / Grothendieck polynomials.** 18 µs/product, CFM-backed, with real
  algorithmic content (transition formula, Monk's rule). Moderate on both axes.

---

## 7. Caveats

- **Single runs on a laptop, coarse timers.** Treat every number here as
  order-of-magnitude. `docs/record/README.md` documents a 2.1x battery/thermal
  artifact on this same machine; the ratios and the *profile shapes* are the
  durable quantities, not absolute times.
- **`tottime` vs `cumtime`.** The percentages in §1 compare `tottime` of the
  marshalling function against `cumtime` of the whole operation, which is the
  right comparison for "how much of this could a better boundary remove" but is
  not a partition of the total.
- **The marshalling fraction is input-dependent.** For lrcalc specifically the
  underlying enumeration is O(tableaux) while marshalling is O(terms), so at very
  large shapes the C algorithm reclaims share. The 52% figure is for a small
  product. Conversely, against a *fast* kernel (where the algorithm no longer
  dominates), marshalling share goes up — see §5.
- **Replacing `CombinatorialFreeModule`'s backing store upstream is a large
  project** with real review and social cost, touching every corner of
  `sage.combinat`. Nothing in Tier 0 requires it. The incremental path already
  proven by `scripts/sage_backend.py` — a drop-in backend for specific
  operations, shipped one at a time, each independently validated — is the
  shippable one.

---

## 8. Bottom line

An efficient kernel can have broad impact, but the leverage does not come from
faster mathematics. It comes from:

1. a **compact term representation** that crosses the Python boundary once, and
2. **interned, lazily-materialized key objects**,

both of which sit on the one substrate — `CombinatorialFreeModule` plus
`Partition`/`Composition`/`Permutation` — shared by every algebra in §2. Do that
first, then port operations behind it, cheapest-and-most-starved first: the
pure-Python symmetric-function operations, then the other Hopf algebras, and
only then the LR paths that already have C underneath them.

### Reproducing

Profiling scripts used for this document are ad hoc and were not committed. The
measurements are: `cProfile` over 20 repetitions for §1, wall-clock over
20 000–200 000 repetitions for §3, 500–2 000 repetitions for §4, and a
distinct-input rerun for the QSym outlier. Sage 10.9, Python 3.14, macOS,
single-threaded.
