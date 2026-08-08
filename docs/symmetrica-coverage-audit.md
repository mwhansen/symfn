# What Sage actually calls in Symmetrica — a coverage audit

*Audit of `sagemath/sage` at `09472ff` (10.10.beta7, 2026-07-26): 3053 sagelib
source files scanned against the 66 entry points exported by
`sage/libs/symmetrica/all.py`. Companion to
[docs/sage-packaging-audit.md](sage-packaging-audit.md), which answered the
packaging question; this one answers the coverage question that
[docs/release-readiness.md](release-readiness.md) Phase 5c called the highest-
value unknown.*

**Headline: Sage reaches 36 of Symmetrica's 66 exported entry points, from six
files. All 36 are now computed by symfn; 29 of them are intercepted by
the adapter today, across five of the six consumer files.** That adapter was
`scripts/sage_backend.py` when this was written and is now `sage/libs/symfn/`
inside Sage; `scripts/check_backend.py` is what drives the comparison either
way.**

That is a claim about what Sage calls, and only that. The other 30 entry points
Symmetrica exports ([§3](#3-the-30-unreached-entry-points)) stay reachable by
anyone who imports `sage.libs.symmetrica` directly, and symfn does not cover
them. "symfn covers every Symmetrica entry point Sage calls" is the accurate
sentence; "symfn replaces Symmetrica" is not.

> ⚠️ **Superseded in four places.** The original audit's gap list was read off
> function *names* rather than off what Sage does with the return value, and two
> of its three claims did not survive implementation. What it called a binding
> gap was a mathematics gap; what it called a one-to-one mapping had one entry
> that maps to a different operation. The fourth correction is this audit's own
> decision not to reimplement that entry, reversed in
> [§2](#the-schubert-seven-do-not-map-one-to-one). The corrections are marked ⚠️
> inline and collected in [§5](#5-what-implementation-corrected); the counts
> above are post-correction.

This is a much smaller displacement target than Symmetrica's own scope
statement suggests. Symmetrica advertises modular and projective representation
theory, classical groups, Hecke algebras and finite group operations
([the record index](record/README.md)'s "Beyond the core" list, near-verbatim). **Sage
calls none of it.**

---

## 1. The consumer surface is six files

Every reference to `sage.libs.symmetrica` in sagelib, by import form:

| file | form | what it uses |
|---|---|---|
| `combinat/sf/classical.py:55` | `import ... all as symmetrica` | the 20 basis conversions, via `getattr` |
| `combinat/sf/sfa.py:5656` | `import ... all as symmetrica` | the 5 `compute_*_with_alphabet`, via `getattr` |
| `combinat/sf/hall_littlewood.py:28` | `from ...symmetrica import (...)` | `hall_littlewood` |
| `combinat/sf/monomial.py:22` | `import ... all as symmetrica` | `mult_monomial_monomial` |
| `combinat/tableau.py:116` | `lazy_import` | `kostka_number`, `kostka_tab` |
| `combinat/schubert_polynomial.py:89` | `lazy_import` | 7 Schubert entry points |

Plus two non-consumers: `libs/all.py` (re-exports the module into Sage's global
namespace) and `misc/citation.pyx` (a one-line citation-registry entry,
`systems['Symmetrica'] = ['sage.libs.symmetrica']`).

Everything else that matches "symmetrica" in sagelib is the word
*symmetrically* in prose. `sagelib` is also the only package in Sage's 441 that
depends on the symmetrica spkg at all.

## 2. The 36 reachable entry points, and symfn's coverage

Two of the six call sites build names dynamically, so they don't appear as
literal strings — `classical.py` does `getattr(symmetrica, f't_{other_name}_{name}')`
over `{Schur, monomial, homogeneous, elementary, powersum}`, and `sfa.py` does
`getattr(symmetrica, f'compute_{basis}_with_alphabet')`. Those account for 25 of
the 36.

"Covered" below means symfn computes it *and* the adapter intercepts the
call; "library only" means symfn computes it but Sage still reaches the C
function.

| group | count | state | notes |
|---|---|---|---|
| `t_{FROM}_{TO}` conversions | 20 | ✅ covered | `convert.rs`, intercepted via the conversion table |
| `compute_*_with_alphabet` | 5 | ✅ covered | ⚠️ **was mis-classified** — needed `Monomial::expand`, not a binding of `eval()`; see §5 |
| `hall_littlewood` | 1 | ✅ covered | `hl.rs`; the adapter now intercepts the direct import too |
| `mult_monomial_monomial` | 1 | ✅ covered | ⚠️ the pre-existing `Monomial::mul` was in `convert.rs` and went through Schur; now a direct rule in `sym.rs` |
| `kostka_number` | 1 | ✅ covered | `kostka.rs` |
| `kostka_tab` | 1 | ✅ covered | was the one true gap; now `kostka::semistandard_tableaux` |
| Schubert | 7 | 7 library only | ⚠️ **not one-to-one** — `scalarproduct_schubert` is a different operation; see below |
| **total** | **36** | **29 covered, 7 library only** | |

### The Schubert seven do *not* map one-to-one

`docs/record/schubert.md` designed the bindings against this list, and six
of the seven land exactly — 1679 comparisons against Symmetrica over `S₁`–`S₄`,
with **zero disagreements on any input Symmetrica answers at all**. The seventh
does not map:

| Symmetrica | symfn | checked |
|---|---|---|
| `mult_schubert_schubert` | `schubert_multiply` | 1089 pairs, exact |
| `mult_schubert_variable` | `schubert_multiply_variable` (⚠️ 1-based; Symmetrica's is 0-based) | 132, exact |
| `t_SCHUBERT_POLYNOM` | `schubert_expand` | 33, exact |
| `t_POLYNOM_SCHUBERT` | `polynomial_to_schubert` | 32, exact |
| `divdiff_schubert` | `schubert_divided_difference` | 132, exact where Symmetrica answers |
| `divdiff_perm_schubert` | `schubert_divided_difference_perm` | 264, exact where Symmetrica answers |
| `scalarproduct_schubert` | `schubert_scalar_product` | 617 single-term pairs over S₁–S₄, exact |

**`scalarproduct_schubert` returns a Schubert polynomial, not a scalar.** The
original table paired it with `schubert_pairing`, which returns an integer — the
Poincaré pairing on `H*(Fl(n))`. They are not the same map, and the smallest
witness takes one line:

```text
  X([2,1]).scalar_product(X([2,1]))  =  X[1,3,2]     (Symmetrica)
  symfn.schubert_pairing([([2,1],1)], [([2,1],1)], 2)  =  0
```

Sage reaches it from `SchubertPolynomial.scalar_product`
(`schubert_polynomial.py:350`), so it is one of the 36, and it was for a while
the only remaining mathematical gap in the whole displacement. It is now
[`Schubert::scalar_product`](../src/schubert.rs) and
`symfn.schubert_scalar_product`.

**What the map is.** `∂_{w₀⁽ⁿ⁾}(S_u · S_v)`, and the Poincaré pairing is its
coefficient at the identity — the two operations agree on nothing else. `n` is
an explicit argument here, as it is for `schubert_pairing`.

Two further obstacles to intercepting this file, recorded so they are not
rediscovered:

- **symfn is more total than Symmetrica here.** Of the 1679 comparisons, 187 are
  cases where Symmetrica or its Sage wrapper raises `ValueError` — `∂_i` with
  `i` past the permutation's length, `∂_w` on the identity — and symfn returns
  the mathematically correct value instead. A faithful drop-in has to
  *reproduce the exceptions*, because Sage's doctests assert them.
- **Symmetrica leaks and then blocks.** After a run of `divdiff_perm_schubert`
  calls, Symmetrica prints `ERROR: permutation memory not freed?:
  mem_counter_perm = 99` at teardown and drops into an **interactive prompt**
  (`enter a to abort with core dump, g to go, …`). On a non-tty that hangs
  forever; the harness has to close stdin. This is a robustness argument for
  displacement, and it is why the Schubert comparison is a standalone probe
  rather than part of `check_backend.py`.

**Decided: `scalarproduct_schubert` will not be reimplemented.** Sage keeps
Symmetrica as an *optional* package rather than losing it, so the operation
stays available to anyone who installs it, and the displacement target becomes
"demote Symmetrica from standard to optional" rather than "remove it". See
[release-readiness](release-readiness.md) Phase 5c for what that changes.

The fact that makes this cheap: **nothing in sagelib calls
`SchubertPolynomial.scalar_product`.** Its only references are its own
definition (`schubert_polynomial.py:323`) and its own three doctest lines. It is
public API a user could invoke with zero internal dependents — structurally the
same position as the 30 unreached entry points in §3, and it wants the same
answer.

⚠️ **That decision is reversed, and the paragraph above is kept because its
premise is still true.** Nothing in sagelib does call the method. What the
argument missed is that the method is not in the same position as the 30
unreached entry points: those are reachable only by writing
`from sage.libs.symmetrica.all import ...`, where a user knows they are calling
Symmetrica, while `scalar_product` is a public method on the user-facing
`SchubertPolynomial` class, one of the six consumer files. Third-party code can
therefore depend on it without its author ever knowing Symmetrica was
underneath — which is the case "no internal dependents" does not cover, and the
one that decides it. The operation is implemented, and no operation a user could
already be calling is left behind.

The adapter does not intercept `schubert_polynomial.py` today. **That is
sequencing, not a second decision.** The argument for waiting was that wiring
six of the seven buys nothing while the file loads Symmetrica for the seventh
regardless; with all seven available that argument is gone, and wiring the file
is what takes Schubert polynomials off Symmetrica entirely. The
exception-fidelity requirement above still has to be met first, and it is the
whole remaining cost: `scalar_product` itself needs the adapter to pass the rank
Symmetrica infers, which is the longest one-line form among the two arguments
after Sage strips trailing fixed points (measured over 1089 pairs; see
[docs/record/schubert.md](record/schubert.md)).

## 3. The 30 unreached entry points

Exported by `all.py`, and therefore public API a user could call directly, but
**no sagelib code path reaches them**:

```
bdg  chartafel  charvalue  compute_schur_with_alphabet_det  dimension_schur
dimension_symmetrization  gupta_nm  gupta_tafel  kostka_tafel  kranztafel
mult_schur_schur  ndg  newtrans  odd_to_strict  odg  outerproduct_schur
part_part_skewschur  plethysm  q_core  random_partition  scalarproduct_schur
schur_schur_plet  sdg  specht_dg  start  strict_to_odd_part
t_POLYNOM_ELMSYM  t_POLYNOM_MONOMIAL  t_POLYNOM_POWER  t_POLYNOM_SCHUR
```

(`start` is the library initializer, called by `all.py` itself.)

Three observations:

- **symfn already covers much of this incidentally** — `plethysm`,
  `mult_schur_schur`, `outerproduct_schur`, `part_part_skewschur`,
  `dimension_schur`, `charvalue`/`chartafel`, `newtrans` all have direct
  equivalents (`plethysm.rs`, `schur_multiply`, `skew_schur`, `dimension`,
  `character_value`/`character_table`, `schubert_to_stanley_schur`). The
  overlap is not the problem.
- **The genuinely uncovered ones are the representation-theory group** — `bdg`,
  `sdg`, `odg`, `ndg`, `specht_dg`, `dimension_symmetrization`, `kranztafel`
  (wreath products), `gupta_nm`/`gupta_tafel`, `q_core`,
  `strict_to_odd_part`/`odd_to_strict`. This is `docs/record/README.md`'s "Beyond the core"
  territory, and it is **exactly the part Sage never calls**.
- Removing them is a **deprecation question, not an implementation question**.
  They are reachable by users via `from sage.libs.symmetrica.all import ...`,
  so displacement means either reimplementing them, or deprecating them through
  Sage's normal cycle. That decision belongs upstream, not here.

## 4. What this means for the plan

The work these findings imply is tracked in
[docs/release-readiness.md](release-readiness.md) under **Phase 5b →
"Becoming a complete drop-in"** — deliberately there and not in Phase 5c,
because every item is ordinary library or binding work needing no upstream
involvement:

| task | Sage call site | state | what it actually took |
|---|---|---|---|
| `compute_*_with_alphabet` | `sf/sfa.py:5656` | done | ⚠️ **new mathematics**, not a binding — `Monomial::expand` + `expand_alphabet` |
| `mult_monomial_monomial` | `sf/monomial.py:129` | done | direct overlay rule in `sym.rs` + `monomial_multiply` |
| `kostka_tab` | `tableau.py:7016,7036` | done | `kostka::semistandard_tableaux` + `semistandard_tableaux` |
| wire `hall_littlewood` | `sf/hall_littlewood.py:28` | done | adapter rebinds the module-level name |
| `scalarproduct_schubert` | `schubert_polynomial.py:350` | done | new gap this audit missed, and a decision it then got wrong — `∂_{w₀⁽ⁿ⁾}(S_u·S_v)`, with `n` explicit; see §2 |

the adapter now displaces **five of the six consumer sites**.
`scripts/check_backend.py` covers the five: **8647 computations at degree 8, 0
mismatches**, up from 4678 covering the conversion table alone.

**The 30 entry points that stay behind** are §3's: exported, never reached from
inside Sage, and reachable by a user only through an explicit
`sage.libs.symmetrica` import. What to do about them — reimplement, deprecate,
or keep Symmetrica installable as an optional package — is a Sage-side decision
and is not settled here. `scalarproduct_schubert` is no longer among them, and
it is the one that could not have been left to that decision, because a caller
reaches it without naming Symmetrica.

**The revised risk ranking for Phase 5c** is therefore: platform reach of the
wheel (from [the packaging audit](sage-packaging-audit.md)) > upstream appetite
for displacing Symmetrica on the paths Sage uses. Coverage is no longer on that
list: every entry point Sage calls is computed here.

## 5. What implementation corrected

Every claim in §2 that implementation touched, and how it moved.

> **Two of this audit's three gap classifications were wrong, and both were
> wrong in the same direction.** The audit matched a Symmetrica name to a symfn
> name and stopped there, where what the caller does with the return value was
> the deciding fact — which is what a name-based audit cannot see.

1. **`compute_*_with_alphabet` was called "a binding gap, not a mathematics
   gap". It was the reverse.** `eval()` evaluates at an alphabet of `Ring`
   elements and returns one value; `sfa._expand` needs a *polynomial* in `n`
   indeterminates, which it feeds to `resPR(...)`. No amount of binding turns
   the first into the second. The fix was a new operation, `Monomial::expand`,
   returning exponent vectors — which then also gave `mult_monomial_monomial`
   its multiset machinery for free.
2. **`Monomial::mul` was said to be in `sym.rs`. It was in `convert.rs`**, and
   it went `m → s → LR → m`, which costs the whole degree because `m → s`
   inverts the Kostka matrix. The classification ("binding only") was right; the
   location and the cost were not. It is now the direct rule in `sym.rs`, with
   the Schur route retained as the reference oracle.
3. **The Schubert seven were said to map one-to-one. Six do.** See §2.
4. **This audit then decided the seventh would not be reimplemented, and that
   was wrong.** The reasoning — no internal dependents — was sound and is still
   true; it was applied to the wrong category. A public method on a user-facing
   class is not in the position of a low-level export, because the caller need
   not know which library answers. See §2 for the reversal, and
   [docs/record/schubert.md](record/schubert.md) for what implementing it
   found.

What the audit got right and implementation confirmed: `kostka_tab` was a real
gap and the only one it identified correctly; the 20 conversions, `kostka_number`
and `hall_littlewood` needed no library work; and the whole
representation-theory half of Symmetrica is genuinely never called.

One contract the audit did not anticipate at all: **`kostka_tab`'s order is part
of the interface.** Sage's `SemistandardTableaux(λ, μ).list()` returns the
backend's list verbatim and its doctests print it, so the enumeration order is
visible to callers. It is increasing lexicographic in the row-major reading
word — which the natural chain-of-horizontal-strips walk does *not* produce, so
`semistandard_tableaux` sorts. Checked against Symmetrica over all 1818 (λ, μ)
pairs through degree 9.

---

### Reproducing

```
git clone --filter=blob:none --sparse --depth 1 https://github.com/sagemath/sage
cd sage && git sparse-checkout set src/sage
grep -rnE "(from|import)[^\n]*sage\.libs\.symmetrica|lazy_import\([^)]*symmetrica" \
     --include="*.py" --include="*.pyx" src/sage | grep -v libs/symmetrica/all.py
```

The dynamic call sites will not appear in that grep; check
`combinat/sf/classical.py:55` and `combinat/sf/sfa.py:5656` for the two
`getattr` constructions by hand.
