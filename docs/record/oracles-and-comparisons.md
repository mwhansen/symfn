# Oracles and comparison harnesses

Sage, lrcalc and Symmetrica are all used as oracles and as baselines, and
the three differ in what they can prove. This file holds the harnesses
themselves and the findings that came out of building them — including
two whole-table deficits that only a like-for-like C comparison could
show.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

## Against Sage (`scripts/compare_sage.py`)

Sage is now a second comparison harness alongside lrcalc, covering the
operations lrcalc cannot. It differs in one respect and that respect is the
whole design: lrcalc runs **out of process**, one invocation per case, so both
sides pay startup and neither reuses a warm cache. Sage's interpreter costs
seconds, so both sides run **in one process** and startup is excluded by
construction — which forfeits the cache isolation a fresh process gave free, and
that has to be bought back explicitly:

* Sage memoizes, so every case uses *distinct inputs of comparable size*, each
  computed once. Repeating one input measures Sage's cache from the second call.
* symfn memoizes too, so `clear_caches()` (newly exposed through PyO3) runs
  before each timed call.
* Sage's first touch of a basis is lazy, so an untimed warm-up runs first or the
  first case absorbs setup belonging to all of them.

Every case is verified, not merely timed. **All cases agree with Sage at every
degree measured**, which independently cross-checks the rewrites above.

### ⚠️ Half these rows are not benchmarks against Python

Sage's conversions between the five classical bases **are Symmetrica**:
`sage.combinat.sf.classical.init` populates `conversion_functions` with
`t_<FROM>_<TO>_symmetrica`, confirmed at runtime. Everything else in the table
is Sage's own Python. The harness now prints a `via` column (`C` / `py`) so a
ratio is never read without its baseline, because the two mean very different
things: **3.7x on a `C` row is a stronger result than 9x on a `py` row.**

| case | via | deg 8 | deg 12 | deg 16 | deg 20 |
|---|---|---|---|---|---|
| p → s | C | 6.3x | 6.8x | 7.5x | 7.6x |
| s → m | C | 3.2x | 3.7x | 4.3x | 5.6x |
| **m → s** | C | 5.7x | 6.4x | 5.5x | **3.7x** *(was 0.004x)* |
| **s → e** | C | 6.6x | 8.1x | 4.9x | **3.7x** *(was 0.4x)* |
| s → h | C | 4.8x | 4.6x | 4.2x | 4.0x |
| **s → p** | C | 3.1x | 3.5x | 2.1x | **2.0x** *(was 1.2–1.6x)* |
| Kostka | py | 1.3x | 8.6x | 7.2x | 6.0x |
| plethysm | py | 9.0x | 8.4x | 8.5x | 9.0x |
| coproduct | py | 2.6x | 2.7x | 3.1x | 2.4x |
| skew | py | 7.1x | 6.4x | 4.5x | 4.5x |
| Hall | py | 2.5x | 1.8x | 1.9x | 1.7x |

Both former scaling failures are fixed (see the commit "Fix the two scaling
failures the Sage ladder found"); the ladder is what exposed them, since each
was *faster than Sage at degree 8* and only lost as degree climbed.

## Against Symmetrica directly (`scripts/compare_symmetrica.py`)

`compare_sage.py` can only see half of what it measures: Sage's five classical-
basis conversions dispatch into Symmetrica's C, but its *other* operations are
pure Python — even where Symmetrica ships a C implementation Sage never calls.
Plethysm read as **9x faster** than Sage while being **0.14x** against
`sage.libs.symmetrica.all.plethysm` on the same inputs. That defect sat inside a
green benchmark.

This harness drives Symmetrica directly through the bindings Sage installs but
mostly leaves unused, and verifies every case rather than only timing it — which
also makes it a second independent oracle for LR alongside lrcalc.

| case | deg 10 | deg 12 | deg 14 |
|---|---|---|---|
| LR product s_μ·s_ν | 3.5x | 3.8x | **4.7x** |
| skew s_{λ/μ} | 3.4x | 3.6x | 4.1x |
| Kostka K_{λμ} (one value) | 3.0x | 4.8x | 5.1x |
| character χ^λ(μ) (one value) | 7.9x | 0.6x | 0.5x |
| Hall ⟨s_λ, s_λ⟩ | 8.0x | 8.2x | 9.4x |
| **character table** | — | 3.6x *(was 0.6x)* | **2.9x** *(was 0.7x)* |
| **Kostka table** | — | 3.6x *(was 0.5x)* | **2.5x** *(was 0.39x)* |
| plethysm (single-row outer) | 1.5x | 1.7x | 2.1x |

At degree 16: character table **3.8x**, Kostka table **2.8x**, LR product 5.3x.

**LR is comfortably ahead of Symmetrica** — and the small-degree number badly
understates it. The ladder above tops out where a whole product is ~0.1 ms,
which says nothing about the regime this library targets. At real sizes:

| product | degree | Symmetrica | symfn | ratio |
|---|---|---|---|---|
| `[4,3,2,1]²` | 20 | 0.0026s | 0.0008s | 3.2x |
| `[6,5,4,3]²` | 36 | 0.259s | 0.0028s | **91x** |
| `[8,7,6,5]²` | 52 | 7.50s | 0.0170s | **442x** |
| `[10,8,6,4]²` | 56 | 530s | 0.0589s | **9002x** |

Every coefficient agrees, which also gives the LR engine a second independent
oracle alongside lrcalc — the first it has had. Symmetrica's
`outerproduct_schur` degrades violently over this range (0.0026s → 530s while
symfn goes 0.0008s → 0.059s), so `run_big_lr` stops as soon as it passes the
budget. At toy sizes this same comparison reads 3.2x, so sizing the benchmark
where the work actually lives changed the answer by three orders of magnitude.

**Two new deficits, both on whole tables.** Our *per-value*
Kostka and character are 3–5x faster, but the *whole table* is 2–2.5x slower and
the Kostka gap widens with degree (1.4x → 0.51x → 0.39x). That is the signature
of Symmetrica computing a table **as a table**, sharing work across entries,
while we answer p(n)² independent memoized queries. It is exactly the "produce
the whole answer in one sweep rather than query it entry by entry" pattern this
library has already applied to Kostka rows, the coproduct, p → s and m → s — and
has not applied to either table.

**Both are fixed, by the same observation: a table is p(n) *sweeps*, not p(n)²
numbers.**

* `kostka` bounds its chain DP by λ and reads one entry out of the final
  layer, discarding everything else that layer holds. Drop the bound and
  the layer at the end of μ's chain **is** the whole column — every λ with its
  K_{λμ} — for barely more than the single-value cost.
* `p_expand` already computed p_μ = Σ_λ χ^λ(μ) s_λ in one Murnaghan–Nakayama
  sweep, which *is* a column of the character table. The machinery was there; it
  had simply never been pointed at the table.

Columns then share with each other. K_{λμ} depends on μ only as a multiset, so
its parts can be consumed in any order, and characters likewise — taking them
descending lets partitions with a common prefix share the whole initial segment
of their chain. One traversal covers every μ. That is the same trie as
`convert::p_expand_shared`, which both now use.

Roughly **6–7x** on each table, turning both from losses into 2.5–3.8x wins.
Verified entry-by-entry against the per-pair functions through degree 12 (Kostka)
and 11 (characters) — a prefix-grouping slip would misattribute a column, and the
β-mask row index would drop rows rather than corrupt them, neither of which a
spot check catches.

⚠️ Single-value character rows above are unreliable and should not be read as a
deficit on their own: at 20–50 µs they are dominated by harness overhead, and
Symmetrica caches internally across calls while symfn calls `clear_caches()`
before each timed one. The *table* rows are the trustworthy comparison, which is
why they exist — a big symmetric unit of work with no caching asymmetry.

### The comparison above is on Symmetrica's home turf

Symmetrica's plethysm **refuses a multi-row outer partition** — it reports "for
the moment only for outer S_n" and computes nothing. Every case in the table is
therefore a single-row outer, the easy case, and possibly a specialized path.

symfn has no such restriction. `s_{21}[s_{21}]` (17 terms, 0.0011s),
`s_{22}[s_2]`, `s_{32}[s_{11}]`, `s_{21}[s_{31}]` (39 terms) and
`s_{31}[s_{22}]` (104 terms) all compute and all agree with Sage — 48/48 cases
across both shapes of input. So the honest summary is that symfn is faster than
Symmetrica on the inputs Symmetrica accepts, and is the only one of the two that
handles the rest.

### s → p: fixed, but not the way expected

`PowerSum::from_schur` calls `character_in(λ, μ)` once per μ, which looked like
the pattern removed from `p → s` — p(n) independent queries where one sweep
would do. Two measurements redirected the work:

1. **The obvious fix is a loss.** `p_expand` produces a *column* of the character
   table (one μ, all λ); `s → p` needs a *row*. Running the batched column sweep
   over every μ of the degree and reading off one row costs 0.0398s at degree 20
   against 0.0117s for the three shapes it would serve — **3.4x worse**, with
   break-even only around ten shapes.
2. **Characters were 84–100% of `s → p`**, so the transition arithmetic was never
   worth touching.

The real cost was representation, not algorithm. `border_strips` runs at *every
node* of the recursion and allocated a `Vec` for β, a heap **`HashSet`** for
membership, then per strip another `Vec` plus a sort — thousands of heap
allocations per character, to do arithmetic that fits in registers. Holding the
β-set in a u64 makes membership a bit test and the height a masked
`count_ones`. Worth **1.5x** at every degree (1.51x, 1.57x, 1.51x, over 6
interleaved rounds), taking `s → p` from 1.2–2.0x to **2.0–3.5x**.

Same β-number mathematics as before, and the general form is retained both as
the fallback past |λ| = 32 and as the reference the masked path is checked
against. That test earned its place immediately: asserting on strip *order*
failed, because the masked form walks β upward and the general one downward —
invisible to callers, which only sum, but noticed rather than assumed harmless.

Still open: `character_cached` clones both partitions into its key and takes a
global `RwLock` at every node, and `character_uncached` allocates a fresh
`Partition` for the μ-suffix at every node. Interning shapes to integer ids
would remove both.

⚠️ The plethysm row is capped at degree 10 regardless of the ladder setting
(its input is the *outer* partition and the result reaches degree 30), so that
row does not vary across the columns above — its four entries are the same
measurement repeated.
