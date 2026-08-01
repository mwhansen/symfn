# Transitions between the classical bases

Kostka numbers, characters, the coproduct, and the Jacobi–Trudi
conversions `s → h` / `s → e`, plus the sixth classical basis. Every
item here began as a benchmark finding: the operations were correct from
Phase 2 onward and exponential where they did not need to be.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

## Non-LR baselines (`examples/bench_ops.rs`)

Every non-LR operation was unmeasured; `scripts/bench_vs_sage.py` covers only
Schur products. The first run of the new harness found the hot spot immediately:

| operation | time | work |
|---|---|---|
| `character_table_n16` | 0.036s | 53 361 values |
| `coproduct_[7,6,5,4,3]` | 0.059s | 8 336 terms |
| `convert_m_to_s` (deg 12) | 0.254s | inverts the 77×77 Kostka matrix |
| `kostka_all_pairs_n12` | 0.264s | 5 929 values |
| `kostka_warm_50x` | 0.00027s | 3 850 cached lookups (66 ns each) |

**The β-sweep's row index was still on SipHash.** `character_table_in` builds
the whole table through one `p_expand` sweep whose frontier is already a
`fasthash::Map` on `u64` β-masks — but the `index` that maps a swept mask back
to a table row was missed, and the sink hits it once per emitted entry. Fixing
it is **1.12×**, on AC power, min-of-3 through the new
`character_beta_sweep_n{24,28}` rows of `bench_ops`:

| case | SipHash | `MixHasher` |
|---|---|---|
| `character_beta_sweep_n24` | 0.2088s | 0.1883s |
| `character_beta_sweep_n28` | 1.5993s | 1.4325s |

Those rows are a different engine from `character_table_n16`/`n18` above,
which loop `character(λ, μ)` pairwise over the memoized MN recursion; the name
collision was already misleading and the new rows are named apart from it.

The sweep of the rest of the tree that found this is in
[llt.md](llt.md) — the short version is that the win tracks whether the key is
already a word, and the `Vec<u32>`-keyed frontiers (`kostka_uncached` 1.08×,
`strip_lr` and `convert::jt_terms` nil) are not worth converting.

**Kostka is the bottleneck, and it is exponential.** Cost of one K_{λμ}:

| λ | degree | per value |
|---|---|---|
| `[3,2,1]` | 6 | 1.1 µs |
| `[4,3,2,1]` | 10 | 12 µs |
| `[5,4,3,2]` | 14 | 527 µs |
| `[5,4,3,2,1]` | 15 | 1.9 ms |
| `[6,5,4,3,2]` | 20 | **638 ms** |

At degree 20 that makes `convert_s_to_m` take **400 seconds** for a single
Schur function, since it needs K_{λμ} for all p(20) = 627 partitions μ. Compare
characters at 0.6 µs per value — Kostka is ~500x slower per number, and it is
the s ↔ m transition, so it also gates `convert_m_to_s` and anything routed
through the monomial basis.

**Fixed.** `kostka_uncached` enumerated SSYT one cell at a time. K_{λμ} instead
counts chains ∅ ⊆ λ¹ ⊆ … ⊆ λ with λⁱ/λⁱ⁻¹ a horizontal strip of size μᵢ, and
chains through the same intermediate shape *merge* — `strip_lr`'s DP without the
lattice condition. Pruning every intermediate to λ bounds the state space by the
partitions inside λ instead of by the tableaux of shape λ.

Measured by interleaved A/B of two `bench_ops` binaries, min of 2 rounds:

| operation | before | after | |
|---|---|---|---|
| `kostka_row_[5,4,3,2,1]` | 0.3246s | 0.0006s | **538x** |
| `kostka_row_[5,4,3,2]` | 0.0704s | 0.0003s | **216x** |
| `convert_s_to_m` | 0.0058s | 0.00013s | **46x** |
| `convert_m_to_s` | 0.2506s | 0.0075s | **34x** |
| `kostka_all_pairs_n12` | 0.2504s | 0.0080s | **31x** |

The speedup grows with degree, as exponential → polynomial should. On the case
that prompted this, `convert_s_to_m` for `s[6,5,4,3,2]` went from **400 s to
11.7 ms** (~34 000x), and degrees that were simply unreachable now run: degree
30 in 154 ms (2 317 terms), degree 40 in 2.70 s (16 306 terms).

Nothing else moved — characters, plethysm and the Hopf operations are unchanged
to within noise, which is the expected result since none routes through Kostka,
and is worth stating because a rewrite that quietly perturbed them would be a
regression hiding behind a headline.

Correctness: the Sage fixture covers only small degrees, so it cannot exercise
the new implementation where it now operates. `K_{λ,1ⁿ}` counts standard
tableaux, so the hook-length formula `n!/∏h(i,j)` checks it independently at
degrees up to 25.

**The coproduct had the same defect, one function away from its own fix.**
`skew_schur`'s doc comment records that it used to sweep all p(n) partitions
running a full LR backtrack per candidate, "the same answer for orders of
magnitude more work". `coproduct` was still doing exactly that: sweeping every
(μ, ν) pair of the right total degree and calling `NaiveLr::lr_coeff` on each.
Computing Δ(s_λ) = Σ_{μ⊆λ} s_μ ⊗ s_{λ/μ} instead gets every ν from one
traversal per μ — **0.057s → 0.0039s, ~15x**, identical 8 336 terms.

Worth noting as a pattern rather than a one-off: the fix already existed in the
same file, and the benchmark is what made the second instance visible. Nothing
about reading the code had surfaced it, in a codebase this well-commented,
because the comment explaining the mistake sat on the function that no longer
made it.

**After the Kostka rewrite the profile changed shape.** Rescaled to degree 20
(`examples/bench_ops.rs`), the leaders are now:

| operation | time | work |
|---|---|---|
| `convert_m_to_s` | **2.03s** | **1 term** |
| `kostka_all_pairs_n20` | 1.97s | 393 129 values |
| `convert_p_to_s` | 0.195s | 1 term |
| `hall_s_p_degree17` | 0.173s | 1 pairing |
| `coproduct_[8,7,6,5,4]` | 0.0092s | 19 758 terms |

`convert_m_to_s` builds the entire p(20)×p(20) Kostka matrix and inverts it to
use **one row**, and the matrix construction is ~97% of that (2.03s against
1.97s for the same number of Kostka values). The inversion itself is nearly
free by comparison.

A dominance early-out (`K_{λμ} ≠ 0` iff λ ⊵ μ) was added and is worth only
**1.1x** — I expected much more, on the reasoning that dominance is far sparser
than the lex triangle the matrix is built in. The chain DP's capacity pruning
was already rejecting those pairs quickly, so the test mostly replaces a fast
zero with a faster one.

**Fixed by computing only the row that is needed**, via the forward recurrence
`w[j] = 1`, `w[jj] = −Σ_{m=j}^{jj−1} w[m]·K[m][jj]`, cached per μ rather than per
degree. Measured **2.02s → 0.527s, 3.8x** — better than the ~2x predicted from
the triangle-versus-square argument alone, because skipping `m` where
`w[m] = 0` avoids the Kostka call entirely rather than computing a value and
multiplying it by zero. Inverse-Kostka rows are sparse enough for that to be
the larger half of the win.

**p → s without computing characters.** `p_μ = Σ_λ χ^λ(μ) s_λ` was implemented
by asking for χ^λ(μ) once per λ — p(n) independent recursions per term. But
Murnaghan–Nakayama is itself a multiplication rule, `p_k·s_λ = Σ (−1)^{ht} s_{λ∪ξ}`
over k-rim-hooks ξ, so multiplying successively by each part of μ builds the whole
expansion in ℓ(μ) passes and never evaluates a character.

| operation | characters | iterated MN | |
|---|---|---|---|
| `convert_p_to_s` | 0.203s | 0.0997s | **2.0x** |
| `hall_s_p_degree17` | 0.179s | 0.0991s | **1.8x** |
| `plethysm_[3,2][[2,1]]` | 0.00464s | 0.00188s | **2.5x** |
| `plethysm_[4][[3]]` | 0.00266s | 0.00109s | **2.4x** |

⚠️ **The first implementation of this measured 0.8x — a regression — and the
algorithm was not what changed.** Keying the frontier on `Vec<i64>` β-numbers
meant hashing a 160-byte key per rim hook, and the character path it was
competing against has a memo cache that shares subproblems across every term.
β values here are below 2n, so for n ≤ 32 the whole set is a `u64` bitmask:
adding a rim hook becomes two shifts and a popcount, and the key is one word.
Same algorithm, 2.5x swing. Degrees past 32 fall back to characters, which stay
exact in `C` for bignum rings.

That is the second time this session a data structure inverted an algorithmic
conclusion — the other was `HashMap` versus generation-stamped dense tables in
`three_row`, also worth 3.3x. Both times the operation counts were unchanged and
only the constants moved. **A structural idea that benchmarks badly on its first
implementation has not been tested yet.**

## The forgotten basis: Macdonald's sixth

`f_λ = ω(m_λ)` (Macdonald I.2). This is the last of the six classical bases and
the only one with no independent combinatorial description — it is *defined* as
the image of the monomial basis under ω, which is where the name comes from.
Adding it makes `convert` total over the standard set.

The implementation is two lines each way, because ω is an involution:

```text
  to_schur:    Σ c_λ f_λ = ω(Σ c_λ m_λ)         →  ω(monomial_to_schur(c))
  from_schur:  x = Σ c_λ f_λ ⟺ ω(x) = Σ c_λ m_λ →  monomial coeffs of ω(x)
```

and ω on the Schur basis is conjugation of every index. So both directions
inherit Muir's rule and the Kostka machinery — and their asymptotics — for one
transpose per term.

Being nearly free makes the *testing* the real work, since a test phrased in
terms of ω would only restate the implementation. Three checks that don't:

* **Endpoints, hand-derived.** `f_(n) = (−1)^{n−1} p_n` and `f_(1^n) = h_n`,
  from `m_(n) = p_n` and `m_(1^n) = e_n` pushed through ω. Neither mentions the
  forgotten basis; both are facts about the other bases.
* **Duality**, the structural characterisation: `⟨f_λ, e_μ⟩ = δ_{λμ}` for every
  pair through degree 7, reached through `hall` and the e → s expansion — code
  the conversion never touches. Wiring `Forgotten` to the wrong involution, or
  to conjugation of the *index* rather than of the Schur expansion, would fail
  this while a round trip still closed.
* **Round trips** out through m, e, and h and back, all partitions to degree 7.

Both directions agree with Sage at degrees 8 and 12.

⚠️ The measured ratios are 146x/989x (s → f) and 379x/888x (f → s), and they
should not be quoted. Symmetrica has **no forgotten basis** — verified, not
assumed: `sage.combinat.sf.classical.conversion_functions` holds exactly 20
entries, the ordered pairs among {Schur, elementary, homogeneous, monomial,
powersum}, and forgotten appears in none. So Sage falls back to a generic
Python basis-change through its own machinery, and this is a `py` row against an
unoptimised path. It says the basis is not a bottleneck; it says nothing more.

## s → h and s → e: a 200x regression, then a 30x win — from reading Symmetrica

**The first thing letting Sage pick the inputs found.** `s_(14) → e` took
**1.54s** against Symmetrica's 0.0094s. Not a constant factor: the cost is
exponential in the Jacobi–Trudi matrix size, which is ℓ(λ) for s → h and **λ₁**
for s → e, and `jt_terms` enumerates permutations of it.

`convert.rs` had *named* wide shapes as the hazard in its own docstring. The
ladder still reported 3.5–7.6x ahead at every degree, because `shapes_of()`
built only balanced shapes of 3–6 rows, so λ₁ never exceeded about 7. **Those
earlier s → e and s → h ladder numbers measured one shape family, not the
conversion.** `shapes_of` now emits the single row and the hook first.

### What Symmetrica actually does

`Symmetrica_2.0/tse.c` conjugates and calls `tsh_jt` + `tsh_eval_jt` — **the
same algorithm as ours**: enumerate permutations of the Jacobi–Trudi matrix,
skip the entries that vanish, sign by the permutation, collect the index
multiset. No better formula, no special-casing. The entire gap was one line of
representation:

| | matrix entry | vanishes when | tight row |
|---|---|---|---|
| Symmetrica `tsh_jt` | λ_i + i − j | j > λ_i + i | **row 0** |
| symfn, before | c_i − i + j | j < i − c_i | last row |

The two are transposes and describe the same determinant. But c is weakly
decreasing, so `i − c_i` *increases* with i: our constraint tightened as the row
index grew, meaning the walk began at the least constrained row — row 0 accepted
any column — and only met the dead ends near the leaves, long after the
branching had happened. Symmetrica's orientation puts its tight row first for
free.

**Assigning rows last-to-first instead of first-to-last is the whole fix.** For
the degree-24 hook that alone is 0.070s → 0.00105s; whole subtrees now die at
depth 1 rather than at depth 12.

### The other two changes, one of which was wrong

* **Take the smaller matrix and flip h ↔ e** (via Newton's identity as a linear
  recursion, symmetric under the swap so one routine serves both ways). Correct,
  but ⚠️ **applying it whenever the other matrix is smaller made things worse**:
  over every partition of degree 20 it gave 0.33x, where never flipping gave
  2.20x. Flipping is not free, and for a shape like (5,5,5,5) — matrices of 4
  and 5 — a 5-wide determinant costs far less than expanding a degree-20
  h-element into e. It is now reserved for matrices of 14 or more, where the
  determinant really is exponential and the conjugate collapses it: a single row
  of degree 24 is 1.66s direct and 0.004s flipped, its h-expansion being the one
  term h_24.
* **A Muir sweep as backstop**, using
  `coefficient of h_μ in s_λ = coefficient of s_λ in m_μ`. ⚠️ **Nearly useless**:
  for the degree-20 hook the determinant takes 0.0048s against the sweep's
  0.37s. p(n) Muir expansions cost more than a wide determinant. Kept only to
  guard shapes where the determinant would genuinely explode.

### Where it lands

Over **all** partitions of a degree, not a sample:

| degree | s → e | s → h |
|---|---|---|
| 12 | 10.41x | 9.69x |
| 14 | 4.85x | 5.30x |
| 16 | 4.82x | 5.87x |
| 18 | 5.21x | 5.74x |
| 20 | 4.79x | 5.53x |

From 0.33–0.87x to 4.8–10.4x, on a representation change plus a threshold. The
4678 Sage-driven computations still agree.

## Profiling all ten transitions at once (`examples/profile_convert.rs`)

`bench_ops` times every transition on **one** shape, `[6,5,4,3,2]`. That is the
input on which the two Jacobi–Trudi directions look alike, and it is how the
200x `s → e` regression above survived a whole degree ladder. The new harness
sweeps all ten against five shape families chosen as the corners of the
ℓ(λ)-against-λ₁ trade — row, column, hook, staircase, rectangle — and has a
second mode that repeats one case to a deadline for `sample`:

    cargo build --release --example profile_convert
    ./target/release/examples/profile_convert sweep 20

    cargo build --profile profiling --example profile_convert
    ./target/profiling/examples/profile_convert loop m2s 20 12 &
    sample $! 9 -mayDie -f /tmp/convert.txt

All ten rows run over ℚ, deliberately: the table compares transitions against
each other, and mixing ℤ and ℚ rows would attribute the ring's cost to the
transition. Read the integral directions against `bench_ops`, not across this
table.

**The forward directions are not where the time is.** At degree 20 every
`s → X` is under 2 ms on every family; the whole cost sits in the reverse
directions, and the h ↔ e flip is doing its job — `s → h`/`h → s` and
`s → e`/`e → s` are clean mirror images across the row and the column, which is
what the fix above predicted and nothing had checked since.

### `muir_expand` was accumulating on a `Partition`, not on the β-mask

Sampling `m → s` on `s_(20)` (0.386 s/iteration, against the sweep row's
0.362 s — the check that the profile is of the same work), self time out of
8269 samples: `muir_rec` 5409, `Partition::new` **765**, SipHash's
`DefaultHasher::write` **425**, malloc/free ~950.

`muir_rec` reached its leaf holding the β-mask — the recursion's own state, one
word — and then built, heap-allocated and SipHashed a fresh `Vec<u32>` from it
to key a `HashMap<Partition, i128>`, once per leaf, to produce a few hundred
distinct terms. `p_expand_shared`'s caller had had exactly this fix for a while
("Accumulate on the β-mask, not on a Partition"), a thousand lines up the same
file; the sibling never got it. Keying on the `u64` and calling
`mask_to_partition` once per surviving term instead:

| case | before | after | |
|---|---|---|---|
| `m → s`, `s_(20)` | 0.3668s | 0.1703s | **2.15x** |
| `m → s`, hook `[10,1^10]` | 0.0793s | 0.0300s | **2.64x** |
| `m → s`, rectangle `[4^5]` | 0.1080s | 0.0413s | **2.62x** |
| `m → s`, staircase | 0.1062s | 0.0405s | **2.62x** |
| `f → s`, column `[1^20]` | 0.3644s | 0.1703s | **2.14x** |
| `bench_ops convert_m_to_s` | 0.2457s | 0.1044s | **2.35x** |

Interleaved A/B of two binaries, min of 4 rounds (3 for `bench_ops`), **on AC
power**. Every other row of both harnesses is flat — 0.93–1.04x across the
whole sweep, and `kostka_all_pairs_n20`, the character sweeps and the coproduct
all within 1% — which is the result to want, since `muir_expand` should touch
nothing but the monomial directions. `f → s` moves exactly as `m → s` does
because it *is* `m → s` plus a transpose.

⚠️ **The prediction from the self-time column was 1.33x and the fix delivered
2.15–2.64x.** Reading `Partition::new` + SipHash as the recoverable share
undercounted it: the allocator traffic those two generate was sitting in
`malloc`/`free`/`memmove` frames that a flat profile attributes to libsystem,
not to the line that caused them. The confirming profile is the cleaner
evidence than the arithmetic — after the change `muir_rec` is 6781 of ~7100
samples (95%), `Partition::new` and `DefaultHasher::write` have left the
top-of-stack list entirely, malloc/free is down about 10x, and
`mask_to_partition` costs 89 samples (1.2%) doing the conversion once per term.

Correctness rests on `muir_agrees_with_the_inverse_kostka_solve`, which checks
every μ through degree 12 against the linear solve — an independent route to
the same row of K⁻¹, and the implementation this replaced — plus
`muir_matches_hand_computation`.

### h → s and e → s: Pieri instead of the LR engine

The open tail above proposed batching `Elementary::to_schur` and
`Homogeneous::to_schur` the way `p_expand_shared` batches p → s. Measuring
before writing it said the batching was the *small* half, and that turned out
to be right — but it also pointed at the large half.

Both build h_λ and e_λ as a product of one-row or one-column Schur functions,
and `Schur::mul` routes to `AutoLr`, which built and expanded a **skew shape**
for each factor. Multiplying by s_{(k)} or s_{(1^k)} is Pieri — add a
horizontal or a vertical k-strip — a direct enumeration with no LR machinery
under it. The frontier was also rebuilt from `Schur::unit()` per term, and
`out.add(&prod.scale(c))` allocated a scaled copy plus a merged map per input
term; `terms()` is a `BTreeMap` keyed by `Partition`, which orders
lexicographically by parts, so shared prefixes are *already contiguous* and the
traversal needs no sort (unlike `p_expand_shared`, which is handed a `Vec`).

Interleaved A/B of three binaries, min of 4 rounds, on AC power. The two
columns separate the changes: Pieri alone, then prefix sharing on top of it.

| case | before | Pieri | + sharing | Pieri | sharing | total |
|---|---|---|---|---|---|---|
| `e → s`, `s_(20)` | 0.2197 | 0.1918 | 0.1822 | 1.15x | 1.05x | **1.21x** |
| `e → s`, hook `[10,1^10]` | 0.0049 | 0.0024 | 0.0023 | 1.99x | 1.08x | **2.14x** |
| `e → s`, staircase | 0.0024 | 0.0009 | 0.0009 | 2.72x | 1.03x | **2.80x** |
| `e → s`, rectangle `[4^5]` | 0.0028 | 0.0012 | 0.0011 | 2.40x | 1.02x | **2.46x** |
| `h → s`, column `[1^20]` | 0.2135 | 0.1931 | 0.1849 | 1.11x | 1.04x | **1.15x** |
| `h → s`, hook | 0.0063 | 0.0041 | 0.0039 | 1.53x | 1.06x | **1.63x** |
| `h → s`, rectangle | 0.0074 | 0.0050 | 0.0048 | 1.49x | 1.04x | **1.55x** |

Every other row of the sweep is flat, 0.96–1.05x.

⚠️ **The prefix sharing is worth 2–8%, not the 16x its p → s analogue carries.**
That was predicted before it was implemented, by counting rather than timing:
over the partitions of 20 the sharing removes 1.71x of the Pieri *steps* but
only 1.30x weighted by the degree of the element each step multiplies into.
What two terms share is a short cheap prefix; the leaves, which are the
expensive steps, are shared by nothing. The 1.30x is an upper bound and the
measured 1.02–1.08x sits under it because copying the frontier is not free
either. Kept because it is nearly free once the traversal is written this way,
not because it carries the win. p → s gets 16x from the same shape because its
step is a rim hook on a `u64` mask, not because sharing is inherently worth
more there.

⚠️ **A methodology note, since it nearly produced a false result.** The first
attempt to separate the two columns built the Pieri-only variant from a source
edit that *failed to compile*; the `cp` that followed copied the previous
binary, so the "sharing" column was the same binary measured twice and read as
a clean 0.93–1.04x null. Interleaving and min-of-N do nothing about this class
of error. Checking that the two binaries differ (`md5`) before believing an A/B
is now the habit.

### Open tail

* **The big cases are ~73% data structure, not mathematics.** Re-profiling
  `e → s` on `s_(20)` after the above, the self time splits: allocator
  (malloc/free/memset) **53.8%**, `memmove` **13.3%**, strip enumeration only
  **11.3%**, `Partition::new` 6.4%, BTreeMap ops 2.1%, unattributed 13.2%. The
  frontier is a `Schur<C>` — a `BTreeMap<Partition, C>` — so every shape a
  Pieri step emits allocates a fresh heap `Vec<u32>`, sorts it in
  `Partition::new` (already weakly decreasing, so the sort is pure waste), and
  is inserted into a B-tree that memmoves on the way in. That is why the two
  biggest cases take 1.11–1.15x from Pieri where the mid-size ones take
  2.4–2.7x: once the frontier holds thousands of terms, emitting and inserting
  them dominates whatever the multiply costs. It is also the same defect, in
  the same file, that the β-mask keying fixed for Muir — degrees here are
  inside `MASK_LIMIT`, so a `Map<u64, C>` frontier with partitions built once
  at the end is the shape of the fix. It is *not* a transcription of the Muir
  change: rim hooks move one β value, whereas a horizontal or vertical strip is
  an interlacing condition on the whole β-set, and getting that wrong yields
  plausible wrong answers rather than errors. Pin it with a doctest against the
  Pieri rule before trusting it.
* **`Partition::new` sorts shapes that are already sorted.** Cheap and
  self-contained: the strip enumerators, `mask_to_partition`, and the Muir leaf
  all produce weakly decreasing parts by construction. A checked
  `from_sorted_desc` constructor would drop the sort from every one of them.
  Unmeasured on its own — it is 6.4% of `e → s` *including* the allocation the
  sort does not cause, so the prize is smaller than that number looks.
