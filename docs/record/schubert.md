# Schubert polynomials

The last subsystem Sage still routed through Symmetrica, so retiring it
completes the displacement the README opens with. Three engines, a ranking
that inverted twice, and a single-coefficient query that answers products no
machine can hold.

Built 2026-07-30 (`src/permutation.rs`, `src/schubert.rs`). Split out of
[the record index](README.md), which carries the phase plan and a summary of
this file.

---

**Three engines, and the ranking inverted twice.** E1 is Symmetrica's route
(expand a factor to monomials, Monk chain per monomial), kept forever as
the oracle; E3 walks the same peel recursion as a merged DAG; E2 is the
memoized Lascoux–Schützenberger transition. A pre-implementation
measurement of *state compression* — pipe dreams ÷ peel-DAG states, 687×
on `stair7`, 125 599× on an S₁₃ element — pointed hard at E3, so E3 was
built first and E2 demoted. E3 duly beat E1 by 11.2× on `stair6²`… and lost
to the C `schubmult` by 4–51×. E2 then beat E3 by up to 97×.

**The compression metric was the mistake.** Both engines cost (nodes) × (size
of the running element); the metric counted only nodes. On `S_11.1` E2 uses
*more* nodes than E3 (1 102 vs 571) and is 97× faster, because E1/E3 expand a
factor into monomials so the running element inflates to answer-size early,
while transition never expands. The model ranked E3 above E2 on the one term
it measured.

**Against the incumbents** (single-threaded; schubmult C re-measured the same
day and reproducible to 5%; its rows carry ~3 ms of process startup, so rows
≤0.03s measure `exec` and are excluded — floor measurement below):

| case | terms | Sage/Symmetrica | schubmult C | symfn E2 | vs C |
|---|---|---|---|---|---|
| stair6² | 247 | >120s | 0.355s | 0.070s | **5.0x** |
| stair7² | 1 111 | >120s | 20.6s | 1.51s | **13.6x** |
| S₁₁.1 ℓ=35,26 | 118 822 | >120s | 0.748s | 0.392s | 1.9x |
| S₁₂.2 ℓ=43,32 | 242 507 | >120s | 19.0s | 0.634s | **30.0x** |
| S₁₃.2 ℓ=40,41 | 3 241 903 | >120s | 105.8s | 5.10s | **20.8x** |

Ahead on every row either incumbent finishes. Sage's own engine finishes none
of these but `stair5²`, where it takes 6.85s against our 0.035s.

**Where it stops: S₁₅.** S₁₄ is routine (2.6s/53s/83s for 0.37M/7.1M/6.7M
terms); the first S₁₅ pair takes 421.9s for 12.4M terms. Peak RSS 1.52 GB on
the 3.2M-term S₁₃ row.

**Coefficient growth: the "≤16" figure was an artefact of the incumbent.** It
came from the products Symmetrica could finish. Measured on our own engine:
**130** at S₁₃, **591** at S₁₄, **863** at S₁₅. Still far from `i64`, but
growing, and `guard` earns its keep.

**The part no engine tuning could deliver: `schubert_coeff`.** E2 with Bruhat
pruning answers structure constants for products that **cannot be
materialized** — `S_13 ℓ=25,36` has a monomial mass of 4.3×10¹⁶ — in ~0.04s
each, 578× faster than reading the coefficient off a whole product where one
exists. Values up to `c = 18` exhibited, checked 0-mismatch against a
computable neighbor (below). No other package has such a query.

**Memory was a bug, not a trade.** The first E2 retained every node's product
and died at 6.56 GB. A parent-count pre-pass over the transition tree (which
builds no elements) plus drop-on-last-use fixed it with **no time regression
on any row** (full account below).

Correctness: 355 lib tests, engines cross-checked exhaustively on S₅ × S₅
(E1 ≡ E2 ≡ E3, both argument orders), `∂` against `∂` by exact division on
expanded polynomials, the `w₀` recursion, Bruhat order laws, and `stanley`
against a hand-computed value and Grassmannian vexillarity. Bindings
separately: **427 checks against Sage, 0 failures**, including `newtrans`,
against which `stanley` is now **1.3× at S₁₆, 4.4× at S₁₄** — from removing
heap allocations, not memoization, which was measured and never written
(below).

## Measured walls: Sage, schubmult C, and three traps along the way

SageMath 10.9, this machine, 2026-07-29. Single runs, per-case alarm of 120s.
⚠️ **On battery power (68%, discharging)** — these are order-of-magnitude
walls, not benchmarks; ratios between the two tables below are safe (adjacent
runs, same conditions), absolute times are not.

**Products through Sage** (`SchubertPolynomialRing(ZZ)`, i.e. Symmetrica):

```text
  case                                   sec      terms   max|coeff|
  stair4² (S_8)                        0.078         16       4
  stair5² (S_10)                       6.85          59      16
  stair6² (S_12)                       >120           —       —
  rand S_10  ℓ=21,22                   47.2        1070       5
  rand S_10  ℓ=19,18                   0.0004        12       1
  rand S_10  ℓ=28,20                   0.0067        25       1
  rand S_11  ℓ=32,27 / 35,26 / 31,21   >120 (all three)
  rand S_12  (three pairs)             >120 (all three)
```

`stairk` is the Grassmannian permutation [2,4,…,2k,1,3,…,2k−1], Schubert
polynomial `s_{(k,k−1,…,1)}(x₁..x_k)` — its stable shadow is the staircase
LR family this crate measures everything on.

**Same cases, lrcalc's C `schubmult`** (out of process, min-of-3 below
0.5s; ~3 ms startup per row, measured exactly below). Re-measured
2026-07-30: 10 of 11 rows within 5% (the outlier, `S_12 ℓ=27,21`, went
1.61 → 2.25s, inside the LR record's ±30% noise
(`docs/record/littlewood-richardson.md`)) — licensing this table as a
baseline, and confirming the machine matched state for the Sage rows above:

```text
  case                                   sec      terms
  stair4² (S_8)                        0.005         16
  stair5² (S_10)                       0.016         59
  stair6² (S_12)                       0.355        247
  stair7² (S_14)                      20.6         1111
  rand S_10  ℓ=21,22                   0.066       1070
  rand S_11  ℓ=32,27                   0.222      30143
  rand S_11  ℓ=35,26                   0.748     118822
  rand S_11  ℓ=31,21                   3.86      185284
  rand S_12  ℓ=37,24                  21.0       230901
  rand S_12  ℓ=27,21                   1.61        8280
  rand S_12  ℓ=43,32                  19.0       242507
  rand S_13  ℓ=25,36                  >120            —
  rand S_13  ℓ=36,41                  20.9       947254
  rand S_13  ℓ=40,41                 105.8      3241903
```

Three methodology traps came out of building this ladder, recorded so they
are not hit again:

- **Rectangle Grassmannians are dominant.** The first ladder used
  `w = [k+1,…,2k,1,…,k]` as a "hard LR case", but `S_w = s_{(kᵏ)}(x₁..x_k)`
  in exactly `k` variables is the single monomial `(x₁⋯x_k)ᵏ`, so every
  product collapsed to one term and measured nothing. Check the *flagged*
  object is non-trivial; do not import intuition from the stable limit.
- **The peel memo key must carry the level, not just `(perm, stufe)`.**
  `level = n − len(p) + 1` pins the level from the permutation length — but
  only *within one top-level call*, because `n` differs between
  permutations. A memo shared across `S_{132}` (n = 3) and `S_{1423}` (n = 4)
  returns a length-3 subtree computed in `x₁, x₂` to a caller expecting
  `x₂, x₃`. **The bug is invisible to the obvious test**: the
  single-permutation `expand → from_polynomial` round trip passed
  exhaustively through S₆, because each call built its own memo — only
  `from_polynomial` on a *sum of permutations of different sizes* shares one
  memo across differing `n`, and it is the only test that failed. Any engine
  that memoizes the peel inherits this, E3 included.
- **Two variable-index conventions circulate inside one C file.**
  Sage's `multiply_variable(i)` and Symmetrica's `mult_schubert_variable`
  are 0-based; Symmetrica's own `divdiff_schubert` is **1-based** — `sb.c`
  disagrees with itself, and the first transition verification **failed
  100% of cases** on exactly this. The crate fixes one convention (1-based,
  matching the mathematics) and tests it at the boundary.

## What Symmetrica's C actually does

Read at function level from `sb.c` / `mss.c` / `perm.c` (public domain: a
legitimate algorithm source, not merely an oracle). Adoptable pieces ✓,
defects ✗:

- **Expansion** (`m_perm_schubert_monom_summe`, `sb.c:123`; worker
  `algorithmus2`, `sb.c:323`): iterated inverse Monk at position 1, peeling
  x₁'s exponent and recursing on covers `w ⋖ w·t_{1,i}` via a
  running-minimum scan (✓ the right cover enumeration, reused in E1/E3). One
  leaf per pipe dream, cost `S_w(1,…,1)`. ✗ no sharing between leaves.
- **Product** (`mult_schubert_schubert`, `sb.c:965`): expand one factor to
  monomials, then run `|α|` single-variable Monk passes per monomial over
  the accumulating sum (`mult_schubert_variable`, `sb.c:1003` — ✓ its
  two-scan signed cover enumeration is exactly the verified signed Monk
  rule). ✗ cost = (#pipe dreams) × (Monk chains re-run per monomial), the
  two blowups behind the S₁₀ wall above. ✗ "expand the smaller factor"
  compares only the head term's vector length.
- **From-polynomial** (`t_POLYNOM_SCHUBERT`, `sb.c:232`): greedy peel on the
  lex-smallest monomial, pad to a valid Lehmer code, subtract `c·S_w`,
  repeat (✓ triangularity verified; adopted). ✗ re-expands `S_w` from
  scratch every iteration, no caching.
- **Divided differences** (`divdiff_schubert`, `sb.c:1334`):
  `∂ᵢS_w = S_{wsᵢ}` or drop — ✓ trivially cheap, adopted. ✗ 1-based,
  disagreeing with the 0-based product routine (the third trap above); ✗
  latent out-of-range read on a stored permutation shorter than the letter
  applied.
- **Pairing** (`scalarproduct_schubert`, `sb.c:1840`): the Poincaré pairing
  as full product then `n(n−1)/2` divided-difference passes, with **n read
  off the stored vector lengths** — the answer depends on how padded the
  inputs happen to be. ✗ semantics a caller cannot predict.
- **Stanley/Schur** (`newtrans`, `mss.c:44` — the one well-engineered
  routine): the Lascoux–Schützenberger transition, iterative with an
  explicit stack, a Grassmannian base case emitting a single `s_λ`, and the
  `F_w = F_{1×w}` shift when the left sum is empty (✓ adopted intact). ✗
  static `char[1000]` state, not reentrant. ⚠️ **Its Sage name
  `t_SCHUBERT_SCHUR` oversells it**: it computes the Stanley symmetric
  function `F_w`, equal to `S_w` only in the stable range —
  `newtrans([2,1,4,3]) = s₂ + s₁₁` (verified) while `S₂₁₄₃` is not even
  symmetric.
- ✗ **Nothing is memoized anywhere in the module**, every routine mutates
  its inputs in place, and the identity permutation is `[1,2]` — stability
  is patched in the comparator rather than by a normal form.

The bill comes due at interpreter exit: after a product sweep, Sage's
process printed Symmetrica's `ERROR: permutation memory not freed?` banner
with **`mem_counter_perm = 215164`** live objects, in a library embedded in
Sage for twenty years — the same argument [lib.rs](../../src/lib.rs) already
makes against `OP`.

## Three engines for the product

E1 (Symmetrica's route, kept forever as the oracle), E3 (the peel recursion
walked as a merged DAG) and E2 (the memoized transition) were the three
candidates. What follows is the measurement history that picked among them.

### E3: compression measured before the code existed

Before any Rust existed, the peel recursion's state-compression ratio — pipe
dreams ÷ merged-DAG states — was measured against Symmetrica's
`algorithmus2` (`sb.c:323`): a state is `(perm, alphabetindex, stufe)`, and
merging on `perm` alone is valid because `len(perm)` determines
`alphabetindex` (only one transition changes either) and `stufe` enters only
as an overall power-of-one-variable factor, so a cached state pays a
Monk-style shift rather than a re-walk. Validated against Sage's `expand()`
exhaustively on S₂–S₆, 872 permutations, 0 mismatches:

```text
  case                     ℓ     pipe dreams      states  edges  peak live  compress
  stair4  (S_8)           10              64          88    113        12      0.7×
  stair5  (S_10)          15           1 024         283    382        28      3.6×
  stair6  (S_12)          21          32 768         923  1 286        69     35.5×
  stair7  (S_14)          28       2 097 152       3 052  4 340       182       687×
  stair8  (S_16)          36     268 435 456      10 192 14 693       499    26 338×
  rand S_11.0u            32         142 818         521    744        27       274×
  rand S_12.1u            27      15 468 012       5 092  7 706       279     3 038×
  rand S_13.0u            25      31 877 880       2 740  4 094       133    11 634×
  rand S_13.0v            36   1 355 962 209      10 796 17 006       546   125 599×
  w0(S_12), dominant      66               1          11     10         2       0.1×
```

`compress = pipe dreams ÷ states`, the direct analogue of the LR record's
(`docs/record/littlewood-richardson.md`) `LR tableaux ÷ states produced`.
Four readings: **the leaf count is exponential and the state count is not**
(staircase pipe dreams are exactly `2^C(k,2)` while states grow ~3.3× per
rung — the gap is the whole engine bet); **the compression is largest
exactly where the incumbents die** (the `ℓ=25,36` pair neither incumbent
finishes, above, peels its smaller factor at 31 877 880 pipe dreams against
2 740 states); **memory is not the trade this time** (`peak live` —
simultaneously-live state values under DFS post-order with refcounting,
each a whole Schubert element — never exceeds 546, even on the
1.4-billion-pipe-dream row, where the LR work bought speed with memory and
said so); and **compression below 1× exists and is bounded** (dominant and
short permutations have more states than pipe dreams, but counts stay ≤ 68
— a small constant, not a regime).

⚠️ Still a node count, not a runtime: per-state work is element-sized, so
the ratio is *indicative* of the speedup and not equal to it.

### E3 built: 11.2× over E1, and what the compression ratio overstated

E3 landed keyed on `(perm, level, stufe)`:

```text
  case        compression   E1        E3        E3 vs E1   schubmult C   vs C
  stair4²           0.7×    0.0047s   0.0049s      1.0×       0.005s      1.1×
  stair5²           3.6×    0.1767s   0.0674s      2.6×       0.016s      0.24×
  stair6²          35.5×   27.7062s   2.4751s     11.2×       0.355s      0.14×
  stair7²           687×    (~8h est) 141.54s        —        20.6s       0.15×
  S_11.1 ℓ=35,26     50×          —    38.05s        —        0.748s      0.02×
```

**The engine bet paid off against E1, at roughly a third of the nominal
compression** (35.5× → 11.2×): the `(perm, level, stufe)` key merges less
than the `perm`-only key the compression table measures, and per-state work
is element-sized where E1's per-leaf work is a monomial chain — the
calibration for the next time a state-count ratio is used to predict a
speedup.

The displacement bar was met — `stair6²`/`stair7²` complete where Sage's
engine completes neither, and `stair5²` went 6.85s → 0.067s, 102×, where
both run. The parity bar (matching schubmult C) was missed two ways.
On the staircase family it was a flat ~7×: a profile of `stair7²` (15 660
samples) put **~47% of the time in `BTreeMap` and the allocator**
(malloc/free 25%, insert/remove 15%, memmove/bzero 7%) — `Schubert` stores
`BTreeMap<Perm, C>` with `Perm` a heap `Vec`, so every Bruhat cover
allocates: `three_row`'s lesson again, a sound idea on the wrong data
structure. On `S_11.1` the gap was 51×, and that was **not** a
data-structure story: with 118 822 output terms, E3's cost is (Monk passes)
× (element size), and the element is answer-sized for most of the DAG;
schubmult produces those terms in 0.748s while E3 does ~10⁹ term updates to
get there.

**Tuning E3, same day: one win, one instructive null result.**

```text
                                     stair5²   stair6²   stair7²   vs schubmult C
  E3 as first written                 0.0674s   2.4751s   141.54s      0.15×
  + inline Perm ([u8; 32], Copy)      0.0345s   1.3862s    81.12s      0.25×
  + in-place add + Rc-keyed memo      0.0565s   1.3540s        —       0.26×
  schubmult C                         0.016s    0.355s     20.6s        —
```

**Inline `Perm` was worth ~1.8×**: with `[u8; 32]` plus length, `Copy` and
zero-filled so `Ord` is unchanged, the allocator fell from 25% to ~8% of
samples and `BTreeMap` insert/remove from 15% to ~4% — below the 3–5× the
profile predicted, since removing the *allocations* does not remove the
*inserts*.

**In-place `add_assign` plus an `Rc`-keyed memo bought ~2%, i.e. nothing.**
The hypothesis — that copying the accumulator per cover, and cloning a
whole element per memo hit, were the residue after the profile's
`clone_subtree` line — was wrong: with both gone the time stayed put,
locating the cost in the *number* of insertions, not their cost. Kept as
strictly better code, but a dead end no one should re-derive expecting a
win.

### E2: the memoized transition, and it is the engine

Implemented as `Schubert::mul_e2`, verified against E1 exhaustively on
S₅ × S₅ in both argument orders (E2 recurses on the *second* factor, so the
orders walk different trees), plus staircases.

**`Schubert::mul` now delegates to it** — and for a while it did not, which is
worth recording because the failure was invisible from inside either half.
`mul` kept calling E3 after E2 superseded it, while `python.rs` called
`mul_e2` directly; so a Rust caller and a Sage caller got engines differing by
up to 97×, and `examples/bench_schubert.rs` — which called `mul` — reported the
*superseded* engine's times against the C `schubmult`, making the crate look
4–51× behind when it is 1.9–27.6× ahead. Nothing was wrong in any one file.

Two changes keep it from recurring: the binding now goes through `mul` rather
than naming an engine, so "best engine" has exactly one definition; and the
three `e1_and_e3_agree_*` tests call `mul_e3` explicitly, so pointing `mul`
somewhere else cannot silently orphan E3's coverage. E3 is kept and kept
tested — it is the historical comparison the compression-metric lesson rests
on.

```text
  case              terms      E3 nodes  E2 nodes  E2 passes   E3 time    E2 time   schubmult C   E2 vs C
  stair5²              59           641       191        333    0.0345s   0.0066s     0.016s        2.4×
  stair6²             247         2 559       633        931    1.3540s   0.0665s     0.355s        5.3×
  stair7²           1 111        10 052     2 085      2 667   81.1190s   1.4821s     20.6s        13.9×
  S_10.0            1 070           307       152        305          —   0.0125s     0.066s        5.3×
  S_11.0           30 143           181        91        204          —   0.0848s     0.222s        2.6×
  S_11.1          118 822           571     1 102      5 477   38.0513s   0.3923s     0.748s        1.9×
  S_11.2          185 284           265     1 156      2 064          —   0.7605s     3.86s         5.1×
  S_12.0          230 901         7 067     1 862      3 154          —   6.4134s    21.0s          3.3×
  S_12.1            8 280         1 787       556      1 117          —   0.6541s     1.61s         2.5×
  S_12.2          242 507           267     1 126      1 718          —   0.6221s    19.0s         30.5×
  S_13.1          947 254           734       351      1 368          —   1.8852s    20.9s         11.1×
  S_13.2        3 241 903         2 288       447        943          —   5.0792s   105.8s         20.8×
```

(Pre-eviction-fix numbers; the fix below moves times by noise only —
`stair7²` 1.4821s → 1.515s, `S_13.2` 5.0792s → 5.099s — matching the
headline table at the top of this file.)

Two-clause parity bar, one clause met: ✅ within 2× of C `schubmult` on
every row it finishes, ahead by 1.9× to 30.5× (`stair4²` excluded at
0.005s, the startup floor). ❌ complete `S_13 ℓ=25,36`, which schubmult
does not — **E2 does not complete it either**, running 428s and dying at
**6.56 GB** peak RSS.

**The bug was not the mathematics**: E2's memo, `HashMap<Perm,
Rc<Schubert<C>>>`, had no eviction, so every node's whole product stayed
retained for the run — 447 nodes against `S_13.2`'s 3.2M-term answer, each
holding hundreds of megabytes. Exactly the failure the performance targets
had already warned about for E3 — *a rising RSS curve means the
implementation is retaining the whole DAG rather than refcounting it, and
that is a bug, not a trade* — shipped in E2 anyway.

**Fixed the same day, for free**: `count_transition_parents` walks the
transition tree using only `Perm` operations, counting each node's
parents; `E2::release` drops a child's cached value once its last parent
has consumed it. No time regression anywhere (`stair7²` 1.482s → 1.515s,
`S_13.2` 5.079s → 5.099s — the pre-pass is `Perm`-only, and the drops
replace work the allocator would do at the end regardless). `S_13.2`, the
largest row that finishes — 3 241 903 terms in 5.27s — now peaks at
**1.52 GB**. `S_13.0` (`ℓ=25,36`) no longer dies but still does not
finish: RSS oscillates 1.2–2.4 GB for 20+ minutes instead of climbing to
6.56 GB and aborting at 428s. The memory bug is gone; what is left on that
row is genuine time.

### The pathological row, and why E2 wins where it wins

`S_13.0` (`ℓ=25,36`) defeats both engines while `ℓ=40,41` — longer, 3.4× the
output — finishes in 5s. Since *both* implementations fail the same row and
no other, the difficulty is a property of the pair, visible without running
a product (`dimension` is `S_w(1,…,1)`; `transition_tree` walks E2's own
recursion building no elements, costing microseconds even where the
product is hopeless):

```text
  case             deg           dim(u)           dim(v)  tree(u)  tree(v)   E2 time
  stair6²           42            32 768           32 768      633      633    0.070s
  stair7²           56         2 097 152        2 097 152    2 085    2 085     1.51s
  S_11.1 ℓ=35,26    61            13 305           11 067      376    1 102     0.40s
  S_12.0 ℓ=37,24    61         3 352 856        1 215 900    1 384    1 862     6.58s
  S_12.2 ℓ=43,32    75             4 228        1 396 206      157    1 126     0.63s
  S_13.1 ℓ=36,41    77           496 776            9 310      856      351     1.90s
  S_13.2 ℓ=40,41    81           585 004          832 723      663      447     5.10s
  S_13.0 ℓ=25,36 *  61        31 877 880    1 355 962 209    1 477    4 909   NEITHER
```

Three conclusions: **not the node count** (1 477 nodes, smaller than
`stair7²`'s 2 085 at 1.51s and comparable to `S_12.0`'s 1 384 at 6.58s — E2
visits few nodes there, each enormous); **not the degree** (`deg = 61`,
shared with `S_11.1` at 0.40s and `S_12.0` at 6.58s, *lower* than
`S_13.2`'s 81 at 5.10s); and **the all-ones specialization singles it out
by 10⁴** — `S_u(1,…,1)·S_v(1,…,1) = Σ_w c^w_{uv} S_w(1,…,1)` is an exact
identity with every term non-negative, so `dim(u)·dim(v)` is the product's
total monomial mass: **4.32×10¹⁶**, against a ladder maximum of 4.40×10¹²
elsewhere — the only quantity that separates the row at all.

⚠️ Two limits, since the quantity tempts over-use: not a runtime predictor
in general (`stair7²` at 4.40×10¹²/1.51s and `S_12.0` at 4.08×10¹²/6.58s
share the mass, not the time), nor a term-count predictor (mass ÷ terms
ranges 1.2×10³–4.0×10⁹ across completed rows). It flags a row as out of
family, in microseconds, and nothing more.

**Why E2 wins is not the node count.** On `S_11.1`
E2 uses *more* nodes than E3 — 1 102 against 571 — and is still 97× faster.
Both engines cost (nodes) × (size of the running element); E1 and E3
expand a factor into monomials and push the other through Monk chains, so
the element inflates toward answer size early and stays there, while the
transition recursion never expands — each node holds a product only as
large as that subproblem needs. The element-size term, not the node count,
was the whole gap to `schubmult`, and the compression measurement — which
only ever counted nodes — was structurally incapable of seeing it.

A profile of E2 on `stair7²` (~13 400 samples) is healthier than E3's ever
was: `mul_variable` 50%, `add_assign` 11%, cover scans 10%, `BTreeMap`
insert/remove/drop 8%, allocator 8%, memmove/memset 8% — the time is in
the mathematics, not the plumbing. Two follow-ups remain, neither needed
to hold the bar: the cover scans still allocate a `Vec<(u32, Perm)>` per
term per pass, and `BTreeMap::remove` alone is 504 samples, signed Monk's
cancellation churning entries in and out.

### Leaf dispatch into LR: verified, and its cost

Products where both factors are Grassmannian with the same descent k reduce
to LR with a k-row flag truncation: `AutoLr::schur_product`, then dropping
output partitions with more than k rows and re-indexing as Grassmannian
`Perm`s. **Verified**: exhaustive over same-descent pairs for k ≤ 4,
|λ|,|μ| ≤ 5 — 760 pairs, 0 failures — nine larger cases through k = 5, and
the shape ↔ permutation convention `code(w) = (λ_k,…,λ_1)`,
`w(i) = λ_{k+1−i} + i` pinned against `s_λ(x₁..x_k)` over 128 shapes. The
surviving term counts (2, 5, 16, 59, 247, 1111 for k = 2..7) independently
match the measured `stair4²`–`stair7²` product sizes above; a
different-descent pair produces non-Grassmannian output, so "same descent"
is a real condition, not a formality.

The identity's truth conceals what it costs to use:

```text
  k       LR terms   kept (≤k rows)   dropped   waste
  2              7                2         5   71.4%
  3             34                5        29   85.3%
  4            206               16       190   92.2%
  5          1 433               59     1 374   95.9%
  6         10 873              247    10 626   97.7%
  7         87 452            1 111    86 341   98.7%
```

`AutoLr::schur_product` has no row bound, so a `stair7` leaf costs 87 452
terms to keep 1 111 — 79× redundant work, eating the point of dispatching
to the fast engine. Leaf dispatch is therefore contingent on a row-bounded
LR product pruning inside the strip recursion rather than filtering the
finished expansion; until it exists, leaves stay on E1/E2.

**Not chosen: cotransition.** Knutson's co-transition formula
(arXiv:1909.13777, the Lemma on p.2) recurses *upward* —
`(xᵢ − y_{π(i)}) P_π = Σ P_σ` over covers — and computing `S_π` needs
dividing a sum of Schubert polynomials by a linear form (`Ring::div_exact`).
Whole products need no division here, so the route was recorded and
skipped, per the standing rule that being able to do a thing is not a
reason to.

## The single-coefficient query

`schubert_coeff(u, v, w)` runs E2 with **Bruhat pruning**: at every node,
terms not `≤ w` are discarded, sound because the signed Monk rule moves
strictly *up* the Bruhat order (both sums run over covers), so nothing
below the cut can reach `w`, directly or by cancelling against something
that does. Two necessary conditions are checked first and answer most
queries with no work: `ℓ(w) = ℓ(u) + ℓ(v)`, and `u ≤ w`, `v ≤ w` — all
three *verified*, not assumed, as named tests
(`product_support_lies_above_both_factors`,
`monk_covers_move_up_the_bruhat_order`,
`bruhat_is_a_partial_order_with_known_bounds`). The query agrees with
reading the coefficient off the full product exhaustively on S₄ × S₄ with
targets over S₆ — **including every zero answer**, since a pruning bug
shows up as a spurious zero rather than a wrong number.

```text
  S_13.2 ℓ=40,41   full product, 3 241 903 terms   5.08s
                   one coefficient (c = 130)       0.0088s      578×
  S_13.0 ℓ=25,36   full product                    IMPOSSIBLE (mass 4.3×10¹⁶)
                   one coefficient                 0.028–0.190s
```

**The row no engine can complete is queryable in a tenth of a second** —
nobody wants 2.9×10¹¹ terms, they want particular ones.

**Nonzero constants exhibited on the impossible pair.** Random `w` are
essentially never in the support, so sampling the group only ever returned
0. Walking the support works instead: **propose** by applying the Monk
chain of `x^{code(v)}` to `S_u` with a beam of 4 000 — candidates come out
automatically of the right degree and `≥ u` — and **dispose** by
re-checking each with the exact pruned query (the beam's own coefficients
are wrong, since truncation discards cancelling terms). `S_13.0` produced
3 545 candidates in 0.06s, 2 238 with `v ≤ w`, including:

```text
  c^[8,16,12,10,4,6,5,7,13,2,3,11,1,9,14,15] = 5
  c^[8,16,12,10,4,6,5,7,14,2,3,9,1,11,13,15] = 5
  c^[8,16,12,13,4,6,5,7,10,2,3,9,1,11,14,15] = 5
  c^[8,16,13,10,4,6,5,7,12,2,3,9,1,11,14,15] = 3
  c^[9,16,12,10,4,6,5,7,13,2,3,8,1,11,14,15] = 18
```

**Refereed on the computable neighbor**, the only way to trust answers on
a pair with no referee by construction: the same pipeline run on `S_13.2`
against its full 3 241 903-term product agrees on **400 candidates with 0
mismatches**, 214 nonzero — the zeros matter as much as the nonzeros, since
a pruning bug shows up as a spurious zero, not a wrong number.

So the claim is the strong one: **structure constants of a product that
cannot exist in memory, computed and checked.** `c^w_{uv} = 18` above is a
number no other package can produce by any route.

## Bindings, and stanley against newtrans

Ten Python entry points — `schubert_multiply`, `schubert_multiply_variable`,
`schubert_divided_difference`, `schubert_divided_difference_perm`,
`schubert_expand`, `polynomial_to_schubert`, `schubert_pairing(a, b, n)`,
`schubert_dimension`, `schubert_coefficient`, `schubert_monomial_mass` — all
through `guarded`/`escalate`, permutations normalized on entry exactly as
`part()` normalizes partitions. Checked against Sage: **415 checks, 0
failures** — products, the 1-based/0-based `multiply_variable` boundary,
divided differences, `expand`'s round trip, `dimension`, single
coefficients over S₅, stability under padded input. Run separately from
the engine tests, so a failure here can only be marshalling: it earned its
keep immediately, since four apparent failures were the *harness* looking
up a padded key in a dict Sage keys unpadded.

`stanley()` and `schubert_to_stanley_schur` closed the list, taking the
check to **427 checks, 0 failures**, including agreement with Symmetrica's
own `newtrans` on 11 permutations.

**Against `newtrans`: parity, then 1.9× from removing allocations.** As
first written `stanley` was 2.7× `newtrans` at S₁₄ and **0.7× at
S₁₆ — slower** — unsurprising, since `newtrans` is the one well-engineered
routine in Symmetrica's module (above). **Memoization would have bought
nothing, and was not written**: the obvious fix — `stanley` walked the
transition tree as a *tree*, so cache the shared subtrees — was tested by
measuring the sharing first (`examples/probe_stanley.rs`): **1.0×–2.1×,
mostly 1.0–1.3×**, the tree really is a tree, recorded as a dead end that
cost one measurement instead of an implementation. **Marshalling was not
the cost either** — 0.00241s in Rust against 0.00287s through the binding,
so the FFI is 16%.

What was left was ~180 ns/node: `transition()` heap-allocated three times
per node — `descents()` built a whole `Vec` for its last element, and
`covers_left` built a `Vec<(u32, Perm)>` re-collected by callers into a
`Vec<Perm>`. `Perm::last_descent` and `Perm::for_each_cover_left` remove
all three, ups pushed straight onto the work stack — **1.9× uniformly**
(S₁₈ 0.00241s → 0.00133s). `stanley` is now ahead of `newtrans` everywhere
measured: **1.3× at S₁₆** (was 0.7×), 4.4× at S₁₄, 2.5× at S₁₂. ⚠️ The
first row's 205× is Sage warmup, not a result — its very next call costs
0.0003s.

## The schubmult startup floor

Every ratio against schubmult C in this file rests on one number: how much
of its reported time is `exec`, not computation. Measured directly, min of
5:

```text
  startup floor (trivial input, 1 term out)        ~3.0–3.3 ms
  printing: pipe vs /dev/null, 30 143 terms              +1.0%
  printing: pipe vs /dev/null, 118 822 terms             −0.6%
  printing: pipe vs /dev/null, 185 284 terms             −0.4%
```

Two consequences, opposite in sign. **Output formatting is not a
confound** — the one worth worrying about, since it would scale with term
count, which is the axis the engine comparison varies; at 185 284 terms the
pipe-vs-`/dev/null` difference is inside noise, so the comparison axis is
clean. **Startup is ~3 ms of every row, and it inflates the incumbent**: the
schubmult rows at or below ~0.03s — `stair4²` (0.005s), `stair5²` (0.016s),
and both fast random S₁₀ rows above — are substantially measuring `exec` and
are **not valid targets**; beating them proves nothing about the engine.
Every ratio quoted elsewhere in this file uses rows ≥ ~0.1s, or excludes the
floor rows explicitly (`stair4²` in the E2 ladder above) — this measurement
is what licenses doing so.

## The 1-based check in `Perm::at` costs nothing

`Perm::at` is the innermost accessor on this path, and its "positions are
1-based" precondition was a `debug_assert` — so in the release profile as it
then stood, `at(0)` computed `0 - 1` as `u32::MAX`, missed the bounds check,
and returned `0` as if it were an answer. The failure policy's R2 does not
allow a public contract to be enforced only in debug, and the fix is an
unconditional `assert!`; the only question was what it costs on the hot path.

⚠️ **The premise changed under this while it was in flight.** The policy's
item 1 shipped `overflow-checks = true` in `[profile.release]`, so the
underflow is now caught and `at(0)` no longer returns anything. What the
`assert!` buys past that is the message: "positions are 1-based" instead of
"attempt to subtract with overflow". The measurement below was taken before
that landed and is unaffected by it — the branch it times is the same branch.

Measured with a throwaway probe (`examples/tmp_at_probe.rs`, deleted
after the run), release profile, best-of-5 spread reported:

```text
                             debug_assert      assert!
  at-loop, 240M calls        91.8–93.0 ms   90.1–91.9 ms
  stair6² product            69.4–71.2 ms   68.9–69.2 ms
```

Free, and marginally on the fast side of free — the branch is never taken
and folds away wherever `i` is a loop index. ⚠️ **The first attempt at this
measurement said 98 ms → 140 ms, a 43% regression, and it was an artifact**:
that probe ran a tiny mixed workload ahead of the tight loop, and the
resulting code layout, not the branch, moved the number. A synthetic loop
around a single `#[inline]` accessor is sensitive to layout at that
amplitude; the stair rows are the ones to trust, and they agree with the
corrected loop. The same promotion was applied to `transpose`,
`covers_right` and `covers_left`, which allocate a `Vec` and were never in
question.

## Offline oracle fixture

617 products `X_u · X_v` for u, v through S_4, against Sage's
`SchubertPolynomialRing`, committed and checked without Sage.

The expansion alone does not check the support, so the test also sweeps the
zeros: for u, v ∈ S_3, every w Sage does not list must come back zero. The
universe that sweep runs over is **measured by the generator** rather than
guessed at the test site — the `schubbound` records say u, v ∈ S_N have support
within one-line length M (1→1, 2→3, 3→5, 4→7). A zero-sweep over too small a
universe silently checks nothing, which is what makes the measured bound part
of the check rather than a convenience.

## Next

- **Memo policy for pair-keyed products.** `product_cached` works on
  `(Partition, Partition)`; permutation pairs are a much bigger key space
  with an unmeasured hit rate. Measure before paying for the table.
- **The hot representation.** E2's own profile on `stair7²` still spends
  ~24% of its time in `BTreeMap<Perm, C>` and the allocator, plus a
  `Vec<(u32, Perm)>` allocated per term per pass in the cover scans. The
  inline `[u8; 32]` `Perm` worth ~1.8× on E3 (above) was never carried to
  E2, the engine actually shipped, though the fix is engine-independent.
- **A hybrid dispatching E3 on predicted small-output/large-DAG cases.**
  Unvalidated: no case measured so far has E3 beating E2 outright, and
  `dim(u)·dim(v)` — cheap, but not a runtime predictor (above) — is not by
  itself a workable trigger.
- **A row-bounded LR product**, `schur_product_bounded(μ, ν, max_rows)`,
  pruning shapes exceeding `max_rows` inside `SkewLr`'s strip recursion
  rather than filtering the finished expansion. Without it, verified
  Grassmannian leaf dispatch wastes 79× the work at `stair7` (above), so
  leaves stay on E1/E2. Open: does the pruning actually cut the recursion,
  and does `rect`/`two_row`/`three_row` dispatch survive a row bound at all.
- **Refuse hopeless products, or attempt and let them die?**
  `monomial_mass` flags an out-of-family input (`S_13.0`'s 4.32×10¹⁶ against
  a ladder maximum of 4.40×10¹²) in microseconds, but the threshold is soft
  (above) and a caller with a big machine may want the attempt anyway.
  Leaning: expose the mass, refuse nothing — not yet a rule the code
  enforces.
- **Double Schubert polynomials** (two alphabets). Symmetrica has them as a
  copy-paste second recursion (`sb.c:1451`); the actively-maintained
  `schubmult` package centers on them. Out of scope until single Schuberts
  hold their bars above; the coefficient story needs a `QtPoly`-shaped ring.
- **Type B/C/D and Grothendieck/K-theory.** Recorded, out of scope. Type B
  lives entirely in expanded-polynomial land in Symmetrica's `bar.c`,
  sharing nothing with the engines here; Grothendieck polynomials (Knutson
  §7) share the same recursions with isobaric operators. Neither should
  force the module layout to move.
