# Schubert polynomials — the last subsystem Sage routed through Symmetrica

The last subsystem Sage still routed through Symmetrica, so retiring it
completes the displacement. Three engines, a ranking that inverted twice,
and a single-coefficient query that answers products no machine can hold.

Split out of [docs/record/README.md](../../docs/record/README.md), which carries the phase plan
and a summary of this file.

---

Spec: `docs/record/schubert-spec.md`. Built 2026-07-30 (`src/permutation.rs`,
`src/schubert.rs`). This is the one Symmetrica subsystem Sage still actively
used, so retiring it completes the displacement the README opens with.

**Three engines, and the ranking inverted twice.** E1 is Symmetrica's route
(expand a factor to monomials, Monk chain per monomial) kept forever as the
oracle. E3 walks the same peel recursion as a merged DAG. E2 is the memoized
Lascoux–Schützenberger transition. A pre-implementation measurement of *state
compression* — pipe dreams ÷ peel-DAG states, 687× on `stair7`, 125 599× on an
S₁₃ element — pointed hard at E3, so E3 was built first and E2 demoted. E3
duly beat E1 by 11.2× on `stair6²`… and lost to the C `schubmult` by 4–51×.
E2 then beat E3 by up to 97×.

**The compression metric was the mistake, and it is the transferable lesson.**
Both engines cost (nodes) × (size of the running element); the metric counted
only nodes. On `S_11.1` E2 uses *more* nodes than E3 (1 102 vs 571) and is 97×
faster, because E1/E3 expand a factor into monomials so the running element
inflates to answer-size early, while transition never expands. **A cost model
that omits a factor will rank engines confidently and wrongly** — the same
shape of error as `three_row`'s "structural idea on the wrong data structure",
but at the level of the measurement rather than the implementation.

**Against the incumbents** (single-threaded; schubmult C re-measured the same
day and reproducible to 5%; its rows carry ~3 ms of process startup, so rows
≤0.03s are measuring `exec` and are excluded):

| case | terms | Sage/Symmetrica | schubmult C | symfn E2 | vs C |
|---|---|---|---|---|---|
| stair6² | 247 | >120s | 0.355s | 0.070s | **5.0x** |
| stair7² | 1 111 | >120s | 20.6s | 1.51s | **13.6x** |
| S₁₁.1 ℓ=35,26 | 118 822 | >120s | 0.748s | 0.392s | 1.9x |
| S₁₂.2 ℓ=43,32 | 242 507 | >120s | 19.0s | 0.634s | **30.0x** |
| S₁₃.2 ℓ=40,41 | 3 241 903 | >120s | 105.8s | 5.10s | **20.8x** |

Ahead on every row either incumbent finishes. Sage's own engine finishes none
of these but `stair5²`, where it takes 6.85s against our 0.035s.

**Where it stops: S₁₅.** S₁₄ is routine (2.6s / 53s / 83s for 0.37M / 7.1M /
6.7M terms); the first S₁₅ pair takes 421.9s for 12.4M terms. Peak RSS 1.52 GB
on the 3.2M-term S₁₃ row.

**Coefficient growth: the "≤16" figure was an artefact of the incumbent.** It
came from the products Symmetrica could finish. Measured on our own engine:
**130** at S₁₃, **591** at S₁₄, **863** at S₁₅. Still far from `i64`, but
growing, and `guard` earns its keep.

**The part no engine tuning could deliver: `schubert_coeff`.** E2 with Bruhat
pruning — signed Monk moves strictly *up* the Bruhat order (verified, along
with `u ≤ w`, `v ≤ w` on every product support), so terms not `≤ w` cannot
reach `w` and are dropped. 578× faster than the whole product on `S_13.2`.
More importantly it answers pairs whose product **cannot be materialised**:
`S_13 ℓ=25,36` has a monomial mass of 4.3×10¹⁶ — the C `schubmult` fails it,
E2 fails it, and no machine holds the answer — yet its structure constants
come back in ~0.04s each. Nonzero ones exhibited (`c = 18` among them) by
walking the support with a beam and refereeing each candidate with the exact
query; the pipeline is validated on the computable neighbour against its full
3 241 903-term product, 400 candidates, **0 mismatches**. No other package has
such a query at all.

**Memory was a bug, not a trade.** The first E2 retained every node's product
and died at 6.56 GB. A parent-count pre-pass over the transition tree (which
builds no elements) plus drop-on-last-use fixed it with **no time regression
on any row**.

Correctness: 355 lib tests — engines cross-checked exhaustively on S₅ × S₅
(E1 ≡ E2 ≡ E3, both argument orders), `∂` on the basis against `∂` by exact
division on expanded polynomials, the defining `w₀` recursion, Bruhat order
laws, and `stanley` against the §3.8 hand value and Grassmannian vexillarity.
Bindings separately: **427 checks against Sage, 0 failures**, including
`newtrans`. `stanley` started at *parity* with `newtrans` — the one routine in
Symmetrica's Schubert module that was well engineered — and is now **1.3× at
S₁₆, 4.4× at S₁₄**, from removing three heap allocations per node
(`Perm::last_descent`, `Perm::for_each_cover_left`). Notably **not** from
memoization: the obvious fix, and the defect §3.1 records for Symmetrica, was
walking the transition tree as a tree — but its sharing measured 1.0–1.3×, so
the memo was never written. One measurement instead of an implementation.

Open, in `docs/record/schubert-spec.md` §7: a row-bounded LR product to make
Grassmannian leaf dispatch pay (the flag discards 98.7% of the expansion at
k=7), the cover scans still allocating per term per pass, and whether `mul`
should refuse out-of-family inputs or merely report their mass (leaning:
report, refuse nothing).
