# Labeled Dyck paths, and both sides of the Delta conjecture

The combinatorial side of the Delta conjecture: the rise version (a
theorem, so a mismatch is our bug) and the valley version (open, so a
mismatch is a result). Both are held against the operator route.

Split out of [the record index](README.md), which carries the phase plan
and a summary of this file.

---

`src/dyck.rs` — the combinatorial side of [HRW]'s Conjecture 1.1, built because
the operators made it the next constraint rather than because the operators
needed it. Once `Δ'_{e_k} e_n` costs 0.1s at degree 8, what stops a search is the
enumeration, and this is it.

Both versions are produced: the **rise** version (a theorem, so a mismatch is our
bug) and the **valley** version (**open**, so a mismatch is a result). Both are
held against `deltaop::delta_prime_e`.

## Sage has the objects but not the statistics

`ParkingFunctions(n)` *is* labeled Dyck paths with distinct labels, and carries
`.area()`, `.dinv()`, `.to_labelling_area_sequence_pair()`. Its conventions were
confirmed to be [HRW]'s before being relied on — `Σ_PF q^dinv t^area` equals
`⟨∇e_n, h_1ⁿ⟩` for n ≤ 6 — so it is a genuine external oracle for the one slice
it covers, and `scripts/check_deltaop.py` now uses it (exact agreement through
n = 7, 262144 paths).

What it does **not** have is everything the conjecture needs on top: no `Rise(P)`
or `Val(P)`, no per-row `d_i(P)`, no `z`-extraction, and no *word* parking
functions — repeated labels — so no whole symmetric function, only the `h_1ⁿ`
coefficient. So Sage is an oracle for the `k = n−1` slice and not a substitute
for this module.

## Val is not `a_i ≤ a_{i−1}`

```text
  Val(P) = { i : a_i < a_{i−1} } ∪ { i : a_i = a_{i−1}, ℓ_i > ℓ_{i−1} }
```

It is **not** `a_i ≤ a_{i−1}`. Taking the easy reading makes the valley side
disagree with the rise side at every k except `k = n−1`, where nothing is chosen
and both collapse to the shuffle theorem — which is exactly the shape a
counterexample to an open conjecture would have. That was the first thing tried
here, and it is why the check reports the two columns separately.

## Two things that keep the inner loop cheap

- **The `z` extraction is `e_m` of the weights**, so it is a knapsack over the
  exponents, not a walk over `C(|Rise|, m)` subsets. One DP gives the whole
  selection polynomial.
- **The valley side visits negative powers of `q`.** `Σ_{i∈S}(d_i + 1)` can
  exceed `dinv(P)` by up to `|S| ≤ n`, and `QtPoly` has unsigned exponents, so the
  accumulation carries a uniform `q^n` and divides it out exactly at the end. A
  remainder there means the offset was too small, not that the conjecture failed.
  The rise side needs no offset: the chosen `a_i` sum to at most `area(P)`.

## Measured

Whole symmetric function — every content, hence every `m_μ` coefficient — both
versions, every k, against the operator:

```text
  n    paths (distinct labels)    per (n,k), rise / valley
  6                     16,807      0.005s / 0.006s
  7                    262,144      0.05s  / 0.06s
  8                  4,782,969      1.1s   / 1.2s
  9                100,000,000      25-49s / 29-51s
```

**Everything agrees, for every k, at every n ≤ 9.** The rise half is a check on
this crate; the valley half is evidence for an open conjecture, and it is the
whole symmetric function rather than the `h_1ⁿ` coefficient — every `m_μ`, so a
statistic that was wrong only on repeated labels could not hide.

⚠️ The growth is `(n+1)^{n−1}` and nothing here changes that. n = 9 is ~1×10⁸
paths and takes about 12 minutes for the full ladder; n = 10 is 2.4×10⁹ and is
out of reach by direct enumeration. Past that the enumeration needs a recursive
decomposition, not a faster loop — and that is now the binding constraint on the
whole line of work, not `deltaop`.

## One walk fills the whole k ladder

`k` never enters the enumeration — only the `z`-extraction, which is `e_j` of the
weights for `j = n−1−k`. The knapsack that computes one `j` fills the whole table
on the way, so `ladder` returns every `k` from a single walk over the labeled
paths. Asking per `k` re-walked everything `n` times: the full n = 9 ladder went
from **~12 minutes to 2.7 minutes** (78.6s rise + 83.2s valley, both versions,
all nine k, whole symmetric function).

## The decomposition that would actually change the exponent — and why it is not here

Grouping by Dyck path is already done; what remains is computing, for a
*fixed* area sequence, the sum over labelings without enumerating them.
That sum is a **vertical-strip LLT polynomial** — `∇e_n = Σ_D t^{area(D)} LLT_D(x;q)`
is the standard decomposition — so this is the same machinery as the LLT item in
`research-gaps.md`, and it belongs with that work rather than bolted on here.

⚠️ Two obstructions were identified before stopping, both worth recording because
the route looks routine until you try it:

- The usual device is **standardization**: expand in fundamental quasisymmetric
  functions indexed by standard labelings, of which there are only
  `(n+1)^{n−1}` — reducing "all word parking functions" to "parking functions"
  *and* yielding the whole symmetric function rather than the `h_1ⁿ` slice.
- But `dinv` is **not** invariant under either tie-breaking convention. A tie
  `ℓ_i = ℓ_j` contributes 0; breaking ties earlier-smaller adds a `dinv` in the
  same-diagonal case, later-smaller adds one in the adjacent-diagonal case.
- And the **valley** statistics depend on the labels directly: `Val(P)`'s tie
  clause (`a_i = a_{i−1}`, `ℓ_i > ℓ_{i−1}`) is decided by an equality that
  standardization destroys. Later-smaller preserves `Val` but breaks `dinv`.

So the rise side would standardize with care, and the valley side — the open
one, and the reason for the exercise — needs an argument that was not
attempted here.
Getting it wrong would silently produce false evidence about an open conjecture,
which is the one failure mode worth being slow about.

### Resolved, 2026-07-30 — `src/llt.rs`

Both obstructions are resolved, and the tuple model removes the first rather
than working around it.

- **The `dinv` tie-break does not arise in the tuple model.** Standardizing
  *labeled paths* has no `dinv`-invariant convention, which is what stopped
  this. But a Dyck path's labelings are in bijection with semistandard
  fillings of a tuple of vertical strips ([DA] Rem 2.2's reversal — components
  in *reverse* row order, verified pointwise per path to n = 6), and on that
  side [HHL] (82) is an **identity**: standardization preserves attacking
  inversions, no convention chosen. So `G_D` costs `#SYT` of the tuple, not
  `#labelings`, and yields the whole symmetric function rather than a slice.
  `llt::nabla_e_by_path` is the decomposition this section wanted; measured,
  all 16 796 paths of n = 10 with their Schur expansions in 84 s, exact
  against `deltaop::nabla_e(10)`, every piece Schur-positive.
- **The valley obstruction was correctly identified and is correctly out of
  scope.** `Val` is not an LLT statistic — no amount of LLT machinery reaches
  it — so the valley side stays here, with its honest enumeration. Nothing in
  `llt.rs` claims otherwise, and the open conjecture is not touched.

## Next

- ~~Vertical-strip LLT polynomials as a first-class object~~ — **done**,
  `src/llt.rs` / `docs/record/llt.md`. What it left behind: Python bindings for
  the module, and the [BHMPS] Catalanimal route for `∇` of a *general* LLT,
  which is the `research-gaps.md` row at line 255 and the open item at
  `docs/record/llt.md`.
- The same treatment for the **compositional** refinements, where the open cases
  are. This is the next item, and the one place the valley obstruction above
  still applies.

### The rise ladder now goes through the LLT engine

`ladder(n, Side::Rise)` dispatches to `rise_ladder_via_llt`: `Rise(P)` and its
weights `t^{−a_i}` are area-only, so the `z`-extraction factors out of the
labeling sum and the whole ladder is `C_n` LLT evaluations plus a knapsack.
Measured 29× at n = 8 and 56× at n = 9 (~2× per degree) against the labeled
walk, which survives as `ladder_at_content` — the oracle, and still the cheaper
route for one coarse content. The **valley** side is unchanged and unchangeable:
`Val` reads the labels. So the open side is now the whole cost of testing the
conjecture.
