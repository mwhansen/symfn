# Clean-room specification: skew Littlewood–Richardson expansion

## Why this document exists

symfn is MIT/Apache-2.0. Its LR engine was originally written after reading
`lrcalc` (GPL) — not by copying it, but the author had *access* to lrcalc's
source, including implementation code that ships inline in its headers.
Copyright infringement is established by access plus substantial similarity, and
clean-room practice exists to eliminate the access prong outright.

This is the specification half of a two-team clean room. It was written by
someone who has read lrcalc, and it deliberately contains **only**:

- textbook mathematics (definitions that predate lrcalc by decades),
- functional requirements (what the code must compute),
- performance requirements (how much work it may do),
- interface requirements (the API the rest of the crate already calls).

It deliberately contains **no** implementation technique: no data structures, no
loop shapes, no pruning tricks, no algebraic reformulations. Those are for the
implementation half to derive independently. If the independent implementation
converges on the same structures anyway, that is evidence those structures are
obvious consequences of the problem — which is the useful result either way.

## 1. Mathematical background

All of this is standard; see Macdonald, *Symmetric Functions and Hall
Polynomials* (2nd ed.), I.9 and Appendix A, or Fulton, *Young Tableaux*, Ch. 5.

A **partition** λ = (λ₁ ≥ λ₂ ≥ … ≥ λ_k > 0) is a weakly decreasing sequence of
positive integers. |λ| = Σλᵢ. Its **Young diagram** is the set of cells (i, j)
with 0 ≤ i < k and 0 ≤ j < λᵢ. We write μ ⊆ λ when μᵢ ≤ λᵢ for all i.

For μ ⊆ λ, the **skew diagram** λ/μ is the set-difference of the two diagrams:
the cells of λ not in μ.

A **semistandard filling** of λ/μ assigns a positive integer to each cell so
that entries are **weakly increasing left-to-right along each row** and
**strictly increasing top-to-bottom down each column**.

The **content** of a filling is the sequence (c₁, c₂, …) where cᵢ is the number
of cells containing i.

The **reading word** of a filling is obtained by reading rows top to bottom, and
each row right to left.

A word is a **lattice word** (equivalently: satisfies the ballot condition) if
in every prefix, the number of occurrences of i is at least the number of
occurrences of i+1, for every i ≥ 1.

A **Littlewood–Richardson (LR) tableau** of shape λ/μ is a semistandard filling
of λ/μ whose reading word is a lattice word.

**The fact this module computes.** The Littlewood–Richardson coefficient
c^λ_{μν} equals the number of LR tableaux of shape λ/μ with content ν. It is the
multiplicity of s_λ in s_μ·s_ν, and equivalently the coefficient of s_ν in the
skew Schur function:

    s_{λ/μ} = Σ_ν c^λ_{μν} · s_ν

Note that the content of any LR tableau is automatically a partition — this is a
consequence of the lattice word condition, not an extra requirement.

## 2. What to implement

A new module `src/skew_lr.rs` providing a backend named `SkewLr`.

**The crate does not currently compile on this branch.** The module was deleted,
but its callers, its memo table, its tests, and its benchmark entry were all
left in place. That scaffolding is not a hint about the algorithm — it is the
interface contract, and it tells you exactly what shape the module must have.
Compile errors are your task list. Nothing in that scaffolding describes *how*
to compute anything.

### 2.1 Primary operation

```rust
pub fn expand_skew(outer: &Partition, inner: &Partition) -> Vec<(Partition, u128)>
```

Returns every ν with c^outer_{inner,ν} ≠ 0, paired with that coefficient.

Requirements:

- **Sorted by ν** (using `Partition`'s `Ord`). This is part of the contract:
  callers compare results across backends for equality.
- **No zero coefficients** in the output.
- `inner ⊄ outer` ⇒ empty vector.
- `outer == inner` ⇒ exactly `[(empty partition, 1)]`, since s_{λ/λ} = 1.
- Results must be memoized. `src/memo.rs` already provides `skew_cached`, keyed
  on the shape; use it rather than inventing a new mechanism. Note what that
  key implies about §2.3(3).

### 2.2 Backend trait

`SkewLr` must implement the existing `LrBackend` trait (`src/lr.rs`), i.e. both
`lr_coeff` and `schur_product`. Read that trait's documentation for the contract
its implementors owe — in particular `schur_product`'s ordering guarantee.

### 2.3 Performance requirements

These are the point of the exercise. The existing backends are too slow for a
specific structural reason, and merely reimplementing them is not success.

1. **`expand_skew` must produce all of its coefficients from a single
   enumeration.** The obvious implementation iterates over all p(n) partitions ν
   of n = |outer| − |inner| and runs a complete independent search for each one
   — `NaiveLr.schur_product` in `src/lr.rs` still does exactly that, and it is
   what you are replacing. The cost of expanding a shape must not scale with the
   number of candidate outputs.

2. **`schur_product` must not be implemented as a loop over candidate λ calling
   `lr_coeff`.** The default method on `LrBackend` does exactly that; override
   it. A full product expansion should cost one enumeration, not p(|μ|+|ν|) of
   them. How to arrange that is for you to work out — it follows from the
   mathematics in §1 and from properties of skew diagrams.

3. **`lr_coeff` must not be more expensive than `expand_skew`.** Given the
   memoization in §2.1, a caller sweeping many ν against one (λ, μ) should pay
   for one enumeration in total.

4. The result must substantially beat both existing backends.
   `examples/bench_lr.rs` already has a column for `SkewLr` and will measure it
   once it compiles. For scale: s[7,6,5,4,3]² takes ~0.95s via `NaiveLr` and
   ~0.66s via `StripLr`. Both are beatable by a wide margin.

Beyond these, you are free — and encouraged — to find your own optimizations.
Pruning, normalization, representation choices, and recursion strategy are all
yours to decide. Measure rather than assume; `examples/bench_lr.rs` is the
harness and should be extended to cover the new backend.

### 2.4 Integration

`src/lib.rs`, `src/strip_lr.rs` (`AutoLr`) and `src/hopf.rs` (`skew_schur`)
already reference the module and will compile once it exists. Check that each
call site still makes sense given what you actually built — in particular
`AutoLr`, which selects the library default. If your measurements do not support
the choice encoded there, change it and say why. Base that on numbers you
produce, not on what the surrounding comments assert.

## 3. Correctness requirements

Correctness is not negotiable and is cheap to establish here, because the crate
already contains an independent implementation and two external oracles.

1. **Exhaustive agreement with `NaiveLr`**, the reference backend, over every
   product with |μ|+|ν| ≤ 7 and — separately — every skew shape λ/μ with
   |λ| ≤ 8. Compare *full expansions*, not sampled coefficients: a
   single-enumeration design fails by dropping or duplicating whole terms, which
   spot-checking a coefficient will not catch.
2. **The committed oracle fixtures must pass**: `tests/sage_oracle.rs` and
   `tests/lrcalc_oracle.rs`. Do not regenerate or edit either fixture — they are
   ground truth computed by independent programs, and editing one to make a test
   pass would defeat its purpose.
3. Cover the degenerate cases explicitly: empty `inner`, `outer == inner`,
   `inner ⊄ outer`, one-row and one-column shapes, and a disconnected skew
   diagram (one where the cells fall into pieces sharing no row or column).
4. The full suite (`cargo test`, and `cargo test --features bignum`) must pass, and
   `cargo doc --no-deps` must not introduce new warnings.

## 4. Constraints on how you work

The value of this exercise depends entirely on these being respected.

**Do not consult lrcalc in any form.** Not its source, headers, documentation,
website, papers describing its internals, or any third-party description of how
it is implemented. A copy exists on this machine under
`/opt/homebrew/anaconda3/pkgs/lrcalc-*`; do not open it. If a web search result
starts describing lrcalc's implementation, stop reading it.

**Do not read the previous implementation.** A module `src/skew_lr.rs` existed
and was deleted on this branch. Do not recover it: not via `git show`,
`git log -p`, `git checkout`, the reflog, `target/` build artifacts, or the
`main` branch. Treat it as not existing.

**Do not read `NOTICE.md`, `ROADMAP.md`, or `README.md`.** These describe the
previous implementation's internals — they name the techniques it used and one
of them spells out a key construction. Reading them would hand you the answers
this exercise is trying to derive independently, which defeats the point. The
same goes for any doc comment elsewhere in the crate that describes what
`SkewLr` does internally (there are a few, in `src/sym.rs`, `src/lr.rs` and
`src/strip_lr.rs`): note the *interface* and the *performance targets*, and
disregard any description of mechanism. Do not go looking for more of them.

**Everything else in the crate is fair game** — `src/lr.rs`, `src/strip_lr.rs`,
`src/partition.rs`, `src/memo.rs`, the tests, the fixtures. That code is symfn's
own. In particular `NaiveLr` in `src/lr.rs` is a correct, readable LR tableau
enumerator and a perfectly good starting point to think from.

If you find yourself having read something off-limits, do not pretend otherwise
— say so in your report. A clean room with an honest disclosure is still useful;
one with a concealed leak is worthless.

**General literature is fair game.** Macdonald, Fulton, Stanley, Sagan, the
original Littlewood–Richardson papers, and standard expositions of the LR rule
are all appropriate references. The mathematics is not anyone's property.

## 5. What to report

When done, report:

- the design you arrived at, and *why* — especially how you satisfied §2.3,
- benchmark numbers against `NaiveLr` and `StripLr` on the shapes in
  `examples/bench_lr.rs`, produced by running it,
- test results, stated plainly, including anything that does not pass,
- any requirement in this spec you could not meet, and what you did instead,
- anything you tried that did not work, and what the measurement showed.

Do not report success you have not observed. If a benchmark was not run, say it
was not run.
