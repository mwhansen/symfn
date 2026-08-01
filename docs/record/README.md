# The record

What was built, what it was measured against, and what was learned — one file
per subsystem, per the genre rules in [../style.md](../style.md). Nothing in
this directory is a plan: the forward-looking layer is
[release-readiness.md](../release-readiness.md) plus each file's open tail.
Until 2026-07-31 this directory was `docs/roadmap/` — a roadmap when nothing
was written, which became the record as the plans were executed; the name now
matches the function.

Five subsystems were built against written specifications
(`docs/spec-*.md`). Those specifications have been **merged into the
subsystem files here** and no longer exist separately: their measurements,
dead ends and verification ledgers are chapters below, their conventions were
copied into the module docs, and their planning scaffolding — proposed API
signatures, crate-fit tables of `file.rs:NNN` pointers, checklists now
realized as named tests — was dropped. Nothing cites a spec section number
any more; see [../style.md](../style.md), "Specs, and how they end", for why
that citation format was retired rather than repaired.

The clean-room LR specification is the one exception and is deliberately
*not* here — it stays at
[../cleanroom-spec-skew-lr.md](../cleanroom-spec-skew-lr.md) as frozen
evidence for the licensing story, not as record.

## Current state

**Phases 0–6 complete.** 407 tests green on the default build — 369 unit,
6 algebra-law, 23 oracle (5 lrcalc-fixture, 8 Sage-fixture, 6 non-field
coefficient ring, 4 in-house), 4 overflow-profile canary, 1 memory-budget and
4 doctests — and the default build still has **no external dependencies**.
`--features bignum` adds 12 more.

The memory-budget test is the one that is not about correctness: it runs all
twelve workloads in `symfn::measure::workloads` and checks peak live bytes and
allocation count against ceilings, in 0.7s. That is assertable where a timing is
not, because those two numbers are bit-reproducible on a single-threaded
workload — see [memory.md](memory.md).

All five bases exist as real types with multiplication; every ordered pair of
bases converts; ω, the Hall inner product, and the full Hopf structure (skew
Schur, coproduct, counit, antipode) are implemented and law-tested. The sixth
classical basis, `f_λ = ω(m_λ)`, was added later and makes `convert` total.

That was the target. What the crate actually covers now goes well past it —
Littlewood–Richardson at scales neither lrcalc nor Symmetrica reaches,
Hall–Littlewood, Macdonald, the (q,t)-Kostka table, the Macdonald operator
algebra, LLT, Jack with the Goulden–Jackson tables, and Schubert polynomials —
together with a PyO3 module that can stand in for Symmetrica underneath Sage.
All of that is recorded in this directory, one file per
subsystem, summarised and linked under
[The record, subsystem by subsystem](#the-record-subsystem-by-subsystem) below.

## The founding roadmap — base functionality to par

**Target ("at par"):** the five classical bases {m, e, h, p, s}, multiplication in
each, all conversions between them, the ω involution, the Hall inner product, the
partition conjugate — **plus the Hopf structure**: skew Schur functions,
comultiplication, counit, and antipode. All tested against Sage as oracle.
Everything from here through Phase 6 is that original roadmap, kept as the
record of the build; the checkboxes are all checked.

## Design decisions (settled)

1. **Introduce ℚ** via a `Ring` / `Field` trait split. Integral bases (s, h, e, m)
   stay exact over ℤ; only the power-sum basis, ω on p, and the inner product pull
   in ℚ. (Later: `num_rational::BigRational` under the `bignum` feature —
   `rug`/GMP was the original plan and was reversed; see Phase 5.)
2. **Hub-and-spoke conversions through Schur** for integral bases (compose to
   cover all 20 ordered pairs) rather than a full pairwise matrix.
3. **Tensor type** `SymTensor` (map `(Partition, Partition) → C` in the Schur
   basis) as the codomain of the coproduct.

---

## Phase 0 — coefficient & partition groundwork
*Small; unblocks everything.*
- [x] Split `Coeff`: `Ring` (add/mul, today's trait) and `Field: Ring` (+ `inv`).
- [x] Pure-Rust `Rational` implementing `Field` (i128 numer/denom, normalized).
- [x] `Partition::conjugate()` (transpose) and `Partition::z()` = ∏ iᵐⁱ·mᵢ!.
- **Done when:** rationals arithmetic tested; `conjugate` involutive; `z` matches
  known values (z of (2,1,1)=4, etc.).

## Phase 1 — the multiplicative bases
*Easy; mostly free.*
- [x] Add `Elementary<C>` (e) and `Homogeneous<C>` (h) basis types.
- [x] Multiplication for e, h = multiset union of parts (same rule as p).
- [x] m-multiplication (routed through Schur, in Phase 2).
- **Done when:** eλ·eμ, hλ·hμ tested; all five basis types exist.

## Phase 2 — the transition engine  *(the core)*
Per-basis-element conversion to/from the Schur hub, extended linearly; compose
for all pairs. A uniform `convert::<Target>()` on top.
- [x] `e ↔ h`: recursion Σ(−1)ⁱ eᵢ h₍ₙ₋ᵢ₎ = 0 (integral).
- [x] `s ↔ h`: Jacobi–Trudi `sλ = det(h₍λᵢ−i+j₎)`; inverse via Pieri (integral).
- [x] `s ↔ e`: dual Jacobi–Trudi (integral).
- [x] `s ↔ m`: Kostka `sλ = Σ Kλμ mμ`; inverse Kostka (integral).
- [x] `s ↔ p`: Murnaghan–Nakayama characters; `sλ = Σ zμ⁻¹ χλ(μ) pμ` (needs ℚ).
- [x] `m` multiplication via convert→Schur→multiply→convert-back.
- **Done when:** every ordered pair converts; round-trips are identity on a sweep
  of small partitions; transition values match Sage.

## Phase 3 — standard operations
*Small once Phase 2 lands.*
- [x] ω involution per basis: `ω(pλ)=(−1)^{|λ|−ℓ(λ)}pλ`, `ω(hλ)=eλ`, `ω(sλ)=s_{λ'}`.
- [x] Hall inner product `⟨·,·⟩` via p-expansion + `z` (check vs Schur orthonormality).
- **Done when:** ω²=id; ⟨sλ,sμ⟩=δλμ; ⟨pλ,pμ⟩=zλ δ across small degrees.

## Phase 4 — Hopf structure
*Leverages the existing LR machinery.*
- [x] `SymTensor<C>` type (element of Sym ⊗ Sym in the Schur basis) with linear ops.
- [x] Skew Schur `s_{λ/μ} = Σν c^λ_{μν} sν` (direct from `LrBackend`).
- [x] Coproduct Δ: on multiplicative bases `Δ(hₙ)=Σ hᵢ⊗h₍ₙ₋ᵢ₎`, `Δ(eₙ)` dual,
      `Δ(pₙ)=pₙ⊗1+1⊗pₙ`; on Schur `Δ(sλ)=Σ c^λ_{μν} sμ⊗sν`.
- [x] Counit ε (project to degree-0 coefficient) and antipode `S(sλ)=(−1)^{|λ|}s_{λ'}`.
- **Done when:** coassociativity on small inputs; `Δ` agrees across bases;
      antipode satisfies the Hopf axiom `m∘(S⊗id)∘Δ = ε·1` on small degrees.

## Phase 5 — validation & interop
*Makes it trustworthy and actually usable from Sage.*
- [x] `Ring::from_u128` injection, so structure constants (LR, Kostka, z_λ) are
      never silently truncated by a cast — the seam a bignum type needs.
- [x] Cross-cutting algebraic-law suite (`tests/algebra_laws.rs`): conversions are
      ring homomorphisms, round-trips are identity, ω is an involutive algebra
      map, Hall pairings ⟨s,s⟩/⟨h,m⟩/⟨p,p⟩ correct, Δ is an algebra map.
- [x] Sage oracle: `scripts/gen_sage_oracle.sage` generates a 743-line fixture
      (Kostka, characters, Schur products, all four conversions, skew Schur);
      `tests/sage_oracle.rs` checks symfn against it. Fixture committed, so
      `cargo test` needs no Sage.
- [x] `gmp` feature: `impl Ring for rug::Integer`, `Ring`+`Field` for
      `rug::Rational`. `from_u128` exact; verified on coefficients beyond i64.
- [x] `python` feature: PyO3 module (abi3, so one artifact serves CPython 3.9+),
      coarse-grained whole-object API. Verified importing into Sage (Py 3.14).
- **Done:** `import symfn` works inside Sage; 49 Schur products cross-checked
      in-process with zero mismatches; benchmarks show 2–20x over the pure-Python
      path (see below).

⚠️ **The `gmp` decision was later reversed and this checklist item describes
what was built, not what ships.** The escalation ring is now `num-bigint`
behind the `bignum` feature: measured coefficients past `i128` are 2–5 limbs,
where every library runs schoolbook and GMP's asymptotic algorithms never
engage, and `rug`/`gmp-mpfr-sys` are LGPL-3.0+ against a crate that stays
MIT OR Apache-2.0. See
[The Python boundary](python-and-sage-interop.md).

## Phase 6 — performance (memoization → plethysm → optimized LR)

- [x] **Memoization** (`src/memo.rs`). Thread-safe caches of *pure* function
      results: partition lists, characters, Kostka numbers, LR coefficients,
      inverse-Kostka data, and whole Schur products. Referentially transparent,
      so unlike Symmetrica's mutable globals nothing is observable in results.
      `clear_caches()` reclaims the memory.
- [x] **Plethysm** (`src/plethysm.rs`). Computed through the power-sum basis,
      where `p_n[g]` is just part-scaling and plethysm is multiplicative in f.
      Cross-checked against Sage on 27 values.
- [x] **Optimized LR backend** (`src/strip_lr.rs`). `StripLr`: a row-level DP
      over horizontal strips with state merging, replacing cell-by-cell
      enumeration. Verified against `NaiveLr` on every product with |μ|+|ν| ≤ 7.
- [x] **Whole-shape LR backend** (`src/skew_lr.rs`). `SkewLr`: expand a skew
      shape in one traversal, advancing a *merged frontier* of partial fillings
      rather than enumerating tableaux individually. Products go through the
      disconnected skew shape whose skew Schur function is s_μ·s_ν. Now the
      default; verified against both other backends exhaustively and against
      `lrcalc` on large shapes. Written clean-room — see `NOTICE.md` and
      `docs/cleanroom-spec-skew-lr.md`.

### Benchmark summary (vs Sage, same machine)

⚠️ These are the figures as Phase 6 closed, kept as the record of that point.
They are undated and uncontrolled, and every one of them has been superseded —
often by more than an order of magnitude — by the per-subsystem files below.
Do not quote from this table.

| workload | speedup |
|---|---|
| 144 small Schur products | **15.2x** |
| repeated large products | **8–10x** |
| single large products, cold | **1.2–4x** |
| 225 Kostka numbers | **62x** |
| plethysm (12 cases, cold) | **13.9x** |

---

## Shipping it — [release-readiness.md](../release-readiness.md)

A separate plan, and the one thing here that is forward-looking rather than a
record. This file tracks what the library computes; that one tracks what stands
between the tree and a package someone else can depend on — CI (there is none
today), the 173 rustdoc warnings, which of the 39 public modules are actually
API, the panic/overflow contract, crate and wheel metadata, and keeping the
wheel free of any Sage dependency with the Sage adapter layered on top.

Two supporting audits of the Sage side, both against 10.10.beta7:
[docs/sage-packaging-audit.md](../sage-packaging-audit.md) — can a *standard*
Sage package be a prebuilt Rust wheel? (yes; `rpds_py` is maturin-built and
standard, and Sage builds no Rust from source at all) — and
[docs/symmetrica-coverage-audit.md](../symmetrica-coverage-audit.md) — what
would displacing Symmetrica actually require? (Sage reaches 36 of its 66 entry
points from six files; symfn covers 34 of the 36 today).

## The record, subsystem by subsystem

Everything below Phase 6 used to live in this file, which had reached 3800 lines.
Each subsystem is now its own document in this directory,
with the prose unchanged; what follows is a summary and a pointer. The order is
roughly the order the work happened in.

Nothing here is a plan. These are records of what was built, what it was
measured against, and — as often as not — which predicted optimisation turned
out to be worth nothing.

### [Littlewood–Richardson](littlewood-richardson.md)

`SkewLr` expands a skew shape in one traversal, advancing a merged frontier of
partial fillings rather than enumerating tableaux, and is the default backend:
7.6–110.6x over `NaiveLr`, 11.1x over lrcalc on the largest shape both finish,
and `[24,20,16,12]²` (5 313 471 terms) completes in 148 s where lrcalc does not
finish at all. On top of that sit byte-packed inline keys, a conjugate-orientation
dispatch worth 7.9x on the big case, a sharded parallel merge worth 2.86x, and
per-output counting routes for two- and three-row factors. Correctness comes from
four independent directions: `NaiveLr`, lrcalc, Symmetrica, and a principal-
specialization checksum that ships with a negative control (412/412 perturbations
detected).

The last deficit — few-row factors, where the frontier compresses 1.0x and so
does a naive enumerator's work plus hashing — closed on 2026-07-31: a ~2x
cheaper fibre count (packed state, window-form transitions) dropped the
counting crossover to n ≥ 48, and every above-floor case in the comparison
sweep now measures ahead of lrcalc, 1.06x to two orders of magnitude —
confirmed on AC the same day (1.02–1.49x interleaved on the former loss band). The same session proved the ballot condition
survives column-by-column scanning (a plactic argument, verified exhaustively)
and measured why it does not help: any one-traversal DP carries partial content
in its state, and at three rows content pins the filling, so enumeration cannot
be beaten by merging in either scan direction — only per-output counting
escapes. This file is also where the project's measurement discipline was
learned, and most of it the hard way: battery versus AC is worth 2x and changes
ratios rather than just times, an undated table with no control silently
understated the library by up to 6x, a frontier-free enumerator that should
have won by the profile's own numbers turned out to be parity at best, and a
dispatch bound calibrated correctly in one machine era quietly inverted in the
next.

### [Transitions between the classical bases](transitions.md)

Kostka was exponential — 638 ms for a single K_{λμ} at degree 20, making
`convert_s_to_m` take 400 seconds for one Schur function. Counting chains of
horizontal strips instead of enumerating SSYT, so that chains through a shared
intermediate merge, made it polynomial: ~34 000x on the case that prompted it,
and degrees 30 and 40 became reachable. The coproduct had the identical defect
with its own fix one function away in the same file (15x), `convert_m_to_s` now
computes only the inverse-Kostka row it needs (3.8x), and `p → s` iterates
Murnaghan–Nakayama as a multiplication rule rather than evaluating p(n)
characters (2.0x).

Also here: `s → h` and `s → e` were a 200x regression on wide shapes, hidden
because the benchmark's shape family never let λ₁ exceed about 7. Reading
Symmetrica showed it runs the same algorithm — the whole gap was one line of
representation, its transposed Jacobi–Trudi orientation putting the tight row
first, so subtrees die at depth 1 rather than depth 12. And the sixth classical
basis `f_λ = ω(m_λ)`, which makes `convert` total over the standard set and whose
real work was finding three tests that do not merely restate ω.

### [Plethysm](plethysm.md)

Sage's plethysm is Python, but Symmetrica ships a C implementation Sage never
calls — so the benchmark that read **9x faster than Sage** was 0.14–0.79x against
the baseline that mattered, and a scaling defect sat inside a green result. Three
fixes in `p → s` (which profiling showed was 98–99% of plethysm's runtime), plus
a common-denominator integer sweep, put symfn ahead on every case and took
plethysm from 9x to ~40x against Sage.

The file also carries a retraction worth reading on its own: a cheap experiment
had concluded the rational leaf arithmetic was *not* the bottleneck, on a single
un-interleaved run. A sampling profile put rational arithmetic at 55%, the
re-run gave a 2.02x ceiling, and the fix that followed was worth 3.0x. The lesson
recorded is that an experiment used to **cancel** work needs the same rigour as
one used to justify it.

### [The Kronecker product, ordinary and reduced](kronecker.md)

The internal product is diagonal in the power-sum basis, so a famously hard
object — Kronecker coefficients have no known positive combinatorial rule — costs
`s → p` on both sides, a coefficientwise multiply weighted by z_λ, and `p → s`
back. Nothing enumerates anything. 434/434 products agree with Sage's `itensor`,
and the implementation is separately checked against the character-theoretic
definition, since it rests entirely on the p-basis identity and could otherwise be
self-consistently wrong.

The reduced (stable) coefficients come from the Orellana–Zabrocki `st` basis,
whose outer-product structure constants *are* the stable Kronecker coefficients.
Reading their Theorem 14 as a **map** rather than as a formula makes the change of
basis a tensor product of univariate triangular matrices, so both directions are
whole-element operations and the multiplication in the middle happens in the
power-sum basis — there is no Littlewood–Richardson coefficient anywhere in a
reduced Kronecker calculation. 3400x on the largest case Sage still answers, and
it reaches sizes Sage cannot: the wall is `z_γ` overflow in the intermediates, not
the answers, which stay under 20 bits.

### [Skewing and evaluation](skew-and-evaluation.md)

`g^⊥`, the adjoint of multiplication under the Hall inner product, generic over
the basis `g` is written in — because three bases have direct rules (Pieri, dual
Pieri, and Murnaghan–Nakayama read backwards) worth 2.2–268x over expanding into
Schur and using Littlewood–Richardson. The power-sum case is the one that most
repays it and not only for speed: rim-hook removal stays in ℤ, so `Schur<i64>` can
be skewed by a power sum, which the generic route could not have offered at all.
32 448 Sage checks, zero mismatches.

Evaluation at an alphabet is the bridge back to concrete numbers, one algorithm
per basis rather than convert-then-evaluate, with Schur going through the
branching rule. The bialternant is deliberately absent: it needs division, so it
is not generic over `Ring`, and it is 0/0 whenever two variables coincide — the
oracle script exercises exactly that hole. The closed forms for `f^λ` and both
principal specializations interleave their divisions with their multiplications,
which is not a micro-optimisation: for the staircase, `f^λ` has 35 digits and
fits `u128` while 55! has 74.

### [The Python boundary, and running as Sage's backend](python-and-sage-interop.md)

`scripts/sage_backend.py` fills Sage's own `conversion_functions` with symfn
shims, so Sage drives and the comparison is against the C library the shim
displaces — 4678 computations agree, and within minutes it found something no
conversion-table check could: Sage has two calling conventions, and every non-QQ
base ring goes through the one a dict-only backend fails on. End to end the honest
figure is **1.84x like-for-like**, 4.37x with a partition cache that Symmetrica's
wrapper could equally adopt; at these sizes a conversion is microseconds of
arithmetic wrapped in milliseconds of object marshalling, and a Cython per-term
loop took the shim from 91% glue to 40%.

The `i128` ceiling is closed by compute-and-escalate: checked arithmetic measures
0–1%, indistinguishable from noise, so the fixed-width path runs first and the
whole call re-runs over `BigInt`/`BigRational` if anything overflowed. Overflow is
recorded by a monotone global counter rather than a flag, because clear-run-test
loses answers under concurrency and the counter's failure direction is the safe
one. Separately, `Field` was replaced by `QAlgebra` in every dividing bound —
every division in this library is by a positive integer — which is what lets ℚ[t]
and ℚ[q,t] through, and those are exactly the rings Hall–Littlewood and Macdonald
need.

### [Oracles and comparison harnesses](oracles-and-comparisons.md)

Three external references with different powers, and knowing which is which is
the point. Half of Sage's ladder dispatches into Symmetrica's C and half is Sage's
own Python, so the harness prints a `via` column: 3.7x on a C row is a stronger
result than 9x on a py row. Sage memoizes, so cases must use distinct cold inputs;
symfn memoizes too, so `clear_caches()` runs before each timed call.

Driving Symmetrica directly — through bindings Sage installs but mostly leaves
unused — turned up two deficits a Sage-only harness could not see: per-value
Kostka and characters were 3–5x ahead while the **whole tables** were 2–2.5x
behind, the signature of Symmetrica computing a table as a table. Both were fixed
by the same observation that a table is p(n) sweeps and not p(n)² numbers, worth
6–7x each and turning both losses into 2.5–3.8x wins. On real-sized LR products
the same harness reads 91x, 442x and 9002x, where the small-degree ladder reads
3.2x — sizing a benchmark where the work lives changed the answer by three orders
of magnitude.

### [Hall–Littlewood and Kostka–Foulkes](hall-littlewood.md)

Reading Symmetrica's C before building changed the plan: its `hall_littlewood`
does not use charge at all, so `Q'_λ` comes from the Morris recursion over
machinery this crate already had and had already made fast, and charge became an
**independent** oracle rather than the engine — the more valuable role, since the
two routes then share no code. 3.36x Symmetrica end to end, from three
profile-driven changes worth 3.05x on the Rust path, none of which was the one
the plan predicted (that one paid 1.15x and is the smallest of the four effects
measured).

Kostka–Foulkes then falls out of the transition, since the polynomials *are* the
coefficients Hall–Littlewood already produces: 40x Sage asked per pair, and 884x
asked by column, which is the right unit. Also here is `QtPoly`, the sparse
bivariate coefficient ring, and a representation premise that was wrong twice
before being measured rather than assumed.

### [Macdonald polynomials](macdonald.md)

`P_λ(x; q, t)` by the branching formula, at ~94x Sage. The blocker was the
coefficient ring: a general fraction field over ℚ(q,t) needs a bivariate gcd,
which would have been the bulk of the work and none of the point. It is also
unnecessary — every denominator Macdonald produces is a product of binomials
`1 − qᵃtᵇ`, and that family is closed under product and lcm, so `Frac` keeps the
denominator factored as a multiset of exponent pairs and never expands it. One
consequence had to be handled rather than assumed away: the factored form is not
canonical, so `PartialEq` cross-multiplies.

`Q`, `J` and Hall–Littlewood `P` follow, the last by inverting the Kostka–Foulkes
matrix rather than by any new enumeration, and `q = 0` turning Macdonald `P` into
Hall–Littlewood `P` is the cross-check the whole layer was built for — two routes
with nothing in common. Two profiling passes were worth 3.8x then a further 5.9x,
the second contradicting three things the first had left in place; the sharpest
was that **100% of `mul` calls had a two-term operand**, so 39% of the profile was
quicksorting a concatenation of two already-sorted runs.

### [(q,t)-Kostka polynomials](qt-kostka.md)

Three routes to one table, all kept, and a benchmark asserts they agree at every
degree it times. The branching route is 2.8x Sage but on a losing curve (4.7x per
degree against Sage's 2.5x, crossover at n = 11). The Lapointe–Lascoux–Morse
eigenvector route fixes the curve — 3.1x per degree, 7.4x the branching route at
degree 12 — but the operator matrix turns out to be 0.01% of its own runtime, and
after the crossing into the Kostka basis the whole table still runs 0.5x. Then a
report on Sage's internals settled what all of it had been guessing at: **Sage has
not used LLM for this table since 2015**, so that comparison was against a third
thing entirely.

The engine Sage does use is Bergeron–Haiman's Pieri recursion, which is not a
per-shape enumeration but a recursion over pairs ordered by containment, where
every value is shared by every μ and ν whose recursion reaches it. That is now
17–26x Sage and is what every entry point reads. The two slow routes stay by the
crate's standing policy — the operator one twice over, since it shares no
*mathematics* with either alternative.

### [The Macdonald operator algebra](macdonald-operators.md)

∇, Δ_f, Δ'_f, Π and Θ_f. Sage has `nabla` and nothing else — measured, not
assumed, with the two near-misses identified numerically as different operators.
Every formula went into a Sage comparison before it went into the specification,
which caught a rendered `(−1)^{|μ|}` in a published eigenvalue that is invisible
on `∇e_2`, the first case anyone checks. The obvious route inverts a 77×77 matrix
over ℚ(q,t); it is not needed, because `H̃` is orthogonal for the star scalar
product and pairing against the *Schur* basis cancels the `z_ρ` outright.

~20–27x Sage — `∇e_13` in 16 s against 5m38s — against a target of 100x that is
recorded as the guess it was. Δ, Δ' and Θ have no oracle anywhere, so they are
tied to the one oracled operator by published identities instead. Three
performance claims in this section were written before the measurements meant to
support them and all three were wrong; they are corrected in place.

### [Labelled Dyck paths and the Delta conjecture](dyck-paths.md)

Both sides of [HRW] Conjecture 1.1: the **rise** version, a theorem, so a mismatch
is our bug; and the **valley** version, open, so a mismatch is a result. Both are
held against `deltaop::delta_prime_e`, and everything agrees for every k at every
n ≤ 9 — as whole symmetric functions, every `m_μ` coefficient, not the `h_1ⁿ`
slice Sage's parking functions can reach. The definition of `Val(P)` is not the
easy reading, and taking the easy one produces disagreement at every k except the
one where both sides collapse to the shuffle theorem, which is exactly the shape a
counterexample would have.

`k` never enters the enumeration, so the whole ladder comes from one walk (n = 9
went from ~12 minutes to 2.7), and the rise half now dispatches through
`src/llt.rs` for a further 29–56x. The valley half reads the labels directly, is
not an LLT statistic, and is therefore unchangeable — so the open side is now the
whole cost of testing the conjecture.

### [LLT polynomials, ribbon and vertical-strip](llt.md)

Three engines sharing nothing but the coefficient ring — SYT enumeration by
descent set, β-set ribbon strips on the abacus, and Fock-space straightening —
so agreement between them is evidence rather than restatement. Sage's is the
only other implementation anywhere and walls at whole-degree n = 10–11;
Symmetrica has no LLT at all, so this is a capability gap rather than a backend
swap. 700–19 100x on mains, and two of its outputs have no Sage entry point at
any speed: the `∇e_n` by-path Schur-positive refinement, and parabolic affine
Kazhdan–Lusztig columns by exact straightening.

The file is also where the convention traps are priced. Four normalizations of
`G` circulate and the literature reuses `G̃` for two different ones, every trap
being silent — a wrong-by-a-twist answer rather than an error. The min-inv floor
is the sharpest: its witness was written down as λ = (2,2) twice before anyone
checked, and the true smallest is (2,2,2), (2,2) having floor 0 — which the
convention gate elsewhere depends on, so the two claims were never consistent.

### [Jack polynomials and the Goulden–Jackson tables](jack.md)

Three engines over `AFrac`, the fourth factored fraction field in the crate and
the only **canonical** one: every scalar in the Jack calculus is a ratio of
integer-linear forms `uα + v`, and normalising the atoms primitive makes them
irreducible and pairwise coprime. The Laplace–Beltrami eigenoperator route wins
because it enumerates nothing at all: 3000–8000x Sage, n = 16 in 0.73 s where
Stembridge's SF ships precomputed archives and Sage cannot reach n = 12. The 200x
target was beaten by more than an order of magnitude, and the file records why —
it priced the arithmetic correctly and the incumbent wrongly.

`src/gj.rs` then computes the Matchings-Jack and b-conjecture coefficients, which
**no package computes**, pinned at both known specializations against objects with
independent definitions (the class algebra at b = 0, matchings counted directly at
b = 1 — deliberately not via zonal polynomials, which would have re-used Jack and
checked nothing). A modular-arithmetic second engine took the ladder from n = 10
to n = 14. Positivity is open for both families and is only *observed*; the file
also reports what fraction of its own output is not already covered by a theorem,
which is the argument that pushing degree further is the wrong next rung.

### [Schubert polynomials](schubert.md)

The one Symmetrica subsystem Sage still actively used, so retiring it completes
the displacement. Three engines, and the ranking inverted twice: a
pre-implementation measurement of state compression pointed hard at the peel-DAG,
which duly beat Symmetrica's route by 11.2x and then lost to the C `schubmult` by
4–51x, after which the memoized Lascoux–Schützenberger transition beat *it* by up
to 97x. The metric was the mistake and it is the transferable lesson — both
engines cost (nodes) × (size of the running element) and it counted only nodes, so
a cost model with a factor missing ranked them confidently and wrongly.

E2 is ahead of both incumbents on every row either finishes (5.0–30.0x the C
`schubmult`; Sage finishes none of them). The part no engine tuning could deliver
is `schubert_coeff`: E2 with Bruhat pruning answers structure constants for pairs
whose product **cannot be materialised** — one has a monomial mass of 4.3×10¹⁶ —
in about 0.04 s each. No other package has such a query at all.

### [Executing the failure policy](failure-and-overflow.md)

The record half of [../policies/failure.md](../policies/failure.md): what
changed in the tree to meet the rulebook, and what it cost. `[profile.release]`
now carries `overflow-checks = true` — measured interleaved against the same
tree without it at 0–6% on every harness the crate has, nothing on the bignum
routes, and no change to peak bytes or allocation counts — with a four-test
canary that fails the day the line leaves `Cargo.toml`, verified against a
build with the flag off.

Turning it on immediately found one: the allocation harness counted live bytes
unsigned, while `measure::reset()` zeroes the counter with earlier allocations
still held, so freeing them underflowed and the high-water mark latched ~2^64.
**A counter that is reset while its subject is still live is signed, whatever
it counts.** It also inverted the premise of an `#[ignore]`d test that had been
asserting the *wrong answer* release wrapping produced — that call now refuses,
in every profile.

The `# Panics` sweep that closed the policy's last documentation gap is the
same file's later chapter: four public functions in the crate documented their
panics, 46 now do, and every bare `.unwrap()` in `src/` is gone. It found a
public accessor that was not panicking at all — `Perm::at` returned
`w(0) = 0`, its 1-based precondition being a `debug_assert` — which the
release-profile flag above had already converted from a wrong answer into a
panic by the time it landed.

### [Memory: measurement and findings](memory.md)

Memory numbers here had been inconsistent because "memory" meant three
quantities that move independently — peak live heap, total bytes allocated, and
peak RSS — and only the last is easy to measure and the only one that drifts,
since it depends on allocation history. `src/measure/` is the shared accounting
(a counting `GlobalAlloc`, no dependencies) behind one catalogue of workloads
that feeds **both** an exploratory report (`examples/heapstat.rs`, with a
size-class histogram that attributes churn to a specific buffer) and a
regression test (`tests/memory.rs`). The rules that came out of it: **a memoized value returned by
clone is a design error** (`expand_skew` deep-copied 164 041 terms and 14.1 MB
per call to hand back what the cache already held — fixed by `expand_skew_shared`,
and `SkewLr::lr_coeff` was copying an entire expansion to read one coefficient);
and **churn is not a memory problem until it is shown to be one** — halving
`QtPoly`'s allocation churn with an in-place merge moved RSS by 6% and cost 16%
of the run time, because uniform, promptly-freed buffers are exactly what an
allocator recycles perfectly. That change is reverted and recorded, next to the
frontier-pooling experiment it rhymes with.

---

## Beyond the core (deferred, but intended)

The target above is the symmetric-function core. Symmetrica — the library this
one succeeds in spirit — also covers, and symfn does not:

- modular and projective representation theory of the symmetric group
- Schubert polynomials, commutative and non-commutative
- Hecke algebras of type A
- finite group operations
- ordinary representation theory of the classical groups

That is the eventual scope, not the current one. Nothing here is scheduled and
none of it should be read as implied by the phases above; the core comes first
and is where all the depth work (LR backends, Kostka, conversions) lives.

Recorded now for one practical reason: **Symmetrica is public domain**, so when
this work does begin it is a legitimate source of algorithms, not merely a
reference point — see `NOTICE.md`. It would also serve as a third test oracle
alongside Sage and lrcalc, with no licensing friction and coverage of
operations neither of those makes convenient.

---

### Dependency order
`Phase 0` → `Phase 1` → `Phase 2` (needs 0,1) → `Phase 3` (needs 2) →
`Phase 4` (needs LR + conjugate + ω) → `Phase 5` (needs it all).

### Notes
- Keep `NaiveLr` as the in-house LR oracle: it is the most obviously-correct
  backend, and the faster ones are held to exhaustive agreement with it.
  (Skew Schur and the Schur coproduct now go through `SkewLr` instead.)
- Every phase ships with tests; Sage supplies ground truth wherever a value is
  non-trivial (Kostka tables, characters, transition matrices).
