# Executing the failure policy

[docs/policies/failure.md](../policies/failure.md) is the rulebook — every way
a computation here may fail, and which mechanism each situation demands. This
file is the other half: what changed in the tree to meet it, what the changes
cost, and what turned up along the way. The policy's "What this changes" list
is the agenda; each item lands here as a chapter when it is executed.

**State.** All seven items are executed. `[profile.release]` carries
`overflow-checks = true`; the guard's `i128::MIN` corners report; every
fixed-width injection checks; the Python boundary raises where it panicked;
`src/` is cast-clean and gated at deny in CI; every `(q,t)` family states its
measured reach; and the crate front page carries the contract. Four of the
seven turned up live defects rather than gaps, and each has its chapter below.

Two things followed from the measurements rather than from the policy. Every
`(q,t)` wall a caller reaches in under a minute now has an escalation ladder —
`llt_h` at μ = 1ⁿ (n = 87), Hall–Littlewood on one shape, and all three
Macdonald forms at λ = (n) (n = 26–30) — while the whole-degree tables keep
none, their walls being 2–4× further out in degree than anything that finishes.
And CI has a release lane, because the tests that pin the profile flag are the
ones that only mean anything there.

A later pass audited the three rules the seven items never covered. **R4**
holds at all 27 wrapping/saturating sites, each for a different reason, and the
one `pub` function that saturates now says so. **R10** holds too, but
`algebra_laws.rs` compares two sides that share `i64` and survives only on an
a-priori bound that nothing stated — measured at 120 against `i64::MAX`, and
now pinned by a test that fails when a raised degree cap spends the headroom.

**R6** is audited too, and held the one live defect of the three: inside a
guarded scope `powersum_scalar` formed z_μ in native `u128`, so `jack_scalar`
past |μ| = 34 **panicked** where the ladder was watching for `None` — a crash
on an input the wide pass answers exactly. The generalizable half is that
`overflow-checks` (item 1) did not remove this failure mode, it changed its
shape: silent-wrong became loud-crash, and inside a fast pass loud-crash is
still wrong, because `escalate` cannot catch a panic.

Still open, with premises recorded in [Open](#open): the two-tier cache (specified, deliberately unbuilt), the `# Panics` sweep across
the whole public surface, clippy's own 155-warning backlog, and the fact that
CI has never actually run.

Every number below is from one machine — macOS arm64, rustc 1.96, on AC. CI now
exists but has never executed, so that caveat still stands
([release-readiness.md](../release-readiness.md), Phase 0).

## The release profile carries `overflow-checks` (policy item 1, R3)

Before this, `Cargo.toml` had no `[profile.release]` at all, so the profile
users ship wrapped: `impl Ring for i64/i128` multiplies with a plain `*`, and
the Hall–Littlewood, Kostka–Foulkes, Macdonald, qt-Kostka, nabla/delta and LLT
pyfunctions instantiate at plain `<i128>` with no escalation. Past their walls
they returned confident nonsense.

### What it costs

Measured by building the same tree twice — once with
`CARGO_PROFILE_RELEASE_OVERFLOW_CHECKS=true`, once without — and running the
two binaries **interleaved**, alternating which goes first each round, per the
discipline `examples/bench_guarded.rs` and
[littlewood-richardson.md](littlewood-richardson.md) both had to learn the hard
way. Minimum per (build, case) over 4 rounds (2 for the bignum harness).

| harness | workload | off | on | on/off |
|---|---|---|---|---|
| `bench_shapes` | `[16,13,10,7]²` (390 075 terms) | 0.3491s | 0.3507s | 1.00x |
| `bench_shapes` | `[8,7,6,5,4,3]²` (164 037 terms) | 0.1272s | 0.1215s | 0.96x |
| `bench_shapes` | `[20,16,12]²` | 0.0578s | 0.0582s | 1.01x |
| `bench_ops` | `kostka_all_pairs_n20` (393 129 values) | 1.9783s | 2.0989s | 1.06x |
| `bench_ops` | `convert_m_to_s` | 0.2374s | 0.2471s | 1.04x |
| `bench_ops` | `convert_p_to_s` | 0.0116s | 0.0126s | 1.09x |
| `bench_ops` | `character_table_n18` (148 225 values) | 0.0706s | 0.0691s | 0.98x |
| `bench_guarded` | combined, plain ring | 2.5369s | 2.5732s | 1.01x |
| `bench_guarded` | dividing (`s → p`, plethysm, kron) | 0.0985s | 0.1022s | 1.04x |
| `bench_kron_coeff` | product route, n = 32, `BigRational` | 42.40s | 42.98s | 1.01x |
| `bench_kron_coeff` | character sum, n = 44, `BigRational` | 1.02s | 1.01s | 0.99x |

Nothing exceeds 1.09x, and the rows that "win" say what the rows that lose
also say: most of this is noise around a real cost of a few percent. The
pattern that is not noise is where the cost concentrates — the rational paths
(`convert_p_to_s`, `bench_guarded`'s dividing block, the `hall_s_p_degree17`
pairing at 1.08x), where a single `Rational` op is four fixed-width multiplies
and a gcd, so there are simply more checks per unit of work. The bignum routes
pay nothing, as expected: `BigInt` arithmetic is function calls into
`num-bigint`, and the flag never reaches inside them.

`tests/memory.rs` was run under both builds to rule out an allocation
confound. Peak live bytes and allocation counts are identical workload for
workload, except for `skew-big`, which fans out across threads and moves ~1%
either way by the shard split it has always moved by.

**The burden of proof sits on turning it off, not on turning it on.** That is
the policy's rule (R3) and the measurement is why it costs nothing to hold: a
hot loop the flag visibly slows gets explicitly checked or proven arithmetic,
not the flag removed.

### The bug the flag found: `measure` counted live bytes unsigned

`cargo test --release --all-features` failed the moment the flag went on, in
`src/measure/mod.rs` — "attempt to add with overflow" in the allocator's own
high-water tracker, and reproducibly, but only when the test harness captured
output.

`measure::reset()` zeroes the live-bytes counter, and the measurement then runs
with memory allocated *before* the reset still held — the harness's capture
buffer is the reliable example. Every one of those blocks is eventually freed
against a counter that no longer counts it, so the counter goes below zero.
Unsigned, that underflowed to ~2^64 and `PEAK.fetch_max` latched it. **A
memory budget could have reported 16 exabytes**, or, more quietly, any garbage
above the true peak, and nothing in the harness would have said so — the
budgets are ceilings, so the failure mode is a red test with an absurd number,
but the number itself is the thing this harness exists to produce.

The fix is to count live bytes as `isize`, which is what they are: the counter
legitimately goes negative, and the peak is floored at 0 on the way out.
`tests/memory.rs` pins it with a region that only frees, ahead of the budgets
so that a regression is not read as a workload growing.

**A counter that is reset while its subject is still live is signed, whatever
it counts.** The unsigned type was not a micro-decision that happened to be
wrong; it was a claim — "this only ever goes up from here" — that `reset()`
had already falsified.

### The test whose premise inverted

`ops.rs` carried `unguarded_fixed_width_is_wrong_where_the_guarded_path_escalates`,
`#[ignore]`d, asserting that at n = 40 the unguarded `Rational` route returns
*the wrong answer* — a fraction, where the Kronecker coefficient is 1. It was
release-only, and the doc comment said why: in debug the same call panics
instead, so the wrong answer it pinned existed only in the profile users ship.

That is exactly the hazard the flag removes, so the test's premise inverted
with it. It is now
`unguarded_fixed_width_refuses_where_the_guarded_path_escalates`, a
`#[should_panic]` test that runs in the ordinary suite in both profiles, paired
with `guarded_path_escalates_where_fixed_width_refuses` for the half a
`should_panic` body cannot check. It doubles as a canary: it is the only test
that exercises the flag on a real computation rather than on a synthetic
product.

The synthetic ones are `tests/overflow_checks.rs` — four `#[should_panic]`
tests through `Ring::mul` / `Ring::add_assign` on `i64`, `i128` and
`Rational`, `black_box`ed so the optimizer cannot const-fold the overflow into
a compile error. Verified to fail as a set against a build with
`overflow-checks = false`, which is the only property that makes them a pin
rather than four tests that happen to pass.

## The `i128::MIN` corners (policy item 2, R5)

`Guarded` is a promise: `guarded(|| …) -> Some(v)` says every intermediate
stayed inside the fixed width. Two operations broke it on one input.

`Guarded::neg` and `GuardedRat::neg` used `wrapping_neg`. `-i128::MIN` is not
`i128::MIN` mathematically — it is the one negation the width does not have —
so on that value the wrap is a **silently wrong sign**, and, uniquely, one that
`guarded` hands back inside a `Some`. It is reachable: `checked_mul` returns
`MIN` happily, since `MIN` is a perfectly good product (`(MIN/2) · 2`), and
that path reports nothing because nothing overflowed. Both now route through
`checked_neg` to `note_overflow`, and `negating_the_width_minimum_is_reported_rather_than_wrapped`
pins both halves — the product that must succeed, and the negation that must
escalate.

Both `gcd`s took `.abs()` of possibly-`MIN` arguments, which is the same
non-existent negation. The guard's now computes over `u128` via
`unsigned_abs`, which is total, and `GuardedRat::new` refuses `MIN` in either
part before narrowing back — a report, not a value, because normalizing needs
`|num|`, `|den|` and possibly a sign flip and `MIN` has none of them. That
narrowing carries its bound proof at the site, per R5. `GuardedRat::from_i128`
refuses `MIN` at the seam for the same reason, rather than storing one for a
later operation to discover.

`Rational` — the unguarded fixed-width field — takes the other branch of R5:
it has no report channel, so its corner is a **documented wall**. `gcd`,
`Rational::new` and `Rational::neg` assert with one shared message naming the
requirement ("a Rational part of i128::MIN has no negation in i128; use the
bignum ring"), and `neg` checks rather than trusting the constructor because
the arithmetic fast paths build the fields directly.

One more narrowing turned up next to them, and it was not a corner case at
all: `Rational::div_u128` did `let n = n as i128`. The divisor there is `z_μ`,
which reaches `|μ|!` — a `u128` past `i128::MAX` becomes a *negative* divisor
and the answer comes back with a flipped sign and no signal. It is now
`try_from` with a panic naming the divisor, and
`dividing_by_a_divisor_past_the_width_refuses` pins it.

**The value a `checked_` operation legitimately returns can still be one the
next operation cannot use.** `checked_mul` did its job on every input here;
what was missing was that `MIN`'s *successor* operations are the partial ones.

## The injection seams refuse (policy item 3, R8)

`Ring::from_u128` / `from_i128` exist precisely because structure constants —
LR coefficients, Kostka numbers, `z_λ`, characters — are naturally wider than
`i64`, and routing them through `from_i64` would truncate. The seam existed;
the fixed-width implementations of it did not honor it. `impl_ring_for_int!`
wrote `n as $t`, `Rational::from_u128` wrote `n as i128`, and the trait's own
*defaults* wrote `Self::from_i64(n as i64)` — so the mechanism built to stop
silent truncation performed one at every fixed-width leaf.

All four now check. The message names the constant and the ring ("the
structure constant 340282366920938463463374607431768211455 does not fit i128;
use the bignum ring"), because a panic is documentation printed at the worst
moment. `Guarded` and `GuardedRat` keep reporting-and-escalating, which is R8's
other branch and was already right.

Cost, `bench_ops` interleaved before/after, min of 3 rounds each: ≤1.03x on
every workload of substance, with the two 1.05x rows being 20 µs cases. That
is the R8 claim — "injection is never a hot loop, so the check costs nothing
that matters" — measured rather than assumed.

The trait defaults are worth the separate note. Nothing in the crate uses them
today, since every concrete ring overrides both; they are a *future*
implementor's silent truncation, and the doc now states the obligation
(check and panic, or report and escalate) next to a default that meets it.

**A seam that exists to prevent a mistake still has to make the mistake
impossible.** This one was documented, universally routed through, and
truncating.

Also inverted here: `tests/bignum.rs` asserted that
`<i64 as Ring>::from_u128` *loses* a value past the width — "which is why the
seam exists". That assertion was pinning the bug. It is now
`from_u128_past_i64_refuses_rather_than_truncating`.

## The panic audit (policy item 4, R1/R2)

[release-readiness.md](../release-readiness.md) Phase 3 counted 138
`panic!`/`unwrap`/`expect` sites in `src/`. Re-grepped at audit time, outside
`#[cfg(test)]`, it was **93**. Two clusters were 43 of them, and only one
cluster was a real defect.

### The Python boundary was panicking on caller input (the defect)

`build_schubert(&a).unwrap()` on the escalation path. `build_schubert`
returns `None` for two unrelated reasons — a coefficient too wide for the
fixed-width pass, and **a one-line word that is not a permutation** — and the
`unwrap` treated both as impossible. The fast pass declined the malformed word
by returning `None`, escalation ran, and the `unwrap` fired: a
`PanicException` in Sage, which R2 names as by definition a bug report and
never an interface. Any coefficient small enough to fit — nearly all of them —
took that route.

The fix is structural rather than a better message. A `Wide` trait marks the
rings that *cannot* decline an input (`BigInt`, `BigRational`), and
`build_wide` / `build_rat_wide` / `build_schubert_wide` build over them with no
`Option` to unwrap; the one-line words are validated once, up front, by
`schub_terms`, so a malformed word is a `ValueError` before either pass runs.
That removed all 21 `unwrap`s in `python.rs`, and the two that were live bugs
with them.

**An `unwrap` is a claim that nothing checks.** Both readings of the `Option`
were true of the same call, and the type system was the only place the
difference could be recorded.

The pin lives in `scripts/check_schubert_bindings.py`, not in the Rust suite:
the `extension-module` build has no interpreter to link, so a unit test that
so much as constructs a `PyErr` aborts the test binary at load
(`symbol not found in flat namespace '_PyExc_BaseException'`). It exercises
both coefficient sizes, because it was the escalation path that panicked.

### The cache locks (not a defect, but a bad failure direction)

22 of the 93 were `RwLock::read().unwrap()` in `memo.rs` — lock poisoning.
Poisoning here means a neighbour panicked while holding a guard; it does not
mean the table is unsound, since every value is a pure function of its key and
`compute` runs *outside* the guard. Propagating it converts one thread's
failure into a permanent crash of every cached path in the process — a caller
who overflowed an `i128` and caught it would find the library dead. Two
helpers (`rd`, `wr`) now clear the flag, with the reasoning at the definition.

### The rest

The remaining ~50 were in better shape than the count suggested: mostly
`expect` with the invariant stated, which is the R2 form. What the audit
changed there was the bare ones — `chain.last().unwrap()`,
`by_degree.remove(&n).unwrap()`, `c.try_into().unwrap()` — which now name the
invariant they rely on ("the chain starts at the empty shape", "n came from
the map's own keys", "chunks_exact(8) yields 8 bytes").

Two reachable walls gained `# Panics`: `expand_skew` (an LR multiplicity past
`u128` — far beyond any shape that fits in memory), and
`reduced_kronecker_product` (`|λ|+|μ| = 24` without `bignum`). One turned up
that nobody had noticed: `class_algebra_coefficient` computes `n!` in `i128`
and so walls at **n = 34** — silent wrapping until item 1, now a panic, now
documented. Its `dimension(..).expect("a partition has a dimension")` was
claiming a proof it did not have; `dimension` declines past `u128` at |λ| ≈ 55.
The claim is true only because the `n!` wall fires two decades earlier, and the
comment now says *that* instead.

### The `# Panics` sweep, and the accessor it found

The audit above closed the reachable walls it found; it did not visit every
`pub fn`, and this is the pass that did — [style.md](../style.md)'s delta 2.
Before it, **four** public functions in the crate documented their panics
while public functions asserted throughout. After it, every `pub fn` in `src/`
outside `python.rs` that can panic carries a `# Panics` section: 46 added, of
which 8 had said "panics if …" in prose without the section rustdoc renders.

```text
                              before   after
  panic-family sites              66      38
    of which bare .unwrap()      ~30       0
  pub fn that can panic,
    with no `# Panics` section    46       0
```

⚠️ **That 66 is not the 93 the audit above quotes, and neither is wrong.** 93
counted `python.rs`; 66 excludes it (its boundary work had already landed) and
counts the panic family alone, with 79 assert-family sites beside it. Any
later comparison needs to say which of the three it means.

Sections **point at the exposed bound rather than restating it** —
`MAX_SUPPORT`, `MAX_CELLS`, `MAX_FREE_EDGES`, `abacus_reach` against
`ABACUS_REACH_LIMIT` — on R11's reasoning for the Python boundary: a bound
copied to a second site is a bound that drifts.

**`Perm::at` was returning `w(0) = 0`.** Its "positions are 1-based"
precondition was a `debug_assert`, so in release `i - 1` wrapped to `u32::MAX`,
the bounds check missed, and the `None` arm handed back `i`. The module doc
already called 1-based indexing a documented source of silently-wrong results
in Symmetrica — `mult_schubert_variable` is 0-based there and
`divdiff_schubert` is 1-based — and it was one here too. ⚠️ **Item 1 changed
this finding's consequence while the work was in flight**: with
`overflow-checks = true` shipped, `at(0)` panics on the subtraction instead of
returning `0`, so the wrong answer is already gone and what the unconditional
`assert!` now buys is a message naming the requirement rather than "attempt to
subtract with overflow". The finding stands as evidence for item 1, not as a
live defect past it.

Measured free even though `at` is the innermost accessor on the Schubert path
— see [schubert.md](schubert.md) for the table, and for the 43% "regression"
that turned out to be code layout. `transpose`, `covers_right` and
`covers_left` took the same promotion and were never in question.

**Four more debug-only public preconditions**, all the same degenerate-atom
condition that `Atom::unit`, `Frac::inv_factor` and `afrac::split` already
asserted hard: `QtPoly::mul_binomial`, `mul_diff`, `Frac::from_factors` and
`mul_factors`. `from_factors` was the one with teeth — `mul_binomial` covers
its numerator branch, but a negative exponent inserts `(0,0)` straight into
the *denominator*, and nothing downstream catches a zero denominator factor.

**A wall that was not one.** `stanley`'s `steps < 1 << 24` reads like a
capacity limit and is not: it is a backstop against a non-terminating
transition recursion, which without it would hang rather than fail. Its
section says so, and says the reach is unmeasured — a first draft claimed the
cap was unreachable because the tree is bounded by `ℓ(w)` under `MAX_SUPPORT`,
which is **not established**; the transition tree is not bounded by `ℓ(w)`.

**Where the fix was to delete the panic rather than document it.** Four sites
asked one question twice — an emptiness test standing next to the `Option`
that answers it. `charge`, `divide_by_factor`, the `s → s̃` pivot loop, and
`convert`'s degree drain now ask once, via `let … else` or `while let`, and
the `unwrap` has nowhere left to live. That is the cheaper end of R2: a site
that needs no panic beats a site with a well-worded one.

## The cast audit, phase 1 (policy item 5, R5)

`overflow-checks` does not reach `as`; nothing does. The three lints that see
narrowing are now on as warnings in `[lints.clippy]`, which required installing
clippy — **it had never run on this codebase** (release-readiness Phase 0). Its
default backlog is 155 warnings and is that phase's problem, not this one; the
three cast lints add **~370 more in `src/`**, and the plan's instruction was to
take the coefficient-adjacent modules first and leave the bounded `u32`
partition bookkeeping.

Audited and clean: `convert.rs`, `eval.rs`, `character.rs`, `llt.rs`. What the
audit found in 117 sites across them was **one** narrowing that carries a value
rather than an index:

```rust
acc += w[m] * kostka(&parts[m], &parts[jj]) as i128;   // inverse_kostka_row
```

`kostka` answers in `u128`. Past `i128::MAX` that `as` produces a *negative*
Kostka number, and the alternating sum it feeds absorbs the sign without a
trace — the inverse-Kostka row would come back wrong rather than absent. It
checks now, and names the constant, which is the same rule R8 applies to the
same numbers at the `Ring` seam.

Everything else in those four modules is β-set arithmetic, cell coordinates,
hook offsets, and q-exponents: bounded by `|λ|`, `ℓ(λ)` or the mask width, all
`u32` where they are stored. Those carry `#[allow]` plus the one-line proof, at
function scope, so the module still warns if a *new* cast appears somewhere
else in it.

`llt.rs` (67 of the 117) takes a module-level allow instead, and the reason is
worth stating because it is a structural argument rather than an inventory:
coefficients there are `QtPoly<C>` over a generic `C: Ring`, and **a generic
parameter cannot be `as`-cast at all**. A value-carrying narrowing cannot be
written in that module without first introducing a concrete integer
coefficient, so the blanket allow cannot hide one.

**A cast audit is mostly an exercise in telling indices from values**, and the
ratio here — 1 in 117 — is the argument for doing it by module rather than by
warning count.

### Phase 2: the rest of `src/`, and the gate

The remaining ~250 sites went the same way, and the ratio held. Triaged by
*kind* rather than by module, which is the faster cut: of the whole tree only
**26** casts were between wide integer types at all (`u128 ↔ i128`, `u128 → u64`),
and those are the only ones that can carry a coefficient. Everything else —
`usize → u32`, `i64 → usize`, `u32 → u8` — is shape indices, DP window bounds,
and byte-packing widths.

Five of the 26 needed a check rather than a proof, and each is a value with no
a-priori bound: the two `h̃`-basis structure constants in `character_basis.rs`,
the common denominator and gcd in `gjmod.rs`, and the tableau count crossing
into the signed difference array in `two_row.rs`. The others had real proofs and
now state them — `guard.rs`'s five are each guarded by an explicit
`n > i128::MAX as u128` immediately above, `modular.rs` is residues mod
`p < 2^31` throughout, `gj.rs` is bounded by the `n!` wall that fires at n = 34
first, and `kostka.rs`'s is a count and so non-negative by definition.

`src/` is clean under all three lints. Every module carries a module-level
`#[allow]` with a proof **in that module's own terms** — not a blanket
suppression, which is what R5 forbids: `permutation.rs` says its casts are
positions in a word `Perm` already bounds, `skew_lr.rs` says every element it
packs is at most the cell count `elem_width` sized the buffer from, `two_row.rs`
points at the two checks it does carry.

The gate is **`cargo clippy --lib`, denied in CI**, rather than `deny` in
`[lints.clippy]`. `[lints]` reaches examples and tests too, and those are
research drivers — 94 of the remaining warnings are theirs, all index
arithmetic in harnesses that never ship. Gating them would buy noise; gating
the library buys the invariant.

## The (q,t) walls, measured (policy item 6, R9)

After item 1 the `(q,t)` families fail loudly instead of wrongly, which makes
them honest but not *stated*. R9 asks for the wall in reproducible terms, or
for the admission that it is unmeasured. `examples/probe_qt_walls.rs` measures
it.

**The instrument is the point, and the first version of it was wrong.** Walking
the degree up until something panics answers "did it overflow by n = N?", which
is not a reach statement — it says nothing about whether N+1 is fine, or
whether the family spent the whole walk one degree from the wall. The probe
reports the **bit width of the widest coefficient** at each degree instead, so
the growth per degree is what extrapolates, in the form `memo.rs` already uses
for its own bounds. The projection was then checked against a measured wall and
came within one degree (below).

### Two walls, and the unit of work decides which binds

| entry point | reached | widest | bits/degree | 127 bits near |
|---|---|---|---|---|
| `hall_littlewood_table` | n = 30 | 45 | +2.2 | n ≈ 67 |
| `kostka_foulkes_table` | n = 30 | 45 | +2.2 | n ≈ 67 |
| `hall_littlewood_p_table` | n = 24 | 14 | +0.9 | n ≈ 156 |
| `qt_kostka_table` | n = 18 | 20 | +1.8 | n ≈ 77 |
| `macdonald_j`, every λ ⊢ n | n = 13 | 39 | +4.3 | n ≈ 33 |
| `nabla_e` | n = 16 | 29 | +2.8 | n ≈ 52 |
| `delta_prime_e` | n = 16 | 30 | +2.8 | n ≈ 51 |
| `llt_h_table` (k = 3) | n = 20 | 39 | +2.6 | n ≈ 54 |
| `llt_gtilde_table` (k = 3) | n = 17 | 31 | +2.5 | n ≈ 55 |
| `hall_littlewood`, λ = 1ⁿ | n = 47 | 87 | +2.5 | n ≈ 63 |
| `hall_littlewood`, λ = (n-1,1) | n = 200 | 1 | +0.0 | never |
| `macdonald_j`, λ = 1ⁿ | n = 96 | 23 | +0.3 | n ≈ 499 |
| `llt_h`, μ = 1ⁿ (k = 3) | **overflows at n = 87** | 127 at n = 86 | +1.6 | measured |

Every *whole-degree* row is stopped by runtime with its arithmetic wall two to
four times further out in degree. For those, "exact as far as you can afford to
compute" is the true reach statement and the `i128` wall is a fact about a
degree nobody reaches.

The single-shape rows are a different statement about the same families, and
that distinction is the chapter's main result: `hall_littlewood_table(n)` is
p(n) polynomials and stops finishing near n = 30, while one `Q'_λ` runs to
n = 47 in the same budget and projects to a wall at n ≈ 63 — hours of compute,
not never. **Reach is a property of the entry point, not of the family.**

### The one wall a caller reaches cheaply

`llt_h` at μ = 1ⁿ, k = 3 **overflows at n = 87**, in about 1.5 seconds. Not a
research campaign — a call anyone could make. The widest coefficient is
`87456140265747761930074962686430006450` (127 bits) at n = 86, against
`i128::MAX ≈ 1.7·10³⁸`, and the panic is `attempt to add with overflow` in
`impl_ring_for_int!`'s `add_assign`: an *answer* coefficient accumulating past
the width, not a transient intermediate, which makes it a harder wall than the
`z_γ`-style ones this crate usually meets.

It moves sharply with the ribbon level — k = 2 at n = 124, k = 3 at n = 87,
k = 4 at n = 71 — so a larger `k` packs more coefficient into the same degree.
None of this is visible from `llt_h_table`, which runtime stops at n ≈ 20.

**A projection is worth what its first check says it is worth.** The
extrapolation put this wall at n = 86; the measured wall is 87. That is the
only calibration the other rows' projections have, and it is why they are
quoted as projections.

### The shape dominates the degree

Two rows above make a reach statement phrased in `n` alone unwritable:
`Q'_{1ⁿ}` walls near n ≈ 63 while `Q'_{(n-1,1)}` has 1-bit coefficients at
*every* degree out to the probe's cap. For Hall–Littlewood the extremal shape
is `1ⁿ` and it is extremal for a reason — the coefficients there are the
Kostka–Foulkes ones at content `1ⁿ`, summing to `K_{μ,1ⁿ} = f^μ`, and
`Σ_μ (f^μ)² = n!` caps them at `√(n!)`, the same bound behind the character
ceiling at n ≈ 58. Measured widths track it about 9 bits below, and the
measured slope (+2.5) matches its derivative `½·log₂ n`.

That bound does *not* transfer. Swept over every λ, `macdonald_j` gains 4.3
bits per degree; at λ = 1ⁿ it gains 0.3 and is 23 bits at n = 96. So `1ⁿ` is
nowhere near extremal for `J`, and the single-shape wall there is recorded as
**unmeasured** — which λ maximizes the numerator has not been established.
`hall_littlewood_p` likewise has no bound, `P` coming from inverting the
Kostka–Foulkes matrix, so its coefficients are signed and are not counts.

### Two ways a degree-walking probe lies to itself

Both were found by running it, and both are now comments in the harness.

*A time budget does not stop a family that is cheap at every degree.*
`Q'_{(n-1,1)}` is: the Morris recursion peels one part and the expansion stays
tiny. The first run reached **n = 85 799 505** that way, having long since
stopped measuring anything a caller would ask for. Fixed with `MAX_DEGREE`.

*A lookahead assumes smooth growth, and a vacuous degree breaks it.* The
staircase λ = (n, n-1, …, 1) has size `n(n+1)/2`, which carries no 3-ribbon
tableaux unless 3 divides it — so two degrees in three return zero instantly.
The walk read n = 13 as free (it is vacuous) and started n = 14, which is not,
and sat inside that single call for over an hour before it was killed. The row
is gone: **a row that alternates empty and enormous measures neither.**

Timings above are upper bounds — two probes shared the machine for part of the
run — which affects the degrees reached, not the bit widths, since those are
deterministic.

## The ladder, where the measurement asked for one (policy item 6, second half)

Item 6 said escalation lands per family when a workload demands it. `llt_h` at
μ = 1ⁿ demanded one: n = 87 in about a second is not a research campaign, it is
an ordinary call. Through the Python boundary the single-shape entry points of
both families now escalate — `hall_littlewood`, and LLT's `llt_gtilde`,
`llt_h`, `llt_h_tilde`, `llt_g_lt`, `llt_schur` — with the fixed-width pass
over `Guarded` reporting and the same generic code re-running over `BigInt`.
The table entry points deliberately do not: their walls are 2–4× further out in
degree than anything that finishes.

Both families memoize locally and generically, so this needed no cache work at
all — the two-tier cache above stays unbuilt, and stays right.

One design note worth keeping. The boundary's `Boundary` trait means "one of
the two escalation passes", and the unescalated `(q,t)` entry points still have
to emit coefficients over plain `i128`, which is emphatically not a pass — it
panics at its wall rather than reporting. Making the output helpers generic by
widening `Boundary` to include `i128` would have quietly made `i128` acceptable
everywhere an escalating pass is meant. The outbound half is split out as
`ToCoeff` instead, so the two cannot be confused.

### Macdonald: a documented claim the measurement contradicted

Item 6 left `macdonald_j`'s single-shape wall unmeasured, because λ = 1ⁿ — the
extremal shape for Hall–Littlewood — turned out to be nowhere near extremal
here. Sweeping every λ at each degree settles it: **the extremal shape is the
single row `λ = (n)`**, at every degree from 4 to 12, with 1ⁿ at the other end.

Walking that shape reaches the wall in about a minute per call:

| at λ = (n) | overflows | last clean |
|---|---|---|
| `macdonald_p` | n = 30 | 122 bits, 66 s |
| `macdonald_q` | n = 26 | 101 bits, 63 s |
| `macdonald_j` | n = 26 | 101 bits, 59 s |

101 bits and then overflow one degree later, against a ~5.5 bits/degree slope,
says the wall is in the **intermediates** of the `Frac` arithmetic rather than
in the answers — the shape the record keeps meeting, and the reason the answer
widths understate it.

**Also corrected:** `macdonald_p`'s rustdoc claimed "the enumeration becomes
impractical long before the width does", citing degree 10, where the widest
coefficient is 25 bits and growth is ~3.5 bits per degree. At the extremal
shape that is false — n = 29 finishes in 66 seconds and n = 30 overflows, and
66 seconds is not impractical. The claim was true of the shapes it was measured
on and was generalized past them. All three entry points now escalate, so the
boundary has no wall at all.

**The pin is in `tests/bignum.rs`, not in `python.rs`**, and the reason is the
one the panic audit already found: an `extension-module` test binary has no
interpreter, so CI builds the `python` feature rather than testing it. The
mechanism — the fixed-width pass *reports* at n = 87, and the wide pass returns
a coefficient that provably could not have fitted — is testable without the
boundary, and that is where CI's release lane can see it. Release-only, like
the memory budgets: 9 s there against 78 s in a debug build.

## The rules with no audit item: R4, R6, R10

The seven-item list was organized around the rules with known gaps, so three
rules never got an item. Absence of an item is not evidence, and two of the
three turned out to hold for reasons nothing wrote down.

### R4 — wrapping only where the ring is modular

27 sites in `src/` use `wrapping_*` / `saturating_*` / `overflowing_*`. Most
are `saturating_sub` on a length, the clamp-at-zero idiom, where no value is
at stake. Five carry something that is not an index, and each is legal for a
*different* reason — which is why a blanket "these are fine" would have been
the wrong note to leave:

| site | why it is legal |
|---|---|
| `shard_of` (`skew_lr.rs`) | hash mixer; the product mod 2^64 **is** the operation |
| `overlap`'s row index (`skew_lr.rs`) | `0.wrapping_sub(1)` is the "no row above" sentinel, and `overlap` rejects it |
| `masks.saturating_mul(width)` (`llt.rs`) | a budget test, where saturating routes exactly as the true value would |
| `total_states` / `total_dimension` (`schubert.rs`) | dispatch magnitudes, compared only against each other |
| `schubert_monomial_mass_of` (`schubert.rs`) | same, but **`pub`** — so it now says so under `# Reach` |

The last one is the only one that was arguably a defect: a public function
returning `u128::MAX` in place of a true mass, with nothing in its docs to say
it saturates. It is a cost signal and saturation does not change how it
classifies, so the fix is the sentence rather than a checked multiply.

### R10 — the exact side never shares the width under test

The oracles are clean by construction: `oracle.rs`, `lrcalc_oracle.rs` and
`sage_oracle.rs` compare against committed Sage/lrcalc fixtures, so the exact
side is external and arbitrary-precision. `oracle.rs` runs `i64` against an
`i128` production path, which is the rule working as designed.

`algebra_laws.rs` is the one that reads like a violation: every law computes
**both sides over `i64`**, which is precisely the shape R10 names. It survives
on R10's own escape clause, an a-priori bound — and until now that bound was
assumed. Measured, the widest value anywhere in those sweeps is **120**, and
it is `z_{1⁵} = 5!` from the Hall pairing rather than any structure constant.
Against `i64::MAX` that is 56 bits of headroom.

The bound belongs to the **degree caps** (5, 3 and 2 across the three sweeps),
not to the laws, so raising a cap spends the headroom without touching
anything that looks like a bound.
`the_laws_run_far_below_the_width_both_sides_share` now pins it: raising the
sweep to degree 7 makes it report `5040 (at z at [1,1,1,1,1,1,1])` and name
what to re-derive. That converts the premise from a comment nobody would check
into a test that fails when it expires.

Worth recording alongside it: since item 1 put `overflow-checks` in the
release profile, **no profile wraps a native integer silently**, so R10's
failure mode there is now a panic rather than two sides agreeing on a wrong
answer. The rule is not retired, because `as` casts, `wrapping_*` and
`Guarded` all sit outside that backstop — but for plain `i64`/`i128` laws it
is a second line of defence rather than the only one.

### R6 — audited, and it held one live defect

Every operation inside a `guarded` scope must report, check, or prove. The
audit went in expecting the quiet failure — `Guarded` does not panic, so a
seam that neither reports nor proves returns `Some(garbage)`. What it found
was the opposite, and the inversion is the lesson.

**Since item 1, the quiet failure is mostly gone and a loud one replaced it.**
`overflow-checks` in the release profile means native `+`/`*` no longer wrap in
*any* profile. Inside a guarded scope that is not the improvement it looks
like: the ladder is watching for `None`, and a panic is not something
`escalate` can catch. So the seam that used to return a wrong answer now kills
the process, on an input where the wide pass has the answer sitting right
there. Silent-to-loud is the right trade everywhere else in this policy; here
it converts one bug into another.

**The instance.** `powersum_scalar` ([jack.rs](../../src/jack.rs)) needs z_μ
and formed it with `Partition::z`, which accumulates in native `u128`:

| \|μ\| | z_μ = \|μ\|! | what happened |
|---|---|---|
| ≤ 33 | fits `i128` | fast pass answers |
| 34 | ≈ 2.95·10³⁸ — fits `u128`, past `i128` | `from_u128` **reports** → escalates |
| ≥ 35 | past `u128` | native multiply **panics** → ladder cannot escalate |

The one-degree window at 34 is why this survived: the wall *appears* to be
handled, because at the first degree that exceeds `i128` only the injection is
too wide and the injection is checked (item 3). The bug starts one degree
later, when the native accumulation itself goes.

Reachable from `jack_scalar` and `jack_structure_constant`, both of which wrap
their fast pass in `guarded` — and cheaply, since `powersum_scalar` is exported
from `lib.rs` and the probe hit it in milliseconds. In Sage it surfaces as a
`PanicException`, which R2 says is a bug report and never an interface.

The fix is one line and was already prescribed: `Partition::z`'s own docs name
`z_in` as the escape for a caller multiplying *by* z_λ, because it accumulates
in the coefficient ring — so every factor becomes a `Guarded` multiply that
reports. After it, |μ| ≥ 35 returns `None` and the wide pass carries \|μ\|!
exactly, verified to 40!. Pinned by
`the_z_wall_reports_inside_a_guarded_scope_rather_than_panicking`.

**What else was checked, and by what argument.** The rest holds:

| surface | why it is compliant |
|---|---|
| `Guarded` / `GuardedRat` impls | every op is `checked_*` → `note_overflow`; the `sub_assign` default routes through `neg` + `add_assign`, both reporting, which is where item 2's `i128::MIN` fix earns its keep |
| `eval::dimension`, `principal_specialization`, `character_uncached` | `checked_*` → `None`, R6's case (b) |
| `gj.rs`, `gjmod.rs` | concrete `i128` and modular; never instantiated at `Guarded`, so outside the rule |
| `schubert::dimension`'s saturation | reached only by the documented cost signal and tests, never by a coefficient |
| LR counts (`lr.rs`, `skew_lr.rs`) | one increment per enumerated tableau, so the enumeration walls first |

**And the rule's own instruction was outstanding.** R6 ends "this precondition
belongs in `guarded`'s rustdoc as part of its contract, not only here" — and it
did not. The doc asserted the promise ("a `Some` is a promise that every
intermediate stayed inside the fixed width") with nothing about what the
closure must do to make it true, which reads as unconditional. It now carries
the three obligations, both failure modes, and the instance.

## Open

- **R4, R6 and R10 were assumed rather than audited.** The seven-item list
  covered R1/R2, R3, R5, R8 and R9; the other three rules had no item and no
  evidence. All three are now checked. R4 and R10 hold, for reasons that were
  unstated and are not the same reason twice. R6 did not: it held a live
  defect, now fixed and pinned — see the chapter below. Nothing here is open
  any longer; the entry stays as the record of what "no item" was hiding.

- **Two-tier caches are the answer for memoized intermediates, and nothing
  needs them yet.** The escalation ladder assumes the wide pass can re-run the
  computation; a memo typed at the narrow width breaks that, because the wide
  pass reaches the same cache. `htilde_table::<C>` is the instance —
  it computes `htilde_table_uncached::<i128>` whatever `C` is, so
  `qt_kostka_table::<BigInt>` walls exactly where `<i128>` does, and a Python
  layer dispatching to a "bignum backend" there would buy nothing. The
  mechanism is written up in
  [../policies/failure.md](../policies/failure.md) ("Two-tier caches"): a
  second static per coefficient regime, the wide tier seeded by widening
  injection from the narrow one so only the entries that overflowed are
  recomputed, and R7's peek/store split on the fill side. The portable half is
  a rule about keys — **cache the unit that overflows, not the unit that is
  asked for**: `bh_pieri_table` is keyed per `(μ,ν)` and would degrade entry by
  entry, while `htilde_table` is keyed by degree and would redo a whole degree
  for one large value. It stays unbuilt because the measurement says so: the
  families blocked on the cache cannot reach their arithmetic wall, and the
  family that can reach one (`hall_littlewood` at λ = 1ⁿ) memoizes locally and
  generically, so it needs no cache work at all.
- **The release lane exists but has never run.** `.github/workflows/ci.yml`
  now carries one, because the canary and the escalation pin only carry
  information under `--release`. This repository has no remote, so the workflow
  is written and unverified; the first push is where its portability claims get
  tested. Making `-D warnings` affordable took clearing nine warnings, four of
  which were functions that are exercised by tests and dead only outside them.
