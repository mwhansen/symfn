# Executing the failure policy

[docs/policies/failure.md](../policies/failure.md) is the rulebook — every way
a computation here may fail, and which mechanism each situation demands. This
file is the other half: what changed in the tree to meet it, what the changes
cost, and what turned up along the way. The policy's "What this changes" list
is the agenda; each item lands here as a chapter when it is executed.

**State.** Item 1 is done: `[profile.release]` carries
`overflow-checks = true`, measured at 0–6% on every harness the crate has and
0% on the bignum routes, with a four-test canary
(`tests/overflow_checks.rs`) that fails the day the line leaves `Cargo.toml`.
Items 2–7 are open; the policy file holds them in execution order.

Every number below is from one machine — macOS arm64, rustc 1.96, on AC —
which is the standing caveat until CI exists
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

## Open

- **The (q,t) walls are now loud but still unmeasured** (policy items 1 and 6).
  Hall–Littlewood, Kostka–Foulkes, Macdonald, qt-Kostka, nabla/delta and LLT
  instantiate at plain `<i128>`; the flag turns their walls from wrong answers
  into panics, which is honest but not yet *stated*. The arithmetic bound is
  that the coefficients of one `H̃_μ` sum to `n!` at `q = t = 1`, and
  `34! ≈ 3·10³⁸` already exceeds `i128::MAX ≈ 1.7·10³⁸`, so a single-coefficient
  wall is at most a few degrees past n = 34. Where each family's wall actually
  sits is an R9 gap.
- **CI has no release-profile lane** (release-readiness Phase 0). The canary
  and the escalation pin only carry information under `--release`; a CI that
  runs the default profile alone will report them green while the flag is gone.
