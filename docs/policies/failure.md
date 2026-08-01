# Failure policy — how this library is allowed to fail

What this file governs: every way a computation in this crate may fail —
overflow, capacity walls, violated preconditions, panics — and which mechanism
each situation demands. What it does not govern: prose and documentation style
([style.md](../style.md)), or performance policy — measurements and their
context belong to the record.

Like the style guide, this is a rulebook, not a description of current
practice. Where the two disagree, the gap is listed in
[What this changes](#what-this-changes) at the bottom, with the reason — and,
because the gaps here are correctness gaps, in execution order.

## The invariant

**Every value that leaves this library is exact, or the call fails loudly — in
every build profile.** No code path may return a wrapped, truncated, or
rounded result.

The invariant is stated this hard because trust is the product. The
researcher this crate serves publishes results under their own name, and the
failure that destroys them is not the crash — it is the plausible wrong value
([style.md](../style.md), "The readers"). A crash costs a session; a wrapped
coefficient costs a retraction.

A computation has exactly **three legal outcomes**:

1. the exact answer;
2. escalation to a wider ring, then the exact answer;
3. a **loud refusal** — `None`, `Err`, a typed Python exception, or a
   documented panic. Loud means the caller cannot mistake it for an answer.

There is no fourth outcome. "Right below some unstated degree", "wrapped but
rarely", and "wrong with a warning printed" are all the fourth outcome.

## The rules

### R1 — Overflow is capacity, not a bug, and never silent

A coefficient outgrowing a fixed width is an ordinary event on a
caller-reachable path, so it is never a wrapped value and never (in an
escalating path) a panic. Entry points that promise exactness escalate
(`escalate` in [python.rs](../../src/python.rs), `kronecker_coeff` in
[ops.rs](../../src/ops.rs)); paths that cannot escalate return
`Option`/`Result` (`try_character`; the numerator products in
[eval.rs](../../src/eval.rs)) or document the wall under `# Panics` with
range stated in reproducible terms (R9).

### R2 — Panics are for violated contracts and proven-unreachable states

"This partition is not a partition" panics: the caller broke a documented
precondition, and the message names the violated requirement in the problem's
terms — "an abacus needs at least ℓ(λ) beads to hold λ"
([partition.rs](../../src/partition.rs)) is the model. A state reachable by a
caller who followed the docs must not panic, except at a wall the docs state.
A `PanicException` surfacing in Sage is by definition a bug report — the
crate's bug, or evidence that a wall needs a real mechanism — never an
interface.

### R3 — The release profile carries `overflow-checks = true`

The escalation ladder handles the overflow that was foreseen; the profile
flag converts the overflow that was not from silent to loud. The profile is
part of the correctness surface, so a release-mode canary test pins the flag:
a `#[should_panic]` overflow that fails the suite the day someone drops the
line from Cargo.toml (`tests/overflow_checks.rs`, verified against a build
with the flag off). The burden of proof sits on *off*, not on — measured at
0–6% across every harness the crate has
([failure-and-overflow.md](../record/failure-and-overflow.md)) — so a hot loop
the flag visibly slows gets explicitly checked or proven arithmetic, not the
flag removed.

`overflow-checks` reaches neither `as` casts (R5) nor `wrapping_*` calls
(R4); those have their own rules.

### R4 — Wrapping arithmetic only where the ring is modular by definition

Hashing ([fasthash.rs](../../src/fasthash.rs), the key mix in
[skew_lr.rs](../../src/skew_lr.rs)) and `ℤ/p`
([modular.rs](../../src/modular.rs)) wrap because their rings do; that is
exactness, not overflow. Such sites are spelled `wrapping_*`, so one grep
finds every place modular arithmetic is claimed. Bare arithmetic relying on
release-mode wrap is never legal. A modular engine used to produce an answer
in ℤ pairs with an a-priori bound or a cross-check — `engines_agree` in
[gjmod.rs](../../src/gjmod.rs) is the shape.

### R5 — Narrowing conversions prove themselves or check themselves

A widening `as` is free. A narrowing or sign-changing conversion carries
either a one-line bound proof at the site or `try_from` with a loud failure.
The record already knows why: "narrowing is where overflow bugs live"
([littlewood-richardson.md](../record/littlewood-richardson.md)). Inside the
guarded types, *everything* is checked, including `neg` and `abs`:
`checked_mul` can legitimately produce `i128::MIN`, after which a
`wrapping_neg` is a silently wrong sign (item 2 below).

### R6 — Inside a `guarded` scope, every operation reports, checks, or proves

`guarded(|| …) -> Some(v)` is a promise that every intermediate stayed inside
the fixed width ([guard.rs](../../src/guard.rs)). The promise is only as good
as the closure: raw arithmetic that wraps without `note_overflow` yields
`Some(garbage)`. So inside a guarded scope, every arithmetic operation either
(a) runs through `Guarded`/`GuardedRat`, (b) is `checked_*` with an explicit
refusal or fallback on `None` — `integral_sweep` in
[convert.rs](../../src/convert.rs) keeps its own flag and bails to the
generic path, which is the model — or (c) carries a bound proof at the site.
This precondition belongs in `guarded`'s rustdoc as part of its contract, not
only here — **it is now there**, with both failure modes and the instance.

⚠️ The second failure mode is the one that reads as safe. Since R3 put
`overflow-checks` in the release profile, native arithmetic no longer wraps in
any profile — but a *panic* inside a fast pass is not an improvement over a
wrong answer, because `escalate` is watching for `None` and cannot catch it.
The caller gets a crash on an input the wide pass answers exactly.
`powersum_scalar` formed z_μ in native `u128` and did this past |μ| = 34; the
fix is `Partition::z_in`, which accumulates in the coefficient ring so every
factor reports ([failure-and-overflow.md](../record/failure-and-overflow.md)).

### R7 — Caches never launder overflow

A value computed by a call that overflowed is garbage; caching it converts a
detected overflow into an undetected one, because a later reader sees a clean
counter and accepts the poisoned entry. So: overflow reports are never cached
as values, values are stored only after their computation is confirmed clean,
and where fill and check cannot be one operation, peek and store split —
`bold_p_peek`/`bold_p_store` in [memo.rs](../../src/memo.rs) is the model,
with the reasoning in its doc.

### R8 — Fixed-width injections fail loudly on truncation

Structure constants — LR coefficients, Kostka numbers, `z_λ`, characters —
are exactly where large values enter, so the injection seams
`Ring::from_u128`/`from_i128` must not truncate. On `i64`/`i128` they check
and panic naming the constant and the ring; injection is never a hot loop, so
the check costs nothing that matters. `Guarded` instead reports and
escalates, which is its job.

### R9 — Every fixed-width public family states its range

Per [style.md](../style.md), "Range and performance in rustdoc": the wall, in
reproducible terms — "the wall is `z_γ` in the intermediates at total degree
24, a fact about `i128`" is the house form — with the measurement in the
record file that owns it. A wall not yet measured is stated as unmeasured. An
unstated wall is how a refusal the design intended becomes a surprise the
caller publishes past.

### R10 — The exact side of a comparison never shares the width under test

Oracles, cross-checks, and both sides of an `engines_agree`-style pairing run
over arbitrary precision, or over a different width with an a-priori bound.
Two sides sharing one fixed width wrap identically and agree on the same
wrong answer. `bench_kron_coeff` runs `BigRational` on both Kronecker routes
for exactly this reason (its manifest note in
[Cargo.toml](../../Cargo.toml) says so).

### R11 — A foreign caller's argument is data, and is validated, not repaired

*That* a Python-reachable precondition raises rather than panics is
[python.md](python.md)'s P8, which owns the Python surface and states the bar.
This rule is the failure-mechanism half P8 defers here: **which** wrong answer
the boundary is guarding against, and how to tell it from a right one. R2 lets
a violated precondition panic because a Rust caller who broke one has a bug in
the same tree; neither half holds across the FFI, where the argument is data
and the panic is unhandleable. `part_arg`, `perm_arg`, `level_arg` and
`variable_arg` in [python.rs](../../src/python.rs) are the models, and
`scripts/check_python_boundary.py` is the pin.

**Repairing the input is not the alternative.** `Partition::new` normalizes —
it sorts and drops zeros — which is right for a Rust caller who built the
vector from a generator and wrong for a foreign caller whose list is data:
`[1, 3]` used to reach the mathematics as `[3, 1]`, and the caller got a
well-formed answer to a question they had not asked. That is the plausible
wrong value this policy ranks below a crash, arriving through the door built
to be convenient. Padding is the one exception, because it carries no other
reading: Sage hands over fixed-width lists, so trailing zeros are dropped and
everything else raises.

**A zero is not automatically a refusal in disguise.** Where the value is a
theorem — `c^λ_{μν} = 0` off-degree, `K_{λμ} = 0` when λ does not dominate μ,
`s_λ(1^n) = 0` when `ℓ(λ) > n` — it is an answer, and a caller sweeping a range
depends on getting it; those say so in rustdoc and stay. Where the zero is a
*convention* standing in for an object that has no value at all — `χ^λ(μ)` with
`|λ| ≠ |μ|`, `g^ν_{λμ}` across degrees — the boundary raises. The core keeps
the convention, because totality is what a composing Rust caller needs and the
routes to one coefficient and to a whole product must agree; only the boundary
is stricter. Distinguishing the two is a mathematical judgment per function,
not a rule that can be applied by grep.

A **capacity** wall on this boundary is the same story with a different cause
— the caller violated nothing, the representation ran out — and gets the same
treatment for the same reason. Where the bound lives in another module, that
module exposes it (`llt::abacus_reach`, `llt::MAX_CELLS`) rather than the
boundary restating it, since a bound copied to a second site is a bound that
drifts.

## Choosing a mechanism

### Two questions before any mechanism

**Can it overflow at all, for inputs the contract accepts?** If a bound
proves it cannot, raw arithmetic plus the one-line proof at the site is the
correct technique — most partition bookkeeping (`u32` parts, lengths bounded
by what fits in memory) and sites like the `i128` sweep in
[convert.rs](../../src/convert.rs) ("for every input it accepts, i128 cannot
overflow"). The proof is stated where it is relied on, and R3 is its safety
net: a proof makes the fast path free, the backstop makes a broken proof loud
instead of silent.

**Can the numbers be made smaller?** Algebraic restructuring moves the wall
for free and is usually also faster, and the record shows it repeatedly:
cancel the gcd against the numerator before growing the denominator
(`div_u128` in [coeff.rs](../../src/coeff.rs) — the divisor is `z_μ`, which
reaches `|μ|!`); keep intermediates near the answer instead of near `|λ|!`
([eval.rs](../../src/eval.rs)); route through the basis whose intermediates
stay small ([kronecker.md](../record/kronecker.md) — the overflow was
entirely in the intermediate rationals). Only after the numbers are as small
as the mathematics allows does a failure mechanism get chosen for what
remains.

### The table

| situation | mechanism | in-tree model |
|---|---|---|
| boundary entry point; exactness at any size; common case fits `i128` | two-pass escalation: `guarded` fast pass, then the same generic code over `BigInt`/`BigRational` | `escalate` in [python.rs](../../src/python.rs); `kronecker_coeff` |
| fixed width fails on most of the intended range | honest absence: the function exists only under `bignum` | `kronecker_coeff`, wall at n ≈ 26 (its doc), gated rather than panicking |
| hot loop accumulating unsigned counts | width-retry: checked ops in the narrow type, `None` → rerun wider; the widest rung refuses loudly | `Acc` (u64 → u128) in [skew_lr.rs](../../src/skew_lr.rs); `expect("LR tableau multiplicity exceeded u128")` |
| single value that may not fit the signature | `try_*` returning `Option`, beside a ring-generic infallible form | `try_character` / `character_in` |
| local, non-generic calculation | plain `checked_*` + `Option`; no global counter | the `u128` numerator products in [eval.rs](../../src/eval.rs) |
| ring modular by definition | `wrapping_*`, spelled | [fasthash.rs](../../src/fasthash.rs), [modular.rs](../../src/modular.rs) |
| cold path; oracle or verification code | arbitrary precision from the start | `bench_kron_coeff` |
| **memoized** intermediate that can outgrow the width while the answers built from it fit | two-tier cache: the fixed-width table, plus a wide one that seeds from it by widening injection and computes only the entries that overflowed | none yet — see below for where it would apply |
| violated precondition | panic naming the requirement; `# Panics` section | the abacus assert in [partition.rs](../../src/partition.rs) |
| the same, reachable from Python | validate at the entry point; `PyValueError` naming the requirement (R11) | `part_arg` / `perm_arg` / `level_arg` in [python.rs](../../src/python.rs) |
| capacity wall reachable from Python | the owning module exposes the bound; the entry point refuses on it (R11) | `abacus_arg` against `llt::abacus_reach` |
| narrowing conversion | `try_from` with loud failure, or a bound proof | R5 |
| everything unforeseen | the profile backstop — a net, never an interface | R3 and its canary |

### Two-tier caches, when what overflows is a memoized intermediate

The escalation ladder assumes the wide pass can *re-run* the computation. A
memo breaks that assumption in one specific way, and it is the crate's most
familiar shape seen from a new angle: the answers fit and the
**intermediates** do not ([kronecker.md](../record/kronecker.md) — the
overflow was entirely in the intermediate rationals; the `st` basis, where the
answers are under 20 bits and `z_γ` is not). When such an intermediate is
memoized, the cache is keyed on a subproblem whose value can leave the width
even though everything built from it fits — and a cache typed at the narrow
width then walls the wide pass too, because the wide pass reaches the same
cache. `htilde_cached` is the in-tree instance: `htilde_table::<C>` computes
`htilde_table_uncached::<i128>` whatever `C` is, so a `BigInt` instantiation
does not escape the `i128` wall at all.

The mechanism, when a workload needs it:

- **Two statics, not one generic cache.** A `static` cannot be generic, but
  the crate ships exactly two coefficient regimes, so the answer is one more
  instantiation rather than a type-keyed registry.
- **The wide tier seeds from the narrow one.** Widening is exact, so every
  entry the fixed-width cache already holds is a free, correct wide entry.
  Only the entries that actually overflowed get computed wide — which is the
  point: past the wall a few values are large and most are not.
- **Cache the unit that overflows, not the unit that is asked for.** This is
  what decides how much the seeding buys. `bh_pieri_table` and `bh_ell_table`
  are keyed per `(μ,ν)` pair and would degrade entry by entry; `htilde_table`
  is keyed by *degree* and holds a whole table, so one overflowing entry costs
  a wide recomputation of the entire degree.
- **R7 still binds, and binds harder.** Today an entry is trustworthy by
  accident: an overflow panics before the `insert`. A fast tier that *reports*
  instead (`Guarded`) can reach the store with a poisoned value, so fill and
  check split — `bold_p_peek`/`bold_p_store` in [memo.rs](../../src/memo.rs) is
  the model, and its doc already carries the reasoning.
- **Not "just use bignum for the cache".** The narrow tier is the common case
  and carries the residency the memory budgets are calibrated against
  ([memory.md](../record/memory.md)).

Like every other mechanism here, it lands when a measured workload demands it
and not before. The `H̃` family is the one blocked on it today and is
explicitly *not* the case that justifies it: its coefficients gain ~1.3
bits/degree against a runtime wall roughly seven times sooner than the
arithmetic one ([failure-and-overflow.md](../record/failure-and-overflow.md)),
so nothing can reach the wall the cache would move.

### The distinctions that get miscalled

- **Counter vs `Option`.** The global counter exists for one reason:
  `Ring::mul` returns `Self`, so generic code has no per-operation `None`
  channel. Use `Guarded` and the counter only where the trait forces it.
  Local concrete code uses `checked_*` with an explicit `Option`, which keeps
  the failure in the signature, needs no global state, and cannot be tripped
  by an unrelated thread.
- **Escalate vs absent.** Both `character` (ceiling `|λ| ≳ 58`,
  [memo.rs](../../src/memo.rs)) and `kronecker_coeff` (n ≈ 26) run the
  ladder, but only the latter is feature-gated: when the fixed rung fails on
  most of the intended range, a build that cannot escalate would ship a
  function that refuses almost everything it exists for, and the honest form
  of that is an absent function ([ops.rs](../../src/ops.rs), its doc).
- **Answer vs intermediates.** An answer that does not fit needs a wider
  return channel (the Python `int` path) or an `Option`; intermediates that
  do not fit ask the shrink question first, and are otherwise the cheap case
  for two-pass escalation — the re-run is rare and the answer itself is
  small.
- **The wide pass is the same code.** Escalation and retry rerun *the same
  generic code* over a wider type — the coefficient ring is a parameter, and
  that is the whole point ([coeff.rs](../../src/coeff.rs)). A fallback that
  is a second algorithm is an untested path exercised only on the inputs
  least understood.
- **Panic vs refusal.** Reachability under the documented contract decides:
  unreachable-unless-the-caller-lied panics (R2); reachable refuses (R1).
- **Panic vs raise.** *Who* the caller is decides, not what they did wrong. The
  same violated precondition panics for a Rust caller (R2) and raises for a
  Python one (R11), because a panic is a bug report to someone who can fix the
  bug and an unhandleable abort to someone who cannot.
- **Validate vs normalize.** A constructor that repairs its input is a
  convenience for a caller who built that input and a trap for one who was
  handed it. `Partition::new` and `Partition::try_new` exist as a pair for this
  reason; the boundary takes the second (R11).

### Defaults when unsure

New public entry point → two-pass escalation (0–1% on the fast path,
`examples/bench_guarded.rs`). New internal helper → `checked_*` + `Option`,
escalate at the caller. New hot loop → prove the bound or retry-wider, and
measure before inventing anything. New cast → `try_from` unless obviously
widening. A site that fits no row of the table is a policy gap: extend this
file in the same change, rather than improvising silently.

## Where the policy surfaces

This rulebook is for maintainers. Callers meet it as:

- the overflow contract on the crate front page
  ([lib.rs](../../src/lib.rs) rustdoc, "The overflow contract");
- `# Panics` / `# Errors` on every public function that can
  ([style.md](../style.md), the pre-ship checklist);
- the scope precondition in `guarded`'s rustdoc (R6);
- a Reach statement per family (R9), with the wall's measurement owned by the
  subsystem's record file.

## What this changes

Current practice already embodies most of this policy — the escalation
ladder, the width-retry, peek/store caching, and honest absence all exist and
are kept as-is. These are the deltas, in execution order; each names its
gate.

1. ~~**Measure `overflow-checks = true`, then flip it (R3).**~~ **Done** —
   [Cargo.toml](../../Cargo.toml) carries `[profile.release]
   overflow-checks = true`, `tests/overflow_checks.rs` is the canary, and the
   `ops.rs` escalation test now pins the refusal rather than the wrapped
   answer it used to assert. Measurements, and the live-bytes underflow the
   flag found in the allocation harness, are in
   [failure-and-overflow.md](../record/failure-and-overflow.md). Two things it
   did not close: the (q,t) walls are loud but still unstated (item 6), and CI
   has no release-profile lane, so the canary and the escalation pin carry
   information only when the suite is run with `--release`.
2. ~~**Close the guard's own `i128::MIN` corners (R5).**~~ **Done** — the
   guard reports through `checked_neg` and a `u128` `gcd`, `Rational` asserts
   with a message naming the requirement, and `Rational::div_u128`'s `z_μ`
   narrowing is `try_from` rather than `as`. Chapter in
   [failure-and-overflow.md](../record/failure-and-overflow.md).
3. ~~**Make fixed-width injections loud (R8).**~~ **Done** — `from_u128` /
   `from_i128` on `i64`, `i128`, `Rational` **and the trait's own defaults**
   check and panic naming the constant and the ring; measured at ≤1.03x on
   `bench_ops`. Chapter in
   [failure-and-overflow.md](../record/failure-and-overflow.md).
4. ~~**Audit the panic sites against R1/R2.**~~ **Done** — 93 sites at audit
   time (release-readiness Phase 3 counted 138 earlier), each ending as a
   documented contract violation, a documented wall, or a `Result`/`Option`.
   This absorbs that phase's first two checklist items.

   The **Python-reachable** subset closed under R11: five clusters over ~30
   entry points, every `unwrap` in [python.rs](../../src/python.rs) removed,
   pinned by `scripts/check_python_boundary.py`
   ([python-and-sage-interop.md](../record/python-and-sage-interop.md)). The
   entry points that returned a plausible `0` are resolved with it — five whose
   zero was a convention over an undefined question now raise, and the rest are
   theorems and say so.

   The live defect the audit found was the boundary panicking on a malformed
   permutation, reachable only along the escalation path; fixed structurally
   with a `Wide` trait and up-front validation, so that path has no `unwrap`
   left to make. Chapter in
   [failure-and-overflow.md](../record/failure-and-overflow.md).

   The `# Panics` sweep across the whole public surface
   ([style.md](../style.md), delta 2) followed and is **done**: 46 sections
   added, every `pub fn` outside `python.rs` that can panic now says so, and
   every bare `.unwrap()` in `src/` is gone. It found one live defect of its
   own — `Perm::at` returned `w(0) = 0`, its 1-based precondition being a
   `debug_assert` — which item 1 had already converted from a wrong answer
   into a panic by the time it landed. Chapter in
   [failure-and-overflow.md](../record/failure-and-overflow.md).
5. ~~**Cast audit, phased (R5).**~~ **Done** — `src/` is clean under all three
   cast lints, and CI gates the **library** at deny
   (`cargo clippy --lib -- -D clippy::cast_*`); examples and tests stay
   advisory, being drivers whose index arithmetic never ships. Of ~370 sites,
   6 carried a value rather than an index and now check; the rest carry a bound
   proof in their own module's terms. Chapters in
   [failure-and-overflow.md](../record/failure-and-overflow.md).
6. **(q,t) escalation on demand; range documented now (R9).** *Range done*:
   every family's module doc now states its wall, measured rather than
   asserted — `examples/probe_qt_walls.rs` reports the widest coefficient per
   degree, and the growth per degree is what extrapolates. Table in
   [failure-and-overflow.md](../record/failure-and-overflow.md). Two results
   shape what is left: every *whole-degree* entry point is stopped by runtime
   with its arithmetic wall two to four times further out, while
   [`llt_h`](../../src/llt.rs) at μ = 1ⁿ **overflows at n = 87 in about a
   second** — so escalation was worth building for the single-shape LLT and
   Hall–Littlewood entry points and for nothing else yet. *That ladder is now
   built*: those six pyfunctions run `QtPoly` over `Guarded` and re-run over
   `BigInt`, pinned in `tests/bignum.rs` where CI's release lane can see it.
   Both families memoize locally and generically, so it needed no cache work;
   the `H̃` family would need the two-tier cache above, for a wall nothing can
   reach. Macdonald followed the same rule and the same
   evidence: its extremal shape is the single row `λ = (n)`, whose wall is
   n = 26–30 in about a minute, so `P`, `Q` and `J` escalate too.
7. ~~**Promote the caller-facing contract into `lib.rs`.**~~ **Done** — the
   crate front page carries the three legal outcomes, what each coefficient
   type does at its wall, and a pointer here. The stale "Roadmap (the marked
   seams)" it replaced listed four features that had all shipped, which is
   [style.md](../style.md)'s own exhibit for why the reference must not carry
   future work.
