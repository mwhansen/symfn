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
reach stated in reproducible terms (R9).

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
line from Cargo.toml. The flag is gated on measurement (item 1 below), but
the burden of proof sits on *off*, not on: a hot loop the flag visibly slows
gets explicitly checked or proven arithmetic, not the flag removed.

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
only here.

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

### R9 — Every fixed-width public family states its reach

Per [style.md](../style.md), "Reach and performance in rustdoc": the wall, in
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

R2 lets a violated precondition panic because a Rust caller who broke one has
a bug in the same tree. Neither half of that holds across the FFI: a Python
caller's list is *data*, often built from a file or a Sage object, and a panic
reaches them as a `PanicException` that R2 already calls a bug report and
never an interface. So every precondition a Python caller can violate —
non-partitions, non-permutations, unknown basis strings, malformed graphs,
indices past a representation's ceiling — is checked at the entry point and
raised as a typed exception naming the requirement. `part_arg`, `perm_arg`,
`level_arg` and `variable_arg` in [python.rs](../../src/python.rs) are the
models, and `scripts/check_python_boundary.py` is the pin.

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
| violated precondition | panic naming the requirement; `# Panics` section | the abacus assert in [partition.rs](../../src/partition.rs) |
| the same, reachable from Python | validate at the entry point; `PyValueError` naming the requirement (R11) | `part_arg` / `perm_arg` / `level_arg` in [python.rs](../../src/python.rs) |
| capacity wall reachable from Python | the owning module exposes the bound; the entry point refuses on it (R11) | `abacus_arg` against `llt::abacus_reach` |
| narrowing conversion | `try_from` with loud failure, or a bound proof | R5 |
| everything unforeseen | the profile backstop — a net, never an interface | R3 and its canary |

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
  ([lib.rs](../../src/lib.rs) rustdoc — pending, item 7 below);
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

1. **Measure `overflow-checks = true`, then flip it (R3).** Today
   [Cargo.toml](../../Cargo.toml) has no `[profile.release]`, so the profile
   users actually ship wraps: `impl Ring for i64/i128` multiplies with a
   plain `*` ([coeff.rs](../../src/coeff.rs)), and the Hall–Littlewood,
   Kostka–Foulkes, Macdonald, qt-Kostka, nabla/delta, and LLT pyfunctions
   instantiate at plain `<i128>` with no escalation — confident nonsense past
   their walls. Those walls sit inside advertised territory: the coefficients
   of one `H̃_μ` sum to `n!` at `q = t = 1`, and `34! ≈ 3·10³⁸` already
   exceeds `i128::MAX ≈ 1.7·10³⁸`, so the single-coefficient walls are at
   most a few degrees past n = 34 — unmeasured, which is itself the R9 gap.
   Gate: run the harnesses with the flag on (`bench_guarded`, the skew
   products, `bench_kron_coeff`, the memory harness to rule out allocation
   confounds), and land the numbers in the record. Then: add the R3 canary,
   and upgrade `unguarded_fixed_width_is_wrong_where_the_guarded_path_escalates`
   ([ops.rs](../../src/ops.rs)) from `#[ignore]` documentation to a
   release-mode CI pin — its premise inverts from "release wraps" to
   "release panics". CI grows a release-profile lane for exactly these tests,
   because the flag only exists in that profile.
2. **Close the guard's own `i128::MIN` corners (R5).** `Guarded::neg` and
   `GuardedRat::neg` use `wrapping_neg`, and both `gcd`s take `.abs()` of
   possibly-`MIN` values ([guard.rs](../../src/guard.rs),
   [coeff.rs](../../src/coeff.rs)); `checked_mul` can legitimately return
   `i128::MIN`, after which negation wraps with no report — `Some(garbage)`
   from a measure-zero input. `checked_neg`/`checked_abs` routing to
   `note_overflow` (guard) or an assert (`Rational`), with a `MIN`-injection
   regression test.
3. **Make fixed-width injections loud (R8).** `from_u128`/`from_i128` on
   `i64`/`i128` become checked, with tests pinning the panic; `Guarded` keeps
   reporting-then-escalating.
4. **Audit the panic sites against R1/R2.**
   [release-readiness.md](../release-readiness.md) Phase 3 counts 138
   `panic!`/`unwrap`/`expect` sites in `src/` (re-grep at audit time). Each
   ends as a documented contract violation, a documented wall, or a
   `Result`/`Option`. This absorbs that phase's first two checklist items.
   The **Python-reachable** subset is done under R11 — five clusters over ~30
   entry points, every `unwrap` in [python.rs](../../src/python.rs) removed,
   pinned by `scripts/check_python_boundary.py`
   ([python-and-sage-interop.md](../record/python-and-sage-interop.md)). What
   The entry points that returned a plausible `0` are also resolved: five
   whose zero was a convention over an undefined question now raise, and the
   rest are theorems and say so. What remains is the Rust-facing sites.
5. **Cast audit, phased (R5).** Roughly 600 `as` sites. Enable
   `clippy::cast_possible_truncation`, `cast_sign_loss`, and
   `cast_possible_wrap` as warnings in `[lints]`; audit coefficient-adjacent
   modules first ([convert.rs](../../src/convert.rs),
   [eval.rs](../../src/eval.rs), [character.rs](../../src/character.rs), the
   value paths of [llt.rs](../../src/llt.rs)) rather than bounded `u32`
   partition bookkeeping; justified sites get `#[allow]` plus the proof;
   flip to deny in CI when clean.
6. **(q,t) escalation on demand; reach documented now (R9).** After item 1
   the (q,t) walls fail loudly, which makes them honest; real escalation
   (`QtPoly` over `Guarded`, rerun over `BigInt`) lands per family when a
   workload demands it — mechanism follows a measured need here, as
   everywhere else in this crate. Until then, each family's module doc states
   its wall, or states that it is unmeasured.
7. **Promote the caller-facing contract into `lib.rs`,** closing the
   release-readiness Phase 3 exit criterion ("the overflow contract is on the
   type, not only in the README") — compressed to the contract plus a pointer
   here, per the promote stage of [style.md](../style.md), "How a learning
   ages".
