# Failure paths: the panic-site audit

What this file records: the audit of every way this crate can die, run against
[../policies/failure.md](../policies/failure.md) item 4 — what was counted,
what each class turned into, and the two things the count itself got wrong.
The rulebook is that policy file; this is the ledger of executing it.

The Python-reachable half was done first and lives in
[python-and-sage-interop.md](python-and-sage-interop.md): five clusters over
~30 entry points, every `unwrap` in `src/python.rs` removed, pinned by
`scripts/check_python_boundary.py`. What follows is the Rust-facing
remainder.

## The count, and why the old one was not comparable

The policy quoted **138 `panic!`/`unwrap`/`expect` sites in `src/`** and told
its executor to re-grep at audit time. That was the right instruction, because
138 was not a like-for-like figure: it counted `#[cfg(test)]` modules,
`python.rs`, and the assert family all together. Separated, at audit time:

```text
                              non-test, outside python.rs
  panic! / unwrap / expect / unreachable          66
  assert! / assert_eq! / debug_assert             79
```

⚠️ **Do not compare a later count against 138.** The comparable baseline is
66 + 79, and the two families are worth counting apart: an `assert!` on a
documented precondition is the policy working (R2), while a bare `.unwrap()`
is usually a place nobody decided anything.

After the audit:

```text
                              before   after
  panic-family                    66      38
    of which bare .unwrap()      ~30       0
  assert-family                   79      81
  pub fn that can panic,
    with no `# Panics` section    46       0
```

The assert family grew by two because four `debug_assert`s on **public**
preconditions became unconditional, and three new 1-based checks were added
where the wrong answer had been silent.

## One function in the crate documented its panics

That was the headline finding, and it is the R2/style.md gap rather than a
correctness one: `character::character` had a `# Panics` section and nothing
else did. Eight more mentioned panicking in prose without the section
(`Rational::new`, `Atom::unit`, `Frac::inv_factor`, `Frac::ratio`, `Md::new`,
`Md::inv`, `kronecker_coeff`, `Ht::mul`); 37 said nothing at all.

Sections **point at the exposed bound rather than restating it** —
`MAX_SUPPORT`, `MAX_CELLS`, `MAX_FREE_EDGES`, `abacus_reach` against
`ABACUS_REACH_LIMIT` — on the same reasoning R11 gives for the Python
boundary: a bound copied to a second site is a bound that drifts.

## What the sites actually were

The `.unwrap()`s split three ways, and only the third was a mechanism
question at all.

**A redundant test standing next to the `Option` that answers it.** `charge`
checked `word.is_empty()` and then unwrapped `word.iter().max()`;
`divide_by_factor` did the same with its degree bound; the `s → s̃` pivot loop
ran `while !rem.is_empty()` and then unwrapped both the `max_by_key` and the
`remove`. Each asks one question twice. A `let … else` or `while let` asks it
once and the unwrap has nowhere to live — these were *deleted*, not
documented.

**Proven-unreachable invariants with messages a reader cannot act on.**
`.expect("non-empty")`, `.expect("pivot present")`, `.expect("just
inserted")`. R2 licenses the panic and asks the message to name the violated
requirement in the problem's terms; these needed only that.

**Walls**, and this is where R9 bit. Several were undocumented, and one was
mis-classified: `stanley`'s `steps < 1 << 24` reads like a capacity limit and
is not one — it is a backstop against a non-terminating transition recursion,
which without it would hang rather than fail. ⚠️ A first draft of its
`# Panics` section claimed the cap was unreachable because the tree is
bounded by `ℓ(w)` under `MAX_SUPPORT`. **That is not established** — the
transition tree is not bounded by `ℓ(w)` — and the section now states the
reach as unmeasured, which is what R9 requires of a wall nobody has measured.

## Two mechanism changes

### Cache locks recover from poisoning instead of unwrapping

All 22 `RwLock` accesses in `src/memo.rs` were `.unwrap()`. A poisoned lock is
neither a violated contract nor a proven-unreachable state, so R2 never
licensed it; worse, the resulting panic is a *different* failure from the one
that poisoned the lock, it repeats at every later cache access for the rest of
the process, and across the FFI it is unhandleable.

`read_table`/`write_table` recover via `into_inner()`. Sound because a table
is a pure function of its keys and is only ever written whole: `lookup`
releases the read guard before `compute` runs and inserts only after it
returns, and `bold_p_store` is called only on a value its caller has confirmed
clean. So no half-built entry is reachable and recovery cannot launder an
overflow into the cache — which is R7's concern, and the reason `bold_p` has
the peek/store split in the first place.

### `Perm::at` returned a wrong answer in release, and the fix is free

`at`'s "positions are 1-based" precondition was a `debug_assert`. In release
`at(0)` computed `i - 1` as `u32::MAX`, the bounds check missed, and the
`None` arm handed back `i` — so `w(0)` was `0`, a plausible wrong value on a
public accessor, which is exactly the outcome the invariant ranks below a
crash. The module doc already called 1-based indexing a documented source of
silently-wrong results in Symmetrica; it was one here too.

The measurement is in [schubert.md](schubert.md): free, and the first attempt
at it was wrong by 43% in the scary direction because of code layout. The
same promotion was applied to `transpose`, `covers_right`, `covers_left`,
`QtPoly::mul_binomial`, `mul_diff`, `Frac::from_factors` and `mul_factors`,
none of which are on a hot path. `from_factors` was the one with real teeth:
`mul_binomial` covers its numerator branch, but a negative exponent inserts
`(0,0)` straight into the denominator, and nothing downstream catches a zero
denominator factor.

## Open

- **The release profile still wraps.** `overflow-checks = true` is item 1 of
  the policy and is untouched: `Cargo.toml` has no `[profile.release]`, so
  everything above documents failure paths in a profile whose arithmetic can
  still go quiet. That item, not this one, is the correctness gap.
- **`stanley`'s step cap is unmeasured.** Where in `S_n` a legitimate
  expansion first reaches `2²⁴` steps is unknown; so is whether the cap is
  above or below the point where the expansion stops fitting in memory.
- **`class_algebra_coefficient`'s `i128` wall is unmeasured.** `n!`, the
  common denominator and the character products all grow with `n` and the
  route holds nothing back.
- **Nothing here about rustdoc warnings.** They predate this work and are
  already owned by [../release-readiness.md](../release-readiness.md), which
  has the breakdown. Recorded only as a control: `cargo doc --no-deps` emitted
  183 before the audit and 183 after, so 46 new `# Panics` sections added
  none. ⚠️ That 183 is the default-feature build and does not match the 173
  release-readiness quotes for `--all-features`; neither number is wrong, they
  are different builds.
