# Coefficient arithmetic: the fixed-width rationals every family shares

`Rational` and `GuardedRat` (`src/coeff.rs`, `src/guard.rs`) are the field
under `s → p`, the Hall and internal products, plethysm, and — through the
`i128::div_exact` the same file defines — the reduction step of Jack's
`AFrac<i128>`. Everything that leaves ℤ passes through them, so a cost here is
paid by every family and a fix here is shared. This file records where that
cost actually was, what the shared fix is, and what it does not buy.

Split out of the family records because the premise that prompted it was
cross-family: profiles of the (q,t)-Kostka, Macdonald, plethysm and Jack work
had each reported 128-bit division (`u128_div_rem`, `__divti3`) at 15-38% of
wall time, and each family had routed around it locally — integer sweeps in
plethysm, the `i128` route for `H̃`, `Frac<i128>` for Macdonald, Barrett
reduction in `gj.rs`. The question was whether one shared fix was left.

---

## Where the 128-bit division still was

Sampled with `sample` on the `profiling` profile (Apple M4), on the tree
before the character recursion moved to β-masks
([transitions.md](transitions.md), "The character recursion on β-masks"):

| workload | ring | 128-bit division share |
|---|---|---|
| `s → p`, one term, `profile_convert loop s2p 24` | `Rational` | 24% (+5% `div_u128` self) |
| (q,t)-Kostka **branching check route**, `profile_qtk 9 loop` | `Frac<Rational>` | 49% (+45% `Rational` self time) |
| Jack `P → m`, `profile_jack 18` | `AFrac<i128>` | 13%, in `i128::div_exact` from `AFrac::reduce_at` |
| plethysm ladder, Macdonald P/J, ∇e₁₁, `H̃` table, LLT, HL, `p → s`, `hall` | various | 0-0.4% |

So the family figures were historical. The shipped (q,t)-Kostka route is BH
over `i128`, the shipped Macdonald ring is `Frac<i128>`, and plethysm's
`p → s` is the integer sweep with the ladder route on `den == 1` fast paths:
each local fix had removed the gcd from its hot loop entirely, and a faster
gcd cannot beat no gcd. What was left over ℚ on shipped paths was `s → p`
(and everything routed through it: `hall`, `internal`, the character bases,
the wheel's `s_to_p`/`h_to_p`/`e_to_p`) and Jack over `AFrac<i128>`.

## Why the remaining cost was call overhead, not arithmetic

The gcd was instrumented (temporary counters, since removed) across those
workloads: **99-100% of calls had both operands below 2³²**, mean 2-3 Euclid
steps. `z_μ` is divided one factor at a time by `Partition::div_by_z`, so
numerators stay reduced and the divisors are at most n. The whole cost was
that `i128 %` and `i128 /` are calls into `compiler_builtins` on aarch64 (and
x86-64) — about 9 ns each against one instruction for the 64-bit forms — and
every `Rational` operation made two to four of them. A standalone
microbenchmark on the observed operand shape: the `i128` Euclid at 21 ns per
call, the same narrowed to `u64` at 12 ns; on random 64-bit operands 177 ns
against 52 ns.

## The shared fix

Four changes, all in `Rational` and mirrored in `GuardedRat`, plus the two
integer helpers the other rings share:

1. **`gcd_u128` narrows to `u64`** when both magnitudes fit
   (`(a | b) >> 64 == 0`), and falls back to a binary gcd above. Plain `u64`
   Euclid measured the same in situ as a Euclid-step-then-binary hybrid, so
   the simple form is the one kept. `afrac.rs`'s own gcd and `convert.rs`'s
   integral-sweep gcd now go through it.
2. **`quo` narrows the exact quotients** (`num / g`, `den / g`) to `i64` the
   same way, and skips them when `g == 1`. This mattered as much as (1): after
   the gcd alone, `__divti3` from the quotients was still 32% of the
   (q,t)-Kostka check loop.
3. **`div_u128` without the second normalization.** `gcd(num, den) = 1` and
   `gcd(num/g, n/g) = 1` imply `gcd(num/g, den·(n/g)) = 1`, so the
   `Rational::new` it called re-ran a full gcd to learn nothing — half of all
   gcds in `s → p`.
4. **Henrici's addition and cross-cancelled multiplication** (Knuth, TAOCP
   4.5.1). `a/b + c/d` with `g = gcd(b, d)`: equal denominators need one gcd
   of the sum against `b`; `g == 1` means `(ad + cb)/(bd)` is already reduced,
   no gcd at all; otherwise `t = a(d/g) + c(b/g)`, `g' = gcd(t, g)`.
   Multiplication cancels `gcd(a, d)` and `gcd(c, b)` first. Beyond the saved
   gcds this shrinks the intermediates by a factor of `g` — with `z_μ`-derived
   denominators that share factors heavily, `GuardedRat` escalates to
   `BigRational` later than it did.

`i128::div_exact` (`div_exact_i128`) narrows to `i64` the same way, which is
what Jack's `AFrac<i128>` reduction runs on.

**`GuardedRat` had fallen behind `Rational`.** It had no `den == 1` fast
paths, its own gcd, and a `div_u128` that renormalized — so the wheel's
`plethysm`, `internal_product`, `s/h/e → p` and Jack over the guard rings had
never received the fast paths the record measured for `Rational` (the
1.07-1.09x on guarded plethysm below is that fast path finally arriving).
They cannot drift again: the arithmetic is written once, in `coeff.rs`
(`rat_normalize`, `rat_add`, `rat_mul`, `rat_div`), generic over an
`Overflow` policy that says only what an integer `+` or `*` does when it
leaves `i128` — `Panics` for `Rational` (native arithmetic under
`overflow-checks`), `Reporting` for `GuardedRat` (`checked_*`, note, continue
on `0`). Each ring keeps its own walls and injections, and `GuardedRat`
builds every result through `reduced`, which makes the two refusals `new`
makes (`den ≤ 0` after a reported overflow, a `MIN` part) without a gcd.
Measured after the extraction, both rings are 0.97-1.02x against the
hand-duplicated version on every case below — the policy monomorphizes to
what was there. The unification changes one `Rational` behavior, for the
better: cross-cancelling against a `MIN` numerator (which `from_int` admits)
used to panic in the gcd and now gives the exact answer, since the gcd is
against a positive denominator and fits; `new` and `neg` still refuse `MIN`.

## Measured

Min of 3 interleaved rounds, cold caches, Apple M4 **on battery** — both arms
in the same state, so the ratios stand and the absolute times read ~1.8x slow.
"Before" is the tree with the β-mask character recursion already in, so this
is the arithmetic change alone; the combined figures against the tree before
both are in the last column.

| case | before | after | | combined vs two commits back |
|---|---|---|---|---|
| `s → p`, one term: [24] / [8,8,8] / [10,5,3,2] / [6,5,4,3,2] / [7,6,5,4,3,2] | 0.68 / 0.95 / 0.28 / 0.25 / 1.61 ms | 0.32 / 0.55 / 0.17 / 0.17 / 1.12 ms | **2.12x / 1.71x / 1.65x / 1.47x / 1.44x** | 3.5x / 5.1x / 5.1x / 6.1x / 6.7x |
| the same over `GuardedRat` (the wheel's ring) | 0.72 / 0.97 / 0.29 / 0.26 / 1.63 ms | 0.34 / 0.53 / 0.17 / 0.17 / 1.10 ms | 2.12x / 1.86x / 1.72x / 1.53x / 1.48x | |
| `convert_s_to_p` (`bench_ops`, [6,5,4,3,2]) | 0.435 ms | 0.311 ms | **1.40x** | 6.5x |
| `internal` [5,3,2]·[4,4,2]; `kronecker` | 56 / 118 µs | 48 / 103 µs | 1.17x / 1.15x | |
| `hall(s, p)` and `p → s`, [8,8,8] | 0.154 / 0.153 s | 0.144 / 0.144 s | 1.07x / 1.06x | |
| plethysm s₆[s₆] / s₅[s₃₂] / s₄₂[s₃₁] over `Rational` | 0.147 / 0.093 / 0.0155 s | 0.147 / 0.092 / 0.0141 s | 1.00x / 1.01x / 1.10x | |
| the same over `GuardedRat` | 0.158 / 0.099 / 0.0152 s | 0.144 / 0.092 / 0.0138 s | 1.09x / 1.08x / 1.10x | |
| `jack_table::<i128>` n = 14 / 16 (from `div_exact` alone) | 0.099 / 0.399 s | 0.088 / 0.349 s | **1.12x / 1.15x** | |
| `jack_table::<Guarded>` n = 14 | 0.098 s | 0.087 s | 1.14x | |
| `jack_table::<Rational>` n = 13 (a check route) | 0.127 s | 0.076 s | 1.68x | |
| `qt_kostka_table_via_branching` n = 8 / 9 (check route, `Frac<Rational>`) | 0.140 / 0.583 s | 0.096 / 0.404 s | 1.46x / 1.44x | |
| `schur_in_j_table` n = 8 | 0.078 s | 0.051 s | 1.52x | |
| characters, Kostka, `s → m/h/e`, `m → s`, coproduct, `bench_ops` LR rows | — | — | 0.97-1.01x (untouched) | |

The `character_beta_sweep_n24/n28` rows of `bench_ops` read 0.93-0.95x
against the immediately preceding build and 1.00-1.03x against the one before
it; the sweep touches no rational arithmetic (it is `i128` through
`p_expand_shared`), and the preceding build had moved the same rows +6% for
equally unrelated reasons. That is code layout around a hot loop that neither
change touched, and it nets to nothing.

Suites: `cargo test`, `cargo test --features bignum` and the doctests all
green; the arithmetic is exact either way and every output is byte-identical.

## What it does not buy, and where the next levers are

* **The shipped headline paths were already off ℚ.** Macdonald, `H̃`, ∇,
  LLT, HL and plethysm's `p → s` gain nothing here because their local fixes
  removed the gcd rather than speeding it. The premise's "15-38% across
  families" was true when each record measured it and false by the time this
  was asked; the family records stand as history, and this file is the
  present.
* **The (q,t)-Kostka branching route** (`Frac<Rational>`, a cross-check, not
  the shipped route) is now 92% `Rational::add_assign`/`mul` self time inside
  `Frac`'s numerator polynomials, because `Frac::div_u128` divides every
  coefficient and the polynomials carry genuinely fractional coefficients
  from then on. `AFrac` solved exactly this with an integer `scale`; a `Frac`
  scale would let that route run over `Frac<i128>`. Low value while it stays
  a check route.
* **`s → p` after both changes** (`sample`, the staircase of 27, harness
  timer excluded): the mask recursion itself 31%, the memo's local inserts and
  the per-call merge into the shared table 21%, malloc/free 16% (the output
  `BTreeMap` and the `partitions_cached` walk), `div_u128` 8% with no 128-bit
  division routine left in the top thirty frames. The next lever is the 21%:
  one `try_character` call per (λ, μ) merges its new entries into the shared
  table p(n) times per row, and a row-batched recursion (one local memo for
  all μ of a λ, one merge) would take most of it.

## Tried and dropped

* **A Euclid-step-then-binary hybrid at 64 bits** (one `%` to equalize sizes,
  then Stein's algorithm): identical to plain Euclid in situ, on every case
  above; the operands are 2-3 steps from done either way, and the branch-free
  loop buys nothing at that length. Kept only for the wide fallback, where it
  avoids the 128-bit `%` altogether.
* **Looking for a shared fix in `Frac`.** `Frac<i128>` has no integer scale
  and no gcd on its hot path (its denominators are factored binomials), so the
  Macdonald and `H̃` profiles have no rational arithmetic to speed up.

## Pinned by

`shortcut_arithmetic_matches_reduce_after_the_fact` (every branch of the
addition, multiplication and `div_u128` shortcuts against "form the cross
product and reduce it", on a grid of both signs, zeros, equal, coprime and
factor-sharing denominators, and 40-bit operands),
`narrowed_gcd_and_quotient_agree_with_the_wide_forms` (the `u64`/binary gcd
against Euclid on operands straddling 64 bits and on products of large primes;
`quo` and `div_exact_i128` at the `i64` edges),
`guarded_shortcut_arithmetic_matches_rational` (the same grid, `GuardedRat`
against `Rational` bit for bit, reporting nothing) and
`shortcut_paths_refuse_what_normalization_refuses` (a product landing on
`MIN` and a coprime-denominator overflow both leave the guarded scope as
`None`, not as a value).
