# The public API, and what a version number promises

[layout.md](layout.md) maps the tree; this file maps the *interface*, which
is smaller. The public module list is decided rather than accumulated, and
the test that sorts it is whether a caller who only wants symmetric functions
would ever name the module. Three tiers, checkable with `cargo doc --no-deps`:

- **API.** Every family and every conversion, plus the types and coefficient
  rings they are written over: `sym`, `convert`, `ops`, `hopf`, `plethysm`,
  `eval`, `kostka`, `character`, `character_basis`, `charge`, `kf`, `hl`,
  `jack`, `macdonald`, `qtkostka`, `llt`, `schubert`, `lr`, `skew_lr`,
  `deltaop`, `dyck`, `gj`, `partition`, `permutation`, `coeff`, `guard`, `qt`,
  `frac`, `afrac`. Every item in them carries a doc comment —
  `#![deny(missing_docs)]` is what keeps that true rather than a habit.
- **Reachable, promised nothing** — all `#[doc(hidden)]`. `bh`, `gjmod` and
  `macop` are cross-check engines that exist to disagree with a primary route;
  `rect`, `two_row`, `three_row` and `strip_lr` are product strategies whose
  *results* are API even though their paths are not (`okada_coeff`,
  `okada_product`, `two_row_coeff`, `two_row_product`, `three_row_product`,
  `AutoLr` and `StripLr` are re-exported at the crate root and documented
  there); `measure` is heap accounting and `python` is the PyO3 bridge, and
  neither is symmetric functions. Naming one of these compiles. It is not a
  promise, and it is where a break lands without a version bump.
- **Private.** `memo`, `modular`, `fasthash`, `candidates`. `clear_caches` is
  re-exported at the crate root, because timing a run means clearing the
  caches between them, and `cache_stats`, `cache_budget` and
  `set_cache_budget` beside it, because a session that can clear the caches
  should be able to see what they hold and bound it.

**What consumers build on is the coefficient-ring layer**: `Ring`, and the
`QAlgebra` and `Plethystic` refinements above it. Generic code bounded on those
three is what survives a basis or backend being rewritten underneath it — the
bound is deliberately weaker than `Field` so that ℚ[t] and ℚ[q,t] qualify,
which is what Macdonald and Hall–Littlewood need. `LrBackend` is the same shape
one level down: three native backends implement it, and swapping one for
another changed no caller.

## Operators on the element types

The basis types, `Schubert` and `SymTensor` implement the `core::ops` traits
over their named methods: `+`, `-`, unary `-`, `*` by a coefficient on the
right, the assigning forms, and `*` between elements where the type has a
product. Each binary operator is implemented for both operands owned, both
borrowed, and each mixed pair. That is decided before the first tag rather
than left additive because a by-value impl added later changes how existing
calls resolve: with `core::ops::Mul` in scope, `a.mul(&b)` resolves to the
operator's by-value method and moves `a` (rustc E0382, checked 2026-09-03).
Shipping the by-value impls in the first release means no later release
changes that resolution. The named methods stay, because generic code bounded
on `SymFn` has no operator bounds to use. A coefficient on the left, `c * f`,
is not implemented: the impl would have to be on the coefficient type, and
the orphan rules forbid that for a type parameter. `Sum` and `Product` over
iterators are not implemented; adding them later is additive.

## Pre-release

What a version number promises is decided at the first release
([release-readiness.md](release-readiness.md)); the crate is at 1.0.0-rc and
until then nothing about the number line is promised. The hidden tier
sits outside whatever is decided — it can move in a patch, which is why it
is a separate tier rather than a naming convention.

Two things break API-tier callers that do not look like breaks, whatever the
policy ends up being, so they are worth naming now:

- **A method added to `Ring`, `SymFn`, `LrBackend` or `SkewBy`** breaks any
  code implementing the trait outside this crate, while breaking no caller.
  The coefficient-ring traits are the seams this library changes behavior
  through ([policies/failure.md](policies/failure.md)), so they are
  the ones most likely to gain a method.
- **The Python surface is frozen harder than the crate**, not in step with it.
  It is the contract nearly every consumer reaches this library through, and
  [policies/python.md](policies/python.md) is its rulebook; a
  crate-internal change is free, and the same change at that boundary is not.
