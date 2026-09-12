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
  `frac`, `afrac`, and `interrupt`, which an embedder names to stop a running
  call. Every item in them carries a doc comment —
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

## The crate root

The root re-exports what a consumer names and nothing else. The test is the
module test applied to items, and it sorts them into three groups:

- **At the root.** The types and traits every call is written with —
  `Partition`, the basis types, `SymFn`, `Ring` and its refinements, the
  coefficient rings, `LrBackend` and the backends that implement it;
  each family's entry points, in every basis and normalization it is
  computed in; the whole-degree tables and columns; the conversions between
  bases; the Macdonald and Delta operators and the pairings; the `try_`
  twins beside their primaries; and the cache controls. A type that appears
  in the signature of a root function is at the root with it — `CacheStat`,
  `GjTables`, `BPoly`, `Atom`, `Ratio`, `DecoratedGraph`, `SkewTuple`.
- **Behind the module path only.** Still public and documented there, and
  promotable later without a break. A second route to a value the root
  already computes: `jack::jack_p_lb`, `jack::jack_p_branching`,
  `jack::jack_j_tableaux`, `qtkostka::qt_kostka_table_via_bh`,
  `qtkostka::qt_kostka_table_via_branching`,
  `qtkostka::qt_kostka_table_via_operator`,
  `character_basis::reduced_kronecker_via_ht`,
  `charge::kostka_foulkes_by_charge`, `llt::htilde_by_llt` and
  `skew_lr::expand_skew_shared`. Element arithmetic over plain maps, which
  exists for the Python bridge: `jack::jack_element_add`,
  `jack::jack_element_scale`, `macdonald::macdonald_element_add`,
  `macdonald::macdonald_element_scale`, `deltaop::htilde_element_add` and
  `deltaop::htilde_element_scale`. Building blocks and enumeration
  primitives: `jack::hook_lower`, `jack::hook_upper`,
  `character_basis::ht_product_terms`, `kostka::semistandard_tableaux`,
  `charge::charge`, `llt::llt_min_inv` and `llt::llt_max_inv`. A
  classification against the literature rather than a computation:
  `gj::matchings_jack_coverage` and its `Coverage`.
- **Not re-exported.** The hidden engines `gjmod` and `macop` had hidden root
  re-exports; those are gone, because a hidden re-export promised nothing
  the module path does not.

The re-exports from the hidden strategy modules — `okada_coeff`,
`okada_product`, `two_row_coeff`, `two_row_product`, `three_row_product`,
`AutoLr`, `StripLr` — are unchanged. For those, leaving the root means
leaving the reference, which is a tier question and was decided above.

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

## Operators on the coefficient types

The coefficient types in the API tier — `Rational`, `Guarded`, `GuardedRat`,
`QtPoly`, `Frac`, `AFrac` and `Ratio` — implement `+`, `-`, unary `-`, `*`,
`+=`, `-=` and `*=` over their `Ring` methods, each binary operator for both
operands owned, both borrowed, and each mixed pair. `Rational` also
implements `/` and `/=` over `Field::div`. One macro, `impl_ring_ops` in
`src/coeff.rs`, writes every impl. Decided 2026-09-11, for 0.9.0.

The reason is the one the element types' section gives, and it is stronger
here, because `Ring` and `Field` use the method names `core::ops` uses:
`mul`, `neg`, `add_assign`, `sub_assign`, `div`. Checked 2026-09-11 on rustc
1.98.0 with a stand-in type carrying both a `Ring`-shaped trait and the
operator impls, in a module that imports the operator trait:

- `a.mul(&b)` and `a.neg()` on an owned `a` resolve to the operator's
  by-value method and move `a`; a later use of `a` is E0382.
- `r.mul(&b)` and `r.neg()` on a borrowed `r`, and `a.add_assign(&b)` on
  either, find two methods at the same step and are E0034.

No impl shape avoids this. Without the by-value impls, the owned case finds
`Ring::mul` and the by-reference operator at the same autoref step and is
E0034 as well. So whichever release first ships the impls changes what such
a call means, for any caller that imports an operator trait, and 0.9.0 is the
only release with no callers to break. A module that imports both traits
calls the named methods through the trait, as `Ring::mul(&a, &b)`; the
module doc of `src/coeff.rs` shows that, and pins the E0034 case with a
`compile_fail` doctest. Code inside this crate is unaffected: no module
imports an operator trait by name (the impls use `core::ops::` paths), and
generic code over `C: Ring` cannot see impls on concrete types, so no
existing call changes what it resolves to.

The named methods stay, because generic code bounded on `Ring` has no
operator bounds to use. Adding the operator traits as supertraits of `Ring`
would break every ring implemented outside this crate (the trait-method rule
below), so it is not done.

`/` goes with `Field`: a type implements `Div` exactly when it implements
`Field`, and both arrive in the same change, so there is never a `.div` call
already written for the new impl to re-resolve. `Frac`, `AFrac` and `Ratio`
are not fields in this sense and have no `/`.

Not covered: `i64` and `i128` are primitives; `BigInt` and `BigRational`
carry their own crate's operators; `bh::Rat` is in the hidden tier, which
may change in a patch. A right operand of another type — `QtPoly<C> * C`,
`Rational + i64` — is not implemented. Adding one later adds an impl of a
trait whose by-value method a call already finds, rather than a trait it did
not find before. `Sum` and `Product` are not implemented; adding them later is
additive. The operators call the `Ring` methods and fail exactly as those do
([policies/failure.md](policies/failure.md)): `Guarded` reports overflow,
`Rational` panics on it.

## The number line

The first release is **0.9.0**, not 1.0.0. Decided 2026-09-04 with the
2026-09-03 review's finding in view: the tree has no external caller yet,
the root re-exports 161 names, and the API tier holds about 400 `pub fn`s.
A 1.0 tag would freeze every one of those under semver on the strength of
six weeks of in-house use, and the cost of a wrong promise is paid forever.
A 0.x first release gathers callers first; 1.0 is cut once the surface has
held still for a few months under them.

What the number promises, from 0.9.0 on:

- **A patch** (`0.9.x`) adds or fixes and breaks nothing: no signature,
  name, or convention in the API tier changes, and no re-export leaves the
  root. Cargo resolves `0.9` as compatible with every `0.9.x`, so this is
  what a consumer's version requirement relies on.
- **A minor** (`0.10.0`) may break the API tier, and the `CHANGELOG.md`
  entry names every break. Under the Cargo convention a minor bump before
  1.0 is the breaking bump, and this tree uses it that way rather than
  promising a stricter rule no consumer holds yet.
- **The hidden tier** can move in a patch, which is why it is a separate
  tier rather than a naming convention.
- **The Python contract layer** changes only in a minor, never a patch —
  the same rule as the crate's API tier, and the one that matters most,
  because the Sage adapter pins a version range and a break there is a Sage
  bug ([policies/python.md](policies/python.md)). The convenience layer is
  held to the same line, since a consumer cannot tell the two apart.

The version lives in `Cargo.toml` alone; `pyproject.toml` and
`symfn.__version__` read it from there, and the release workflow refuses a
tag that disagrees with it. The tree carries a `-rc.N` suffix until the tag is
cut, so a wheel built from the working tree cannot be mistaken for the
release.

Two things break API-tier callers that do not look like breaks, so they
count as breaks under the rule above:

- **A method added to `Ring`, `SymFn`, `LrBackend` or `SkewBy`** breaks any
  code implementing the trait outside this crate, while breaking no caller.
  The coefficient-ring traits are the seams this library changes behavior
  through ([policies/failure.md](policies/failure.md)), so they are
  the ones most likely to gain a method.
- **The Python surface is frozen harder than the crate**, not in step with it.
  It is the contract nearly every consumer reaches this library through, and
  [policies/python.md](policies/python.md) is its rulebook; a
  crate-internal change is free, and the same change at that boundary is not.
