//! A kernel for computing with symmetric functions: the six classical bases,
//! every transition between them, and the Hall–Littlewood, Macdonald, LLT and
//! Jack families above them, together with Schubert polynomials.
//!
//! Each basis is a distinct Rust type over a coefficient ring the caller
//! chooses, so a basis mix-up is a compile error and `ℚ[q,t]` is as ordinary a
//! coefficient ring as `ℤ`. Three products — ordinary, plethystic and internal
//! (Kronecker) — the full Hopf structure, symmetric-group characters, and
//! evaluation at a finite alphabet sit on top. It is a clean-room successor
//! *in spirit* to Symmetrica, sharing no code with that library, and the
//! default build has no dependencies.
//!
//! ## What this opens
//!
//! - **A single structure constant for a product that cannot be
//!   materialized.** `schubert::schubert_coeff` answers `c^w_{uv}` by Bruhat
//!   pruning for pairs whose product no machine holds, and
//!   `ops::kronecker_coeff` (under `bignum`) answers one `g^ν_{λμ}` where the
//!   whole internal product does not fit. The query itself is the capability,
//!   not a faster route to the full product.
//! - **Reduced (stable) Kronecker coefficients as an outer product** in the
//!   Orellana–Zabrocki bases ([`character_basis`]). The calculation never
//!   leaves the power-sum basis, so no Littlewood–Richardson coefficient
//!   enters it at all.
//! - **Matchings–Jack and b-conjecture coefficients** ([`gj`]), which no other
//!   package computes, and the two [`llt`] outputs with no Sage entry point at
//!   any speed: the `∇e_n` by-path Schur-positive refinement, and parabolic
//!   affine Kazhdan–Lusztig columns.
//!
//! `docs/research-gaps.md` is the measured survey of the incumbents behind
//! those claims, including what does exist elsewhere; `docs/record/` carries
//! the degrees and shapes each engine reaches.
//!
//! ## Exactness
//!
//! **Every value that leaves this library is exact, or the call fails
//! loudly.** A computation has three legal outcomes and no fourth: the exact
//! answer; escalation to a wider ring and then the exact answer; or a refusal
//! the caller cannot mistake for an answer — `None`, `Err`, a typed Python
//! exception, or a panic the item documents. No path returns a wrapped,
//! truncated or rounded value.
//!
//! Fixed width is the fast path, not the promise. The default coefficient
//! rings are `i64`, `i128` and [`Rational`], and a coefficient outgrowing one
//! of them is an ordinary event rather than a bug: entry points that promise
//! exactness compute over [`Guarded`] / [`GuardedRat`], which report leaving
//! the fixed width instead of wrapping, and [`guarded`] turns that report into
//! a re-run of the same generic code over `BigInt` / `BigRational` — the
//! `bignum` feature, which `python` always enables. A path that cannot
//! escalate returns `Option` or `Result` ([`try_character`]) or documents its
//! wall under `# Panics`, and a function whose fixed-width path fails across
//! most of its intended range exists only under `bignum` rather than existing
//! and refusing.
//!
//! - **`i64` / `i128` / [`Rational`]** are exact until a value leaves the
//!   width, and then they panic. The release profile carries
//!   `overflow-checks = true`, so this holds in the profile you ship, not only
//!   in debug. Structure constants injected through [`Ring::from_u128`] check
//!   the same way, naming the constant and the ring.
//! - **`BigInt` / `BigRational`** (the `bignum` feature) have no wall.
//! - **[`Guarded`] / [`GuardedRat`]**, inside a [`guarded`] scope, *report*
//!   instead of panicking: `None` means "an intermediate left the width", and
//!   the caller re-runs the same generic code over a bignum ring. That two-pass
//!   escalation is what the Python boundary and `ops::kronecker_coeff` do.
//!
//! Where a fixed-width family has a wall a caller can reach, its own docs
//! state that wall in reproducible terms. `docs/policies/failure.md` is the
//! rulebook this compresses — which mechanism each situation demands — and
//! `docs/record/failure-and-overflow.md` is what executing it cost and turned
//! up.
//!
//! ## The model
//!
//! - **A type per basis, not a tagged union.** The six classical bases are
//!   distinct types — [`Schur`], [`Homogeneous`], [`Elementary`],
//!   [`Monomial`], [`PowerSum`], [`Forgotten`] — unified by the [`SymFn`]
//!   trait, with [`convert()`] between every ordered pair through the Schur
//!   hub. The modified Macdonald and Orellana–Zabrocki bases, [`Ht`] and
//!   [`St`], are types on the same footing. Basis confusion is a compile
//!   error rather than a wrong answer; contrast Symmetrica's single untyped
//!   `OP` object.
//! - **The coefficient ring is a parameter.** Everything is generic over
//!   [`Ring`]; the paths that divide ask only for [`QAlgebra`], a ring
//!   containing `ℚ`, because every division in the library is by `z_μ` — an
//!   *integer*. That is weaker than a field on purpose: `ℚ[t]` and `ℚ[q,t]`
//!   are not fields, and they are exactly the rings Hall–Littlewood and
//!   Macdonald need. Plethysm asks one thing more, [`Plethystic`], since
//!   `p_n` acts on the coefficients too.
//! - **Littlewood–Richardson behind a trait, native in Rust.** [`LrBackend`]
//!   has three implementations and no external C library: [`NaiveLr`], the
//!   in-house oracle; [`StripLr`], a row-strip DP over horizontal strips; and
//!   [`SkewLr`], which expands a whole skew shape in one traversal of a
//!   merged layer and carries the general case. [`AutoLr`] is where
//!   dispatch lives, taking the closed-form and counting routes for
//!   rectangles and few-row factors first. The three backends agree
//!   exhaustively on every product with `|μ| + |ν| ≤ 7`, so the choice is
//!   unobservable except in timing.
//! - **Oracles are committed, not assumed.** `tests/sage_oracle.rs` and
//!   `tests/lrcalc_oracle.rs` check against fixtures Sage and `lrcalc`
//!   produced, with `scripts/gen_sage_oracle.sage` in the tree so an auditor
//!   can regenerate rather than trust; `tests/algebra_laws.rs` checks the laws
//!   a value pin cannot — conversions are ring homomorphisms, `ω` is an
//!   involutive algebra map, `Δ` is an algebra map; and
//!   `scripts/check_backend.py` runs this crate in place of Symmetrica
//!   underneath Sage, so the inputs are chosen by Sage rather than by these
//!   tests.
//!
//! Each module's docs carry its conventions, its references by equation
//! number, and its traps: the circulating conventions in [`llt`] are the
//! fullest example, and naming which normalization ships is what keeps a
//! wrong-by-a-twist answer from passing for a right one.
//!
//! ## Example
//!
//! ```
//! use symfn::{Schur, SymFn, Partition};
//!
//! // s_2 · s_1 = s_3 + s_{21}
//! let s2: Schur<i64> = Schur::monomial(Partition::new([2]), 1);
//! let s1: Schur<i64> = Schur::monomial(Partition::new([1]), 1);
//! let prod = s2.mul(&s1);
//! assert_eq!(prod.coeff(&Partition::new([3])), 1);
//! assert_eq!(prod.coeff(&Partition::new([2, 1])), 1);
//! ```

// Affordable only because the module sort below shrank the surface it applies
// to: the lint skips `#[doc(hidden)]` items, so the nine hidden modules and the
// two private ones account for 17 of the 53 undocumented items this found, and
// the 36 that remained were accessors and trait-method signatures rather than
// mathematics (`docs/release-readiness.md`, Phase 2).
#![deny(missing_docs)]

// The public module list is a decided list, not an accumulated one
// (`docs/release-readiness.md`, Phase 2). Three tiers, and the test that sorts
// them is whether a caller who only wants symmetric functions would ever name
// the module:
//
// * **API** — the families, the types, the coefficient rings. Documented,
//   semver-stable, and what the policy in `README.md` promises.
// * **`#[doc(hidden)]`** — reachable and compiled, absent from the reference,
//   promised nothing. Two kinds live here: the cross-check engines that exist
//   to disagree with a primary route (`bh`, `gjmod`, `macop`), and the
//   strategy modules whose *results* are API but whose paths are not (`rect`,
//   `two_row`, `three_row`, `strip_lr` — their re-exports below stay
//   documented). `measure` is a heap-accounting harness and `python` is a
//   PyO3 bridge; neither is symmetric functions.
// * **`pub(crate)`** — `memo` and `modular`, named from nowhere outside
//   `src/`. `clear_caches` is re-exported below because the measurement
//   discipline needs it (`CLAUDE.md`).
pub mod afrac;
#[doc(hidden)]
pub mod bh;
pub mod character;
pub mod character_basis;
pub mod charge;
pub mod coeff;
pub mod convert;
pub mod deltaop;
pub mod dyck;
pub mod eval;
mod fasthash;
pub mod frac;
pub mod gj;
#[doc(hidden)]
pub mod gjmod;
pub mod guard;
pub mod hl;
pub mod hopf;
pub mod jack;
pub mod kf;
pub mod kostka;
pub mod llt;
pub mod lr;
pub mod macdonald;
#[doc(hidden)]
pub mod macop;
#[doc(hidden)]
pub mod measure;
pub(crate) mod memo;
pub(crate) mod modular;
pub mod ops;
pub mod partition;
pub mod permutation;
pub mod plethysm;
#[cfg(feature = "python")]
#[doc(hidden)]
pub mod python;
pub mod qt;
pub mod qtkostka;
#[doc(hidden)]
pub mod rect;
pub mod schubert;
pub mod skew_lr;
#[doc(hidden)]
pub mod strip_lr;
pub mod sym;
#[doc(hidden)]
pub mod three_row;
#[doc(hidden)]
pub mod two_row;

pub use afrac::AFrac;
pub use character::{character, character_in, try_character};
pub use character_basis::{
    ht_product_terms, reduced_kronecker, reduced_kronecker_product, reduced_kronecker_via_ht,
};
pub use charge::{charge, kostka_foulkes_by_charge};
pub use coeff::{Field, Plethystic, QAlgebra, Rational, Ring};
pub use convert::{convert, FromSchur, ToSchur};
pub use deltaop::{
    big_pi, big_pi_inverse, delta, delta_prime, delta_prime_e, nabla, nabla_e, nabla_power, theta,
    Atom, Ratio,
};
pub use dyck::{ladder, ladder_at_content, side, side_at_content, Side};
pub use eval::{dimension, principal_specialization, principal_specialization_q};
pub use frac::Frac;
pub use gj::{
    class_algebra_coefficient, double_coset_coefficient, double_coset_table, gj_connection_tables,
    matchings_jack_coverage, BPoly, Coverage, GjTables,
};
#[doc(hidden)]
pub use gjmod::{engines_agree, gj_connection_tables_modular};
pub use guard::{guarded, Guarded, GuardedRat};
pub use hl::{hall_littlewood, hall_littlewood_p, hall_littlewood_p_table, hall_littlewood_table};
pub use hopf::{antipode, coproduct, counit, skew_schur, SkewBy, SymTensor};
pub use jack::{
    hook_lower, hook_upper, jack_j, jack_j_powersum, jack_j_table, jack_j_tableaux, jack_norm_j,
    jack_norm_p, jack_p, jack_p_branching, jack_p_lb, jack_powersum_table, jack_q, jack_scalar,
    jack_structure_constant, jack_table, omega_alpha, powersum_scalar, stanley_table, zonal_j,
    zonal_p,
};
pub use kf::{kostka_foulkes, kostka_foulkes_column, kostka_foulkes_table};
pub use kostka::{kostka, semistandard_tableaux};
pub use llt::{
    chromatic_from_llt, htilde_by_llt, llt_e_expansion, llt_fundamental, llt_g, llt_g_lt,
    llt_graph, llt_gtilde, llt_gtilde_table, llt_h, llt_h_table, llt_h_tilde, llt_kl_column,
    llt_max_inv, llt_min_inv, llt_schur, nabla_e_by_path, DecoratedGraph, SkewTuple,
};
pub use lr::{LrBackend, NaiveLr};
pub use macdonald::{macdonald_j, macdonald_p, macdonald_q};
#[doc(hidden)]
pub use macop::{eigenvector, eigenvectors, operator_matrix};
pub use memo::clear_caches;
pub use ops::{hall, internal, kronecker, omega};
pub use partition::{partitions_of, Partition, PartitionError};
pub use plethysm::plethysm;
pub use qt::QtPoly;
pub use qtkostka::{
    macdonald_ht, modified_qt_kostka, qt_kostka, qt_kostka_column, qt_kostka_table,
    qt_kostka_table_via_bh, qt_kostka_table_via_branching, qt_kostka_table_via_operator,
    schur_in_j_table,
};
pub use rect::{okada_coeff, okada_product};
pub use skew_lr::{expand_skew, expand_skew_shared, SkewLr};
pub use strip_lr::{AutoLr, StripLr};
pub use sym::{
    Elementary, Forgotten, Homogeneous, Ht, Monomial, PowerSum, Schur, St, SymAlgebra, SymFn,
};
pub use three_row::three_row_product;
pub use two_row::{two_row_coeff, two_row_product};
