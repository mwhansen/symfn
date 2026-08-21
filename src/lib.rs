//! A kernel for computing with symmetric functions: the six classical bases,
//! every transition between them, and the Hall–Littlewood, Macdonald, LLT and
//! Jack families above them, together with Schubert polynomials.
//!
//! Each basis is a distinct Rust type over a coefficient ring the caller
//! chooses, so a basis mix-up is a compile error and `ℚ[q,t]` is as ordinary a
//! coefficient ring as `ℤ`. Every value that leaves the library is exact, or
//! the call fails loudly.
//!
//! ## What it computes
//!
//! **The classical core.** The six bases — [`Schur`], [`Homogeneous`],
//! [`Elementary`], [`Monomial`], [`PowerSum`], [`Forgotten`] — as distinct
//! types unified by the [`SymFn`] trait, with [`convert()`] between every
//! ordered pair. Three products: ordinary — Littlewood–Richardson behind the
//! [`LrBackend`] trait, native in Rust with no external C library —
//! plethystic ([`plethysm()`]) and internal ([`kronecker`]). The full Hopf
//! structure: [`coproduct`], [`counit`], [`antipode`], and skewing by an
//! *arbitrary* symmetric function ([`SkewBy`]). Evaluation at a finite
//! alphabet, the principal specializations and their q-analogue ([`eval`]),
//! symmetric-group characters ([`character()`]) and Kostka numbers
//! ([`kostka()`]) — as single values and as whole tables.
//!
//! **The parametrized families.** Hall–Littlewood `Q'_λ(x;t)` and `P_λ(x;t)`
//! ([`hl`]) and the Kostka–Foulkes polynomials ([`kf`]); Macdonald `P_λ`,
//! `Q_λ` and `J_λ` ([`macdonald`]); the `(q,t)`-Kostka polynomials and the
//! modified Macdonald basis [`Ht`] ([`qtkostka`]); Jack `P/Q/J_λ(x;α)` and
//! the zonal specialization ([`jack`]); LLT polynomials in both the ribbon
//! and tuple models ([`llt`]); the delta-operator tower `∇`, `Δ'_f`, `Θ`
//! ([`deltaop`]); Schubert polynomials and their structure constants
//! ([`schubert`]).
//!
//! Everything is generic over [`Ring`], and the paths that divide ask only
//! for [`QAlgebra`], a ring containing `ℚ`, because every division in the
//! library is by `z_μ` — an *integer*. That is weaker than a field on
//! purpose: `ℚ[t]` and `ℚ[q,t]` are not fields, and they are exactly the
//! rings Hall–Littlewood and Macdonald need. Plethysm asks one thing more,
//! [`Plethystic`], since `p_n` acts on the coefficients too.
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
//! ## Conventions
//!
//! Every `K_{λμ}` variant, LLT `G`, and `H̃` in circulation differs from its
//! rivals by a normalization twist that yields plausible wrong answers, not
//! errors. Each family's module doc therefore opens with the conventions in
//! circulation, names the one that ships with references by equation number,
//! and pins it with a doctest whose value distinguishes it from its rivals —
//! "The conventions in circulation" in [`llt`] is the fullest example.
//!
//! ## Exactness
//!
//! **Every value that leaves this library is exact, or the call fails
//! loudly — in every build profile.** A computation has three legal outcomes
//! and no fourth: the exact answer; escalation to a wider ring and then the
//! exact answer; or a refusal the caller cannot mistake for an answer —
//! `None`, `Err`, a typed Python exception, or a panic the item documents. No
//! path returns a wrapped, truncated or rounded value.
//!
//! The default coefficient rings `i64`, `i128` and [`Rational`] are the fast
//! path: exact until a value leaves the width, and then a panic rather than a
//! wrap. `BigInt` and `BigRational` (the `bignum` feature) have no wall. In
//! between, [`Guarded`] and [`GuardedRat`] inside a [`guarded`] scope
//! *report* an overflow instead of panicking, so a caller can re-run the same
//! generic code over a bignum ring — which is what the Python boundary does
//! on every call. Where a fixed-width family has a wall a caller can reach,
//! its own module doc states the wall in reproducible terms, next to the
//! degrees and shapes the engine is measured to reach.
//! `docs/policies/failure.md` is the rulebook this compresses.
//!
//! ## Verification
//!
//! The test suite checks computed values against reference output from
//! independent software, committed to the tree: `tests/sage_oracle.rs` and
//! `tests/lrcalc_oracle.rs` compare against fixtures Sage and `lrcalc`
//! produced, with `scripts/gen_sage_oracle.sage` alongside so an auditor can
//! regenerate the fixtures rather than trust them; `tests/algebra_laws.rs`
//! checks the laws a fixed value cannot — conversions are ring homomorphisms,
//! `ω` is an involutive algebra map, `Δ` is an algebra map; and
//! `scripts/check_backend.py` runs this crate underneath Sage, so the inputs
//! are chosen by Sage rather than by these tests. Where no independent
//! software computes a family at all, the cross-check is internal: the
//! `(q,t)`-Kostka table is computed by three algorithms sharing nothing above
//! [`Partition`], and the three [`LrBackend`] implementations agree
//! exhaustively on every product with `|μ| + |ν| ≤ 7`.
//! `docs/policies/validation.md` is the standard each family is held to.
//!
//! ## Feature flags
//!
//! - **`bignum`** — `BigInt` / `BigRational` coefficient rings and the
//!   [`guarded`] escalation they complete. An entry point whose fixed-width
//!   path fails across most of its intended range — `ops::kronecker_coeff`
//!   is one — exists only under this feature.
//! - **`python`** — the PyO3 bindings the wheel is built from; enables
//!   `bignum`. The Python surface has its own reference at
//!   <https://symfn.readthedocs.io>.
//!
//! ## Examples
//!
//! Littlewood–Richardson over `ℤ`:
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
//!
//! Hall–Littlewood over `ℤ[q,t]` — the same [`Schur`] type, with [`QtPoly`]
//! coefficients:
//!
//! ```
//! use symfn::{hall_littlewood, Partition, QtPoly, SymFn};
//!
//! // Q'_{21}(x; t) = s_{21} + t·s_3
//! let qp = hall_littlewood::<i64>(&Partition::new([2, 1]));
//! assert_eq!(qp.coeff(&Partition::new([3])), QtPoly::t());
//! assert_eq!(qp.coeff(&Partition::new([2, 1])).coeff(0, 0), 1);
//! assert!(qp.coeff(&Partition::new([1, 1, 1])).is_empty());
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
// * **`pub(crate)`** — `memo`, `modular` and `candidates`, named from nowhere
//   outside `src/`. `clear_caches` is re-exported below because the
//   measurement discipline needs it (`CLAUDE.md`).
pub mod afrac;
#[doc(hidden)]
pub mod bh;
pub(crate) mod candidates;
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
pub mod interrupt;
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
    big_pi, big_pi_inverse, delta, delta_prime, delta_prime_e, nabla, nabla_e, nabla_power,
    schur_to_macdonald_ht, theta, Atom, Ratio,
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
pub use hl::{
    hall_littlewood, hall_littlewood_p, hall_littlewood_p_table, hall_littlewood_table,
    schur_to_hall_littlewood_p, schur_to_hall_littlewood_qp,
};
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
