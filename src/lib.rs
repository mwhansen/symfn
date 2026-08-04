//! # symfn — a modern kernel for symmetric functions
//!
//! A clean-room Rust library for computing with symmetric functions, designed
//! as a fast, testable, Sage-interoperable successor *in spirit* to Symmetrica
//! (it shares no code with the old C library).
//!
//! ## Design principles
//!
//! - **Real types, not a tagged union.** Each basis is its own type
//!   ([`Schur`], [`PowerSum`], [`Monomial`]) unified by the [`SymFn`] trait, so
//!   basis confusion is a compile error rather than a runtime bug. Contrast
//!   Symmetrica's single untyped `OP` object.
//! - **Coefficient ring is a parameter.** Everything is generic over [`Ring`];
//!    the paths that divide ask only for [`QAlgebra`] (a ring containing ℚ),
//!    because every division in the library is by z_μ — an *integer*. That is
//!    weaker than a field on purpose: `ℚ[t]` and `ℚ[q,t]` are not fields, and
//!    they are exactly the rings Macdonald/Hall–Littlewood need. Plethysm asks
//!    for one thing more, [`Plethystic`], since `p_n` acts on the coefficients
//!    too. The scaffold uses `i64` and [`Rational`]; the `bignum` feature swaps
//!    in `BigInt` / `BigRational`.
//! - **Swappable backends behind traits.** Littlewood–Richardson lives behind
//!   [`LrBackend`] and is computed **natively in Rust** ([`NaiveLr`]) — no
//!   external C library. The trait lets a future optimized backend (memoized /
//!   DP) drop in without touching callers, cross-checked against this one.
//! - **Correct by construction, tested against an oracle.** Sage computes all
//!    of this correctly (if slowly); those values are the test oracle. Unit
//!    tests here pin known expansions; the `tests/` integration suite checks
//!    algebraic laws, and property tests vs. Sage slot in once dependencies are
//!    available.
//!
//! ## The overflow contract
//!
//! **Every value that leaves this library is exact, or the call fails loudly —
//! in every build profile.** No path returns a wrapped, truncated, or rounded
//! result. A computation has three legal outcomes and no fourth:
//!
//! 1. the exact answer;
//! 2. escalation to a wider ring, then the exact answer;
//! 3. a loud refusal — `None`, `Err`, a typed Python exception, or a documented
//!    panic.
//!
//! What that means in practice, for a caller choosing a coefficient type:
//!
//! - **`i64` / `i128` / [`Rational`]** are exact until a value leaves the
//!   width, and then they panic. The release profile carries
//!   `overflow-checks = true` — measured at 0–6% across every harness the
//!   crate has — so this holds in the profile you ship, not only in debug.
//!   Structure constants injected through [`Ring::from_u128`] check the same
//!   way, naming the constant and the ring.
//! - **`BigInt` / `BigRational`** (the `bignum` feature) have no wall.
//! - **[`Guarded`] / [`GuardedRat`]**, inside a [`guarded`] scope, *report*
//!   instead of panicking: `None` means "an intermediate left the width", and
//!   the caller re-runs the same generic code over a bignum ring. That two-pass
//!   escalation is what the Python boundary and `ops::kronecker_coeff` do,
//!   and it costs 0–1% on the fast path.
//!
//! Where a fixed-width family has a wall a caller can reach, its own docs state
//! that wall in reproducible terms. The rules behind all of this, and which
//! mechanism each situation demands, are in `docs/policies/failure.md`; what
//! executing them cost and turned up is in
//! `docs/record/failure-and-overflow.md`.
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

pub mod afrac;
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
pub mod macop;
pub mod measure;
pub mod memo;
pub mod modular;
pub mod ops;
pub mod partition;
pub mod permutation;
pub mod plethysm;
#[cfg(feature = "python")]
pub mod python;
pub mod qt;
pub mod qtkostka;
pub mod rect;
pub mod schubert;
pub mod skew_lr;
pub mod strip_lr;
pub mod sym;
pub mod three_row;
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
