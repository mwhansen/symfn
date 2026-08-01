//! The canary for `[profile.release] overflow-checks = true`.
//!
//! ```text
//!   cargo test --release --test overflow_checks
//! ```
//!
//! The escalation ladder handles the overflow that was foreseen; the profile
//! flag is what converts the overflow that was *not* from silent to loud
//! (`docs/policies/failure.md`, R3). That makes the flag part of the
//! correctness surface rather than a build detail — and a build detail is
//! exactly what a stray edit to `Cargo.toml` would treat it as. These tests
//! fail the day the line goes, which is the only way that edit becomes visible.
//!
//! They pass in debug too, where the flag is on by default: the assertion is
//! "arithmetic in this crate is checked", not "this profile differs from that
//! one". `--release` is the run that carries information.
//!
//! Everything goes through [`Ring`], not through bare operators, because the
//! `Ring` impls are the seam the whole library's generic code multiplies at
//! (`impl_ring_for_int!` in `src/coeff.rs` uses a plain `*`), and `black_box`
//! keeps the optimiser from const-folding the overflow into a compile error.

use std::hint::black_box;

use symfn::{Rational, Ring};

#[test]
#[should_panic(expected = "overflow")]
fn i64_ring_mul_is_checked() {
    let big = black_box(i64::MAX / 2 + 1);
    let _ = black_box(Ring::mul(&big, &black_box(3i64)));
}

#[test]
#[should_panic(expected = "overflow")]
fn i128_ring_mul_is_checked() {
    let big = black_box(i128::MAX / 2 + 1);
    let _ = black_box(Ring::mul(&big, &black_box(3i128)));
}

#[test]
#[should_panic(expected = "overflow")]
fn i128_ring_add_is_checked() {
    let mut big = black_box(i128::MAX);
    Ring::add_assign(&mut big, &black_box(1i128));
    black_box(big);
}

/// `Rational`'s integer fast path — the one the (q,t) families and every
/// `s → p` run through, and the site the Kronecker wall shows up at.
#[test]
#[should_panic(expected = "overflow")]
fn rational_ring_mul_is_checked() {
    let big = black_box(Rational::new(i128::MAX / 2 + 1, 1));
    let _ = black_box(Ring::mul(&big, &black_box(Rational::new(3, 1))));
}
