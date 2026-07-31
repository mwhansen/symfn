#!/bin/sh
# The local gate: what must be green before a commit, in the order that fails
# fastest. There is no CI yet (docs/release-readiness.md, Phase 0); until
# there is, a green run of this script is what "the suite passes" means.
#
#     scripts/preflight.sh
#
# Needs nothing beyond the Rust toolchain: the oracle tests read committed
# fixtures, and nothing here touches Sage. The Sage-side harnesses
# (scripts/check_*.py) are separate, and run when the subsystem they oracle
# changes.

set -e

step() { printf '\n== preflight: %s\n' "$1"; }

step "cargo fmt --all --check"
cargo fmt --all --check

step "cargo test (default features)"
cargo test --quiet

step "cargo test --features bignum"
cargo test --quiet --features bignum

step "clean"
