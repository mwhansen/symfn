#!/bin/sh
# The local gate: what must be green before a commit, in the order that fails
# fastest. There is no CI yet (docs/release-readiness.md, Phase 0); until
# there is, a green run of this script is what "the suite passes" means.
#
#     scripts/preflight.sh
#
# Needs nothing beyond the Rust toolchain and python3: the oracle tests read
# committed fixtures, and nothing here touches Sage. The Sage-side harnesses
# (scripts/check_*.py) are separate, and run when the subsystem they oracle
# changes — except check_panics_documented.py, which reads sources only and is
# run here because the invariant it pins decays on any commit adding a `pub fn`.

set -e

step() { printf '\n== preflight: %s\n' "$1"; }

step "cargo fmt --all --check"
cargo fmt --all --check

step "public panics documented"
python3 "$(dirname "$0")/check_panics_documented.py"

step "cargo test (default features)"
cargo test --quiet

step "cargo test --features bignum"
cargo test --quiet --features bignum

step "clean"
