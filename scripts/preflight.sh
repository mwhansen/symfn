#!/bin/sh
# The local gate: what must be green before a commit, in the order that fails
# fastest. CI runs the same suites and more on every push
# (.github/workflows/ci.yml); a green run of this script is what "the suite
# passes" means locally, and it is what the pre-commit hook's formatting check
# sits in front of.
#
#     scripts/preflight.sh
#
# Needs nothing beyond the Rust toolchain and python3: the oracle tests read
# committed fixtures, and nothing here touches Sage. The Sage-side harnesses
# (scripts/check_*.py) are separate, and run when the subsystem they oracle
# changes — except check_panics_documented.py, check_doc_sentences.py,
# check_spelling.py, check_figures.py, check_links.py and
# check_sage_guards.py, which read sources only and are run here because the
# invariants they pin decay on any commit that adds a `pub fn`, writes a doc
# comment, moves a file, or adds a script.

set -e

step() { printf '\n== preflight: %s\n' "$1"; }

step "cargo fmt --all --check"
cargo fmt --all --check

step "public panics documented"
python3 "$(dirname "$0")/check_panics_documented.py"

step "item-doc sentence discipline"
python3 "$(dirname "$0")/check_doc_sentences.py"

step "one spelling"
python3 "$(dirname "$0")/check_spelling.py"

step "no aphorisms, no metaphors"
python3 "$(dirname "$0")/check_figures.py"

step "markdown links resolve"
python3 "$(dirname "$0")/check_links.py"

step "Sage harnesses refuse a live backend"
python3 "$(dirname "$0")/check_sage_guards.py"

step "cargo test (default features)"
cargo test --quiet

# Its own target directory, which is worth 5x. A feature set change
# invalidates everything in a shared one: measured 9:44 for this suite after
# a default-feature run and 0:56 for the same suite warm in its own
# directory, so a shared directory also left the next `cargo test` or
# `cargo build --features python` to pay another nine minutes
# (docs/record/README.md, "Shipping it").
step "cargo test --features bignum"
CARGO_TARGET_DIR=${CARGO_TARGET_DIR:-target}/bignum cargo test --quiet --features bignum

step "clean"
