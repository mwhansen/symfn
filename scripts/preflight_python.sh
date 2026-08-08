#!/bin/sh
# The Python-side gate, in the order that fails fastest.
#
#     scripts/preflight_python.sh
#
# Separate from scripts/preflight.sh because the bars are different. That script
# needs nothing beyond the Rust toolchain and python3, and builds the default
# features; everything here needs `cargo build --features python`, which pulls
# pyo3 and its dependencies, and the last two steps need ruff and Sphinx. Run
# this when the Python surface changes; CI runs it on every push.
#
# What each step pins:
#
#   stubs        symfn.pyi and the module agree, name for name and arity for
#                arity (docs/policies/python.md, P10)
#   boundary     every precondition a caller can violate is a typed exception,
#                never a panic (P8)
#   contract     every entry point's docstring example runs and is true (P11)
#   convenience  every convenience method equals the contract calls it claims
#                to be, and the families hit their classical limits (P4, P7)
#   examples     the convenience layer's own doctests, and that every public
#                item has one (P11)
#   ruff         the lint configured in pyproject.toml, `ANN` included, so a
#                public signature cannot go back to being unannotated
#   mypy         --strict over the layer, reading symfn.pyi for the compiled
#                half: the check that the annotations are *right*, where ruff
#                only checks they are *there*
#   docs         every supported name reaches a rendered page

set -e

here=$(dirname "$0")
root=$here/..
step() { printf '\n== python preflight: %s\n' "$1"; }

case $(uname -s) in
Darwin) lib=$root/target/debug/libsymfn.dylib ;;
*) lib=$root/target/debug/libsymfn.so ;;
esac

step "cargo build --features python"
cargo build --features python --manifest-path "$root/Cargo.toml"

step "stubs agree with the module"
python3 "$here/check_python_stubs.py" "$lib"

step "the boundary raises rather than panicking"
python3 "$here/check_python_boundary.py" "$lib"

step "contract-layer docstring examples"
python3 "$here/check_python_docs.py" "$lib"

step "convenience layer equals the contract layer"
python3 "$here/check_convenience.py"

step "convenience-layer docstring examples"
python3 "$here/check_convenience_docs.py"

if command -v ruff >/dev/null 2>&1 || python3 -c "import ruff" 2>/dev/null; then
	step "ruff"
	(cd "$root" && python3 -m ruff check .)
else
	printf '\n== python preflight: ruff not installed, skipped\n'
fi

if python3 -c "import mypy" 2>/dev/null; then
	step "mypy --strict"
	(cd "$root" && python3 -m mypy)
else
	printf '\n== python preflight: mypy not installed, skipped\n'
fi

if python3 -c "import sphinx, myst_parser, furo" 2>/dev/null; then
	step "the rendered documentation is complete"
	python3 -m sphinx -b html -W "$root/docsite" "$root/docsite/_build/html" >/dev/null
	python3 "$here/check_docs_complete.py"
else
	printf '\n== python preflight: sphinx not installed, skipped\n'
fi

step "clean"
