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
#   pointers     no docstring under python/symfn/, and none on an exported
#                item of src/python.rs, names a path in this tree, which a
#                reader of the wheel cannot open (P11)
#   boundary     every precondition a caller can violate is a typed exception,
#                never a panic (P8)
#   marshalling  a value handed in comes back out intact, at the widths and
#                in the shapes the boundary promises (P1)
#   interrupt    a call that runs for seconds stops when Ctrl-C arrives, rather
#                than when it would have finished anyway (docs/policies/
#                failure.md, the cancellation row)
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

step "docstrings point only where a Python reader can go"
python3 "$here/check_python_pointers.py"

step "the boundary raises rather than panicking"
python3 "$here/check_python_boundary.py" "$lib"

step "values cross intact, in the promised shapes"
python3 "$here/check_python_marshalling.py" "$lib"

step "a long call answers Ctrl-C while it runs"
python3 "$here/check_python_interrupt.py" "$lib"

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

# The floor is 2.0, and it is checked rather than assumed: mypy 1.x reads
# `check_basis`'s membership test without narrowing `str` to `Basis` and reports
# the missing `cast` as a return-type error, while 2.x reports that same cast as
# redundant. Nothing satisfies both, and CI runs the newer one. An old checker
# says so here instead of printing a failure that is about itself.
if python3 -c "import mypy" 2>/dev/null; then
	mypy_major=$(python3 -c "import mypy.version; print(mypy.version.__version__.split('.')[0])")
	if [ "$mypy_major" -ge 2 ] 2>/dev/null; then
		step "mypy --strict"
		(cd "$root" && python3 -m mypy)
	else
		printf '\n== python preflight: mypy %s is below the 2.0 this gate needs, skipped\n' \
			"$(python3 -c 'import mypy.version; print(mypy.version.__version__)')"
	fi
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
