"""Execute every example in the extension module's docstrings.

    cargo build --features python
    python3 scripts/check_python_docs.py target/debug/libsymfn.dylib

The proposition: **an example in a Python docstring is a pin, not an
illustration.** `docs/policies/python.md` P11 requires each entry point to
carry a convention-distinguishing example, and an example nothing runs stops
being true silently — the module docstring's own `schur_multiply` line showed
lists where the boundary returns tuples, which is the exact claim P1 makes
about the outbound type, printed wrong on the surface's front page.

Examples live in ` ```text ` fences, the house form: cargo's doctest runner
does not compile them as Rust, and this script does run them. Each fence that
contains a `>>>` line is one doctest, run with `symfn` bound to the module
loaded here, so a fence is self-contained and the prose between fences is not
parsed as expected output.

Needs no Sage and no maturin — it loads the artifact `cargo build` leaves
behind, as `scripts/check_python_boundary.py` and
`scripts/check_python_stubs.py` do. `scripts/preflight.sh` cannot run it,
since it builds only the default features.

Completeness is checked too: an entry point with no example fails, because
P11 requires one and a surface that is half-pinned decays back to none. What
is **not** checked is whether an example distinguishes anything — that a value
separates the shipped convention from its rivals is a judgment, and the
reviewer makes it.
"""

import doctest
import importlib.machinery
import importlib.util
import pathlib
import re
import sys

FENCE = re.compile(r"^ *```text\n(.*?)^ *```", re.MULTILINE | re.DOTALL)
DEFAULT = "target/debug/libsymfn.dylib"


def load(path):
    """Import the built cdylib under the name its `#[pymodule]` declares."""
    loader = importlib.machinery.ExtensionFileLoader("symfn", path)
    spec = importlib.util.spec_from_loader("symfn", loader)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["symfn"] = mod
    loader.exec_module(mod)
    return mod


def documented(mod):
    """Yield (name, docstring) for the module and every exported callable."""
    yield "symfn", mod.__doc__ or ""
    for name in sorted(n for n in dir(mod) if not n.startswith("_")):
        yield name, getattr(mod, name).__doc__ or ""


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT
    if not pathlib.Path(path).exists():
        print(f"no extension module at {path}")
        print("build one: cargo build --features python")
        return 1

    mod = load(path)
    parser, runner = doctest.DocTestParser(), doctest.DocTestRunner()
    examples = total = covered = 0
    bare = []
    for name, doc in documented(mod):
        fences = [m.group(1) for m in FENCE.finditer(doc) if ">>>" in m.group(1)]
        if name != "symfn":
            total += 1
            covered += bool(fences)
            if not fences:
                bare.append(name)
        for i, block in enumerate(fences):
            test = parser.get_doctest(block, {"symfn": mod}, f"{name}[{i}]", None, 0)
            examples += len(test.examples)
            runner.run(test)

    failed = runner.summarize(verbose=False).failed
    if failed:
        print(f"\npython docs: {failed} of {examples} examples failed")
        return 1
    if bare:
        print(f"python docs: {len(bare)} of {total} entry points carry no "
              "example, which `docs/policies/python.md` P11 requires")
        for name in bare:
            print(f"  {name}")
        return 1
    print(f"python docs: {examples} examples pass, over all {total} entry "
          "points")
    return 0


if __name__ == "__main__":
    sys.exit(main())
