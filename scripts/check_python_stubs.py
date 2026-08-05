"""Hold `symfn.pyi` to the module it describes, name for name and arity for arity.

    cargo build --features python
    python3 scripts/check_python_stubs.py target/debug/libsymfn.dylib

The proposition: **the stub file is the supported surface, not a description
of it.** `docs/policies/python.md` P10 makes membership a decided list rather
than an accumulated one, and a stub file is only that list if nothing can be
exported without appearing here — so a new `#[pyfunction]` that nobody stubs
fails this script rather than shipping unlisted.

Three ways the two can disagree, and this checks all three:

  * **exported but not stubbed** — the failure P10 exists against. A caller's
    editor would not know the function exists, and the "deliberate list" claim
    would be false the moment it happened.
  * **stubbed but not exported** — a name a caller can write, autocomplete,
    and type-check successfully, that raises `AttributeError` at run time.
  * **parameter names disagree** — the quieter one, and the reason this
    compares signatures rather than just names. PyO3 exports every argument as
    keyword-callable, so a parameter name is contract, not decoration: a stub
    saying `mu` where the module says `nu` type-checks a call that fails.

Needs no Sage, and no maturin — it takes the artifact `cargo build` leaves
behind, the way `scripts/check_python_boundary.py` does.

Not checked here: whether the *types* in the stub are right. Nothing at this
boundary carries a Python type annotation to compare against — PyO3 erases
them into extraction code — so the types are hand-maintained and reviewed,
while the names and arities are mechanical and therefore checkable.
"""

import ast
import importlib.machinery
import importlib.util
import pathlib
import sys

STUBS = pathlib.Path(__file__).resolve().parent.parent / "symfn.pyi"


def load(path):
    """Import the built cdylib under the name its `#[pymodule]` declares."""
    loader = importlib.machinery.ExtensionFileLoader("symfn", path)
    spec = importlib.util.spec_from_loader("symfn", loader)
    module = importlib.util.module_from_spec(spec)
    loader.exec_module(module)
    return module


def stub_signatures():
    """`{name: [parameter, ...]}` for every `def` at the top level of the stub."""
    tree = ast.parse(STUBS.read_text())
    out = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef):
            args = node.args
            out[node.name] = [a.arg for a in args.posonlyargs + args.args + args.kwonlyargs]
    return out


def module_signatures(mod):
    """The same, read from PyO3's `__text_signature__`.

    Every entry point is a plain function, so a missing text signature means
    PyO3 stopped emitting one rather than that the function takes no
    arguments — worth failing on rather than reading as `[]`.
    """
    out = {}
    for name in dir(mod):
        if name.startswith("_"):
            continue
        obj = getattr(mod, name)
        if not callable(obj):
            continue
        sig = getattr(obj, "__text_signature__", None)
        if sig is None:
            out[name] = None
            continue
        params = []
        for part in sig.strip("()").split(","):
            part = part.strip().split("=")[0].strip()
            if part and part not in ("/", "*"):
                params.append(part)
        out[name] = params
    return out


def main():
    if len(sys.argv) != 2:
        raise SystemExit(f"usage: {sys.argv[0]} <path to built symfn extension>")
    mod = load(sys.argv[1])
    stubs, exports = stub_signatures(), module_signatures(mod)

    failures = []
    for name in sorted(set(exports) - set(stubs)):
        failures.append(f"{name} is exported but has no stub in symfn.pyi")
    for name in sorted(set(stubs) - set(exports)):
        failures.append(f"{name} is stubbed in symfn.pyi but not exported")
    for name in sorted(set(stubs) & set(exports)):
        if exports[name] is None:
            failures.append(f"{name} has no __text_signature__; PyO3 stopped emitting one")
        elif stubs[name] != exports[name]:
            failures.append(
                f"{name} parameters disagree: stub says "
                f"({', '.join(stubs[name])}), module says ({', '.join(exports[name])})"
            )

    if not hasattr(mod, "__version__"):
        failures.append("the module exports no __version__")
    elif "__version__: str" not in STUBS.read_text():
        failures.append("symfn.pyi does not declare __version__")

    if not (mod.__doc__ or "").strip():
        failures.append("the #[pymodule] carries no docstring (docs/policies/python.md, P11)")

    if failures:
        print(f"FAIL: {len(failures)} problem(s)\n")
        for f in failures:
            print(f"  {f}")
        raise SystemExit(1)
    print(
        f"ok: symfn.pyi and the module agree on {len(stubs)} functions, "
        f"names and parameters; __version__ = {mod.__version__}"
    )


if __name__ == "__main__":
    main()
