"""Execute every example in the convenience layer's docstrings.

    cargo build --features python
    python3 scripts/check_convenience_docs.py

The proposition is `scripts/check_python_docs.py`'s, one layer up: **an example
in a docstring is a pin, not an illustration** (``docs/policies/python.md``,
P11). The difference is the runner. The contract layer's examples sit in
` ```text ` fences so cargo does not compile them as Rust, and that script pulls
them out by hand; the convenience layer is Python all the way down, so its
examples are ordinary doctests — but `doctest`'s stock finder must be steered.
Every class in the layer says `__module__ = "symfn"` so `help()` reads well,
and that lie fails the finder's ownership tests: the module scan drops the
class, the package scan drops its methods, and an instance — the basis
factories, the family namespaces — is never recursed into at all. The layer's
examples sat unexecuted behind exactly that gap once, and rotted there, so the
collection below walks the supported classes explicitly.

Completeness is checked the same way and for the same reason: a public class or
method with no example fails, because a surface that is half-pinned decays back
to none. What is not checked is whether an example *distinguishes* anything —
that a value separates the shipped convention from its rivals is a judgment,
and the reviewer makes it.

Needs no Sage, no maturin and no installed wheel: it puts the built cdylib
where the package expects its compiled half and imports `symfn` from the source
tree, so it runs on a stock interpreter with the Rust artifact alone.
"""

import doctest
import importlib
import pathlib
import shutil
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "python" / "symfn"
CDYLIB = ROOT / "target" / "debug" / "libsymfn.dylib"
CDYLIB_LINUX = ROOT / "target" / "debug" / "libsymfn.so"

#: The pure-Python modules, in the order a reader meets them.
MODULES = ["_bases", "_param", "_sym", "_families", "_schubert"]


class TrustingFinder(doctest.DocTestFinder):
    """A finder that takes every member of the object it is handed.

    Used only on classes whose `__module__` is spoofed to `symfn`: the stock
    ownership test would reject their methods for carrying the defining
    module's globals, and every member of such a class is defined beside it.
    """

    def _from_module(self, module, object):  # noqa: A002  doctest's signature
        return True


def stage():
    """Put the built extension where `symfn/__init__.py` imports it from.

    Returns the path staged, or `None` if no build was found.
    """
    for built in (CDYLIB, CDYLIB_LINUX):
        if built.exists():
            target = PACKAGE / "symfn.so"
            if not target.exists() or target.stat().st_mtime < built.stat().st_mtime:
                shutil.copy2(built, target)
            return target
    return None


def public(package, contract):
    """Yield `(label, object)` for the convenience layer's supported surface.

    The list is `symfn.__all__` minus the contract layer, which
    `scripts/check_python_docs.py` already covers, plus the public methods and
    properties of every class and namespace on it. A helper in a private
    module is not on it: reachability from `symfn` is what "supported" means
    (``docs/policies/python.md``, P10).
    """
    for name in package.__all__:
        if name in contract or name == "__version__":
            continue
        obj = getattr(package, name)
        holder = obj if isinstance(obj, type) else type(obj)
        if isinstance(obj, type) or holder.__module__ == "symfn":
            yield name, obj
        # A constructor is exempt when its class's own docstring shows one.
        # The class doc owns the type's representation and how it is built
        # (`docs/policies/python.md`, P11); a second example under `__init__`
        # would restate it, and the rule is here to catch surfaces nothing
        # demonstrates, not to count.
        built = ">>>" in (holder.__doc__ or "")
        for attr in sorted(vars(holder)):
            if attr.startswith("_") and not (attr == "__init__" and not built):
                continue
            member = vars(holder)[attr]
            if callable(member) or isinstance(member, property):
                yield f"{name}.{attr}", member


def documented(obj):
    """The docstring of a function, method or property."""
    if isinstance(obj, property):
        return obj.fget.__doc__ or ""
    return getattr(obj, "__doc__", "") or ""


def main():
    """Run the doctests, then the completeness check. Returns an exit code."""
    if stage() is None:
        print(f"no extension module at {CDYLIB}")
        print("build one: cargo build --features python")
        return 1

    sys.path.insert(0, str(ROOT / "python"))
    package = importlib.import_module("symfn")

    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    tests, seen = [], set()

    def collect(found):
        for test in found:
            if test.examples and test.name not in seen:
                seen.add(test.name)
                tests.append(test)

    # `extraglobs`, not `globs`: the latter replaces the module's own globals,
    # and a helper's example that calls its neighbor by bare name then fails
    # with a `NameError` that says nothing about the helper.
    finder = doctest.DocTestFinder(exclude_empty=False)
    for name in MODULES + [""]:
        module = package if not name else importlib.import_module(f"symfn.{name}")
        label = "symfn" if not name else f"symfn.{name}"
        collect(finder.find(module, label, extraglobs={"symfn": package}))

    # The stock finder's ownership test loses the layer's classes twice over.
    # Every class here says `__module__ = "symfn"` so `help()` reads well,
    # which makes the module scan reject the class; the package scan accepts
    # it but rejects its methods, whose globals are the defining module's; and
    # a class reachable only as the type of an instance — the basis factories,
    # the four family namespaces — is never reached at all, because the finder
    # does not recurse into instances. So the supported classes are walked
    # directly, with a finder that waives the ownership test: every member of
    # such a class is defined beside it.
    trusting = TrustingFinder(exclude_empty=False)
    contract = {n for n in dir(package.symfn) if not n.startswith("_")}
    holders = {}  # class -> the first public name that reaches it
    for name in package.__all__:
        if name in contract or name == "__version__":
            continue
        obj = getattr(package, name)
        holder = obj if isinstance(obj, type) else type(obj)
        if holder.__module__ == "symfn":
            holders.setdefault(holder, name)
    for holder, name in holders.items():
        collect(
            trusting.find(
                holder, f"symfn.{name}", module=package, extraglobs={"symfn": package}
            )
        )

    examples = sum(len(test.examples) for test in tests)
    for test in tests:
        runner.run(test)

    failed = runner.summarize(verbose=False).failed
    if failed:
        print(f"\nconvenience docs: {failed} of {examples} examples failed")
        return 1

    contract = {n for n in dir(package.symfn) if not n.startswith("_")}
    bare, total = [], 0
    for label, obj in public(package, contract):
        total += 1
        if ">>>" not in documented(obj):
            bare.append(label)
    if bare:
        print(
            f"convenience docs: {len(bare)} of {total} public items carry no "
            "example, which `docs/policies/python.md` P11 requires"
        )
        for label in bare:
            print(f"  {label}")
        return 1

    print(
        f"convenience docs: {examples} examples pass, over all {total} public "
        "items"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
