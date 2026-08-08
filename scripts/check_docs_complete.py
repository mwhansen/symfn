"""Fail when a supported name is missing from the rendered documentation.

    cargo build --features python
    python3 -m sphinx -b html docsite docsite/_build/html
    python3 scripts/check_docs_complete.py

The proposition: **the supported surface is a deliberate list
(``docs/policies/python.md``, P10), and the website is where an outside reader
meets it — so a name on the list that no page documents is a hole, not an
omission someone will notice.** Sphinx does not check this on its own. It warns
about a reference that does not resolve; it says nothing about an entry point
nobody wrote an `automodule` for, which is the failure that actually happens
when the surface grows.

What is compared is `symfn.__all__` against the objects Sphinx recorded in
`objects.inv` — the same inventory `intersphinx` publishes, so a name that
passes here is also a name another project can link to. Public methods and
properties of the convenience layer's classes are checked too: a class that is
documented with half its methods missing passes an `autoclass` and fails a
reader.

The check reads a built site rather than building one, so it stays cheap enough
for `scripts/preflight.sh` and does not need Sphinx installed to run — the
inventory is a zlib stream with a four-line plain header.
"""

import pathlib
import re
import sys
import zlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
INVENTORY = ROOT / "docsite" / "_build" / "html" / "objects.inv"

#: What counts as the surface: the names a reader can type. Operators are
#: rendered when they carry a docstring — `Sym.__mul__` states the routing and
#: the Littlewood-Richardson convention, so it should be on the page — but they
#: are not *required* here, because the rule this file enforces has to be the
#: same one `scripts/check_convenience_docs.py` enforces for examples, and that
#: one is "the names a reader can type" too. Requiring dunders in one gate and
#: not the other is how the two drift.
def public_attrs(holder):
    """The non-underscore methods and properties `holder` defines."""
    for attr, member in vars(holder).items():
        if attr.startswith("_"):
            continue
        if callable(member) or isinstance(member, property):
            yield attr


def inventory(path):
    """The set of object names in a Sphinx `objects.inv`.

    The format is four plain-text header lines followed by one zlib stream of
    `name domain:role priority uri dispname` records.
    """
    raw = path.read_bytes()
    start = 0
    for _ in range(4):
        start = raw.index(b"\n", start) + 1
    body = zlib.decompress(raw[start:]).decode("utf-8")
    names = set()
    for line in body.splitlines():
        match = re.match(r"(\S+)\s+\S+:\S+\s+\S+\s+\S+\s+.*", line)
        if match:
            names.add(match.group(1))
    return names


def expected(package):
    """Every name the documentation must carry.

    Each entry is the set of paths any one of which satisfies it. A method of
    an exported *instance* — `symfn.macdonald.P`, `symfn.s.__call__` — is
    documented under the class that defines it, whose own path is what
    `autoclass` puts in the inventory, so both spellings count.
    """
    wanted = []
    for name in package.__all__:
        if name == "__version__":
            continue
        obj = getattr(package, name)
        holder = obj if isinstance(obj, type) else type(obj)
        # The convenience classes set `__module__ = "symfn"` so a traceback
        # prints the public name, which means the attribute does not say where
        # `autoclass` filed them. The module that actually holds the class does.
        paths = {
            f"{holder.__module__}.{holder.__qualname__}",
            f"{_defining_module(holder)}.{holder.__qualname__}",
        }
        wanted.append({f"symfn.{name}"} | paths)
        if not _is_ours(holder):
            continue
        for attr in public_attrs(holder):
            wanted.append(
                {f"symfn.{name}.{attr}"} | {f"{p}.{attr}" for p in paths}
            )
    return wanted


def _defining_module(holder):
    """The name of the module whose namespace actually holds `holder`."""
    for name, module in list(sys.modules.items()):
        if name.startswith("symfn.") and (
            getattr(module, holder.__qualname__, None) is holder
        ):
            return name
    return holder.__module__


def _is_ours(holder):
    """Whether `holder` is a convenience-layer type rather than a builtin."""
    return _defining_module(holder).startswith("symfn.")


def main():
    """Compare the inventory to the surface. Returns an exit code."""
    if not INVENTORY.exists():
        print(f"no built site at {INVENTORY}")
        print("build one: python3 -m sphinx -b html docsite docsite/_build/html")
        return 1

    sys.path.insert(0, str(ROOT / "scripts"))
    from check_convenience_docs import stage

    if stage() is None:
        print("no extension module; build one: cargo build --features python")
        return 1
    sys.path.insert(0, str(ROOT / "python"))
    import symfn

    documented = inventory(INVENTORY)
    wanted = expected(symfn)

    # The contract layer is documented under its own module path, since that is
    # where `automodule` puts it; the flat alias is the same object.
    missing = []
    for paths in wanted:
        spellings = paths | {p.replace("symfn.", "symfn.symfn.", 1) for p in paths}
        if not spellings & documented:
            missing.append(min(paths))
    if missing:
        print(
            f"docs: {len(missing)} of {len(wanted)} supported names appear on "
            "no page"
        )
        for name in sorted(missing):
            print(f"  {name}")
        return 1
    print(f"docs: all {len(wanted)} supported names are documented")
    return 0


if __name__ == "__main__":
    sys.exit(main())
