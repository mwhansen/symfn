"""Fail when a docstring on the Python surface names a path in this tree.

    python3 scripts/check_python_pointers.py

The proposition: **a pointer in a Python docstring leads somewhere the Python
reader can go** (`docs/policies/python.md`, P11). That reader has the wheel,
`help()`, a stub tooltip and the two published sites — not a checkout — so
`docs/policies/python.md`, `docs/record/llt.md` and `src/llt.rs` are dead ends
for them, and P11 says what to write instead: the fact itself, the docs.rs
module path, or the docsite's Conventions page.

Two surfaces are scanned, because both are what `help()` prints:

  * every string in `python/symfn/*.py` and `symfn.pyi` — the docstrings, and
    the `#:` attribute comments Sphinx picks up;
  * in `src/python.rs`, the `///` block on every `#[pyfunction]` and on the
    `#[pymodule]`, which PyO3 ships as `__doc__`.

What is looked for is `docs/`, `scripts/`, `examples/`, and a `.rs` file name.
Not scanned: plain `#` and `//` comments, and the `///` on private Rust items,
which address a maintainer with a checkout, where the tree-path form is the
right one (`docs/style.md`, "Pointers").

Needs nothing built and no Sage; it reads source.
"""

import pathlib
import re
import sys
import tokenize

ROOT = pathlib.Path(__file__).resolve().parent.parent
PACKAGE = ROOT / "python" / "symfn"
CONTRACT = ROOT / "src" / "python.rs"

# A tree path in any of the forms it has been written: bare, backticked, or
# double-backticked. A `.rs` name is caught with or without a directory —
# `dyck.rs` and `src/llt.rs` both — while `src` alone is not, since several
# entry points have an argument of that name, and `docs.rs` is the site the
# rule points at, not a file.
TREE_PATH = re.compile(
    r"\b(?:docs/[\w./-]+|scripts/[\w./-]+|examples/[\w./-]+|(?!docs\.rs\b)[\w/]+\.rs\b)"
)


def python_offenders(path):
    """`(line, match)` for every tree path in a docstring or `#:` comment."""
    out = []
    with tokenize.open(path) as f:
        for tok in tokenize.generate_tokens(f.readline):
            if tok.type == tokenize.STRING or (
                tok.type == tokenize.COMMENT and tok.string.startswith("#:")
            ):
                for line_offset, line in enumerate(tok.string.splitlines()):
                    for m in TREE_PATH.finditer(line):
                        out.append((tok.start[0] + line_offset, m.group(0)))
    return out


def rust_offenders(path):
    """`(line, match)` for every tree path in a `///` block on an exported item.

    A block is exported when the attributes between it and its item include
    `#[pyfunction` or `#[pymodule`; every other `///` is a private item's and
    is left alone.
    """
    lines = path.read_text().splitlines()
    out = []
    i = 0
    while i < len(lines):
        if not lines[i].lstrip().startswith("///"):
            i += 1
            continue
        start = i
        while i < len(lines) and lines[i].lstrip().startswith("///"):
            i += 1
        j = i
        exported = False
        while j < len(lines) and lines[j].lstrip().startswith("#["):
            attr = lines[j].strip()
            exported |= attr.startswith("#[pyfunction") or attr.startswith("#[pymodule")
            j += 1
        if exported:
            for k in range(start, i):
                for m in TREE_PATH.finditer(lines[k]):
                    out.append((k + 1, m.group(0)))
    return out


def main():
    failures = []
    files = sorted(PACKAGE.glob("*.py")) + [PACKAGE / "symfn.pyi"]
    for path in files:
        for line, hit in python_offenders(path):
            failures.append(f"{path.relative_to(ROOT)}:{line}: names {hit!r}")
    for line, hit in rust_offenders(CONTRACT):
        failures.append(f"{CONTRACT.relative_to(ROOT)}:{line}: names {hit!r}")
    if failures:
        print(f"FAIL: {len(failures)} tree path(s) in Python-facing docstrings\n")
        for f in failures:
            print(f"  {f}")
        print(
            "\nA reader of the wheel cannot open these. State the fact instead, or\n"
            "point at the docs.rs module path or the docsite's Conventions page\n"
            "(docs/policies/python.md, P11)."
        )
        raise SystemExit(1)
    print(
        f"ok: no tree paths in the docstrings of {len(files)} files under "
        f"python/symfn/ or on the exported items of src/python.rs"
    )


if __name__ == "__main__":
    main()
