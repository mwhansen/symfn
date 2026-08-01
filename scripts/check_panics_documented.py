"""Demand a `# Panics` section from every public function that can panic.

    python3 scripts/check_panics_documented.py

Two propositions, both from `docs/policies/failure.md`:

  * **R2, and `docs/style.md`'s pre-ship checklist** — a `pub fn` that can
    panic says so, in the section rustdoc renders, not in prose a reader has to
    notice. An undocumented panic is how a refusal the design intended becomes
    a surprise the caller publishes past (R9).
  * **No bare `.unwrap()` in `src/`.** An `unwrap` is a claim that nothing
    checks; where the claim is true it costs one `expect` to say what it is,
    and where it is false it was a defect. Both readings were found live during
    the audit — see `docs/record/failure-and-overflow.md`.

Needs no toolchain: it reads the sources. `scripts/preflight.sh` runs it, which
is the point — this invariant was established once and silently broken by the
next merge that added a `pub fn`, twice, before it was pinned.

**`src/python.rs` is excluded** and has its own, stronger pin
(`scripts/check_python_boundary.py`): nothing there may panic *at all*, so
"documented its panic" would be the wrong bar.

This is a lint, not a proof. It cannot tell a reachable wall from a
proven-unreachable invariant — that judgment is per function and lives in the
prose. What it can tell is that nobody wrote the prose.
"""

import pathlib
import re
import sys

SKIP = {"python.rs"}

PANICS = re.compile(
    r"panic!|\.unwrap\(\)|\.expect\(|unreachable!|todo!|unimplemented!"
    r"|\bassert!|\bassert_eq!|\bassert_ne!"
)
# `debug_assert*` is deliberately not in PANICS: it does not fire in release, so
# it documents an invariant rather than declaring a failure mode. A *public*
# precondition guarded only that way is a different defect, and the audit
# closed the four that existed; this script does not police it.
DEBUG_ASSERT = re.compile(r"\bdebug_assert")
BARE_UNWRAP = re.compile(r"\.unwrap\(\)")
FN = re.compile(
    r"^\s*(pub(?:\(crate\))?\s+)?"
    r"(?:const\s+|unsafe\s+|extern\s+\"[^\"]*\"\s+)*fn\s+(\w+)"
)


def live_lines(lines):
    """Index -> is this line production code (outside every `#[cfg(test)]` item)?

    ⚠️ Masks each `#[cfg(test)]` item individually rather than truncating the
    file at the first one. Truncating is what this did originally, and it is
    wrong twice over: a `#[cfg(test)] fn` *helper* sitting mid-file hid the 396
    lines of `convert.rs` after it — including a `pub fn` — and `frac.rs` put
    390 lines of production code after its test module, none of it scanned.
    The `frac.rs` half is now also held by clippy's `items_after_test_module`
    at deny; the helper half is only held here, because no lint objects to a
    test-only function living beside the code it tests.
    """
    live = [True] * len(lines)
    i = 0
    while i < len(lines):
        if not lines[i].strip().startswith("#[cfg(test)]"):
            i += 1
            continue
        # The attribute plus the item it applies to. An item is either braced
        # (`mod`, `fn`, `impl`) or a one-liner ending in `;` (`use`).
        live[i] = False
        j, depth, opened = i + 1, 0, False
        while j < len(lines):
            live[j] = False
            code = lines[j].split("//")[0]
            depth += code.count("{") - code.count("}")
            if "{" in code:
                opened = True
            if opened and depth <= 0:
                break
            if not opened and code.strip().endswith(";"):
                break
            j += 1
        i = j + 1
    return live


def scan(path):
    """(undocumented pub fns, bare unwraps) outside every `#[cfg(test)]` item."""
    lines = path.read_text().split("\n")
    live = live_lines(lines)

    starts = [
        (i, m.group(1), m.group(2))
        for i, line in enumerate(lines)
        if live[i] and (m := FN.match(line))
    ]

    undocumented, bare = [], []
    seen = set()
    for i, line in enumerate(lines):
        stripped = line.strip()
        if not live[i] or stripped.startswith("//"):
            continue
        if BARE_UNWRAP.search(line):
            bare.append((i + 1, stripped))
        if DEBUG_ASSERT.search(line) or not PANICS.search(line):
            continue
        enclosing = [s for s in starts if s[0] <= i]
        if not enclosing:
            continue
        start, vis, name = enclosing[-1]
        if (vis or "").strip() != "pub" or name in seen:
            continue
        seen.add(name)
        # The doc block is the run of `///` and attribute lines above `fn`.
        j, doc = start - 1, []
        while j >= 0 and (
            lines[j].strip().startswith("///") or lines[j].strip().startswith("#[")
        ):
            doc.append(lines[j])
            j -= 1
        if "# Panics" not in "\n".join(doc):
            undocumented.append((start + 1, name))
    return undocumented, bare


def main():
    root = pathlib.Path(__file__).resolve().parent.parent / "src"
    undocumented, bare, files = [], [], 0
    for path in sorted(root.rglob("*.rs")):
        if path.name in SKIP:
            continue
        files += 1
        rel = path.relative_to(root.parent)
        u, b = scan(path)
        undocumented += [f"{rel}:{ln} `{n}` can panic with no `# Panics`" for ln, n in u]
        bare += [f"{rel}:{ln} bare `.unwrap()`: {src}" for ln, src in b]

    problems = undocumented + bare
    if problems:
        print(f"FAIL: {len(problems)} problem(s) across {files} files\n")
        for p in problems:
            print(f"  {p}")
        print(
            "\nEither document the panic (name the requirement or the wall), "
            "restructure so\nit cannot arise, or give the `unwrap` an `expect` "
            "naming the invariant it rests on.\ndocs/policies/failure.md R2 and "
            "R9 are the rules; the audit chapter in\ndocs/record/"
            "failure-and-overflow.md has the worked examples."
        )
        return 1

    print(
        f"ok: {files} files, every public panic documented, no bare `.unwrap()` "
        "(python.rs has its own pin)"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
