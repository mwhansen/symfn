"""Every markdown link to a path in the tree points at something that exists.

    python3 scripts/check_links.py

The house rule is that pointers are greppable file paths, because section
numbers drift and nothing checks them (``docs/style.md``, "Specs, and how
they end"). A greppable path only makes drift findable; this is what makes it
found. The scope is every ``.md`` file the spelling gate reads — ``README.md``,
``CLAUDE.md``, ``docs/`` and ``docsite/`` — and the check is existence: a
relative link target, resolved from the file that carries it, is a file or
directory in the tree. External URLs and same-page anchors are out of scope,
and a fragment on a path link (``file.md#section``) is stripped rather than
verified, since markdown anchors are a rendering detail no two hosts agree
on.

Fenced blocks are skipped: the layout tree and quoted material are data, not
links. `sphinx -b html -W` also fails on a broken cross-page reference inside
``docsite/``, but only there, and only for references Sphinx resolves; this
check covers the tree's markdown uniformly, links into ``src/`` and
``scripts/`` included.
"""

import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

DOC_ROOTS = ("docs", "docsite")
DOC_FILES = ("README.md", "CLAUDE.md", "CONTRIBUTING.md", "SECURITY.md", "CODE_OF_CONDUCT.md")

#: An inline markdown link or image, `[text](target)`. Reference-style links
#: are not in the tree's markdown, so they are not parsed. The lookbehind
#: rejects an indexing expression glued to an identifier — `p[2](q*p[1])` in
#: a plethysm example is not a link, and a real link never follows a word
#: character directly.
LINK = re.compile(r"(?<![\w\]`])(?:!\[|\[)[^\]]*\]\(([^)\s]+)\)")
#: Code spans are masked before parsing: `s[6](s[6])` quoted in prose is an
#: expression, not a link. Double-backtick first, as in check_spelling.py.
CODE_SPAN = re.compile(r"``.*?``|`[^`\n]*`")
FENCE = re.compile(r"^\s*(```|~~~)")
#: Targets that are not paths in this tree.
EXTERNAL = re.compile(r"^(https?:|mailto:)|^#")


def sources():
    """Every scanned file, sorted for stable output."""
    files = [ROOT / name for name in DOC_FILES]
    for root in DOC_ROOTS:
        files += (ROOT / root).rglob("*.md")
    return sorted(p for p in files if "_build" not in p.parts)


def main():
    """Check every link in every scanned file. Returns an exit code."""
    broken = []
    checked = 0
    for path in sources():
        in_fence = False
        for lineno, line in enumerate(path.read_text().splitlines(), 1):
            if FENCE.match(line):
                in_fence = not in_fence
                continue
            if in_fence:
                continue
            line = CODE_SPAN.sub(lambda m: " " * len(m.group()), line)
            for match in LINK.finditer(line):
                target = match.group(1)
                if EXTERNAL.match(target):
                    continue
                checked += 1
                plain = target.split("#")[0]
                if not plain or not (path.parent / plain).exists():
                    rel = path.relative_to(ROOT)
                    broken.append(f"{rel}:{lineno}: {target}")

    if broken:
        print("\n".join(broken))
        print(f"\nlinks: {len(broken)} of {checked} targets do not resolve")
        return 1
    print(f"links: {checked} targets resolve, over {len(sources())} files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
