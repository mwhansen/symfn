"""Execute every Python example on the docsite's pages and the README.

    cargo build --features python
    python3 scripts/check_docsite_docs.py

The proposition is `scripts/check_convenience_docs.py`'s, one surface over:
**an example on a rendered page is a pin, not an illustration**
(``docs/policies/python.md``, P11). The narrative pages repeat values whose
home is a docstring — that is what lets a reader stay on the page — and the
repetition is safe only while both copies run against the library.

A page is one interpreter session: a name imported in an earlier ` ```pycon `
fence is in scope for every later one, which is how the pages read and how a
reader would type them. Sphinx is not involved — the fences are read from the
Markdown source, so this needs the built cdylib and nothing else. The README
is on the page list because its usage block is the same genre: a value a
registry visitor will read as a promise. Its Rust twin runs under
`cargo test` (`src/lib.rs`, the `readme` doctest module).
"""

import doctest
import pathlib
import re
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

sys.path.insert(0, str(ROOT / "scripts"))
from check_convenience_docs import stage  # noqa: E402

#: A ```pycon fence and its body, on the narrative pages' flush-left form.
FENCE = re.compile(r"^```pycon\n(.*?)^```", re.MULTILINE | re.DOTALL)


def main():
    """Run every fence on every page. Returns an exit code."""
    if stage() is None:
        print("no extension module built; run: cargo build --features python")
        return 1
    sys.path.insert(0, str(ROOT / "python"))

    pages = (
        [ROOT / "README.md"]
        + sorted((ROOT / "docsite").glob("*.md"))
        + sorted((ROOT / "docsite" / "api").glob("*.md"))
    )
    parser = doctest.DocTestParser()
    runner = doctest.DocTestRunner(optionflags=doctest.ELLIPSIS)
    examples = fenced_pages = 0
    for page in pages:
        text = page.read_text()
        globs = {}
        found = False
        for index, match in enumerate(FENCE.finditer(text)):
            found = True
            test = parser.get_doctest(
                match.group(1),
                globs,
                f"{page.name}[{index}]",
                str(page),
                text.count("\n", 0, match.start(1)),
            )
            examples += len(test.examples)
            runner.run(test, clear_globs=False)
            # `get_doctest` copies its globs, so the session the next fence
            # continues is the one this run just built.
            globs = test.globs
        fenced_pages += found

    failed = runner.summarize(verbose=False).failed
    if failed:
        print(f"\npage examples: {failed} of {examples} failed")
        return 1
    print(
        f"page examples: {examples} pass, over {fenced_pages} pages "
        f"of {len(pages)}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
