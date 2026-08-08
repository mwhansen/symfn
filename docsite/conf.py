"""Sphinx configuration for the rendered reference at symfn.readthedocs.io.

The site documents both layers of the Python surface from the objects
themselves — the compiled module's `__doc__` is what PyO3 ships from the `///`
on each `#[pyfunction]`, and the convenience layer's is ordinary Python — so
there is one copy of every sentence and `help()` and the website cannot drift.

That decision has one consequence worth stating, because it is the reason for
`markdown_docstrings` below. The docstrings are **Markdown**: ``# Raises``
headings, single-backtick code spans, ` ```text ` fences holding doctests
(``docs/policies/python.md``, P11). Sphinx hands autodoc's output to the RST
parser, which reads a `#` line as a comment, a single backtick as a title
reference, and a fence as a paragraph of literal backticks. The hook below
translates the house form to RST before the parser sees it, so the docstrings
stay readable as text in `help()` — which is where most callers meet them —
rather than being rewritten into a markup language for the website's benefit.

The narrative pages are Markdown too, parsed by MyST.
"""

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Read the Docs installs the wheel, so `symfn` imports with its compiled half in
# place and the source tree must stay off the path — `python/symfn` without a
# built extension imports as far as `from .symfn import *` and then fails. A
# local build that has staged the cdylib into the package is the other case, and
# the presence of that file is what tells the two apart.
if list((ROOT / "python" / "symfn").glob("symfn*.so")):
    sys.path.insert(0, str(ROOT / "python"))

project = "symfn"
author = "Mike Hansen"
copyright = "Mike Hansen"  # noqa: A001  Sphinx names this option

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "myst_parser",
]

exclude_patterns = ["_build"]

html_theme = "furo"
html_title = "symfn"

intersphinx_mapping = {"python": ("https://docs.python.org/3", None)}

myst_enable_extensions = ["deflist", "fieldlist"]

autodoc_member_order = "bysource"
autoclass_content = "both"
# Set once here rather than per directive: `scripts/check_docs_complete.py`
# requires every public method and operator to reach a page, and a per-page
# list drifts out of agreement with it one class at a time.
autodoc_default_options = {
    "members": True,
    "special-members": (
        "__call__,__getitem__,__repr__,__eq__,__len__,__bool__,__iter__,"
        "__add__,__sub__,__mul__,__pow__"
    ),
}

# --- the house Markdown, translated for the RST parser -----------------------

#: A fenced block, with the language tag captured.
FENCE = re.compile(
    r"^(?P<indent> *)```(?P<lang>\w*)\n(?P<body>.*?)^(?P=indent)```",
    re.MULTILINE | re.DOTALL,
)
#: An ATX heading, which the docstrings use for `# Raises`.
#
# The trailing class is `[ \t]*`, not `\s*`: `\s` matches newlines, so a greedy
# `\s*$` swallows the blank line after the heading and the rubric then runs
# straight into its paragraph — which docutils reports as "explicit markup ends
# without a blank line", once per entry point.
HEADING = re.compile(r"^(?P<indent> *)#+ +(?P<title>.+?)[ \t]*$", re.MULTILINE)
#: A single-backtick code span that is not part of a double-backtick one.
#
# The span may contain a newline, because prose wraps at 80 columns and a long
# one lands across the wrap: `stanley_table` opens with `|λ| = |μ| = k` broken
# over two lines. An earlier `[^`\n]+` did not match those, so the pair kept
# its single backticks while every span around it became double, and docutils
# reported an inline literal start-string with no end. What it must still not
# cross is a blank line — a span that ran to the next paragraph would pair with
# whatever backtick it found there — so the newline is admitted only when the
# line after it has something on it.
CODE_SPAN = re.compile(r"(?<!`)`((?:[^`\n]|\n(?![ \t]*\n))+)`(?!`)")


def _fence_to_block(match):
    """Rewrite a fenced block as an RST literal block.

    A fence holding a doctest becomes a `pycon` block, which is what makes the
    `>>>` lines render as a session rather than as a paragraph.
    """
    indent, lang, body = match["indent"], match["lang"], match["body"]
    language = "pycon" if ">>>" in body else (lang or "text")
    inner = "".join(
        (indent + "   " + line if line.strip() else line)
        for line in body.splitlines(keepends=True)
    )
    return f"{indent}.. code-block:: {language}\n\n{inner}"


def _literal(match):
    """Rewrite a Markdown code span as an RST inline literal.

    RST requires whitespace or punctuation after a literal's closing quotes, so
    a span a word runs into — ``` `int`s ``` is the one in the module doc —
    needs the escaped space that separates them without printing anything.
    """
    body, after = match.group(1), match.string[match.end() : match.end() + 1]
    return f"``{body}``" + ("\\ " if after.isalnum() else "")


def markdown_docstrings(app, what, name, obj, options, lines):
    """Translate the house docstring Markdown into RST, in place.

    Four rules cover the whole form:

    * a fenced block becomes a literal block, `pycon` when it holds a doctest;
    * an ATX heading — `# Raises` is the only one in use — becomes a rubric,
      which renders as a bold run-in title rather than splitting the page's
      section tree;
    * a single-backtick span becomes a double-backtick literal, so `s_lambda`
      does not become a dangling title reference;
    * a bare vertical bar is escaped, because the mathematics writes `|lambda|`
      for the size of a partition and RST reads that as a substitution.

    The order matters. Fences and code spans are set aside before anything else
    runs, so a `>>>` line full of backticks, hashes and bars is left alone, and
    so the bar inside a literal stays a bar rather than gaining a backslash.
    """
    text = "\n".join(lines)
    stashed = []

    def stash(rendered):
        stashed.append(rendered)
        return f"\x00{len(stashed) - 1}\x00"

    text = FENCE.sub(lambda m: stash(_fence_to_block(m)), text)
    text = HEADING.sub(lambda m: f"{m['indent']}.. rubric:: {m['title']}", text)
    text = CODE_SPAN.sub(lambda m: stash(_literal(m)), text)
    text = text.replace("|", "\\|")
    text = re.sub(r"\x00(\d+)\x00", lambda m: stashed[int(m.group(1))], text)

    lines[:] = text.split("\n")


def setup(app):
    """Register the docstring translation."""
    app.connect("autodoc-process-docstring", markdown_docstrings)
    return {"parallel_read_safe": True}
