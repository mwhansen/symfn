"""Hold the tree to one spelling, the American one.

    python3 scripts/check_spelling.py            # report, exit 1 if any
    python3 scripts/check_spelling.py --fix      # rewrite them in place

`docs/style.md` has asked for American English since it was written (delta 6,
"One spelling, where today there are two"), and nothing checked. Two years of
commits later `src/` held 70 British forms across 22 files — `labelled` 22
times, `normalisation` 11, `optimisation` 12 — and the sweep that was supposed
to fix them caught 4, because it fixed the ones somebody happened to read.
A rule with no gate decays at the rate the tree grows.

**Every line is scanned, not only the doc comments.** The words turn up in
panic messages, `assert!` messages, and test names, which `docs/style.md`
governs as prose exactly like rustdoc. No public identifier in `src/` carries
one; if one ever must — an external API spelled the other way — this script
has no escape hatch and wants a real one rather than a skipped file.

**Markdown is scanned too, with code masked.** For its first year this script
read `.rs` only, on the reasoning that a prose linter for markdown is a
different tool with a different false-positive profile. It is not: the stems
are the same and the rule is the same, and the exemption let 30 British forms
accumulate in `docs/` while `src/` held none. In `.md` files, fenced blocks and
inline code spans are masked the way URLs are, which is what makes the two
genuine exemptions work without a skip list — Sage's
`to_labelling_area_sequence_pair` is an API name and stays in backticks, and
`docs/style.md` quotes `normalise` as the form to avoid. Prose in a fence is
invisible to the lint; the tree has one such line, the `dyck.rs` entry in the
README's layout tree, and it was fixed by hand.

Matching is by stem, and each stem deliberately omits the word's first letter
(`ormalis`, not `normalis`), so `Normalisation` and `normalisation` are one
rule and no case table is needed. An all-caps `NORMALISATION` would slip
through; none exists, and inventing the machinery for none is how a lint gets
too clever to trust.

This is a lint, not a proof. It knows the stems it was told about. A British
form outside the table is invisible to it, so a new one arrives the same way
these did — which is why `--fix` prints what it changed rather than working
silently.
"""

import pathlib
import re
import sys

# Stem -> replacement, with the leading letter of each word omitted so that
# capitalization needs no second table. `erialis` covers both `serialise` and
# `materialise`; the shared tail is the point.
STEMS = {
    "abell": "abel",  # labelled, labelling
    "arallelis": "aralleliz",
    "arameteris": "arameteriz",
    "atalogue": "atalog",
    "avour": "avor",
    "ecognis": "ecogniz",
    "efence": "efense",
    "ehaviour": "ehavior",
    "eighbour": "eighbor",
    "emois": "emoiz",  # memoise
    "enalis": "enaliz",
    "eneralis": "eneraliz",
    "erialis": "erializ",  # serialise, materialise
    "haracteris": "haracteriz",
    "icence": "icense",
    "mortis": "mortiz",  # amortise
    "ocalis": "ocaliz",
    "olour": "olor",
    "ormalis": "ormaliz",
    "pecialis": "pecializ",
    "ptimis": "ptimiz",
    "riticis": "riticiz",
    "rioritis": "rioritiz",
    "tabilis": "tabiliz",
    "tandardis": "tandardiz",
    "ummaris": "ummariz",
}

# ⚠️ An `-is` stem fires only before a British inflection. Without this, the
# stem `arallelis` swallows the `-ism` of `available_parallelism` — a std API
# call, in code, renamed to nothing by a spelling lint. `-ise` is a suffix, not
# a substring, and the table has to say so.
# Every prose surface in the tree: Rust sources, then the markdown that
# `docs/style.md` governs by the same rule.
ROOTS = ("src", "tests", "examples", "benches")
DOC_ROOTS = ("docs",)
DOC_FILES = ("README.md", "CLAUDE.md")

PATTERN = re.compile(
    "|".join(
        stem + ("(?=e|ing|ation|able)" if stem.endswith("is") else "")
        for stem in sorted(STEMS, key=len, reverse=True)
    )
)
# A URL is somebody else's spelling and not ours to correct.
URL = re.compile(r"https?://\S+")
# So is anything in a code span: an API name, or a form quoted to be rejected.
# The double-backtick alternative comes first so it is not read as two empty
# spans — `docs/style.md` uses it for headings that contain a backtick.
CODE_SPAN = re.compile(r"``.*?``|`[^`\n]*`")
FENCE = re.compile(r"^\s*(```|~~~)")


def word_at(line, pos):
    """The whole word containing `pos`, for a message that names it."""
    lo = pos
    while lo > 0 and (line[lo - 1].isalpha() or line[lo - 1] == "_"):
        lo -= 1
    hi = pos
    while hi < len(line) and (line[hi].isalpha() or line[hi] == "_"):
        hi += 1
    return line[lo:hi]


def fix_line(line, code_spans=False):
    """Return `(new_line, [(british, american), ...])`.

    URLs are masked rather than skipped, because a line can hold a link and a
    sentence, and the sentence still has to obey the rule. In markdown, code
    spans are masked the same way and for the same reason.
    """
    holes = [m.span() for m in URL.finditer(line)]
    if code_spans:
        holes += [m.span() for m in CODE_SPAN.finditer(line)]

    def masked(pos):
        return any(lo <= pos < hi for lo, hi in holes)

    changes = []
    out, last = [], 0
    for m in PATTERN.finditer(line):
        if masked(m.start()):
            continue
        before = word_at(line, m.start())
        after = before.replace(m.group(), STEMS[m.group()])
        changes.append((before, after))
        out.append(line[last : m.start()])
        out.append(STEMS[m.group()])
        last = m.end()
    out.append(line[last:])
    return "".join(out), changes


def sources(repo):
    """Every scanned file, paired with whether it is markdown."""
    for path in sorted(p for r in ROOTS for p in (repo / r).rglob("*.rs")):
        yield path, False
    docs = [p for r in DOC_ROOTS for p in (repo / r).rglob("*.md")]
    docs += [repo / name for name in DOC_FILES]
    for path in sorted(docs):
        yield path, True


def main():
    fix = "--fix" in sys.argv[1:]
    repo = pathlib.Path(__file__).resolve().parent.parent
    total, files = 0, 0
    for path, markdown in sources(repo):
        lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
        touched, fenced = False, False
        for i, line in enumerate(lines):
            if markdown and FENCE.match(line):
                fenced = not fenced
                continue
            if markdown and fenced:
                continue
            new, changes = fix_line(line, code_spans=markdown)
            if not changes:
                continue
            total += len(changes)
            touched = True
            rel = path.relative_to(repo)
            for before, after in changes:
                print(f"{rel}:{i + 1}: {before} -> {after}")
            lines[i] = new
        if touched:
            files += 1
            if fix:
                path.write_text("".join(lines), encoding="utf-8")

    if not total:
        print("spelling: one spelling in the tree, and it is American")
        return 0
    verb = "fixed" if fix else "found"
    print(f"\nspelling: {verb} {total} British spellings in {files} files")
    if not fix:
        print("re-run with --fix to rewrite them")
    return 0 if fix else 1


if __name__ == "__main__":
    sys.exit(main())
