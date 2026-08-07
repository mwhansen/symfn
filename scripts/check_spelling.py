"""Hold `src/` to one spelling, the American one.

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
    "atalogue": "atalog",
    "avour": "avor",
    "ecognis": "ecogniz",
    "efence": "efense",
    "ehaviour": "ehavior",
    "eighbour": "eighbor",
    "emois": "emoiz",  # memoise
    "eneralis": "eneraliz",
    "erialis": "erializ",  # serialise, materialise
    "haracteris": "haracteriz",
    "icence": "icense",
    "ocalis": "ocaliz",
    "olour": "olor",
    "ormalis": "ormaliz",
    "pecialis": "pecializ",
    "ptimis": "ptimiz",
    "tabilis": "tabiliz",
    "tandardis": "tandardiz",
    "ummaris": "ummariz",
}

# ⚠️ An `-is` stem fires only before a British inflection. Without this, the
# stem `arallelis` swallows the `-ism` of `available_parallelism` — a std API
# call, in code, renamed to nothing by a spelling lint. `-ise` is a suffix, not
# a substring, and the table has to say so.
# Every Rust prose surface in the tree. `docs/` obeys the same rule and is not
# scanned here: this script reads sources, like the rest of what
# `scripts/preflight.sh` runs, and a prose linter for markdown is a different
# tool with a different false-positive profile.
ROOTS = ("src", "tests", "examples", "benches")

PATTERN = re.compile(
    "|".join(
        stem + ("(?=e|ing|ation|able)" if stem.endswith("is") else "")
        for stem in sorted(STEMS, key=len, reverse=True)
    )
)
# A URL is somebody else's spelling and not ours to correct.
URL = re.compile(r"https?://\S+")


def word_at(line, pos):
    """The whole word containing `pos`, for a message that names it."""
    lo = pos
    while lo > 0 and (line[lo - 1].isalpha() or line[lo - 1] == "_"):
        lo -= 1
    hi = pos
    while hi < len(line) and (line[hi].isalpha() or line[hi] == "_"):
        hi += 1
    return line[lo:hi]


def fix_line(line):
    """Return `(new_line, [(british, american), ...])`.

    URLs are masked rather than skipped, because a line can hold a link and a
    sentence, and the sentence still has to obey the rule.
    """
    holes = [m.span() for m in URL.finditer(line)]

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


def main():
    fix = "--fix" in sys.argv[1:]
    repo = pathlib.Path(__file__).resolve().parent.parent
    total, files = 0, 0
    for path in sorted(p for r in ROOTS for p in (repo / r).rglob("*.rs")):
        lines = path.read_text(encoding="utf-8").splitlines(keepends=True)
        touched = False
        for i, line in enumerate(lines):
            new, changes = fix_line(line)
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
        print("spelling: one spelling in the Rust tree, and it is American")
        return 0
    verb = "fixed" if fix else "found"
    print(f"\nspelling: {verb} {total} British spellings in {files} files")
    if not fix:
        print("re-run with --fix to rewrite them")
    return 0 if fix else 1


if __name__ == "__main__":
    sys.exit(main())
