"""Hold every prose surface to the ban on aphorisms and metaphors.

    python3 scripts/check_figures.py             # report, exit 1 if any
    python3 scripts/check_figures.py --measure REV   # hit rate against REV

`docs/style.md`, "No aphorisms, no metaphors", bans two things: a maxim stated
as if self-evident, and naming a thing as something it is not. Neither is
greppable in general — that is why the rule's own checklist entry called for a
gate only "if one is cheap". This is the cheap part of it.

**The vocabulary is measured, not imagined.** A 2026-08-08 sweep read every
prose surface in the tree and replaced 130 figures. The patterns below are the
ones that recurred, with the count each had when the sweep found it. A word
that appeared once is not here: a lint that fires on prose nobody writes twice
costs more attention than it saves.

Two classes, and the split is what makes the gate usable:

* **Announcements.** "The lesson is", "worth stating as a rule", "the portable
  lesson", "generalizes past". These introduce a maxim, and the maxim is
  almost always redundant with the paragraph above it — in the sweep every
  single one was. They are also nearly impossible to write by accident, so the
  false-positive rate is the lowest of anything here.
* **Figures with no other use in this domain.** `posture`, `launder`,
  `minefield`, `papering over`, `in disguise`, `shop window`. These are
  ordinary English with ordinary meanings, and none of those meanings comes up
  in a symmetric-functions library. `posture` reached 12 sites before anyone
  noticed it was one word doing the job of "rule".

⚠️ **What this cannot check, which is most of the rule.** "One edge per tableau
is enumeration wearing a hash map" was found by reading, not by grep, and so
were the four bolded maxims in `llt.md`. This gate holds the vocabulary that
already decayed once; it says nothing about the next figure someone invents.
`--measure` exists so the claim "the false-positive rate is low" stays a
measurement rather than an assertion: it runs the patterns against an older
revision and reports what fraction of the hits that sweep actually changed.
"""

import pathlib
import re
import subprocess
import sys

# Pattern -> what to write instead. Each entry recurred in the 2026-08-08
# sweep; the comment is the count at the time it was found.
FIGURES = {
    r"\bpostures?\b": "rule, or the specific behavior",  # 12
    r"\blaunder(s|ed|ing)?\b": "hide, or state what is hidden",  # 4
    r"\bminefields?\b": "the conventions in circulation",  # 3
    r"\bpapering over\b": "hiding, or leaving unstated",  # 3
    r"\bin disguise\b": "name the thing it actually is",  # 3
    r"\bshop window\b": "what the README states",  # 3
    r"rather than a gamble\b": "rather than a guess",  # 3
    r"\bthe durable half\b": "name the half",  # 2
    r"\bthe payoff\b": "what it buys, stated",  # 2
    r"\bthe whole (point|reason|value) (of|is)\b": "that is why …",  # 22
    r"\bwearing (a|an|its|the)\b": "name what it is",  # 2
    r"\brhymes with\b": "failed the same way",  # 2
    r"\bshelf life\b": "state what goes stale and when",  # 1, kept: it recurs
    r"\bairtight\b": "sound",  # 1, kept: it recurs
}

# The maxim announcements. These introduce a general rule; the rule is what the
# ban is about, and the finding it generalizes is nearly always right above it.
ANNOUNCEMENTS = {
    r"\b[Tt]he (lesson|moral) (is|generalizes|here)\b": "state the finding",
    r"\bthe (portable|transferable|general) lesson\b": "state the finding",
    r"\bworth stating (as a rule|portably)\b": "state the finding",
    r"\b(lesson|decision|pattern|result) generalizes past\b": "state the case",
    r"\bthat is the finding\b": "delete; the finding follows",
    r"\bWorth stating as a rule\b": "state the finding",
}

PATTERNS = {**FIGURES, **ANNOUNCEMENTS}
COMPILED = [(re.compile(p), p, fix) for p, fix in PATTERNS.items()]

# A word inside a code span is being *named*, not used — this file's own
# vocabulary list and `docs/release-readiness.md`'s correction both quote the
# banned words that way. Same mechanism `scripts/check_spelling.py` uses, and
# the same reason.
CODE_SPAN = re.compile(r"``.*?``|`[^`\n]*`")
FENCE = re.compile(r"^\s*(```|~~~)")

ROOTS = ("src", "tests", "examples", "benches", "docs", "scripts")
DOC_FILES = ("README.md", "CLAUDE.md")
SUFFIXES = (".rs", ".md")

# `python` and `docsite` joined when the convenience layer did: a docstring on
# the surface a user actually reads is a prose surface `docs/style.md` governs,
# and the rendered site is what an outside reader meets first.
#
# `scripts/*.py` is deliberately *not* here, and that is a debt rather than a
# judgment: adding it surfaces four figures in `check_bindings.py`,
# `check_jack.py` and `check_llt.py` that predate this gate. They should be
# rewritten and the root added; doing it here would have buried the convenience
# layer's own review under unrelated edits.
PY_ROOTS = ("python",)
MD_ROOTS = ("docsite",)

# This file quotes every pattern it bans, and `docs/style.md` quotes the
# figures it retires as exhibits. Both are the rule, not violations of it.
EXEMPT = ("scripts/check_figures.py", "docs/style.md")


def sources(repo):
    for root in ROOTS:
        for suffix in SUFFIXES:
            yield from sorted((repo / root).rglob(f"*{suffix}"))
    for root in PY_ROOTS:
        yield from sorted((repo / root).rglob("*.py"))
    for root in MD_ROOTS:
        yield from sorted((repo / root).rglob("*.md"))
    for name in DOC_FILES:
        yield repo / name


def scan(text):
    """Yield (lineno, matched text, suggested fix)."""
    fenced = False
    for i, line in enumerate(text.splitlines(), 1):
        if FENCE.match(line):
            fenced = not fenced
            continue
        if fenced:
            continue
        line = CODE_SPAN.sub("", line)
        for rx, _, fix in COMPILED:
            m = rx.search(line)
            if m:
                yield i, m.group(0), fix


def measure(repo, rev):
    """Report what fraction of an older revision's hits the sweep changed.

    A hit is a true positive if the line it sits on is not in the tree today —
    that is, if someone came along and rewrote it. This undercounts (a figure
    could survive on a line that changed for another reason) and overcounts
    (a line could be deleted wholesale), so it is a rate and not a proof.
    """
    listing = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", rev],
        cwd=repo, capture_output=True, text=True, check=True,
    ).stdout.split()
    today = {}
    hits = changed = 0
    for path in listing:
        if not path.endswith(SUFFIXES) or path in EXEMPT:
            continue
        old = subprocess.run(
            ["git", "show", f"{rev}:{path}"],
            cwd=repo, capture_output=True, text=True, check=True,
        ).stdout
        if path not in today:
            here = repo / path
            today[path] = here.read_text(encoding="utf-8") if here.exists() else ""
        current = today[path]
        for lineno, found, _ in scan(old):
            hits += 1
            line = old.splitlines()[lineno - 1]
            if line not in current:
                changed += 1
            else:
                print(f"  survived: {path}:{lineno}: {found}")
    if not hits:
        print(f"no hits at {rev}")
        return 0
    print(f"\n{rev}: {hits} hits, {changed} on lines the sweep rewrote "
          f"({100 * changed // hits}%)")
    return 0


def main():
    repo = pathlib.Path(__file__).resolve().parent.parent
    args = sys.argv[1:]
    if args and args[0] == "--measure":
        return measure(repo, args[1] if len(args) > 1 else "HEAD~20")

    total = 0
    for path in sources(repo):
        rel = str(path.relative_to(repo))
        if rel in EXEMPT:
            continue
        for lineno, found, fix in scan(path.read_text(encoding="utf-8")):
            print(f"{rel}:{lineno}: {found!r} — {fix}")
            total += 1

    if not total:
        print("figures: no banned figure or maxim announcement in the tree")
        return 0
    print(f"\nfigures: {total} sites")
    print('see docs/style.md, "No aphorisms, no metaphors"')
    return 1


if __name__ == "__main__":
    sys.exit(main())
