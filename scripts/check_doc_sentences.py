"""Hold item docs to the sentence discipline in `docs/style.md`.

    python3 scripts/check_doc_sentences.py            # the two gates
    python3 scripts/check_doc_sentences.py --report    # ranked long sentences

Two checks, both on `///` docs only. Module docs (`//!`) are the treatise and
are deliberately out of scope: the rules below serve a reader who landed on a
search hit and needs the contract, which is what an item doc is.

  * **The contract opener states its verb.** A `# Panics` or `# Errors`
    section that opens `If the value exceeds i128, around n ≈ 58.` has elided
    the verb, not the subject. `docs/style.md` permits the headless form the
    Rust ecosystem writes — `Panics if …`, `Returns every ν with …` — because
    the elided subject is the item the doc sits on and no reader resolves it
    wrongly. A missing *verb* leaves the reader to supply the relationship.
  * **A long sentence is reported, and the worst are a gate.** The cap is 25
    words for one sentence of an item doc. Sentences above `HARD` are the
    ones where a reader is asked to hold four facts at once; they fail the
    build so the tail does not grow back.

Needs no toolchain: it reads the sources. `scripts/preflight.sh` runs it.

**The cap is a ratchet, not a wall.** `HARD` is set just under the longest
sentence surviving the pass that introduced this script, so the check is green
on the tree it shipped with and red on anything longer. Lowering it is the
work; raising it is a decision, and belongs in a commit body that says why.

This is a lint, not a proof. Word count is a proxy for the thing that matters —
how many facts one sentence asks a reader to hold — and a 30-word list of
partitions is easier than a 20-word clause chain. `--report` exists because the
judgment is per sentence, and this script only ranks the candidates.
"""

import pathlib
import re
import sys

# The cap `docs/style.md` states, and the ratchet this script gates on.
CAP = 25
HARD = 50

# Sentence-final periods this tree does not have: a citation coordinate, an
# abbreviation, a section sign. Splitting on them reports phantom fragments.
ABBREV = re.compile(
    r"(?:e\.g|i\.e|cf|vs|ch|approx|Def|Thm|Prop|Cor|Lem|Eq|Ex|Fig|Ch|§\s*\d+)\.$"
)
# A sentence ends at `.`, `?` or `!` followed by space and a capital, a digit,
# an opening bracket, or a backtick — never mid-decimal, never inside `3.2`.
# Any non-ASCII glyph counts: this tree opens sentences with `ℓ(λ)`, `λ ⊢ n`
# and `Σ_ν`, and without them two sentences are measured as one.
SPLIT = re.compile(r"(?<=[.?!])\s+(?=[A-Z\[`(\d]|[^\x00-\x7f])")
CODE_SPAN = re.compile(r"`[^`]*`")
FENCE = re.compile(r"^\s*```")
# A table row is layout, not prose, and its pipes defeat sentence splitting.
TABLE = re.compile(r"^\s*\|")
# A list item is its own unit: run them together and one bullet's sentence
# swallows the next bullet's.
BULLET = re.compile(r"^\s*(?:[-*+]\s|\d+\.\s)")
# `# Panics` opens with a verb; a subordinator or a preposition means it did
# not. A heuristic, and the docstring's "lint, not a proof" covers the rest:
# it reads the first word only, so `Never, for any λ: the hooks are …` passes
# on a verb four clauses away.
VERBLESS = re.compile(
    r"^(?:If|When|Whenever|Where|Unless|On|For|In|At|As|With|Without"
    r"|Only|Past|After|Before|Above|Beyond)\b"
)
SECTION = re.compile(r"^#+\s+(Panics|Errors)\s*$")


def doc_blocks(path):
    """Yield (line number, list of doc lines) for each run of `///` in a file.

    A run is one item's documentation: consecutive `///` lines, broken by any
    line that is not one. Attributes between the doc and the item do not
    interrupt a run, because they come after it.
    """
    lines = path.read_text().splitlines()
    run, start = [], 0
    for i, raw in enumerate(lines, 1):
        stripped = raw.strip()
        if stripped.startswith("///"):
            if not run:
                start = i
            run.append(stripped[3:].strip())
        elif run:
            yield start, run
            run = []
    if run:
        yield start, run


def paragraphs(block):
    """Split a doc block into (line offset, section heading, text) paragraphs.

    Fenced blocks are dropped whole: a display formula is not prose, and its
    line breaks are load-bearing.
    """
    out, buf, offset, section, fenced = [], [], 0, None, False

    def flush():
        # `section` is cleared on the way out, so only the *first* paragraph
        # under `# Panics` is the opener the verb rule is about.
        nonlocal buf, section
        if buf:
            out.append((offset, section, " ".join(buf)))
            buf, section = [], None

    for i, line in enumerate(block):
        if FENCE.match(line):
            fenced = not fenced
            flush()
            continue
        if fenced or TABLE.match(line):
            continue
        heading = SECTION.match(line)
        if heading:
            flush()
            section = heading.group(1)
            continue
        if line.startswith("#"):
            flush()
            section = None
            continue
        if not line:
            flush()
            continue
        if BULLET.match(line):
            flush()
        if not buf:
            offset = i
        buf.append(BULLET.sub("", line))
    flush()
    return out


def sentences(text):
    """Split prose into sentences, keeping each code span as one word."""
    # Bold runs to a sentence end more often than not — `**Range.** The wall…`
    # — and the markers hide the boundary from the splitter.
    text = CODE_SPAN.sub("`x`", text).replace("**", "")
    parts, current = [], ""
    for piece in SPLIT.split(text):
        current = f"{current} {piece}".strip() if current else piece
        if ABBREV.search(current):
            continue
        parts.append(current)
        current = ""
    if current:
        parts.append(current)
    return parts


def words(sentence):
    return len(re.sub(r"[*_>]", "", sentence).split())


def scan(root):
    """Return (verbless openers, every sentence) across `src/*.rs`."""
    verbless, measured = [], []
    for path in sorted(root.glob("*.rs")):
        for start, block in doc_blocks(path):
            for offset, section, text in paragraphs(block):
                line = start + offset
                for n, sentence in enumerate(sentences(text)):
                    measured.append((words(sentence), path.name, line, sentence))
                    if n == 0 and section and VERBLESS.match(sentence):
                        verbless.append((path.name, line, section, sentence))
    return verbless, measured


def main():
    root = pathlib.Path(__file__).resolve().parent.parent / "src"
    verbless, measured = scan(root)
    report = "--report" in sys.argv

    if report:
        for count, name, line, sentence in sorted(measured, reverse=True):
            if count <= CAP:
                break
            print(f"{count:4d}  {name}:{line}  {sentence[:110]}")
        over = sum(1 for c, *_ in measured if c > CAP)
        print(
            f"\n{len(measured)} sentences, {over} over {CAP} words "
            f"({100 * over // len(measured)}%), longest {max(measured)[0]}"
        )
        return 0

    failures = 0
    for name, line, section, sentence in verbless:
        print(f"{name}:{line}: `# {section}` opens without a verb: {sentence[:70]}")
        failures += 1
    for count, name, line, sentence in sorted(measured, reverse=True):
        if count <= HARD:
            break
        print(f"{name}:{line}: {count} words, over the {HARD} ratchet: {sentence[:70]}")
        failures += 1

    if failures:
        print(f"\n{failures} item-doc sentences to fix (docs/style.md).")
        return 1
    over = sum(1 for c, *_ in measured if c > CAP)
    print(f"item docs: {over}/{len(measured)} sentences over {CAP} words, none over {HARD}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
