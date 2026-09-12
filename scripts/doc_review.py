#!/usr/bin/env python3
"""Split every rustdoc paragraph in `src/` into a reviewable worksheet.

One paragraph on the screen at a time, and enter means "nothing to say":

    python3 scripts/doc_review.py review --skip-tests
    python3 scripts/doc_review.py collect --skip-tests  # -> only what you flagged

`review` writes an answer to `doc-review.txt` before drawing the next
paragraph, so quitting or Ctrl-C never costs more than the one on screen, and
re-running resumes at the first *unanswered* paragraph. Saying nothing is an
answer and is recorded as one — a `~ reviewed` line — so a resumed pass does
not ask again about the ones you were happy with.

The same file can be filled in by hand instead, which is better for a long
sitting with a lot to say:

    python3 scripts/doc_review.py extract              # -> doc-review.txt
    $EDITOR doc-review.txt                             # type under `>` lines
    python3 scripts/doc_review.py collect              # -> only what you marked

Either way, one block per paragraph with a stable id, its `file:line`, and the
item the doc is attached to; `collect` prints back only the blocks you
commented on, which is the worklist.

Fenced code blocks (```text ... ```) stay whole rather than splitting on the
blank lines inside them, and a `# Panics` / `# Examples` heading rides with the
paragraph it introduces.

Scope it down when a full pass is too much at once:

    python3 scripts/doc_review.py extract --path src/character.rs src/kostka.rs
    python3 scripts/doc_review.py extract --module-only     # `//!` blocks only
    python3 scripts/doc_review.py extract --min-words 12    # skip one-liners
    python3 scripts/doc_review.py extract --skip-tests      # drop `#[cfg(test)]`

Re-running `extract` preserves comments already typed into the output file: a
block is matched by its `file:line:item` plus paragraph text, so edits to the
tree only orphan the paragraphs that actually changed (reported at the top).
"""

import argparse
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DEFAULT_OUT = REPO / "doc-review.txt"

# A block header, the paragraph, then the comment slot. `@@@` never begins a
# line of rustdoc, so the split is unambiguous even when a paragraph contains
# markdown headings, tables, or a fenced block.
SEP = "@@@"
SLOT = ">"
SEEN = "~"

HEADER = """\
# Rustdoc review worksheet — {count} paragraphs from {files} files
#
# HOW TO COMMENT. Each block is a paragraph followed by a bare `{slot}` line.
# Type after the `{slot}`. Leave the rest of the block alone. A bare `{slot}` means
# no comment, so skipping is the default and you only pay for what you flag:
#
#     {sep}
#     @ src/character.rs:20 [item] pub fn character(...) -> i128
#     Zero when |λ| ≠ |μ|, so a size mismatch is not an error.
#
#     {slot} shorter -- and say what happens at |λ| = |μ| = 0
#
# Multiple lines are fine, with or without repeating the `{slot}`; everything
# from the first `{slot}` to the end of the block is the comment.
#
# A `{seen} reviewed` line under the header means you looked at that paragraph and
# had nothing to say. `review` writes it for you; delete it to be asked again.
#
# Then:  python3 scripts/doc_review.py collect --skip-tests
#
# Shorthand that reads fine as a comment and needs no parsing on my end:
#   cut          this paragraph does not belong in the reference
#   -> record    move it to docs/record/, keep the shape here
#   shorter      right content, too many words
#   unclear      say what you did not follow
#   wrong        say what is actually true
"""


def rust_files(paths):
    """Every `.rs` file under the given roots, or under `src/` by default."""
    roots = [Path(p) for p in paths] if paths else [REPO / "src"]
    out = []
    for root in roots:
        if root.is_file():
            out.append(root)
        else:
            out.extend(sorted(root.rglob("*.rs")))
    return out


def doc_blocks(path, module_only, skip_tests):
    """Yield `(kind, start_line, item, [doc_lines])` for each run of doc lines.

    `item` is the signature the doc is attached to, trimmed to one line — it is
    what makes a block identifiable in the worksheet without opening the file.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    depth_of_test_mod = None
    depth = 0
    i = 0
    while i < len(lines):
        raw = lines[i]
        stripped = raw.strip()

        # Track `#[cfg(test)]` blocks by brace depth so `--skip-tests` can drop
        # them without a parser: test-module docs are a different genre and are
        # allowed history the reference is not.
        if skip_tests and re.match(r"#\[cfg\(test\)\]", stripped):
            depth_of_test_mod = depth
        depth += raw.count("{") - raw.count("}")
        if depth_of_test_mod is not None and depth <= depth_of_test_mod:
            depth_of_test_mod = None

        m = re.match(r"^\s*//(!|/)(?!/)(.*)$", raw)
        if not m:
            i += 1
            continue
        kind = "module" if m.group(1) == "!" else "item"
        if module_only and kind != "module":
            i += 1
            continue
        in_test = depth_of_test_mod is not None

        start = i + 1
        body = []
        while i < len(lines):
            m2 = re.match(r"^\s*//(!|/)(?!/)(.*)$", lines[i])
            if not m2 or ("module" if m2.group(1) == "!" else "item") != kind:
                break
            # Exactly one space, not `lstrip`: list indentation and fenced
            # alignment are the content in a good half of these blocks.
            text = m2.group(2)
            body.append(text[1:] if text.startswith(" ") else text)
            i += 1

        item = ""
        if kind == "item":
            j = i
            while j < len(lines) and (
                lines[j].strip().startswith("#[") or lines[j].strip().startswith("#!")
            ):
                j += 1
            if j < len(lines):
                item = lines[j].strip().rstrip("{").strip()
        yield kind, start, item, body, in_test


def paragraphs(body, start_line):
    """Split a doc block into paragraphs, keeping fenced code whole.

    A `#` heading is glued to what follows it: reviewing "# Panics" on its own
    tells you nothing, and the contract is in the sentence under it.
    """
    out = []
    cur, cur_line, fence = [], start_line, False
    for offset, line in enumerate(body):
        if line.startswith("```"):
            fence = not fence
        if not fence and not line.strip() and cur:
            out.append((cur_line, "\n".join(cur).rstrip()))
            cur, cur_line = [], start_line + offset + 1
            continue
        if not line.strip() and not cur:
            cur_line = start_line + offset + 1
            continue
        cur.append(line)
    if cur:
        out.append((cur_line, "\n".join(cur).rstrip()))

    merged = []
    for line_no, text in out:
        prev = merged[-1][1].strip() if merged else ""
        # A heading alone is not reviewable, and neither is a display block on
        # its own — both belong to the prose that introduces them.
        if merged and re.fullmatch(r"#{1,6} .*", prev):
            merged[-1][1] += "\n\n" + text
        elif merged and text.strip().startswith("```"):
            merged[-1][1] += "\n\n" + text
        else:
            merged.append([line_no, text])
    return [(ln, t) for ln, t in merged if not is_link_defs(t)]


def is_link_defs(text):
    """Markdown reference definitions carry no prose to comment on."""
    lines = [l for l in text.splitlines() if l.strip()]
    return bool(lines) and all(re.match(r"^\[[^\]]+\]:\s*\S+$", l.strip()) for l in lines)


def collect_paragraphs(args):
    items = []
    for path in rust_files(args.path):
        rel = path.relative_to(REPO)
        for kind, start, item, body, in_test in doc_blocks(
            path, args.module_only, args.skip_tests
        ):
            if args.skip_tests and in_test:
                continue
            for line_no, text in paragraphs(body, start):
                if len(text.split()) < args.min_words:
                    continue
                items.append(
                    {
                        "file": str(rel),
                        "line": line_no,
                        "kind": kind,
                        "item": item,
                        "text": text,
                    }
                )
    return items


def key_of(rec):
    """Identity used to carry comments across a re-extract.

    Deliberately *not* the line number: paragraphs move every time anything
    above them is edited, and a review file that loses its comments on the
    first unrelated edit is a review file nobody uses twice.
    """
    return (rec["file"], rec["item"], " ".join(rec["text"].split()))


def parse_worksheet(path):
    """Read an existing worksheet back into `{key: (comment, reviewed)}`.

    `reviewed` and `comment` are independent: a paragraph you looked at and had
    nothing to say about is *answered*, and a resumed
    pass must not show it again. Silence and not-yet-asked are different
    states, and only the file can tell them apart.
    """
    if not path.exists():
        return {}
    out = {}
    for block in path.read_text(encoding="utf-8").split("\n" + SEP + "\n"):
        m = re.search(r"^@ (\S+):(\d+) \[(\w+)\](?: (.*))?$", block, re.M)
        if not m:
            continue
        lines = block.splitlines()
        try:
            head = next(i for i, l in enumerate(lines) if l.startswith("@ "))
        except StopIteration:
            continue
        # Everything from the first `>` to the end of the block is the comment,
        # whether or not the continuation lines repeat the marker — typing
        # freely under the prompt is the obvious thing to do, and no paragraph
        # in the tree begins a line with `>`, so the split stays unambiguous.
        text, comment, in_slot, seen = [], [], False, False
        for line in lines[head + 1 :]:
            if line.startswith(SEEN) and not text:
                seen = True
            elif line.startswith(SLOT):
                in_slot = True
                comment.append(line[len(SLOT) :].strip())
            elif in_slot:
                comment.append(line.strip())
            else:
                text.append(line)
        body = " ".join("\n".join(text).split())
        note = "\n".join(c for c in comment if c).strip()
        if note or seen:
            out[(m.group(1), (m.group(4) or "").strip(), body)] = (note, True)
    return out


def render(items, existing):
    chunks = [
        HEADER.format(
            count=len(items),
            files=len({i["file"] for i in items}),
            slot=SLOT,
            sep=SEP,
            seen=SEEN,
        )
    ]
    orphans = dict(existing)
    for rec in items:
        note, seen = existing.get(key_of(rec), ("", False))
        orphans.pop(key_of(rec), None)
        head = f"@ {rec['file']}:{rec['line']} [{rec['kind']}]"
        if rec["item"]:
            head += " " + rec["item"]
        # A comment implies the paragraph was answered, so the marker is only
        # written where it carries information: looked at, nothing to say.
        mark = f"{SEEN} reviewed\n" if seen and not note else ""
        slot = "\n".join(f"{SLOT} {l}" for l in note.splitlines()) if note else SLOT
        chunks.append(f"{head}\n{mark}{rec['text']}\n\n{slot}\n")
    body = ("\n" + SEP + "\n").join(chunks)
    if orphans:
        lost = "\n".join(
            f"#   {f} [{i}] {(c or '(reviewed, no comment)').splitlines()[0][:60]}"
            for (f, i, _), (c, _s) in orphans.items()
        )
        body = (
            f"# ⚠️ {len(orphans)} comment(s) had no paragraph to attach to — the\n"
            f"# text they were written against changed. They are listed here and\n"
            f"# nowhere else:\n{lost}\n#\n" + body
        )
    return body


def cmd_extract(args):
    items = collect_paragraphs(args)
    out = Path(args.out)
    existing = parse_worksheet(out)
    out.write_text(render(items, existing), encoding="utf-8")
    kept = sum(1 for i in items if key_of(i) in existing)
    print(f"{len(items)} paragraphs from {len({i['file'] for i in items})} files")
    if existing:
        print(f"{kept}/{len(existing)} existing comments carried over")
    print(f"wrote {out}")


def cmd_collect(args):
    out = Path(args.out)
    if not out.exists():
        sys.exit(f"{out} does not exist — run `extract` first")
    items = collect_paragraphs(args)
    state = parse_worksheet(out)
    marked = {k: note for k, (note, _seen) in state.items() if note}
    if not marked:
        slots = sum(1 for l in out.read_text(encoding="utf-8").splitlines() if l.startswith(SLOT))
        answered = len(state)
        print(
            f"{out}: {slots} paragraphs, {answered} reviewed, none commented on"
            " — nothing to collect.\n"
            f"Type under a `{SLOT}` line, or run `review`, and re-run this."
        )
        return
    by_key = {key_of(i): i for i in items}
    dest = Path(args.into) if args.into else None
    lines = [f"# {len(marked)} commented paragraphs\n"]
    for key, note in marked.items():
        rec = by_key.get(key)
        where = f"{rec['file']}:{rec['line']}" if rec else f"{key[0]} (moved)"
        item = (rec["item"] if rec else key[1]) or ""
        text = rec["text"] if rec else key[2]
        lines.append(f"## {where}  {item}\n\n{text}\n\nCOMMENT: {note}\n")
    body = "\n".join(lines)
    if dest:
        dest.write_text(body, encoding="utf-8")
        print(f"{len(marked)} comments -> {dest}")
    else:
        print(body)


REVIEW_HELP = f"""\
  <enter>   nothing to say — marks it reviewed so it does not come back
  <text>    your comment (it is saved before the next paragraph is drawn)
  {SLOT}         start a multi-line comment; end it with a blank line
  b         back one paragraph, to change an answer
  s         skip without answering — it returns next session
  f         skip the rest of this file
  q         save and quit
  ?         this list

Shorthand that needs no parsing on my end: cut / -> record / shorter /
unclear / wrong. Say what you did not follow, or what is actually true."""


def clear():
    sys.stdout.write("\033[2J\033[3J\033[H")
    sys.stdout.flush()


def cmd_review(args):
    """One paragraph at a time, with the answers saved after every entry."""
    out = Path(args.out)
    items = collect_paragraphs(args)
    state = parse_worksheet(out)
    if not sys.stdin.isatty():
        sys.exit("`review` needs a terminal; use `extract` and edit the file instead")

    def save():
        out.write_text(render(items, state), encoding="utf-8")

    def answered(rec):
        return key_of(rec) in state

    total = len(items)
    done_at_start = sum(1 for r in items if answered(r))
    todo = [n for n, r in enumerate(items) if args.redo or not answered(r)]
    if not todo:
        print(f"all {total} paragraphs answered — nothing left in this scope")
        return

    print(
        f"{len(todo)} to review, {done_at_start} already answered, {total} in scope.\n"
        f"{REVIEW_HELP}\n"
    )
    input("enter to start ")

    pos = 0
    while 0 <= pos < len(todo):
        rec = items[todo[pos]]
        key = key_of(rec)
        note, _seen = state.get(key, ("", False))
        clear()
        left = len(todo) - pos
        print(f"\033[2m{len(todo) - left + 1}/{len(todo)}  ({left - 1} after this)\033[0m")
        print(f"\033[1m{rec['file']}:{rec['line']}\033[0m  [{rec['kind']}]")
        if rec["item"]:
            print(f"\033[2m{rec['item']}\033[0m")
        print()
        print(rec["text"])
        print()
        if note:
            print(f"\033[2mcurrent comment: {note}\033[0m")
        try:
            reply = input(f"{SLOT} ").strip()
        except (EOFError, KeyboardInterrupt):
            print()
            break

        if reply == "?":
            clear()
            print(REVIEW_HELP)
            input("\nenter to continue ")
            continue
        if reply == "q":
            break
        if reply == "b":
            pos = max(0, pos - 1)
            continue
        if reply == "s":
            pos += 1
            continue
        if reply == "f":
            here = rec["file"]
            while pos < len(todo) and items[todo[pos]]["file"] == here:
                pos += 1
            continue
        if reply == SLOT:
            print("(blank line ends the comment)")
            lines = []
            while True:
                try:
                    line = input("  ")
                except (EOFError, KeyboardInterrupt):
                    break
                if not line.strip():
                    break
                lines.append(line.rstrip())
            reply = "\n".join(lines)

        state[key] = (reply, True)
        save()
        pos += 1

    save()
    remaining = sum(1 for r in items if not answered(r))
    commented = sum(1 for _k, (n, _s) in state.items() if n)
    clear()
    print(f"saved {out}")
    print(f"{total - remaining}/{total} answered, {commented} with a comment")
    if remaining:
        print(f"{remaining} left — re-run `review` to pick up where you stopped")
    if commented:
        tail = " --skip-tests" if args.skip_tests else ""
        print(f"\nnext:  python3 scripts/doc_review.py collect{tail}")


def cmd_files(args):
    """Paragraph counts per file, so a pass can be scoped to one sitting."""
    counts = {}
    for rec in collect_paragraphs(args):
        counts[rec["file"]] = counts.get(rec["file"], 0) + 1
    for name, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        print(f"{n:5d}  {name}")
    print(f"{sum(counts.values()):5d}  total")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", choices=["review", "extract", "collect", "files"])
    ap.add_argument("--out", default=str(DEFAULT_OUT), help="worksheet path")
    ap.add_argument("--into", help="collect: write the worklist here instead of stdout")
    ap.add_argument("--path", nargs="*", default=[], help="files or dirs (default: src/)")
    ap.add_argument("--module-only", action="store_true", help="only `//!` blocks")
    ap.add_argument("--skip-tests", action="store_true", help="drop `#[cfg(test)]` docs")
    ap.add_argument("--min-words", type=int, default=1, help="skip shorter paragraphs")
    ap.add_argument("--redo", action="store_true", help="review: revisit answered paragraphs too")
    args = ap.parse_args()
    {
        "review": cmd_review,
        "extract": cmd_extract,
        "collect": cmd_collect,
        "files": cmd_files,
    }[args.mode](args)


if __name__ == "__main__":
    main()
