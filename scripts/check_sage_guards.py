"""Fail when a script that uses Sage could be comparing symfn against itself.

    python3 scripts/check_sage_guards.py

Sage reaches this library through its optional backend whenever symfn is
installed (README.md, "Under Sage"), so a script that treats Sage as an oracle,
or as a benchmark's control arm, has to refuse to run unless
`SAGE_DISABLE_SYMFN` is set. `scripts/sage_guard.py` states that rule and
`require_own_sage` enforces it; this file is what makes a *new* script obey it,
by failing when one imports Sage and neither calls the guard nor appears below.

It is a static scan: it imports nothing, needs no Sage, and runs inside
`scripts/preflight.sh` so it is checked on every commit.

The hole this closes was real and open. Four oracle scripts -- `check_jack.py`,
`check_macdonald.py`, `check_hl.py`, `check_hl_p.py` -- had no guard while the
backend was live in the development environment, so the families that already
dispatch were being held to themselves (`docs/record/python-and-sage-interop.md`).
"""

import ast
import pathlib
import re
import sys

HERE = pathlib.Path(__file__).resolve().parent

#: Scripts that use Sage and must NOT call the guard, each with the reason.
#: Every one of these manages `SAGE_DISABLE_SYMFN` itself, per arm, because
#: each measures what the backend does and needs it on for at least one arm.
MANAGES_ITS_OWN = {
    "check_backend.py": "A/Bs the backend and asserts which side answered",
    "bench_backend.py": "times both arms, setting the variable per arm",
    "bench_inverse.py": "sets it in the Sage arm's environment and refuses without it",
    "bench_macdonald_cache.py": "times the cache with the backend on and off",
    "bench_sf_candidates.py": "sets it in the LLT control arm's environment and refuses without it",
    "spec_llt_inverse.py": "needs the backend on in both arms and refuses without it",
    "check_qt_kostka.py": "carries the guard inline, predating sage_guard.py",
    "gen_sage_oracle.sage": "carries the guard inline; not a .py file",
}

SAGE_IMPORT = re.compile(r"\bsage\b")


def uses_sage(tree):
    """Whether the module imports anything from Sage, at any nesting depth."""
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            if any(SAGE_IMPORT.match(a.name) for a in node.names):
                return True
        elif isinstance(node, ast.ImportFrom):
            if node.module and SAGE_IMPORT.match(node.module):
                return True
    return False


def calls_guard(tree):
    """Whether the module calls `require_own_sage`, at any nesting depth."""
    return any(
        isinstance(node, ast.Call)
        and isinstance(node.func, ast.Name)
        and node.func.id == "require_own_sage"
        for node in ast.walk(tree)
    )


def main():
    problems = []
    guarded_count = 0
    allowed_count = 0
    for path in sorted(HERE.glob("*.py")):
        if path.name in ("sage_guard.py", pathlib.Path(__file__).name):
            continue
        tree = ast.parse(path.read_text(), filename=str(path))
        if not uses_sage(tree):
            continue
        allowed = path.name in MANAGES_ITS_OWN
        guarded = calls_guard(tree)
        if guarded:
            guarded_count += 1
        if allowed:
            allowed_count += 1
        if guarded and allowed:
            problems.append(
                f"{path.name}: calls require_own_sage but is listed in "
                "MANAGES_ITS_OWN; remove it from one or the other"
            )
        elif not guarded and not allowed:
            problems.append(
                f"{path.name}: imports Sage without calling require_own_sage. "
                "Add it (see scripts/sage_guard.py), or list the script in "
                "MANAGES_ITS_OWN here with the reason it measures the backend"
            )

    # A name that no longer exists would make the allowlist read as coverage it
    # does not have.
    for name in MANAGES_ITS_OWN:
        if not (HERE / name).exists():
            problems.append(f"MANAGES_ITS_OWN names {name}, which does not exist")

    if problems:
        print(f"FAIL: {len(problems)} problem(s)\n", file=sys.stderr)
        for p in problems:
            print(f"  {p}", file=sys.stderr)
        return 1
    print(
        f"ok: {guarded_count} Sage-using scripts refuse to run against a live "
        f"backend, {allowed_count} measure it by design"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
