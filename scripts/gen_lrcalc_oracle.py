"""Generate the lrcalc oracle fixture consumed by tests/lrcalc_oracle.rs.

    python3 scripts/gen_lrcalc_oracle.py > tests/fixtures/lrcalc_oracle.txt

lrcalc (Anders S. Buch, https://bitbucket.org/asbuch/lrcalc) is a second
independent oracle alongside Sage, and it earns its place by reaching shapes
Sage is too slow to enumerate: this fixture pins symfn's single-traversal skew
expansion on inputs an order of magnitude larger than the Sage fixture covers.

Only lrcalc's *output* is used here. symfn shares no code with lrcalc and is
not a derivative work of it; the binary is invoked as an external program, so
its GPL does not reach this crate.

Point LRCALC at the binary if it is not on PATH, e.g.

    LRCALC=/path/to/bin/lrcalc \\
    DYLD_LIBRARY_PATH=/path/to/lib \\
    python3 scripts/gen_lrcalc_oracle.py > tests/fixtures/lrcalc_oracle.txt

Line format (partitions comma-separated, empty partition = empty string):
    smul  MU|NU  LAM:COEFF ...
    skew  LAM|MU NU:COEFF ...
"""

import os
import re
import subprocess
import sys

LRCALC = os.environ.get("LRCALC", "lrcalc")

TERM = re.compile(r"^\s*(\d+)\s+\((.*)\)\s*$")


def run(*args):
    """Run lrcalc and parse its `COEFF  (p1,p2,...)` lines into a sorted list."""
    try:
        out = subprocess.run(
            [LRCALC, *[str(a) for a in args]],
            check=True,
            capture_output=True,
            text=True,
        ).stdout
    except FileNotFoundError:
        sys.exit(f"lrcalc not found at {LRCALC!r}; set $LRCALC to its path")
    except subprocess.CalledProcessError as e:
        sys.exit(f"lrcalc {' '.join(map(str, args))} failed:\n{e.stderr}")

    terms = []
    for line in out.splitlines():
        if not line.strip():
            continue
        m = TERM.match(line)
        if not m:
            sys.exit(f"unparsed lrcalc output line: {line!r}")
        coeff, parts = m.group(1), m.group(2)
        # "()" is the empty partition; lrcalc never emits a zero coefficient.
        key = [int(x) for x in parts.split(",")] if parts else []
        terms.append((key, int(coeff)))
    terms.sort()
    return terms


def enc(lam):
    return ",".join(str(x) for x in lam)


def emit(tag, a, b, terms):
    body = " ".join(f"{enc(k)}:{c}" for k, c in terms)
    print(f"{tag} {enc(a)}|{enc(b)} {body}")


# Products. The tail entries are deliberately far beyond what the Sage fixture
# covers -- that is where the single-traversal backend behaves differently from
# the coefficient-at-a-time one, so that is where cross-validation pays.
PRODUCTS = [
    ([1], [1]),
    ([2, 1], [2, 1]),
    ([3, 1], [2, 2]),
    ([3, 2, 1], [3, 2, 1]),
    ([4, 2], [3, 3, 1]),
    ([4, 3, 2, 1], [2, 2, 1]),
    ([5, 4, 3, 2, 1], [3, 2, 1]),
    ([5, 5, 5], [4, 4, 4]),
    ([6, 5, 4, 3, 2], [3, 3, 2]),
    ([2, 2, 2, 2, 2, 2], [3, 3, 3]),  # tall x wide: exercises the transpose path
    ([6, 5, 4, 3, 2, 1], [6, 5, 4, 3, 2, 1]),
    ([7, 6, 5, 4, 3], [4, 3, 2, 1]),
]

# Skew expansions, including shapes whose skew diagram is disconnected.
SKEWS = [
    ([3, 2, 1], [2, 1]),
    ([4, 3, 2], [2, 1]),
    ([5, 4, 3, 2], [3, 1]),
    ([4, 4, 4, 4], [2, 2]),
    ([6, 5, 4, 3, 2], [3, 2, 1]),
    ([7, 6, 5, 4, 3], [3, 2, 1]),
    ([8, 7, 6, 5, 4], [3, 2, 1]),
    ([2, 2, 2, 2, 2, 2, 2], [1, 1]),  # tall: exercises the transpose path
    ([5, 5, 5, 5, 5], [5, 5, 5, 5, 5]),  # empty skew: the unit
    ([6, 6, 4, 4, 2, 2], [4, 4, 2, 2]),
]

for mu, nu in PRODUCTS:
    emit("smul", mu, nu, run("mult", *mu, "-", *nu))

for lam, mu in SKEWS:
    emit("skew", lam, mu, run("skew", *lam, "/", *mu))
