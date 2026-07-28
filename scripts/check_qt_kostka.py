"""Hold `qt_kostka` to Sage.

    sage -python scripts/check_qt_kostka.py [max_degree]

Timing lives in `bench_qt_kostka.py`, not here -- see the note in the loop.

Sage's `sage.combinat.sf.macdonald.qt_kostka` is the comparison, used strictly
as a black box: run, timed, and compared, never read (see NOTICE.md). Symmetrica
has no Macdonald polynomials at all, so there is no C reference for any of this
layer.

Every (lambda, mu) pair of each degree is checked. That is not the usual
belt-and-braces about zeros -- the (q,t)-Kostka matrix has **no** zero entries,
unlike the Kostka and Kostka-Foulkes tables it sits above -- but the two indices
enter asymmetrically and a transposed answer is otherwise a plausible-looking
matrix. The `q = 0` slice being Kostka-Foulkes is checked in Rust; what is
checked here is the whole polynomial, in both variables.
"""

import sys

sys.path.insert(0, "pybuild")

from sage.all import QQ, PolynomialRing, Partitions  # noqa: E402
from sage.combinat.sf.macdonald import qt_kostka  # noqa: E402

import symfn  # noqa: E402

TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 7

R = PolynomialRing(QQ, "q,t")
q, t = R.gens()


def mine(lam, mu):
    return R(sum(int(c) * q ** int(a) * t ** int(b)
                 for a, b, c in symfn.qt_kostka(list(lam), list(mu))))


print(f"{'n':>3} {'pairs':>7}")
for n in range(1, TOP + 1):
    shapes = [tuple(p) for p in Partitions(n)]
    pairs = [(lam, mu) for mu in shapes for lam in shapes]

    bad = 0
    for lam, mu in pairs:
        got = mine(lam, mu)
        want = R(qt_kostka(list(lam), list(mu)))
        if got != want:
            bad += 1
            print(f"  MISMATCH {lam} {mu}: symfn {got}  sage {want}")
    if bad:
        sys.exit(1)

    # The loop above already asked for every pair, and Sage caches the transition
    # matrices behind `qt_kostka` -- timing it again here would measure a
    # dictionary. That cache is also why per-pair timing is the wrong unit: it
    # charges symfn p(n) recomputations of J_mu against one of Sage's. The whole
    # table, cold on both sides, is the comparison both libraries can answer, and
    # it runs in `bench_qt_kostka.py` where a fresh process per degree is cheap
    # to arrange. Nothing is timed here.
    print(f"{n:>3} {len(pairs):>7}   ok")

# The table must be the columns, in the orientation `partitions(n)` fixes.
for n in range(1, TOP + 1):
    parts = symfn.partitions(n)
    table = symfn.qt_kostka_table(n)
    for j, mu in enumerate(parts):
        column = dict((tuple(lam), k) for lam, k in symfn.qt_kostka_column(mu))
        for i, lam in enumerate(parts):
            if table[i][j] != column[tuple(lam)]:
                print(f"  ORIENTATION n={n} [{i}][{j}] = K_{{{lam},{mu}}}")
                sys.exit(1)
print(f"qt_kostka_table: agrees with the columns through degree {TOP}")

