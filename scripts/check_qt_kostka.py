"""Hold `qt_kostka` to Sage.

    SAGE_DISABLE_SYMFN=1 sage -python scripts/check_qt_kostka.py [max_degree]

`SAGE_DISABLE_SYMFN` is not optional. Sage's Macdonald bases reach this library
through the backend whenever it is installed, and every comparison below would
then be symfn against symfn -- passing, and proving nothing. It has to be in the
environment before the process starts, because `Feature.is_present` caches. The
`symfn` module imported here is loaded directly and is unaffected.

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

from sage.all import QQ, PolynomialRing, Partitions, SymmetricFunctions  # noqa: E402
from sage.combinat.sf.macdonald import qt_kostka  # noqa: E402

import symfn  # noqa: E402

try:
    from sage.libs.symfn import is_available  # noqa: E402
except ImportError:
    pass  # stock Sage, with no backend to disable
else:
    if is_available():
        sys.exit("set SAGE_DISABLE_SYMFN=1: Sage's Macdonald bases would be symfn")

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

# The modified form, against Sage's own Ht basis rather than against a reflection
# of the check above -- otherwise the t^{n(mu)} would only ever be compared with
# itself. n(mu) is the one place a shape-dependent power could go wrong while
# still yielding polynomials.
QT = R.fraction_field()
schur = SymmetricFunctions(QT).schur()
ht = SymmetricFunctions(QT).macdonald().Ht()
count = 0
for n in range(1, TOP + 1):
    for mu in Partitions(n):
        got = {tuple(lam): QT(sum(int(c) * q ** int(a) * t ** int(b) for a, b, c in k))
               for lam, k in symfn.macdonald_ht(list(mu))}
        want = {tuple(lam): QT(c)
                for lam, c in schur(ht[mu]).monomial_coefficients().items()}
        if got != want:
            print(f"  MISMATCH Ht{list(mu)}: symfn {got}  sage {want}")
            sys.exit(1)
        count += len(got)
print(f"macdonald_ht: {count} coefficients against Sage's Ht through degree {TOP}")

# `schur_in_macdonald_j` is the inverse of the J -> s transition, and Sage
# reaches the same matrix by a triangular solve over QQ(q,t) -- so the two share
# the mathematics of J and nothing of how the inverse is obtained. Compared as
# elements of the fraction field, for the reason `check_macdonald.py` gives:
# the factored denominator symfn hands over and the expanded one Sage prints are
# the same object, and a string comparison would call them different.
J = SymmetricFunctions(QT).macdonald().J()
s = SymmetricFunctions(QT).schur()
count = 0
for n in range(0, TOP + 1):
    for lam, row in symfn.schur_in_macdonald_j(n):
        got = {tuple(mu): QT(sum(int(c) * q ** int(a) * t ** int(b) for a, b, c in num))
               / QT(R.prod((1 - q ** int(a) * t ** int(b)) ** int(e) for a, b, e in den))
               for mu, num, den in row}
        want = {tuple(mu): QT(c)
                for mu, c in J(s[list(lam)]).monomial_coefficients().items()}
        if got != want:
            print(f"  MISMATCH s{list(lam)} in J: symfn {got}  sage {want}")
            sys.exit(1)
        count += len(got)
print(f"schur_in_macdonald_j: {count} coefficients against Sage's J "
      f"through degree {TOP}")

