"""Hold `kostka_foulkes` to Sage, and time it against Sage's own `kfpoly`.

    sage -python scripts/check_kf.py [max_degree]

Symmetrica has no Kostka-Foulkes entry point -- `hall_littlewood` is the only
way to reach these from it, and the transition has to be read off by hand -- so
unlike the Hall-Littlewood benchmark there is no C reference here. Sage's
`sage.combinat.sf.kfpoly` is the comparison, used strictly as a black box: run,
timed, and compared, never read (see NOTICE.md).

Every (lambda, mu) pair of each degree is checked, including the zeros, since a
transition that is right on its support and wrong about where the support *is*
would pass a nonzero-only comparison.
"""

import sys
import time

sys.path.insert(0, "pybuild")

from sage.all import QQ, PolynomialRing, Partitions  # noqa: E402
from sage.combinat.sf.kfpoly import KostkaFoulkesPolynomial  # noqa: E402

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's Kostka-Foulkes polynomials")

import symfn  # noqa: E402

TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 9

R = PolynomialRing(QQ, "t")
t = R.gen()


def mine(lam, mu):
    return sum(int(c) * t ** int(e) for e, c in symfn.kostka_foulkes(list(lam), list(mu)))


print(f"{'n':>3} {'pairs':>7} {'symfn':>10} {'sage':>10} {'speedup':>8}")
for n in range(1, TOP + 1):
    shapes = [tuple(p) for p in Partitions(n)]
    pairs = [(lam, mu) for mu in shapes for lam in shapes]

    bad = 0
    for lam, mu in pairs:
        got = mine(lam, mu)
        want = R(KostkaFoulkesPolynomial(list(lam), list(mu), t))
        if got != want:
            bad += 1
            print(f"  MISMATCH {lam} {mu}: symfn {got}  sage {want}")
    if bad:
        sys.exit(1)

    # symfn is asked per pair here, the same way Sage is, even though a whole
    # column comes out of one Q'_mu -- otherwise the two sides would not be
    # answering the same question. `kostka_foulkes_column` is the honest way to
    # ask for more than one, and is measured separately below.
    start = time.perf_counter()
    for lam, mu in pairs:
        mine(lam, mu)
    ours = time.perf_counter() - start

    start = time.perf_counter()
    for lam, mu in pairs:
        KostkaFoulkesPolynomial(list(lam), list(mu), t)
    theirs = time.perf_counter() - start

    print(f"{n:>3} {len(pairs):>7} {ours:>10.4f} {theirs:>10.4f} {theirs / ours:>7.2f}x")

# The same numbers, asked for the way the recursion wants to produce them.
n = TOP
shapes = [tuple(p) for p in Partitions(n)]
start = time.perf_counter()
for mu in shapes:
    symfn.kostka_foulkes_column(list(mu))
by_column = time.perf_counter() - start
print(f"\ndegree {n}: {len(shapes)**2} values by column in {by_column:.4f}s")
