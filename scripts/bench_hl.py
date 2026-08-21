"""A/B `Q'_λ` against Symmetrica's own `hall_littlewood`, in one process.

    sage -python scripts/bench_hl.py [max_degree]

Both sides are asked for the same thing -- Q'_λ expanded in the Schur basis with
polynomial coefficients -- and both are charged for building the Sage object, so
the comparison is end to end rather than kernel to kernel.

Symmetrica's C is reached through `sage.libs.symmetrica`, which needs `sage.all`
imported first; without it the binding raises inside `_py_polynom`. It is used
purely as a black box here: called, timed, and compared.

Correctness is checked before timing. That matters more than usual because this
*is* a port of Symmetrica's recursion, so agreement is a port check rather than
independent evidence -- the independent evidence is `charge::kostka_foulkes` in
the crate's own tests and Sage's `hall_littlewood().Qp()` in `check_hl.py`.
"""

import sys
import time

sys.path.insert(0, "pybuild")

from sage.all import QQ, PolynomialRing, Partitions, SymmetricFunctions  # noqa: E402
from sage.libs.symmetrica import all as symmetrica  # noqa: E402

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("the Symmetrica arm, which is reached through Sage")

import symfn  # noqa: E402

TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 13

R = PolynomialRing(QQ, "x")
x = R.gen()
s = SymmetricFunctions(R).schur()


def via_symfn(lam):
    return s._from_dict(
        {
            Partitions()(mu): R({e: c for e, c in poly})
            for mu, poly in symfn.hall_littlewood(list(lam))
        }
    )


def via_symmetrica(lam):
    return symmetrica.hall_littlewood(list(lam))


def clock(f, shapes):
    start = time.perf_counter()
    for lam in shapes:
        f(lam)
    return time.perf_counter() - start


print(f"{'n':>3} {'p(n)':>6} {'symfn':>10} {'symmetrica':>11} {'speedup':>8}")
for n in range(4, TOP + 1):
    shapes = [tuple(p) for p in Partitions(n)]

    for lam in shapes:
        a, b = via_symfn(lam), via_symmetrica(lam)
        if a != b:
            print(f"MISMATCH at {lam}:\n  symfn      {a}\n  symmetrica {b}")
            sys.exit(1)

    # Rotated, not interleaved: whichever pass runs second inherits the other's
    # cache footprint, and that alone was worth 36% in an earlier benchmark here.
    mine = min(clock(via_symfn, shapes), clock(via_symfn, shapes[::-1]))
    theirs = min(clock(via_symmetrica, shapes), clock(via_symmetrica, shapes[::-1]))
    print(f"{n:>3} {len(shapes):>6} {mine:>10.4f} {theirs:>11.4f} {theirs / mine:>7.2f}x")
