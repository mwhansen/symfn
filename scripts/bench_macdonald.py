"""Time Sage's Macdonald P against symfn's, on the same work.

    sage -python scripts/bench_macdonald.py 8

Both sides expand `P_lambda(x; q, t)` in the monomial basis for **every**
partition of each degree, which is the unit of work either library is actually
asked for. Sage is timed as a black box: only the call, not the import, and the
`SymmetricFunctions` objects are built once outside the timed region so the
comparison is of the expansion and not of Sage's start-up.

symfn is timed from a separate `--example bench_mac` process, so its number
includes process start; that charge is against symfn, not for it.
"""

import subprocess
import sys
import time

from sage.all import QQ, PolynomialRing, SymmetricFunctions

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("the Sage arm")

top = int(sys.argv[1]) if len(sys.argv) > 1 else 7

R = PolynomialRing(QQ, "q,t").fraction_field()
Sym = SymmetricFunctions(R)
P = Sym.macdonald().P()
m = Sym.monomial()
Partitions = __import__("sage.all", fromlist=["Partitions"]).Partitions

# Warm the parent coercion machinery so the first degree is not charged for it.
m(P[1])

print(f"{'n':>3} {'sage':>10} {'symfn':>10} {'ratio':>8}")
for n in range(1, top + 1):
    t0 = time.perf_counter()
    for lam in Partitions(n):
        m(P[lam])
    sage_t = time.perf_counter() - t0

    out = subprocess.run(
        ["cargo", "run", "--release", "--quiet", "--example", "bench_mac", "--", str(n)],
        capture_output=True,
        text=True,
    ).stdout
    line = [l for l in out.splitlines() if f"mac P   n={n:2}" in l]
    mine = float(line[0].split()[4].rstrip("s"))

    ratio = f"{sage_t / mine:>7.1f}x" if mine > 0 else "      --"
    print(f"{n:>3} {sage_t:>9.4f}s {mine:>9.4f}s {ratio}")
