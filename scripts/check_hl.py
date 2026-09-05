"""Hold `hall_littlewood` to Sage: Q'_lambda expanded in the Schur basis.

    cargo run --release --example hldump -- 9 > /tmp/hl.txt
    sage -python scripts/check_hl.py /tmp/hl.txt

Sage is used as a black-box oracle only -- it is run as a separate program and
its output compared, never read for algorithms (see NOTICE.md).

The in-crate tests already check this against `charge::kostka_foulkes`, which is
independent of the recursion. Sage is a *third* opinion, and the one that would
catch a convention shared by both of ours -- t vs 1/t, Q' vs P, or a transposed
index -- since those are exactly the mistakes an internal cross-check cannot see.
"""

import sys
from collections import defaultdict

from sage.all import QQ, PolynomialRing, SymmetricFunctions

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's Hall-Littlewood Q'")

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/hl.txt"

R = PolynomialRing(QQ, "t")
t = R.gen()
Sym = SymmetricFunctions(R.fraction_field())
Qp = Sym.hall_littlewood().Qp()
s = Sym.schur()

mine = defaultdict(dict)
for line in open(path):
    lam, mu, body = (f.strip() for f in line.split("|"))
    lam = tuple(int(x) for x in lam.split(","))
    mu = tuple(int(x) for x in mu.split(","))
    poly = sum(int(c) * t**int(e) for e, c in (p.split(":") for p in body.split(",")))
    mine[lam][mu] = poly

bad = 0
for lam in sorted(mine, key=lambda p: (sum(p), p)):
    want = {}
    for mu, c in s(Qp[list(lam)]).monomial_coefficients().items():
        want[tuple(mu)] = R(c)
    got = mine[lam]
    if got != want:
        bad += 1
        keys = set(got) | set(want)
        for mu in sorted(keys):
            if got.get(mu, R(0)) != want.get(mu, R(0)):
                print(f"  {lam} -> {mu}: symfn {got.get(mu, 0)}  sage {want.get(mu, 0)}")

total = sum(len(v) for v in mine.values())
print(f"{len(mine)} partitions, {total} coefficients, {bad} disagreeing shapes")
sys.exit(1 if bad else 0)
