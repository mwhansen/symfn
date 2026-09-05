"""Hold `hall_littlewood_p` to Sage: P_lambda(x;t) in the Schur basis.

    cargo run --release --example hldump -- 12 P > /tmp/hlp.txt
    sage -python scripts/check_hl_p.py /tmp/hlp.txt

P is obtained here by inverting the Kostka-Foulkes matrix rather than by its own
enumeration, so this checks the inversion as much as the recursion feeding it.
Sage is a black-box oracle only (see NOTICE.md).
"""

import sys
from collections import defaultdict

from sage.all import QQ, PolynomialRing, SymmetricFunctions

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's Hall-Littlewood P")

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/hlp.txt"

R = PolynomialRing(QQ, "t")
t = R.gen()
Sym = SymmetricFunctions(R.fraction_field())
P = Sym.hall_littlewood().P()
s = Sym.schur()

mine = defaultdict(dict)
for line in open(path):
    lam, mu, body = (f.strip() for f in line.split("|"))
    lam = tuple(int(x) for x in lam.split(","))
    mu = tuple(int(x) for x in mu.split(","))
    mine[lam][mu] = sum(
        int(c) * t ** int(e) for e, c in (p.split(":") for p in body.split(","))
    )

bad = 0
for lam in sorted(mine, key=lambda p: (sum(p), p)):
    want = {tuple(mu): R(c) for mu, c in s(P[list(lam)]).monomial_coefficients().items()}
    got = {k: v for k, v in mine[lam].items() if v != 0}
    if got != want:
        bad += 1
        for mu in sorted(set(got) | set(want)):
            a, b = got.get(mu, R(0)), want.get(mu, R(0))
            if a != b:
                print(f"  {lam} -> {mu}: symfn {a}  sage {b}")

total = sum(len(v) for v in mine.values())
print(f"{len(mine)} partitions, {total} coefficients, {bad} disagreeing shapes")
sys.exit(1 if bad else 0)
