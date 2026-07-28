"""Hold `macdonald_p` to Sage.

    cargo run --release --example macdump -- 7 > /tmp/mac.txt
    sage -python scripts/check_macdonald.py /tmp/mac.txt

Symmetrica has no Macdonald polynomials, so Sage is the *only* external oracle
here -- unlike Hall-Littlewood, where Symmetrica's C could be A/B'd as well.
It is used as a black box: run as a separate program, output compared, source
not read (see NOTICE.md).

Coefficients are compared as elements of the fraction field, so the factored
denominator symfn keeps and the expanded one Sage returns are the same object
regardless of how either side chose to write it. That matters: `(1+q)/(1-q^2)`
and `1/(1-q)` are equal and structurally different, and a string comparison
would call them different.
"""

import sys
from collections import defaultdict

from sage.all import QQ, PolynomialRing, SymmetricFunctions

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/mac.txt"
which = sys.argv[2] if len(sys.argv) > 2 else "P"

R = PolynomialRing(QQ, "q,t").fraction_field()
q, t = R.gens()
Sym = SymmetricFunctions(R)
P = {"P": Sym.macdonald().P, "Q": Sym.macdonald().Q, "J": Sym.macdonald().J}[which]()
m = Sym.monomial()


def parse(field):
    total = R(0)
    for piece in field.split(";"):
        exps, c = piece.split(":")
        a, b = (int(x) for x in exps.split(","))
        total += int(c) * q**a * t**b
    return total


mine = defaultdict(dict)
for line in open(path):
    lam, mu, num, den = (f.strip() for f in line.split("|"))
    lam = tuple(int(x) for x in lam.split(","))
    mu = tuple(int(x) for x in mu.split(","))
    mine[lam][mu] = parse(num) / parse(den)

bad = 0
for lam in sorted(mine, key=lambda p: (sum(p), p)):
    want = {tuple(mu): R(c) for mu, c in m(P[list(lam)]).monomial_coefficients().items()}
    got = mine[lam]
    if got != want:
        bad += 1
        for mu in sorted(set(got) | set(want)):
            a, b = got.get(mu, R(0)), want.get(mu, R(0))
            if a != b:
                print(f"  {lam} -> {mu}:\n    symfn {a}\n    sage  {b}")

total = sum(len(v) for v in mine.values())
print(f"{which}: {len(mine)} partitions, {total} coefficients, {bad} disagreeing shapes")
sys.exit(1 if bad else 0)
