"""Verify symfn's evaluation and specializations against Sage.

    sage -python scripts/check_eval.py

Correctness only, no timings. These are new surface rather than a faster route
to something we already had, so the question is "does it agree", and Sage is
the oracle for all three:

* **evaluation** at an explicit alphabet, against `f.expand(n)` substituted;
* **s_λ(1^n)**, against `principal_specialization(n, q=1)`;
* **s_λ(1,q,…,q^{n-1})** as a polynomial in q, against
  `principal_specialization(n, q=q)`.

The alphabets deliberately include a **repeat and a zero**. That is the case the
bialternant s_λ = a_{λ+δ}/a_δ cannot express — a_δ vanishes when two variables
coincide, so the textbook formula is 0/0 there — and symfn uses the branching
rule precisely so that it has no such hole. An oracle check that only ever saw
distinct nonzero alphabets would not touch it.

Sage is invoked only as a separate program whose output is used; see NOTICE.md.
"""

import sys

sys.path.insert(0, "pybuild")

import symfn  # noqa: E402
from sage.all import QQ, PolynomialRing, SymmetricFunctions, Partitions  # noqa: E402

Sym = SymmetricFunctions(QQ)
s = Sym.schur()
R = PolynomialRing(QQ, "q")
q = R.gen()

MAX_DEGREE = int(sys.argv[1]) if len(sys.argv) > 1 else 7
ALPHABETS = [
    [2, -1, 3],
    [2, 2, 2, 0, 5],  # repeated and zero: the bialternant's blind spot
    [1, 1, 1, 1],
    [3, -2, 1, 4, -1],
]

bad = 0
checked = 0

for deg in range(1, MAX_DEGREE + 1):
    for lam in Partitions(deg):
        lam = list(lam)
        terms = [(lam, 1)]

        for xs in ALPHABETS:
            n = len(xs)
            want = QQ(s[lam].expand(n)(*xs))
            got = QQ(symfn.evaluate_schur(terms, xs))
            checked += 1
            if want != got:
                bad += 1
                print(f"MISMATCH eval s_{lam} at {xs}: sage {want} symfn {got}")

        for n in range(1, 6):
            want = QQ(s[lam].principal_specialization(n, q=QQ(1)))
            got = QQ(symfn.principal_specialization(lam, n))
            checked += 1
            if want != got:
                bad += 1
                print(f"MISMATCH s_{lam}(1^{n}): sage {want} symfn {got}")

            want_poly = R(s[lam].principal_specialization(n, q=q))
            got_poly = R(list(symfn.principal_specialization_q(lam, n)))
            checked += 1
            if want_poly != got_poly:
                bad += 1
                print(f"MISMATCH s_{lam}(1,q..q^{n-1}): sage {want_poly} symfn {got_poly}")

        # f^λ, against Sage's own standard-tableaux count.
        want = int(Partitions(deg)([*lam]).dimension())
        got = symfn.dimension(lam)
        checked += 1
        if want != got:
            bad += 1
            print(f"MISMATCH dim {lam}: sage {want} symfn {got}")

print(f"\n{checked} checks through degree {MAX_DEGREE}, {bad} mismatches")
sys.exit(1 if bad else 0)
