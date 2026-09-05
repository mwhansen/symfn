"""Verify symfn's general skewing against Sage's `skew_by`.

    sage -python scripts/check_skew.py [max_degree]

Correctness only. Skewing by g is the adjoint of multiplication by g, and the
point of the implementation is that the *basis g is written in* selects the
algorithm: h, e, and p take native Pieri / dual-Pieri / Murnaghan-Nakayama
paths that never touch Littlewood-Richardson, while s, m, and f go through it.

So each g is put to symfn **twice** — once in its own basis and once expanded
into Schur — and both are compared. Agreeing with Sage says the answer is
right; agreeing with each other says the fast path is a shortcut and not a
different operation.

The power-sum row is the one to watch. symfn computes it over the integers,
because rim-hook removal never divides; expanding p_mu into Schur needs 1/z_mu
and is rational. The boundary carries integers, so the via-Schur leg clears
denominators by the lcm first and the comparison is made against the scaled
result — which is also a check that the scaling is the only difference.

Sage is invoked only as a separate program whose output is used; see NOTICE.md.
"""

import sys
from functools import reduce
from math import lcm

sys.path.insert(0, "pybuild")

import symfn  # noqa: E402
from sage.all import QQ, Partitions, SymmetricFunctions  # noqa: E402

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's basis conversions")

Sym = SymmetricFunctions(QQ)
s = Sym.schur()
BASES = {
    "s": Sym.schur(),
    "h": Sym.homogeneous(),
    "e": Sym.elementary(),
    "p": Sym.power(),
    "m": Sym.monomial(),
    "f": Sym.forgotten(),
}

MAX_DEGREE = int(sys.argv[1]) if len(sys.argv) > 1 else 7


def norm(d):
    return sorted((tuple(k), QQ(v)) for k, v in d.items())


def norm_symfn(terms, scale=1):
    return sorted((tuple(k), QQ(c) / scale) for k, c in terms)


bad = 0
checked = 0

for df in range(0, MAX_DEGREE + 1):
    for lam in Partitions(df):
        lam = list(lam)
        f_terms = [(lam, 1)]
        for dg in range(0, df + 1):
            for mu in Partitions(dg):
                mu = list(mu)
                for name, basis in BASES.items():
                    g = basis[mu] if mu else basis.one()

                    want = norm(s[lam].skew_by(g).monomial_coefficients())
                    got = norm_symfn(symfn.skew_by(f_terms, [(mu, 1)], name))
                    checked += 1
                    if want != got:
                        bad += 1
                        print(f"MISMATCH {name}_{mu}^perp s_{lam}: sage {want} symfn {got}")

                    # The same g expanded into Schur, forced down the LR route.
                    # Clear denominators so it fits the integer boundary.
                    coeffs = s(g).monomial_coefficients()
                    den = reduce(lcm, (QQ(v).denominator() for v in coeffs.values()), 1)
                    scaled = [(list(k), int(QQ(v) * den)) for k, v in coeffs.items()]
                    via_s = norm_symfn(symfn.skew_by(f_terms, scaled, "s"), den)
                    checked += 1
                    if via_s != got:
                        bad += 1
                        print(f"ROUTE SPLIT {name}_{mu}^perp s_{lam}: native {got} via-s {via_s}")

print(f"\n{checked} checks through degree {MAX_DEGREE}, {bad} mismatches")
sys.exit(1 if bad else 0)
