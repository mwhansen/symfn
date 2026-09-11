"""Check the LLT `m -> H` inverse by back-substitution in the Schur basis.

    python scripts/spec_llt_inverse.py

Run with the sage-dev environment's `python`. A pre-implementation check: it
decides whether a symfn entry point for this inverse is worth building, by
doing the same algebra in the Sage adapter's language first.

The `sage` arm is the route on `mwhansen/sage` branch `symfn`: symfn fills
`H -> m` and Sage's `_invert_morphism` produces `m -> H` by a generic inverse
over `QQ(t)`, because the tables are not triangular in `m`. The `proto` arm
produces `m -> H` as `(m -> s) . (s -> H)`. `s -> H` inverts the Schur
expansion of `H`, which is lower triangular in `symfn.partitions(n)` order --
supported on shapes dominating mu -- with diagonal 1 for spin and a power of
`t` for cospin. The substitution runs in `ZZ[t, 1/t]` and divides only by that
monomial. The triangularity is asserted on every row, not assumed: a
counterexample stops the run.

Each arm runs in its own process, because Sage keeps the LLT caches at module
level. The two `m -> H` caches are compared as strings, exactly.

The columns are not the same quantity. Sage's is `_m_cache(n)`, both
directions. The prototype's is `m -> H` only, including the `H -> s` build
it needs; the forward table `llt_m_table` would be added to it in the adapter.

Needs the backend on in both arms, so it does not call the guard in
`sage_guard.py`; the `sage` arm refuses if symfn is not reached.

⚠️ Record the power state (`pmset -g batt`).
"""

import subprocess
import sys
import time

CASES = [
    ("spin", 2, 8), ("spin", 2, 10), ("spin", 2, 11), ("spin", 2, 12),
    ("spin", 3, 8), ("spin", 3, 9), ("spin", 3, 10),
    ("cospin", 2, 10), ("cospin", 2, 12), ("cospin", 3, 9), ("cospin", 3, 10),
]


def dump(cache):
    return "\n".join(
        f"{list(la)} "
        + repr(sorted((list(mu), str(c)) for mu, c in cache[la].items() if c))
        for la in sorted(cache))


def sage_arm(fam, k, n):
    import symfn
    from sage.all import QQ, SymmetricFunctions
    from sage.libs.symfn import is_available

    # Without the backend the fill is the ribbon-tableau enumeration, and the
    # timing would be of a different route.
    assert is_available(), "the sage arm measures the symfn-filled route"
    symfn.clear_caches()
    L = SymmetricFunctions(QQ["t"].fraction_field()).llt(k)
    B = L.hspin() if fam == "spin" else L.hcospin()
    start = time.perf_counter()
    B._m_cache(n)
    return time.perf_counter() - start, dump(B._m_to_self_cache[n])


def proto_arm(fam, k, n):
    import symfn
    from sage.all import QQ, ZZ, LaurentPolynomialRing
    from sage.combinat.partition import _Partitions

    symfn.clear_caches()
    QQt = QQ["t"].fraction_field()
    Lt = LaurentPolynomialRing(ZZ, "t")
    start = time.perf_counter()

    parts = symfn.partitions(n)
    index = {tuple(p): i for i, p in enumerate(parts)}
    size = len(parts)

    def lpoly(cs):
        return Lt({int(e): int(c) for e, _, c in cs})

    # m -> s, integral: kinv[nu] = {lambda index: coefficient}.
    kinv = [{index[tuple(la)]: c
             for la, c in symfn.monomial_to_schur([(list(nu), 1)])}
            for nu in parts]

    # H -> s, as rows M[mu] = {lambda index: Laurent polynomial}.
    if fam == "cospin":
        M = [{index[tuple(la)]: lpoly(cs)
              for la, cs in symfn.llt_schur([k * p for p in mu], k)}
             for mu in parts]
    else:
        M = [None] * size
        for mu, expansion in symfn.llt_h_table(n, k):
            row = {}
            for nu, cs in expansion:
                a = lpoly(cs)
                for la, c in kinv[index[tuple(nu)]].items():
                    row[la] = row.get(la, 0) + c * a
            M[index[tuple(mu)]] = {la: v for la, v in row.items() if v}
    fill = time.perf_counter() - start

    for i, row in enumerate(M):
        assert all(j <= i for j in row), (fam, k, n, parts[i], "not triangular")
        d = row[i]
        assert d.number_of_terms() == 1 and d.coefficients()[0] == 1, \
            (fam, k, n, parts[i], "diagonal is not a monomial", d)

    # N = M^-1, lower triangular:
    # N[i][j] = -(1/M[i][i]) * sum over j <= mid < i of M[i][mid] N[mid][j].
    N = [None] * size
    for i in range(size):
        inv_d = ~M[i][i]
        acc = {}
        for mid, c in M[i].items():
            if mid != i:
                for j, v in N[mid].items():
                    acc[j] = acc.get(j, 0) + c * v
        row = {i: inv_d}
        row.update({j: -inv_d * v for j, v in acc.items() if v})
        N[i] = row

    # m -> H = (m -> s) . (s -> H).
    sage_parts = [_Partitions.from_parts(p) for p in parts]
    out = {}
    for a in range(size):
        acc = {}
        for la, c in kinv[a].items():
            for mu, v in N[la].items():
                acc[mu] = acc.get(mu, 0) + c * v
        out[sage_parts[a]] = {sage_parts[mu]: QQt(v)
                              for mu, v in acc.items() if v}
    return time.perf_counter() - start, fill, dump(out)


def main():
    if sys.argv[1:2] == ["--arm"]:
        arm, fam, k, n = sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5])
        if arm == "sage":
            total, d = sage_arm(fam, k, n)
            print(f"{total}\n{d}")
        else:
            total, fill, d = proto_arm(fam, k, n)
            print(f"{total} {fill}\n{d}")
        return
    print(f"{'family':<7}{'k':>2}{'n':>4}  {'Sage _m_cache':>14}  "
          f"{'proto m->H':>11} {'(H->s part)':>12}  identical")
    for fam, k, n in CASES:
        res = {}
        for arm in ("sage", "proto"):
            r = subprocess.run(
                [sys.executable, __file__, "--arm", arm, fam, str(k), str(n)],
                capture_output=True, text=True, check=True)
            head, _, body = r.stdout.partition("\n")
            res[arm] = ([float(x) for x in head.split()], body)
        (st,), sd = res["sage"]
        (pt, pf), pd = res["proto"]
        print(f"{fam:<7}{k:>2}{n:>4}  {st:>13.3f}s  {pt:>10.3f}s "
              f"{pf:>11.3f}s  {sd == pd}")


if __name__ == "__main__":
    main()
