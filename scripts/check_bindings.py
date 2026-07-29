"""Hold the *Python bindings* to Sage, for the non-classical families.

    maturin build --release --features python   # then unpack into pybuild/
    sage -python scripts/check_bindings.py 7

The other check scripts compare a dump written by a Rust example; this one calls
`symfn.macdonald_p` and friends the way Sage would, and rebuilds the answers as
Sage objects. It is the boundary that is under test — the marshalling of a
factored denominator, the exponent conventions, and the table orientations —
not the mathematics, which the dumps already cover.

Sage is a black-box oracle here, run as a separate program with its source not
read (see NOTICE.md).
"""

import sys

sys.path.insert(0, "pybuild")

import symfn  # noqa: E402
from sage.all import QQ, PolynomialRing, Partitions, SymmetricFunctions  # noqa: E402

top = int(sys.argv[1]) if len(sys.argv) > 1 else 6

QT = PolynomialRing(QQ, "q,t").fraction_field()
q, t = QT.gens()
MacSym = SymmetricFunctions(QT)
m = MacSym.monomial()
mac = {"P": MacSym.macdonald().P(), "Q": MacSym.macdonald().Q(), "J": MacSym.macdonald().J()}

T = PolynomialRing(QQ, "t")
(tt,) = T.gens()
HlSym = SymmetricFunctions(T.fraction_field())
hl_p = HlSym.hall_littlewood().P()
schur = HlSym.schur()

failures = 0


def fail(what, got, want):
    global failures
    failures += 1
    print(f"  {what}:\n    symfn {got}\n    sage  {want}")


def as_frac(num, den):
    """Rebuild a fraction-field element from the binding's two lists.

    This is the whole point of handing the denominator over *factored*: the
    caller multiplies the binomials it names, rather than being given an
    expanded polynomial it would have to factor again.
    """
    n = sum(c * q**a * t**b for a, b, c in num)
    d = QT(1)
    for a, b, mult in den:
        d *= (1 - q**a * t**b) ** mult
    return n / d


# --- Macdonald P, Q, J ------------------------------------------------------
for basis in ("P", "Q", "J"):
    count = 0
    call = {"P": symfn.macdonald_p, "Q": symfn.macdonald_q, "J": symfn.macdonald_j}[basis]
    for n in range(1, top + 1):
        for lam in Partitions(n):
            got = {tuple(mu): as_frac(num, den) for mu, num, den in call(list(lam))}
            want = {
                tuple(mu): QT(c)
                for mu, c in m(mac[basis][lam]).monomial_coefficients().items()
            }
            if got != want:
                fail(f"macdonald_{basis.lower()}({list(lam)})", got, want)
            count += len(got)
    print(f"macdonald_{basis.lower()}: {count} coefficients through degree {top}")

# J is the integral form: its coefficients are polynomials, so the binding must
# hand back an empty denominator list. Checking the value alone would not catch
# a denominator that failed to cancel.
for n in range(1, top + 1):
    for lam in Partitions(n):
        for mu, _num, den in symfn.macdonald_j(list(lam)):
            if den:
                fail(f"macdonald_j({list(lam)}) at {mu} has denominator", den, [])
print(f"macdonald_j: denominators empty through degree {top}")


# --- Hall-Littlewood P ------------------------------------------------------
def as_t_poly(terms):
    # Coerced, not left as a Python int: an empty term list is the zero
    # polynomial, and `sum` of nothing is `0`, which has no `.subs`.
    return T(sum(c * tt**e for e, c in terms))


count = 0
for n in range(1, top + 1):
    for lam in Partitions(n):
        got = {tuple(mu): as_t_poly(p) for mu, p in symfn.hall_littlewood_p(list(lam))}
        want = {
            tuple(mu): T(c)
            for mu, c in schur(hl_p[lam]).monomial_coefficients().items()
        }
        if {k: v for k, v in got.items() if v != 0} != want:
            fail(f"hall_littlewood_p({list(lam)})", got, want)
        count += len(got)
print(f"hall_littlewood_p: {count} coefficients through degree {top}")

# The table must agree with the one-shape call, shape for shape -- the two go
# through different code (`hall_littlewood_p` reads its answer out of the table,
# but the binding could still mismarshal one and not the other).
for n in range(1, top + 1):
    table = dict((tuple(lam), rows) for lam, rows in symfn.hall_littlewood_p_table(n))
    for lam in Partitions(n):
        solo = symfn.hall_littlewood_p(list(lam))
        if table[tuple(lam)] != solo:
            fail(f"hall_littlewood_p_table({n})[{list(lam)}]", table[tuple(lam)], solo)
print(f"hall_littlewood_p_table: agrees with the per-shape call through degree {top}")


# --- Kostka-Foulkes table ---------------------------------------------------
# Orientation is the thing under test: table[i][j] = K_{lambda^i lambda^j}(t),
# indexed as `partitions(n)` is, matching `kostka_table`. Transposing it would
# still produce a plausible-looking triangular matrix.
count = 0
for n in range(1, top + 1):
    parts = symfn.partitions(n)
    table = symfn.kostka_foulkes_table(n)
    kostka = symfn.kostka_table(n)
    for i, lam in enumerate(parts):
        for j, mu in enumerate(parts):
            got = as_t_poly(table[i][j])
            want = as_t_poly(symfn.kostka_foulkes(lam, mu))
            if got != want:
                fail(f"kostka_foulkes_table({n})[{i}][{j}] = K_{{{lam},{mu}}}", got, want)
            # t = 1 must be the ordinary Kostka number, in the same slot.
            if got(1) != kostka[i][j]:
                fail(
                    f"kostka_foulkes_table({n})[{i}][{j}] at t=1",
                    got(1),
                    kostka[i][j],
                )
            count += 1
# ...and the orientation check only has teeth if the matrix is not symmetric.
asym = symfn.kostka_foulkes_table(4)
if not any(asym[i][j] != asym[j][i] for i in range(len(asym)) for j in range(len(asym))):
    fail("kostka_foulkes_table(4) is symmetric, so the orientation check proves nothing",
         asym, "an asymmetric matrix")

print(f"kostka_foulkes_table: {count} entries through degree {top}, orientation and t=1 checked")

# --- the Macdonald operator algebra -----------------------------------------
#
# Only nabla has an external oracle; Delta', Theta and Pi exist nowhere else, so
# they are tied to it through the identities `deltaop` is tested on. What is
# under test here is the *boundary* -- exponent order in the (q_exp, t_exp, c)
# triples, and which slot a partition lands in -- not the mathematics.
def as_qt_schur(rows):
    return {tuple(lam): sum(c * q**a * t**b for a, b, c in poly) for lam, poly in rows}


s_qt = MacSym.schur()
e_qt = MacSym.elementary()

count = 0
for n in range(1, top + 1):
    got = as_qt_schur(symfn.nabla_e(n))
    want = {tuple(lam): QT(c)
            for lam, c in s_qt(e_qt[n].nabla()).monomial_coefficients().items()}
    if got != want:
        fail(f"nabla_e({n})", got, want)
    count += len(got)

    for lam in Partitions(n):
        arg = [(list(lam), [(0, 0, 1)])]
        got = as_qt_schur(symfn.nabla(arg))
        want = {tuple(mu): QT(c)
                for mu, c in s_qt(s_qt[lam].nabla()).monomial_coefficients().items()}
        if got != want:
            fail(f"nabla(s{list(lam)})", got, want)
        count += len(got)
print(f"nabla / nabla_e: {count} coefficients through degree {top}, vs Sage")

# Theta_{e_k} nabla e_{n-k} = Delta'_{e_{n-k-1}} e_n, across the boundary.
count = 0
for n in range(2, top + 1):
    for k in range(1, n):
        lhs = as_qt_schur(symfn.theta_ek(k, symfn.nabla_e(n - k)))
        rhs = as_qt_schur(symfn.delta_prime_e(n - k - 1, n))
        if lhs != rhs:
            fail(f"theta_ek({k}, nabla_e({n - k})) vs delta_prime_e({n - k - 1}, {n})",
                 lhs, rhs)
        count += 1
print(f"theta / delta_prime identity: {count} cases through degree {top}")

# Both combinatorial sides, across the boundary, at the degrees that are quick.
count = 0
for n in range(1, min(top, 6) + 1):
    want = [as_qt_schur(symfn.delta_prime_e(k, n)) for k in range(n)]
    for side in ("rise", "valley"):
        ladder = symfn.delta_conjecture_side(n, side)
        for k in range(n):
            # the binding returns the monomial basis; convert with Sage
            got_m = sum(c * m(list(mu)) for mu, poly in ladder[k]
                        for c in [sum(cc * q**a * t**b for a, b, cc in poly)])
            got = {tuple(lam): QT(c)
                   for lam, c in s_qt(got_m).monomial_coefficients().items()}
            if got != want[k]:
                fail(f"delta_conjecture_side({n},{side!r})[{k}]", got, want[k])
            count += 1
print(f"delta_conjecture_side: {count} (n,k,side) cases through degree "
      f"{min(top, 6)}, both versions")

print("FAILURES:", failures)
sys.exit(1 if failures else 0)
