"""Numerical verification of every formula in docs/record/jack.md,
against Sage.

Sage is an oracle here, never a source: nothing below was learned from
Sage's implementation, and each formula is stated from the papers cited in
the spec ([KS], [MOPS], Macdonald VI.10 via the q = t^alpha limit of the
Macdonald branching rule the crate already implements) and then CHECKED
against Sage's public API.

Run:  sage -python spec_jack_verify.py
Prints one PASS line per family and "all formulas verified" at the end.
"""

import sys
from itertools import product as iproduct

from sage.all import QQ, Partitions, Partition, SymmetricFunctions, prod, factorial

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's Jack bases")

R = QQ["t"]
ALPHA = R.gen()  # Sage calls the Jack parameter t; it is alpha throughout.
F = R.fraction_field()
Sym = SymmetricFunctions(F)
jack = Sym.jack()
P, J, Q, Qp = jack.P(), jack.J(), jack.Q(), jack.Qp()
m, s, p = Sym.monomial(), Sym.schur(), Sym.powersum()

MAXN_CHAIN = 6  # branching / LB recursion sweeps
MAXN_TAB = 5  # Knop-Sahi tableau enumeration (n^|lambda| labelings)


# ---------------------------------------------------------------- shapes

def parts(la):
    return list(la)


def conj(la):
    la = parts(la)
    if not la:
        return []
    return [sum(1 for x in la if x > j) for j in range(la[0])]


def cells(la):
    for i, row in enumerate(parts(la)):
        for j in range(row):
            yield (i, j)


def arm(la, i, j):
    return parts(la)[i] - j - 1


def leg(la, i, j):
    return conj(la)[j] - i - 1


def h_low(la, i, j):
    """Lower hook  alpha*a + l + 1  — the alpha-limit of (1 - q^a t^{l+1})."""
    return ALPHA * arm(la, i, j) + leg(la, i, j) + 1


def h_up(la, i, j):
    """Upper hook  alpha*(a+1) + l  — the alpha-limit of (1 - q^{a+1} t^l)."""
    return ALPHA * (arm(la, i, j) + 1) + leg(la, i, j)


def nfn(la):
    return sum(i * x for i, x in enumerate(parts(la)))


def dominates(la, mu):
    la, mu = parts(la), parts(mu)
    ca = cb = 0
    for i in range(max(len(la), len(mu))):
        ca += la[i] if i < len(la) else 0
        cb += mu[i] if i < len(mu) else 0
        if ca < cb:
            return False
    return True


# ------------------------------------------------- 1. branching formula

def horizontal_strips(base, size, cap):
    """All nu with base <= nu <= cap (containment), nu/base a horizontal
    strip of `size` cells.  Strip condition: nu_i <= base_{i-1} for i >= 1
    (at most one new cell per column); weak decrease then comes free."""
    base = parts(base)
    cap = parts(cap)
    rows = len(cap)
    b = base + [0] * (rows - len(base))
    out = []

    def rec(i, left, cur):
        if left < 0:
            return
        if i == rows:
            if left == 0:
                out.append([x for x in cur if x])
            return
        lo = b[i]
        hi = min(cap[i], lo + left)
        if i > 0:
            hi = min(hi, b[i - 1])
        for v in range(lo, hi + 1):
            rec(i + 1, left - (v - lo), cur + [v])

    rec(0, size, [])
    return out


def psi_alpha(lam, mu):
    """psi_{lam/mu}(alpha): product over cells of mu in rows that meet the
    strip lam/mu but not in columns that meet it, of
        [h_low_mu / h_up_mu] * [h_up_lam / h_low_lam]  at that cell
    — the per-atom q = t^alpha, t -> 1 limit of the Macdonald psi the crate
    already computes (macdonald.rs psi_factors)."""
    lam, mu = parts(lam), parts(mu)
    mup = mu + [0] * (len(lam) - len(mu))
    width = lam[0] if lam else 0
    colmeet = [sum(1 for x in lam if x > j) > sum(1 for x in mu if x > j) for j in range(width)]
    val = F(1)
    for i in range(len(mu)):
        if i < len(lam) and lam[i] == mup[i]:
            continue  # row does not meet the strip
        for j in range(mu[i]):
            if colmeet[j]:
                continue  # column meets the strip
            val *= h_low(mu, i, j) * h_up(lam, i, j)
            val /= h_up(mu, i, j) * h_low(lam, i, j)
    return val


def jack_p_via_branching(la):
    """P_la in the monomial basis, coefficient of m_mu = sum over chains of
    horizontal strips with sizes mu_1, mu_2, ... of the product of psi's."""
    la = parts(la)
    out = {}
    n = sum(la)
    for mu in Partitions(n):
        mu_l = list(mu)
        total = F(0)

        def rec(k, cur, acc):
            nonlocal total
            if k == len(mu_l):
                if cur == la:
                    total += acc
                return
            for nu in horizontal_strips(cur, mu_l[k], la):
                w = psi_alpha(nu, cur)
                rec(k + 1, nu, acc * w)

        rec(0, [], F(1))
        if total:
            out[tuple(mu_l)] = total
    return out


def check_branching():
    for n in range(1, MAXN_CHAIN + 1):
        for la in Partitions(n):
            mine = jack_p_via_branching(la)
            sage_ = m(P[la])
            got = {tuple(k): v for k, v in sage_.monomial_coefficients().items()}
            assert set(mine) == set(got), (la, set(mine) ^ set(got))
            for mu, v in mine.items():
                assert v == got[mu], (la, mu, v, got[mu])
    print(f"PASS  branching: P_la -> m via chains-of-strips psi^alpha == Sage, n <= {MAXN_CHAIN}")


# ---------------------------------------- 2. Laplace-Beltrami recursion

def jack_p_via_lb(kappa):
    """P_kappa in the monomial basis by the [MOPS] eigenoperator recursion:

        c_{kappa,la} = [ sum weights * c_{kappa,mu} ] / (E(kappa) - E(la)),
        E(nu) = alpha*n(nu') - n(nu),

    mu running over partitions obtained from la by moving t boxes from part
    j to an EARLIER part i (positions i < j in la, 1 <= t), sorted, with
    la < mu <= kappa in dominance; weight (la_i - la_j + 2t).
    """
    kappa = parts(kappa)
    n = sum(kappa)
    Ek = ALPHA * nfn(conj(kappa)) - nfn(kappa)
    coeffs = {tuple(kappa): F(1)}
    # dominance-descending order so every mu needed is already present
    order = sorted(Partitions(n), key=lambda x: nfn(x))
    for la in order:
        la_l = list(la)
        if tuple(la_l) == tuple(kappa):
            continue
        if not dominates(kappa, la_l):
            continue
        acc = F(0)
        L = len(la_l)
        for i in range(L):
            for j in range(i + 1, L):
                for t in range(1, la_l[j] + 1):
                    mu = la_l[:]
                    mu[i] += t
                    mu[j] -= t
                    mu = sorted((x for x in mu if x), reverse=True)
                    if not dominates(kappa, mu):
                        continue
                    c = coeffs.get(tuple(mu))
                    if c is None:
                        continue
                    acc += (la_l[i] - la_l[j] + 2 * t) * c
        den = Ek - (ALPHA * nfn(conj(la_l)) - nfn(la_l))
        if acc:
            coeffs[tuple(la_l)] = acc / den
    return coeffs


def check_lb():
    for n in range(1, MAXN_CHAIN + 1):
        for la in Partitions(n):
            mine = jack_p_via_lb(la)
            got = {tuple(k): v for k, v in m(P[la]).monomial_coefficients().items()}
            mine = {k: v for k, v in mine.items() if v}
            assert set(mine) == set(got), (la, set(mine) ^ set(got))
            for mu, v in mine.items():
                assert v == got[mu], (la, mu, v, got[mu])
    print(f"PASS  LB recursion: [MOPS] moving-box recursion == Sage, n <= {MAXN_CHAIN}")


# ------------------------------------------------ 3. Knop-Sahi tableaux

def jack_j_via_ks(la):
    """J_la as sum over admissible generalized tableaux [KS] Theorem 5.1.

    T labels cells with 1..n (n = |la|).  Admissible:
      (a) T(i,j) != T(i',j)   for i' > i          (below in the column)
      (b) T(i,j) != T(i',j-1) for i' < i, j > 1   (above in the column left)
    Critical cell: j > 1 and T(i,j) = T(i,j-1); weight prod of
      d(s) = alpha*(a(s)+1) + l(s) + 1  over critical s.
    """
    la = parts(la)
    k = sum(la)
    cs = list(cells(la))
    out = {}
    for labels in iproduct(range(1, k + 1), repeat=len(cs)):
        T = {c: v for c, v in zip(cs, labels)}
        ok = True
        for (i, j), v in T.items():
            for (i2, j2), v2 in T.items():
                if j2 == j and i2 > i and v2 == v:
                    ok = False
                    break
                if j > 0 and j2 == j - 1 and i2 < i and v2 == v:
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            continue
        w = F(1)
        for (i, j) in cs:
            if j > 0 and T[(i, j)] == T[(i, j - 1)]:
                w *= ALPHA * (arm(la, i, j) + 1) + leg(la, i, j) + 1
        content = [0] * k
        for c, v in T.items():
            content[v - 1] += 1
        out.setdefault(tuple(content), F(0))
        out[tuple(content)] = out[tuple(content)] + w
    # coefficient of m_mu = coefficient of x_1^{mu_1} x_2^{mu_2} ...
    res = {}
    for mu_key in {tuple(sorted((x for x in c if x), reverse=True)) for c in out}:
        exact = list(mu_key) + [0] * (k - len(mu_key))
        v = out.get(tuple(exact), F(0))
        if v:
            res[mu_key] = v
    return res


def check_ks():
    for n in range(1, MAXN_TAB + 1):
        for la in Partitions(n):
            mine = jack_j_via_ks(la)
            got = {tuple(k): v for k, v in m(J[la]).monomial_coefficients().items()}
            assert set(mine) == set(got), (la, set(mine) ^ set(got))
            for mu, v in mine.items():
                assert v == got[mu], (la, mu, v, got[mu])
    print(f"PASS  Knop-Sahi: J_la = sum over admissible tableaux of d_T(alpha) x^T, n <= {MAXN_TAB}")


# --------------------------------------------- 4. norms and normalizers

def check_norms_and_scalars():
    for n in range(1, 7):
        for la in Partitions(n):
            H = prod(h_low(la, i, j) for (i, j) in cells(la))
            Hp = prod(h_up(la, i, j) for (i, j) in cells(la))
            # J = H_la * P  (leading coefficient of J in m is H_la)
            assert m(J[la]).coefficient(la) == H, ("J leading", la)
            # <J,J> = H * H'
            assert J[la].scalar_jack(J[la]) == H * Hp, ("norm J", la)
            # <P,P> = H'/H
            assert P[la].scalar_jack(P[la]) == Hp / H, ("norm P", la)
            # Q = P / <P,P>, i.e. <P,Q> = 1
            assert P[la].scalar_jack(Q[la]) == 1, ("P-Q duality", la)
    # the alpha-deformed power-sum pairing.  ⚠️ Sage oracle wart, recorded:
    # scalar_jack crashes with TypeError when the scalar is the plain int 0
    # (jack.py normalize_coefficients calls c.denominator() and Python ints
    # have denominator as a property) — so the vanishing off-diagonal case
    # is exactly the case the oracle cannot state directly.
    for n in range(1, 7):
        for la in Partitions(n):
            for mu in Partitions(n):
                want = (
                    Partition(la).centralizer_size() * ALPHA ** len(la) if la == mu else 0
                )
                try:
                    got = P(p[la]).scalar_jack(P(p[mu]))
                except TypeError:
                    got = 0  # the int-zero crash path: Sage cannot RETURN a
                    # zero scalar_jack, so the crash set must coincide exactly
                    # with the pairs where want == 0 — asserted next line.
                assert got == want, ("p pairing", la, mu, got)
    print("PASS  norms: [m_la]J = H_la, <J,J> = H H', <P,P> = H'/H, <p,p> = z alpha^l, n <= 6")


def check_ks_positivity_and_u():
    """[KS] Theorem 1.1: coefficients of J in m, divided by u_mu = prod m_i(mu)!,
    are in N[alpha]."""
    for n in range(1, 8):
        for la in Partitions(n):
            for mu, c in m(J[la]).monomial_coefficients().items():
                u = prod(factorial(mult) for mult in Partition(mu).to_exp() if mult)
                q = R(c) / u
                assert q.denominator() == 1, ("u-divisibility", la, mu)
                assert all(x >= 0 and x in QQ and QQ(x).denominator() == 1 for x in R(q).coefficients()), (
                    "positivity",
                    la,
                    mu,
                )
    print("PASS  [KS] Thm 1.1: v_{la,mu}/u_mu in N[alpha] on Sage's J -> m, n <= 7")


# ------------------------------------------------- 5. specializations

def check_specializations():
    Sym1 = SymmetricFunctions(QQ)
    s1, m1 = Sym1.schur(), Sym1.monomial()
    for n in range(1, 7):
        for la in Partitions(n):
            # alpha = 1: J = hookprod * s_la
            hook1 = prod(arm(la, i, j) + leg(la, i, j) + 1 for (i, j) in cells(la))
            exp = m(J[la])
            got = {tuple(k): R(v).subs({ALPHA: 1}) for k, v in exp.monomial_coefficients().items()}
            want = m1(s1[la]) * hook1
            wantd = {tuple(k): v for k, v in want.monomial_coefficients().items()}
            assert got == wantd, ("alpha=1", la)
    # alpha = 2: which normalization is Sage's zonal Z?  Measure the ratio
    # Z_la / J_la^{(2)} on the leading m_la coefficient and report it.
    SymZ = SymmetricFunctions(QQ)
    Z = SymZ.zonal()
    mz = SymZ.monomial()
    ratios = []
    exact = True
    for n in range(1, 6):
        for la in Partitions(n):
            gotJ = {tuple(k): R(v).subs({ALPHA: 2}) for k, v in m(J[la]).monomial_coefficients().items()}
            wantZ = {tuple(k): v for k, v in mz(Z[la]).monomial_coefficients().items()}
            r0 = wantZ[tuple(la)] / gotJ[tuple(la)]
            ratios.append((list(la), r0))
            # whatever the leading ratio is, it must be a GLOBAL scalar
            for k, v in gotJ.items():
                if wantZ.get(k, 0) != r0 * v:
                    exact = False
    scal = "Z = J^(2) exactly" if all(r == 1 for _, r in ratios) else (
        f"Z_la = c_la * J_la^(2) with c_la a scalar per shape; c on n<=3: "
        + str([(la, r) for la, r in ratios if sum(la) <= 3])
    )
    assert exact, "Sage Z is not even a per-shape scalar multiple of J^(2)"
    print(f"PASS  alpha = 1: J_la = (prod hooks) s_la, n <= 6;  alpha = 2: {scal}")

    # duality: omega_alpha P_la^(alpha) = Q_{la'}^(1/alpha) in the p basis,
    # where omega_alpha := omega . (p_r -> alpha p_r), i.e. the p_mu
    # coefficient picks up alpha^{l(mu)} on top of omega's sign.  Plain
    # omega FAILS already at la = (1) — recorded as the trap.
    inv = R.fraction_field().hom([1 / ALPHA])
    for n in range(1, 6):
        for la in Partitions(n):
            left = p(P[la]).omega()
            left = {
                tuple(k): inv(ALPHA ** len(k) * v)
                for k, v in left.monomial_coefficients().items()
            }
            right = {tuple(k): v for k, v in p(Q[Partition(la).conjugate()]).monomial_coefficients().items()}
            assert left == right, ("duality", la, left, right)
    print("PASS  duality: omega_alpha P_la^(a) |_{a -> 1/a} = Q_{la'}^(1/a), n <= 5")

    # principal specialization: J_la(1^N) = prod over cells (N + alpha a'(s) - l'(s))
    for n in range(1, 6):
        for la in Partitions(n):
            for N in range(len(parts(la)), 5):
                want = prod(N + ALPHA * j - i for (i, j) in cells(la))
                poly = J[la].expand(N)
                got = poly(*([1] * N)) if N else 0
                assert got == want, ("principal", la, N)
    print("PASS  principal specialization: J_la(1^N) = prod (N + a'(s) alpha - l'(s)), n <= 5")


# ------------------------------------------------- 6. Stanley's g^nu

def check_stanley_g():
    """Stanley's 1989 conjecture is about the UNDIVIDED pairing:
    <J_la J_mu, J_nu>_alpha in N[alpha].  Observing it on small triples pins
    our normalization conventions; it is not evidence for the conjecture."""
    for na, nb in [(2, 2), (2, 3), (3, 3)]:
        for la in Partitions(na):
            for mu in Partitions(nb):
                prodJ = J[la] * J[mu]
                for nu in Partitions(na + nb):
                    g = prodJ.scalar_jack(J[nu])
                    gp = R(g)
                    assert all(c >= 0 and c.denominator() == 1 for c in gp.coefficients()), (
                        "g not in N[alpha]",
                        la,
                        mu,
                        nu,
                        g,
                    )
    print("PASS  Stanley: <J_la J_mu, J_nu> observed in N[alpha] on all triples, |la|,|mu| <= 3")


if __name__ == "__main__":
    # The convention gate: if Sage's parameter is not our alpha this dies
    # immediately and loudly.
    assert m(J[[2]]) == (ALPHA + 1) * m[[2]] + 2 * m[[1, 1]], "parameter convention"
    check_branching()
    check_lb()
    check_ks()
    check_norms_and_scalars()
    check_ks_positivity_and_u()
    check_specializations()
    check_stanley_g()
    print("all formulas verified")
