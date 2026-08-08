"""Verify every formula `docs/record/macdonald-operators.md` asserts,
numerically.

    sage -python scripts/verify_deltaop_formulas.py

This is not a `check_*.py`: there is no symfn code to dump yet. It is the step
that has to pass *before* a formula is written into the spec, per the rule the
`st` spec was built under -- a formula goes into a Sage comparison before it
goes into a document.

Sage is an ORACLE and nothing more. Its `Ht` basis and its `nabla` are called as
a black box and compared against formulas implemented here from the papers;
Sage's source is not read (see `../NOTICE.md`).

Sources, all read from the arXiv PDFs on 2026-07-29:

  [DM]  D'Adderio, Mellit, *A proof of the compositional Delta conjecture*,
        arXiv:2011.11467 -- (5)-(8) M, B_mu, T_mu, Pi_mu; (10) f*; (11) nabla;
        (12) Delta_f and Delta'_f; (21) Pi; (22) Theta_f.
  [DIV] D'Adderio, Iraci, Vanden Wyngaerd, *Theta operators, refined Delta
        conjectures, and coinvariants*, arXiv:1906.02623 -- (9) w_mu.
  [IR]  Iraci, Romero, *Delta and Theta operator expansions*, arXiv:2203.10342
        -- the star scalar product <F,G>_* = <F, (omega G)[MX]>.
  [HRW] Haglund, Remmel, Wilson, *The Delta Conjecture*, arXiv:1509.07058 --
        the rise and valley versions, Val(P), d_i(P), dinv, area.
"""

from sage.all import (SymmetricFunctions, QQ, PolynomialRing, FractionField,
                      Partitions, Partition, Permutations, StandardTableaux,
                      Subsets, prod)

R = FractionField(PolynomialRing(QQ, ['q', 't']))
q, t = R.gens()
Sym = SymmetricFunctions(R)
p, s, e, h, m = Sym.p(), Sym.s(), Sym.e(), Sym.h(), Sym.m()
Ht = Sym.macdonald().Ht()

M = (1 - q) * (1 - t)                                                   # [DM] (5)

FAILURES = []


def ok(label, cond):
    if not cond:
        FAILURES.append(label)
    print(f"  {'PASS' if cond else '*** FAIL ***'}  {label}")


# --------------------------------------------------------------- cell statistics
def cells(mu):
    """Every cell as (coleg, coarm) = (row, column), both 0-based."""
    return [(i, j) for i, r in enumerate(mu) for j in range(r)]


def B(mu):                                                              # [DM] (6)
    return sum(q**j * t**i for (i, j) in cells(mu))


def T(mu):                                                              # [DM] (7)
    return prod([q**j * t**i for (i, j) in cells(mu)], R.one())


def Pi(mu):                                                             # [DM] (8)
    return prod([1 - q**j * t**i for (i, j) in cells(mu) if (i, j) != (0, 0)],
                R.one())


def w(mu):                                                             # [DIV] (9)
    conj = Partition(mu).conjugate()
    out = R.one()
    for (i, j) in cells(mu):
        a, l = mu[i] - j - 1, conj[j] - i - 1
        out *= (q**a - t**(l + 1)) * (t**l - q**(a + 1))
    return out


# ------------------------------------------------------------ the star structure
def star_pair(F, G):
    """<F,G>_* diagonalised in the power sums:

        <p_rho, p_rho>_* = z_rho (-1)^{|rho|-l(rho)} prod_i (1-q^{rho_i})(1-t^{rho_i})
    """
    pf, pg = p(F), p(G)
    total = R.zero()
    for rho, c in pf:
        d = pg.coefficient(rho)
        if d == 0:
            continue
        total += (c * d * rho.centralizer_size() * (-1)**(sum(rho) - len(rho))
                  * prod([(1 - q**r) * (1 - t**r) for r in rho], R.one()))
    return total


def star_pair_via_omega(F, G):
    """[IR]'s spelling: <F, (omega G)[MX]>."""
    sub = p.zero()
    for rho, c in p(G.omega()):
        sub += c * prod([(1 - q**r) * (1 - t**r) for r in rho], R.one()) * p(rho)
    return p(F).scalar(sub)


def star_sub(F):
    """f* = f[X/M], [DM] (10): p_k -> p_k / ((1-q^k)(1-t^k)).

    A LINEAR map on the alphabet. Coefficients are not raised -- the same
    distinction `qtkostka.rs` documents for phi_t, and the same trap.
    """
    out = p.zero()
    for rho, c in p(F):
        out += c / prod([(1 - q**r) * (1 - t**r) for r in rho], R.one()) * p(rho)
    return out


def ht_coeffs(F, n):
    """F = sum_mu c_mu H~_mu, with c_mu = <F,H~_mu>_* / w_mu."""
    return {mu: star_pair(F, Ht[mu]) / w(mu) for mu in Partitions(n)}


def from_ht(coeffs):
    return sum((c * s(Ht[mu]) for mu, c in coeffs.items()), s.zero())


# ------------------------------------------------------------------- the operators
def nabla(F, n):
    return from_ht({mu: c * T(mu) for mu, c in ht_coeffs(F, n).items()})


def delta(f, F, n, prime=False):
    """Delta_f / Delta'_f, [DM] (12).

    f[B_mu] is f evaluated at the multiset of cell monomials; f[B_mu - 1] is the
    same with the cell (0,0) dropped, since that cell contributes exactly 1. No
    virtual alphabet is needed.
    """
    def ev(mu):
        alphabet = [q**j * t**i for (i, j) in cells(mu)]
        if prime:
            alphabet.remove(R.one())
        out = R.zero()
        for rho, c in p(f):
            out += c * prod([sum(x**r for x in alphabet) for r in rho], R.one())
        return out
    return from_ht({mu: c * ev(mu) for mu, c in ht_coeffs(F, n).items()})


def theta(f, F, n):
    """Theta_f F = Pi f* Pi^{-1} F, [DM] (21)-(22). F homogeneous of degree n."""
    if n == 0:
        return s(f) if f.degree() == 0 else s.zero()
    inner = from_ht({mu: c / Pi(mu) for mu, c in ht_coeffs(F, n).items()})
    prod_ = s(star_sub(f)) * inner
    return from_ht({mu: c * Pi(mu)
                    for mu, c in ht_coeffs(prod_, n + f.degree()).items()})


# ------------------------------------------------- labeled Dyck paths, [HRW]
def area_sequences(n):
    def rec(seq):
        if len(seq) == n:
            yield tuple(seq)
            return
        for a in range(1 if not seq else seq[-1] + 2):
            yield from rec(seq + [a])
    yield from rec([])


def labelings(a):
    """Labels strictly increasing up a rise."""
    n = len(a)
    for wd in Permutations(n):
        if all(a[i] != a[i - 1] + 1 or wd[i] > wd[i - 1] for i in range(1, n)):
            yield wd


def d_row(a, wd, i):
    """d_i(P): the dinv pairs beginning in row i."""
    return sum(1 for j in range(i + 1, len(a))
               if (a[i] == a[j] and wd[i] < wd[j])
               or (a[i] == a[j] + 1 and wd[i] > wd[j]))


def rise_side(n, k):
    want = n - k - 1
    total = R.zero()
    for a in area_sequences(n):
        rs = [i for i in range(1, n) if a[i] == a[i - 1] + 1]
        if len(rs) < want:
            continue
        for wd in labelings(a):
            base = q**sum(d_row(a, wd, i) for i in range(n)) * t**sum(a)
            for S in Subsets(rs, want):
                total += base * prod([t**(-a[i]) for i in S], R.one())
    return total


def valley_side(n, k):
    want = n - k - 1
    total = R.zero()
    for a in area_sequences(n):
        for wd in labelings(a):
            val = [i for i in range(1, n)
                   if a[i] < a[i - 1] or (a[i] == a[i - 1] and wd[i] > wd[i - 1])]
            if len(val) < want:
                continue
            dv = [d_row(a, wd, i) for i in range(n)]
            base = q**sum(dv) * t**sum(a)
            for S in Subsets(val, want):
                total += base * prod([q**(-(dv[i] + 1)) for i in S], R.one())
    return total


# =============================================================================
print("=== 1. the two spellings of <,>_* agree (Schur pairs, n <= 4) ===")
bad = [(mu, nu) for n in range(1, 5) for mu in Partitions(n) for nu in Partitions(n)
       if star_pair(s(mu), s(nu)) != star_pair_via_omega(s(mu), s(nu))]
ok("<F,G>_* = <F,(omega G)[MX]>", not bad)

print()
print("=== 2. H~ is *-orthogonal with <H~_mu,H~_mu>_* = w_mu ===")
for n in range(1, 6):
    diag = all(star_pair(Ht[mu], Ht[mu]) == w(mu) for mu in Partitions(n))
    off = all(star_pair(Ht[mu], Ht[nu]) == 0
              for mu in Partitions(n) for nu in Partitions(n) if mu != nu)
    ok(f"n={n} diagonal is w_mu", diag)
    ok(f"n={n} off-diagonal vanishes", off)

print()
print("=== 3. the H~-expansion round trips ===")
for n in range(1, 6):
    for F in [s(Partitions(n)[0]), e[n], h[n], m(Partitions(n)[-1])]:
        ok(f"n={n} round trip {F}", from_ht(ht_coeffs(F, n)) == s(F))

print()
print("=== 4. the closed forms for e_n ===")
for n in range(1, 7):
    ok(f"e_{n} = sum_mu M B_mu Pi_mu H~_mu / w_mu",
       s(e[n]) == sum((M * B(mu) * Pi(mu) / w(mu) * s(Ht[mu])
                       for mu in Partitions(n)), s.zero()))
    ok(f"e_{n}[X/M] = sum_mu H~_mu / w_mu",
       s(star_sub(e[n])) == sum((s(Ht[mu]) / w(mu) for mu in Partitions(n)),
                                s.zero()))

print()
print("=== 5. nabla has eigenvalue T_mu, and agrees with Sage's nabla ===")
for n in range(1, 6):
    ok(f"n={n} Sage's nabla H~_mu = T_mu H~_mu, every mu",
       all(s(Ht[mu].nabla()) == T(mu) * s(Ht[mu]) for mu in Partitions(n)))
    for F in [e[n], h[n], s(Partitions(n)[0]), m(Partitions(n)[-1])]:
        ok(f"n={n} our nabla = Sage's nabla on {F}", nabla(F, n) == s(F.nabla()))

print()
print("=== 6. dim of diagonal harmonics: <nabla e_n, h_1^n> = (n+1)^(n-1) ===")
for n in range(1, 7):
    got = s(e[n].nabla()).scalar(h[[1] * n])
    ok(f"n={n}: {got.subs(q=1, t=1)} = {(n + 1)**(n - 1)}",
       got.subs(q=1, t=1) == (n + 1)**(n - 1))

print()
print("=== 7. Delta_{e_n} = nabla on degree n ===")
for n in range(1, 6):
    for F in [e[n], h[n], m(Partitions(n)[-1])]:
        ok(f"n={n} Delta_e{n} = nabla on {F}", delta(e[n], F, n) == s(F.nabla()))

print()
print("=== 8. Delta'_{e_{n-1}} e_n = nabla e_n ===")
for n in range(1, 6):
    ok(f"n={n}", delta(e[n - 1], e[n], n, prime=True) == s(e[n].nabla()))

print()
print("=== 9. Theta_{e_k} nabla e_{n-k} = Delta'_{e_{n-k-1}} e_n ===")
for n in range(2, 6):
    for k in range(1, n):
        ok(f"n={n}, k={k}",
           theta(e[k], e[n - k].nabla(), n - k)
           == delta(e[n - k - 1], e[n], n, prime=True))

print()
print("=== 10. Delta conjecture, RISE version (a theorem) ===")
for n in range(1, 6):
    for k in range(0, n):
        ok(f"n={n}, k={k}  <Delta'_e{k} e_{n}, h_1^{n}> = Rise_{{n,k}}",
           delta(e[k], e[n], n, prime=True).scalar(h[[1] * n]) == rise_side(n, k))

print()
print("=== 11. Delta conjecture, VALLEY version (open in general) ===")
for n in range(1, 6):
    for k in range(0, n):
        ok(f"n={n}, k={k}  valley side",
           delta(e[k], e[n], n, prime=True).scalar(h[[1] * n]) == valley_side(n, k))

print()
print("=== 12. Delta'_{e_k} e_n has coefficients in N[q,t] ===")
for n in range(2, 6):
    for k in range(0, n):
        f = delta(e[k], e[n], n, prime=True)
        ok(f"n={n}, k={k}",
           all(c.denominator() == 1 and all(co > 0
                                            for co in c.numerator().coefficients())
               for _, c in f))

print()
print("=== 13. what Sage does and does not have ===")
for f in [e[2], s[2, 1], h[3]]:
    want = p.zero()
    for rho, c in p(f):
        want += c * prod([(1 - q**r) / (1 - t**r) for r in rho], R.one()) * p(rho)
    ok(f"theta_qt({f}) is the classical f[X(1-q)/(1-t)], not Theta_f",
       s(f.theta_qt(q, t)) == s(want))
for lam in [[2], [1, 1]]:
    got = p(lam).scalar_qt(p(lam))
    want = (Partition(lam).centralizer_size()
            * prod([(1 - q**r) / (1 - t**r) for r in lam], R.one()))
    ok(f"scalar_qt on p{lam} is Macdonald's <,>_{{q,t}}, not <,>_*", got == want)
names = [a for a in dir(s[2, 1]) if not a.startswith('_')]
ok("no Delta / Delta' / Theta_f / Pi anywhere on a symmetric function",
   not any(a in names for a in
           ('delta', 'Delta', 'delta_prime', 'Theta', 'delta_qt', 'big_pi')))

print()
if FAILURES:
    print(f"{len(FAILURES)} FAILURE(S):")
    for f in FAILURES:
        print("  -", f)
    raise SystemExit(1)
print("all formulas verified")
