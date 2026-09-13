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
from sage.all import (QQ, Partition, Partitions, PolynomialRing,  # noqa: E402
                      SymmetricFunctions, catalan_number, prod)

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

    The denominator crosses *factored* so the caller multiplies the binomials
    it names, rather than receiving an expanded polynomial it would have to
    factor again.
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


# --- Jack ---------------------------------------------------------------------

# Sage calls the Jack parameter `t`; it is alpha throughout symfn, and this is
# the same ring Hall-Littlewood already uses above.  Naming it anything else
# breaks `scalar_jack`, which reaches for the ring's own generator.
AF = T.fraction_field()
alpha = AF(tt)
JackSym = SymmetricFunctions(AF)
jm = JackSym.monomial()
jp_basis = JackSym.powersum()
jack_bases = {
    "P": JackSym.jack().P(),
    "Q": JackSym.jack().Q(),
    "J": JackSym.jack().J(),
}


def as_jack(rows):
    """A JackTerms payload -> {partition: element of Q(alpha)}.

    This is the marshalling under test: a DENSE numerator indexed by the
    alpha-exponent, a FACTORED denominator of primitive atoms (u,v,mult)
    meaning (u*alpha+v)^mult, an integer scalar, and a dense TAIL polynomial
    that divides as well. The tail is the one denominator that is not a
    product of linear forms; it is empty except after a plethysm, which
    raises alpha to a power and so leaves the linear class. The row carried
    four fields until the 2026-08-25 widening added it
    (docs/record/python-and-sage-interop.md); the Sage adapter's _jack_cell
    is the same arithmetic on the same payload.
    """
    out = {}
    for mu, num, den, scale, tail in rows:
        v = sum(AF(c) * alpha**k for k, c in enumerate(num)) / AF(scale)
        for u, w, mult in den:
            v /= (u * alpha + w) ** mult
        if tail:
            v /= sum(AF(c) * alpha**k for k, c in enumerate(tail) if c)
        out[tuple(mu)] = v
    return out


def as_cell(cell):
    """A JackCell payload -> one element of Q(alpha).

    The same fields as a row of as_jack without the partition, and the same
    arithmetic; the two exist separately only because the entry points
    return one shape or the other.
    """
    num, den, scale, tail = cell
    v = sum(AF(c) * alpha**k for k, c in enumerate(num)) / AF(scale)
    for u, w, mult in den:
        v /= (u * alpha + w) ** mult
    if tail:
        v /= sum(AF(c) * alpha**k for k, c in enumerate(tail) if c)
    return v


count = 0
for name, call in (("P", symfn.jack_p), ("Q", symfn.jack_q), ("J", symfn.jack_j)):
    for n in range(1, top + 1):
        for la in Partitions(n):
            got = as_jack(call(list(la)))
            want = {tuple(k): v
                    for k, v in jm(jack_bases[name][la]).monomial_coefficients().items()}
            if got != want:
                fail(f"jack_{name.lower()}({list(la)})", got, want)
            count += 1
print(f"jack_p / jack_q / jack_j: {count} expansions through degree {top}, vs Sage")

# J is integral, so the denominator list, the scalar and the tail must all
# come back empty/1 -- the payload saying so is a separate claim from the
# value being right, and a wrong answer in the right slot is a different
# failure.
for n in range(1, top + 1):
    for la in Partitions(n):
        for mu, num, den, scale, tail in symfn.jack_j(list(la)):
            if den or scale != 1 or tail:
                fail(f"jack_j({list(la)}) at {mu}: J must be a polynomial",
                     (den, scale, tail), ([], 1, []))
print(f"jack_j: denominator, scalar and tail all empty through degree {top}")

count = 0
for n in range(1, top + 1):
    for la in Partitions(n):
        got = as_jack(symfn.jack_j_powersum(list(la)))
        want = {tuple(k): v
                for k, v in jp_basis(jack_bases["J"][la]).monomial_coefficients().items()}
        if got != want:
            fail(f"jack_j_powersum({list(la)})", got, want)
        count += 1
print(f"jack_j_powersum: {count} Jack character rows through degree {top}, vs Sage")

# The table must agree with the per-shape call -- the orientation check.
for n in range(1, top + 1):
    per_shape = {tuple(la): as_jack(symfn.jack_p(list(la))) for la in Partitions(n)}
    for la, rows in symfn.jack_table(n):
        if as_jack(rows) != per_shape[tuple(la)]:
            fail(f"jack_table({n}) at {la}", as_jack(rows), per_shape[tuple(la)])
print(f"jack_table: agrees with the per-shape call through degree {top}")

# The closed-form norm, handed over FACTORED, against Sage's scalar_jack.
count = 0
for n in range(1, top + 1):
    for la in Partitions(n):
        got = prod((u * alpha + w) ** mult for u, w, mult in symfn.jack_norm_j(list(la)))
        want = jack_bases["J"][la].scalar_jack(jack_bases["J"][la])
        if AF(got) != AF(want):
            fail(f"jack_norm_j({list(la)})", got, want)
        count += 1
print(f"jack_norm_j: {count} factored norms through degree {top}, vs scalar_jack")

# Stanley's structure constants.  OPEN conjecture: positivity is observed and
# never asserted -- only agreement with Sage is.
count = 0
for na, nb in ((2, 2), (2, 3)):
    for la in Partitions(na):
        for mu in Partitions(nb):
            prodJ = jack_bases["J"][la] * jack_bases["J"][mu]
            for nu in Partitions(na + nb):
                got = as_cell(
                    symfn.jack_structure_constant(list(la), list(mu), list(nu)))
                want = prodJ.scalar_jack(jack_bases["J"][nu])
                if got != want:
                    fail(f"jack_structure_constant({list(la)},{list(mu)},{list(nu)})",
                         got, want)
                count += 1
print(f"jack_structure_constant: {count} Stanley triples, vs Sage")

# jack_scalar takes INTEGRAL alpha-polynomial coefficients, so J is exactly the
# shape it accepts -- jack_j's rows go in unchanged, since the boundary refuses
# a row that is not the full five fields; <J_la, J_la> must be the closed-form
# norm.
count = 0
for n in range(1, min(top, 5) + 1):
    for la in Partitions(n):
        rows = symfn.jack_j(list(la))
        got = as_cell(symfn.jack_scalar(rows, rows))
        want = prod((u * alpha + w) ** mult for u, w, mult in symfn.jack_norm_j(list(la)))
        if got != AF(want):
            fail(f"jack_scalar(J_{list(la)}, J_{list(la)})", got, want)
        count += 1
print(f"jack_scalar: {count} norms via the general pairing, through degree {min(top, 5)}")

# The batch Stanley table must agree with the single-shot binding, including
# on the entries it omits -- a hoist that reused the wrong expansion would
# still produce a plausible, positive table.
#
# The two shapes differ by one field and that is the contract, not drift: a
# stanley_table row is (la, mu, nu, num, atoms, scale) with no tail, because
# Stanley's object is a polynomial in alpha, while jack_structure_constant
# returns a general JackCell and so carries the tail slot. The comparison is
# against the cell's first three fields, with the tail asserted empty rather
# than ignored -- dropping a field silently is how the comparison would stop
# checking anything.
for k in range(1, min(top, 3) + 1):
    batch = {(tuple(la), tuple(mu), tuple(nu)): (num, den, scale)
             for la, mu, nu, num, den, scale in symfn.stanley_table(k)}
    for la in Partitions(k):
        for mu in Partitions(k):
            for nu in Partitions(2 * k):
                key = (tuple(la), tuple(mu), tuple(nu))
                one = symfn.jack_structure_constant(list(la), list(mu), list(nu))
                if one[3]:
                    fail(f"jack_structure_constant{key}: a Stanley cell has no tail",
                         one[3], [])
                if key in batch:
                    if batch[key] != one[:3]:
                        fail(f"stanley_table({k}) at {key}", batch[key], one[:3])
                elif one[0]:
                    fail(f"stanley_table({k}) omits {key}", "missing", one[:3])
print(f"stanley_table: agrees with the per-triple call through k = {min(top, 3)}")

# Both zonal normalizations.  Sage's zonal() is P^(2), NOT J^(2); returning
# only one under an ambiguous name is how a caller gets plausible garbage.
Z = SymmetricFunctions(QQ).zonal()
mz = SymmetricFunctions(QQ).monomial()
count = 0
for n in range(1, min(top, 5) + 1):
    for la in Partitions(n):
        want = {tuple(k): v for k, v in mz(Z[la]).monomial_coefficients().items()}
        got_p = {tuple(mu): QQ(a) / QQ(b) for mu, a, b in symfn.zonal(list(la), False)}
        if got_p != want:
            fail(f"zonal({list(la)}, integral_form=False) vs Sage zonal()", got_p, want)
        # J = H_lambda * P, and P is monic at lambda, so the integral form's own
        # leading coefficient IS H_lambda(2) -- the whole difference between the
        # two normalizations, checked term by term rather than asserted.
        got_j = {tuple(mu): QQ(a) / QQ(b) for mu, a, b in symfn.zonal(list(la), True)}
        h2 = got_j[tuple(la)]
        if got_j != {k: v * h2 for k, v in got_p.items()}:
            fail(f"zonal({list(la)}): J^(2) is not H_lambda(2) * P^(2)", got_j, got_p)
        count += 1
print(f"zonal: {count} shapes, both normalizations, through degree {min(top, 5)}")

# The Goulden-Jackson tables.  b = 0 must be the class algebra of S_n, which
# the binding also exposes -- computed from characters alone, no Jack anywhere.
count = 0
for n in range(1, min(top, 6) + 1):
    c_tab, h_tab = symfn.gj_connection_tables(n)
    seen = {(tuple(la), tuple(mu), tuple(nu)): num
            for la, mu, nu, num, den in c_tab if den == 1}
    if len(seen) != len(c_tab):
        fail(f"gj_connection_tables({n}): c must be integral ([BD])", len(seen), len(c_tab))
    for la in Partitions(n):
        for mu in Partitions(n):
            for nu in Partitions(n):
                key = (tuple(la), tuple(mu), tuple(nu))
                num = seen.get(key, [])
                got = num[0] if num else 0
                want = symfn.class_algebra_coefficient(list(la), list(mu), list(nu))
                if got != want:
                    fail(f"c^{list(la)}_{{{list(mu)},{list(nu)}}}(0)", got, want)
                count += 1
    if not h_tab:
        fail(f"gj_connection_tables({n})", "empty h table", "nonempty")
print(f"gj_connection_tables: {count} c-coefficients at b=0 vs the S_n class algebra, "
      f"through degree {min(top, 6)}")


# --- LLT polynomials --------------------------------------------------------
# The mathematics is covered by check_llt.py (which compares a dump). What is
# under test here is the *boundary*: the q/t slot convention, the table
# orientations, the floor arriving as a separate call, and the graph edge
# conventions -- every one of which would type-check in Python while being wrong.
Q1 = PolynomialRing(QQ, "q")
(qq,) = Q1.gens()
LltSym = SymmetricFunctions(Q1.fraction_field())
llt_m = LltSym.monomial()
llt_s = LltSym.schur()
LLT = {k: LltSym.llt(k, t=qq) for k in (1, 2, 3)}


def as_q_poly(terms, what):
    """A binding's [(q_exp, t_exp, coeff)] -> a polynomial in q.

    The `t` slot must be zero: the LLT families proper live in q alone, and a
    stray t exponent here would otherwise be silently dropped.
    """
    out = Q1.zero()
    for a, b, c in terms:
        if b != 0:
            fail(f"{what} has a nonzero t exponent", (a, b, c), "t_exp == 0")
            return out
        out += c * qq**a
    return out


def as_m(rows, what):
    return {tuple(mu): as_q_poly(p, what) for mu, p in rows}


def sage_m(el):
    return {
        tuple(mu): Q1(c) for mu, c in llt_m(el).monomial_coefficients().items()
    }


count = 0
for k in (1, 2, 3):
    for n in range(1, min(top, 4) + 1):
        for mu in Partitions(n):
            lam = list(mu)
            got = as_m(symfn.llt_h(lam, k), f"llt_h({lam},{k})")
            want = sage_m(LLT[k].hspin()[mu])
            if {a: b for a, b in got.items() if b != 0} != want:
                fail(f"llt_h({lam}, {k})", got, want)
            got = as_m(symfn.llt_h_tilde(lam, k), f"llt_h_tilde({lam},{k})")
            want = sage_m(LLT[k].hcospin()[mu])
            if {a: b for a, b in got.items() if b != 0} != want:
                fail(f"llt_h_tilde({lam}, {k})", got, want)
            count += 2
print(f"llt_h / llt_h_tilde: {count} calls vs Sage hspin/hcospin, k = 1..3")

# cospin on a plain shape, and the Schur binding against it -- two different
# marshalling paths (monomial and Schur) onto one object.
count = 0
for k in (2, 3):
    for r in range(1, min(top, 4) + 1):
        for lam in Partitions(k * r):
            if list(Partition(list(lam)).core(k)):
                continue  # no k-ribbon tableaux
            l = list(lam)
            got = as_m(symfn.llt_gtilde(l, k), f"llt_gtilde({l},{k})")
            want = sage_m(LLT[k].cospin(Partition(l)))
            if {a: b for a, b in got.items() if b != 0} != want:
                fail(f"llt_gtilde({l}, {k})", got, want)
            got_s = {
                tuple(mu): as_q_poly(p, "llt_schur") for mu, p in symfn.llt_schur(l, k)
            }
            want_s = {
                tuple(nu): Q1(c)
                for nu, c in llt_s(LLT[k].cospin(Partition(l)))
                .monomial_coefficients()
                .items()
            }
            if {a: b for a, b in got_s.items() if b != 0} != want_s:
                fail(f"llt_schur({l}, {k})", got_s, want_s)
            count += 2
print(f"llt_gtilde / llt_schur: {count} calls vs Sage cospin, in m and in s")

# The tuple model, and the floor. Sage floors its tuple entry point and the
# binding deliberately does not, so the dictionary is only testable if BOTH
# calls cross correctly -- a floor of 0 would make a broken pair look fine, so
# the sweep has to contain a nonzero one, and this asserts that it does.
tuples = [[[1], [1]], [[2], [1]], [[1, 1], [1]], [[2], [2]], [[2, 1], [1]],
          [[1], [1], [1]], [[2], [1], [1]], [[1], [1, 1]], [[2, 2], [2, 1]]]
count, floored = 0, 0
for shapes in tuples:
    k = len(shapes)
    raw = as_m(symfn.llt_g(shapes), f"llt_g({shapes})")
    floor = symfn.llt_min_inv(shapes)
    floored += floor > 0
    got = {mu: Q1(p / qq**floor) for mu, p in raw.items()}
    want = sage_m(LLT[k].cospin(shapes))
    if {a: b for a, b in got.items() if b != 0} != want:
        fail(f"llt_g({shapes}) / q^{floor}", got, want)
    count += 1
if not floored:
    fail("llt_min_inv sweep", "every floor is 0", "at least one nonzero floor")
print(f"llt_g + llt_min_inv: {count} tuples vs Sage cospin ({floored} floored)")

# Table orientation: the whole-degree call must agree with the per-shape one.
for k in (2, 3):
    for n in range(1, min(top, 4) + 1):
        table = {tuple(mu): rows for mu, rows in symfn.llt_h_table(n, k)}
        if sorted(table) != sorted(tuple(mu) for mu in Partitions(n)):
            fail(f"llt_h_table({n}, {k}) index set", sorted(table),
                 sorted(tuple(mu) for mu in Partitions(n)))
        for mu in Partitions(n):
            solo = symfn.llt_h(list(mu), k)
            if table[tuple(mu)] != solo:
                fail(f"llt_h_table({n}, {k})[{list(mu)}]", table[tuple(mu)], solo)
print(f"llt_h_table: index set and per-shape agreement through degree "
      f"{min(top, 4)}, k = 2..3")

# k-core / k-quotient, including the component ORDER, which `G_nu` is not
# symmetric in and which nothing else at this boundary would catch.
count = 0
for k in (2, 3, 4):
    for n in range(1, min(top, 6) + 1):
        for lam in Partitions(n):
            core, quot = symfn.k_core_quotient(list(lam), k)
            p = Partition(list(lam))
            if tuple(core) != tuple(int(x) for x in p.core(k)):
                fail(f"k_core_quotient({list(lam)}, {k}) core", core, list(p.core(k)))
            want = [tuple(int(x) for x in c) for c in p.quotient(k)]
            if [tuple(c) for c in quot] != want:
                fail(f"k_core_quotient({list(lam)}, {k}) quotient", quot, want)
            count += 1
print(f"k_core_quotient: {count} shapes vs Sage core/quotient, k = 2..4")

# The t slot, which everything above asserted was zero. `nabla_e_by_path` is
# where it matters: t carries the area grading, and the pieces must sum
# to Sage's own nabla.
QT2 = PolynomialRing(QQ, "q,t").fraction_field()
q2, t2 = QT2.gens()
NabSym = SymmetricFunctions(QT2)
nab_m, nab_e = NabSym.monomial(), NabSym.elementary()
for n in range(1, min(top, 5) + 1):
    total = nab_m.zero()
    pieces = symfn.nabla_e_by_path(n)
    for _area, rows in pieces:
        for mu, terms in rows:
            total += nab_m(Partition(mu)) * sum(c * q2**a * t2**b for a, b, c in terms)
    want = nab_m(nab_e[n].nabla())
    if total != want:
        fail(f"nabla_e_by_path({n}) sums to nabla e_{n}", total, want)
    if len(pieces) != catalan_number(n):
        fail(f"nabla_e_by_path({n}) piece count", len(pieces), catalan_number(n))
print(f"nabla_e_by_path: sums to Sage nabla e_n and has C_n pieces, n <= "
      f"{min(top, 5)}")

# The graph bindings. Edge conventions are the whole risk: weak edges are
# ORDERED (an ascent is kappa(u) < kappa(v)) and strict edges constrain without
# scoring, so a swapped orientation or a double-counted strict edge is a
# plausible wrong answer rather than an error.
from sage.graphs.graph import Graph  # noqa: E402

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's parametric bases")

for n in range(2, min(top, 5) + 1):
    # the path P_n as a unit interval graph: edges (i, i+1), natural orientation
    weak = [(i, i + 1) for i in range(n - 1)]
    got = as_m(symfn.chromatic_from_llt(n, weak, []), "chromatic_from_llt")
    G = Graph([list(range(1, n + 1)), [(a + 1, b + 1) for a, b in weak]],
              format="vertices_and_edges")
    X = llt_m(G.chromatic_quasisymmetric_function(t=qq).to_symmetric_function())
    want = {tuple(mu): Q1(c) for mu, c in X.monomial_coefficients().items()}
    if {a: b for a, b in got.items() if b != 0} != want:
        fail(f"chromatic_from_llt(P_{n})", got, want)
    # and the LLT of the same graph must be the coloring generating function
    g = as_m(symfn.llt_graph(n, weak, []), "llt_graph")
    if not g:
        fail(f"llt_graph(P_{n})", "empty", "nonempty")
print(f"chromatic_from_llt: path graphs P_2..P_{min(top, 5)} vs Sage chromatic QSF")

# Isolated vertices must survive the boundary: a 2-vertex edgeless graph is not
# the empty graph, and Sage's own Graph([...]) drops them unless the vertex set
# is given explicitly -- the trap this fixture exists for.
iso = as_m(symfn.chromatic_from_llt(2, [], []), "chromatic_from_llt(edgeless)")
if iso.get((1, 1)) != 2 or iso.get((2,)) != 1:
    fail("chromatic_from_llt on 2 isolated vertices", iso, {(1, 1): 2, (2,): 1})
print("chromatic_from_llt: isolated vertices survive")

# The two error paths, which are the reason the graph binding validates at all.
for bad, why in (
    ((2, [], [(1, 0)]), "a strict edge oriented u > v"),
    ((2, [(0, 1)], [(0, 1)]), "an edge that is both weak and strict"),
    ((2, [(0, 5)], []), "an edge outside 0..n"),
):
    try:
        symfn.llt_graph(*bad)
        fail(f"llt_graph{bad} must raise", "accepted", why)
    except ValueError:
        pass
print("llt_graph: rejects bad edge sets instead of computing a wrong statistic")

# --- kronecker_coefficient -------------------------------------------------
#
# The binding is checked separately from the library because the failure it can
# have is its own: three partition arguments of the same shape, marshalled in
# the wrong order, give a *plausible* number rather than an error. g is symmetric
# in its three indices, so a permutation bug hides on symmetric inputs -- which
# is why the cases below are deliberately asymmetric.
SchurQQ = SymmetricFunctions(QQ).schur()
for lam, mu, nu in (
    ([3, 1], [2, 2], [2, 1, 1]),
    ([4, 2], [3, 2, 1], [3, 3]),
    ([5, 2, 1], [4, 3, 1], [3, 3, 2]),
    ([6, 3], [5, 4], [4, 4, 1]),
):
    want = (
        SchurQQ(Partition(lam))
        .itensor(SchurQQ(Partition(mu)))
        .coefficient(Partition(nu))
    )
    got = symfn.kronecker_coefficient(lam, mu, nu)
    if got != want:
        fail(f"kronecker_coefficient({lam},{mu},{nu})", got, want)
print("kronecker_coefficient: 4 asymmetric triples vs Sage itensor")

# Past where Sage's itensor is usable, the identities are the check: tensoring
# with the trivial character is the identity, so g^nu_{lam,(n)} = delta.
lam, n = [38, 2], 40
if symfn.kronecker_coefficient(lam, [n], lam) != 1:
    fail("kronecker_coefficient trivial-tensor at n=40", "!=1", 1)
if symfn.kronecker_coefficient(lam, [n], [37, 3]) != 0:
    fail("kronecker_coefficient trivial-tensor off-diagonal at n=40", "!=0", 0)
print("kronecker_coefficient: trivial-character identity holds at n = 40")

# The return type must be a Python int at both widths: Coeff::Big stringifies
# nothing on the way out.
if not isinstance(symfn.kronecker_coefficient([3, 1], [2, 2], [2, 1, 1]), int):
    fail("kronecker_coefficient return type", "not int", "int")
print("kronecker_coefficient: returns a Python int")

print("FAILURES:", failures)
sys.exit(1 if failures else 0)
