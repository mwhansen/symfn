"""Numerical verification of every formula in docs/record/llt.md,
against Sage.

Sage is an oracle here, never a source: nothing below was learned from
Sage's implementation. Each model is stated from the papers cited in the
spec ([HHL] math/0409538, [LLT] q-alg/9512031, [LT] math/9809122,
[BHMPS] 2112.07063, [DA] 1906.02633, [AS] 2004.09198, [KMS] q-alg/9508006)
and then CHECKED against Sage's public API and against the papers' own
worked examples.

Run:  sage -python spec_llt_verify.py
Prints one PASS line per family and "all formulas verified" at the end.
"""

import itertools
import sys
from fractions import Fraction

from sage.all import (QQ, FractionField, Partition, Partitions,
                      SymmetricFunctions, prod)

# ---------------------------------------------------------------- rings

Rq = QQ['q']
q = Rq.gen()
Rqt = QQ['q,t']
qq, tt = Rqt.gens()

Sym = SymmetricFunctions(FractionField(Rqt))
s_b = Sym.schur()
m_b = Sym.monomial()
p_b = Sym.powersum()
e_b = Sym.elementary()
h_b = Sym.homogeneous()


def parts_tuple(mu):
    return tuple(int(x) for x in mu)


# ================================================================ models
# --- the tuple (attacking-inversion) model, [BHMPS] Def 2.2.1 == [HHL] §3.
#
# A component is a list of raw cells (row, col), 0-based, forming a skew
# shape (or any subset closed under the adjacency we test).  Content of a
# cell is col - row + offset(component)  (classical Macdonald convention;
# [HHL]'s printed c(u) = i - j is row - col, the same numbers negated —
# their (ii) condition maps onto ours, checked by the gate fixture).
# Attacking pairs (a, b), a before b in reading order (content ascending,
# component index breaking ties):
#   same content, comp(a) < comp(b);   content(b) = content(a)+1, comp(b) < comp(a).
# inv(T) = #{attacking pairs with T(a) > T(b)}.   G = sum q^inv x^T.


def skew_cells(outer, inner=()):
    """Cells (row, col), 0-based, of outer/inner."""
    inner = tuple(inner) + (0,) * (len(outer) - len(inner))
    return [(r, c) for r, row in enumerate(outer) for c in range(inner[r], row)]


class TupleShape:
    """A tuple of components; each component = (cells, content_offset)."""

    def __init__(self, comps):
        # comps: list of (cells, offset)
        self.comps = comps
        self.cells = []  # (comp, row, col, content)
        for k, (cells, off) in enumerate(comps):
            for (r, c) in cells:
                self.cells.append((k, r, c, c - r + off))
        self.n = len(self.cells)
        # reading order: content ascending, then component ascending, then
        # position along the diagonal (row ascending is fine within a comp).
        self.order = sorted(range(self.n),
                            key=lambda i: (self.cells[i][3], self.cells[i][0],
                                           self.cells[i][1]))
        self.pos = {i: p for p, i in enumerate(self.order)}
        # attacking pairs as (earlier, later) index pairs
        self.attack = []
        for a in range(self.n):
            ka, ra, ca, da = self.cells[a]
            for b in range(self.n):
                kb, rb, cb, db = self.cells[b]
                if db == da and ka < kb:
                    self.attack.append((a, b))
                elif db == da + 1 and kb < ka:
                    self.attack.append((a, b))
        # in-component covers for SSYT constraints: right-neighbor weak,
        # up-neighbor strict (French: row above)
        self.weak = []   # T[a] <= T[b]
        self.strict = []  # T[a] < T[b]
        index = {(k, r, c): i for i, (k, r, c, _) in enumerate(self.cells)}
        for i, (k, r, c, _) in enumerate(self.cells):
            j = index.get((k, r, c + 1))
            if j is not None:
                self.weak.append((i, j))
            j = index.get((k, r + 1, c))
            if j is not None:
                self.strict.append((i, j))

    @staticmethod
    def of_partitions(shapes, offsets=None):
        offsets = offsets or [0] * len(shapes)
        return TupleShape([(skew_cells(sh), off)
                           for sh, off in zip(shapes, offsets)])

    @staticmethod
    def of_skews(skews, offsets=None):
        offsets = offsets or [0] * len(skews)
        return TupleShape([(skew_cells(outer, inner), off)
                           for (outer, inner), off in zip(skews, offsets)])

    # ---- standard tableaux (linear extensions), with inv and descent set
    def syt_buckets(self):
        """dict: frozenset(descent set) -> poly in q, summing q^inv."""
        n = self.n
        succ = [[] for _ in range(n)]
        indeg = [0] * n
        for (a, b) in self.weak + self.strict:
            succ[a].append(b)
            indeg[b] += 1
        buckets = {}
        assign = [0] * n  # cell -> value (1..n)
        avail = [i for i in range(n) if indeg[i] == 0]

        def rec(step, avail, indeg):
            if step > n:
                inv = 0
                for (a, b) in self.attack:
                    if assign[a] > assign[b]:
                        inv += 1
                des = frozenset(
                    i for i in range(1, n)
                    if self.pos[assign.index(i + 1)] < self.pos[assign.index(i)])
                buckets[des] = buckets.get(des, Rq(0)) + q ** inv
                return
            for idx in range(len(avail)):
                cell = avail[idx]
                assign[cell] = step
                new_avail = avail[:idx] + avail[idx + 1:]
                new_indeg = list(indeg)
                for t in succ[cell]:
                    new_indeg[t] -= 1
                    if new_indeg[t] == 0:
                        new_avail = new_avail + [t]
                rec(step + 1, new_avail, new_indeg)
                assign[cell] = 0

        rec(1, avail, indeg)
        return buckets

    def g_monomial(self):
        """m-expansion of G_nu as dict {mu (desc tuple): poly in q}.

        [HHL] (82): G = sum_{S in SYT} q^{inv(S)} Q_{n,D(S)}, and
        [x^mu] Q_{n,D} = 1 iff D is contained in the partial sums of mu.
        """
        buckets = self.syt_buckets()
        out = {}
        for mu in Partitions(self.n):
            mu_t = parts_tuple(mu)
            sums = set(itertools.accumulate(mu_t[:-1]))
            coeff = Rq(0)
            for des, poly in buckets.items():
                if des <= sums:
                    coeff += poly
            if coeff:
                out[mu_t] = coeff
        return out

    # ---- direct semistandard enumeration (independent slow route)
    def g_monomial_direct(self):
        """Same m-expansion by enumerating semistandard fillings with
        entries bounded by n, reading off x^mu coefficients directly."""
        n = self.n
        out = {}
        for mu in Partitions(n):
            mu_t = parts_tuple(mu)
            nvals = len(mu_t)
            total = Rq(0)
            # enumerate fillings with content exactly mu
            counts = [0] * (nvals + 1)

            def rec(k, filling):
                nonlocal total
                if k == n:
                    inv = 0
                    for (a, b) in self.attack:
                        if filling[a] > filling[b]:
                            inv += 1
                    total += q ** inv
                    return
                cell = k  # fill cells in index order; constraints checked partially
                for v in range(1, nvals + 1):
                    if counts[v] == mu_t[v - 1]:
                        continue
                    ok = True
                    for (x, y) in self.weak:
                        if y == cell and filling[x] > v:
                            ok = False
                            break
                        if x == cell and y < cell and v > filling[y]:
                            ok = False
                            break
                    if ok:
                        for (x, y) in self.strict:
                            if y == cell and filling[x] >= v:
                                ok = False
                                break
                            if x == cell and y < cell and v >= filling[y]:
                                ok = False
                                break
                    if ok:
                        counts[v] += 1
                        filling.append(v)
                        rec(k + 1, filling)
                        filling.pop()
                        counts[v] -= 1

            rec(0, [])
            if total:
                out[mu_t] = total
        return out

    def max_inv(self):
        """max inv over semistandard fillings = attained on standard ones
        (standardization preserves inv), so read it off the buckets."""
        best = 0
        for poly in self.syt_buckets().values():
            best = max(best, poly.degree())
        return best

    def n_attack(self):
        return len(self.attack)


# --- the ribbon model, [LLT] §4 + [LT] §6 (beta-number strips).
#
# Everything integer-spin: spin_LT(R) = h(R) - 1, spin_LT(T) = sum.
# ([LLT]'s s(T) is half of this; their cospin s~ = (spinmax - spin_LT)/2.)
# Horizontal k-ribbon strips of lambda/nu of weight m <-> add k to m
# distinct beta numbers of nu' to reach beta(lambda'), with
#   spin_LT(strip) = (k-1)*m - #crossings   [LT Lemma 6.5]
# where #crossings = number of unmoved beads jumped over.


def beta_set(la, nrows):
    la = tuple(la) + (0,) * (nrows - len(la))
    return [la[j] + nrows - 1 - j for j in range(nrows)]


def conj(la):
    la = [x for x in la if x > 0]
    if not la:
        return ()
    return tuple(sum(1 for x in la if x > j) for j in range(la[0]))


def strips_down(la, m, k, nrows):
    """All nu with la/nu a horizontal k-ribbon strip of weight m,
    yielding (nu, spin_LT). Works on conjugates via beta moves."""
    lc = conj(la)
    B = beta_set(lc, nrows)
    out = []
    for subset in itertools.combinations(range(nrows), m):
        newB = list(B)
        ok = True
        for i in subset:
            newB[i] -= k
            if newB[i] < 0:
                ok = False
                break
        if not ok or len(set(newB)) != nrows:
            continue
        # crossings: unmoved beads passed over while sliding down by k
        cross = 0
        unmoved = [B[i] for i in range(nrows) if i not in subset]
        for i in subset:
            cross += sum(1 for b in unmoved if B[i] - k < b < B[i])
        srt = sorted(newB, reverse=True)
        nu = tuple(srt[j] - (nrows - 1 - j) for j in range(nrows))
        nu = tuple(x for x in nu if x > 0)
        if any(nu[i] < nu[i + 1] for i in range(len(nu) - 1)):
            continue
        out.append((conj(nu), (k - 1) * m - cross))
    return out


def ribbon_g(la, k):
    """G_LT(lambda; q) = sum over k-ribbon tableaux of shape lambda of
    q^{spin_LT(T)} x^{weight}, as dict {mu: poly}. Requires |la| % k == 0
    and empty k-core (else returns {} when no tableaux exist)."""
    n = sum(la)
    assert n % k == 0
    r = n // k
    nrows = max(len(conj(la)), len(la), r * k) + k
    out = {}
    for mu in Partitions(r):
        mu_t = parts_tuple(mu)
        total = Rq(0)
        # chains of horizontal strips, peeling weights mu_r, ..., mu_1
        def rec(shape, i, spin):
            nonlocal total
            if i < 0:
                if not shape:
                    total += q ** spin
                return
            for (nu, sp) in strips_down(shape, mu_t[i], k, nrows):
                rec(nu, i - 1, spin + sp)

        rec(tuple(la), len(mu_t) - 1, 0)
        if total:
            out[mu_t] = total
    return out


def dict_subs_q(d, sub):
    """Apply q |-> sub (a rational function) to each coefficient, clearing
    to polynomials; sub given as callable poly -> new value in Frac(Rq)."""
    out = {}
    F = FractionField(Rq)
    for k_, v in d.items():
        out[k_] = F(sub(v))
    return out


def normalize_frac_dict(d):
    """Convert Frac(Rq) values that are actually polynomials to Rq."""
    out = {}
    for k_, v in d.items():
        num, den = v.numerator(), v.denominator()
        assert den.is_unit() or num % den == 0, (k_, v)
        out[k_] = num // den if not den.is_unit() else num * den.inverse_of_unit()
    return out


def spin_max(la, k):
    g = ribbon_g(la, k)
    return max(v.degree() for v in g.values())


def ribbon_g_cospin(la, k):
    """LLT97 (26): cospin-normalized  G~ = q^{s*} G_LT(1/q) with exponents
    halved... in spin_LT units: exponent = (spinmax - spin_LT)/2, always
    an integer. Returns dict {mu: poly in q}."""
    g = ribbon_g(la, k)
    smax = max(v.degree() for v in g.values())
    out = {}
    for mu, poly in g.items():
        acc = Rq(0)
        for d in range(poly.degree() + 1):
            c = poly[d]
            if c:
                diff = smax - d
                assert diff % 2 == 0, (la, k, mu, poly)
                acc += c * q ** (diff // 2)
        out[mu] = acc
    return out


def ribbon_g_spin(la, k):
    """LLT97 (28) H-normalization: exponents spin_LT/2 (integral for kmu
    shapes)."""
    g = ribbon_g(la, k)
    out = {}
    for mu, poly in g.items():
        acc = Rq(0)
        for d in range(poly.degree() + 1):
            c = poly[d]
            if c:
                assert d % 2 == 0, (la, k, mu, poly)
                acc += c * q ** (d // 2)
        out[mu] = acc
    return out


# --- k-quotient with content offsets (Stanton-White via the abacus).


def k_quotient_with_offsets(la, k):
    """Return list of (partition, offset) for components r = 0..k-1.

    Beta numbers on k runners; runner r holds beta values b = k*p + r.
    Component r's partition read off runner positions. The content offset
    is pinned empirically by the D-sweep in section D (hypothesis: the
    offset that makes cell (i,j) of component r sit on tuple-diagonal
    m = j - i + off_r, matching the D_m diagonals of the ribbon picture).
    """
    n = sum(la)
    assert n % k == 0
    nrows = k * ((max(len(la), 1) + k) // k)  # multiple of k
    B = beta_set(la, nrows)
    comps = []
    for r in range(k):
        pos = sorted(((b - r) // k for b in B if b % k == r), reverse=True)
        cnt = len(pos)
        part = tuple(pos[i] - (cnt - 1 - i) for i in range(cnt))
        part = tuple(x for x in part if x > 0)
        # offset hypothesis: derived in section D; see spec.
        # count of beads on runners < r ... start with 0, fixed by sweep.
        comps.append([part, 0])
    return comps, nrows


# --- labelled Dyck paths ([HRW] statistics, dyck.rs conventions).


def dyck_paths(n):
    """Area sequences a_1..a_n, a_1 = 0, 0 <= a_i <= a_{i-1}+1."""
    def rec(prefix):
        if len(prefix) == n:
            yield tuple(prefix)
            return
        top = prefix[-1] + 1 if prefix else 0
        for a in range(0, top + 1):
            yield from rec(prefix + [a])
    yield from rec([])


def dyck_g(area, nvals):
    """G_D as dict {mu: poly}: sum over labelings l_1..l_n (l_i > l_{i-1}
    if a_i = a_{i-1}+1) of q^{dinv} x^l, restricted to entries <= nvals,
    coefficient of x^mu."""
    n = len(area)
    out = {}
    for mu in Partitions(n):
        mu_t = parts_tuple(mu)
        if len(mu_t) > nvals:
            continue
        total = Rq(0)
        counts = [0] * (len(mu_t) + 1)

        def rec(i, labels):
            nonlocal total
            if i == n:
                dinv = 0
                for x in range(n):
                    for y in range(x + 1, n):
                        if area[x] == area[y] and labels[x] < labels[y]:
                            dinv += 1
                        if area[x] == area[y] + 1 and labels[x] > labels[y]:
                            dinv += 1
                total += q ** dinv
                return
            for v in range(1, len(mu_t) + 1):
                if counts[v] == mu_t[v - 1]:
                    continue
                if i > 0 and area[i] == area[i - 1] + 1 and v <= labels[-1]:
                    continue
                counts[v] += 1
                labels.append(v)
                rec(i + 1, labels)
                labels.pop()
                counts[v] -= 1

        rec(0, [])
        if total:
            out[mu_t] = total
    return out


def dyck_tuple(area):
    """The vertical-strip tuple of a Dyck path: components = maximal runs
    of consecutive +1 rises; row i sits at content -a_i; a run is a
    vertical strip (cells stacked, one per content). Components in
    REVERSE row order — [DA] Remark 2.2's reversal between the [HHL]
    convention (this script's TupleShape) and the shuffle-world one; with
    it, attacking inversions coincide with [HRW] dinv pairs pointwise
    (verified in section I)."""
    n = len(area)
    comps = []
    start = 0
    for i in range(1, n + 1):
        if i == n or area[i] != area[i - 1] + 1:
            cells = [(r - start, 0) for r in range(start, i)]
            off = -area[start]
            comps.append((cells, off))
            start = i
    return TupleShape(comps[::-1])


# --- unit-interval-graph coloring model ([AS] §2.6-2.8, [DA] §2.2).


def graph_from_area(area):
    """Unit interval graph by the cell-below-path rule ([AS] picture,
    [CM] §3): vertices 1..n along the diagonal, edge {j, i}, j < i, iff
    j >= i - a_i (the cell in column j, row i is between path and
    diagonal). This is a GENERATOR of unit-interval graphs for the
    graph<->graph chromatic bridge (section K); it is NOT the dictionary
    from the path's dinv model to a graph — that is dinv_graph below."""
    n = len(area)
    edges = set()
    for i in range(n):
        for j in range(i - area[i], i):
            edges.add((j, i))
    return n, edges


def dinv_graph(area):
    """Decorated graph faithful to the dinv model of the path: vertices =
    rows 0..n-1; ORIENTED weak edges = dinv pairs — primary (a_x = a_y,
    x < y) emitted as (x, y) so asc counts kappa_x < kappa_y, secondary
    (a_x = a_y + 1, x < y) emitted as (y, x) so asc counts kappa_x >
    kappa_y; strict edges = consecutive same-run rows (i-1, i) with
    a_i = a_{i-1} + 1 (labels strictly increase up a column). Weak and
    strict sets are disjoint: same-run pairs are never dinv pairs."""
    n = len(area)
    edges = set()
    strict = set()
    for x_ in range(n):
        for y_ in range(x_ + 1, n):
            if area[x_] == area[y_]:
                edges.add((x_, y_))
            elif area[x_] == area[y_] + 1:
                edges.add((y_, x_))
    for i in range(1, n):
        if area[i] == area[i - 1] + 1:
            strict.add((i - 1, i))
    return n, edges, strict


def coloring_g(n, edges, strict, nvals_mu):
    """G = sum over colorings kappa (kappa(u) < kappa(v) on strict edges
    u<v) of q^{asc} x^kappa; asc = non-strict edges u<v with
    kappa(u) < kappa(v). m-expansion dict."""
    out = {}
    for mu in Partitions(n):
        mu_t = parts_tuple(mu)
        total = Rq(0)
        counts = [0] * (len(mu_t) + 1)
        kappa = []

        def rec(v):
            nonlocal total
            if v == n:
                asc = sum(1 for (a, b) in edges
                          if (a, b) not in strict and kappa[a] < kappa[b])
                total += q ** asc
                return
            for c in range(1, len(mu_t) + 1):
                if counts[c] == mu_t[c - 1]:
                    continue
                ok = all(not (a == v or b == v) or not ((a, b) in strict)
                         or (kappa[a] < c if b == v else True)
                         for (a, b) in strict)
                if ok:
                    counts[c] += 1
                    kappa.append(c)
                    rec(v + 1)
                    kappa.pop()
                    counts[c] -= 1

        rec(0)
        if total:
            out[mu_t] = total
    return out


# ---------------------------------------------------------------- helpers


def sage_m_dict(el):
    """Sage symmetric function -> {mu: value in Rqt fraction field}."""
    el = m_b(el)
    out = {}
    for mu, c in el.monomial_coefficients().items():
        out[parts_tuple(mu)] = c
    return out


def to_sage_m(d, var=qq):
    """my dict {mu: poly in q} -> Sage m-basis element with q -> var."""
    el = m_b.zero()
    for mu, poly in d.items():
        el += m_b(Partition(list(mu))) * poly(var)
    return el


def dicts_equal(d1, d2):
    keys = set(d1) | set(d2)
    for k_ in keys:
        v1 = d1.get(k_, 0)
        v2 = d2.get(k_, 0)
        if Rqt(str(v1)) != Rqt(str(v2)):
            return False
    return True


def check(name, cond):
    if not cond:
        print(f"FAIL  {name}")
        sys.exit(1)
    print(f"PASS  {name}")


# ================================================================ sections


def gate():
    # hand fixture: nu = ((1),(1)) both at content 0:
    # G = m2 + (1+q) m11
    T = TupleShape.of_partitions([(1,), (1,)])
    g = T.g_monomial()
    ok = g == {(2,): Rq(1), (1, 1): 1 + q}
    check("gate: G_((1),(1)) = m2 + (1+q)m11", ok)
    # and the direct route agrees
    check("gate: standardization == direct on ((1),(1))",
          g == T.g_monomial_direct())


def section_A():
    tuples = [
        [(2,), (1,)],
        [(1, 1), (1,)],
        [(2, 1), (1,)],
        [(1,), (1,), (1,)],
        [(2,), (2,)],
        [(2, 1), (2,)],
    ]
    skews = [
        [((2, 1), (1,)), ((1,), ())],
        [((2, 2), (1,)), ((2,), (1,))],
    ]
    for shapes in tuples:
        T = TupleShape.of_partitions(shapes)
        assert T.g_monomial() == T.g_monomial_direct(), shapes
    for sk in skews:
        T = TupleShape.of_skews(sk)
        assert T.g_monomial() == T.g_monomial_direct(), sk
    check("A: SYT/descent route == direct SSYT route (8 tuples)", True)
    # offsets change G (contents matter): shifted second component
    T0 = TupleShape.of_partitions([(1,), (1,)], [0, 0])
    T1 = TupleShape.of_partitions([(1,), (1,)], [0, 5])
    check("A: content offsets are load-bearing (G differs)",
          T0.g_monomial() != T1.g_monomial())


def section_B():
    # Sage's tuple entry point: level-k cospin on a list of partitions.
    # Dictionary to pin: Sage cospin(tuple) vs my G_nu (all offsets 0).
    LL = {}
    for k in (2, 3, 4):
        LL[k] = Sym.llt(k, t=qq)
    cases = [
        [(1,), (1,)],
        [(2,), (1,)],
        [(1, 1), (1,)],
        [(2,), (2,)],
        [(2, 1), (1,)],
        [(1,), (1,), (1,)],
        [(2,), (1,), (1,)],
        [(2, 1), (2,), (1,)],
        [(2, 2), (2, 1)],
    ]
    for shapes in cases:
        k = len(shapes)
        T = TupleShape.of_partitions(shapes)
        mine = T.g_monomial()
        sage_el = LL[k].cospin([list(sh) for sh in shapes])
        sage_d = sage_m_dict(sage_el)
        # the one dictionary: Sage cospin = q^{-min inv} G_nu
        minv = min(min(e for e in range(v.degree() + 1) if v[e])
                   for v in mine.values())
        scaled = {k_: v // q ** minv for k_, v in mine.items()}
        assert dicts_equal({k_: v(qq) for k_, v in scaled.items()}, sage_d), \
            (shapes, minv, mine, sage_d)
        if minv:
            print(f"      [B] {shapes}: min inv = {minv}")
    check("B: Sage llt.cospin(tuple) = q^(-min inv) * G_nu (9 tuples)", True)


def section_C():
    # LLT97 Example 6.8(i): shape (3,3,3,2,1), k=3, cospin normalization
    want_m = {(3, 1): Rq(1), (2, 2): 1 + q, (2, 1, 1): 2 + 2 * q + q ** 2,
              (1, 1, 1, 1): 3 + 5 * q + 3 * q ** 2 + q ** 3}
    got = ribbon_g_cospin((3, 3, 3, 2, 1), 3)
    check("C: LLT97 Ex 6.8(i) G~_33321 monomial fixture", got == want_m)
    got_s = sage_m_dict(s_b(to_sage_m(got)))  # convert via Sage for schur
    el = s_b(to_sage_m(got))
    d = {parts_tuple(mu): c for mu, c in el.monomial_coefficients().items()}
    want_s = {(3, 1): Rqt(1), (2, 2): qq, (2, 1, 1): qq + qq ** 2,
              (1, 1, 1, 1): qq ** 3}
    check("C: LLT97 Ex 6.8(i) Schur expansion s31+qs22+(q+q^2)s211+q^3s1111",
          d == want_s)

    # LLT97 Example 6.8(ii): H^(k)_3211 for k = 2,3,4; k=4 equals Qp_3211.
    def h_spin_schur(mu, k):
        la = tuple(k * x for x in mu)
        g = ribbon_g_spin(la, k)
        el = s_b(to_sage_m(g))
        return {parts_tuple(nu): c
                for nu, c in el.monomial_coefficients().items()}

    h2 = h_spin_schur((3, 2, 1, 1), 2)
    want_h2 = {(3, 2, 1, 1): Rqt(1), (3, 2, 2): qq, (3, 3, 1): qq,
               (4, 1, 1, 1): qq, (4, 2, 1): qq + qq ** 2, (4, 3): qq ** 2,
               (5, 1, 1): qq ** 2, (5, 2): qq ** 3}
    check("C: LLT97 Ex 6.8(ii) H^(2)_3211 Schur fixture", h2 == want_h2)

    h4 = h_spin_schur((3, 2, 1, 1), 4)
    HLQp = Sym.hall_littlewood(t=qq).Qp()
    qp = s_b(HLQp[3, 2, 1, 1])
    want_h4 = {parts_tuple(nu): c
               for nu, c in qp.monomial_coefficients().items()}
    check("C: LLT97 Thm 6.6: H^(4)_3211 = Q'_3211 (Sage HL Qp oracle)",
          h4 == want_h4)

    # Sage basis dictionaries, k = 2 and 3 sweeps
    for k in (2, 3):
        L = Sym.llt(k, t=qq)
        HSp, HCosp = L.hspin(), L.hcospin()
        for n in range(1, 5):
            for mu in Partitions(n):
                mu_t = parts_tuple(mu)
                la = tuple(k * x for x in mu_t)
                mine_sp = to_sage_m(ribbon_g_spin(la, k))
                mine_co = to_sage_m(ribbon_g_cospin(la, k))
                assert m_b(HSp[mu]) == mine_sp, (k, mu)
                assert m_b(HCosp[mu]) == mine_co, (k, mu)
    check("C: HSp_k[mu] = spin/2-normalized, HCosp_k[mu] = cospin ribbon G "
          "(k=2,3, |mu| <= 4)", True)

    # Sage cospin on a plain partition = my cospin ribbon function
    L3 = Sym.llt(3, t=qq)
    for la in [(3, 3, 3, 2, 1), (4, 3, 2), (6, 3), (2, 2, 2)]:
        mine = to_sage_m(ribbon_g_cospin(la, 3))
        assert m_b(L3.cospin(Partition(list(la)))) == mine, la
    check("C: Sage cospin(lambda) == cospin ribbon model (k=3 shapes)", True)


def section_D():
    # quotient dictionary: my tuple G on the k-quotient (offsets 0) equals
    # the cospin ribbon G~ DIRECTLY after dividing by q^{min inv}. (First
    # guess was [HHL]'s remark-shaped q^e G(1/q) — wrong: their G~ symbol
    # is the spin-flavored function, not [LLT]'s cospin G~. Recorded trap.)
    F = FractionField(Rq)
    todo = []
    for k in (2, 3):
        for n in (k, 2 * k, 3 * k):
            for la in Partitions(n):
                la_t = parts_tuple(la)
                comps, _ = k_quotient_with_offsets(la_t, k)
                if all(not p for p, _ in comps):
                    continue
                # empty k-core test: quotient sizes sum to n/k
                if sum(sum(p) for p, _ in comps) * k != n:
                    continue
                todo.append((k, la_t, comps))
    fails = []
    nontrivial = 0
    for (k, la_t, comps) in todo:
        cos = ribbon_g_cospin(la_t, k)
        T = TupleShape.of_partitions([p for p, _ in comps],
                                     [o for _, o in comps])
        mine = T.g_monomial()
        minv = min(min(e for e in range(v.degree() + 1) if v[e])
                   for v in mine.values())
        nontrivial += bool(minv)
        scaled = {k_: v // q ** minv for k_, v in mine.items()}
        if scaled != cos:
            fails.append((k, la_t, comps, mine, cos))
    if fails:
        k, la_t, comps, mine, cos = fails[0]
        print(f"      [D] first failure: k={k} la={la_t} comps={comps}")
        print(f"          tuple G: {mine}")
        print(f"          ribbon cospin: {cos}")
    check(f"D: q^(-min inv) G_(k-quotient, offsets 0) == cospin ribbon G~ "
          f"({len(todo)} empty-core shapes, k=2,3; {nontrivial} with "
          f"min inv > 0)", not fails)
    # LLT97 (28) internal relation H = q^{s*} G~(1/q), plus the recorded
    # trap: the tuple model only sees RELATIVE spin. The absolute grading
    # q^{s* - (maxinv - minv)} is shape data lost by the quotient: e.g.
    # lambda = (1,1,1,1), k=2 has a single tableau (spin_LT 2), tuple
    # G = m11 with maxinv = minv = 0, so no tuple statistic recovers s*=1.
    for (k, la_t, comps) in todo[:20]:
        try:
            hsp = ribbon_g_spin(la_t, k)
        except AssertionError:
            continue  # odd spin_LT: (28)'s q^{1/2} shapes, skipped
        cos = ribbon_g_cospin(la_t, k)
        smax = max(v.degree() for v in ribbon_g(la_t, k).values())
        assert smax % 2 == 0
        sstar = smax // 2
        flip = {k_: Rq((q ** sstar * v(F(1) / q)).numerator())
                for k_, v in cos.items()}
        assert flip == hsp, (k, la_t)
    check("D: LLT97 (28) H = q^{s*} G~(1/q) (even-spin shapes)", True)

    # LT Example 4.1: G_LT((1),(11),(1); q) for mu = (3,3,3,2,1)
    g = ribbon_g((3, 3, 3, 2, 1), 3)
    el = s_b(to_sage_m(g))
    d = {parts_tuple(nu): c for nu, c in el.monomial_coefficients().items()}
    want = {(3, 1): qq ** 7, (2, 2): qq ** 5, (2, 1, 1): qq ** 5 + qq ** 3,
            (1, 1, 1, 1): qq}
    check("D: LT Ex 4.1 spin-generating G matches (q^7 s31 + q^5 s22 + ...)",
          d == want)


def section_E():
    # q = 1: G_nu(x;1) = prod of skew Schurs (Sage oracle)
    cases = [
        [((2, 1), ()), ((1,), ())],
        [((2, 2), (1,)), ((2, 1), ())],
        [((3, 1), (1,)), ((2,), ()), ((1, 1), ())],
    ]
    for sk in cases:
        T = TupleShape.of_skews(sk)
        mine = to_sage_m(T.g_monomial(), var=Rqt(1))
        want = prod(s_b(Partition(list(o))).skew_by(s_b(Partition(list(i))))
                    for (o, i) in sk)
        assert m_b(want) == mine, sk
    check("E: G_nu(x;1) = product of skew Schurs (3 cases)", True)


def section_F():
    # omega duality [HHL] Lemma 10.1 / (83):
    # G_{nu'}(x;q) = q^a omega G_nu(x; 1/q), a = #attacking pairs,
    # nu' = transpose each component and reverse the tuple.
    F = FractionField(Rq)
    cases = [
        [(1,), (1,)],
        [(2,), (1,)],
        [(2, 1), (1,)],
        [(2,), (2,)],
    ]
    for shapes in cases:
        T = TupleShape.of_partitions(shapes)
        a = T.n_attack()
        # transpose components, reverse order; offsets: negate (transpose
        # flips content sign) — pinned by this very check.
        shapes_t = [conj(sh) for sh in shapes][::-1]
        Tt = TupleShape.of_partitions(shapes_t)
        lhs = to_sage_m(Tt.g_monomial())
        g = T.g_monomial()
        rhs = m_b.zero()
        for mu, poly in g.items():
            val = q ** a * poly(F(1) / q)
            rhs += m_b(Sym.monomial()(Partition(list(mu))).omega()) \
                * Rq(val.numerator())(qq)
        assert lhs == m_b(rhs), shapes
    check("F: omega duality G_{nu'} = q^{#attack} omega G_nu(1/q) (4 cases)",
          True)


def section_G():
    # BHMPS (11) coproduct: G_nu[X+Y] = sum q^{A} G_{nu'}[X] G_{nu''}[Y]
    # over order-ideal splits. Verified as: coefficient extraction over
    # two finite alphabets x_1, x_2 vs y_1 — small case.
    # (implemented as: check G(x1,x2,y1) expansion equality)
    shapes = [(2,), (1,)]
    T = TupleShape.of_partitions(shapes)
    # brute force both sides over alphabet split {1,2} + {3}:
    # LHS: fillings with entries in 1..3, x-vars for 1,2, y for 3
    # RHS: split cells into ideal (gets 1..2) / filter (gets 3), count
    # attacking pairs (a in upper, b in lower), q^A, and SSYT each part.
    n = T.n
    from collections import defaultdict
    lhs = defaultdict(lambda: Rq(0))

    def all_ssyt(maxv):
        res = []

        def rec(k, filling):
            if k == n:
                res.append(tuple(filling))
                return
            for v in range(1, maxv + 1):
                ok = True
                for (x, y) in T.weak:
                    if y == k and filling[x] > v:
                        ok = False
                for (x, y) in T.strict:
                    if y == k and filling[x] >= v:
                        ok = False
                if ok:
                    filling.append(v)
                    rec(k + 1, filling)
                    filling.pop()

        rec(0, [])
        return res

    for filling in all_ssyt(3):
        inv = sum(1 for (a, b) in T.attack if filling[a] > filling[b])
        key = (tuple(sorted((v for v in filling if v <= 2), reverse=True)),
               sum(1 for v in filling if v == 3))
        lhs[(key, )] = lhs[(key,)] + q ** inv
    # RHS over ideal splits: lower ideal = downward-closed cell set per comp
    rhs = defaultdict(lambda: Rq(0))
    cellidx = list(range(n))
    import itertools as it
    for mask in range(1 << n):
        lower = [i for i in cellidx if mask >> i & 1]
        upper = [i for i in cellidx if not (mask >> i & 1)]
        # lower must be downward closed under weak+strict covers
        okl = all(not (b in lower and a not in lower)
                  for (a, b) in T.weak + T.strict)
        if not okl:
            continue
        A = sum(1 for (a, b) in T.attack if a in upper and b in lower)
        # G_lower over {1,2}, G_upper over {3} (single variable: weight only)
        def sub_ssyt(cells, maxv, offset=0):
            res = []

            def rec(kk, filling):
                if kk == len(cells):
                    res.append(dict(zip(cells, filling)))
                    return
                cell = cells[kk]
                for v in range(1 + offset, maxv + 1 + offset):
                    ok = True
                    for (x, y) in T.weak:
                        if y == cell and x in cells[:kk + 1] and \
                                dict(zip(cells, filling + [v])).get(x, 0) > v:
                            ok = False
                    for (x, y) in T.strict:
                        if y == cell and x in cells[:kk + 1] and \
                                dict(zip(cells, filling + [v])).get(x, -1) >= v:
                            ok = False
                    if ok:
                        rec(kk + 1, filling + [v])

            rec(0, [])
            return res

        for fl in sub_ssyt(sorted(lower), 2):
            invl = sum(1 for (a, b) in T.attack
                       if a in lower and b in lower
                       and fl[a] > fl[b])
            for fu in sub_ssyt(sorted(upper), 1, offset=2):
                invu = sum(1 for (a, b) in T.attack
                           if a in upper and b in upper and fu[a] > fu[b])
                key = (tuple(sorted(fl.values(), reverse=True)),
                       len(upper))
                rhs[(key,)] = rhs[(key,)] + q ** (A + invl + invu)
    check("G: BHMPS coproduct law on ((2),(1)) with X={x1,x2}, Y={y1}",
          dict(lhs) == dict(rhs))


def section_H():
    # HHL: H~_mu(x;q,t) = sum_D q^{-a(D)} t^{maj(D)} G_{nu(mu,D)}(x;q)
    # nu(mu,D): k = mu_1 ribbons, nu^(j) has size mu'_j, HHL-contents
    # 1..mu'_j (row-col), descent set from column j of D.
    Mac = Sym.macdonald(q=qq, t=tt)
    Ht = Mac.Ht()

    def ribbon_from_descents(size, des):
        """cells (row, col) 0-based with MY content (col-row) equal to
        -1, -2, ..., -size read upward; HHL content c_H = row-col = 1..size.
        Build: start cell (row 0, col -1+0)... place cell for c_H = 1 at
        (1, 0): row-col = 1. c_H = c+1: if c+1 in des: directly above;
        else: directly left."""
        cells = [(1, 0)]
        for c in range(2, size + 1):
            (r, cc) = cells[-1]
            if c in des:
                cells.append((r + 1, cc))
            else:
                cells.append((r + 1, cc + 1))  # up-left keeps skew shape:
                # row-col increases by 1 going up... (r+1)-(cc+1) = r-cc ✗
        # fix: moving "left" in HHL means next content cell is left
        # neighbor: (r, cc-1): row-col = r-cc+1 ✓
        cells = [(1, 0)]
        for c in range(2, size + 1):
            (r, cc) = cells[-1]
            if c in des:
                cells.append((r + 1, cc))
            else:
                cells.append((r, cc - 1))
        return cells

    for n in range(1, 6):
        for mu in Partitions(n):
            mu_t = parts_tuple(mu)
            k = mu_t[0]
            muc = conj(mu_t)
            # descent candidate cells: (i,j) in mu with i > 1 (1-based rows)
            cand = [(i, j) for j in range(1, k + 1)
                    for i in range(2, muc[j - 1] + 1)]
            total = m_b.zero()
            for rmask in range(1 << len(cand)):
                D = [cand[i] for i in range(len(cand)) if rmask >> i & 1]
                # arm/leg in mu (French): arm = cells right in row,
                # leg = cells above in column
                def arm(cell):
                    (i, j) = cell
                    return mu_t[i - 1] - j

                def leg(cell):
                    (i, j) = cell
                    return muc[j - 1] - i

                aD = sum(arm(c) for c in D)
                majD = sum(leg(c) + 1 for c in D)
                comps = []
                for j in range(1, k + 1):
                    size = muc[j - 1]
                    des = set(i for (i, jj) in D if jj == j)
                    cells = ribbon_from_descents(size, des)
                    comps.append((cells, 0))
                T = TupleShape(comps)
                g = T.g_monomial()
                el = to_sage_m(g)
                total += el * qq ** (-aD) * tt ** majD
            want = m_b(Ht[mu])
            assert m_b(total) == want, (mu_t, m_b(total), want)
    check("H: HHL decomposition sum_D q^-a t^maj G_nu(mu,D) = H~_mu "
          "(all mu |- n <= 5)", True)


def section_I():
    # shuffle: nabla e_n = sum_D t^{area} G_D(x;q), G_D from [HRW] dinv
    for n in range(1, 6):
        total = m_b.zero()
        for area in dyck_paths(n):
            g = dyck_g(area, n)
            total += to_sage_m(g) * tt ** sum(area)
        want = m_b(e_b[n].nabla())
        assert m_b(total) == want, n
    check("I: sum_D t^area G_D(dinv model) = nabla e_n (n <= 5, Sage nabla)",
          True)
    # and the tuple construction reproduces the per-path dinv model after
    # the cospin normalization q^{-min inv} — the same q-floor that turned
    # up in sections B and D. The dinv statistic starts at 0 per path; the
    # anchored tuple carries forced inversions.
    for n in range(1, 7):
        for area in dyck_paths(n):
            g1 = dyck_tuple(area).g_monomial()
            g2 = dyck_g(area, n)
            assert g1 == g2, (area, g1, g2)
    check("I: reverse-run tuple G == dinv model pointwise (all D, n <= 6)",
          True)


def section_J():
    # coloring model == dinv model, via the dinv-faithful decorated graph
    # (weak edges = oriented dinv pairs, strict edges = same-run pairs).
    # NOT via graph_from_area: the cell-below-path graph of a is a
    # different object (e.g. a=(0,0) has no cell below the path but rows
    # 1,2 form a primary dinv pair).
    for n in range(1, 7):
        for area in dyck_paths(n):
            nn, edges, strict = dinv_graph(area)
            g1 = coloring_g(nn, edges, strict, None)
            g2 = dyck_g(area, n)
            assert g1 == g2, (area, g1, g2)
    check("J: coloring model on dinv graph (oriented ascents, strict "
          "runs) == dinv model (all D, n <= 6)", True)


def section_K():
    # chromatic bridge: X_P(x;q) = (q-1)^{-n} G_P[x(q-1); q]  [CM Prop 3.5]
    # G -> p-basis; p_r |-> p_r * (q^r - 1); divide by (q-1)^n; compare to
    # Sage's chromatic_quasisymmetric_function of the unit interval graph.
    from sage.graphs.graph import Graph
    for n in range(2, 6):
        for area in dyck_paths(n):
            nn, edges = graph_from_area(area)
            g = coloring_g(nn, edges, set(), None)
            el = p_b(to_sage_m(g))
            F = el.parent().base_ring()
            tot = p_b.zero()
            for rho, c in el.monomial_coefficients().items():
                fac = prod(qq ** r - 1 for r in rho)
                tot += p_b(rho) * c * F(fac)
            tot = tot / F((qq - 1) ** n)
            # vertex set explicit: isolated vertices must survive
            G = Graph([list(range(1, nn + 1)),
                       [(a + 1, b + 1) for (a, b) in edges]],
                      format='vertices_and_edges')
            X = G.chromatic_quasisymmetric_function(t=qq)
            Xs = p_b(X.to_symmetric_function())
            assert p_b(tot) == Xs, (area,)
    check("K: X_Gamma = (q-1)^{-n} G[x(q-1)] vs Sage chromatic QSF "
          "(all unit-interval graphs, n <= 5)", True)


def section_L():
    # AS fixtures. Path words in {n,d,e}; vertices = diagonal cells.
    # nndee (n=3): G = q^2 s111 + q s21  [AS Ex 6.1]
    # graph: from the Schroder path; strict edge from the diagonal step.
    # Decorated graph of nndee: n=3 vertices; edges/strictness derived
    # from the path per [AS] §2.6 (cells below path; diagonal-step
    # endpoints are strict).
    def graph_from_word(word):
        # follow [AS]: vertices 1..n numbered along the diagonal; for
        # u < v edge iff cell (column u, row v) lies below the path;
        # diagonal step ending at (u, v) makes uv strict.
        # We implement via the path's column heights.
        steps = list(word)
        x = y = 0
        below = set()
        strict = set()
        # record the path as the set of lattice points visited
        colheight = {}  # for column c (1..n): y of path when x goes c-1->c
        diag_end = []
        for st in steps:
            if st == 'n':
                y += 1
            elif st == 'e':
                colheight[x + 1] = y
                x += 1
            else:  # d
                colheight[x + 1] = y + 1  # diagonal covers the step
                x += 1
                y += 1
                diag_end.append((x, y))
        n = x
        edges = set()
        for u in range(1, n + 1):
            for v in range(u + 1, n + 1):
                # cell in column u, row v below the path: the path at
                # column u has height colheight[u]; cell row v is below
                # path iff v <= colheight[u]
                if v <= colheight[u]:
                    edges.add((u - 1, v - 1))
        for (u, v) in diag_end:
            # endpoint of diagonal step at lattice point (u, v): strict
            # edge {u, v}
            strict.add((u - 1, v - 1))
            edges.add((u - 1, v - 1))
        return n, edges, strict

    n, edges, strict = graph_from_word("nndee")
    g = coloring_g(n, edges, strict, None)
    el = s_b(to_sage_m(g))
    want = qq ** 2 * s_b[1, 1, 1] + qq * s_b[2, 1]
    check("L: AS Ex 6.1  G_nndee = q^2 s111 + q s21", s_b(el) == s_b(want))

    # D'Adderio Ex 5.6: word (-,0,-,0,+,+) i.e. n d n d e e ->
    # d_P(1) = q e1 e3 + q(q-1) e4
    n, edges, strict = graph_from_word("ndndee")
    g = coloring_g(n, edges, strict, None)
    el = e_b(to_sage_m(g))
    want = qq * e_b[3, 1] + qq * (qq - 1) * e_b[4]
    check("L: D'Adderio Ex 5.6  G_ndndee = q e1e3 + q(q-1)e4",
          e_b(el) == e_b(want))

    # AS main theorem on these fixtures: orientation formula
    def orientation_e(n, edges, strict):
        edges_l = sorted(edges)
        free = [e for e in edges_l if e not in strict]
        acc = {}
        for mask in range(1 << len(free)):
            # orientation: strict edges natural u->v; free edge oriented
            # u->v if bit set else v->u; ascending = natural-oriented free
            asc = 0
            adj = {v: [] for v in range(n)}
            for (u, v) in strict:
                adj[u].append(v)
            for i, (u, v) in enumerate(free):
                if mask >> i & 1:
                    adj[u].append(v)
                    asc += 1
                else:
                    adj[v].append(u)
            # hrv via directed reachability on strict+ascending edges:
            # [AS] uses strict and ascending edges only — here all edges
            # oriented; paths use edges u->v with v>u?? — no: "directed
            # path from u to v using only strict and ascending edges";
            # ascending = u->v with u<v, uv not strict. So drop descending.
            radj = {v: [] for v in range(n)}
            for (u, v) in strict:
                radj[u].append(v)
            for i, (u, v) in enumerate(free):
                if mask >> i & 1:
                    radj[u].append(v)
            hrv = [0] * n
            for u in range(n):
                seen = {u}
                stack = [u]
                while stack:
                    w = stack.pop()
                    for z in radj[w]:
                        if z not in seen:
                            seen.add(z)
                            stack.append(z)
                hrv[u] = max(seen)
            blocks = {}
            for u in range(n):
                blocks.setdefault(hrv[u], []).append(u)
            lam = tuple(sorted((len(b) for b in blocks.values()),
                               reverse=True))
            acc[lam] = acc.get(lam, Rq(0)) + q ** asc
        return acc

    n, edges, strict = graph_from_word("nndee")
    ori = orientation_e(n, edges, strict)
    # \hat G(x; q+1) = sum q^{asc} e_{lambda}: compare after q -> q-1
    el2 = e_b.zero()
    for lam, poly in ori.items():
        el2 += e_b(Partition(list(lam))) * poly(qq - 1)
    check("L: AS orientation formula reproduces G_nndee",
          e_b(el2) == e_b(to_sage_m(coloring_g(*graph_from_word('nndee'),
                                               None))))


def section_M():
    # Fock straightening [KMS] (43)/(45): compute V_beta route coefficients
    # and check S_lambda|rho> coefficient = c^lambda_mu vs ribbon model.
    # Wedge = tuple of strictly decreasing ints, coefficient in Frac(Rq).
    # Normal ordering with the KMS rules; q_KMS here is written v.
    Rv = QQ['v']
    v = Rv.gen()
    Fv = FractionField(Rv)

    def straighten(wedge, coeff, k, out):
        # find first ascent (w[i] <= w[i+1])
        w = list(wedge)
        for i in range(len(w) - 1):
            if w[i] < w[i + 1]:
                l_, m_ = w[i], w[i + 1]
                diff = m_ - l_
                if diff % k == 0:
                    nw = w[:i] + [m_, l_] + w[i + 2:]
                    straighten(nw, -coeff, k, out)
                else:
                    im = diff % k
                    # KMS (45): u_l ^ u_m = -v u_m ^ u_l + (v^2-1)(...)
                    nw = w[:i] + [m_, l_] + w[i + 2:]
                    straighten(nw, -v * coeff, k, out)
                    sign = Fv(v ** 2 - 1)
                    # offsets alternate: i, k, k+i, 2k, 2k+i, ...
                    cand = sorted(set(
                        [j * k for j in range(1, diff // k + 1)] +
                        [j * k + im for j in range(0, diff // k + 1)]))
                    cand = [a for a in cand if 0 < a and 2 * a < diff]
                    for t_, a in enumerate(cand):
                        c2 = sign * (-v) ** t_ * coeff
                        nw2 = w[:i] + [m_ - a, l_ + a] + w[i + 2:]
                        straighten(nw2, c2, k, out)
                    return
                return
            if w[i] == w[i + 1]:
                return  # zero
        key = tuple(w)
        out[key] = out.get(key, Fv(0)) + coeff
        return

    def v_op(wedges, m_, k):
        """apply V_m = h_m(Y^{-1}): add k to entries, all weak comps."""
        out = {}
        for wedge, coeff in wedges.items():
            r = len(wedge)
            for comp in itertools.combinations_with_replacement(range(r), m_):
                nw = list(wedge)
                for j in comp:
                    nw[j] += k
                straighten(nw, coeff, k, out)
        return {k_: c for k_, c in out.items() if c}

    # check c^lambda coefficients for k=2, deg <= 4 (deg 4 has multi-term
    # KL entries, e.g. q^5+q^3 on s_211 — pins the variable dictionary
    # beyond what single monomials can).
    # S_lambda|rho> = sum_nu kappa_{lambda,nu} V_nu |rho>, kappa: s in h.
    r = 10
    rho = tuple(range(r - 1, -1, -1))
    for lam in [(2,), (1, 1), (3,), (2, 1),
                (4,), (3, 1), (2, 2), (2, 1, 1), (1, 1, 1, 1)]:
        deg = sum(lam)
        el = h_b(s_b(Partition(list(lam))))
        acc = {}
        for nu, c in el.monomial_coefficients().items():
            wedges = {rho: Fv(QQ(c))}
            for part in nu:
                wedges = v_op(wedges, part, 2)
            for k_, cc in wedges.items():
                acc[k_] = acc.get(k_, Fv(0)) + cc
        # coefficient of |mu + rho| = c^lam_{2-quotient(mu)}(-v^{-1}):
        # compare with ribbon model: G_LT s_lam coefficient c^lam_mu(q).
        for mu in Partitions(2 * deg):
            mu_t = parts_tuple(mu)
            key = tuple(x + y for x, y in
                        zip(tuple(mu_t) + (0,) * (r - len(mu_t)), rho))
            got = acc.get(key, Fv(0))
            # ribbon side
            g = ribbon_g(mu_t, 2) if sum(mu_t) % 2 == 0 else {}
            el_s = s_b(to_sage_m(g)) if g else s_b.zero()
            cmu = el_s.coefficient(Partition(list(lam))) if g else 0
            # dictionary: fock coefficient (KMS variable v) equals the
            # ribbon-model c^lam_mu in the spin_LT grading at q = -v.
            # (Measured from the deg <= 4 table; the multi-term deg-4
            # entries rule out the mirrored laws.)
            from sage.all import SR
            # SR.var (not var()): var() injects into globals and would
            # clobber this module's Rq generator q.
            v_sr, q_sr = SR.var('v'), SR.var('q')
            lhs = SR(str(got))
            rhs = SR(str(cmu)).subs({q_sr: -v_sr})
            assert (lhs - rhs).simplify_rational() == 0, \
                (lam, mu_t, got, cmu)
    check("M: KMS straightening S_lam|rho> == ribbon c^lam_mu(q) at "
          "q = -v (k=2, deg <= 4, incl. multi-term KL entries)", True)


if __name__ == "__main__":
    gate()
    section_A()
    section_B()
    section_C()
    section_D()
    section_E()
    section_F()
    section_G()
    section_H()
    section_I()
    section_J()
    section_K()
    section_L()
    section_M()
    print("all formulas verified")
