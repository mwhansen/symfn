# Generate the Sage oracle fixture consumed by tests/sage_oracle.rs.
#
#   sage scripts/gen_sage_oracle.sage > tests/fixtures/sage_oracle.txt
#
# Sage is the independent implementation: it computes the same quantities by its
# own algorithms, so agreement is real cross-validation rather than symfn
# checking itself. The fixture is committed so `cargo test` needs no Sage.
#
# Line format (partitions are comma-separated, empty partition = empty string):
#   kostka   LAM|MU VALUE
#   char     LAM|MU VALUE
#   smul     MU|NU LAM:COEFF ...
#   s2h      LAM PART:COEFF ...
#   s2e      LAM PART:COEFF ...
#   s2m      LAM PART:COEFF ...
#   s2p      LAM PART:NUM/DEN ...
#   skew     LAM|MU NU:COEFF ...
#   jackp    LAM MU:NUM|DEN ...   (Jack P -> m, coefficients in Q(alpha))
#   jackj    LAM MU:NUM|DEN ...   (Jack J -> m; DEN is always 1)
#   jackjp   LAM MU:NUM|DEN ...   (Jack J -> p, the Jack character table)
#   pleth    F|G PART:COEFF ...
#   hlqp     LAM MU:QTPOLY ...    (Hall-Littlewood Q' -> s)
#   hlp      LAM MU:QTPOLY ...    (Hall-Littlewood P -> s)
#   kf       LAM|MU QTPOLY        (Kostka-Foulkes, zeros included)
#   qtk      LAM|MU QTPOLY        ((q,t)-Kostka)
#   macht    MU PART:QTPOLY ...   (H~ -> s)
#   nablae   N PART:QTPOLY ...    (nabla e_N -> s)
#   kron     LAM|MU|NU VALUE      (Kronecker, zeros included)
#   macp     LAM MU:QTNUM|QTDEN ...  (Macdonald P -> m; likewise macq, macj)
#   lltspin  K|MU PART:QTPOLY ...    (H^(k); likewise lltcospin, lltgtilde)
#   schub    U|V W:COEFF ...      (Schubert structure constants)
#   schubbound N M                (measured: u,v in S_N have support in S_M)
#   cop      LAM MU/NU:COEFF ...  (coproduct, into the tensor square)
#   anti     LAM PART:COEFF ...   (antipode)
#   counit   LAM VALUE
#   dim      LAM VALUE            (f^lambda, counted as standard tableaux)
#   psp      LAM|N VALUE          (s_lambda(1^N))
#   pspq     LAM|N c0,c1,...      (s_lambda(1,q,..,q^{N-1}), dense; Z when zero)
#   redkron  LAM|MU NU:COEFF ... (reduced Kronecker, in the s~ basis)
#
# Jack coefficients are rational FUNCTIONS of alpha, so NUM and DEN are each a
# comma-separated dense list of integer coefficients, index = power of alpha.
# Sage calls the parameter t; it is alpha throughout symfn.
#
# QTPOLY is a polynomial in q and t as 'QE.TE.COEFF;...' sorted, or 'Z' when
# zero. ⚠️ The LLT families are graded in q here, because that is where symfn
# puts them, while Sage names the same parameter t -- see the LLT block.

Sym = SymmetricFunctions(QQ)
s = Sym.schur()
h = Sym.homogeneous()
e = Sym.elementary()
m = Sym.monomial()
p = Sym.powersum()

MAX_N = 6        # kostka / character sweep
MAX_PROD = 6     # |mu| + |nu| for Schur products
MAX_CONV = 5     # degree for basis conversions and skew shapes
MAX_JACK = 7     # degree for the Jack sweep


def enc(lam):
    return ",".join(str(x) for x in lam)


def expansion(elt):
    """Render a symmetric-function element as 'PART:COEFF ...' sorted by part."""
    items = []
    for part, coeff in elt.monomial_coefficients().items():
        items.append((list(part), coeff))
    items.sort()
    out = []
    for part, coeff in items:
        if coeff.denominator() == 1:
            out.append("%s:%s" % (enc(part), coeff.numerator()))
        else:
            out.append("%s:%s/%s" % (enc(part), coeff.numerator(), coeff.denominator()))
    return " ".join(out)


lines = []

# --- Kostka numbers, counted independently as semistandard tableaux ----------
for n in range(0, MAX_N + 1):
    for lam in Partitions(n):
        for mu in Partitions(n):
            k = SemistandardTableaux(list(lam), list(mu)).cardinality()
            lines.append("kostka %s|%s %s" % (enc(lam), enc(mu), k))

# --- Symmetric group characters chi^lam(mu) ---------------------------------
for n in range(0, MAX_N + 1):
    for mu in Partitions(n):
        expanded = s(p[list(mu)])
        for lam in Partitions(n):
            chi = expanded.coefficient(list(lam))
            lines.append("char %s|%s %s" % (enc(lam), enc(mu), chi))

# --- Schur products (encodes every nonzero LR coefficient) ------------------
for a in range(0, MAX_PROD + 1):
    for b in range(0, MAX_PROD + 1 - a):
        for mu in Partitions(a):
            for nu in Partitions(b):
                prod = s(s[list(mu)] * s[list(nu)])
                lines.append("smul %s|%s %s" % (enc(mu), enc(nu), expansion(prod)))

# --- Basis conversions out of Schur ----------------------------------------
for n in range(0, MAX_CONV + 1):
    for lam in Partitions(n):
        sl = s[list(lam)]
        lines.append("s2h %s %s" % (enc(lam), expansion(h(sl))))
        lines.append("s2e %s %s" % (enc(lam), expansion(e(sl))))
        lines.append("s2m %s %s" % (enc(lam), expansion(m(sl))))
        lines.append("s2p %s %s" % (enc(lam), expansion(p(sl))))

# --- Skew Schur functions ---------------------------------------------------
for n in range(0, MAX_CONV + 1):
    for lam in Partitions(n):
        for k in range(0, n + 1):
            for mu in Partitions(k):
                if not Partition(list(lam)).contains(Partition(list(mu))):
                    continue
                sk = s(s[list(lam)].skew_by(s[list(mu)]))
                lines.append("skew %s|%s %s" % (enc(lam), enc(mu), expansion(sk)))

# --- plethysms (degree |f|*|g| grows fast, so bound the product) ------------
for a in range(0, 4):
    for b in range(0, 4):
        if a * b > 8:
            continue
        for f in Partitions(a):
            for g in Partitions(b):
                pl = s(s[list(f)].plethysm(s[list(g)]))
                lines.append("pleth %s|%s %s" % (enc(f), enc(g), expansion(pl)))

print("\n".join(lines))


# --- Jack -------------------------------------------------------------------
#
# Sage has all three normalizations, so this is a real external oracle rather
# than an identity chain.  It is committed as a fixture so `cargo test` checks
# it with no Sage installed -- the live, wider check is scripts/check_jack.py.

JR = QQ["t"]
ALPHA = JR.gen()
JF = JR.fraction_field()
JSym = SymmetricFunctions(JF)
jack = JSym.jack()
jP, jJ = jack.P(), jack.J()
jm, jp = JSym.monomial(), JSym.powersum()


def ratfun(c):
    """A rational function of alpha as 'NUM|DEN', each a dense INTEGER list.

    Sage's numerator()/denominator() over Frac(QQ[t]) are polynomials with
    RATIONAL coefficients -- 1/2*t + 3 is a perfectly ordinary numerator -- so
    both sides are cleared by the lcm of those denominators and then divided by
    the gcd of the result.  Without that the fixture carries '1/2' tokens into
    an integer parser.
    """
    c = JF(c)
    num = JR(c.numerator()).list()
    den = JR(c.denominator()).list()
    scale = lcm([QQ(x).denominator() for x in num + den] + [1])
    num = [ZZ(QQ(x) * scale) for x in num]
    den = [ZZ(QQ(x) * scale) for x in den]
    g = gcd([abs(x) for x in num + den] + [0])
    if g > 1:
        num = [x // g for x in num]
        den = [x // g for x in den]
    # Normalize the sign so the fixture is stable: denominator leading > 0.
    if den and den[-1] < 0:
        num = [-x for x in num]
        den = [-x for x in den]
    enc_num = ",".join(str(x) for x in num) if num else "0"
    enc_den = ",".join(str(x) for x in den) if den else "1"
    return f"{enc_num}|{enc_den}"


def jack_expansion(elt):
    items = sorted((list(part), coeff) for part, coeff in elt.monomial_coefficients().items())
    return " ".join(f"{enc(part)}:{ratfun(coeff)}" for part, coeff in items)


for n in range(0, MAX_JACK + 1):
    for lam in Partitions(n):
        print(f"jackp {enc(lam)} {jack_expansion(jm(jP[lam]))}")
        print(f"jackj {enc(lam)} {jack_expansion(jm(jJ[lam]))}")
        print(f"jackjp {enc(lam)} {jack_expansion(jp(jJ[lam]))}")


# --- The (q,t) layer: Hall-Littlewood, Kostka-Foulkes, (q,t)-Kostka ---------
#
# These are the durable half of each family's oracle: the scripts/check_*.py
# harnesses run wider, this runs on every `cargo test`
# (docs/policies/validation.md V4).
#
# ⚠️ Convention, and the reason all four tags are carried rather than one:
# Sage's Qp() is Q'(x;t), the *modified* Hall-Littlewood whose Schur expansion
# has the Kostka-Foulkes polynomials as its coefficients. P is a different
# object (`Qp[2,1]` is `s[2,1] + t*s[3]`; `P[2,1]` is `s[2,1] - (t^2+t)*s[111]`),
# and neither is Q. Carrying Q', P and K(t) separately means a swap between them
# disagrees on a value rather than producing a plausible table.

MAX_HL = 6       # Hall-Littlewood Q' and P, in the Schur basis
MAX_KF = 6       # Kostka-Foulkes, every (lambda, mu) pair including the zeros
MAX_QTK = 5      # the (q,t)-Kostka table, every pair

QTR = PolynomialRing(QQ, "q,t")


def qtpoly(c):
    """A polynomial in q and t as 'QE.TE.COEFF;...' sorted, or 'Z' when zero.

    The coefficients of every family here are integers, so a rational one is a
    generator bug and must not reach the fixture as a silently-truncated value:
    the assert is the strict-parser rule applied on the writing side.
    """
    items = sorted(QTR(c).dict().items())
    if not items:
        return "Z"
    out = []
    for (a, b), coeff in items:
        assert QQ(coeff).denominator() == 1, f"non-integral coefficient {coeff}"
        out.append(f"{a}.{b}.{QQ(coeff).numerator()}")
    return ";".join(out)


HLR = PolynomialRing(QQ, "t")
HLSym = SymmetricFunctions(HLR.fraction_field())
hlQp, hlP = HLSym.hall_littlewood().Qp(), HLSym.hall_littlewood().P()
hls = HLSym.schur()


def qt_expansion(elt):
    items = sorted((list(part), coeff) for part, coeff in elt.monomial_coefficients().items())
    return " ".join(f"{enc(part)}:{qtpoly(coeff)}" for part, coeff in items)


for n in range(0, MAX_HL + 1):
    for lam in Partitions(n):
        print(f"hlqp {enc(lam)} {qt_expansion(hls(hlQp[list(lam)]))}")
        print(f"hlp {enc(lam)} {qt_expansion(hls(hlP[list(lam)]))}")

# Every pair, including the zeros: a transition that is right on its support and
# wrong about where the support *is* would pass a nonzero-only comparison.
from sage.combinat.sf.kfpoly import KostkaFoulkesPolynomial  # noqa: E402

for n in range(0, MAX_KF + 1):
    for lam in Partitions(n):
        for mu in Partitions(n):
            k = KostkaFoulkesPolynomial(list(lam), list(mu), HLR.gen())
            print(f"kf {enc(lam)}|{enc(mu)} {qtpoly(k)}")

# The (q,t)-Kostka matrix has no zero entries, but its two indices enter
# asymmetrically, so a transposed answer is otherwise a plausible matrix.
from sage.combinat.sf.macdonald import qt_kostka  # noqa: E402

for n in range(0, MAX_QTK + 1):
    for lam in Partitions(n):
        for mu in Partitions(n):
            print(f"qtk {enc(lam)}|{enc(mu)} {qtpoly(qt_kostka(list(lam), list(mu)))}")


# --- H-tilde, nabla, and the Kronecker product -------------------------------
#
# ⚠️ Four normalizations of the modified Macdonald polynomial circulate and the
# fixture is chosen to tell them apart: `Ht[2] = s_2 + q*s_11` against
# `Ht[1,1] = s_2 + t*s_11` is the smallest pair that distinguishes H~ from H,
# from J, and from a q<->t transpose, all of which agree on `Ht[1]`.
#
# nabla is the one Macdonald operator with an external oracle, and Delta, Delta'
# and Theta are all tied back to it (docs/record/macdonald-operators.md), so a
# durable fixture here is load-bearing for that whole family.

MAX_HT = 6       # H~_mu in the Schur basis
MAX_NABLA = 6    # nabla e_n in the Schur basis
MAX_KRON = 5     # every (lambda, mu, nu) triple, including the zeros

QTF = QTR.fraction_field()
QTSym = SymmetricFunctions(QTF)
mHt = QTSym.macdonald().Ht()
qts, qte = QTSym.schur(), QTSym.elementary()

for n in range(0, MAX_HT + 1):
    for mu in Partitions(n):
        print(f"macht {enc(mu)} {qt_expansion(qts(mHt[list(mu)]))}")

for n in range(0, MAX_NABLA + 1):
    en = qte[[n]] if n > 0 else qte[[]]
    print(f"nablae {n} {qt_expansion(qts(en.nabla()))}")

# Every triple, zeros included: the Kronecker product's support is the hard
# part, and `kronecker_via_characters` and the p-basis route would have to
# agree about it independently for a nonzero-only sweep to mean anything.
KSym = SymmetricFunctions(QQ)
ks = KSym.schur()

for n in range(0, MAX_KRON + 1):
    for lam in Partitions(n):
        for mu in Partitions(n):
            prod = ks(ks[list(lam)].itensor(ks[list(mu)]))
            for nu in Partitions(n):
                g = prod.coefficient(list(nu))
                print(f"kron {enc(lam)}|{enc(mu)}|{enc(nu)} {g}")


# --- Macdonald P, Q and J ----------------------------------------------------
#
# Coefficients are rational FUNCTIONS of q and t, so each is emitted as
# 'NUM|DEN' with both sides qtpoly-encoded integer polynomials. symfn keeps its
# denominators factored into binomial atoms and Sage returns them expanded;
# those are different normal forms for the same element, so tests/sage_oracle.rs
# compares by evaluating at generic (q,t) rather than structurally -- the same
# reason check_macdonald.py compares in the fraction field.
#
# ⚠️ P, Q and J differ by b_lambda(q,t) and by c_lambda(q,t) respectively, and
# all three are called "the Macdonald polynomial". All three are carried.

MAX_MAC = 5

MSym = SymmetricFunctions(QTF)
mP, mQ, mJ = MSym.macdonald().P(), MSym.macdonald().Q(), MSym.macdonald().J()
mm = MSym.monomial()


def qtratfun(c):
    """A rational function of q,t as 'NUM|DEN', both qtpoly-encoded.

    Sage's numerator()/denominator() over Frac(QQ[q,t]) carry RATIONAL
    coefficients, so both sides are cleared by the lcm of those denominators
    and reduced by the content -- otherwise a '1/2' token reaches an integer
    parser, which V4 requires be a parse error rather than an absorbed value.
    """
    c = QTF(c)
    num, den = QTR(c.numerator()), QTR(c.denominator())
    scale = lcm(
        [QQ(x).denominator() for x in list(num.coefficients()) + list(den.coefficients())] + [1]
    )
    num, den = num * scale, den * scale
    g = gcd([ZZ(QQ(x)) for x in list(num.coefficients()) + list(den.coefficients())] + [0])
    if g > 1:
        num, den = num // g, den // g
    # Stable sign so a regeneration is a no-op diff: the lex-largest term of the
    # denominator is positive.
    if den.dict() and QQ(den.dict()[max(den.dict())]) < 0:
        num, den = -num, -den
    return f"{qtpoly(num)}|{qtpoly(den)}"


def mac_expansion(elt):
    items = sorted((list(part), coeff) for part, coeff in elt.monomial_coefficients().items())
    return " ".join(f"{enc(part)}:{qtratfun(coeff)}" for part, coeff in items)


for n in range(0, MAX_MAC + 1):
    for lam in Partitions(n):
        print(f"macp {enc(lam)} {mac_expansion(mm(mP[list(lam)]))}")
        print(f"macq {enc(lam)} {mac_expansion(mm(mQ[list(lam)]))}")
        print(f"macj {enc(lam)} {mac_expansion(mm(mJ[list(lam)]))}")


# --- LLT: the ribbon dictionaries -------------------------------------------
#
# ⚠️ Four normalizations of G~ circulate and they agree on the easy cases, so
# the fixture carries three separate dictionaries rather than deriving any from
# another: hspin (H^(k)), hcospin (H~^(k)) and cospin (G~^(k)).
#
# ⚠️ The q/t translation is the trap this block exists to pin. Sage names the
# LLT parameter **t**; symfn grades the LLT families in **q**, with the t slot
# always zero (`llt_dump.rs` asserts it). So the exponent is written into the q
# slot below, and a fixture that put it in t would disagree with symfn on
# everything -- which is the point of writing it once, here, rather than in
# each consumer.
#
# k = 1 is included because hspin at level 1 is the Schur function, which is the
# cheapest place a spin/cospin swap shows.

MAX_LLT = 4
LLT_LEVELS = (1, 2, 3)

LSym = SymmetricFunctions(HLR.fraction_field())
lm = LSym.monomial()
LFAM = {}


def lfamily(k):
    # Built once per level: the constructor is where Sage's coercion wart lives
    # (`llt(k)` over a field whose variable is named q fails to coerce).
    if k not in LFAM:
        LFAM[k] = LSym.llt(k, t=HLR.gen())
    return LFAM[k]


def llt_poly(c):
    """A Sage LLT coefficient (a polynomial in Sage's t) in symfn's q slot."""
    coeffs = HLR(c).list()
    out = []
    for i, co in enumerate(coeffs):
        if co == 0:
            continue
        assert QQ(co).denominator() == 1, f"non-integral LLT coefficient {co}"
        out.append(f"{i}.0.{QQ(co).numerator()}")
    return ";".join(out) if out else "Z"


def llt_expansion(elt):
    items = sorted((list(part), coeff) for part, coeff in lm(elt).monomial_coefficients().items())
    return " ".join(f"{enc(part)}:{llt_poly(coeff)}" for part, coeff in items)


for k in LLT_LEVELS:
    fam = lfamily(k)
    for n in range(0, MAX_LLT + 1):
        for mu in Partitions(n):
            pmu = Partition(list(mu))
            print(f"lltspin {k}|{enc(mu)} {llt_expansion(fam.hspin()[pmu])}")
            print(f"lltcospin {k}|{enc(mu)} {llt_expansion(fam.hcospin()[pmu])}")

# cospin is indexed by the k-ribbon shape itself, so it exists only when k
# divides |lambda|, and Sage rejects the empty shape outright (ValueError) --
# a floor above zero with its reason, per V6.
for k in LLT_LEVELS:
    fam = lfamily(k)
    for n in range(k, MAX_LLT * 2 + 1):
        if n % k:
            continue
        for lam in Partitions(n):
            print(f"lltgtilde {k}|{enc(lam)} {llt_expansion(fam.cospin(Partition(list(lam))))}")


# --- Schubert structure constants -------------------------------------------
#
# c^w_{uv} for Schubert polynomials, as the full expansion of each product.
# schubmult is the other external implementation; Sage's SchubertPolynomialRing
# is a second one, and neither shares authorship with symfn's transition tree.
#
# The `schubbound N M` records are measured, not assumed: for u,v in S_N the
# support never reaches past a one-line length of M. tests/sage_oracle.rs sweeps
# S_M for the zeros with that bound, so the universe the zeros are checked over
# comes from the generator rather than from a guess at the test site.

SX = SchubertPolynomialRing(ZZ)


def enc_perm(w):
    return ",".join(str(x) for x in w)


for n in (1, 2, 3, 4):
    bound = 0
    for u in Permutations(n):
        for v in Permutations(n):
            prod = SX(list(u)) * SX(list(v))
            items = sorted((list(w), c) for w, c in prod.monomial_coefficients().items())
            for w, _ in items:
                bound = max(bound, len(w))
            body = " ".join(f"{enc_perm(w)}:{c}" for w, c in items)
            print(f"schub {enc_perm(u)}|{enc_perm(v)} {body}")
    print(f"schubbound {n} {bound}")


# --- The Hopf structure: coproduct, antipode, counit -------------------------
#
# Sage computes all three, so this is an oracle rather than an identity chain.
# The in-tree evidence is the Hopf axioms plus agreement with the LR route, and
# laws are convention-blind: every self-consistent normalization satisfies them
# (docs/policies/validation.md, evidence class 4 vs 5). What the fixture adds is
# the values.
#
# ⚠️ The antipode's sign is the pin. S(s_λ) = (−1)^{|λ|} s_{λ'} couples a sign
# to a conjugation, and both halves are silent when wrong on a self-conjugate
# shape of even size -- so the sweep runs over whole degrees rather than over
# hand-picked shapes.

MAX_HOPF = 6


def tensor_expansion(elt):
    items = sorted(
        ((list(a), list(b)), c) for (a, b), c in elt.monomial_coefficients().items()
    )
    return " ".join(f"{enc(a)}/{enc(b)}:{c}" for (a, b), c in items)


for n in range(0, MAX_HOPF + 1):
    for lam in Partitions(n):
        sl = s[list(lam)]
        print(f"cop {enc(lam)} {tensor_expansion(sl.coproduct())}")
        print(f"anti {enc(lam)} {expansion(s(sl.antipode()))}")
        print(f"counit {enc(lam)} {sl.counit()}")


# --- Principal specialization, and the dimension -----------------------------
#
# The alphabet size runs from 0, where s_lambda(1^0) is 1 for the empty shape
# and 0 for every other -- a convention over a well-posed question, so it is
# swept rather than skipped (V6).
#
# The dimension comes from StandardTableaux(...).cardinality(), a direct count,
# rather than from Sage's own hook-length formula: symfn computes f^lambda by
# hooks, and an oracle that used the same formula would check the arithmetic
# and nothing else.

MAX_EVAL = 6
EVAL_ALPHABET = 6

PSR = PolynomialRing(QQ, "q")

for n in range(0, MAX_EVAL + 1):
    for lam in Partitions(n):
        plam = Partition(list(lam))
        print(f"dim {enc(lam)} {StandardTableaux(plam).cardinality()}")
        sl = s[list(lam)]
        for k in range(0, EVAL_ALPHABET + 1):
            print(f"psp {enc(lam)}|{k} {sl.principal_specialization(k, q=1)}")
            coeffs = PSR(sl.principal_specialization(k, q=PSR.gen())).list()
            body = ",".join(str(ZZ(c)) for c in coeffs) if coeffs else "Z"
            print(f"pspq {enc(lam)}|{k} {body}")


# --- Reduced (stable) Kronecker coefficients ---------------------------------
#
# Orellana-Zabrocki Theorem 7: an ordinary product in the s~ basis has the
# reduced Kronecker coefficients as its structure constants, so Sage's st()
# basis is a direct oracle for symfn's reduced_kronecker. The in-tree evidence
# is published expansions, the LR top degree, and stability against the
# ordinary Kronecker product -- all real, none of them an independent
# implementation of g-bar itself.
#
# The expansion is inhomogeneous, which is the whole point of the basis: terms
# of every degree up to |lambda|+|mu| appear, and a route that dropped the
# lower-degree tail would still look like a plausible product.

MAX_REDKRON = 3

stb = Sym.st()

for a in range(0, MAX_REDKRON + 1):
    for b in range(0, MAX_REDKRON + 1):
        for lam in Partitions(a):
            for mu in Partitions(b):
                prod = stb(stb[list(lam)] * stb[list(mu)])
                items = sorted((list(nu), c) for nu, c in prod.monomial_coefficients().items())
                body = " ".join(f"{enc(nu)}:{c}" for nu, c in items)
                print(f"redkron {enc(lam)}|{enc(mu)} {body}")
