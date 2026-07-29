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
#
# Jack coefficients are rational FUNCTIONS of alpha, so NUM and DEN are each a
# comma-separated dense list of integer coefficients, index = power of alpha.
# Sage calls the parameter t; it is alpha throughout symfn.

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
for n in range(1, MAX_N + 1):
    for lam in Partitions(n):
        for mu in Partitions(n):
            k = SemistandardTableaux(list(lam), list(mu)).cardinality()
            lines.append("kostka %s|%s %s" % (enc(lam), enc(mu), k))

# --- Symmetric group characters chi^lam(mu) ---------------------------------
for n in range(1, MAX_N + 1):
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
for a in range(1, 4):
    for b in range(1, 4):
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


for n in range(1, MAX_JACK + 1):
    for lam in Partitions(n):
        print(f"jackp {enc(lam)} {jack_expansion(jm(jP[lam]))}")
        print(f"jackj {enc(lam)} {jack_expansion(jm(jJ[lam]))}")
        print(f"jackjp {enc(lam)} {jack_expansion(jp(jJ[lam]))}")
