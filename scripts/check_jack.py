"""Hold `jack` and `afrac` to Sage, and to the laws Sage cannot state.

    cargo run --release --example jack_dump -- 7 5 > /tmp/jack.txt
    sage -python scripts/check_jack.py /tmp/jack.txt

Sage is used as a black box: run as a separate program, its output compared,
its source not read (see ../NOTICE.md).  Sage calls the Jack parameter `t`;
it is alpha throughout here and in `docs/record/jack.md`.

Unlike the Delta-operator check, Sage *does* have all three normalizations
(P, Q, J), `scalar_jack`, and `zonal` — so most of this is a real external
oracle rather than an identity chain.  Two oracle warts were measured before
relying on it (`spec_jack_verify.py`, and re-checked here):

  * `scalar_jack` cannot RETURN zero — a vanishing pairing raises TypeError
    from inside Sage.  So the crash set is treated as the zero set and is
    checked to coincide exactly with the pairs that must vanish.
  * `zonal()` is `P^(2)`, NOT `J^(2)`.  The two differ by H_lambda(2), and a
    fixture that gets this backwards still looks plausible.

Three things here have no oracle anywhere and are checked as laws:

  * [KS] Thm 1.1 — every [m_mu]J_lambda lies in N[alpha] and is divisible by
    u_mu = prod m_i(mu)!.  Our coefficients arrive through fraction
    arithmetic, so this is a check on the whole route.
  * The closed-form norms <J,J> = H H' — compared against Sage's own
    scalar_jack, which is a different computation, not the same formula.
  * Stanley's 1989 conjecture, <J_la J_mu, J_nu> in N[alpha].  OPEN.  A
    violation is a RESULT TO REPORT, never a bug to fix; check_deltaop.py
    treats the valley Delta conjecture the same way.
"""

import sys
from collections import defaultdict

from sage.all import QQ, Partition, Partitions, SymmetricFunctions, factorial, prod

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's Jack bases")

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/jack.txt"

R = QQ["t"]
ALPHA = R.gen()  # Sage's Jack parameter, our alpha
F = R.fraction_field()
Sym = SymmetricFunctions(F)
jack = Sym.jack()
P, Q, J = jack.P(), jack.Q(), jack.J()
m, p = Sym.monomial(), Sym.powersum()

fail = 0
checked = 0


def bad(msg):
    global fail
    fail += 1
    if fail <= 25:
        print(f"  FAIL {msg}")


def shape(s):
    return tuple(int(x) for x in s.split(",") if x)


def value(num, den, scale):
    """Rebuild one AFrac cell: (sum c_k alpha^k) / (scale * prod (u a + v)^m)."""
    v = sum(F(int(c)) * ALPHA**k for k, c in enumerate(num.split(",")) if c)
    v /= F(int(scale))
    if den:
        for atom in den.split(";"):
            u, w, mult = (int(x) for x in atom.split(":"))
            v /= (u * ALPHA + w) ** mult
    return v


rows = defaultdict(dict)
with open(path) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        kind, lam, mu, num, den, scale = line.split("|")
        rows[(kind, lam)][shape(mu)] = value(num, den, scale)

print(f"read {sum(len(v) for v in rows.values())} coefficients "
      f"in {len(rows)} expansions from {path}")


def compare(label, got, want_elt):
    """`got` is dict mu -> value; `want_elt` a Sage symmetric function."""
    global checked
    want = {tuple(k): v for k, v in want_elt.monomial_coefficients().items()}
    if set(got) != set(want):
        bad(f"{label}: shapes differ, ours-only {set(got) - set(want)}, "
            f"Sage-only {set(want) - set(got)}")
        return
    for mu, v in got.items():
        checked += 1
        if v != want[mu]:
            bad(f"{label} at {mu}: ours {v}, Sage {want[mu]}")


# ----------------------------------------------------- 1. the three bases

for kind, basis, into in (("p", P, m), ("q", Q, m), ("j", J, m), ("jp", J, p)):
    n_rows = 0
    for (k, lam), got in rows.items():
        if k != kind:
            continue
        n_rows += 1
        la = Partition(list(shape(lam)))
        compare(f"{kind}_{la}", got, into(basis[la]))
    print(f"  {kind}: {n_rows} expansions vs Sage")


# ------------------------------------------------ 1b. the inverse direction
#
# `m_lambda` written in each normalization, against Sage's own conversion.
# The round trip in `cargo test` cannot see an error the forward direction
# shares; this can, because Sage solves the same triangular system from its
# own P.

for kind, basis in (("mp", P), ("mq", Q), ("mj", J)):
    n_rows = 0
    for (k, lam), got in rows.items():
        if k != kind:
            continue
        n_rows += 1
        la = Partition(list(shape(lam)))
        compare(f"{kind}_{la}", got, basis(m[la]))
    print(f"  {kind}: {n_rows} expansions vs Sage")


# -------------------------------------------- 2. norms, the closed form

n_norm = 0
for (k, lam), got in rows.items():
    if k != "normj":
        continue
    n_norm += 1
    la = Partition(list(shape(lam)))
    ours = got[()]
    theirs = J[la].scalar_jack(J[la])
    checked += 1
    if ours != theirs:
        bad(f"<J_{la},J_{la}>: ours {ours}, Sage {theirs}")
print(f"  normj: {n_norm} closed-form norms vs Sage's scalar_jack")


# ------------------------------------- 3. Stanley's structure constants

n_g = 0
neg = []
for (k, pair), got in rows.items():
    if k != "g":
        continue
    a, b = pair.split("+")
    la, mu = Partition(list(shape(a))), Partition(list(shape(b)))
    prodJ = J[la] * J[mu]
    for nu_t, ours in got.items():
        n_g += 1
        checked += 1
        theirs = prodJ.scalar_jack(J[Partition(list(nu_t))])
        if ours != theirs:
            bad(f"<J_{la} J_{mu}, J_{list(nu_t)}>: ours {ours}, Sage {theirs}")
        poly = R(ours)
        if any(c < 0 or c.denominator() != 1 for c in poly.coefficients()):
            neg.append((list(la), list(mu), list(nu_t), poly))
print(f"  g: {n_g} Stanley structure constants vs Sage")


# --------------------------------- 4. [KS] Thm 1.1, a law with no oracle

n_ks = 0
for (k, lam), got in rows.items():
    if k != "j":
        continue
    for mu_t, v in got.items():
        n_ks += 1
        checked += 1
        poly = R(v)  # raises if a denominator survived, and that raise is the check
        u = prod(factorial(e) for e in Partition(list(mu_t)).to_exp() if e)
        quo = poly / u
        if quo.denominator() != 1:
            bad(f"[KS] Thm 1.1: u_{list(mu_t)} = {u} does not divide "
                f"[m]J_{list(shape(lam))} = {poly}")
        elif any(c < 0 for c in poly.coefficients()):
            bad(f"[KS] Thm 1.1: [m_{list(mu_t)}]J_{list(shape(lam))} = {poly} "
                f"left N[alpha]")
print(f"  [KS] Thm 1.1: {n_ks} coefficients in N[alpha] and divisible by u_mu")


# --------------------------------------- 5. the two zonal normalizations

# Measured, not assumed: Sage's zonal() is P^(2).  Re-checked here so a change
# on either side is caught rather than silently retested against the wrong
# object.
Z = SymmetricFunctions(QQ).zonal()
mz = SymmetricFunctions(QQ).monomial()
n_z = 0
for (k, lam), got in rows.items():
    if k != "p":
        continue
    la = Partition(list(shape(lam)))
    if la.size() > 5:
        continue
    n_z += 1
    want = {tuple(kk): v for kk, v in mz(Z[la]).monomial_coefficients().items()}
    for mu, v in got.items():
        checked += 1
        if R(v.numerator())(2) / R(v.denominator())(2) != want.get(mu, 0):
            bad(f"zonal: P_{la}(alpha=2) at {mu} is not Sage's zonal")
print(f"  zonal: {n_z} shapes confirm Sage's zonal() = P^(2), not J^(2)")


# -------------------------------- 6. the scalar_jack zero-crash wart

# Sage cannot return a zero scalar_jack (jack.py calls .denominator() on the
# plain int 0).  Confirm the crash set is EXACTLY the vanishing set, so the
# oracle's silence above is silence about zeros and nothing else.
crashes = 0
for n in range(1, 5):
    for la in Partitions(n):
        for mu in Partitions(n):
            try:
                got = P(p[la]).scalar_jack(P(p[mu]))
            except TypeError:
                crashes += 1
                checked += 1
                if la == mu:
                    bad(f"scalar_jack crashed on a NONZERO pairing <p_{la},p_{la}>")
                continue
            checked += 1
            want = Partition(la).centralizer_size() * ALPHA ** len(la) if la == mu else 0
            if got != want:
                bad(f"<p_{la},p_{mu}>_alpha: Sage {got}, expected {want}")
print(f"  oracle wart: scalar_jack raised on {crashes} vanishing pairings, "
      f"none nonzero")


print()
print(f"compared {checked} values")
if neg:
    print()
    print("*** STANLEY POSITIVITY VIOLATED — THIS IS A RESULT, NOT A BUG ***")
    for la, mu, nu, poly in neg:
        print(f"    <J_{la} J_{mu}, J_{nu}> = {poly}")
    print("Check the normalization against docs/record/jack.md before")
    print("believing it; then report it.")
else:
    print("Stanley: every structure constant computed lies in N[alpha] "
          "(the conjecture is OPEN; this is evidence)")
print(f"FAILURES: {fail}")
sys.exit(1 if fail else 0)
