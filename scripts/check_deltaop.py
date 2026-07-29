"""Hold `deltaop` to Sage, and to the identities Sage cannot check.

    cargo run --release --example deltaop_dump -- 7 > /tmp/dop.txt
    sage -python scripts/check_deltaop.py /tmp/dop.txt

Sage is used as a black box: run as a separate program, output compared, source
not read (see ../NOTICE.md).

**Sage can only check one of the five operators.** It has `nabla` and nothing
else -- no Δ, no Δ', no Θ, no Π (measured; see
`../docs/spec-macdonald-operators.md` §2.1). So the script does two different
things:

  * `nabla` and `nabla2` rows are compared against Sage's own `nabla`, term for
    term. That is a genuine external oracle.
  * `deltaprime` and `theta` rows have no oracle anywhere, and are checked
    against published *identities* instead -- Θ_{e_k} ∇ e_{n-k} = Δ'_{e_{n-k-1}}
    e_n, Δ'_{e_{n-1}} e_n = ∇ e_n, and the rise version of the Delta conjecture
    against a direct enumeration of labelled Dyck paths. The first two tie the
    unoracled operators to the oracled one.

Coefficients are compared as exact polynomials in ℤ[q,t].
"""

import sys
from collections import defaultdict

from sage.all import (SymmetricFunctions, QQ, PolynomialRing, FractionField,
                      Partition, Permutations, Subsets, prod)

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/dop.txt"

R = FractionField(PolynomialRing(QQ, ['q', 't']))
q, t = R.gens()
Sym = SymmetricFunctions(R)
s, e, h = Sym.s(), Sym.e(), Sym.h()

rows = defaultdict(dict)
with open(path) as fh:
    for line in fh:
        line = line.strip()
        if not line:
            continue
        op, arg, lam, terms = (x.strip() for x in line.split("|"))
        poly = R.zero()
        for term in terms.split(";"):
            ab, c = term.split(":")
            a, b = ab.split(",")
            poly += R(int(c)) * q**int(a) * t**int(b)
        rows[(op, arg)][Partition([int(x) for x in lam.split(",") if x])] = poly

print(f"read {sum(len(v) for v in rows.values())} coefficients "
      f"in {len(rows)} expansions from {path}")

fail = 0


def compare(label, got, want):
    """`got` is a dict λ -> poly; `want` is a Sage symmetric function."""
    global fail
    want_s = s(want)
    keys = set(got) | {lam for lam, _ in want_s}
    for lam in keys:
        a = got.get(lam, R.zero())
        b = want_s.coefficient(lam)
        if a != b:
            print(f"  MISMATCH {label} at {lam}:\n    ours {a}\n    want {b}")
            fail += 1
            return
    print(f"  ok  {label}")


# ---------------------------------------------------------------- 1. vs Sage
print()
print("=== nabla, against Sage's own (the only external oracle available) ===")
for (op, arg), got in sorted(rows.items()):
    if op not in ("nabla", "nabla2"):
        continue
    if arg.startswith("e"):
        f = e[int(arg[1:])]
    elif arg.startswith("s"):
        f = s[[int(x) for x in arg[1:].split(",")]]
    else:
        continue
    want = f.nabla()
    if op == "nabla2":
        want = want.nabla()
    compare(f"{op} {arg}", got, want)

# --------------------------------------------------- 2. identities, no oracle
print()
print("=== Delta' and Theta: no oracle exists, so identities instead ===")

for (op, arg), got in sorted(rows.items()):
    if op != "theta":
        continue
    # arg is "e{k},nabla_e{m}" and the theorem says this equals
    # Delta'_{e_{m-1}} e_{k+m}.
    k = int(arg.split(",")[0][1:])
    m = int(arg.split("nabla_e")[1])
    n = k + m
    other = rows.get(("deltaprime", f"e{m - 1},e{n}"))
    if other is None:
        continue
    if got == other:
        print(f"  ok  Theta_e{k} nabla e_{m} = Delta'_e{m - 1} e_{n}")
    else:
        print(f"  MISMATCH Theta_e{k} nabla e_{m} vs Delta'_e{m - 1} e_{n}")
        fail += 1

for (op, arg), got in sorted(rows.items()):
    if op != "deltaprime":
        continue
    k, n = (int(x[1:]) for x in arg.split(","))
    if k == n - 1:
        nab = rows.get(("nabla", f"e{n}"))
        if nab is not None:
            if got == nab:
                print(f"  ok  Delta'_e{n - 1} e_{n} = nabla e_{n}")
            else:
                print(f"  MISMATCH Delta'_e{n - 1} e_{n} vs nabla e_{n}")
                fail += 1


# ------------------------------ 3. the Dyck enumeration, against parking functions
print()
print("=== labelled Dyck paths, against Sage's ParkingFunctions ===")
print("    (sum q^dinv t^area over distinct labellings; Sage's dinv/area were")
print("     confirmed to be HRW's by the shuffle theorem before being used here)")
from sage.all import ParkingFunctions

for (op, arg), got in sorted(rows.items()):
    if op != "dyck":
        continue
    n = int(arg[2:])
    ours = got[Partition([])]
    want = sum(q**p.dinv() * t**p.area() for p in ParkingFunctions(n))
    if ours == want:
        print(f"  ok  n={n}  ({ParkingFunctions(n).cardinality()} parking functions)")
    else:
        print(f"  MISMATCH n={n}\n    ours {ours}\n    want {want}")
        fail += 1


# ------------------------------------- 4. the Delta conjecture, rise version
def area_sequences(n):
    def rec(seq):
        if len(seq) == n:
            yield tuple(seq)
            return
        for a in range(1 if not seq else seq[-1] + 2):
            yield from rec(seq + [a])
    yield from rec([])


def labellings(a):
    n = len(a)
    for wd in Permutations(n):
        if all(a[i] != a[i - 1] + 1 or wd[i] > wd[i - 1] for i in range(1, n)):
            yield wd


def d_row(a, wd, i):
    return sum(1 for j in range(i + 1, len(a))
               if (a[i] == a[j] and wd[i] < wd[j])
               or (a[i] == a[j] + 1 and wd[i] > wd[j]))


def rise_side(n, k):
    """<Rise_{n,k}, h_1^n>, Haglund-Remmel-Wilson arXiv:1509.07058."""
    want = n - k - 1
    total = R.zero()
    for a in area_sequences(n):
        rs = [i for i in range(1, n) if a[i] == a[i - 1] + 1]
        if len(rs) < want:
            continue
        for wd in labellings(a):
            base = q**sum(d_row(a, wd, i) for i in range(n)) * t**sum(a)
            for S in Subsets(rs, want):
                total += base * prod([t**(-a[i]) for i in S], R.one())
    return total


print()
print("=== the Delta conjecture, rise version, against labelled Dyck paths ===")
for (op, arg), got in sorted(rows.items()):
    if op != "deltaprime":
        continue
    k, n = (int(x[1:]) for x in arg.split(","))
    if n > 6:
        continue                       # the enumeration is the slow side here
    f = sum((c * s(lam) for lam, c in got.items()), s.zero())
    lhs = f.scalar(h[[1] * n])
    rhs = rise_side(n, k)
    if lhs == rhs:
        print(f"  ok  <Delta'_e{k} e_{n}, h_1^{n}> = Rise_{{{n},{k}}}")
    else:
        print(f"  MISMATCH Delta'_e{k} e_{n} vs the rise side")
        fail += 1

print()
if fail:
    print(f"{fail} MISMATCH(ES)")
    raise SystemExit(1)
print("all checks passed")
