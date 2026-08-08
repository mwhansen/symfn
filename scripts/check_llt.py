"""Hold `llt` to Sage: the four ribbon normalizations and the tuple model.

    cargo run --release --example llt_dump -- 5 > /tmp/llt.txt
    sage -python scripts/check_llt.py /tmp/llt.txt

Sage is used as a black-box oracle only — run as a separate program and its
output compared, never read for algorithms (see ../NOTICE.md). This subsystem
was built clean-room: `sage/combinat/llt.py` and `ribbon_tableau.py` were not
read and must not be.

The in-crate tests already hold every model against every other model, plus
`hl`, `qtkostka`, `deltaop` and `skew_lr`. Sage is the *third* opinion, and the
only one that can catch a convention the crate shares with itself — spin vs
cospin, q vs 1/q, a transposed index, or a permuted tuple.

Four oracle warts, all measured before relying on this:

  * Sage's LLT parameter is named **t**, not q, and `llt(k)` over a fraction
    field in a variable named q fails to coerce. Construct over QQ['t'] or pass
    `t=` explicitly.
  * `cospin` rejects plain int lists for shapes ('float' object has no attribute
    'floor'); wrap in `Partition`.
  * `cospin(tuple)` is **floored**: it returns `q^{-min inv} G_nu`, where our
    `llt_g` returns the raw inv grading. The dump ships the floor separately so
    the dictionary is checked rather than assumed.
  * A **global cache** survives fresh `SymmetricFunctions` parents, so timings
    here mean nothing; this script checks values only.

`hspin` at k = 1 is the Schur function and `cospin` at k = 1 likewise, which is
why the dump includes k = 1: it is the cheapest place a spin/cospin swap shows.
"""

import sys
from collections import defaultdict

from sage.all import Partition, PolynomialRing, QQ, SymmetricFunctions

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/llt.txt"

R = PolynomialRing(QQ, "t")
t = R.gen()
Sym = SymmetricFunctions(R.fraction_field())
m = Sym.monomial()

# One family object per level, built once: the constructor is where the
# coercion wart lives, and re-entering it per record would also re-enter the
# global cache and hide nothing useful.
FAMILY = {}


def family(k):
    if k not in FAMILY:
        FAMILY[k] = Sym.llt(k, t=t)
    return FAMILY[k]


def parse(field):
    return tuple(int(x) for x in field.split(",") if x != "")


def poly(field):
    out = R.zero()
    for term in field.split(","):
        e, c = term.split(":")
        out += R(int(c)) * t ** int(e)
    return out


# ------------------------------------------------------------------ the dump

# kind -> key -> {nu: poly}
mine = defaultdict(lambda: defaultdict(dict))
floors = {}
quotients = {}

for line in open(path):
    fields = [f.strip() for f in line.split("|")]
    kind = fields[0]
    if kind == "floor":
        floors[fields[1]] = int(fields[2])
    elif kind == "quot":
        quotients[(int(fields[1]), parse(fields[2]))] = [
            parse(s) for s in fields[3].split(";")
        ]
    elif kind == "tuple":
        mine[kind][fields[1]][parse(fields[2])] = poly(fields[3])
    else:
        key = (int(fields[1]), parse(fields[2]))
        mine[kind][key][parse(fields[3])] = poly(fields[4])

fail = 0
checked = 0


def compare(name, got, want):
    """`got` is a {nu: poly} dict of ours; `want` a Sage symmetric function."""
    global fail, checked
    checked += 1
    want_d = {
        tuple(int(x) for x in nu): c
        for nu, c in m(want).monomial_coefficients().items()
    }
    keys = set(got) | set(want_d)
    for nu in sorted(keys):
        a = got.get(nu, R.zero())
        b = want_d.get(nu, R.zero())
        if a != b:
            print(f"FAIL {name}  m_{list(nu)}: mine {a}, Sage {b}")
            fail += 1
            return


# -------------------------------------------------- the ribbon dictionaries

for (k, mu), got in sorted(mine["hspin"].items()):
    compare(f"hspin k={k} mu={list(mu)}", got, family(k).hspin()[Partition(list(mu))])

for (k, mu), got in sorted(mine["hcospin"].items()):
    compare(
        f"hcospin k={k} mu={list(mu)}", got, family(k).hcospin()[Partition(list(mu))]
    )

for (k, lam), got in sorted(mine["cospin"].items()):
    # `Partition`, not a list: the int-list path raises inside Sage.
    compare(
        f"cospin k={k} lambda={list(lam)}", got, family(k).cospin(Partition(list(lam)))
    )

# ------------------------------------------------------ the tuple dictionary
#
# Sage floors, we do not.  Checking `sage == mine / t^floor` tests the floor and
# the model in one line; checking `sage == mine` would pass on every tuple whose
# floor happens to be zero and fail on the interesting ones.

for key, got in sorted(mine["tuple"].items()):
    shapes = [list(parse(s)) for s in key.split(";")]
    k = len(shapes)
    floor = floors[key]
    scaled = {}
    for nu, c in got.items():
        num = c / t**floor
        assert num.denominator() == 1, (key, nu, c, floor)
        scaled[nu] = R(num)
    compare(f"tuple {shapes} (floor {floor})", scaled, family(k).cospin(shapes))

# ---------------------------------------------- the quotient, against Sage's
#
# Not the dictionary — that is checked in-crate — but the component *order*,
# which nothing else here would catch and which `G_nu` is not symmetric in.

for (k, lam), got in sorted(quotients.items()):
    want = [tuple(int(x) for x in p) for p in Partition(list(lam)).quotient(k)]
    checked += 1
    if [tuple(g) for g in got] != want:
        print(f"FAIL quotient k={k} lambda={list(lam)}: mine {got}, Sage {want}")
        fail += 1

# ------------------------------------------------------------------- verdict

print(f"{checked} comparisons, {fail} failures")
if fail:
    sys.exit(1)
print("llt agrees with Sage")
