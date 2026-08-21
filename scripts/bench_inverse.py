"""Time the seven inverse expansions against Sage: the four Macdonald ones
`s -> Htilde`, `s -> J`, `m -> P`, `m -> Q`, and the three Jack ones
`m -> P`, `m -> Q`, `m -> J`.

    python scripts/bench_inverse.py [top_degree]

The unit of work is the **degree**: every partition of `n`, expanded out of the
Schur basis into the parametric one. That is what a caller asking whether some
family of functions is `Htilde`- or `J`-positive does, and it is the unit both
sides amortize their transition matrix over.

Each degree and each arm runs in its **own process**, and the arm part of that
is not optional: Sage caches the transition matrices behind these conversions,
the cache is not confined to one degree, and it is *shared between the
normalizations of one family*. Two Jack arms in one process read 20x faster
than they do alone, because the first one built the matrix the other two then
reused (`docs/record/jack.md`). symfn memoizes too; `clear_caches` in the
example drops it between workloads.

The Sage arm needs `SAGE_DISABLE_SYMFN` in its *environment*, and this script
refuses to run without it taking effect: Sage's Macdonald bases reach this
library through the optional backend, so an arm with it enabled compares symfn
to symfn and reports about 1.0x.

The two families print as two tables rather than one wide row, and they run
over different base rings: the Macdonald arms over the fraction field of
`QQ[q,t]`, the Jack ones over `QQ[alpha]`'s. Building the ring is outside the
timed region on both sides.

Sage is a black box here: run and timed, never read (see NOTICE.md).

WARNING: Record the power state with the numbers. This machine drifts about
1.8x on battery, and it changes ratios and not only times.
"""

import os
import subprocess
import sys
import time


def one(n, arm):
    from sage.all import Partitions, PolynomialRing, QQ, SymmetricFunctions

    from sage.libs.symfn import is_available

    if is_available():
        print(
            "the symfn backend is enabled; run with SAGE_DISABLE_SYMFN=1 in the "
            "environment or this compares symfn to symfn",
            file=sys.stderr,
        )
        sys.exit(1)

    if arm < 4:
        R = PolynomialRing(QQ, "q,t").fraction_field()
        Sym = SymmetricFunctions(R)
        mac = Sym.macdonald()
        basis = (mac.Ht(), mac.J(), mac.P(), mac.Q())[arm]
        source = Sym.schur() if arm < 2 else Sym.monomial()
    else:
        A = PolynomialRing(QQ, "a").fraction_field()
        Sym = SymmetricFunctions(A)
        jack = Sym.jack(t=A.gen())
        basis = (jack.P(), jack.Q(), jack.J())[arm - 4]
        source = Sym.monomial()

    # One degree-1 conversion first, untimed: Sage builds a family's coercion
    # machinery on its first use, and without this the arm that ran first paid
    # for the whole family. It touches no transition matrix of the degree
    # being measured.
    basis(source(Partitions(1)[0]))

    shapes = list(Partitions(n))
    start = time.perf_counter()
    for la in shapes:
        basis(source(la))
    print(len(shapes), f"{time.perf_counter() - start:.4f}")


if len(sys.argv) > 3 and sys.argv[1] == "--one":
    one(int(sys.argv[2]), int(sys.argv[3]))
    sys.exit(0)

TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 6
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

env = dict(os.environ)
env["SAGE_DISABLE_SYMFN"] = "1"

#: The Sage arm's seven timings, in order, against the example's workload
#: names, split into the two tables they print as.
MACDONALD = (("s->Ht", "ht"), ("s->J", "j"), ("m->P", "p"), ("m->Q", "q"))
JACK = (("m->P", "jack_p"), ("m->Q", "jack_q"), ("m->J", "jack_j"))
ARMS = MACDONALD + JACK


def ratio(theirs, ours):
    return f"{theirs / ours:>6.1f}x" if ours > 0 else "      --"


def header(arms):
    out = f"{'n':>3} {'p(n)':>5}"
    for label, _ in arms:
        out += f" {label + ' sage':>12} {'symfn':>9} {'ratio':>7}  "
    return out


rows = []
for n in range(1, TOP + 1):
    sage = []
    for arm in range(len(ARMS)):
        theirs = subprocess.run(
            [sys.executable, os.path.abspath(__file__), "--one", str(n), str(arm)],
            capture_output=True,
            text=True,
            env=env,
            cwd=ROOT,
        )
        if theirs.returncode:
            sys.stderr.write(theirs.stderr)
            sys.exit(1)
        count, secs = theirs.stdout.split()
        sage.append(float(secs))

    cargo = ["cargo", "run", "--release", "--quiet", "--example", "bench_inverse"]
    mine = subprocess.run(
        [*cargo, "--", str(n)], capture_output=True, text=True, cwd=ROOT
    ).stdout.split()
    ours = {mine[i]: float(mine[i + 1]) for i in range(0, len(mine), 3)}

    rows.append((n, count, sage, ours))

for title, arms in (("Macdonald", MACDONALD), ("Jack", JACK)):
    print(f"\n{title}")
    print(header(arms))
    for n, count, sage, ours in rows:
        row = f"{n:>3} {count:>5}"
        for label, name in arms:
            theirs_secs = sage[[a[1] for a in ARMS].index(name)]
            row += (
                f" {theirs_secs:>11.4f}s {ours[name]:>8.4f}s "
                f"{ratio(theirs_secs, ours[name]):>7}  "
            )
        print(row)
