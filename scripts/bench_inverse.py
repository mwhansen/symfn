"""Time the inverse expansions `s -> Htilde`, `s -> J`, `m -> P` and `m -> Q`
against Sage.

    python scripts/bench_inverse.py [top_degree]

The unit of work is the **degree**: every partition of `n`, expanded out of the
Schur basis into the parametric one. That is what a caller asking whether some
family of functions is `Htilde`- or `J`-positive does, and it is the unit both
sides amortize their transition matrix over.

Each degree and each arm runs in its **own process**. Sage caches the
transition matrices behind these conversions and the cache is not confined to
one degree, so a sequential run charges the first degree for machinery all of
them share (`docs/record/qt-kostka.md` measures this). symfn memoizes too;
`clear_caches` in the example drops it between workloads.

The Sage arm needs `SAGE_DISABLE_SYMFN` in its *environment*, and this script
refuses to run without it taking effect: Sage's Macdonald bases reach this
library through the optional backend, so an arm with it enabled compares symfn
to symfn and reports about 1.0x.

Sage is a black box here: run and timed, never read (see NOTICE.md).

WARNING: Record the power state with the numbers. This machine drifts about
1.8x on battery, and it changes ratios and not only times.
"""

import os
import subprocess
import sys
import time


def one(n):
    from sage.all import Partitions, PolynomialRing, QQ, SymmetricFunctions

    from sage.libs.symfn import is_available

    if is_available():
        print(
            "the symfn backend is enabled; run with SAGE_DISABLE_SYMFN=1 in the "
            "environment or this compares symfn to symfn",
            file=sys.stderr,
        )
        sys.exit(1)

    R = PolynomialRing(QQ, "q,t").fraction_field()
    Sym = SymmetricFunctions(R)
    s = Sym.schur()
    m = Sym.monomial()
    Ht = Sym.macdonald().Ht()
    J = Sym.macdonald().J()
    P = Sym.macdonald().P()
    Q = Sym.macdonald().Q()
    shapes = list(Partitions(n))

    times = []
    for basis, source in ((Ht, s), (J, s), (P, m), (Q, m)):
        start = time.perf_counter()
        for la in shapes:
            basis(source(la))
        times.append(time.perf_counter() - start)

    print(len(shapes), " ".join(f"{v:.4f}" for v in times))


if len(sys.argv) > 2 and sys.argv[1] == "--one":
    one(int(sys.argv[2]))
    sys.exit(0)

TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 6
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

env = dict(os.environ)
env["SAGE_DISABLE_SYMFN"] = "1"

#: The Sage arm's four timings, in order, against the example's workload names.
ARMS = (("s->Ht", "ht"), ("s->J", "j"), ("m->P", "p"), ("m->Q", "q"))


def ratio(theirs, ours):
    return f"{theirs / ours:>6.1f}x" if ours > 0 else "      --"


header = f"{'n':>3} {'p(n)':>5}"
for label, _ in ARMS:
    header += f" {label + ' sage':>12} {'symfn':>9} {'ratio':>7}  "
print(header)

for n in range(1, TOP + 1):
    theirs = subprocess.run(
        [sys.executable, os.path.abspath(__file__), "--one", str(n)],
        capture_output=True,
        text=True,
        env=env,
        cwd=ROOT,
    )
    if theirs.returncode:
        sys.stderr.write(theirs.stderr)
        sys.exit(1)
    count, *sage = theirs.stdout.split()
    sage = [float(v) for v in sage]

    cargo = ["cargo", "run", "--release", "--quiet", "--example", "bench_inverse"]
    mine = subprocess.run(
        [*cargo, "--", str(n)], capture_output=True, text=True, cwd=ROOT
    ).stdout.split()
    ours = {mine[i]: float(mine[i + 1]) for i in range(0, len(mine), 3)}

    row = f"{n:>3} {count:>5}"
    for theirs_secs, (_, name) in zip(sage, ARMS):
        row += (
            f" {theirs_secs:>11.4f}s {ours[name]:>8.4f}s "
            f"{ratio(theirs_secs, ours[name]):>7}  "
        )
    print(row)
