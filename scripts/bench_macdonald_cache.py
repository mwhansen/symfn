"""Time Sage's Macdonald `J` change of basis with symfn against without it.

    python scripts/bench_macdonald_cache.py [top_degree] [rounds]

`J._s_cache(n)` fills both directions of the `J` <-> `s` transition for a whole
degree, and it is the only thing anything Macdonald-shaped in Sage waits for --
one call per degree, then every conversion at that degree is a table lookup. So
the degree is the unit of work here, as it is for the backend itself.

The two arms run in **separate processes**, alternating, because both Sage and
symfn memoise and whichever ran first would otherwise be charged for the cold
cache. The control arm needs `SAGE_DISABLE_SYMFN` in its *environment*: Sage
fills its conversion table at import, before any line here runs, and
`Feature.is_present` caches, so setting it from inside the child is too late and
would silently compare symfn to symfn.

⚠️ Record the power state with the numbers. This machine drifts about 1.8x on
battery, and it changes ratios and not only times.
"""

import os
import subprocess
import sys
import time

TOP = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "--run" else 10
ROUNDS = int(sys.argv[2]) if len(sys.argv) > 2 and sys.argv[1] != "--run" else 3


def child_env(mode):
    env = dict(os.environ)
    if mode == "symmetrica":
        env["SAGE_DISABLE_SYMFN"] = "1"
    else:
        env.pop("SAGE_DISABLE_SYMFN", None)
    return env


def run(mode, top):
    from sage.all import QQ, PolynomialRing, SymmetricFunctions

    from sage.libs.symfn import is_available
    if (mode == "symmetrica") == is_available():
        print(f"the {mode} arm has is_available() = {is_available()}", file=sys.stderr)
        sys.exit(1)

    QQqt = PolynomialRing(QQ, "q,t").fraction_field()
    J = SymmetricFunctions(QQqt).macdonald().J()
    for n in range(1, top + 1):
        t = time.perf_counter()
        J._s_cache(n)
        elapsed = time.perf_counter() - t
        # The cell count is the comparability check: the two arms must fill the
        # same table, and a differing count means they did not.
        cells = sum(len(row) for row in J._s_to_self_cache[n].values())
        print(f"{n}\t{elapsed:.6f}\t{cells}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--run":
        run(sys.argv[2], int(sys.argv[3]))
        sys.exit(0)

    totals = {"symmetrica": {}, "symfn": {}}
    cells = {}
    for r in range(ROUNDS):
        order = ("symmetrica", "symfn") if r % 2 == 0 else ("symfn", "symmetrica")
        for mode in order:
            proc = subprocess.run(
                [sys.executable, __file__, "--run", mode, str(TOP)],
                capture_output=True,
                text=True,
                env=child_env(mode),
            )
            if proc.returncode != 0:
                print(f"{mode} failed:\n{proc.stderr[-3000:]}")
                sys.exit(1)
            for line in proc.stdout.splitlines():
                n, secs, count = line.split("\t")
                totals[mode][int(n)] = totals[mode].get(int(n), 0.0) + float(secs)
                cells.setdefault(int(n), set()).add(count)

    bad = [n for n, v in cells.items() if len(v) != 1]
    if bad:
        print(f"table sizes differ, not comparable: {bad}")
        sys.exit(1)

    print(f"{ROUNDS} rounds, alternating, separate processes\n")
    print(f"{'n':>3} {'cells':>7} {'symmetrica':>11} {'symfn':>11} {'ratio':>8}")
    for n in sorted(totals["symmetrica"]):
        a, b = totals["symmetrica"][n] / ROUNDS, totals["symfn"][n] / ROUNDS
        print(f"{n:>3} {cells[n].pop():>7} {a:>10.4f}s {b:>10.4f}s {a / b:>7.2f}x")
