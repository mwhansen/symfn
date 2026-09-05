"""Time the whole (q,t)-Kostka table of each degree, cold on both sides.

    sage -python scripts/bench_qt_kostka.py 9

The unit is the **table**, not the pair. Sage caches the transition matrices
behind `qt_kostka`, so a second ask for any pair of a degree it has already
touched is a dictionary lookup; timing pair-by-pair would compare symfn's
arithmetic against Sage's `dict.__getitem__`. symfn caches too (partitions,
characters), which `clear_caches` drops.

Worse, the caching is not even confined to one degree: running degrees 5, 6, 7
in one Sage process gives 0.106s, 0.056s, 0.144s, where cold they are 0.107s,
0.150s, 0.250s -- degree 6 comes out *faster* than degree 5 because degree 5 paid
for machinery both share. So each degree is timed in its own process, which is
what the `--one` mode is for; the top-level run just drives them.

Sage is a black box here: run and timed, never read (see NOTICE.md).
"""

import subprocess
import sys
import time


def one(n):
    sys.path.insert(0, "pybuild")
    from sage.all import Partitions
    from sage.combinat.sf.macdonald import qt_kostka

    from sage_guard import require_own_sage

    require_own_sage("the Sage arm")

    import symfn

    shapes = [list(p) for p in Partitions(n)]

    start = time.perf_counter()
    for mu in shapes:
        for lam in shapes:
            qt_kostka(lam, mu)
    theirs = time.perf_counter() - start

    symfn.clear_caches()
    start = time.perf_counter()
    symfn.qt_kostka_table(n)
    ours = time.perf_counter() - start

    print(f"{n} {len(shapes) ** 2} {ours:.4f} {theirs:.4f}")


if len(sys.argv) > 2 and sys.argv[1] == "--one":
    one(int(sys.argv[2]))
    sys.exit(0)

TOP = int(sys.argv[1]) if len(sys.argv) > 1 else 9
print(f"{'n':>3} {'values':>8} {'symfn':>10} {'sage':>10} {'ratio':>8}")
for n in range(1, TOP + 1):
    out = subprocess.run(
        ["sage", "-python", __file__, "--one", str(n)],
        capture_output=True,
        text=True,
    ).stdout.split()
    _, count, ours, theirs = out[-4:]
    ratio = float(theirs) / float(ours) if float(ours) > 0 else float("inf")
    print(f"{n:>3} {count:>8} {ours:>10} {theirs:>10} {ratio:>7.1f}x")
