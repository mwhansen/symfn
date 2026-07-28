"""Time Sage with symfn substituted for Symmetrica, against Symmetrica itself.

    sage -python scripts/bench_backend.py [rounds]

`scripts/check_backend.py` establishes that the substitution is *correct*. This
asks whether it is worth making, and the honest answer has to include the shim:
symfn's Rust core being faster does not by itself make Sage faster, because the
shim pays Python costs Symmetrica's Cython path does not -- building `Partition`
objects, a dict per call, and marshalling every coefficient across the FFI
boundary twice.

So what is timed here is Sage-level operations end to end, with the backend as
the only difference. Both regimes appear, because they are different code paths:
over QQ the whole element crosses in one call, and over QQ['t'] Sage calls the
backend once per partition and recombines, which multiplies the per-call cost by
the size of the support and is the harder case for a shim.

The two backends run in **separate processes**, alternating, for the reason the
rest of this project has learned repeatedly: Sage memoises, this machine drifts
as it warms, and whichever side runs first pays for a cold cache.
"""

import subprocess
import sys
import time

ROUNDS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "--run" else 5


def run(mode):
    sys.path.insert(0, "scripts")
    from sage.all import QQ, Partitions, SymmetricFunctions
    from sage.combinat.sf import classical

    import sage_backend

    if mode == "symfn":
        sage_backend.install()
    else:
        classical.init()

    cases = []

    def timed(label, fn):
        t = time.perf_counter()
        n = fn()
        cases.append((label, time.perf_counter() - t, n))

    for ring, label in ((QQ, "QQ"), (QQ["t"], "QQ[t]")):
        Sym = SymmetricFunctions(ring)
        s, m, h, e, p = (
            Sym.schur(),
            Sym.monomial(),
            Sym.homogeneous(),
            Sym.elementary(),
            Sym.power(),
        )
        for deg in (10, 14, 18):
            parts = [list(x) for x in Partitions(deg)][:12]

            def sweep(src, dst, parts=parts):
                return sum(len(dst(src[L]).monomial_coefficients()) for L in parts)

            timed(f"[{label}] s -> m  deg {deg}", lambda: sweep(s, m))
            timed(f"[{label}] m -> s  deg {deg}", lambda: sweep(m, s))
            timed(f"[{label}] s -> p  deg {deg}", lambda: sweep(s, p))
            timed(f"[{label}] p -> s  deg {deg}", lambda: sweep(p, s))
            timed(f"[{label}] s -> h  deg {deg}", lambda: sweep(s, h))
            timed(f"[{label}] s -> e  deg {deg}", lambda: sweep(s, e))

    for label, elapsed, n in cases:
        print(f"{label}\t{elapsed:.6f}\t{n}")


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--run":
        run(sys.argv[2])
        sys.exit(0)

    totals = {"symmetrica": {}, "symfn": {}}
    counts = {}
    for r in range(ROUNDS):
        # Rotate which backend goes first.
        order = ("symmetrica", "symfn") if r % 2 == 0 else ("symfn", "symmetrica")
        for mode in order:
            proc = subprocess.run(
                [sys.executable, __file__, "--run", mode],
                capture_output=True,
                text=True,
            )
            if proc.returncode != 0:
                print(f"{mode} failed:\n{proc.stderr[-3000:]}")
                sys.exit(1)
            for line in proc.stdout.splitlines():
                label, secs, n = line.split("\t")
                totals[mode][label] = totals[mode].get(label, 0.0) + float(secs)
                counts.setdefault(label, set()).add(n)

    bad = [k for k, v in counts.items() if len(v) != 1]
    if bad:
        print(f"term counts differ, not comparable: {bad}")
        sys.exit(1)

    print(f"{ROUNDS} rounds, alternating, separate processes\n")
    print(f"{'case':<24} {'symmetrica':>11} {'symfn':>11} {'ratio':>8}")
    for label in totals["symmetrica"]:
        a = totals["symmetrica"][label]
        b = totals["symfn"][label]
        print(f"{label:<24} {a:>10.4f}s {b:>10.4f}s {a / b:>7.2f}x")
    ta = sum(totals["symmetrica"].values())
    tb = sum(totals["symfn"].values())
    print(f"\n{'total':<24} {ta:>10.4f}s {tb:>10.4f}s {ta / tb:>7.2f}x")
