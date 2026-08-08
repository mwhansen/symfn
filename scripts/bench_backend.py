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
rest of this project has learned repeatedly: Sage memoizes, this machine drifts
as it warms, and whichever side runs first pays for a cold cache.
"""

import os
import subprocess
import sys
import time

ROUNDS = int(sys.argv[1]) if len(sys.argv) > 1 and sys.argv[1] != "--run" else 5


def child_env(mode):
    """The environment the `mode` arm has to start with.

    `SAGE_DISABLE_SYMFN` cannot be set from inside the child: Sage fills its
    conversion table when `sage.combinat.sf.classical` is imported, which is
    before any line of the harness runs.
    """
    env = dict(os.environ)
    if mode == "symmetrica":
        env["SAGE_DISABLE_SYMFN"] = "1"
    else:
        env.pop("SAGE_DISABLE_SYMFN", None)
    return env


def run(mode):
    from sage.all import QQ, Partitions, SymmetricFunctions
    from sage.libs.symfn import is_available

    # Checked rather than installed -- see the same block in
    # scripts/check_backend.py for why an unverified control is worse than none.
    want = mode == "symfn"
    if is_available() != want:
        raise SystemExit(
            f"the {mode} arm wanted is_available() == {want}; Sage says "
            f"{is_available()}."
        )

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
            # The six pairs among h, e and p, which reach each other without
            # the Schur hub. Going through it inflates: p_lambda has a
            # handful of terms and its Schur expansion has p(n) of them.
            timed(f"[{label}] p -> h  deg {deg}", lambda: sweep(p, h))
            timed(f"[{label}] p -> e  deg {deg}", lambda: sweep(p, e))
            timed(f"[{label}] h -> p  deg {deg}", lambda: sweep(h, p))
            timed(f"[{label}] e -> p  deg {deg}", lambda: sweep(e, p))
            timed(f"[{label}] h -> e  deg {deg}", lambda: sweep(h, e))
            timed(f"[{label}] e -> h  deg {deg}", lambda: sweep(e, h))

    # The character bases reach a conversion by *peeling* -- one small
    # conversion per term removed -- and each peel step is an h -> p and a
    # p -> h. Nothing else in this file makes thousands of small conversions,
    # and that is the shape a missing direct route punishes hardest.
    Sym = SymmetricFunctions(QQ)
    h, ht = Sym.homogeneous(), Sym.ht()
    for shape in ([4, 3], [5, 3]):
        timed(
            f"[QQ] h -> ht  {shape}",
            lambda shape=shape: len(ht(h[shape] * h[shape]).monomial_coefficients()),
        )

    # And the same conversion *again*, at a degree already visited. Sage's peel
    # caches its own expansions, but the conversions underneath it are what
    # repeat; whichever backend memoizes those answers the second call for
    # nothing. One-call-per-process timings hide this entirely, which is the
    # cold case and not the one a session spends its time in.
    repeats = [[6, 2], [4, 4], [7, 1]]
    timed(
        "[QQ] h -> ht  again",
        lambda: sum(
            len(ht(h[a] * h[5, 3]).monomial_coefficients()) for a in repeats
        ),
    )

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
                env=child_env(mode),
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
