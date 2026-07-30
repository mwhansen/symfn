"""The Sage side of the Jack ladder, one fresh process per measurement.

    sage -python -u scripts/bench_jack.py 11          # the sweep
    sage -python -u scripts/bench_jack.py --one pm 9  # one cell, used internally

`docs/record/jack-spec.md` §2 measured these walls before any Rust existed, but
that sweep straddled a power-state change (AC detached mid-run). This script
exists so both sides can be re-measured in ONE session and ONE power state — the
rule `docs/record/macdonald-operators-spec.md` arrived at the hard way, after a
1.8x battery drift on this machine inflated every ratio in an earlier table.

⚠️ THE ISOLATION IS NOT OPTIONAL, and finding that out is what this docstring
is for.  A first version ran every unit in one process and produced nonsense
in two different ways:

  * Sage MEMOIZES its Jack transition matrices, so a `J -> m` sweep run after a
    `P -> m` sweep reads a warm cache.  Measured 0.152 s for the whole of
    degree 10 against 44.1 s for the same work cold — a 290x difference that is
    entirely cache reuse and says nothing about either operation.
  * SIGALRM fires INSIDE Sage's cache construction and leaves it half-built.
    The very next degree then died with `KeyError: [11]` out of
    `sfa.py:_from_cache`.  So a timeout does not merely abandon an item; it
    poisons the process for everything after it.

Hence: one subprocess per (unit, degree), killed by the parent on timeout.
Slower to run, and the only version whose numbers mean anything.

Sage is a black-box oracle here: run as a separate program, timed from the
outside, source not read (see ../NOTICE.md).

The units are the ones a user actually pays, and the ones
`examples/bench_jack.rs` measures on the Rust side:

  pm     whole-degree P -> m   (every lambda |- n expanded)
  jm     whole-degree J -> m
  jp     whole-degree J -> p   (the [GJ] unit)
  norm   the norms table <J_la, J_la> for every lambda |- n
"""

import subprocess
import sys
import time

UNITS = ("pm", "jm", "jp", "norm")


def run_one(unit, n):
    """The measurement itself, in a process of its own."""
    from sage.all import QQ, Partitions, PolynomialRing, SymmetricFunctions

    R = PolynomialRing(QQ, "t")
    Sym = SymmetricFunctions(R.fraction_field())
    jack = Sym.jack()
    P, J = jack.P(), jack.J()
    m, p = Sym.monomial(), Sym.powersum()
    work = {
        "pm": lambda la: m(P[la]),
        "jm": lambda la: m(J[la]),
        "jp": lambda la: p(J[la]),
        "norm": lambda la: J[la].scalar_jack(J[la]),
    }[unit]

    total, worst, worst_shape = 0.0, 0.0, None
    for la in Partitions(n):
        t0 = time.time()
        work(la)
        secs = time.time() - t0
        total += secs
        if secs > worst:
            worst, worst_shape = secs, list(la)
    print(f"RESULT {total:.3f} {worst:.3f} {worst_shape}")


def sweep(top, limit):
    print(f"SageMath side, one fresh process per cell, per-cell timeout {limit}s")
    print("⚠️ record the power state: this machine drifts ~1.8x on battery")
    print()
    for unit in UNITS:
        print(f"== {unit} ==")
        print(f"{'n':>3} {'total(s)':>10} {'max(s)':>10}  worst")
        for n in range(1, top + 1):
            t0 = time.time()
            try:
                out = subprocess.run(
                    [sys.executable, "-u", __file__, "--one", unit, str(n)],
                    capture_output=True,
                    text=True,
                    timeout=limit,
                )
            except subprocess.TimeoutExpired:
                print(f"{n:>3} {'>' + str(limit):>10} {'—':>10}  killed")
                sys.stdout.flush()
                break
            line = next(
                (x for x in out.stdout.splitlines() if x.startswith("RESULT")), None
            )
            if line is None:
                print(f"{n:>3} {'ERR':>10} {'—':>10}  {out.stderr.strip()[-90:]}")
                sys.stdout.flush()
                break
            _, total, worst, shape = line.split(maxsplit=3)
            # Also report the wall time including Sage startup, so the overhead
            # this isolation costs is visible rather than hidden.
            print(
                f"{n:>3} {float(total):>10.3f} {float(worst):>10.3f}  {shape}"
                f"   (+{time.time() - t0 - float(total):.1f}s startup)"
            )
            sys.stdout.flush()
        print()


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--one":
        run_one(sys.argv[2], int(sys.argv[3]))
    else:
        sweep(
            int(sys.argv[1]) if len(sys.argv) > 1 else 11,
            int(sys.argv[2]) if len(sys.argv) > 2 else 180,
        )
