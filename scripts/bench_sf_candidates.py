"""Time the `sage.combinat.sf` paths that do not reach symfn, against symfn.

    python scripts/bench_sf_candidates.py            # every case
    python scripts/bench_sf_candidates.py llt_spin   # one case

Run with the sage-dev environment's `python`; this Sage's CLI has no
`sage -python`.

Every row times one Sage-level call both ways and includes everything the
caller pays, marshalling included -- which is the cost this backend keeps
finding once the algorithm is gone. Each arm is a child process, and the Sage
arm's has `SAGE_DISABLE_SYMFN=1` in its environment, so the two arms differ in
whether the call dispatches. The rows are `qt_kostka` for a whole degree,
`kfpoly` over every pair of that degree, `nabla` of every Schur function,
`reduced_kronecker` for the reduced Kronecker square of their sum, and
`llt_spin` and `llt_cospin` for `_m_cache(n)`, both directions of the level-`k`
basis at degree `n`.

The `qt_kostka` and `kfpoly` rows read as kernel-only figures until
2026-09-12, and that reading overstated what routing them was worth by two
orders of magnitude; `docs/record/python-and-sage-interop.md` records what the
difference turned out to be.

Every timing is one cold run in a fresh process: Sage memoizes changes of basis
at module level, and `symfn.clear_caches()` clears the kernel's own. The
parameter ring and the first dispatch are warmed before the clock starts, so
the row is the work and not the import.

⚠️ Record the power state (`pmset -g batt`). This machine drifts about 1.8x on
battery.
"""

import os
import subprocess
import sys
import time

# (case, degree, level) -- the level is used by the LLT rows only.
CASES = [
    ("qt_kostka", 7, None), ("qt_kostka", 8, None), ("qt_kostka", 9, None),
    ("kfpoly", 8, None), ("kfpoly", 10, None),
    ("nabla", 7, None), ("nabla", 8, None),
    ("reduced_kronecker", 4, None), ("reduced_kronecker", 5, None),
    ("reduced_kronecker", 6, None),
    ("llt_spin", 6, 3), ("llt_spin", 7, 3), ("llt_spin", 8, 3),
    ("llt_spin", 9, 2),
    ("llt_cospin", 7, 3), ("llt_cospin", 8, 3),
]


def child_env(case, arm):
    """`SAGE_DISABLE_SYMFN` is read when Sage imports, so it is set here.

    Turning the backend off also moves the classical conversions underneath
    the call back onto Symmetrica and Sage's own Python, which is what the
    Sage arm is: the whole route as it stands without this library.
    """
    env = dict(os.environ)
    if arm == "sage":
        env["SAGE_DISABLE_SYMFN"] = "1"
    else:
        env.pop("SAGE_DISABLE_SYMFN", None)
    return env


def run(case, n, k, arm):
    import symfn
    from sage.all import QQ, Partitions, SymmetricFunctions

    from sage.libs.symfn import is_available

    symfn.clear_caches()
    # A control arm that silently reached symfn would print a ratio near 1.0x
    # and prove nothing.
    assert is_available() == (arm == "symfn"), (case, arm)
    if case in ("llt_spin", "llt_cospin"):
        L = SymmetricFunctions(QQ["t"].fraction_field()).llt(k)
        B = L.hspin() if case == "llt_spin" else L.hcospin()
        start = time.perf_counter()
        B._m_cache(n)
    elif case == "nabla":
        s = SymmetricFunctions(QQ["q,t"].fraction_field()).s()
        f = sum(s(la) for la in Partitions(n))
        start = time.perf_counter()
        f.nabla()
    elif case == "reduced_kronecker":
        s = SymmetricFunctions(QQ).schur()
        f = sum(s(la) for la in Partitions(n))
        start = time.perf_counter()
        f.reduced_kronecker_product(f)
    elif case == "qt_kostka":
        from sage.combinat.sf.macdonald import qt_kostka
        qt_kostka([2], [1, 1])
        start = time.perf_counter()
        # One pair fills the whole degree, both ways.
        qt_kostka([n], [1] * n)
    elif case == "kfpoly":
        from sage.combinat.sf.kfpoly import kfpoly
        P = Partitions(n).list()
        kfpoly([2], [1, 1])
        start = time.perf_counter()
        [kfpoly(a, b) for a in P for b in P]
    return time.perf_counter() - start


def main():
    if sys.argv[1:2] == ["--run"]:
        case, n, k, arm = sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
        print(run(case, n, None if k == "-" else int(k), arm))
        return
    wanted = set(sys.argv[1:])
    print(f"{'case':<17} {'n':>2} {'k':>2} {'Sage':>9} {'symfn':>9} {'ratio':>8}")
    for case, n, k in CASES:
        if wanted and case not in wanted:
            continue
        times = {}
        for arm in ("sage", "symfn"):
            out = subprocess.run(
                [sys.executable, __file__, "--run", case, str(n),
                 "-" if k is None else str(k), arm],
                env=child_env(case, arm), capture_output=True, text=True,
                check=True)
            times[arm] = float(out.stdout.strip().splitlines()[-1])
        print(f"{case:<17} {n:>2} {k or '':>2} {times['sage']:>8.3f}s "
              f"{times['symfn']:>8.4f}s {times['sage'] / times['symfn']:>7.0f}x")


if __name__ == "__main__":
    main()
