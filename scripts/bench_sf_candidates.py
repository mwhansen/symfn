"""Time the `sage.combinat.sf` paths that do not reach symfn, against symfn.

    python scripts/bench_sf_candidates.py            # every case
    python scripts/bench_sf_candidates.py llt_spin   # one case

Run with the sage-dev environment's `python`; this Sage's CLI has no
`sage -python`.

Two kinds of row, and they are not the same claim:

- **Kernel rows** (`qt_kostka`, `kfpoly`) time Sage's own function against the
  bare symfn entry point, with no conversion of the answer into Sage objects.
  They bound what routing the path through symfn could save; they are not an
  end-to-end speedup, because the marshalling is exactly the cost this backend
  keeps finding once the algorithm is gone.
- **End-to-end rows** (`nabla`, `llt_spin`, `llt_cospin`) time one Sage-level
  operation both ways and include everything the caller pays. `nabla` compares
  `f.nabla()` with `sage.libs.symfn.extras.nabla`. The two LLT rows time
  `_m_cache(n)` -- both directions of the level-`k` basis at degree `n` -- in a
  child with `SAGE_DISABLE_SYMFN=1` against one without, which is only an A/B
  once `llt.py` dispatches to symfn; before that both arms are Sage.

Sage's arm of `qt_kostka` runs with the backend on, since its Macdonald `H`
expansion already comes from symfn: the row measures the wall that remains,
the same reading as the standing list in
`docs/record/python-and-sage-interop.md`.

Every timing is one cold run in a fresh process: Sage memoizes changes of basis
at module level, and `symfn.clear_caches()` clears the kernel's own.

⚠️ Record the power state (`pmset -g batt`). This machine drifts about 1.8x on
battery.
"""

import os
import subprocess
import sys
import time

# (case, degree, level) -- the level is used by the LLT rows only.
CASES = [
    ("qt_kostka", 5, None), ("qt_kostka", 7, None),
    ("kfpoly", 8, None), ("kfpoly", 10, None),
    ("nabla", 5, None), ("nabla", 7, None),
    ("llt_spin", 6, 3), ("llt_spin", 7, 3), ("llt_spin", 8, 3),
    ("llt_spin", 9, 2),
    ("llt_cospin", 7, 3), ("llt_cospin", 8, 3),
]


def child_env(case, arm):
    """`SAGE_DISABLE_SYMFN` is read when Sage imports, so it is set here.

    Only the LLT rows turn the backend off: the others call a Sage function
    that does not dispatch, and turning it off there would also move the
    classical conversions underneath them back onto Symmetrica.
    """
    env = dict(os.environ)
    if arm == "sage" and case in ("llt_spin", "llt_cospin"):
        env["SAGE_DISABLE_SYMFN"] = "1"
    else:
        env.pop("SAGE_DISABLE_SYMFN", None)
    return env


def run(case, n, k, arm):
    import symfn
    from sage.all import QQ, Partitions, SymmetricFunctions

    symfn.clear_caches()
    if case in ("llt_spin", "llt_cospin"):
        from sage.libs.symfn import is_available
        # A control arm that silently reached symfn would print a ratio near
        # 1.0x and prove nothing.
        assert is_available() == (arm == "symfn"), (case, arm)
        L = SymmetricFunctions(QQ["t"].fraction_field()).llt(k)
        B = L.hspin() if case == "llt_spin" else L.hcospin()
        start = time.perf_counter()
        B._m_cache(n)
    elif case == "nabla":
        s = SymmetricFunctions(QQ["q,t"].fraction_field()).s()
        f = sum(s(la) for la in Partitions(n))
        start = time.perf_counter()
        if arm == "sage":
            f.nabla()
        else:
            from sage.libs.symfn.extras import nabla
            nabla(f)
    elif case == "qt_kostka":
        start = time.perf_counter()
        if arm == "sage":
            from sage.combinat.sf.macdonald import qt_kostka
            P = Partitions(n).list()
            [qt_kostka(a, b) for a in P for b in P]
        else:
            symfn.qt_kostka_table(n)
    elif case == "kfpoly":
        start = time.perf_counter()
        if arm == "sage":
            from sage.combinat.sf.kfpoly import kfpoly
            P = Partitions(n).list()
            [kfpoly(a, b) for a in P for b in P]
        else:
            symfn.kostka_foulkes_table(n)
    return time.perf_counter() - start


def main():
    if sys.argv[1:2] == ["--run"]:
        case, n, k, arm = sys.argv[2], int(sys.argv[3]), sys.argv[4], sys.argv[5]
        print(run(case, n, None if k == "-" else int(k), arm))
        return
    wanted = set(sys.argv[1:])
    print(f"{'case':<11} {'n':>2} {'k':>2} {'Sage':>9} {'symfn':>9} {'ratio':>8}")
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
        print(f"{case:<11} {n:>2} {k or '':>2} {times['sage']:>8.3f}s "
              f"{times['symfn']:>8.4f}s {times['sage'] / times['symfn']:>7.0f}x")


if __name__ == "__main__":
    main()
