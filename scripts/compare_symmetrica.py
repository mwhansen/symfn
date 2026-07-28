"""Head-to-head timing of symfn against **Symmetrica's C**, on identical inputs.

    maturin build --release --features python     # then unpack the wheel into pybuild/
    python scripts/compare_symmetrica.py [REPEATS] [BUDGET] [DEGREES]

This exists because `compare_sage.py` cannot see half of what it is measuring.

Sage's five classical-basis conversions dispatch into Symmetrica's C, but its
*other* symmetric-function operations are pure Python — even where Symmetrica
ships a C implementation Sage does not call. So those rows compared against the
weaker of two available baselines, and the gap was not small: plethysm read as
**9x faster** than Sage while being **0.14x** against
`sage.libs.symmetrica.all.plethysm` on the same inputs. That defect had been
sitting in a green benchmark.

This harness therefore drives Symmetrica *directly*, through the bindings Sage
installs but mostly leaves unused:

    outerproduct_schur    Littlewood-Richardson products  <- our flagship
    part_part_skewschur   skew Schur functions
    kostka_number         Kostka numbers
    charvalue             symmetric-group characters
    scalarproduct_schur   Hall inner product
    plethysm              plethysm (single-row outer only, see below)

Symmetrica is public domain (see NOTICE.md) and is used here only as a separate
program whose output is read. Every case is **verified**, not merely timed: a
disagreement is reported instead of a ratio, which also makes this a second
independent oracle for LR alongside lrcalc.

Two things to keep in mind when reading the output:

* **Symmetrica refuses some inputs.** Its plethysm handles only a single-row
  outer partition ("for the moment only for outer S_n"). A row restricted that
  way is being measured on Symmetrica's home turf and is noted as such; it is
  not a like-for-like comparison of generality.
* **It reports errors on stdout rather than raising**, and can prompt
  interactively, so every call goes through `guard()` and stdin is closed.
"""

import sys
import time

sys.path.insert(0, "pybuild")

import symfn  # noqa: E402
import sage.libs.symmetrica.all as sym  # noqa: E402
from sage.all import QQ, Partition, SymmetricFunctions  # noqa: E402

Sym = SymmetricFunctions(QQ)
s = Sym.schur()

REPEATS = int(sys.argv[1]) if len(sys.argv) > 1 else 1
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0
DEGREES = [int(x) for x in (sys.argv[3].split(",") if len(sys.argv) > 3 else "8,12,16")]


def guard(fn, *a):
    """Call a Symmetrica binding, turning its refusals into None.

    It signals some unsupported inputs by printing rather than raising, so a
    bare call can return a wrong-looking answer or block on a prompt.
    """
    try:
        return fn(*a)
    except Exception:
        return None


def norm(x):
    """A Schur-basis element from either side -> sorted (parts, coeff)."""
    if x is None:
        return None
    if hasattr(x, "monomial_coefficients"):
        return sorted((tuple(k), QQ(v)) for k, v in x.monomial_coefficients().items())
    return sorted((tuple(k), QQ(c)) for k, c in x)


def shapes_of(n, count=3):
    """A few distinct partitions of n, biased to several rows."""
    out, seen = [], set()
    for rows in (5, 4, 3, 6, 2):
        base, rem = divmod(n, rows)
        lam = sorted(
            (x for x in (base + (1 if i < rem else 0) for i in range(rows)) if x > 0),
            reverse=True,
        )
        if sum(lam) == n and tuple(lam) not in seen:
            seen.add(tuple(lam))
            out.append(lam)
        if len(out) == count:
            break
    return out


def inner_of(lam):
    """A non-trivial partition contained in `lam` (halve the leading rows)."""
    return [x for x in (p // 2 for p in lam[: max(1, len(lam) - 1)]) if x > 0] or [1]


def halves_of(n):
    """Pairs (mu, nu) with |mu| + |nu| = n, for products."""
    a = n // 2
    b = n - a
    return [(m, v) for m in shapes_of(a, 2) for v in shapes_of(b, 2)][:3]


# --- cases -------------------------------------------------------------------
# (label, note, build inputs from degree, symmetrica call, symfn call)

CASES = [
    (
        "LR product s_mu * s_nu",
        "",
        halves_of,
        lambda mu, nu: norm(guard(sym.outerproduct_schur, Partition(mu), Partition(nu))),
        lambda mu, nu: norm(symfn.schur_multiply([(mu, 1)], [(nu, 1)])),
    ),
    (
        "skew s_{lam/mu}",
        "",
        # The inner shape is derived from the outer so that mu subset lam always
        # holds. A fixed inner like [3,2,1] is not contained in every partition
        # of n, and Symmetrica does not return an error for that -- it prints one
        # and aborts the process.
        lambda n: [(lam, inner_of(lam)) for lam in shapes_of(n)],
        lambda lam, mu: norm(guard(sym.part_part_skewschur, Partition(lam), Partition(mu))),
        lambda lam, mu: norm(symfn.skew_schur(lam, mu)),
    ),
    (
        "Kostka K_{lam,mu}",
        "",
        lambda n: [(lam, [2, 2] + [1] * (n - 4)) for lam in shapes_of(n)],
        lambda lam, mu: QQ(guard(sym.kostka_number, Partition(lam), Partition(mu))),
        lambda lam, mu: QQ(symfn.kostka_number(lam, mu)),
    ),
    (
        "character chi^lam(mu)",
        "",
        lambda n: [(lam, [2, 2] + [1] * (n - 4)) for lam in shapes_of(n)],
        lambda lam, mu: QQ(guard(sym.charvalue, Partition(lam), Partition(mu))),
        lambda lam, mu: QQ(symfn.character_value(lam, mu)),
    ),
    (
        "Hall <s_lam, s_lam>",
        "",
        lambda n: [(lam, lam) for lam in shapes_of(n)],
        lambda a, b: QQ(guard(sym.scalarproduct_schur, s[a], s[b])),
        lambda a, b: QQ(symfn.hall_inner_product([(a, 1)], [(b, 1)])),
    ),
    (
        "character table S_n",
        "whole p(n)^2 table",
        lambda n: [(n,)],
        lambda n: sorted(guard(sym.chartafel, n).list()),
        lambda n: sorted(v for row in symfn.character_table(n) for v in row),
    ),
    (
        "Kostka table degree n",
        "whole p(n)^2 table",
        lambda n: [(n,)],
        lambda n: sorted(guard(sym.kostka_tafel, n).list()),
        lambda n: sorted(v for row in symfn.kostka_table(n) for v in row),
    ),
    (
        "plethysm s_a[s_b]",
        "single-row outer only (Symmetrica refuses the rest)",
        lambda n: [([k], b) for k, b in [(3, [2, 1]), (4, [2, 1]), (4, [2, 2])]],
        lambda a, b: norm(guard(sym.plethysm, s[a], s[b])),
        lambda a, b: norm(symfn.plethysm([(a, 1)], [(b, 1)])),
    ),
]


def warm_up():
    """Untimed: Sage's lazy basis setup and Symmetrica's first allocation."""
    guard(sym.outerproduct_schur, Partition([2, 1]), Partition([1]))
    guard(sym.part_part_skewschur, Partition([3, 2]), Partition([1]))
    guard(sym.kostka_number, Partition([2, 1]), Partition([1, 1, 1]))
    guard(sym.charvalue, Partition([2, 1]), Partition([1, 1, 1]))
    guard(sym.scalarproduct_schur, s[[2, 1]], s[[2, 1]])
    guard(sym.plethysm, s[[2]], s[[1]])
    symfn.schur_multiply([([1], 1)], [([1], 1)])


def main():
    print("symfn vs Symmetrica's C, called directly through Sage's bindings")
    print("(Sage's own Python is NOT involved except as a data carrier)\n")
    warm_up()
    mismatches, refused = [], []
    for deg in DEGREES:
        run_degree(deg, mismatches, refused)

    print()
    for label, inp in refused:
        print(f"REFUSED by Symmetrica: {label} on {inp}")
    if mismatches:
        for label, inp in mismatches:
            print(f"MISMATCH: {label} on {inp}")
        sys.exit(1)
    print(f"all comparable cases agree with Symmetrica across degrees {DEGREES}")


def run_degree(deg, mismatches, refused):
    print(f"\n=== degree {deg} ===", flush=True)
    print(f"{'case':<26} {'Symmetrica':>11} {'symfn':>10} {'ratio':>9}  {'n':>3}", flush=True)
    for label, note, build, theirs, mine in CASES:
        inputs = build(deg)
        t_them = t_me = 0.0
        used = 0
        for _ in range(REPEATS):
            for inp in inputs:
                if t_them > BUDGET:
                    break
                print(f"  .. {label} {inp}".ljust(70), end="\r", flush=True)
                t0 = time.perf_counter()
                want = theirs(*inp)
                t_them += time.perf_counter() - t0

                symfn.clear_caches()
                t0 = time.perf_counter()
                got = mine(*inp)
                t_me += time.perf_counter() - t0

                if want is None:
                    refused.append((label, inp))
                    continue
                used += 1
                if want != got:
                    mismatches.append((label, inp))
        ratio = f"{t_them / t_me:>8.2f}x" if t_me and used else f"{'n/a':>9}"
        flag = "   <-- MISMATCH" if any(l == label for l, _ in mismatches) else ""
        if note:
            flag += f"   [{note}]"
        print(" " * 70, end="\r")
        print(
            f"{label:<26} {t_them:>10.5f}s {t_me:>9.5f}s {ratio}  {used:>3}{flag}",
            flush=True,
        )


if __name__ == "__main__":
    main()
