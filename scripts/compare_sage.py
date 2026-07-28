"""Head-to-head timing of symfn against Sage, on identical inputs.

    maturin build --release --features python     # then unpack the wheel into pybuild/
    sage -python scripts/compare_sage.py

This is the Sage counterpart of `compare_lrcalc.py`, and it differs in exactly
one respect: lrcalc is driven **out of process**, one invocation per case, so
both sides pay process startup and neither reuses a warm cache. Sage cannot be
measured that way — its interpreter takes seconds to start, which would swamp
every operation here — so both sides run **in one process** and startup is
excluded by construction. What that costs is the cache isolation a fresh process
gave for free, so it has to be bought back explicitly:

* **Sage memoizes.** Timing the same computation twice measures Sage's cache on
  the second run, not its algorithm. Every case therefore supplies *several
  distinct inputs of comparable size* and each is computed exactly once.
* **symfn memoizes too**, so `symfn.clear_caches()` runs before each timed call.
* **Sage's first touch of a basis is lazy** (building `SymmetricFunctions`,
  coercions, conversion morphisms), so one untimed warm-up runs before any
  measurement. Without it the first case absorbs setup belonging to all of them.

Every case is also *verified*: both sides' outputs are normalised and compared,
and a mismatch is reported instead of a timing. A fast wrong answer is not a
result.

**Not every row has the same baseline, and the difference is large.** Sage's
conversions between the five classical bases are not Python at all: they
dispatch straight into *Symmetrica*'s C (`sage.combinat.sf.classical.init` fills
`conversion_functions` with `t_<FROM>_<TO>_symmetrica`, verifiable at runtime).
The remaining rows are Sage's own Python. So the `C` rows below are a comparison
against optimised C and a 3x there is a strong result, while a 9x on a `py` row
is against an interpreter and means much less.

That distinction is not cosmetic. Symmetrica *also* ships C implementations of
plethysm and Schur products which Sage does **not** use, so those rows were
measuring the weaker of two available baselines. Checked directly against
`sage.libs.symmetrica.all.plethysm`, symfn's plethysm is 0.14-0.79x -- i.e.
mostly slower -- on the cases Symmetrica supports, against 9x here. See
ROADMAP.md.

Sage and Symmetrica are invoked only as separate programs whose output is used;
see NOTICE.md.
"""

import sys
import time

sys.path.insert(0, "pybuild")

import symfn  # noqa: E402
from sage.all import SemistandardTableaux  # noqa: E402
from sage.all import QQ, SymmetricFunctions  # noqa: E402

Sym = SymmetricFunctions(QQ)
s, m, p, e, h = Sym.schur(), Sym.monomial(), Sym.power(), Sym.elementary(), Sym.homogeneous()
f = Sym.forgotten()

REPEATS = int(sys.argv[1]) if len(sys.argv) > 1 else 1
# Per-case budget on Sage's side. Some operations are minutes for Sage at
# degree 20 while symfn is milliseconds, and one of those must not be able to
# wedge the sweep -- the same reason compare_lrcalc.py has a per-invocation
# timeout. A case that exceeds it reports the inputs it managed.
BUDGET = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0


def norm_terms(d):
    """Sage `monomial_coefficients()` -> sorted (parts, coeff) with ints."""
    return sorted((tuple(k), QQ(v)) for k, v in d.items())


def norm_symfn(terms):
    """symfn `[(parts, coeff)]` -> the same normal form."""
    return sorted((tuple(k), QQ(c)) for k, c in terms)


def norm_symfn_rat(terms):
    return sorted((tuple(k), QQ(n) / QQ(d)) for k, (n, d) in terms)


# --- cases -------------------------------------------------------------------
# Each entry: (label, [inputs], sage_fn, symfn_fn). The input lists are distinct
# shapes of comparable difficulty, so neither side is ever asked the same
# question twice.

DEGREES = [int(x) for x in (sys.argv[3].split(",") if len(sys.argv) > 3 else "8,12,16,20")]


def shapes_of(n, count=3):
    """A few distinct partitions of `n`, spanning the shape families.

    Distinct inputs are the whole point: repeating one measures Sage's cache.

    **The families matter as much as the degree.** This function used to emit
    only balanced shapes of 3-6 rows, which kept both l(lambda) and lambda_1
    under about 7 -- and the two Jacobi-Trudi conversions are exponential in
    exactly those, s -> h in l(lambda) and s -> e in lambda_1. So the ladder
    reported s -> e as 3.5-7.6x ahead at every degree while s_(14) took 1.5
    seconds against Symmetrica's 0.02, a case it never generated. The single
    row and the hook below are the two families that were missing; both are
    now first, so a regression in them cannot hide behind an average.
    """
    out, seen = [], set()
    # The extremes first: one row (worst for s -> e), one hook (worst for both,
    # since l and lambda_1 are each about n/2).
    for extreme in ([n], [max(1, n - n // 2)] + [1] * (n // 2)):
        extreme = sorted((x for x in extreme if x > 0), reverse=True)
        if sum(extreme) == n and tuple(extreme) not in seen:
            seen.add(tuple(extreme))
            out.append(extreme)
    for rows in (5, 4, 3, 6):
        base, rem = divmod(n, rows)
        lam = [base + (1 if i < rem else 0) for i in range(rows)]
        lam = [x for x in lam if x > 0]
        # nudge into a strictly different shape when the flat one repeats
        if tuple(lam) in seen and len(lam) > 1 and lam[0] + 1 <= n:
            lam = [lam[0] + 1] + lam[1:-1] + ([lam[-1] - 1] if lam[-1] > 1 else [])
        lam = sorted((x for x in lam if x > 0), reverse=True)
        if sum(lam) == n and tuple(lam) not in seen:
            seen.add(tuple(lam))
            out.append(lam)
        if len(out) == count:
            break
    return out

CASES = [
    (
        "kostka K_{λ,1^n}",
        None,
        # SemistandardTableaux(...).cardinality() is Sage computing a Kostka
        # number. `Partition.dimension()` was wrong here: it is the hook-length
        # formula, a closed form for the 1^n case only, so it measured a
        # different algorithm and flattered Sage by orders of magnitude.
        lambda lam: QQ(SemistandardTableaux(lam, [2, 2] + [1] * (sum(lam) - 4)).cardinality()),
        lambda lam: QQ(symfn.kostka_number(lam, [2, 2] + [1] * (sum(lam) - 4))),
    ),
    (
        "s -> m  (Kostka row)",
        None,
        lambda lam: norm_terms(m(s[lam]).monomial_coefficients()),
        lambda lam: norm_symfn(symfn.schur_to_monomial([(lam, 1)])),
    ),
    (
        "m -> s  (inverse Kostka)",
        None,
        lambda mu: norm_terms(s(m[mu]).monomial_coefficients()),
        lambda mu: norm_symfn(symfn.monomial_to_schur([(mu, 1)])),
    ),
    (
        "p -> s  (Murnaghan-Nakayama)",
        None,
        lambda mu: norm_terms(s(p[mu]).monomial_coefficients()),
        lambda mu: norm_symfn(symfn.power_to_schur([(mu, 1)])),
    ),
    (
        "s -> p",
        None,
        lambda lam: norm_terms(p(s[lam]).monomial_coefficients()),
        lambda lam: norm_symfn_rat(symfn.schur_to_power([(lam, 1)])),
    ),
    (
        "s -> h  (Jacobi-Trudi)",
        None,
        lambda lam: norm_terms(h(s[lam]).monomial_coefficients()),
        lambda lam: norm_symfn(symfn.schur_to_homogeneous([(lam, 1)])),
    ),
    (
        "s -> e  (dual Jacobi-Trudi)",
        None,
        lambda lam: norm_terms(e(s[lam]).monomial_coefficients()),
        lambda lam: norm_symfn(symfn.schur_to_elementary([(lam, 1)])),
    ),
    (
        "s -> f  (forgotten)",
        None,
        lambda lam: norm_terms(f(s[lam]).monomial_coefficients()),
        lambda lam: norm_symfn(symfn.schur_to_forgotten([(lam, 1)])),
    ),
    (
        "f -> s",
        None,
        lambda lam: norm_terms(s(f[lam]).monomial_coefficients()),
        lambda lam: norm_symfn(symfn.forgotten_to_schur([(lam, 1)])),
    ),
    (
        "coproduct",
        None,
        lambda lam: sorted(
            ((tuple(a), tuple(b)), QQ(c))
            for (a, b), c in s[lam].coproduct().monomial_coefficients().items()
        ),
        lambda lam: sorted(
            ((tuple(a), tuple(b)), QQ(c)) for (a, b), c in symfn.coproduct([(lam, 1)])
        ),
    ),
    (
        "skew s_{λ/μ}",
        None,
        lambda lam: norm_terms(s[lam].skew_by(s[[3, 2, 1]]).monomial_coefficients()),
        lambda lam: norm_symfn(symfn.skew_schur(lam, [3, 2, 1])),
    ),
    (
        "plethysm s_a[s_b]",
        "small",
        lambda a: norm_terms(s[a](s[[2, 1]]).monomial_coefficients()),
        lambda a: norm_symfn(symfn.plethysm([(a, 1)], [([2, 1], 1)])),
    ),
    (
        "Kronecker s_λ * s_λ",
        None,
        lambda lam: norm_terms(s[lam].itensor(s[lam]).monomial_coefficients()),
        lambda lam: norm_symfn(symfn.internal_product([(lam, 1)], [(lam, 1)])),
    ),
    (
        "Hall <s_λ, h_λ>",
        None,
        lambda lam: QQ(s[lam].scalar(h[lam])),
        lambda lam: QQ(symfn.hall_inner_product([(lam, 1)], symfn.homogeneous_to_schur([(lam, 1)]))),
    ),
]


# Which side of Sage each case actually lands on. The five classical-basis
# conversions go through Symmetrica's C; everything else is Sage's Python.
# Printed per row so a ratio is never read without its baseline.
#
# The **forgotten** basis is deliberately absent from this set, and that was
# checked rather than assumed: `conversion_functions` holds exactly 20 entries,
# the ordered pairs among {Schur, elementary, homogeneous, monomial, powersum},
# and forgotten appears in none of them. Symmetrica has no forgotten basis at
# all, so `s -> f` and `f -> s` are Sage's own Python and read as `py`.
SYMMETRICA_BACKED = {
    "s -> m  (Kostka row)",
    "m -> s  (inverse Kostka)",
    "p -> s  (Murnaghan-Nakayama)",
    "s -> p",
    "s -> h  (Jacobi-Trudi)",
    "s -> e  (dual Jacobi-Trudi)",
}


def backend_of(label):
    return "C " if label in SYMMETRICA_BACKED else "py"


def warm_up():
    """Untimed: force Sage's lazy basis and coercion setup."""
    lam = [2, 1]
    for basis in (m, p, e, h, f):
        basis(s[lam])
    s(m[lam]), s(p[lam]), s(f[lam])
    s[lam].coproduct()
    s[[3, 2]].skew_by(s[[1]])
    s[[2]](s[[1, 1]])
    s[lam].scalar(h[lam])
    s[[2, 1]].itensor(s[[2, 1]])
    SemistandardTableaux([2, 1], [1, 1, 1]).cardinality()


def main():
    print(f"in-process under Sage, so interpreter startup is excluded entirely")
    print(f"{REPEATS} pass(es); each case uses distinct inputs so neither side answers twice\n")
    warm_up()

    mismatches = []
    for deg in DEGREES:
        run_degree(deg, mismatches)

    print()
    if mismatches:
        for label, inp in mismatches:
            print(f"MISMATCH: {label} on {inp}")
        sys.exit(1)
    print(f"all cases agree with Sage across degrees {DEGREES}")


def run_degree(deg, mismatches):
    print(f"\n=== degree {deg} ===", flush=True)
    print(
        f"{'case':<32} {'via':>3} {'Sage':>10} {'symfn':>10} {'ratio':>9}  {'inputs':>7}",
        flush=True,
    )
    for label, spec, sage_fn, mine_fn in CASES:
        inputs = shapes_of(min(deg, 10)) if spec == "small" else shapes_of(deg)
        t_sage = t_mine = 0.0
        used = 0
        for _ in range(REPEATS):
            for inp in inputs:
                if t_sage > BUDGET:
                    break
                used += 1
                print(f"  .. {label} {inp}".ljust(72), end="\r", flush=True)
                t0 = time.perf_counter()
                want = sage_fn(inp)
                t_sage += time.perf_counter() - t0

                symfn.clear_caches()
                t0 = time.perf_counter()
                got = mine_fn(inp)
                t_mine += time.perf_counter() - t0

                if want != got:
                    mismatches.append((label, inp))
        ratio = f"{t_sage / t_mine:>8.1f}x" if t_mine else f"{'inf':>9}"
        flag = "   <-- MISMATCH" if any(l == label for l, _ in mismatches) else ""
        if t_sage > BUDGET:
            flag += f"  (Sage over {BUDGET:.0f}s budget)"
        print(" " * 72, end="\r")
        print(
            f"{label:<32} {backend_of(label):>3} {t_sage:>9.4f}s {t_mine:>9.4f}s"
            f" {ratio}  {used:>7}{flag}",
            flush=True,
        )



if __name__ == "__main__":
    main()
