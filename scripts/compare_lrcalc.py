"""Head-to-head timing of symfn against lrcalc, on identical inputs.

    cargo build --release --example lr_cli
    LRCALC=/path/to/lrcalc python3 scripts/compare_lrcalc.py

Both sides run as a separate process per case, so both pay process startup and
both format their full output to stdout. Neither reuses a warm memo cache across
cases. That is the fair comparison; anything measured in-process on our side and
out-of-process on theirs would flatter us.

Every case is also *verified*: the two outputs are parsed, sorted and compared,
and a mismatch is reported as a failure rather than a timing. A fast wrong
answer is not a result.

lrcalc is invoked as an external program and only its output is used; see
NOTICE.md.
"""

import os
import re
import subprocess
import sys
import time

LRCALC = os.environ.get("LRCALC", "lrcalc")
LR_CLI = os.environ.get(
    "LR_CLI", os.path.join(os.path.dirname(__file__), "..", "target", "release", "examples", "lr_cli")
)
REPEATS = int(os.environ.get("REPEATS", "3"))
# Per-invocation ceiling. Some shapes are pathological for one side or the
# other -- lrcalc needs well over 20 minutes on [24,20,16,12]^2 -- and a sweep
# that blocks on one case reports nothing at all.
TIMEOUT = float(os.environ.get("TIMEOUT", "120"))

TERM = re.compile(r"^\s*(\d+)\s+\((.*)\)\s*$")


def parse_terms(out):
    """Parse `COEFF  (parts)` lines into a sorted canonical form."""
    terms = []
    for line in out.splitlines():
        if not line.strip():
            continue
        m = TERM.match(line)
        if not m:
            return None  # `coef` output: a bare integer
        parts = m.group(2)
        terms.append((tuple(int(x) for x in parts.split(",")) if parts else (), int(m.group(1))))
    return sorted(terms)


def run(cmd):
    """Run once, returning (best-of-REPEATS seconds, stdout).

    Returns (None, None) if the command exceeds TIMEOUT. A single pathological
    case must not be able to wedge the whole sweep -- and "did not finish in
    TIMEOUT" is itself a result worth reporting, so it is recorded rather than
    treated as an error.
    """
    best, out = None, None
    for _ in range(REPEATS):
        t0 = time.perf_counter()
        try:
            r = subprocess.run(cmd, capture_output=True, text=True, timeout=TIMEOUT)
        except subprocess.TimeoutExpired:
            return None, None
        dt = time.perf_counter() - t0
        if r.returncode != 0:
            sys.exit(f"failed: {' '.join(cmd)}\n{r.stderr}")
        best = dt if best is None else min(best, dt)
        out = r.stdout
    return best, out


def canonical(out):
    terms = parse_terms(out)
    return terms if terms is not None else out.strip()


CASES = [
    # (label, argv-tail shared by both CLIs)
    # --- products: staircases, the classic stress case ---
    ("mult s[4,3,2,1]^2", ["mult", "4", "3", "2", "1", "-", "4", "3", "2", "1"]),
    ("mult s[5,4,3,2,1]^2", ["mult", "5", "4", "3", "2", "1", "-", "5", "4", "3", "2", "1"]),
    ("mult s[6,5,4,3,2,1]^2", ["mult", "6", "5", "4", "3", "2", "1", "-", "6", "5", "4", "3", "2", "1"]),
    ("mult s[7,6,5,4,3]^2", ["mult", "7", "6", "5", "4", "3", "-", "7", "6", "5", "4", "3"]),
    ("mult s[8,7,6,5,4,3]^2", ["mult", "8", "7", "6", "5", "4", "3", "-", "8", "7", "6", "5", "4", "3"]),
    ("mult s[8,7,6,5,4]^2", ["mult", "8", "7", "6", "5", "4", "-", "8", "7", "6", "5", "4"]),
    ("mult s[9,8,7,6,5]^2", ["mult", "9", "8", "7", "6", "5", "-", "9", "8", "7", "6", "5"]),
    # --- products: other regimes, where the shape is not a staircase.
    #     Sized to clear the ~6ms process-startup floor, or the row measures
    #     exec() rather than either algorithm. ---
    ("mult rectangle [5^5]^2", ["mult", "5", "5", "5", "5", "5", "-", "5", "5", "5", "5", "5"]),
    ("mult rectangle [7^7]^2", ["mult"] + ["7"] * 7 + ["-"] + ["7"] * 7),
    # Rectangles are the one shape with a closed form (see `src/rect.rs`), so
    # the sweep needs one large enough for the algorithm rather than `exec` to
    # dominate: [5^5]^2 and [7^7]^2 above are both inside the startup floor.
    ("mult rectangle [12^6]^2", ["mult"] + ["12"] * 6 + ["-"] + ["12"] * 6),
    ("mult rectangle [14^7]^2", ["mult"] + ["14"] * 7 + ["-"] + ["14"] * 7),
    ("mult wide [12,10,8]^2", ["mult", "12", "10", "8", "-", "12", "10", "8"]),
    ("mult wide [14,12,10]^2", ["mult", "14", "12", "10", "-", "14", "12", "10"]),
    ("mult wide [20,16,12]^2", ["mult", "20", "16", "12", "-", "20", "16", "12"]),
    # Four wide rows: seconds rather than milliseconds, so it separates the two
    # implementations well clear of the startup floor and of run-to-run noise.
    ("mult wide [16,13,10,7]^2", ["mult", "16", "13", "10", "7", "-", "16", "13", "10", "7"]),
    ("mult wide [24,20,16,12]^2", ["mult", "24", "20", "16", "12", "-", "24", "20", "16", "12"]),
    ("mult tall [2^8]^2", ["mult"] + ["2"] * 8 + ["-"] + ["2"] * 8),
    ("mult tall [3^12]^2", ["mult"] + ["3"] * 12 + ["-"] + ["3"] * 12),
    ("mult two-row [20,10]^2", ["mult", "20", "10", "-", "20", "10"]),
    ("mult lopsided [9,8,7,6]x[3,2,1]", ["mult", "9", "8", "7", "6", "-", "3", "2", "1"]),
    # --- skew expansions ---
    ("skew [8,7,6,5,4]/[3,2,1]", ["skew", "8", "7", "6", "5", "4", "/", "3", "2", "1"]),
    ("skew [9,9,8,8,7]/[2,1]", ["skew", "9", "9", "8", "8", "7", "/", "2", "1"]),
    ("skew [10,9,8,7,6,5]/[4,3,2,1]", ["skew", "10", "9", "8", "7", "6", "5", "/", "4", "3", "2", "1"]),
    ("skew [12,11,10,9,8,7]/[5,4,3]", ["skew", "12", "11", "10", "9", "8", "7", "/", "5", "4", "3"]),
    # --- single coefficient: the regime where computing a whole expansion
    #     to answer one question could plausibly lose ---
    (
        "coef c^[8,7,6,5,4,3]_[4,3,2,1],[4,3,2,1]",
        ["coef", "8", "7", "6", "5", "4", "3", "-", "4", "3", "2", "1", "-", "4", "3", "2", "1"],
    ),
    (
        "coef c^[12,10,8,6]_[6,5,4,3],[6,5,4,3]",
        ["coef", "12", "10", "8", "6", "-", "6", "5", "4", "3", "-", "6", "5", "4", "3"],
    ),
]


def fmt(t):
    return f">{TIMEOUT:.0f}s".rjust(10) if t is None else f"{t:>9.4f}s"


print(f"best of {REPEATS}, both as separate processes, full output formatted")
print(f"per-invocation timeout {TIMEOUT:.0f}s\n")
print(f"{'case':<38} {'lrcalc':>10} {'symfn':>10} {'ratio':>9}  {'terms':>8}", flush=True)
mismatches, compared = [], 0
for label, argv in CASES:
    t_lr, out_lr = run([LRCALC] + argv)
    t_sy, out_sy = run([LR_CLI] + argv)

    if out_lr is None or out_sy is None:
        # One side timed out, so there is nothing to verify and no ratio.
        who = "lrcalc" if out_lr is None else "symfn"
        print(f"{label:<38} {fmt(t_lr)} {fmt(t_sy)} {'--':>9}  {'--':>8}   ({who} timed out)", flush=True)
        continue

    ok = canonical(out_lr) == canonical(out_sy)
    compared += 1
    if not ok:
        mismatches.append(label)
    n = len([l for l in out_lr.splitlines() if l.strip()])
    ratio = f"{t_lr / t_sy:>8.2f}x" if t_sy else f"{'inf':>9}"
    flag = "" if ok else "   <-- MISMATCH"
    print(f"{label:<38} {fmt(t_lr)} {fmt(t_sy)} {ratio}  {n:>8}{flag}", flush=True)

print()
if mismatches:
    print(f"{len(mismatches)} MISMATCH(es): {', '.join(mismatches)}")
    sys.exit(1)
print(f"all {compared} completed cases agree with lrcalc")
