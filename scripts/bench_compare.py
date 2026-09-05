"""Compare two `bench_suite` runs, workload by workload.

    cargo run --release --example bench_suite > after.tsv
    python3 scripts/bench_compare.py docs/record/bench_suite.tsv after.tsv

Prints one line per workload with both times and the ratio `after / before`,
and marks the ratios outside the noise band. The band is ±20%, the line the
LR record draws for a single run on this machine
(`docs/record/littlewood-richardson.md`: "treat anything inside ±20% as a
tie"). A ratio outside it on one workload is a lead to measure properly with
that subsystem's own `bench_*` example; it is not a result by itself.

A workload under 5 ms is printed and never flagged: at that scale the timer
and the scheduler are what is being compared, which is why the LR record's
out-of-process sweep discards its rows under the startup floor. `schubert` in
the catalog is such a workload — its product is one term, sized for the
memory budget rather than for time.

Exits nonzero when any workload is outside the band, so a shell loop can
notice, and prints nothing else: the interpretation belongs to whoever ran it.
Both files must come from the same machine and power state, which the header
lines of a committed run record.
"""

import sys

BAND = 0.20
FLOOR = 0.005


def read(path):
    rows = {}
    header = []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("#"):
                header.append(line.lstrip("# "))
                continue
            if not line or line.startswith("workload\t"):
                continue
            name, seconds, *_ = line.split("\t")
            rows[name] = float(seconds)
    return header, rows


def main():
    if len(sys.argv) != 3:
        sys.exit(f"usage: {sys.argv[0]} before.tsv after.tsv")
    before_header, before = read(sys.argv[1])
    after_header, after = read(sys.argv[2])
    print("before: " + " | ".join(before_header))
    print("after:  " + " | ".join(after_header))
    print(f"{'workload':<16} {'before':>9} {'after':>9} {'ratio':>7}")
    outside = 0
    for name, b in before.items():
        a = after.get(name)
        if a is None:
            print(f"{name:<16} {b:>9.4f} {'missing':>9}")
            continue
        if b < FLOOR and a < FLOOR:
            print(f"{name:<16} {b:>9.4f} {a:>9.4f}   under {FLOOR * 1000:.0f} ms")
            continue
        ratio = a / b if b > 0 else float("inf")
        flag = ""
        if abs(ratio - 1.0) > BAND:
            flag = "  outside the band"
            outside += 1
        print(f"{name:<16} {b:>9.4f} {a:>9.4f} {ratio:>6.2f}x{flag}")
    for name in after:
        if name not in before:
            print(f"{name:<16} {'missing':>9} {after[name]:>9.4f}")
    sys.exit(1 if outside else 0)


if __name__ == "__main__":
    main()
