"""Hold the character bases to Sage: transitions and reduced Kronecker products.

    cargo run --release --example stdump -- 6 4 > /tmp/st.txt
    sage -python scripts/check_st.py /tmp/st.txt

Sage is used as a black-box oracle only -- it is run as a separate program and
its output compared, never read for algorithms (see NOTICE.md).

The in-crate tests check the two published expansions from the paper (OZ Eq 20
and Eq 21), the round trips, and the agreement of two independent product
routes. Sage is the third opinion, and the one that would catch a *convention*
shared by both of ours -- the implicit first row, a transposed index, or the
`ht`/`st` pair swapped -- which is exactly what an internal cross-check cannot
see.

Products are compared only where Sage will finish: a per-item SIGALRM skips the
cases past its wall and says so, rather than hanging or silently narrowing the
claim.
"""

import signal
import sys
import time
from collections import defaultdict

from sage.all import SymmetricFunctions, QQ

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's character bases")

path = sys.argv[1] if len(sys.argv) > 1 else "/tmp/st.txt"
LIMIT = int(sys.argv[2]) if len(sys.argv) > 2 else 60

Sym = SymmetricFunctions(QQ)
s = Sym.s()
st = Sym.st()


class Timeout(Exception):
    pass


signal.signal(signal.SIGALRM, lambda *_: (_ for _ in ()).throw(Timeout()))


def with_limit(fn, limit=LIMIT):
    signal.alarm(limit)
    try:
        return fn()
    except Timeout:
        return None
    finally:
        signal.alarm(0)


def parse_partition(text):
    return () if text == "0" else tuple(int(x) for x in text.split(","))


def parse_body(text):
    out = {}
    if text:
        for item in text.split():
            key, val = item.rsplit(":", 1)
            out[parse_partition(key)] = int(val)
    return out


def sage_terms(elt):
    return {tuple(k): int(v) for k, v in elt.monomial_coefficients().items() if v}


sections = defaultdict(dict)
for line in open(path):
    kind, key, rest = (f.strip() for f in line.split("|", 2))
    if kind == "prod":
        a, b = key.split(";")
        sections[kind][(parse_partition(a), parse_partition(b))] = parse_body(rest)
    else:
        sections[kind][parse_partition(key)] = parse_body(rest)

bad = skipped = checked = 0


def report(label, got, want):
    global bad
    if got == want:
        return True
    bad += 1
    print(f"MISMATCH {label}")
    for k in sorted(set(got) | set(want)):
        if got.get(k, 0) != want.get(k, 0):
            print(f"    {k}: ours={got.get(k, 0)} sage={want.get(k, 0)}")
    return False


print(f"=== st -> s  ({len(sections['st2s'])} partitions) ===")
for lam, got in sorted(sections["st2s"].items(), key=lambda kv: (sum(kv[0]), kv[0])):
    want = sage_terms(s(st[list(lam)] if lam else st.one()))
    checked += 1
    report(f"st{list(lam)} -> s", got, want)

print(f"=== s -> st  ({len(sections['s2st'])} partitions) ===")
for lam, got in sorted(sections["s2st"].items(), key=lambda kv: (sum(kv[0]), kv[0])):
    want = sage_terms(st(s[list(lam)] if lam else s.one()))
    checked += 1
    report(f"s{list(lam)} -> st", got, want)

print(f"=== products ({len(sections['prod'])} pairs) ===")
for (a, b), got in sorted(sections["prod"].items(), key=lambda kv: (sum(kv[0][0]) + sum(kv[0][1]), kv[0])):
    left = st[list(a)] if a else st.one()
    right = st[list(b)] if b else st.one()
    t0 = time.time()
    prod = with_limit(lambda: left * right)
    dt = time.time() - t0
    if prod is None:
        skipped += 1
        print(f"  SKIP st{list(a)} * st{list(b)} — sage exceeded {LIMIT}s "
              f"(ours: {len(got)} terms)")
        continue
    checked += 1
    if report(f"st{list(a)} * st{list(b)}", got, sage_terms(prod)) and sum(a) + sum(b) >= 7:
        print(f"  ok   st{list(a)} * st{list(b)}  {len(got):>5} terms, sage took {dt:.2f}s")

print()
print(f"checked {checked}, mismatches {bad}, skipped (sage too slow) {skipped}")
sys.exit(1 if bad else 0)
