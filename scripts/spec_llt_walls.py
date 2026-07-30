"""Sage LLT walls for docs/spec-llt.md sections 2.2 and 2.3.

Round 2 of 2: exact wall degrees, per-shape split, level 4, single-shape
cache-wart demonstration, tuple walls, cospin_polynomial coefficients.

Run twice, and both runs are in the spec because the difference is the point:

  * 2026-07-29, **battery ~96%**  -> section 2.2
  * 2026-07-30, **mains**         -> section 2.3, and this is the one the
    section 6 ratios are computed against

Mains came out ~2x faster across the board, which is the drift figure the
provenance notes had been assuming, now measured on this very workload. Two
items that TIMED OUT on battery complete on mains — level 4 at n=8 (101.45s)
and the n=12 tuple (66.99s) — so section 2.3 has real numbers where 2.2 had
bounds. Anything quoted as a speedup must come from one power state on both
sides."""
import signal
import time

from sage.all import QQ, Partitions, SymmetricFunctions, FractionField


class Timeout(Exception):
    pass


signal.signal(signal.SIGALRM, lambda s, f: (_ for _ in ()).throw(Timeout()))

LIMIT = 120


def timed(label, thunk, limit=LIMIT):
    signal.alarm(limit)
    t0 = time.time()
    try:
        val = thunk()
        dt = time.time() - t0
        signal.alarm(0)
        print(f"  {label:52s} {dt:9.2f}s", flush=True)
        return dt, val
    except Timeout:
        print(f"  {label:52s}   >{limit}s TIMEOUT", flush=True)
        return None, None
    except Exception as e:
        signal.alarm(0)
        print(f"  {label:52s}   ERROR {type(e).__name__}: {e}", flush=True)
        return None, None


print("== whole-degree walls, levels 2/3/4, fresh parent per (k,n) ==", flush=True)
for k in (2, 3, 4):
    for n in range(7, 12):
        Rn = FractionField(QQ['t'])
        Symn = SymmetricFunctions(Rn)
        sn = Symn.schur()
        H = Symn.llt(k).hspin()

        def whole_degree(H=H, sn=sn, n=n):
            per_shape = []
            for mu in Partitions(n):
                t0 = time.time()
                sn(H[mu])
                per_shape.append((time.time() - t0, mu))
            per_shape.sort(reverse=True)
            return per_shape

        dt, per = timed(f"level {k}: whole degree n={n}", whole_degree)
        if dt is not None and per:
            worst = per[0]
            print(f"      worst shape {list(worst[1])}: {worst[0]:.2f}s "
                  f"({100 * worst[0] / max(dt, 1e-9):.0f}% of degree)", flush=True)
        if dt is None:
            break

print("== single-shape wall, fresh parent each, level 3, mu=(n) ==", flush=True)
for n in range(6, 11):
    Rn = FractionField(QQ['t'])
    Symn = SymmetricFunctions(Rn)
    sn = Symn.schur()
    H = Symn.llt(3).hspin()
    dt, _ = timed(f"level 3: single HSp3[[{n}]] -> s", lambda H=H, sn=sn, n=n: sn(H[[n]]))
    if dt is None:
        break

print("== tuple LLT walls: 3-tuples growing n ==", flush=True)
for tup in ([[2, 2], [2, 1], [2]], [[2, 2], [2, 2], [2]], [[3, 2], [2, 2], [1]],
            [[2, 2], [2, 2], [2, 1]], [[3, 2], [2, 2], [2, 1]]):
    Rn = FractionField(QQ['t'])
    Symn = SymmetricFunctions(Rn)
    lltk = Symn.llt(len(tup))
    total = sum(sum(p) for p in tup)
    dt, _ = timed(f"cospin{tup} (n={total})", lambda lltk=lltk, tup=tup: lltk.cospin(tup))

print("== D: cospin_polynomial single coefficients (fixed) ==", flush=True)
from sage.combinat.ribbon_tableau import cospin_polynomial
cases = [
    ([[6, 6, 4, 2], []], [1] * 9, 2),
    ([[6, 6, 6, 4, 2], []], [1] * 12, 2),
    ([[8, 8, 6, 4, 2], []], [1] * 14, 2),
    ([[9, 6, 3], []], [1] * 6, 3),
    ([[9, 9, 6, 3], []], [1] * 9, 3),
    ([[12, 9, 6, 3], []], [1] * 10, 3),
]
for shape, weight, k in cases:
    dt, val = timed(f"cospin_poly {shape[0]} wt 1^{len(weight)} k={k}",
                    lambda shape=shape, weight=weight, k=k: cospin_polynomial(shape, weight, k))
    if val is not None:
        print(f"      = {val}", flush=True)

print("done", flush=True)
