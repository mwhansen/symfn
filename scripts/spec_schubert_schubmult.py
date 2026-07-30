"""Time Buch's C schubmult (lrcalc suite, bundled with Sage) on the SAME cases
as schubert_walls2.py. Out of process, one timing per case (walls, not benchmarks);
small cases get min-of-3. sage -python schubmult_walls.py"""
import subprocess, time
from sage.all import *

SB = "/var/tmp/sage-10.9-current/local/bin/schubmult"

def run(u, v, limit=120):
    args = [SB] + [str(a) for a in u] + ["-"] + [str(a) for a in v]
    t0 = time.perf_counter()
    try:
        r = subprocess.run(args, capture_output=True, text=True, timeout=limit)
        dt = time.perf_counter() - t0
        lines = [l for l in r.stdout.splitlines() if l.strip()]
        return dt, len(lines)
    except subprocess.TimeoutExpired:
        return None, None

def bench(label, u, v):
    dt, nt = run(u, v)
    if dt is None:
        print(f"  {label}: TIMEOUT >120s", flush=True)
        return
    if dt < 0.5:  # min-of-3 for small ones
        for _ in range(2):
            d2, _ = run(u, v)
            if d2 is not None: dt = min(dt, d2)
    print(f"  {label}: {dt:.4f}s  terms={nt}", flush=True)

def stair(k):
    return [2*i for i in range(1, k+1)] + [2*i-1 for i in range(1, k+1)]

print("== schubmult (lrcalc C): staircase Grassmannian squares ==", flush=True)
for k in [4, 5, 6, 7]:
    w = stair(k)
    bench(f"stair{k}^2 (S_{2*k})", w, w)

print("\n== schubmult: same seed-1 random pairs ==", flush=True)
set_random_seed(1)
pairs = {}
for n in [10, 11, 12, 13]:
    pairs[n] = [(Permutations(n).random_element(),
                 Permutations(n).random_element()) for _ in range(3)]
for n in [10, 11, 12, 13]:
    for i, (u, v) in enumerate(pairs[n]):
        bench(f"S_{n}.{i} l(u)={u.length()} l(v)={v.length()}", list(u), list(v))

print("DONE", flush=True)
