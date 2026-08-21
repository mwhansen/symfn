"""Sharper Schubert product ladder: S_10..S_13, staircase Grassmannians,
long random pairs. sage -python schubert_walls2.py"""
import time
from sage.all import *
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt

X = SchubertPolynomialRing(ZZ)

def timed(label, f, limit=120):
    t0 = time.perf_counter()
    try:
        alarm(limit)
        r = f()
        cancel_alarm()
        return r, time.perf_counter() - t0
    except AlarmInterrupt:
        print(f"  {label}: TIMEOUT >{limit}s", flush=True)
        return None, None

def stats(el):
    mc = el.monomial_coefficients()
    return len(mc), max(abs(c) for c in mc.values())

# staircase Grassmannian: shape (k,k-1,...,1), descent k -> w(i)=2i, in S_2k
def stair(k):
    return [2*i for i in range(1, k+1)] + [2*i-1 for i in range(1, k+1)]

print("== staircase Grassmannian squares (S_w = s_{(k..1)}(x_1..x_k)) ==", flush=True)
for k in [4, 5, 6]:
    w = stair(k)
    a = X(w)
    r, dt = timed(f"stair{k}^2 (S_{2*k})", lambda: a*a)
    if r is not None:
        nt, mx = stats(r)
        print(f"  stair{k}^2 (S_{2*k}): {dt:.4f}s  terms={nt}  maxcoeff={mx}", flush=True)

print("\n== random long pairs, ladder ==", flush=True)
set_random_seed(1)
for n in [10, 11, 12, 13]:
    hit_wall = False
    for i in range(3):
        u = Permutations(n).random_element()
        v = Permutations(n).random_element()
        a, b = X(list(u)), X(list(v))
        label = f"S_{n}.{i} l(u)={u.length()} l(v)={v.length()}"
        r, dt = timed(label, lambda: a*b)
        if r is None:
            hit_wall = True
            continue
        nt, mx = stats(r)
        print(f"  {label}: {dt:.4f}s  terms={nt}  maxcoeff={mx}", flush=True)
    if hit_wall and n >= 12:
        break

print("\n== expand, long non-dominant ==", flush=True)
for k in [5, 6]:
    w = stair(k)
    el = X(w)
    r, dt = timed(f"expand stair{k}", lambda: el.expand())
    if r is not None:
        print(f"  expand stair{k}: {dt:.4f}s  monomials={len(r.monomials())}", flush=True)
set_random_seed(2)
for n in [11, 12]:
    w = Permutations(n).random_element()
    el = X(list(w))
    r, dt = timed(f"expand rand S_{n} l={w.length()}", lambda: el.expand())
    if r is not None:
        print(f"  expand rand S_{n} l={w.length()}: {dt:.4f}s  monomials={len(r.monomials())}", flush=True)

print("\n== newtrans (Stanley symmetric function) scaling ==", flush=True)
import sage.libs.symmetrica.all as symca

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's Schubert polynomials")
set_random_seed(3)
for n in [8, 10, 12]:
    w = Permutations(n).random_element()
    r, dt = timed(f"newtrans S_{n} l={w.length()}", lambda: symca.newtrans(Permutation(list(w))), limit=60)
    if r is not None:
        print(f"  newtrans S_{n} l={w.length()}: {dt:.4f}s  terms={len(r)}", flush=True)

print("\nDONE", flush=True)
