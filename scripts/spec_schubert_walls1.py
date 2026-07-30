"""Measure Sage/Symmetrica Schubert-polynomial walls. Run: sage -python schubert_walls.py
Single process; per-case cysignals alarm; incremental flushed output so a hard
crash preserves earlier rows. Order-of-magnitude walls, not benchmarks.
"""
import time, sys
from sage.all import *
from cysignals.alarm import alarm, cancel_alarm, AlarmInterrupt

X = SchubertPolynomialRing(ZZ)
P = Permutation

def timed(label, f, limit=90):
    t0 = time.perf_counter()
    try:
        alarm(limit)
        r = f()
        cancel_alarm()
        dt = time.perf_counter() - t0
        return r, dt
    except AlarmInterrupt:
        print(f"  {label}: TIMEOUT >{limit}s", flush=True)
        return None, None
    except Exception as e:
        cancel_alarm()
        print(f"  {label}: ERROR {type(e).__name__}: {e}", flush=True)
        return None, None

def stats(el):
    mc = el.monomial_coefficients()
    return len(mc), max(abs(c) for c in mc.values())

print("== sage version:", version(), flush=True)

# --- A. products ---------------------------------------------------------
print("\n== A. products u*v (through SchubertPolynomialRing) ==", flush=True)

# box Grassmannian in S_{2k}: S_w = s_{(k^k)}(x_1..x_k); product = LR rectangle^2
def box(k):
    return list(range(k+1, 2*k+1)) + list(range(1, k+1))

# staircase-ish 'middle' permutation families, deterministic:
set_random_seed(0)
rand_pairs = {}
for n in [6, 7, 8, 9, 10]:
    rand_pairs[n] = [(Permutations(n).random_element(),
                      Permutations(n).random_element()) for _ in range(2)]

cases = []
for k in [3, 4, 5, 6]:
    w = box(k); cases.append((f"box{k}xbox{k} (S_{2*k}, rect LR {k}^{k})", w, w))
for n in [6, 7, 8, 9, 10]:
    for i, (u, v) in enumerate(rand_pairs[n]):
        cases.append((f"rand{n}.{i} l(u)={u.length()} l(v)={v.length()}",
                      list(u), list(v)))
# dominant sanity: w0 * w0 (should be instant, 1 term)
cases.append(("w0(8)*w0(8)", list(Permutation([8,7,6,5,4,3,2,1])),
              list(Permutation([8,7,6,5,4,3,2,1]))))

for label, u, v in cases:
    a, b = X(u), X(v)
    r, dt = timed(label, lambda: a*b)
    if r is not None:
        nt, mx = stats(r)
        print(f"  {label}: {dt:.4f}s  terms={nt}  maxcoeff={mx}", flush=True)

# --- B. expand (t_SCHUBERT_POLYNOM) --------------------------------------
print("\n== B. expand to monomials ==", flush=True)
for label, w in [("box4 (S_8)", box(4)), ("box5 (S_10)", box(5)),
                 ("rand9.0u", list(rand_pairs[9][0][0])),
                 ("rand10.0u", list(rand_pairs[10][0][0])),
                 ("w0(9)", [9,8,7,6,5,4,3,2,1])]:
    el = X(w)
    r, dt = timed(label, lambda: el.expand())
    if r is not None:
        print(f"  {label}: {dt:.4f}s  monomials={len(r.monomials())}", flush=True)

# --- C. polynomial -> Schubert (t_POLYNOM_SCHUBERT) ----------------------
print("\n== C. from-polynomial roundtrip ==", flush=True)
import sage.libs.symmetrica.all as symca
for label, w in [("box4", box(4)), ("rand8.0u", list(rand_pairs[8][0][0])),
                 ("rand9.0u", list(rand_pairs[9][0][0]))]:
    pol = X(w).expand()
    r, dt = timed(label, lambda: symca.t_POLYNOM_SCHUBERT(pol))
    if r is not None:
        print(f"  {label}: {dt:.4f}s  (roundtrip ok: {r == X(w)})", flush=True)

# --- D. divided differences ----------------------------------------------
print("\n== D. divided_difference on a product ==", flush=True)
big = None
try:
    alarm(60); big = X(box(4)) * X(box(4)); cancel_alarm()
except AlarmInterrupt:
    print("  (skipped: box4^2 timed out)", flush=True)
if big is not None:
    r, dt = timed("d_1(box4^2)", lambda: big.divided_difference(1))
    if r is not None:
        print(f"  d_1(box4^2): {dt:.4f}s", flush=True)

# --- E. API surface (black box) -------------------------------------------
print("\n== E. element API surface ==", flush=True)
el = X([2,1,4,3])
print("  element methods:", sorted(m for m in dir(el)
      if not m.startswith('_') and m not in dir(X)), flush=True)

# --- F. newtrans probe (LAST: may crash) ----------------------------------
print("\n== F. newtrans probe ==", flush=True)
for w in [[2,1,4,3], [3,2,1], [1,4,2,3]]:
    r, dt = timed(f"newtrans{w}", lambda: symca.newtrans(Permutation(w)), limit=30)
    if r is not None:
        print(f"  newtrans({w}) = {r}  [{dt:.4f}s]", flush=True)

print("\nDONE", flush=True)
