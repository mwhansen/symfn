"""Verify Schubert formulas against Sage (oracle) BEFORE they enter the spec.
sage -python verify_formulas.py"""
from sage.all import *

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's Schubert polynomials")

X = SchubertPolynomialRing(ZZ)

def perm_mul_t(w, a, b):
    w = list(w); w[a-1], w[b-1] = w[b-1], w[a-1]
    return Permutation(w)

def length(w): return Permutation(list(w)).length()

# --- 1. LS transition:  S_w = x_r * S_v + sum_{q<r cover} S_{v t_{qr}} ------
# r = LAST descent of w; s = max s>r with w(s) < w(r); v = w t_{rs}.
def transition_check(w):
    w = Permutation(list(w))
    if w.length() == 0: return True
    des = w.descents()          # sage: positions i with w(i)>w(i+1), 1-based?
    r = max(des)
    n = len(w)
    # ensure indexing convention: w.descents() returns i where w(i)>w(i+1), 1-based
    assert w[r-1] > w[r], f"descent convention mismatch {w} {r}"
    cand = [s for s in range(r+1, n+1) if w(s) < w(r)]
    s = max(cand)
    v = perm_mul_t(w, r, s)
    assert v.length() == w.length() - 1
    rhs = X(v).multiply_variable(Integer(r-1))   # x_r * S_v (1-based per probe)
    for q in range(1, r):
        u = perm_mul_t(v, q, r)
        if u.length() == v.length() + 1:
            rhs = rhs + X(u)
    return rhs == X(w)

# pin multiply_variable indexing: does multiply_variable(i) mean x_i with
# sage gens x0,x1,... (0-based) or math x_1.. (1-based)?
w0 = Permutation([2,1,3])
e = X(w0).expand(); R0 = e.parent(); g = R0.gens()
mv1 = X(w0).multiply_variable(Integer(1)).expand()
print("multiply_variable(1) == x0*S?", R0(mv1) == g[0]*e,
      " == x1*S?", R0(mv1) == g[1]*e, flush=True)

ok, bad = 0, []
for w in Permutations(5):
    if w.length() == 0: continue
    try:
        if transition_check(w): ok += 1
        else: bad.append(list(w))
    except AssertionError as e:
        bad.append((list(w), str(e)))
print(f"transition on S_5: {ok} ok, bad: {bad[:5]}", flush=True)

set_random_seed(0)
bad7 = []
for _ in range(20):
    w = Permutations(8).random_element()
    if w.length() == 0: continue
    if not transition_check(w): bad7.append(list(w))
print(f"transition on 20 random S_8: bad: {bad7}", flush=True)

# --- 2. signed Monk: x_i S_w = sum_{b>i cover} S_{w t_{ib}} - sum_{a<i cover} S_{w t_{ai}}
def monk_check(w, i):
    w = Permutation(list(w))
    n = max(len(w), i+1) + 6      # headroom for covers past n
    wl = list(w) + list(range(len(w)+1, n+1))
    lhs = X(w).multiply_variable(Integer(i-1))
    rhs = X.zero()
    for b in range(i+1, n+1):
        u = perm_mul_t(wl, i, b)
        if u.length() == w.length() + 1: rhs += X(u)
    for a in range(1, i):
        u = perm_mul_t(wl, a, i)
        if u.length() == w.length() + 1: rhs -= X(u)
    return lhs == rhs

bad = []
for w in Permutations(4):
    for i in [1, 2, 3, 4]:
        if not monk_check(w, i): bad.append((list(w), i))
print(f"signed Monk on S_4 x i in 1..4: bad: {bad[:5]} ({len(bad)} total)", flush=True)

# --- 3. BJS / pipe dream example from Knutson p.6: C_1423 = x2^2+x1x2+x1^2 (y=0)
e = X([1,4,2,3]).expand()
R = e.parent(); x = R.gens()
print("BJS example 1423:", e == x[0]**2 + x[0]*x[1] + x[1]**2, "->", e, flush=True)

# --- 4. dominant base case: S_w = x^code for dominant w (code (3,1,1): w?) ---
# code (2,1) -> w = [3,2,1]; S = x1^2 x2
e = X([3,2,1]).expand()
print("dominant [3,2,1]:", e == x[0]**2 * x[1], "->", e, flush=True)

# --- 5. stability: S_w == S_{w x 1} (append fixed point) --------------------
w = [3,1,4,2]
print("stability:", X(w) == X(w + [5]), flush=True)

# --- 6. Macdonald reduced-word checksum: S_w(1,...,1) = (1/l!) sum_{red words} prod a_i
def macdonald_count(w):
    w = Permutation(list(w))
    l = w.length()
    tot = 0
    for rw in w.reduced_words():
        p = prod(rw)
        tot += p
    return tot / factorial(l)

for w in [[2,1,4,3], [1,4,2,3], [3,2,1], [2,4,1,3]]:
    ww = Permutation(w)
    spec = X(w).expand()(**{str(g): 1 for g in X(w).expand().parent().gens()})
    mc = macdonald_count(w)
    print(f"Macdonald checksum {w}: S_w(1..1)={spec} formula={mc} match={spec==mc}", flush=True)

print("DONE", flush=True)
