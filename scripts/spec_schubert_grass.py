"""Settle the ONE unverified formula in docs/record/schubert-spec.md (§3.5 leaf
dispatch, §7 open question 2): the Grassmannian flag-truncation identity.

Claim under test.  If u, v are Grassmannian with the SAME descent k, with
shapes lambda, mu (at most k rows), then

    S_u * S_v = sum_{nu : ell(nu) <= k} c^nu_{lambda,mu} * S_{w(nu,k)}

i.e. the ordinary LR expansion with every output partition of more than k rows
dropped, each surviving nu re-indexed as the Grassmannian permutation of
descent k and shape nu.  If this holds, transition-tree leaves that are
Grassmannian pairs can be handed to AutoLr.

Also pins the two index conventions the claim rides on:
  (a) shape <-> Grassmannian permutation, and
  (b) S_w = s_lambda(x_1..x_k) for such w.

Run: sage -python scripts/spec_schubert_grass.py
"""
from sage.all import *

X = SchubertPolynomialRing(ZZ)
s = SymmetricFunctions(QQ).schur()


def grass_perm(lam, k):
    """Grassmannian permutation, descent <= k, shape lam (at most k rows).

    code(w) = (lam_k, ..., lam_1) weakly increasing, so w(i) = lam_{k+1-i} + i
    for i <= k and the unused values follow in increasing order.
    """
    lam = list(lam) + [0] * (k - len(lam))
    assert len(lam) == k and all(lam[i] >= lam[i + 1] for i in range(k - 1))
    head = [lam[k - i] + i for i in range(1, k + 1)]
    used = set(head)
    n = max(head) if head else 1
    rest = [j for j in range(1, n + 1) if j not in used]
    return head + rest


# --- (a)+(b): the indexing conventions --------------------------------------
# Both sides are expanded into ONE common ring: Sage's Schubert `expand()` and
# Sym's `expand()` build their own rings with their own generator counts, so
# comparing them needs an explicit hom, not coercion.
NV = 24
Rbig = PolynomialRing(QQ, NV, 'z')
ZS = [str(g) for g in Rbig.gens()]


def lift(poly):
    P = poly.parent()
    m = P.ngens()
    assert m <= NV, f"{m} generators > {NV}"
    return P.hom([Rbig.gen(i) for i in range(m)], Rbig)(poly)


print("== convention check: S_{w(lam,k)} == s_lam(x_1..x_k) ==", flush=True)
bad_conv, n_conv = [], 0
for k in [1, 2, 3, 4]:
    for size in range(0, 9):
        for lam in Partitions(size, max_length=k):
            w = grass_perm(lam, k)
            got = lift(X(w).expand())
            want = lift(s(lam).expand(k, alphabet=ZS[:k]))
            des = Permutation(w).descents()
            ok_des = des == ([] if max(list(lam) + [0]) == 0 else [k])
            n_conv += 1
            if got != want or not ok_des:
                bad_conv.append((k, list(lam), w, des, str(got)[:50],
                                 str(want)[:50]))
print(f"  checked k=1..4, |lam|<=8 ({n_conv} shapes): "
      f"{len(bad_conv)} mismatches {bad_conv[:3]}", flush=True)

# --- the truncation identity ------------------------------------------------
print("\n== flag truncation: S_u * S_v for same-descent Grassmannians ==",
      flush=True)


def truncation_check(lam, mu, k, verbose=False):
    u, v = grass_perm(lam, k), grass_perm(mu, k)
    lhs = X(u) * X(v)
    prod = s(lam) * s(mu)
    rhs = X.zero()
    dropped = 0
    for nu, c in prod.monomial_coefficients().items():
        if len(nu) > k:
            dropped += 1
            continue
        rhs = rhs + Integer(c) * X(grass_perm(nu, k))
    if verbose:
        print(f"    lam={list(lam)} mu={list(mu)} k={k}: "
              f"{len(prod.monomial_coefficients())} LR terms, "
              f"{dropped} dropped by the >{k}-row flag", flush=True)
    return lhs == rhs, dropped


tot, bad, tot_dropped, nontrivial = 0, [], 0, 0
for k in [1, 2, 3, 4]:
    for a in range(0, 6):
        for b in range(0, 6):
            for lam in Partitions(a, max_length=k):
                for mu in Partitions(b, max_length=k):
                    ok, dropped = truncation_check(lam, mu, k)
                    tot += 1
                    tot_dropped += dropped
                    if dropped:
                        nontrivial += 1
                    if not ok:
                        bad.append((k, list(lam), list(mu)))
print(f"  exhaustive k=1..4, |lam|,|mu| <= 5: {tot} pairs, {len(bad)} FAIL "
      f"{bad[:5]}", flush=True)
print(f"  pairs where the flag actually dropped something: {nontrivial}"
      f"  (total terms dropped: {tot_dropped})", flush=True)

# larger, where the truncation really bites
print("\n== larger cases (where dropping matters most) ==", flush=True)
bad2 = []
for k, lam, mu in [
        (2, [3, 1], [3, 1]), (2, [4, 2], [3, 3]),
        (3, [3, 2, 1], [3, 2, 1]), (3, [4, 2, 1], [3, 3, 2]),
        (4, [4, 3, 2, 1], [4, 3, 2, 1]), (4, [5, 3, 2, 1], [4, 4, 3, 1]),
        (5, [5, 4, 3, 2, 1], [5, 4, 3, 2, 1]),
        (3, [5, 5, 5], [5, 5, 5]), (2, [6, 3], [5, 5]),
]:
    ok, dropped = truncation_check(lam, mu, k, verbose=True)
    if not ok:
        bad2.append((k, lam, mu))
print(f"  {len(bad2)} FAIL {bad2}", flush=True)

# --- how much of the LR expansion the flag throws away ----------------------
# Matters because AutoLr::schur_product (strip_lr.rs:281) has no row bound: as
# specified, leaf dispatch computes the full product and discards the excess.
print("\n== flag waste on the staircase family (k, k-1, ..., 1) ==", flush=True)
print(f"  {'k':>2} {'LR terms':>9} {'kept':>6} {'dropped':>8} {'waste':>7}",
      flush=True)
for k in range(2, 8):
    lam = list(range(k, 0, -1))
    prod = s(lam) * s(lam)
    mc = prod.monomial_coefficients()
    kept = sum(1 for nu in mc if len(nu) <= k)
    print(f"  {k:>2} {len(mc):>9} {kept:>6} {len(mc) - kept:>8} "
          f"{100 * (1 - kept / len(mc)):>6.1f}%", flush=True)

# --- negative control: DIFFERENT descents must NOT obey the rule ------------
print("\n== negative control: different descents ==", flush=True)
u = grass_perm([2, 1], 2)      # descent 2
v = grass_perm([2, 1, 1], 3)   # descent 3
lhs = X(u) * X(v)
allg = all(len(Permutation(w).descents()) <= 1 for w in lhs.monomial_coefficients())
print(f"  u={u} (k=2) * v={v} (k=3): {len(lhs.monomial_coefficients())} terms;"
      f" all outputs Grassmannian? {allg}", flush=True)

print("\nDONE", flush=True)
