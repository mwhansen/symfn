"""Check the Schubert Python bindings against Sage (which is Symmetrica).

Separate from the maths: every engine test in src/schubert.rs already passes,
so anything failing here is a *marshalling* fault — a correct answer in the
wrong slot, a 0-vs-1-based index, a permutation that lost its normalization
crossing the FFI. Those are a different failure mode and deserve their own
harness.

  maturin build --release --features python   # then unpack into pybuild/
  sage -python scripts/check_schubert_bindings.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "pybuild"))
import symfn
from sage.all import *

X = SchubertPolynomialRing(ZZ)
P = Permutation
fails, checks = [], 0


def sch(el):
    """Sage Schubert element -> our (one-line, coeff) wire format."""
    return [(list(w), int(c)) for w, c in el.monomial_coefficients().items()]


def norm(terms):
    """Compare ignoring order and trailing fixed points."""
    out = {}
    for w, c in terms:
        w = tuple(w)
        while w and w[-1] == len(w):
            w = w[:-1]
        out[w] = out.get(w, 0) + c
    return {k: v for k, v in out.items() if v}


def check(label, got, want):
    global checks
    checks += 1
    if norm(got) != norm(want):
        fails.append((label, norm(got), norm(want)))


print("== products ==", flush=True)
set_random_seed(7)
cases = [([1, 3, 2], [1, 3, 2]), ([3, 2, 1], [2, 1, 3]), ([2, 4, 1, 3], [2, 4, 1, 3]),
         ([1, 4, 2, 3], [3, 1, 4, 2]), ([2, 4, 6, 1, 3, 5], [2, 4, 6, 1, 3, 5])]
for _ in range(6):
    cases.append((list(Permutations(7).random_element()),
                  list(Permutations(7).random_element())))
for u, v in cases:
    check(f"mul {u}x{v}", symfn.schubert_multiply(sch(X(u)), sch(X(v))), sch(X(u) * X(v)))

print("== multiply_variable (ours 1-based, Sage's 0-based) ==", flush=True)
for w in [[1, 3, 2], [3, 2, 1], [2, 4, 1, 3], [1, 4, 2, 3]]:
    for i in range(1, 5):
        check(f"x_{i}·S_{w}", symfn.schubert_multiply_variable(sch(X(w)), i),
              sch(X(w).multiply_variable(Integer(i - 1))))

print("== divided differences (both 1-based here) ==", flush=True)
for w in [[3, 2, 1], [2, 4, 1, 3], [1, 4, 2, 3], [4, 2, 3, 1]]:
    for i in range(1, 4):
        check(f"d_{i} S_{w}", symfn.schubert_divided_difference(sch(X(w)), i),
              sch(X(w).divided_difference(Integer(i))))

print("== expand / from-polynomial round trip ==", flush=True)
for w in [[1, 3, 2], [3, 2, 1], [1, 4, 2, 3], [2, 4, 1, 3], [3, 1, 4, 2]]:
    e = symfn.schubert_expand(sch(X(w)))
    # against Sage's own expansion
    se = X(w).expand()
    R = se.parent()
    got = R(sum(R.monomial(*(list(v) + [0] * (R.ngens() - len(v)))) * c for v, c in e))
    checks += 1
    if got != se:
        fails.append((f"expand {w}", str(got)[:60], str(se)[:60]))
    check(f"roundtrip {w}", symfn.polynomial_to_schubert(e), sch(X(w)))

print("== dimension == S_w(1,...,1) ==", flush=True)
for w in [[1, 3, 2], [3, 2, 1], [1, 4, 2, 3], [2, 4, 6, 1, 3, 5], [3, 1, 4, 2]]:
    checks += 1
    want = sum(int(c) for c in X(w).expand().coefficients())
    if symfn.schubert_dimension(w) != want:
        fails.append((f"dimension {w}", symfn.schubert_dimension(w), want))

print("== single coefficient vs the full Sage product ==", flush=True)
for u, v in [([1, 3, 2], [1, 3, 2]), ([2, 4, 1, 3], [2, 4, 1, 3]), ([1, 4, 2, 3], [3, 1, 4, 2])]:
    # Sage keys its dict by unpadded permutations, so normalize BOTH sides --
    # looking up [1,4,2,3,5] in a dict keyed [1,4,2,3] silently returns 0 and
    # reads as a binding bug. (It did, for four rows, until this was fixed.)
    full = {}
    for ww, cc in (X(u) * X(v)).monomial_coefficients().items():
        t = tuple(ww)
        while t and t[-1] == len(t):
            t = t[:-1]
        full[t] = int(cc)
    for w in Permutations(5):
        checks += 1
        got = symfn.schubert_coefficient(u, v, list(w))
        t = tuple(w)
        while t and t[-1] == len(t):
            t = t[:-1]
        want = full.get(t, 0)
        if got != want:
            fails.append((f"c^{list(w)}_{u},{v}", got, want))

print("== stanley symmetric function vs Symmetrica's newtrans ==", flush=True)
import sage.libs.symmetrica.all as symca
for w in [[2, 1, 4, 3], [1, 3, 2], [3, 2, 1], [2, 4, 1, 3], [1, 4, 2, 3],
          [3, 1, 4, 2], [4, 2, 3, 1], [2, 4, 6, 1, 3, 5], [5, 3, 1, 4, 2],
          [1, 5, 2, 4, 3], [4, 1, 5, 2, 3]]:
    checks += 1
    got = {tuple(l): c for l, c in symfn.schubert_to_stanley_schur(w)}
    want = {}
    for lam, c in symca.newtrans(P(w)).monomial_coefficients().items():
        want[tuple(lam)] = int(c)
    got = {k: v for k, v in got.items() if v}
    want = {k: v for k, v in want.items() if v}
    if got != want:
        fails.append((f"stanley {w}", got, want))
# the name is a lie in Symmetrica: F_w != S_w outside the stable range
checks += 1
if {tuple(l): c for l, c in symfn.schubert_to_stanley_schur([2, 1, 4, 3])} != {(2,): 1, (1, 1): 1}:
    fails.append(("stanley 2143 hand value", None, None))

print("== stability: padded input is the same element ==", flush=True)
a = symfn.schubert_multiply(sch(X([3, 1, 4, 2])), sch(X([1, 3, 2])))
b = symfn.schubert_multiply([([3, 1, 4, 2, 5, 6], 1)], [([1, 3, 2, 4], 1)])
checks += 1
if norm(a) != norm(b):
    fails.append(("stability under padding", norm(a), norm(b)))

print("== a malformed word raises, and does not abort ==", flush=True)
# A repeated entry is not a permutation. It must come back as a ValueError:
# a PanicException here is a bug report, never an interface
# (docs/policies/failure.md, R2). Both coefficient sizes are exercised on
# purpose -- the small one runs the fixed-width pass, the wide one forces the
# escalation path, and it was the escalation path that used to panic.
for coeff, tag in [(1, "small"), (2 ** 130, "wide")]:
    checks += 1
    try:
        symfn.schubert_expand([([1, 1, 3], coeff)])
        fails.append((f"malformed word ({tag}) did not raise", None, None))
    except ValueError:
        pass
    except BaseException as e:  # PanicException included
        fails.append((f"malformed word ({tag}) raised {type(e).__name__}", None, None))

print(f"\n{checks} checks, {len(fails)} failures", flush=True)
for f in fails[:10]:
    print("  FAIL", f, flush=True)
print("DONE", flush=True)
