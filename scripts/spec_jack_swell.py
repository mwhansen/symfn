"""Pre-implementation swell measurement for spec-jack.md.

Simulates, in Python, the exact arithmetic the Rust implementation would
use for Jack coefficients over Q(alpha): a numerator in Z[alpha] over a
denominator held as a MULTISET OF PRIMITIVE INTEGER-LINEAR ATOMS
(u*alpha + v, gcd(u,v) = 1), plus one rational content.  Cancellation is
detected by the integer root test

    (u*alpha + v) | p  in Q[alpha]   <=>   sum_k p_k (-v)^k u^(d-k) = 0,

exact for linear atoms, so no failed trial division ever runs — unlike the
q,t case, where divide_exact cannot see failure early and needed a
separate necessary-condition pass (spec-macdonald-operators.md 6.2).
Primitivity matters: E-differences are NOT primitive in general
(kappa = (2,2), la = (1,1,1,1) gives 2*alpha + 4), and dividing Z[alpha]
by a non-primitive linear form leaves Z[alpha].  Gauss's lemma makes the
primitive quotient integral, which the synthetic division asserts.

Reports, per degree, for the Laplace-Beltrami recursion at the worst
shape kappa = (n): peak numerator degree, peak numerator coefficient
bits, peak denominator atom count, peak content bits.  Small degrees are
cross-checked against plain rational-function arithmetic at alpha = 5,
so the simulated arithmetic cannot silently be wrong.

This is the measurement that licenses (or kills) the factored design
before any Rust exists — the move spec-macdonald-operators.md 7.1 records.

Run:  sage -python spec_jack_swell.py
(sage only for Partitions; the arithmetic is deliberately hand-rolled.)
"""

import sys
from fractions import Fraction
from math import gcd

from sage.all import Partitions

# ---------------------------------------------------------- Z[alpha]


def pnorm(p):
    while p and p[-1] == 0:
        p.pop()
    return p


def padd(a, b):
    n = max(len(a), len(b))
    return pnorm([(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)])


def pscale(a, c):
    return [x * c for x in a] if c else []


def pmul_atom(a, u, v):
    """p * (u*alpha + v)"""
    return padd(pscale([0] + a, u), pscale(a, v))


def pdiv_atom(p, u, v):
    """Exact division by primitive (u*alpha + v); root test must have passed."""
    q = [0] * (len(p) - 1)
    r = list(p)
    for i in range(len(p) - 1, 0, -1):
        c = r[i]
        assert c % u == 0, "Gauss's lemma violated — atom not primitive?"
        q[i - 1] = c // u
        r[i] = 0
        r[i - 1] -= q[i - 1] * v
    assert not pnorm(r), (p, u, v, r)
    return pnorm(q)


def root_test(p, u, v):
    """sum_k p_k (-v)^k u^(d-k) == 0, i.e. u^d p(-v/u) == 0, in Z."""
    d = len(p) - 1
    acc = 0
    for k, c in enumerate(p):
        acc += c * ((-v) ** k) * (u ** (d - k))
    return acc == 0


class AFrac:
    """content * num(alpha) / prod (u alpha + v)^mult, num in Z[alpha] primitive-ish."""

    def __init__(self, num, den=None, content=Fraction(1)):
        self.num = pnorm(list(num))
        self.den = dict(den or {})
        self.content = content

    @staticmethod
    def zero():
        return AFrac([])

    @staticmethod
    def integer(k):
        return AFrac([k])

    def is_zero(self):
        return not self.num or self.content == 0

    def clone(self):
        return AFrac(list(self.num), dict(self.den), self.content)

    def mul_scalar(self, k):
        out = self.clone()
        out.content *= k
        return out

    def div_atom(self, u, v):
        g = gcd(u, v)
        u, v = u // g, v // g
        out = self.clone()
        out.content /= g
        out.den[(u, v)] = out.den.get((u, v), 0) + 1
        out.reduce_at((u, v))
        return out

    def reduce_at(self, atom):
        u, v = atom
        while self.den.get(atom, 0) > 0 and self.num and root_test(self.num, u, v):
            self.num = pdiv_atom(self.num, u, v)
            self.den[atom] -= 1
            if not self.den[atom]:
                del self.den[atom]
        self.normalize_content()

    def normalize_content(self):
        if not self.num:
            self.content = Fraction(1)
            return
        g = 0
        for c in self.num:
            g = gcd(g, abs(c))
        if g > 1:
            self.num = [c // g for c in self.num]
            self.content *= g

    def add(self, other):
        if self.is_zero():
            return other.clone()
        if other.is_zero():
            return self.clone()
        # distinct primitive linear atoms are pairwise coprime, so
        # atom-wise max multiplicity IS the lcm — exactly, not merely a
        # common multiple (better than the q,t family, which is not coprime).
        den = dict(self.den)
        for a, m in other.den.items():
            den[a] = max(den.get(a, 0), m)

        def lifted(f):
            num = list(f.num)
            for a, m in den.items():
                for _ in range(m - f.den.get(a, 0)):
                    num = pmul_atom(num, a[0], a[1])
            return num

        n1, n2 = lifted(self), lifted(other)
        c1, c2 = self.content, other.content
        q = c1.denominator * c2.denominator // gcd(c1.denominator, c2.denominator)
        num = padd(
            pscale(n1, c1.numerator * (q // c1.denominator)),
            pscale(n2, c2.numerator * (q // c2.denominator)),
        )
        out = AFrac(num, den, Fraction(1, q))
        out.normalize_content()
        for atom in list(out.den):
            out.reduce_at(atom)
        return out

    def eval_at(self, x):
        """Value at alpha = x (Fraction), for the cross-check."""
        num = Fraction(0)
        for c in reversed(self.num):
            num = num * x + c
        den = Fraction(1)
        for (u, v), m in self.den.items():
            den *= (u * x + v) ** m
        return self.content * num / den

    def stats(self):
        deg = len(self.num) - 1 if self.num else 0
        bits = max((abs(c).bit_length() for c in self.num), default=0)
        cbits = max(self.content.numerator.bit_length(), self.content.denominator.bit_length())
        return deg, bits, sum(self.den.values()), cbits


# ------------------------------------------------ the LB recursion


def nfn(la):
    return sum(i * x for i, x in enumerate(la))


def conj(la):
    if not la:
        return []
    return [sum(1 for x in la if x > j) for j in range(la[0])]


def dominates(la, mu):
    ca = cb = 0
    for i in range(max(len(la), len(mu))):
        ca += la[i] if i < len(la) else 0
        cb += mu[i] if i < len(mu) else 0
        if ca < cb:
            return False
    return True


def lb_row(kappa, coeff_zero, coeff_one, add, scale, divide):
    """The [MOPS] recursion, generic over the coefficient arithmetic."""
    kappa = [int(x) for x in kappa]
    n = sum(kappa)
    Ekp, Ekm = nfn(conj(kappa)), nfn(kappa)
    coeffs = {tuple(kappa): coeff_one()}
    peak = (0, 0, 0, 0)
    order = sorted(Partitions(n), key=lambda x: nfn(list(x)))
    for la in order:
        la_l = [int(x) for x in la]
        if tuple(la_l) == tuple(kappa) or not dominates(kappa, la_l):
            continue
        acc = coeff_zero()
        L = len(la_l)
        for i in range(L):
            for j in range(i + 1, L):
                for t in range(1, la_l[j] + 1):
                    mu = la_l[:]
                    mu[i] += t
                    mu[j] -= t
                    mu = tuple(sorted((x for x in mu if x), reverse=True))
                    if not dominates(kappa, list(mu)):
                        continue
                    c = coeffs.get(mu)
                    if c is None:
                        continue
                    acc = add(acc, scale(c, la_l[i] - la_l[j] + 2 * t))
                    if hasattr(acc, "stats"):
                        peak = tuple(max(a, b) for a, b in zip(peak, acc.stats()))
        u = Ekp - nfn(conj(la_l))
        v = nfn(la_l) - Ekm
        assert u > 0 and v > 0, (kappa, la_l, u, v)
        if not (hasattr(acc, "is_zero") and acc.is_zero() or acc == 0):
            coeffs[tuple(la_l)] = divide(acc, u, v)
    return coeffs, peak


def afrac_row(kappa):
    return lb_row(
        kappa,
        AFrac.zero,
        lambda: AFrac.integer(1),
        lambda a, b: a.add(b),
        lambda c, k: c.mul_scalar(k),
        lambda c, u, v: c.div_atom(u, v),
    )


def rational_row(kappa, x):
    """Same recursion over plain Fractions at alpha = x — the cross-check."""
    return lb_row(
        kappa,
        lambda: Fraction(0),
        lambda: Fraction(1),
        lambda a, b: a + b,
        lambda c, k: c * k,
        lambda c, u, v: c / (u * x + v),
    )


def main():
    x = Fraction(5)
    for n in range(3, 8):
        for kappa in Partitions(n):
            got, _ = afrac_row(list(kappa))
            want, _ = rational_row(list(kappa), x)
            assert set(got) == set(want), kappa
            for la, c in got.items():
                assert c.eval_at(x) == want[la], (kappa, la)
    print("PASS  AFrac arithmetic == plain rationals at alpha = 5, all kappa, n <= 7")
    print()
    print("# LB recursion at kappa = (n), factored-atom simulation, reduce after every add")
    print(f"{'n':>3} {'p(n)':>5} {'num deg':>8} {'num bits':>9} {'den atoms':>10} {'content bits':>13}")
    for n in range(4, 13):
        _, peak = afrac_row([n])
        print(f"{n:>3} {Partitions(n).cardinality():>5} {peak[0]:>8} {peak[1]:>9} {peak[2]:>10} {peak[3]:>13}")
        sys.stdout.flush()
    print()
    print("# worst peak over ALL kappa |- n (is (n) really the worst shape?)")
    print(f"{'n':>3} {'num deg':>8} {'num bits':>9} {'den atoms':>10} {'content bits':>13}  worst kappa per column")
    for n in range(6, 12):
        worst = [(0, None)] * 4
        for kappa in Partitions(n):
            _, peak = afrac_row(list(kappa))
            for i in range(4):
                if peak[i] > worst[i][0]:
                    worst[i] = (peak[i], list(kappa))
        cols = "  ".join(str(w[1]) for w in worst)
        print(f"{n:>3} {worst[0][0]:>8} {worst[1][0]:>9} {worst[2][0]:>10} {worst[3][0]:>13}  {cols}")
        sys.stdout.flush()


if __name__ == "__main__":
    main()
