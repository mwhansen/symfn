"""`Sym`, a symmetric function tagged with the basis it is written in.

Every method here is a composition of contract-layer calls with bookkeeping
around them — picking a route, clearing a denominator so an integer entry
point can be reached, and putting the basis tag back on the result — so a
`Sym` operation returns the same values as the flat `symfn.*` functions it
is built from. Where a method would need mathematics the contract layer does
not expose, it raises rather than implementing it here.

The tag is the point. No Python type distinguishes `s_λ` from `h_λ` — both are
a partition and a number — so a basis mix-up is a plausible wrong answer rather
than an error. `Sym` carries the basis and
refuses to combine two elements that disagree, which is the crate's "basis
confusion is a compile error" in a language with no compiler.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Sequence
from fractions import Fraction
from typing import TYPE_CHECKING, Callable, Union

from . import symfn as _c
from ._bases import BasisError, check_basis, clear_denominators, exact, restore
from ._types import Basis, Coefficient, Partition, PartitionArg, TermsArg

if TYPE_CHECKING:
    from ._param import Poly

#: What a binary operation accepts beside another element: a scalar is the
#: multiple of the unit, which is the empty partition in every basis.
Operand = Union["Sym", int, Fraction]

__all__ = ["Sym", "s", "h", "e", "p", "m", "f", "skew"]


def _partition(la: PartitionArg) -> Partition:
    """Normalize a partition argument to a zero-free tuple.

        >>> _partition([3, 1, 0])
        (3, 1)
        >>> _partition(4)
        (4,)

    # Raises

    Raises `ValueError` unless the parts are weakly decreasing and positive
    once trailing zeros are dropped. Trailing zeros are tolerated because that
    is the fixed-width form callers hand over.
    """
    if isinstance(la, int):
        la = (la,)
    parts = tuple(int(x) for x in la)
    while parts and parts[-1] == 0:
        parts = parts[:-1]
    if any(x <= 0 for x in parts):
        raise ValueError(f"partition parts must be positive: {tuple(la)!r}")
    if any(a < b for a, b in zip(parts, parts[1:])):
        raise ValueError(f"partition must be weakly decreasing: {tuple(la)!r}")
    return parts


class Sym:
    """A symmetric function, as a basis tag and a zero-free term dictionary.

    Construct one through the basis factories `s`, `h`, `e`, `p`, `m`, `f`
    rather than directly; the constructor is what those call.

        >>> from symfn import s
        >>> s([2, 1])
        s[2,1]
        >>> s([2, 1]) * s([1])
        s[2,1,1] + s[2,2] + s[3,1]

    Coefficients are `int` or `Fraction` and never anything else, so every
    value is exact. Terms are held zero-free and in the contract layer's
    element order — increasing lexicographic by partition — so two equal
    elements have equal `terms`.

    Instances are immutable and hashable. `repr` is the readable form and is
    **not** `eval`-able when a coefficient is rational; `terms` is the
    machine-readable form.

    # Raises

    Raises `BasisError` from `+`, `-`, `*` and comparison when the two operands
    are in different bases; convert one with `to` first.
    """

    __slots__ = ("_basis", "_terms")
    __module__ = "symfn"

    def __init__(self, basis: str, terms: TermsArg) -> None:
        """Build an element from a basis code and a `{partition: coefficient}`
        mapping or a `(partition, coefficient)` sequence.

            >>> Sym("h", {(2,): 1, (1, 1): -1})
            -h[1,1] + h[2]
            >>> Sym("s", [((2,), 3)])
            3*s[2]

        # Raises

        Raises `ValueError` unless `basis` names a basis and every key is a
        partition, and `TypeError` unless every coefficient is an `int` or a
        `Fraction`.
        """
        self._basis = check_basis(basis)
        items = terms.items() if hasattr(terms, "items") else terms
        collected: dict[Partition, Coefficient] = {}
        for la, c in items:
            la = _partition(la)
            c = exact(c) + collected.get(la, 0)
            if c:
                collected[la] = exact(c)
            else:
                collected.pop(la, None)
        self._terms: dict[Partition, Coefficient] = dict(sorted(collected.items()))

    # --- what it is ---------------------------------------------------------

    @property
    def basis(self) -> Basis:
        """The one-letter basis code: one of `s`, `h`, `e`, `p`, `m`, `f`.

        >>> from symfn import p
        >>> p([2]).basis
        'p'
        """
        return self._basis

    @property
    def terms(self) -> dict[Partition, Coefficient]:
        """A copy of the `{partition: coefficient}` mapping, zero-free.

        Partitions are tuples, in increasing lexicographic order.

            >>> from symfn import s
            >>> (s([2]) + 3 * s([1, 1])).terms
            {(1, 1): 3, (2,): 1}
        """
        return dict(self._terms)

    def coefficient(self, la: PartitionArg) -> Coefficient:
        """The coefficient of `la`, or `0` if it does not appear.

        >>> from symfn import s
        >>> (s([2]) * s([1])).coefficient([2, 1])
        1
        >>> s([2]).coefficient([5])
        0
        """
        return self._terms.get(_partition(la), 0)

    def support(self) -> list[Partition]:
        """The partitions carrying a nonzero coefficient, in element order.

        >>> from symfn import s
        >>> (s([2]) * s([1])).support()
        [(2, 1), (3,)]
        """
        return list(self._terms)

    def degree(self) -> int | None:
        """The common degree of every term, or `None` if the element is not
        homogeneous — and `None` for the zero element, which has no degree.

            >>> from symfn import s
            >>> (s([2]) + s([1, 1])).degree()
            2
            >>> (s([2]) + s([1])).degree() is None
            True
        """
        degrees = {sum(la) for la in self._terms}
        return degrees.pop() if len(degrees) == 1 else None

    def is_homogeneous(self) -> bool:
        """Whether every term has the same degree. The zero element is
        homogeneous.

            >>> from symfn import s
            >>> (s([2]) + s([1])).is_homogeneous()
            False
        """
        return len({sum(la) for la in self._terms}) <= 1

    def __len__(self) -> int:
        """The number of nonzero terms.

        >>> from symfn import s
        >>> len(s([2]) * s([1]))
        2
        """
        return len(self._terms)

    def __bool__(self) -> bool:
        """False exactly for the zero element.

        >>> from symfn import s
        >>> bool(s([2]) - s([2]))
        False
        """
        return bool(self._terms)

    def __iter__(self) -> Iterator[tuple[Partition, Coefficient]]:
        """Iterate `(partition, coefficient)` pairs in element order.

        >>> from symfn import s
        >>> list(s([2]) * s([1]))
        [((2, 1), 1), ((3,), 1)]
        """
        return iter(self._terms.items())

    def __eq__(self, other: object) -> bool:
        """Equality within a basis; an element never equals one in another
        basis, and an `int` compares against the multiple of the unit.

            >>> from symfn import s, h
            >>> s([]) * 3 == 3
            True
            >>> s([2]) == h([2])
            False
        """
        if isinstance(other, Sym):
            return self._basis == other._basis and self._terms == other._terms
        if isinstance(other, (int, Fraction)):
            return self._terms == ({(): exact(other)} if other else {})
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self._basis, frozenset(self._terms.items())))

    def __repr__(self) -> str:
        """The readable form: `s[2,1] + 2*s[3]`, with the unit term written as
        a bare coefficient.

            >>> from symfn import s
            >>> s([2, 1]) - 2 * s([3]) + 1
            1 + s[2,1] - 2*s[3]
            >>> s([2]) - s([2])
            0
        """
        if not self._terms:
            return "0"
        pieces = []
        for la, c in self._terms.items():
            negative = c < 0
            mag = -c if negative else c
            atom = f"{self._basis}[{','.join(map(str, la))}]" if la else ""
            if not atom:
                body = _fmt(mag)
            elif mag == 1:
                body = atom
            else:
                body = f"{_fmt(mag)}*{atom}"
            pieces.append(("-" if negative else "+", body))
        sign, body = pieces[0]
        out = ("-" if sign == "-" else "") + body
        for sign, body in pieces[1:]:
            out += f" {sign} {body}"
        return out

    # --- arithmetic ---------------------------------------------------------

    def _same(self, other: object, op: str) -> Sym:
        """Coerce `other` to this basis, or raise. An `int` or `Fraction` is
        the multiple of the unit, which is the empty partition in every basis.
        """
        if isinstance(other, Sym):
            if other._basis != self._basis:
                raise BasisError(
                    f"cannot combine {self._basis} with {other._basis}; "
                    "convert one with .to()"
                )
            return other
        if isinstance(other, (int, Fraction)):
            unit: dict[Partition, Coefficient] = {(): other} if other else {}
            return Sym(self._basis, unit)
        raise TypeError(
            f"cannot {op} {type(other).__name__} with a symmetric function"
        )

    def __add__(self, other: Operand) -> Sym:
        """Add within a basis.

        >>> from symfn import h
        >>> h([2]) + h([2]) + h([1])
        h[1] + 2*h[2]
        """
        other = self._same(other, "add")
        out = dict(self._terms)
        for la, c in other._terms.items():
            out[la] = out.get(la, 0) + c
        return Sym(self._basis, out)

    __radd__ = __add__

    def __neg__(self) -> Sym:
        """Negate every coefficient.

        >>> from symfn import s
        >>> -s([2])
        -s[2]
        """
        return Sym(self._basis, {la: -c for la, c in self._terms.items()})

    def __sub__(self, other: Operand) -> Sym:
        """Subtract within a basis.

        >>> from symfn import s
        >>> s([1]) * s([1]) - s([1, 1])
        s[2]
        """
        return self + (-self._same(other, "subtract"))

    def __rsub__(self, other: Operand) -> Sym:
        return (-self) + other

    def __mul__(self, other: Operand) -> Sym:
        """Multiply, in the basis both operands are written in.

        A scalar scales. Two elements multiply through the contract layer:
        directly in the Schur and monomial bases, which have entry points of
        their own, and through Schur otherwise. The routing is decided here so
        it can be changed in one place.

            >>> from symfn import s, m
            >>> s([1]) * s([1])
            s[1,1] + s[2]
            >>> m([1]) * m([1])
            2*m[1,1] + m[2]
            >>> 2 * s([1])
            2*s[1]

        The Schur values are `c^λ_{μν}` in the Littlewood-Richardson rule; the
        pair above distinguishes it from the Kronecker product, where
        `s_1 * s_1 = s_1`.
        """
        if isinstance(other, (int, Fraction)):
            k = exact(other)
            return Sym(self._basis, {la: c * k for la, c in self._terms.items()})
        other = self._same(other, "multiply")
        if not self._terms or not other._terms:
            return Sym(self._basis, {})
        a, sa = clear_denominators(self._terms)
        b, sb = clear_denominators(other._terms)
        if self._basis == "s":
            return Sym("s", restore(_c.schur_multiply(a, b), sa * sb))
        if self._basis == "m":
            return Sym("m", restore(_c.monomial_multiply(a, b), sa * sb))
        via = self._basis
        sa_, sb_ = _c.convert_terms(a, via, "s"), _c.convert_terms(b, via, "s")
        product = _c.schur_multiply(sa_, sb_)
        return Sym("s", restore(product, sa * sb)).to(self._basis)

    __rmul__ = __mul__

    def __pow__(self, n: int) -> Sym:
        """A non-negative integer power, by repeated squaring.

            >>> from symfn import s
            >>> s([1]) ** 3
            s[1,1,1] + 2*s[2,1] + s[3]

        # Raises

        Raises `ValueError` unless `n` is a non-negative `int`; there is no
        inverse in this ring.
        """
        if not isinstance(n, int) or n < 0:
            raise ValueError(f"exponent must be a non-negative int, not {n!r}")
        one: dict[Partition, Coefficient] = {(): 1}
        out, base = Sym(self._basis, one), self
        while n:
            if n & 1:
                out = out * base
            n >>= 1
            if n:
                base = base * base
        return out

    # --- the ring's own operations -----------------------------------------

    def to(self, basis: str) -> Sym:
        """The element rewritten in `basis`, exactly.

        Conversions into the power-sum basis are rational, so coefficients come
        back as `Fraction`; every other pair lands in ℤ.

            >>> from symfn import s, h
            >>> h([2]).to("s")
            s[2]
            >>> s([1, 1]).to("p")
            1/2*p[1,1] - 1/2*p[2]

        The `p` values pin the convention: `s_11 = (p_1² − p_2)/2`, the second
        elementary symmetric function, which distinguishes it from `s_2`, where
        the sign is `+`.

        # Raises

        Raises `ValueError` unless `basis` names a basis.
        """
        basis = check_basis(basis)
        if basis == self._basis or not self._terms:
            return Sym(basis, self._terms)
        pairs, scale = clear_denominators(self._terms)
        src = self._basis
        if basis == "p":
            rows = [(la, Fraction(n, d)) for la, (n, d) in _c.to_power(pairs, src)]
            return Sym("p", restore(rows, scale))
        return Sym(basis, restore(_c.convert_terms(pairs, src, basis), scale))

    def omega(self) -> Sym:
        """The ω involution, returned in this element's basis.

        ω is its own inverse and exchanges `e` with `h`.

            >>> from symfn import s
            >>> s([2, 1, 1]).omega()
            s[3,1]

        The value pins the convention: ω conjugates the shape, so ω is not the
        identity and not the antipode, which carries a sign.
        """
        return self._through_schur(lambda a: _c.omega(a))

    def antipode(self) -> Sym:
        """The antipode S of the Hopf algebra, in this element's basis.

        `S(s_λ) = (−1)^|λ| s_{λ'}`, which is ω up to that sign — the value
        below is what separates them.

            >>> from symfn import s
            >>> s([2, 1]).antipode()
            -s[2,1]
        """
        return self._through_schur(lambda a: _c.antipode(a))

    def plethysm(self, g: Sym) -> Sym:
        """The plethysm `f[g]`, with `f` this element, in the Schur basis.

            >>> from symfn import s
            >>> s([2]).plethysm(s([2]))
            s[2,2] + s[4]

        The value distinguishes plethysm from the ordinary product, where
        `s_2 · s_2` also carries `s[3,1]`.

        # Raises

        Raises `ValueError` if `g` has a rational coefficient. Plethysm is not
        linear in `g`, so the denominator-clearing every other method here uses
        does not apply, and the contract layer takes integers.
        """
        g = self._same(g, "compose")
        gs = g.to("s")
        if any(isinstance(c, Fraction) for c in gs._terms.values()):
            raise ValueError(
                "plethysm needs integer coefficients in its inner argument"
            )
        pairs, scale = clear_denominators(self.to("s")._terms)
        inner = [(la, int(c)) for la, c in gs._terms.items()]
        return Sym("s", restore(_c.plethysm(pairs, inner), scale)).to(self._basis)

    def internal_product(self, other: Sym) -> Sym:
        """The internal (Kronecker) product, in this element's basis.

            >>> from symfn import s
            >>> s([2, 1]).internal_product(s([2, 1]))
            s[1,1,1] + s[2,1] + s[3]

        The value distinguishes it from the outer product, which is homogeneous
        of degree 6 rather than 3.
        """
        return self._through_schur_pair(other, _c.internal_product)

    def scalar(self, other: Sym) -> Coefficient:
        """The Hall inner product `⟨self, other⟩`, an `int` or `Fraction`.

        The Schur basis is orthonormal for it, which the values pin:

            >>> from symfn import s, p
            >>> s([2, 1]).scalar(s([2, 1])), s([2, 1]).scalar(s([3]))
            (1, 0)
            >>> p([2]).scalar(p([2]))
            2

        `⟨p_λ, p_λ⟩ = z_λ`, so the power-sum basis is orthogonal but not
        orthonormal — the second value is what says which.
        """
        other = self._same(other, "pair")
        a, sa = clear_denominators(self.to("s")._terms)
        b, sb = clear_denominators(other.to("s")._terms)
        return exact(Fraction(_c.hall_inner_product(a, b), sa * sb))

    def skew_by(self, g: Sym) -> Sym:
        """The element skewed by `g`, the adjoint of multiplication by `g`
        under the Hall inner product, in this element's basis.

            >>> from symfn import s
            >>> s([2, 1]).skew_by(s([1]))
            s[1,1] + s[2]

        `s_{λ/μ}` for a partition `μ` is `skew(la, mu)`; this is the general
        form, linear in `g`.
        """
        g = self._same(g, "skew")
        a, sa = clear_denominators(self.to("s")._terms)
        b, sb = clear_denominators(g._terms)
        out = _c.skew_by(a, b, g._basis)
        return Sym("s", restore(out, sa * sb)).to(self._basis)

    def coproduct(self) -> dict[tuple[Partition, Partition], Coefficient]:
        """The coproduct Δ, as a `{(mu, nu): coefficient}` mapping over the
        Schur basis of each factor.

            >>> from symfn import s
            >>> s([2]).coproduct()
            {((), (2,)): 1, ((1,), (1,)): 1, ((2,), ()): 1}

        The keys are partition pairs rather than `Sym` objects, because the
        result lives in a tensor square this type does not model.
        """
        pairs, scale = clear_denominators(self.to("s")._terms)
        rows = _c.coproduct(pairs)
        return {k: exact(Fraction(c, scale)) for k, c in rows}

    # --- leaving the ring ---------------------------------------------------

    def expand(self, n: int) -> dict[tuple[int, ...], Coefficient]:
        """The expansion in `n` variables, as a `{exponent vector: coefficient}`
        mapping with every vector of length `n`.

            >>> from symfn import s
            >>> s([2]).expand(2)
            {(1, 1): 1, (2, 0): 1, (0, 2): 1}

        The order is the contract layer's: grouped by monomial term, and
        unspecified within a group.
        """
        pairs, scale = clear_denominators(self._terms)
        rows = _c.expand_alphabet(pairs, self._basis, n)
        return {v: exact(Fraction(c, scale)) for v, c in rows}

    def evaluate(self, xs: Sequence[int]) -> Coefficient:
        """The value at the alphabet `xs`, a sequence of integers.

            >>> from symfn import s
            >>> s([2, 1]).evaluate([1, 1, 1])
            8

        `s_21(1,1,1) = 8` is the dimension of the `GL_3` irreducible, which
        `principal_specialization(3)` gives without an alphabet.
        """
        pairs, scale = clear_denominators(self.to("s")._terms)
        return exact(Fraction(_c.evaluate_schur(pairs, list(xs)), scale))

    def principal_specialization(self, n: int) -> Coefficient:
        """The value at `1^n`, summed over the Schur expansion.

            >>> from symfn import s
            >>> (s([2, 1]) + s([3])).principal_specialization(3)
            18

        # Raises

        Raises `OverflowError` if a term's value exceeds what the contract
        layer's `principal_specialization` can represent, which reports the
        wall rather than wrapping.
        """
        total: Coefficient = 0
        for la, c in self.to("s")._terms.items():
            v = _c.principal_specialization(la, n)
            if v is None:
                raise OverflowError(
                    f"s_{la} at 1^{n} exceeds the fixed-width specialization"
                )
            total += c * v
        return exact(total)

    def principal_specialization_q(self, n: int) -> Poly:
        """The value at `1, q, …, q^{n−1}`, as a `Poly` in `q`.

            >>> from symfn import s
            >>> s([2, 1]).principal_specialization_q(3)
            q + 2*q^2 + 2*q^3 + 2*q^4 + q^5
            >>> s([2, 1]).principal_specialization_q(3).at(1)
            8

        At `q = 1` this is `principal_specialization`, which is the check that
        the grading is the one that sums to it. The lowest power is `q^{n(λ)}`
        rather than `q^0`, which is what distinguishes this normalization from
        the one that divides the leading power out.

        # Raises

        Raises `ValueError` unless `n` is non-negative.
        """
        from ._param import Poly

        total: dict[int, Coefficient] = {}
        for la, c in self.to("s")._terms.items():
            for k, v in enumerate(_c.principal_specialization_q(la, n)):
                if v:
                    total[k] = total.get(k, 0) + c * v
        return Poly("q", total)

    def dimension(self) -> Coefficient:
        """The dimension `Σ c_λ f^λ`, with `f^λ` the standard-tableaux count.

            >>> from symfn import s
            >>> (s([2, 1]) + s([1, 1, 1])).dimension()
            3

        # Raises

        Raises `OverflowError` if a term's `f^λ` exceeds `u128`.
        """
        total: Coefficient = 0
        for la, c in self.to("s")._terms.items():
            v = _c.dimension(la)
            if v is None:
                raise OverflowError(f"f^{la} exceeds u128")
            total += c * v
        return exact(total)

    # --- routing ------------------------------------------------------------

    def _through_schur(
        self,
        call: Callable[
            [list[tuple[Partition, int]]], Iterable[tuple[Partition, int]]
        ],
    ) -> Sym:
        """Apply a Schur-basis contract call and come back to this basis."""
        pairs, scale = clear_denominators(self.to("s")._terms)
        return Sym("s", restore(call(pairs), scale)).to(self._basis)

    def _through_schur_pair(
        self,
        other: Sym,
        call: Callable[..., Iterable[tuple[Partition, int]]],
    ) -> Sym:
        other = self._same(other, "combine")
        a, sa = clear_denominators(self.to("s")._terms)
        b, sb = clear_denominators(other.to("s")._terms)
        return Sym("s", restore(call(a, b), sa * sb)).to(self._basis)


def _fmt(c: Coefficient) -> str:
    """A coefficient as it appears in a `repr`: `3`, or `1/2` for a rational."""
    return f"{c.numerator}/{c.denominator}" if isinstance(c, Fraction) else str(c)


class _Factory:
    """One basis, callable and subscriptable, producing single terms.

    `s([2, 1])` and `s[2, 1]` are the same element; the subscript reads like
    the `repr` and the call takes a partition a caller already has.

        >>> from symfn import s
        >>> s[2, 1] == s([2, 1])
        True
        >>> s()
        1
    """

    __slots__ = ("_basis",)
    __module__ = "symfn"

    def __init__(self, basis: str) -> None:
        self._basis = check_basis(basis)

    def __call__(self, la: PartitionArg = ()) -> Sym:
        """The basis element indexed by `la`; with no argument, the unit.

        # Raises

        Raises `ValueError` unless `la` is a partition.
        """
        return Sym(self._basis, {_partition(la): 1})

    def __getitem__(self, la: PartitionArg) -> Sym:
        """The basis element indexed by `la`, written as a subscript.

        `s[2, 1]` arrives here as a tuple and `s[2]` as an `int`; both mean the
        partition they read as.
        """
        return Sym(self._basis, {_partition(la): 1})

    def __repr__(self) -> str:
        return self._basis


#: The Schur basis.
s = _Factory("s")
#: The complete homogeneous basis.
h = _Factory("h")
#: The elementary basis.
e = _Factory("e")
#: The power-sum basis.
p = _Factory("p")
#: The monomial basis.
m = _Factory("m")
#: The forgotten basis.
f = _Factory("f")


def skew(la: PartitionArg, mu: PartitionArg) -> Sym:
    """The skew Schur function `s_{λ/μ}`, in the Schur basis.

        >>> from symfn import skew
        >>> skew([2, 1], [1])
        s[1,1] + s[2]
        >>> skew([2, 1], [3])
        0

    A `μ` not contained in `λ` gives zero, which is the theorem rather than a
    convention: the skew shape is empty.

    # Raises

    Raises `ValueError` unless both arguments are partitions.
    """
    return Sym("s", _c.skew_schur(_partition(la), _partition(mu)))
