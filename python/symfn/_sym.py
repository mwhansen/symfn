"""`Sym`, a symmetric function tagged with the basis it is written in.

**One class for every element**, whether or not its coefficients carry
parameters. The difference between `s([2])` and `macdonald.P([2])` is what a
coefficient is and which basis the tag names, not what kind of thing the
element is. Each method below picks its route from `parameters`: empty is the
classical one, straight to the integer entry points; anything else goes through
the per-ring entry points in `_families`.

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
from typing import TYPE_CHECKING, Callable, Union, cast

from . import symfn as _c
from ._bases import (
    BASES,
    BasisError,
    check_basis,
    check_param_basis,
    clear_denominators,
    exact,
    restore,
)
from ._types import (
    AnyCoefficient,
    Coefficient,
    ParamBasis,
    Partition,
    PartitionArg,
    TermsArg,
)

if TYPE_CHECKING:
    from ._families import Scalar
    from ._param import ParamCoefficient, Poly, QtPoly

#: What a binary operation accepts beside another element: a scalar is the
#: multiple of the unit, which is the empty partition in every basis.
Operand = Union["Sym", int, Fraction]

__all__ = ["Sym", "s", "h", "e", "p", "m", "f", "skew"]


def _scalars() -> tuple[type, ...]:
    """What may stand in for an element in the parametric arithmetic: a scalar
    of the base ring, or an integer or rational that injects into it.

    A function rather than a constant because the coefficient classes live in
    `_param`, which imports this module; naming them at import time would close
    the cycle.
    """
    from ._param import _SCALARS

    return _SCALARS


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


if TYPE_CHECKING:  # the parameter types, for annotations only
    from ._param import Poly, QtPoly


def _check_coefficient(c: object, params: tuple[str, ...]) -> None:
    """Refuse a parameter-carrying coefficient that does not fit `params`.

    A mismatch here used to construct and die later — `Sym("m", {(2,): q})`
    surfaced as a `'<' not supported` comparison from inside `repr` — so the
    constructor is where the promise is checked (P7: misuse raises rather
    than computing garbage).

    # Raises

    Raises `TypeError` for a coefficient of no known class or one whose
    parameters the element does not declare, and `ValueError` when the
    class's variables are not the ones `params` names.
    """
    from ._param import AlphaFrac, Poly, QtFrac, QtPoly, QtRatio

    if not isinstance(c, (Poly, QtPoly, QtFrac, QtRatio, AlphaFrac)):
        raise TypeError(
            "a coefficient must be an int, a Fraction, or one of the "
            f"parameter-carrying coefficient classes, not {type(c).__name__}"
        )
    if not params:
        raise TypeError(
            f"a {type(c).__name__} coefficient carries parameters the "
            "element does not declare; name them in the 'parameters' "
            "argument"
        )
    if isinstance(c, Poly):
        if params != (c.variable,):
            raise ValueError(
                f"a Poly coefficient in {c.variable!r} does not match "
                f"parameters {params}"
            )
    elif isinstance(c, AlphaFrac):
        if params != ("alpha",):
            raise ValueError(
                "an AlphaFrac coefficient is in alpha, which does not match "
                f"parameters {params}"
            )
    elif params != ("q", "t"):
        raise ValueError(
            f"a {type(c).__name__} coefficient is in q and t, which does "
            f"not match parameters {params}"
        )


class Sym:
    """A symmetric function, as a basis tag and a zero-free term dictionary.

    Construct one through the basis factories `s`, `h`, `e`, `p`, `m`, `f`, or
    through a family — `macdonald`, `jack`, `hl`, `llt` — rather than directly;
    the constructor is what those call.

        >>> from symfn import s, jack
        >>> s([2, 1])
        s[2,1]
        >>> s([2, 1]) * s([1])
        s[2,1,1] + s[2,2] + s[3,1]
        >>> jack.P([2]).to("m")
        2/(alpha + 1)*m[1,1] + m[2]

    **There is one kind of element here.** A coefficient is an `int` or a
    `Fraction`, or one of the classes that hold a parameter — `Poly`,
    `QtPoly`, `QtFrac`, `QtRatio`, `AlphaFrac` — and `parameters` names which,
    empty for the first two. The basis tag is one of the six classical codes or
    one of the nine parametric ones.

    Every value is exact whichever coefficient it carries. Terms are held
    zero-free and in the contract layer's element order — increasing
    lexicographic by partition — so two equal elements have equal `terms`.

    Instances are immutable and hashable. `repr` is the readable form and is
    **not** `eval`-able when a coefficient is rational; `terms` is the
    machine-readable form.

    # Raises

    Raises `BasisError` from `+`, `-` and `*` when the two operands are in
    different bases; convert one with `to` first. Raises `BaseRingError` when
    they are over different base rings, which `to` cannot fix. Comparison
    never raises: elements the operations refuse to combine are unequal.
    """

    __slots__ = ("_basis", "_terms", "_params")
    __module__ = "symfn"

    def __init__(
        self,
        basis: str,
        terms: TermsArg,
        parameters: Sequence[str] = (),
    ) -> None:
        """Build an element from a basis code, a `{partition: coefficient}`
        mapping or a `(partition, coefficient)` sequence, and the parameter
        names `at` accepts.

            >>> Sym("h", {(2,): 1, (1, 1): -1})
            -h[1,1] + h[2]
            >>> Sym("s", [((2,), 3)])
            3*s[2]

        `parameters` is empty for an element whose coefficients are numbers,
        and naming it is what tells the operations which route to take. The
        families supply it; a caller building an element by hand rarely does.

        # Raises

        Raises `ValueError` unless `basis` names a basis and every key is a
        partition, and when a repeated shape carries a parameter-carrying
        coefficient — numbers accumulate, but summing two ring coefficients
        is the ring's business, and refusing beats dropping one. Raises
        `TypeError` unless every coefficient is an `int`, a `Fraction`, or
        one of the parameter-carrying coefficient classes, and unless the
        coefficients fit `parameters`: numbers when it is empty, one
        coefficient class in the named parameters when it is not.
        """
        self._basis = check_param_basis(basis)
        self._params = tuple(parameters)
        items = terms.items() if hasattr(terms, "items") else terms
        collected: dict[Partition, AnyCoefficient] = {}
        kind: type | None = None
        for la, c in items:
            la = _partition(la)
            if isinstance(c, (int, Fraction)):
                if self._params:
                    raise TypeError(
                        f"an element with parameters {self._params} takes "
                        "coefficients in one parameter-carrying class, not "
                        f"{type(c).__name__}"
                    )
                # Numbers accumulate, so a repeated shape sums rather than the
                # last one winning.
                c = exact(exact(c) + cast("Coefficient", collected.get(la, 0)))
            else:
                _check_coefficient(c, self._params)
                if kind is None:
                    kind = type(c)
                elif type(c) is not kind:
                    raise TypeError(
                        f"coefficients mix {kind.__name__} and "
                        f"{type(c).__name__}; an element holds one "
                        "coefficient class"
                    )
                if la in collected:
                    # A parameter-carrying coefficient cannot accumulate:
                    # adding two of those is the ring's business, and every
                    # caller that builds one has already done it.
                    raise ValueError(
                        f"the shape {la} is repeated and its coefficient "
                        "carries a parameter; sum the coefficients before "
                        "constructing"
                    )
            if c:
                # `_check_coefficient` has vouched for the class; the cast
                # only makes that legible to the checker.
                collected[la] = cast("AnyCoefficient", c)
            else:
                collected.pop(la, None)
        self._terms: dict[Partition, AnyCoefficient] = dict(sorted(collected.items()))

    # --- what it is ---------------------------------------------------------

    @property
    def basis(self) -> ParamBasis:
        """The basis the terms are indexed by: one of the six classical codes
        `s`, `h`, `e`, `p`, `m`, `f`, or one of the nine parametric tags for an
        element written in a family's own basis.

        >>> from symfn import hl, macdonald, p, s
        >>> p([2]).basis
        'p'
        >>> macdonald.P([2]).basis, hl.Qp([2]).basis, hl.to_Qp(s([2])).basis
        ('McdP', 'HLQp', 'HLQp')
        >>> macdonald.P([2]).to("m").basis
        'm'
        """
        return self._basis

    @property
    def parameters(self) -> tuple[str, ...]:
        """The parameter names `at` takes, in order — empty for an element
        whose coefficients are numbers.

        >>> from symfn import jack, macdonald, s
        >>> macdonald.P([2]).parameters, jack.P([2]).parameters
        (('q', 't'), ('alpha',))
        >>> s([2]).parameters
        ()

        This is what every operation below reads to pick its route, so it is
        the one place the two kinds of element are still told apart.
        """
        return self._params

    @property
    def terms(self) -> dict[Partition, AnyCoefficient]:
        """A copy of the `{partition: coefficient}` mapping, zero-free.

        Partitions are tuples, in increasing lexicographic order.

            >>> from symfn import hl, s
            >>> (s([2]) + 3 * s([1, 1])).terms
            {(1, 1): 3, (2,): 1}
            >>> hl.Qp([2]).to("s").terms
            {(2,): 1}
        """
        return dict(self._terms)

    def coefficient(self, la: PartitionArg) -> AnyCoefficient:
        """The coefficient of `la`, or `0` if it does not appear.

        >>> from symfn import macdonald, s
        >>> (s([2]) * s([1])).coefficient([2, 1])
        1
        >>> s([2]).coefficient([5])
        0
        >>> macdonald.P([2]).to("m").coefficient([5])
        0
        """
        return self._terms.get(_partition(la), 0)

    def support(self) -> list[Partition]:
        """The partitions carrying a nonzero coefficient, in element order.

        >>> from symfn import jack, s
        >>> (s([2]) * s([1])).support()
        [(2, 1), (3,)]
        >>> jack.P([2]).to("m").support()
        [(1, 1), (2,)]
        """
        return list(self._terms)

    def degree(self) -> int | None:
        """The common degree of every term, or `None` if the element is not
        homogeneous — and `None` for the zero element, which has no degree.

            >>> from symfn import m, macdonald, q, s
            >>> (s([2]) + s([1, 1])).degree()
            2
            >>> (s([2]) + s([1])).degree() is None
            True
            >>> macdonald.P([2]).degree()
            2

        The degree is a fact about the partitions alone, so it needs no
        expansion out of a parametric basis and no arithmetic on the
        coefficients.
        """
        degrees = {sum(la) for la in self._terms}
        return degrees.pop() if len(degrees) == 1 else None

    def is_homogeneous(self) -> bool:
        """Whether every term has the same degree. The zero element is
        homogeneous.

            >>> from symfn import hl, s
            >>> (s([2]) + s([1])).is_homogeneous()
            False
            >>> (hl.Qp([2]) + hl.Qp([1, 1])).is_homogeneous()
            True
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

    def __iter__(self) -> Iterator[tuple[Partition, AnyCoefficient]]:
        """Iterate `(partition, coefficient)` pairs in element order.

        >>> from symfn import s
        >>> list(s([2]) * s([1]))
        [((2, 1), 1), ((3,), 1)]
        """
        return iter(self._terms.items())

    def __eq__(self, other: object) -> bool:
        """Equality within a basis, except that zero and the constants are
        equal wherever they are written; an `int` or a `Fraction` compares
        against the multiple of the unit.

            >>> from symfn import h, hl, jack, m, macdonald, q, s
            >>> q * m([2]) - q * m([2]) == 0
            True
            >>> hl.P([1]) ** 0 == 1
            True
            >>> jack.P([]) == 1
            True
            >>> s([]) * 3 == h([]) * 3
            True
            >>> s([2]) == h([2])
            False
            >>> macdonald.P([2]) == s([2])
            False

        An element with no terms is the zero of every ring here, and an
        element supported on the empty partition alone is the constant its
        one coefficient names — the empty partition indexes 1 in all fifteen
        bases, so a constant is the same element wherever it is written.
        Everything else compares within its basis and base ring: the last
        two values are `False` on those grounds, never an error.
        """
        if isinstance(other, Sym):
            if not self._terms or not other._terms:
                return not self._terms and not other._terms
            mine, theirs = self._terms, other._terms
            if len(mine) == 1 and () in mine and len(theirs) == 1 and () in theirs:
                a, b = mine[()], theirs[()]
                if isinstance(a, (int, Fraction)):
                    return bool(a == b)
                if isinstance(b, (int, Fraction)):
                    return bool(b == a)
                # Coefficient classes answer `==` against numbers but not
                # against each other — and a constant `Poly` compares its
                # variable, which a constant does not depend on — so two
                # constants meet on the numbers they equal.
                from ._param import _constant_value

                av, bv = _constant_value(a), _constant_value(b)
                if av is not None and bv is not None:
                    return av == bv
                return bool(a == b)
            return (
                self._basis == other._basis
                and self._params == other._params
                and self._terms == other._terms
            )
        if isinstance(other, (int, Fraction)):
            if not self._terms:
                return not other
            if len(self._terms) == 1 and () in self._terms:
                return bool(self._terms[()] == exact(other))
            return False
        return NotImplemented

    def __hash__(self) -> int:
        """Hash follows `==`: the empty element hashes as `0` and a constant
        as its one coefficient, which itself hashes as the number it equals.
        """
        if not self._terms:
            return hash(0)
        if len(self._terms) == 1 and () in self._terms:
            return hash(self._terms[()])
        return hash((self._basis, self._params, frozenset(self._terms.items())))

    def __repr__(self) -> str:
        """The readable form: `s[2,1] + 2*s[3]`, with the unit term written as
        a bare coefficient.

            >>> from symfn import jack, s
            >>> s([2, 1]) - 2 * s([3]) + 1
            1 + s[2,1] - 2*s[3]
            >>> s([2]) - s([2])
            0
            >>> jack.P([2]).to("m")
            2/(alpha + 1)*m[1,1] + m[2]

        A parameter-carrying coefficient is rendered by its own class and
        parenthesized where it adds at the top level, so the sum stays
        readable as a sum.
        """
        if not self._terms:
            return "0"
        if self._params:
            return self._param_repr()
        pieces = []
        for la, c in self._numbers().items():
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

    def _param_repr(self) -> str:
        """`__repr__` for the coefficient classes that carry a parameter.

        Separate because the classical path formats numbers itself — a
        `Fraction` prints as `1/2` and not as `Fraction(1, 2)` — while these
        classes render themselves and only need a sign lifted out and, where
        they add at the top level, a pair of parentheses.
        """
        from ._param import _factor

        pieces = []
        for la, c in self._terms.items():
            atom = f"{self._basis}[{','.join(map(str, la))}]" if la else ""
            sign, body = _factor(c)
            if not atom:
                pass
            elif body == "1":
                body = atom
            else:
                body = f"{body}*{atom}"
            pieces.append((sign, body))
        sign, body = pieces[0]
        out = ("-" if sign == "-" else "") + body
        for sign, body in pieces[1:]:
            out += f" {sign} {body}"
        return out

    def at(self, *args: Coefficient, **kwargs: Coefficient) -> Sym:
        """The element with its parameters set, in the same basis.

            >>> from symfn import hl, macdonald
            >>> macdonald.P([2]).at(q=7, t=7)
            m[1,1] + m[2]
            >>> hl.Qp([2, 1]).at(t=0)
            s[2,1]

        The first is `P_λ(x; q, q) = s_λ` written in the monomial basis; the
        second is `Q'_λ(x; 0) = s_λ`. Both are theorems, and both fail under a
        `q ↔ t` or `t → 1/t` twist of the convention.

        The classical bases are the only ones a parameter-free element can be
        written in, so an element in a parametric basis is expanded first —
        through `to`, into the basis its family is written in — and specialized
        there.

        # Raises

        Raises `TypeError` unless exactly this element's `parameters` are
        supplied, by name or in that order, and for an element that has none.
        Raises `ValueError` on what `to` raises for, when the element is in a
        parametric basis.
        """
        from ._families import EXPANDS_IN

        if self._basis not in BASES:
            return self.to(EXPANDS_IN[self._basis]).at(*args, **kwargs)
        if not self._params:
            raise TypeError(
                "this element carries no parameters, so there is nothing to "
                "set; it is already the value at every point"
            )
        if args and kwargs:
            raise TypeError("give the parameters by name or by position, not both")
        if args:
            if len(args) != len(self._params):
                raise TypeError(
                    f"{len(self._params)} parameters expected "
                    f"{self._params}, got {len(args)}"
                )
            values = dict(zip(self._params, args))
        else:
            if set(kwargs) != set(self._params):
                raise TypeError(
                    f"parameters {self._params} expected, got {tuple(kwargs)}"
                )
            values = kwargs
        ordered = [values[name] for name in self._params]
        # `parameters` being non-empty is the promise that every coefficient
        # is one of the five classes, each of which has `at`.
        return Sym(
            self._basis,
            {
                la: cast("ParamCoefficient", c).at(*ordered)
                for la, c in self._terms.items()
            },
        )

    # --- arithmetic ---------------------------------------------------------

    def _numbers(self) -> dict[Partition, Coefficient]:
        """The terms, typed as numbers.

        Every classical route below runs only when `parameters` is empty, which
        is the promise that each coefficient is an `int` or a `Fraction`. The
        cast makes that promise legible to a checker rather than widening
        anything: no value changes hands.
        """
        return cast("dict[Partition, Coefficient]", self._terms)

    def _needs_ring(self, other: object) -> bool:
        """Whether an operation with `other` has to run over a base ring that
        carries parameters, rather than over the integers.

        Either operand can put it there: `s([2]) + q * s([2])` is an element of
        `ℚ[q,t]` even though the first term is not, and ℚ sits inside every
        base ring here, so the parameter-free side is lifted rather than
        refused.

        A bare parameter is not one of those operands. `q * m([2])` is a
        parameter-free element scaled by a coefficient, which `_lift` already
        answers by moving the value to the class that can hold it — routing it
        here instead would ask the ring machinery to scale integers.
        """
        if self._params:
            return True
        return isinstance(other, Sym) and bool(other._params)

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

            >>> from symfn import h, jack, macdonald, q, s, t
            >>> h([2]) + h([2]) + h([1])
            h[1] + 2*h[2]
            >>> q * macdonald.Htilde([2, 1]) + t * macdonald.Htilde([3])
            q*McdHt[2,1] + t*McdHt[3]
            >>> s([2]) + q * s([2])
            (1 + q)*s[2]

        Two different bases raise rather than one being converted:
        `McdP[2] + McdQ[2]` names no element. Two different base rings raise as
        well, and say so separately, because `.to()` cannot fix that one — but
        ℚ sits inside every base ring here, so an operand without parameters is
        **lifted** rather than refused, which is the third value.

        **A scalar adds**, as the constant it names times the unit, which the
        empty partition indexes in every basis:

            >>> 2 + jack.P([1])
            2 + JackP[1]
            >>> sum([jack.P([1]), jack.P([2])])
            JackP[1] + JackP[2]

        The second is why: `sum` starts from `0`, so without this it raises on
        the first term.

        # Raises

        Raises `BasisError` unless both elements are in the same basis, and
        `BaseRingError` unless both are over the same base ring. Raises
        `ValueError` for `McdHt` alone if two coefficients at one shape have
        different denominators, which the boundary encoding cannot put over a
        common one.
        """
        if self._needs_ring(other):
            from ._families import _add, _constant, _ring_pair

            if isinstance(other, Sym):
                return _add(*_ring_pair(self, other, "add"))
            if isinstance(other, _scalars()):
                return _add(self, _constant(self, other))
            return NotImplemented
        other = self._same(other, "add")
        out = self._numbers()
        for la, c in other._numbers().items():
            out[la] = out.get(la, 0) + c
        return Sym(self._basis, out)

    __radd__ = __add__

    def __neg__(self) -> Sym:
        """Negate every coefficient.

        >>> from symfn import s
        >>> -s([2])
        -s[2]
        """
        if self._params:
            return self * -1
        return Sym(self._basis, {la: -c for la, c in self._numbers().items()})

    def __sub__(self, other: Operand) -> Sym:
        """Subtract within a basis.

        >>> from symfn import s
        >>> s([1]) * s([1]) - s([1, 1])
        s[2]
        """
        if self._needs_ring(other):
            from ._families import _add, _constant, _ring_pair

            if isinstance(other, Sym):
                f, g = _ring_pair(self, other, "subtract")
                return _add(f, -g)
            if isinstance(other, _scalars()):
                return _add(self, -_constant(self, other))
            return NotImplemented
        return self + (-self._same(other, "subtract"))

    def __rsub__(self, other: Operand) -> Sym:
        if self._needs_ring(other):
            from ._families import _add, _constant

            if isinstance(other, _scalars()):
                return _add(-self, _constant(self, cast("Scalar", other)))
            return NotImplemented
        return (-self) + other

    def __mul__(self, other: Operand | ParamCoefficient) -> Sym:
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

        A **parameter** scales too, and the element then carries it:

            >>> from symfn import alpha, jack, q, macdonald
            >>> q * m([2])
            q*m[2]
            >>> macdonald.to_P(q * m([2])).coefficient([2])
            q
            >>> (1 - alpha) * jack.P([2])
            (1 - alpha)*JackP[2]

        which is how a scaled classical element reaches the inverse
        expansions. A scalar is an `int`, a `Fraction`, a polynomial in this
        element's parameters, or a coefficient of the kind it carries.

        **A parametric basis multiplies like any other.** The product is taken
        where the family expands and rewritten back by the inverse expansion,
        so all nine tags have it:

            >>> from symfn import hl
            >>> hl.P([1]) * hl.P([1])
            (1 + t)*HLP[1,1] + HLP[2]
            >>> jack.P([1]) * jack.P([1])
            2*alpha/(alpha + 1)*JackP[1,1] + JackP[2]
            >>> (hl.P([2, 1]) * hl.P([2, 1])).coefficient([3, 1, 1, 1])
            1 + t - t^3 - t^4

        The last value is what pins the Hall-Littlewood normalization, and the
        first is not: `P_1² = P_2 + (1 + t)·P_11` is what every convention in
        circulation gives. It has **negative** coefficients, so these structure
        constants are in ℤ[t] rather than the ℕ[t] of the classical Hall
        polynomials counting subgroups of abelian p-groups, which differ from
        these by a twist. Two independent readings confirm it: the constant
        term is the Littlewood-Richardson coefficient `c^{3111}_{21,21} = 1`,
        since `P_λ(x; 0) = s_λ`, and the value at `t = 1` is 0, matching
        `m[2,1]²` having no `m[3,1,1,1]` term, since `P_λ(x; 1) = m_λ`.

        # Raises

        Raises `TypeError` if the scalar is in the wrong parameters,
        `BasisError` if the two elements are in different bases, and
        `BaseRingError` if they are over different base rings.
        """
        if self._needs_ring(other):
            from ._families import _product, _ring_pair, _scale

            if isinstance(other, Sym):
                return _product(*_ring_pair(self, other, "multiply"))
            if not isinstance(other, _scalars()):
                return NotImplemented
            return _scale(self, other)
        from ._param import AlphaFrac as _AlphaFrac
        from ._param import Poly as _Poly
        from ._param import QtFrac as _QtFrac
        from ._param import QtPoly as _QtPoly
        from ._param import QtRatio as _QtRatio

        if isinstance(other, (_Poly, _QtPoly)):
            return self._lift(other)
        if isinstance(other, (_QtFrac, _QtRatio, _AlphaFrac)):
            # A fraction coefficient names its ring, so the element is lifted
            # into it and scaled there — the same move `_ring_pair` makes when
            # one operand of a binary operation carries no parameters.
            from ._families import _COEFF_RING, _lift_to, _scale

            ring = _COEFF_RING[type(other)]
            if not other:
                return Sym(self._basis, {}, ring)
            like = Sym(self._basis, [((), other)], ring)
            return _scale(_lift_to(self, like, "multiply"), other)
        return self._times(other)

    __rmul__ = __mul__

    def __truediv__(self, other: object) -> Sym:
        """Division by a nonzero `int` or `Fraction` scalar.

            >>> from fractions import Fraction
            >>> from symfn import hl, s
            >>> s([2]) / 2
            1/2*s[2]
            >>> s([2]) / 2 == Fraction(1, 2) * s([2])
            True
            >>> (2 * hl.P([2])) / 2
            HLP[2]

        Only a scalar divisor is taken: nothing in this ring is invertible
        but the scalars, so dividing by an element names no operation.

        # Raises

        Raises `ZeroDivisionError` on zero, and `TypeError` for a divisor
        that is not an `int` or a `Fraction`.
        """
        if not isinstance(other, (int, Fraction)):
            raise TypeError(
                "cannot divide a symmetric function by "
                f"{type(other).__name__}; only an int or Fraction divisor "
                "is taken"
            )
        if not other:
            raise ZeroDivisionError("cannot divide an element by zero")
        return self * (Fraction(1) / other)

    def _times(self, other: Operand) -> Sym:
        """`__mul__` with the parameter case already handled, so the result is
        a `Sym`.

        `__pow__` composes products and calls this rather than `*`: a
        parameter has no place in the middle of one, and the narrower return
        is what keeps that loop typed.
        """
        if isinstance(other, (int, Fraction)):
            k = exact(other)
            return Sym(self._basis, {la: c * k for la, c in self._numbers().items()})
        other = self._same(other, "multiply")
        if not self._terms or not other._terms:
            return Sym(self._basis, {})
        a, sa = clear_denominators(self._numbers())
        b, sb = clear_denominators(other._numbers())
        if self._basis == "s":
            return Sym("s", restore(_c.schur_multiply(a, b), sa * sb))
        if self._basis == "m":
            return Sym("m", restore(_c.monomial_multiply(a, b), sa * sb))
        via = self._basis
        sa_, sb_ = _c.convert_terms(a, via, "s"), _c.convert_terms(b, via, "s")
        product = _c.schur_multiply(sa_, sb_)
        return Sym("s", restore(product, sa * sb)).to(self._basis)

    def _lift(self, c: Poly | QtPoly) -> Sym:
        """The element scaled by a parameter, in the same basis.

        The coefficients move from numbers to the class that can hold the
        parameter, and `parameters` names it, so every operation takes the ring
        route afterwards. That is what makes `q * m([2])` an argument the
        inverse expansions accept.
        """
        from ._param import QtPoly

        params = ("q", "t") if isinstance(c, QtPoly) else (c.variable,)
        terms: dict[Partition, ParamCoefficient] = {
            la: c * v for la, v in self._numbers().items() if v
        }
        return Sym(self._basis, terms, params)

    def __pow__(self, n: int) -> Sym:
        """A non-negative integer power, by repeated squaring.

            >>> from symfn import hl, s
            >>> s([1]) ** 3
            s[1,1,1] + 2*s[2,1] + s[3]
            >>> hl.P([1]) ** 2
            (1 + t)*HLP[1,1] + HLP[2]
            >>> hl.P([1]) ** 0
            1

        # Raises

        Raises `ValueError` unless `n` is a non-negative `int`; there is no
        inverse in this ring. Raises what `*` raises for the coefficient
        classes it does not carry.
        """
        if self._params:
            from ._families import _product, _unit_like

            if not isinstance(n, int) or n < 0:
                raise ValueError(f"exponent must be a non-negative int, not {n!r}")
            if n == 0:
                return _unit_like(self)
            out, base = self, self
            n -= 1
            while n:
                if n & 1:
                    out = _product(out, base)
                n >>= 1
                if n:
                    base = _product(base, base)
            return out
        if not isinstance(n, int) or n < 0:
            raise ValueError(f"exponent must be a non-negative int, not {n!r}")
        one: dict[Partition, Coefficient] = {(): 1}
        out, base = Sym(self._basis, one), self
        while n:
            if n & 1:
                out = out._times(base)
            n >>= 1
            if n:
                base = base._times(base)
        return out

    # --- the ring's own operations -----------------------------------------

    def to(self, basis: str) -> Sym:
        """The element rewritten in `basis`, any of the fifteen codes.

        Conversions into the power-sum basis are rational, so coefficients come
        back as `Fraction`; every other classical pair lands in ℤ.

            >>> from symfn import s, h, hl, jack, m, macdonald, q
            >>> h([2]).to("s")
            s[2]
            >>> s([1, 1]).to("p")
            1/2*p[1,1] - 1/2*p[2]
            >>> jack.P([2]).to("m")
            2/(alpha + 1)*m[1,1] + m[2]
            >>> jack.P([2]).to("s")
            (1 - alpha)/(alpha + 1)*s[1,1] + s[2]
            >>> hl.Qp([1, 1]).to("m")
            (1 + t)*m[1,1] + t*m[2]
            >>> (q * m([2])).to("s")
            -q*s[1,1] + q*s[2]
            >>> s([2]).to("HLP")
            t*HLP[1,1] + HLP[2]
            >>> s([2]).to("HLP") == hl.to_P(s([2]))
            True
            >>> macdonald.P([2]).to("m").to("McdP")
            McdP[2]
            >>> jack.P([2]).to("p")
            1/(alpha + 1)*p[1,1] + alpha/(alpha + 1)*p[2]
            >>> hl.Qp([1, 1]).to("p")
            (1/2 + 1/2*t)*p[1,1] + (-1/2 + 1/2*t)*p[2]

        The `p` values pin the convention: `s_11 = (p_1² − p_2)/2`, the second
        elementary symmetric function, which distinguishes it from `s_2`, where
        the sign is `+`. The power-sum target is the one conversion that
        divides — by z_μ — so it is the one whose coefficients leave ℤ (and,
        with parameters, leave the polynomial ring); the Jack value above is
        `J_2 = p_11 + α·p_2` scaled by `1/H_2`, and its coefficients still
        live over ℚ(α).

        A parametric basis expands into one classical basis — monomial for the
        Macdonald and Jack normalizations, Schur for Hall-Littlewood and `H̃` —
        because that is the basis its family's forward direction is defined in.
        Any other classical basis is that expansion followed by an ordinary
        change of basis, which is a ℤ-linear map on the partitions and so
        carries the coefficients through untouched. The Jack value vanishes at
        α = 1, where `P_λ` is the Schur function, which distinguishes this
        convention from its `α → 1/α` mirror.

        A parametric *target* runs the family's inverse expansion from the
        classical basis it reads, as the last two examples show. Each
        family's `to_P`, `to_Q`, `to_J`, `to_Qp` and `to_Htilde` method is
        the home of that expansion's convention; `to` delegates to them and
        adds nothing. Converting between two parametric tags goes through
        the classical basis their expansions share, and a pair over
        different coefficient rings refuses the same way the family method
        would.

        The coefficients are unchanged by the move: an element over ℚ(α) is
        still over ℚ(α) in whichever basis it is written.

        # Raises

        Raises `ValueError` unless `basis` names a basis.
        """
        basis = check_param_basis(basis)
        if basis == self._basis:
            return Sym(self._basis, self._terms, self._params)
        if basis not in BASES:
            from ._families import _INVERSE, _back_to, _expand

            if self._params:
                f = self if self._basis in BASES else _expand(self)
                return _back_to(f, basis, "to")
            src, inverse = _INVERSE[basis]
            return inverse(self.to(src))
        if self._params:
            from ._families import _convert, _expand

            if basis == "p":
                from ._families import _to_power

                return _to_power(self)
            if self._basis in BASES:
                return _convert(self, basis)
            expanded = _expand(self)
            return (
                expanded
                if expanded.basis == basis
                else _convert(expanded, basis)
            )
        if not self._terms:
            return Sym(basis, self._terms)
        pairs, scale = clear_denominators(self._numbers())
        src = self._basis
        if basis == "p":
            rows = [(la, Fraction(n, d)) for la, (n, d) in _c.to_power(pairs, src)]
            return Sym("p", restore(rows, scale))
        return Sym(basis, restore(_c.convert_terms(pairs, src, basis), scale))

    def omega(self) -> Sym:
        """The ω involution, returned in this element's basis.

        ω is its own inverse and exchanges `e` with `h`. It is defined on the
        Schur basis, where it conjugates the shape, and reaches any other by a
        change of basis on each side — so an element in a family's own basis
        comes back in it.

            >>> from symfn import hl, jack, m, q, s
            >>> s([2, 1, 1]).omega()
            s[3,1]
            >>> hl.Qp([2]).omega()
            HLQp[1,1] - t*HLQp[2]
            >>> (q * m([2, 1])).omega()
            -q*m[2,1] - 2*q*m[3]
            >>> jack.P([2]).omega()
            4*alpha/(alpha + 1)^2*JackP[1,1] + (1 - alpha)/(alpha + 1)*JackP[2]

        The first value pins the convention: ω conjugates the shape, so it is
        not the identity and not the antipode, which carries a sign. The third
        is `m([2, 1]).omega()` scaled by `q`, which is the check that a
        parameter is carried and not acted on.

        The Jack value pins which ω this is. It is the plain involution,
        `p_r ↦ (−1)^{r−1} p_r`, so α is carried and not touched; the
        α-deformed one sends `P_λ^{(α)}` to `Q_{λ'}^{(1/α)}`, which inverts
        the parameter. Sage's `.omega()` gives this same value.

        # Raises

        Raises `ValueError` if the element carries a coefficient class this
        layer does not know, and for a parametric basis whose inverse expansion
        runs over one it does not.
        """
        if self._params:
            from ._families import _hopf

            return _hopf(self, "omega")
        return self._through_schur(lambda a: _c.omega(a))

    def antipode(self) -> Sym:
        """The antipode S of the Hopf algebra, in this element's basis.

        `S(s_λ) = (−1)^|λ| s_{λ'}`, which is ω up to that sign — the values
        below are what separate them.

            >>> from symfn import macdonald, m, q, s
            >>> s([2, 1]).antipode()
            -s[2,1]
            >>> (q * m([2, 1])).antipode()
            q*m[2,1] + 2*q*m[3]
            >>> macdonald.P([1]).antipode()
            -McdP[1]

        Each shape here has odd degree, so each value is the negative of
        `omega`'s and the two are told apart; at even degree every convention
        agrees.

        # Raises

        Raises `ValueError` for the same reasons `omega` does.
        """
        if self._params:
            from ._families import _hopf

            return _hopf(self, "antipode")
        return self._through_schur(lambda a: _c.antipode(a))

    def plethysm(self, g: Sym) -> Sym:
        """The plethysm `f[g]`, with `f` this element, in the Schur basis.

            >>> from symfn import hl, jack, s, Poly
            >>> s([2]).plethysm(s([2]))
            s[2,2] + s[4]
            >>> hl.P([2]).plethysm(Poly("t", {1: 1}) * hl.P([1]))
            t^2*HLP[2]
            >>> hl.P([2]).plethysm(hl.P([1, 1]))
            (1 - t^3)*HLP[1,1,1,1] + HLP[2,2]
            >>> jack.P([2]).plethysm(jack.P([1, 1]))
            6*alpha/((alpha + 1)*(alpha + 2))*JackP[1,1,1,1] + JackP[2,2]

        The first distinguishes plethysm from the ordinary product, where
        `s_2 · s_2` also carries `s[3,1]`. The rest are Sage's values.

        **The parameters are part of the alphabet, so `p_n` raises them**: the
        second is `t²`, not `t`, and that is the value separating this
        convention from the one that holds `t` fixed — Sage's `exclude=`, which
        has no counterpart here. Over ℚ(α) the raising is α ↦ α^n, so a
        denominator `α + 1` can become `α² + 1` — irreducible, and carried in
        the coefficient's `tail` rather than among its atoms.

        # Raises

        Raises `ValueError` if `g` has a rational coefficient. Plethysm is not
        linear in `g`, so the denominator-clearing every other method here uses
        does not apply, and the contract layer takes integers. Raises
        `BaseRingError` unless both are over the same base ring.
        """
        if self._needs_ring(g):
            from ._families import _plethysm

            return _plethysm(self, g)
        g = self._same(g, "compose")
        gs = g.to("s")
        if any(isinstance(c, Fraction) for c in gs._numbers().values()):
            raise ValueError(
                "plethysm needs integer coefficients in its inner argument"
            )
        pairs, scale = clear_denominators(self.to("s")._numbers())
        inner = [(la, int(c)) for la, c in gs._numbers().items()]
        return Sym("s", restore(_c.plethysm(pairs, inner), scale)).to(self._basis)

    def internal_product(self, other: Sym) -> Sym:
        """The internal (Kronecker) product, in this element's basis.

            >>> from symfn import hl, jack, macdonald, s
            >>> s([2, 1]).internal_product(s([2, 1]))
            s[1,1,1] + s[2,1] + s[3]
            >>> jack.P([1]).internal_product(jack.P([1]))
            JackP[1]
            >>> macdonald.P([2]).internal_product(macdonald.P([1, 1]))
            (1 - t^2 - q^2 + q^2*t^2)/(1 - q*t)^2*McdP[1,1] + (-t + q)/(1 - q*t)*McdP[2]
            >>> hl.P([2, 1]).internal_product(hl.P([2, 1])).coefficient([3])
            1 + t^2 + 2*t^3 + t^4

        The first distinguishes this from the outer product, which is
        homogeneous of degree 6 rather than 3. The other three are Sage's. The
        structure constants are Kronecker coefficients, which are integers
        carrying no parameter, so the parameters are multiplied through rather
        than acted on. At `t = 0` the last is 1, one term of the first, since
        `P_λ(x; 0) = s_λ`.

        Two elements combine here, unlike `scalar` and `skew_by`, so they must
        be in the same basis, on the same grounds `*` refuses. An `other`
        without parameters is lifted into this element's base ring rather than
        refused.

        # Raises

        Raises `BasisError` unless both are in the same basis, `BaseRingError`
        unless both are over the same base ring, and `ValueError` if the
        power-sum route produces a non-integral coefficient.
        """
        if self._needs_ring(other):
            from ._families import _internal

            return _internal(self, other)
        return self._through_schur_pair(other, _c.internal_product)

    def scalar(self, other: Sym) -> AnyCoefficient:
        """The classical Hall inner product `⟨self, other⟩`, as one
        coefficient.

        This is the pairing with `⟨p_λ, p_μ⟩ = z_λ δ_{λμ}`, under which the
        Schur basis is orthonormal — which the values pin:

            >>> from symfn import s, p, h, m, macdonald, jack
            >>> s([2, 1]).scalar(s([2, 1])), s([2, 1]).scalar(s([3]))
            (1, 0)
            >>> p([2]).scalar(p([2]))
            2
            >>> h([2]).scalar(m([2])), h([2]).scalar(m([1, 1]))
            (1, 0)
            >>> macdonald.P([2]).scalar(macdonald.P([1, 1]))
            (-t + q)/(1 - q*t)
            >>> jack.P([2]).scalar(jack.P([1, 1]))
            (1 - alpha)/(alpha + 1)
            >>> jack.P([2, 1]).scalar(jack.P([2, 1]))
            (8 - 4*alpha + 5*alpha^2)/(alpha + 2)^2

        `⟨p_λ, p_λ⟩ = z_λ`, so the power-sum basis is orthogonal but not
        orthonormal — the second value is what says which. The parametric
        families are not orthogonal under this pairing: Macdonald and Jack
        `P` are orthogonal under their own deformed pairings — `scalar_qt`
        and `scalar_jack`, with `scalar_t` the Hall-Littlewood one — where
        both off-diagonal values above are 0. Under the Hall pairing they
        are the nonzero values shown — correct, not defects. The Macdonald
        value and the Jack norm are Sage's. The pairing is bilinear over
        whatever ring the coefficients live in, so the parameters are
        carried and never acted on.

        The `h`/`m` pair is the duality of those two bases, and it is why
        `other` may be in any basis: the pairing is defined on the ring, so two
        spellings of one argument give one number and there is nothing to
        refuse. That is unlike `+` and `*`, which combine two elements and do
        refuse. An `other` without parameters is lifted into this element's
        base ring rather than refused.

        # Raises

        Raises `ValueError` if this element carries a coefficient class this
        layer does not know, and `BaseRingError` if `other` is over a different
        base ring.
        """
        if self._needs_ring(other):
            from ._families import _scalar

            return _scalar(self, other)
        other = other if isinstance(other, Sym) else self._same(other, "pair")
        a, sa = clear_denominators(self.to("s")._numbers())
        b, sb = clear_denominators(other.to("s")._numbers())
        return exact(Fraction(_c.hall_inner_product(a, b), sa * sb))

    def scalar_t(self, other: Sym) -> AnyCoefficient:
        """The t-deformed Hall pairing `⟨self, other⟩_t` — Sage's `scalar_t`.

        This is the pairing with `⟨p_λ, p_μ⟩_t = δ_{λμ} z_λ Π 1/(1 − t^{λ_i})`,
        under which the Hall-Littlewood bases are orthogonal — the property
        `scalar` lacks for them — with `⟨P_λ, P_λ⟩_t = 1/b_λ(t)`. At t = 0 it
        degenerates to `scalar`.

            >>> from symfn import s, hl
            >>> s([1]).scalar_t(s([1]))
            1/(1 - t)
            >>> hl.P([2]).scalar_t(hl.P([1, 1]))
            0
            >>> hl.P([2]).scalar_t(hl.P([2]))
            (1 + t)/(1 - t^2)

        The first value pins the convention: `scalar` gives 1 there, and the
        reciprocal convention `z_λ Π (1 − t^{λ_i})` gives `1 − t`. The last is
        `1/(1 − t)` — the norm `1/b_(2)` — in the representation the sum
        arrives in; `==` compares values, not spellings.

        The value is a `QtFrac` whatever went in, because the pairing itself
        introduces the binomial denominators — so even a pair of parameter-free
        elements lands in the fraction field, and an element in `t` alone (or
        `q` alone, for LLT) is read into it. Both sides convert to the Schur
        basis first, so `other` may be in any basis, and a side without
        parameters is lifted rather than refused.

        # Raises

        Raises `ValueError` for coefficients in α — that ring pairs under
        `scalar_jack` — and for the `H̃` ring, whose denominators only
        `scalar_qt` reaches; `BaseRingError` if `other` is over a different
        base ring.
        """
        from ._families import _scalar_t

        other = other if isinstance(other, Sym) else self._same(other, "pair")
        return _scalar_t(self, other)

    def scalar_qt(self, other: Sym) -> AnyCoefficient:
        """The `(q,t)`-deformed Hall pairing `⟨self, other⟩_{q,t}` — Sage's
        `scalar_qt`.

        This is the pairing with
        `⟨p_λ, p_μ⟩_{q,t} = δ_{λμ} z_λ Π (1 − q^{λ_i})/(1 − t^{λ_i})`, under
        which Macdonald's `P` and `Q` are dual bases — the property `scalar`
        lacks for them. At q = t it degenerates to `scalar`, and at q = 0 to
        `scalar_t`.

            >>> from symfn import s, macdonald
            >>> s([1]).scalar_qt(s([1]))
            (1 - q)/(1 - t)
            >>> macdonald.P([2]).scalar_qt(macdonald.Q([1, 1]))
            0
            >>> macdonald.P([2]).scalar_qt(macdonald.Q([2]))
            1

        The first value pins the convention against its `q ↔ t` twist, which
        gives the reciprocal, and against `scalar_t`, which gives
        `1/(1 − t)`. ⚠️ This is not the pairing `H̃` is orthogonal under —
        that is the star product, whose weight carries an extra sign and
        `Π (1 − q^{λ_i})(1 − t^{λ_i})` — so `⟨H̃_μ, H̃_ν⟩_{q,t} ≠ 0` for
        `μ ≠ ν` is a value, not a defect.

        The value is a `QtFrac`, as in `scalar_t` — except over the `H̃`
        ring, whose `q^a − t^b` denominators come back in that ring's own
        `QtRatio`. Both sides convert to the Schur basis first, so `other`
        may be in any basis, and a side without parameters is lifted rather
        than refused.

        # Raises

        Raises `ValueError` for coefficients in α, which pair under
        `scalar_jack`, and `BaseRingError` if `other` is over a different
        base ring.
        """
        from ._families import _scalar_qt

        other = other if isinstance(other, Sym) else self._same(other, "pair")
        return _scalar_qt(self, other)

    def scalar_jack(self, other: Sym) -> AnyCoefficient:
        """The α-deformed Hall pairing `⟨self, other⟩_α` — Sage's
        `scalar_jack`.

        This is the pairing with `⟨p_λ, p_μ⟩_α = δ_{λμ} z_λ α^{ℓ(λ)}`, under
        which Jack's `P` and `Q` are dual bases — the property `scalar` lacks
        for them, and the pairing whose absence made the `scalar` values of
        Jack elements read as defects. At α = 1 it degenerates to `scalar`.

            >>> from symfn import s, jack
            >>> s([1]).scalar_jack(s([1]))
            alpha
            >>> jack.P([2]).scalar_jack(jack.P([1, 1]))
            0
            >>> jack.P([2]).scalar_jack(jack.Q([2]))
            1
            >>> jack.P([2]).scalar_jack(jack.P([2]))
            2*alpha^2/(alpha + 1)

        The first value pins the convention — `scalar` gives 1 — and the
        last is the norm `⟨P_λ, P_λ⟩_α = H'_λ/H_λ`, Sage's value. The value
        is an `AlphaFrac` whatever went in. Both sides convert to the
        monomial basis — the one Jack's expansions read — so `other` may be
        in any basis, and a side without parameters is lifted rather than
        refused.

        # Raises

        Raises `ValueError` for coefficients in `q` and `t`, which pair
        under `scalar_t` and `scalar_qt`, and `BaseRingError` if `other` is
        over a different base ring.
        """
        from ._families import _scalar_jack

        other = other if isinstance(other, Sym) else self._same(other, "pair")
        return _scalar_jack(self, other)

    def skew_by(self, g: Sym) -> Sym:
        """The element skewed by `g`, the adjoint of multiplication by `g`
        under the Hall inner product, in this element's basis.

            >>> from symfn import s, h, hl, macdonald
            >>> s([2, 1]).skew_by(s([1]))
            s[1,1] + s[2]
            >>> s([2, 1]).skew_by(h([1]))
            s[1,1] + s[2]
            >>> macdonald.P([2, 1]).skew_by(s([1]))
            (1 - t^2 - q^2*t + q^2*t^3)/((1 - q*t)*(1 - q*t^2))*McdP[1,1] + McdP[2]
            >>> hl.P([2, 1]).skew_by(h([1]))
            (1 - t^2)*HLP[1,1] + HLP[2]

        `s_{λ/μ}` for a partition `μ` is `skew(la, mu)`; this is the general
        form, linear in `g`. The last two are Sage's values.

        `g` keeps the basis it is written in, and that basis selects which rule
        runs rather than merely how `g` is read: `h`, `e` and `p` take the
        Pieri, dual-Pieri and Murnaghan-Nakayama paths, and `s`, `m` and `f`
        go through Littlewood-Richardson. So the first two values are the same
        answer by two different algorithms. All six have integer structure
        constants, so the parameters are carried and never acted on.

        Unlike `+` and `*`, a `g` in another basis is not a mismatch — this is
        not a combination of two elements of one ring but an operator built
        from `g` and applied here. A `g` without parameters is lifted into this
        element's base ring rather than refused.

        # Raises

        Raises `ValueError` if this element carries a coefficient class this
        layer does not know, and `BaseRingError` if `g` is over a different
        base ring.
        """
        if self._needs_ring(g):
            from ._families import _skew

            return _skew(self, g)
        g = g if isinstance(g, Sym) else self._same(g, "skew")
        a, sa = clear_denominators(self.to("s")._numbers())
        b, sb = clear_denominators(g._numbers())
        out = _c.skew_by(a, b, g._basis)
        return Sym("s", restore(out, sa * sb)).to(self._basis)

    def coproduct(self) -> dict[tuple[Partition, Partition], Coefficient]:
        """The coproduct Δ, as a `{(mu, nu): coefficient}` mapping over the
        Schur basis of each factor.

            >>> from symfn import s, hl, jack
            >>> s([2]).coproduct()
            {((), (2,)): 1, ((1,), (1,)): 1, ((2,), ()): 1}
            >>> hl.P([1]).coproduct()
            {((), (1,)): 1, ((1,), ()): 1}
            >>> jack.P([2]).coproduct()[(1,), (1,)]
            2/(alpha + 1)
            >>> jack.P([2]).coproduct()[(), (1, 1)]
            (1 - alpha)/(alpha + 1)

        The last two are Sage's, which writes them `2/(a+1)` and `(-a+1)/(a+1)`
        after expanding `s(JackP[2])`. The second pins the α convention: the
        `α → 1/α` mirror gives `(α − 1)/(α + 1)`, the negative of it. At α = 1
        both agree with `s([2]).coproduct()`, where they are 1 and 0.

        The keys are partition pairs rather than elements, because the result
        lives in a tensor square this type does not model — which is also why
        both factors come back in the Schur basis rather than in the basis this
        element is written in.

        `Δ(s_λ) = Σ c^λ_{μν} s_μ ⊗ s_ν` with Littlewood-Richardson
        coefficients, so the parameters are carried and never acted on.

        # Raises

        Raises `ValueError` if this element carries a coefficient class this
        layer does not know.
        """
        if self._params:
            from ._families import _coproduct

            return _coproduct(self)
        pairs, scale = clear_denominators(self.to("s")._numbers())
        rows = _c.coproduct(pairs)
        return {k: exact(Fraction(c, scale)) for k, c in rows}

    # --- leaving the ring ---------------------------------------------------

    def expand(self, n: int) -> dict[tuple[int, ...], Coefficient]:
        """The expansion in `n` variables, as a `{exponent vector: coefficient}`
        mapping with every vector of length `n`.

            >>> from symfn import macdonald, s
            >>> s([2]).expand(2)
            {(1, 1): 1, (2, 0): 1, (0, 2): 1}
            >>> macdonald.P([2]).expand(2)[1, 1]
            (1 - t + q - q*t)/(1 - q*t)

        Sage writes that coefficient of `x0*x1` as `(q*t-q+t-1)/(q*t-1)`, the
        same after clearing signs. The `q ↔ t` swap gives
        `(1 - q + t - q*t)/(1 - q*t)`, which is the twist to check.

        The order is the contract layer's: grouped by monomial term, and
        unspecified within a group.

        # Raises

        Raises `ValueError` if this element carries a coefficient class this
        layer does not know.
        """
        if self._params:
            from ._families import _expand_alphabet

            return _expand_alphabet(self, n)
        pairs, scale = clear_denominators(self._numbers())
        rows = _c.expand_alphabet(pairs, self._basis, n)
        return {v: exact(Fraction(c, scale)) for v, c in rows}

    def evaluate(self, xs: Sequence[int]) -> AnyCoefficient:
        """The value at the alphabet `xs`, a sequence of integers.

            >>> from symfn import hl, jack, s
            >>> s([2, 1]).evaluate([1, 1, 1])
            8
            >>> hl.P([2, 1]).evaluate([1, 1, 1])
            8 - t - t^2
            >>> jack.P([2, 1]).evaluate([1, 1, 1])
            (18 + 6*alpha)/(alpha + 2)

        `s_21(1,1,1) = 8` is the dimension of the `GL_3` irreducible, which
        `principal_specialization(3)` gives without an alphabet. The other two
        are Sage's values: the alphabet injects into the base ring, so it meets
        the shape and not the parameters, which ride through. Both degenerate
        to 8 — the first at `t = 0` and the second at α = 1, which are where
        each family is the Schur basis.

        # Raises

        Raises `TypeError` unless every alphabet entry is an `int` — the
        boundary crosses integer alphabets only — and `ValueError` if this
        element carries a coefficient class this layer does not know.
        """
        xs = list(xs)
        for x in xs:
            if not isinstance(x, int):
                raise TypeError(
                    "evaluate takes an alphabet of integers, not "
                    f"{type(x).__name__}"
                )
        if self._params:
            from ._families import _evaluate

            return _evaluate(self, xs)
        pairs, scale = clear_denominators(self.to("s")._numbers())
        return exact(Fraction(_c.evaluate_schur(pairs, xs), scale))


    def principal_specialization(
        self, n: int, q: AnyCoefficient | None = None
    ) -> AnyCoefficient:
        """The value at `1, q, …, q^{n−1}`, summed over the Schur expansion,
        with `q = 1` by default.

            >>> from symfn import hl, jack, s, Poly, t_hl
            >>> (s([2, 1]) + s([3])).principal_specialization(3)
            18
            >>> s([2, 1]).principal_specialization(3, q=2)
            90
            >>> hl.P([2, 1]).principal_specialization(3, q=Poly("t", {1: 1}))
            t + 2*t^2 + 2*t^3 + t^4
            >>> jack.P([2]).principal_specialization(3, q=2)
            (49 + 21*alpha)/(alpha + 1)
            >>> s([2]).principal_specialization(2, q=t_hl)
            1 + t + t^2

        Without `q` this is the value at `1^n`, the same `evaluate([1] * n)`
        gives, by a different route: this one weighs each shape by `s_λ(1^n)`
        and never lays out an alphabet.

        **`q` is an element of this element's own base ring** — ℚ for an
        element without parameters, so any integer or `Fraction` is one. That
        is how Sage spells it, and every ring here has it, including ℚ(α),
        where `principal_specialization_q` has no variable to introduce. The
        third value substitutes the `t` the element already carries, so it is
        `s_21(1,t,t²)` weighted by `P_21`'s own coefficients rather than a
        two-variable answer; the second is `s_21(1,2,4)`, which
        `principal_specialization_q(3).at(2)` also gives. A classical element
        also takes a `q` from a wider ring, and the sum lands there — the
        last value is `s_2(1,t)` in `ℤ[t]`.

        # Raises

        Raises `ValueError` if this element carries a coefficient class this
        layer does not know, or if `q` is not integral in that class's
        encoding, and `OverflowError` if a term's value exceeds what the
        contract layer's `principal_specialization` can represent, which
        reports the wall rather than wrapping.
        """
        if self._params:
            from ._families import _PRINCIPAL, _functional, _principal_at

            if q is None:
                return _functional(self, _PRINCIPAL, "the value at 1^n", n)
            return _principal_at(self, n, q)
        value: AnyCoefficient = q if q is not None else 1
        powers: list[AnyCoefficient] = [1]
        total: AnyCoefficient = 0
        for la, c in self.to("s")._numbers().items():
            if q is None:
                v = _c.principal_specialization(la, n)
                if v is None:
                    raise OverflowError(
                        f"s_{la} at 1^{n} exceeds the fixed-width specialization"
                    )
                total += c * v
            else:
                # The q-analogue's coefficients are integers, so substituting
                # is arithmetic in the base ring and needs no second route —
                # the powers are built by multiplication because the fraction
                # classes have `*` but no `**`.
                for k, w in enumerate(_c.principal_specialization_q(la, n)):
                    while len(powers) <= k:
                        powers.append(powers[-1] * value)
                    total += c * w * powers[k]
        return exact(total) if isinstance(total, (int, Fraction)) else total

    def principal_specialization_q(self, n: int) -> Poly | QtPoly:
        """The value at `1, q, …, q^{n−1}`, as a `Poly` in `q`.

            >>> from symfn import hl, s
            >>> s([2, 1]).principal_specialization_q(3)
            q + 2*q^2 + 2*q^3 + 2*q^4 + q^5
            >>> s([2, 1]).principal_specialization_q(3).at(1)
            8
            >>> hl.P([2, 1]).principal_specialization_q(3)
            q + 2*q^2 + 2*q^3 - q^3*t - q^3*t^2 + 2*q^4 + q^5
            >>> hl.P([2, 1]).principal_specialization_q(3).at(1, 0)
            8

        At `q = 1` this is `principal_specialization`, which is the check that
        the grading is the one that sums to it. The lowest power is `q^{n(λ)}`
        rather than `q^0`, which is what distinguishes this normalization from
        the one that divides the leading power out. Sage writes the third value
        as `q^5 + 2q^4 + (−t² − t + 2)q³ + 2q² + q`; setting `t = 0` there
        gives `s_21(1,q,q²)`, and the fourth does both at once, which is the
        check that the `q` this introduces and the `t` already present stayed
        apart.

        **Only where the base ring leaves room for `q`.** This layer's
        coefficient classes carry at most two variables, so an element already
        over ℚ(q,t) or ℚ(α) has nowhere to put the new one and is refused — the
        wall Sage reports as "the variable q is in the base ring, pass it
        explicitly". `principal_specialization(n, q=...)` is that explicit
        form: it substitutes an element of the ring rather than introducing a
        variable, and every ring here has it.

        # Raises

        Raises `ValueError` unless `n` is non-negative, if the base ring
        already carries `q`, and for the coefficient classes with no free
        variable at all.
        """
        if self._params:
            from ._families import _principal_q

            return _principal_q(self, n)
        from ._param import Poly

        total: dict[int, Coefficient] = {}
        for la, c in self.to("s")._numbers().items():
            for k, v in enumerate(_c.principal_specialization_q(la, n)):
                if v:
                    total[k] = total.get(k, 0) + c * v
        return Poly("q", total)

    def dimension(self) -> AnyCoefficient:
        """The dimension `Σ c_λ f^λ`, with `f^λ` the standard-tableaux count.

            >>> from symfn import hl, jack, s
            >>> (s([2, 1]) + s([1, 1, 1])).dimension()
            3
            >>> hl.P([2, 1]).dimension()
            2 - t - t^2
            >>> jack.P([2, 1]).dimension()
            6/(alpha + 2)

        `f^λ` depends on the shape alone, so the parameters are carried and
        never acted on. The last two degenerate to `f^{21} = 2`, the first at
        `t = 0` and the second at α = 1.

        # Raises

        Raises `ValueError` if this element carries a coefficient class this
        layer does not know, and `OverflowError` if a term's `f^λ` exceeds
        `u128`.
        """
        if self._params:
            from ._families import _DIMENSION, _functional

            return _functional(self, _DIMENSION, "the dimension")
        total: Coefficient = 0
        for la, c in self.to("s")._numbers().items():
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
        pairs, scale = clear_denominators(self.to("s")._numbers())
        return Sym("s", restore(call(pairs), scale)).to(self._basis)

    def _through_schur_pair(
        self,
        other: Sym,
        call: Callable[..., Iterable[tuple[Partition, int]]],
    ) -> Sym:
        other = self._same(other, "combine")
        a, sa = clear_denominators(self.to("s")._numbers())
        b, sb = clear_denominators(other.to("s")._numbers())
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
