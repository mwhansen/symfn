"""The coefficient types the parameter families return.

The contract layer hands a parameter family back as exponent-keyed rows: a
Hall-Littlewood coefficient is `[(t_exponent, coefficient), ...]`, a Macdonald
one is a numerator in `q` and `t` over a *factored* denominator, a Jack one is a
dense numerator in α over factored atoms. That encoding is deliberate — it
assumes no coefficient ring on the far side, which is what lets Sage, SymPy and
a bare interpreter each rebuild it in their own — and it is not what a person
wants to read.

The types here are that encoding with a `repr` and an evaluation map. They hold
the rows they were given, unchanged; `at` substitutes numbers into them. That
is arithmetic on plain data rather than symmetric-function mathematics:
nothing here can produce a coefficient the contract layer did not.

Evaluation is also what makes the families checkable from Python. Every one of
them degenerates to a classical basis at a particular parameter value —
`P_λ(x; q, q) = s_λ`, `P_λ(x; α = 1) = s_λ`, `Q'_λ(x; 0) = s_λ` — so `at` turns
a convention claim into a value a test can compare against a contract call.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from fractions import Fraction
from typing import (
    TYPE_CHECKING,
    Union,
    cast,
)

from ._bases import exact
from ._types import Coefficient

if TYPE_CHECKING:
    pass

#: A coefficient that carries parameters, as an element holds them.
ParamCoefficient = Union["Poly", "QtPoly", "QtFrac", "QtRatio", "AlphaFrac"]

__all__ = ["Poly", "QtPoly", "QtFrac", "QtRatio", "AlphaFrac"]


class Poly:
    """A polynomial in one variable with exact coefficients.

    This is what the `t`-indexed families return: Kostka-Foulkes polynomials,
    the coefficients of `Q'_λ`, and the `q`-analogue of a principal
    specialization.

        >>> from symfn import hl
        >>> c = hl.kostka_foulkes([2, 2], [1, 1, 1, 1])
        >>> c
        t^2 + t^4
        >>> c.at(1)
        2

    `K_{λμ}(1)` is the Kostka number, which is the check that pins this against
    `symfn.kostka_number`, and `K_{λμ}(t)` having no constant term is what
    distinguishes the charge convention from its transpose-statistic rival.
    """

    __slots__ = ("_var", "_terms")
    __module__ = "symfn"

    def __init__(
        self,
        var: str,
        terms: Mapping[int, Coefficient] | Iterable[tuple[int, Coefficient]],
    ) -> None:
        """Build from a variable name and `{exponent: coefficient}` or
        `(exponent, coefficient)` rows.
        """
        self._var = var
        items = terms.items() if hasattr(terms, "items") else terms
        self._terms = {int(k): exact(c) for k, c in items if c}

    @property
    def variable(self) -> str:
        """The variable's name, as it appears in the `repr`.

        >>> from symfn import hl
        >>> hl.kostka_foulkes([2], [1, 1]).variable
        't'
        """
        return self._var

    def coefficients(self) -> dict[int, Coefficient]:
        """A copy of the `{exponent: coefficient}` mapping, zero-free.

        >>> from symfn import hl
        >>> hl.kostka_foulkes([2], [1, 1]).coefficients()
        {1: 1}
        """
        return dict(self._terms)

    def degree(self) -> int | None:
        """The largest exponent present, or `None` for the zero polynomial.

        >>> from symfn import hl
        >>> hl.kostka_foulkes([2, 2], [1, 1, 1, 1]).degree()
        4
        """
        return max(self._terms) if self._terms else None

    def at(self, value: Coefficient) -> Coefficient:
        """The value at `value`, exactly.

        >>> from symfn import hl
        >>> hl.kostka_foulkes([2, 2], [1, 1, 1, 1]).at(0)
        0

        # Raises

        Raises `TypeError` unless `value` is an `int` or a `Fraction`.
        """
        value = exact(value)
        return exact(sum(c * value**k for k, c in self._terms.items()))

    __call__ = at

    def _same(self, other: Poly) -> None:
        """Refuse two polynomials in different variables.

        `t` and α are both `Poly`, and adding one to the other would build a
        value in neither variable rather than raising.
        """
        if self._var != other._var:
            raise ValueError(
                f"cannot combine a polynomial in {self._var} with one in "
                f"{other._var}"
            )

    def __neg__(self) -> Poly:
        return Poly(self._var, {k: -c for k, c in self._terms.items()})

    def __add__(self, other: object) -> Poly:
        if isinstance(other, (int, Fraction)):
            other = Poly(self._var, {0: other})
        if not isinstance(other, Poly):
            return NotImplemented
        self._same(other)
        out = dict(self._terms)
        for k, c in other._terms.items():
            out[k] = out.get(k, 0) + c
        return Poly(self._var, out)

    __radd__ = __add__

    def __sub__(self, other: object) -> Poly:
        neg = -other if isinstance(other, (int, Fraction, Poly)) else NotImplemented
        return NotImplemented if neg is NotImplemented else self.__add__(neg)

    def __rsub__(self, other: object) -> Poly:
        return (-self).__add__(other)

    def __mul__(self, other: object) -> Poly:
        if isinstance(other, (int, Fraction)):
            return Poly(self._var, {k: c * other for k, c in self._terms.items()})
        if not isinstance(other, Poly):
            return NotImplemented
        self._same(other)
        out: dict[int, Coefficient] = {}
        for j, a in self._terms.items():
            for k, b in other._terms.items():
                out[j + k] = out.get(j + k, 0) + a * b
        return Poly(self._var, out)

    __rmul__ = __mul__

    def __pow__(self, n: int) -> Poly:
        if not isinstance(n, int) or n < 0:
            raise ValueError(
                f"a polynomial power must be a non-negative int, got {n!r}"
            )
        out = Poly(self._var, {0: 1})
        for _ in range(n):
            out = out * self
        return out

    def __truediv__(self, other: object) -> Poly:
        if not isinstance(other, (int, Fraction)):
            return NotImplemented
        if not other:
            raise ZeroDivisionError("cannot divide a polynomial by zero")
        return self * (Fraction(1) / other)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Poly):
            return self._var == other._var and self._terms == other._terms
        if isinstance(other, (int, Fraction)):
            return self._terms == ({0: exact(other)} if other else {})
        return NotImplemented

    def __hash__(self) -> int:
        c = _constant_value(self)
        if c is not None:
            return hash(c)
        return hash((self._var, frozenset(self._terms.items())))

    def __bool__(self) -> bool:
        return bool(self._terms)

    def __repr__(self) -> str:
        return _sum(
            (c, _power(self._var, k)) for k, c in sorted(self._terms.items())
        )


class QtPoly:
    """A polynomial in `q` and `t` with exact coefficients.

    This is the numerator of a Macdonald coefficient, and the whole of a
    `(q,t)`-Kostka polynomial or a `∇` output.

        >>> from symfn import macdonald
        >>> macdonald.qt_kostka([2], [1, 1])
        t
        >>> macdonald.qt_kostka([2], [1, 1]).at(q=3, t=5)
        5

    `K̃_{(2),(11)}(q, t) = t` rather than `q` is the orientation of the
    Garsia-Haiman convention, and the value that separates it from its
    `q ↔ t` mirror.
    """

    __slots__ = ("_terms",)
    __module__ = "symfn"

    def __init__(
        self,
        terms: Mapping[tuple[int, int], Coefficient]
        | Iterable[tuple[int, int, Coefficient]],
    ) -> None:
        """Build from `{(q_exponent, t_exponent): coefficient}` or
        `(q_exponent, t_exponent, coefficient)` rows.
        """
        items = terms.items() if hasattr(terms, "items") else (
            ((a, b), c) for a, b, c in terms
        )
        self._terms = {(int(a), int(b)): exact(c) for (a, b), c in items if c}

    def coefficients(self) -> dict[tuple[int, int], Coefficient]:
        """A copy of the `{(q_exponent, t_exponent): coefficient}` mapping.

        >>> from symfn import macdonald
        >>> macdonald.qt_kostka([2], [1, 1]).coefficients()
        {(0, 1): 1}
        """
        return dict(self._terms)

    def at(self, q: Coefficient, t: Coefficient) -> Coefficient:
        """The value at `q` and `t`, exactly.

        >>> from symfn import macdonald
        >>> macdonald.qt_kostka([3, 1], [1, 1, 1, 1]).at(q=1, t=1)
        3

        # Raises

        Raises `TypeError` unless both are an `int` or a `Fraction`.
        """
        q, t = exact(q), exact(t)
        return exact(sum(c * q**a * t**b for (a, b), c in self._terms.items()))

    __call__ = at

    def __neg__(self) -> QtPoly:
        return QtPoly({k: -c for k, c in self._terms.items()})

    def __add__(self, other: object) -> QtPoly:
        if isinstance(other, (int, Fraction)):
            other = QtPoly({(0, 0): other})
        if not isinstance(other, QtPoly):
            return NotImplemented
        out = dict(self._terms)
        for k, c in other._terms.items():
            out[k] = out.get(k, 0) + c
        return QtPoly(out)

    __radd__ = __add__

    def __sub__(self, other: object) -> QtPoly:
        neg = -other if isinstance(other, (int, Fraction, QtPoly)) else NotImplemented
        return NotImplemented if neg is NotImplemented else self.__add__(neg)

    def __rsub__(self, other: object) -> QtPoly:
        return (-self).__add__(other)

    def __mul__(self, other: object) -> QtPoly:
        if isinstance(other, (int, Fraction)):
            return QtPoly({k: c * other for k, c in self._terms.items()})
        if not isinstance(other, QtPoly):
            return NotImplemented
        out: dict[tuple[int, int], Coefficient] = {}
        for (a, b), x in self._terms.items():
            for (c, d), y in other._terms.items():
                out[(a + c, b + d)] = out.get((a + c, b + d), 0) + x * y
        return QtPoly(out)

    __rmul__ = __mul__

    def __pow__(self, n: int) -> QtPoly:
        if not isinstance(n, int) or n < 0:
            raise ValueError(
                f"a polynomial power must be a non-negative int, got {n!r}"
            )
        out = QtPoly({(0, 0): 1})
        for _ in range(n):
            out = out * self
        return out

    def __truediv__(self, other: object) -> QtPoly:
        if not isinstance(other, (int, Fraction)):
            return NotImplemented
        if not other:
            raise ZeroDivisionError("cannot divide a polynomial by zero")
        return self * (Fraction(1) / other)

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QtPoly):
            return self._terms == other._terms
        if isinstance(other, (int, Fraction)):
            return self._terms == ({(0, 0): exact(other)} if other else {})
        return NotImplemented

    def __hash__(self) -> int:
        c = _constant_value(self)
        if c is not None:
            return hash(c)
        return hash(frozenset(self._terms.items()))

    def __bool__(self) -> bool:
        return bool(self._terms)

    def __repr__(self) -> str:
        return _sum(
            (c, _monomial(_power("q", a), _power("t", b)))
            for (a, b), c in sorted(self._terms.items())
        )


class QtFrac:
    """A Macdonald coefficient: a `QtPoly` over a factored denominator.

    The denominator is a list of `(a, b, multiplicity)` standing for
    `(1 − q^a t^b)^multiplicity`, and it is kept factored because that is the
    representation rather than an optimization — a caller rebuilding this in a
    fraction field wants the factors, and expanding here to re-factor there is
    work done twice.

        >>> from symfn import macdonald
        >>> macdonald.Q([1]).to("m").coefficient([1])
        (1 - t)/(1 - q)
        >>> macdonald.P([1]).to("m").coefficient([1])
        1

    `Q_(1) = (1 − t)/(1 − q) · m_1` where `P_(1) = m_1` is what separates the
    two normalizations at the smallest shape.
    """

    __slots__ = ("_num", "_den")
    __module__ = "symfn"

    def __init__(
        self,
        numerator: QtPoly | Iterable[tuple[int, int, Coefficient]],
        denominator: Iterable[tuple[int, int, int]] = (),
    ) -> None:
        """Build from `(q, t, coefficient)` numerator rows and
        `(q, t, multiplicity)` denominator factors.
        """
        self._num = numerator if isinstance(numerator, QtPoly) else QtPoly(numerator)
        self._den = tuple(sorted((int(a), int(b), int(k)) for a, b, k in denominator))

    @property
    def numerator(self) -> QtPoly:
        """The numerator, as a `QtPoly`.

        >>> from symfn import macdonald
        >>> macdonald.Q([1]).to("m").coefficient([1]).numerator
        1 - t
        """
        return self._num

    @property
    def denominator(self) -> tuple[tuple[int, int, int], ...]:
        """The denominator's factors, as `(q, t, multiplicity)` triples for
        `(1 − q^a t^b)^multiplicity`. Empty when the coefficient is a
        polynomial, which is the integral form's signature.

            >>> from symfn import macdonald
            >>> macdonald.Q([1]).to("m").coefficient([1]).denominator
            ((1, 0, 1),)
            >>> macdonald.J([1, 1]).to("m").coefficient([1, 1]).denominator
            ()
        """
        return self._den

    def at(self, q: Coefficient, t: Coefficient) -> Coefficient:
        """The value at `q` and `t`, exactly.

        >>> from fractions import Fraction
        >>> from symfn import macdonald
        >>> macdonald.Q([1]).to("m").coefficient([1]).at(q=0, t=Fraction(1, 2))
        Fraction(1, 2)

        # Raises

        Raises `ZeroDivisionError`, naming the factor, when a denominator
        factor vanishes at the given values.
        """
        q, t = exact(q), exact(t)
        value = Fraction(self._num.at(q, t))
        for a, b, k in self._den:
            factor = 1 - q**a * t**b
            if factor == 0:
                raise ZeroDivisionError(
                    f"the denominator factor (1 - q^{a}*t^{b}) vanishes here"
                )
            value /= Fraction(factor) ** k
        return exact(value)

    __call__ = at

    def __neg__(self) -> QtFrac:
        from ._families import _coeff_scale

        return cast("QtFrac", _coeff_scale(self, -1))

    def __add__(self, other: object) -> QtFrac:
        if not isinstance(other, (int, Fraction, QtPoly, QtFrac)):
            return NotImplemented
        from ._families import _coeff_add

        return cast("QtFrac", _coeff_add(self, other))

    __radd__ = __add__

    def __sub__(self, other: object) -> QtFrac:
        if not isinstance(other, (int, Fraction, QtPoly, QtFrac)):
            return NotImplemented
        return self.__add__(-other)

    def __rsub__(self, other: object) -> QtFrac:
        return (-self).__add__(other)

    def __mul__(self, other: object) -> QtFrac:
        if not isinstance(other, (int, Fraction, QtPoly, QtFrac)):
            return NotImplemented
        from ._families import _coeff_scale

        return cast("QtFrac", _coeff_scale(self, other))

    __rmul__ = __mul__

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QtFrac):
            return self._num == other._num and self._den == other._den
        if isinstance(other, (int, Fraction)):
            return not self._den and self._num == other
        return NotImplemented

    def __hash__(self) -> int:
        c = _constant_value(self)
        if c is not None:
            return hash(c)
        return hash((self._num, self._den))

    def __bool__(self) -> bool:
        return bool(self._num)

    def __repr__(self) -> str:
        if not self._den:
            return repr(self._num)
        factors = [
            f"(1 - {_monomial(_power('q', a), _power('t', b)) or '1'})"
            + (f"^{k}" if k > 1 else "")
            for a, b, k in self._den
        ]
        above = repr(self._num)
        if _has_top_level_sum(above):
            above = f"({above})"
        # A denominator of several factors is parenthesized as a whole:
        # `n/(1 - q)*(1 - t)` would read as a product with the quotient.
        below = factors[0] if len(factors) == 1 else "(" + "*".join(factors) + ")"
        return f"{above}/{below}"


class QtRatio:
    """A modified-Macdonald coefficient: a `QtPoly` over factored atoms.

    `QtFrac` holds the denominators the `P`, `Q` and `J` expansions produce,
    every factor a `1 − q^a t^b`. Expanding *into* `H̃` divides by `w_μ`
    instead, whose factors are `q^a − t^b`, and a product of those is not a
    product of the first kind — so this type carries **two** families, tagged
    `0` for `1 − q^a t^b` and `1` for `q^a − t^b`.

        >>> from symfn import macdonald, s
        >>> macdonald.to_Htilde(s([2])).coefficient([1, 1])
        q/(q - t)
        >>> macdonald.to_Htilde(s([2])).coefficient([1, 1]).denominator
        ((1, 1, 1, 1),)

    `s_2 = q/(q−t)·H̃_11 − t/(q−t)·H̃_2`. `H̃` is not symmetric in `q` and `t`,
    so the swap gives a different answer rather than an error; which variable
    sits upstairs on the column shape is the convention.

    ⚠️ The two families are why a kind tag is needed at all, and why `q^0 − t^b`
    is refused as kind `1`: it *is* `1 − t^b`, and letting it in under the
    other name would leave two spellings of one polynomial that never cancel.
    """

    __slots__ = ("_num", "_den")
    __module__ = "symfn"

    def __init__(
        self,
        numerator: QtPoly | Iterable[tuple[int, int, Coefficient]],
        denominator: Iterable[tuple[int, int, int, int]] = (),
    ) -> None:
        """Build from `(q, t, coefficient)` rows and `(kind, a, b,
        multiplicity)` atoms.
        """
        self._num = numerator if isinstance(numerator, QtPoly) else QtPoly(numerator)
        atoms: dict[tuple[int, int, int], int] = {}
        for kind, a, b, k in denominator:
            if kind not in (0, 1):
                raise ValueError(
                    f"unknown atom kind {kind}; 0 is 1 - q^a t^b and 1 is "
                    "q^a - t^b"
                )
            if kind == 0 and a == 0 and b == 0:
                raise ValueError("not a denominator atom: 1 - q^0 t^0 is zero")
            if kind == 1 and (a == 0 or b == 0):
                raise ValueError(
                    f"q^{a} - t^{b} is not of kind 1: with a zero exponent it "
                    "is 1 - q^a t^b up to sign, which is kind 0"
                )
            if k:
                atoms[(kind, a, b)] = atoms.get((kind, a, b), 0) + k
        self._den = tuple(
            (kind, a, b, k) for (kind, a, b), k in sorted(atoms.items())
        )

    @property
    def numerator(self) -> QtPoly:
        """The numerator, as a `QtPoly`.

        >>> from symfn import macdonald, s
        >>> macdonald.to_Htilde(s([2])).coefficient([1, 1]).numerator
        q
        """
        return self._num

    @property
    def denominator(self) -> tuple[tuple[int, int, int, int], ...]:
        """The denominator's atoms with their multiplicities, as
        `(kind, a, b, multiplicity)` in ascending order — kind `0` standing for
        `1 − q^a t^b` and kind `1` for `q^a − t^b`.

            >>> from symfn import macdonald, s
            >>> macdonald.to_Htilde(s([2])).coefficient([1, 1]).denominator
            ((1, 1, 1, 1),)
            >>> macdonald.to_Htilde(s([1])).coefficient([1]).denominator
            ()

        The empty tuple is a denominator of 1, which is what a polynomial
        coefficient has.
        """
        return self._den

    def _atom_at(
        self, kind: int, a: int, b: int, q: Coefficient, t: Coefficient
    ) -> Coefficient:
        """One atom's value, by kind."""
        return 1 - q**a * t**b if kind == 0 else q**a - t**b

    def at(self, q: Coefficient, t: Coefficient) -> Coefficient:
        """The value at `q` and `t`, exactly.

        >>> from symfn import macdonald, s
        >>> macdonald.to_Htilde(s([2])).coefficient([1, 1]).at(q=1, t=0)
        1

        # Raises

        Raises `ZeroDivisionError`, naming the atom, when one vanishes at the
        given values — which `H̃` does on the diagonal `q = t`.
        """
        q, t = exact(q), exact(t)
        value = Fraction(self._num.at(q, t))
        for kind, a, b, k in self._den:
            factor = self._atom_at(kind, a, b, q, t)
            if factor == 0:
                raise ZeroDivisionError(
                    f"the denominator atom ({_atom_repr(kind, a, b)}) vanishes "
                    f"at q = {q}, t = {t}"
                )
            value /= Fraction(factor) ** k
        return exact(value)

    __call__ = at

    def __neg__(self) -> QtRatio:
        from ._families import _coeff_scale

        return cast("QtRatio", _coeff_scale(self, -1))

    def __add__(self, other: object) -> QtRatio:
        if not isinstance(other, (int, Fraction, QtPoly, QtRatio)):
            return NotImplemented
        from ._families import _coeff_add

        return cast("QtRatio", _coeff_add(self, other))

    __radd__ = __add__

    def __sub__(self, other: object) -> QtRatio:
        if not isinstance(other, (int, Fraction, QtPoly, QtRatio)):
            return NotImplemented
        return self.__add__(-other)

    def __rsub__(self, other: object) -> QtRatio:
        return (-self).__add__(other)

    def __mul__(self, other: object) -> QtRatio:
        if not isinstance(other, (int, Fraction, QtPoly, QtRatio)):
            return NotImplemented
        from ._families import _coeff_scale

        return cast("QtRatio", _coeff_scale(self, other))

    __rmul__ = __mul__

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QtRatio):
            return self._num == other._num and self._den == other._den
        if isinstance(other, (int, Fraction)):
            return not self._den and self._num == other
        return NotImplemented

    def __hash__(self) -> int:
        c = _constant_value(self)
        if c is not None:
            return hash(c)
        return hash((self._num, self._den))

    def __bool__(self) -> bool:
        return bool(self._num)

    def __repr__(self) -> str:
        if not self._den:
            return repr(self._num)
        factors = [
            f"({_atom_repr(kind, a, b)})" + (f"^{k}" if k > 1 else "")
            for kind, a, b, k in self._den
        ]
        above = repr(self._num)
        if _has_top_level_sum(above):
            above = f"({above})"
        below = factors[0] if len(factors) == 1 else "(" + "*".join(factors) + ")"
        return f"{above}/{below}"


def _atom_repr(kind: int, a: int, b: int) -> str:
    """One atom as it prints: `1 - q^a*t^b` or `q^a - t^b`."""
    if kind == 0:
        return f"1 - {_monomial(_power('q', a), _power('t', b)) or '1'}"
    return f"{_power('q', a)} - {_power('t', b)}"


class AlphaFrac:
    """A Jack coefficient: a polynomial in α over factored linear atoms.

    The value is
    `(Σ_k numerator[k]·α^k) / (scale · ∏ (u·α + v)^multiplicity · tail(α))`,
    with the numerator dense in the α-exponent and the atoms primitive
    (`gcd(u, v) = 1`), so the factorization is canonical.

        >>> from symfn import jack
        >>> jack.P([2]).to("m").coefficient([1, 1])
        2/(alpha + 1)
        >>> jack.P([2]).to("m").coefficient([1, 1]).at(1)
        1

    `P_(2) = 2/(α+1)·m_11 + m_2` is monic in `m_λ`, which is what separates `P`
    from `Q` and `J`; at α = 1 it becomes `s_(2) = m_11 + m_2`.

    **`tail` is empty except after a plethysm.** It is a further denominator
    factor, dense in α like the numerator, that is not a product of linear
    forms — `p_n` raises the variable, so an atom `α + 1` becomes `α² + 1`,
    which is irreducible over ℚ. Every other operation leaves it empty.
    """

    __slots__ = ("_num", "_atoms", "_scale", "_tail")
    __module__ = "symfn"

    def __init__(
        self,
        numerator: Sequence[Coefficient],
        atoms: Iterable[tuple[int, int, int]] = (),
        scale: int = 1,
        tail: Sequence[Coefficient] = (),
    ) -> None:
        """Build from a dense numerator, `(u, v, multiplicity)` atoms, an
        integer scale, and a dense general denominator factor.
        """
        self._num = tuple(exact(c) for c in numerator)
        self._atoms = tuple(sorted((int(u), int(v), int(k)) for u, v, k in atoms))
        self._scale = int(scale)
        self._tail = tuple(exact(c) for c in tail)

    @property
    def numerator(self) -> tuple[Coefficient, ...]:
        """The numerator, dense in the α-exponent: index `k` is the
        coefficient of `α^k`.

            >>> from symfn import jack
            >>> jack.P([2]).to("m").coefficient([1, 1]).numerator
            (2,)
        """
        return self._num

    @property
    def atoms(self) -> tuple[tuple[int, int, int], ...]:
        """The denominator's atoms, as `(u, v, multiplicity)` for
        `(u·α + v)^multiplicity`.

            >>> from symfn import jack
            >>> jack.P([2]).to("m").coefficient([1, 1]).atoms
            ((1, 1, 1),)
        """
        return self._atoms

    @property
    def scale(self) -> int:
        """The integer the denominator also carries.

        >>> from symfn import jack
        >>> jack.P([2]).to("m").coefficient([1, 1]).scale
        1
        """
        return self._scale

    @property
    def tail(self) -> tuple[Coefficient, ...]:
        """The denominator's general factor, dense in the α-exponent — empty
        unless a plethysm put something there.

            >>> from symfn import jack
            >>> jack.P([2]).to("m").coefficient([1, 1]).tail
            ()
            >>> jack.P([2]).plethysm(jack.P([2])).coefficient([2, 2]).tail
            (1, 0, 1)

        The second is `α² + 1`, which `p_2` produced from `α + 1` by raising
        the variable. It is irreducible over ℚ, which is why it cannot join
        the atoms.
        """
        return self._tail

    def at(self, alpha: Coefficient) -> Coefficient:
        """The value at `alpha`, exactly.

        >>> from symfn import jack
        >>> jack.P([3, 1]).to("m").coefficient([2, 1, 1]).at(2)
        Fraction(11, 9)

        # Raises

        Raises `ZeroDivisionError`, naming the atom, when one vanishes at
        `alpha`.
        """
        alpha = exact(alpha)
        value = Fraction(sum(c * alpha**k for k, c in enumerate(self._num)))
        value /= self._scale
        for u, v, k in self._atoms:
            factor = u * alpha + v
            if factor == 0:
                raise ZeroDivisionError(
                    f"the denominator atom ({u}*alpha + {v}) vanishes at "
                    f"alpha = {alpha}"
                )
            value /= Fraction(factor) ** k
        if self._tail:
            below = sum(c * alpha**k for k, c in enumerate(self._tail))
            if below == 0:
                raise ZeroDivisionError(
                    f"the denominator factor {self._tail} vanishes at "
                    f"alpha = {alpha}"
                )
            value /= Fraction(below)
        return exact(value)

    __call__ = at

    def __neg__(self) -> AlphaFrac:
        from ._families import _coeff_scale

        return cast("AlphaFrac", _coeff_scale(self, -1))

    def __add__(self, other: object) -> AlphaFrac:
        if not isinstance(other, (int, Fraction, Poly, AlphaFrac)):
            return NotImplemented
        from ._families import _coeff_add

        return cast("AlphaFrac", _coeff_add(self, other))

    __radd__ = __add__

    def __sub__(self, other: object) -> AlphaFrac:
        if not isinstance(other, (int, Fraction, Poly, AlphaFrac)):
            return NotImplemented
        return self.__add__(-other)

    def __rsub__(self, other: object) -> AlphaFrac:
        return (-self).__add__(other)

    def __mul__(self, other: object) -> AlphaFrac:
        if not isinstance(other, (int, Fraction, Poly, AlphaFrac)):
            return NotImplemented
        from ._families import _coeff_scale

        return cast("AlphaFrac", _coeff_scale(self, other))

    __rmul__ = __mul__

    def __eq__(self, other: object) -> bool:
        if isinstance(other, AlphaFrac):
            return (self._num, self._atoms, self._scale, self._tail) == (
                other._num,
                other._atoms,
                other._scale,
                other._tail,
            )
        if isinstance(other, (int, Fraction)):
            # `num[:1]` alone would call `3 + 5α` equal to 3: the higher
            # α-coefficients have to be zero before the constant term decides.
            return (
                not self._atoms
                and not self._tail
                and self._scale == 1
                and not any(self._num[1:])
                and (self._num[0] if self._num else 0) == exact(other)
            )
        return NotImplemented

    def __hash__(self) -> int:
        c = _constant_value(self)
        if c is not None:
            return hash(c)
        return hash((self._num, self._atoms, self._scale, self._tail))

    def __bool__(self) -> bool:
        return any(self._num)

    def __repr__(self) -> str:
        above = _sum(
            (c, _power("alpha", k)) for k, c in enumerate(self._num) if c
        )
        factors = []
        if self._scale != 1:
            factors.append(str(self._scale))
        for u, v, k in self._atoms:
            atom = _sum([(u, "alpha"), (v, "")])
            if _has_top_level_sum(atom):
                atom = f"({atom})"
            factors.append(atom + (f"^{k}" if k > 1 else ""))
        if self._tail:
            below = _sum(
                (c, _power("alpha", k)) for k, c in enumerate(self._tail) if c
            )
            factors.append(f"({below})" if _has_top_level_sum(below) else below)
        if not factors:
            return above
        if _has_top_level_sum(above):
            above = f"({above})"
        below = factors[0] if len(factors) == 1 else "(" + "*".join(factors) + ")"
        return f"{above}/{below}"


def _constant_value(c: object) -> Coefficient | None:
    """The number a coefficient equals, or `None` when it is not constant.

    Every class here answers `==` against an `int` or a `Fraction`, so hash
    can only follow equality if a constant hashes as the number it equals —
    each `__hash__` above routes through this, and `Sym` reads it to compare
    constant elements across coefficient classes.
    """
    if isinstance(c, (int, Fraction)):
        return c
    if isinstance(c, Poly):
        return c._terms.get(0, 0) if set(c._terms) <= {0} else None
    if isinstance(c, QtPoly):
        return c._terms.get((0, 0), 0) if set(c._terms) <= {(0, 0)} else None
    if isinstance(c, (QtFrac, QtRatio)):
        return None if c._den else _constant_value(c._num)
    if isinstance(c, AlphaFrac):
        if c._atoms or c._tail or c._scale != 1 or any(c._num[1:]):
            return None
        return c._num[0] if c._num else 0
    return None


#: The Macdonald and LLT parameter `q`, as a `QtPoly`, so a coefficient can be
#: written as `q` rather than as the rows that encode it.
q: QtPoly = QtPoly({(1, 0): 1})

#: The parameter `t`, as a `QtPoly` — the one Macdonald and `H̃` take. The
#: Hall-Littlewood families are in `t` alone and take `t_hl` instead, because
#: their coefficients are a `Poly` in one variable and the two do not mix.
t: QtPoly = QtPoly({(0, 1): 1})

#: The Hall-Littlewood parameter `t`, as a one-variable `Poly`.
t_hl: Poly = Poly("t", {1: 1})

#: The LLT parameter `q`, as a one-variable `Poly` — the LLT family is in `q`
#: alone, so like the Hall-Littlewood families it is over a one-variable ring
#: rather than the `q`, `t` one the `QtPoly` `q` names.
q_llt: Poly = Poly("q", {1: 1})

#: The Jack parameter α, as a one-variable `Poly`.
alpha: Poly = Poly("alpha", {1: 1})


def _as_one_variable(c: QtPoly, var: str) -> Poly | None:
    """`c` as a one-variable `Poly` in `var`, or `None` when `c` involves the
    other variable.

    The exported `q` and `t` are `QtPoly`, and the Hall-Littlewood and LLT
    elements are over one-variable rings — so a `QtPoly` supported on `var`
    alone is the same scalar in a wider encoding, and the operations demote
    it rather than refuse it for its class.
    """
    if var not in ("q", "t"):
        return None
    out: dict[int, Coefficient] = {}
    for (a, b), v in c._terms.items():
        if b if var == "q" else a:
            return None
        out[a if var == "q" else b] = v
    return Poly(var, out)


#: What may stand in for an element in the parametric arithmetic: a scalar of
#: the base ring, or an integer or rational that injects into it. The same set
#: for `+`, `-` and `*`, so a scalar that can multiply an element can add to
#: it.
_SCALARS = (int, Fraction, Poly, QtPoly, QtFrac, QtRatio, AlphaFrac)


def _power(var: str, k: int) -> str:
    """`""` for exponent 0, `t` for 1, `t^3` otherwise."""
    return "" if k == 0 else var if k == 1 else f"{var}^{k}"


def _monomial(*powers: str) -> str:
    """Join rendered powers with `*`, dropping the empty ones: `q^2*t`."""
    return "*".join(x for x in powers if x)


def _has_top_level_sum(text: str) -> bool:
    """Whether `text` adds at depth zero, so using it as a factor needs
    parentheses.

        >>> _has_top_level_sum("1 - t")
        True
        >>> _has_top_level_sum("(1 - t)/(1 - q)")
        False
        >>> _has_top_level_sum("-t")
        False
    """
    depth = 0
    for i, ch in enumerate(text):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        elif ch in "+-" and depth == 0 and i and text[i - 1] == " ":
            return True
    return False


def _factor(coefficient: object) -> tuple[str, str]:
    """A coefficient rendered for use in front of a basis element, as
    `(sign, body)` with the body ready to prefix `*basis[...]`.

    A coefficient that adds at the top level is parenthesized; a leading minus
    on an atomic one is lifted out so the sum reads `a - b` rather than
    `a + -b`.
    """
    text = repr(coefficient)
    if _has_top_level_sum(text):
        return "+", f"({text})"
    if text.startswith("-"):
        return "-", text[1:]
    return "+", text


def _sum(pieces: Iterable[tuple[Coefficient, str]]) -> str:
    """Render `(coefficient, monomial)` pairs as `a*x + b`, signs folded in."""
    out = ""
    for c, mono in pieces:
        if not c:
            continue
        sign = "-" if c < 0 else "+"
        mag = -c if c < 0 else c
        if not mono:
            body = _num(mag)
        elif mag == 1:
            body = mono
        else:
            body = f"{_num(mag)}*{mono}"
        out += body if not out and sign == "+" else (
            f"-{body}" if not out else f" {sign} {body}"
        )
    return out or "0"


def _num(c: Coefficient) -> str:
    return f"{c.numerator}/{c.denominator}" if isinstance(c, Fraction) else str(c)
