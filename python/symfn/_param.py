"""The coefficient types the parameter families return, and `Param` over them.

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

from collections.abc import Iterable, Iterator, Mapping, Sequence
from fractions import Fraction
from typing import (
    TYPE_CHECKING,
    Union,
)

from ._bases import BASES, check_basis, check_param_basis, exact
from ._sym import _partition
from ._types import Coefficient, ParamBasis, Partition, PartitionArg

if TYPE_CHECKING:
    from ._sym import Sym

#: A coefficient that carries parameters, as `Param` holds them.
ParamCoefficient = Union["Poly", "QtPoly", "QtFrac", "QtRatio", "AlphaFrac"]

__all__ = ["Poly", "QtPoly", "QtFrac", "QtRatio", "AlphaFrac", "Param"]


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

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Poly):
            return self._var == other._var and self._terms == other._terms
        if isinstance(other, (int, Fraction)):
            return self._terms == ({0: exact(other)} if other else {})
        return NotImplemented

    def __hash__(self) -> int:
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

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QtPoly):
            return self._terms == other._terms
        if isinstance(other, (int, Fraction)):
            return self._terms == ({(0, 0): exact(other)} if other else {})
        return NotImplemented

    def __hash__(self) -> int:
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

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QtFrac):
            return self._num == other._num and self._den == other._den
        if isinstance(other, (int, Fraction)):
            return not self._den and self._num == other
        return NotImplemented

    def __hash__(self) -> int:
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

    def __eq__(self, other: object) -> bool:
        if isinstance(other, QtRatio):
            return self._num == other._num and self._den == other._den
        if isinstance(other, (int, Fraction)):
            return not self._den and self._num == other
        return NotImplemented

    def __hash__(self) -> int:
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

    The value is `(Σ_k numerator[k]·α^k) / (scale · ∏ (u·α + v)^multiplicity)`,
    with the numerator dense in the α-exponent and the atoms primitive
    (`gcd(u, v) = 1`), so the factorization is canonical.

        >>> from symfn import jack
        >>> jack.P([2]).to("m").coefficient([1, 1])
        2/(alpha + 1)
        >>> jack.P([2]).to("m").coefficient([1, 1]).at(1)
        1

    `P_(2) = 2/(α+1)·m_11 + m_2` is monic in `m_λ`, which is what separates `P`
    from `Q` and `J`; at α = 1 it becomes `s_(2) = m_11 + m_2`.
    """

    __slots__ = ("_num", "_atoms", "_scale")
    __module__ = "symfn"

    def __init__(
        self,
        numerator: Sequence[Coefficient],
        atoms: Iterable[tuple[int, int, int]] = (),
        scale: int = 1,
    ) -> None:
        """Build from a dense numerator, `(u, v, multiplicity)` atoms, and an
        integer scale.
        """
        self._num = tuple(exact(c) for c in numerator)
        self._atoms = tuple(sorted((int(u), int(v), int(k)) for u, v, k in atoms))
        self._scale = int(scale)

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
        return exact(value)

    __call__ = at

    def __eq__(self, other: object) -> bool:
        if isinstance(other, AlphaFrac):
            return (self._num, self._atoms, self._scale) == (
                other._num,
                other._atoms,
                other._scale,
            )
        if isinstance(other, (int, Fraction)):
            return not self._atoms and self._scale == 1 and self._num[:1] == (
                (exact(other),) if other else ()
            )
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self._num, self._atoms, self._scale))

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
        if not factors:
            return above
        if _has_top_level_sum(above):
            above = f"({above})"
        below = factors[0] if len(factors) == 1 else "(" + "*".join(factors) + ")"
        return f"{above}/{below}"


#: The Macdonald and LLT parameter `q`, as a `QtPoly`, so a coefficient can be
#: written as `q` rather than as the rows that encode it.
q: QtPoly = QtPoly({(1, 0): 1})

#: The parameter `t`, as a `QtPoly` — the one Macdonald and `H̃` take. The
#: Hall-Littlewood families are in `t` alone and take `t_hl` instead, because
#: their coefficients are a `Poly` in one variable and the two do not mix.
t: QtPoly = QtPoly({(0, 1): 1})

#: The Hall-Littlewood parameter `t`, as a one-variable `Poly`.
t_hl: Poly = Poly("t", {1: 1})

#: The Jack parameter α, as a one-variable `Poly`.
alpha: Poly = Poly("alpha", {1: 1})


class Param:
    """An element whose coefficients carry parameters, tagged with its basis.

    This is what the Macdonald, Jack, Hall-Littlewood and LLT families return.
    It is the parameter-carrying sibling of `Sym`: same basis tag, same term
    dictionary, but each coefficient is a `Poly`, `QtPoly`, `QtFrac` or
    `AlphaFrac` rather than a number.

        >>> from symfn import jack
        >>> jack.P([2])
        JackP[2]
        >>> jack.P([2]).to("m")
        2/(alpha + 1)*m[1,1] + m[2]
        >>> jack.P([2]).at(alpha=1)
        m[1,1] + m[2]

    `at` returns a `Sym`, so a specialization rejoins the ordinary arithmetic:
    the value above is `s_(2)` written in the monomial basis, which is the
    theorem `P_λ(x; 1) = s_λ` and the check that the α convention is the one
    in the literature rather than its `α → 1/α` mirror.

    The basis tag is one of the six classical codes or one of the nine
    parametric bases — `HLP` and `HLQp`; `McdHt`, `McdJ`, `McdP` and `McdQ`;
    `JackP`, `JackQ` and `JackJ`, each spelled as Sage prints it. A family's
    forward constructor returns a shape in its own basis and the inverse
    expansions land there too, so both directions speak of the same objects:

        >>> from symfn import hl, s
        >>> hl.to_P(s([2]))
        t*HLP[1,1] + HLP[2]
        >>> hl.to_P(s([2])).to("s")
        s[2]
        >>> hl.to_P(s([2])).at(t=0)
        s[2]

    `at` on a parametric basis expands through `to` first, since a `Sym`
    carries only the six classical codes; at `t = 0` the `P` basis *is* the
    Schur basis, which is why the last two agree here and not in general.
    """

    __slots__ = ("_basis", "_terms", "_params")
    __module__ = "symfn"

    def __init__(
        self,
        basis: str,
        terms: Mapping[Partition, ParamCoefficient]
        | Iterable[tuple[Partition, ParamCoefficient]],
        parameters: Sequence[str],
    ) -> None:
        """Build from a basis code, `{partition: coefficient}` rows, and the
        parameter names `at` accepts.
        """
        self._basis = check_param_basis(basis)
        items = terms.items() if hasattr(terms, "items") else terms
        self._terms = {_partition(la): c for la, c in items if c}
        self._params = tuple(parameters)

    @property
    def basis(self) -> ParamBasis:
        """The basis code the terms are indexed by: a classical one-letter
        code, or one of the nine parametric tags for an element written in a
        family's own basis.

        >>> from symfn import hl, macdonald, s
        >>> macdonald.P([2]).basis, hl.Qp([2]).basis, hl.to_Qp(s([2])).basis
        ('McdP', 'HLQp', 'HLQp')
        >>> macdonald.P([2]).to("m").basis
        'm'
        """
        return self._basis

    @property
    def parameters(self) -> tuple[str, ...]:
        """The parameter names `at` takes, in order.

        >>> from symfn import macdonald, jack
        >>> macdonald.P([2]).parameters, jack.P([2]).parameters
        (('q', 't'), ('alpha',))
        """
        return self._params

    @property
    def terms(self) -> dict[Partition, ParamCoefficient]:
        """A copy of the `{partition: coefficient}` mapping, zero-free and in
        the contract layer's element order.

            >>> from symfn import hl
            >>> hl.Qp([2]).to("s").terms
            {(2,): 1}
        """
        return dict(self._terms)

    def coefficient(self, la: PartitionArg) -> ParamCoefficient | int:
        """The coefficient of `la`, or `0` if it does not appear.

        >>> from symfn import macdonald
        >>> macdonald.P([2]).to("m").coefficient([5])
        0
        """
        return self._terms.get(_partition(la), 0)

    def support(self) -> list[Partition]:
        """The partitions carrying a nonzero coefficient, in element order.

        >>> from symfn import jack
        >>> jack.P([2]).to("m").support()
        [(1, 1), (2,)]
        """
        return list(self._terms)

    def degree(self) -> int | None:
        """The common degree of every term, or `None` if the element is not
        homogeneous — and `None` for the zero element, which has no degree.

            >>> from symfn import macdonald, m, q
            >>> macdonald.P([2]).degree()
            2
            >>> (q * m([2]) + q * m([1])).degree() is None
            True

        The degree is a fact about the partitions alone, so it needs no
        expansion out of a parametric basis and no arithmetic on the
        coefficients. `Sym.degree` returns the same value on the same shapes.
        """
        degrees = {sum(la) for la in self._terms}
        return degrees.pop() if len(degrees) == 1 else None

    def is_homogeneous(self) -> bool:
        """Whether every term has the same degree. The zero element is
        homogeneous.

            >>> from symfn import hl
            >>> (hl.Qp([2]) + hl.Qp([1, 1])).is_homogeneous()
            True
            >>> (hl.Qp([2]) + hl.Qp([1])).is_homogeneous()
            False
        """
        return len({sum(la) for la in self._terms}) <= 1

    def __neg__(self) -> Param:
        return self * -1

    def __add__(self, other: object) -> Param:
        """`f + g`, termwise, for two elements written in the same basis.

            >>> from symfn import macdonald, q, t
            >>> q * macdonald.Htilde([2, 1]) + t * macdonald.Htilde([3])
            q*McdHt[2,1] + t*McdHt[3]

        Two different bases raise rather than one being converted, on the same
        grounds `Sym` refuses: `McdP[2] + McdQ[2]` names no element.

        # Raises

        Raises `ValueError` unless both elements are in the same basis, and —
        for `McdHt` alone — if two coefficients at one shape have different
        denominators, which the boundary encoding cannot put over a common
        one.
        """
        from ._families import _add

        return _add(self, other) if isinstance(other, Param) else NotImplemented

    def __sub__(self, other: object) -> Param:
        return self + (-other) if isinstance(other, Param) else NotImplemented

    def __mul__(self, other: object) -> Param:
        """`c*f`, `c` a scalar in this element's own parameters.

            >>> from symfn import jack, alpha
            >>> alpha * jack.P([2])
            alpha*JackP[2]
            >>> (1 - alpha) * jack.P([2])
            (1 - alpha)*JackP[2]

        A scalar is an `int`, a `Fraction`, a polynomial in this element's
        parameters, or a coefficient of the kind this element carries. Two
        elements cannot be multiplied: that is a product in the ring, and a
        parametric basis has structure constants this does not compute.

        # Raises

        Raises `TypeError` if the scalar is in the wrong parameters, or if
        `other` is an element rather than a scalar.
        """
        from ._families import _scale

        if isinstance(other, Param):
            raise TypeError(
                "two elements cannot be multiplied; a product in a parametric "
                "basis needs its structure constants"
            )
        if not isinstance(other, (int, Fraction, Poly, QtPoly, QtFrac, AlphaFrac)):
            return NotImplemented
        return _scale(self, other)

    __rmul__ = __mul__

    def to(self, basis: str) -> Param:
        """The element rewritten in `basis`, one of the six classical codes.

            >>> from symfn import jack, macdonald, hl, q, m
            >>> jack.P([2]).to("m")
            2/(alpha + 1)*m[1,1] + m[2]
            >>> jack.P([2]).to("s")
            (1 - alpha)/(alpha + 1)*s[1,1] + s[2]
            >>> hl.Qp([1, 1]).to("m")
            (1 + t)*m[1,1] + t*m[2]
            >>> (q * m([2])).to("s")
            -q*s[1,1] + q*s[2]
            >>> macdonald.to_P(macdonald.P([2]).to("m")) == macdonald.P([2])
            True

        A parametric basis expands into one classical basis — monomial for the
        Macdonald and Jack normalizations, Schur for Hall-Littlewood and `H̃` —
        because that is the basis its family's forward direction is defined in.
        Any other classical basis is that expansion followed by an ordinary
        change of basis, which is a ℤ-linear map on the partitions and so
        carries the coefficients through untouched.

        The Jack value vanishes at α = 1, where `P_λ` is the Schur function,
        which distinguishes this convention from its `α → 1/α` mirror.

        The result is a `Param`, not a `Sym`: the coefficients still carry `q`,
        `t` or α.

        # Raises

        Raises `ValueError` unless `basis` is one of `s`, `h`, `e`, `m`, `f`.
        The power-sum basis is not reachable: that conversion divides by z_μ,
        and the rings these coefficients live in are not all closed under it.
        """
        from ._families import _convert, _expand

        basis = check_basis(basis)
        if basis == "p":
            raise ValueError(
                "the power-sum basis is not reachable from a coefficient that "
                "carries a parameter; that conversion divides by z_mu"
            )
        if self._basis in BASES:
            return _convert(self, basis)
        expanded = _expand(self)
        return expanded if expanded.basis == basis else _convert(expanded, basis)

    def omega(self) -> Param:
        """The ω involution, returned in this element's basis.

            >>> from symfn import hl, q, m
            >>> hl.Qp([2]).omega()
            HLQp[1,1] - t*HLQp[2]
            >>> (q * m([2, 1])).omega()
            -q*m[2,1] - 2*q*m[3]

        ω is its own inverse and exchanges `e` with `h`. It is defined on the
        Schur basis, where it conjugates the shape, and reaches any other by a
        change of basis on each side — so an element in a family's own basis
        comes back in it. The second value is `m([2, 1]).omega()` scaled by
        `q`, which is the check that the parameter is carried and not acted on.

        # Raises

        Raises `ValueError` for the coefficient classes `to` declines, and for
        a parametric basis whose inverse expansion runs over them.
        """
        from ._families import _hopf

        return _hopf(self, "omega")

    def antipode(self) -> Param:
        """The antipode S of the Hopf algebra, in this element's basis.

            >>> from symfn import q, m
            >>> (q * m([2, 1])).antipode()
            q*m[2,1] + 2*q*m[3]

        `S(s_λ) = (−1)^|λ| s_{λ'}`, which is ω up to that sign. The shape here
        has odd degree, so this value is the negative of `omega`'s and the two
        are told apart; at even degree every convention agrees.

        # Raises

        Raises `ValueError` for the same reasons `omega` does.
        """
        from ._families import _hopf

        return _hopf(self, "antipode")

    def at(self, *args: Coefficient, **kwargs: Coefficient) -> Sym:
        """The element with its parameters set, as a `Sym` in the same basis.

            >>> from symfn import hl, macdonald
            >>> macdonald.P([2]).at(q=7, t=7)
            m[1,1] + m[2]
            >>> hl.Qp([2, 1]).at(t=0)
            s[2,1]

        The first is `P_λ(x; q, q) = s_λ` written in the monomial basis; the
        second is `Q'_λ(x; 0) = s_λ`. Both are theorems, and both fail under a
        `q ↔ t` or `t → 1/t` twist of the convention.

        A `Sym` carries only the six classical bases, so an element in a
        parametric one is expanded first — through `to`, into the basis its
        family is written in — and specialized there.

        # Raises

        Raises `TypeError` unless exactly this element's `parameters` are
        supplied, by name or in that order. Raises `ValueError` on what `to`
        raises for, when the element is in a parametric basis.
        """
        from ._families import EXPANDS_IN
        from ._sym import Sym

        if self._basis not in BASES:
            return self.to(EXPANDS_IN[self._basis]).at(*args, **kwargs)
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
        return Sym(self._basis, {la: c.at(*ordered) for la, c in self._terms.items()})

    def __len__(self) -> int:
        return len(self._terms)

    def __bool__(self) -> bool:
        return bool(self._terms)

    def __iter__(self) -> Iterator[tuple[Partition, ParamCoefficient]]:
        return iter(self._terms.items())

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Param):
            return (
                self._basis == other._basis
                and self._params == other._params
                and self._terms == other._terms
            )
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self._basis, self._params, frozenset(self._terms.items())))

    def __repr__(self) -> str:
        if not self._terms:
            return "0"
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
