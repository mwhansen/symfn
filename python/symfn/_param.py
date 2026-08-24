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


#: What may stand in for an element in `Param`'s arithmetic: a scalar of the
#: base ring, or an integer or rational that injects into it. The same set for
#: `+`, `-` and `*`, so a scalar that can multiply an element can add to it.
_SCALARS = (int, Fraction, Poly, QtPoly, QtFrac, AlphaFrac)


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
        # Sorted, as `Sym` sorts: the order is what `repr` and iteration show,
        # and leaving it as the caller inserted made `f + g` and `g + f` print
        # differently for the two coefficient classes this layer adds itself.
        self._terms = dict(
            sorted((_partition(la), c) for la, c in items if c)
        )
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
        grounds `Sym` refuses: `McdP[2] + McdQ[2]` names no element. Two
        different base rings raise as well, and say so separately, because
        `.to()` cannot fix that one.

        **A scalar adds**, as the constant it names times the unit, which the
        empty partition indexes in every basis here:

            >>> from symfn import jack
            >>> 2 + jack.P([1])
            2 + JackP[1]
            >>> sum([jack.P([1]), jack.P([2])])
            JackP[1] + JackP[2]

        The second is why: `sum` starts from `0`, so without this it raises on
        the first term.

        # Raises

        Raises `BasisError` unless both elements are in the same basis, and
        `BaseRingError` unless both are over the same base ring — the same two
        exceptions `Sym` raises for the same two questions. Raises
        `ValueError` for `McdHt` alone if two coefficients at one shape have
        different denominators, which the boundary encoding cannot put over a
        common one.
        """
        from ._families import _add, _constant

        if isinstance(other, Param):
            return _add(self, other)
        if isinstance(other, _SCALARS):
            return _add(self, _constant(self, other))
        return NotImplemented

    __radd__ = __add__

    def __sub__(self, other: object) -> Param:
        from ._families import _add, _constant

        if isinstance(other, Param):
            return _add(self, -other)
        if isinstance(other, _SCALARS):
            return _add(self, -_constant(self, other))
        return NotImplemented

    def __rsub__(self, other: object) -> Param:
        from ._families import _add, _constant

        if isinstance(other, _SCALARS):
            return _add(-self, _constant(self, other))
        return NotImplemented

    def __mul__(self, other: object) -> Param:
        """`c*f`, `c` a scalar in this element's own parameters.

            >>> from symfn import jack, alpha
            >>> alpha * jack.P([2])
            alpha*JackP[2]
            >>> (1 - alpha) * jack.P([2])
            (1 - alpha)*JackP[2]

        A scalar is an `int`, a `Fraction`, a polynomial in this element's
        parameters, or a coefficient of the kind this element carries.

        **Two elements multiply**, in the basis both are written in. A
        parametric basis is a basis of the ring like any other: the product is
        taken where the family expands and rewritten back by the inverse
        expansion. All nine tags multiply.

            >>> from symfn import hl, jack
            >>> hl.P([1]) * hl.P([1])
            (1 + t)*HLP[1,1] + HLP[2]
            >>> jack.P([1]) * jack.P([1])
            2*alpha/(alpha + 1)*JackP[1,1] + JackP[2]
            >>> (hl.P([2, 1]) * hl.P([2, 1])).coefficient([3, 1, 1, 1])
            1 + t - t^3 - t^4

        The second value is what pins the normalization, and the first is not:
        `P_1² = P_2 + (1 + t)·P_11` is what every convention in circulation
        gives. The second has **negative** coefficients, so these structure
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
        from ._families import _product, _scale

        if isinstance(other, Param):
            return _product(self, other)
        if not isinstance(other, (int, Fraction, Poly, QtPoly, QtFrac, AlphaFrac)):
            return NotImplemented
        return _scale(self, other)

    __rmul__ = __mul__

    def __pow__(self, n: int) -> Param:
        """A non-negative integer power, by repeated squaring.

            >>> from symfn import hl
            >>> hl.P([1]) ** 2
            (1 + t)*HLP[1,1] + HLP[2]
            >>> hl.P([1]) ** 0
            1

        # Raises

        Raises `ValueError` unless `n` is a non-negative `int`; there is no
        inverse in this ring. Raises what `*` raises for the coefficient
        classes it does not carry.
        """
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

            >>> from symfn import hl, jack, q, m
            >>> hl.Qp([2]).omega()
            HLQp[1,1] - t*HLQp[2]
            >>> (q * m([2, 1])).omega()
            -q*m[2,1] - 2*q*m[3]
            >>> jack.P([2]).omega()
            4*alpha/(alpha + 1)^2*JackP[1,1] + (1 - alpha)/(alpha + 1)*JackP[2]

        ω is its own inverse and exchanges `e` with `h`. It is defined on the
        Schur basis, where it conjugates the shape, and reaches any other by a
        change of basis on each side — so an element in a family's own basis
        comes back in it. The second value is `m([2, 1]).omega()` scaled by
        `q`, which is the check that the parameter is carried and not acted on.

        The Jack value pins which ω this is. It is the plain involution,
        `p_r ↦ (−1)^{r−1} p_r`, so α is carried and not touched; the
        α-deformed one sends `P_λ^{(α)}` to `Q_{λ'}^{(1/α)}`, which inverts
        the parameter. Sage's `.omega()` gives this same value.

        # Raises

        Raises `ValueError` if the element carries no coefficient class this
        layer knows, and for a parametric basis whose inverse expansion runs
        over one it does not.
        """
        from ._families import _hopf

        return _hopf(self, "omega")

    def internal_product(self, g: Sym | Param) -> Param:
        """The internal (Kronecker) product, in this element's basis.

            >>> from symfn import hl, jack, macdonald
            >>> jack.P([1]).internal_product(jack.P([1]))
            JackP[1]
            >>> macdonald.P([2]).internal_product(macdonald.P([1, 1]))
            (1 - t^2 - q^2 + q^2*t^2)/(1 - q*t)^2*McdP[1,1] + (-t + q)/(1 - q*t)*McdP[2]
            >>> hl.P([2, 1]).internal_product(hl.P([2, 1])).coefficient([3])
            1 + t^2 + 2*t^3 + t^4

        All three are Sage's values. The structure constants are Kronecker
        coefficients, which are integers carrying no parameter, so the
        parameters are multiplied through rather than acted on. At `t = 0` the
        third is 1, one term of `s_21 ∗ s_21 = s_111 + s_21 + s_3`, since
        `P_λ(x; 0) = s_λ`.

        Two elements combine here, unlike `scalar` and `skew_by`, so they must
        be in the same basis, on the same grounds `*` refuses. A `Sym` is
        lifted into this element's base ring rather than refused.

        # Raises

        Raises `BasisError` unless both are in the same basis, `BaseRingError`
        unless both are over the same base ring, and `ValueError` if the
        power-sum route produces a non-integral coefficient.
        """
        from ._families import _internal

        return _internal(self, g)

    def dimension(self) -> ParamCoefficient:
        """The dimension `Σ c_λ f^λ`, with `f^λ` the standard-tableaux count.

            >>> from symfn import hl, jack
            >>> hl.P([2, 1]).dimension()
            2 - t - t^2
            >>> jack.P([2, 1]).dimension()
            6/(alpha + 2)

        `f^λ` depends on the shape alone, so the parameters are carried and
        never acted on. Both degenerate to `f^{21} = 2`, the first at `t = 0`
        and the second at α = 1.

        # Raises

        Raises `ValueError` if this element carries no coefficient class this
        layer knows, and `OverflowError` if a term's `f^λ` exceeds `u128`.
        """
        from ._families import _DIMENSION, _functional

        return _functional(self, _DIMENSION, "the dimension")

    def principal_specialization(self, n: int) -> ParamCoefficient:
        """The value at `1^n`, summed over the Schur expansion.

            >>> from symfn import hl
            >>> hl.P([2, 1]).principal_specialization(3)
            8 - t - t^2

        The same value `evaluate([1, 1, 1])` gives, by a different route: this
        one weighs each shape by `s_λ(1^n)` and never lays out an alphabet.

        # Raises

        Raises `ValueError` if this element carries no coefficient class this
        layer knows, and `OverflowError` if a term's value exceeds the
        fixed-width specialization.
        """
        from ._families import _PRINCIPAL, _functional

        return _functional(self, _PRINCIPAL, "the value at 1^n", n)

    def principal_specialization_q(self, n: int) -> QtPoly:
        """The value at `1, q, …, q^{n−1}`, as a `QtPoly`.

            >>> from symfn import hl
            >>> hl.P([2, 1]).principal_specialization_q(3)
            q + 2*q^2 + 2*q^3 - q^3*t - q^3*t^2 + 2*q^4 + q^5
            >>> hl.P([2, 1]).principal_specialization_q(3).at(1, 0)
            8

        Sage writes the first as `q^5 + 2q^4 + (−t² − t + 2)q³ + 2q² + q`, the
        same value. Setting `q = 1` gives `principal_specialization(3)`, and
        setting `t = 0` gives `s_21(1,q,q²)`; the second value does both at
        once and is `s_21(1,1,1) = 8`. That the two agree is the check that
        the `q` this introduces and the `t` already there stayed apart.

        **Only where the base ring leaves room for `q`.** This layer's
        coefficient classes carry at most two variables, so an element already
        over `ℚ(q,t)` or `ℚ(α)` has nowhere to put the new one and is refused —
        the wall Sage reports as "the variable q is in the base ring, pass it
        explicitly". Use `evaluate` with an alphabet you name yourself instead.

        # Raises

        Raises `ValueError` if the base ring already carries `q`, and for the
        coefficient classes with no free variable at all.
        """
        from ._families import _principal_q

        return _principal_q(self, n)

    def expand(self, n: int) -> dict[tuple[int, ...], ParamCoefficient]:
        """The expansion in `n` variables, as a `{exponent vector: coefficient}`
        mapping with every vector of length `n`.

            >>> from symfn import macdonald
            >>> macdonald.P([2]).expand(2)[1, 1]
            (1 - t + q - q*t)/(1 - q*t)

        Sage writes that coefficient of `x0*x1` as `(q*t-q+t-1)/(q*t-1)`, the
        same after clearing signs. The `q ↔ t` swap gives
        `(1 - q + t - q*t)/(1 - q*t)`, which is the twist to check.

        The order is the contract layer's: grouped by monomial term, and
        unspecified within a group.

        # Raises

        Raises `ValueError` if this element carries no coefficient class this
        layer knows.
        """
        from ._families import _expand_alphabet

        return _expand_alphabet(self, n)

    def evaluate(self, xs: Sequence[Coefficient]) -> ParamCoefficient:
        """The value at the alphabet `xs`, a sequence of integers.

            >>> from symfn import hl, jack
            >>> hl.P([2, 1]).evaluate([1, 1, 1])
            8 - t - t^2
            >>> jack.P([2, 1]).evaluate([1, 1, 1])
            (18 + 6*alpha)/(alpha + 2)

        Both are Sage's values. The alphabet injects into the base ring, so it
        meets the shape and not the parameters, which ride through. Both
        degenerate to `s_21(1,1,1) = 8` — the first at `t = 0` and the second
        at α = 1, which are where each family is the Schur basis.

        # Raises

        Raises `ValueError` if this element carries no coefficient class this
        layer knows.
        """
        from ._families import _evaluate

        return _evaluate(self, xs)

    def coproduct(self) -> dict[tuple[Partition, Partition], ParamCoefficient]:
        """The coproduct Δ, as a `{(mu, nu): coefficient}` mapping over the
        Schur basis of each factor.

            >>> from symfn import hl, jack
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
        both factors come back in the Schur basis rather than the basis this
        element is written in.

        `Δ(s_λ) = Σ c^λ_{μν} s_μ ⊗ s_ν` with Littlewood-Richardson
        coefficients, so the parameters are carried and never acted on.

        # Raises

        Raises `ValueError` if this element carries no coefficient class this
        layer knows.
        """
        from ._families import _coproduct

        return _coproduct(self)

    def scalar(self, g: Sym | Param) -> ParamCoefficient | Coefficient:
        """The Hall inner product `⟨self, g⟩`, as one coefficient.

            >>> from symfn import macdonald, jack, s
            >>> macdonald.P([2]).scalar(macdonald.P([1, 1]))
            (-t + q)/(1 - q*t)
            >>> macdonald.P([2]).scalar(s([2]))
            1
            >>> jack.P([2, 1]).scalar(jack.P([2, 1]))
            (8 - 4*alpha + 5*alpha^2)/(alpha + 2)^2

        Both are Sage's values. The Schur basis is orthonormal for this
        pairing, so the value is the sum of the products of matching
        coefficients — bilinear over whatever ring they live in, which is why
        the parameters are carried and never acted on.

        `g` may be in any basis and may be a `Sym`. The pairing is defined on
        the ring, so two spellings of one argument give one number and there is
        nothing to refuse; `h([2]).scalar(m([2]))` is 1, a value rather than a
        mismatch.

        # Raises

        Raises `ValueError` if this element carries no coefficient class this
        layer knows, and `BaseRingError` if `g` is over a different base ring.
        """
        from ._families import _scalar

        return _scalar(self, g)

    def skew_by(self, g: Sym | Param) -> Param:
        """The element skewed by `g`, the adjoint of multiplication by `g`
        under the Hall inner product, in this element's basis.

            >>> from symfn import macdonald, hl, s, h
            >>> macdonald.P([2, 1]).skew_by(s([1]))
            (1 - t^2 - q^2*t + q^2*t^3)/((1 - q*t)*(1 - q*t^2))*McdP[1,1] + McdP[2]
            >>> hl.P([2, 1]).skew_by(h([1]))
            (1 - t^2)*HLP[1,1] + HLP[2]

        Both are Sage's values. `g` keeps the basis it is written in, because
        that basis selects which rule runs and not merely how `g` is read:
        `h`, `e` and `p` take the Pieri, dual-Pieri and Murnaghan-Nakayama
        paths, and `s`, `m` and `f` go through Littlewood-Richardson. All six
        have integer structure constants, so the parameters are carried and
        never acted on.

        `g` may be a `Sym`, which is lifted into this element's base ring
        rather than refused.

        # Raises

        Raises `ValueError` if this element carries no coefficient class this
        layer knows, and `BaseRingError` if `g` is over a different base ring.
        """
        from ._families import _skew

        return _skew(self, g)

    def antipode(self) -> Param:
        """The antipode S of the Hopf algebra, in this element's basis.

            >>> from symfn import macdonald, q, m
            >>> (q * m([2, 1])).antipode()
            q*m[2,1] + 2*q*m[3]
            >>> macdonald.P([1]).antipode()
            -McdP[1]

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
