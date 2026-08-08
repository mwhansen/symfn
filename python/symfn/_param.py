"""The coefficient types the parameter families return, and `Param` over them.

The contract layer hands a parameter family back as exponent-keyed rows: a
Hall-Littlewood coefficient is `[(t_exponent, coefficient), ...]`, a Macdonald
one is a numerator in `q` and `t` over a *factored* denominator, a Jack one is a
dense numerator in α over factored atoms. That encoding is deliberate — it
assumes no coefficient ring on the far side, which is what lets Sage, SymPy and
a bare interpreter each rebuild it in their own
(``docs/policies/python.md``, P1) — and it is not what a person wants to read.

The types here are that encoding with a `repr` and an evaluation map. They hold
the rows they were given, unchanged; `at` substitutes numbers into them. That
is arithmetic on plain data rather than symmetric-function mathematics, so the
"convenience computes nothing" rule (``docs/policies/python.md``, P4) is not
bent by it: nothing here can produce a coefficient the contract layer did not.

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

from ._bases import check_basis, exact
from ._sym import _partition
from ._types import Basis, Coefficient, Partition, PartitionArg

if TYPE_CHECKING:
    from ._sym import Sym

#: A coefficient that carries parameters, as `Param` holds them.
ParamCoefficient = Union["Poly", "QtPoly", "QtFrac", "AlphaFrac"]

__all__ = ["Poly", "QtPoly", "QtFrac", "AlphaFrac", "Param"]


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
        >>> macdonald.Q([1]).coefficient([1])
        (1 - t)/(1 - q)
        >>> macdonald.P([1]).coefficient([1])
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
        >>> macdonald.Q([1]).coefficient([1]).numerator
        1 - t
        """
        return self._num

    @property
    def denominator(self) -> tuple[tuple[int, int, int], ...]:
        """The denominator's factors, as `(q, t, multiplicity)` triples for
        `(1 − q^a t^b)^multiplicity`. Empty when the coefficient is a
        polynomial, which is the integral form's signature.

            >>> from symfn import macdonald
            >>> macdonald.Q([1]).coefficient([1]).denominator
            ((1, 0, 1),)
            >>> macdonald.J([1, 1]).coefficient([1, 1]).denominator
            ()
        """
        return self._den

    def at(self, q: Coefficient, t: Coefficient) -> Coefficient:
        """The value at `q` and `t`, exactly.

        >>> from fractions import Fraction
        >>> from symfn import macdonald
        >>> macdonald.Q([1]).coefficient([1]).at(q=0, t=Fraction(1, 2))
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


class AlphaFrac:
    """A Jack coefficient: a polynomial in α over factored linear atoms.

    The value is `(Σ_k numerator[k]·α^k) / (scale · ∏ (u·α + v)^multiplicity)`,
    with the numerator dense in the α-exponent and the atoms primitive
    (`gcd(u, v) = 1`), so the factorization is canonical.

        >>> from symfn import jack
        >>> jack.P([2]).coefficient([1, 1])
        2/(alpha + 1)
        >>> jack.P([2]).coefficient([1, 1]).at(1)
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
            >>> jack.P([2]).coefficient([1, 1]).numerator
            (2,)
        """
        return self._num

    @property
    def atoms(self) -> tuple[tuple[int, int, int], ...]:
        """The denominator's atoms, as `(u, v, multiplicity)` for
        `(u·α + v)^multiplicity`.

            >>> from symfn import jack
            >>> jack.P([2]).coefficient([1, 1]).atoms
            ((1, 1, 1),)
        """
        return self._atoms

    @property
    def scale(self) -> int:
        """The integer the denominator also carries.

        >>> from symfn import jack
        >>> jack.P([2]).coefficient([1, 1]).scale
        1
        """
        return self._scale

    def at(self, alpha: Coefficient) -> Coefficient:
        """The value at `alpha`, exactly.

        >>> from symfn import jack
        >>> jack.P([3, 1]).coefficient([2, 1, 1]).at(2)
        Fraction(9, 5)

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


class Param:
    """An element whose coefficients carry parameters, tagged with its basis.

    This is what the Macdonald, Jack, Hall-Littlewood and LLT families return.
    It is the parameter-carrying sibling of `Sym`: same basis tag, same term
    dictionary, but each coefficient is a `Poly`, `QtPoly`, `QtFrac` or
    `AlphaFrac` rather than a number.

        >>> from symfn import jack
        >>> jack.P([2])
        2/(alpha + 1)*m[1,1] + m[2]
        >>> jack.P([2]).at(alpha=1)
        m[1,1] + m[2]

    `at` returns a `Sym`, so a specialization rejoins the ordinary arithmetic:
    the value above is `s_(2)` written in the monomial basis, which is the
    theorem `P_λ(x; 1) = s_λ` and the check that the α convention is the one
    in the literature rather than its `α → 1/α` mirror.
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
        self._basis = check_basis(basis)
        items = terms.items() if hasattr(terms, "items") else terms
        self._terms = {_partition(la): c for la, c in items if c}
        self._params = tuple(parameters)

    @property
    def basis(self) -> Basis:
        """The one-letter basis code the terms are indexed by.

        >>> from symfn import hl, macdonald
        >>> macdonald.P([2]).basis, hl.Qp([2]).basis
        ('m', 's')
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
            >>> hl.Qp([2]).terms
            {(2,): 1}
        """
        return dict(self._terms)

    def coefficient(self, la: PartitionArg) -> ParamCoefficient | int:
        """The coefficient of `la`, or `0` if it does not appear.

        >>> from symfn import macdonald
        >>> macdonald.P([2]).coefficient([5])
        0
        """
        return self._terms.get(_partition(la), 0)

    def support(self) -> list[Partition]:
        """The partitions carrying a nonzero coefficient, in element order.

        >>> from symfn import jack
        >>> jack.P([2]).support()
        [(1, 1), (2,)]
        """
        return list(self._terms)

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

        # Raises

        Raises `TypeError` unless exactly this element's `parameters` are
        supplied, by name or in that order.
        """
        from ._sym import Sym

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
