"""`Schub`, a Schubert polynomial keyed by permutations rather than partitions.

Same shape as `Sym` — an immutable term dictionary, arithmetic through the
contract layer, nothing computed here — over a different index set. A support
is a permutation in one-line notation, so the two types never mix: `Sym` and
`Schub` have no common basis and no operation takes one of each.

Sage's equivalent is `SchubertPolynomialRing(ZZ)`.
"""

from __future__ import annotations

from collections.abc import Iterator
from typing import Union

from . import symfn as _c
from ._bases import exact
from ._types import Permutation, PermutationArg, SchubertTermsArg

#: What a binary operation accepts beside another Schubert polynomial.
Operand = Union["Schub", int]

__all__ = ["Schub", "X"]


def _permutation(w: PermutationArg) -> Permutation:
    """Normalize a permutation argument to a tuple with trailing fixed points
    dropped.

        >>> _permutation([1, 3, 2, 4])
        (1, 3, 2)
        >>> _permutation([1, 2, 3])
        ()
        >>> _permutation(1)
        ()

    # Raises

    Raises `ValueError` unless the argument is a permutation of `1..n`.
    """
    if isinstance(w, int):
        w = (w,)
    w = tuple(int(x) for x in w)
    if sorted(w) != list(range(1, len(w) + 1)):
        raise ValueError(f"not a permutation of 1..{len(w)}: {w!r}")
    while w and w[-1] == len(w):
        w = w[:-1]
    return w


class Schub:
    """A Schubert polynomial, as permutation-keyed terms.

        >>> from symfn import X
        >>> X[1, 3, 2] * X[1, 3, 2]
        X[1,4,2,3] + X[2,3,1]

    The identity permutation indexes 1, and is written as a bare coefficient in
    the `repr`, as `Sym` writes the empty partition. Terms are held zero-free
    and in the contract layer's order.

    # Raises

    Raises `TypeError` from `+`, `-` and `*` when the other operand is neither
    a `Schub` nor an integer.
    """

    __slots__ = ("_terms",)
    __module__ = "symfn"

    def __init__(self, terms: SchubertTermsArg) -> None:
        """Build from a `{permutation: coefficient}` mapping or
        `(permutation, coefficient)` rows.

        # Raises

        Raises `ValueError` unless every key is a permutation, and `TypeError`
        unless every coefficient is an `int`.
        """
        items = terms.items() if hasattr(terms, "items") else terms
        collected: dict[Permutation, int] = {}
        for key, value in items:
            w = _permutation(key)
            c = exact(value)
            if not isinstance(c, int):
                raise TypeError(
                    "a Schubert coefficient is an integer; the basis has no "
                    f"denominators, and {c!r} has one"
                )
            c += collected.get(w, 0)
            if c:
                collected[w] = c
            else:
                collected.pop(w, None)
        self._terms: dict[Permutation, int] = dict(sorted(collected.items()))

    @property
    def terms(self) -> dict[Permutation, int]:
        """A copy of the `{permutation: coefficient}` mapping, zero-free.

        >>> from symfn import X
        >>> X[1, 3, 2].terms
        {(1, 3, 2): 1}
        """
        return dict(self._terms)

    def coefficient(self, w: PermutationArg) -> int:
        """The coefficient of `w`, or `0` if it does not appear.

        >>> from symfn import X
        >>> (X[1, 3, 2] * X[1, 3, 2]).coefficient([2, 3, 1])
        1
        """
        return self._terms.get(_permutation(w), 0)

    def support(self) -> list[Permutation]:
        """The permutations carrying a nonzero coefficient, in order.

        >>> from symfn import X
        >>> (X[1, 3, 2] * X[1, 3, 2]).support()
        [(1, 4, 2, 3), (2, 3, 1)]
        """
        return list(self._terms)

    def __len__(self) -> int:
        return len(self._terms)

    def __bool__(self) -> bool:
        return bool(self._terms)

    def __iter__(self) -> Iterator[tuple[Permutation, int]]:
        return iter(self._terms.items())

    def _same(self, other: object, op: str) -> Schub:
        if isinstance(other, Schub):
            return other
        if isinstance(other, int):
            unit: dict[Permutation, int] = {(): other} if other else {}
            return Schub(unit)
        raise TypeError(
            f"cannot {op} {type(other).__name__} with a Schubert polynomial"
        )

    def __add__(self, other: Operand) -> Schub:
        """Add two Schubert polynomials.

        >>> from symfn import X
        >>> X[1, 3, 2] + X[1, 3, 2]
        2*X[1,3,2]
        """
        other = self._same(other, "add")
        out = dict(self._terms)
        for w, c in other._terms.items():
            out[w] = out.get(w, 0) + c
        return Schub(out)

    __radd__ = __add__

    def __neg__(self) -> Schub:
        return Schub({w: -c for w, c in self._terms.items()})

    def __sub__(self, other: Operand) -> Schub:
        """Subtract two Schubert polynomials.

        >>> from symfn import X
        >>> X[2, 1] - X[2, 1]
        0
        """
        return self + (-self._same(other, "subtract"))

    def __rsub__(self, other: Operand) -> Schub:
        return (-self) + other

    def __mul__(self, other: Operand) -> Schub:
        """Multiply, expanding back into the Schubert basis.

            >>> from symfn import X
            >>> X[2, 1] * X[2, 1]
            X[3,1,2]

        The structure constants are the generalized Littlewood-Richardson
        coefficients; `𝔖_21 · 𝔖_21 = 𝔖_312` rather than `𝔖_231` is what
        distinguishes this from the convention indexed by inverse
        permutations.
        """
        if isinstance(other, int):
            return Schub({w: c * other for w, c in self._terms.items()})
        other = self._same(other, "multiply")
        if not self._terms or not other._terms:
            return Schub({})
        return Schub(_c.schubert_multiply(list(self), list(other)))

    __rmul__ = __mul__

    def __pow__(self, n: int) -> Schub:
        """A non-negative integer power.

        >>> from symfn import X
        >>> X[2, 1] ** 2
        X[3,1,2]

        # Raises

        Raises `ValueError` unless `n` is a non-negative `int`.
        """
        if not isinstance(n, int) or n < 0:
            raise ValueError(f"exponent must be a non-negative int, not {n!r}")
        one: dict[Permutation, int] = {(): 1}
        out, base = Schub(one), self
        while n:
            if n & 1:
                out = out * base
            n >>= 1
            if n:
                base = base * base
        return out

    def divided_difference(self, i: int) -> Schub:
        """The divided difference operator `∂_i`.

            >>> from symfn import X
            >>> X[1, 3, 2].divided_difference(2)
            1

        `∂_i 𝔖_w` is `𝔖_{w s_i}` when `w` descends at `i` and 0 otherwise,
        which the value above pins.

        # Raises

        Raises `ValueError` unless `i` is a positive integer.
        """
        return Schub(_c.schubert_divided_difference(list(self), i))

    def expand(self) -> dict[tuple[int, ...], int]:
        """The expansion as a polynomial, as a `{exponent vector: coefficient}`
        mapping.

            >>> from symfn import X
            >>> X[1, 3, 2].expand()
            {(0, 1): 1, (1,): 1}

        The vectors are ragged: trailing zeros are dropped, so `(1,)` means
        `x_1`.
        """
        return dict(_c.schubert_expand(list(self)))

    def __eq__(self, other: object) -> bool:
        if isinstance(other, Schub):
            return self._terms == other._terms
        if isinstance(other, int):
            return self._terms == ({(): other} if other else {})
        return NotImplemented

    def __hash__(self) -> int:
        return hash(frozenset(self._terms.items()))

    def __repr__(self) -> str:
        if not self._terms:
            return "0"
        pieces = []
        for w, c in self._terms.items():
            negative = c < 0
            mag = -c if negative else c
            atom = f"X[{','.join(map(str, w))}]" if w else ""
            if not atom:
                body = str(mag)
            elif mag == 1:
                body = atom
            else:
                body = f"{mag}*{atom}"
            pieces.append(("-" if negative else "+", body))
        sign, body = pieces[0]
        out = ("-" if sign == "-" else "") + body
        for sign, body in pieces[1:]:
            out += f" {sign} {body}"
        return out


class _SchubFactory:
    """The Schubert basis, callable and subscriptable.

    >>> from symfn import X
    >>> X[1, 3, 2] == X([1, 3, 2])
    True
    >>> X()
    1
    """

    __slots__ = ()
    __module__ = "symfn"

    def __call__(self, w: PermutationArg = ()) -> Schub:
        """`𝔖_w`, the Schubert polynomial of `w`; with no argument, 1.

        # Raises

        Raises `ValueError` unless `w` is a permutation of `1..n`.
        """
        return Schub({_permutation(w): 1})

    def __getitem__(self, w: PermutationArg) -> Schub:
        """`𝔖_w`, written as a subscript."""
        return Schub({_permutation(w): 1})

    def __repr__(self) -> str:
        return "X"


#: The Schubert basis.
X = _SchubFactory()
