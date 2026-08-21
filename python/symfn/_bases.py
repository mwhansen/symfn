"""Basis codes and coefficient scaling for the convenience layer.

The convenience layer identifies the six classical bases by one-letter code:
``s``, ``h``, ``e``, ``p``, ``m``, ``f``. The contract layer accepts the same
codes wherever it takes a basis argument, so the codes pass straight through;
``check_basis`` validates one and ``BASES`` lists them.

The contract layer accepts only integer coefficients, but the convenience
layer holds rational ones. ``clear_denominators`` multiplies an element by
the least common denominator of its coefficients so it can be passed across,
and ``restore`` divides the result back. Every contract entry point the
convenience layer calls is ℚ-linear, so this round trip is exact.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from fractions import Fraction
from math import gcd

from ._types import Basis, Coefficient, ParamBasis, Partition

__all__ = [
    "BasisError",
    "BASES",
    "PARAM_BASES",
    "basis_name",
    "check_basis",
    "check_param_basis",
]

#: The six basis codes, in the order error messages list them.
BASES: tuple[Basis, ...] = ("s", "h", "e", "p", "m", "f")

#: The codes a `Param` accepts: the six, then the parametric bases the inverse
#: expansions return in, spelled as Sage prints them.
PARAM_BASES: tuple[ParamBasis, ...] = (*BASES, "HLP", "HLQp", "McdHt")

#: One-letter code to the name a human reads.
LONG: dict[Basis, str] = {
    "s": "Schur",
    "h": "homogeneous",
    "e": "elementary",
    "p": "power-sum",
    "m": "monomial",
    "f": "forgotten",
}


class BasisError(TypeError):
    """Raised when two elements in different bases are combined.

    The convenience layer refuses rather than converting, because a silent
    conversion picks a basis for the result that the caller did not choose and
    hides its cost. Convert explicitly with `Sym.to`:

        >>> from symfn import s, h
        >>> s([2, 1]) + h([2])
        Traceback (most recent call last):
          ...
        symfn.BasisError: cannot combine s with h; convert one with .to()
        >>> s([2, 1]) + h([2]).to("s")
        s[2] + s[2,1]

    It subclasses `TypeError` so a caller who wraps arithmetic in
    ``except TypeError`` still catches it.
    """

    # The public name is `symfn.BasisError`, and that is what a traceback
    # should print; without this it reads `symfn._bases.BasisError`.
    __module__ = "symfn"


def check_basis(code: str) -> Basis:
    """Return `code` if it names a basis, raising `ValueError` otherwise.

        >>> check_basis("s")
        's'
        >>> check_basis("Schur")
        Traceback (most recent call last):
          ...
        ValueError: unknown basis 'Schur'; expected one of s, h, e, p, m, f

    This is the one place a plain `str` becomes a `Basis`: the membership
    test against `BASES` is what narrows the type, so every string from
    outside the package passes through here first.

    # Raises

    Raises `ValueError` unless `code` is one of `s`, `h`, `e`, `p`, `m`, `f`.
    """
    if code not in BASES:
        raise ValueError(
            f"unknown basis {code!r}; expected one of " + ", ".join(BASES)
        )
    return code


def check_param_basis(code: str) -> ParamBasis:
    """Return `code` if it names a basis a `Param` can carry, raising
    `ValueError` otherwise.

        >>> check_param_basis("HLP")
        'HLP'
        >>> check_param_basis("P")
        Traceback (most recent call last):
          ...
        ValueError: unknown basis 'P'; expected one of s, h, e, ... HLQp, McdHt

    The parametric codes are spelled as Sage prints them, so `HLP[2,1]` in a
    `repr` here and `HLP[2, 1]` in Sage name the same element.

    # Raises

    Raises `ValueError` unless `code` is one of `PARAM_BASES`.
    """
    if code not in PARAM_BASES:
        raise ValueError(
            f"unknown basis {code!r}; expected one of " + ", ".join(PARAM_BASES)
        )
    return code


def basis_name(code: str) -> str:
    """The human-readable name of a basis code.

    >>> basis_name("p")
    'power-sum'
    """
    return LONG[check_basis(code)]


def clear_denominators(
    terms: Mapping[Partition, Coefficient],
) -> tuple[list[tuple[Partition, int]], int]:
    """Scale rational coefficients to integers, returning `(pairs, scale)`.

    `pairs` is the ``(partition, integer)`` list a contract entry point
    accepts, and the element it stands for is `pairs` divided by `scale`. A
    caller passes `pairs` across the boundary and hands the result and `scale`
    to `restore`.

        >>> from fractions import Fraction
        >>> clear_denominators({(2,): Fraction(1, 2), (1, 1): Fraction(1, 3)})
        ([((2,), 3), ((1, 1), 2)], 6)

    An all-integer element comes back with `scale` 1 and no copying of values:

        >>> clear_denominators({(2,): 5})
        ([((2,), 5)], 1)
    """
    scale = 1
    for c in terms.values():
        if isinstance(c, Fraction):
            d = c.denominator
            scale = scale * d // gcd(scale, d)
    return [(la, int(c * scale)) for la, c in terms.items()], scale


def restore(
    pairs: Iterable[tuple[Partition, Coefficient]], scale: int
) -> dict[Partition, Coefficient]:
    """Undo `clear_denominators`: divide `pairs` by `scale`, exactly.

    Coefficients that come out whole come out as `int`, so an element that
    happened to cross the boundary scaled is indistinguishable from one that
    did not.

        >>> restore([((2,), 3), ((1, 1), 2)], 6)
        {(2,): Fraction(1, 2), (1, 1): Fraction(1, 3)}
        >>> restore([((2,), 6)], 3)
        {(2,): 2}

    A coefficient already rational is divided as one, so this also undoes the
    scaling of a call that returned `(numerator, denominator)` pairs.
    """
    if scale == 1:
        return {la: exact(c) for la, c in pairs}
    return {la: exact(Fraction(c, scale)) for la, c in pairs}


def exact(value: object) -> Coefficient:
    """Normalize a coefficient: a whole `Fraction` becomes an `int`.

        >>> from fractions import Fraction
        >>> exact(Fraction(4, 2)), exact(Fraction(1, 2)), exact(3)
        (2, Fraction(1, 2), 3)

    # Raises

    Raises `TypeError` unless `value` is an `int` or a `Fraction`. Floats and
    other inexact types are rejected here so that every coefficient this
    library returns is exact.
    """
    if isinstance(value, Fraction):
        return int(value) if value.denominator == 1 else value
    if isinstance(value, int):
        return value
    raise TypeError(
        f"coefficient must be an int or a Fraction, not {type(value).__name__}"
    )
