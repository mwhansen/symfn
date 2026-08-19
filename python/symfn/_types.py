"""The type vocabulary the convenience layer is annotated in.

The aliases live in one module because the distinction between what the
boundary *promises* and what it *accepts* is part of the contract, and a
structural type restated at every signature loses it. A partition comes back
as `Partition` and goes in as `PartitionArg`, and the asymmetry is deliberate
at both layers.

`Basis` is a `Literal` rather than `str` on purpose. The convenience layer
carries basis identity — the crate's "basis confusion is a compile error"
translated into a language with no compiler — and a checked `Literal` is
as close to that as Python gets: `element.to("Schur")` is now a type error
before it is a `ValueError`. The narrowing from `str` is earned in exactly one
place, `check_basis`, which is the function that validates it.

The 3.9 floor is why `Union` and `Optional` appear in the alias assignments
below rather than `|`: PEP 585 generics (`tuple[int, ...]`) evaluate at runtime
on 3.9, PEP 604 unions do not. Inside annotations either spelling is fine,
because every module here imports `annotations` from `__future__`.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from fractions import Fraction
from typing import Literal, Union

__all__ = [
    "Basis",
    "Coefficient",
    "Partition",
    "PartitionArg",
    "Permutation",
    "PermutationArg",
]

#: A basis of the ring of symmetric functions, by one-letter code. The contract
#: layer accepts the same codes, so a `Basis` passes through unchanged.
Basis = Literal["s", "h", "e", "p", "m", "f"]

#: A partition as it crosses *out*: a zero-free tuple, weakly decreasing.
Partition = tuple[int, ...]

#: A partition as it crosses *in*. A bare `int` is the one-row shape, which is
#: what makes `s[2]` mean `s_(2)` rather than an error; trailing zeros are
#: tolerated because that is the fixed-width form callers hand over.
PartitionArg = Union[int, Sequence[int]]

#: A permutation in one-line notation, out and in. Trailing fixed points are
#: dropped on the way out.
Permutation = tuple[int, ...]
PermutationArg = Union[int, Sequence[int]]

#: Every coefficient in this layer, and there are no others: exact integers of
#: any size, and rationals where a conversion divides. Nothing inexact ever
#: enters.
Coefficient = Union[int, Fraction]

#: The two shapes a term collection is accepted in — a mapping, or the pairs
#: themselves. Both constructors take either.
#
# The mapping branches name `Partition`/`Permutation` and `int` rather than the
# permissive `*Arg` aliases, and that is not a tightening: a mapping key has to
# be hashable, so a support that arrives as a key is already a tuple or a bare
# `int`. A list reaches the constructor only through the pairs form, which
# `Iterable` and `tuple` being covariant let it do. Spelling the mapping branch
# permissively instead would be *rejected* by a checker, since `Mapping` is
# invariant in its key.
TermsArg = Union[
    Mapping[Partition, Coefficient],
    Mapping[int, Coefficient],
    Iterable[tuple[PartitionArg, Coefficient]],
]
#: A polynomial keyed by exponent vector, which is what `Schub.expand` returns
#: and `from_polynomial` reads back. The key is a tuple for the same reason the
#: term mappings above use one.
PolynomialArg = Union[
    Mapping[tuple[int, ...], int],
    Iterable[tuple[Sequence[int], int]],
]
SchubertTermsArg = Union[
    Mapping[Permutation, int],
    Mapping[int, int],
    Iterable[tuple[PermutationArg, int]],
]
