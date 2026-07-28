"""Install symfn as Sage's classical-basis conversion backend, replacing Symmetrica.

    sage -python scripts/check_backend.py      # A/B the two backends

Sage routes every conversion between the five classical bases through
`sage.combinat.sf.classical.conversion_functions`, a table of 20 entries keyed
by (from_basis, to_basis) and filled at import with Symmetrica's C functions.
That table is the entire integration surface: replacing its values swaps the
backend for all of Sage's symmetric-function machinery at once, including the
paths that only reach a conversion indirectly (products in a non-Schur basis,
`scalar`, `expand`, plethysm, Hall-Littlewood and Macdonald bases, ...).

The contract each entry must honour, read off Symmetrica's own behaviour:

* input is a **nonempty** dict {Partition: coefficient}; Sage guards the empty
  case itself, and Symmetrica aborts the process rather than raising on it;
* coefficients may be rational, and the empty partition is a legal key;
* the return value is anything with `._monomial_coefficients`, a dict
  {Partition: coefficient}. Symmetrica hands back an element over **ZZ** when
  every coefficient is integral and over QQ otherwise, and this matches that:
  in the non-QQ path Sage calls `_from_dict` *without* coercion, so the value
  types are not merely cosmetic.

Denominators are cleared before the call and restored after. Every conversion
here is linear over the coefficient ring, so scaling by the common denominator
D and dividing the result by D is exact -- and it is what lets an integer-only
boundary serve a rational input.

symfn is invoked as a library from this process; Sage is not modified on disk.
"""

import sys
from functools import reduce
from math import lcm

sys.path.insert(0, "pybuild")

import symfn  # noqa: E402
from sage.all import QQ, ZZ  # noqa: E402
from sage.combinat.partition import _Partitions  # noqa: E402
from sage.combinat.sf import classical  # noqa: E402
from sage.combinat.sf.sf import SymmetricFunctions  # noqa: E402

# Sage's basis names, as they appear as keys in `conversion_functions`.
NAMES = ["Schur", "monomial", "homogeneous", "elementary", "powersum"]

_TO_SCHUR = {
    "monomial": symfn.monomial_to_schur,
    "homogeneous": symfn.homogeneous_to_schur,
    "elementary": symfn.elementary_to_schur,
    "powersum": symfn.power_to_schur,
}
_FROM_SCHUR = {
    "monomial": symfn.schur_to_monomial,
    "homogeneous": symfn.schur_to_homogeneous,
    "elementary": symfn.schur_to_elementary,
}


def _basis(ring, name):
    sym = SymmetricFunctions(ring)
    return {
        "Schur": sym.schur,
        "monomial": sym.monomial,
        "homogeneous": sym.homogeneous,
        "elementary": sym.elementary,
        "powersum": sym.power,
    }[name]()


def _items(x):
    """The input, whichever of Sage's two calling conventions produced it.

    `sage.combinat.sf.classical` passes a plain `{Partition: coeff}` dict, but
    `sage.combinat.sf.sf.SymmetricaConversionOnBasis` — the wrapper that builds
    conversion *morphisms*, and therefore the one every non-QQ base ring goes
    through — passes a `CombinatorialFreeModule` element instead. Its docstring
    describes the argument as "a monomial in CombinatorialFreeModule(QQ,
    Partitions())", and calls `dict(...)` on whatever comes back.

    Only the first convention is visible from the conversion table, so this was
    invisible until Hall-Littlewood was exercised. A backend that handles only
    dicts appears to work and then fails the moment a caller reaches for
    Macdonald, Jack, HL, or any base ring other than QQ.
    """
    mc = getattr(x, "monomial_coefficients", None)
    return mc().items() if mc is not None else x.items()


def _convert(d, src, dst):
    """One entry of the table: `src` basis -> element in `dst`.

    The return value is a Sage element, which satisfies both consumers: the
    table path reads `._monomial_coefficients`, and the morphism path calls
    `dict(...)`, for which iterating as (partition, coefficient) pairs is
    exactly right.
    """
    items = list(_items(d))
    # Clear denominators so the integer boundary can carry the input.
    den = reduce(lcm, (int(QQ(v).denominator()) for _, v in items), 1)
    terms = [(list(k), int(QQ(v) * den)) for k, v in items]

    if src != "Schur":
        terms = _TO_SCHUR[src](terms)

    if dst == "powersum":
        out = [(k, QQ(n) / QQ(dd) / den) for k, (n, dd) in symfn.schur_to_power(terms)]
    elif dst == "Schur":
        out = [(k, QQ(c) / den) for k, c in terms]
    else:
        out = [(k, QQ(c) / den) for k, c in _FROM_SCHUR[dst](terms)]

    out = [(k, v) for k, v in out if v]
    ring = ZZ if all(QQ(v).denominator() == 1 for _, v in out) else QQ
    target = _basis(ring, dst)
    return target._from_dict({_Partitions(list(k)): ring(v) for k, v in out})


def _make(src, dst):
    def entry(d):
        return _convert(d, src, dst)

    entry.__name__ = f"t_{src}_{dst}_symfn"
    return entry


_ORIGINAL = {}


def install():
    """Point Sage's conversion table at symfn. Returns the previous table."""
    if not classical.conversion_functions:
        classical.init()
    if not _ORIGINAL:
        _ORIGINAL.update(classical.conversion_functions)
    for src in NAMES:
        for dst in NAMES:
            if src != dst:
                classical.conversion_functions[(src, dst)] = _make(src, dst)
    return _ORIGINAL


def restore():
    """Put Symmetrica back."""
    classical.conversion_functions.update(_ORIGINAL)


def original():
    """The Symmetrica table, for A/B comparison."""
    if not _ORIGINAL:
        install()
        restore()
    return dict(_ORIGINAL)
