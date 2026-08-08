"""Exact symmetric functions, for any CPython 3.9 or later.

The package has two layers, and both are supported.

The **contract layer** is the compiled module: 108 whole-object entry points
over plain data — lists of `(partition, coefficient)` pairs, `int`
coefficients of any size, nothing to import in order to unpack a result. Its
names are re-exported here, so `symfn.schur_multiply` and
`symfn.symfn.schur_multiply` are the same function, and `help(symfn.symfn)`
carries the data model in full.

The **convenience layer** is this package's pure Python: `Sym`, the basis
factories `s`, `h`, `e`, `p`, `m`, `f`, and the parameter families under
`macdonald`, `jack`, `hl` and `llt`. It computes nothing of its
own — every method is a contract call with the basis bookkeeping done for you
(``docs/policies/python.md``, P4).

    >>> import symfn
    >>> from symfn import s, h
    >>> s([2, 1]) * s([1])
    s[2,1,1] + s[2,2] + s[3,1]
    >>> h([2]).to("s")
    s[2]
    >>> symfn.schur_multiply([([2, 1], 1)], [([1], 1)])
    [((2, 1, 1), 1), ((2, 2), 1), ((3, 1), 1)]

Reach for the contract layer when you are marshalling in bulk or building
another library on top; reach for the convenience layer when a person is
reading the code.

Sage is not imported here, depended on, or required. The adapter that lets
Sage use this library lives on the Sage side of the boundary.
"""

from . import symfn
from ._bases import BasisError
from ._families import hl, jack, llt, macdonald
from ._param import AlphaFrac, Param, Poly, QtFrac, QtPoly
from ._schubert import Schub, X
from ._sym import Sym, e, f, h, m, p, s, skew
from .symfn import *  # noqa: F403  the contract layer, re-exported flat

__version__ = symfn.__version__

#: Every name this package supports, contract layer first.
__all__ = (
    [n for n in dir(symfn) if not n.startswith("_")]
    + [
        "Sym",
        "s",
        "h",
        "e",
        "p",
        "m",
        "f",
        "skew",
        "BasisError",
        "macdonald",
        "jack",
        "hl",
        "llt",
        "Param",
        "Poly",
        "QtPoly",
        "QtFrac",
        "AlphaFrac",
        "Schub",
        "X",
    ]
    + ["__version__"]
)
