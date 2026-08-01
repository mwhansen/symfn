"""Install symfn in place of Symmetrica everywhere Sage calls it.

    sage -python scripts/check_backend.py      # A/B the two backends

Sage reaches `sage.libs.symmetrica` from exactly six files
(docs/symmetrica-coverage-audit.md). One of them is the conversion table, which
is the large surface and is documented at length below; the other five call
Symmetrica directly and are handled in "the five direct call sites" further
down. `install()` displaces all six.

Two patching mechanisms are needed, because the six sites bind differently.
Five of them — `sf/classical.py`, `sf/sfa.py`, `sf/monomial.py`,
`combinat/tableau.py` — reach a function through the `sage.libs.symmetrica.all`
*module object* at call time, whether by `getattr`, by attribute access, or by
`lazy_import` of the module; rebinding attributes on that module reaches all of
them. `sf/hall_littlewood.py` is the exception: it does
`from ...symmetrica import hall_littlewood_symmetrica as hall_littlewood` at
import, so the name has to be rebound in *that* module's namespace instead.

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

import os
import sys
from functools import reduce
from math import lcm

sys.path.insert(0, "pybuild")

import symfn  # noqa: E402
from sage.all import QQ, ZZ  # noqa: E402

try:
    # The per-term loop, compiled. Optional: the pure-Python fallback below is
    # the same computation, and `scripts/setup_cy.py` builds this when wanted.
    import symfn_cy
except ImportError:  # pragma: no cover
    symfn_cy = None
from sage.combinat.partition import _Partitions  # noqa: E402
from sage.combinat.sf import classical  # noqa: E402
from sage.combinat.sf.sf import SymmetricFunctions  # noqa: E402
from sage.combinat.tableau import Tableau  # noqa: E402
from sage.libs.symmetrica import all as symmetrica_all  # noqa: E402
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing  # noqa: E402

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


# Sage `Partition` objects, keyed by their part tuple.
#
# Rebuilding them is the single largest cost of this shim: `_Partitions(list)`
# validates and interns on every call, and a conversion produces one per output
# term. Measured on 40 shapes of degree 14, the symfn call itself is 0.033s
# while the shim around it is 0.389s -- **91% Python glue**, of which the
# Partition rebuild alone is 9.3x the computation.
#
# They repeat almost perfectly: every conversion of degree n draws from the same
# p(n) partitions, so a plain dict turns the rebuild into a lookup and takes
# that step from 0.279s to 0.018s. `element_class` (skipping validation) only
# reaches 0.154s, so the win is avoiding construction, not avoiding checks.
_PART_CACHE = {}


def _part(key):
    k = tuple(key)
    p = _PART_CACHE.get(k)
    if p is None:
        p = _PART_CACHE[k] = _Partitions(list(k))
    return p


# Sage `Partition` objects for a whole degree, in symfn's own order.
#
# The natural key for a cached partition is its tuple of parts, but building
# that tuple from the list symfn returns, then hashing it, was ~40% of
# everything outside symfn itself. `convert_indexed` sidesteps the question:
# it returns the *position* of each output partition in `symfn.partitions(n)`,
# so the lookup is a list index rather than a hash. No key is faster than no
# key.
_PARTS_BY_DEGREE = {}

# Set SYMFN_NO_PARTITION_CACHE=1 to rebuild the table on every call. This is not
# a tuning knob -- it is how the cache's contribution is measured, since
# Symmetrica's wrapper has no equivalent and a comparison that quietly assumes
# one is not a fair one.
_NO_CACHE = bool(os.environ.get("SYMFN_NO_PARTITION_CACHE"))


def _parts(n):
    if _NO_CACHE:
        return [_Partitions(list(p)) for p in symfn.partitions(n)]
    t = _PARTS_BY_DEGREE.get(n)
    if t is None:
        t = _PARTS_BY_DEGREE[n] = [_Partitions(list(p)) for p in symfn.partitions(n)]
    return t


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

    # Fast path: an integral input converting to an integral basis, which is
    # nearly every call. Everything below stays in Python ints and ZZ, and the
    # output is walked **once**.
    #
    # The general path walks it four times -- build with QQ(c)/den, drop zeros,
    # scan the denominators to choose ZZ or QQ, then build the dict -- and
    # profiling put that, plus the QQ round-trip it implies, above the symfn
    # call itself. Denominators only ever arise from a rational input or from
    # s -> p, so the common case should not pay for them.
    if den == 1 and dst != "powersum":
        raw = symfn.convert_indexed(terms, "Schur", dst)
        # Degrees come from the *input*: a basis change preserves degree, and
        # the input has a handful of terms where the output has thousands.
        # Deriving them from `raw` instead put a full Python pass back over the
        # output and cancelled the compiled loop exactly.
        by_degree = {sum(k): None for k, _ in terms}
        for n in by_degree:
            by_degree[n] = _parts(n)
        if symfn_cy is not None:
            d = symfn_cy.build_terms(raw, by_degree)
        else:
            d = {}
            for n, i, c in raw:
                if c:
                    d[_parts(n)[i]] = ZZ(c)
        return _basis(ZZ, dst)._from_dict(d)

    if dst == "powersum":
        out = [(k, QQ(n) / QQ(dd) / den) for k, (n, dd) in symfn.schur_to_power(terms)]
    elif dst == "Schur":
        out = [(k, QQ(c) / den) for k, c in terms]
    else:
        out = [(k, QQ(c) / den) for k, c in _FROM_SCHUR[dst](terms)]

    out = [(k, v) for k, v in out if v]
    ring = ZZ if all(QQ(v).denominator() == 1 for _, v in out) else QQ
    target = _basis(ring, dst)
    return target._from_dict({_part(k): ring(v) for k, v in out})


def _make(src, dst):
    def entry(d):
        return _convert(d, src, dst)

    entry.__name__ = f"t_{src}_{dst}_symfn"
    return entry


# --- the five direct call sites --------------------------------------------
#
# Everything below replaces a function Sage calls on `sage.libs.symmetrica.all`
# by name, rather than a table entry. Each mirrors its Symmetrica counterpart's
# return *type*, not just its value: Sage consumes these without coercion --
# `sfa._expand` feeds the result to `resPR(...)`, `monomial.py` reads
# `.monomial_coefficients()`, `tableau.py` wraps each element in
# `self.element_class`, and `hall_littlewood.py` calls `.coefficient(part).subs`
# -- so a bare dict or a list of lists would fail at the caller, not here.


def _with_alphabet(basis):
    """`compute_{basis}_with_alphabet`: one basis element as a polynomial.

    This is the site the coverage audit mis-read as a binding gap. `eval()`
    evaluates at an alphabet of *ring elements* and cannot serve it: what Sage
    wants back is a polynomial in `n` indeterminates, which is why symfn grew
    `expand_alphabet` rather than a wrapper around `evaluate_schur`.
    """

    def entry(part, n, alphabet="x"):
        n = int(n)
        ring = PolynomialRing(ZZ, n, alphabet)
        terms = symfn.expand_alphabet([(list(part), 1)], basis, n)
        # symfn guarantees distinct exponent vectors and no zero coefficients,
        # so this dict needs no merging pass.
        return ring({tuple(a): ZZ(c) for a, c in terms})

    entry.__name__ = f"compute_{basis}_with_alphabet_symfn"
    return entry


def _mult_monomial_monomial(left, right):
    """`mult_monomial_monomial`: the product of two monomial-basis elements.

    Sage only ever calls this with two single-term, coefficient-1 dicts (see
    `sf/monomial.py:129`, which does its own outer double loop and its own
    empty-partition special case -- the latter because Symmetrica *aborts the
    process* on two empty partitions). Denominators are cleared anyway, on the
    same reasoning as `_convert`: the caller's contract is a dict of
    coefficients, and nothing in it promises they are integral.
    """
    items = [(list(k), QQ(v)) for k, v in _items(left)], [
        (list(k), QQ(v)) for k, v in _items(right)
    ]
    den = 1
    for side in items:
        for _, v in side:
            den = lcm(den, int(v.denominator()))
    a, b = ([(k, int(v * den)) for k, v in side] for side in items)
    out = [(k, QQ(c) / QQ(den * den)) for k, c in symfn.monomial_multiply(a, b)]
    out = [(k, v) for k, v in out if v]
    ring = ZZ if all(v.denominator() == 1 for _, v in out) else QQ
    return _basis(ring, "monomial")._from_dict({_part(k): ring(v) for k, v in out})


def _kostka_number(shape, weight):
    """`kostka_number`: K_{shape,weight}, an `Integer` as Symmetrica returns."""
    return ZZ(symfn.kostka_number(list(shape), list(weight)))


def _kostka_tab(shape, weight):
    """`kostka_tab`: the SSYT themselves, as Sage `Tableau` objects.

    The order is Symmetrica's -- increasing lexicographic in the row-major
    reading word -- and it is part of the interface, not a formatting choice:
    `SemistandardTableaux(shape, weight).list()` returns this list verbatim and
    Sage's doctests print it. symfn sorts to that order; see
    `semistandard_tableaux` in `src/kostka.rs`.
    """
    return [Tableau(t) for t in symfn.semistandard_tableaux(list(shape), list(weight))]


def _hall_littlewood(part):
    """`hall_littlewood`: `Q'_part` in the Schur basis, over `ZZ[x]`.

    The parameter is called `x` here and not `t` because that is the ring
    Symmetrica hands back and `sf/hall_littlewood.py:989` substitutes for it by
    name (`res.coefficient(part2).subs(x=t)`). Returning `ZZ[t]` would raise
    there rather than differ quietly.
    """
    ring = ZZ["x"]
    rows = symfn.hall_littlewood(list(part))
    return _basis(ring, "Schur")._from_dict(
        {_part(mu): ring({int(e): ZZ(c) for e, c in poly}) for mu, poly in rows}
    )


# The `sage.libs.symmetrica.all` attributes to rebind, and what to rebind them
# to. `sfa._expand` builds these names with `getattr`, so a typo here is a
# silent no-op that leaves Symmetrica installed -- `check_backend.py` is what
# catches that, by comparing against Symmetrica rather than against nothing.
def _direct_entries():
    entries = {
        "mult_monomial_monomial": _mult_monomial_monomial,
        "kostka_number": _kostka_number,
        "kostka_tab": _kostka_tab,
    }
    for sage_name, ours in (
        ("schur", "Schur"),
        ("monomial", "monomial"),
        ("homsym", "homogeneous"),
        ("elmsym", "elementary"),
        ("powsym", "powersum"),
    ):
        entries[f"compute_{sage_name}_with_alphabet"] = _with_alphabet(ours)
    return entries


_ORIGINAL = {}
_ORIGINAL_DIRECT = {}


def install():
    """Point every Sage call site at symfn. Returns the previous table."""
    if not classical.conversion_functions:
        classical.init()
    if not _ORIGINAL:
        _ORIGINAL.update(classical.conversion_functions)
    for src in NAMES:
        for dst in NAMES:
            if src != dst:
                classical.conversion_functions[(src, dst)] = _make(src, dst)

    # Importing here rather than at module scope so that merely importing this
    # module stays side-effect-free on Sage's symmetric-function machinery.
    from sage.combinat.sf import hall_littlewood as hl_module

    if not _ORIGINAL_DIRECT:
        for name in _direct_entries():
            _ORIGINAL_DIRECT[name] = getattr(symmetrica_all, name)
        _ORIGINAL_DIRECT["hall_littlewood"] = hl_module.hall_littlewood
    for name, fn in _direct_entries().items():
        setattr(symmetrica_all, name, fn)
    hl_module.hall_littlewood = _hall_littlewood
    return _ORIGINAL


def restore():
    """Put Symmetrica back, at all six sites."""
    classical.conversion_functions.update(_ORIGINAL)
    if not _ORIGINAL_DIRECT:
        return
    from sage.combinat.sf import hall_littlewood as hl_module

    for name, fn in _ORIGINAL_DIRECT.items():
        if name == "hall_littlewood":
            hl_module.hall_littlewood = fn
        else:
            setattr(symmetrica_all, name, fn)


def original():
    """The Symmetrica table, for A/B comparison."""
    if not _ORIGINAL:
        install()
        restore()
    return dict(_ORIGINAL)
