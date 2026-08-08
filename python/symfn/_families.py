"""The parameter families, as namespaces over the contract layer.

Four objects — `macdonald`, `jack`, `hall_littlewood`, `llt` — group the entry
points that share a parameter and a set of conventions, and wrap their rows in
the types of `_param.py`. Nothing here computes: each method is one contract
call and one constructor (``docs/policies/python.md``, P4).

Every family in this file has a normalization that plausible rivals disagree
with, so each method's doc names the convention and each carries an example
whose value rules the rivals out. The Rust module docs are the depth —
``src/macdonald.rs``, ``src/jack.rs``, ``src/hl.rs``, ``src/llt.rs`` — and this
layer states the convention rather than delegating it.
"""

from . import symfn as _c
from ._param import AlphaFrac, Param, Poly, QtFrac, QtPoly
from ._sym import _partition

__all__ = ["macdonald", "jack", "hl", "llt"]


def _qt_element(rows, basis):
    """Wrap `(partition, [(a, b, coefficient)])` rows as a `Param` in q, t."""
    return Param(basis, [(la, QtPoly(c)) for la, c in rows], ("q", "t"))


def _q_element(rows, basis):
    """Wrap `(partition, [(a, b, coefficient)])` rows as a `Param` in q alone.

    The LLT entry points share the `(q_exponent, t_exponent, coefficient)`
    encoding with the Macdonald operators and use only the first, because LLT
    is a one-parameter family. Reading the rows as `(q, t)` would make
    `parameters` say so and force every caller to pass a `t` that appears
    nowhere.

    # Raises

    Raises `ValueError` if a row carries a nonzero `t` exponent, which would
    mean the family had grown a second parameter and this projection was
    dropping it.
    """
    terms = []
    for la, cells in rows:
        poly = {}
        for a, b, c in cells:
            if b:
                raise ValueError(
                    f"LLT row for {tuple(la)} carries t^{b}; this family is in "
                    "q alone"
                )
            poly[a] = c
        terms.append((la, Poly("q", poly)))
    return Param(basis, terms, ("q",))


def _mac_element(rows):
    """Wrap `(partition, numerator, denominator)` rows as a `Param` in q, t."""
    return Param("m", [(la, QtFrac(n, d)) for la, n, d in rows], ("q", "t"))


def _jack_element(rows, basis="m"):
    """Wrap `(partition, numerator, atoms, scale)` rows as a `Param` in α."""
    return Param(
        basis, [(la, AlphaFrac(n, d, k)) for la, n, d, k in rows], ("alpha",)
    )


class _Macdonald:
    """The Macdonald family in `q` and `t`, and the operators built on it.

    `P` is monic in the monomial basis, `Q = b_λ·P`, and `J = c_λ·P` is the
    integral form whose coefficients are polynomials. `Htilde` is the modified
    form `H̃_μ`, the one the `(q,t)`-Kostka polynomials expand.

        >>> from symfn import macdonald
        >>> macdonald.P([1])
        m[1]
        >>> macdonald.Q([1])
        (1 - t)/(1 - q)*m[1]

    Sage's equivalents are `Sym.macdonald().P()`, `.Q()`, `.J()` and `.Ht()`.
    """

    __module__ = "symfn"

    def P(self, la):
        """`P_λ(x; q, t)` in the monomial basis, monic in `m_λ`.

            >>> from symfn import macdonald
            >>> macdonald.P([2])
            (1 - q)*(1 - t^2)/(1 - t)*(1 - q*t)*m[1,1] + m[2]
            >>> macdonald.P([2]).at(q=5, t=5)
            m[1,1] + m[2]

        The leading coefficient 1 is the normalization; at `q = t` the whole
        family collapses to the Schur function, which is the check that this is
        `P` and not `Q`.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _mac_element(_c.macdonald_p(_partition(la)))

    def Q(self, la):
        """`Q_λ = b_λ · P_λ`, in the monomial basis.

            >>> from symfn import macdonald
            >>> macdonald.Q([1]).coefficient([1])
            (1 - t)/(1 - q)

        `Q_(1) = (1 − t)/(1 − q)·m_1` where `P_(1) = m_1`: the value that
        separates the two normalizations at the smallest shape.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _mac_element(_c.macdonald_q(_partition(la)))

    def J(self, la):
        """`J_λ = c_λ · P_λ`, the integral form, in the monomial basis.

            >>> from symfn import macdonald
            >>> macdonald.J([1, 1]).coefficient([1, 1])
            1 - t - t^2 + t^3

        Every coefficient is a polynomial — the empty denominator is the
        integral form's signature, since `J` clears what `Q` carries.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _mac_element(_c.macdonald_j(_partition(la)))

    def Htilde(self, mu):
        """The modified Macdonald polynomial `H̃_μ`, in the Schur basis.

            >>> from symfn import macdonald
            >>> macdonald.Htilde([2])
            q*s[1,1] + s[2]

        The Schur coefficients are the `(q,t)`-Kostka polynomials `K̃_{λμ}`,
        and `H̃_(2)` carrying `q` rather than `t` on `s_11` is the orientation
        of the Garsia-Haiman convention.

        # Raises

        Raises `ValueError` unless μ is a partition.
        """
        return _qt_element(_c.macdonald_ht(_partition(mu)), "s")

    def qt_kostka(self, la, mu):
        """The `(q,t)`-Kostka polynomial `K̃_{λμ}(q, t)`, as a `QtPoly`.

            >>> from symfn import macdonald
            >>> macdonald.qt_kostka([2], [1, 1])
            q
            >>> macdonald.qt_kostka([1, 1], [2]).at(q=1, t=1)
            1

        # Raises

        Raises `ValueError` unless both arguments are partitions of the same
        integer.
        """
        return QtPoly(_c.qt_kostka(_partition(la), _partition(mu)))

    def nabla_e(self, n):
        """`∇ e_n`, in the Schur basis — the Shuffle Theorem's left side.

            >>> from symfn import macdonald
            >>> macdonald.nabla_e(2)
            (q + t)*s[1,1] + s[2]

        # Raises

        Raises `ValueError` unless `n` is positive.
        """
        return _qt_element(_c.nabla_e(n), "s")

    def nabla(self, f):
        """`∇` applied to a `Param` in the Schur basis, or to contract rows.

            >>> from symfn import macdonald, s
            >>> macdonald.nabla(s([1, 1]) + s([2]))
            (q + t)*s[1,1] + s[2]

        # Raises

        Raises `ValueError` unless the argument is homogeneous in the Schur
        basis, which `∇` requires: it acts by a scalar on each `H̃_μ`, and a
        sum across degrees has no single one.
        """
        rows = _schur_rows(f)
        return _qt_element(_c.nabla(rows), "s")

    def __repr__(self):
        return "symfn.macdonald"


class _Jack:
    """The Jack family in α: `P` monic in `m_λ`, `Q` and `J` its rescalings.

    At `α = 1` every one of them degenerates to a Schur function, and at
    `α = 2` to a zonal polynomial; both are checks a caller can run.

        >>> from symfn import jack
        >>> jack.P([2]).at(alpha=1)
        m[1,1] + m[2]

    Sage's equivalents are `Sym.jack().P()`, `.Q()` and `.J()`.
    """

    __module__ = "symfn"

    def P(self, la):
        """`P_λ(x; α)` in the monomial basis, monic in `m_λ`.

            >>> from symfn import jack
            >>> jack.P([2])
            2/(alpha + 1)*m[1,1] + m[2]

        The leading 1 is what separates `P` from `Q` and `J`, both of which
        scale it.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _jack_element(_c.jack_p(_partition(la)))

    def Q(self, la):
        """`Q_λ(x; α)` in the monomial basis.

            >>> from symfn import jack
            >>> jack.Q([2]).coefficient([2])
            (1 + alpha)/2*alpha

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _jack_element(_c.jack_q(_partition(la)))

    def J(self, la):
        """`J_λ(x; α)`, the integral form, in the monomial basis.

            >>> from symfn import jack
            >>> jack.J([1, 1])
            2*m[1,1]

        The α-free coefficients are the integral form's signature.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _jack_element(_c.jack_j(_partition(la)))

    def zonal(self, la, integral_form=False):
        """The zonal polynomial, `α = 2`, in the monomial basis, as a `Sym`.

            >>> from symfn import jack
            >>> jack.zonal([2])
            2/3*m[1,1] + m[2]

        The parameter is already substituted, so this returns a `Sym` rather
        than a `Param`; `jack.P(la).at(alpha=2)` is the same element.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        from fractions import Fraction

        from ._sym import Sym

        rows = _c.zonal(_partition(la), bool(integral_form))
        return Sym("m", {la_: Fraction(n, d) for la_, n, d in rows})

    def __repr__(self):
        return "symfn.jack"


class _HallLittlewood:
    """The Hall-Littlewood family in `t`, and the Kostka-Foulkes polynomials.

    `Qp` is `Q'_λ`, the one whose Schur coefficients are the Kostka-Foulkes
    polynomials `K_{μλ}(t)`; `P` is `P_λ`, monic in the monomial basis.

        >>> from symfn import hl
        >>> hl.Qp([2, 1]).at(t=0)
        s[2,1]
        >>> hl.Qp([2, 1]).at(t=1)
        s[1,1,1] + 2*s[2,1] + s[3]

    `Q'_λ(x; 0) = s_λ` and `Q'_λ(x; 1) = h_λ` expanded in Schur functions,
    whose coefficients are the Kostka numbers. Both are theorems, and both
    fail under the `t → 1/t` twist that the other convention in circulation
    uses. Sage's equivalent is `Sym.hall_littlewood().Qp()`.
    """

    __module__ = "symfn"

    def Qp(self, la):
        """`Q'_λ(x; t) = Σ_μ K_{μλ}(t) s_μ`, in the Schur basis.

            >>> from symfn import hl
            >>> hl.Qp([1, 1])
            t*s[1,1] + s[2]

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _t_element(_c.hall_littlewood(_partition(la)), "s")

    def P(self, la):
        """`P_λ(x; t)` in the Schur basis.

            >>> from symfn import hl
            >>> hl.P([2])
            -t*s[1,1] + s[2]

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _t_element(_c.hall_littlewood_p(_partition(la)), "s")

    def kostka_foulkes(self, la, mu):
        """The Kostka-Foulkes polynomial `K_{λμ}(t)`, as a `Poly`.

            >>> from symfn import hl
            >>> hl.kostka_foulkes([2, 2], [1, 1, 1, 1])
            t^2 + t^3 + t^4

        `K_{λμ}(1)` is the Kostka number, and the absence of a constant term
        distinguishes the charge convention from its cocharge rival.

        # Raises

        Raises `ValueError` unless both arguments are partitions of the same
        integer.
        """
        return Poly("t", _c.kostka_foulkes(_partition(la), _partition(mu)))

    def __repr__(self):
        return "symfn.hl"


class _LLT:
    """The LLT family in `q`, over a tuple of shapes.

    The conventions here differ by more than a twist — see the module doc of
    ``src/llt.rs``, "The conventions in circulation" — so each method names
    which `G` or `H` it computes.

        >>> from symfn import llt
        >>> llt.H([2], 2)
        (q + 1)*s[2]
    """

    __module__ = "symfn"

    def Gtilde(self, la, k):
        """`G̃` for the `k`-quotient of λ, in the Schur basis.

            >>> from symfn import llt
            >>> llt.Gtilde([2, 1], 2)
            q*s[1,1,1] + s[2,1]

        # Raises

        Raises `ValueError` unless λ is a partition and `k` is positive.
        """
        return _q_element(_c.llt_gtilde(_partition(la), k), "m")

    def H(self, mu, k):
        """The LLT `H` for μ at level `k`, in the Schur basis.

            >>> from symfn import llt
            >>> llt.H([1, 1], 2)
            (q + 1)*s[1,1]

        # Raises

        Raises `ValueError` unless μ is a partition and `k` is positive.
        """
        return _q_element(_c.llt_h(_partition(mu), k), "m")

    def G(self, shapes, offsets=None):
        """The LLT product `G` over a tuple of shapes, in the Schur basis.

            >>> from symfn import llt
            >>> llt.G([[1], [1]])
            q*s[1,1] + s[2]

        `q` on `s_11` rather than on `s_2` is the inversion statistic's
        orientation, and what separates this from the `q → 1/q` convention.

        # Raises

        Raises `ValueError` unless every shape is a partition, and unless
        `offsets` — when given — has one entry per shape.
        """
        rows = _c.llt_g([_partition(sh) for sh in shapes], offsets)
        return _q_element(rows, "m")

    def __repr__(self):
        return "symfn.llt"


def _t_element(rows, basis):
    """Wrap `(partition, [(t_exponent, coefficient)])` rows as a `Param`."""
    return Param(basis, [(la, Poly("t", c)) for la, c in rows], ("t",))


def _schur_rows(f):
    """The `(partition, [(a, b, coefficient)])` rows `nabla` takes.

    Accepts a `Sym` in the Schur basis, a `Param` in `q` and `t`, or the rows
    themselves.
    """
    from ._sym import Sym

    if isinstance(f, Sym):
        if f.basis != "s":
            raise ValueError(f"nabla needs a Schur-basis element, not {f.basis}")
        return [(la, [(0, 0, int(c))]) for la, c in f]
    if isinstance(f, Param):
        if f.basis != "s":
            raise ValueError(f"nabla needs a Schur-basis element, not {f.basis}")
        return [
            (la, [(a, b, c) for (a, b), c in coeff.coefficients().items()])
            for la, coeff in f
        ]
    return list(f)


#: The Macdonald family.
macdonald = _Macdonald()
#: The Jack family.
jack = _Jack()
#: The Hall-Littlewood family. Named `hl` because `hall_littlewood`
#: is a contract entry point and the flat surface keeps its names (P10).
hl = _HallLittlewood()
#: The LLT family.
llt = _LLT()
