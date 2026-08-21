"""The parameter families: `macdonald`, `jack`, `hl` and `llt`.

Each is a namespace object whose methods call one contract entry point and
wrap the returned rows in the coefficient types of `_param.py`, so a method
here returns the same values as the flat `symfn.*` function it names, in a
form with a `repr` and an `at` for substituting parameter values.

Every family has rival normalizations in the literature that differ by a
twist (`q ↔ t`, `t → 1/t`, `α → 1/α`), and a wrong one returns a plausible
answer rather than an error. So each method's docstring states its
convention and gives an example whose value distinguishes it from the
rivals, and the Conventions page of the rendered documentation collects
these per family.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from fractions import Fraction
from math import gcd
from typing import Any, Union

from . import symfn as _c
from ._param import AlphaFrac, Param, Poly, QtFrac, QtPoly, QtRatio
from ._sym import Sym, _partition
from ._types import Coefficient, Partition, PartitionArg

#: What `nabla` accepts: an element in the Schur basis in either of the two
#: types that can carry one, or the contract layer's rows themselves.
NablaArg = Union["Sym", Param, Iterable[Any]]

#: What the contract layer's `(q, t)`-graded rows look like on arrival.
QtRows = Iterable[tuple[Partition, Iterable[tuple[int, int, Coefficient]]]]

#: What the Hall-Littlewood inverse expansions accept: a Schur-basis element
#: as a `Sym`, as a `Param` in `t`, or as the contract layer's `t`-rows.
TSchurArg = Union["Sym", Param, Iterable[Any]]

__all__ = ["macdonald", "jack", "hl", "llt"]


def _qt_element(rows: QtRows, basis: str) -> Param:
    """Wrap `(partition, [(a, b, coefficient)])` rows as a `Param` in q, t."""
    return Param(basis, [(la, QtPoly(c)) for la, c in rows], ("q", "t"))


def _q_element(rows: QtRows, basis: str) -> Param:
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


def _ht_element(rows: Iterable[Any]) -> Param:
    """Wrap `(partition, numerator, denominator)` rows as a `Param` in the
    modified Macdonald basis, tagged as Sage prints it.
    """
    return Param("McdHt", [(la, QtRatio(n, d)) for la, n, d in rows], ("q", "t"))


def _mac_element(rows: Iterable[Any]) -> Param:
    """Wrap `(partition, numerator, denominator)` rows as a `Param` in q, t."""
    return Param("m", [(la, QtFrac(n, d)) for la, n, d in rows], ("q", "t"))


def _jack_element(rows: Iterable[Any], basis: str = "m") -> Param:
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

    def P(self, la: PartitionArg) -> Param:
        """`P_λ(x; q, t)` in the monomial basis, monic in `m_λ`.

            >>> from symfn import macdonald
            >>> macdonald.P([2])
            (1 - t + q - q*t)/(1 - q*t)*m[1,1] + m[2]
            >>> macdonald.P([2]).at(q=5, t=5)
            m[1,1] + m[2]

        The `m_11` coefficient is `(1 + q)(1 − t)/(1 − q·t)` with its
        numerator expanded. The leading coefficient 1 is the normalization; at
        `q = t` the whole family collapses to the Schur function, which is the
        check that this is `P` and not `Q`.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _mac_element(_c.macdonald_p(_partition(la)))

    def Q(self, la: PartitionArg) -> Param:
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

    def J(self, la: PartitionArg) -> Param:
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

    def Htilde(self, mu: PartitionArg) -> Param:
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

    def to_Htilde(self, f: NablaArg) -> Param:
        """`f`, given in the Schur basis, rewritten in the `H̃` basis — the
        direction `Htilde` does not go.

            >>> from symfn import macdonald, s
            >>> macdonald.to_Htilde(s([2]))
            q/(-t + q)*McdHt[1,1] - t/(-t + q)*McdHt[2]
            >>> macdonald.to_Htilde(macdonald.Htilde([2, 1]))
            McdHt[2,1]

        `s_2 = q/(q−t)·H̃_11 − t/(q−t)·H̃_2`; the `q` upstairs on the column
        shape is the orientation, and the `q ↔ t` swap gives a different
        answer rather than an error. The coefficients are `QtRatio`s and
        genuinely not polynomials: this divides by `w_μ`, whose factors are
        `q^a − t^b` and do not cancel.

        Accepts what `nabla` accepts — a `Sym` or a `Param` in the Schur
        basis, or the contract rows — and unlike `nabla` it takes mixed
        degrees, expanding each degree on its own.

        # Raises

        Raises `ValueError` unless the argument is in the Schur basis with
        integer coefficients, or coefficients in `q` and `t`.
        """
        return _ht_element(_c.schur_to_macdonald_ht(_schur_rows(f, "to_Htilde")))

    def qt_kostka(self, la: PartitionArg, mu: PartitionArg) -> QtPoly:
        """The `(q,t)`-Kostka polynomial `K̃_{λμ}(q, t)`, as a `QtPoly`.

            >>> from symfn import macdonald
            >>> macdonald.qt_kostka([2], [1, 1])
            t
            >>> macdonald.qt_kostka([1, 1], [2]).at(q=1, t=1)
            1

        `K̃_{(2),(11)} = t` rather than `q` is the Garsia-Haiman orientation;
        the `q ↔ t` mirror swaps this value with `Htilde`'s example.

        # Raises

        Raises `ValueError` unless both arguments are partitions of the same
        integer.
        """
        return QtPoly(_c.qt_kostka(_partition(la), _partition(mu)))

    def nabla_e(self, n: int) -> Param:
        """`∇ e_n`, in the Schur basis — the Shuffle Theorem's left side.

            >>> from symfn import macdonald
            >>> macdonald.nabla_e(2)
            (t + q)*s[1,1] + s[2]

        # Raises

        Raises `ValueError` unless `n` is positive.
        """
        return _qt_element(_c.nabla_e(n), "s")

    def nabla(self, f: NablaArg) -> Param:
        """`∇` applied to a `Param` in the Schur basis, or to contract rows.

            >>> from symfn import macdonald, s
            >>> macdonald.nabla(s([1, 1]) + s([2]))
            (t + q - q*t)*s[1,1] + s[2]

        The `−q·t` against `nabla_e(2)`'s value is `∇s_2 = −q·t·s_11`.

        # Raises

        Raises `ValueError` unless the argument is homogeneous in the Schur
        basis, which `∇` requires: it acts by a scalar on each `H̃_μ`, and a
        sum across degrees has no single one. Coefficients must be integers
        or polynomials in `q` and `t`.
        """
        rows = _schur_rows(f)
        return _qt_element(_c.nabla(rows), "s")

    def nabla_power(self, f: NablaArg, r: int) -> Param:
        """`∇^r F`, sharing one change of basis across the powers.

            >>> from symfn import macdonald, s
            >>> macdonald.nabla_power(s([1, 1]), 1) == macdonald.nabla_e(2)
            True

        `r = 1` is `nabla`, which is the check that the shared change of basis
        computes the same operator as applying it once.

        # Raises

        Raises `ValueError` unless the argument is homogeneous in the Schur
        basis, and unless `r` is non-negative.
        """
        return _qt_element(_c.nabla_power(_schur_rows(f), r), "s")

    def delta_prime_e(self, k: int, n: int) -> Param:
        """`Δ'_{e_k} e_n` in the Schur basis — the Delta conjecture's object.

            >>> from symfn import macdonald
            >>> macdonald.delta_prime_e(1, 2)
            (t + q)*s[1,1] + s[2]

        At `k = n − 1` this is `∇e_n`, which is the Shuffle Theorem's object
        and the value above at `n = 2`.

        # Raises

        Raises `ValueError` unless `0 < k < n`.
        """
        return _qt_element(_c.delta_prime_e(k, n), "s")

    def delta_ek(self, k: int, f: NablaArg) -> Param:
        """`Δ_{e_k} F`, with eigenvalue `e_k[B_μ]`, in the Schur basis.

            >>> from symfn import macdonald, s
            >>> macdonald.delta_ek(1, s([1, 1]))
            (1 + t + q)*s[1,1] + s[2]

        The eigenvalue is `e_k[B_μ]` where `delta_prime_ek`'s is
        `e_k[B_μ − 1]`, and the constant term above is exactly that
        difference — the value that separates the primed operator from the
        unprimed one.

        # Raises

        Raises `ValueError` unless `F` is homogeneous in the Schur basis,
        with integer coefficients or coefficients in `q` and `t`.
        """
        return _qt_element(_c.delta_ek(k, _schur_rows(f)), "s")

    def delta_prime_ek(self, k: int, f: NablaArg) -> Param:
        """`Δ'_{e_k} F`, with eigenvalue `e_k[B_μ − 1]`, in the Schur basis.

            >>> from symfn import macdonald, s
            >>> macdonald.delta_prime_ek(1, s([1, 1])) == macdonald.delta_prime_e(1, 2)
            True

        `e_2 = s_{11}`, so this is the general form of `delta_prime_e` and
        agreeing with it is what pins the eigenvalue.

        # Raises

        Raises `ValueError` unless `F` is homogeneous in the Schur basis,
        with integer coefficients or coefficients in `q` and `t`.
        """
        return _qt_element(_c.delta_prime_ek(k, _schur_rows(f)), "s")

    def theta_ek(self, k: int, f: NablaArg) -> Param:
        """`Θ_{e_k} F`, which raises the degree by `k`, in the Schur basis.

            >>> from symfn import macdonald, s
            >>> macdonald.theta_ek(1, s([1, 1])).support()
            [(1, 1, 1), (2, 1)]

        The support sits one box higher than the argument's, which is what
        distinguishes `Θ` from the `Δ` operators — those preserve degree.

        # Raises

        Raises `ValueError` unless `F` is homogeneous in the Schur basis,
        with integer coefficients or coefficients in `q` and `t`.
        """
        return _qt_element(_c.theta_ek(k, _schur_rows(f)), "s")

    def big_pi(self, f: NablaArg) -> Param:
        """`Π F`, with eigenvalue `Π_μ`, in the Schur basis.

            >>> from symfn import macdonald, s
            >>> macdonald.big_pi(s([1, 1])).coefficient([2])
            -1

        # Raises

        Raises `ValueError` unless `F` is homogeneous in the Schur basis,
        with integer coefficients or coefficients in `q` and `t`.
        """
        return _qt_element(_c.big_pi(_schur_rows(f)), "s")

    def __repr__(self) -> str:
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

    def P(self, la: PartitionArg) -> Param:
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

    def Q(self, la: PartitionArg) -> Param:
        """`Q_λ = b_λ·P_λ`, in the monomial basis.

            >>> from symfn import jack
            >>> jack.Q([2]).coefficient([2])
            (1 + alpha)/(2*alpha^2)

        `b_(2) = (1 + α)/(2α²)` is `1/⟨P_(2), P_(2)⟩_α`, the norm `Q` divides
        out; `P_(2)` carries 1 in this position.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _jack_element(_c.jack_q(_partition(la)))

    def J(self, la: PartitionArg) -> Param:
        """`J_λ(x; α)`, the integral form, in the monomial basis.

            >>> from symfn import jack
            >>> jack.J([1, 1])
            2*m[1,1]

        The α-free coefficients are the integral form's signature.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _jack_element(_c.jack_j(_partition(la)))

    def zonal(self, la: PartitionArg, integral_form: bool = False) -> Sym:
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

    def structure_constant(
        self, la: PartitionArg, mu: PartitionArg, nu: PartitionArg
    ) -> AlphaFrac:
        """`⟨J_λ J_μ, J_ν⟩_α`, Stanley's object, as an `AlphaFrac`.

            >>> from symfn import jack
            >>> jack.structure_constant([1], [1], [2])
            2*alpha^2

        Whether every such value lies in `ℕ[α]` is Stanley's 1989 conjecture
        and is open; a coefficient that came back negative would be a result
        to report rather than a defect to fix.

        # Raises

        Raises `ValueError` unless all three are partitions and
        `|λ| + |μ| = |ν|`.
        """
        num, atoms, scale = _c.jack_structure_constant(
            _partition(la), _partition(mu), _partition(nu)
        )
        return AlphaFrac(num, atoms, scale)

    def __repr__(self) -> str:
        return "symfn.jack"


class _HallLittlewood:
    """The Hall-Littlewood family in `t`, and the Kostka-Foulkes polynomials.

    `Qp` is `Q'_λ`, the one whose Schur coefficients are the Kostka-Foulkes
    polynomials `K_{μλ}(t)`; `P` is `P_λ`, monic in the monomial basis. Both
    come back expanded in Schur functions; `to_P` and `to_Qp` go the other
    way, rewriting a Schur-basis element in the `P` or `Q'` basis.

        >>> from symfn import hl
        >>> hl.Qp([2, 1]).at(t=0)
        s[2,1]
        >>> hl.Qp([2, 1]).at(t=1)
        s[2,1] + s[3]

    `Q'_λ(x; 0) = s_λ` and `Q'_λ(x; 1) = h_λ` expanded in Schur functions,
    whose coefficients are the Kostka numbers. Both are theorems, and both
    fail under the `t → 1/t` twist that the other convention in circulation
    uses. Sage's equivalent is `Sym.hall_littlewood().Qp()`.
    """

    __module__ = "symfn"

    def Qp(self, la: PartitionArg) -> Param:
        """`Q'_λ(x; t) = Σ_μ K_{μλ}(t) s_μ`, in the Schur basis.

            >>> from symfn import hl
            >>> hl.Qp([1, 1])
            s[1,1] + t*s[2]

        `t` on `s_2` rather than on `s_11` is the charge convention:
        `K_{(2),(11)}(t) = t` where the cocharge rival puts the `t` on the
        diagonal term.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _t_element(_c.hall_littlewood(_partition(la)), "s")

    def P(self, la: PartitionArg) -> Param:
        """`P_λ(x; t)` in the Schur basis.

            >>> from symfn import hl
            >>> hl.P([2])
            -t*s[1,1] + s[2]

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _t_element(_c.hall_littlewood_p(_partition(la)), "s")

    def to_P(self, f: TSchurArg) -> Param:
        """`f`, given in the Schur basis, rewritten in the `P` basis, as a
        `Param` tagged `HLP`.

            >>> from symfn import hl, s
            >>> hl.to_P(s([2]))
            t*HLP[1,1] + HLP[2]
            >>> hl.to_P(hl.P([3, 1]))
            HLP[3,1]

        `s_2 = P_2 + t·P_11`: the coefficient of `P_λ` in `s_μ` is the
        Kostka-Foulkes polynomial `K_{μλ}(t)`, so the `t` sits on the
        dominance-smaller shape. Sage's `HLP(s[2])` prints the same value.
        `f` may be a `Sym`, a `Param` in `t` (so a `P` or `Qp` value feeds
        back in, as the second example does), or the contract layer's rows;
        rational coefficients are scaled through the boundary and restored.

        # Raises

        Raises `ValueError` unless `f` is in the Schur basis with coefficients
        in `t` alone, and unless every support is a partition.
        """
        rows, scale = _t_schur_rows(f, "to_P")
        return _t_element(_c.schur_to_hall_littlewood_p(rows), "HLP", scale)

    def to_Qp(self, f: TSchurArg) -> Param:
        """`f`, given in the Schur basis, rewritten in the `Q'` basis, as a
        `Param` tagged `HLQp`.

            >>> from symfn import hl, s
            >>> hl.to_Qp(s([1, 1]))
            HLQp[1,1] - t*HLQp[2]
            >>> hl.to_Qp(hl.Qp([2, 1]))
            HLQp[2,1]

        `s_11 = Q'_11 − t·Q'_2`, which is `Q'_11 = s_11 + t·s_2` read
        backwards, and Sage's `HLQp(s[1,1])` prints the same value. Against
        `to_P`: there the `t` lands on the smaller shape with a plus sign,
        here on the larger one with a minus, so a swap of the two
        normalizations shows in the sign. Accepts what `to_P` accepts.

        # Raises

        Raises `ValueError` on the same conditions as `to_P`.
        """
        rows, scale = _t_schur_rows(f, "to_Qp")
        return _t_element(_c.schur_to_hall_littlewood_qp(rows), "HLQp", scale)

    def kostka_foulkes(self, la: PartitionArg, mu: PartitionArg) -> Poly:
        """The Kostka-Foulkes polynomial `K_{λμ}(t)`, as a `Poly`.

            >>> from symfn import hl
            >>> hl.kostka_foulkes([2, 2], [1, 1, 1, 1])
            t^2 + t^4

        `K_{λμ}(1)` is the Kostka number, and the absence of a constant term
        distinguishes the charge convention from its cocharge rival.

        # Raises

        Raises `ValueError` unless both arguments are partitions of the same
        integer.
        """
        return Poly("t", _c.kostka_foulkes(_partition(la), _partition(mu)))

    def __repr__(self) -> str:
        return "symfn.hl"


class _LLT:
    """The LLT family in `q`: ribbon tableaux on one side, tuples on the other.

    One family, two presentations. `G` takes a tuple of skew shapes and sums
    `q^{inv(T)} x^T` over semistandard fillings, `inv` counting attacking
    pairs that are out of order. `Gtilde`, `Htilde` and `H` take a partition
    and a level `k` and sum over `k`-ribbon tableaux: `H` grades by spin,
    `Gtilde` and `Htilde` by cospin, and `Htilde(mu, k)` is
    `Gtilde(k·mu, k)`. Every entry point here returns the **monomial** basis;
    `schur` is the one that converts.

        >>> from symfn import llt
        >>> llt.H([1, 1], 2)
        (1 + q)*m[1,1] + q*m[2]

    The conventions in circulation differ by more than a twist, so each
    method states which function it computes. Sage's dictionary, with the
    grading variable named `t` there and `q` here: `Sym.llt(k).hspin()` is
    `H`, `hcospin()` is `Htilde`, and `cospin()` on a partition is `Gtilde`.
    On a tuple, Sage's `cospin()` divides out the floor `q^{min_inv(...)}`
    that `G` deliberately keeps.
    """

    __module__ = "symfn"

    def Gtilde(self, la: PartitionArg, k: int) -> Param:
        """`G̃^(k)_λ`, cospin-graded over `k`-ribbon tableaux of λ, in the
        monomial basis.

            >>> from symfn import llt
            >>> llt.Gtilde([2, 2], 2)
            (1 + q)*m[1,1] + m[2]
            >>> llt.Gtilde([2, 1], 2)
            0

        A shape with a nonempty `k`-core admits no ribbon tiling at all, so
        the sum is empty — `(2, 1)` above — rather than an error.

        # Raises

        Raises `ValueError` unless λ is a partition and `k` is positive.
        """
        return _q_element(_c.llt_gtilde(_partition(la), k), "m")

    def H(self, mu: PartitionArg, k: int) -> Param:
        """`H^(k)_μ`, spin-graded over `k`-ribbon tableaux of `k·μ`, in the
        monomial basis.

            >>> from symfn import llt
            >>> llt.H([1, 1], 2)
            (1 + q)*m[1,1] + q*m[2]

        `H = q^{s*}·H̃(x; 1/q)`: spin where `Htilde` is cospin, and the
        reversal is visible above — the `q` on `m_2` sits where
        `Htilde([1, 1], 2)` carries the constant.

        # Raises

        Raises `ValueError` unless μ is a partition and `k` is positive.
        """
        return _q_element(_c.llt_h(_partition(mu), k), "m")

    def G(
        self,
        shapes: Sequence[PartitionArg],
        offsets: Sequence[int] | None = None,
    ) -> Param:
        """`G_ν(x; q) = Σ_T q^{inv(T)} x^T` over a tuple of shapes, in the
        monomial basis.

            >>> from symfn import llt
            >>> llt.G([[1], [1]])
            (1 + q)*m[1,1] + m[2]

        Read in the Schur basis this is `s_2 + q·s_11`: `q` on `s_11` rather
        than on `s_2` is the inversion statistic's orientation, and what
        separates this from the `q → 1/q` convention. The raw `inv` grading
        is kept — `min_inv` gives the floor a comparison with Sage's
        `cospin()` divides out.

        # Raises

        Raises `ValueError` unless every shape is a partition, and unless
        `offsets` — when given — has one entry per shape.
        """
        rows = _c.llt_g([_partition(sh) for sh in shapes], offsets)
        return _q_element(rows, "m")

    def schur(self, la: PartitionArg, k: int) -> Param:
        """`G̃^(k)_λ` in the **Schur** basis, where `Gtilde` gives monomial.

            >>> from symfn import llt
            >>> llt.schur([1, 1], 2).basis
            's'

        # Raises

        Raises `ValueError` unless λ is a partition and `k ≥ 1`.
        """
        return _q_element(_c.llt_schur(_partition(la), k), "s")

    def Htilde(self, mu: PartitionArg, k: int) -> Param:
        """`H̃^(k)_μ = G̃^(k)_{kμ}`, the **cospin** family, monomial basis.

            >>> from symfn import llt
            >>> llt.Htilde([1], 2)
            m[1]

        `H` is spin and this is cospin; the two differ by `q^{s*} → q^{−s}`
        and disagree on any shape with a positive spin range.

        # Raises

        Raises `ValueError` unless μ is a partition and `k ≥ 1`.
        """
        return _q_element(_c.llt_h_tilde(_partition(mu), k), "m")

    def min_inv(
        self,
        shapes: Sequence[PartitionArg],
        offsets: Sequence[int] | None = None,
    ) -> int:
        """`min_T inv(T)` over the fillings of a tuple — the floor `G` leaves in.

            >>> from symfn import llt
            >>> llt.min_inv([[1], [1]])
            0
            >>> llt.min_inv([[1], [1, 1]])
            1

        `G` is deliberately not divided by `q` to this power, because the floor
        is real data about the shape tuple; Sage's `cospin` divides it out, so
        this is what a comparison needs. The second value shows the floor is
        forced: `((1), (11))` is the 2-quotient of `(2, 2, 2)`, and no offset
        choice brings its floor to zero.

        # Raises

        Raises `ValueError` unless every shape is a partition and `offsets`,
        when given, has one entry per shape.
        """
        return _c.llt_min_inv([_partition(sh) for sh in shapes], offsets)

    def __repr__(self) -> str:
        return "symfn.llt"


def _t_element(
    rows: Iterable[tuple[Partition, Iterable[tuple[int, Coefficient]]]],
    basis: str,
    scale: int = 1,
) -> Param:
    """Wrap `(partition, [(t_exponent, coefficient)])` rows as a `Param`,
    dividing every coefficient by `scale` — the `restore` half of the
    denominator round trip `_t_schur_rows` begins.
    """
    if scale == 1:
        return Param(basis, [(la, Poly("t", c)) for la, c in rows], ("t",))
    return Param(
        basis,
        [
            (la, Poly("t", [(k, Fraction(v, scale)) for k, v in c]))
            for la, c in rows
        ],
        ("t",),
    )


def _t_schur_rows(f: TSchurArg, what: str) -> tuple[list[Any], int]:
    """The `(partition, [(t_exponent, coefficient)])` rows the Hall-Littlewood
    inverse expansions take, and the integer the rows were scaled by.

    Accepts a `Sym` in the Schur basis, a `Param` in the Schur basis whose
    coefficients are polynomials in `t`, or the rows themselves. The contract
    layer takes integers, so rational coefficients are multiplied up by their
    least common denominator here and divided back out in `_t_element` —
    the same round trip `Sym.to` makes, exact because the expansion is
    linear.
    """
    if isinstance(f, Sym):
        if f.basis != "s":
            raise ValueError(f"{what} needs a Schur-basis element, not {f.basis}")
        polys = [(la, {0: c}) for la, c in f]
    elif isinstance(f, Param):
        if f.basis != "s":
            raise ValueError(f"{what} needs a Schur-basis element, not {f.basis}")
        polys = []
        for la, coeff in f:
            # The expansion is over ℤ[t]; a coefficient in `q` and `t`, or in
            # α, has the same term structure and a different meaning, so it is
            # refused rather than read through whichever accessor exists.
            if not isinstance(coeff, Poly) or coeff.variable != "t":
                raise ValueError(
                    f"{what} needs coefficients in t, not "
                    f"{type(coeff).__name__}"
                )
            polys.append((la, coeff.coefficients()))
    else:
        return list(f), 1
    scale = 1
    for _, terms in polys:
        for c in terms.values():
            if isinstance(c, Fraction):
                d = c.denominator
                scale = scale * d // gcd(scale, d)
    return [
        (la, [(k, int(c * scale)) for k, c in terms.items()]) for la, terms in polys
    ], scale


def _schur_rows(f: NablaArg, what: str = "nabla") -> list[Any]:
    """The `(partition, [(a, b, coefficient)])` rows `nabla` takes.

    Accepts a `Sym` in the Schur basis, a `Param` in `q` and `t`, or the rows
    themselves.
    """
    if isinstance(f, Sym):
        if f.basis != "s":
            raise ValueError(f"{what} needs a Schur-basis element, not {f.basis}")
        out: list[Any] = []
        for la, c in f:
            # `int()` on a `Fraction` truncates, so a rational coefficient
            # would cross as a different element rather than as an error.
            # These rows are ℤ[q,t]; refusing is the boundary's contract
            # (`docs/policies/failure.md`, P8).
            if c != int(c):
                raise ValueError(
                    f"{what} needs integer coefficients; {tuple(la)} carries {c}"
                )
            out.append((la, [(0, 0, int(c))]))
        return out
    if isinstance(f, Param):
        if f.basis != "s":
            raise ValueError(f"{what} needs a Schur-basis element, not {f.basis}")
        rows = []
        for la, coeff in f:
            # `∇` takes `(q, t)`-graded rows, and only a `QtPoly` coefficient
            # has them. A Schur-basis `Param` over any other coefficient type
            # has the same term structure and different coefficients, so it is
            # refused rather than read through whichever accessor exists.
            if not isinstance(coeff, QtPoly):
                raise ValueError(
                    f"{what} needs coefficients in q and t, not "
                    f"{type(coeff).__name__}"
                )
            rows.append(
                (la, [(a, b, c) for (a, b), c in coeff.coefficients().items()])
            )
        return rows
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
