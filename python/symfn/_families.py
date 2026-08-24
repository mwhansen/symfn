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

from collections.abc import Callable, Iterable, Sequence
from fractions import Fraction
from math import gcd
from typing import Any, Union

from . import symfn as _c
from ._bases import BASES
from ._param import (
    AlphaFrac,
    Param,
    ParamCoefficient,
    Poly,
    QtFrac,
    QtPoly,
    QtRatio,
)
from ._sym import Sym, _partition
from ._types import Coefficient, Partition, PartitionArg

#: What `nabla` accepts: an element in the Schur basis in either of the two
#: types that can carry one, or the contract layer's rows themselves.
NablaArg = Union["Sym", Param, Iterable[Any]]

#: What the contract layer's `(q, t)`-graded rows look like on arrival.
QtRows = Iterable[tuple[Partition, Iterable[tuple[int, int, Coefficient]]]]

#: What a `Param` may be multiplied by: an integer, a rational, a polynomial in
#: its own parameters, or a coefficient of the kind it carries.
Scalar = Union[int, Fraction, "Poly", "QtPoly", "QtFrac", "AlphaFrac"]

#: What the Hall-Littlewood inverse expansions accept: a Schur-basis element
#: as a `Sym`, as a `Param` in `t`, or as the contract layer's `t`-rows.
TSchurArg = Union["Sym", Param, Iterable[Any]]

#: What the Macdonald inverse expansions take: a `Sym` or a `Param` in the
#: monomial basis, or the contract layer's `(partition, numerator,
#: denominator)` rows.
MacArg = Union["Sym", Param, Iterable[Any]]

#: What the Jack inverse expansions take: a `Sym` or a `Param` in the monomial
#: basis, or the contract layer's `(partition, numerator, atoms, scale)` rows.
JackArg = Union["Sym", Param, Iterable[Any]]

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


def _ht_element(rows: Iterable[Any], basis: str = "McdHt") -> Param:
    """Wrap `(partition, numerator, denominator atoms)` rows as a `Param` in
    `q` and `t`, tagged as Sage prints it.

    The expansion out of `H̃` lands in the Schur basis with the same
    coefficients, so this builds both ends of the pair.
    """
    return Param(basis, [(la, QtRatio(n, d)) for la, n, d in rows], ("q", "t"))


def _mac_element(rows: Iterable[Any], basis: str = "m", scale: int = 1) -> Param:
    """Wrap `(partition, numerator, denominator)` rows as a `Param` in q, t,
    dividing every numerator coefficient by `scale` — the `restore` half of
    the denominator round trip `_mac_rows` begins.
    """
    if scale == 1:
        return Param(basis, [(la, QtFrac(n, d)) for la, n, d in rows], ("q", "t"))
    return Param(
        basis,
        [
            (la, QtFrac([(a, b, Fraction(v, scale)) for a, b, v in n], d))
            for la, n, d in rows
        ],
        ("q", "t"),
    )


def _jack_element(rows: Iterable[Any], basis: str = "m") -> Param:
    """Wrap `(partition, numerator, atoms, scale)` rows as a `Param` in α."""
    return Param(
        basis, [(la, AlphaFrac(n, d, k)) for la, n, d, k in rows], ("alpha",)
    )


def _unit(basis: str, la: PartitionArg) -> Param:
    """The single basis element `X_λ` in the parametric basis `basis`, with
    coefficient 1.

    Every family's forward constructor returns one of these. The expansion is
    a `to` away and is not computed here, which is the whole point: a shape is
    what the caller named, and `P_λ` is the object, not its coordinates.
    """
    key = _partition(la)
    if basis.startswith("Mcd") and basis != "McdHt":
        return _mac_element([(key, [(0, 0, 1)], [])], basis)
    if basis == "McdHt":
        return _ht_element([(key, [(0, 0, 1)], [])])
    if basis.startswith("Jack"):
        return _jack_element([(key, [1], [], 1)], basis)
    return _t_element([(key, [(0, 1)])], basis)


class _Macdonald:
    """The Macdonald family in `q` and `t`, and the operators built on it.

    `P` is monic in the monomial basis, `Q = b_λ·P`, and `J = c_λ·P` is the
    integral form whose coefficients are polynomials. `Htilde` is the modified
    form `H̃_μ`, the one the `(q,t)`-Kostka polynomials expand.

        >>> from symfn import macdonald
        >>> macdonald.P([1])
        McdP[1]
        >>> macdonald.Q([1]).to("m")
        (1 - t)/(1 - q)*m[1]

    Each of the four names a shape in its own basis, and `to` expands it —
    into the monomial basis for `P`, `Q` and `J`, into Schur for `H̃`, because
    that is the basis each is defined in. `to_P`, `to_Q`, `to_J` and
    `to_Htilde` run the other way, rewriting a classical element *into* one of
    the four, which is the direction a positivity question asks in. Each takes
    the basis its forward sibling expands in.

    Sage's equivalents are `Sym.macdonald().P()`, `.Q()`, `.J()` and `.Ht()`,
    called on a shape or on an element respectively.
    """

    __module__ = "symfn"

    def P(self, la: PartitionArg) -> Param:
        """`P_λ(x; q, t)` in the monomial basis, monic in `m_λ`.

            >>> from symfn import macdonald
            >>> macdonald.P([2])
            McdP[2]
            >>> macdonald.P([2]).to("m")
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
        return _unit("McdP", la)

    def Q(self, la: PartitionArg) -> Param:
        """`Q_λ = b_λ · P_λ`, in the monomial basis.

            >>> from symfn import macdonald
            >>> macdonald.Q([1]).to("m").coefficient([1])
            (1 - t)/(1 - q)

        `Q_(1) = (1 − t)/(1 − q)·m_1` where `P_(1) = m_1`: the value that
        separates the two normalizations at the smallest shape.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _unit("McdQ", la)

    def J(self, la: PartitionArg) -> Param:
        """`J_λ = c_λ · P_λ`, the integral form, in the monomial basis.

            >>> from symfn import macdonald
            >>> macdonald.J([1, 1]).to("m").coefficient([1, 1])
            1 - t - t^2 + t^3

        Every coefficient is a polynomial — the empty denominator is the
        integral form's signature, since `J` clears what `Q` carries.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _unit("McdJ", la)

    def Htilde(self, mu: PartitionArg) -> Param:
        """The modified Macdonald polynomial `H̃_μ`, in the Schur basis.

            >>> from symfn import macdonald
            >>> macdonald.Htilde([2])
            McdHt[2]
            >>> macdonald.Htilde([2]).to("s")
            q*s[1,1] + s[2]

        The Schur coefficients are the `(q,t)`-Kostka polynomials `K̃_{λμ}`,
        and `H̃_(2)` carrying `q` rather than `t` on `s_11` is the orientation
        of the Garsia-Haiman convention.

        # Raises

        Raises `ValueError` unless μ is a partition.
        """
        return _unit("McdHt", mu)

    def to_Htilde(self, f: NablaArg) -> Param:
        """`f`, given in the Schur basis, rewritten in the `H̃` basis — the
        direction `Htilde` does not go.

            >>> from symfn import macdonald, s
            >>> macdonald.to_Htilde(s([2]))
            q/(q - t)*McdHt[1,1] - t/(q - t)*McdHt[2]
            >>> macdonald.to_Htilde(s([2])).to("s")
            s[2]
            >>> macdonald.to_Htilde(macdonald.Htilde([2, 1]).to("s"))
            McdHt[2,1]

        `s_2 = q/(q−t)·H̃_11 − t/(q−t)·H̃_2`; the `q` upstairs on the column
        shape is the orientation, and the `q ↔ t` swap gives a different
        answer rather than an error. The coefficients are `QtRatio`s and
        genuinely not polynomials: this divides by `w_μ`, whose factors are
        `q^a − t^b` and do not cancel. They cross **factored**, so `to` runs
        the expansion back and the round trip closes, as the second example
        shows.

        Accepts what `nabla` accepts — a `Sym` or a `Param` in the Schur
        basis, or the contract rows — and unlike `nabla` it takes mixed
        degrees, expanding each degree on its own.

        # Raises

        Raises `ValueError` unless the argument is in the Schur basis with
        integer coefficients, or coefficients in `q` and `t`.
        """
        return _ht_element(_c.schur_to_macdonald_ht(_schur_rows(f, "to_Htilde")))

    def to_J(self, f: NablaArg) -> Param:
        """`f`, given in the Schur basis, rewritten in the `J` basis.

            >>> from symfn import macdonald, s
            >>> macdonald.to_J(s([1, 1]))
            1/((1 - t)*(1 - t^2))*McdJ[1,1]
            >>> macdonald.to_J(s([2])).support()
            [(1, 1), (2,)]
            >>> macdonald.to_J(s([2])).coefficient([2])
            1/((1 - t)*(1 - q*t))

        `s_11 = J_11/((1−t)(1−t²))` and nothing else, where `s_2` reaches both
        shapes — the triangularity, which runs the opposite way from `J → s`.
        The denominator is the hook product `c_μ`, not `c'_μ = (1−q)(1−q²)`,
        which is the twist to check. `J` is the integral form, so the forward
        direction has polynomial coefficients and this one does not.

        Accepts what `nabla` accepts, and unlike `nabla` it takes mixed
        degrees, expanding each degree on its own. `schur_in_macdonald_j`
        is the whole-degree table this reads.

        # Raises

        Raises `ValueError` unless the argument is in the Schur basis with
        integer coefficients, or coefficients in `q` and `t`.
        """
        return _mac_element(
            _c.schur_to_macdonald_j(_schur_rows(f, "to_J")), "McdJ"
        )

    def to_P(self, f: MacArg) -> Param:
        """`f`, given in the monomial basis, rewritten in the `P` basis, as a
        `Param` tagged `McdP`.

            >>> from symfn import macdonald, m
            >>> macdonald.to_P(m([2])).coefficient([1, 1])
            (-1 + t - q + q*t)/(1 - q*t)
            >>> macdonald.to_P(macdonald.P([2, 1]).to("m"))
            McdP[2,1]

        `m_2 = P_2 − [(1−t)(1+q)/(1−q·t)] P_11`: the coefficient `P → m` puts
        on the dominance-smaller shape, negated. The `q ↔ t` swap gives
        `(1−q)(1+t)/(1−q·t)` instead, which is the twist to check.

        `f` may be a `Sym` or a `Param` in the monomial basis — so a `P`, `Q`
        or `J` value feeds back in, as the second example does — or the
        contract layer's rows; rational coefficients are scaled through the
        boundary and restored. An element in another classical basis is
        refused rather than converted, on the same grounds as `BasisError`:
        write `macdonald.to_P(f.to("m"))` and the conversion is the caller's,
        with its cost visible.

        # Raises

        Raises `ValueError` unless `f` is in the monomial basis with
        coefficients in `q` and `t`, and unless every support is a partition.
        """
        rows, scale = _mac_rows(f, "to_P")
        return _mac_element(_c.monomial_to_macdonald_p(rows), "McdP", scale)

    def to_Q(self, f: MacArg) -> Param:
        """`f`, given in the monomial basis, rewritten in the `Q` basis, as a
        `Param` tagged `McdQ`.

            >>> from symfn import macdonald, m
            >>> macdonald.to_Q(m([1, 1]))
            (1 - q - q*t + q^2*t)/((1 - t)*(1 - t^2))*McdQ[1,1]
            >>> macdonald.to_P(m([1, 1]))
            McdP[1,1]

        `Q_λ = b_λ P_λ`, so this is `to_P` with each coefficient divided by
        that shape's `b_λ`. `m_11 = P_11` outright where the `Q` coefficient
        is `(1−q·t)(1−q)/((1−t)(1−t²))` — the value that separates the two
        normalizations at the smallest shape where they differ. Accepts what
        `to_P` accepts.

        # Raises

        Raises `ValueError` on the same conditions as `to_P`.
        """
        rows, scale = _mac_rows(f, "to_Q")
        return _mac_element(_c.monomial_to_macdonald_q(rows), "McdQ", scale)

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

    Each names a shape in its own basis and `to("m")` expands it, since the
    monomial basis is where all three are defined. At `α = 1` every one of
    them degenerates to a Schur function, and at `α = 2` to a zonal
    polynomial; both are checks a caller can run, and `at` expands first so
    neither needs a conversion written out.

        >>> from symfn import jack
        >>> jack.P([2])
        JackP[2]
        >>> jack.P([2]).at(alpha=1)
        m[1,1] + m[2]

    `to_P`, `to_Q` and `to_J` run the other way, rewriting a monomial-basis
    element into one of the three. ⚠️ The `α = 1` check above cannot tell them
    apart on the way back — `m_11` is `P_11` outright and `[α(α+1)/2] Q_11`,
    and both are 1 there — so each of those methods carries a value at free
    `α` instead.

    Sage's equivalents are `Sym.jack().P()`, `.Q()` and `.J()`, called on a
    shape or on an element respectively.
    """

    __module__ = "symfn"

    def P(self, la: PartitionArg) -> Param:
        """`P_λ(x; α)` in the monomial basis, monic in `m_λ`.

            >>> from symfn import jack
            >>> jack.P([2])
            JackP[2]
            >>> jack.P([2]).to("m")
            2/(alpha + 1)*m[1,1] + m[2]

        The leading 1 is what separates `P` from `Q` and `J`, both of which
        scale it.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _unit("JackP", la)

    def Q(self, la: PartitionArg) -> Param:
        """`Q_λ = b_λ·P_λ`, in the monomial basis.

            >>> from symfn import jack
            >>> jack.Q([2]).to("m").coefficient([2])
            (1 + alpha)/(2*alpha^2)

        `b_(2) = (1 + α)/(2α²)` is `1/⟨P_(2), P_(2)⟩_α`, the norm `Q` divides
        out; `P_(2)` carries 1 in this position.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _unit("JackQ", la)

    def J(self, la: PartitionArg) -> Param:
        """`J_λ(x; α)`, the integral form, in the monomial basis.

            >>> from symfn import jack
            >>> jack.J([1, 1])
            JackJ[1,1]
            >>> jack.J([1, 1]).to("m")
            2*m[1,1]

        The α-free coefficients are the integral form's signature.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _unit("JackJ", la)

    def to_P(self, f: JackArg) -> Param:
        """`f`, given in the monomial basis, rewritten in the `P` basis, as a
        `Param` tagged `JackP`.

            >>> from symfn import jack, m
            >>> jack.to_P(m([2])).coefficient([1, 1])
            -2/(alpha + 1)
            >>> jack.to_P(jack.P([2, 1]).to("m"))
            JackP[2,1]

        `m_2 = P_2 − [2/(α+1)] P_11`: the coefficient `P → m` puts on the
        dominance-smaller shape, negated. Sending `α → 1/α` would give
        `−2α/(α+1)` instead, which is the twist to check; the two agree at
        `α = 1`, so that specialization cannot see it.

        `f` may be a `Sym` or a `Param` in the monomial basis — so a `P`, `Q`
        or `J` value feeds back in, as the second example does — or the
        contract layer's rows. An element in another classical basis is
        refused rather than converted, on the same grounds as `BasisError`:
        write `jack.to_P(f.to("m"))` and the conversion is the caller's, with
        its cost visible.

        # Raises

        Raises `ValueError` unless `f` is in the monomial basis with
        coefficients in α, and unless every support is a partition.
        """
        return _jack_element(_c.monomial_to_jack_p(_jack_rows(f, "to_P")), "JackP")

    def to_Q(self, f: JackArg) -> Param:
        """`f`, given in the monomial basis, rewritten in the `Q` basis, as a
        `Param` tagged `JackQ`.

            >>> from symfn import jack, m
            >>> jack.to_Q(m([1, 1]))
            (alpha + alpha^2)/2*JackQ[1,1]
            >>> jack.to_P(m([1, 1]))
            JackP[1,1]

        `Q_λ = (H_λ/H'_λ)·P_λ`, so this is `to_P` with each coefficient
        multiplied by that shape's `⟨P_λ, P_λ⟩_α`. `m_11 = P_11` outright
        where the `Q` coefficient is `α(α+1)/2` — the smallest shape at which
        the two normalizations differ. At `α = 1` both are 1, so setting the
        parameter cannot tell them apart either. Accepts what `to_P` accepts.

        # Raises

        Raises `ValueError` on the same conditions as `to_P`.
        """
        return _jack_element(_c.monomial_to_jack_q(_jack_rows(f, "to_Q")), "JackQ")

    def to_J(self, f: JackArg) -> Param:
        """`f`, given in the monomial basis, rewritten in the `J` basis, as a
        `Param` tagged `JackJ`.

            >>> from symfn import jack, m
            >>> jack.to_J(m([2])).coefficient([1, 1])
            -1/(alpha + 1)

        `J_λ = H_λ·P_λ`, so this is `to_P` with each coefficient divided by
        that shape's lower hooks. `m_2 = [J_2 − J_11]/(α+1)`, which is
        `J_(2) = (α+1)·m_2 + 2·m_11` and `J_(1,1) = 2·m_11` read backwards.
        Unlike `J` itself the coefficients are not polynomials in α: only the
        forward direction is integral. Accepts what `to_P` accepts.

        # Raises

        Raises `ValueError` on the same conditions as `to_P`.
        """
        return _jack_element(_c.monomial_to_jack_j(_jack_rows(f, "to_J")), "JackJ")

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
    polynomials `K_{μλ}(t)`; `P` is `P_λ`, monic in the monomial basis. Each
    names a shape in its own basis and `to("s")` expands it, since Schur is
    where both are defined here; `to_P` and `to_Qp` go the other way,
    rewriting a Schur-basis element in the `P` or `Q'` basis.

        >>> from symfn import hl
        >>> hl.Qp([2, 1])
        HLQp[2,1]
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
            HLQp[1,1]
            >>> hl.Qp([1, 1]).to("s")
            s[1,1] + t*s[2]

        `t` on `s_2` rather than on `s_11` is the charge convention:
        `K_{(2),(11)}(t) = t` where the cocharge rival puts the `t` on the
        diagonal term.

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _unit("HLQp", la)

    def P(self, la: PartitionArg) -> Param:
        """`P_λ(x; t)` in the Schur basis.

            >>> from symfn import hl
            >>> hl.P([2])
            HLP[2]
            >>> hl.P([2]).to("s")
            -t*s[1,1] + s[2]

        # Raises

        Raises `ValueError` unless λ is a partition.
        """
        return _unit("HLP", la)

    def to_P(self, f: TSchurArg) -> Param:
        """`f`, given in the Schur basis, rewritten in the `P` basis, as a
        `Param` tagged `HLP`.

            >>> from symfn import hl, s
            >>> hl.to_P(s([2]))
            t*HLP[1,1] + HLP[2]
            >>> hl.to_P(hl.P([3, 1]).to("s"))
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
            >>> hl.to_Qp(hl.Qp([2, 1]).to("s"))
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
    `Gtilde(k·mu, k)`. Every entry point here returns the **monomial** basis
    directly, unlike the other three families, which name a shape in a basis
    of their own; `schur` is the one that converts.

        >>> from symfn import llt
        >>> llt.H([1, 1], 2)
        (1 + q)*m[1,1] + q*m[2]

    The conventions in circulation differ by more than a twist, so each
    method states which function it computes. Sage's dictionary, with the
    grading variable named `t` there and `q` here: `Sym.llt(k).hspin()` is
    `H`, `hcospin()` is `Htilde`, and `cospin()` on a partition is `Gtilde`.
    On a tuple, Sage's `cospin()` divides out the floor `q^{min_inv(...)}`
    that `G` deliberately keeps.

    Unlike the other three families this one has no `to_*` method, and cannot:
    the `G̃^{(k)}_λ` are indexed by tuples and by the level `k`, and are
    linearly dependent across `k`, so they are not a basis of Λ and "expand in
    the LLT basis" names no unique answer. `schur` is the conversion that does
    exist.
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


#: Where each parametric basis expands: the classical basis its family's
#: forward direction is defined in. `Param.to` reads this, and `Param.at`
#: goes through it, since a `Sym` carries only the six classical bases.
EXPANDS_IN: dict[str, str] = {
    "McdP": "m",
    "McdQ": "m",
    "McdJ": "m",
    "JackP": "m",
    "JackQ": "m",
    "JackJ": "m",
    "HLP": "s",
    "HLQp": "s",
    "McdHt": "s",
}


def _expand(f: Param) -> Param:
    """A `Param` in a parametric basis, expanded in the basis `EXPANDS_IN`
    names for it.

    One contract call and no arithmetic here, which is what the two directions
    sharing a row encoding buys: the coefficients go back out in the shape they
    arrived in.
    """
    tag = f.basis
    if tag in ("McdP", "McdQ", "McdJ"):
        rows, scale = _mac_rows(f, "to", tag)
        out = {
            "McdP": _c.macdonald_p_to_monomial,
            "McdQ": _c.macdonald_q_to_monomial,
            "McdJ": _c.macdonald_j_to_monomial,
        }[tag](rows)
        return _mac_element(out, "m", scale)
    if tag in ("JackP", "JackQ", "JackJ"):
        jack_rows = _jack_rows(f, "to", tag)
        jack_out = {
            "JackP": _c.jack_p_to_monomial,
            "JackQ": _c.jack_q_to_monomial,
            "JackJ": _c.jack_j_to_monomial,
        }[tag](jack_rows)
        return _jack_element(jack_out, "m")
    if tag in ("HLP", "HLQp"):
        t_rows, t_scale = _t_schur_rows(f, "to", tag)
        t_out = {
            "HLP": _c.hall_littlewood_p_to_schur,
            "HLQp": _c.hall_littlewood_qp_to_schur,
        }[tag](t_rows)
        return _t_element(t_out, "s", t_scale)
    rows = _c.macdonald_ht_to_schur(_ht_rows(f, "to", tag))
    # `K̃_{λμ}` is a polynomial, so an `H̃` element with polynomial
    # coefficients expands to one — and the Schur-basis element it becomes
    # should be the same kind every other Schur-basis value in `q` and `t` is,
    # so that `nabla`, `to_J` and `to_Htilde` take it without a conversion. A
    # genuine ratio survives only when one went in.
    if all(not atoms for _, _, atoms in rows):
        return _qt_element([(la, num) for la, num, _ in rows], "s")
    return _ht_element(rows, "s")


def _qt_pack(f: Param) -> tuple[list[Any], int]:
    """An element's coefficients as integer `(a, b, value)` exponent rows, and
    the integer they were scaled by.

    The polynomial encoding is what the `(q,t)` entry points read, and a
    one-variable `Poly` uses the second slot for whatever its variable is —
    the entry points never ask what the exponents count.

    # Raises

    Raises `ValueError` for a coefficient class that is not a polynomial.
    """
    kind = _kind(f)
    if kind is Poly:
        return _qt_int_rows([(la, _one_variable_row(c)) for la, c in f])
    if kind is QtPoly:
        return _qt_int_rows([(la, _two_variable_row(c)) for la, c in f])
    raise ValueError(
        f"this route reads the polynomial encoding, which "
        f"{kind.__name__ if kind else 'an empty element'} is not"
    )


def _qt_unpack(rows: list[Any], like: Param, dst: str, scale: int) -> Param:
    """Exponent rows back as a `Param` in `dst`, in the coefficient class and
    parameters `like` carries, dividing out the scale `_qt_pack` cleared.
    """
    if _kind(like) is Poly:
        var = next(c.variable for _, c in like if isinstance(c, Poly))
        return Param(
            dst,
            [
                (la, Poly(var, [(b, _unscale(v, scale)) for _, b, v in c]))
                for la, c in rows
            ],
            like.parameters,
        )
    return Param(
        dst,
        [
            (la, QtPoly([(a, b, _unscale(v, scale)) for a, b, v in c]))
            for la, c in rows
        ],
        like.parameters,
    )


def _carry_qt(f: Param, dst: str, call: Callable[[list[Any]], list[Any]]) -> Param:
    """An element's exponent rows through one contract call, rebuilt in `dst`
    in the coefficient class they went in as.

    The shared half of every operation this layer sends through the polynomial
    encoding: pack, clear denominators, call, restore, rebuild. Nothing here is
    arithmetic on the element — `call` is the whole computation.
    """
    rows, scale = _qt_pack(f)
    return _qt_unpack(call(rows), f, dst, scale)


def _convert(f: Param, dst: str) -> Param:
    """A `Param` in a classical basis, rewritten in another classical basis.

    One contract call and no arithmetic here. A basis change is a ℤ-linear map
    on the partitions, so a coefficient that carries a parameter crosses it
    intact and comes back in the class it went in as.

    Which entry point runs is decided by that class, because each coefficient
    ring has its own encoding: polynomials cross as exponent rows, the
    Macdonald and Jack families as a numerator over factored denominator
    atoms, and `H̃` as its own atoms. The four are the same routing over four
    rings — `crate::convert_named` is generic — so nothing about the basis pair
    differs between them.

    # Raises

    Raises `ValueError` if the element carries no coefficient class this layer
    knows.
    """
    src = f.basis
    if not len(f) or src == dst:
        return Param(dst, dict(f), f.parameters)
    kind = _kind(f)
    if kind in (Poly, QtPoly):
        return _carry_qt(f, dst, lambda rows: _c.convert_qt_terms(rows, src, dst))
    if kind is QtFrac:
        rows, scale = _mac_rows(f, "to", src)
        return _mac_element(_c.convert_macdonald_terms(rows, src, dst), dst, scale)
    if kind is AlphaFrac:
        jack_rows = _jack_rows(f, "to", src)
        return _jack_element(_c.convert_jack_terms(jack_rows, src, dst), dst)
    if kind is QtRatio:
        ht_rows = _ht_rows(f, "to", src)
        return _ht_element(_c.convert_ht_terms(ht_rows, src, dst), dst)
    raise ValueError(
        f"converting {src!r} to {dst!r} is not written for "
        f"{kind.__name__ if kind else 'these'} coefficients"
    )


#: The two Hopf operations this layer sends through the Schur basis, by name.
#: Both are defined there — ω conjugates the index, the antipode conjugates and
#: signs — and reach any other basis by a change of basis on each side.
_HOPF = {"omega": _c.omega_qt_terms, "antipode": _c.antipode_qt_terms}


def _hopf(f: Param, op: str) -> Param:
    """`omega` or `antipode`, returned in the basis it was given in.

    Three legs where `Sym` needs one: a parametric basis expands into its
    pivot, the pivot converts to Schur, and the result travels back the same
    way — so an element in a family's own basis comes back in it, as
    `Sym.omega` comes back in the basis it was handed.

    # Raises

    Raises `ValueError` for the coefficient classes the polynomial encoding
    does not carry, and for a parametric basis whose inverse expansion runs
    over them.
    """
    tag = f.basis
    schur = _convert(f if tag in BASES else _expand(f), "s")
    if _kind(schur) not in (Poly, QtPoly):
        raise ValueError(
            f"{op} is not written for {_kind(schur).__name__} coefficients "  # type: ignore[union-attr]
            "yet; it acts in the Schur basis, and the entry point that does so "
            "reads the polynomial encoding"
        )
    return _back_to(_carry_qt(schur, "s", _HOPF[op]), tag, op)


def _back_to(acted: Param, tag: str, what: str) -> Param:
    """A Schur-basis result rewritten in the basis the operand was written in.

    The last leg of every operation that leaves a parametric basis to compute:
    a classical tag is one more change of basis, and a parametric one is its
    family's inverse expansion.

    # Raises

    Raises `ValueError` for a parametric basis whose inverse expansion runs
    over coefficients this route does not carry.
    """
    if tag in BASES:
        return _convert(acted, tag)
    if tag == "HLP":
        return hl.to_P(acted)
    if tag == "HLQp":
        return hl.to_Qp(acted)
    if tag == "McdHt":
        return macdonald.to_Htilde(acted)
    raise ValueError(
        f"{what} in the {tag} basis needs an inverse expansion over "
        "rational-function coefficients, which is not written yet"
    )


def _product(f: Param, g: Param) -> Param:
    """`f*g`, in the basis both are written in.

    A product in a parametric basis is the product in the basis its family
    expands in, rewritten back: expand, convert to Schur, multiply there, and
    return the same way. The structure constants of the Schur product are the
    Littlewood-Richardson coefficients, which are integers and carry no
    parameter, so the coefficient ring is multiplied through rather than acted
    on.

    # Raises

    Raises `ValueError` if the two are in different bases or over different
    base rings, and for the coefficient classes the polynomial encoding does
    not carry.
    """
    if f.basis != g.basis:
        raise ValueError(
            f"cannot multiply an element in {f.basis} by one in {g.basis}; "
            "convert one with .to()"
        )
    if f.parameters != g.parameters:
        raise ValueError(
            f"cannot multiply an element in {_ring(f.parameters)} by one in "
            f"{_ring(g.parameters)}"
        )
    tag = f.basis
    if not len(f) or not len(g):
        return Param(tag, {}, f.parameters)
    left = _convert(f if tag in BASES else _expand(f), "s")
    right = _convert(g if tag in BASES else _expand(g), "s")
    rows_a, scale_a = _qt_pack(left)
    rows_b, scale_b = _qt_pack(right)
    product = _c.schur_multiply_qt(rows_a, rows_b)
    return _back_to(_qt_unpack(product, left, "s", scale_a * scale_b), tag, "a product")


def _unit_like(f: Param) -> Param:
    """The multiplicative identity in `f`'s basis, over `f`'s coefficient ring.

    The empty partition indexes 1 in every basis here — `P_∅`, `s_∅` and `m_∅`
    are all the constant 1 — so the unit differs only in which coefficient
    class carries it.
    """
    kind = _kind(f)
    one: ParamCoefficient
    if kind is Poly:
        one = Poly(next(c.variable for _, c in f if isinstance(c, Poly)), {0: 1})
    elif kind is QtPoly:
        one = QtPoly([(0, 0, 1)])
    else:
        raise ValueError(
            "the unit is not written for "
            f"{kind.__name__ if kind else 'an empty element'}"
        )
    # A pair rather than a mapping: `Mapping` is invariant in its value
    # type, so a dict of one coefficient class is not a dict of the union.
    return Param(f.basis, [((), one)], f.parameters)


def _ring(parameters: tuple[str, ...]) -> str:
    """The base ring a set of parameter names denotes, as an error names it."""
    return f"Q({', '.join(parameters)})" if parameters else "Q"


def _one_variable_row(c: ParamCoefficient) -> list[tuple[int, int, Coefficient]]:
    """A `Poly` in the two-exponent row encoding, its variable in the second
    slot — the one `t` occupies everywhere else at this boundary.

    Every term of an element shares a coefficient class, so the refusal is an
    invariant of `_convert`'s dispatch rather than a branch a caller reaches.
    """
    if not isinstance(c, Poly):
        raise TypeError(f"expected a polynomial, not {type(c).__name__}")
    return [(0, k, v) for k, v in c.coefficients().items()]


def _two_variable_row(c: ParamCoefficient) -> list[tuple[int, int, Coefficient]]:
    """A `QtPoly` in the same encoding, both exponents used. Same invariant."""
    if not isinstance(c, QtPoly):
        raise TypeError(f"expected a (q,t)-polynomial, not {type(c).__name__}")
    return _qt_rows(c)


def _qt_int_rows(rows: list[Any]) -> tuple[list[Any], int]:
    """Exponent rows with every coefficient an integer, and what they were
    scaled by.

    The boundary's coefficients are integers, and a `Poly` may hold a
    `Fraction` — `_t_element` divides by the denominator the Hall-Littlewood
    row builder cleared. Same clear-then-restore the `Sym` conversions use, one
    level further in: the scale is common to the whole element, so it survives
    a linear map and divides back out at the end.
    """
    scale = 1
    for _, c in rows:
        for row in c:
            d = getattr(row[-1], "denominator", 1)
            scale = scale * d // gcd(scale, d)
    if scale == 1:
        return rows, 1
    return [
        (la, [(*row[:-1], int(row[-1] * scale)) for row in c]) for la, c in rows
    ], scale


def _unscale(v: int, scale: int) -> Coefficient:
    """`v/scale`, exactly, and as an `int` when it divides."""
    return v if scale == 1 else Fraction(v, scale)


def _ht_rows(f: Param, what: str, expect: str) -> list[Any]:
    """The `(partition, numerator, denominator atoms)` rows the `H̃` entry
    points take, read off an element.

    Both directions use this encoding, so a value from either feeds back into
    the other.
    """
    if f.basis != expect:
        raise _needs(what, expect, f.basis, expect in BASES)
    rows = []
    for la, coeff in f:
        if not isinstance(coeff, QtRatio):
            raise ValueError(
                f"{what} needs coefficients in q and t, not "
                f"{type(coeff).__name__}"
            )
        rows.append((la, _qt_rows(coeff.numerator), list(coeff.denominator)))
    return rows


def _kind(f: Param) -> type | None:
    """The coefficient class an element's terms are written in, or `None` if it
    has no terms.

    Every term of an element shares one, because a family's entry point builds
    them all through the same wrapper.
    """
    for _, c in f:
        return type(c)
    return None


def _add(f: Param, g: Param) -> Param:
    """`f + g`, termwise in the basis both are written in.

    The fraction kinds go through the contract layer rather than adding here,
    because the coefficient classes compare **structurally**: an unreduced sum
    is the right value in a representation nothing else produces, and `==`
    against a value built another way would then be false. `Poly` and `QtPoly`
    need no such care — a polynomial sum is already canonical.
    """
    if f.basis != g.basis:
        raise ValueError(
            f"cannot add an element in {f.basis} to one in {g.basis}"
        )
    kind = _kind(f) or _kind(g)
    if kind is None or not len(f):
        return g if kind is not None or len(g) else f
    if not len(g):
        return f
    if kind is QtFrac:
        rows_f, sf = _mac_rows(f, "+", f.basis)
        rows_g, sg = _mac_rows(g, "+", g.basis)
        lcm = sf * sg // gcd(sf, sg)
        out = _c.macdonald_element_add(
            _mac_restale(rows_f, lcm // sf), _mac_restale(rows_g, lcm // sg)
        )
        return _mac_element(out, f.basis, lcm)
    if kind is AlphaFrac:
        jack_out = _c.jack_element_add(
            _jack_rows(f, "+", f.basis), _jack_rows(g, "+", g.basis)
        )
        return _jack_element(jack_out, f.basis)
    if kind is QtRatio:
        ht_out = _c.macdonald_ht_element_add(
            _ht_rows(f, "+", f.basis), _ht_rows(g, "+", g.basis)
        )
        return _ht_element(ht_out, f.basis)
    return Param(f.basis, _add_terms(f, g), f.parameters)


def _mac_restale(rows: list[Any], by: int) -> list[Any]:
    """Multiply every numerator in `(partition, numerator, denominator)` rows
    by an integer, so two elements scaled by different amounts can be added
    over one common scale.
    """
    if by == 1:
        return rows
    return [
        (la, [(a, b, c * by) for a, b, c in num], den) for la, num, den in rows
    ]


def _add_terms(f: Param, g: Param) -> dict[Partition, Any]:
    """The termwise sum for the two polynomial coefficient kinds, `Poly` and
    `QtPoly`.

    A polynomial sum is already canonical, so unlike the fraction kinds these
    need no reduction and no contract call.
    """
    out = dict(f.terms)
    for la, c in g:
        if la not in out:
            out[la] = c
            continue
        have = out[la]
        if not isinstance(have, (Poly, QtPoly)) or not isinstance(c, (Poly, QtPoly)):
            raise TypeError(
                f"cannot add a {type(c).__name__} coefficient to a "
                f"{type(have).__name__} one"
            )
        total = have + c
        if total:
            out[la] = total
        else:
            del out[la]
    return out


def _qt_rows(p: QtPoly) -> list[tuple[int, int, Coefficient]]:
    """A `QtPoly` as the `(a, b, coefficient)` rows its constructor takes."""
    return [(a, b, c) for (a, b), c in p.coefficients().items()]


def _scale(f: Param, c: Scalar) -> Param:
    """`c*f`, `c` a scalar in this element's parameters.

    The fraction kinds go through the contract layer for the reason `_add`
    gives, and with one more of their own: multiplying by `1 - q*t` when that
    factor sits in a denominator has to cancel it.
    """
    kind = _kind(f)
    if kind is None:
        return f
    if kind is QtFrac:
        rows, scale = _mac_rows(f, "*", f.basis)
        num, den, over = _qt_scalar(c)
        return _mac_element(
            _c.macdonald_element_scale(rows, num, den), f.basis, scale * over
        )
    if kind is AlphaFrac:
        num, atoms, over = _alpha_scalar(c)
        out = _c.jack_element_scale(
            _jack_rows(f, "*", f.basis), num, atoms, over
        )
        return _jack_element(out, f.basis)
    if kind is QtRatio:
        num, atoms, over = _ht_scalar(c)
        ht_out = _c.macdonald_ht_element_scale(
            _ht_rows(f, "*", f.basis), num, atoms
        )
        rows = ht_out if over == 1 else _ht_divide(ht_out, over)
        return _ht_element(rows, f.basis)
    terms: dict[Partition, Any] = {}
    w: Any
    for la, v in f:
        if isinstance(v, QtPoly):
            if not isinstance(c, (int, Fraction, QtPoly)):
                raise TypeError(
                    f"cannot scale an element in q and t by {type(c).__name__}"
                )
            w = v * c
        elif isinstance(v, Poly):
            if not isinstance(c, (int, Fraction)) and not (
                isinstance(c, Poly) and c.variable == v.variable
            ):
                raise TypeError(
                    f"cannot scale an element in {v.variable} by "
                    f"{type(c).__name__}"
                )
            w = v * c
        else:
            raise TypeError(f"cannot scale a {type(v).__name__} coefficient")
        if w:
            terms[la] = w
    return Param(f.basis, terms, f.parameters)


def _qt_scalar(c: Scalar) -> tuple[list[Any], list[Any], int]:
    """A scalar as `(numerator rows, denominator factors, divisor)` for
    `macdonald_element_scale`.

    The contract layer takes integer numerators, so a rational one is
    multiplied up here and the divisor hands it back in `_mac_element` — the
    same round trip `_mac_rows` makes, exact because scaling is linear.
    """
    if isinstance(c, (int, Fraction)):
        c = QtPoly({(0, 0): c})
    if isinstance(c, QtPoly):
        num, den = c.coefficients(), []
    elif isinstance(c, QtFrac):
        num, den = c.numerator.coefficients(), list(c.denominator)
    else:
        raise TypeError(f"cannot scale an element in q and t by {type(c).__name__}")
    over = 1
    for v in num.values():
        if isinstance(v, Fraction):
            over = over * v.denominator // gcd(over, v.denominator)
    return [(a, b, int(v * over)) for (a, b), v in num.items()], den, over


def _ht_scalar(c: Scalar) -> tuple[list[Any], list[Any], int]:
    """A scalar as `(numerator rows, denominator atoms, divisor)` for
    `macdonald_ht_element_scale`.

    The divisor is `_qt_scalar`'s: the contract layer takes integer
    numerators, so a rational one is multiplied up here and divided back out
    by `_ht_divide`.
    """
    if isinstance(c, (int, Fraction)):
        c = QtPoly({(0, 0): c})
    den: list[Any] = []
    if isinstance(c, QtPoly):
        num = c.coefficients()
    elif isinstance(c, QtRatio):
        num, den = c.numerator.coefficients(), list(c.denominator)
    else:
        raise TypeError(f"cannot scale an McdHt element by {type(c).__name__}")
    over = 1
    for v in num.values():
        if isinstance(v, Fraction):
            over = over * v.denominator // gcd(over, v.denominator)
    return [(a, b, int(v * over)) for (a, b), v in num.items()], den, over


def _ht_divide(rows: list[Any], by: int) -> list[Any]:
    """Divide every numerator in `(partition, numerator, atoms)` rows by an
    integer — the `restore` half of the round trip `_ht_scalar` begins.

    The atoms are untouched: they are the two binomial families, and an
    integer is neither.
    """
    return [
        (la, [(a, b, Fraction(v, by)) for a, b, v in num], den)
        for la, num, den in rows
    ]


def _alpha_scalar(c: Scalar) -> tuple[list[int], list[Any], int]:
    """A scalar as `(dense numerator, atoms, scale)` for `jack_element_scale`.

    Unlike `_qt_scalar` there is no divisor to hand back: an `AlphaFrac` row
    carries its own integer `scale`, so a rational coefficient goes there and
    crosses unchanged.
    """
    if isinstance(c, int):
        return [c], [], 1
    if isinstance(c, Fraction):
        return [c.numerator], [], c.denominator
    if isinstance(c, Poly) and c.variable == "alpha":
        num, atoms, scale = _dense(c), [], 1
    elif isinstance(c, AlphaFrac):
        num, atoms, scale = list(c.numerator), list(c.atoms), c.scale
    else:
        raise TypeError(f"cannot scale an element in alpha by {type(c).__name__}")
    lcm = 1
    for v in num:
        if isinstance(v, Fraction):
            lcm = lcm * v.denominator // gcd(lcm, v.denominator)
    return [int(v * lcm) for v in num], atoms, scale * lcm


def _needs(what: str, expect: str, got: str, convert: bool) -> ValueError:
    """The refusal every row builder raises when the basis is wrong.

    `convert` is for the classical bases alone: `.to()` is the way out of one
    of those and there is no way out of a parametric tag, so pointing at it
    from a `McdP` refusal would name a method that cannot help.
    """
    name = {"m": "a monomial-basis", "s": "a Schur-basis"}.get(expect, f"a {expect}")
    tail = f"; convert with .to({expect!r})" if convert else ""
    return ValueError(f"{what} needs {name} element, not {got}{tail}")


def _t_schur_rows(
    f: TSchurArg, what: str, expect: str = "s"
) -> tuple[list[Any], int]:
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
        if f.basis != expect:
            raise _needs(what, expect, f.basis, True)
        polys = [(la, {0: c}) for la, c in f]
    elif isinstance(f, Param):
        if f.basis != expect:
            raise _needs(what, expect, f.basis, expect in BASES)
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


def _mac_rows(
    f: MacArg, what: str, expect: str = "m"
) -> tuple[list[Any], int]:
    """The `(partition, numerator, denominator)` rows the Macdonald inverse
    expansions take, and the integer the numerators were scaled by.

    Accepts a `Sym` in the monomial basis, a `Param` in the monomial basis
    whose coefficients are in `q` and `t`, or the rows themselves. The
    contract layer takes integer numerators, so rational ones are multiplied
    up by their least common denominator here and divided back out in
    `_mac_element`; the expansion is linear, so the round trip is exact and
    the denominators are untouched by it.
    """
    if isinstance(f, Sym):
        if f.basis != expect:
            raise _needs(what, expect, f.basis, True)
        cells: list[Any] = [(la, {(0, 0): c}, ()) for la, c in f]
    elif isinstance(f, Param):
        if f.basis != expect:
            raise _needs(what, expect, f.basis, expect in BASES)
        cells = []
        for la, coeff in f:
            # The expansion is over ℚ(q,t); a coefficient in `t` alone, or in
            # α, has the same term structure and a different meaning, so it is
            # refused rather than read through whichever accessor exists.
            if isinstance(coeff, QtFrac):
                cells.append(
                    (la, coeff.numerator.coefficients(), coeff.denominator)
                )
            elif isinstance(coeff, QtPoly):
                cells.append((la, coeff.coefficients(), ()))
            else:
                raise ValueError(
                    f"{what} needs coefficients in q and t, not "
                    f"{type(coeff).__name__}"
                )
    else:
        return list(f), 1
    scale = 1
    for _, terms, _ in cells:
        for c in terms.values():
            if isinstance(c, Fraction):
                d = c.denominator
                scale = scale * d // gcd(scale, d)
    return [
        (
            la,
            [(a, b, int(c * scale)) for (a, b), c in terms.items()],
            list(den),
        )
        for la, terms, den in cells
    ], scale


def _jack_rows(f: JackArg, what: str, expect: str = "m") -> list[Any]:
    """The `(partition, numerator, atoms, scale)` rows the Jack inverse
    expansions take.

    Accepts a `Sym` in the monomial basis, a `Param` in the monomial basis
    whose coefficients are in α, or the rows themselves. A rational numerator
    coefficient needs no round trip the way `_mac_rows` does: every row
    already carries its own integer `scale`, so the row's least common
    denominator goes there and the value crosses unchanged.
    """
    if isinstance(f, Sym):
        if f.basis != expect:
            raise _needs(what, expect, f.basis, True)
        cells: list[Any] = [(la, [c], (), 1) for la, c in f]
    elif isinstance(f, Param):
        if f.basis != expect:
            raise _needs(what, expect, f.basis, expect in BASES)
        cells = []
        for la, coeff in f:
            # The expansion is over ℚ(α); a coefficient in `t` or in `q` and
            # `t` has the same term structure and a different meaning, so it
            # is refused rather than read through whichever accessor exists.
            if isinstance(coeff, AlphaFrac):
                cells.append((la, coeff.numerator, coeff.atoms, coeff.scale))
            elif isinstance(coeff, Poly) and coeff.variable == "alpha":
                cells.append((la, _dense(coeff), (), 1))
            else:
                raise ValueError(
                    f"{what} needs coefficients in alpha, not "
                    f"{type(coeff).__name__}"
                )
    else:
        return list(f)
    rows = []
    for la, num, atoms, scale in cells:
        lcm = 1
        for c in num:
            if isinstance(c, Fraction):
                lcm = lcm * c.denominator // gcd(lcm, c.denominator)
        rows.append(
            (la, [int(c * lcm) for c in num], list(atoms), scale * lcm)
        )
    return rows


def _dense(p: Poly) -> list[Coefficient]:
    """A `Poly` as a dense coefficient list, index the exponent."""
    terms = p.coefficients()
    if not terms:
        return []
    return [terms.get(k, 0) for k in range(max(terms) + 1)]


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
