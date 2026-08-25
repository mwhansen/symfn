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
from typing import Any, Union, cast

from . import symfn as _c
from ._bases import BASES, BaseRingError, BasisError
from ._param import (
    AlphaFrac,
    ParamCoefficient,
    Poly,
    QtFrac,
    QtPoly,
    QtRatio,
    _as_one_variable,
)
from ._sym import Sym, _partition
from ._types import Coefficient, Partition, PartitionArg

#: What `nabla` accepts: an element in the Schur basis, or the contract layer's
#: rows themselves.
NablaArg = Union["Sym", Iterable[Any]]

#: What the contract layer's `(q, t)`-graded rows look like on arrival.
QtRows = Iterable[tuple[Partition, Iterable[tuple[int, int, Coefficient]]]]

#: What a `Sym` may be multiplied by: an integer, a rational, a polynomial in
#: its own parameters, or a coefficient of the kind it carries.
Scalar = Union[int, Fraction, "Poly", "QtPoly", "QtFrac", "QtRatio", "AlphaFrac"]

#: What the Hall-Littlewood inverse expansions accept: a Schur-basis element,
#: with or without a parameter in its coefficients, or the contract layer's
#: `t`-rows.
TSchurArg = Union["Sym", Iterable[Any]]

#: What the Macdonald inverse expansions take: an element in the monomial
#: basis, or the contract layer's `(partition, numerator, denominator)` rows.
MacArg = Union["Sym", Iterable[Any]]

#: What the Jack inverse expansions take: an element in the monomial basis, or
#: the contract layer's `(partition, numerator, atoms, scale, tail)` rows.
JackArg = Union["Sym", Iterable[Any]]

__all__ = ["macdonald", "jack", "hl", "llt"]


def _qt_element(rows: QtRows, basis: str) -> Sym:
    """Wrap `(partition, [(a, b, coefficient)])` rows as a `Sym` in q, t."""
    return Sym(basis, [(la, QtPoly(c)) for la, c in rows], ("q", "t"))


def _q_element(rows: QtRows, basis: str) -> Sym:
    """Wrap `(partition, [(a, b, coefficient)])` rows as a `Sym` in q alone.

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
    return Sym(basis, terms, ("q",))


def _ht_element(rows: Iterable[Any], basis: str = "McdHt") -> Sym:
    """Wrap `(partition, numerator, denominator atoms)` rows as a `Sym` in
    `q` and `t`, tagged as Sage prints it.

    The expansion out of `H̃` lands in the Schur basis with the same
    coefficients, so this builds both ends of the pair.
    """
    return Sym(basis, [(la, QtRatio(n, d)) for la, n, d in rows], ("q", "t"))


def _mac_element(rows: Iterable[Any], basis: str = "m", scale: int = 1) -> Sym:
    """Wrap `(partition, numerator, denominator)` rows as a `Sym` in q, t,
    dividing every numerator coefficient by `scale` — the `restore` half of
    the denominator round trip `_mac_rows` begins.
    """
    if scale == 1:
        return Sym(basis, [(la, QtFrac(n, d)) for la, n, d in rows], ("q", "t"))
    return Sym(
        basis,
        [
            (la, QtFrac([(a, b, Fraction(v, scale)) for a, b, v in n], d))
            for la, n, d in rows
        ],
        ("q", "t"),
    )


def _jack_element(rows: Iterable[Any], basis: str = "m") -> Sym:
    """Wrap `(partition, numerator, atoms, scale, tail)` rows as a Jack
    `Sym`.

    The tail is the general denominator factor only a plethysm produces; every
    other row carries an empty one.
    """
    return Sym(
        basis,
        [(la, AlphaFrac(n, d, k, t)) for la, n, d, k, t in rows],
        ("alpha",),
    )


def _unit(basis: str, la: PartitionArg) -> Sym:
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
        return _jack_element([(key, [1], [], 1, [])], basis)
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

    def P(self, la: PartitionArg) -> Sym:
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

    def Q(self, la: PartitionArg) -> Sym:
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

    def J(self, la: PartitionArg) -> Sym:
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

    def Htilde(self, mu: PartitionArg) -> Sym:
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

    def to_Htilde(self, f: NablaArg) -> Sym:
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

        Accepts what `nabla` accepts — a `Sym` in the Schur basis,
        parameter-free or with coefficients in `q` and `t`, or the contract
        rows — and unlike `nabla` it takes mixed degrees, expanding each
        degree on its own.

        # Raises

        Raises `ValueError` unless the argument is in the Schur basis with
        integer coefficients, or coefficients in `q` and `t`.
        """
        return _ht_element(_c.schur_to_macdonald_ht(_schur_rows(f, "to_Htilde")))

    def to_J(self, f: NablaArg) -> Sym:
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

    def to_P(self, f: MacArg) -> Sym:
        """`f`, given in the monomial basis, rewritten in the `P` basis, as a
        `Sym` tagged `McdP`.

            >>> from symfn import macdonald, m
            >>> macdonald.to_P(m([2])).coefficient([1, 1])
            (-1 + t - q + q*t)/(1 - q*t)
            >>> macdonald.to_P(macdonald.P([2, 1]).to("m"))
            McdP[2,1]

        `m_2 = P_2 − [(1−t)(1+q)/(1−q·t)] P_11`: the coefficient `P → m` puts
        on the dominance-smaller shape, negated. The `q ↔ t` swap gives
        `(1−q)(1+t)/(1−q·t)` instead, which is the twist to check.

        `f` may be a `Sym` in the monomial basis — with parameters set, so
        a `P`, `Q` or `J` value feeds back in, as the second example does,
        or parameter-free — or the
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

    def to_Q(self, f: MacArg) -> Sym:
        """`f`, given in the monomial basis, rewritten in the `Q` basis, as a
        `Sym` tagged `McdQ`.

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

    def nabla_e(self, n: int) -> Sym:
        """`∇ e_n`, in the Schur basis — the Shuffle Theorem's left side.

            >>> from symfn import macdonald
            >>> macdonald.nabla_e(2)
            (t + q)*s[1,1] + s[2]

        # Raises

        Raises `ValueError` unless `n` is positive.
        """
        return _qt_element(_c.nabla_e(n), "s")

    def nabla(self, f: NablaArg) -> Sym:
        """`∇` applied to a `Sym` in the Schur basis, or to contract rows.

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

    def nabla_power(self, f: NablaArg, r: int) -> Sym:
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

    def delta_prime_e(self, k: int, n: int) -> Sym:
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

    def delta_ek(self, k: int, f: NablaArg) -> Sym:
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

    def delta_prime_ek(self, k: int, f: NablaArg) -> Sym:
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

    def theta_ek(self, k: int, f: NablaArg) -> Sym:
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

    def big_pi(self, f: NablaArg) -> Sym:
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

    def P(self, la: PartitionArg) -> Sym:
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

    def Q(self, la: PartitionArg) -> Sym:
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

    def J(self, la: PartitionArg) -> Sym:
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

    def to_P(self, f: JackArg) -> Sym:
        """`f`, given in the monomial basis, rewritten in the `P` basis, as a
        `Sym` tagged `JackP`.

            >>> from symfn import jack, m
            >>> jack.to_P(m([2])).coefficient([1, 1])
            -2/(alpha + 1)
            >>> jack.to_P(jack.P([2, 1]).to("m"))
            JackP[2,1]

        `m_2 = P_2 − [2/(α+1)] P_11`: the coefficient `P → m` puts on the
        dominance-smaller shape, negated. Sending `α → 1/α` would give
        `−2α/(α+1)` instead, which is the twist to check; the two agree at
        `α = 1`, so that specialization cannot see it.

        `f` may be a `Sym` in the monomial basis — with parameters set, so
        a `P`, `Q` or `J` value feeds back in, as the second example does,
        or parameter-free — or the
        contract layer's rows. An element in another classical basis is
        refused rather than converted, on the same grounds as `BasisError`:
        write `jack.to_P(f.to("m"))` and the conversion is the caller's, with
        its cost visible.

        # Raises

        Raises `ValueError` unless `f` is in the monomial basis with
        coefficients in α, and unless every support is a partition.
        """
        return _jack_element(_c.monomial_to_jack_p(_jack_rows(f, "to_P")), "JackP")

    def to_Q(self, f: JackArg) -> Sym:
        """`f`, given in the monomial basis, rewritten in the `Q` basis, as a
        `Sym` tagged `JackQ`.

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

    def to_J(self, f: JackArg) -> Sym:
        """`f`, given in the monomial basis, rewritten in the `J` basis, as a
        `Sym` tagged `JackJ`.

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

        The parameter is already substituted, so this returns a
        parameter-free element; `jack.P(la).at(alpha=2)` is the same element.

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
        num, atoms, scale, tail = _c.jack_structure_constant(
            _partition(la), _partition(mu), _partition(nu)
        )
        return AlphaFrac(num, atoms, scale, tail)

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

    def Qp(self, la: PartitionArg) -> Sym:
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

    def P(self, la: PartitionArg) -> Sym:
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

    def to_P(self, f: TSchurArg) -> Sym:
        """`f`, given in the Schur basis, rewritten in the `P` basis, as a
        `Sym` tagged `HLP`.

            >>> from symfn import hl, s
            >>> hl.to_P(s([2]))
            t*HLP[1,1] + HLP[2]
            >>> hl.to_P(hl.P([3, 1]).to("s"))
            HLP[3,1]

        `s_2 = P_2 + t·P_11`: the coefficient of `P_λ` in `s_μ` is the
        Kostka-Foulkes polynomial `K_{μλ}(t)`, so the `t` sits on the
        dominance-smaller shape. Sage's `HLP(s[2])` prints the same value.
        `f` may be a parameter-free `Sym`, a `Sym` in `t` (so a `P` or `Qp`
        value feeds back in, as the second example does), or the contract
        layer's rows;
        rational coefficients are scaled through the boundary and restored.

        # Raises

        Raises `ValueError` unless `f` is in the Schur basis with coefficients
        in `t` alone, and unless every support is a partition.
        """
        rows, scale = _t_schur_rows(f, "to_P")
        return _t_element(_c.schur_to_hall_littlewood_p(rows), "HLP", scale)

    def to_Qp(self, f: TSchurArg) -> Sym:
        """`f`, given in the Schur basis, rewritten in the `Q'` basis, as a
        `Sym` tagged `HLQp`.

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

    One family, two presentations. `G` takes a tuple of partitions, with an
    optional content offset per shape, and sums `q^{inv(T)} x^T` over
    semistandard fillings, `inv` counting attacking pairs that are out of
    order. `Gtilde`, `Htilde` and `H` take a partition
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

    def Gtilde(self, la: PartitionArg, k: int) -> Sym:
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

    def H(self, mu: PartitionArg, k: int) -> Sym:
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
    ) -> Sym:
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

    def schur(self, la: PartitionArg, k: int) -> Sym:
        """`G̃^(k)_λ` in the **Schur** basis, where `Gtilde` gives monomial.

            >>> from symfn import llt
            >>> llt.schur([1, 1], 2).basis
            's'

        # Raises

        Raises `ValueError` unless λ is a partition and `k ≥ 1`.
        """
        return _q_element(_c.llt_schur(_partition(la), k), "s")

    def Htilde(self, mu: PartitionArg, k: int) -> Sym:
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
) -> Sym:
    """Wrap `(partition, [(t_exponent, coefficient)])` rows as a `Sym`,
    dividing every coefficient by `scale` — the `restore` half of the
    denominator round trip `_t_schur_rows` begins.
    """
    if scale == 1:
        return Sym(basis, [(la, Poly("t", c)) for la, c in rows], ("t",))
    return Sym(
        basis,
        [
            (la, Poly("t", [(k, Fraction(v, scale)) for k, v in c]))
            for la, c in rows
        ],
        ("t",),
    )


#: Where each parametric basis expands: the classical basis its family's
#: forward direction is defined in. `Sym.to` reads this, and `Sym.at`
#: goes through it, since evaluation is defined on the classical bases.
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


def _expand(f: Sym) -> Sym:
    """A `Sym` in a parametric basis, expanded in the basis `EXPANDS_IN`
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


def _qt_pack(f: Sym) -> tuple[list[Any], int]:
    """An element's coefficients as integer `(a, b, value)` exponent rows, and
    the integer they were scaled by.

    The polynomial encoding is what the `(q,t)` entry points read, and a
    one-variable `Poly` uses the second slot for whatever its variable is —
    the entry points never ask what the exponents count.

    # Raises

    Raises `ValueError` for a coefficient class that is not a polynomial.
    """
    kind = _kind(f)
    # `_kind` naming a class is the promise that every coefficient is one.
    coeffs = cast("Iterable[tuple[Any, ParamCoefficient]]", f)
    if kind is Poly:
        return _qt_int_rows([(la, _one_variable_row(c)) for la, c in coeffs])
    if kind is QtPoly:
        return _qt_int_rows([(la, _two_variable_row(c)) for la, c in coeffs])
    raise ValueError(
        f"this route reads the polynomial encoding, which "
        f"{kind.__name__ if kind else 'an empty element'} is not"
    )


def _qt_unpack(rows: list[Any], like: Sym, dst: str, scale: int) -> Sym:
    """Exponent rows back as a `Sym` in `dst`, in the coefficient class and
    parameters `like` carries, dividing out the scale `_qt_pack` cleared.
    """
    if _kind(like) is Poly:
        var = next(c.variable for _, c in like if isinstance(c, Poly))
        return Sym(
            dst,
            [
                (la, Poly(var, [(b, _unscale(v, scale)) for _, b, v in c]))
                for la, c in rows
            ],
            like.parameters,
        )
    return Sym(
        dst,
        [
            (la, QtPoly([(a, b, _unscale(v, scale)) for a, b, v in c]))
            for la, c in rows
        ],
        like.parameters,
    )


def _carry_qt(f: Sym, dst: str, call: Callable[[list[Any]], list[Any]]) -> Sym:
    """An element's exponent rows through one contract call, rebuilt in `dst`
    in the coefficient class they went in as.

    The shared half of every operation this layer sends through the polynomial
    encoding: pack, clear denominators, call, restore, rebuild. Nothing here is
    arithmetic on the element — `call` is the whole computation.
    """
    rows, scale = _qt_pack(f)
    return _qt_unpack(call(rows), f, dst, scale)


def _rat(n: int, d: int, scale: int = 1) -> Coefficient:
    """`n/(d·scale)`, exactly, and as an `int` when it divides."""
    v = Fraction(n, d * scale)
    return int(v) if v.denominator == 1 else v


def _to_power(f: Sym) -> Sym:
    """A `Sym` with parameters, rewritten in the power-sum basis over its own
    coefficient ring.

    The conversion divides by z_μ — the one classical target that does — so
    it crosses through the `to_power_*` entry points, which return the
    division: a Jack row carries it in its own `scale` slot, and the `(q,t)`
    encodings hand each numerator coefficient back as an explicit
    `(numerator, denominator)` pair. A parametric basis expands first, as it
    does for every classical target.
    """
    g = f if f.basis in BASES else _expand(f)
    if not len(g):
        return Sym("p", [], g.parameters)
    if g.basis == "p":
        return g
    kind = _kind(g)
    if kind is AlphaFrac or (
        kind is Poly
        and next(c.variable for _, c in g if isinstance(c, Poly)) == "alpha"
    ):
        rows = _jack_rows(g, "to", g.basis)
        return _jack_element(_c.to_power_jack(rows, g.basis), "p")
    if kind in (Poly, QtPoly):
        rows, scale = _qt_pack(g)
        out = _c.to_power_qt(rows, g.basis)
        if kind is Poly:
            var = next(c.variable for _, c in g if isinstance(c, Poly))
            return Sym(
                "p",
                [
                    (la, Poly(var, [(b, _rat(n, d, scale)) for _, b, n, d in c]))
                    for la, c in out
                ],
                g.parameters,
            )
        return Sym(
            "p",
            [
                (la, QtPoly([(a, b, _rat(n, d, scale)) for a, b, n, d in c]))
                for la, c in out
            ],
            g.parameters,
        )
    if kind is QtFrac:
        rows, scale = _mac_rows(g, "to", g.basis)
        mac_out = _c.to_power_macdonald(rows, g.basis)
        return Sym(
            "p",
            [
                (
                    la,
                    QtFrac(
                        [(a, b, _rat(n, d, scale)) for a, b, n, d in num], den
                    ),
                )
                for la, num, den in mac_out
            ],
            ("q", "t"),
        )
    ht_out = _c.to_power_ht(_ht_rows(g, "to", g.basis), g.basis)
    return Sym(
        "p",
        [
            (la, QtRatio([(a, b, _rat(n, d)) for a, b, n, d in num], den))
            for la, num, den in ht_out
        ],
        ("q", "t"),
    )


def _convert(f: Sym, dst: str) -> Sym:
    """A `Sym` in a classical basis, rewritten in another classical basis.

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
        return Sym(dst, dict(f), f.parameters)
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


#: The entry point for each Hopf operation over each coefficient ring. Both
#: operations are defined on the Schur basis — ω conjugates the index, the
#: antipode conjugates and signs — and reach any other basis by a change of
#: basis on each side. Each ring has its own encoding, so the pair is restated
#: once per ring exactly as the converters are.
_HOPF: dict[str, dict[type, Callable[[list[Any]], list[Any]]]] = {
    "omega": {
        Poly: _c.omega_qt_terms,
        QtPoly: _c.omega_qt_terms,
        QtFrac: _c.omega_macdonald_terms,
        AlphaFrac: _c.omega_jack_terms,
        QtRatio: _c.omega_ht_terms,
    },
    "antipode": {
        Poly: _c.antipode_qt_terms,
        QtPoly: _c.antipode_qt_terms,
        QtFrac: _c.antipode_macdonald_terms,
        AlphaFrac: _c.antipode_jack_terms,
        QtRatio: _c.antipode_ht_terms,
    },
}


def _hopf(f: Sym, op: str) -> Sym:
    """`omega` or `antipode`, returned in the basis it was given in.

    Three legs where a classical basis needs one: a parametric basis expands
    into its pivot, the pivot converts to Schur, and the result travels back
    the same way — so an element in a family's own basis comes back in it, as
    `Sym.omega` comes back in the basis it was handed.

    Which entry point runs is decided by the coefficient class, for the reason
    `_convert` gives: the act itself is a relabeling of the index, and only the
    encoding the coefficients cross in differs between the rings.

    # Raises

    Raises `ValueError` if the element carries no coefficient class this layer
    knows, and for a parametric basis whose inverse expansion runs over one it
    does not.
    """
    tag = f.basis
    if not len(f):
        return Sym(tag, {}, f.parameters)
    schur = _convert(f if tag in BASES else _expand(f), "s")
    kind = _kind(schur)
    call = _HOPF[op].get(kind)  # type: ignore[arg-type]
    if call is None:
        raise ValueError(
            f"{op} is not written for "
            f"{kind.__name__ if kind else 'these'} coefficients"
        )
    acted: Sym
    if kind in (Poly, QtPoly):
        acted = _carry_qt(schur, "s", call)
    elif kind is QtFrac:
        mac_rows, mac_scale = _mac_rows(schur, op, "s")
        acted = _mac_element(call(mac_rows), "s", mac_scale)
    elif kind is AlphaFrac:
        acted = _jack_element(call(_jack_rows(schur, op, "s")), "s")
    else:
        acted = _ht_element(call(_ht_rows(schur, op, "s")), "s")
    return _back_to(acted, tag, op)


#: The skew entry point for each coefficient ring, on the model of `_HOPF`.
_SKEW: dict[type, Callable[..., list[Any]]] = {
    Poly: _c.skew_by_qt,
    QtPoly: _c.skew_by_qt,
    QtFrac: _c.skew_by_macdonald,
    AlphaFrac: _c.skew_by_jack,
    QtRatio: _c.skew_by_ht,
}


def _skew(f: Sym, g: Sym) -> Sym:
    """`g^⊥ f`, the adjoint of multiplication by `g`, in `f`'s own basis.

    The same three legs `_hopf` takes, and for the same reason: the rule is
    written on the Schur basis, so a parametric basis expands into its pivot,
    converts, acts, and travels back.

    `g` keeps its own basis, because that basis selects which rule runs and not
    merely how `g` is read — `h`, `e` and `p` take the Pieri, dual-Pieri and
    Murnaghan-Nakayama paths. A parametric `g` has no such rule and is expanded
    into a classical basis first. A `g` without parameters is lifted into `f`'s
    ring rather than refused, which is what Sage does.

    # Raises

    Raises `ValueError` if the element carries no coefficient class this layer
    knows, and `BaseRingError` if `g` is over a different base ring.
    """
    tag = f.basis
    if not len(f):
        return Sym(tag, {}, f.parameters)
    schur = _convert(f if tag in BASES else _expand(f), "s")
    kind = _kind(schur)
    call = _SKEW.get(kind)  # type: ignore[arg-type]
    if call is None:
        raise ValueError(
            f"skew_by is not written for "
            f"{kind.__name__ if kind else 'these'} coefficients"
        )
    other = _lift_to(g, schur, "skew_by")
    if other.basis not in BASES:
        other = _convert(_expand(other), "s")
    code = other.basis
    if kind in (Poly, QtPoly):
        rows_f, scale_f = _qt_pack(schur)
        rows_g, scale_g = _qt_pack(other)
        acted = _qt_unpack(
            call(rows_f, rows_g, code), schur, "s", scale_f * scale_g
        )
    elif kind is QtFrac:
        mac_f, mac_sf = _mac_rows(schur, "skew_by", "s")
        mac_g, mac_sg = _mac_rows(other, "skew_by", code)
        acted = _mac_element(
            call(mac_f, mac_g, code), "s", mac_sf * mac_sg
        )
    elif kind is AlphaFrac:
        acted = _jack_element(
            call(
                _jack_rows(schur, "skew_by", "s"),
                _jack_rows(other, "skew_by", code),
                code,
            ),
            "s",
        )
    else:
        acted = _ht_element(
            call(
                _ht_rows(schur, "skew_by", "s"),
                _ht_rows(other, "skew_by", code),
                code,
            ),
            "s",
        )
    return _back_to(acted, tag, "skew_by")


def _ring_rows(
    f: Sym, table: dict[type, Callable[..., Any]], what: str, expect: str
) -> tuple[Callable[..., Any], list[Any], int]:
    """The entry point for `f`'s coefficient ring, `f`'s rows in that ring's
    encoding, and the integer the numerators were scaled by.

    The shared first half of every operation this layer sends over a ring:
    which of the four entry points runs is decided by the coefficient class,
    because each ring has its own encoding, and the routing is otherwise the
    same.

    # Raises

    Raises `ValueError` if the element carries no coefficient class this layer
    knows.
    """
    kind = _kind(f)
    call = table.get(kind)  # type: ignore[arg-type]
    if call is None:
        raise ValueError(
            f"{what} is not written for "
            f"{kind.__name__ if kind else 'these'} coefficients"
        )
    if kind in (Poly, QtPoly):
        rows, scale = _qt_pack(f)
    elif kind is QtFrac:
        rows, scale = _mac_rows(f, what, expect)
    elif kind is AlphaFrac:
        rows, scale = _jack_rows(f, what, expect), 1
    else:
        rows, scale = _ht_rows(f, what, expect), 1
    return call, rows, scale


def _cell_coeff(
    cell: list[Any] | tuple[Any, ...], like: Sym, scale: int
) -> ParamCoefficient:
    """One coefficient in `like`'s class, from the cell its ring encodes it as,
    dividing out the scale `_ring_rows` cleared.

    The inverse of the encoding half of `_ring_rows`, for the operations whose
    answer is a coefficient rather than an element.
    """
    kind = _kind(like)
    if kind in (Poly, QtPoly):
        return _qt_coeff(list(cell), like, scale)
    if kind is QtFrac:
        num, den = cell
        return QtFrac([(a, b, _unscale(v, scale)) for a, b, v in num], den)
    if kind is AlphaFrac:
        num, den, over, tail = cell
        return AlphaFrac(num, den, over, tail)
    num, den = cell
    return QtRatio(num, den)


#: The expansion entry point for each coefficient ring.
_EXPAND: dict[type, Callable[..., Any]] = {
    Poly: _c.expand_qt,
    QtPoly: _c.expand_qt,
    QtFrac: _c.expand_macdonald,
    AlphaFrac: _c.expand_jack,
    QtRatio: _c.expand_ht,
}

#: The evaluation entry point for each coefficient ring.
_EVALUATE: dict[type, Callable[..., Any]] = {
    Poly: _c.evaluate_qt,
    QtPoly: _c.evaluate_qt,
    QtFrac: _c.evaluate_macdonald,
    AlphaFrac: _c.evaluate_jack,
    QtRatio: _c.evaluate_ht,
}


def _expand_alphabet(f: Sym, n: int) -> dict[tuple[int, ...], Any]:
    """`f` laid out over `n` variables, as a `{exponent vector: coefficient}`
    mapping.

    Every basis reaches the layout through the monomial basis, where the
    expansion is definitional, so the conversion happens here and the entry
    points take only that basis. Laying the exponents out copies coefficients
    and does no arithmetic, so nothing about the ring is exercised past the
    conversion.
    """
    if not len(f):
        return {}
    mono = _convert(f if f.basis in BASES else _expand(f), "m")
    call, rows, scale = _ring_rows(mono, _EXPAND, "an expansion", "m")
    return {
        v: _cell_coeff(c[0] if len(c) == 1 else tuple(c), mono, scale)
        for v, *c in call(rows, n)
    }


def _evaluate(f: Sym, xs: Sequence[Coefficient]) -> ParamCoefficient:
    """`f` at the integer alphabet `xs`, as one coefficient.

    The alphabet injects into the coefficient ring, so this answers over the
    parameters what `Sym.evaluate` answers over ℚ, with them carried.
    """
    schur = _convert(f if f.basis in BASES else _expand(f), "s")
    if not len(schur):
        return 0  # type: ignore[return-value]
    call, rows, scale = _ring_rows(schur, _EVALUATE, "an evaluation", "s")
    return _cell_coeff(call(rows, list(xs)), schur, scale)


#: The dimension entry point for each coefficient ring.
_DIMENSION: dict[type, Callable[..., Any]] = {
    Poly: _c.dimension_qt,
    QtPoly: _c.dimension_qt,
    QtFrac: _c.dimension_macdonald,
    AlphaFrac: _c.dimension_jack,
    QtRatio: _c.dimension_ht,
}

#: The principal-specialization entry point for each coefficient ring.
_PRINCIPAL: dict[type, Callable[..., Any]] = {
    Poly: _c.principal_specialization_qt,
    QtPoly: _c.principal_specialization_qt,
    QtFrac: _c.principal_specialization_macdonald,
    AlphaFrac: _c.principal_specialization_jack,
    QtRatio: _c.principal_specialization_ht,
}


def _functional(
    f: Sym, table: dict[type, Callable[..., Any]], what: str, *args: int
) -> ParamCoefficient:
    """`Σ_λ c_λ w(λ)` for a weight the shape alone decides, as one coefficient.

    The dimension and the principal specialization differ only in that weight,
    so they differ only in `table` here. Both are computed on the Schur
    expansion, where the weight is defined.
    """
    schur = _convert(f if f.basis in BASES else _expand(f), "s")
    if not len(schur):
        return 0  # type: ignore[return-value]
    call, rows, scale = _ring_rows(schur, table, what, "s")
    return _cell_coeff(call(rows, *args), schur, scale)


def _principal_q(f: Sym, n: int) -> QtPoly:
    """The value at `1, q, …, q^{n−1}`, as a `(q,t)`-polynomial.

    The specialization introduces `q`, and this layer's coefficient classes
    carry at most two variables, so it is written only where the base ring is
    a single variable other than `q` — Hall-Littlewood and LLT, over `ℚ[t]`.
    Anything else is refused for want of a place to put `q`, which is the wall
    Sage reports as "the variable q is in the base ring, pass it explicitly" —
    and `_principal_at` is the way past it, in every ring.

    # Raises

    Raises `ValueError` if the base ring already carries `q`, and for the
    coefficient classes that have no free variable at all.
    """
    schur = _convert(f if f.basis in BASES else _expand(f), "s")
    if not len(schur):
        return QtPoly([])
    kind = _kind(schur)
    if "q" in f.parameters:
        raise ValueError(
            f"the variable q is in the base ring {_ring(f.parameters)}, and "
            "this specialization introduces it; pass the value explicitly, as "
            "principal_specialization(n, q=...)"
        )
    if kind is not Poly:
        raise ValueError(
            "this specialization introduces q, and the "
            f"{kind.__name__ if kind else 'these'} coefficients here have no "
            "free variable for it; pass the value explicitly, as "
            "principal_specialization(n, q=...)"
        )
    rows, scale = _qt_pack(schur)
    out = _c.principal_specialization_q_qt(rows, n)
    return QtPoly([(a, b, _unscale(v, scale)) for a, b, v in out])


#: The specialization at a base-ring alphabet, for each coefficient ring.
_PRINCIPAL_AT: dict[type, Callable[..., Any]] = {
    Poly: _c.principal_specialization_at_qt,
    QtPoly: _c.principal_specialization_at_qt,
    QtFrac: _c.principal_specialization_at_macdonald,
    AlphaFrac: _c.principal_specialization_at_jack,
    QtRatio: _c.principal_specialization_at_ht,
}

#: The zero coefficient in each ring's encoding, for a `_coeff_cell` whose
#: term dropped out.
_ZERO_CELL: dict[type, tuple[Any, ...] | list[Any]] = {
    Poly: [],
    QtPoly: [],
    QtFrac: ([], []),
    AlphaFrac: ([], [], 1, []),
    QtRatio: ([], []),
}


def _coeff_cell(
    c: ParamCoefficient | int | Fraction, like: Sym, what: str
) -> tuple[Any, ...] | list[Any]:
    """One coefficient in its ring's encoding: the inverse of `_cell_coeff`.

    Built by handing the coefficient through the same path a row of an element
    takes — a one-term element at the empty partition — so a value the ring
    cannot hold raises here rather than crossing malformed.

    # Raises

    Raises `ValueError` if the coefficient is not integral in the encoding,
    since `_ring_rows` clears such a denominator by scaling and nothing
    downstream can restore it: the alphabet enters at every power from 0 to
    the degree, not linearly.
    """
    kind = _kind(like)
    term = (
        Sym(like.basis, [((), c)], like.parameters)
        if isinstance(c, QtRatio)
        else _constant(like, c)
    )
    if not len(term):
        return _ZERO_CELL[kind]  # type: ignore[index]
    if _kind(term) is not kind:
        raise ValueError(
            f"{what}: the alphabet is over a different coefficient ring than "
            "the element"
        )
    _, rows, scale = _ring_rows(term, _PRINCIPAL_AT, what, like.basis)
    if scale != 1:
        raise ValueError(
            f"{what}: the alphabet {c} is not integral in this encoding, and "
            "it enters at every power rather than linearly, so a cleared "
            "denominator cannot be restored"
        )
    cell = tuple(rows[0][1:])
    return cell[0] if kind in (Poly, QtPoly) else cell


def _principal_at(
    f: Sym, n: int, q: ParamCoefficient | int | Fraction
) -> ParamCoefficient:
    """The value at `1, q, …, q^{n−1}` with `q` a coefficient of `f`'s own
    ring, rather than a variable the ring must have room for.

    This is what `_principal_q` refuses to guess at, and it is how Sage spells
    the same question: `P[2].principal_specialization(3, q=q)` substitutes an
    element of the base ring. Every ring here has it, including ℚ(α), which
    has no free variable at all.

    # Raises

    Raises `ValueError` if `f` carries no coefficient class this layer knows,
    and if the alphabet is not integral in that class's encoding.
    """
    schur = _convert(f if f.basis in BASES else _expand(f), "s")
    if not len(schur):
        return 0  # type: ignore[return-value]
    cell = _coeff_cell(q, schur, "the value at 1, q, ...")
    call, rows, scale = _ring_rows(
        schur, _PRINCIPAL_AT, "the value at 1, q, ...", "s"
    )
    return _cell_coeff(call(rows, n, cell), schur, scale)


#: The internal (Kronecker) product entry point for each coefficient ring.
_INTERNAL: dict[type, Callable[..., Any]] = {
    Poly: _c.internal_product_qt,
    QtPoly: _c.internal_product_qt,
    QtFrac: _c.internal_product_macdonald,
    AlphaFrac: _c.internal_product_jack,
    QtRatio: _c.internal_product_ht,
}


def _internal(f: Sym, g: Sym) -> Sym:
    """`f * g` under the internal (Kronecker) product, in `f`'s own basis.

    The same three legs the other Schur-basis operations take. The structure
    constants are Kronecker coefficients, which are integers carrying no
    parameter, so the coefficient ring is multiplied through.

    Unlike the Hall pairing this combines two elements of the ring, so the
    two must be in the same basis, on the same grounds `*` refuses. A
    parameter-free element is lifted into `f`'s ring rather than refused.

    # Raises

    Raises `BasisError` unless both are in the same basis, `BaseRingError`
    unless both are over the same base ring, and `ValueError` if the element
    carries no coefficient class this layer knows.
    """
    if f.basis != g.basis:
        raise BasisError(
            f"cannot combine {f.basis} with {g.basis}; convert one with .to()"
        )
    if g.parameters:
        _same_ring(f, g)
    tag = f.basis
    if not len(f) or not len(g):
        return Sym(tag, {}, f.parameters)
    schur = _convert(f if tag in BASES else _expand(f), "s")
    other = _lift_to(g, schur, "the internal product")
    other = _convert(other if other.basis in BASES else _expand(other), "s")
    call, rows_f, scale_f = _ring_rows(
        schur, _INTERNAL, "the internal product", "s"
    )
    _, rows_g, scale_g = _ring_rows(other, _INTERNAL, "the internal product", "s")
    raw = call(rows_f, rows_g)
    scale = scale_f * scale_g
    kind = _kind(schur)
    acted: Sym
    if kind in (Poly, QtPoly):
        acted = _qt_unpack(raw, schur, "s", scale)
    elif kind is QtFrac:
        acted = _mac_element(raw, "s", scale)
    elif kind is AlphaFrac:
        acted = _jack_element(raw, "s")
    else:
        acted = _ht_element(raw, "s")
    return _back_to(acted, tag, "the internal product")


#: The plethysm entry point for each coefficient ring.
#:
#: Jack's is the one that can produce a coefficient with a `tail`: `p_n` raises
#: the variable, so over ℚ(α) it is α ↦ α^n, and an atom `α + 1` becomes
#: `α² + 1`, which is irreducible and so belongs to no linear factorization.
_PLETHYSM: dict[type, Callable[..., Any]] = {
    Poly: _c.plethysm_qt,
    QtPoly: _c.plethysm_qt,
    QtFrac: _c.plethysm_macdonald,
    AlphaFrac: _c.plethysm_jack,
    QtRatio: _c.plethysm_ht,
}


def _plethysm(f: Sym, g: Sym) -> Sym:
    """The plethysm `f[g]`, with `f` this element, in `f`'s own basis.

    The same three legs the other Schur-basis operations take. The bases of
    `f` and `g` need not agree — a plethysm composes two elements rather than
    combining two elements of one basis — but the base ring must, and a
    parameter-free element is lifted into `f`'s ring rather than refused.

    The parameters are part of the alphabet, so `p_n` raises them:
    `s_2[q·s_1]` is `q²·s_2`. That is Sage's default, and Sage's `exclude=`,
    which holds a variable constant instead, has no counterpart here.

    # Raises

    Raises `BaseRingError` unless both are over the same base ring, and
    `ValueError` if `g` has a rational coefficient, which `Sym.plethysm`
    refuses for the same reason: plethysm is linear in `f`, so its cleared
    denominator is restored afterwards, and not in `g`, so `g`'s cannot be.
    """
    tag = f.basis
    if not len(f):
        return Sym(tag, {}, f.parameters)
    schur = _convert(f if tag in BASES else _expand(f), "s")
    other = _lift_to(g, schur, "plethysm")
    other = _convert(other if other.basis in BASES else _expand(other), "s")
    call, rows_f, scale_f = _ring_rows(schur, _PLETHYSM, "plethysm", "s")
    _, rows_g, scale_g = _ring_rows(other, _PLETHYSM, "plethysm", "s")
    if scale_g != 1:
        raise ValueError(
            "plethysm needs integer coefficients in its inner argument; it is "
            "not linear there, so a cleared denominator cannot be restored"
        )
    raw = call(rows_f, rows_g)
    kind = _kind(schur)
    acted: Sym
    if kind in (Poly, QtPoly):
        acted = _qt_unpack(raw, schur, "s", scale_f)
    elif kind is QtFrac:
        acted = _mac_element(raw, "s", scale_f)
    elif kind is AlphaFrac:
        acted = _jack_element(raw, "s")
    else:
        acted = _ht_element(raw, "s")
    return _back_to(acted, tag, "plethysm")


#: The Hall inner product entry point for each coefficient ring.
_HALL: dict[type, Callable[..., Any]] = {
    Poly: _c.hall_inner_product_qt,
    QtPoly: _c.hall_inner_product_qt,
    QtFrac: _c.hall_inner_product_macdonald,
    AlphaFrac: _c.hall_inner_product_jack,
    QtRatio: _c.hall_inner_product_ht,
}


def _scalar_pair(f: Sym, g: Sym, what: str, pivot: str = "s") -> tuple[Sym, Sym]:
    """Both sides of a pairing in the `pivot` basis over one coefficient ring.

    Either side may be parameter-free: it is lifted into the other's ring, so
    the two orders of a mixed pairing accept the same pairs rather than only
    the one whose left side already carries the parameters.
    """
    lhs = f.to(pivot)
    if not lhs.parameters and g.parameters:
        lhs = _lift_to(lhs, g, what)
    other = _lift_to(g, lhs, what) if lhs.parameters else g
    return lhs, other.to(pivot)


def _scalar(f: Sym, g: Sym) -> ParamCoefficient | Coefficient:
    """`⟨f, g⟩`, the Hall inner product, as one coefficient.

    Both sides convert to the Schur basis, which is orthonormal for this
    pairing, and the value is the sum of the products of matching coefficients.
    That is bilinear over whatever ring they live in, so the parameters are
    carried and never acted on.

    `g` may be in any basis and may be a `Sym`: the pairing is defined on the
    ring, so two spellings of the same argument give the same number, and there
    is nothing to refuse. `⟨h_2, m_2⟩` is 1, which is a value and not a
    mismatch.

    # Raises

    Raises `ValueError` if the element carries no coefficient class this layer
    knows, and `BaseRingError` if `g` is over a different base ring.
    """
    if not len(f) or not len(g):
        return 0
    schur, other = _scalar_pair(f, g, "the Hall inner product")
    kind = _kind(schur)
    call = _HALL.get(kind)  # type: ignore[arg-type]
    if call is None:
        raise ValueError(
            f"the Hall inner product is not written for "
            f"{kind.__name__ if kind else 'these'} coefficients"
        )
    if kind in (Poly, QtPoly):
        rows_f, scale_f = _qt_pack(schur)
        rows_g, scale_g = _qt_pack(other)
        return _qt_coeff(call(rows_f, rows_g), schur, scale_f * scale_g)
    if kind is QtFrac:
        mac_f, mac_sf = _mac_rows(schur, "the Hall inner product", "s")
        mac_g, mac_sg = _mac_rows(other, "the Hall inner product", "s")
        num, den = call(mac_f, mac_g)
        over = mac_sf * mac_sg
        return QtFrac([(a, b, _unscale(v, over)) for a, b, v in num], den)
    if kind is AlphaFrac:
        num, den, k, tail = call(
            _jack_rows(schur, "the Hall inner product", "s"),
            _jack_rows(other, "the Hall inner product", "s"),
        )
        return AlphaFrac(num, den, k, tail)
    ht_num, ht_den = call(
        _ht_rows(schur, "the Hall inner product", "s"),
        _ht_rows(other, "the Hall inner product", "s"),
    )
    return QtRatio(ht_num, ht_den)


def _deformed_rows(f: Sym, what: str) -> tuple[list[Any], int]:
    """`f`'s terms in the ℚ(q,t) row encoding the deformed pairings read.

    A one-variable polynomial coefficient is read into ℤ[q,t] first: the
    pairings deform in `q` and `t`, so the Hall-Littlewood ring in `t` and the
    LLT ring in `q` sit inside their coefficient field.

    # Raises

    Raises `ValueError` for coefficients in α, which pair under `scalar_jack`.
    """
    if _kind(f) is Poly:
        var = next(c.variable for _, c in f if isinstance(c, Poly))
        if var == "alpha":
            raise ValueError(
                f"{what} needs coefficients in q and t, not alpha; "
                "the alpha ring pairs under scalar_jack"
            )
        cells = []
        for la, c in f:
            assert isinstance(c, Poly)
            items = c.coefficients().items()
            triples = (
                [(0, e, v) for e, v in items]
                if var == "t"
                else [(e, 0, v) for e, v in items]
            )
            cells.append((la, QtPoly(triples)))
        f = Sym(f.basis, cells, ("q", "t"))
    return _mac_rows(f, what, "s")


def _deformed_scalar(
    f: Sym, g: Sym, what: str, call: Callable[..., Any]
) -> ParamCoefficient | Coefficient:
    """The engine of `Sym.scalar_t` and `Sym.scalar_qt`: both sides to the
    Schur basis over one ring, then one contract call.

    The value is a `QtFrac` — the pairings introduce binomial denominators, so
    even a `ℤ[t]` pair lands in the fraction field — except over the `H̃`
    ring, where `scalar_qt` answers in that ring's own `QtRatio`.

    # Raises

    Raises `ValueError` for coefficients the pairing is not defined over: α
    pairs under `scalar_jack`, and the `H̃` ring carries `q^a − t^b` atoms
    that only `scalar_qt` has an entry point for.
    """
    if not len(f) or not len(g):
        return 0
    schur, other = _scalar_pair(f, g, what)
    kind = _kind(schur)
    if kind is QtRatio:
        if call is not _c.scalar_qt:
            raise ValueError(
                f"{what} is not written for QtRatio coefficients; the "
                "Htilde ring pairs under scalar_qt"
            )
        ht_num, ht_den = _c.scalar_qt_ht(
            _ht_rows(schur, what, "s"), _ht_rows(other, what, "s")
        )
        return QtRatio(ht_num, ht_den)
    if kind is AlphaFrac:
        raise ValueError(
            f"{what} needs coefficients in q and t, not AlphaFrac; "
            "the alpha ring pairs under scalar_jack"
        )
    rows_f, sf = _deformed_rows(schur, what)
    rows_g, sg = _deformed_rows(other, what)
    num, den = call(rows_f, rows_g)
    over = sf * sg
    return QtFrac([(a, b, _unscale(v, over)) for a, b, v in num], den)


def _scalar_t(f: Sym, g: Sym) -> ParamCoefficient | Coefficient:
    """`⟨f, g⟩_t`, as one coefficient; the contract is `Sym.scalar_t`'s."""
    return _deformed_scalar(f, g, "scalar_t", _c.scalar_t)


def _scalar_qt(f: Sym, g: Sym) -> ParamCoefficient | Coefficient:
    """`⟨f, g⟩_{q,t}`, as one coefficient; the contract is `Sym.scalar_qt`'s."""
    return _deformed_scalar(f, g, "scalar_qt", _c.scalar_qt)


def _scalar_jack(f: Sym, g: Sym) -> ParamCoefficient | Coefficient:
    """`⟨f, g⟩_α`, as one coefficient; the contract is `Sym.scalar_jack`'s.

    The pivot is the monomial basis — the one Jack's family expansions read —
    rather than `_scalar`'s Schur, so the rows reach `jack_scalar` in the
    encoding it takes.

    # Raises

    Raises `ValueError` for coefficients not in α, naming the pairing that
    does take them.
    """
    if not len(f) or not len(g):
        return 0
    what = "scalar_jack"
    mono, other = _scalar_pair(f, g, what, "m")
    kind = _kind(mono)
    if kind in (QtPoly, QtFrac, QtRatio):
        raise ValueError(
            f"{what} needs coefficients in alpha, not {kind.__name__}; "
            "the q,t rings pair under scalar_t and scalar_qt"
        )
    num, den, k, tail = _c.jack_scalar(
        _jack_rows(mono, what, "m"), _jack_rows(other, what, "m")
    )
    return AlphaFrac(num, den, k, tail)


#: The coproduct entry point for each coefficient ring.
_COPRODUCT: dict[type, Callable[..., Any]] = {
    Poly: _c.coproduct_qt,
    QtPoly: _c.coproduct_qt,
    QtFrac: _c.coproduct_macdonald,
    AlphaFrac: _c.coproduct_jack,
    QtRatio: _c.coproduct_ht,
}


def _coproduct(f: Sym) -> dict[tuple[Partition, Partition], Any]:
    """Δf, as a `{(mu, nu): coefficient}` mapping over the Schur basis of each
    factor.

    Both factors are Schur-basis, as `Sym.coproduct` returns them, and not the
    basis `f` was written in — the result lives in a tensor square this layer
    does not model, so there is nothing to return it in. Sage does write it in
    the element's own basis; that would need the inverse expansion applied to
    both factors of a tensor, which is not an operation here.

    `Δ(s_λ) = Σ c^λ_{μν} s_μ ⊗ s_ν` and those are Littlewood-Richardson
    coefficients, so the parameters are carried and never acted on.

    # Raises

    Raises `ValueError` if the element carries no coefficient class this layer
    knows.
    """
    if not len(f):
        return {}
    schur = _convert(f if f.basis in BASES else _expand(f), "s")
    kind = _kind(schur)
    call = _COPRODUCT.get(kind)  # type: ignore[arg-type]
    if call is None:
        raise ValueError(
            f"the coproduct is not written for "
            f"{kind.__name__ if kind else 'these'} coefficients"
        )
    if kind in (Poly, QtPoly):
        rows, scale = _qt_pack(schur)
        return {k: _qt_coeff(c, schur, scale) for k, c in call(rows)}
    if kind is QtFrac:
        mac, over = _mac_rows(schur, "the coproduct", "s")
        return {
            k: QtFrac([(a, b, _unscale(v, over)) for a, b, v in num], den)
            for k, num, den in call(mac)
        }
    if kind is AlphaFrac:
        return {
            k: AlphaFrac(n, d, j, t)
            for k, n, d, j, t in call(_jack_rows(schur, "the coproduct", "s"))
        }
    return {
        k: QtRatio(n, d)
        for k, n, d in call(_ht_rows(schur, "the coproduct", "s"))
    }


def _qt_coeff(
    rows: list[Any], like: Sym, scale: int
) -> ParamCoefficient:
    """One exponent row back as a coefficient of `like`'s polynomial class,
    dividing out the scale `_qt_pack` cleared — the single-coefficient half of
    `_qt_unpack`.
    """
    if _kind(like) is Poly:
        var = next(c.variable for _, c in like if isinstance(c, Poly))
        return Poly(var, [(b, _unscale(v, scale)) for _, b, v in rows])
    return QtPoly([(a, b, _unscale(v, scale)) for a, b, v in rows])


def _back_to(acted: Sym, tag: str, what: str) -> Sym:
    """A Schur-basis result rewritten in the basis the operand was written in.

    The last leg of every operation that leaves a parametric basis to compute.
    A classical tag is one more change of basis; a parametric one is a change
    of basis into the classical basis its inverse expansion reads, followed by
    that expansion. `_INVERSE` names that basis, which is not always the one
    the family expands in.

    `what` names the operation, so a refusal says which one could not return.
    """
    if tag in BASES:
        return _convert(acted, tag)
    row = _INVERSE.get(tag)
    if row is None:
        raise ValueError(f"{what} cannot be returned in the {tag} basis")
    src, inverse = row
    return inverse(_demote(_convert(acted, src)))


def _demote(f: Sym) -> Sym:
    """An element whose coefficients are fractions with no denominator,
    rewritten over the polynomial class those numerators already are.

    A change of encoding, not of value: `QtFrac.denominator` is the factored
    denominator, and it is empty exactly when there is none. Anything else is
    returned untouched.

    `McdJ` is why this exists. `J` is the integral form, so a Schur-basis
    element on its way back into it has polynomial coefficients — but the
    route there passes through the monomial basis over `ℚ(q,t)`, and comes out
    in that ring's class rather than the one `to_J` reads.
    """
    out: list[tuple[Any, ParamCoefficient]] = []
    for la, c in f:
        if not isinstance(c, QtFrac) or c.denominator:
            return f
        out.append((la, c.numerator))
    return f if not out else Sym(f.basis, out, f.parameters)


def _product(f: Sym, g: Sym) -> Sym:
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
        raise BasisError(
            f"cannot combine {f.basis} with {g.basis}; convert one with .to()"
        )
    _same_ring(f, g)
    tag = f.basis
    if not len(f) or not len(g):
        return Sym(tag, {}, f.parameters)
    left = _convert(f if tag in BASES else _expand(f), "s")
    right = _convert(g if tag in BASES else _expand(g), "s")
    kind = _kind(left)
    schur: Sym
    if kind in (Poly, QtPoly):
        qt_a, scale_a = _qt_pack(left)
        qt_b, scale_b = _qt_pack(right)
        schur = _qt_unpack(
            _c.schur_multiply_qt(qt_a, qt_b), left, "s", scale_a * scale_b
        )
    elif kind is QtFrac:
        mac_a, mac_scale_a = _mac_rows(left, "a product", "s")
        mac_b, mac_scale_b = _mac_rows(right, "a product", "s")
        schur = _mac_element(
            _c.schur_multiply_macdonald(mac_a, mac_b),
            "s",
            mac_scale_a * mac_scale_b,
        )
    elif kind is AlphaFrac:
        schur = _jack_element(
            _c.schur_multiply_jack(
                _jack_rows(left, "a product", "s"),
                _jack_rows(right, "a product", "s"),
            ),
            "s",
        )
    elif kind is QtRatio:
        schur = _ht_element(
            _c.schur_multiply_ht(
                _ht_rows(left, "a product", "s"),
                _ht_rows(right, "a product", "s"),
            ),
            "s",
        )
    else:
        raise ValueError(
            f"a product is not written for "
            f"{kind.__name__ if kind else 'these'} coefficients"
        )
    return _back_to(schur, tag, "a product")


def _unit_like(f: Sym) -> Sym:
    """The multiplicative identity in `f`'s basis, over `f`'s coefficient ring.

    The empty partition indexes 1 in every basis here — `P_∅`, `s_∅` and `m_∅`
    are all the constant 1 — so the unit differs only in which coefficient
    class carries it.

    An empty element records no class, only its parameters, and the unit is
    built over the polynomial ring in those — the smallest of the five that
    contains both 1 and them.
    """
    kind = _kind(f)
    one: ParamCoefficient
    if kind is Poly:
        one = Poly(next(c.variable for _, c in f if isinstance(c, Poly)), {0: 1})
    elif kind is QtPoly:
        one = QtPoly([(0, 0, 1)])
    elif kind is QtFrac:
        one = QtFrac([(0, 0, 1)], [])
    elif kind is AlphaFrac:
        one = AlphaFrac([1], [], 1)
    elif kind is QtRatio:
        one = QtRatio([(0, 0, 1)], [])
    elif len(f.parameters) == 1:
        one = Poly(f.parameters[0], {0: 1})
    elif f.parameters == ("q", "t"):
        one = QtPoly([(0, 0, 1)])
    else:
        raise ValueError(
            f"the unit is not written over {_ring(f.parameters)}"
        )
    # A pair rather than a mapping: `Mapping` is invariant in its value
    # type, so a dict of one coefficient class is not a dict of the union.
    return Sym(f.basis, [((), one)], f.parameters)


def _lift_to(g: Sym, like: Sym, what: str) -> Sym:
    """`g` rewritten over `like`'s coefficient ring, in `g`'s own basis.

    A change of encoding, not of value: an integer is the constant polynomial,
    the fraction with that numerator, or the α-rational with that numerator
    over 1. It is what lets an operation take one operand with parameters and
    one without, the way Sage's `P[2,1].skew_by(s[1])` does.

    A `Sym` is returned unchanged once its base ring matches, since it is
    already over the ring.

    # Raises

    Raises `BaseRingError` if `g` carries a different base ring, and
    `ValueError` for a rational coefficient the target encoding cannot hold —
    `H̃`'s takes integer numerators only.
    """
    if g.parameters:
        _same_ring(g, like)
        return g
    kind = _kind(like)
    cells: list[tuple[Any, ParamCoefficient]] = []
    # `g` carries no parameters, checked above, so its coefficients are numbers.
    for la, c in cast("Iterable[tuple[Any, Coefficient]]", g):
        num, den = Fraction(c).numerator, Fraction(c).denominator
        if kind is Poly:
            var = next(x.variable for _, x in like if isinstance(x, Poly))
            cells.append((la, Poly(var, {0: c})))
        elif kind is QtPoly:
            cells.append((la, QtPoly([(0, 0, c)])))
        elif kind is QtFrac:
            cells.append((la, QtFrac([(0, 0, c)], [])))
        elif kind is AlphaFrac:
            cells.append((la, AlphaFrac([num], [], den)))
        elif kind is QtRatio:
            if den != 1:
                raise ValueError(
                    f"{what}: {c} is not an integer, and this encoding takes "
                    "integer numerators only"
                )
            cells.append((la, QtRatio([(0, 0, num)], [])))
        else:
            raise ValueError(f"{what} is not written for an empty element")
    return Sym(g.basis, cells, like.parameters)


def _constant(f: Sym, c: Scalar) -> Sym:
    """The scalar `c` as an element of `f`'s basis over `f`'s coefficient ring.

    A scalar is the constant it names times the unit, and the unit is indexed
    by the empty partition in every basis here. The scaling goes through
    `_scale`, so a fraction ring reduces the product the way it would for any
    other coefficient rather than being multiplied out here.
    """
    return _scale(_unit_like(f), c)


#: The ring each fraction coefficient class is over, as an element declares it.
_COEFF_RING: dict[type, tuple[str, ...]] = {
    QtFrac: ("q", "t"),
    QtRatio: ("q", "t"),
    AlphaFrac: ("alpha",),
}


def _coeff_element(c: ParamCoefficient) -> Sym:
    """`c` as a one-term element at the empty partition.

    The fraction classes route their own `+`, `-` and `*` through here, so
    the arithmetic runs on the element entry points — the same reduction
    `_add` and `_scale` use. That keeps the result in the encoding
    structural `==` expects, and computes nothing in Python (P4).
    """
    return Sym("m", [((), c)], _COEFF_RING[type(c)])


def _coeff_zero(a: ParamCoefficient) -> ParamCoefficient:
    """The zero of `a`'s class, for the sum an empty element encodes."""
    if isinstance(a, AlphaFrac):
        return AlphaFrac([])
    return QtFrac([]) if isinstance(a, QtFrac) else QtRatio([])


def _coeff_add(a: ParamCoefficient, b: Scalar) -> ParamCoefficient:
    """`a + b` for a fraction coefficient, through the ring's element add."""
    f = _coeff_element(a)
    g = _coeff_element(b) if isinstance(b, type(a)) else _constant(f, b)
    total = _add(f, g)
    if not len(total):
        return _coeff_zero(a)
    return cast("ParamCoefficient", total.coefficient(()))


def _coeff_scale(a: ParamCoefficient, b: Scalar) -> ParamCoefficient:
    """`a * b` for a fraction coefficient, through the ring's element scale."""
    scaled = _scale(_coeff_element(a), b)
    if not len(scaled):
        return _coeff_zero(a)
    return cast("ParamCoefficient", scaled.coefficient(()))


def _same_ring(f: Sym, g: Sym) -> None:
    """Refuse two elements whose coefficients are over different base rings.

    An empty element is over any of them, so it is not refused: it is the zero
    of whichever ring the other one names.

    Without this the mismatch surfaced as `TypeError: unsupported operand
    type(s) for +: 'Poly' and 'QtPoly'` — the coefficient classes' own failure,
    leaking through what is a question about the elements. Both operands can be
    in the same basis, so `BasisError` is not the same refusal and `.to()` is
    not the fix.

    # Raises

    Raises `BaseRingError` unless the two carry the same parameters.
    """
    if len(f) and len(g) and f.parameters != g.parameters:
        raise BaseRingError(
            f"cannot combine an element over {_ring(f.parameters)} with one "
            f"over {_ring(g.parameters)}"
        )


def _ring_pair(f: Sym, g: Sym, what: str) -> tuple[Sym, Sym]:
    """Both operands over one base ring, lifting whichever carries none.

    ℚ sits inside every base ring here, so an element without parameters is an
    element of the other one's ring — written in a different encoding, which is
    all `_lift_to` changes. That is what lets `s([2]) + q * s([2])` mean what it
    says now that both are the same class; before the merge the first operand's
    type refused it.

    An empty element is left alone: it is the zero of whichever ring the other
    one names, and `_lift_to` has no coefficient to read a kind from.

    # Raises

    Raises `BaseRingError` when both carry parameters and they differ, which no
    lifting can fix.
    """
    if f.parameters == g.parameters:
        return f, g
    if not g.parameters and len(g) and len(f):
        return f, _lift_to(g, f, what)
    if not f.parameters and len(f) and len(g):
        return _lift_to(f, g, what), g
    _same_ring(f, g)
    return f, g


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


def _ht_rows(f: Sym, what: str, expect: str) -> list[Any]:
    """The `(partition, numerator, denominator atoms)` rows the `H̃` entry
    points take, read off an element.

    Both directions use this encoding, so a value from either feeds back into
    the other.

    # Raises

    Raises `ValueError` for the wrong basis or coefficient class, and for a
    numerator coefficient that is not an integer. `H̃`'s encoding puts its
    denominator in factored `q^a − t^b` atoms and nothing else, so there is no
    slot for a rational numerator — scaling by `1/2` produces one, and without
    this check it reached the boundary and came back as `TypeError: 'Fraction'
    object cannot be interpreted as an integer`.
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
        num = _qt_rows(coeff.numerator)
        for a, b, v in num:
            if getattr(v, "denominator", 1) != 1:
                raise ValueError(
                    f"{what}: the coefficient of q^{a}t^{b} at {list(la)} is "
                    f"{v}, and this encoding takes integer numerators only"
                )
        rows.append((la, num, list(coeff.denominator)))
    return rows


def _kind(f: Sym) -> type | None:
    """The coefficient class an element's terms are written in, or `None` if it
    has no terms.

    Every term of an element shares one, because a family's entry point builds
    them all through the same wrapper.
    """
    for _, c in f:
        return type(c)
    return None


def _add(f: Sym, g: Sym) -> Sym:
    """`f + g`, termwise in the basis both are written in.

    The fraction kinds go through the contract layer rather than adding here,
    because the coefficient classes compare **structurally**: an unreduced sum
    is the right value in a representation nothing else produces, and `==`
    against a value built another way would then be false. `Poly` and `QtPoly`
    need no such care — a polynomial sum is already canonical.
    """
    if f.basis != g.basis:
        raise BasisError(
            f"cannot combine {f.basis} with {g.basis}; convert one with .to()"
        )
    _same_ring(f, g)
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
    return Sym(f.basis, _add_terms(f, g), f.parameters)


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


def _add_terms(f: Sym, g: Sym) -> dict[Partition, Any]:
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


def _scale(f: Sym, c: Scalar) -> Sym:
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
        num, atoms, over, tail = _alpha_scalar(c)
        out = _c.jack_element_scale(
            _jack_rows(f, "*", f.basis), num, atoms, over, tail
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
            scalar: Any = None
            if isinstance(c, (int, Fraction)) or (
                isinstance(c, Poly) and c.variable == v.variable
            ):
                scalar = c
            elif isinstance(c, QtPoly):
                # The exported `q` and `t` are `QtPoly`; one supported on
                # this element's variable alone is the same scalar in a
                # wider encoding, so it is demoted rather than refused.
                scalar = _as_one_variable(c, v.variable)
            if scalar is None:
                symbol = {"t": "t_hl", "q": "q_llt", "alpha": "alpha"}.get(
                    v.variable, v.variable
                )
                what = (
                    "a polynomial in q and t"
                    if isinstance(c, QtPoly)
                    else type(c).__name__
                )
                raise TypeError(
                    f"cannot scale an element in {v.variable} by {what}; "
                    f"the one-variable {v.variable} is exported as {symbol}"
                )
            w = v * scalar
        else:
            raise TypeError(f"cannot scale a {type(v).__name__} coefficient")
        if w:
            terms[la] = w
    return Sym(f.basis, terms, f.parameters)


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


def _alpha_scalar(c: Scalar) -> tuple[list[int], list[Any], int, list[Any]]:
    """A scalar as `(dense numerator, atoms, scale, tail)` for
    `jack_element_scale`.

    Unlike `_qt_scalar` there is no divisor to hand back: an `AlphaFrac` row
    carries its own integer `scale`, so a rational coefficient goes there and
    crosses unchanged. The tail — the denominator factor only a plethysm
    produces — is carried too: dropping it would scale by a different value
    than the coefficient names.
    """
    if isinstance(c, int):
        return [c], [], 1, []
    if isinstance(c, Fraction):
        return [c.numerator], [], c.denominator, []
    if isinstance(c, Poly) and c.variable == "alpha":
        num, atoms, scale, tail = _dense(c), [], 1, []
    elif isinstance(c, AlphaFrac):
        num, atoms, scale, tail = (
            list(c.numerator),
            list(c.atoms),
            c.scale,
            list(c.tail),
        )
    else:
        raise TypeError(f"cannot scale an element in alpha by {type(c).__name__}")
    lcm = 1
    for v in num:
        if isinstance(v, Fraction):
            lcm = lcm * v.denominator // gcd(lcm, v.denominator)
    return [int(v * lcm) for v in num], atoms, scale * lcm, tail


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

    Accepts a `Sym` in the Schur basis — parameter-free, or with
    coefficients polynomial in `t` — or the rows themselves. The contract
    layer takes integers, so rational coefficients are multiplied up by their
    least common denominator here and divided back out in `_t_element` —
    the same round trip `Sym.to` makes, exact because the expansion is
    linear.
    """
    if isinstance(f, Sym) and not f.parameters:
        if f.basis != expect:
            raise _needs(what, expect, f.basis, True)
        polys = [
            (la, {0: c})
            for la, c in cast("Iterable[tuple[Any, Coefficient]]", f)
        ]
    elif isinstance(f, Sym):
        if f.basis != expect:
            raise _needs(what, expect, f.basis, expect in BASES)
        polys = []
        for la, coeff in f:
            # The expansion is over ℤ[t]. The exported `t` is a `QtPoly`, so
            # a two-variable coefficient supported on `t` alone is demoted;
            # one genuinely in `q`, or in α, has the same term structure and
            # a different meaning, so it is refused rather than read through
            # whichever accessor exists.
            if isinstance(coeff, QtPoly):
                narrowed = _as_one_variable(coeff, "t")
                if narrowed is not None:
                    coeff = narrowed
            if not isinstance(coeff, Poly) or coeff.variable != "t":
                raise ValueError(
                    f"{what} needs coefficients in t alone — the exported "
                    f"t_hl — not {type(coeff).__name__}"
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

    Accepts a `Sym` in the monomial basis — parameter-free, or with
    coefficients in `q` and `t` — or the rows themselves. The
    contract layer takes integer numerators, so rational ones are multiplied
    up by their least common denominator here and divided back out in
    `_mac_element`; the expansion is linear, so the round trip is exact and
    the denominators are untouched by it.
    """
    if isinstance(f, Sym) and not f.parameters:
        if f.basis != expect:
            raise _needs(what, expect, f.basis, True)
        cells: list[Any] = [(la, {(0, 0): c}, ()) for la, c in f]
    elif isinstance(f, Sym):
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
    """The `(partition, numerator, atoms, scale, tail)` rows the Jack inverse
    expansions take.

    Accepts a `Sym` in the monomial basis — parameter-free, or with
    coefficients in α — or the rows themselves. A rational numerator
    coefficient needs no round trip the way `_mac_rows` does: every row
    already carries its own integer `scale`, so the row's least common
    denominator goes there and the value crosses unchanged.
    """
    if isinstance(f, Sym) and not f.parameters:
        if f.basis != expect:
            raise _needs(what, expect, f.basis, True)
        cells: list[Any] = [(la, [c], (), 1, ()) for la, c in f]
    elif isinstance(f, Sym):
        if f.basis != expect:
            raise _needs(what, expect, f.basis, expect in BASES)
        cells = []
        for la, coeff in f:
            # The expansion is over ℚ(α); a coefficient in `t` or in `q` and
            # `t` has the same term structure and a different meaning, so it
            # is refused rather than read through whichever accessor exists.
            if isinstance(coeff, AlphaFrac):
                cells.append(
                    (la, coeff.numerator, coeff.atoms, coeff.scale, coeff.tail)
                )
            elif isinstance(coeff, Poly) and coeff.variable == "alpha":
                cells.append((la, _dense(coeff), (), 1, ()))
            else:
                raise ValueError(
                    f"{what} needs coefficients in alpha, not "
                    f"{type(coeff).__name__}"
                )
    else:
        return list(f)
    rows = []
    for la, num, atoms, scale, tail in cells:
        lcm = 1
        for c in num:
            if isinstance(c, Fraction):
                lcm = lcm * c.denominator // gcd(lcm, c.denominator)
        rows.append(
            (la, [int(c * lcm) for c in num], list(atoms), scale * lcm, list(tail))
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

    Accepts a `Sym` in the Schur basis — parameter-free, or in `q` and `t` —
    or the rows themselves.
    """
    if isinstance(f, Sym) and not f.parameters:
        if f.basis != "s":
            raise ValueError(f"{what} needs a Schur-basis element, not {f.basis}")
        out: list[Any] = []
        for la, c in f:
            # `int()` on a `Fraction` truncates, so a rational coefficient
            # would cross as a different element rather than as an error.
            # These rows are ℤ[q,t]; refusing is the boundary's contract
            # (`docs/policies/failure.md`, P8).
            number = cast("Coefficient", c)
            if number != int(number):
                raise ValueError(
                    f"{what} needs integer coefficients; {tuple(la)} carries {c}"
                )
            out.append((la, [(0, 0, int(number))]))
        return out
    if isinstance(f, Sym):
        if f.basis != "s":
            raise ValueError(f"{what} needs a Schur-basis element, not {f.basis}")
        rows = []
        for la, coeff in f:
            # `∇` takes `(q, t)`-graded rows, and only a `QtPoly` coefficient
            # has them. A Schur-basis `Sym` over any other coefficient type
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


#: How each parametric basis is re-entered: the family's inverse expansion,
#: which takes an element in the basis `EXPANDS_IN` names for that tag.
#: `_back_to` reads this, so every operation that leaves a parametric basis to
#: compute returns to it the same way.
#:
#: Defined here rather than beside `EXPANDS_IN` because the family singletons
#: it names do not exist until this point in the module.
#: For each parametric tag, the classical basis its inverse expansion reads and
#: the function that runs it. That is the basis the family expands in for eight
#: of the nine; `McdJ` is the exception, because `J` is triangular against the
#: Schur basis in the direction the inverse needs, so `EXPANDS_IN` is not this
#: table.
_INVERSE: dict[str, tuple[str, Callable[[Any], Sym]]] = {
    "HLP": ("s", hl.to_P),
    "HLQp": ("s", hl.to_Qp),
    "McdHt": ("s", macdonald.to_Htilde),
    "McdP": ("m", macdonald.to_P),
    "McdQ": ("m", macdonald.to_Q),
    "McdJ": ("s", macdonald.to_J),
    "JackP": ("m", jack.to_P),
    "JackQ": ("m", jack.to_Q),
    "JackJ": ("m", jack.to_J),
}
