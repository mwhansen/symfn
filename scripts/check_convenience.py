"""Hold the convenience layer to the contract layer it is defined over.

    cargo build --features python
    python3 scripts/check_convenience.py

`docs/policies/python.md` P4 says the convenience layer computes nothing: every
method is a composition of contract calls. That is a claim about equality, and
this is where it is checked — each method run against the contract calls it
claims to be, over a sweep rather than at one shape.

Four propositions, in the order they catch things:

1. **Every method equals its contract composition.** A `Sym` operation and the
   hand-written contract sequence agree, term for term, across every partition
   to degree 6 and every basis. This is the one that fails if a convenience
   method ever grows mathematics of its own.
2. **The parameter families degenerate correctly.** `P_λ(x; q, q) = s_λ`,
   `P_λ(x; α = 1) = s_λ`, `Q'_λ(x; 0) = s_λ` and `Q'_λ(x; 1) = h_λ` are
   theorems, and each fails under a `q ↔ t`, `α → 1/α` or `t → 1/t` twist of
   the convention. Checking them from Python is what makes the wrappers' basis
   tags and parameter order honest — the LLT wrappers were tagged Schur when
   the entry points return monomial, and a check of this shape is what finds
   that class of defect.
3. **No convenience name shadows a contract name.** The supported surface is
   flat at `symfn.*` and frozen (P10), so a convenience export that reuses an
   entry point's name silently removes it. `hall_littlewood` did exactly that
   before this check existed, which is why the namespace is `hl`.
4. **Basis identity holds.** Mixing bases raises rather than converting (P7),
   in every binary operation, and `to` round-trips.

**What this does not check is the kernel**, deliberately. Catching a defect in
the mathematics is `docs/policies/validation.md`'s job, discharged by the Rust
suites and the committed fixtures under `tests/fixtures/`; a second oracle
reached through PyO3 would be a slower, narrower copy of one that exists, with
the boundary in the way of every failure it reported. What is owed here is the
boundary and the layer above it, and comparing the two Python layers to each
other is the right instrument for that.

Needs no Sage and no installed wheel: it stages the built cdylib into the
source package, as `scripts/check_convenience_docs.py` does.
"""

import itertools
import sys
from collections import Counter
from fractions import Fraction

from check_convenience_docs import ROOT, stage

BASES = "shepmf"

#: Degree bound for the sweeps. Six is where the LR products stay small enough
#: to run in a second and large enough that a basis with a sign — `f` and `e` —
#: has canceled at least once.
DEGREE = 6


def partitions(n):
    """The partitions of `n`, generated here rather than read from the module,
    so the sweep does not depend on the surface it is checking.
    """
    if n == 0:
        return [()]
    out = []

    def build(rest, cap, acc):
        if rest == 0:
            out.append(tuple(acc))
            return
        for part in range(min(rest, cap), 0, -1):
            build(rest - part, part, acc + [part])

    build(n, n, [])
    return out


def every_shape(bound=DEGREE):
    """Every partition up to `bound`, smallest first."""
    return [la for n in range(bound + 1) for la in partitions(n)]


class Check:
    """A running tally, so one failure does not hide the rest."""

    def __init__(self):
        self.failures = []
        self.checks = 0

    def equal(self, got, want, what):
        """Record `got == want`, naming `what` when it does not hold."""
        self.checks += 1
        if got != want:
            self.failures.append(f"{what}: got {got!r}, want {want!r}")

    def raises(self, exception, call, what):
        """Record that `call` raises `exception`."""
        self.checks += 1
        try:
            call()
        except exception:
            return
        except Exception as other:  # noqa: BLE001  the point is which type
            self.failures.append(f"{what}: raised {other!r}, want {exception}")
            return
        self.failures.append(f"{what}: did not raise {exception}")


def check_compositions(sf, c, check):
    """Proposition 1 — each `Sym` method equals its contract composition."""
    shapes = every_shape()
    for la in shapes:
        a = [(la, 1)]
        check.equal(
            sf.s(la).omega().terms, dict(c.omega(a)), f"omega {la}"
        )
        check.equal(
            sf.s(la).antipode().terms, dict(c.antipode(a)), f"antipode {la}"
        )
        check.equal(
            sf.s(la).coproduct(), dict(c.coproduct(a)), f"coproduct {la}"
        )
        for n in (1, 3):
            check.equal(
                sf.s(la).expand(n),
                dict(c.expand_alphabet(a, "Schur", n)),
                f"expand {la} in {n}",
            )
        for dst in BASES:
            want = (
                {p: Fraction(x, y) for p, (x, y) in c.to_power(a, "Schur")}
                if dst == "p"
                else dict(c.convert_terms(a, "Schur", _long(dst)))
            )
            check.equal(
                sf.s(la).to(dst).terms,
                {k: _exact(v) for k, v in want.items()},
                f"s{la}.to({dst})",
            )

    for la, mu in itertools.combinations_with_replacement(shapes[:12], 2):
        a, b = [(la, 1)], [(mu, 1)]
        check.equal(
            (sf.s(la) * sf.s(mu)).terms,
            dict(c.schur_multiply(a, b)),
            f"s{la}*s{mu}",
        )
        check.equal(
            (sf.m(la) * sf.m(mu)).terms,
            dict(c.monomial_multiply(a, b)),
            f"m{la}*m{mu}",
        )
        check.equal(
            sf.s(la).scalar(sf.s(mu)),
            c.hall_inner_product(a, b),
            f"<s{la},s{mu}>",
        )
        check.equal(
            sf.s(la).skew_by(sf.s(mu)).terms,
            dict(c.skew_by(a, b, "s")),
            f"s{la}/s{mu}",
        )
        check.equal(
            sf.s(la).internal_product(sf.s(mu)).terms,
            dict(c.internal_product(a, b)),
            f"s{la}*_int s{mu}",
        )
        if sum(la) and sum(mu):
            check.equal(
                sf.s(la).plethysm(sf.s(mu)).terms,
                dict(c.plethysm(a, b)),
                f"s{la}[s{mu}]",
            )


def check_products_in_every_basis(sf, c, check):
    """Multiplication routed through Schur agrees with multiplying in Schur.

    The `h`, `e`, `p` and `f` bases have no product entry point of their own,
    so `Sym.__mul__` converts, multiplies and converts back. The claim is that
    the route changes no value, which is what makes it a routing decision
    rather than a second implementation.
    """
    shapes = every_shape(4)
    for basis in BASES:
        for la, mu in itertools.combinations_with_replacement(shapes[:8], 2):
            got = sf.Sym(basis, {la: 1}) * sf.Sym(basis, {mu: 1})
            want = (
                sf.Sym(basis, {la: 1}).to("s") * sf.Sym(basis, {mu: 1}).to("s")
            ).to(basis)
            check.equal(got, want, f"{basis}{la}*{basis}{mu} routed via Schur")


def check_round_trips(sf, check):
    """`to` is exact in both directions, over every ordered pair of bases."""
    for la in every_shape(5):
        start = sf.s(la)
        for first, second in itertools.permutations(BASES, 2):
            check.equal(
                start.to(first).to(second).to("s"),
                start,
                f"s{la} through {first} then {second}",
            )


def check_parametric_conversions(sf, check):
    """Setting the parameter commutes with changing the basis.

    The two sides share no route: `Sym.to` carries `(q,t)`-polynomial
    coefficients through `convert_qt_terms`, while `Sym.to` clears
    denominators and calls the integer conversions. So an agreement here is
    evidence about the parametric route, not about one implementation
    compared with itself.
    """
    targets = [b for b in BASES if b != "p"]
    for la in every_shape(4):
        if not la:
            continue
        one = (
            ("HLP", sf.hl.P(la)),
            ("HLQp", sf.hl.Qp(la)),
            ("JackP", sf.jack.P(la)),
            ("JackQ", sf.jack.Q(la)),
            ("JackJ", sf.jack.J(la)),
        )
        for name, f in one:
            param = "t" if name.startswith("HL") else "alpha"
            for dst in targets:
                moved = f.to(dst)
                check.equal(moved.basis, dst, f"{name}{la}.to({dst!r}) basis")
                for value in (0, 1, 3, Fraction(1, 5)):
                    if param == "alpha" and value == 0:
                        continue  # Jack's denominators vanish at alpha = 0
                    check.equal(
                        moved.at(**{param: value}),
                        f.at(**{param: value}).to(dst),
                        f"{name}{la} -> {dst} at {param}={value}",
                    )
        two = (
            ("McdP", sf.macdonald.P(la)),
            ("McdQ", sf.macdonald.Q(la)),
            ("McdJ", sf.macdonald.J(la)),
            ("McdHt", sf.macdonald.Htilde(la)),
        )
        for name, f in two:
            for dst in targets:
                moved = f.to(dst)
                check.equal(moved.basis, dst, f"{name}{la}.to({dst!r}) basis")
                for qv, tv in ((2, 3), (5, 7)):
                    check.equal(
                        moved.at(q=qv, t=tv),
                        f.at(q=qv, t=tv).to(dst),
                        f"{name}{la} -> {dst} at q={qv}, t={tv}",
                    )
        # The classical limits, reached through the new route rather than
        # through the pivot the family expands in: P_lambda(x; 1) = s_lambda
        # for Jack and P_lambda(x; q, q) = s_lambda for Macdonald. Both fail
        # under the alpha -> 1/alpha and q <-> t twists of the convention.
        check.equal(
            sf.jack.P(la).to("s").at(alpha=1), sf.s(la), f"JackP{la} at alpha=1"
        )
        check.equal(
            sf.macdonald.P(la).to("s").at(q=4, t=4),
            sf.s(la),
            f"McdP{la} at q=t=4",
        )
        # A coefficient that keeps a factored denominator through the change of
        # basis, which the families above reach only when one survives.
        ratio = sf.Sym("s", {tuple(la): sf.QtRatio([(1, 0, 1)], [(1, 1, 1, 1)])}, ("q", "t"))
        for dst in targets:
            check.equal(
                ratio.to(dst).at(q=2, t=3),
                ratio.at(q=2, t=3).to(dst),
                f"q/(q-t)*s{la} -> {dst}",
            )
        # A classical element whose coefficients carry parameters is the case
        # with no expansion in front of it, and the one that used to lose every
        # method to the class it was in.
        scaled = sf.q * sf.m(la) + sf.t * sf.m([1] * sum(la))
        for dst in targets:
            for qv, tv in ((1, 1), (2, 3), (0, 5)):
                check.equal(
                    scaled.to(dst).at(q=qv, t=tv),
                    scaled.at(q=qv, t=tv).to(dst),
                    f"q*m{la} + t*m(1^n) -> {dst} at q={qv}, t={tv}",
                )


def check_hall_littlewood_products(sf, check):
    """The Hall-Littlewood product, against three things that do not share its
    route.

    `P_λ(x; 0) = s_λ` sends the structure constants to the
    Littlewood-Richardson coefficients and `P_λ(x; 1) = m_λ` sends them to the
    monomial ones, so the product at those two values must equal a product this
    library computes over integers. At a generic `t` the check is against
    multiplying the two specialized expansions, which runs the same integer
    Schur product on numbers rather than on polynomials.

    The negative coefficients are the reason this needs pinning at all: they
    put these constants in Z[t], where the classical Hall polynomials counting
    subgroups of abelian p-groups are in N[t], and the two differ by a twist.
    """
    negative = 0
    for n in range(1, 5):
        for mu in partitions(n):
            for nu in partitions(n):
                product = sf.hl.P(list(mu)) * sf.hl.P(list(nu))
                check.equal(product.basis, "HLP", f"P{mu}*P{nu} basis")
                for _, c in product:
                    negative += sum(1 for v in c.coefficients().values() if v < 0)
                check.equal(
                    product.at(t=0),
                    (sf.s(list(mu)) * sf.s(list(nu))).to("s"),
                    f"P{mu}*P{nu} at t=0 is Littlewood-Richardson",
                )
                check.equal(
                    product.at(t=1).to("m"),
                    sf.m(list(mu)) * sf.m(list(nu)),
                    f"P{mu}*P{nu} at t=1 is the monomial product",
                )
                for value in (2, Fraction(1, 3)):
                    left = sf.hl.P(list(mu)).at(t=value).to("s")
                    right = sf.hl.P(list(nu)).at(t=value).to("s")
                    check.equal(
                        product.at(t=value),
                        (left * right).to("s"),
                        f"P{mu}*P{nu} at t={value}",
                    )
    # The sweep has to reach the negative coefficients or it pins nothing.
    check.equal(negative > 0, True, "the sweep saw a negative coefficient")
    # Repeated squaring against repeated multiplication.
    cube = sf.hl.P([2, 1]) ** 3
    check.equal(
        cube,
        sf.hl.P([2, 1]) * sf.hl.P([2, 1]) * sf.hl.P([2, 1]),
        "P[2,1]**3 is three multiplications",
    )
    check.equal(len(sf.hl.P([1]) ** 0), 1, "the zeroth power is the unit")
    # Q' multiplies too, and its own degeneration is at t = 0 rather than 1.
    for la in ([1], [2], [2, 1]):
        for mu in ([1], [2], [1, 1]):
            if sum(la) != sum(mu):
                continue
            got = sf.hl.Qp(la) * sf.hl.Qp(mu)
            check.equal(got.basis, "HLQp", f"Qp{la}*Qp{mu} basis")
            check.equal(
                got.at(t=0),
                (sf.s(la) * sf.s(mu)).to("s"),
                f"Qp{la}*Qp{mu} at t=0 is Littlewood-Richardson",
            )


def check_hall_littlewood_products_against_sage(sf, check):
    """The convenience layer's products against the committed Sage fixture.

    `tests/sage_oracle.rs` already holds the crate's route to the same values.
    This is the same fixture read through the Python side, so what it adds is
    the boundary: the packing into exponent rows, the denominator clearing, and
    the rebuild. Sage multiplies by coercing into the Schur basis and inverting
    the transition matrix; symfn expands through its own forward polynomials
    and back-substitutes, so the two share the definition of P and Q' and
    nothing about how the product is reached.
    """
    fixture = ROOT / "tests" / "fixtures" / "sage_oracle.txt"
    if not fixture.exists():
        check.equal(True, False, "the Sage oracle fixture is missing")
        return
    seen = 0
    for line in fixture.read_text().splitlines():
        tag, _, rest = line.partition(" ")
        if tag not in ("hlpmul", "hlqpmul"):
            continue
        arg, _, body = rest.partition(" ")
        mu, nu = (
            [int(x) for x in side.split(",")] if side else []
            for side in arg.split("|")
        )
        ctor = sf.hl.P if tag == "hlpmul" else sf.hl.Qp
        want = {}
        for token in body.split():
            part, poly = token.split(":")
            key = tuple(int(x) for x in part.split(",")) if part else ()
            want[key] = {
                int(b): int(v) for _, b, v in (x.split(".") for x in poly.split(";"))
            }
        got = {la: c.coefficients() for la, c in ctor(mu) * ctor(nu)}
        check.equal(got, want, f"{tag} {mu} * {nu} against Sage")
        seen += 1
    check.equal(seen > 30, True, f"expected a real sweep, got {seen}")


def check_parametric_products_degenerate(sf, check):
    """The rational-function families' products, at the values where the family
    is classical.

    `P_λ(x; 1) = s_λ` for Jack and `P_λ(x; q, q) = s_λ` for Macdonald, so a
    product of two of them must specialize to the Schur product — computed here
    over integers by a route that shares nothing with the rational-function
    one. `tests/sage_oracle.rs` holds the generic-parameter agreement; this is
    the half that runs without a fixture, and it is the pin that fails under
    the `alpha -> 1/alpha` and `q <-> t` twists.
    """
    shapes = [la for la in every_shape(3) if la]
    for mu in shapes:
        for nu in shapes:
            if sum(mu) != sum(nu):
                continue
            classical = (sf.s(list(mu)) * sf.s(list(nu))).to("m")
            jack = sf.jack.P(list(mu)) * sf.jack.P(list(nu))
            check.equal(jack.basis, "JackP", f"JackP{mu}*P{nu} basis")
            check.equal(
                jack.at(alpha=1), classical, f"JackP{mu}*P{nu} at alpha=1"
            )
            mac = sf.macdonald.P(list(mu)) * sf.macdonald.P(list(nu))
            check.equal(mac.basis, "McdP", f"McdP{mu}*P{nu} basis")
            check.equal(mac.at(q=5, t=5), classical, f"McdP{mu}*P{nu} at q=t=5")


def check_parametric_scalars(sf, check):
    """A scalar adds to a parametric element the way it adds to a `Sym`, and
    the zeroth power is the unit over every coefficient ring.

    `sum` is the case that motivated it: it starts from `0`, so without scalar
    addition the first term raises. The independent side is the
    specialization — adding over ℚ after setting the parameter runs `Sym`'s
    integer arithmetic and shares no entry point with the parametric route.
    """
    for la in every_shape(3):
        if not la:
            continue
        for name, f, kw in (
            ("HLP", sf.hl.P(la), {"t": 3}),
            ("McdP", sf.macdonald.P(la), {"q": 2, "t": 3}),
            ("McdHt", sf.macdonald.Htilde(la), {"q": 2, "t": 3}),
            ("JackP", sf.jack.P(la), {"alpha": 3}),
            ("q*m", sf.q * sf.m(la), {"q": 2, "t": 3}),
        ):
            # `H̃`'s encoding takes integer numerators only, so a rational
            # scalar is refused there rather than answered; that refusal is
            # checked below.
            scalars = (2,) if name == "McdHt" else (2, Fraction(1, 2))
            for c in scalars:
                check.equal(
                    (c + f) - c, f, f"{name}{la}: (c + f) - c is f, c={c}"
                )
                check.equal(c + f, f + c, f"{name}{la}: c + f is f + c, c={c}")
                check.equal(
                    (c - f).at(**kw),
                    c - f.at(**kw),
                    f"{name}{la}: c - f at {kw}, c={c}",
                )
                check.equal(
                    (c + f).at(**kw),
                    c + f.at(**kw),
                    f"{name}{la}: c + f at {kw}, c={c}",
                )
            check.equal(sum([f]), f, f"{name}{la}: sum of one term")
            check.equal(f**0 * f, f, f"{name}{la}: the zeroth power is a unit")
            check.equal(
                (f**0).coefficient([]), 1, f"{name}{la}: the unit is 1 at ()"
            )
            check.equal((f**0).basis, f.basis, f"{name}{la}: unit basis")
        check.raises(
            ValueError,
            lambda la=la: Fraction(1, 2) + sf.macdonald.Htilde(la),
            f"McdHt{la}: a rational scalar is refused, not truncated",
        )


def check_parametric_skew(sf, check):
    """`skew_by` over parameters agrees with the integer route after
    specializing, and the six bases of `g` agree with each other.

    Two independent readings. The specialization runs `Sym.skew_by` over ℚ and
    shares no entry point with the parametric route. The basis sweep runs six
    different rules — Pieri, dual Pieri, Murnaghan-Nakayama and three through
    Littlewood-Richardson — on the same `g`, so an error in any one of them
    shows as a disagreement rather than a wrong answer everywhere.
    """
    for la in every_shape(4):
        if len(la) < 1:
            continue
        for name, f, kw in (
            ("HLP", sf.hl.P(la), {"t": 3}),
            ("McdP", sf.macdonald.P(la), {"q": 2, "t": 3}),
            ("McdHt", sf.macdonald.Htilde(la), {"q": 2, "t": 3}),
            ("JackP", sf.jack.P(la), {"alpha": 3}),
            ("q*m", sf.q * sf.m(la), {"q": 2, "t": 3}),
        ):
            for mu in ([1], [2], [1, 1]):
                if sum(mu) > sum(la):
                    continue
                got = f.skew_by(sf.s(mu))
                check.equal(got.basis, f.basis, f"{name}{la}.skew_by(s{mu}) basis")
                check.equal(
                    got.at(**kw).to("m"),
                    f.at(**kw).to("s").skew_by(sf.s(mu)).to("m"),
                    f"{name}{la}.skew_by(s{mu}) at {kw}",
                )
                # Same g, six rules. `h` and `e` take the Pieri paths and `p`
                # the Murnaghan-Nakayama one, none of which touches
                # Littlewood-Richardson.
                # Compared by subtracting, not by `==`. The fraction
                # coefficient classes compare structurally, and two of these
                # rules reach the same value over a different denominator —
                # `(1-t+q-qt)/(1-qt)` and its multiple by `(1+qt)/(1+qt)`.
                # The difference goes through the contract layer, which
                # reduces.
                for code in ("h", "e", "p", "m", "f"):
                    check.equal(
                        len(f.skew_by(sf.s(mu).to(code)) - got),
                        0,
                        f"{name}{la}.skew_by(s{mu}) via {code}",
                    )
            # The Schur basis is orthonormal for the Hall pairing, so pairing
            # against `s_mu` reads off the coefficient of `s_mu`. The two sides
            # share no entry point: one is the pairing, the other a change of
            # basis.
            for mu in every_shape(sum(la)):
                check.equal(
                    f.scalar(sf.s(mu)),
                    f.to("s").coefficient(mu),
                    f"{name}{la}.scalar(s{mu}) reads the Schur coefficient",
                )
            for code in ("h", "e", "p", "m", "f"):
                check.equal(
                    f.scalar(sf.s(la).to(code)),
                    f.scalar(sf.s(la)),
                    f"{name}{la}.scalar(s{la}) via {code}",
                )


def check_parametric_internal(sf, check):
    """The internal product over parameters agrees with the integer route
    after specializing, is symmetric, and has `h_n` as its identity.

    `s_lambda * h_n = s_lambda` in degree n for the Kronecker product, which is
    a fact about the operation and not about any coefficient, so it holds over
    every ring — and reads the answer without computing a second one.
    """
    for la in every_shape(4):
        if not la:
            continue
        for name, f, kw in (
            ("HLP", sf.hl.P(la), {"t": 3}),
            ("McdP", sf.macdonald.P(la), {"q": 2, "t": 3}),
            ("McdHt", sf.macdonald.Htilde(la), {"q": 2, "t": 3}),
            ("JackP", sf.jack.P(la), {"alpha": 3}),
        ):
            # `h_n` is the Kronecker identity in degree n, which is a fact
            # about the operation and not about any coefficient, so it holds
            # over every ring. Taken in the Schur basis, where both sides can
            # be written.
            schur = f.to("s")
            check.equal(
                len(schur.internal_product(sf.h([sum(la)]).to("s")) - schur),
                0,
                f"{name}{la}: h_{sum(la)} is the Kronecker identity",
            )
            for mu in every_shape(sum(la)):
                other = {
                    "HLP": sf.hl.P,
                    "McdP": sf.macdonald.P,
                    "McdHt": sf.macdonald.Htilde,
                    "JackP": sf.jack.P,
                }[name](mu)
                got = f.internal_product(other)
                check.equal(got.basis, f.basis, f"{name}{la}*{mu} basis")
                check.equal(
                    len(got - other.internal_product(f)),
                    0,
                    f"{name}: {la}*{mu} is symmetric",
                )
                check.equal(
                    got.at(**kw).to("s"),
                    f.at(**kw)
                    .to("s")
                    .internal_product(other.at(**kw).to("s")),
                    f"{name}: {la}*{mu} at {kw}",
                )


def check_parametric_principal_at(sf, check):
    """The specialization at an alphabet drawn from the base ring agrees with
    laying that alphabet out, degenerates to the value at `1^n`, and survives
    specializing the parameters.

    `principal_specialization(n, q=c)` weighs each shape by `s_lambda`'s
    q-analog and substitutes; `evaluate([1, c, ..., c^{n-1}])` lays the
    alphabet out and expands in the monomial basis. Neither route knows the
    other, and at `c = 1` the first must also meet `principal_specialization`,
    which takes a third.
    """
    for la in every_shape(4):
        if not la:
            continue
        for name, f, kw in (
            ("HLP", sf.hl.P(la), {"t": 3}),
            ("McdP", sf.macdonald.P(la), {"q": 2, "t": 3}),
            ("McdHt", sf.macdonald.Htilde(la), {"q": 2, "t": 3}),
            ("JackP", sf.jack.P(la), {"alpha": 3}),
        ):

            def value(v, f=f, kw=kw):
                """One coefficient with the parameters set, as the constant
                term of a one-term element — the coefficient classes take
                their `at` arguments differently.
                """
                one = sf.Sym("m", [((), v)], f.parameters)
                return one.at(**kw).coefficient([])

            check.equal(
                f.principal_specialization(3, q=1),
                f.principal_specialization(3),
                f"{name}{la}: at q = 1 it is the value at 1^n",
            )
            for c in (2, 3):
                check.equal(
                    f.principal_specialization(3, q=c),
                    f.evaluate([c**k for k in range(3)]),
                    f"{name}{la}: at q = {c} it is the value on that alphabet",
                )
            check.equal(
                value(f.principal_specialization(3, q=2)),
                f.at(**kw).principal_specialization(3, q=2),
                f"{name}{la}: at q = 2 agrees with the integer route at {kw}",
            )
    # The parameters themselves are alphabet values: ring elements the
    # q-introducing form has nowhere to put. Setting q = 2
    # afterwards must agree with having specialized first, since `at` is a ring
    # homomorphism and the alphabet is a ring element like any other.
    for la in every_shape(4):
        if not la:
            continue
        f = sf.macdonald.P(la)
        one = sf.Sym("m", [((), f.principal_specialization(3, q=sf.q))], f.parameters)
        check.equal(
            one.at(q=2, t=3).coefficient([]),
            f.at(q=2, t=3).principal_specialization(3, q=2),
            f"McdP{la}: the ring's own q as the alphabet, then q = 2",
        )
    check.raises(
        ValueError,
        lambda: sf.macdonald.P([2]).principal_specialization_q(3),
        "the q-introducing form still refuses a ring that carries q",
    )


def check_parametric_plethysm(sf, check):
    """Plethysm over parameters raises them, is linear and multiplicative in
    its outer argument, and agrees with the integer route where the two
    questions coincide.

    **Specializing does not commute with plethysm in general**, which is the
    whole content of the convention: `t·s_1` composed into `p_2` gives `t²p_2`,
    and setting `t = 3` first gives `3p_2` instead of `9p_2`. So the crossing
    check runs with an inner argument whose coefficients carry no parameter,
    where there is nothing to raise and the two orders do agree.

    The raising itself is pinned by a law rather than by a table:
    `f[t^k·g] = t^{kd}·f[g]` for `f` homogeneous of degree `d`, since
    `p_n` sends the monomial `t^k` to `t^{kn}` and the parts of every `mu`
    contributing to `f` sum to `d`.
    """
    from symfn import Poly

    for name, fam, scalar, kw in (
        ("HLP", sf.hl.P, Poly("t", {1: 1}), {"t": 3}),
        ("McdP", sf.macdonald.P, sf.q, {"q": 2, "t": 3}),
        ("McdHt", sf.macdonald.Htilde, sf.q, {"q": 2, "t": 3}),
    ):
        for la in every_shape(3):
            if not la:
                continue
            f, g = fam(la), fam([2])
            check.equal(
                len(f.plethysm(scalar * g) - scalar ** sum(la) * f.plethysm(g)),
                0,
                f"{name}{la}: p_n raises the scalar to the degree",
            )
            # The inner argument's coefficients are integers here, so nothing
            # is raised and the two orders agree.
            check.equal(
                f.plethysm(sf.s([2])).at(**kw).to("s"),
                f.at(**kw).to("s").plethysm(sf.s([2])),
                f"{name}{la}: agrees with the integer route at {kw}",
            )
            check.equal(
                f.plethysm(g).basis, f.basis, f"{name}{la}: answers in its basis"
            )
        for mu in every_shape(2):
            if not mu:
                continue
            f1, f2, g = fam([2]), fam(mu), fam([2])
            check.equal(
                len((f1 + f2).plethysm(g) - (f1.plethysm(g) + f2.plethysm(g))),
                0,
                f"{name}: linear in f at {mu}",
            )
            check.equal(
                len((f1 * f2).plethysm(g) - f1.plethysm(g) * f2.plethysm(g)),
                0,
                f"{name}: multiplicative in f at {mu}",
            )
    # Jack is the one family whose plethysm puts a tail in its coefficients:
    # alpha + 1 raised is alpha^2 + 1, which no product of linear forms holds.
    # These five coefficients are Sage's, exactly, read as num/den at three
    # values of alpha because a ratio of `Poly`s is not a coefficient class.
    got = sf.jack.P([2]).plethysm(sf.jack.P([2]))
    al = sf.alpha
    one = sf.Poly("alpha", {0: 1})
    square = al**2 + one
    want = {
        (1, 1, 1, 1): (
            24 * al**2 - 24 * al**3,
            (al + one) ** 2 * (al + 2 * one) * (al + 3 * one) * square,
        ),
        (2, 1, 1): (4 * al**3 - 4 * al**2, (al + one) ** 3 * square),
        (2, 2): (
            2 * al**5 + 12 * al**4 + 14 * al**3 + 16 * al**2 + 4 * al,
            (al + one) ** 3 * (2 * al + one) * square,
        ),
        (3, 1): (4 * al - 4 * al**2, (al + one) ** 2 * (3 * al + one)),
        (4,): (one, one),
    }
    for mu, (num, den) in want.items():
        for a in (2, 3, 5):
            check.equal(
                got.coefficient(mu).at(a) * den.at(a),
                num.at(a),
                f"JackP[2][JackP[2]] at {mu}, alpha = {a}",
            )
    check.equal(
        got.coefficient([2, 2]).tail,
        (1, 0, 1),
        "the plethysm's coefficients carry alpha^2 + 1 in the tail",
    )
    check.equal(
        got.coefficient([3, 1]).tail,
        (),
        "and not where the raised factor canceled",
    )


def check_parametric_coproduct(sf, check):
    """The coproduct over parameters agrees with the integer route after
    specializing, and with the identity that defines it.

    `⟨Δf, g ⊗ h⟩ = ⟨f, gh⟩`, and the Schur basis of each factor is orthonormal,
    so the coefficient at `(mu, nu)` is `⟨f, s_mu · s_nu⟩`. That reads the same
    number through the product and the pairing instead — three entry points,
    none of them the coproduct's.
    """
    for la in every_shape(4):
        if not la:
            continue
        for name, f, kw in (
            ("HLP", sf.hl.P(la), {"t": 3}),
            ("McdP", sf.macdonald.P(la), {"q": 2, "t": 3}),
            ("JackP", sf.jack.P(la), {"alpha": 3}),
            ("q*m", sf.q * sf.m(la), {"q": 2, "t": 3}),
        ):
            rows = f.coproduct()

            def value(v, f=f, kw=kw):
                """One coefficient with the parameters set. The coefficient
                classes take their arguments differently, so it is specialized
                as the constant term of a one-term element instead.
                """
                one = sf.Sym("m", [((), v)], f.parameters)
                return one.at(**kw).coefficient([])

            # Zeros are dropped: a coefficient can be a nonzero rational
            # function that vanishes at the point it is specialized to, and
            # the integer route never builds a term for it.
            check.equal(
                {k: c for k, v in rows.items() if (c := value(v))},
                f.at(**kw).to("s").coproduct(),
                f"{name}{la}.coproduct() at {kw}",
            )
            for mu, nu in rows:
                check.equal(
                    rows[mu, nu],
                    f.scalar(sf.s(mu) * sf.s(nu)),
                    f"{name}{la}.coproduct()[{mu},{nu}] is <f, s{mu} s{nu}>",
                )


def check_parametric_alphabet(sf, check):
    """`expand` and `evaluate` over parameters agree with the integer route
    after specializing, and with each other at the all-ones alphabet.

    `f(1, …, 1)` is the sum of the coefficients of the expansion, and the two
    reach it by different engines: the expansion lays out the monomial basis,
    the evaluation runs the Schur one. Summed over the specialized values, so
    the addition is ℚ's and not the coefficient classes'.
    """
    for la in every_shape(4):
        if not la:
            continue
        for name, f, kw in (
            ("HLP", sf.hl.P(la), {"t": 3}),
            ("McdP", sf.macdonald.P(la), {"q": 2, "t": 3}),
            ("McdHt", sf.macdonald.Htilde(la), {"q": 2, "t": 3}),
            ("JackP", sf.jack.P(la), {"alpha": 3}),
            ("q*m", sf.q * sf.m(la), {"q": 2, "t": 3}),
        ):

            def value(v, f=f, kw=kw):
                one = sf.Sym("m", [((), v)], f.parameters)
                return one.at(**kw).coefficient([])

            for n in (len(la), len(la) + 1):
                rows = f.expand(n)
                check.equal(
                    {k: c for k, v in rows.items() if (c := value(v))},
                    f.at(**kw).to("m").expand(n),
                    f"{name}{la}.expand({n}) at {kw}",
                )
                ones = [1] * n
                check.equal(
                    value(f.evaluate(ones)),
                    f.at(**kw).to("s").evaluate(ones),
                    f"{name}{la}.evaluate(1^{n}) at {kw}",
                )
                check.equal(
                    sum(value(v) for v in rows.values()),
                    value(f.evaluate(ones)),
                    f"{name}{la}: expand({n}) sums to the value at 1^{n}",
                )
                # The same number by a third route: weigh each shape by
                # s_lambda(1^n) and never lay out an alphabet.
                check.equal(
                    f.principal_specialization(n),
                    f.evaluate(ones),
                    f"{name}{la}.principal_specialization({n}) is the value at 1^{n}",
                )
            check.equal(
                value(f.dimension()),
                f.at(**kw).to("s").dimension(),
                f"{name}{la}.dimension() at {kw}",
            )
            check.equal(
                value(f.principal_specialization(3)),
                f.at(**kw).to("s").principal_specialization(3),
                f"{name}{la}.principal_specialization(3) at {kw}",
            )
            if name == "HLP":
                # The one family with a free variable for the q this
                # introduces. Setting q = 1 has to give the value at 1^n, and
                # setting t = 0 the classical q-analog, since P(x; 0) = s.
                psq = f.principal_specialization_q(3)
                check.equal(
                    psq.at(1, 3),
                    f.principal_specialization(3).at(3),
                    f"HLP{la}.principal_specialization_q(3) at q=1",
                )
                check.equal(
                    psq.at(2, 0),
                    sf.s(la).principal_specialization_q(3).at(2),
                    f"HLP{la}.principal_specialization_q(3) at t=0",
                )
            else:
                # Everything else has no room for it, and says which of the
                # two reasons applies.
                check.raises(
                    ValueError,
                    lambda f=f: f.principal_specialization_q(3),
                    f"{name}{la}.principal_specialization_q(3) is refused",
                )


def check_parametric_hopf(sf, check):
    """ω and the antipode agree with the classical ones after specializing, and
    ω is still an involution.

    Same independence as the conversions: `Sym.omega` runs over integer
    coefficients through the integer entry points, and the parametric route
    goes out to Schur over polynomial coefficients and back. An element in a
    family's own basis has to come back in it, which the basis assertion pins —
    that leg is an inverse expansion the classical route never runs.
    """
    for la in every_shape(4):
        if not la:
            continue
        cases = [
            ("HLP", sf.hl.P(la), {"t"}),
            ("HLQp", sf.hl.Qp(la), {"t"}),
            ("McdHt", sf.macdonald.Htilde(la), {"q", "t"}),
        ]
        for name, f, _ in cases:
            check.equal(f.omega().omega(), f, f"{name}{la}: omega is an involution")
            for op in ("omega", "antipode"):
                moved = getattr(f, op)()
                check.equal(moved.basis, name, f"{name}{la}.{op}() basis")
                for value in (2, 3, Fraction(1, 2)):
                    if name == "McdHt":
                        got, want = value, value + 1
                        check.equal(
                            moved.at(q=got, t=want),
                            getattr(f.at(q=got, t=want), op)(),
                            f"{name}{la} {op} at q={got}, t={want}",
                        )
                    else:
                        check.equal(
                            moved.at(t=value),
                            getattr(f.at(t=value), op)(),
                            f"{name}{la} {op} at t={value}",
                        )
        scaled = sf.q * sf.m(la)
        for op in ("omega", "antipode"):
            check.equal(
                getattr(scaled, op)().at(q=2, t=3),
                getattr(scaled.at(q=2, t=3), op)(),
                f"q*m{la} {op}",
            )
        # The six families over a rational-function ring, which reach the Schur
        # basis through their own encoding rather than the polynomial one. The
        # specialization is the independent check: setting the parameter first
        # and acting over ℚ runs the integer entry points end to end.
        for name, f, kw in (
            ("McdP", sf.macdonald.P(la), {"q": 2, "t": 3}),
            ("McdQ", sf.macdonald.Q(la), {"q": 2, "t": 3}),
            ("McdJ", sf.macdonald.J(la), {"q": 2, "t": 3}),
            ("JackP", sf.jack.P(la), {"alpha": 3}),
            ("JackQ", sf.jack.Q(la), {"alpha": 3}),
            ("JackJ", sf.jack.J(la), {"alpha": 3}),
        ):
            check.equal(
                f.omega().omega(), f, f"{name}{la}: omega is an involution"
            )
            for op in ("omega", "antipode"):
                moved = getattr(f, op)()
                check.equal(moved.basis, name, f"{name}{la}.{op}() basis")
                check.equal(
                    moved.at(**kw),
                    getattr(f.at(**kw), op)(),
                    f"{name}{la} {op} at {kw}",
                )
        # At α = 1 every Jack basis is the Schur basis, so ω there is the
        # classical involution on `s_lambda` — computed by `Sym.omega`, which
        # shares no entry point with the α route.
        check.equal(
            sf.jack.P(la).omega().at(alpha=1),
            sf.s(la).omega().to("m"),
            f"JackP{la}.omega() at alpha=1 is s{la}.omega()",
        )
        # At q = t every Macdonald P is a Schur function, the same way.
        check.equal(
            sf.macdonald.P(la).omega().at(q=5, t=5),
            sf.s(la).omega().to("m"),
            f"McdP{la}.omega() at q=t is s{la}.omega()",
        )


def check_degenerations(sf, c, check):
    """Proposition 2 — the parameter families hit their classical limits."""
    for la in every_shape(5):
        if not la:
            continue
        schur_in_m = sf.s(la).to("m")

        # P_λ(x; q, q) = s_λ, at two unrelated values so a coincidence at one
        # of them is not mistaken for the identity.
        for value in (3, Fraction(2, 7)):
            check.equal(
                sf.macdonald.P(la).at(q=value, t=value),
                schur_in_m,
                f"macdonald.P({la}).at(q=t={value})",
            )

        # P_λ(x; α = 1) = s_λ, and α = 2 is the zonal polynomial.
        check.equal(
            sf.jack.P(la).at(alpha=1), schur_in_m, f"jack.P({la}).at(alpha=1)"
        )
        check.equal(
            sf.jack.P(la).at(alpha=2), sf.jack.zonal(la), f"jack.zonal({la})"
        )

        # Q'_λ(x; 0) = s_λ, and Q'_λ(x; 1) = h_λ.
        check.equal(sf.hl.Qp(la).at(t=0), sf.s(la), f"hl.Qp({la}).at(t=0)")
        check.equal(
            sf.hl.Qp(la).at(t=1), sf.h(la).to("s"), f"hl.Qp({la}).at(t=1)"
        )

        # K_{λμ}(1) is the Kostka number, and K̃_{λμ}(1,1) agrees with it too.
        for mu in every_shape(sum(la)):
            if sum(mu) != sum(la):
                continue
            check.equal(
                sf.hl.kostka_foulkes(la, mu).at(1),
                c.kostka_number(la, mu),
                f"K_{la}{mu}(1)",
            )

        # The inverse expansions undo the forward ones, term for term, and
        # equal the contract calls they wrap. `s_μ = Σ_λ K_{μλ}(t) P_λ` ties
        # `to_P` to `kostka_foulkes`, which reaches the same matrix through a
        # different entry point.
        in_p = sf.hl.to_P(sf.hl.P(la).to("s"))
        check.equal(in_p.terms, {tuple(la): 1}, f'hl.to_P(hl.P({la}).to("s"))')
        check.equal(in_p.basis, "HLP", f'hl.to_P(hl.P({la}).to("s")).basis')
        in_qp = sf.hl.to_Qp(sf.hl.Qp(la).to("s"))
        check.equal(in_qp.terms, {tuple(la): 1}, f'hl.to_Qp(hl.Qp({la}).to("s"))')
        check.equal(in_qp.basis, "HLQp", f'hl.to_Qp(hl.Qp({la}).to("s")).basis')
        rows = [(tuple(la), [(0, 1)])]
        check.equal(
            sf.hl.to_P(sf.s(la)).terms,
            {mu: sf.Poly("t", c) for mu, c in c.schur_to_hall_littlewood_p(rows)},
            f"hl.to_P(s({la})) against the contract rows",
        )
        check.equal(
            sf.hl.to_Qp(sf.s(la)).terms,
            {mu: sf.Poly("t", c) for mu, c in c.schur_to_hall_littlewood_qp(rows)},
            f"hl.to_Qp(s({la})) against the contract rows",
        )
        for mu, coeff in sf.hl.to_P(sf.s(la)):
            check.equal(
                coeff, sf.hl.kostka_foulkes(la, mu), f"[P_{mu}] s_{la} = K_{la}{mu}(t)"
            )

        # `to_Htilde` undoes `Htilde`, carries the tag Sage prints, and equals
        # the contract rows. The independent check is by value: specializing
        # the coefficients and recombining with the `H̃_μ` at the same point
        # must rebuild `s_λ`, which reads the expansion rather than repeating
        # the solve. `q = 3`, `t = 2/7` keeps every `q^a − t^b` away from zero.
        in_ht = sf.macdonald.to_Htilde(sf.macdonald.Htilde(la).to("s"))
        check.equal(in_ht.terms, {tuple(la): 1}, f"macdonald.to_Htilde(Htilde({la}))")
        check.equal(
            in_ht.basis, "McdHt", f"macdonald.to_Htilde(Htilde({la})).basis"
        )
        qt_rows = [(tuple(la), [(0, 0, 1)])]
        check.equal(
            sf.macdonald.to_Htilde(sf.s(la)).terms,
            {mu: sf.QtRatio(n, d) for mu, n, d in c.schur_to_macdonald_ht(qt_rows)},
            f"macdonald.to_Htilde(s({la})) against the contract rows",
        )
        q, t = 3, Fraction(2, 7)
        rebuilt = sf.Sym("s", {})
        for mu, coeff in sf.macdonald.to_Htilde(sf.s(la)):
            rebuilt += sf.macdonald.Htilde(mu).at(q=q, t=t) * coeff.at(q=q, t=t)
        check.equal(rebuilt, sf.s(la), f"Sigma_mu c_mu Ht_mu at (3, 2/7) = s{la}")

        # `to_J` is the element-wise form of the `schur_in_macdonald_j`
        # table, so the check is the table row it must equal, read through a
        # different entry point.
        in_j = sf.macdonald.to_J(sf.s(la))
        check.equal(in_j.basis, "McdJ", f"macdonald.to_J(s({la})).basis")
        table = dict(c.schur_in_macdonald_j(sum(la)))[tuple(la)]
        check.equal(
            in_j.terms,
            {mu: sf.QtFrac(n, d) for mu, n, d in table},
            f"macdonald.to_J(s({la})) against the whole-degree table",
        )

        # `to_P` and `to_Q` undo `P` and `Q`, carry the tags Sage prints, and
        # equal the contract rows. The independent check is again by value:
        # specializing the coefficients and recombining with `P_mu` at the
        # same point must rebuild `m_lambda`. `q = 3`, `t = 2/7` keeps every
        # `1 - q^a t^b` away from zero.
        in_mp = sf.macdonald.to_P(sf.macdonald.P(la).to("m"))
        check.equal(in_mp.terms, {tuple(la): 1}, f"macdonald.to_P(P({la}))")
        check.equal(in_mp.basis, "McdP", f"macdonald.to_P(P({la})).basis")
        in_mq = sf.macdonald.to_Q(sf.macdonald.Q(la).to("m"))
        check.equal(in_mq.terms, {tuple(la): 1}, f"macdonald.to_Q(Q({la}))")
        check.equal(in_mq.basis, "McdQ", f"macdonald.to_Q(Q({la})).basis")
        mac_rows = [(tuple(la), [(0, 0, 1)], [])]
        check.equal(
            sf.macdonald.to_P(sf.m(la)).terms,
            {mu: sf.QtFrac(n, d) for mu, n, d in c.monomial_to_macdonald_p(mac_rows)},
            f"macdonald.to_P(m({la})) against the contract rows",
        )
        check.equal(
            sf.macdonald.to_Q(sf.m(la)).terms,
            {mu: sf.QtFrac(n, d) for mu, n, d in c.monomial_to_macdonald_q(mac_rows)},
            f"macdonald.to_Q(m({la})) against the contract rows",
        )
        rebuilt = sf.Sym("m", {})
        for mu, coeff in sf.macdonald.to_P(sf.m(la)):
            rebuilt += sf.macdonald.P(mu).at(q=q, t=t) * coeff.at(q=q, t=t)
        check.equal(rebuilt, sf.m(la), f"Sigma_mu c_mu P_mu at (3, 2/7) = m{la}")

        # `to_P`, `to_Q` and `to_J` undo `P`, `Q` and `J`, carry the tags
        # Sage prints, and equal the contract rows. The independent check is
        # again by value: specializing the coefficients at `alpha = 5` and
        # recombining with `P_mu` there must rebuild `m_lambda`. A whole
        # number keeps every hook `u*alpha + v` away from zero, and 5 is not
        # 1 or 2, where the family degenerates.
        for name, forward, inverse, tag in (
            ("P", sf.jack.P, sf.jack.to_P, "JackP"),
            ("Q", sf.jack.Q, sf.jack.to_Q, "JackQ"),
            ("J", sf.jack.J, sf.jack.to_J, "JackJ"),
        ):
            back = inverse(forward(la).to("m"))
            check.equal(
                back.terms, {tuple(la): 1}, f'jack.to_{name}({name}({la}).to("m"))'
            )
            check.equal(
                back.basis, tag, f'jack.to_{name}({name}({la}).to("m")).basis'
            )
        jack_rows = [(tuple(la), [1], [], 1, [])]
        for name, inverse, entry in (
            ("P", sf.jack.to_P, c.monomial_to_jack_p),
            ("Q", sf.jack.to_Q, c.monomial_to_jack_q),
            ("J", sf.jack.to_J, c.monomial_to_jack_j),
        ):
            check.equal(
                inverse(sf.m(la)).terms,
                {mu: sf.AlphaFrac(n, d, k, t) for mu, n, d, k, t in entry(jack_rows)},
                f"jack.to_{name}(m({la})) against the contract rows",
            )
        alpha = 5
        rebuilt = sf.Sym("m", {})
        for mu, coeff in sf.jack.to_P(sf.m(la)):
            rebuilt += sf.jack.P(mu).at(alpha=alpha) * coeff.at(alpha)
        check.equal(rebuilt, sf.m(la), f"Sigma_mu c_mu P_mu at alpha = 5 = m{la}")

        # The Macdonald and Jack wrappers must carry the basis their entry
        # points return; a wrong tag survives every value check but this one.
        # `H̃` runs both ways now that its denominators cross factored, so
        # the round trip closes on the shape it started from rather than only
        # on the tag. `q = 3, t = 2/7` keeps every `q^a - t^b` away from zero.
        ht = sf.macdonald.Htilde(la)
        check.equal(
            sf.macdonald.to_Htilde(ht.to("s")).terms,
            {tuple(la): 1},
            f"to_Htilde(Htilde({la}).to('s'))",
        )
        # The expansion of a genuinely rational `H̃` element: every atom
        # cancels, so `s_lambda` comes back exactly and not merely at a point.
        general = sf.macdonald.to_Htilde(sf.s(la))
        check.equal(
            general.to("s").at(q=3, t=Fraction(2, 7)),
            sf.s(la),
            f"to_Htilde(s({la})).to('s') rebuilds s_{la}",
        )
        check.equal(
            len(general - general), 0, f"to_Htilde(s({la})) minus itself is zero"
        )
        check.equal(
            (sf.q * general).to("s").at(q=3, t=Fraction(2, 7)),
            sf.s(la) * 3,
            f"q*to_Htilde(s({la})) expanded and evaluated",
        )

        # Arithmetic in a parametric basis, against the same arithmetic done
        # after expanding. Scaling and addition go through the contract layer
        # so the coefficients come back reduced; expanding is a separate path,
        # so agreement here is not the two sharing an implementation.
        # `q = 3, t = 2/7` keeps every `1 - q^a t^b` away from zero, and
        # `alpha = 5` every hook `u*alpha + v`; 5 is neither 1 nor 2, where the
        # Jack family degenerates.
        for name, unit, gen, basis, point in (
            ("macdonald", sf.macdonald.P, sf.q, "m", {"q": 3, "t": Fraction(2, 7)}),
            ("jack", sf.jack.P, sf.alpha, "m", {"alpha": 5}),
            ("hl", sf.hl.Qp, sf.t_hl, "s", {"t": 3}),
        ):
            f = unit(la)
            check.equal(
                (gen * f).to(basis),
                gen * f.to(basis),
                f"scaling {name}.P({la}) commutes with expanding",
            )
            check.equal(
                (f + f).to(basis),
                2 * f.to(basis),
                f"{name}.P({la}) + itself commutes with expanding",
            )
            check.equal(len(f - f), 0, f"{name}.P({la}) - itself is zero")
            scalar = gen.at(*(point[v] for v in f.parameters)) if isinstance(
                gen, sf.QtPoly
            ) else gen.at(point[f.parameters[0]])
            check.equal(
                (gen * f).at(**point),
                f.at(**point) * scalar,
                f"{name}.P({la}) scaled by a generator, evaluated",
            )

        check.equal(sf.macdonald.P(la).basis, "McdP", f"macdonald.P({la}).basis")
        check.equal(sf.jack.P(la).basis, "JackP", f"jack.P({la}).basis")
        check.equal(sf.hl.Qp(la).basis, "HLQp", f"hl.Qp({la}).basis")
        check.equal(
            sf.macdonald.P(la).to("m").basis, "m", f'macdonald.P({la}).to("m").basis'
        )
        check.equal(sf.jack.P(la).to("m").basis, "m", f'jack.P({la}).to("m").basis')
        check.equal(sf.hl.Qp(la).to("s").basis, "s", f'hl.Qp({la}).to("s").basis')

    # ∇ of an element equals ∇e_n when that element is e_n = s_{1^n}.
    for n in range(1, 5):
        check.equal(
            sf.macdonald.nabla(sf.s([1] * n)),
            sf.macdonald.nabla_e(n),
            f"nabla(s[1^{n}]) = nabla_e({n})",
        )

    # The LLT wrappers return in the basis their entry points document, and
    # carry the one parameter the family has. Both were wrong at first, and
    # both are invisible to a value check: LLT rows share the `(q, t, c)`
    # encoding with the Macdonald operators and use only the first slot, so
    # reading them as `(q, t)` produces right values under a signature that
    # demands a `t` appearing nowhere.
    for label, element in (
        ("llt.G", sf.llt.G([[1], [1]])),
        ("llt.H", sf.llt.H([1, 1], 2)),
        ("llt.Gtilde", sf.llt.Gtilde([2, 1], 2)),
    ):
        check.equal(element.basis, "m", f"{label} basis")
        check.equal(element.parameters, ("q",), f"{label} parameters")
        check.equal(element.at(q=1).basis, "m", f"{label}.at(q=1)")


def check_new_wrappers(sf, c, check):
    """The wrappers added after the first pass, each against what defines it.

    Every one of these is a composition like the rest, but they are the ones
    whose *identity* is easy to get wrong: an operator that agrees with its
    specialized form, a specialization that agrees at `q = 1`, and a pair of
    functions that are supposed to be inverse.
    """
    for n in range(1, 5):
        e_n = sf.s([1] * n)
        check.equal(
            sf.macdonald.nabla_power(e_n, 1),
            sf.macdonald.nabla(e_n),
            f"nabla_power(e_{n}, 1) = nabla(e_{n})",
        )
        for k in range(1, n):
            check.equal(
                sf.macdonald.delta_prime_ek(k, e_n),
                sf.macdonald.delta_prime_e(k, n),
                f"delta_prime_ek({k}, e_{n}) = delta_prime_e({k}, {n})",
            )
        # Theta raises the degree by k; the Deltas preserve it.
        check.equal(
            {sum(la) for la in sf.macdonald.theta_ek(1, e_n).support()},
            {n + 1},
            f"theta_ek(1, e_{n}) raises the degree",
        )
        check.equal(
            {sum(la) for la in sf.macdonald.delta_ek(1, e_n).support()},
            {n},
            f"delta_ek(1, e_{n}) preserves the degree",
        )

    for la in every_shape(5):
        if not la:
            continue
        for n in (1, 3):
            # The q-analog sums to the plain specialization at q = 1.
            check.equal(
                sf.s(la).principal_specialization_q(n).at(1),
                sf.s(la).principal_specialization(n),
                f"s{la} at 1^{n}, graded then evaluated",
            )
        check.equal(
            sf.jack.structure_constant([1], [1], [2]).at(1),
            2,
            "jack.structure_constant is Stanley's object",
        )

    # The LLT wrappers added later carry the bases their entry points state.
    check.equal(sf.llt.schur([1, 1], 2).basis, "s", "llt.schur basis")
    check.equal(sf.llt.Htilde([1], 2).basis, "m", "llt.Htilde basis")
    check.equal(sf.llt.min_inv([[1], [1]]), c.llt_min_inv([[1], [1]]), "llt.min_inv")

    for w in itertools.permutations(range(1, 5)):
        # `expand` and `from_polynomial` are inverse, and `dimension` counts
        # the monomials `expand` would produce without producing them.
        check.equal(sf.from_polynomial(sf.X(w).expand()), sf.X(w), f"X{w} round trip")
        check.equal(
            sf.X(w).dimension(), len(sf.X(w).expand()), f"X{w}.dimension()"
        )
        check.equal(
            sf.stanley_schur(w).terms,
            dict(c.schubert_to_stanley_schur(w)),
            f"stanley_schur{w}",
        )
        for v in itertools.permutations(range(1, 5)):
            check.equal(
                sf.X(w).pairing(sf.X(v), 4),
                c.schubert_pairing([(w, 1)], [(v, 1)], 4),
                f"<X{w}, X{v}>",
            )
            check.equal(
                dict(sf.X(w).scalar_product(sf.X(v), 4)),
                dict(c.schubert_scalar_product([(w, 1)], [(v, 1)], 4)),
                f"X{w}.scalar_product(X{v}, 4)",
            )


def check_no_shadowing(sf, check):
    """Proposition 3 — no convenience export hides a contract entry point."""
    contract = {n for n in dir(sf.symfn) if not n.startswith("_")}
    for name in contract:
        check.equal(
            getattr(sf, name) is getattr(sf.symfn, name),
            True,
            f"symfn.{name} is the contract entry point, not a convenience "
            "export of the same name",
        )


def check_basis_identity(sf, check):
    """Proposition 4 — mixing bases raises rather than converting.

    `skew_by` and `scalar` are not in the list, and that is deliberate. Neither
    combines two elements of one ring: one is an operator built from `g` and
    applied here, the other is a pairing defined on the ring, and both give one
    answer whatever basis the argument is spelled in. Sage accepts any basis
    for both. What is checked instead is that the six spellings agree.
    """
    for first, second in itertools.permutations(BASES, 2):
        a, b = sf.Sym(first, {(2,): 1}), sf.Sym(second, {(2,): 1})
        for label, op in (
            ("+", lambda: a + b),
            ("-", lambda: a - b),
            ("*", lambda: a * b),
        ):
            check.raises(sf.BasisError, op, f"{first} {label} {second}")
        check.equal(a == b, False, f"{first}[2] == {second}[2]")
        check.equal(
            sf.s([3, 1]).skew_by(sf.s([2]).to(first)),
            sf.s([3, 1]).skew_by(sf.s([2]).to(second)),
            f"s[3,1].skew_by(s[2]) in {first} and in {second}",
        )
        check.equal(
            sf.s([3, 1]).scalar(sf.s([3, 1]).to(first)),
            sf.s([3, 1]).scalar(sf.s([3, 1]).to(second)),
            f"s[3,1].scalar(s[3,1]) in {first} and in {second}",
        )


def check_fifteen_way_to(sf, check):
    """`to` takes all fifteen codes: a parametric target is the family's
    inverse expansion from the classical basis it reads, and `to(x.basis)`
    is the identity for every tag.
    """
    inverses = {
        "HLP": ("s", sf.hl.to_P),
        "HLQp": ("s", sf.hl.to_Qp),
        "McdHt": ("s", sf.macdonald.to_Htilde),
        "McdJ": ("s", sf.macdonald.to_J),
        "McdP": ("m", sf.macdonald.to_P),
        "McdQ": ("m", sf.macdonald.to_Q),
        "JackP": ("m", sf.jack.to_P),
        "JackQ": ("m", sf.jack.to_Q),
        "JackJ": ("m", sf.jack.to_J),
    }
    for la in every_shape(4):
        x = sf.s(la)
        for tag, (pivot, inverse) in inverses.items():
            check.equal(
                x.to(tag),
                inverse(x.to(pivot)),
                f"s{list(la)}.to({tag!r}) is the family method's value",
            )
    identity_cases = [sf.Sym(b, {(2, 1): 1}) for b in BASES] + [
        sf.macdonald.P([2, 1]),
        sf.macdonald.Q([2]),
        sf.macdonald.J([2]),
        sf.macdonald.Htilde([2]),
        sf.hl.P([2, 1]),
        sf.hl.Qp([2]),
        sf.jack.P([2, 1]),
        sf.jack.Q([2]),
        sf.jack.J([2]),
    ]
    for x in identity_cases:
        check.equal(x.to(x.basis), x, f"to({x.basis!r}) is the identity")
    check.raises(
        ValueError,
        lambda: (sf.q * sf.s([2])).to("HLP"),
        "to('HLP') from coefficients in q and t",
    )


def check_parameter_symbols(sf, check):
    """The exported two-variable `q` and `t` scale the one-variable families:
    a `QtPoly` supported on the element's own variable demotes rather than
    refusing.
    """
    for la in every_shape(4):
        if not la:
            continue
        check.equal(
            sf.t * sf.hl.P(la),
            sf.t_hl * sf.hl.P(la),
            f"t * hl.P({list(la)}) is t_hl * hl.P({list(la)})",
        )
        check.equal(
            sf.q * sf.llt.H(la, 2),
            sf.q_llt * sf.llt.H(la, 2),
            f"q * llt.H({list(la)}, 2) is q_llt * llt.H({list(la)}, 2)",
        )
    check.equal(
        sf.hl.to_P(sf.t * sf.s([2])),
        sf.t_hl * sf.hl.to_P(sf.s([2])),
        "hl.to_P accepts the element the exported t builds",
    )
    check.raises(
        TypeError,
        lambda: sf.q * sf.hl.P([2]),
        "q * hl.P([2]) refuses: q is not in the Hall-Littlewood ring",
    )


def check_scalar_division(sf, check):
    """`x / k` is `1/k * x` for a nonzero `int` or `Fraction`, on both the
    classical and the parametric route.
    """
    for la in every_shape(4):
        check.equal(
            sf.s(la) / 2, Fraction(1, 2) * sf.s(la), f"s{list(la)} / 2"
        )
    check.equal(
        sf.macdonald.P([2, 1]) / 2,
        Fraction(1, 2) * sf.macdonald.P([2, 1]),
        "McdP[2,1] / 2",
    )
    check.equal(
        sf.hl.P([2]) / Fraction(3, 2),
        Fraction(2, 3) * sf.hl.P([2]),
        "HLP[2] / (3/2)",
    )
    check.raises(ZeroDivisionError, lambda: sf.s([2]) / 0, "s[2] / 0")
    check.raises(TypeError, lambda: sf.s([2]) / sf.s([1]), "s[2] / s[1]")


def check_coefficient_arithmetic(sf, check):
    """Extracted coefficients have `+`, `-` and `*` that keep canonical form:
    a value built along two routes compares equal structurally, and `at`
    distributes over the arithmetic.
    """
    half = Fraction(1, 2)
    c = sf.macdonald.Q([1]).to("m").coefficient([1])
    a = sf.jack.P([2]).to("m").coefficient([1, 1])
    r = sf.macdonald.to_Htilde(sf.s([2])).coefficient([1, 1])
    for name, x in (("QtFrac", c), ("AlphaFrac", a), ("QtRatio", r)):
        check.equal(x + x, 2 * x, f"{name}: c + c is 2*c")
        check.equal(x - x == 0, True, f"{name}: c - c is 0")
        check.equal(-(-x), x, f"{name}: -(-c) is c")
        check.equal(x + 1 - 1, x, f"{name}: c + 1 - 1 is c")
    check.equal(
        (c + c).at(q=0, t=half),
        2 * c.at(q=0, t=half),
        "QtFrac: at distributes over +",
    )
    check.equal(
        (a * a).at(3), a.at(3) ** 2, "AlphaFrac: at distributes over *"
    )
    check.equal(
        (r + r).at(q=2, t=3), 2 * r.at(q=2, t=3), "QtRatio: at distributes over +"
    )
    check.equal(
        c * sf.m([1]),
        c * sf.macdonald.P([1]).to("m"),
        "a QtFrac scales a classical element by lifting it",
    )
    check.equal(
        (a * sf.s([2])).parameters,
        ("alpha",),
        "an AlphaFrac lifts a classical element into its ring",
    )


def check_principal_specialization_polynomial_q(sf, check):
    """The classical route takes a polynomial `q` and agrees with
    `principal_specialization_q` coefficient by coefficient.
    """
    variable = sf.Poly("q", {1: 1})
    for la in every_shape(4):
        got = sf.s(la).principal_specialization(3, q=variable)
        want = sf.s(la).principal_specialization_q(3)
        check.equal(
            got.coefficients() if isinstance(got, sf.Poly) else (
                {0: got} if got else {}
            ),
            want.coefficients(),
            f"s{list(la)}.principal_specialization(3, q) against the q-analog",
        )


def check_constants_hash_like_their_values(sf, check):
    """Zero and the constants compare and hash as the values they are, in
    every basis and over every coefficient ring.

    The sweep requires `hash(x) == hash(y)` wherever `x == y`, across the
    constant and zero cases of all five coefficient classes and of elements
    built over them.
    """
    from symfn._param import AlphaFrac, Poly, QtFrac, QtPoly, QtRatio

    threes = {
        "int": 3,
        "Poly(t)": Poly("t", {0: 3}),
        "Poly(alpha)": Poly("alpha", {0: 3}),
        "QtPoly": QtPoly({(0, 0): 3}),
        "QtFrac": QtFrac([(0, 0, 3)]),
        "QtRatio": QtRatio([(0, 0, 3)]),
        "AlphaFrac": AlphaFrac([3]),
    }
    zeros = {
        "int": 0,
        "Poly(t)": Poly("t", {}),
        "QtPoly": QtPoly({}),
        "QtFrac": QtFrac([]),
        "QtRatio": QtRatio([]),
        "AlphaFrac": AlphaFrac([]),
    }
    for name, c in threes.items():
        check.equal(c == 3, True, f"{name} constant 3 == 3")
        check.equal(hash(c), hash(3), f"hash of {name} constant 3")
    for name, c in zeros.items():
        check.equal(c == 0, True, f"{name} zero == 0")
        check.equal(hash(c), hash(0), f"hash of {name} zero")
    check.equal(
        AlphaFrac([3, 5]) == 3,
        False,
        "3 + 5*alpha == 3 (the num[:1] equality defect)",
    )

    q = sf.q
    constants = {
        "s": sf.s([]) * 3,
        "h": sf.h([]) * 3,
        "hl": 3 * sf.hl.P([]),
        "qt": 3 * (q**0 * sf.s([])),
        "mcd": 3 * sf.macdonald.P([]),
        "ht": 3 * sf.macdonald.Htilde([]),
        "jack": 3 * sf.jack.P([]),
    }
    for name, x in constants.items():
        check.equal(x == 3, True, f"constant 3 over {name} == 3")
        check.equal(hash(x), hash(3), f"hash of constant 3 over {name}")
    for (na, a), (nb, b) in itertools.combinations(constants.items(), 2):
        check.equal(a == b, True, f"constant 3 over {na} == over {nb}")
    empties = {
        name: x - x
        for name, x in {
            "s": sf.s([2]),
            "qt": q * sf.m([2]),
            "hl": sf.hl.P([1]),
            "mcd": sf.macdonald.P([2]),
            "ht": sf.macdonald.Htilde([2]),
            "jack": sf.jack.P([2]),
        }.items()
    }
    for name, x in empties.items():
        check.equal(x == 0, True, f"zero over {name} == 0")
        check.equal(hash(x), hash(0), f"hash of zero over {name}")
    for (na, a), (nb, b) in itertools.combinations(empties.items(), 2):
        check.equal(a == b, True, f"zero over {na} == zero over {nb}")
    check.equal(sf.hl.P([1]) ** 0 == 1, True, "hl.P([1])**0 == 1")
    check.equal(sf.jack.P([]) == 1, True, "jack.P([]) == 1")
    check.equal(sf.s([2]) == 3, False, "s[2] == 3")
    check.equal(sf.s([]) * 3 == Fraction(1, 2), False, "3 == 1/2")


def check_malformed_construction(sf, check):
    """The constructor refuses what it cannot hold, with the typed exception
    each malformed shape names.
    """
    from symfn._param import QtFrac

    q, t_hl = sf.q, sf.t_hl
    cases = [
        (
            ValueError,
            lambda: sf.Sym("m", [((2,), q), ((2,), q)], ("q", "t")),
            "a repeated shape with a parametric coefficient",
        ),
        (
            TypeError,
            lambda: sf.Sym("m", {(2,): q}),
            "a QtPoly coefficient with no parameters declared",
        ),
        (
            TypeError,
            lambda: sf.Sym("m", {(2,): 1}, ("q", "t")),
            "an int coefficient under declared parameters",
        ),
        (
            TypeError,
            lambda: sf.Sym(
                "m", [((2,), q), ((1, 1), QtFrac([(0, 0, 1)]))], ("q", "t")
            ),
            "a mix of QtPoly and QtFrac coefficients",
        ),
        (
            ValueError,
            lambda: sf.Sym("m", {(2,): t_hl}, ("alpha",)),
            "a Poly in t under parameters ('alpha',)",
        ),
        (
            ValueError,
            lambda: sf.Sym("m", {(2,): t_hl}, ("q", "t")),
            "a Poly in t under parameters ('q', 't')",
        ),
        (
            TypeError,
            lambda: sf.Sym("m", {(2,): 1.5}),
            "a float coefficient",
        ),
    ]
    for exception, call, what in cases:
        check.raises(exception, call, what)
    check.equal(
        sf.Sym("m", [((2,), 1), ((2,), Fraction(1, 2))]),
        Fraction(3, 2) * sf.m([2]),
        "numeric coefficients on a repeated shape accumulate",
    )


def check_evaluate_refuses_rational_alphabets(sf, check):
    """`evaluate` refuses a non-integer alphabet with a typed error on both
    routes, rather than leaking PyO3's conversion failure.
    """
    for name, element in (
        ("s[2]", sf.s([2])),
        ("McdP[2]", sf.macdonald.P([2])),
        ("HLP[2]", sf.hl.P([2])),
        ("JackP[2]", sf.jack.P([2])),
    ):
        check.raises(
            TypeError,
            lambda element=element: element.evaluate([Fraction(1, 2), 1]),
            f"{name}.evaluate over a rational alphabet",
        )


def check_deformed_pairings(sf, check):
    """`scalar_t`, `scalar_qt` and `scalar_jack` are the pairings the families
    are orthogonal under, degenerate to `scalar` at their classical points,
    and take a mixed pair in either order.
    """
    shapes = [la for la in every_shape(4) if la]
    for la in shapes:
        for mu in shapes:
            want = 1 if la == mu else 0
            check.equal(
                sf.macdonald.P(la).scalar_qt(sf.macdonald.Q(mu)),
                want,
                f"<McdP{list(la)}, McdQ{list(mu)}>_qt is {want}",
            )
            check.equal(
                sf.jack.P(la).scalar_jack(sf.jack.Q(mu)),
                want,
                f"<JackP{list(la)}, JackQ{list(mu)}>_alpha is {want}",
            )
            if la != mu:
                check.equal(
                    sf.hl.P(la).scalar_t(sf.hl.P(mu)),
                    0,
                    f"<HLP{list(la)}, HLP{list(mu)}>_t is 0",
                )
    # The Hall-Littlewood norm: `<P_λ, P_λ>_t · b_λ(t) = 1` with
    # `b_λ = Π_i Π_{j ≤ m_i} (1 − t^j)` — the value `scalar` gets wrong.
    for la in shapes:
        prod = sf.hl.P(la).scalar_t(sf.hl.P(la))
        for mult in Counter(la).values():
            power = 1
            for j in range(1, mult + 1):
                power = power * sf.t
                prod = prod * (1 - power)
        check.equal(prod, 1, f"<HLP{list(la)}, HLP{list(la)}>_t inverts b")
    # Each pairing degenerates to `scalar` at its classical point.
    third, half = Fraction(1, 3), Fraction(1, 2)
    for la in shapes:
        for mu in shapes:
            f, g = sf.s(la), sf.s(mu)
            hall = f.scalar(g)
            check.equal(
                f.scalar_t(g).at(q=half, t=0),
                hall,
                f"<s{list(la)}, s{list(mu)}>_t at t = 0 is the Hall value",
            )
            check.equal(
                f.scalar_qt(g).at(q=third, t=third),
                hall,
                f"<s{list(la)}, s{list(mu)}>_qt at q = t is the Hall value",
            )
            check.equal(
                f.scalar_jack(g).at(alpha=1),
                hall,
                f"<s{list(la)}, s{list(mu)}>_alpha at alpha = 1 is the Hall value",
            )
    # A mixed pair works in both orders — the parameter-free side lifts.
    for pairing, parametric in [
        ("scalar", sf.jack.P([2])),
        ("scalar_t", sf.hl.P([2])),
        ("scalar_qt", sf.macdonald.P([2])),
        ("scalar_jack", sf.jack.P([2])),
    ]:
        classical = sf.s([2])
        check.equal(
            getattr(classical, pairing)(parametric),
            getattr(parametric, pairing)(classical),
            f"{pairing} takes a mixed pair in either order",
        )
    # The `H̃` ring's `q^a − t^b` denominators answer in `QtRatio`; pinned
    # against the polynomial route by evaluation, since the two fraction
    # classes compare within themselves.
    ratio = sf.Sym(
        "s", [((1,), sf.QtRatio([(0, 0, 1)], [(1, 1, 1, 1)]))], ("q", "t")
    )
    got = ratio.scalar_qt(ratio).at(q=half, t=third)
    base = sf.s([1]).scalar_qt(sf.s([1])).at(q=half, t=third)
    check.equal(
        got,
        base / (half - third) ** 2,
        "scalar_qt over the Htilde ring matches the Frac route by evaluation",
    )
    check.raises(
        ValueError,
        lambda: sf.jack.P([2]).scalar_qt(sf.jack.P([2])),
        "scalar_qt refuses alpha coefficients, naming scalar_jack",
    )
    check.raises(
        ValueError,
        lambda: sf.hl.P([2]).scalar_jack(sf.hl.P([2])),
        "scalar_jack refuses t coefficients, naming the q,t pairings",
    )
    check.raises(
        sf.BaseRingError,
        lambda: sf.hl.P([2]).scalar_t(sf.jack.P([2])),
        "scalar_t refuses a cross-ring pair",
    )


def check_to_power_parametric(sf, check):
    """`to("p")` reaches every coefficient ring: it commutes with `at`, the
    family inverse takes the value back, and a scaled coefficient rides along.
    """
    third, half = Fraction(1, 3), Fraction(1, 2)
    families = [
        (sf.macdonald.P, "McdP", {"q": half, "t": third}),
        (sf.macdonald.Htilde, "McdHt", {"q": half, "t": third}),
        (sf.hl.P, "HLP", {"t": third}),
        (sf.hl.Qp, "HLQp", {"t": third}),
        (sf.jack.P, "JackP", {"alpha": half}),
        (sf.jack.J, "JackJ", {"alpha": half}),
    ]
    for ctor, tag, point in families:
        for la in every_shape(4):
            if not la:
                continue
            x = ctor(la)
            check.equal(
                x.to("p").at(**point),
                x.at(**point).to("p"),
                f"{tag}{list(la)}.to('p') commutes with at",
            )
            check.equal(
                x.to("p").to(tag),
                x,
                f"{tag}{list(la)} comes back from the power sums",
            )
    check.equal(
        (sf.q * sf.m([2])).to("p"),
        sf.q * sf.m([2]).to("p"),
        "a coefficient rides through to('p')",
    )
    check.equal(
        (sf.t_hl * sf.hl.P([2])).to("p").at(t=third),
        third * sf.hl.P([2]).at(t=third).to("p"),
        "the one-variable ring reaches p and keeps its scalar",
    )
    check.equal(
        sf.Sym("p", [], ("q", "t")).to("p"),
        0,
        "the empty parametric element reaches p",
    )


def check_partial_at(sf, check):
    """`at` takes a nonempty subset of the parameters by name, substituting
    those and keeping the rest.
    """
    third, half = Fraction(1, 3), Fraction(1, 2)
    # The q = 0 degeneration of Macdonald P is Hall-Littlewood P — a theorem
    # that fails under the `q ↔ t` twist. The guard is for shapes whose
    # partial value is already constant, where the element comes back
    # parameter-free.
    for la in every_shape(5):
        if not la:
            continue
        x = sf.macdonald.P(la).at(q=0)
        lhs = x.at(t=third) if x.parameters else x
        check.equal(
            lhs,
            sf.hl.P(la).at(t=third).to("m"),
            f"McdP{list(la)} at q = 0 degenerates to HLP",
        )
    # Partial then the rest equals both at once, along both orders.
    for la in every_shape(4):
        if not la:
            continue
        f = sf.macdonald.P(la)
        want = f.at(q=0, t=third)
        x = f.at(q=0)
        check.equal(
            x.at(t=third) if x.parameters else x,
            want,
            f"McdP{list(la)}: q then t equals both at once",
        )
        y = f.at(t=0)
        check.equal(
            y.at(q=half) if y.parameters else y,
            f.at(q=half, t=0),
            f"McdP{list(la)}: t then q equals both at once",
        )
    # The (q,t)-Kostka value at t = 1, the one-variable polynomial that was
    # unreachable while `at` demanded every parameter.
    check.equal(
        sf.macdonald.qt_kostka([3, 1], [2, 1, 1]).at(t=1),
        2 + sf.q,
        "qt_kostka at t = 1 is the one-variable polynomial",
    )
    # The refusals: a value that takes a denominator out of its factored
    # class, a name the element does not carry, and no value at all.
    check.raises(
        ValueError,
        lambda: sf.macdonald.P([2]).to("m").coefficient([1, 1]).at(q=half),
        "a partial value off the atom lattice refuses",
    )
    check.raises(
        ValueError,
        lambda: sf.macdonald.to_Htilde(sf.s([2])).coefficient([1, 1]).at(t=0),
        "a monomial denominator refuses",
    )
    check.raises(
        TypeError,
        lambda: sf.jack.P([2]).at(t=1),
        "a name outside the element's parameters refuses",
    )
    check.raises(
        TypeError,
        lambda: sf.macdonald.P([2]).at(),
        "at with nothing to set refuses",
    )


def check_llt_skew_tuples(sf, check):
    """`llt.G` and `llt.min_inv` take `(outer, inner)` pairs beside plain
    shapes. The values on skew tuples are pinned against Sage in
    `tests/sage_oracle.rs`; what the layer owes is the spellings.
    """
    check.equal(
        sf.llt.G([([1], []), ([1], [])]),
        sf.llt.G([[1], [1]]),
        "an empty inner is the straight shape",
    )
    check.equal(
        sf.llt.G([([2, 1], [1]), ([1], [])]),
        sf.llt.G([([2, 1], [1]), [1]]),
        "pairs and plain shapes mix in one tuple",
    )
    # `(2)/(1)` is `(1)` moved one column right, and content shifts move
    # nothing `inv` sees.
    check.equal(
        sf.llt.G([([2, 1], [1]), ([2], [1])]),
        sf.llt.G([([2, 1], [1]), [1]]),
        "a content translation of a component changes nothing",
    )
    check.equal(
        sf.llt.min_inv([([2, 1], [1]), ([2], [1])]),
        sf.llt.min_inv([([2, 1], [1]), [1]]),
        "min_inv agrees across the translation",
    )
    check.raises(
        ValueError,
        lambda: sf.llt.G([([1], [2])]),
        "an inner not contained in its outer refuses",
    )


def check_tailed_coefficients_scale_exactly(sf, check):
    """A coefficient whose denominator carries a tail — the factor only a
    plethysm produces — scales an element by its whole value, tail included.
    """
    c = sf.jack.P([2]).plethysm(sf.jack.P([2])).coefficient([2, 2])
    check.equal(bool(c.tail), True, "the plethysm coefficient carries a tail")
    check.equal(
        (c * sf.jack.P([1])).coefficient([1]),
        c,
        "scaling by a tailed coefficient keeps the tail",
    )


def check_schubert(sf, c, check):
    """The Schubert type against its contract calls, over `S_4`."""
    perms = [w for w in itertools.permutations(range(1, 5))]
    for w in perms:
        check.equal(
            sf.X(w).expand(), dict(c.schubert_expand([(w, 1)])), f"X{w}.expand()"
        )
    for u, v in itertools.combinations_with_replacement(perms[:8], 2):
        check.equal(
            (sf.X(u) * sf.X(v)).terms,
            dict(c.schubert_multiply([(u, 1)], [(v, 1)])),
            f"X{u}*X{v}",
        )


def _long(code):
    """The contract layer's name for a one-letter basis code."""
    return {
        "s": "Schur",
        "h": "homogeneous",
        "e": "elementary",
        "p": "powersum",
        "m": "monomial",
        "f": "forgotten",
    }[code]


def _exact(v):
    return int(v) if isinstance(v, Fraction) and v.denominator == 1 else v


def main():
    """Run every proposition. Returns an exit code."""
    if stage() is None:
        print("no extension module; build one: cargo build --features python")
        return 1
    sys.path.insert(0, str(ROOT / "python"))
    import symfn as sf

    check = Check()
    check_compositions(sf, sf.symfn, check)
    check_products_in_every_basis(sf, sf.symfn, check)
    check_round_trips(sf, check)
    check_parametric_conversions(sf, check)
    check_parametric_scalars(sf, check)
    check_parametric_skew(sf, check)
    check_parametric_coproduct(sf, check)
    check_parametric_internal(sf, check)
    check_parametric_principal_at(sf, check)
    check_parametric_plethysm(sf, check)
    check_parametric_alphabet(sf, check)
    check_parametric_hopf(sf, check)
    check_hall_littlewood_products(sf, check)
    check_hall_littlewood_products_against_sage(sf, check)
    check_parametric_products_degenerate(sf, check)
    check_degenerations(sf, sf.symfn, check)
    check_new_wrappers(sf, sf.symfn, check)
    check_no_shadowing(sf, check)
    check_basis_identity(sf, check)
    check_fifteen_way_to(sf, check)
    check_parameter_symbols(sf, check)
    check_scalar_division(sf, check)
    check_coefficient_arithmetic(sf, check)
    check_principal_specialization_polynomial_q(sf, check)
    check_constants_hash_like_their_values(sf, check)
    check_malformed_construction(sf, check)
    check_evaluate_refuses_rational_alphabets(sf, check)
    check_deformed_pairings(sf, check)
    check_to_power_parametric(sf, check)
    check_partial_at(sf, check)
    check_llt_skew_tuples(sf, check)
    check_tailed_coefficients_scale_exactly(sf, check)
    check_schubert(sf, sf.symfn, check)

    if check.failures:
        print(f"convenience: {len(check.failures)} of {check.checks} checks failed")
        for line in check.failures[:40]:
            print(f"  {line}")
        if len(check.failures) > 40:
            print(f"  ... and {len(check.failures) - 40} more")
        return 1
    print(f"convenience: {check.checks} checks pass against the contract layer")
    return 0


if __name__ == "__main__":
    sys.exit(main())
