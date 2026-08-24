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
from fractions import Fraction

from check_convenience_docs import ROOT, stage

BASES = "shepmf"

#: Degree bound for the sweeps. Six is where the LR products stay small enough
#: to run in a second and large enough that a basis with a sign — `f` and `e` —
#: has cancelled at least once.
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

    The two sides share no route: `Param.to` carries `(q,t)`-polynomial
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
        ratio = sf.Param("s", {tuple(la): sf.QtRatio([(1, 0, 1)], [(1, 1, 1, 1)])}, ("q", "t"))
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
        jack_rows = [(tuple(la), [1], [], 1)]
        for name, inverse, entry in (
            ("P", sf.jack.to_P, c.monomial_to_jack_p),
            ("Q", sf.jack.to_Q, c.monomial_to_jack_q),
            ("J", sf.jack.to_J, c.monomial_to_jack_j),
        ):
            check.equal(
                inverse(sf.m(la)).terms,
                {mu: sf.AlphaFrac(n, d, k) for mu, n, d, k in entry(jack_rows)},
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
            # The q-analogue sums to the plain specialization at q = 1.
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
    check_parametric_hopf(sf, check)
    check_hall_littlewood_products(sf, check)
    check_hall_littlewood_products_against_sage(sf, check)
    check_parametric_products_degenerate(sf, check)
    check_degenerations(sf, sf.symfn, check)
    check_new_wrappers(sf, sf.symfn, check)
    check_no_shadowing(sf, check)
    check_basis_identity(sf, check)
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
