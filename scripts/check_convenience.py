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

        # The Macdonald and Jack wrappers must carry the basis their entry
        # points return; a wrong tag survives every value check but this one.
        check.equal(sf.macdonald.P(la).basis, "m", f"macdonald.P({la}).basis")
        check.equal(sf.jack.P(la).basis, "m", f"jack.P({la}).basis")
        check.equal(sf.hl.Qp(la).basis, "s", f"hl.Qp({la}).basis")

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
    """Proposition 4 — mixing bases raises rather than converting."""
    for first, second in itertools.permutations(BASES, 2):
        a, b = sf.Sym(first, {(2,): 1}), sf.Sym(second, {(2,): 1})
        for label, op in (
            ("+", lambda: a + b),
            ("-", lambda: a - b),
            ("*", lambda: a * b),
            ("scalar", lambda: a.scalar(b)),
            ("skew_by", lambda: a.skew_by(b)),
        ):
            check.raises(sf.BasisError, op, f"{first} {label} {second}")
        check.equal(a == b, False, f"{first}[2] == {second}[2]")


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
    check_degenerations(sf, sf.symfn, check)
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
