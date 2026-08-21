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

        # The inverse expansions undo the forward ones, term for term, and
        # equal the contract calls they wrap. `s_μ = Σ_λ K_{μλ}(t) P_λ` ties
        # `to_P` to `kostka_foulkes`, which reaches the same matrix through a
        # different entry point.
        in_p = sf.hl.to_P(sf.hl.P(la))
        check.equal(in_p.terms, {tuple(la): 1}, f"hl.to_P(hl.P({la}))")
        check.equal(in_p.basis, "HLP", f"hl.to_P(hl.P({la})).basis")
        in_qp = sf.hl.to_Qp(sf.hl.Qp(la))
        check.equal(in_qp.terms, {tuple(la): 1}, f"hl.to_Qp(hl.Qp({la}))")
        check.equal(in_qp.basis, "HLQp", f"hl.to_Qp(hl.Qp({la})).basis")
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
        in_ht = sf.macdonald.to_Htilde(sf.macdonald.Htilde(la))
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
