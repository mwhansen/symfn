"""A/B Sage with symfn substituted for Symmetrica as the conversion backend.

    sage -python scripts/check_backend.py [max_degree]

This is the first test that *replaces* Symmetrica rather than consulting Sage as
an oracle, and the difference matters: every earlier script chose its own inputs,
so it could only find bugs in cases we thought of. Here Sage drives, through its
own dispatch, and the comparison is against the C library the shim displaces --
on identical inputs, in identical code.

Coverage is deliberately indirect as well as direct. The 20 table entries are
checked head-on, but so are the operations that only *reach* a conversion on the
way to something else: products in a non-Schur basis, `scalar`, `expand`,
plethysm, `itensor`, `skew_by`, the antipode, and the Hall-Littlewood, Jack and
Macdonald bases, which are built on top of the classical ones and are where an
integration bug is most likely to surface first.

Both coefficient regimes are exercised, because Sage takes **different code
paths** for them (see `sage.combinat.sf.classical`): over QQ it hands the whole
dict to the backend in one call, and over any other ring it calls once per
partition with coefficient 1 and recombines itself. QQ['t'] therefore tests the
per-monomial path that a Macdonald or Hall-Littlewood user would take.

The two halves run in **separate processes**. Sage memoises conversion morphisms
and basis elements aggressively, so computing one answer and then swapping the
backend in the same process risks comparing a cached value against a fresh one
and calling it agreement.
"""

import subprocess
import sys

_ARGS = sys.argv[1:]
MAX_DEGREE = int(_ARGS[0]) if _ARGS and _ARGS[0] != "--dump" else 6


def dump(mode, max_degree):
    """Print a canonical transcript of many Sage computations, one per line."""
    sys.path.insert(0, "scripts")
    from sage.all import QQ, Partitions, SymmetricFunctions  # noqa: E402
    from sage.combinat.partition import Partition  # noqa: E402
    from sage.combinat.sf import classical  # noqa: E402

    import sage_backend  # noqa: E402

    if mode == "symfn":
        sage_backend.install()
    else:
        classical.init()

    out = []

    def say(tag, value):
        out.append(f"{tag}\t{value}")

    def canon(elt):
        mc = elt.monomial_coefficients()
        return "{" + ", ".join(f"{list(k)}:{v}" for k, v in sorted(mc.items())) + "}"

    NAMES = ["Schur", "monomial", "homogeneous", "elementary", "powersum"]

    # --- the 20 table entries, head-on ------------------------------------
    for src in NAMES:
        for dst in NAMES:
            if src == dst:
                continue
            fn = classical.conversion_functions[(src, dst)]
            for deg in range(0, max_degree + 1):
                for lam in Partitions(deg):
                    d = {Partition(list(lam)): QQ(1)}
                    say(f"table {src}->{dst} {list(lam)}", sorted(fn(d)._monomial_coefficients.items()))
            # a multi-term, rational-coefficient input
            d = {}
            for i, lam in enumerate(Partitions(min(4, max_degree))):
                d[Partition(list(lam))] = QQ(i + 1) / QQ(i + 3)
            if d:
                say(f"table {src}->{dst} mixed", sorted(fn(d)._monomial_coefficients.items()))

    # --- through Sage's own dispatch, in both coefficient regimes ----------
    for ring, label in ((QQ, "QQ"), (QQ["t"], "QQ[t]")):
        Sym = SymmetricFunctions(ring)
        bases = {
            "s": Sym.schur(),
            "m": Sym.monomial(),
            "h": Sym.homogeneous(),
            "e": Sym.elementary(),
            "p": Sym.power(),
        }
        for sn, src in bases.items():
            for dn, dst in bases.items():
                if sn == dn:
                    continue
                for deg in range(0, max_degree + 1):
                    for lam in Partitions(deg):
                        say(f"conv[{label}] {sn}->{dn} {list(lam)}", canon(dst(src[list(lam)])))

    # --- operations that reach a conversion indirectly --------------------
    Sym = SymmetricFunctions(QQ)
    s, m, h, e, p = Sym.schur(), Sym.monomial(), Sym.homogeneous(), Sym.elementary(), Sym.power()
    for deg in range(1, max_degree + 1):
        for lam in Partitions(deg):
            L = list(lam)
            say(f"omega {L}", canon(s[L].omega()))
            say(f"antipode {L}", canon(s[L].antipode()))
            say(f"expand {L}", str(s[L].expand(3)))
            say(f"scalar {L}", str(s[L].scalar(h[L])))
            say(f"skew_by {L}", canon(s[L].skew_by(s[[1]])))
            say(f"h-prod {L}", canon(s(h[L] * h[[1]])))
            say(f"p-prod {L}", canon(s(p[L] * p[[1]])))
            say(f"e-prod {L}", canon(m(e[L] * e[[1]])))
            if deg <= 5:
                say(f"itensor {L}", canon(s[L].itensor(s[L])))
                say(f"coproduct {L}", str(sorted(
                    ((list(a), list(b)), c)
                    for (a, b), c in s[L].coproduct().monomial_coefficients().items())))
    for a in Partitions(min(3, max_degree)):
        for b in Partitions(min(3, max_degree)):
            if sum(a) * sum(b) <= 9:
                say(f"plethysm {list(a)}[{list(b)}]", canon(s[list(a)](s[list(b)])))

    # --- bases built on top of the classical ones -------------------------
    # These are the real integration surface: each is defined by a transition
    # from a classical basis, so a shim that is subtly wrong shows up here even
    # when the direct table checks pass.
    QQt = QQ["t"].fraction_field()
    SymT = SymmetricFunctions(QQt)
    hl = SymT.hall_littlewood()
    jack = SymT.jack()
    QQqt = QQ["q", "t"].fraction_field()
    mac = SymmetricFunctions(QQqt).macdonald()
    for deg in range(1, min(max_degree, 5) + 1):
        for lam in Partitions(deg):
            L = list(lam)
            say(f"HL P {L}", str(SymT.schur()(hl.P()[L])))
            say(f"HL Q {L}", str(SymT.schur()(hl.Q()[L])))
            say(f"Jack P {L}", str(jack.P()[L].expand(2)))
            if deg <= 4:
                say(f"Mac P {L}", str(mac.P()[L].expand(2)))

    print("\n".join(out))


if __name__ == "__main__":
    if len(sys.argv) > 2 and sys.argv[1] == "--dump":
        dump(sys.argv[2], int(sys.argv[3]))
        sys.exit(0)

    results = {}
    for mode in ("symmetrica", "symfn"):
        print(f"running {mode} ...", file=sys.stderr, flush=True)
        r = subprocess.run(
            [sys.executable, __file__, "--dump", mode, str(MAX_DEGREE)],
            capture_output=True,
            text=True,
        )
        if r.returncode != 0:
            print(f"{mode} failed:\n{r.stderr[-4000:]}")
            sys.exit(1)
        results[mode] = r.stdout.splitlines()

    a, b = results["symmetrica"], results["symfn"]
    if len(a) != len(b):
        print(f"transcript lengths differ: {len(a)} vs {len(b)}")
        sys.exit(1)

    bad = [(x, y) for x, y in zip(a, b) if x != y]
    for x, y in bad[:20]:
        print(f"MISMATCH\n  symmetrica: {x}\n  symfn     : {y}")
    print(f"\n{len(a)} computations through Sage, {len(bad)} mismatches")
    sys.exit(1 if bad else 0)
