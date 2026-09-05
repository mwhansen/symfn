"""Measured walls for docs/record/jack.md: where Sage's Jack
implementation stops.

Run:  sage -python spec_jack_walls.py <mode>
Modes: pm  jm  jp  norms  products  ps  api

Each mode runs in its own Sage process so per-family caches don't leak
across families; within a family the whole-degree sweep DOES share Sage's
caches, which is the honest unit (that is how a user computes a table).
Single runs, per-item SIGALRM.  Order-of-magnitude walls, not benchmarks.
"""

import signal
import sys
import time

from sage.all import SymmetricFunctions, Partitions, QQ

from sage_guard import require_own_sage  # noqa: E402

require_own_sage("Sage's Jack bases")

TIMEOUT = 120


class Timeout(Exception):
    pass


def _handler(signum, frame):
    raise Timeout()


signal.signal(signal.SIGALRM, _handler)


def timed(thunk):
    """(seconds, result) or (None, None) on timeout."""
    t0 = time.time()
    try:
        signal.alarm(TIMEOUT)
        r = thunk()
        return time.time() - t0, r
    except Timeout:
        return None, None
    finally:
        signal.alarm(0)


def setup():
    R = QQ["t"].fraction_field()
    Sym = SymmetricFunctions(R)
    jack = Sym.jack()
    return Sym, jack


def degree_sweep(target_basis_of, source, label, degrees):
    """Time the whole-degree table: every lambda |- n expanded in the target."""
    print(f"# {label}: whole-degree table, per-item timeout {TIMEOUT}s")
    print(f"{'n':>4} {'p(n)':>6} {'total(s)':>10} {'max(s)':>9}  worst shape / terms")
    for n in degrees:
        total = 0.0
        worst = (0.0, None, 0)
        dead = None
        for la in Partitions(n):
            dt, val = timed(lambda la=la: target_basis_of(source[la]))
            if dt is None:
                dead = la
                break
            total += dt
            if dt > worst[0]:
                worst = (dt, la, len(val))
        if dead is not None:
            print(f"{n:>4} {Partitions(n).cardinality():>6} {'>120 at ' + str(list(dead)):>10}")
            break
        print(
            f"{n:>4} {Partitions(n).cardinality():>6} {total:>10.3f} {worst[0]:>9.3f}  "
            f"{list(worst[1])} / {worst[2]}"
        )
        sys.stdout.flush()
    print()


def main():
    mode = sys.argv[1]
    Sym, jack = setup()
    P, J, Q, Qp = jack.P(), jack.J(), jack.Q(), jack.Qp()
    m, s, p = Sym.monomial(), Sym.schur(), Sym.powersum()

    if mode == "api":
        # What Sage has, by introspection (never by reading source).
        import sage.combinat.sf.sfa as sfa

        el = P.an_element()
        names = sorted(
            n for n in dir(el) if "jack" in n.lower() or "zonal" in n.lower() or "scalar" in n.lower()
        )
        print("element methods mentioning jack/zonal/scalar:", names)
        print("bases:", [str(b) for b in (P, Q, J, Qp)])
        z = Sym.base_ring()
        print("base ring:", z)
        try:
            Z = SymmetricFunctions(QQ).zonal()
            print("zonal basis exists over QQ:", Z)
        except Exception as e:
            print("zonal:", e)

    elif mode == "pm":
        degree_sweep(m, P, "P -> m", range(6, 16))
    elif mode == "jm":
        degree_sweep(m, J, "J -> m", range(6, 16))
    elif mode == "jp":
        degree_sweep(p, J, "J -> p  (the Goulden-Jackson pipeline unit)", range(6, 16))
    elif mode == "ps":
        # The research-gaps.md row was P -> s; recalibrate it on mains.
        degree_sweep(s, P, "P -> s", range(6, 16))
    elif mode == "qpm":
        degree_sweep(m, Qp, "Qp -> m", range(6, 14))
    elif mode == "norms":
        print("# <J_la, J_la>_alpha for every la |- n (scalar_jack)")
        for n in range(6, 16):
            t0 = time.time()
            try:
                signal.alarm(3 * TIMEOUT)
                for la in Partitions(n):
                    J[la].scalar_jack(J[la])
                signal.alarm(0)
            except Timeout:
                print(f"{n:>4}  >360")
                break
            print(f"{n:>4}  {time.time() - t0:>8.3f}s")
            sys.stdout.flush()
    elif mode == "products":
        # Stanley's g^nu_{la,mu}: one J*J product expanded back in J.
        print("# J[la]*J[mu] re-expanded in J  (structure constants g^nu)")
        cases = [
            ([2, 1], [2, 1]),
            ([2, 2], [2, 1]),
            ([3, 2], [2, 2]),
            ([3, 2, 1], [3, 2, 1]),
            ([4, 3], [4, 3]),
            ([4, 3, 2], [4, 3, 2]),
        ]
        for la, mu in cases:
            dt, val = timed(lambda la=la, mu=mu: J(J[la] * J[mu]))
            shown = f"{dt:>8.3f}s  {len(val)} terms" if dt is not None else "   >120s"
            print(f"J{la} * J{mu}: {shown}")
            sys.stdout.flush()
    else:
        raise SystemExit(f"unknown mode {mode}")


if __name__ == "__main__":
    main()
