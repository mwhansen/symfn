"""Send SIGINT into a running call and demand it stops.

    cargo build --features python
    python3 scripts/check_python_interrupt.py target/debug/libsymfn.dylib

The proposition: **a long call raises `KeyboardInterrupt` when the signal
arrives, not when the call would have finished anyway.** Both halves matter,
and only the first is obvious. An extension that never checks signals still
raises `KeyboardInterrupt` eventually — CPython recorded it and runs the
handler at the first bytecode boundary after the call returns — so a test that
only asserts the exception type passes against the defect it is meant to catch.
What separates them is *when*: this measures the delay and requires it to be a
fraction of the call it interrupted.

That is the shape the defect had. Every entry point held the GIL for its whole
duration with no `check_signals` anywhere in the tree, so Ctrl-C during a
multi-minute conversion under Sage did nothing until the conversion ended
(`docs/record/python-and-sage-interop.md`).

The call is chosen to run for a few seconds and to spend that time inside the
kernel rather than marshalling, so the signal lands somewhere a poll site has
to catch it. Timing is a wall-clock comparison against the same call left
alone, so a slow machine moves both numbers and the ratio still holds.

Needs no Sage. `SIGALRM` delivers the interrupt to this process, which is what
a terminal's Ctrl-C does to a foreground one.
"""

import importlib.machinery
import importlib.util
import os
import signal
import sys
import time

# The interrupt must land well inside the call, and the call must be long
# enough that "stopped early" and "ran to completion" are far apart. Degree 30
# takes a few seconds in a debug build and rather less in release, so the
# assertion is a fraction of the measured baseline rather than a fixed second
# count.
OUTER, INNER = [([5], 1)], [([6], 1)]
FIRE_AFTER = 0.5

# The delay allowed between the signal and the exception, as a share of the
# uninterrupted call. Generous on purpose: what it has to separate is "a poll
# site saw it" from "the call finished first", and those differ by the whole
# length of the call.
BUDGET = 0.5


def load(path):
    """Import the built cdylib under the name its `#[pymodule]` declares."""
    loader = importlib.machinery.ExtensionFileLoader("symfn", path)
    spec = importlib.util.spec_from_loader("symfn", loader)
    module = importlib.util.module_from_spec(spec)
    sys.modules["symfn"] = module
    loader.exec_module(module)
    return module


def time_uninterrupted(symfn):
    start = time.monotonic()
    symfn.plethysm(OUTER, INNER)
    return time.monotonic() - start


def time_to_interrupt(symfn):
    """Fire SIGINT `FIRE_AFTER` seconds in; return the delay until it lands.

    Returns `None` if the call completed without raising, which is the defect
    this script exists to catch.
    """
    signal.signal(signal.SIGALRM, lambda *_: os.kill(os.getpid(), signal.SIGINT))
    signal.setitimer(signal.ITIMER_REAL, FIRE_AFTER)
    start = time.monotonic()
    try:
        symfn.plethysm(OUTER, INNER)
    except KeyboardInterrupt:
        return time.monotonic() - start - FIRE_AFTER
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0)
    return None


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    symfn = load(sys.argv[1])

    baseline = time_uninterrupted(symfn)
    if baseline < 2 * FIRE_AFTER:
        raise SystemExit(
            f"the probe call takes {baseline:.2f}s, too short to interrupt "
            f"{FIRE_AFTER}s in; raise the degree in OUTER/INNER"
        )

    delay = time_to_interrupt(symfn)
    if delay is None:
        raise SystemExit(
            "FAIL: the call ran to completion with SIGINT delivered "
            f"{FIRE_AFTER}s in — nothing is checking signals"
        )

    allowed = BUDGET * baseline
    if delay > allowed:
        raise SystemExit(
            f"FAIL: KeyboardInterrupt arrived {delay:.2f}s after the signal, "
            f"over the {allowed:.2f}s allowed against a {baseline:.2f}s call — "
            "the poll sites are too far apart, or the interrupt only landed "
            "when the call ended"
        )

    print(
        f"ok: {baseline:.2f}s call, interrupted {delay:.3f}s after the signal "
        f"(allowed {allowed:.2f}s)"
    )


if __name__ == "__main__":
    main()
