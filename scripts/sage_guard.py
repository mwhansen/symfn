"""The one statement of when a Sage comparison is honest.

Sage reaches this library through its optional backend whenever symfn is
installed (README.md, "Under Sage"). A script that uses Sage as an oracle, or
as the control arm of a benchmark, then compares symfn against symfn: it
passes, it proves nothing, and the symptom is a run of ratios near 1.0x. The
failure is silent in both directions, which is why it is refused here rather
than left to the invoker to remember.

`SAGE_DISABLE_SYMFN` has to be in the process environment *before* Sage
starts, because `Feature.is_present` caches -- so this cannot set it, only
check it and exit. Call it once, right after the Sage import:

    from sage_guard import require_own_sage

    require_own_sage("Sage's Jack bases")

`scripts/check_sage_guards.py` fails when a Sage-importing script under
`scripts/` neither calls this nor names itself in that file's allowlist, so a
new script cannot quietly reopen the hole.

Sibling import with no path setup: Python puts the running script's own
directory on `sys.path[0]`, and `sage -python scripts/foo.py` is plain Python.
"""

import sys


def require_own_sage(what):
    """Exit unless Sage will answer out of its own code.

    `what` names the part of Sage that would otherwise be symfn, so the
    message says which comparison was about to be vacuous.

    A stock Sage with no backend installed needs nothing and passes straight
    through.
    """
    try:
        from sage.libs.symfn import is_available
    except ImportError:
        return  # stock Sage, with no backend to disable
    if is_available():
        sys.exit(
            f"set SAGE_DISABLE_SYMFN=1 in the environment: {what} would be "
            "symfn, so this would compare symfn against itself"
        )


def describe():
    """Say which Sage is answering, and exit if it would be symfn.

    Printed at the top of the Sage job in `.github/workflows/sage.yml`, so a
    run's log states what its oracle was: a stock Sage with no backend
    module, or one that carries the backend and has it disabled. The third
    state -- backend present and live -- is the vacuous one, and
    `require_own_sage` exits on it.
    """
    from importlib.util import find_spec

    from sage.version import version

    if find_spec("sage.libs.symfn") is None:
        print(f"Sage {version}: stock, no symfn backend module")
        return
    require_own_sage("every family the backend covers")
    print(f"Sage {version}: symfn backend present and disabled by SAGE_DISABLE_SYMFN")


if __name__ == "__main__":
    describe()
