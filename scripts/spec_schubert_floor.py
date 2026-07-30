"""How much of a `schubmult` timing is NOT schubmult.

The §2 and §6 comparisons drive Buch's C `schubmult` out of process and time
the whole `subprocess.run`, which bundles three things:

  1. process startup (fork/exec/dyld),
  2. schubmult's actual computation,
  3. formatting every output term to stdout, and the parent draining the pipe.

(3) is the dangerous one, because it **scales with the number of output
terms** — the very axis the engine comparison is about. A win measured against
a competitor that is also printing 3.2M lines is partly a win against printf.

Measures: the startup floor, and pipe-vs-/dev/null on cases spanning four
orders of magnitude in output size.

Run: sage -python scripts/spec_schubert_floor.py
"""
import subprocess
import time
from sage.all import Permutations, set_random_seed

SB = "/var/tmp/sage-10.9-current/local/bin/schubmult"
REPS = 5


def run(u, v, sink, limit=180):
    args = [SB] + [str(a) for a in u] + ["-"] + [str(a) for a in v]
    best, n = None, None
    for _ in range(REPS):
        t0 = time.perf_counter()
        if sink == "devnull":
            with open("/dev/null", "wb") as fh:
                subprocess.run(args, stdout=fh, stderr=subprocess.DEVNULL,
                               timeout=limit)
        else:
            r = subprocess.run(args, capture_output=True, text=True,
                               timeout=limit)
            n = sum(1 for l in r.stdout.splitlines() if l.strip())
        dt = time.perf_counter() - t0
        best = dt if best is None else min(best, dt)
    return best, n


def stair(k):
    return [2 * i for i in range(1, k + 1)] + [2 * i - 1 for i in range(1, k + 1)]


print("== startup floor: trivial inputs, essentially no work, ~1 line out ==",
      flush=True)
for label, u, v in [("id x id", [1, 2], [1, 2]),
                    ("s1 x s1", [2, 1], [2, 1]),
                    ("stair2^2", stair(2), stair(2))]:
    d_pipe, n = run(u, v, "pipe")
    d_null, _ = run(u, v, "devnull")
    print(f"  {label:<12} pipe={d_pipe * 1000:7.2f}ms  devnull={d_null * 1000:7.2f}ms"
          f"  terms={n}", flush=True)

print("\n== pipe vs /dev/null, across output sizes ==", flush=True)
print(f"  {'case':<26} {'terms':>8} {'pipe':>9} {'devnull':>9} {'printing':>9}"
      f" {'share':>7}", flush=True)

cases = [(f"stair{k}^2", stair(k), stair(k)) for k in [4, 5, 6]]
set_random_seed(1)
pairs = {}
for n in [10, 11, 12]:
    pairs[n] = [(Permutations(n).random_element(),
                 Permutations(n).random_element()) for _ in range(3)]
# chosen for output size, not runtime: 1070 / 30143 / 118822 / 185284 terms
for n, i in [(10, 0), (11, 0), (11, 1), (11, 2)]:
    u, v = pairs[n][i]
    cases.append((f"S_{n}.{i} l={u.length()},{v.length()}", list(u), list(v)))

for label, u, v in cases:
    d_pipe, terms = run(u, v, "pipe")
    d_null, _ = run(u, v, "devnull")
    diff = d_pipe - d_null
    share = 100.0 * diff / d_pipe if d_pipe else 0.0
    print(f"  {label:<26} {terms:>8} {d_pipe:>8.4f}s {d_null:>8.4f}s"
          f" {diff:>8.4f}s {share:>6.1f}%", flush=True)

print("\nDONE", flush=True)
