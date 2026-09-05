"""The E3 compression measurement the Schubert design called for, on the same
cases the incumbents were measured on — it predicts whether E3 beats E2 before
either is tuned. Recorded in docs/record/schubert.md.

What is measured.  Symmetrica's expansion recursion (`algorithmus2`, sb.c:323,
public domain, read at function level) is a tree:

    state = (perm, alphabetindex, stufe)
    len(perm) == 2            -> emit one monomial                    [LEAF]
    perm[0] == len(perm)      -> x_alphabetindex^stufe * recurse on perm[1:]
                                 at alphabetindex+1, stufe = len-2    [DESCEND]
    otherwise                 -> sum over Bruhat covers w -> w t_{1,i}
                                 (running minimum of values > perm[0]),
                                 same alphabetindex, stufe-1          [BRANCH]

Leaf count = S_w(1,...,1) = #pipe dreams = exactly Symmetrica's cost, and E1's.
E3's cost is instead the number of DISTINCT states, because the subtree value
at a state is a function of the state alone.  Two merge granularities:

  states_ps  -- merge on (perm, stufe): sound with no bookkeeping at all.
  states_p   -- merge on perm alone.  Sound because the subtree value depends
                on stufe only through an overall factor x_alphabetindex^stufe
                (stufe is read only at the DESCEND branch), so a node reusing a
                cached perm pays a shift by a power of one variable.  Note
                len(perm) determines alphabetindex, so perm alone pins the
                level.  This is the compression an implementation can reach.

compression = leaves / states. the LR record's analogue
(`docs/record/littlewood-richardson.md`) is "LR tableaux / states produced".

Run: sage -python scripts/spec_schubert_peel.py     (Sage only for the seeded
random permutations of the incumbent ladder and the small-case cross-check; the DP is plain
Python and exact.)
"""
import sys
import time

sys.setrecursionlimit(100000)


# --------------------------------------------------------------------------
# the peel DAG
# --------------------------------------------------------------------------
def covers_at_1(p):
    """perm -> [p t_{1,i}] for the Bruhat covers raising p[0], in the running-
    minimum order of Symmetrica's scan."""
    n = len(p)
    out = []
    maximal = n + 1
    p0 = p[0]
    for i in range(1, n):
        if p0 < p[i] < maximal:
            maximal = p[i]
            q = list(p)
            q[0], q[i] = q[i], q[0]
            out.append(tuple(q))
    return out


def peel_stats(w):
    """Return (leaves, states_p, states_ps, max_depth)."""
    w = tuple(int(x) for x in w)
    while len(w) < 2:
        w = w + (len(w) + 1,)
    memo_p = {}
    stats = {"depth": 0}

    def rec(p, stufe, depth):
        """Leaf count, memoized on perm alone. Sound because a subtree's leaf
        count cannot depend on stufe or level -- those change which monomials
        the leaves carry, not how many there are."""
        if depth > stats["depth"]:
            stats["depth"] = depth
        if p in memo_p:
            return memo_p[p]
        n = len(p)
        if n <= 2:
            r = 1
        elif p[0] == n:
            r = rec(p[1:], len(p) - 2, depth + 1)
        else:
            r = 0
            for q in covers_at_1(p):
                r += rec(q, stufe - 1, depth + 1)
        memo_p[p] = r
        return r

    leaves = rec(w, len(w) - 1, 0)

    # The (perm, stufe) state space, counted in its OWN traversal.
    #
    # It must not share the perm-keyed memo above: adding (p, stufe) to a set
    # while pruning on p alone never explores the subtree below an
    # already-seen perm at a *different* stufe, and so undercounts. That is
    # exactly what the first version of this script did -- it reported 114 for
    # stair4 where the true count is 158 -- and it was caught only when the
    # Rust implementation, which walks the real key, disagreed.
    seen_ps = set()
    stack = [(w, len(w) - 1)]
    while stack:
        st = stack.pop()
        if st in seen_ps:
            continue
        seen_ps.add(st)
        p, stufe = st
        n = len(p)
        if n <= 2:
            continue
        if p[0] == n:
            stack.append((p[1:], len(p) - 2))
        else:
            for q in covers_at_1(p):
                stack.append((q, stufe - 1))

    return leaves, len(memo_p), len(seen_ps), stats["depth"]


def peel_dag(w):
    """Build the merged-on-perm DAG.  Returns (children, roots) where
    children[p] is the list of child perms (with multiplicity = 1 here, since
    the covers of a given perm are distinct)."""
    w = tuple(int(x) for x in w)
    children = {}
    stack = [w]
    while stack:
        p = stack.pop()
        if p in children:
            continue
        n = len(p)
        if n <= 2:
            kids = []
        elif p[0] == n:
            kids = [p[1:]]
        else:
            kids = covers_at_1(p)
        children[p] = kids
        stack.extend(k for k in kids if k not in children)
    return children, w


def peel_memory(w):
    """Peak number of simultaneously-live state values under the natural DFS
    post-order evaluation with refcounting -- the memory bound for E3, where
    every live state holds a whole Schubert element.

    Returns (states, edges, peak_live).
    """
    children, root = peel_dag(w)
    parents = {p: 0 for p in children}
    for p, kids in children.items():
        for k in kids:
            parents[k] += 1
    parents[root] += 1  # the caller holds the root
    edges = sum(len(k) for k in children.values())

    live, peak, refs = 0, 0, {}
    done = set()
    # iterative DFS post-order
    stack = [(root, False)]
    while stack:
        p, expanded = stack.pop()
        if expanded:
            if p in done:
                continue
            done.add(p)
            refs[p] = parents[p]
            live += 1
            peak = max(peak, live)
            for k in children[p]:
                refs[k] -= 1
                if refs[k] == 0:
                    live -= 1
            continue
        if p in done:
            continue
        stack.append((p, True))
        for k in children[p]:
            if k not in done:
                stack.append((k, False))
    return len(children), edges, peak


def peel_stats_nomemo(w, cap=6_000_000):
    """states without ANY merging = tree node count (sanity: >= leaves)."""
    w = tuple(int(x) for x in w)
    cnt = [0]

    def rec(p, stufe):
        cnt[0] += 1
        if cnt[0] > cap:
            raise OverflowError
        n = len(p)
        if n <= 2:
            return
        if p[0] == n:
            rec(p[1:], len(p) - 2)
            return
        for q in covers_at_1(p):
            rec(q, stufe - 1)

    try:
        rec(w, len(w) - 1)
        return cnt[0]
    except OverflowError:
        return None


def stair(k):
    return [2 * i for i in range(1, k + 1)] + [2 * i - 1 for i in range(1, k + 1)]


def row(label, w):
    leaves, sp, sps, _ = peel_stats(w)
    states, edges, peak = peel_memory(w)
    assert states == sp
    print(f"  {label:<20} n={len(w):>2} l={inv(w):>3}  "
          f"pipe_dreams={leaves:<12} states={sp:<7} edges={edges:<7} "
          f"peak_live={peak:<6} compress={leaves / sp:>10.1f}x", flush=True)
    return leaves, sp, sps


def inv(w):
    return sum(1 for i in range(len(w)) for j in range(i + 1, len(w)) if w[i] > w[j])


# --------------------------------------------------------------------------
print("== cross-check: peel leaf count == S_w(1,...,1) ==", flush=True)
try:
    from sage.all import (SchubertPolynomialRing, ZZ, Permutations,
                          set_random_seed, Permutation)
    from sage_guard import require_own_sage

    require_own_sage("Sage's Schubert polynomials")
    X = SchubertPolynomialRing(ZZ)
    HAVE_SAGE = True
except ImportError:
    HAVE_SAGE = False
    print("  (no sage: skipping cross-check and random rows)", flush=True)

if HAVE_SAGE:
    bad = []
    for n in [2, 3, 4, 5, 6]:
        for w in Permutations(n):
            leaves, _, _, _ = peel_stats(list(w))
            e = X(list(w)).expand()
            spec = sum(int(c) for c in e.coefficients())
            if leaves != spec:
                bad.append((list(w), leaves, spec))
    print(f"  exhaustive S_2..S_6 ({sum(1 for n in range(2,7) for _ in Permutations(n))}"
          f" permutations): {len(bad)} mismatches {bad[:3]}", flush=True)

    print("\n== no-merge tree size vs merged states (small, sanity) ==", flush=True)
    for w in [[1, 4, 2, 3], stair(3), [4, 3, 2, 1], stair(4)]:
        leaves, sp, sps, _ = peel_stats(w)
        tree = peel_stats_nomemo(w)
        print(f"  {str(w):<26} leaves={leaves:<8} tree_nodes={tree:<9} "
              f"states_p={sp:<7} states_ps={sps}", flush=True)

print("\n== staircase Grassmannian ladder (the benchmark family) ==",
      flush=True)
for k in range(3, 9):
    row(f"stair{k} (S_{2*k})", stair(k))

print("\n== dominant + w0 controls ==", flush=True)
for n in [6, 8, 10, 12]:
    row(f"w0(S_{n}) dominant", list(range(n, 0, -1)))

if HAVE_SAGE:
    print("\n== random ladder (same seed/order as spec_schubert_walls2.py) ==",
          flush=True)
    set_random_seed(1)
    for n in [10, 11, 12, 13]:
        for i in range(3):
            u = Permutations(n).random_element()
            v = Permutations(n).random_element()
            for tag, x in (("u", u), ("v", v)):
                row(f"S_{n}.{i}{tag}", list(x))

    print("\n== expansion wall (rand S_12 l=33, 84084 monomials) ==",
          flush=True)
    set_random_seed(2)
    for n in [11, 12]:
        w = Permutations(n).random_element()
        row(f"expand-case S_{n}", list(w))

print("\n== bigger, past every incumbent's wall ==", flush=True)
for n in [13, 14, 15, 16]:
    row(f"w0-1 S_{n}", list(range(n, 2, -1)) + [1, 2])

print("\nDONE", flush=True)
