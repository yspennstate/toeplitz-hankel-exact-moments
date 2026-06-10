"""
Temperley-Lieb cross-check, adapted from the owner's "Quantum_Gossip" R code
(references/Quantum_Gossip_original.R).

The R code builds all planar (non-crossing) perfect matchings between 2N boundary
points (N bottom "inputs", N top "outputs"), verifies their count is the Catalan
number Cat_N = (2N)!/((N+1)! N!), and represents the Temperley-Lieb generators e_i
as 0/1-weighted matrices on the diagram basis with loop factor delta = 2 (each closed
loop formed when stacking diagrams contributes a factor 2), then checks the TL axioms
    e_i^2 = 2 e_i,    e_i e_j = e_j e_i (|i-j| >= 2),    e_i e_{i+-1} e_i = e_i.

Relevance to this project: the non-crossing pairings are EXACTLY the volume-1 class
of our structural theorem (paper Thm 3.1), and the Catalan counts 2, 5, 14, 42, 132,
429 appearing in our volume spectra are the dimensions of the same TL diagram bases.
This module re-derives those counts and the TL structure independently in Python
(port + simplification of the R logic), and verifies against our results files.

Diagram encoding here: a non-crossing perfect matching of {0,..,2N-1} as a partner
tuple (same convention as tk_core).
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from fractions import Fraction
from math import comb
import json
import os

import tk_core


def noncrossing_pairings(n2k):
    """All non-crossing perfect matchings, by direct recursion (the R code grows them
    edge-by-edge from the first uncovered point; this is the same tree, condensed)."""
    if n2k == 0:
        yield ()
        return
    # point 0 pairs with an odd offset partner j so both sides are even-sized
    for j in range(1, n2k, 2):
        for left in noncrossing_pairings(j - 1):
            for right in noncrossing_pairings(n2k - j - 1):
                pm = [0] * n2k
                pm[0] = j
                pm[j] = 0
                for a in range(j - 1):
                    pm[1 + a] = 1 + left[a]
                for b in range(n2k - j - 1):
                    pm[j + 1 + b] = j + 1 + right[b]
                yield tuple(pm)


def catalan(k):
    return comb(2 * k, k) // (k + 1)


# ---- Temperley-Lieb structure on planar (N bottom, N top) diagrams ----------
# A TL diagram on N strands = non-crossing perfect matching of 2N points where
# bottom points are 0..N-1 (left to right) and top points are N..2N-1 (right to left
# along the cycle, matching the standard "rectangle" boundary order).

def tl_diagrams(N):
    return list(noncrossing_pairings(2 * N))


def compose(d1, d2, N):
    """Stack diagram d2 on top of d1 (both partner tuples on 2N cycle points).
    Returns (diagram, loops). Cycle convention: bottom 0..N-1 left->right, then top
    N..2N-1 right->left; so bottom point i meets d2's bottom-mirrored top point."""
    # union-find over 3 layers: bottom(0..N-1), middle(N..2N-1), top(2N..3N-1)
    parent = list(range(3 * N))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(x, y):
        rx, ry = find(x), find(y)
        if rx != ry:
            parent[rx] = ry
            return False
        return True  # closed a loop

    loops = 0
    # d1 edges: bottom layer uses points 0..N-1; its top boundary cycle point N+t
    # (t = 0..N-1, right-to-left) is the middle strand position N + (2N-1 - (N+t)) =
    # middle index (N-1-t) counted left-to-right. Map cycle point p of d1:
    #   p < N  -> node p          (bottom)
    #   p >= N -> node N + (2N-1-p)   (middle, left-to-right)
    def n1(p):
        return p if p < N else N + (2 * N - 1 - p)

    # d2: bottom boundary attaches to the same middle strands, top to the top layer:
    #   p < N  -> node N + p      (middle)
    #   p >= N -> node 2N + (2N-1-p)  (top, left-to-right)
    def n2(p):
        return N + p if p < N else 2 * N + (2 * N - 1 - p)

    for i in range(2 * N):
        if i < d1[i]:
            if union(n1(i), n1(d1[i])):
                loops += 1
        if i < d2[i]:
            if union(n2(i), n2(d2[i])):
                loops += 1
    # read off the composed pairing on bottom+top boundary
    boundary = list(range(N)) + [2 * N + t for t in range(N)]
    groups = {}
    for b in boundary:
        groups.setdefault(find(b), []).append(b)
    pm = [0] * (2 * N)

    def cyc(node):
        # map node back to cycle coordinates of the composed diagram
        return node if node < N else N + (2 * N - 1 - (node - N))  # node-N in 0..N-1 top l->r

    for g in groups.values():
        assert len(g) == 2, "boundary points must pair up"
        a, b = g
        ca = a if a < N else N + (3 * N - 1 - a)   # top node 2N+t -> cycle N + (N-1-t)
        cb = b if b < N else N + (3 * N - 1 - b)
        pm[ca] = cb
        pm[cb] = ca
    return tuple(pm), loops


def tl_generator(i, N):
    """e_i: cup-cap between strands i, i+1 (1-indexed i in 1..N-1)."""
    pm = [0] * (2 * N)
    # bottom cup between bottom points i-1, i ; top cap between the matching top points
    b1, b2 = i - 1, i
    t1, t2 = N + (N - 1 - (i - 1)), N + (N - 1 - i)  # cycle coords of top points above b1,b2
    pm[b1] = b2
    pm[b2] = b1
    pm[t1] = t2
    pm[t2] = t1
    # all other strands go straight through
    for s in range(N):
        if s in (b1, b2):
            continue
        top = N + (N - 1 - s)
        pm[s] = top
        pm[top] = s
    return tuple(pm)


def check_all(N=5):
    diags = tl_diagrams(N)
    cN = catalan(N)
    print(f"N={N}: TL diagrams found {len(diags)}, Catalan {cN}", "OK" if len(diags) == cN else "FAIL")
    assert len(diags) == cN
    index = {d: i for i, d in enumerate(diags)}

    # generator matrices with delta = 2 (the R code's 2^loops entries)
    import numpy as np
    mats = {}
    for i in range(1, N):
        E = np.zeros((cN, cN), dtype=object)
        gi = tl_generator(i, N)
        for j, d in enumerate(diags):
            comp, loops = compose(d, gi, N)
            E[index[comp], j] += Fraction(2) ** loops
        mats[i] = E
    ok = True
    for i in range(1, N):
        ok &= (mats[i] @ mats[i] == 2 * mats[i]).all()
    for i in range(1, N):
        for j in range(1, N):
            if abs(i - j) >= 2:
                ok &= (mats[i] @ mats[j] == mats[j] @ mats[i]).all()
    for i in range(1, N - 1):
        ok &= (mats[i] @ mats[i + 1] @ mats[i] == mats[i]).all()
        ok &= (mats[i + 1] @ mats[i] @ mats[i + 1] == mats[i + 1]).all()
    print(f"N={N}: Temperley-Lieb axioms (delta=2): {'ALL HOLD' if ok else 'FAIL'}")
    assert ok
    return diags


def crosscheck_with_results(n2k):
    """Volume-1 pairings in our results == non-crossing pairings (as sets)."""
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", f"volumes_T{n2k}.json")
    if not os.path.exists(path):
        print(f"results for T{n2k} not found, skipping")
        return
    with open(path) as f:
        data = json.load(f)
    vol1 = set()
    for key, rec in data["orbits"].items():
        if Fraction(rec["volume"]) == 1:
            rep = tuple(json.loads(key))
            vol1.update(tk_core.dihedral_images(rep))
    nc = set(noncrossing_pairings(n2k))
    print(f"T{n2k}: volume-1 set == non-crossing set: {vol1 == nc} "
          f"({len(nc)} = Cat_{n2k//2} = {catalan(n2k//2)})")
    assert vol1 == nc


if __name__ == "__main__":
    for N in (3, 4, 5):
        check_all(N)
    for n2k in (6, 8, 10, 12):
        crosscheck_with_results(n2k)
    print("\nAll Temperley-Lieb cross-checks passed: the volume-1 class of the paper's")
    print("Theorem 3.1 is exactly the TL diagram basis of the Quantum_Gossip code.")
