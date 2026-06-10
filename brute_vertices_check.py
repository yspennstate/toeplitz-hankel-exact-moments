"""Independent per-orbit validation at the disputed order T8 (and H8 control):
1. Enumerate vertices by BRUTE FORCE over all d-subsets of constraints (own Fraction
   solver; zero shared code with the DD path); check vertex-set equality with DD.
2. Compute volume of the vertex hull with scipy ConvexHull (float, independent
   triangulation); compare to exact volume to ~1e-9.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from fractions import Fraction
from itertools import combinations
import numpy as np
from scipy.spatial import ConvexHull
import tk_core
import exact_volume


def solve_square(rows, rhs):
    """Solve d x d Fraction system; return None if singular."""
    d = len(rows)
    M = [list(map(Fraction, rows[i])) + [Fraction(rhs[i])] for i in range(d)]
    for c in range(d):
        piv = None
        for r in range(c, d):
            if M[r][c] != 0:
                piv = r
                break
        if piv is None:
            return None
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        for r in range(d):
            if r != c and M[r][c] != 0:
                f = M[r][c] / pv
                for cc in range(c, d + 1):
                    M[r][cc] -= f * M[c][cc]
    return tuple(M[i][d] / M[i][i] for i in range(d))


def brute_vertices(halfspaces, d):
    """All feasible basic solutions of the halfspace system (a.v + b >= 0)."""
    verts = set()
    m = len(halfspaces)
    for combo in combinations(range(m), d):
        rows = [halfspaces[i][0] for i in combo]
        rhs = [-halfspaces[i][1] for i in combo]
        v = solve_square(rows, rhs)
        if v is None:
            continue
        ok = True
        for (a, b) in halfspaces:
            s = Fraction(b)
            for ai, vi in zip(a, v):
                if ai:
                    s += ai * vi
            if s < 0:
                ok = False
                break
        if ok:
            verts.add(v)
    return verts


def halfspaces_for(pm, kind):
    if kind == "T":
        rows = tk_core.dedupe_rows(tk_core.toeplitz_forms(pm))
        k = len(pm) // 2
        bounds = [(0, 1)] + [(-1, 1)] * k
    else:
        rows = tk_core.dedupe_rows(tk_core.hankel_forms(pm))
        k = len(pm) // 2
        bounds = [(0, 1)] + [(0, 2)] * k
    hs = []
    d = k + 1
    for i, (lo, hi) in enumerate(bounds):
        a = [0] * d; a[i] = 1
        hs.append((tuple(a), -lo))
        a = [0] * d; a[i] = -1
        hs.append((tuple(a), hi))
    for r in rows:
        hs.append((tuple(r), 0))
        hs.append((tuple(-c for c in r), 1))
    # dedupe
    out = []
    seen = set()
    for h in hs:
        if h not in seen:
            seen.add(h)
            out.append(h)
    return out


def check(n2k, kind):
    orbits = tk_core.orbit_reps(
        tk_core.all_pairings(n2k) if kind == "T" else tk_core.parity_pairings(n2k))
    total_exact = Fraction(0)
    total_hull = 0.0
    print(f"=== {kind}{n2k}: {len(orbits)} orbits ===")
    for rep, size in sorted(orbits.items()):
        d = n2k // 2 + 1
        hs = halfspaces_for(rep, kind)
        bv = brute_vertices(hs, d)
        verts, all_hs, masks = exact_volume.dd_vertices(
            [h for h in hs], [(0, 1)] + ([(-1, 1)] if kind == "T" else [(0, 2)]) * (n2k // 2))
        dd_set = set(verts)
        same = (bv == dd_set)
        exact = exact_volume.polytope_volume(verts, masks, len(all_hs), d)
        pts = np.array([[float(c) for c in v] for v in sorted(bv)])
        try:
            hull = ConvexHull(pts, qhull_options="QJ")
            hv = hull.volume
        except Exception as e:
            hv = float("nan")
        diff = abs(float(exact) - hv)
        total_exact += size * exact
        total_hull += size * hv
        flag = "OK" if (same and diff < 1e-6) else "PROBLEM"
        print(f"  size {size:2d}  |V|={len(bv):3d} vertsEqual={same}  exact={float(exact):.6f} hull={hv:.6f} diff={diff:.2e}  {flag}")
        assert same, f"vertex sets differ for {rep}"
        assert diff < 1e-6, f"volume mismatch for {rep}"
    print(f"TOTAL {kind}{n2k}: exact = {total_exact} = {float(total_exact):.6f}; hull-sum = {total_hull:.6f}")


if __name__ == "__main__":
    check(8, "T")
    check(8, "H")
