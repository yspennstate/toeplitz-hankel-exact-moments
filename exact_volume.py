"""
Exact rational polytope volume engine (pure Python, fractions only).

Pipeline per polytope:
  1. Halfspaces (a, b) meaning a.v + b >= 0, integer data.
  2. Incremental double-description from a bounding box: exact vertex enumeration.
  3. Combinatorial face triangulation (anchor coning) using exact tight-sets + rank tests.
  4. Volume = sum |det| / d! over top simplices, exact Fractions.

Built for the small (dim <= 8), highly degenerate 0/1-coefficient polytopes arising from
Toeplitz/Hankel moment computations, but fully generic.
"""
from __future__ import annotations
from fractions import Fraction
from itertools import combinations
from math import factorial

F0 = Fraction(0)


# ----------------------------------------------------------------- lin alg

def mat_rank(rows):
    """Exact rank of a list of tuples of Fractions/ints."""
    m = [list(map(Fraction, r)) for r in rows]
    if not m:
        return 0
    nr, nc = len(m), len(m[0])
    rank = 0
    col = 0
    for col in range(nc):
        piv = None
        for r in range(rank, nr):
            if m[r][col] != 0:
                piv = r
                break
        if piv is None:
            continue
        m[rank], m[piv] = m[piv], m[rank]
        pv = m[rank][col]
        for r in range(nr):
            if r != rank and m[r][col] != 0:
                f = m[r][col] / pv
                mr = m[r]
                mk = m[rank]
                for c in range(col, nc):
                    mr[c] -= f * mk[c]
        rank += 1
        if rank == nr:
            break
    return rank


def det_frac(rows):
    """Exact determinant of a square matrix of Fractions."""
    m = [list(map(Fraction, r)) for r in rows]
    n = len(m)
    det = Fraction(1)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if m[r][c] != 0:
                piv = r
                break
        if piv is None:
            return F0
        if piv != c:
            m[c], m[piv] = m[piv], m[c]
            det = -det
        pv = m[c][c]
        det *= pv
        for r in range(c + 1, n):
            if m[r][c] != 0:
                f = m[r][c] / pv
                mr = m[r]
                mc = m[c]
                for cc in range(c, n):
                    mr[cc] -= f * mc[cc]
    return det


# ----------------------------------------------------- double description

def box_vertices(bounds):
    """bounds: list of (lo, hi) ints per coordinate. Return vertices + box halfspaces.
    Halfspace list: for coord i: (e_i, -lo) [v_i >= lo] and (-e_i, hi) [v_i <= hi]."""
    d = len(bounds)
    halfspaces = []
    for i, (lo, hi) in enumerate(bounds):
        a = [0] * d
        a[i] = 1
        halfspaces.append((tuple(a), -lo))
        a = [0] * d
        a[i] = -1
        halfspaces.append((tuple(a), hi))
    verts = []
    for bits in range(1 << d):
        v = tuple(Fraction(bounds[i][1] if (bits >> i) & 1 else bounds[i][0]) for i in range(d))
        verts.append(v)
    return verts, halfspaces


def hval(h, v):
    a, b = h
    s = Fraction(b)
    for ai, vi in zip(a, v):
        if ai:
            s += ai * vi
    return s


def dd_vertices(extra_halfspaces, bounds):
    """Exact vertex enumeration of {v in box(bounds): a.v + b >= 0 for all extra}.
    Returns (verts, halfspaces, masks): vertices, full halfspace list, tight bitmasks."""
    verts, halfspaces = box_vertices(bounds)
    d = len(bounds)
    # dedupe extra halfspaces, drop any equal to box ones
    seen = set(halfspaces)
    extras = []
    for h in extra_halfspaces:
        h = (tuple(h[0]), h[1])
        if h not in seen:
            seen.add(h)
            extras.append(h)
    all_hs = list(halfspaces)

    def tight_mask(v, hs):
        m = 0
        for idx, h in enumerate(hs):
            if hval(h, v) == 0:
                m |= 1 << idx
        return m

    masks = [tight_mask(v, all_hs) for v in verts]

    for h in extras:
        vals = [hval(h, v) for v in verts]
        pos = [i for i, x in enumerate(vals) if x > 0]
        neg = [i for i, x in enumerate(vals) if x < 0]
        zer = [i for i, x in enumerate(vals) if x == 0]
        hidx = len(all_hs)
        all_hs.append(h)
        if not neg:
            for i in zer:
                masks[i] |= 1 << hidx
            continue
        if not pos and not zer:
            # empty polytope
            return [], all_hs, []
        new_verts = []
        dm1 = d - 1
        for i in pos:
            mi = masks[i]
            vi = vals[i]
            vvi = verts[i]
            for j in neg:
                common = mi & masks[j]
                if common.bit_count() < dm1:
                    continue
                # combinatorial adjacency: no other vertex's tight set contains common
                adjacent = True
                for t in range(len(verts)):
                    if t != i and t != j and (masks[t] & common) == common:
                        adjacent = False
                        break
                if not adjacent:
                    continue
                vj = vals[j]
                vvj = verts[j]
                lam = vi / (vi - vj)  # in (0,1): x = vi*? -> x = v_i + lam (v_j - v_i)
                x = tuple(a + lam * (b - a) for a, b in zip(vvi, vvj))
                new_verts.append(x)
        keep = [verts[i] for i in pos + zer]
        keep_masks = [masks[i] for i in pos + zer]
        # add new vertices, dedupe coordinates
        existing = set(keep)
        for x in new_verts:
            if x not in existing:
                existing.add(x)
                keep.append(x)
                keep_masks.append(0)  # recomputed below
        # recompute masks fresh against all halfspaces so far (safe & simple)
        verts = keep
        masks = [tight_mask(v, all_hs) for v in verts]

    return verts, all_hs, masks


# ----------------------------------------------------------- triangulation

def polytope_volume(verts, masks, n_halfspaces, d):
    """Exact volume of conv(verts) (vertices with tight masks over the halfspaces).
    Returns Fraction. Assumes verts are exactly the polytope's vertices."""
    if len(verts) <= d:
        return F0
    # affine rank check
    v0 = verts[0]
    if mat_rank([tuple(a - b for a, b in zip(v, v0)) for v in verts[1:]]) < d:
        return F0

    all_idx = frozenset(range(len(verts)))
    rank_memo = {}

    def affine_dim(idxset):
        r = rank_memo.get(idxset)
        if r is None:
            it = iter(idxset)
            base = verts[next(it)]
            rows = [tuple(a - b for a, b in zip(verts[i], base)) for i in it]
            r = mat_rank(rows) if rows else 0
            rank_memo[idxset] = r
        return r

    tri_memo = {}

    def triangulate(face, j):
        """face: frozenset of vertex indices with affine dim j. Returns list of
        (j+1)-tuples of indices triangulating it."""
        key = face
        got = tri_memo.get(key)
        if got is not None:
            return got
        if j == 0:
            res = [tuple(face)]
            tri_memo[key] = res
            return res
        if len(face) == j + 1:
            res = [tuple(sorted(face))]
            tri_memo[key] = res
            return res
        anchor = min(face)
        # common tight bits of the whole face
        common = ~0
        for i in face:
            common &= masks[i]
        simplices = []
        seen_facets = set()
        for hbit_idx in range(n_halfspaces):
            bit = 1 << hbit_idx
            if common & bit:
                continue  # tight on whole face -> not a proper subface
            G = frozenset(i for i in face if masks[i] & bit)
            if len(G) < j or G == face or G in seen_facets:
                continue
            if affine_dim(G) != j - 1:
                continue
            seen_facets.add(G)
            if anchor in G:
                continue
            for s in triangulate(G, j - 1):
                simplices.append(s + (anchor,))
        tri_memo[key] = simplices
        return simplices

    total = F0
    for s in triangulate(all_idx, d):
        base = verts[s[0]]
        rows = [tuple(a - b for a, b in zip(verts[i], base)) for i in s[1:]]
        total += abs(det_frac(rows))
    return total / factorial(d)


def volume_from_halfspaces(extra_halfspaces, bounds):
    verts, hs, masks = dd_vertices(extra_halfspaces, bounds)
    if not verts:
        return F0
    return polytope_volume(verts, masks, len(hs), len(bounds))


# ------------------------------------------------------- pairing volumes

def pairing_volume_exact(pm, kind):
    """Exact volume p(pi) for a pairing pm, kind 'T' (Toeplitz) or 'H' (Hankel)."""
    import tk_core
    if kind == "T":
        rows = tk_core.dedupe_rows(tk_core.toeplitz_forms(pm))
        k = len(pm) // 2
        bounds = [(0, 1)] + [(-1, 1)] * k
    else:
        rows = tk_core.dedupe_rows(tk_core.hankel_forms(pm))
        k = len(pm) // 2
        bounds = [(0, 1)] + [(0, 2)] * k
    extra = []
    for r in rows:
        # 0 <= r.v <= 1  -> (r, 0) and (-r, 1); skip rows identical to box bounds
        extra.append((tuple(r), 0))
        extra.append((tuple(-c for c in r), 1))
    return volume_from_halfspaces(extra, bounds)


# ------------------------------------------------------------- self tests

if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    import time
    import tk_core

    t0 = time.time()
    # 1. boxes
    for d in range(2, 7):
        v = volume_from_halfspaces([], [(0, 1)] * d)
        assert v == 1, (d, v)
    # 2. standard simplex x_i >= 0, sum x_i <= 1
    for d in range(2, 7):
        hs = [(tuple(-1 for _ in range(d)), 1)]
        v = volume_from_halfspaces(hs, [(0, 1)] * d)
        assert v == Fraction(1, factorial(d)), (d, v)
    # 3. cross polytope |x|_1 <= 1 in [-1,1]^d: vol = 2^d/d!
    for d in range(2, 6):
        hs = []
        for signs in range(1 << d):
            a = tuple(1 if (signs >> i) & 1 else -1 for i in range(d))
            hs.append((a, 1))
        v = volume_from_halfspaces(hs, [(-1, 1)] * d)
        assert v == Fraction(2 ** d, factorial(d)), (d, v)
    # 4. half-open slice: x0+x1 <= 1/2 in unit square: triangle area 1/8
    v = volume_from_halfspaces([((-2, -2), 1)], [(0, 1)] * 2)
    assert v == Fraction(1, 8), v
    print(f"generic engine tests OK ({time.time()-t0:.1f}s)")

    # 5. Toeplitz known: crossing pairing 2k=4 -> 2/3; m4 = 8/3; m6 = 11 (Hammond-Miller)
    t0 = time.time()
    for n2k, kind, expect in [(4, "T", Fraction(8, 3)), (6, "T", Fraction(11)),
                              (4, "H", Fraction(2)), (6, "H", Fraction(11, 2)),
                              (8, "H", Fraction(281, 15))]:
        if kind == "T":
            orbits = tk_core.orbit_reps(tk_core.all_pairings(n2k))
        else:
            orbits = tk_core.orbit_reps(tk_core.parity_pairings(n2k))
        tot = F0
        for rep, size in orbits.items():
            tot += size * pairing_volume_exact(rep, kind)
        status = "OK" if tot == expect else f"MISMATCH expected {expect}"
        print(f"m_{n2k}^{kind} exact = {tot}   {status}   ({time.time()-t0:.1f}s cum)")
        assert tot == expect, (n2k, kind, tot, expect)
    print("ALL EXACT CHECKS PASSED")
