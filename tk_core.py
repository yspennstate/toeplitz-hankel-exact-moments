"""
Core combinatorics for exact moments of the Toeplitz / Hankel limiting spectral
distributions.

Moment-volume formulas (see NOTES_derivation.md):

Toeplitz:  m_{2k}^T = sum over all pairings pi of [0..2k-1] of vol(P_pi),
  P_pi = { (x0, t_1..t_k) : 0 <= x0 + S_l(t) <= 1, l = 0..2k-1 },
  S_l = sum_{r<l} sigma(r) t_{v(r)},  sigma=+1 at the first slot of a pair, -1 at the second.

Hankel:    m_{2k}^H = sum over parity-mixing pairings (even<->odd slots only) of vol(P^H_pi),
  P^H_pi = { (x0, y_1..y_k) : 0 <= A_l <= 1, l = 0..2k-1 },
  A_l = (-1)^l x0 + sum_{r<l} (-1)^(l-1-r) y_{v(r)}.

A pairing on 2k slots is stored as a tuple pm of length 2k with pm[i] = partner of i.
"""
from __future__ import annotations
import itertools
from fractions import Fraction


# ---------------------------------------------------------------- pairings

def all_pairings(n2k: int):
    """Yield all perfect matchings of range(n2k) as partner tuples."""
    slots = list(range(n2k))

    def rec(avail):
        if not avail:
            yield ()
            return
        a = avail[0]
        for jx in range(1, len(avail)):
            b = avail[jx]
            rest = avail[1:jx] + avail[jx + 1:]
            for sub in rec(rest):
                yield ((a, b),) + sub

    for pairs in rec(slots):
        pm = [0] * n2k
        for a, b in pairs:
            pm[a] = b
            pm[b] = a
        yield tuple(pm)


def is_parity_mixing(pm) -> bool:
    return all((i + pm[i]) % 2 == 1 for i in range(len(pm)))


def parity_pairings(n2k: int):
    """All pairings joining even slots to odd slots (k! of them)."""
    evens = list(range(0, n2k, 2))
    odds = list(range(1, n2k, 2))
    for perm in itertools.permutations(odds):
        pm = [0] * n2k
        for e, o in zip(evens, perm):
            pm[e] = o
            pm[o] = e
        yield tuple(pm)


def is_crossing(pm) -> bool:
    n = len(pm)
    pairs = [(i, pm[i]) for i in range(n) if i < pm[i]]
    for (a, b), (c, d) in itertools.combinations(pairs, 2):
        if a < c < b < d or c < a < d < b:
            return True
    return False


# ------------------------------------------------- dihedral canonicalization

def dihedral_images(pm):
    """All images of the pairing under the dihedral group D_{2k} (order 4k)."""
    n = len(pm)
    out = []
    for s in range(n):
        # rotation i -> (i+s) % n
        im = [0] * n
        for i in range(n):
            im[(i + s) % n] = (pm[i] + s) % n
        out.append(tuple(im))
        # reflection i -> (s - i) % n
        im2 = [0] * n
        for i in range(n):
            im2[(s - i) % n] = (s - pm[i]) % n
        out.append(tuple(im2))
    return out


def canonical(pm):
    return min(dihedral_images(pm))


def orbit_reps(pairings_iter):
    """Group pairings into dihedral orbits. Returns dict canonical_rep -> orbit size."""
    orbits = {}
    for pm in pairings_iter:
        orbits[canonical(pm)] = orbits.get(canonical(pm), 0) + 1
    return orbits


# ------------------------------------------------------------- constraints

def toeplitz_forms(pm):
    """Integer matrix F (rows l=0..2k-1) of forms x0 + S_l over vars (x0, t_1..t_k).
    Constraint: 0 <= F . v <= 1 for each row. Rows deduped preserve full system."""
    n = len(pm)
    k = n // 2
    # assign pair indices and sigma
    pair_idx = {}
    sigma = [0] * n
    nxt = 0
    for i in range(n):
        j = pm[i]
        if i < j:
            pair_idx[i] = nxt
            pair_idx[j] = nxt
            sigma[i] = 1
            sigma[j] = -1
            nxt += 1
    rows = []
    cur = [0] * (k + 1)
    cur[0] = 1  # x0 coefficient
    rows.append(tuple(cur))  # l = 0: x0
    for l in range(1, n):
        r = l - 1
        cur = list(rows[-1])
        cur[1 + pair_idx[r]] += sigma[r]
        rows.append(tuple(cur))
    return rows


def hankel_forms(pm):
    """Integer matrix rows A_l over vars (x0, y_1..y_k); constraints 0 <= A_l <= 1."""
    n = len(pm)
    k = n // 2
    pair_idx = {}
    nxt = 0
    for i in range(n):
        j = pm[i]
        if i < j:
            pair_idx[i] = nxt
            pair_idx[j] = nxt
            nxt += 1
    rows = []
    prev = [0] * (k + 1)
    prev[0] = 1  # A_0 = x0
    rows.append(tuple(prev))
    for l in range(1, n):
        # A_l = s_{l-1} - A_{l-1} where s_{l-1} = y_{pair_idx[l-1]}
        cur = [-c for c in rows[-1]]
        cur[1 + pair_idx[l - 1]] += 1
        rows.append(tuple(cur))
    return rows


def dedupe_rows(rows):
    seen = []
    out = []
    for r in rows:
        if r not in seen:
            seen.append(r)
            out.append(r)
    return out


# ------------------------------------------------------------- Monte Carlo

def mc_volume(rows, kind: str, n_samples: int = 2_000_000, seed: int = 0, batch: int = 500_000):
    """Monte Carlo estimate of vol{v in box : 0 <= F v <= 1}.
    kind 'T': box = [0,1] x [-1,1]^k (vol 2^k); kind 'H': [0,1] x [0,2]^k (vol 2^k)."""
    import numpy as np
    F = np.array(rows, dtype=np.float64)
    k = F.shape[1] - 1
    rng = np.random.default_rng(seed)
    hits = 0
    tot = 0
    while tot < n_samples:
        b = min(batch, n_samples - tot)
        x0 = rng.random(b)
        if kind == "T":
            t = rng.random((b, k)) * 2.0 - 1.0
        else:
            t = rng.random((b, k)) * 2.0
        V = np.column_stack([x0, t])
        W = V @ F.T
        ok = ((W >= 0.0) & (W <= 1.0)).all(axis=1)
        hits += int(ok.sum())
        tot += b
    boxvol = 2.0 ** k
    return hits / tot * boxvol


# ---------------------------------------------------------------- summaries

def moment_mc(n2k: int, kind: str, n_samples=400_000, seed=1):
    """Quick MC estimate of m_{2k} using orbit reduction."""
    if kind == "T":
        orbits = orbit_reps(all_pairings(n2k))
        formsf = toeplitz_forms
    else:
        orbits = orbit_reps(parity_pairings(n2k))
        formsf = hankel_forms
    total = 0.0
    for rep, size in orbits.items():
        rows = dedupe_rows(formsf(rep))
        total += size * mc_volume(rows, kind, n_samples=n_samples, seed=seed)
    return total, len(orbits)


if __name__ == "__main__":
    import sys
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    from math import factorial

    def dfact(n):  # (n)!! for odd n
        r = 1
        while n > 1:
            r *= n
            n -= 2
        return r

    # sanity of counts
    for k2 in (2, 4, 6, 8):
        ps = list(all_pairings(k2))
        assert len(ps) == dfact(k2 - 1), (k2, len(ps))
        pp = [p for p in ps if is_parity_mixing(p)]
        assert len(pp) == factorial(k2 // 2)
        assert set(pp) == set(parity_pairings(k2))
    print("pairing counts OK")

    # quick MC moments: expect m2=1, m4T=8/3=2.667, m4H=2, m6H=5.5, m8H=18.733=281/15
    for kind, k2, expect in [("T", 2, 1.0), ("T", 4, 8 / 3), ("H", 4, 2.0),
                             ("H", 6, 11 / 2), ("T", 6, None), ("H", 8, 281 / 15)]:
        val, norb = moment_mc(k2, kind, n_samples=300_000)
        print(f"m_{k2}^{kind} MC ~ {val:.4f}   (orbits {norb})"
              + (f"   expected {expect:.4f}" if expect else "   [NEW]"))
