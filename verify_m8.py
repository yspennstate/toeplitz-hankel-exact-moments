"""Resolve the m8 Toeplitz discrepancy vs Hammond-Miller.
1. Brute-force exact sum over ALL 105 pairings (no orbit reduction).
2. Direct simulation: empirical 8th moment of large random Toeplitz matrices,
   with Richardson extrapolation in 1/n.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from fractions import Fraction
import numpy as np
import tk_core
import exact_volume

# ---- 1. brute force exact, all pairings individually
tot = Fraction(0)
vols = {}
for pm in tk_core.all_pairings(8):
    v = exact_volume.pairing_volume_exact(pm, "T")
    vols[pm] = v
    tot += v
print("brute-force exact m8_T over all 105 pairings:", tot, "=", float(tot))

# orbit-reduced again for comparison
orb = tk_core.orbit_reps(tk_core.all_pairings(8))
tot2 = sum(exact_volume.pairing_volume_exact(r, "T") * s for r, s in orb.items())
print("orbit-reduced m8_T:", tot2, "=", float(tot2))

# within-orbit constancy check (validates Lemma 1 computationally)
bad = 0
for rep in orb:
    imgs = set(tk_core.dihedral_images(rep))
    vals = {vols[i] for i in imgs}
    if len(vals) != 1:
        bad += 1
        print("ORBIT NOT CONSTANT:", rep, vals)
print("orbits with non-constant volume:", bad)

# ---- 2. matrix simulation
rng = np.random.default_rng(7)
print("\nmatrix-simulation check of m8_T (and m6, m4):")
for n in (500, 1000, 2000):
    reps = max(6, 24000 // n)
    m8s, m6s, m4s = [], [], []
    for _ in range(reps):
        a = rng.standard_normal(2 * n)
        c = a[:n]
        T = c[np.abs(np.subtract.outer(np.arange(n), np.arange(n)))] / np.sqrt(n)
        ev = np.linalg.eigvalsh(T)
        m8s.append(np.mean(ev ** 8)); m6s.append(np.mean(ev ** 6)); m4s.append(np.mean(ev ** 4))
    print(f" n={n:5d} reps={reps:3d}  m4={np.mean(m4s):8.4f}  m6={np.mean(m6s):8.4f}  m8={np.mean(m8s):9.4f}")
print("targets: m4=8/3=2.6667, m6=11, m8: mine=908/15=60.5333 vs HM-text 644/15=42.93 or 64+4/15=64.2667")
