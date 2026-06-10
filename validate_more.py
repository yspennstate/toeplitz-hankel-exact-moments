"""Extended validation on the production results:
(a) orbit-constancy (Lemma 1) spot checks: random orbit members recomputed and compared;
(b) non-crossing <=> volume 1 check across all computed orbits;
(c) summary tables of the volume spectra per order (for the paper).
Usage: python validate_more.py T 10 [n_orbit_samples]
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import json
import random
from fractions import Fraction
import tk_core
import exact_volume

kind = sys.argv[1]
n2k = int(sys.argv[2])
nsamp = int(sys.argv[3]) if len(sys.argv) > 3 else 12

with open(f"results/volumes_{kind}{n2k}.json") as f:
    data = json.load(f)
orbits = {tuple(json.loads(k)): rec for k, rec in data["orbits"].items()}

# (a) orbit constancy on random members
rng = random.Random(99)
reps = list(orbits.keys())
rng.shuffle(reps)
bad = 0
for rep in reps[:nsamp]:
    imgs = list(set(tk_core.dihedral_images(rep)))
    pm = rng.choice(imgs)
    v = exact_volume.pairing_volume_exact(pm, kind)
    v0 = Fraction(orbits[rep]["volume"])
    ok = (v == v0)
    bad += (not ok)
    print(f"orbit-constancy: rep vol {v0} random-member vol {v}  {'OK' if ok else 'FAIL'}")
assert bad == 0

# (b) non-crossing <=> volume 1
nc_ok = True
n_nc_pairings = 0
n_vol1_pairings = 0
for rep, rec in orbits.items():
    v = Fraction(rec["volume"])
    cross = tk_core.is_crossing(rep)
    if not cross:
        n_nc_pairings += rec["size"]
        if v != 1:
            nc_ok = False
            print("NONCROSSING WITH VOL != 1:", rep, v)
    if v == 1:
        n_vol1_pairings += rec["size"]
        if cross:
            nc_ok = False
            print("CROSSING WITH VOL == 1:", rep, v)
print(f"non-crossing pairings: {n_nc_pairings}; volume-1 pairings: {n_vol1_pairings}; "
      f"equivalence {'HOLDS' if nc_ok and n_nc_pairings == n_vol1_pairings else 'FAILS'}")

# (c) volume spectrum
from collections import Counter
spec = Counter()
for rep, rec in orbits.items():
    spec[Fraction(rec["volume"])] += rec["size"]
tot = Fraction(0)
print("\nvolume spectrum (volume: #pairings):")
for v in sorted(spec, reverse=True):
    print(f"  {str(v):>10}  x {spec[v]}")
    tot += v * spec[v]
print(f"moment m_{n2k}^{kind} = {tot} = {float(tot):.6f}")
