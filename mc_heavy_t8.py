"""10^8-sample MC per T8 orbit: second independent numeric pillar for m8_T."""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from fractions import Fraction
import json
import tk_core

orb = tk_core.orbit_reps(tk_core.all_pairings(8))
with open("results/volumes_T8.json") as f:
    saved = json.load(f)
ex = {tuple(json.loads(k)): Fraction(rec["volume"]) for k, rec in saved["orbits"].items()}

tot_mc = 0.0
tot_var = 0.0
N = 100_000_000
print("orbit-by-orbit: exact vs heavy MC (1e8 samples)")
for rep, size in sorted(orb.items()):
    rows = tk_core.dedupe_rows(tk_core.toeplitz_forms(rep))
    mc = tk_core.mc_volume(rows, "T", n_samples=N, seed=hash(rep) % (2**31), batch=2_000_000)
    exv = float(ex[rep])
    boxvol = 2.0 ** 4
    p_acc = exv / boxvol
    sigma = boxvol * (p_acc * (1 - p_acc) / N) ** 0.5
    z = (mc - exv) / sigma if sigma > 0 else 0.0
    tot_mc += size * mc
    tot_var += (size * sigma) ** 2
    print(f"  size {size:2d}  exact {exv:8.5f}  mc {mc:8.5f}  z={z:+5.2f}")
tot_sig = tot_var ** 0.5
print(f"\nTOTAL m8_T: MC = {tot_mc:.5f} +- {tot_sig:.5f}")
print(f"            exact mine = {float(sum(ex[r]*s for r, s in orb.items())):.5f} (=908/15)")
print(f"            HM published = 964/15 = {964/15:.5f}  ->  z = {(tot_mc - 964/15)/tot_sig:+.1f} sigma")
