"""Emit the exact volume spectrum of a results file as a LaTeX-ready line + checks."""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import json
from collections import Counter
from fractions import Fraction

kind, n2k = sys.argv[1], int(sys.argv[2])
with open(f"results/volumes_{kind}{n2k}.json") as f:
    data = json.load(f)
spec = Counter()
tot = Fraction(0)
npair = 0
for k, rec in data["orbits"].items():
    v = Fraction(rec["volume"])
    spec[v] += rec["size"]
    tot += v * rec["size"]
    npair += rec["size"]
print(f"{kind}{n2k}: {npair} pairings, {len(data['orbits'])} orbits, {len(spec)} distinct volumes")
print(f"moment = {tot} = {float(tot):.6f}")
items = sorted(spec.items(), key=lambda kv: -kv[0])
line = ";\\ ".join(
    (f"$1\\times{c}$" if v == 1 else f"$\\tfrac{{{v.numerator}}}{{{v.denominator}}}\\times{c}$")
    for v, c in items)
print(line)
