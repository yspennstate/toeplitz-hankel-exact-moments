"""
Fully independent exact verification of LSD moments by finite-n circuit counting.

For Rademacher entries, E[(1/n) tr (A_n/sqrt n)^{2k}] = count(n)/n^{k+1}, where count(n)
counts circuits (i_0..i_{2k-1}) in [n]^{2k} whose multiset of slot-values (|i_l - i_{l+1}|
for Toeplitz, i_l + i_{l+1} for Hankel) has every value with EVEN multiplicity.

count(n) is an Ehrhart-type quasi-polynomial of degree k+1; its leading coefficient IS the
limit moment m_{2k}. We compute count(n) exactly by meet-in-the-middle (value-multiset
parity tracked by XOR of 128-bit random keys), interpolate the even-n branch with exact
rationals, verify by predicting held-out n, and extract the leading coefficient.

This shares NO code or mathematics with the polytope-volume pipeline.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from fractions import Fraction
import random


def count_even_multiplicity_circuits(n, n2k, kind, key2):
    """Exact count via meet-in-the-middle over the cut {i_0, i_h} with h = n2k//2."""
    rnd = random.Random(12345 + key2)
    if kind == "T":
        keys = {d: rnd.getrandbits(127) for d in range(0, n)}
        def K(a, b):
            return keys[abs(a - b)]
    else:
        keys = {s: rnd.getrandbits(127) for s in range(0, 2 * n - 1)}
        def K(a, b):
            return keys[a + b]
    h = n2k // 2  # slots 0..h-1 use i_0..i_h ; slots h..2k-1 use i_h..i_{2k-1}, i_0
    rng_n = range(n)

    total = 0
    # A-half: paths i_0 -> i_1 -> ... -> i_h : XOR of K over consecutive steps
    # B-half: paths i_h -> i_{h+1} -> ... -> i_{2k-1} -> i_0
    # for each (i_0, i_h): hash maps of XOR over middle indices
    # build A-map: for fixed i_0: dict over i_h of {xor: count}? we need (i_0, i_h) pairs:
    # do DP: paths of length h from i_0 to i_h with xor x: dict[(end, xor)] -> count
    for i0 in rng_n:
        # DP forward over h steps
        cur = {(i0, 0): 1}
        for _ in range(h):
            nxt = {}
            for (a, x), c in cur.items():
                for b in rng_n:
                    key = (b, x ^ K(a, b))
                    nxt[key] = nxt.get(key, 0) + c
            cur = nxt
        amap = cur  # (i_h, xor) -> count
        # B DP backward: paths from i_h (free) through h steps ending back at i0
        # do reverse DP from i0: paths of length h from i0 backwards
        curb = {(i0, 0): 1}
        for _ in range(h):
            nxt = {}
            for (a, x), c in curb.items():
                for b in rng_n:
                    key = (b, x ^ K(b, a))
                    nxt[key] = nxt.get(key, 0) + c
            curb = nxt
        bmap = curb  # (i_h, xor) -> count of B-half paths i_h -> ... -> i0
        # combine: same i_h, xors equal (xor total == 0)
        for (a, x), c in amap.items():
            c2 = bmap.get((a, x))
            if c2:
                total += c * c2
    return total


def lagrange_leading_coeff(xs, ys, deg):
    """Exact leading coefficient of the unique degree<=deg polynomial through points."""
    assert len(xs) == deg + 1
    lead = Fraction(0)
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        denom = 1
        for j, xj in enumerate(xs):
            if j != i:
                denom *= (xi - xj)
        lead += Fraction(yi, denom)
    return lead


def eval_lagrange(xs, ys, x):
    tot = Fraction(0)
    for i, (xi, yi) in enumerate(zip(xs, ys)):
        num = Fraction(yi)
        for j, xj in enumerate(xs):
            if j != i:
                num *= Fraction(x - xj, xi - xj)
        tot += num
    return tot


def run(kind, n2k, fit_ns, test_ns):
    k = n2k // 2
    print(f"=== {kind}{n2k}: counting circuits ===", flush=True)
    counts = {}
    for n in sorted(set(fit_ns + test_ns)):
        c1 = count_even_multiplicity_circuits(n, n2k, kind, key2=0)
        c2 = count_even_multiplicity_circuits(n, n2k, kind, key2=1)
        assert c1 == c2, f"XOR collision suspected at n={n}: {c1} vs {c2}"
        counts[n] = c1
        print(f"  n={n:2d}: count = {c1}   count/n^{k+1} = {c1 / n**(k+1):.5f}", flush=True)
    lead = lagrange_leading_coeff(fit_ns, [counts[n] for n in fit_ns], k + 1)
    print(f"  leading coeff from even-branch interpolation: {lead} = {float(lead):.6f}")
    ok = True
    for n in test_ns:
        pred = eval_lagrange(fit_ns, [counts[x] for x in fit_ns], n)
        match = (pred == counts[n])
        ok = ok and match
        print(f"  predict n={n}: {pred} vs actual {counts[n]}  {'MATCH' if match else 'MISMATCH'}")
    print(f"  => m_{n2k}^{kind} = {lead} ({'verified by held-out prediction' if ok else 'PREDICTION FAILED'})\n", flush=True)
    return lead, ok


if __name__ == "__main__":
    # controls first: T6 (=11, Hammond-Miller & ours), H6 (=11/2 BDJ), H8 (=281/15 BDJ)
    run("T", 6, fit_ns=[4, 6, 8, 10, 12], test_ns=[14, 16])
    run("H", 6, fit_ns=[4, 6, 8, 10, 12], test_ns=[14, 16])
    run("H", 8, fit_ns=[4, 6, 8, 10, 12, 14], test_ns=[16, 18])
    # the disputed one
    run("T", 8, fit_ns=[4, 6, 8, 10, 12, 14], test_ns=[16, 18])
