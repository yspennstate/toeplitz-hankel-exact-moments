"""Heavy matrix simulation: empirical LSD moments of random Toeplitz & Hankel matrices
with SE bars and 1/n extrapolation. Verifies exact moment computations end-to-end."""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np

rng = np.random.default_rng(20260610)
KS = [4, 6, 8, 10, 12]

def run(n, reps):
    msT = {k: [] for k in KS}
    msH = {k: [] for k in KS}
    idx = np.arange(n)
    D = np.abs(np.subtract.outer(idx, idx))
    S = np.add.outer(idx, idx)
    for _ in range(reps):
        a = rng.standard_normal(2 * n) ; aH = rng.standard_normal(2 * n)
        T = a[D] / np.sqrt(n)
        H = aH[S] / np.sqrt(n)
        evT = np.linalg.eigvalsh(T)
        evH = np.linalg.eigvalsh(H)
        for k in KS:
            msT[k].append(np.mean(evT ** k))
            msH[k].append(np.mean(evH ** k))
    out = {}
    for k in KS:
        for kind, ms in (("T", msT), ("H", msH)):
            arr = np.array(ms[k])
            out[(kind, k)] = (arr.mean(), arr.std(ddof=1) / np.sqrt(len(arr)))
    return out

print("n, reps, then  kind k: mean +- SE")
results = {}
for n, reps in [(1000, 320), (2000, 96), (4000, 28)]:
    res = run(n, reps)
    results[n] = res
    line = f"n={n} reps={reps}\n"
    for (kind, k), (m, se) in sorted(res.items()):
        line += f"  {kind}{k:2d}: {m:10.4f} +- {se:.4f}\n"
    print(line, flush=True)

# 1/n extrapolation: fit m(n) = m_inf + c/n by weighted least squares over the 3 n's
print("\n1/n-extrapolated limits (m_inf +- SE):")
ns = np.array(sorted(results.keys()), dtype=float)
for kind in ("T", "H"):
    for k in KS:
        y = np.array([results[int(n)][(kind, k)][0] for n in ns])
        se = np.array([results[int(n)][(kind, k)][1] for n in ns])
        X = np.column_stack([np.ones_like(ns), 1.0 / ns])
        W = np.diag(1.0 / se ** 2)
        beta, cov = np.linalg.solve(X.T @ W @ X, X.T @ W @ y), np.linalg.inv(X.T @ W @ X)
        print(f"  {kind}{k:2d}: m_inf = {beta[0]:10.4f} +- {np.sqrt(cov[0,0]):.4f}")
