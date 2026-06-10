"""Figure: simulated spectral densities of gamma_T, gamma_H with candidate overlays;
plus a crude simulation check of the predicted Hankel norm constant ~0.57."""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
import numpy as np
from math import gamma, sqrt, log

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

rng = np.random.default_rng(31)

def eigs(kind, n, reps):
    idx = np.arange(n)
    M = np.abs(np.subtract.outer(idx, idx)) if kind == "T" else np.add.outer(idx, idx)
    out = []
    for _ in range(reps):
        a = rng.standard_normal(2 * n)
        out.append(np.linalg.eigvalsh(a[M] / np.sqrt(n)))
    return np.concatenate(out)

n, reps = 1500, 60
evT = eigs("T", n, reps)
evH = eigs("H", n, reps)

def ged_pdf(x, beta):
    a = sqrt(gamma(1.0 / beta) / gamma(3.0 / beta))  # variance 1
    c = beta / (2 * a * gamma(1.0 / beta))
    return c * np.exp(-np.abs(x / a) ** beta)

xs = np.linspace(-4, 4, 801)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
axes[0].hist(evT, bins=160, density=True, alpha=0.45, color="steelblue", label="Toeplitz eigenvalues ($n{=}1500$)")
axes[0].plot(xs, ged_pdf(xs, 2.4367), "r-", lw=1.6, label=r"GED $\beta=2.4367$ (certified $\leq0.4\%$ to $m_{10}$)")
axes[0].plot(xs, np.exp(-xs**2/2)/np.sqrt(2*np.pi), "k--", lw=1.2, label=r"$\mathcal{N}(0,1)$")
axes[0].set_title(r"$\gamma_T$"); axes[0].legend(fontsize=8); axes[0].set_xlim(-4, 4)
axes[1].hist(evH, bins=160, density=True, alpha=0.45, color="seagreen", label="Hankel eigenvalues ($n{=}1500$)")
axes[1].plot(xs, np.abs(xs)*np.exp(-xs**2), "r-", lw=1.6, label=r"reverse-circulant $|x|e^{-x^2}$")
axes[1].plot(xs, ged_pdf(xs, 6.0), "k--", lw=1.2, label=r"GED $\beta=6$ (excluded)")
axes[1].set_title(r"$\gamma_H$"); axes[1].legend(fontsize=8); axes[1].set_xlim(-4, 4)
fig.tight_layout()
fig.savefig("paper/fig_density.pdf"); fig.savefig("paper/fig_density.png", dpi=150)
print("figure saved")

# empirical check: ||H_n|| / sqrt(2 n log n)
print("\nHankel norm-constant check (prediction ~0.56-0.57):")
for nn, rr in [(1000, 12), (2000, 8), (4000, 4), (8000, 2)]:
    idx = np.arange(nn)
    S = np.add.outer(idx, idx)
    vals = []
    for _ in range(rr):
        a = rng.standard_normal(2 * nn)
        H = a[S] / np.sqrt(nn)
        vals.append(np.abs(np.linalg.eigvalsh(H)).max())
    c = np.mean(vals) / sqrt(2 * log(nn))
    print(f"  n={nn:5d}: ||H||/sqrt(2 n log n) = {c:.4f} +- {np.std(vals)/sqrt(2*log(nn))/sqrt(rr):.4f}")
print("(Toeplitz comparison, Sen-Virag 0.8289):")
for nn, rr in [(2000, 6), (4000, 3)]:
    idx = np.arange(nn)
    D = np.abs(np.subtract.outer(idx, idx))
    vals = []
    for _ in range(rr):
        a = rng.standard_normal(2 * nn)
        T = a[D] / np.sqrt(nn)
        vals.append(np.abs(np.linalg.eigvalsh(T)).max())
    c = np.mean(vals) / sqrt(2 * log(nn))
    print(f"  n={nn:5d}: ||T||/sqrt(2 n log n) = {c:.4f}")
