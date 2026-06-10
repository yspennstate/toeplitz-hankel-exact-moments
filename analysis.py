"""Quantitative consequences of the exact moments:
(1) deficit sequences + ratio extrapolations vs the Sen-Virag constant;
(2) exclusion of candidate closed forms (Gaussian scale mixtures; generalized error
    distributions; symmetrized gamma) using the new exact moments.
"""
import sys
sys.stdout.reconfigure(encoding="utf-8", errors="replace")
from fractions import Fraction as F
import numpy as np
from math import factorial, gamma, sqrt

def dfact(n):
    r = 1
    while n > 1:
        r *= n; n -= 2
    return r

mT = {1: F(1), 2: F(8,3), 3: F(11), 4: F(908,15), 5: F(415), 6: F(23840,7), 7: F(325719,10)}
mH = {1: F(1), 2: F(2), 3: F(11,2), 4: F(281,15), 5: F(2717,36), 6: F(1052,3), 7: F(1331087,720)}

print("== deficits and ratios ==")
print("Toeplitz delta_k = m2k/(2k-1)!!:")
dT = {k: mT[k] / dfact(2*k-1) for k in mT}
for k in sorted(dT): print(f"  k={k}: {dT[k]} = {float(dT[k]):.6f}")
rT = {k: float(dT[k]/dT[k-1]) for k in sorted(dT) if k-1 in dT}
print("ratios:", {k: round(v,5) for k,v in rT.items()})

print("Hankel delta_k = m2k/k!:")
dH = {k: mH[k] / factorial(k) for k in mH}
for k in sorted(dH): print(f"  k={k}: {dH[k]} = {float(dH[k]):.6f}")
rH = {k: float(dH[k]/dH[k-1]) for k in sorted(dH) if k-1 in dH}
print("ratios (-> theta):", {k: round(v,5) for k,v in rH.items()})
print("Hankel Gaussian-base ratios m2k/((2k-1) m2k-2):",
      {k: round(float(mH[k]/( (2*k-1)*mH[k-1])),5) for k in sorted(mH) if k-1 in mH})

def extrap(ks, vs, name):
    # fit v = a + b/k (last two) and a + b/k + c/k^2 (last three)
    k1,k2 = ks[-2],ks[-1]; v1,v2 = vs[-2],vs[-1]
    b = (v1-v2)/(1/k1-1/k2); a_lin = v2 - b/k2
    k0,v0 = ks[-3],vs[-3]
    A = np.array([[1,1/k0,1/k0**2],[1,1/k1,1/k1**2],[1,1/k2,1/k2**2]])
    a_quad = np.linalg.solve(A, np.array([v0,v1,v2]))[0]
    print(f"  {name}: linear-in-1/k -> {a_lin:.4f};  quadratic -> {a_quad:.4f}")
    return a_lin, a_quad

print("\n== extrapolations ==")
ks = sorted(rT); vs = [rT[k] for k in ks]
extrap(ks, vs, "Toeplitz deficit-ratio (-> sigma*^2; Sen-Virag^2 = 0.687074)")
ks = sorted(rH); vs = [rH[k] for k in ks]
lH = extrap(ks, vs, "Hankel k!-deficit-ratio (-> theta)")
th = lH[1]
print(f"  => Hankel tail ~ |x| exp(-x^2/theta), theta ~ {th:.3f}; "
      f"predicted ||H||/sqrt(2 n log n) -> sqrt(theta/2) ~ {sqrt(th/2):.3f}")

print("\n== candidate exclusions ==")
# (1) Gaussian scale mixtures: m4 >= 3 m2^2 always; Toeplitz 8/3 < 3, Hankel 2 < 3.
print("Gaussian scale mixture X = W*G: requires m4 >= 3*m2^2 = 3; "
      f"m4(T) = {float(mT[2]):.4f}, m4(H) = {float(mH[2]):.4f}  -> both EXCLUDED")

# (2) Generalized error distribution f ~ exp(-|x/alpha|^beta):
# m_{2k} = alpha^{2k} * Gamma((2k+1)/beta) / Gamma(1/beta). Fit m2=1, m4, predict m6, m8.
from scipy.optimize import brentq
def ged_fit(m4_target):
    def m2k(alpha, beta, k):
        return alpha**(2*k) * gamma((2*k+1)/beta) / gamma(1/beta)
    def eq(beta):
        # alpha fixed by m2 = 1: alpha^2 = Gamma(1/beta)/Gamma(3/beta)
        a2 = gamma(1/beta)/gamma(3/beta)
        return a2**2 * gamma(5/beta)/gamma(1/beta) - m4_target
    beta = brentq(eq, 0.7, 60.0)
    a2 = gamma(1/beta)/gamma(3/beta)
    pred = {k: a2**k * gamma((2*k+1)/beta)/gamma(1/beta) for k in (3,4,5,6)}
    return beta, pred

for name, m4t, m_true in (("Toeplitz", float(mT[2]), mT), ("Hankel", float(mH[2]), mH)):
    beta, pred = ged_fit(m4t)
    print(f"GED fit to {name} (m2=1, m4={m4t:.4f}): beta = {beta:.4f}")
    for k in (3,4,5,6):
        print(f"   m{2*k}: GED {pred[k]:.4f} vs exact {float(m_true[k]):.4f}  "
              f"rel.err {abs(pred[k]-float(m_true[k]))/float(m_true[k])*100:.1f}%")

# (3) symmetrized gamma: f ~ |x|^(a-1) exp(-|x|/s) symmetrized:
# m_{2k} = s^{2k} * Gamma(a+2k)/Gamma(a). Fit m2=1, m4; predict m6, m8.
def symgamma_fit(m4_target):
    def eq(a):
        s2 = 1.0/((a)*(a+1))
        return s2**2 * (a)*(a+1)*(a+2)*(a+3) - m4_target
    a = brentq(eq, 1e-6, 1e6)
    s2 = 1.0/(a*(a+1))
    pred = {}
    for k in (3,4,5,6):
        prod = 1.0
        for j in range(2*k):
            prod *= (a+j)
        pred[k] = s2**k * prod
    return a, pred

for name, m4t, m_true in (("Toeplitz", float(mT[2]), mT), ("Hankel", float(mH[2]), mH)):
    a, pred = symgamma_fit(m4t)
    print(f"sym-gamma fit to {name}: a = {a:.4f}")
    for k in (3,4,5,6):
        print(f"   m{2*k}: symGam {pred[k]:.4f} vs exact {float(m_true[k]):.4f}  "
              f"rel.err {abs(pred[k]-float(m_true[k]))/float(m_true[k])*100:.1f}%")

# tail bounds
print("\n== tail bounds: gamma([x,inf)) <= m_{2k}/(2 x^{2k}) ==")
for x in (2.5, 3.0, 4.0):
    bT = min(float(mT[k])/(2*x**(2*k)) for k in mT)
    bH = min(float(mH[k])/(2*x**(2*k)) for k in mH)
    print(f"  x={x}: gamma_T <= {bT:.3e}, gamma_H <= {bH:.3e}")
