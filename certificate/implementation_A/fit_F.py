# Fit F/N = c + beta * P(N)^2/N over the exact_F.py output; report residual stats, primorial subset, and max P^2/N
import sys, numpy as np
d = np.loadtxt(sys.argv[1], comments="#")
N, FN, P2N, R, TN, rg, sp = d.T
X = np.vstack([np.ones_like(P2N), P2N]).T
coef, *_ = np.linalg.lstsq(X, FN, rcond=None)
res = FN - X @ coef
print(f"n={len(N)}  fit F/N = {coef[0]:.6f} + {coef[1]:.6f} * P^2/N ; resid sd={res.std():.2e} max|res|={abs(res).max():.2e}")
print(f"-4 log(2 pi) = {-4*np.log(2*np.pi):.6f}")
print(f"(F-P^2)/N: mean={R.mean():.6f} sd={R.std():.2e} min={R.min():.6f} max={R.max():.6f}")
print(f"F/N: sd={FN.std():.3e} min={FN.min():.6f} max={FN.max():.6f};  P^2/N max={P2N.max():.4f}")
s = sp == 1
print(f"primorial-type N (n={s.sum()}): (F-P^2)/N mean={R[s].mean():.6f} sd={R[s].std():.2e}; rG(2N)/2N max={rg[s].max():.3f}")
print(f"corr(F/N, P^2/N)={np.corrcoef(FN,P2N)[0,1]:.6f};  corr((F-P^2)/N, rG(2N)/2N)={np.corrcoef(R,rg)[0,1]:.4f}")
