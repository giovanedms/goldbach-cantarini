import numpy as np, math, sys
d = np.loadtxt(sys.argv[1], skiprows=1)
N, FN, P2N, TN, rg = d[:,0], d[:,1], d[:,2], d[:,3], d[:,4]
X = np.c_[np.ones_like(N), P2N]
b, *_ = np.linalg.lstsq(X, FN, rcond=None); res = FN - X@b
print(f"n={len(N)} N in [{N.min():.0f},{N.max():.0f}]")
print(f"fit F/N = {b[0]:.6f} + {b[1]:.6f}*P^2/N   resid sd={res.std():.2e} max={abs(res).max():.2e}")
print(f"-4 log(2pi) = {-4*math.log(2*math.pi):.6f}")
G = FN - P2N
print(f"(F-P^2)/N: mean={G.mean():.6f} sd={G.std():.2e} range=[{G.min():.6f},{G.max():.6f}]")
print(f"F/N range=[{FN.min():.4f},{FN.max():.4f}]  P^2/N max={P2N.max():.4f}  T/N range=[{TN.min():.3f},{TN.max():.3f}]")
print(f"corr(F-P^2, rG/2N)={np.corrcoef(G, rg)[0,1]:.3f}  corr(F-P^2, P^2/N)={np.corrcoef(G,P2N)[0,1]:.3f}")
# alternative hypotheses: is F + k P^2 constant?
for k in (-1, 0, 1, 2):
    H = FN + (k-1)*P2N
    print(f"  sd of (F/N {'+' if k-1>=0 else '-'} {abs(k-1)}*P^2/N) = {H.std():.3e}")
hi = rg > 3
print(f"2N highly composite (rG/2N>3, n={hi.sum()}): (F-P^2)/N mean={G[hi].mean():.6f} sd={G[hi].std():.2e}")
