# Validates (a) A(x) (Ingham, Thm 28) against the sum over zeros; (b) T(N) in Stieltjes form against the truncated double series
# T_K = 2 sum_{rho1} int_2^N P_K(t) (2N-t)^{rho1-1} dt  (identity checked in beta_identity.py), K zeros (RH assumed for these zeros).
# Uniform grid in u=log t (frequencies <= gamma_max), midpoint rule.
import numpy as np, sys
from core import Ctx
g = np.loadtxt(sys.argv[1]); rho_all = 0.5+1j*g
c = Ctx(1000)
for x in (10, 100, 500):
    za = 2*np.real(np.sum(x**(rho_all+1)/(rho_all*(rho_all+1))))
    print(f"A({x}) Ingham={c.A(x):.6f} zeros(1e5)={za:.6f}", flush=True)
for N in [int(a) for a in sys.argv[2].split(',')]:
    Ts = c.T_stieltjes(N)
    P2 = float(c.P(np.array([N]))[0])**2
    print(f"N={N} T_stieltjes={Ts:.4f}  P(N)^2={P2:.4f}  T_stieltjes-2P^2={Ts-2*P2:.4f} T_stieltjes+2P^2={Ts+2*P2:.4f}", flush=True)
    for K in [int(a) for a in sys.argv[3].split(',')]:
        rho = rho_all[:K]; gm = g[K-1]
        du = 0.2/gm
        u = np.arange(np.log(2), np.log(N), du) + du/2
        t = np.exp(u); l2 = np.log(2*N-t)
        PK = np.zeros(len(t)); KK = np.zeros(len(t))
        for i in range(0, K, 25):
            r = rho[i:i+25, None]
            PK += 2*np.real(np.sum(np.exp(r*u)/r, axis=0))
            KK += 2*np.real(np.sum(np.exp((r-1)*l2), axis=0))
        TK = 2*np.sum(PK*KK*t)*du
        print(f"   K={K} gamma_max={gm:.0f} T_K={TK:.4f}  diff_vs_stieltjes={TK-Ts:+.4f}", flush=True)
