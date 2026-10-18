# Check the outer zero-sum: sum_{rho, |gamma|<=G} int_2^N f(t)(2N-t)^(rho-1) dt  ->  int_2^N f(t) dP(2N-t)-form
# with f(t) = t - psi(t) - c0 (piecewise linear, jumps at prime powers), P'(s) = 1 - 1/(s^3-s) - sum_m Lambda(m) delta(s-m)
# usage: python3 -I outer_sum.py ZEROS_FILE
import sys, numpy as np
from scipy.integrate import quad
zeros = np.loadtxt(sys.argv[1])
c0 = np.log(2*np.pi)
def lam_table(M):
    L = np.zeros(M+1)
    isp = np.ones(M+1, bool); isp[:2] = False
    for p in range(2, M+1):
        if isp[p]:
            isp[p*p::p] = False
            q = p
            while q <= M: L[q] = np.log(p); q *= p
    return L
for N in (30, 100, 300):
    L = lam_table(2*N); psi = np.cumsum(L)
    f_int = lambda k: k - psi[k] - c0      # f on (k,k+1) is t - psi(k) - c0
    # RHS (distributional limit)
    rhs = 0.0
    for k in range(2, N):
        ck = psi[k] + c0
        rhs += quad(lambda t: (t-ck)*(1-1/((2*N-t)**3-(2*N-t))), k, k+1)[0]
    for m in range(N, 2*N-1):
        if L[m] == 0: continue
        t = 2*N-m
        if m == N:   val, wgt = N - psi[N-1] - c0, 0.5
        elif t == 2: val, wgt = 2 - psi[2] - c0, 0.5
        else:        val, wgt = t - (psi[t]+psi[t-1])/2 - c0, 1.0
        rhs -= wgt*L[m]*val
    # LHS truncated at G zeros (rho and conj): int_{k}^{k+1} (t-ck)(2N-t)^(rho-1) dt analytic
    ks = np.arange(2, N); ck = psi[ks] + c0
    a = 2*N - ks - 1.0; b = 2*N - ks + 0.0      # s-range for t in (k,k+1)
    for G in (100, 1000, 10000, 100000):
        g = zeros[:G]; rho = 0.5 + 1j*g
        tot = 0.0
        for chunk in np.array_split(rho, max(1, G//2000)):
            r = chunk[:, None]
            # int_a^b (2N - s - ck) s^(r-1) ds
            I = (2*N - ck)[None, :]*(b**r - a**r)/r - (b**(r+1) - a**(r+1))/(r+1)
            tot += 2*np.real(I.sum())
        print(f"N={N} G={G:6d} lhs={tot:.6f} rhs={rhs:.6f} diff={tot-rhs:.3e}", flush=True)
