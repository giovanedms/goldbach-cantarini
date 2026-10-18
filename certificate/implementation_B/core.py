# Exact core (own sieve). c0 = zeta'/zeta(0) = log(2 pi); z1 = zeta'/zeta(-1).
import numpy as np, math, mpmath as mp
C0 = math.log(2*math.pi)
Z1 = float(mp.diff(mp.zeta,-1)/mp.zeta(-1))
def lam_table(M):
    lam = np.zeros(M+1)
    s = np.ones(M+1, bool); s[:2] = False
    for p in range(2, int(M**.5)+1):
        if s[p]: s[p*p::p] = False
    for p in np.nonzero(s)[0]:
        lp = math.log(p); q = int(p)
        while q <= M: lam[q] = lp; q *= p
    return lam
class Ctx:
    def __init__(s, M):
        s.lam = lam_table(M)
        s.psi = np.cumsum(s.lam)                       # psi(n), n integer
        n = np.arange(M+1, dtype=float)
        s.nlam = np.cumsum(n*s.lam)                    # sum_{k<=n} k Lambda(k)
    def psi0(s, y): return s.psi[y] - s.lam[y]/2       # y integer
    def psi1(s, x):                                     # int_0^x psi, x integer
        return x*s.psi[x] - s.nlam[x]
    def A(s, x):  # A(x)=sum_rho x^{rho+1}/(rho(rho+1)) via Ingham Thm 28
        tr = sum(x**(1-2*r)/(2*r*(2*r-1)) for r in range(1,30))
        return x*x/2 - C0*x + Z1 - tr - s.psi1(x)
    def P(s, y):  # P(y)=sum_rho y^rho/rho = y - psi0(y) - c0 - 1/2 log(1-y^-2), vectorised
        y = np.asarray(y)
        return y - (s.psi[y] - s.lam[y]/2) - C0 - 0.5*np.log1p(-1.0/y.astype(float)**2)
    def T_stieltjes(s, N):
        # T = 2 int_N^{2N-2} P(2N-s) dP(s); dP = ds - dpsi - d(1/2 log(1-s^-2))
        # ds term: 2(A(N)-A(2)); dpsi term: -2 sum' Lambda(m) P(2N-m), weight 1/2 at the endpoints
        m = np.arange(N, 2*N-1)
        w = np.ones(len(m)); w[0] = w[-1] = 0.5
        jump = np.sum(w*s.lam[m]*s.P(2*N-m))
        # log term: -int P(2N-s) d log(1-s^-2) ; O(1/N), coarse quadrature
        ss = np.arange(N, 2*N-2) + 0.5
        tl = -np.sum(s.P((2*N-ss).astype(int)) * 2/(ss**3-ss))
        return 2*(s.A(N)-s.A(2)) - 2*jump + tl
    def S_direct(s, N):   # sum_{n<=2N} r_G(n) = sum_{m1} Lambda(m1) psi(2N-m1)  (direct definition)
        m = np.arange(1, 2*N)
        return math.fsum(s.lam[m]*s.psi[2*N-m])
    def rG(s, n):
        m = np.arange(1, n)
        return math.fsum(s.lam[m]*s.lam[n-m])
