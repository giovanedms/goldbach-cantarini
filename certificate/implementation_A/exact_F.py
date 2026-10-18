# Exact evaluation of Cantarini's F(N) = S'(N) - 2N^2 + 2A(2N-2) - T(N)   (Thm 1, arXiv 1801.08475v1)
#   S'(N) = sum_{n<=2N} r_G(n) - r_G(2N)/2,  via sum = 2 sum_{a<=N} L(a) psi(2N-a) - psi(N)^2   (paper eq. 18)
#   A(x)  = sum_rho x^(rho+1)/(rho(rho+1)) = x^2/2 - psi_1(x) - c0 x + O(1)          (Ingham Thm 28; O(1) dropped)
#   T(N)  = 2 sum_rho1 int_2^N P(t)(2N-t)^(rho1-1) dt   (Beta identity: beta_identity.py; inner sum over rho2 by Lemma 5)
#         = 2 int_2^N P(t) dt - 2 sum'_{N<=m<=2N-2} L(m) P(2N-m) + O(1/N)   (outer sum limit: checked in outer_sum.py)
#   P(x)  = sum_rho x^rho/rho = x - psi0(x) - c0 - (1/2)log(1-x^-2)
# Prediction (see the paper): F(N) = P(N)^2 + O(N), O(N) part smooth; i.e. F is NOT O(N) (Littlewood).
# usage: python3 -I exact_F.py NMAX NSAMPLES
import sys, numpy as np
NMAX = int(sys.argv[1]); NS = int(sys.argv[2])
M = 2*NMAX + 2
c0 = np.log(2*np.pi)

# --- sieve Lambda up to M
isp = np.ones(M+1, bool); isp[:2] = False
for p in range(2, int(M**0.5)+1):
    if isp[p]: isp[p*p::p] = False
L = np.zeros(M+1)
L[isp] = np.log(np.nonzero(isp)[0])
for p in np.nonzero(isp[:int(M**0.5)+1])[0]:
    q = int(p)*int(p)
    while q <= M: L[q] = np.log(p); q *= int(p)
psi = np.cumsum(L)                                   # psi[k] = psi(k)
psi0 = psi - L/2
psi1 = np.concatenate(([0.0], np.cumsum(psi.astype(np.longdouble))[:-1]))  # psi1[x] = int_0^x psi = sum_{k<x} psi(k)

def P(x):    # x integer >= 2, midpoint convention
    return x - psi0[x] - c0 - 0.5*np.log1p(-1.0/x**2)

def int_logterm(a, b):   # int_a^b log(1-1/t^2) dt
    G = lambda t: t*np.log1p(-1/t**2) + np.log((t+1)/(t-1))
    return G(b) - G(a)

def F_of(N):
    a = np.arange(1, N+1)
    S = 2*np.dot(L[1:N+1], psi[2*N-a]) - psi[N]**2
    rG2N = np.dot(L[1:2*N], L[2*N-1:0:-1])
    Sp = S - rG2N/2
    x = 2*N-2
    A = x*x/2.0 - float(psi1[x]) - c0*x
    # int_2^N P(t) dt, P(t) = t - psi(k) - c0 - .5 log(1-t^-2) on (k,k+1)
    k = np.arange(2, N)
    IP = np.sum(k + 0.5 - psi[k] - c0) - 0.5*int_logterm(2.0, float(N))
    # sum' over m in [N, 2N-2] of L(m) P(2N-m); endpoints half weight, one-sided values
    m = np.arange(N+1, 2*N-2)
    t = 2*N - m
    Sm = np.dot(L[m], P(t))
    Sm += 0.5*L[N]*(N - psi[N-1] - c0 - 0.5*np.log1p(-1.0/N**2))
    Sm += 0.5*L[2*N-2]*(2 - psi[2] - c0 - 0.5*np.log(0.75))
    T = 2*IP - 2*Sm
    F = Sp - 2.0*N*N + 2*A - T
    return F, P(N), T, rG2N

rng = np.random.default_rng(1)
Ns = sorted(set(rng.integers(NMAX//5, NMAX, NS).tolist()))
# N with 2N = multiples of primorials (many small factors), and neighbours
special = []
for pr in (30030, 510510, 9699690):
    for c in range(1, 2*NMAX//pr + 1):
        if pr*c//2 <= NMAX and pr*c >= 2*(NMAX//5): special += [pr*c//2 - 1, pr*c//2, pr*c//2 + 1]
special = sorted(set(n for n in special if NMAX//5 <= n <= NMAX))[:60]
print("# N  F/N  P(N)^2/N  (F-P^2)/N  T/N  rG(2N)/(2N)  special")
for N in sorted(set(Ns) | set(special)):
    F, PN, T, r = F_of(N)
    print(N, f"{F/N:.6f} {PN*PN/N:.6f} {(F-PN*PN)/N:.6f} {T/N:.6f} {r/(2*N):.6f}", int(N in special), flush=True)
