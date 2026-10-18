# F(N) = S'(N) - 2N^2 + 2A(2N-2) - T(N), S' counted DIRECTLY from the definition (sum_{m1} Lambda(m1) psi(2N-m1)),
# A by Ingham Thm 28, T in Stieltjes form (validated in zeros_check.py). Output: N, F/N, P^2/N, T/N, rG(2N)/2N
import numpy as np, math, sys
from core import Ctx, C0
lo, hi, n, seed = int(float(sys.argv[1])), int(float(sys.argv[2])), int(sys.argv[3]), int(sys.argv[4])
rng = np.random.default_rng(seed)
Ns = sorted(set(rng.integers(lo, hi, n).tolist()))
# N with 2N highly composite (large r_G(2N)) and neighbours
for base in (30030*19, 510510, 9699690//2):
    k = 1
    while base*k//2 < hi:
        if base*k % 2 == 0 and base*k//2 > lo: Ns += [base*k//2 - 1, base*k//2, base*k//2 + 1]
        k += 1
Ns = sorted(set(Ns))
c = Ctx(2*hi + 4)
def _bf(x):  # psi by brute force (trial division)
    t=0.0
    for k in range(2,x+1):
        m=k; p=2
        while p*p<=m and m%p: p+=1
        p = p if p*p<=m else m
        while m%p==0: m//=p
        if m==1: t+=math.log(p)
    return t
assert abs(c.psi[3000]-_bf(3000))<1e-9
print("N\tF_over_N\tP2_over_N\tT_over_N\trG2N_over_2N\tP")
for N in Ns:
    S = c.S_direct(N); r = c.rG(2*N)
    Sp = S - r/2
    T = c.T_stieltjes(N)
    F = Sp - 2*N*N + 2*c.A(2*N-2) - T
    P = float(c.P(np.array([N]))[0])
    print(f"{N}\t{F/N:.7f}\t{P*P/N:.7f}\t{T/N:.5f}\t{r/(2*N):.4f}\t{P:.3f}", flush=True)
