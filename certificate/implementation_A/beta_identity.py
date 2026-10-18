# Check: Gamma(r1)*(2N)^r2*Gamma(r2)/Gamma(r1+r2+1) - (2N)^r2/r2*(B_{1/N}(r2+1,r1)+B_{1/2}(r1,r2+1))
#      = (2N)^r2/r2 * int_{1/N}^{1/2} u^r2 (1-u)^(r1-1) du      (bracket of T(N), per pair)
import mpmath as mp
mp.mp.dps = 30
Z = [mp.mpf(s) for s in "14.134725142 21.022039639 25.010857580 30.424876126".split()]
rho = lambda g: mp.mpc(0.5, g)
for N in (10, 50):
    for r1 in (rho(Z[0]), rho(-Z[1]), rho(Z[3])):
        for r2 in (rho(Z[0]), rho(Z[2]), rho(-Z[3])):
            X = 2*N
            lhs = mp.gamma(r1)*X**r2*mp.gamma(r2)/mp.gamma(r1+r2+1) \
                - X**r2/r2*(mp.betainc(r2+1, r1, 0, mp.mpf(1)/N) + mp.betainc(r1, r2+1, 0, mp.mpf(1)/2))
            rhs = X**r2/r2*mp.quad(lambda u: u**r2*(1-u)**(r1-1), [mp.mpf(1)/N, 0.1, 0.2, 0.3, 0.4, 0.5])
            print(N, mp.nstr(r1.imag, 6), mp.nstr(r2.imag, 6), mp.nstr(abs(lhs-rhs)/abs(rhs), 3))
