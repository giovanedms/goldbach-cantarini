# Check: Gamma(r1)Gamma(r2)/Gamma(r1+r2+1)*(2N)^r2 - (2N)^r2/r2*(B_{1/N}(r2+1,r1)+B_{1/2}(r1,r2+1))
#      = (2N)^r2/r2 * int_{1/N}^{1/2} v^r2 (1-v)^(r1-1) dv
import mpmath as mp
mp.mp.dps=30
g=[mp.mpf('14.134725141734693790'),mp.mpf('21.022039638771554993'),mp.mpf('25.010857580145688763')]
for N in (10,37):
  for a in g:
    for b in g:
      for s1 in (1,-1):
        r1=mp.mpc(.5,s1*a); r2=mp.mpc(.5,b)
        lhs=mp.gamma(r1)*mp.gamma(r2)/mp.gamma(r1+r2+1)*(2*N)**r2 - (2*N)**r2/r2*(mp.betainc(r2+1,r1,0,mp.mpf(1)/N)+mp.betainc(r1,r2+1,0,mp.mpf(1)/2))
        rhs=(2*N)**r2/r2*mp.quad(lambda v: v**r2*(1-v)**(r1-1),[mp.mpf(1)/N,mp.mpf(1)/2])
        print(N, mp.nstr(abs(lhs-rhs)/abs(rhs),3))
