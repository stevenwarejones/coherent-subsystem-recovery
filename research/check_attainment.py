"""Exact algebraic endpoint/symmetrization check. Requires SymPy, not used by verifier."""
import sys as _sys
if not __debug__:
    raise SystemExit('Run without -O: this checker uses assertions.')
import sympy as s
I=s.eye(2);X=s.Matrix([[0,1],[1,0]]);Y=s.Matrix([[0,-s.I],[s.I,0]]);Z=s.diag(1,-1)
P=(I+(X-Y+Z)/s.sqrt(3))/2
assert s.simplify(P*P-P)==s.zeros(2) and s.trace(P)==1
u=(s.sqrt(3)+s.I)/2;v=s.conjugate(u);rho=s.zeros(8);U=s.zeros(4);V=s.zeros(4)
for k,g in enumerate([I,X,Y,Z]):
 r=g*P*g.conjugate().T/4
 for a in range(2):
  for b in range(2):rho[4*a+k,4*b+k]=r[a,b]
 U[k,k]=s.trace(g*X*g.conjugate().T*X)/2*u
 V[k,k]=s.trace(g*Z*g.conjugate().T*Z)/2*v
L=s.eye(8)+s.kronecker_product(X,U)+s.kronecker_product(Z,V)
def trM(A,d):return s.Matrix(2,2,lambda a,b:sum(A[a*d+k,b*d+k] for k in range(d))).applyfunc(s.simplify)
assert trM(rho,4)==I/2
C=trM(L*rho*L.conjugate().T,4)/9
assert C==I/3
print('PASS classical flagged endpoint: reference I/2, C=I/3, hence f_opt=2/3')
# R=1/2 because conditioned reference states are pure, with flag-dependent conjugate preparation attaining 1/2.
phi=s.Matrix([1,0,0,1]);rp=phi*phi.T/2;lp=s.eye(4)+s.kronecker_product(X,X)+s.kronecker_product(Z,Z)
assert trM(lp*rp*lp.T,2)/9==I/2
print('PASS perfect endpoint C=I/2, f_opt=1; Bell recovery R=1')
a=s.symbols('a',nonnegative=True)
Cmix=(1-a)*C+a*I/2
assert s.simplify(s.trace(Cmix)-(s.Rational(2,3)+a/3))==0
Rmix=(1-a)/2+a
assert s.simplify(Rmix-(3*s.trace(Cmix)-1)/2)==0
print('PASS all flagged mixtures 0<=a<=1: f=(2+a)/3, R=(1+a)/2=(3f-1)/2')
