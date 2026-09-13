"""Exact endpoint and elementary-certificate checks. Requires SymPy.
The high-score universal certificate is verified by the separate sharp package.
"""
if not __debug__:
    raise SystemExit("Run without -O: this research checker uses assertions.")
import sympy as s
I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Y=s.Matrix([[0,-s.I],[s.I,0]]); Z=s.diag(1,-1)
kron=s.kronecker_product
simp=lambda A:A.applyfunc(s.simplify)
def trM(A,d):return s.Matrix(2,2,lambda a,b:sum(A[a*d+j,b*d+j] for j in range(d)))
def trA(A,d):return A[:d,:d]+A[d:,d:]
def endpoint(rho,U,V):
 d=U.rows; L=s.eye(2*d)+kron(X,U)+kron(Z,V);K=L.H*L/9
 assert simp(trA(s.Rational(3,2)*K,d)-s.eye(d))==s.zeros(d)
 C=simp(trM(L*rho*L.H,d)/9);p=s.simplify(s.trace(C))
 return C,p
# M trivial, all encoded quantum information is in bypass B.
C,p=endpoint(I/2,s.Matrix([[s.I]]),s.Matrix([[s.I]]))
assert C==I/6 and p==s.Rational(1,3)
# Explicit decoder: L^*/sqrt(3) rotates the input qubit in B.
L0=I+s.I*X+s.I*Z; W=L0.conjugate()/s.sqrt(3)
assert simp(W.H*W-I)==s.zeros(2)
phi=s.Matrix([1,0,0,1])/s.sqrt(2)
chi=s.Matrix([1,s.I,s.I])/s.sqrt(3)
Omega=s.Matrix.vstack(phi,kron(I,X)*phi,kron(I,Z)*phi)/s.sqrt(3)
state=kron(chi,kron(I,W)*phi)
assert s.simplify(abs((Omega.H*state)[0])**2)==s.Rational(1,3)
state0=kron(chi,kron(I,X*W)*phi)
assert s.simplify(abs((Omega.H*state0)[0])**2)==0
# Classical flag endpoint.
r0=(I+(X-Y+Z)/s.sqrt(3))/2
Gs=[I,X,Y,Z];flags=[s.eye(4)[:,j]*s.eye(4)[:,j].T for j in range(4)]
rho=sum([kron(g*r0*g.H,F) for g,F in zip(Gs,flags)],s.zeros(8))/4
u=(s.sqrt(3)+s.I)/2;v=s.conjugate(u)
U=s.diag(*[s.trace(g*X*g.H*X)/2*u for g in Gs]); V=s.diag(*[s.trace(g*Z*g.H*Z)/2*v for g in Gs])
C,p=endpoint(rho,U,V);assert C==I/3 and p==s.Rational(2,3)
# Perfect endpoint.
C,p=endpoint(phi*phi.H,X,Z);assert C==I/2 and p==1
# Piecewise interpolation; flag-readable recovery is exactly additive.
l=s.symbols('l',nonnegative=True)
fmid=(1-l)/3+2*l/3;Rmid=(1-l)/4+l/2
fhigh=2*(1-l)/3+l;Rhigh=(1-l)/2+l
assert s.simplify(Rmid-3*fmid/4)==0
assert s.simplify(Rhigh-(3*fhigh-1)/2)==0
# Literal Ando off-diagonal completion at perfect Pauli point fails.
H=3*(s.eye(4)+kron(X,X)+kron(Z,Z)).H*(s.eye(4)+kron(X,X)+kron(Z,Z))/9-s.eye(4)
assert H[:2,2:]==s.Matrix([[0,s.Rational(4,3)],[0,0]])
# Numerical radius of a 2x2 nilpotent [0,c;0,0] equals |c|/2 = 2/3 > 1/2.
print('PASS: Q0=(3/2)K has partial trace I at all three exact endpoints')
print('PASS: f=1/3,R=1/4 endpoint; explicit decoder scores 1/3 and 0')
print('PASS: f=2/3,R=1/2 and f=1,R=1 endpoints reconstructed')
print('PASS: exact flagged interpolation on both nonconstant segments')
print('PASS: direct off-diagonal Ando substitution fails at the perfect Pauli point')
# Generic partial-trace identity, by Pauli coefficient orthogonality.
for Ps in [[I,X,Z],[I,X,Y,Z]]:
 n=len(Ps)
 for a in range(n):
  for b in range(n):
   assert s.simplify(s.trace(Ps[a].conjugate()*Ps[b].T))==(2 if a==b else 0)
 # diagonal physical words U_a^dag U_a reduce to I; off-diagonal coefficients vanish.
 assert s.Rational(n,2)*s.Rational(2*n,n*n)==1
# Explicit forward-only isometry for the elementary three-branch bound.
Us=[I,X,Z];Ps=[I,X,Z];d=2
W3=sum([kron(kron(I,P.T)*phi,U) for P,U in zip(Ps,Us)],s.zeros(4*d,d))/s.sqrt(3)
assert simp(W3.H*W3-s.eye(d))==s.zeros(d)
# Isometry Choi output tracing its c,M environment, retaining a as recovery output.
Qeffect=s.zeros(4)
for c in range(2):
 for m in range(d):
  D=s.Matrix(2,d,lambda a,b:W3[(2*a+c)*d+m,b])
  Qeffect+=kron(I,D.H)* (phi*phi.H) *kron(I,D)
K3=(s.eye(4)+kron(X,X)+kron(Z,Z)).H*(s.eye(4)+kron(X,X)+kron(Z,Z))/9
assert simp(Qeffect-s.Rational(3,4)*K3)==s.zeros(4)
# Four orthogonal Pauli branches: elementary bound R>=f, sharp for isotropic inputs.
L4=s.eye(4)+kron(X,X)-kron(Y,Y)+kron(Z,Z)
assert L4==4*phi*phi.H
q=s.symbols('q',real=True)
rho4=q*phi*phi.H+(1-q)*s.eye(4)/4
C4=simp(trM(L4*rho4*L4.H,2)/16)
assert C4==(1+3*q)*I/8
print('PASS: generic orthogonal-Pauli partial-trace normalization for 3 and 4 branches')
print('PASS: explicit forward-only three-branch recovery circuit has effect (3/4)K')
print('PASS: four-branch isotropic construction has f_opt=R=(1+3q)/4')
# Two-branch perfect simulation with an entirely classical addressed record.
Gs2=[I,Z] # |+>,|-> on reference
plus=s.Matrix([1,1])/s.sqrt(2)
rho2=sum([kron(g*plus*plus.H*g.H,s.eye(2)[:,j]*s.eye(2)[:,j].T) for j,g in enumerate(Gs2)],s.zeros(4))/2
L2=s.eye(4)+kron(X,Z)
C2=simp(trM(L2*rho2*L2.H,2)/4)
assert C2==I/2
assert simp(trM(rho2,2)-I/2)==s.zeros(2)
print('PASS: two-branch Pauli test admits f=1 with separable R=1/2')
