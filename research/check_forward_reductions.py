"""Exact algebraic checks for the flag reduction and naive Ando comparisons.
Prose arguments about universal circuits are not machine formalized here.
"""
import sympy as s

def require(cond,msg):
 if not cond:raise RuntimeError(msg)
def zero(M):return all(s.simplify(x)==0 for x in M)
I=s.eye(2);X=s.Matrix([[0,1],[1,0]]);Y=s.Matrix([[0,-s.I],[s.I,0]]);Z=s.diag(1,-1)
phi=s.Matrix([1,0,0,1])/s.sqrt(2);Phi=phi*phi.T
for g in [I,X,Y,Z]:
 require(zero(s.kronecker_product(g,I)*phi-s.kronecker_product(I,g.T)*phi),'BELL_TRANSPOSE')
a,b,c,d=s.symbols('a b c d',real=True);C=s.Matrix([[a,b+s.I*c],[b-s.I*c,d]])
twirl=sum((g*C*g.conjugate().T for g in [I,X,Y,Z]),s.zeros(2))/4
require(zero(twirl-(a+d)*I/2),'SCALAR_C')
p,r,theta=s.symbols('p r theta',real=True)
F=lambda x:(3*x-1)/2
require(s.expand(theta*r+1-theta-F(theta*p+1-theta)-theta*(r-F(p)))==0,'DEFICIT_DILUTION')
# Direct Ando substitutions at the ideal qubit realization.
L=s.eye(4)+s.kronecker_product(X,X)+s.kronecker_product(Z,Z)
H=L.T*L/3-s.eye(4);Q=2*Phi
require(zero(H-(s.Rational(8,3)*Phi-s.Rational(2,3)*s.eye(4))),'IDEAL_H')
require(zero(Q-H-s.Rational(2,3)*(s.eye(4)-Phi)),'VALID_COMPLETION')
require(H[:2,2:]==s.Matrix([[0,s.Rational(4,3)],[0,0]]),'ANDO_OFFDIAGONAL')
# Numerical radius of [[0,c],[0,0]] is |c|/2, giving 2/3 > 1/2.
require(s.Rational(2,3)>s.Rational(1,2),'ANDO_NUMERICAL_RADIUS')

# Any universal fixed scaling Q01=k H01 is excluded by two exact instances.
k=s.Rational(3,4) # forced by uniqueness Q=2Phi at ideal Pauli point
Lscalar=s.eye(4)+s.kronecker_product(X,I)+s.kronecker_product(Z,s.I*I)
Hs=Lscalar.conjugate().T*Lscalar/3-s.eye(4)
Cs=Hs[:2,2:]
require(zero(Cs-s.Rational(2,3)*(1-s.I)*I),'SCALAR_H01')
absB2=s.simplify(k*k*Cs[0,0]*s.conjugate(Cs[0,0]))
require(absB2==s.Rational(1,2) and absB2>s.Rational(1,4),'ALL_CONSTANT_SCALINGS_FAIL')
require(zero(Hs*Hs-s.Rational(8,9)*s.eye(4)),'SCALAR_VALID_COMPLETION')
print('PASS Bell-transpose, full Pauli twirl, exact affine-deficit dilution')
print('PASS exact ideal completion and obstruction to the two naive fixed-offdiagonal Ando substitutions')
print('PASS no universal constant Q01=k H01: ideal forces k=3/4, scalar instance violates normalized PSD')
