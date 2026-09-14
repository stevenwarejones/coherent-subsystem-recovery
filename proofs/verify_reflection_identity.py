"""Exact reflection block identity; arbitrary blocks enter only linearly."""
import sympy as s
if not __debug__: raise SystemExit("Run without -O")
X=s.Matrix([[0,1],[1,0]]);Z=s.diag(1,-1);I=s.eye(2)
a,b,c,d,aa,bb,cc,dd=s.symbols('a b c d aa bb cc dd')
U=s.Matrix([[a,b],[c,d]]);Ud=s.Matrix([[aa,cc],[bb,dd]])
V=Z;Pp=(I+V)/2;Pm=(I-V)/2
E=(Pp/2).row_join(Pp*(U+Ud)*Pm/4).col_join((Pm*(U+Ud)*Pp/4).row_join(Pm/2))
# Expanded using U^dag U=V^2=I; remaining expression is linear in U,U^dag.
K9=3*s.eye(4)+s.kronecker_product(X,U+Ud)+2*s.kronecker_product(Z,V)+s.kronecker_product(X*Z,Ud*V-V*U)
D=Pp*U*Pp+Pm*Ud*Pm;Dd=Pp*Ud*Pp+Pm*U*Pm
rhs=-s.zeros(2).row_join(D).col_join(Dd.row_join(s.zeros(2)))/4
assert s.simplify(E-(K9-s.eye(4))/8-rhs)==s.zeros(4)
print('PASS exact reflection block identity, with independent formal adjoint entries')
