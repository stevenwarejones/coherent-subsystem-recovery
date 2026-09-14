"""Exact Laurent-polynomial check for the restricted Bell/phase family."""
if not __debug__: raise SystemExit("Run without -O")
import sympy as s
z,w=s.symbols('z w',nonzero=True)
X=s.Matrix([[0,1],[1,0]]);Z=s.diag(1,-1);I=s.eye(2);Had=(X+Z)/s.sqrt(2)
U=z*X;V=w*Z
R0=(I+V).col_join(U*(I-V))/2
Rv=(V*(I+V)).col_join(U*(I-V))/2
Ru=s.kronecker_product(Had,I)*(U*(I+U)).col_join(V*(I-U))/2
# Conjugate unit-circle variables by their inverses.
def adj(M):return M.conjugate().T.subs({s.conjugate(z):1/z,s.conjugate(w):1/w})
phi=s.Matrix([1,0,0,1])/s.sqrt(2)
a=(z+1/z)/2;b=(w+1/w)/2;c=(z/w+w/z)/2
expected=[(1+a)*(1+b)/4,(1+b)*(1+c)/4,(1+a)*(1+c)/4]
for R,pred in zip((R0,Rv,Ru),expected):
 assert s.simplify(adj(R)*R-s.eye(2))==s.zeros(2)
 A,B=R[:2,:],R[2:,:]
 E=(adj(A)*A).row_join(adj(A)*B).col_join((adj(B)*A).row_join(adj(B)*B))/2
 score=(adj(phi)*E*phi)[0]
 assert s.simplify(score-pred)==0
margin=((1-a)*(1-b)+(1-a)*(1-c)+(1-b)*(1-c))/12
assert s.simplify(sum(expected)/3-(a+b+c)/3-margin)==0
f=(3+2*(a+b+c))/9
assert s.simplify((3*f-1)/2-(a+b+c)/3)==0
print('PASS exact phase-family scores, isometries, and nonnegative-product margin identity')
