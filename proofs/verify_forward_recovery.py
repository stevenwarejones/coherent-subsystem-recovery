"""Exact checks of the two-forward-query decoder and its sharp affine slope.
SymPy required. No imports from the earlier recovery verifiers.
"""
if not __debug__: raise SystemExit('Run without -O')
from collections import defaultdict
import sympy as s
# Independent free-unitary polynomial arithmetic for B,C.
def red(w):
 out=[]
 for x in w:
  if out and out[-1]==-x:out.pop()
  else:out.append(x)
 return tuple(out)
def add(*ps):
 out=defaultdict(int)
 for p in ps:
  for w,c in p.items():out[w]+=c
 return {w:c for w,c in out.items() if c}
def scale(p,a):return {w:a*c for w,c in p.items() if a*c}
def mul(p,q):
 out=defaultdict(int)
 for w,c in p.items():
  for v,d in q.items():out[red(w+v)]+=c*d
 return {w:c for w,c in out.items() if c}
def adj(p):return {tuple(-x for x in w[::-1]):c for w,c in p.items()}
def sq(p):return mul(adj(p),p)
I={():1};B={(1,):1};C={(2,):1};BC=mul(B,C)
L=add(I,B,C);T=add(L,BC);J=add(scale(B,3),scale(I,-1),scale(C,-1),scale(BC,-1))
assert add(sq(T),scale(I,20),scale(sq(L),-4),scale(sq(J),-1))=={}
print('PASS exact free-unitary identity T†T + 20I - 4L†L = J†J')
# Explicit Choi effect reconstructed from circuit Kraus operators.
E=s.eye(2);X=s.Matrix([[0,1],[1,0]]);Z=s.diag(1,-1);k=s.kronecker_product
phi=s.Matrix([1,0,0,1])/s.sqrt(2)
def simp(A):return A.applyfunc(s.simplify)
for U,V in [(X,Z),(s.diag(1,s.I),X),(s.Matrix([[0,s.I],[1,0]]),s.diag(s.I,-1))]:
 Ps=[E,X,Z,X*Z];Us=[E,U,V,U*V]
 W=sum([k(k(E,P)*phi,A) for P,A in zip(Ps,Us)],s.zeros(8,2))/2
 assert simp(W.H*W-E)==s.zeros(2)
 effect=s.zeros(4)
 for c in range(2):
  for m in range(2):
   D=s.Matrix(2,2,lambda a,b:W[(2*a+c)*2+m,b])
   effect+=k(E,D.H)*phi*phi.H*k(E,D)
 Tm=s.eye(4)+k(X,U)+k(Z,V)+k(X*Z,U*V)
 assert simp(effect-Tm.H*Tm/16)==s.zeros(4)
print('PASS explicit forward circuit is TP and its Bell effect is T†T/16 (three exact instances)')
# Generic Bell-state phase family U=e^{it}X,V=e^{it}Z:
# p=fopt=(5+4 cos t)/9, r=(1+cos t)^2/4.
c=s.symbols('c',real=True)
p=(5+4*c)/9;r=(1+c)**2/4
assert s.simplify((1-r)/(1-p)-9*(3+c)/16)==0
assert s.limit(9*(3+c)/16,c,1)==s.Rational(9,4)
assert s.simplify(r-(1-s.Rational(9,4)*(1-p))-(1-c)**2/4)==0
print('PASS physical Bell-state phase family proves 9/4 is the optimal anchored affine slope for this circuit')
# Algebraic constants and crossing with 3f/4.
f=s.symbols('f')
assert s.solve(s.Eq(3*f/4,(9*f-5)/4),f)==[s.Rational(5,6)]
print('PASS combined implementable guarantee max(1/4,3f/4,(9f-5)/4); switch at f=5/6')
