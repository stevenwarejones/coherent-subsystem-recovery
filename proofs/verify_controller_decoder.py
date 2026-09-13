"""Exact controller-only impossibility certificate, with no repository imports.
Dependency: sympy. Prose inequalities and quantifiers are explained in RESEARCH.md.
"""
import sympy as s

def require(p,label):
    if not p: raise ValueError(label)
def zero(M): return all(s.simplify(v)==0 for v in M)
I=s.eye(2); X=s.Matrix([[0,1],[1,0]]); Z=s.diag(1,-1); Y=s.Matrix([[0,-s.I],[s.I,0]])
phi=s.Matrix([1,0,0,1])/s.sqrt(2); Phi=phi*phi.H
F=s.Matrix.hstack(*[s.kronecker_product(I,P)*phi for P in [I,X,Z,X*Z]])
require(zero(F.H*F-s.eye(4)),'BELL_BASIS')
def controller_state(z):
    U=z*X; V=z*Z
    W=sum((s.kronecker_product(F[:,j],w) for j,w in enumerate([I,U,V,U*V])),s.zeros(8,2))/2
    require(zero(W.H*W-I),'QUERY_ISOMETRY')
    v=s.kronecker_product(I,W)*phi
    # axes A,a,c,M. Trace final M.
    return s.Matrix(8,8,lambda j,k:sum(v[2*j+m]*s.conjugate(v[2*k+m]) for m in range(2)))
rho0=s.kronecker_product(Phi,I/2)
require(zero(controller_state(s.Integer(1))-rho0),'IDEAL_STATE')
J0=s.kronecker_product(2*Phi,I) # Choi ordering output, input a, input c
require(zero(J0[:4,:4]+J0[4:,4:]-s.eye(4)),'TP')
require(s.trace(rho0.T*J0)/2==1,'PERFECT_RECOVERY')
D0=s.eye(8)/4-rho0.T/2
require(D0.rank()==6 and zero(D0*J0),'ENDPOINT_KERNEL')
# Parameterize every supported J by Phi_(output,a) tensor S_c.
a,b,d,e=s.symbols('a b d e')
S=s.Matrix([[a,b],[d,e]]); Js=s.kronecker_product(Phi,S)
TP=Js[:4,:4]+Js[4:,4:]-s.eye(4)
require(s.solve(list(TP),[a,b,d,e])=={a:2,b:0,d:0,e:2},'UNIQUE_CHANNEL')
print('PASS exact ideal controller state, endpoint support and unique TP Choi')

c=s.symbols('c',real=True)
YY=s.kronecker_product(Y,Y); O=s.kronecker_product(X,X)+s.kronecker_product(Z,Z)
E=(1+c)/2*s.eye(4)+(1-c)/2*YY
EE=s.kronecker_product(I,E); OO=s.kronecker_product(I,O)
rhobar=s.simplify(EE*rho0*EE+(1-c*c)/4*OO*rho0*OO)
Q=(1+c)**2/16*s.eye(4)+(1-c*c)/16*YY
D=s.simplify(s.kronecker_product(I,Q)-rhobar.T/2)
require(zero(D-D.H),'HERMITIAN_SLACK')
l=s.symbols('lambda')
expected=l**2*(l-(1+c)/8)**2*(l-c/4)**2*(l-(2*c-1)*(c+1)/8)**2
require(s.expand(D.charpoly(l).as_expr()-expected)==0,'EXACT_SPECTRUM')
require(s.simplify(s.trace(Q)-(1+c)**2/4)==0,'DUAL_OBJECTIVE')
require(s.simplify(s.trace(rhobar.T*J0)/2-(1+c)**2/4)==0,'MATCHING_PRIMAL')
# Arithmetic sign proof after c=(1+u)/2, 0<=u<=1.
u=s.symbols('u',nonnegative=True)
for ev in [(1+c)/8,c/4,(2*c-1)*(c+1)/8]:
    require(all(v>=0 for v in s.Poly(s.expand(ev.subs(c,(1+u)/2)),u).all_coeffs()),'PSD_INTERVAL')
# Wrong sign/normalization controls, with exact witnesses of failure.
require((D.subs(c,1)-s.eye(8)/100).det()!=0,'TIGHTENING_CONTROL_SETUP')
require((D.subs(c,1)-s.eye(8)/100).eigenvals().get(-s.Rational(1,100),0)>0,'TIGHTENING_REJECTS')
require(s.factor(((2*c-1)*(c+1)/8).subs(c,s.Rational(1,4)))<0,'OUTSIDE_RANGE_REJECTS')
print('PASS exact dual spectrum, PSD interval, trace, primal match and two rejection controls')
# Independent direct isometry reconstructions at rational phases; no floating arithmetic.
for q in [s.Integer(0),s.Rational(1,4),s.Rational(1,3),s.Rational(1,2),s.Integer(1)]:
    z=(1+s.I*q)**2/(1+q*q); cv=(1-q*q)/(1+q*q)
    rp,rm=controller_state(z),controller_state(s.conjugate(z))
    require(zero((rp+rm)/2-rhobar.subs(c,cv)),'DIRECT_PHASE_AVERAGE')
    for rr in [rp,rm]:
        require(s.simplify(s.trace(rr.T*J0)/2-(1+cv)**2/4)==0,'EACH_SIGN_ATTAINS')
print('PASS five exact rational-phase circuit reconstructions, both signs')
f=s.symbols('f',real=True)
require(s.simplify(((1+c)**2/4).subs(c,(9*f-5)/4)-(9*f-1)**2/64)==0,'SCORE_TRANSLATION')
require(s.simplify((3*f-1)/2-(9*f-1)**2/64+3*(f-1)*(27*f-11)/64)==0,'SHARP_GAP')
print('CERTIFIED ARITHMETIC: controller-only minimax h(f) for 7/9 <= f <= 1')
print('No claim about all two-query circuits, phase-informed maps, or maps retaining M.')
