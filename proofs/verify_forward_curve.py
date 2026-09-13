"""Exact algebra plus numerical physical corroboration. No repository imports.
Run with Python 3, sympy and numpy. This checks arithmetic, not the prose norm proof.
"""
from collections import defaultdict
import sympy as s
import numpy as np

def require(x, message):
    if not x: raise ValueError(message)

def reduce(w):
    out=[]
    for a in w:
        if out and out[-1] == -a: out.pop()
        else: out.append(a)
    return tuple(out)
def add(*ps):
    q=defaultdict(int)
    for p in ps:
        for w,c in p.items(): q[w]+=c
    return {w:c for w,c in q.items() if c}
def scale(p,c): return {w:c*v for w,v in p.items() if c*v}
def mul(p,q):
    ans=defaultdict(int)
    for w,c in p.items():
        for v,d in q.items(): ans[reduce(w+v)]+=c*d
    return {w:c for w,c in ans.items() if c}
def adj(p): return {tuple(-a for a in reversed(w)):c for w,c in p.items()}
I={():1}; B={(1,):1}; C={(2,):1}
L=add(I,B,C); T=mul(add(I,B),add(I,C)); J=add(scale(B,4),scale(T,-1))
require(not add(mul(adj(T),T),scale(I,20),scale(mul(adj(L),L),-4),scale(mul(adj(J),J),-1)), 'FREE_IDENTITY')
require(bool(add(mul(adj(T),T),scale(I,21),scale(mul(adj(L),L),-4),scale(mul(adj(J),J),-1))), 'MUTATION_CONTROL')
require(not add(mul(adj(L),L),scale(I,-1),scale(mul(adj(B),T),-1),scale(mul(adj(T),B),-1)), 'OVERLAP_IDENTITY')
print('PASS exact free-unitary identities; changed constant rejected')
x=s.symbols('x',real=True)
p=(5+4*x)/9; r=(1+x)**2/4; h=lambda f:(9*f-1)**2/64
require(s.expand(r-h(p))==0,'ATTAINMENT')
require(s.expand(h(x)-(9*x-5)/4-81*(1-x)**2/64)==0,'IMPROVEMENT')
c=(11+4*s.sqrt(7))/27
require(s.simplify(h(c)-3*c/4)==0,'SWITCH')
require(s.simplify((3*c-2)/4-(4*s.sqrt(7)-7)/36)==0,'GAP')
print('PASS symbolic phase equality, improvement, switch and maximum gap')
print('switch =',s.N(c,16),'maximum selected-guarantee gap =',s.N((4*s.sqrt(7)-7)/36,16))

# Direct physical circuit: W maps M to output(a), discarded qubit(c), M.
X=np.array([[0,1],[1,0]],complex); Z=np.diag([1,-1]).astype(complex); eye=np.eye(2)
phi=np.array([1,0,0,1],complex)/np.sqrt(2)
Ps=[eye,X,Z,X@Z]
betas=[np.kron(eye,P)@phi for P in Ps]
def effect(U,V):
    d=len(U); Ws=[np.eye(d),U,V,U@V]
    W=sum(np.kron(beta[:,None],w) for beta,w in zip(betas,Ws))/2
    require(np.linalg.norm(W.conj().T@W-np.eye(d))<1e-11,'ISOMETRY')
    # Contract protected A and retained a with Bell bra; environment axes c,M.
    total=np.kron(eye,W).reshape(2,2,2,d,2*d)
    E=np.einsum('aacki->cki',total)/np.sqrt(2)
    E=E.reshape(2*d,2*d)
    return E.conj().T@E
rng=np.random.default_rng(20260913)
def unitary(d): return np.linalg.qr(rng.normal(size=(d,d))+1j*rng.normal(size=(d,d)))[0]
max_error=0.
for d in [1,2,3,4,6]:
    for trial in range(100):
        U,V=unitary(d),unitary(d)
        b=np.kron(X,U); cc=np.kron(Z,V); ident=np.eye(2*d)
        ll=ident+b+cc; tt=(ident+b)@(ident+cc)
        ef=effect(U,V)
        error=np.linalg.norm(ef-tt.conj().T@tt/16); max_error=max(max_error,error)
        require(error<1e-11,'BORN_EFFECT')
        v=rng.normal(size=2*d)+1j*rng.normal(size=2*d); v/=np.linalg.norm(v)
        pp=float(np.linalg.norm(ll@v)**2/9); rr=float(np.vdot(v,ef@v).real)
        require(rr+1e-12>=max(0,9*pp-1)**2/64,'CURVE')
# Physical phase family, including direct original test amplitudes.
for t in np.linspace(0,np.pi,101):
    U=np.exp(1j*t)*X; V=np.exp(1j*t)*Z
    psi=np.concatenate([phi,np.kron(eye,U)@phi,np.kron(eye,V)@phi])/np.sqrt(3)
    target=np.concatenate([np.kron(eye,P)@phi for P in [eye,X,Z]])/np.sqrt(3)
    f=abs(np.vdot(target,psi))**2
    rr=float(np.vdot(phi,effect(U,V)@phi).real)
    require(abs(rr-max(0,9*f-1)**2/64)<1e-12,'PHYSICAL_ATTAINMENT')
    # Reversed query order with paired fourth target ZX has the same value.
    b=np.kron(X,U); cc=np.kron(Z,V); rev=(np.eye(4)+cc)@(np.eye(4)+b)
    require(abs(np.linalg.norm(rev@phi)**2/16-rr)<1e-12,'REVERSE_ORDER')
# At t=pi, fixed forward recovery is zero; a later original-test decoder X makes f=0.
U=-X; V=-Z
psi=np.concatenate([phi,np.kron(eye,U)@phi,np.kron(eye,V)@phi])/np.sqrt(3)
target=np.concatenate([np.kron(eye,P)@phi for P in [eye,X,Z]])/np.sqrt(3)
changed=np.kron(np.eye(3),np.kron(eye,X))@psi
require(abs(np.vdot(target,changed))<1e-12,'LOW_SCORE_ENDPOINT')
print('PASS 500 random matrix/circuit cases; maximum Born-effect discrepancy',max_error)
print('PASS 101 exact-family numerical instances, reversed order, low-score endpoint')
print('Scope: universal proof is in docs/FORWARD_RECOVERY.md; random checks are corroboration only.')
