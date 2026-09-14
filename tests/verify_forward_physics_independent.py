"""Independent gate-level and score-decoder reconstruction; numerical corroboration.
No imports from the supplied verifier/checker modules. Uses an eigensquare-root
and QR unitary extension rather than their Cholesky/SVD construction.
"""
from pathlib import Path
from itertools import product
from fractions import Fraction
import json
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
if not __debug__: raise SystemExit('Run without -O')
rng=np.random.default_rng(13092026)
I=np.eye(2,dtype=complex);X=np.array([[0,1],[1,0]],complex);Z=np.diag([1,-1]).astype(complex);Y=1j*X@Z;H=(X+Z)/np.sqrt(2)
flags=[(I,1,1),(X,1,-1),(Z,-1,1),(Y,-1,-1)]
worst=0.
def check(a,b,label,tol=2e-8):
    global worst
    err=float(np.max(abs(np.asarray(a)-np.asarray(b))))
    if not np.isfinite(err) or err>tol: raise RuntimeError(f'{label}: {err}')
    worst=max(worst,err)
def load(path): return json.loads((ROOT/path).read_text())
def numeric(rows): return np.array([[float(Fraction(x)) for x in row] for row in rows])
def histories(order):
    return [tuple(g for g,bit in reversed(list(zip(order,bits))) if bit) for bits in product((0,1),repeat=len(order))]
def evaluate(w,U,V):
    out=np.eye(len(U),dtype=complex)
    for letter in w: out=out@{'U':U,'V':V}[letter]
    return out
def realize(M,words):
    eigen,basis=np.linalg.eigh(M)
    if min(eigen)<=0: raise RuntimeError('positive moment required')
    vectors=np.sqrt(eigen)[:,None]*basis.T
    ix={w:i for i,w in enumerate(words)};m=len(words);dim=len(M);oracles=[]
    for letter in 'UV':
        domain=[w for w in words if (letter,)+w in ix]
        col=[a*m+ix[w] for a in range(2) for w in domain]
        image=[a*m+ix[(letter,)+w] for a in range(2) for w in domain]
        source=vectors[:,col];target=vectors[:,image]
        Q,R=np.linalg.qr(source,mode='complete');rank=len(col)
        image_basis=np.linalg.solve(R[:rank,:].T,target.T).T
        extended,_=np.linalg.qr(image_basis,mode='complete')
        oracle=image_basis@Q[:,:rank].T+extended[:,rank:]@Q[:,rank:].T
        check(oracle.T@oracle,np.eye(dim),'QR unitary')
        check(oracle@source,target,'partial isometry extension')
        oracles.append(oracle.astype(complex))
    psi=np.stack((vectors[:,ix[()]],vectors[:,m+ix[()]]))
    columns=np.column_stack([evaluate(w,*oracles)@psi[a] for a in range(2) for w in words])
    check(columns.conj().T@columns,M,'all realized words')
    return *oracles,psi

def instrument(psi,V):
    # Actual prepare |+>, controlled V, H; reference/controller/target axes.
    state=np.repeat(psi[:,None,:],2,axis=1)/np.sqrt(2)
    state[:,1,:]=state[:,1,:]@V.T
    return np.einsum('bc,acm->abm',H,state)
def balanced(psi,U,V):
    successes=[]
    for bit in (0,1):
        state=instrument(psi,V)
        state[:,bit,:]=state[:,bit,:]@U.T
        amplitude=(state[0,0,:]+state[1,1,:])/np.sqrt(2)
        successes.append(np.vdot(amplitude,amplitude).real)
    return sum(successes)/2

def calibration(psi,U,V):
    first=instrument(psi,V);prob=np.zeros((2,2))
    for initial in range(2):
        second=instrument(first[:,initial,:]@U.T,V)
        for final in range(2): prob[initial,final]=np.sum(abs(second[:,final,:])**2)
    check(prob.sum(),np.sum(abs(psi)**2),'four outcome normalization')
    return prob[0,0]+prob[1,1]

def score_decoder(components):
    # Orthogonal purifying flags are compressed into disjoint column blocks.
    # These columns belong to M B, so the constructed common decoder acts
    # only there. J is the ideal CA support isometry (six rows, two columns).
    J=np.vstack((I,X.T,Z.T))/np.sqrt(3)
    projected=[];reference=np.zeros((2,2),complex)
    for weight,psi,U,V in components:
        actual=np.vstack((psi,psi@U.T,psi@V.T))/np.sqrt(3)
        projected.append(np.sqrt(weight)*J.conj().T@actual)
        reference+=weight*psi@psi.conj().T
    check(reference,I/2,'Bell-input reference marginal')
    A=np.hstack(projected);C=A@A.conj().T;p=np.trace(C).real
    check(C,p*I/2,'scalar projected marginal')
    left,singular,right=np.linalg.svd(A,full_matrices=False)
    D0=(left@right).conj()
    check(D0@D0.conj().T,I,'decoder isometry rows')
    # Complete D0 to a unitary with arbitrary orthogonal rows. A annihilates
    # their transpose, so the full unconditional score is this contribution.
    check(A-(A@D0.T)@D0.conj(),np.zeros_like(A),'decoder complement contributes zero')
    f=abs(np.trace(A@D0.T))**2/2
    check(f,p,'actual common-decoder success')
    return f

def random_unitary(n):
    Q,R=np.linalg.qr(rng.normal(size=(n,n))+1j*rng.normal(size=(n,n)))
    return Q@np.diag(np.diag(R)/abs(np.diag(R)))

def controller(order,size):
    initial=random_unitary(size)[:,0];gates=[random_unitary(size) for _ in order]
    coeff=initial[:,None]
    mask=np.r_[np.zeros(size//2),np.ones(size//2)]
    for gate in gates:
        coeff=np.stack((gate@((1-mask)[:,None]*coeff),gate@(mask[:,None]*coeff)),axis=-1).reshape(size,-1)
    return initial,gates,coeff.reshape(2,size//2,-1)

def run_controller(psi,U,V,order,initial,gates):
    state=psi[:,None,:]*initial[None,:,None]
    for letter,gate in zip(order,gates):
        state[:,len(initial)//2:,:]=state[:,len(initial)//2:,:]@({'U':U,'V':V}[letter]).T
        state=np.einsum('bc,acm->abm',gate,state)
    state=state.reshape(2,2,len(initial)//2,len(U))
    accepted=(state[0,0]+state[1,1])/np.sqrt(2)
    return np.sum(abs(accepted)**2)

def causal(T):
    n=len(T)//2;G=T[:n,:n]+T[n:,n:]
    while len(G)>1:
        check(G[::2,1::2],np.zeros_like(G[::2,1::2]),'causal cross-bit block')
        G=G[::2,::2]+G[1::2,1::2]
    check(G,[[1]],'causal normalization')

joint=load('certificates/joint_adversary.json')
words=[tuple({1:'U',2:'V'}[x] for x in w) for w in joint['words']]
U,V,psi=realize(numeric(joint['moment']),words)
components=[(.25,g@psi,sx*U,sz*V) for g,sx,sz in flags]
f=score_decoder(components)
r=sum(w*balanced(state,A,B) for w,state,A,B in components)
q=sum(w*calibration(state,A,B) for w,state,A,B in components)
check([f,r,q],[float(Fraction(joint[k])) for k in ('p','r','q')],'joint direct gates')
if not r<(3*f-1)/2 or not q<(1-f)/8: raise RuntimeError('joint physical failure')

rob=load('certificates/forward_queries/robust_VUVUV.json');order='VUVUV';hs=histories(order);words=list(dict.fromkeys(hs))
M=numeric(rob['moment']);U,V,psi=realize(M,words);theta=float(Fraction(rob['theta']))
ideal=I/np.sqrt(2)
components=[(theta/4,g@psi,sx*U,sz*V) for g,sx,sz in flags]+[((1-theta)/4,g@ideal,sx*X,sz*Z) for g,sx,sz in flags]
check(score_decoder(components),.9,'robust measured score')
ix={w:i for i,w in enumerate(words)};indices=[a*len(words)+ix[w] for a in range(2) for w in hs]
Ca=M[np.ix_(indices,indices)]/2
A=np.column_stack([evaluate(w,X,Z)[:,a] for a in range(2) for w in hs]);C0=A.conj().T@A/4
Cm=theta*Ca+(1-theta)*C0
for size in (2,4,8):
    for repeat in range(6):
        initial,gates,coeff=controller(order,size);Ts=[]
        for g,sx,sz in flags:
            sign=np.array([sx**w.count('U')*sz**w.count('V') for w in hs])
            rotated=np.einsum('ab,beh,h->aeh',g.T,coeff,sign).reshape(2*(size//2),32)
            rotated=rotated.reshape(2,size//2,32).transpose(0,2,1).reshape(64,size//2)
            T=rotated@rotated.conj().T;causal(T);Ts.append(T)
        actual=sum(weight*run_controller(state,A,B,order,initial,gates) for weight,state,A,B in components)
        predicted=np.trace(Cm@sum(Ts)/4).real
        check(actual,predicted,'independent actual flagged controller versus objective')
        if actual>float(Fraction(rob['alpha']))+1e-8: raise RuntimeError('robust upper bound')
# Rotated, complex reflections, rather than only diagonal V examples.
for size in (2,3,5):
    for repeat in range(8):
        change=random_unitary(size);V=change@np.diag([1]+[-1]*(size-1))@change.conj().T;U=random_unitary(size)
        I_m=np.eye(size);P=(I_m+V)/2;N=(I_m-V)/2
        effects=[]
        for blocks in ((P,U@N),(U@P,N)):
            effects.append(np.block([[blocks[a].conj().T@blocks[b]/2 for b in range(2)] for a in range(2)]))
        E=sum(effects)/2;L=np.eye(2*size)+np.kron(X,U)+np.kron(Z,V)
        diff=E-(L.conj().T@L-np.eye(2*size))/8
        check(np.linalg.norm(diff,2),np.linalg.norm(U@V+V@U,2)/8,'rotated reflection norm')
print(json.dumps({'independent_gate_and_decoder_checks':'PASS','robust_complex_controllers':18,'rotated_reflection_instances':24,'max_residual':worst},indent=2))
