"""Numerical check of both flagged blocks and the Bell-effect indexing.
Random controllers are generated as actual circuits, not arbitrary PSD matrices.
"""
from pathlib import Path
import itertools,json,sys
import numpy as np
from scipy.linalg import null_space
from scipy.stats import unitary_group
if not __debug__: raise SystemExit('Run without -O')
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'proofs'))
from verify_robust_forward_obstruction import F,words,build
P=Path(__file__).resolve().parents[1]/'certificates'/'forward_queries'
X=np.array([[0,1],[1,0]],complex);Z=np.diag([1,-1]);Y=1j*X@Z;paulis=[np.eye(2),X,Y,Z]

def evaluate(w,U,V):
 a=np.eye(len(U),dtype=complex)
 for x in w:a=a@{1:U,2:V}[x]
 return a

def realize(G,W):
 Fm=np.linalg.cholesky(G).conj().T;m=len(W);ix={w:i for i,w in enumerate(W)};out=[]
 for x in (1,2):
  dom=[w for w in W if (x,)+w in ix]
  A=Fm[:,[a*m+ix[w]for a in range(2)for w in dom]]
  B=Fm[:,[a*m+ix[(x,)+w]for a in range(2)for w in dom]]
  q,s,v=np.linalg.svd(A,full_matrices=False);r=np.count_nonzero(s>1e-9);q=q[:,:r];s=s[:r];v=v[:r]
  t=(B@v.conj().T)/s
  out.append(t@q.conj().T+null_space(t.conj().T)@null_space(q.conj().T).conj().T)
 U,V=out;psi=np.stack((Fm[:,0],Fm[:,m]));reb=np.column_stack([evaluate(w,U,V)@psi[a]for a in range(2)for w in W])
 assert np.max(np.abs(reb.conj().T@reb-G))<1e-10
 assert max(np.linalg.norm(u.conj().T@u-np.eye(len(U)),2)for u in (U,V))<1e-10
 return U,V,psi

def controller(n,k,rng):
 v=unitary_group.rvs(k,random_state=rng)[:,0];hist=[v];proj=[np.r_[np.ones(k//2),np.zeros(k//2)],np.r_[np.zeros(k//2),np.ones(k//2)]]
 for _ in range(n):
  q=unitary_group.rvs(k,random_state=rng);hist=[q@(p*v)for v in hist for p in proj]
 return np.stack(hist,axis=-1).reshape(2,k//2,2**n)

def gram(c):
 A=c.transpose(0,2,1).reshape(2*c.shape[-1],c.shape[1]);return A@A.conj().T

def score(c,ws,psi):
 amplitudes=np.einsum('aeh,hmj,aj->em',c,ws,psi)
 return np.vdot(amplitudes,amplitudes).real/2

def main():
 d=json.loads((P/'robust_VUVUV.json').read_text());hs=words(d['order']);W=list(dict.fromkeys(hs));m=len(W)
 G,C,p,theta=build(d);G=np.array(G,float);C=np.array(C,float);theta=float(theta);U,V,psi=realize(G,W);dim=len(U)
 ws=np.array([evaluate(w,U,V)for w in hs]);ideal=np.array([evaluate(w,X,Z)for w in hs]);phi=np.eye(2)/np.sqrt(2)
 rng=np.random.default_rng(9814);maxerr=0.;maxCerr=0.;idealcollapse=0.;wrongreference=0.
 for k in (2,4,8):
  for repeat in range(10):
   coeff=controller(len(d['order']),k,rng);Tbar=np.zeros((2*len(hs),2*len(hs)),complex);physical=0.;collapsed=0.;wrong=0.;Ctest=np.zeros((2,2),complex)
   for g in paulis:
    su=int(round(np.trace(g@X@g.conj().T@X).real/2));sv=int(round(np.trace(g@Z@g.conj().T@Z).real/2))
    signs=np.array([np.prod([su if x==1 else sv for x in w])for w in hs])
    signed=coeff*signs
    transformed=np.einsum('ab,beh->aeh',g.T,signed);Tbar+=gram(transformed)/4
    radv=score(signed,ws,g@psi);rid=score(signed,ideal,g@phi)
    physical+=(theta*radv+(1-theta)*rid)/4
    collapsed+=(theta*radv+(1-theta)*score(coeff,ideal,phi))/4
    wrong+=(theta*score(signed,ws,psi)+(1-theta)*rid)/4
    L=np.eye(2*dim)+np.kron(X,su*U)+np.kron(Z,sv*V)
    a=(L@(g@psi).reshape(-1)).reshape(2,dim)
    Li=np.eye(4)+np.kron(X,su*X)+np.kron(Z,sv*Z)
    b=(Li@(g@phi).reshape(-1)).reshape(2,2)
    Ctest+=(theta*a@a.conj().T+(1-theta)*b@b.conj().T)/36
   poly=np.trace(C@Tbar).real;maxerr=max(maxerr,abs(physical-poly));idealcollapse=max(idealcollapse,abs(collapsed-poly));wrongreference=max(wrongreference,abs(wrong-poly))
   maxCerr=max(maxCerr,np.max(np.abs(Ctest-.45*np.eye(2))))
   assert physical<=float(F(d['alpha']))+1e-9
 assert maxerr<1e-10 and maxCerr<1e-10
 assert idealcollapse<1e-10,'Ideal Pauli-frame covariance failed'
 assert wrongreference>1e-4,'Missing reference rotation was not distinguished'
 result=dict(actual_random_circuits=30,controller_dimensions=[2,4,8],addressed_dimension=4*dim+2,max_Born_vs_Gram_error=maxerr,max_scalar_C_error=maxCerr,ideal_flag_collapse_error=idealcollapse,missing_reference_rotation_negative_control_difference=wrongreference)
 print(json.dumps(result,indent=2))
if __name__=='__main__':main()
