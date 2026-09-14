"""Numerical circuit and moment-realization corroboration; not the exact proof."""
import json,sys,numpy as np
from pathlib import Path
from fractions import Fraction
from scipy.linalg import null_space
from scipy.stats import unitary_group
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'proofs'))
from verify_joint_recovery import WORDS
root=Path(__file__).resolve().parents[1] / 'certificates'
if not __debug__:raise SystemExit('Run without -O')
X=np.array([[0,1],[1,0]],complex);Z=np.diag([1,-1]).astype(complex);Y=1j*X@Z

def objects(U,V):
 n=len(U);I=np.eye(n);P=(I+V)/2;N=(I-V)/2
 Rs=[np.vstack((P,U@N)),np.vstack((U@P,N))]
 E=sum(np.block([[R[a*n:(a+1)*n].conj().T@R[b*n:(b+1)*n]/2 for b in range(2)]for a in range(2)])for R in Rs)/2
 L=np.eye(2*n)+np.kron(X,U)+np.kron(Z,V);K=L.conj().T@L/9
 Q=sum((A@U@A).conj().T@(A@U@A)for A in [P,N]);return E,K,Q,L

def born(U,V,rho):
 n=len(U);E,K,Q,L=objects(U,V);Cfull=L@rho@L.conj().T/9
 C=np.array([[np.trace(Cfull[a*n:(a+1)*n,b*n:(b+1)*n])for b in range(2)]for a in range(2)])
 return np.trace(K@rho).real,np.trace(E@rho).real,np.trace(np.kron(np.eye(2),Q)@rho).real,C

cert=json.loads((root/'joint_certificate.json').read_text());G=np.array([[float(Fraction(x))for x in row]for row in cert['gram']]);rng=np.random.default_rng(9038)
maxid=0.;min_eval=1.;count=0
for n in [2,3,4,7]:
 for j in range(30):
  U=unitary_group.rvs(n,random_state=rng);V=unitary_group.rvs(n,random_state=rng);E,K,Q,L=objects(U,V)
  operators=[]
  for a in range(2):
   for w in WORDS:
    O=np.eye(n,dtype=complex)
    for x in w:O=O@{1:U,2:V}[x]
    J=np.zeros((n,2*n),complex);J[:,a*n:(a+1)*n]=O;operators.append(J)
  S=sum(G[i,j]*(operators[i].conj().T@operators[j])for i in range(14)for j in range(14))
  T=E-1.25*K+3*np.kron(np.eye(2),Q)+.25*np.eye(2*n)
  maxid=max(maxid,np.max(abs(S-T)));min_eval=min(min_eval,np.linalg.eigvalsh(T)[0]);count+=1
assert maxid<1e-12 and min_eval>-1e-12

ad=json.loads((root/'joint_adversary.json').read_text());H=np.array([[float(Fraction(x))for x in row]for row in ad['moment']]);m=7;n=14
vectors=np.linalg.cholesky(H).T;ix={w:i for i,w in enumerate(WORDS)}
def extend(x):
 dom=[w for w in WORDS if (x,)+w in ix]
 D=np.column_stack([vectors[:,a*m+ix[w]]for a in range(2)for w in dom]);B=np.column_stack([vectors[:,a*m+ix[(x,)+w]]for a in range(2)for w in dom])
 Qd,sv,R=np.linalg.svd(D,full_matrices=False)
 Qi=B@R.conj().T@np.diag(1/sv)
 return Qi@Qd.conj().T+null_space(Qi.conj().T)@null_space(Qd.conj().T).conj().T
U,V=extend(1),extend(2);unitarity=max(np.max(abs(O.conj().T@O-np.eye(n)))for O in [U,V]);assert unitarity<1e-7
psi=np.concatenate([vectors[:,0],vectors[:,m]]);rho=np.outer(psi,psi.conj())
values=[]
for g,sx,sz in [(np.eye(2),1,1),(X,1,-1),(Z,-1,1),(Y,-1,-1)]:
 gg=np.kron(g,np.eye(n));values.append(born(sx*U,sz*V,gg@rho@gg.conj().T))
p,r,q=[sum(v[k]for v in values)/4 for k in range(3)];C=sum(v[3]for v in values)/4
errors=[abs(p-float(Fraction(ad['p']))),abs(r-float(Fraction(ad['r']))),abs(q-float(Fraction(ad['q'])))]
assert max(errors)<1e-7 and np.max(abs(C-p*np.eye(2)/2))<1e-7
f=(np.sqrt(np.maximum(0,np.linalg.eigvalsh(C))).sum())**2/2
assert abs(f-p)<1e-7 and r<(3*f-1)/2 and q<(1-f)/8
result=dict(random_operator_instantiations=count,SOS_identity_residual=maxid,min_random_operator_eigenvalue=min_eval,
 adversary_unitarity_residual=unitarity,adversary_Born_errors=errors,adversary_f_opt=f,adversary_r=r,adversary_q=q,
 flag_C_scalar_residual=float(np.max(abs(C-p*np.eye(2)/2))))
print(json.dumps(result,indent=2))
