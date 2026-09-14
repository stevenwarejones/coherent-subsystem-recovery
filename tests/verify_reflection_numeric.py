"""Numerical reflection identity and norm corroboration, not a universal proof."""
import numpy as np
from scipy.stats import unitary_group
import json
if not __debug__: raise SystemExit("Run without -O")
X=np.array([[0,1],[1,0]],float);Z=np.diag([1,-1])
rng=np.random.default_rng(4871);worst=0.;worst_norm=0.
for n in [2,3,4,6]:
 for j in range(30):
  U=unitary_group.rvs(n,random_state=rng);k=1+j%(n-1);V=np.diag([1]*k+[-1]*(n-k));Pp=(np.eye(n)+V)/2;Pm=(np.eye(n)-V)/2
  R0=np.vstack((Pp,U@Pm));R1=np.vstack((U@Pp,Pm))
  def effect(R):
   B=[R[:n],R[n:]]
   return np.block([[B[a].conj().T@B[b]/2 for b in range(2)]for a in range(2)])
  E=(effect(R0)+effect(R1))/2;L=np.eye(2*n)+np.kron(X,U)+np.kron(Z,V)
  diff=E-(L.conj().T@L-np.eye(2*n))/8;D=Pp@U@Pp+Pm@U.conj().T@Pm
  right=-np.block([[np.zeros((n,n)),D],[D.conj().T,np.zeros((n,n))]])/4
  worst=max(worst,np.max(np.abs(diff-right)))
  worst_norm=max(worst_norm,abs(np.linalg.norm(diff,2)-np.linalg.norm(U@V+V@U,2)/8))
assert worst<1e-12 and worst_norm<1e-12
result=dict(random_tests=120,dimensions=[2,3,4,6],identity_residual=worst,norm_identity_residual=worst_norm)
print(json.dumps(result,indent=2))
