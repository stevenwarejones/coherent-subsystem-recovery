import numpy as np
import cvxpy as cp
I=np.eye(2); X=np.array([[0,1],[1,0]]); Z=np.diag([1,-1]); phi=np.array([1,0,0,1])/np.sqrt(2)
F=np.column_stack([np.kron(I,P)@phi for P in [I,X,Z,X@Z]])
def state(t):
 U=np.exp(1j*t)*X; V=np.exp(1j*t)*Z
 W=sum(np.kron(F[:,j,None],w) for j,w in enumerate([I,U,V,U@V]))/2
 v=(np.kron(I,W)@phi).reshape(8,2)
 return v@v.conj().T
# Choi ordering output o, input controller i, J[o,i;q,j].
def objective_matrix(rho):
 # Bell expectation = 1/2 sum rho[o,i;q,j] J[o,i;q,j].
 return rho.T/2
J=cp.Variable((8,8),hermitian=True)
constraints=[J >> 0,J[:4,:4]+J[4:,4:]==np.eye(4)]
r0=state(0)
for t in [.1,.4,.8,1.2]:
 rp,rm=state(t),state(-t); A=objective_matrix((rp+rm)/2)
 for pinned in [False,True]:
  c=constraints+([cp.real(cp.trace(objective_matrix(r0)@J))==1] if pinned else [])
  prob=cp.Problem(cp.Maximize(cp.real(cp.trace(A@J))),c)
  val=prob.solve(solver='CLARABEL',tol_gap_abs=1e-9,tol_feas=1e-9,tol_gap_rel=1e-9)
  print(t,pinned,prob.status,val,'h=',(1+np.cos(t))**2/4,flush=True)
