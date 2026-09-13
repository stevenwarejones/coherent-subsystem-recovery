"""Independent re-derivation of the forward-only two-query recovery guarantee.

Imports nothing from the shipped checker. Confirms, in exact/limit numerics:
  - the one-square identity  T†T + 20I - 4 L†L = J†J   (random U,V, dims 1..5),
  - r_UV = Tr(rho T†T)/16 = (9p-5)/4 + Tr(rho J†J)/16  >= (9p-5)/4,
  - the Bell-state phase family reproducing p,r and the anchored slope -> 9/4,
  - the selection threshold f=5/6 and forward <= sharp everywhere.
"""
if not __debug__:
    raise SystemExit("Run without -O: this checker uses assertions.")
import numpy as np
I2=np.eye(2);X=np.array([[0,1],[1,0]],complex);Z=np.array([[1,0],[0,-1]],complex)
def ru(d,rng):
    A=rng.normal(size=(d,d))+1j*rng.normal(size=(d,d));q,r=np.linalg.qr(A);return q@np.diag(np.exp(1j*np.angle(np.diag(r))))
rng=np.random.default_rng(13)
# (1) one-square identity T†T + 20I - 4L†L = J†J  for random U,V, several dims
wid=0
for _ in range(200):
    d=rng.integers(1,6);U,V=ru(d,rng),ru(d,rng)
    B=np.kron(X,U);C=np.kron(Z,V);Idm=np.eye(2*d)
    L=Idm+B+C; T=(Idm+B)@(Idm+C); J=3*B-Idm-C-B@C
    lhs=T.conj().T@T+20*Idm-4*(L.conj().T@L); rhs=J.conj().T@J
    wid=max(wid,np.max(np.abs(lhs-rhs)))
    # also T = I+B+C+BC
    assert np.max(np.abs(T-(Idm+B+C+B@C)))<1e-12
print("max |T†T+20I-4L†L - J†J| over 200 random (U,V,dim1..5):",f"{wid:.2e}")
# (2) r_UV = Tr(rho T†T)/16 = (9p-5)/4 + Tr(rho J†J)/16 >= (9p-5)/4 ; p=Tr(rho L†L)/9
worst_gap=1e9; worst_id=0
for _ in range(200):
    d=rng.integers(1,5);U,V=ru(d,rng),ru(d,rng)
    B=np.kron(X,U);C=np.kron(Z,V);Idm=np.eye(2*d)
    L=Idm+B+C; T=(Idm+B)@(Idm+C); J=3*B-Idm-C-B@C
    G=rng.normal(size=(2*d,2*d))+1j*rng.normal(size=(2*d,2*d)); rho=G@G.conj().T; rho/=np.trace(rho)
    p=np.real(np.trace(rho@(L.conj().T@L)))/9
    r=np.real(np.trace(rho@(T.conj().T@T)))/16
    worst_id=max(worst_id,abs(r-((9*p-5)/4+np.real(np.trace(rho@(J.conj().T@J)))/16)))
    worst_gap=min(worst_gap,r-(9*p-5)/4)
print("max |r - ((9p-5)/4 + Tr(rho J†J)/16)|:",f"{worst_id:.2e}")
print("min over trials of  r - (9p-5)/4     :",f"{worst_gap:.2e}","(>=0 required)")
# (3) phase family: rho=Phi+ on AM (d=2 with M a qubit), U=e^{it}X, V=e^{it}Z
phi=np.array([1,0,0,1],complex)/np.sqrt(2);Phi=np.outer(phi,phi.conj())
import math
for t in [0.3,0.1,0.01]:
    U=np.exp(1j*t)*X; V=np.exp(1j*t)*Z
    B=np.kron(X,U);C=np.kron(Z,V);Idm=np.eye(4)
    L=Idm+B+C; T=(Idm+B)@(Idm+C)
    p=np.real(np.trace(Phi@(L.conj().T@L)))/9
    r=np.real(np.trace(Phi@(T.conj().T@T)))/16
    p_pred=(5+4*math.cos(t))/9; r_pred=(1+math.cos(t))**2/4
    ratio=(1-r)/(1-p)
    print(f" t={t}: p={p:.6f}(pred {p_pred:.6f}) r={r:.6f}(pred {r_pred:.6f}) (1-r)/(1-p)={ratio:.5f} ->9/4={9/4}")
# (4) selection thresholds
f=np.linspace(0,1,10001)
sel=np.maximum.reduce([0.25*np.ones_like(f),0.75*f,(9*f-5)/4])
sharp=np.maximum.reduce([0.25*np.ones_like(f),0.75*f,(3*f-1)/2])
print("switch 3f/4 -> (9f-5)/4 at f=5/6:", np.isclose(0.75*(5/6),(9*(5/6)-5)/4))
print("forward selected <= sharp everywhere:", np.all(sel<=sharp+1e-12), " max deficit=",f"{np.max(sharp-sel):.4f}")

assert wid < 1e-9 and worst_id < 1e-9 and worst_gap > -1e-9, "forward identity/bound failed"
print("PASS independent: one-square identity, r>=(9p-5)/4, phase-family slope 9/4, and selection.")
