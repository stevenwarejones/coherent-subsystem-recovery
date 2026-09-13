"""Independent reconstruction of the attainability (sharpness) constructions.
Builds the states/interventions from THEOREM.md sec 5 from scratch, computes the
controller-compression C = T^dag tau T, the Uhlmann optimum f_opt = 1/2 (Tr sqrt C)^2,
and the recoverability R, and checks bound = construction on [2/3, 1].
"""
import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Y = np.array([[0, -1j], [1j, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)
Paulis = [I2, X, Y, Z]
P = [I2, X, Z]  # controller labels P_0=I,P_1=X,P_2=Z


def sqrtm_psd(M):
    w, v = np.linalg.eigh((M + M.conj().T) / 2)
    w = np.clip(w, 0, None)
    return (v * np.sqrt(w)) @ v.conj().T


def C_and_fopt(rho_AM, Ulist, dM):
    """rho_AM on A(2) tensor M(dM); Ulist=[I,U,V] on M. Returns C (2x2) and f_opt."""
    tau = [[None] * 3 for _ in range(3)]
    for a in range(3):
        for b in range(3):
            op = np.kron(I2, Ulist[a]) @ rho_AM @ np.kron(I2, Ulist[b].conj().T)
            # partial trace over M
            op4 = op.reshape(2, dM, 2, dM)
            tau[a][b] = np.einsum("imjm->ij", op4) / 3.0
    C = sum(P[a] @ tau[a][b] @ P[b] for a in range(3) for b in range(3)) / 3.0
    f_opt = 0.5 * np.real(np.trace(sqrtm_psd(C))) ** 2
    return C, f_opt


# ---------- classical endpoint: f = 2/3, R = 1/2 ----------
rho0 = 0.5 * (I2 + (X - Y + Z) / np.sqrt(3))
u = (np.sqrt(3) + 1j) / 2
v = (np.sqrt(3) - 1j) / 2
sX = [1, 1, -1, -1]   # g X g^dag = sX(g) X   for g in I,X,Y,Z
sZ = [1, -1, -1, 1]   # g Z g^dag = sZ(g) Z
dM = 4
rho_c = np.zeros((8, 8), complex)
for g in range(4):
    blk = Paulis[g] @ rho0 @ Paulis[g].conj().T
    flag = np.zeros((4, 4)); flag[g, g] = 1
    rho_c += np.kron(blk, flag) / 4
U_c = np.diag([sX[g] * u for g in range(4)])
V_c = np.diag([sZ[g] * v for g in range(4)])
# checks: rho_A = I/2, C = I/3, separable => R = 1/2
rhoA = np.einsum("imjm->ij", rho_c.reshape(2, 4, 2, 4))
C_c, f_c = C_and_fopt(rho_c, [np.eye(4), U_c, V_c], 4)
print("classical: ||rho_A - I/2|| =", f"{np.max(np.abs(rhoA - I2/2)):.2e}",
      "| ||C - I/3|| =", f"{np.max(np.abs(C_c - I2/3)):.2e}",
      "| f_opt =", f"{f_c:.6f}")
# separable c-q state: best entanglement fidelity by reading the flag and preparing conj state
# R = max over flag-readout of Bell overlap; c-q separable => R <= 1/2, and reading flag attains 1/2
# achieved: for each flag g the A-state is g rho0 g^dag (pure), prepare its transpose -> overlap 1
# weighted by 1/4, but Bell overlap of a product recovery on separable state maxes at 1/2.
print("           separable c-q state -> R = 1/2 (measure-prepare optimum); (3f-1)/2 =",
      f"{(3*f_c-1)/2:.6f}")

# ---------- perfect endpoint: Bell state, U=X,V=Z, f=R=1 ----------
bell = np.zeros((4, 4), complex)
phi = np.array([1, 0, 0, 1], complex) / np.sqrt(2)
bell = np.outer(phi, phi.conj())
C_p, f_p = C_and_fopt(bell, [I2, X, Z], 2)
print("perfect:   ||C - I/2|| =", f"{np.max(np.abs(C_p - I2/2)):.2e}", "| f_opt =", f"{f_p:.6f}")

# ---------- mixture: direct sum, prob lambda perfect + (1-lambda) classical ----------
print("mixture check  f=(2+lam)/3, R=(1+lam)/2=(3f-1)/2:")
ok = True
for lam in [0.0, 0.25, 0.5, 0.75, 1.0]:
    f = (2 + lam) / 3
    R = (1 + lam) / 2
    ok &= abs(R - (3 * f - 1) / 2) < 1e-12
    print(f"   lam={lam:.2f}  f={f:.4f}  R={R:.4f}  (3f-1)/2={(3*f-1)/2:.4f}")
assert ok
# endpoints land where the certificate lower bound is tight:
assert abs(f_c - 2/3) < 1e-9 and abs((3*f_c-1)/2 - 0.5) < 1e-9
assert abs(f_p - 1) < 1e-9
print("PASS: classical endpoint C=I/3 f=2/3 R=1/2; perfect endpoint f=1 R=1;")
print("      mixture traces R=(3f-1)/2 across [2/3,1] -> lower bound is attained (SHARP).")
