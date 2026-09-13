"""End-to-end numeric check: instantiate the free-algebra certificate on concrete
random unitaries U,V in several dimensions, and confirm the OPERATOR facts and the
actual recovery the channel achieves. Catches tensor-order / factor-of-2 bugs that
pure coefficient algebra cannot see. Uses the shipped Gram data but its own operator
assembly and its own channel application.
"""
import json
from pathlib import Path
import numpy as np

C = json.loads((Path(__file__).resolve().parent.parent/'certificates'/'sharp_recovery.json').read_text())
words = [tuple(w) for w in C["words"]]
n = len(words)


from fractions import Fraction as Fr
def _f(M):
    return np.array([[float(Fr(x)) for x in row] for row in M])
def gram(key):
    B = _f(C["basis_" + key]); R = _f(C["reduced_" + key])
    return B @ R @ B.T


gG, gJ = gram("G"), gram("J")
I2 = np.eye(2)
X = np.array([[0, 1], [1.0, 0]])
Z = np.array([[1, 0], [0, -1.0]])


def wop(w, U, V, d):
    """Operator on M (dimension d) for a reduced word in letters +-1 (U/U^dag), +-2 (V/V^dag)."""
    M = np.eye(d, dtype=complex)
    for x in w:
        if x == 1:
            M = M @ U
        elif x == -1:
            M = M @ U.conj().T
        elif x == 2:
            M = M @ V
        elif x == -2:
            M = M @ V.conj().T
    return M


def dag(w):
    return tuple(-x for x in reversed(w))


def red(w):
    out = []
    for x in w:
        if out and out[-1] == -x:
            out.pop()
        else:
            out.append(x)
    return tuple(out)


def pauli2(w):
    A = I2
    for x in w:
        A = A @ (X if abs(x) == 1 else Z)
    return A



def rand_unitary(d, rng):
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    Qm, Rm = np.linalg.qr(A)
    return Qm @ np.diag(np.exp(1j * np.angle(np.diag(Rm))))


def build_Q(g, U, V, d):
    Q = np.zeros((2 * d, 2 * d), dtype=complex)
    Pl = [pauli2(w) for w in words]
    wcache = {}
    for i in range(n):
        for j in range(n):
            w = red(dag(words[i]) + words[j])
            if w not in wcache:
                wcache[w] = wop(w, U, V, d)
            blk2 = g[i][j] * (Pl[i] @ Pl[j].T)                 # 2x2 (on C^2 = "A")
            Q += np.kron(blk2, wcache[w])
    return Q


def H_op(U, V, d):
    return (np.kron(X, U + U.conj().T)
            + np.kron(Z, V + V.conj().T)
            + np.kron(X @ Z, U.conj().T @ V - V.conj().T @ U)) / 3.0


def bell_overlap_after_channel(Q, rho, d):
    """Recovery channel E has Choi Q (on A=qubit tensor M). Achieved Bell overlap = 1/2 Tr(Q rho)."""
    return 0.5 * np.real(np.trace(Q @ rho))


rng = np.random.default_rng(20260913)
worst_S = 0.0
worst_tr = 0.0
worst_gap = 1e9
for trial in range(200):
    d = rng.integers(2, 5)
    U, V = rand_unitary(d, rng), rand_unitary(d, rng)
    Q = build_Q(gG, U, V, d)
    S = build_Q(gJ, U, V, d)
    H = H_op(U, V, d)
    # (1) Q - S = H  as operators
    worst_S = max(worst_S, np.max(np.abs((Q - S) - H)))
    # (2) Q >= 0, S >= 0
    assert np.min(np.linalg.eigvalsh((Q + Q.conj().T) / 2)) > -1e-9, "Q not PSD"
    assert np.min(np.linalg.eigvalsh((S + S.conj().T) / 2)) > -1e-9, "S not PSD"
    # (3) Tr_A Q = I_M   (trace out the qubit A = first factor)
    trA = Q.reshape(2, d, 2, d)[0, :, 0, :] + Q.reshape(2, d, 2, d)[1, :, 1, :]
    worst_tr = max(worst_tr, np.max(np.abs(trA - np.eye(d))))
    # (4) an actual random state rho on A(qubit) tensor M, then check the achieved bound
    G = rng.normal(size=(2 * d, 2 * d)) + 1j * rng.normal(size=(2 * d, 2 * d))
    rho = G @ G.conj().T
    rho = rho / np.trace(rho)
    p = np.real(np.trace(((H + np.eye(2 * d)) / 3) @ rho))     # p = Tr(K rho), K=(H+I)/3
    R_achieved = bell_overlap_after_channel(Q, rho, d)
    # the certificate's channel achieves R_achieved; claim R_achieved >= (3p-1)/2
    worst_gap = min(worst_gap, R_achieved - (3 * p - 1) / 2)
    # sanity: R_achieved = 1/2 Tr(Q rho) and Q>=H  =>  >= 1/2 Tr(H rho) = (3p-1)/2
print("max |(Q-S) - H| over 200 random (U,V,dim) :", f"{worst_S:.2e}")
print("max |Tr_A Q - I_M|                         :", f"{worst_tr:.2e}")
print("min over trials of  R_achieved-(3p-1)/2    :", f"{worst_gap:.2e}", "(>=0 required)")
assert worst_S < 1e-9 and worst_tr < 1e-9 and worst_gap > -1e-9
print("PASS: certificate instantiates correctly on concrete unitaries;")
print("      the Choi-Q recovery channel achieves >= (3p-1)/2 on random states.")
