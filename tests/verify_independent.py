"""Fully independent re-verification. Imports NOTHING from the supplied scripts.

Re-establishes, in exact rational arithmetic:
  (a) reduced_G, reduced_J are positive definite  (my own LDL, not their elimination)
  (b) the assembled 106x106 block Gram matrices G_full, J_full are PSD directly
      (so the operator-SOS conclusion Q,S >= 0 does not rest on trusting the
       congruence prose -- I check the Gram witness itself)
  (c) g = B_G reduced_G B_G^T  and  j = B_J reduced_J B_J^T
  (d) Q - S = H as free-unitary polynomials, with H built from the operator
      formula from scratch, words reduced only by adjacent inverse cancellation
  (e) Tr_A Q = I_M  and  Tr_A S = 0 on every nonempty word
"""
import json
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import sys as _sys
if not __debug__:
    raise SystemExit('Run without -O: this checker uses assertions.')
_CERT = Path(__file__).resolve().parent.parent / 'certificates' / 'sharp_recovery.json'

C = json.loads(_CERT.read_text())


def matmul(A, B):
    Bt = list(zip(*B))
    return [[sum((a * b for a, b in zip(r, c)), F(0)) for c in Bt] for r in A]


def T(A):
    return [list(r) for r in zip(*A)]


def frac(M):
    return [[F(x) for x in row] for row in M]


# ---- my own LDL^T; returns pivots, requires symmetric; PD iff all pivots > 0 ----
def ldl_pivots(A):
    n = len(A)
    assert all(len(r) == n for r in A) and A == T(A), "not symmetric square"
    L = [[F(0)] * n for _ in range(n)]
    D = [F(0)] * n
    for j in range(n):
        d = A[j][j] - sum(L[j][k] ** 2 * D[k] for k in range(j))
        D[j] = d
        L[j][j] = F(1)
        for i in range(j + 1, n):
            if d == 0:
                # off-diagonal must also vanish for PSD with a zero pivot
                num = A[i][j] - sum(L[i][k] * L[j][k] * D[k] for k in range(j))
                assert num == 0, "zero pivot with nonzero off-diagonal: not PSD"
                L[i][j] = F(0)
            else:
                L[i][j] = (A[i][j] - sum(L[i][k] * L[j][k] * D[k] for k in range(j))) / d
    # exact reconstruction check
    R = matmul(matmul(L, [[D[i] if i == k else F(0) for k in range(n)] for i in range(n)]), T(L))
    assert R == A, "LDL reconstruction mismatch"
    return D


def red(w):
    out = []
    for x in w:
        if out and out[-1] == -x:
            out.pop()
        else:
            out.append(x)
    return tuple(out)


def dag(w):
    return tuple(-x for x in reversed(w))


# Pauli image: letter +-1 -> X, +-2 -> Z (sign irrelevant, both Hermitian), in order.
XI, XX, ZZ = [[1, 0], [0, 1]], [[0, 1], [1, 0]], [[1, 0], [0, -1]]


def pauli(w):
    A = XI
    for x in w:
        A = matmul(A, XX if abs(x) == 1 else ZZ)
    return A


words = [tuple(w) for w in C["words"]]
n = len(words)
assert n == 53 and len(set(words)) == 53, "word basis"
assert all(len(w) <= 3 and red(w) == w and all(x in (-2, -1, 1, 2) for x in w) for w in words), "words"
assert F(C["t"]) == 1, "t must be 1"
P = [pauli(w) for w in words]

poly = {}
for key in ("G", "J"):
    B = frac(C["basis_" + key])
    Rm = frac(C["reduced_" + key])
    m = len(Rm)
    assert len(B) == n and all(len(r) == m for r in B), "basis shape"
    # (a) reduced matrix positive definite
    piv = ldl_pivots(Rm)
    assert all(p > 0 for p in piv), (key, "reduced Gram not PD")
    # (c) g = B Rm B^T
    g = matmul(matmul(B, Rm), T(B))
    # (b) assemble the full 2n x 2n block Gram: block (i,j) = g_ij * P_i P_j^T, and check PSD
    full = [[F(0)] * (2 * n) for _ in range(2 * n)]
    for i in range(n):
        for j in range(n):
            blk = [[g[i][j] * v for v in row] for row in matmul(P[i], T(P[j]))]
            for r in range(2):
                for c in range(2):
                    full[2 * i + r][2 * j + c] = blk[r][c]
    fpiv = ldl_pivots(full)  # will assert PSD (zero pivots allowed, off-diagonals must vanish)
    assert all(p >= 0 for p in fpiv), (key, "assembled block Gram not PSD")
    rank = sum(1 for p in fpiv if p != 0)
    print(f"PASS {key}: reduced Gram ({m}x{m}) PD; assembled block Gram (106x106) PSD, rank {rank}")
    # (d/e) build the word polynomial for Q (or S)
    pd_ = defaultdict(lambda: [[F(0), F(0)], [F(0), F(0)]])
    for i in range(n):
        for j in range(n):
            a = matmul(P[i], T(P[j]))
            w = red(dag(words[i]) + words[j])
            for r in range(2):
                for c in range(2):
                    pd_[w][r][c] += g[i][j] * a[r][c]
    poly[key] = pd_

# H from scratch: (1/3)[ X⊗(U+U†) + Z⊗(V+V†) + XZ⊗(U†V − V†U) ]
th = F(1, 3)
Hb = {}
Hb[(1,)] = Hb[(-1,)] = [[F(0), th], [th, F(0)]]        # (1/3) X   for U and U†
Hb[(2,)] = Hb[(-2,)] = [[th, F(0)], [F(0), -th]]        # (1/3) Z   for V and V†
XZmat = matmul(XX, ZZ)                                   # = [[0,-1],[1,0]]
Hb[(-1, 2)] = [[th * v for v in row] for row in XZmat]   # (1/3) XZ  for U†V
Hb[(-2, 1)] = [[-th * v for v in row] for row in XZmat]  # -(1/3)XZ  for V†U
zero = [[F(0), F(0)], [F(0), F(0)]]

allw = set(poly["G"]) | set(poly["J"]) | set(Hb)
bad = 0
for w in allw:
    diff = [[poly["G"][w][r][c] - poly["J"][w][r][c] for c in range(2)] for r in range(2)]
    if diff != Hb.get(w, zero):
        bad += 1
    tr = poly["G"][w][0][0] + poly["G"][w][1][1]
    if tr != (F(1) if not w else F(0)):
        bad += 1
    trS = poly["J"][w][0][0] + poly["J"][w][1][1]
    if trS != F(0):  # Tr_A S has no identity term? check: S = Q - H, Tr_A H = 0, Tr_A Q = I
        # identity word: Tr_A S should be Tr_A Q - Tr_A H = 1 - 0 = 1
        if not w and trS != F(1):
            bad += 1
        elif w and trS != F(0):
            bad += 1
assert bad == 0, f"{bad} coefficient/trace mismatches"
print("PASS Q - S = H exactly (H rebuilt from the operator formula), over all", len(allw), "words")
print("PASS Tr_A Q = I_M (identity word 1, every other word 0)")
print("INDEPENDENT VERIFICATION COMPLETE: the t=1 Gram witnesses are PSD and the identities hold.")
