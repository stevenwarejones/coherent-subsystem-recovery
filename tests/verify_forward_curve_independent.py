"""Independent re-derivation of the exact forward-circuit curve and the controller-only
minimax obstruction. Imports nothing from the shipped checkers.

Forward curve: the identity L^dag L = I + B^dag T + T^dag B over random unitaries; the
resulting r >= h(p) = max(0, 9p-1)^2/64; that h - (9f-5)/4 = 81(f-1)^2/64 (the affine
bound is its tangent); and the elementary->two-query switch f_c = (11+4 sqrt7)/27.

Controller obstruction: the 8x8 dual slack D_x = I_o (x) Q_x - rhobar^T/2 has the exact
spectrum {0, (x+1)/8, x/4, (x+1)(2x-1)/8} (each doubled), PSD on [1/2,1]; Tr(Q_x)=(1+x)^2/4;
and (3f-1)/2 - h(f) = 3(1-f)(27f-11)/64 > 0 on [7/9,1).
"""
import numpy as np
import sympy as sp

if not __debug__:
    raise SystemExit("Run without -O: this checker uses assertions.")

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], complex)
Z = np.array([[1, 0], [0, -1]], complex)


def ru(d, rng):
    A = rng.normal(size=(d, d)) + 1j * rng.normal(size=(d, d))
    q, r = np.linalg.qr(A)
    return q @ np.diag(np.exp(1j * np.angle(np.diag(r))))


rng = np.random.default_rng(31)
worst_id = 0.0
worst_gap = 1e9
for _ in range(300):
    d = int(rng.integers(1, 6))
    U, V = ru(d, rng), ru(d, rng)
    B = np.kron(X, U); Cc = np.kron(Z, V); Idm = np.eye(2 * d)
    L = Idm + B + Cc; T = (Idm + B) @ (Idm + Cc)
    worst_id = max(worst_id, np.max(np.abs(L.conj().T @ L - (Idm + B.conj().T @ T + T.conj().T @ B))))
    G = rng.normal(size=(2 * d, 2 * d)) + 1j * rng.normal(size=(2 * d, 2 * d))
    rho = G @ G.conj().T; rho /= np.trace(rho)
    p = np.real(np.trace(rho @ (L.conj().T @ L))) / 9
    r = np.real(np.trace(rho @ (T.conj().T @ T))) / 16
    worst_gap = min(worst_gap, r - max(0, 9 * p - 1) ** 2 / 64)
assert worst_id < 1e-9, ("E-FWDCURVE-IDENTITY", worst_id)
assert worst_gap > -1e-9, ("E-FWDCURVE-BOUND", worst_gap)
print(f"FWD: max |L†L-(I+B†T+T†B)| = {worst_id:.2e}; min (r - h(p)) = {worst_gap:.2e} (>=0)")

f = sp.symbols('f', real=True)
h = (9 * f - 1) ** 2 / 64
assert sp.factor(h - (9 * f - 5) / 4) == sp.Rational(81, 64) * (f - 1) ** 2, "E-FWDCURVE-TANGENT"
# Only the root in [1/9,1] governs decoder selection; the smaller root of 3f/4=(9f-1)^2/64
# lies below 1/9, where the true (truncated) h is 0 and there is no crossing.
fc = [sp.nsimplify(sv) for sv in sp.solve(sp.Eq(3 * f / 4, h), f) if sv >= sp.Rational(1, 9)]
assert fc == [sp.nsimplify((11 + 4 * sp.sqrt(7)) / 27)], ("E-FWDCURVE-SWITCH", fc)
print(f"FWD: h-affine = 81(f-1)^2/64; decoder switch f_c = {fc[0]}")

# ---- controller obstruction ----
x = sp.symbols('x', real=True)
kron = lambda a, b: sp.Matrix(sp.kronecker_product(a, b))
sI, sX, sZ, sY = sp.eye(2), sp.Matrix([[0, 1], [1, 0]]), sp.Matrix([[1, 0], [0, -1]]), sp.Matrix([[0, -sp.I], [sp.I, 0]])
E = (1 + x) / 2 * sp.eye(4) + (1 - x) / 2 * kron(sY, sY)
O = kron(sX, sX) + kron(sZ, sZ)
ket = lambda b: sp.Matrix([1 if i == b else 0 for i in range(2)])
phi = (kron(ket(0), ket(0)) + kron(ket(1), ket(1))) / sp.sqrt(2)
Phi_Aa = phi * phi.T
rho0 = sp.Matrix(sp.kronecker_product(Phi_Aa, sp.eye(2) / 2))          # A,a,c (8-dim)
IAE = sp.Matrix(sp.kronecker_product(sp.eye(2), E))
IAO = sp.Matrix(sp.kronecker_product(sp.eye(2), O))
rhobar = IAE * rho0 * IAE + (1 - x ** 2) / 4 * (IAO * rho0 * IAO)
Qx = (1 + x) ** 2 / 16 * sp.eye(4) + (1 - x ** 2) / 16 * kron(sY, sY)
Dx = sp.simplify(sp.Matrix(sp.kronecker_product(sp.eye(2), Qx)) - rhobar.T / 2)
assert Dx == Dx.conjugate().T, "E-CTRL-HERMITIAN"
want = {sp.expand((x + 1) / 8): 2, sp.expand(x / 4): 2,
        sp.expand((x + 1) * (2 * x - 1) / 8): 2, sp.Integer(0): 2}
got = {sp.expand(k): v for k, v in Dx.eigenvals().items()}
assert got == want, ("E-CTRL-SPECTRUM", got)
assert sp.simplify(sp.trace(Qx) - (1 + x) ** 2 / 4) == 0, "E-CTRL-TRACE"
f_of_x = (5 + 4 * x) / 9
assert sp.simplify(((9 * f_of_x - 1) ** 2 / 64) - (1 + x) ** 2 / 4) == 0, "E-CTRL-HF"
gap = sp.simplify((3 * f_of_x - 1) / 2 - (1 + x) ** 2 / 4)
assert sp.simplify(gap - 3 * (1 - f_of_x) * (27 * f_of_x - 11) / 64) == 0, "E-CTRL-GAP"
# on [7/9,1): 27f-11 >= 10 > 0 and 1-f > 0, so gap > 0
assert (27 * sp.Rational(7, 9) - 11) == 10
print("CTRL: D_x spectrum {0,(x+1)/8,x/4,(x+1)(2x-1)/8} (x2), PSD on [1/2,1]; gap = 3(1-f)(27f-11)/64 > 0 on [7/9,1)")
print("PASS independent: forward-curve identity+bound+tangent+switch, and controller-only minimax obstruction.")
