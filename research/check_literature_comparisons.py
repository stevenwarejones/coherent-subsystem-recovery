"""Exact literature-comparison arithmetic for docs/PRIOR_ART.md.

Every claim below is ASSERTED, not merely printed: tampering with the numbers makes the
script fail (see tests/test_supporting_mutations.py). Scope: this checks the comparison
arithmetic and the two counterexamples exactly; it does NOT verify the sharp certificate
(that is proofs/verify_sharp_recovery.py) and it is not a novelty audit.
"""
if not __debug__:
    raise SystemExit("Run without -O: this checker uses assertions.")
import sympy as sp

X = sp.Matrix([[0, 1], [1, 0]])
Z = sp.Matrix([[1, 0], [0, -1]])
Y = sp.Matrix([[0, -sp.I], [sp.I, 0]])
I2 = sp.eye(2)


def ptrace_M(M4):
    """Trace out the second qubit (M) of a 4x4 operator, returning the 2x2 on A."""
    return sp.Matrix(2, 2, lambda i, j: sum(M4[2 * i + m, 2 * j + m] for m in range(2)))


# ---- 1. Berta-Coles-Wehner failed substitution: f <= P_pg is FALSE -------------------
ket = lambda b: sp.Matrix([1 if i == b else 0 for i in range(2)])
phi = (sp.kronecker_product(ket(0), ket(0)) + sp.kronecker_product(ket(1), ket(1))) / sp.sqrt(2)
Phi = phi * phi.T
rho = sp.Rational(3, 4) * Phi + sp.Rational(1, 4) * sp.eye(4) / 4
L = sp.kronecker_product(I2, I2) + sp.kronecker_product(X, X) + sp.kronecker_product(Z, Z)
p = sp.nsimplify(sp.trace(L.T * L * rho) / 9)
C = ptrace_M((L * rho * L.T) / 9)
Rpg = sp.nsimplify(sp.trace(rho * rho))
Ppg = sp.nsimplify((2 * Rpg + 1) / 3)
fopt = sp.nsimplify(sp.Rational(1, 2) * (sp.trace(sp.sqrt(C))) ** 2)   # C scalar => Uhlmann

assert p == sp.Rational(5, 6), ("p", p)
assert C == sp.Rational(5, 12) * I2, ("C should be (5/12) I", C)
assert fopt == sp.Rational(5, 6), ("f_opt", fopt)
assert Rpg == sp.Rational(43, 64) and Ppg == sp.Rational(25, 32), ("E-BCW-VALUES", Rpg, Ppg)
assert fopt - Ppg == sp.Rational(5, 6) - sp.Rational(25, 32) > 0, "f_opt must strictly exceed P_pg"
print("PASS BCW: p=f_opt=5/6, C=(5/12)I, R_pg=43/64, P_pg=25/32 < 5/6, so f <= P_pg is FALSE")

# ---- 2. sharp bound strictly beats BCW+Renes2017 on [2/3,1), proven algebraically -----
f = sp.symbols('f', real=True)
s = (1 + 3 * f) / 4
bcw = 3 * s ** 2 - 3 * s + 1
sharp = (3 * f - 1) / 2
gap = sp.simplify(sharp - bcw)
# Prove positivity on (2/3, 1): substitute f = 1 - e, e in (0, 1/3); gap = 3 e (4 - 9 e)/16.
e = sp.symbols('e', positive=True)
gap_e = sp.simplify(gap.subs(f, 1 - e))
assert sp.simplify(gap_e - 3 * e * (4 - 9 * e) / 16) == 0, ("gap(1-e) form", gap_e)
# For 0 < e <= 1/3: 4 - 9 e >= 1 > 0, so gap_e > 0 strictly.
assert sp.simplify((4 - 9 * sp.Rational(1, 3))) == 1 > 0
# leading-order coefficients differ (3/2 vs 9/4): NOT "agreement to leading order"
lin_sharp = sp.series(sharp.subs(f, 1 - e), e, 0, 2).removeO().coeff(e, 1)
lin_bcw = sp.series(bcw.subs(f, 1 - e), e, 0, 2).removeO().coeff(e, 1)
assert lin_sharp == sp.Rational(-3, 2) and lin_bcw == sp.Rational(-9, 4), (lin_sharp, lin_bcw)
# endpoints of the printed table
for fv, a, b in [(sp.Rational(2, 3), sp.Rational(1, 2), sp.Rational(7, 16)),
                 (1, sp.Integer(1), sp.Integer(1))]:
    assert sharp.subs(f, fv) == a and bcw.subs(f, fv) == b, (fv,)
print("PASS comparison: sharp - (BCW+Renes) = 3e(4-9e)/16 > 0 on (2/3,1); slopes 3/2 vs 9/4 differ")

# ---- 3. R >= S/3 FAILS for three independent observables (fully constructed) ----------
# c-q separable state: reference qubit rho0, flag g in {I,X,Y,Z} stored classically in M.
rho0 = sp.Rational(1, 2) * (I2 + (X + Y + Z) / sp.sqrt(3))
paulis = [I2, X, Y, Z]
# rho_AM = (1/4) sum_g (g rho0 g^dag) (x) |g><g|_M ; A is the reference qubit, M the 4-flag.
rhoAM = sp.zeros(8, 8)
for g in range(4):
    blk = paulis[g] * rho0 * paulis[g].conjugate().T
    flag = sp.zeros(4, 4); flag[g, g] = 1
    rhoAM += sp.kronecker_product(blk, flag) / 4
# reference marginal must be maximally mixed (Pauli twirl of any state)
rhoA = sp.Matrix(2, 2, lambda i, j: sum(rhoAM[i * 4 + m, j * 4 + m] for m in range(4)))
assert sp.simplify(rhoA - I2 / 2) == sp.zeros(2, 2), ("rho_A", rhoA)
# Three diagonal memory observables returning each Pauli's conjugation sign under the flag:
# g P g^dag = s_P(g) P, with P in {X, Z, Y}. Correlation <P (x) O_P> where O_P = diag(s_P(g)).
signs = {"X": [1, 1, -1, -1], "Z": [1, -1, -1, 1], "Y": [1, -1, 1, -1]}
refP = {"X": X, "Z": Z, "Y": Y}
corrs = {}
for name in ("X", "Z", "Y"):
    O = sp.diag(*signs[name])                       # diagonal memory observable on the 4-flag
    obs = sp.kronecker_product(refP[name], O)       # reference Pauli (x) memory observable
    corrs[name] = sp.nsimplify(sp.trace(rhoAM * obs))
    assert corrs[name] == 1 / sp.sqrt(3), ("E-STEERING-CORRELATION", name, corrs[name])
Ssum = sp.nsimplify(sum(corrs.values()))
assert Ssum == sp.sqrt(3), ("S", Ssum)
# The state is separable (block-diagonal c-q), so best Bell overlap (singlet fraction) is 1/2.
# Hence R = 1/2 < S/3 = 1/sqrt(3) ~ 0.577, so no universal R >= S/3 for independent observables.
assert sp.Rational(1, 2) < Ssum / 3, "R=1/2 must be below S/3 here"
print("PASS steering c/e: constructed flagged model gives each corr = 1/sqrt3, S = sqrt3,")
print("     rho_A = I/2, separable so R = 1/2 < S/3 -- R >= S/3 fails for independent observables")

# ---- 4. Coherent control: relative phases exceed the individual channel descriptions --
# Branches (I, X, Z) vs (I, -X, Z) induce identical individual channels but differ coherently.
def C_from(U, V, state4):
    Lop = sp.kronecker_product(I2, I2) + sp.kronecker_product(X, U) + sp.kronecker_product(Z, V)
    return ptrace_M((Lop * state4 * Lop.conjugate().T) / 9)


C_plus = C_from(X, Z, Phi)
C_minus = C_from(-X, Z, Phi)
assert C_plus == I2 / 2, ("C(+X)", C_plus)
assert C_minus == I2 / 18, ("C(-X)", C_minus)
fopt_plus = sp.nsimplify(sp.Rational(1, 2) * sp.trace(sp.sqrt(C_plus)) ** 2)
fopt_minus = sp.nsimplify(sp.Rational(1, 2) * sp.trace(sp.sqrt(C_minus)) ** 2)
assert fopt_plus == 1 and fopt_minus == sp.Rational(1, 9), (fopt_plus, fopt_minus)
print("PASS coherent control: (I,X,Z) vs (I,-X,Z) give C = I/2 vs I/18, f_opt = 1 vs 1/9,")
print("     though the individual channels coincide -- relative phases are physical")

print("\nALL literature-comparison assertions passed (exact arithmetic; not a novelty audit).")
