import sympy as sp
# ---- BCW failed-substitution counterexample ----
s2=sp.sqrt(2)
ket=lambda b:sp.Matrix([1 if i==b else 0 for i in range(2)])
phi=(sp.kronecker_product(ket(0),ket(0))+sp.kronecker_product(ket(1),ket(1)))/s2
Phi=phi*phi.T                                   # |Phi+><Phi+| (real)
I4=sp.eye(4)
rho=sp.Rational(3,4)*Phi+sp.Rational(1,4)*I4/4
X=sp.Matrix([[0,1],[1,0]]); Z=sp.Matrix([[1,0],[0,-1]]); I2=sp.eye(2)
L=sp.kronecker_product(I2,I2)+sp.kronecker_product(X,X)+sp.kronecker_product(Z,Z)
p=sp.trace(L.T*L*rho)/9
# C = Tr_M( L rho L^dag )/9 : trace out second qubit (M)
LrL=(L*rho*L.T)/9
C=sp.zeros(2,2)
for i in range(2):
  for j in range(2):
    C[i,j]=sum(LrL[2*i+m,2*j+m] for m in range(2))
Rpg=sp.trace(rho*rho)
Ppg=(2*Rpg+1)/3
fopt=sp.Rational(1,2)*(sp.trace(sp.sqrt(C)))**2   # C scalar => Uhlmann optimum
print("p      =",sp.nsimplify(p),"=",float(p))
print("C      =",C.tolist(),"  (should be (5/12) I)")
print("f_opt  =",sp.nsimplify(fopt),"=",float(fopt))
print("R_pg   =",Rpg,"  P_pg =",sp.nsimplify(Ppg),"=",float(Ppg))
print("f_opt > P_pg ?", sp.simplify(fopt-Ppg)>0, " => 'f <= P_pg' is FALSE as claimed")
# ---- second counterexample: R>=S/3 fails for independent observables ----
Y=sp.Matrix([[0,-sp.I],[sp.I,0]])
rho0=sp.Rational(1,2)*(I2+(X+Y+Z)/sp.sqrt(3))
corr=sp.Rational(1,1)/sp.sqrt(3)
S=3*corr
print("\nsteering c/e: each correlation =",corr,"  S =",sp.nsimplify(S),
      " S/3 =",float(S/3)," vs R=1/2 -> R < S/3 ?", float(S/3)>0.5)
# ---- comparison table: sharp vs BCW+Renes2017 ----
f=sp.symbols('f'); ss=(1+3*f)/4; bcw=3*ss**2-3*ss+1; sharp=(3*f-1)/2
print("\n f      sharp      BCW+Renes   diff")
for fv in [sp.Rational(2,3),sp.Rational(4,5),sp.Rational(9,10),sp.Rational(95,100),sp.Rational(99,100),1]:
  a=sharp.subs(f,fv); b=bcw.subs(f,fv)
  print(f"{float(fv):.4f}  {float(a):.6f}  {float(b):.6f}   {float(a-b):+.6f}")
eps=sp.symbols('epsilon',positive=True)
diff=sp.series(sharp.subs(f,1-eps),eps,0,3).removeO()-sp.series(bcw.subs(f,1-eps),eps,0,3).removeO()
print("\nsharp - (BCW+Renes) at f=1-eps :",sp.simplify(diff)," ; audit says 3*eps*(4-9*eps)/16 =",sp.expand(3*eps*(4-9*eps)/16))
