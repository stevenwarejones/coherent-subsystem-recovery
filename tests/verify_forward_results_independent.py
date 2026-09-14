"""Independent FreeGroup reconstruction and FLINT characteristic-polynomial checks.
Imports no supplied verifier or discovery code. Reads only certificate data.
Requires the optional SymPy/python-flint verification dependencies.
"""
from pathlib import Path
from itertools import product
import json
import sympy as s
from sympy.combinatorics.free_groups import free_group
from flint import fmpq_mat
ROOT=Path(__file__).resolve().parents[1]
if not __debug__: raise SystemExit('Run without -O')
def require(ok,msg):
    if not ok: raise RuntimeError(msg)
FG,u,v=free_group('u,v');e=FG.identity;one={e:s.Rational(1)};zero={}
def add(*polys):
    out={}
    for p in polys:
        for w,c in p.items(): out[w]=out.get(w,0)+c
    return {w:s.expand(c) for w,c in out.items() if s.expand(c)!=0}
def scale(a,p): return {w:a*c for w,c in p.items() if a*c!=0}
def mul(*polys):
    out=one
    for p in polys: out=add(*({w*t:c*d} for w,c in out.items() for t,d in p.items()))
    return out
def adj(p): return {w**-1:s.conjugate(c) for w,c in p.items()}
def positive(M,strict=False):
    # Real symmetric M is PSD iff det(tI+M) has nonnegative coefficients:
    # any negative eigenvalue would produce a positive root, impossible here.
    require(M==M.T,'independent symmetry')
    n=M.rows
    cp=fmpq_mat([[str(x) for x in M.row(i)] for i in range(n)]).charpoly()
    coeff=[(-1)**(n-k)*c for k,c in enumerate(cp)]
    require(all(c>=0 for c in coeff),'characteristic positivity')
    zeros=next(k for k,c in enumerate(coeff) if c!=0)
    require(not strict or zeros==0,'strict positivity')
    return n-zeros
# Negative controls exercise singular and indefinite characteristic cases.
require(positive(s.diag(0,1))==1,'characteristic nullity control')
for bad in (s.diag(-1,0),s.Matrix([[0,1],[1,0]])):
    try: positive(bad)
    except RuntimeError as exc: require(str(exc)=='characteristic positivity','negative control reason')
    else: raise RuntimeError('indefinite matrix accepted')

def load(path): return json.loads((ROOT/path).read_text())
def matrix(rows): return s.Matrix([[s.Rational(x) for x in row] for row in rows])
def word(letters):
    out=e
    for letter in letters: out*={1:u,2:v}[letter]
    return out
P=scale(s.Rational(1,2),add(one,{v:1}));N=scale(s.Rational(1,2),add(one,{v:-1}))
Rs=[(P,mul({u:1},N)),(mul({u:1},P),N)]
E=[[scale(s.Rational(1,4),add(*(mul(adj(R[a]),R[b]) for R in Rs))) for b in range(2)] for a in range(2)]
L=[[add(one,{v:1}),{u:1}],[{u:1},add(one,{v:-1})]]
K=[[scale(s.Rational(1,9),add(*(mul(adj(L[t][a]),L[t][b]) for t in range(2)))) for b in range(2)] for a in range(2)]
outcomes=[mul(last,{u:1},first) for first,last in product((P,N),repeat=2)]
Q=add(*(mul(adj(outcomes[i]),outcomes[i]) for i in (0,3)))
require(add(*(mul(adj(B),B) for B in outcomes))==one,'four calibration outcomes sum to I')
for R in Rs: require(add(*(mul(adj(B),B) for B in R))==one,'recovery trace preservation')
joint=load('certificates/joint_certificate.json');W=[e,v,u,u*v,v*u,v*v,v*u*v]
require([word(w) for w in joint['words']]==W,'joint words')
G=matrix(joint['gram']);require(positive(G)==12,'joint rank')
for a,b in product(range(2),repeat=2):
    rebuilt=add(*({W[i]**-1*W[j]:G[a*7+i,b*7+j]} for i,j in product(range(7),repeat=2)))
    target=add(E[a][b],scale(-s.Rational(5,4),K[a][b]),scale(3,Q) if a==b else zero,scale(s.Rational(1,4),one) if a==b else zero)
    require(rebuilt==target,'independent free-group joint SOS')
def moments(M,words):
    require(positive(M,True)==2*len(words),'moment rank')
    ix={w:i for i,w in enumerate(words)};m=len(words)
    require(M[ix[e],ix[e]]+M[m+ix[e],m+ix[e]]==1,'moment normalization')
    for g in (u,v):
        domain=[w for w in words if g*w in ix]
        for a,b in product(range(2),repeat=2):
            for x,y in product(domain,repeat=2):
                require(M[a*m+ix[x],b*m+ix[y]]==M[a*m+ix[g*x],b*m+ix[g*y]],'unitary Gram consistency')
    return ix
def expect(blocks,M,words):
    ix={w:i for i,w in enumerate(words)};m=len(words);value=0
    for a,b in product(range(2),repeat=2):
        for w,c in blocks[a][b].items():
            left,right=e,e;switched=False
            for symbol,power in w.array_form:
                gen=u if str(symbol)=='u' else v
                for _ in range(abs(power)):
                    if power<0:
                        require(not switched,'hereditary support');left=gen*left
                    else: switched=True;right*=gen
            require(left in ix and right in ix,'moment support')
            value+=c*M[a*m+ix[left],b*m+ix[right]]
    return s.factor(value)
ad=load('certificates/joint_adversary.json');M=matrix(ad['moment']);moments(M,W)
p,r,q=expect(K,M,W),expect(E,M,W),expect([[Q,zero],[zero,Q]],M,W)
require((p,r,q)==(s.Rational(1265617,1406250),s.Rational(424947,500000),s.Rational(58019,5000000)),'adversary values')
require((3*p-1)/2-r==s.Rational(731,7500000),'joint gap')
kappa=q/(1-p);require(kappa==s.Rational(522171,4500256) and kappa<s.Rational(1,8),'joint threshold')
X=s.Matrix([[0,1],[1,0]]);Z=s.diag(1,-1)
def histories(order):
    out=[]
    for bits in product((0,1),repeat=len(order)):
        w=e
        for gate,bit in zip(order,bits):
            if bit: w={'U':u,'V':v}[gate]*w
        out.append(w)
    return out
def pauli(w):
    out=s.eye(2)
    for symbol,power in w.array_form: out=out*({'u':X,'v':Z}[str(symbol)]**power)
    return out
def objectives(order,M):
    hs=histories(order);words=list(dict.fromkeys(hs));ix=moments(M,words)
    A=s.Matrix.hstack(*(pauli(w)[:,a] for a in range(2) for w in hs))
    repeated=[a*len(words)+ix[w] for a in range(2) for w in hs]
    return hs,words,A.T*A/4,M.extract(repeated,repeated)/2

def dual_operator(levels,n):
    # Adjoint of repeated equal-bit sums: tensor each cross-bit multiplier
    # with identity on the already-summed suffix bits and output system.
    H=s.zeros(n)
    for depth,level in enumerate(levels):
        B=matrix(level);off=s.zeros(2*B.rows)
        for i,j in product(range(B.rows),repeat=2):
            off[2*i,2*j+1]=B[i,j]/2;off[2*j+1,2*i]=B[i,j]/2
        H+=s.kronecker_product(off,s.eye(2**depth))
    return s.kronecker_product(s.eye(2),H)
rob=load('certificates/forward_queries/robust_VUVUV.json');M=matrix(rob['moment'])
hs,words,C0,Ca=objectives('VUVUV',M);pstar=expect(K,M,words)
require(pstar==s.Rational(1507007,2254500),'robust support score')
theta=s.Rational(1,10)/(1-pstar);alpha=s.Rational(1692197,2000000)
residual=alpha*s.eye(64)+dual_operator(rob['levels'],32)-((1-theta)*C0+theta*Ca)
positive(residual,True)
require(theta==s.Rational(225450,747493),'robust mixture')
require(10*(1-alpha)==s.Rational(307803,200000),'robust slope')
require(s.Rational(17,20)-alpha==s.Rational(7803,2000000),'robust deficit')
for count in range(1,4):
    for order in product('UV',repeat=count):
        if len(set(order))<2: continue
        cursor=iter('VUVUV')
        require(all(any(g==wanted for g in cursor) for wanted in order),'three-call padding')
# Independent wrong-Bell kernel, with dual lifted before kernel restriction.
for order in ('VUU','VUV','VUVUV','VVU'):
    data=load(f'history/endpoint_forward_queries/certificate_{order}.json')
    hs,words,C0,Ca=objectives(order,matrix(data['moment']));n=len(hs);rows=[]
    for B in (X,Z,X*Z):
        for target in range(2):
            rows.append([sum(B[a,b]*pauli(w)[target,a] for a in range(2)) for b in range(2) for w in hs])
    basis=s.Matrix.hstack(*s.Matrix(rows).nullspace());require(basis.cols>0,'endpoint kernel')
    vals=[s.Rational(x) for x in data['dual']];levels=[];offset=0;size=n
    while size>1:
        count=(size//2)**2
        levels.append(s.Matrix(size//2,size//2,vals[offset:offset+count]).tolist())
        offset+=count;size//=2
    require(offset+1==len(vals),'historical dual length')
    H=vals[-1]*s.eye(2*n)+dual_operator(levels,n)-2*Ca
    positive(basis.T*H*basis,True)
    gap=vals[-1]-(3*expect(K,matrix(data['moment']),words)-1)
    require(gap==s.Rational(data['upper_2r_minus_H']) and gap<0,'historical gap')
# Independently check the reflection identity in the quotient u unitary, v^2=I,
# instead of scalar/block substitution. No commuting-variable simplification.
def reflection(p):
    result={}
    for w,c in p.items():
        stack=[]
        for symbol,power in w.array_form:
            for _ in range(abs(power)):
                letter=2 if str(symbol)=='v' else (1 if power>0 else -1)
                if stack and ((letter==2 and stack[-1]==2) or (letter!=2 and stack[-1]==-letter)):
                    stack.pop()
                else: stack.append(letter)
        key=tuple(stack);result[key]=result.get(key,0)+c
    return {w:s.expand(c) for w,c in result.items() if s.expand(c)!=0}
D=add(mul(P,{u:1},P),mul(N,{u**-1:1},N))
for a,b in product(range(2),repeat=2):
    diff=add(E[a][b],scale(-s.Rational(9,8),K[a][b]),scale(s.Rational(1,8),one) if a==b else zero)
    correction=scale(s.Rational(1,4),D if a==0 else adj(D)) if a!=b else zero
    require(not reflection(add(diff,correction)),'independent reflection quotient identity')
# Full free-group Pauli covariance. This checks also general non-reflection V.
def signs(p,sx,sz):
    return {w:c*s.prod((sx if str(g)=='u' else sz)**power for g,power in w.array_form) for w,c in p.items()}
for g,sx,sz in [(s.eye(2),1,1),(X,1,-1),(Z,-1,1),(s.I*X*Z,-1,-1)]:
    for B in (E,K):
        for a,b in product(range(2),repeat=2):
            transformed=add(*(scale(g[a,i]*s.conjugate(g[b,j]),B[i][j]) for i,j in product(range(2),repeat=2)))
            require(signs(B[a][b],s.Integer(sx),s.Integer(sz))==transformed,'independent Pauli covariance')
    require(signs(Q,s.Integer(sx),s.Integer(sz))==Q,'calibration sign invariance')

# Independent exact coherent-gate calculation of the restricted phase family.
zeta,omega=s.symbols('zeta omega')
Had=(X+Z)/s.sqrt(2)
initial=s.Matrix([1,0,1,0,0,1,0,1])/2
def controlled(A,bit):
    return s.kronecker_product(s.eye(2),s.diag(A,s.eye(2)) if bit==0 else s.diag(s.eye(2),A))
hgate=s.kronecker_product(s.eye(2),Had,s.eye(2))
Uphase=zeta*X;Vphase=omega*Z
base=controlled(Uphase,1)*hgate*controlled(Vphase,1)*initial
echo=controlled(Vphase,0)*base
swapped=hgate*controlled(Uphase,0)*controlled(Vphase,1)*hgate*controlled(Uphase,1)*initial
scores=[]
for state in (base,echo,swapped):
    amplitudes=[(state[j]+state[6+j])/s.sqrt(2) for j in range(2)]
    value=sum(s.conjugate(a)*a for a in amplitudes).subs({s.conjugate(zeta):1/zeta,s.conjugate(omega):1/omega})
    scores.append(s.simplify(value))
a=(zeta+1/zeta)/2;b=(omega+1/omega)/2;c=(zeta/omega+omega/zeta)/2
for actual,expected in zip(scores,[(1+a)*(1+b)/4,(1+b)*(1+c)/4,(1+a)*(1+c)/4]):
    require(s.simplify(actual-expected)==0,'independent phase gate fidelity')
require(s.simplify(sum(scores)/3-(a+b+c)/3-((1-a)*(1-b)+(1-a)*(1-c)+(1-b)*(1-c))/12)==0,'phase margin')
# Direct completion examples, independent of check_forward_reductions.py.
phi=s.Matrix([1,0,0,1])/s.sqrt(2);Phi=phi*phi.T
Lp=s.eye(4)+s.kronecker_product(X,X)+s.kronecker_product(Z,Z)
Hp=Lp.T*Lp/3-s.eye(4)
require(Hp==s.Rational(8,3)*Phi-s.Rational(2,3)*s.eye(4),'ideal completion premise')
require(2*Phi-Hp==s.Rational(2,3)*(s.eye(4)-Phi),'ideal completion slack')
Ls=s.eye(4)+s.kronecker_product(X,s.eye(2))+s.kronecker_product(Z,s.I*s.eye(2))
Hs=s.simplify(Ls.conjugate().T*Ls/3-s.eye(4))
require(Hs[:2,2:]==2*(1-s.I)*s.eye(2)/3,'scalar offdiagonal')
forced=s.Rational(3,4)*Hs[0,2]
require(s.simplify(forced*s.conjugate(forced))==s.Rational(1,2),'normalized PSD violation')
require(s.simplify(Hs*Hs-s.Rational(8,9)*s.eye(4))==s.zeros(4),'scalar feasible completion')

# Validate shot-count rounding with exact integer powers at both cutoffs.
for den,count in ((120,358),(1200,3594),(12000,35948)):
    require(20*pow(den-1,count)<=pow(den,count),'sufficient zero-event count')
    require(20*pow(den-1,count-1)>pow(den,count-1),'minimal zero-event count')
num,den=map(int,s.fraction(kappa/100));count=2537
require(19*pow(den-num,count)<=pow(den,count),'necessary sample cutoff')
require(19*pow(den-num,count-1)>pow(den,count-1),'necessary sample preceding integer')
print('PASS independent FreeGroup recovery/calibration normalization and joint SOS; characteristic rank 12')
print('PASS independent moments, joint values, robust full dual and four wrong-Bell endpoint duals')
print('PASS exhaustive mixed three-call padding and exact integer sample-count cutoffs')
print('PASS independent reflection quotient identity and Pauli covariance for E, K, Q')
print('PASS independent phase-family gate amplitudes and fixed-scaling completion examples')
