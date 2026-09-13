"""Exact robust recovery upper-bound checker. No endpoint restriction or solver.
The numerical certificate generator is not imported.
"""
from fractions import Fraction as F
from fractions import Fraction as Q
from pathlib import Path
import itertools,json,copy,sys
if not __debug__:raise SystemExit("Run without -O")
def need(test,code):
 if not test: raise ValueError(code)
def zero(n,m=None):return [[Q(0)for _ in range(n if m is None else m)]for _ in range(n)]
def eye(n):
 z=zero(n)
 for i in range(n):z[i][i]=Q(1)
 return z
def tr(a):return list(map(list,zip(*a)))
def mul(a,b):
 out=zero(len(a),len(b[0])); nz=[[(j,x)for j,x in enumerate(row)if x]for row in b]
 for i,row in enumerate(a):
  for k,v in enumerate(row):
   if v:
    for j,x in nz[k]:out[i][j]+=v*x
 return out
def pd(a,code):
 n=len(a);need(all(a[i][j]==a[j][i]for i in range(n)for j in range(n)),code)
 a=copy.deepcopy(a)
 for k in range(n):
  p=a[k][k];need(p>0,code)
  for i in range(k+1,n):
   for j in range(i,n):a[j][i]=a[i][j]=a[i][j]-a[i][k]*a[j][k]/p

def words(order):
 return [tuple({'U':1,'V':2}[c]for c,b in list(zip(order,bits))[::-1]if b)for bits in itertools.product([0,1],repeat=len(order))]
def build(d,f=F(9,10)):
 hs=words(d['order']);W=list(dict.fromkeys(hs));N=len(hs);m=len(W);ix={w:i for i,w in enumerate(W)}
 G=[[F(x)for x in row]for row in d['moment']];X=[[F(0),F(1)],[F(1),F(0)]];Z=[[F(1),F(0)],[F(0),F(-1)]];XZ=mul(X,Z)
 A=zero(2,2*N);C=zero(2*N)
 for a in range(2):
  for h,w in enumerate(hs):
   M=eye(2)
   for x in w:M=mul(M,{1:X,2:Z}[x])
   for i in range(2):A[i][a*N+h]=M[i][a]
   for b in range(2):
    for k,v in enumerate(hs):C[b*N+h][a*N+k]=G[a*m+ix[v]][b*m+ix[w]]/2
 C0=[[x/4 for x in row]for row in mul(tr(A),A)]
 terms=[((),(1,),X,1),((1,),(),X,1),((),(2,),Z,1),((2,),(),Z,1),((1,),(2,),XZ,1),((2,),(1,),XZ,-1)]
 H=sum(s*M[a][b]*G[a*m+ix[w]][b*m+ix[v]]/3 for w,v,M,s in terms for a,b in itertools.product(range(2),repeat=2));p=(H+1)/3;theta=(1-f)/(1-p)
 need(0<theta<1,"MIXTURE_RANGE")
 Cmix=[[(1-theta)*x+theta*y for x,y in zip(r,s)]for r,s in zip(C0,C)]
 return G,Cmix,p,theta

def residual(C,levels,alpha):
 n=len(C)//2;N=n;stride=1;R=[[-x for x in row]for row in C]
 for i in range(2*N):R[i][i]+=alpha
 for dual in levels:
  need(len(dual)==n//2 and all(len(r)==n//2 for r in dual),"DUAL_SHAPE")
  for i in range(n//2):
   for j in range(n//2):
    v=dual[i][j]/2
    for a in range(2):
     for s in range(stride):
      h=a*N+2*i*stride+s;k=a*N+(2*j+1)*stride+s
      R[h][k]+=v;R[k][h]+=v
  n//=2;stride*=2
 need(n==1,"DUAL_DEPTH")
 return R


def verify(d):
 need(d['order']=='VUVUV','ORDER');need(F(d['score'])==F(9,10),'SCORE')
 hs=words(d['order']);W=list(dict.fromkeys(hs));m=len(W);ix={w:i for i,w in enumerate(W)}
 G=[[F(x)for x in r]for r in d['moment']]
 need(len(G)==2*m and all(len(r)==2*m for r in G),'MOMENT_SHAPE');pd(G,'MOMENT_PD')
 need(G[0][0]+G[m][m]==1,'NORMALIZATION')
 for x in (1,2):
  dom=[w for w in W if (x,)+w in ix]
  for a,b in itertools.product(range(2),repeat=2):
   for w,v in itertools.product(dom,repeat=2):
    need(G[a*m+ix[(x,)+w]][b*m+ix[(x,)+v]]==G[a*m+ix[w]][b*m+ix[v]],'UNITARY_MOMENTS')
 G,C,p,theta=build(d,F(d['score']))
 need(p==F(d['p']) and theta==F(d['theta']),'MIXTURE_IDENTITY')
 levels=[[[F(x)for x in r]for r in lev]for lev in d['levels']];alpha=F(d['alpha'])
 pd(residual(C,levels,alpha),'DUAL_PD')
 deficit=F(17,20)-alpha
 need(deficit==F(d['deficit']),'DEFICIT_IDENTITY');need(deficit>0,'POSITIVE_GAP')
 return alpha,(1-alpha)*10

def main():
 root=Path(__file__).resolve().parents[1]/'certificates'/'forward_queries';d=json.loads((root/'robust_VUVUV.json').read_text());alpha,beta=verify(d)
 tests=[]
 q=copy.deepcopy(d);q['moment'][0][0]='-1';tests.append((q,'MOMENT_PD'))
 q=copy.deepcopy(d);q['moment'][1][1]=str(F(q['moment'][1][1])+1);tests.append((q,'UNITARY_MOMENTS'))
 q=copy.deepcopy(d);q['theta']='1/2';tests.append((q,'MIXTURE_IDENTITY'))
 q=copy.deepcopy(d);q['alpha']=str(F(q['alpha'])-1);tests.append((q,'DUAL_PD'))
 q=copy.deepcopy(d);q['deficit']='1/10';tests.append((q,'DEFICIT_IDENTITY'))
 for q,reason in tests:
  try:verify(q)
  except ValueError as e:need(str(e)==reason,'WRONG_REJECTION_REASON')
  else:raise ValueError('MUTATION_ACCEPTED')
 print('PASS exact 40x40 moment realization hypotheses')
 print('PASS exact 64x64 full-controller dual is positive definite; no endpoint kernel')
 print(f'PASS at f=9/10, every allowed controller has r <= {alpha} = {float(alpha)}')
 print(f'PASS family consequence: r <= 1 - ({beta})(1-f), for 9/10 <= f <= 1')
 print('PASS all five mutations reject for intended reasons')
if __name__=='__main__':main()
