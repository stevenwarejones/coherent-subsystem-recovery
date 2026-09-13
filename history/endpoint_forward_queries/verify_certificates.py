"""Exact certificate checker. Standard library only; no solver or researcher imports.
Checks arithmetic, moment realizability hypotheses and the endpoint-kernel dual.
The unitary-extension and physical symmetrization arguments require prose review.
"""
from fractions import Fraction as Q
from pathlib import Path
import itertools,json,sys,copy
if not __debug__: raise SystemExit('Run without -O')

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

def nullspace(a):
 a=copy.deepcopy(a);n=len(a[0]);piv=[];r=0
 for c in range(n):
  p=next((i for i in range(r,len(a))if a[i][c]),None)
  if p is None:continue
  a[r],a[p]=a[p],a[r];v=a[r][c];a[r]=[x/v for x in a[r]]
  for i in range(len(a)):
   if i!=r:
    v=a[i][c];a[i]=[x-v*y for x,y in zip(a[i],a[r])]
  piv.append(c);r+=1
  if r==len(a):break
 free=[c for c in range(n)if c not in piv];b=zero(n,len(free))
 for k,c in enumerate(free):
  b[c][k]=Q(1)
  for i,p in enumerate(piv):b[p][k]=-a[i][c]
 return b

def words(order):
 return [tuple({'U':1,'V':2}[c]for c,b in list(zip(order,bits))[::-1]if b)for bits in itertools.product([0,1],repeat=len(order))]
def model(order):
 hs=words(order);W=list(dict.fromkeys(hs));N=len(hs)
 X=[[Q(0),Q(1)],[Q(1),Q(0)]];Z=[[Q(1),Q(0)],[Q(0),Q(-1)]]
 mats=[]
 for w in hs:
  a=eye(2)
  for x in w:a=mul(a,{1:X,2:Z}[x])
  mats.append(a)
 L=[[mul(v,p)[j][a]for a in range(2)for v in mats]for p in [X,Z,mul(X,Z)]for j in range(2)]
 B=nullspace(L);As=[];n=N;stride=1
 while n>1:
  for i in range(0,n,2):
   for j in range(1,n,2):
    A=zero(2*N)
    for a in range(2):
     for s in range(stride):
      h=a*N+i*stride+s;k=a*N+j*stride+s
      A[h][k]+=Q(1,2);A[k][h]+=Q(1,2)
    As.append(mul(mul(tr(B),A),B))
  n//=2;stride*=2
 As.append(mul(tr(B),B))
 return hs,W,B,As,X,Z

def verify(data):
 order=data['order'];need(order in ('VUV','VVU','VUU','VUVUV'),'ORDER')
 hs,W,B,As,X,Z=model(order);N=len(hs);m=len(W);ix={w:i for i,w in enumerate(W)}
 G=[[Q(x)for x in r]for r in data['moment']]
 need(len(G)==2*m and all(len(r)==2*m for r in G),'SHAPE');pd(G,'MOMENT_PD')
 need(G[0][0]+G[m][m]==1,'NORMALIZATION')
 for w in W:need(not w or w[1:]in ix,'SUFFIX_CLOSURE')
 for x in (1,2):
  dom=[w for w in W if (x,)+w in ix]
  for a,b in itertools.product(range(2),repeat=2):
   for w,v in itertools.product(dom,repeat=2):
    need(G[a*m+ix[(x,)+w]][b*m+ix[(x,)+v]]==G[a*m+ix[w]][b*m+ix[v]],'UNITARY_MOMENTS')
 D=zero(2*N)
 for a,b in itertools.product(range(2),repeat=2):
  for h,k in itertools.product(range(N),repeat=2):D[b*N+h][a*N+k]=G[a*m+ix[hs[k]]][b*m+ix[hs[h]]]
 XZ=mul(X,Z)
 terms=[((),(1,),X,1),((1,),(),X,1),((),(2,),Z,1),((2,),(),Z,1),((1,),(2,),XZ,1),((2,),(1,),XZ,-1)]
 H=sum(s*C[a][b]*G[a*m+ix[w]][b*m+ix[v]]/3 for w,v,C,s in terms for a,b in itertools.product(range(2),repeat=2))
 C=mul(mul(tr(B),D),B);lam=[Q(x)for x in data['dual']];need(len(lam)==len(As),'DUAL_SHAPE')
 R=[[sum(t*A[i][j]for t,A in zip(lam,As))-C[i][j]for j in range(len(C))]for i in range(len(C))]
 pd(R,'DUAL_PD');gap=lam[-1]-H
 need(gap==Q(data['upper_2r_minus_H']),'GAP_IDENTITY');need(gap<0,'NEGATIVE_BOUND')
 return gap

def main():
 root=Path(__file__).resolve().parent
 paths=[root/f'certificate_{order}.json' for order in ('VUU','VUV','VUVUV','VVU')]
 need(all(p.is_file() for p in paths),'MISSING_CERTIFICATE')
 messages=[]
 for p in paths:
  d=json.loads(p.read_text());g=verify(d);messages.append(f'PASS {d["order"]}: 2r-<H> <= {g} < 0; dimension <= {len(d["moment"])}')
  tests=[]
  x=copy.deepcopy(d);x['moment'][0][0]='-1';tests.append((x,'MOMENT_PD'))
  x=copy.deepcopy(d);x['moment'][1][1]=str(Q(x['moment'][1][1])+1);tests.append((x,'UNITARY_MOMENTS'))
  x=copy.deepcopy(d);x['dual'][-1]=str(Q(x['dual'][-1])-10);tests.append((x,'DUAL_PD'))
  x=copy.deepcopy(d);x['upper_2r_minus_H']='-1';tests.append((x,'GAP_IDENTITY'))
  for x,reason in tests:
   try:verify(x)
   except ValueError as e:need(str(e)==reason,'WRONG_REJECTION_REASON')
   else:raise ValueError('MUTATION_ACCEPTED')
 for message in messages: print(message)
 print('PASS all four mutation types reject for the intended reason, for each certificate')
if __name__=='__main__':main()
