"""Different exact reconstruction route using python-flint. Same author; not external review."""
import json
from pathlib import Path
from flint import fmpq_mat,fmpq
from collections import defaultdict
D=json.loads((Path(__file__).resolve().parent.parent/'certificates'/'sharp_recovery.json').read_text());assert D['t']=='1'
words=D['words'];n=len(words)
# Build the full 106x106 congruence, then accumulate scalar coefficients.
P=[]
for w in words:
 a=fmpq_mat([[1,0],[0,1]])
 for t in w:a=a*fmpq_mat([[0,1],[1,0]] if abs(t)==1 else [[1,0],[0,-1]])
 P.append(a)
poly={}
for k in ['G','J']:
 B=fmpq_mat(D['basis_'+k]);R=fmpq_mat(D['reduced_'+k]);assert R==R.transpose()
 # Sylvester criterion via exact leading principal determinants.
 for size in range(1,R.nrows()+1):assert fmpq_mat([[R[i,j] for j in range(size)] for i in range(size)]).det()>0
 T=fmpq_mat(2*n,2*B.ncols())
 for i in range(n):
  for j in range(B.ncols()):
   for a in range(2):
    for b in range(2):T[2*i+a,2*j+b]=B[i,j]*P[i][a,b]
 RR=fmpq_mat(2*R.nrows(),2*R.ncols())
 for i in range(R.nrows()):
  for j in range(R.ncols()):
   for a in range(2):RR[2*i+a,2*j+a]=R[i,j]
 A=T*RR*T.transpose();out=defaultdict(lambda:[fmpq(0)]*4)
 for i,u in enumerate(words):
  for j,v in enumerate(words):
   w=[-x for x in reversed(u)]+v
   changed=True
   while changed:
    changed=False
    for z in range(len(w)-1):
     if w[z]==-w[z+1]:del w[z:z+2];changed=True;break
   for a in range(2):
    for b in range(2):out[tuple(w)][2*a+b]+=A[2*i+a,2*j+b]
 poly[k]=out
 print('PASS',k,'Sylvester determinants and full block reconstruction',flush=True)
H={(1,):[0,1,1,0],(-1,):[0,1,1,0],(2,):[1,0,0,-1],(-2,):[1,0,0,-1],(-1,2):[0,-1,1,0],(-2,1):[0,1,-1,0]}
for w in set(poly['G'])|set(poly['J'])|set(H):
 g,j=poly['G'][w],poly['J'][w]
 assert [3*(a-b) for a,b in zip(g,j)]==H.get(w,[0]*4)
 assert g[0]+g[3]==(1 if not w else 0)
print('PASS alternate exact t=1 certificate verification')
