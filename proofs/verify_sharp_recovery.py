"""Exact certificate verification, standard library only. No SDP or numerical tolerances.
The operational interpretation is stated separately in THEOREM.md.
"""
if not __debug__:raise RuntimeError('Verification requires assertions; do not use -O.')
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import json,sys

def mm(A,B):return [[sum((a*b for a,b in zip(row,col)),F(0)) for col in zip(*B)] for row in A]
def transpose(A):return [list(row) for row in zip(*A)]
def pd(A):
 n=len(A);assert n and all(len(row)==n for row in A),'SHAPE'
 assert A==transpose(A),'SYMMETRY'
 # Exact symmetric elimination; positive pivots prove positive definiteness.
 T=[row[:] for row in A]
 for i in range(n):
  p=T[i][i];assert p>0,'POSITIVE_DEFINITE'
  for j in range(i+1,n):
   for k in range(j,n):T[j][k]=T[k][j]=T[j][k]-T[j][i]*T[i][k]/p

def red(w):
 out=[]
 for x in w:
  if out and out[-1]==-x:out.pop()
  else:out.append(x)
 return tuple(out)
def dag(w):return tuple(-x for x in reversed(w))
I=[[1,0],[0,1]];X=[[0,1],[1,0]];Z=[[1,0],[0,-1]]
def pauli(w):
 A=I
 for x in w:A=mm(A,X if abs(x)==1 else Z)
 return A

def verify(path):
 d=json.loads(Path(path).read_text());assert F(d['t'])==1,'CLAIMED_CONSTANT'
 words=[tuple(w) for w in d['words']];n=len(words)
 assert n==53 and len(set(words))==n,'WORD_BASIS'
 assert all(len(w)<=3 and red(w)==w and all(x in [-2,-1,1,2] for x in w) for w in words),'WORD_BASIS'
 P=[pauli(w) for w in words];poly={}
 for key in ['G','J']:
  B=[[F(x) for x in row] for row in d['basis_'+key]];R=[[F(x) for x in row] for row in d['reduced_'+key]]
  assert len(B)==n and all(len(row)==len(R) for row in B),'SHAPE'
  pd(R);g=mm(mm(B,R),transpose(B));p=defaultdict(lambda:[[F(0),F(0)],[F(0),F(0)]])
  for i,u in enumerate(words):
   for j,v in enumerate(words):
    a=mm(P[i],transpose(P[j]));w=red(dag(u)+v)
    for r in range(2):
     for c in range(2):p[w][r][c]+=g[i][j]*a[r][c]
  poly[key]=p
  print('PASS',key,'exact reduced-Gram positivity and polynomial reconstruction',flush=True)
 x=[[F(0),F(1,3)],[F(1,3),F(0)]];z=[[F(1,3),F(0)],[F(0),F(-1,3)]];xz=[[F(0),F(-1,3)],[F(1,3),F(0)]]
 H={(1,):x,(-1,):x,(2,):z,(-2,):z,(-1,2):xz,(-2,1):[[-v for v in row] for row in xz]}
 zero=[[F(0),F(0)],[F(0),F(0)]]
 for w in set(poly['G'])|set(poly['J'])|set(H):
  g=poly['G'][w];j=poly['J'][w]
  assert [[g[r][c]-j[r][c] for c in range(2)] for r in range(2)]==H.get(w,zero),'SOS_IDENTITY'
  assert g[0][0]+g[1][1]==(1 if not w else 0),'PARTIAL_TRACE'
 print('PASS Q-S=H and Tr_A Q=I, exactly in the free-unitary algebra',flush=True)
 return True
if __name__=='__main__':
 verify(sys.argv[1] if len(sys.argv)>1 else Path(__file__).resolve().parent.parent/'certificates'/'sharp_recovery.json')
 print('CERTIFICATE VALID: R >= (3p-1)/2 >= (3f-1)/2 under the stated protocol assumptions.')
