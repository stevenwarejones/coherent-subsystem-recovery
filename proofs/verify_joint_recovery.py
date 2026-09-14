"""Exact rational SOS checker. Standard library; no discovery-code imports."""
from fractions import Fraction as F
from pathlib import Path
from collections import defaultdict
import json,copy
if not __debug__:raise SystemExit('Run without -O')
WORDS=[(),(2,),(1,),(1,2),(2,1),(2,2),(2,1,2)]
def need(v,code):
 if not v:raise ValueError(code)
def red(w):
 out=[]
 for x in w:
  if out and out[-1]==-x:out.pop()
  else:out.append(x)
 return tuple(out)
def dag(w):return tuple(-x for x in reversed(w))
def add_outer(out,a,b,left,right,scale):
 for w,c in left.items():
  for v,d in right.items():out[a,b,red(dag(w)+v)]+=scale*c*d
def clean(D):return {k:v for k,v in D.items() if v}
def target(alpha,beta):
 D=defaultdict(F);P={():F(1,2),(2,):F(1,2)};N={():F(1,2),(2,):F(-1,2)}
 UP={(1,)+w:c for w,c in P.items()};UN={(1,)+w:c for w,c in N.items()}
 for B in [[P,UN],[UP,N]]:
  for a in range(2):
   for b in range(2):add_outer(D,a,b,B[a],B[b],F(1,4))
 L=[[{():F(1),(2,):F(1)},{(1,):F(1)}],[{(1,):F(1)},{():F(1),(2,):F(-1)}]]
 for t in range(2):
  for a in range(2):
   for b in range(2):add_outer(D,a,b,L[t][a],L[t][b],-alpha/9)
 # q: Kraus PUP and NUN, expanded independently of the moment search.
 for a in range(2):
  for B in [{(1,):F(1,4),(1,2):F(1,4),(2,1):F(1,4),(2,1,2):F(1,4)},
            {(1,):F(1,4),(1,2):F(-1,4),(2,1):F(-1,4),(2,1,2):F(1,4)}]:add_outer(D,a,a,B,B,beta)
  D[a,a,()]+=alpha-1
 return clean(D)
def psd(G):
 A=copy.deepcopy(G);n=len(A);rank=0
 need(all(A[i][j]==A[j][i]for i in range(n)for j in range(n)),'GRAM_SYMMETRY')
 for k in range(n):
  p=A[k][k];need(p>=0,'GRAM_PSD')
  if p==0:
   need(all(A[k][j]==0 for j in range(k+1,n)),'GRAM_PSD');continue
  rank+=1
  for i in range(k+1,n):
   for j in range(i,n):A[j][i]=A[i][j]=A[i][j]-A[i][k]*A[j][k]/p
 return rank

def verify(d):
 need(F(d['alpha'])==F(5,4),'ALPHA');need(F(d['beta'])==3,'BETA')
 need([tuple(w)for w in d['words']]==WORDS,'WORDS')
 G=[[F(x)for x in row]for row in d['gram']];n=14
 need(len(G)==n and all(len(row)==n for row in G),'GRAM_SHAPE')
 rank=psd(G);out=defaultdict(F)
 for i in range(n):
  for j in range(n):out[i//7,j//7,red(dag(WORDS[i%7])+WORDS[j%7])]+=G[i][j]
 need(clean(out)==target(F(5,4),F(3)),'POLYNOMIAL_IDENTITY')
 return rank,len(clean(out))

def main():
 d=json.loads((Path(__file__).resolve().parents[1]/'certificates'/'joint_certificate.json').read_text());rank,n=verify(d)
 tests=[]
 t=copy.deepcopy(d);t['alpha']='6/5';tests.append((t,'ALPHA'))
 t=copy.deepcopy(d);t['beta']='2';tests.append((t,'BETA'))
 t=copy.deepcopy(d);t['words'][1]=[1];tests.append((t,'WORDS'))
 t=copy.deepcopy(d);t['gram'][0][0]='-1';tests.append((t,'GRAM_PSD'))
 t=copy.deepcopy(d);t['gram'][0][0]=str(F(t['gram'][0][0])+1);tests.append((t,'POLYNOMIAL_IDENTITY'))
 for t,expected in tests:
  try:verify(t)
  except ValueError as e:need(str(e)==expected,'WRONG_MUTATION_REASON')
  else:raise ValueError('MUTATION_ACCEPTED')
 print('PASS exact 14x14 Gram PSD, rank',rank)
 print('PASS free-unitary identity on',n,'matrix-word coefficients')
 print('CERTIFIED: E - (5/4) K + 3 Q + I/4 >= 0')
 print('CONCLUSION: r >= 1 - (5/4)(1-f) - 3q, using the base bridge f <= p')
 print('PASS five mutations reject for their intended reasons')
if __name__=='__main__':main()
