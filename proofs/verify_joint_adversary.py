"""Exact moment counterexample for the fixed balanced recovery, not all controllers."""
from fractions import Fraction as F
from pathlib import Path
import json,copy
from verify_joint_recovery import WORDS,psd,need

def verify(d):
 need([tuple(w)for w in d['words']]==WORDS,'WORDS')
 H=[[F(x)for x in row]for row in d['moment']];m=len(WORDS);ix={w:i for i,w in enumerate(WORDS)}
 need(len(H)==2*m and all(len(r)==2*m for r in H),'MOMENT_SHAPE')
 need(psd(H)==2*m,'MOMENT_PD');need(H[0][0]+H[m][m]==1,'NORMALIZATION')
 for x in [1,2]:
  dom=[w for w in WORDS if (x,)+w in ix]
  for a in range(2):
   for b in range(2):
    for w in dom:
     for v in dom:need(H[a*m+ix[(x,)+w]][b*m+ix[(x,)+v]]==H[a*m+ix[w]][b*m+ix[v]],'UNITARY_MOMENTS')
 def inner(a,b,D,E):return sum(c*e*H[a*m+ix[w]][b*m+ix[v]]for w,c in D.items()for v,e in E.items())
 P={():F(1,2),(2,):F(1,2)};N={():F(1,2),(2,):F(-1,2)}
 UP={(1,)+w:c for w,c in P.items()};UN={(1,)+w:c for w,c in N.items()}
 r=sum(inner(a,b,B[a],B[b])/4 for B in [[P,UN],[UP,N]]for a in range(2)for b in range(2))
 L=[[{():F(1),(2,):F(1)},{(1,):F(1)}],[{(1,):F(1)},{():F(1),(2,):F(-1)}]]
 p=sum(inner(a,b,L[t][a],L[t][b])/9 for t in range(2)for a in range(2)for b in range(2))
 q=F(0)
 for a in range(2):
  for B in [{(1,):F(1,4),(1,2):F(1,4),(2,1):F(1,4),(2,1,2):F(1,4)},
            {(1,):F(1,4),(1,2):F(-1,4),(2,1):F(-1,4),(2,1,2):F(1,4)}]:q+=inner(a,a,B,B)
 vals=dict(p=p,r=r,q=q,gap=(3*p-1)/2-r,ratio=q/(1-p))
 for key,val in vals.items():need(val==F(d[key]),'STORED_VALUE')
 need(F(2,3)<p<1 and vals['gap']>0,'RECOVERY_FAILURE')
 need(0<vals['ratio']<F(1,8),'RATIO')
 need(r>=(5*p-1)/4-3*q,'LOWER_BOUND_CONSISTENCY')
 return vals

def main():
 d=json.loads((Path(__file__).resolve().parents[1]/'certificates'/'joint_adversary.json').read_text());vals=verify(d)
 for key in ['q','gap']:
  t=copy.deepcopy(d);t[key]=str(F(t[key])+1)
  try:verify(t)
  except ValueError as e:need(str(e)=='STORED_VALUE','WRONG_MUTATION_REASON')
  else:raise ValueError('MUTATION_ACCEPTED')
 # Structural mutations exercising the moment guards themselves (review-recommended):
 # scaling preserves positivity/consistency but breaks normalization; a positive
 # diagonal bump preserves normalization but breaks unitary consistency; a negative
 # diagonal breaks Gram positivity. Each must reject for its own reason.
 def scaled(t):
  t=copy.deepcopy(t);t['moment']=[[str(F(x)*2)for x in row]for row in t['moment']];return t
 def bumped(t):
  t=copy.deepcopy(t);t['moment'][1][1]=str(F(t['moment'][1][1])+1);return t
 def negated(t):
  t=copy.deepcopy(t);t['moment'][0][0]='-1';return t
 for mutate,reason in [(scaled,'NORMALIZATION'),(bumped,'UNITARY_MOMENTS'),(negated,'GRAM_PSD')]:
  try:verify(mutate(d))
  except ValueError as e:need(str(e)==reason,'WRONG_MUTATION_REASON')
  else:raise ValueError('MUTATION_ACCEPTED')
 print('PASS exact 14x14 positive-definite moment matrix and all unitary equalities')
 for k,v in vals.items():print(k,'=',v,'=',float(v))
 print('CONCLUSION: this balanced recovery fails although q < (1-f)/8 after the stated physical flag construction.')
 print('This is NOT an impossibility result for other recovery circuits.')
 print('PASS two stored-claim and three structural mutations rejected for their intended reasons')
if __name__=='__main__':main()
