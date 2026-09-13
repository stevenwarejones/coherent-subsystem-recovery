"""Exact, solver-free check. Run with Python 3; uses only standard library."""
if not __debug__:
    raise RuntimeError("Verification requires assertions enabled; do not use python -O.")
from fractions import Fraction as F
from collections import defaultdict
from pathlib import Path
import json,sys

def reduce_word(word):
    out=[]
    for letter in word:
        if out and out[-1]==-letter:out.pop()
        else:out.append(letter)
    return tuple(out)
def dagger(word):return tuple(-x for x in reversed(word))
def positive_definite(a):
    n=len(a)
    assert all(a[i][j]==a[j][i] for i in range(n) for j in range(n)), 'SYMMETRY'
    l=[[F(0) for _ in range(n)] for _ in range(n)];d=[]
    for i in range(n):
        pivot=a[i][i]-sum(l[i][k]**2*d[k] for k in range(i))
        assert pivot>0, 'POSITIVE_DEFINITE'
        d.append(pivot);l[i][i]=F(1)
        for j in range(i+1,n):
            l[j][i]=(a[j][i]-sum(l[j][k]*l[i][k]*d[k] for k in range(i)))/pivot
    assert all(a[i][j]==sum(l[i][k]*d[k]*l[j][k] for k in range(n)) for i in range(n) for j in range(n)), 'LDL_IDENTITY'
def polynomial(a,words):
    p=defaultdict(lambda:[[F(0),F(0)],[F(0),F(0)]])
    for i,u in enumerate(words):
        for j,v in enumerate(words):
            w=reduce_word(dagger(u)+v)
            for x in range(2):
                for y in range(2):p[w][x][y]+=a[2*i+x][2*j+y]
    return p

def verify(path):
    z=json.loads(Path(path).read_text());words=[tuple(w) for w in z['words']]
    assert len(words)==17 and len(set(words))==17, 'WORD_BASIS'
    assert all(all(x in (-2,-1,1,2) for x in w) and reduce_word(w)==w for w in words), 'WORD_BASIS'
    n=2*len(words);a={k:[[F(v) for v in row] for row in z[k]] for k in ('G','J')}
    assert all(len(m)==n and all(len(row)==n for row in m) for m in a.values()), 'SHAPE'
    t=F(z['t']);assert t==F(1041,1000), 'CLAIMED_CONSTANT'
    for m in a.values():positive_definite(m)
    g,j=(polynomial(a[k],words) for k in ('G','J'))
    x=[[F(0),F(1,3)],[F(1,3),F(0)]];zz=[[F(1,3),F(0)],[F(0),F(-1,3)]]
    xz=[[F(0),F(-1,3)],[F(1,3),F(0)]]
    h={(1,):x,(-1,):x,(2,):zz,(-2,):zz,(-1,2):xz,(-2,1):[[-v for v in row] for row in xz]}
    zero=[[F(0),F(0)],[F(0),F(0)]]
    for w in set(g)|set(j)|set(h):
        assert [[g[w][r][s]-j[w][r][s] for s in range(2)] for r in range(2)]==h.get(w,zero), 'SOS_IDENTITY'
        assert g[w][0][0]+g[w][1][1]==(t if not w else 0), 'PARTIAL_TRACE'
    return True
if __name__=='__main__':
    verify(sys.argv[1] if len(sys.argv)>1 else Path(__file__).with_name('certificate.json'))
    print('PASS: exact positive definiteness, free-unitary polynomial identity, partial trace, and fixed constant.')
    print('CONCLUSION: R >= (500/1041)(3p-1) >= (500/1041)(3f-1), under the stated protocol assumptions.')
