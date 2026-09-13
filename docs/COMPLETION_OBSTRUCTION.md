# A fixed off-diagonal scaling cannot supply the sharp completion

This is a narrow algebraic obstruction, not a historical-novelty claim or an
exclusion of general operator-completion reductions. The exact example identities
are checked in `research/check_forward_reductions.py`; the uniqueness and
positivity steps below are prose mathematics.

## C. New calculation: no fixed scalar off-diagonal prescription works

Suppose one tries to obtain the sharp completion by requiring

```math
Q_{01}=kH_{01}
```

with one universal constant k, leaving the diagonals otherwise free. Two allowed pairs of unitaries contradict this prescription.

### Instance 1 fixes k

At U=X,V=Z on a qubit,

```math
H=\frac83\Phi-\frac23I_4,\qquad H_{01}=\frac43|0\rangle\langle1|.
```

If Q>=0, Tr_A Q=I, and Q>=H, then Tr Q=2 and <Phi|Q|Phi>>=2. A positive matrix of trace2 cannot have expectation greater than2; attaining2 exhausts its trace on that vector. Hence necessarily

```math
Q=2\Phi,\qquad Q_{01}=|0\rangle\langle1|,
```

forcing k=3/4. This also excludes the two initially tempting choices k=1 and k=1/2. The valid slack is Q-H=(2/3)(I-Phi), so this instance is feasible without the wrong prescription.

### Instance 2 excludes k=3/4

Take U=I_2,V=iI_2. Direct expansion gives

```math
H_{01}=\frac{2(1-i)}3I_2.
```

The forced k=3/4 prescription would give Q_01=(1-i)I_2/2. Compress Q to any common unit target vector. Its 2-by-2 scalar matrix must be positive, with diagonal entries summing to1, so its off-diagonal magnitude cannot exceed1/2. But |(1-i)/2|=1/sqrt(2)>1/2. Contradiction.

This second instance itself has a valid completion: H²=(8/9)I_4, so Q=(I_4+H)/2 is positive, has the required partial trace, and Q-H=(I_4-H)/2 is positive. The failure belongs to the fixed-scaling ansatz, not to the sharp theorem.

Therefore **no universal constant k can make the direct fixed-off-diagonal Ando prescription solve this problem**. `research/check_forward_reductions.py` checks both examples in exact arithmetic. The uniqueness/positivity arguments are elementary prose proofs alongside those identities.

## D. What remains of the operator-completion question

A nonlinear prescription for Q_01, a larger block embedding, or a coupled completion might still derive the sharp bound from established results. This pass does not exclude those routes.

A concrete proposed reduction must provide: its operator-algebra domain; the exact matrices to which the cited theorem is applied; verification of that theorem's positivity/numerical-radius hypotheses; the mapping back to both Q>=0 and Q-H>=0 with exact partial trace; and justification at the non-strict boundary. Merely citing complete positivity or WEP leaves the main inequality unproved.

The distinction between a finite-dimensional completion for each pair and one universal polynomial or finite-query implementation also matters. Even a successful abstract completion theorem may answer existence without answering the restricted-access implementation question.
