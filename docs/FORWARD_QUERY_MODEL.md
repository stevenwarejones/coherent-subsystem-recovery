# 1. General forward-only question: an exact reduction, not a no-go theorem

Independent AI derivation with separate exact and numerical checks (see [VERIFICATION.md](VERIFICATION.md)). Human expert review and
historical novelty remain outstanding.


13 September 2026. New research derivation; independent review required.

## Fixed model

A single finite query-history circuit is used in every finite target dimension. Its only actions on M are controlled forward U,V. It may have arbitrary fixed finite controller workspace, controller-only channels, randomness and deferred measurements; no postselection, target operations, or reference/bypass access. Coefficients are independent of the unknown oracles, input state, and observed score. Query order is fixed. The target is r >= (3f-1)/2 for every physical protocol instance with f >= 2/3 (or, equivalently for existence below, any fixed high-score interval [f_0,1] with f_0<1).

The initial proof does not cover circuits chosen separately for each score. A compactness extension below covers dimension-dependent coefficient families when the query count has a uniform finite bound. Classical random mixtures of finitely many fixed orders can be padded into one fixed supersequence using controlled skipping, so they fit after a finite query-budget increase. Indefinite causal order is excluded.

## Result: the operator search is complete for existence after symmetrization

Write L=I+X⊗U+Z⊗V, K=L†L/9, H=3K-I. Every recovery D has Bell acceptance effect E_D=(id⊗D†)(Phi).

**Proposition.** A universal sharp recovery exists in the fixed forward-only model above if and only if one exists in that model with

```math
E_D(U,V)\succeq H(U,V)/2
```

for every pair of finite-dimensional unitaries, with no input-state restriction. The converse follows immediately from f <= Tr(K rho). Necessity uses an implementable Pauli symmetrization and direct-sum flags; proof below.

This removes one uncertainty from the search: after this symmetrization, imposing the operator inequality is not excluding all possible sharp forward-only solutions merely because it also tests states with a nonmaximally mixed reference. It does **not** prove feasibility or infeasibility of the resulting operator problem.

## Proof of necessity

Let g range over I,X,Y,Z and let signs s_X(g),s_Z(g) satisfy gXg†=s_X(g)X and gZg†=s_Z(g)Z. Given a hypothetical universal algorithm D, define

```math
\widetilde D_{U,V}(\sigma)
 =\frac14\sum_g g^T\,D_{s_X(g)U,s_Z(g)V}(\sigma)\,\overline g.
```

This is allowed: multiplying an oracle by -1 is implemented by a Z phase on its control, and output conjugation acts only on the controller output. It needs classical randomness but no new target operation. For fixed query order the query count is unchanged.

Fix arbitrary normalized rho_AM, without imposing rho_A=I/2. Let p=Tr(K rho). Add a four-valued classical flag F **inside M** and set

```math
\rho'=\frac14\bigoplus_g (g\otimes I)\rho(g^\dagger\otimes I),\qquad
U'=\bigoplus_g s_X(g)U,\quad V'=\bigoplus_g s_Z(g)V.
```

Then rho'_A=I/2. The fixed circuit D acts blockwise on these oracles; discarding the flag averages its outputs. The Bell identity (g⊗I)|Phi>=(I⊗g^T)|Phi> shows that its recovery score on this flagged instance is exactly Tr(E_Dtilde rho).

For the original coherent test, define

```math
C=\frac19\operatorname{Tr}_M(L\rho L^\dagger).
```

The flag construction gives C'=(p/2)I. The optimized score with an unrestricted common decoder on the purified encoding output is f_opt=(Tr sqrt(C'))²/2=p. The flagged state has a valid encoding purification because rho'_A=I/2. Thus, whenever p>=2/3, the assumed physical guarantee directly gives Tr(E_Dtilde rho)>=(3p-1)/2.

For p<2/3, append an additional direct-sum block consisting of the ideal Bell/Pauli instance. Mix it with weight 1-theta, retaining the first flagged instance with weight theta>0 chosen so p_theta=theta p+1-theta>=2/3. Both block C matrices are scalar, so f_opt=p_theta; this is an actual score attained by a common decoder, not a replacement of f by an unattainable support score. The ideal instance must have recovery 1 by the original universal hypothesis at f=1. Fixed-circuit direct-sum consistency therefore gives

```math
\theta\operatorname{Tr}(E_{\widetilde D}\rho)+(1-\theta)
\geq\frac{3(\theta p+1-\theta)-1}{2}.
```

Cancel the ideal contribution and divide by theta to obtain the same inequality for the original arbitrary rho. Since this holds for every density matrix, E_Dtilde-H/2 is positive semidefinite.

The proof uses the previously derived optimized-score identity and fixed-circuit direct-sum consistency. Neither may be silently omitted.

## Stronger consequence: restricting to a neighborhood of f=1 does not evade the operator target

The same ideal-block mixture can raise any p to any fixed threshold f_0<1 while preserving a strict deficit by multiplication with theta. Therefore a fixed universal algorithm satisfying the exact affine sharp bound merely on [f_0,1] already has a symmetrized version satisfying the full operator inequality. The algorithm must include perfect recovery at f=1, remain the same algorithm throughout that interval, and be consistent on direct sums.

This is not a robustness claim: the diluted violation can be arbitrarily small. It is an exact-existence statement. It also does not cover score-selected algorithms, since mixing can select a different circuit there.

## Complete finite-query representation

For n fixed controlled calls, histories h∈{0,1}^n give positive words w_h. Purifying the controller processing and final decoder yields coefficients c_(a,e,h) and a PSD matrix T_(ah,bk)=sum_e c_(a,e,h) conjugate(c_(b,e,k)).

Its output trace G_n has zero blocks between histories with different last bits. Summing its two equal-last-bit blocks gives G_(n-1); repeat this condition at every earlier bit, ending in scalar1. Conversely, factor these Grams and use equal-Gram isometries between parent vectors and orthogonal child families to reconstruct a finite controller circuit. This is the standard Gram/comb method specialized to this target interface.

The Bell-effect blocks are

```math
(E_T)_{ab}=\frac12\sum_{h,k}T_{bh,ak}w_k^\dagger w_h.
```

Therefore the remaining problem is precise: find a finite history length and a causal PSD T such that E_T-H/2 is universally positive, or show that no finite history length admits such T. The existing forbidden-word argument rules out one shipped Q, not every admissible E_T.

For arbitrary chosen supports, a finite-degree SOS search is only a sufficient positivity test. Its failure does not establish this general no-go. The suffix-closed history support has a stronger, proposed complete formulation described below; that refinement needs its finite-realization and duality argument. A solver failure or inaccurate optimum is likewise not evidence of nonexistence.

## Compactness extension: dimension-dependent controllers do not evade bounded-query existence

Suppose that for every finite target dimension D there is a controller-only forward algorithm A_D with at most n calls, satisfying the same universal high-score guarantee. Its controller size and fixed gates may depend on D, but not on the unknown U,V, state or observed score. Each circuit has a fixed definite order of U,V calls.

There are finitely many call patterns of length at most n. Therefore an unbounded sequence of target dimensions D_j uses the same pattern. Represent those circuits by history matrices T_j of the same finite size. The causal PSD constraints and Tr T_j=1 define a closed bounded set. Passing to a subsequence gives T_j converging to a causal PSD T_infinity. Controller workspaces need not have a uniform dimension bound: the history matrix already compresses them.

Fix any test dimension d and any U,V,rho on that dimension. For D_j>=4d+2, embed the four-flag construction and a two-dimensional ideal block inside dimension D_j; leave the remaining subspace unoccupied and choose arbitrary unitary blocks there. Apply A_Dj to that physical instance. Although A_Dj is not the circuit used by the family in dimension d, its history coefficients still act blockwise on the chosen d-dimensional oracle blocks. The preceding flag-and-mixture proof thus implies the required operator inequality for the Pauli-symmetrized T_j on that fixed d, U,V, rho.

The Pauli averaging is a continuous linear operation on T. Take the limit. Because the evaluation is linear in its finitely many coefficients, the inequality survives for every fixed d,U,V,rho. The symmetrized T_infinity is causal and PSD, so the Gram realization produces one fixed finite controller circuit of the same query pattern satisfying the operator inequality in every dimension.

**Consequently:** for existence in the ideal arbitrary-controller-gate model, allowing dimension-dependent controller designs with a dimension-uniform finite query budget does not enlarge the class of universally sharp recoveries. A fixed universal history-matrix solution exists if such a bounded-query family exists. This strengthens the search formulation, not the answer to its feasibility question.

Limitations: this is a compactness existence argument, not an effective algorithm for computing T_infinity or its gate coefficients. It does not cover query budgets increasing without bound with dimension, score-selected circuits, extra target operations, or indefinite order. The history realization and limiting argument are prose mathematics and require independent review.

## What should change in the research plan

Stop treating each small circuit variation as a separate foundational question. Use the common T formulation, implement the Pauli symmetry, and demand a universal certificate for a candidate. The present reduction justifies that formulation as an existence search for the stated model. The unrestricted finite-forward-query question is still open.

Literature context: Gram-matrix quantum-query formulations and quantum combs are established techniques, e.g. [Lee et al., arXiv:1011.3020](https://arxiv.org/abs/1011.3020) and [Chiribella et al., arXiv:0904.4483](https://arxiv.org/abs/0904.4483). The protocol-specific step here is the flag-plus-ideal-mixture equivalence. No historical-priority claim is made.

## Support-specific completeness refinement

The following argument from the subsequent classification handoff specializes the
preceding general SOS warning. It is prose mathematics for review, not a theorem
formalized by the numeric or certificate checks. The robust upper bound only
needs necessity of controller constraints and does not depend on their converse.

## 2. Exact realization of the relevant truncated moments

Fix an order of n queries. Each binary history h produces a positive word w_h: take the subsequence of calls that fired, in reverse chronological order for matrix multiplication. Let W be the set of distinct such words, including the empty word. W is closed under deleting the leftmost letter. It has at most 2^n elements. Assume both single-letter words occur; otherwise add them for the general moment statement, or treat the one-generator recovery separately.

For a joint vector psi=sum_a |a> psi_a, form the matrix

\[
\Gamma_{(a,w),(b,v)}=\langle w\psi_a,v\psi_b\rangle,
\qquad a,b\in\{0,1\},\quad w,v\in W.
\]

It satisfies:

1. Gamma is positive semidefinite.
2. Gamma_(0,e),(0,e)+Gamma_(1,e),(1,e)=1.
3. Whenever xw,xv belong to W, for x=U or V,
   Gamma_(a,xw),(b,xv)=Gamma_(a,w),(b,v).

**Realization lemma.** These constraints are sufficient, not just necessary. Every such Gamma is realized by two unitaries and a normalized joint vector in target dimension at most rank(Gamma)<=2|W|.

**Proof.** Factor Gamma as the Gram matrix of vectors v_(a,w) in their finite-dimensional span S. For each x, define the partial map v_(a,w) -> v_(a,xw) on the span of all vectors for which xw lies in W. Condition 3 preserves every inner product between the generators of this domain. It therefore makes the map well-defined, isometric, and onto its image, even when the original vectors are linearly dependent. Domain and image have equal dimension. Their orthogonal complements in the same finite S have equal dimension, so the partial map extends to a unitary U_x on S. U and V extend independently; no commutation condition is imposed. Suffix closure and induction now give v_(a,w)=w(U,V)v_(a,e). Condition 2 normalizes sum_a |a>v_(a,e). This proves the assertion, including singular Gram matrices.

This is a finite Gram/unitary-extension argument in an established mathematical tradition, not a new general positivity principle. Closely related constructions appear in the literature discussed in Section 9.

## 3. A complete SDP for each fixed pattern

After purifying the controller, its history amplitudes have coefficients c_(a,e,h), where a is the output qubit and e the discarded controller environment. Define

\[
T_{ah,bk}=\sum_e c_{a,e,h}\overline{c_{b,e,k}}\succeq0.
\]

Let G_n be the output trace of T. Histories differing in their final query bit are orthogonal before the last controller operation, so the cross-bit block of G_n vanishes. Sum the two equal-bit diagonal blocks to obtain G_(n-1), impose the same rule, and continue to the final scalar 1. These are the causal history constraints. Conversely, successive equal-Gram isometries reconstruct a controller circuit from these constraints; one extends each controller isometry to a unitary after adding workspace. This is the usual Gram/comb realization specialized to controlled calls on an otherwise inaccessible target.

The Bell-effect blocks are

\[
(E_T)_{ab}=\frac12\sum_{h,k}T_{bh,ak}w_k^\dagger w_h.
\]

Every expectation of 2E_T-H depends only on Gamma from Section 2. Consequently, minimizing that expectation over the finite moment spectrahedron is exactly the minimization over all finite-dimensional unitaries and input vectors. It is not an NPA relaxation requiring a higher level.

The moment problem has a strictly positive feasible matrix Gamma=I_(2|W|)/2. Its feasible set is compact: repeated prefix cancellation fixes every diagonal to one of the two empty-word diagonals, and positivity bounds all off-diagonal entries. Thus finite-dimensional SDP strong duality applies. Dualizing the inner moment minimization gives a single joint SDP for T and a positive Gram representation of 2E_T-H. The only zero identities used are the unitary cancellations already imposed on Gamma. Equivalently, a positive-word SOS supported on W suffices for this particular hereditary polynomial problem; there is no need to guess a larger inverse-word support.

**Correction to the preceding research note:** its generic warning that failure of a finite-support SOS need not prove impossibility remains true in general. Under this suffix-closed history support and the realization lemma, the corresponding fixed-pattern search is complete. A floating-point negative optimum still is not a proof: exact dual certificates are required.

The conclusion is an exact finite formulation, not a claim that an ordinary numerical SDP solver decides arbitrary boundary feasibility without numerical difficulties.
