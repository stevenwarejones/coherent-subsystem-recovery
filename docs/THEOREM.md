# Sharp recovery from the three-branch coherent test

13 September 2026. Exact certificate and analytic protocol proof. Produced by one AI system and independently reconstructed from scratch by a second (see `VERIFICATION.md`); **not reviewed by a human domain expert**, and historical novelty remains unresolved. Within this repository the theorem is stated for the full score range in `FULL_CURVE.md`; this document proves the nontrivial high-score segment and its attainability. File names below refer to their repository paths: the certificate is `certificates/sharp_recovery.json`, and the checkers are `proofs/verify_sharp_recovery.py`, `proofs/verify_alternate.py`, `tests/test_certificate_mutations.py`, and `research/check_attainment.py`.

## Main result

Under the protocol assumptions below, for arbitrary finite-dimensional unitaries U,V on the addressed subsystem M and arbitrary finite-dimensional bypass B,

\[
\boxed{R\ge (3f-1)/2.}
\]

For every f in [2/3,1], an allowed construction attains equality. Consequently the minimum recoverability compatible with an exact score f in this interval is exactly (3f-1)/2. The same minimum holds if the constraint is score at least f. The construction uses addressed dimension at most six, together with a finite bypass.

The earlier universal certificate had R >= (500/1041)(3f-1). The new certificate has the exact normalization t=1 rather than t=1041/1000. No numerical tolerance or limiting argument enters the new bound. This is not only the Hermitian or commuting subclass result.

For f=1-epsilon, the sharp guarantee is R >= 1-(3/2)epsilon, reaching perfect recovery at perfect test score. (An earlier draft compared this against a Renes-2016-derived 1-6epsilon guarantee; the current, stronger literature comparison is Berta-Coles-Wehner 2014 + Renes 2017, giving 1-(9/4)epsilon+(27/16)epsilon^2 -- see `PRIOR_ART.md`, which is the canonical comparison. Both are established prior art; neither establishes historical priority for this protocol-specific sharp bound.)

## 1. Protocol and normalization

A trusted input qubit is initially maximally entangled with protected reference A. An encoding isometry maps the input into M tensor B. Let rho_AM denote the reduced state before the intervention. Thus rho_A=I/2. B is an unrestricted finite-dimensional quantum bypass; its size is not assumed bounded by M.

A protected three-dimensional controller starts in the uniform coherent superposition. Its branches apply I,U,V to M alone. A single subsequent trace-preserving decoder may act on MB, but not on A, the controller, or an extra controller-branch label. There is no postselection. Coherent implementation of the controlled unitaries, including their phases, is assumed; specifying their isolated reduced channels does not suffice.

Let P_0=I,P_1=X,P_2=Z, and

\[
|\Omega\rangle={1\over\sqrt3}\sum_{a=0}^2|a\rangle_C
 (I_A\otimes P_a)|\Phi^+\rangle_{AO}.
\]

The score f is the probability of projection onto this full CAO pure target, with decoder garbage discarded. Define

\[
R=\max_{\mathcal D:M\to\mathbb C^2\;\mathrm{CPTP}}
\langle\Phi^+|(\mathrm{id}_A\otimes\mathcal D)(\rho_{AM})|\Phi^+\rangle.
\]

R is squared entanglement fidelity for the maximally mixed qubit input. It is not a worst-case fidelity or a diamond-norm guarantee. The corresponding average pure-input qubit fidelity is (2R+1)/3. The theorem concerns recoverability from M at the tested cut, not continuous retention, semantic records, consciousness, or a unique internal computation.

## 2. From the observed score to the operator polynomial

Define the isometry T from a qubit into CA by blocks T_a=P_a/sqrt(3). Since these three Paulis are real symmetric, T†T=I, and the ideal CA marginal is TT†/2. Put P=TT†.

Writing U_0=I,U_1=U,U_2=V, the actual CA marginal has blocks

\[
\tau_{ab}={1\over3}\operatorname{Tr}_M[(I\otimes U_a)\rho(I\otimes U_b^\dagger)].
\]

The common decoder leaves this marginal unchanged. The ideal target projector is supported on P tensor I_O, hence f <= p:=Tr(P tau).

For

\[
L=I+X\otimes U+Z\otimes V,\qquad K=L^\dagger L/9,
\]

direct block multiplication gives

\[
T^\dagger\tau T={1\over9}\operatorname{Tr}_M(L\rho L^\dagger),\qquad
p=\operatorname{Tr}(K\rho).
\]

Therefore H:=3K-I equals

\[
H={1\over3}\left[X\otimes(U+U^\dagger)+Z\otimes(V+V^\dagger)
+XZ\otimes(U^\dagger V-V^\dagger U)\right].
\]

Only unitarity was used, not commutation or Hermiticity.

## 3. The exact free-unitary certificate

The supplied certificate `certificates/sharp_recovery.json` contains:

- all 53 reduced words of length at most three in U,U†,V,V†;
- rational basis matrices B_G of size 53 by 46 and B_J of size 53 by 48;
- symmetric rational positive-definite matrices R_G and R_J of sizes 46 and 48;
- the claimed constant t=1.

For a word w_i, let P_i be the real orthogonal 2-by-2 matrix obtained by replacing either sign of the U letter by X and either sign of the V letter by Z, retaining word order. This P_i is an auxiliary numerical matrix; it does not impose a Pauli relation on the physical unitaries.

Set g=B_G R_G B_G^T and j=B_J R_J B_J^T. Form 106-by-106 block Gram matrices

\[
\mathsf G_{ij}=g_{ij}P_iP_j^T,\qquad
\mathsf J_{ij}=j_{ij}P_iP_j^T.
\]

Both are positive semidefinite: for example mathsf G is the congruence of g tensor I_2 by diag(P_i). The reduced matrices are positive definite; the full matrices need only be semidefinite.

For every physical pair U,V, define

\[
Q=\sum_{ijab}(\mathsf G_{ij})_{ab}|a\rangle\langle b|\otimes w_i^\dagger w_j,
\quad
S=\sum_{ijab}(\mathsf J_{ij})_{ab}|a\rangle\langle b|\otimes w_i^\dagger w_j.
\]

Positive Gram matrices imply Q,S>=0 by substituting the block column with entries I_A tensor w_i. The exact checker verifies, coefficient by coefficient after adjacent inverse cancellation only,

\[
\boxed{Q-S=H,\qquad \operatorname{Tr}_A Q=I_M.}
\]

Thus Q>=H. These identities hold for every pair of unitaries in every finite dimension. The checker does not infer their truth from floating-point eigenvalues, sampled matrices, or equality-point constraints.

## 4. The recovery channel and the lower bound

Define a map E from qubit operators to operators on M by E(|a><b|)=Q_ab. Its standard Choi matrix is Q, so E is completely positive. The trace identity gives E(I)=I_M. Its adjoint D is therefore a completely positive trace-preserving map from M to a qubit.

The Bell acceptance effect of this recovery is

\[
(\mathrm{id}\otimes E)(|\Phi^+\rangle\langle\Phi^+|)=Q/2.
\]

Consequently

\[
R\ge {1\over2}\operatorname{Tr}(Q\rho)
\ge {1\over2}\operatorname{Tr}(H\rho)
={3p-1\over2}\ge {3f-1\over2}.
\]

This supplies an actual algebraically specified recovery channel on M alone. It does not prove an efficient implementation using only forward black-box calls to U,V. Adjoints appear in the polynomial certificate. The operator inequality itself does not need rho_A=I/2; that marginal condition is used for the physical encoding interpretation.

## 5. Exact attainability

### Optimizing the allowed common decoder

Let C=T†tau T. For a purification on AMB and unrestricted common decoder on MB, Uhlmann's theorem gives

\[
f_{\rm opt}=F_{\rm root}(\tau,P/2)^2
={1\over2}(\operatorname{Tr}\sqrt C)^2.
\]

The decoder can include a discarded auxiliary output. It has no access to the controller. In particular, if C=(p/2)I, then f_opt=p.

### Classical endpoint at f=2/3, R=1/2

Begin with the pure reference state

\[
\rho_0={1\over2}\left[I+{X-Y+Z\over\sqrt3}\right]
\]

and scalar interventions u=(sqrt(3)+i)/2, v=(sqrt(3)-i)/2. To restore the required reference marginal, introduce a four-valued flag as the addressed subsystem. With g ranging over I,X,Y,Z, set

\[
\rho_{AM}={1\over4}\bigoplus_g g\rho_0g^\dagger,
\quad U=\bigoplus_g s_X(g)u,
\quad V=\bigoplus_g s_Z(g)v,
\]

where gXg†=s_X(g)X and gZg†=s_Z(g)Z. Direct calculation gives rho_A=I/2 and C=I/3, hence f_opt=2/3. The flag is part of the pre-intervention addressed subsystem, not an extra controller-branch label.

The state is separable, so no local recovery exceeds Bell overlap 1/2. Reading the flag and preparing the conjugate of its associated pure reference state attains 1/2. Thus R=1/2 exactly. The maximally mixed A marginal guarantees a valid encoding isometry into M and a sufficiently large purifying bypass.

### Perfect endpoint and mixtures

The Bell state on AM with U=X,V=Z has f=R=1 and C=I/2. Take an addressed direct sum of this two-dimensional construction and the preceding four-dimensional construction. Put probability lambda on the perfect block and 1-lambda on the classical block.

Recovery on a flagged direct sum is exactly the weighted sum of the block optima. The compressed C matrices are scalar, so the optimized coherent score is also the weighted sum. Therefore, for 0<=lambda<=1,

\[
f={2+\lambda\over3},\qquad R={1+\lambda\over2}={3f-1\over2}.
\]

This verifies attainability throughout [2/3,1], with addressed dimension at most six. `research/check_attainment.py` verifies the reference marginal, both compressed C matrices, and the mixture identity using exact symbolic arithmetic. The recovery and Uhlmann interpretations are the prose arguments above, not outputs of that script.

## 6. What was verified

- `proofs/verify_sharp_recovery.py`: Python standard library; exact positive pivots for the reduced Gram matrices; exact reconstruction of all free-word coefficients; exact partial-trace normalization; fixed claimed constant. Passes.
- `proofs/verify_alternate.py`: python-flint; a different full-block reconstruction and exact leading principal determinants (Sylvester criterion). Passes. Written by the same AI, not independent external review.
- `tests/test_certificate_mutations.py`: unmodified certificate passes; four corruptions fail for their stated reasons; rejected inputs print no final success claim; -O is refused. One mutation adds the same positive constant summand to Q and S, preserving positivity and their difference while violating only normalization.
- `research/check_attainment.py`: exact symbolic endpoint and mixture checks. Passes. An independent from-scratch reconstruction is in `tests/verify_attainment_independent.py`.

The discovery calculation was numerical. It reached about 0.99999978 at word length three, which alone was not a proof. Equality constraints and Pauli sign symmetry reduced the unknown Gram matrices. A reduced numerical interior point had positivity margin about 0.00012277. Rational rounding followed by exact linear correction produced the shipped certificate. The two exact verification routes, not the numerical optimum, justify the result.

## 7. Literature status and what this does not establish

This pass searched recovery, numerical-radius matrix completions, and unitary operator systems, including:

- Renes 2016, https://arxiv.org/abs/1605.01420 : established complementary-observable recovery principle; the full theorem was inspected in the preceding pass.
- Farenick, Kavruk, Paulsen, https://arxiv.org/abs/1107.0418 : primary abstract inspected this pass. Matrix-completion theorems related to Ando's theorem are relevant established mathematics. No theorem from this source is used to justify the certificate.
- Farenick, Ojo, Plosker, https://arxiv.org/abs/2101.00129 : primary abstract inspected this pass. Its completely positive universality statements assume Weyl relations, which cannot simply be assumed for arbitrary U,V.
- Noncommutative Gram/SOS certificates and Choi positivity are standard methods; neither is claimed as new.

Additional exact-phrase searches did not identify this protocol-specific curve, but they are not a comprehensive novelty audit. The earlier recovery, quantum-memory verification, and supermap-self-testing comparators still require a precise theorem-level comparison.

What is now supplied is an exact proof and matching constructions for the stated protocol. What remains unresolved is historical priority, external mathematical review, efficient physical recovery implementation, robustness to violations of trusted intervention support/coherence, and experimental feasibility. Infinite-dimensional domain extensions are not claimed.

## 8. Instructions for the next reviewer

1. Run the standard-library verifier before installing any solver.
2. Rebuild the full Gram matrices and polynomial identity independently if possible; do not treat the reported numerical optimization value as evidence of proof.
3. Check the congruence argument converting the stored reduced matrices into positive operator polynomials.
4. Reconstruct the controller support-projector inequality f<=p and the Choi factor of two.
5. Check the common-decoder Uhlmann argument, the classical flagged endpoint, and direct-sum recovery additivity. These establish sharpness rather than only validity.
6. Run the mutations and verify their intended failure reasons.
7. Audit novelty separately. Do not infer historical priority from successful arithmetic checks.
8. If these survive, replace the old 'sharp curve conjectural' statements in the new-repo handoff with this theorem and its precise scope. Keep earlier artifacts as dated research history rather than competing current claims. A small separate research repo is appropriate; external publicity should use the reviewed theorem, not this session's confidence.
