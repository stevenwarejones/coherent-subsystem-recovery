# A two-query recovery identity with an algebraic error parameter

Independent AI derivation and separate exact/numerical implementations:
[verification audit](FORWARD_RESULTS_AUDIT.md). Human expert review and
historical novelty remain outstanding.


## Result

Assume V is an exact reflection, meaning V=V-adjoint and V^2=I. U can be any unitary in any finite dimension. Use the balanced two-query extraction defined below, and write

\[
\delta=\|UV+VU\|.
\]

Then its achieved Bell recovery fidelity satisfies

\[
\boxed{r\ge\max\left\{0,\frac{9f-1-\delta}{8}\right\}.}
\]

If U and V anticommute exactly, the stronger identity r=(9p-1)/8 holds, where p=Tr(K rho)>=f. In particular,

\[
r\ge\frac{9f-1}{8}\ge\frac{3f-1}{2}
\quad\text{for }f\le1.
\]

This is a sufficient condition for a simple forward recovery to achieve or exceed the general existential target. It is not an assertion that the unrestricted target is achievable with two calls. The extra algebraic promise restricts the model class and its possible score/recovery pairs.

The balanced circuit is also used in [joint calibration recovery](JOINT_CALIBRATION_RECOVERY.md); the reflection promise here is additional.

The extraction construction and the general use of operator inequalities in recovery/self-testing have substantial prior art. The candidate contribution is the explicit identity and error coefficient for this coherent score, subject to independent review and a more specific novelty check.

## Circuit and operator calculation

Put P_+=(I+V)/2 and P_-=(I-V)/2. With an output controller qubit, define

\[
R_0=\binom{P_+}{UP_-},\qquad
R_1=\binom{UP_+}{P_-}.
\]

Each is an isometry. Start the controller in |+>, apply controlled V, then a controller Hadamard. This produces the two branches P_+,P_-. A controlled U on controller value one gives R_0, while a controlled U on value zero gives R_1. Choose the two circuits with probability 1/2, independently of the state and observed score. Discard M and the random bit. Resources are one controller qubit, one fair classical bit, one controlled V and one controlled U, plus controller gates.

For an isometry R=(B_0;B_1), its Bell acceptance effect has blocks E_ab=B_a-adjoint B_b/2. Thus the balanced map has

\[
E_{00}=P_+/2,\quad E_{11}=P_-/2,\quad
E_{01}=P_+(U+U^\dagger)P_-/4,\quad E_{10}=E_{01}^\dagger.
\]

Let

\[
L=I+X\otimes U+Z\otimes V,\qquad K=L^\dagger L/9,
\quad D=P_+UP_++P_-U^\dagger P_-.
\]

Expanding L-adjoint L using unitarity and V^2=I yields the exact identity

\[
\boxed{
E-\frac{9K-I}{8}
=-\frac14\begin{pmatrix}0&D\\D^\dagger&0\end{pmatrix}.
}
\]

It holds as an operator on reference-qubit tensor M, for every state; no reference-marginal condition enters the algebra.

## Exact defect norm

Relative to the two eigenspaces of V, write

\[
U=\begin{pmatrix}A&B\\C&D_0\end{pmatrix}.
\]

The operator D in the identity is block diagonal with blocks A and D_0-adjoint. Meanwhile UV+VU is block diagonal with blocks 2A and -2D_0. Therefore

\[
\|D\|=\frac12\|UV+VU\|=\delta/2.
\]

A Hermitian off-diagonal block matrix (0,D;D-adjoint,0) has norm ||D||. Consequently,

\[
\left\|E-\frac{9K-I}{8}\right\|=\frac{\delta}{8}.
\]

Taking an expectation in rho gives the lower bound in the result, using f<=p. Equality of this operator norm is exact; it does not claim that the lower fidelity bound is simultaneously tight at every fixed measured f and delta.

If delta=0, D=0, proving the stated recovery identity directly. U need not be Hermitian: the promise is one reflection and exact anticommutation, not two reflections.

## Practical interpretation and limitations

At f=0.99 and delta=0.01, the bound gives r>=0.9875. The unrestricted existential target is 0.985. Thus this is a sufficient condition under which a two-call physical recovery achieves the target with an explicit allowance for an algebraic control defect.

More generally the sufficient condition to meet that target is delta<=3(1-f). This is a sufficient condition, not a necessity or a complete classification. A device failing it may still have an excellent recovery.

The quantity delta is an operator norm on the addressed subsystem. It is additional information about the interventions; it is not inferred from f alone and not automatically estimated by a few control calibrations on selected states. In an unknown large subsystem, establishing a uniform norm bound may be difficult. No new experimentally accessible estimator of delta is proved here. The exact-reflection condition on V is also substantive. If V is only approximately a reflection, the present formula cannot be used unchanged.

There is no contradiction with the adversarial three-query obstruction: that adversary does not satisfy these sufficient conditions. Nor does this establish that the anticommutator is the only obstruction in the general case.

## Verification

`proofs/verify_reflection_identity.py` proves the block-polynomial identity in exact SymPy algebra with independent formal entries for U and its adjoint. The diagonal simplification uses the explicitly stated unitarity relations. Because the remaining expression is linear in the blocks, this calculation represents arbitrary-dimensional blocks rather than a dimension-specific polynomial identity.

`tests/verify_reflection_numeric.py` separately rebuilds the actual two isometries and their Born effects for 120 random unitary/reflection pairs in dimensions 2,3,4,6. The operator identity residual was approximately 2.3e-16, and the norm equality residual approximately 1.7e-16. These numerical tests check signs and normalization; the block proof supplies arbitrary-dimensional validity.

Review priorities: check the factor 1/2 in the Bell effect, the averaging of R_0 and R_1, the off-diagonal D expression, and the precise meaning of the norm bound. Do not market an unmeasured delta as a device-independent conclusion.
