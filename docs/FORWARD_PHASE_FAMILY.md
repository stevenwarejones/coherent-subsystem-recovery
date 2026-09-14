# A restricted phase family does not prove the general obstruction

This diagnostic is distinct from the universal robust obstruction.

## 7. Why phase-only adversaries cannot settle the problem

Consider rho=Phi, U=e^(iu)X and V=e^(iv)Z. Put a=cos(u), b=cos(v), c=cos(u-v). The test score attainable with the identity decoder is

\[
f=(3+2a+2b+2c)/9,
\qquad (3f-1)/2=(a+b+c)/3.
\]

The following isometries transfer M to an output qubit plus a discarded copy of M:

\[
R_0=\frac12\binom{I+V}{U(I-V)},\quad
R_V=\frac12\binom{V(I+V)}{U(I-V)},\quad
R_U=(H_{\rm Had}\otimes I)\frac12\binom{U(I+U)}{V(I-U)}.
\]

R_0 is implemented by preparing the controller in |+>, applying controlled V, then a controller Hadamard, then controlled U. R_V adds V controlled on output value zero; R_U swaps U,V and adds a final output Hadamard. R_0 is the two-call coherent extraction. R_V is its three-call echo, and R_U is the swapped echo with a fixed output Hadamard. Their implementation uses only controlled forward calls and controller gates. Their Bell fidelities on this family are respectively

\[
r_0=(1+a)(1+b)/4,\quad r_V=(1+b)(1+c)/4,\quad r_U=(1+a)(1+c)/4.
\]

Choose these three channels with equal probability, independently of u,v and f. Then

\[
\bar r-\frac{3f-1}{2}
=\frac{(1-a)(1-b)+(1-a)(1-c)+(1-b)(1-c)}{12}\ge0.
\]

The nonnegative factors prove the statement for every phase pair, not just a numerical grid. `research/check_forward_phase_family.py` verifies each channel's phase-family isometry, Born fidelity, and the Laurent-polynomial margin identity exactly. The full universal guarantee still fails, as the five-slot obstruction demonstrates. This is a diagnostic for choosing adversaries, not a competing universal recovery theorem.
