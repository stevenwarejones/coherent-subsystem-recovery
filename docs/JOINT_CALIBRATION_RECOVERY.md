# Joint score–calibration recovery: linear tolerance and a matching scaling obstruction

Integrated from the 2026-09-13 research handoff. Exact certificate checks and
numerical physical corroborations are included here; see the separate exact and
numerical checks documented in [VERIFICATION.md](VERIFICATION.md). The combined
independent reconstruction is supplied in the third PR of this stack. Human
expert review and historical novelty remain outstanding.

## Integration status

The rational lower certificate and adversary are checked by the repository verifiers.
The physical reduction and statistical arguments below remain prose proofs for review;
no human review or novelty clearance is claimed.

**Deferred comparison.** A preceding quadratic calibration bound is referenced in
places below as handoff provenance only. Its derivation is not shipped in this
checkout, so that comparison — its pointwise envelope and its comparative sample
counts — is excluded until its proof is supplied, and is not a verified claim here.
The new affine bound `r ≥ max(0, (5f−1)/4 − 3q)` and its own zero-event sample
counts stand on their own and do not depend on that predecessor.

## Main result

For arbitrary finite-dimensional unitary interventions U,V, the previously defined balanced two-forward-query recovery satisfies

\[
\boxed{r\ge 1-\frac54(1-f)-3q.}
\tag{1}
\]

Here f is the existing coherent score, r is the Bell overlap probability of the actual recovery output, and q is a controller-only calibration probability on separate copies. Negative lower bounds can be replaced by zero. The same fixed recovery is used for every f,q; its coefficients need not be optimized from the data.

The sufficient condition for this physical recovery to achieve the sharp existential target is therefore

\[
\boxed{q\le\frac{1-f}{12}\quad\Longrightarrow\quad r\ge\frac{3f-1}{2}.}
\tag{2}
\]

This does not prove that (1) is the optimal joint bound or that 1/12 is the
optimal threshold constant. The handoff compares it to a preceding quadratic
calibration bound whose derivation is absent from the supplied bundles and
current repository; that comparison is deferred until its proof is integrated.

## Definitions and protocol

Let rho_AM be normalized, A a reference qubit, M the addressed subsystem. U,V are arbitrary unitaries on M. No reflection or anticommutation assumption is made. The state may be purified by an inaccessible bypass. In the Bell-input application rho_A=I/2, though the operator bound does not require that condition.

Write

\[
L=I+X_A\otimes U+Z_A\otimes V,\qquad K=L^\dagger L/9,
\qquad p=\operatorname{Tr}(K\rho).
\]

Use the base theorem's bridge f<=p. This bridge is an explicit inherited premise of the operational interpretation, not a new consequence of the Gram certificate.

Set P=(I+V)/2 and N=(I-V)/2. These are generally Kraus operators, not projectors. Choose with equal probability the two isometries

\[
R_0=\binom{P}{UN},\qquad R_1=\binom{UP}{N}.
\]

Each requires a controller qubit in |+>, controlled V, a controller Hadamard, then controlled U on controller value one or zero, respectively. Discard M; retain the controller qubit. The random choice is independent of the input. Since P-dagger P+N-dagger N=I, both are isometries without requiring V to be Hermitian.

For calibration, on fresh copies of the same state, implement the (P,N) instrument, measure its controller bit, apply U, implement (P,N) again, and measure its second bit. q is the unconditional probability that the two bits agree:

\[
q=\operatorname{Tr}(Q_M\rho_M),\quad
Q_M=(PUP)^\dagger(PUP)+(NUN)^\dagger(NUN).
\]

This costs two controlled V calls and one U call per calibration copy, with controller preparation/reset/readout. No reference or bypass access is needed. For general V the bits are instrument outcomes, not V eigenvalues. There is no postselection in either the calibration probability or the recovery.

Let E be the Bell acceptance effect of the balanced recovery. Its blocks are

\[
E_{00}=P^\dagger P/2,\quad E_{11}=N^\dagger N/2,\quad
E_{01}=P^\dagger(U+U^\dagger)N/4,
\]

with E_10=E_01-dagger. Then r=Tr(E rho).

## Exact proof of the joint inequality

The shipped rational matrix G is 14 by 14, positive semidefinite, with rank 12. Its word list is

\[
\mathcal W=(I,V,U,UV,VU,V^2,VUV).
\]

For a in {0,1} and w in W, define J_(a,w)=<a| tensor w, an operator from A tensor M to M. In the listed order the certificate proves

\[
\sum_{i,j}G_{ij}J_i^\dagger J_j
=E-\frac54K+3(I_A\otimes Q_M)+\frac14I.
\tag{3}
\]

All coefficient arithmetic is rational. Words are reduced only by adjacent inverse cancellation, using unitarity. Neither a dimension-specific matrix identity nor Hermiticity of U,V is used. The appearances of V^2 and VUV in the proof do not add calls to the two-query recovery; they are verification words, and VUV already occurs in the calibration effect.

Positive semidefiniteness of G makes the left side of (3) positive semidefinite in every dimension. Taking its expectation gives

\[
r\ge \frac54p-3q-\frac14\ge\frac54f-3q-\frac14,
\]

which is (1). Subtracting (3f-1)/2 leaves (1-f)/4-3q, proving (2).

`proofs/verify_joint_recovery.py` independently rebuilds E, K, and Q from the formulas above, verifies the exact free-word identity on 30 nonzero matrix-word coefficients, and checks G by exact semidefinite elimination. Five mutations test coefficient changes, a word-list change, loss of positivity, and a positive Gram perturbation that breaks the identity. No optimizer or numerical tolerance is involved in acceptance.

The discovery SDP reported some `optimal_inaccurate` statuses. Those numerical values are not the proof and are not claimed as optima. The final rational certificate is valid regardless of those discovery statuses.

## An exact counterexample: the linear order cannot be improved for this recovery

A separate positive-definite 14 by 14 rational moment matrix yields

\[
p_*=\frac{1265617}{1406250},\quad
r_*=\frac{424947}{500000},\quad
q_*=\frac{58019}{5000000}.
\]

Numerically these are p_* approximately 0.8999943111, r_*=0.849894, q_*=0.0116038. Its deficit is exactly

\[
\frac{3p_*-1}{2}-r_*=\frac{731}{7500000}>0.
\]

Moreover

\[
\kappa:=\frac{q_*}{1-p_*}
=\frac{522171}{4500256}\approx0.1160313991<\frac18.
\tag{4}
\]

Thus replacing the sufficient coefficient 1/12 by 1/8 would be false for the balanced recovery. This counterexample is **not** an upper bound on arbitrary alternative recovery circuits.

### Why the moments correspond to a physical instance

The word list is suffix closed. The certificate is the Gram matrix of vectors v_(a,w), with normalized empty-word components and the exact unitary consistency equations

\[
\langle v_{a,xw},v_{b,xv}\rangle
=\langle v_{a,w},v_{b,v}\rangle
\]

whenever the displayed words belong to W, for x=U,V. These equations make the prescription v_(a,w) -> v_(a,xw) a well-defined isometry on its domain span. Extend each isometry to a unitary on the 14-dimensional Gram span. Suffix closure then proves that every supplied vector equals w v_(a,I). The resulting state |0>v_(0,I)+|1>v_(1,I) realizes all scored moments exactly.

This initial state's reference marginal need not be I/2. Add a four-valued classical flag within M. For each reference Pauli g, rotate the reference by g and change U,V signs according to gXg-dagger=s_X X and gZg-dagger=s_Z Z. Average the flagged blocks. The reference is then maximally mixed.

K and E transform covariantly under these sign/reference operations: changing U's sign flips the off-diagonal reference blocks of E, while changing V's sign swaps P,N and the reference blocks. q is unchanged by either sign operation. Thus p,r,q remain the displayed values.

For C=Tr_M(L rho L-dagger)/9, the same twirl gives C=(p_*/2)I. The base common-decoder formula f_opt=(Tr sqrt(C))^2/2 therefore gives f_opt=p_*. Purifying the state into the bypass gives a legitimate Bell-input encoding. This uses the existing decoder model; it does not assert attainability under a more restricted decoder interface. The flagged target dimension can be 56.

The code checks the moment consistency and scored values exactly. `tests/verify_joint_born.py` also constructs concrete unitary extensions numerically and recomputes the flag-averaged Born values and scalar C. Those numerical checks corroborate the bridge; the finite extension and covariance arguments above establish it mathematically.

### A counterexample at every sufficiently high score

Direct-sum mix this physical instance with an ideal Bell/Pauli block. At any f in [p_*,1), use its weight lambda=(1-f)/(1-p_*). Both C matrices are scalar, so the common-decoder score is exactly f. The actual balanced recovery obeys

\[
q_f=\kappa(1-f),\qquad
\frac{3f-1}{2}-r_f
=\frac{731/7500000}{1-p_*}(1-f)>0.
\tag{5}
\]

Therefore a universally sufficient calibration ceiling growing faster than linearly in 1-f cannot hold for this fixed recovery. For example, a ceiling proportional to sqrt(1-f) would eventually include these failing instances.

Define c_bal as the supremum of constants c for which q<=c(1-f) guarantees the target for this balanced recovery on every instance with f in [0.9,1]. Then

\[
\frac1{12}\le c_{\rm bal}\le\frac{522171}{4500256}.
\]

The endpoints do not match; this is a bracket, not an exact optimum. The counterexample family only obstructs this fixed recovery. The earlier all-three-query obstruction is a different statement and must not be inferred from this one.

## Sampling benefit and its precise limits

With independent score and calibration confidence limits f_L<=f and q<=q_U, (1) gives

\[
r\ge(5f_L-1)/4-3q_U.
\]

Combine failure probabilities for both limits. This certifies the target evaluated at f_L if q_U<=(1-f_L)/12. To certify the target at the unknown true score using an interval [f_L,f_U], a sufficient test is instead

\[
(5f_L-1)/4-3q_U\ge(3f_U-1)/2.
\]

For illustration only, suppose f is known exactly, calibration trials are independent and stable, and zero agreeing outcomes are observed. A one-sided 95% binomial upper limit is q_U=1-0.05^(1/N). Required counts are:

| f | Sufficient ceiling | Zero-event trials |
|---|---:|---:|
| 0.9 | 0.00833333... | 358 |
| 0.99 | 0.00083333... | 3,594 |
| 0.999 | 0.0000833333... | 35,948 |

The calibration component scales as O((1-f)^(-1)) in this favorable zero-event
setting. These are not total experimental costs: score estimation, failures
actually observed, drift, imperfect control/readout, and repetition resources
are excluded. No apparatus pass rate is guaranteed.

### Restricted sample-complexity lower bound

There is also a matching order obstruction for **tests using only the binary agreement records at known fixed f**. A Pauli-block noisy Bell state can have the same f, q=0, and recovery above the target. Equation (5) supplies a failing instance at that f with Bernoulli agreement probability q= kappa(1-f).

Suppose a test accepts the good q=0 instance with probability at least 1-eta and falsely accepts the failing instance with probability at most eta. The all-zero record has probability one in the former and (1-kappa(1-f))^N in the latter. Even for a randomized decision rule this forces

\[
N\ge\frac{\log((1-\eta)/\eta)}{-\log(1-\kappa(1-f))}.
\]

For fixed eta<1/2 this is Omega((1-f)^(-1)). At f=0.99 and eta=0.05 the necessary count is at least 2,537, compared with the sufficient zero-event count of 3,594 above.

This lower bound assumes a fixed sample size and access only to the agreement bits plus known f. It is not a lower bound on every experimental protocol, richer controller records, alternative recovery maps, or quantum-coherent calibration. It establishes the order of the calibration bottleneck under the specified statistical interface, not global experimental optimality.

## Interpretation and remaining work

The conceptual gain is that the score constrains the same errors that appeared in the earlier calibration correction. Bounding the correction from q alone discarded that joint information and gave a square-root penalty. The exact SOS certifies a joint affine tradeoff, enabling a linear tolerance. The counterexample shows that linear order is the right target for this map and statistic.

Historical novelty is not established by use of an SOS, a small Gram matrix, or multiple observed quantities; these techniques have extensive prior art. The literature handoff identifies that overlap. The precise inequality, restricted access model, and paired upper/lower calibration scaling are the candidate contributions.

Independent review should prioritize: the free-word identity and PSD check; the interpretation of E and Q; the physical moment extension and score-attainment flags; and the decision-rule quantifiers in the sampling lower bound. The all-finite-query recovery problem, exact c_bal, and a direct experimental advantage over reference-accessible fidelity estimation remain open.
