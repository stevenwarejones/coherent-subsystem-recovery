# Independent reconstruction of the forward-query additions

Scope: the joint score/calibration result, robust controller obstruction, finite
history formulation, reflection identity, restricted phase diagnostic, completion
examples, and historical endpoint certificates in this change. This is an AI
mathematical/code review, not human expert review, peer review, or a Lean proof.
“Independent” here means a new derivation and implementation with no imports from
the supplied verifier or discovery modules. It is not a second research team.

## Separate computational routes

`tests/verify_forward_results_independent.py` reads only the certificate data and
uses SymPy's free-group operations to multiply channel/calibration operators. It
constructs the joint Bell effect from the two isometries, the calibration from all
four sequential outcomes, and K by matrix multiplication of L. It checks trace
preservation before comparing the resulting operator with the supplied Gram sum.
This does not reuse the tuple-word simplifier in the supplied checker.

For positivity it uses FLINT characteristic polynomials, rather than rational
pivot elimination. For a real symmetric M, all eigenvalues are real. The
coefficients of det(tI+M) are nonnegative iff M is positive semidefinite: necessity
follows by multiplying (t+lambda_i), and sufficiency follows because a negative
lambda_i would give a positive root of a monic polynomial with nonnegative
coefficients. The multiplicity of zero determines the nullity. This checks the
joint Gram's rank 12, all moment matrices, the full robust residual, and all four
historical endpoint residuals via a different exact algorithm.

The robust objective is reconstructed as a repeated moment Gram after Bell
contraction. The ideal objective is a Gram of actual Pauli-word amplitudes. Dual
multipliers are lifted as adjoints of recursive partial sums using Kronecker
products. Historical endpoint kernels are reconstructed from all three wrong-Bell
amplitudes, before restricting the lifted dual. The supplied kernel and dual
builder are not imported.

The same checker independently verifies the reflection identity in the quotient
algebra with V^2=I, Pauli covariance for E/K and sign invariance of Q in the full
free-unitary algebra, the phase-family fidelities by explicit coherent gates,
and the two completion examples. Sample-count rounding uses integer powers on
both sides of each proposed cutoff, not a floating-point ceiling.

`tests/verify_forward_physics_independent.py` uses an eigensquare-root of the
moment matrix and QR completion of partial isometries. The supplied code instead
uses Cholesky/SVD. It applies prepare/controlled-V/H/controlled-U gates directly,
computes all four calibration outcomes, and constructs a common-decoder success
subspace. It tests 18 independently generated complex controllers against the
flagged robust objective and 24 complex, rotated reflections. These numerical
checks corroborate the derivations below; they are not universal proofs.

## Physical moment realization

Let M be a Gram matrix of vectors x_(a,w), with suffix-closed positive words W.
For each generator g and each w with gw in W, map x_(a,w) to x_(a,gw). The exact
consistency equations state that all inner products on this domain are preserved.
If a linear combination of domain vectors is zero, its squared norm is zero;
the image combination has the same squared norm and is zero. Thus the map is
well defined even for singular Gram matrices, and is an isometry of its domain
onto its image. Equal-dimensional orthogonal complements in the same finite
span allow a unitary extension. The two generators can be extended independently.
Induction on word length, deleting the leftmost letter, realizes every stored
word. The two empty-word squared norms sum to one, normalizing the joint vector.

The new numerical implementation realizes the shipped nonsingular adversaries.
The preceding argument also covers singular matrices in the general formulation;
singular cases are a prose conclusion, not a numerical test of every matrix.

## Why the support score becomes an actually attainable coherent score

Write the pure encoded state as a matrix Psi with reference rows and MB columns.
After the coherent query, use CA rows and MB columns. The ideal CA support
isometry is J=(I,X^T,Z^T)^T/sqrt(3). Projecting onto it gives
A=J†Psi_after, with C=A A† and Tr C=p. Any later decoder has f<=p because it cannot
change the CA support weight.

For the four flagged reference Pauli rotations and matching intervention signs,
L_g(g tensor I)=(g tensor I)L. The averaged reference marginal is I/2 and the
averaged C is (p/2)I. An encoding purification therefore exists. When p>0, the
two rows D0=sqrt(2/p) conjugate(A) are orthonormal. They can be completed to an
isometry acting on MB alone with output qubit and discarded environment. Moreover
A D0^T=sqrt(p/2)I, so its Bell/target success amplitude is sqrt(p). Orthogonal
completion rows contribute zero. Hence f=p is attained with a common decoder;
no branch label or reference access is given to it. For p=0 the support bound
already forces zero. This supplies the needed physical bridge without simply
identifying an arbitrary unflagged support score with a measured score.

Direct-sum mixtures of flagged instances retain scalar C, so their measured
scores are the weighted averages. Off-diagonal coherence is removed by an
orthogonal purifying flag in the bypass as in the stated construction. Flags
inside M are not made accessible to the restricted recovery.

## Controller constraints: necessity and converse

Purify all controller-only channels, including stored classical outcomes. For
history h, let c_(a,e,h) be the final output/environment coefficient. Its outer
product gives a positive history matrix T with trace one. Before the last
controller isometry, histories with different final query bits occupy orthogonal
control subspaces. Isometries preserve this orthogonality. Tracing output and
environment therefore makes the cross-bit blocks zero. Summing equal-bit blocks
recovers the preceding histories' Gram. Repeating proves every causal equality
used by the robust dual. The Bell projection produces the coefficient objective
with the factor 1/2; the real rational adversary makes taking the real part of a
complex feasible T harmless to its objective and positivity.

For the converse, factor the successive causal Grams into vector families
z_h^(t). At each stage the two child spans are orthogonal, and the sum of their
Grams is the parent Gram. Therefore the map
z_h^(t-1) -> z_(h0)^(t)+z_(h1)^(t) preserves inner products and is a well-defined
isometry. Embed the child spans in the zero/one control subspaces, using enough
controller workspace, and extend this isometry to controller gates. The next
controlled oracle splits these components as required. After the final call,
an equal-Gram isometry maps to the specified output/environment coefficients.
No operation on M other than the specified controlled query is introduced.

Only necessity is needed for the robust upper bound. The converse is used for
the complete finite-pattern search and the compactness existence argument.

## Finite support and dimension-dependent designs

For a fixed pattern, its distinct positive history words are suffix closed.
Every feasible moment matrix has the finite unitary realization above. Thus
minimizing the relevant hereditary polynomial over these matrices is an exact
physical minimization, not merely a hierarchy relaxation. The feasible moment
set is compact: normalization plus cancellation bounds all diagonals and PSD
bounds every entry. M=I/2 is strictly positive and satisfies the affine moment
equalities, so finite-dimensional SDP duality applies in their affine space.
This justifies the support-specific complete formulation. It is not a claim
that an arbitrary chosen SOS support is complete or that numerical solvers
reliably decide boundary feasibility.

For fixed score-independent circuits, the allowed Pauli output/controlled-sign
symmetrization and ideal-block mixing turn a universal high-score sharp guarantee
into the operator inequality. The same circuit must be used across the mixture;
this step does not cover score-selected families. At f=1 the sharp guarantee
forces ideal recovery one, permitting exact cancellation of the ideal term.
For dimension-dependent designs with uniformly bounded query count, choose an
unbounded sequence of dimensions using one pattern. The finite history matrices
have a convergent subsequence in a compact PSD trace-one set. Embedding each
fixed test plus flags and ideal block into those dimensions, then taking the
limit, gives the same inequality. The converse realization turns the limit into
a finite controller. This is an existence argument, not a numerical compiler.

## Robust high-score interval and active-call coverage

At f=0.9 the independent full residual is positive definite and bounds every
causal T by alpha=1692197/2000000. At higher f, dilute that fixed physical
instance with the ideal block by lambda=10(1-f). Every controller has ideal
success at most one, so r_f <= 1-lambda+lambda alpha. This yields the stated
slope 307803/200000 and positive gap for f<1. It holds pointwise for every
controller, allowing coefficients selected using f without endpoint perfection.

Every mixed order of at most three letters is a subsequence of VUVUV; this was
also exhaustively checked. A controller using only one unitary is a
measure-and-prepare channel from M: diagonalize the unitary, and after tracing M
only its eigenvalue label can affect the output. Thus it is entanglement breaking
and its Bell overlap is at most 1/2. Preselected classical mixtures obey the same
upper bound by linearity. No assertion about arbitrary outcome-dependent or
indefinite order follows from the scoped argument used here.

## Constructive bounds and statistical quantifiers

The exact joint identity is E-5K/4+3Q+I/4 >= 0. Taking an expectation and using
p>=f gives the claimed affine lower bound and q<=(1-f)/12 sufficient condition.
The independent moment values give a positive deficit at q/(1-p)=522171/4500256.
The same flag and ideal-mixture construction preserves these values and yields
the entire high-score counterexample family. This limits the balanced recovery,
not arbitrary alternative circuits. A failing example at that ratio bounds the
supremum of sufficient constants; it does not identify that supremum exactly.

For the reflection identity, relative to V's eigenspaces write U=[A B;C D]. The
correction block is diag(A,D†), while UV+VU=diag(2A,-2D). Their norm ratio is
exactly two, including zero-dimensional eigenspaces. This proves the operator
norm coefficient after the independently checked quotient identity. Approximate
reflections require an additional argument and are not covered.

With iid, stable calibration trials and known f, zero agreeing outcomes have
probability (1-q)^N. The one-sided upper limit is 1-alpha^(1/N). The sufficient
counts 358, 3594, and 35948 at f=.9,.99,.999 were checked by exact integer
inequalities, including failure at N-1. For the lower bound, a randomized test
accepting the q=0 instance with probability >=1-eta must accept its all-zero
record with that probability. On the failing instance its false acceptance is
at least (1-eta)(1-q)^N. Requiring this <=eta gives the stated bound and the
necessary 2537-trial cutoff at f=.99, eta=.05. This concerns fixed-N,
agreement-only records at known f; it is not a total experimental lower bound.
A Pauli/Bell mixture supplies the good q=0 example with the same score because
r=(9f-1)/8 and f is attained by the scalar-C construction.

When f is estimated, a lower recovery bound uses f_L and q_U. Certification
against the unknown true affine target requires comparison with f_U, not f_L.
Coverage must be simultaneous, for example by allocating failure probabilities
and applying the union bound. No independence between confidence procedures is
needed for that union bound.

## Disposition

The independent exact and physical reconstructions substantiate the included
bounds within their stated models. They do not formalize all prose, establish
human agreement, prove historical novelty, solve every finite-query count, or
identify optimal threshold constants. The LCU/adjoint-access compiler and its
control-error guarantee remain excluded because their compiler prerequisite is
not part of the current repository. Comparisons with the missing predecessor
calibration derivation remain deferred.

## Recorded local outcome

Both independent entry points passed. The numerical route's maximum residual
was approximately 6.97e-12. The supplied core suite, historical witness suite,
new physical corroborations, existing independent theorem/forward-curve routes,
and existing supporting mutation suite also passed. Optimized-Python refusal
was checked for the new assertion-dependent entry points. Certificate data
was preserved. Remote CI status must be assessed separately on the PR commit.
