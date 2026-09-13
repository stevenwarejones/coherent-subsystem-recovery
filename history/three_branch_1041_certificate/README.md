# Exact recovery bound from a three-branch coherent test

12 September 2026. Research result with an exact arithmetic certificate; independently written checker within the same AI session, not independent external or human review. No historical priority claim. The sharp linear curve remains conjectural.

## Main result

Under the protocol assumptions below, write f for the full coherent-test score and R for optimal entanglement fidelity after recovery from the addressed memory alone. The enclosed certificate proves

$$R\geq\frac{500}{1041}(3f-1).$$

A simpler, slightly weaker consequence is

$$R\geq\frac{12}{25}(3f-1).$$

These are dimension-independent statements, valid for every finite memory dimension and with no bound on the finite-dimensional quantum bypass. They do not require the physical interventions to be Pauli operators, Hermitian, commuting, or anticommuting. Arbitrary unitaries U,V are allowed. Infinite-dimensional domain extensions have not been audited here.

The same round's literature-based derivation adds

$$R\geq \left[2\left(\frac{1+3f}{4}\right)^2-1\right]_+^2,$$

where [x]_+ = max(0,x). This second inequality follows from Renes's published complementary-observable recovery theorem, not from a new general recovery principle. It approaches one with leading infidelity 6(1-f).

Together with the previous explicit recovery channel, we now have

$$R\geq\max\left\{\frac14,\frac{3f}{4},\frac{500}{1041}(3f-1),\left[2\left(\frac{1+3f}{4}\right)^2-1\right]_+^2\right\}.$$

The previously constructed attainable candidate is R=(3f-1)/2 for 2/3 <= f <= 1. This package does not reverify that construction. Conditional on that earlier construction, the new rational bound reaches 1000/1041, or approximately 96.06%, of the candidate value throughout that interval. This is not a statement that the exact curve is known.

| Coherent score f | Earlier bound 3f/4 | New rational bound | Renes-derived bound | Attainable candidate |
|---|---:|---:|---:|---:|
| 0.80 | 0.6000 | 0.6724 | 0.1980 | 0.7000 |
| 0.90 | 0.6750 | 0.8165 | 0.5059 | 0.8500 |
| 0.95 | 0.7125 | 0.8886 | 0.7273 | 0.9250 |
| 0.99 | 0.7425 | 0.9462 | 0.9411 | 0.9850 |

Use the maximum of the lower bounds, not one formula everywhere. None of these decimals is the certificate; the certificate is rational.

## Verification

Requires ordinary Python 3, with no third-party packages and no solver:

```bash
python3 verify.py
python3 test_mutations.py
```

Run from the extracted folder, or use an absolute path to verify.py. Do not use Python's -O option: the checker explicitly refuses it.

Files:

- `certificate.json`: two 34-by-34 rational Gram matrices, 17 free-unitary words, and the exact constant t=1041/1000.
- `verify.py`: independent of the discovery solver; exact Fraction arithmetic, positive LDL pivots with exact reconstruction, all polynomial coefficients, partial trace, and fixed claimed constant.
- `test_mutations.py`: verifies the original, then confirms rejection for four intended reasons: changed constant, failed positivity, broken SOS identity, and broken partial trace. The last mutation preserves the SOS difference and positivity, so it specifically exercises normalization.
- `README.md`: definitions, proof, comparison, and literature status.

## Protocol and scope

A trusted input qubit Q is maximally entangled with protected reference A. An unknown isometry W encodes Q into M tensor B. M is the addressed internal region; B is an arbitrary quantum bypass. The input is initially independent of any machine ancilla.

A protected qutrit controller C begins in the equal superposition. Its branches coherently apply I, U, V on M alone. The common subsequent decoder can act on MB but cannot access A, C, or an extra branch label. There is no postselection. Coherent implementations of the unitaries are part of the assumption; calibrating their separate reduced channels would not establish it.

The ideal state on CAO is

$$|\Omega\rangle=\frac{1}{\sqrt3}\sum_{a=0}^2 |a\rangle_C (I_A\otimes P_a)|\Phi^+\rangle_{AO},\qquad(P_0,P_1,P_2)=(I,X,Z).$$

Here f is the probability of projecting onto this full pure target, after discarding decoder garbage. Define the pre-probe state rho=rho_AM. R is the maximum Bell overlap achievable by a trace-preserving channel from M to a qubit. This is squared entanglement fidelity at the maximally mixed input. It is not worst-case input fidelity. Average qubit-state fidelity is (2R+1)/3.

The conclusion concerns recoverable information in M at this cut. It does not establish an observer, meaningful record, unique computation, or continuous retention between tested times. The trusted restriction that the probe acts only on M is essential.

## Exact algebraic proof

Let

$$L=I_2\otimes I_M+X\otimes U+Z\otimes V,\quad K=L^\dagger L/9,\quad H=3K-I.$$

The ideal reduced controller-reference state has support projector P=2 Tr_O(|Omega><Omega|). If tau_CA is the actual reduced state after the probe, put p=Tr(P tau). The common decoder does not alter that marginal, and |Omega><Omega| <= P tensor I_O. Expanding the controlled unitary gives

$$f\leq p=\operatorname{Tr}(K\rho).$$

Explicitly,

$$H=\frac13\left[X\otimes(U+U^\dagger)+Z\otimes(V+V^\dagger)+XZ\otimes(U^\dagger V-V^\dagger U)\right].$$

Let w_i range over the 17 reduced words of length at most two in U,U†,V,V†. With the word-major, qubit-minor matrix index (i,a), define

$$\mathcal Q=\sum_{ijab}G_{ia,jb}|a\rangle\langle b|\otimes w_i^\dagger w_j,$$

and define S in the same way from J. The certificate verifies, exactly in the free unitary algebra,

$$G\succ0,\quad J\succ0,\quad \mathcal Q-S=H,\quad \operatorname{Tr}_A\mathcal Q=t I_M,\qquad t=1041/1000.$$

Consequently Q and S are positive for every pair of unitaries of every finite dimension. To see this, apply the positive Gram matrix G tensor I_M to the column with blocks (I_A tensor w_i)|psi>; its quadratic form equals <psi|Q|psi>. The same argument applies to J. The checker reduces words only by adjacent inverse cancellation, so it assumes no extra relations between U and V.

The standard recovery SDP reads

$$2R(\rho)=\min_{\sigma}\{\operatorname{Tr}\sigma:I_A\otimes\sigma\succeq\rho\}.$$

For any feasible sigma, positivity implies

$$3p-1=\operatorname{Tr}(H\rho)\leq\operatorname{Tr}(\mathcal Q\rho)\leq\operatorname{Tr}[\mathcal Q(I_A\otimes\sigma)]=t\operatorname{Tr}\sigma.$$

Minimizing proves 3p-1 <= 2tR, and f<=p gives the stated bound. The polynomial inequality actually needs no maximally mixed marginal assumption; that marginal is needed for the stated input-channel interpretation and fidelity normalization of this protocol.

The certificate also specifies a recovery channel algebraically. Its adjoint is defined on matrix units by

$$\mathcal D^\dagger(|a\rangle\langle b|)=\mathcal Q_{ab}/t.$$

Positive Choi matrix Q/t makes this map completely positive, and Tr_A Q=tI makes the adjoint unital. Thus D is trace-preserving, and its Bell acceptance effect is Q/(2t). It has no access to B. This gives an actual channel, not only an SDP existence argument.

Implementation caveat: this algebraic channel depends on U,V and includes adjoint words. The certificate does not show a low-cost black-box implementation using only forward calls to the three original interventions. Recoverability and efficient implementability are distinct claims.

## Derivation using the published Renes inequality

[Renes, *Uncertainty relations and approximate quantum error correction*, PRA 94, 032314 (2016)](https://arxiv.org/html/1605.01420v2), Eq. (3), gives

$$\arccos\sqrt R\leq\arccos P_j+\arccos P_k$$

for optimal guessing probabilities of two conjugate reference observables from M. His fidelity convention is unsquared, hence sqrt(R) here. Theorem 1 is a stronger variant. I inspected the primary theorem, definitions, and construction; the 2010 precursor is explicitly superseded by this paper.

For our application define real correlations

$$a=\langle X\otimes\operatorname{Re}U\rangle,\quad b=\langle Z\otimes\operatorname{Re}V\rangle,\quad c=\langle Y\otimes\operatorname{Im}(U^\dagger V)\rangle.$$

These obey a+b+c=(9p-3)/2. Each memory observable is a Hermitian contraction, hence defines a binary POVM that guesses its reference Pauli outcome with probability (1+a)/2, (1+b)/2, or (1+c)/2. Optimal guessing is at least as good.

Choose the two largest correlations. Their sum is at least two-thirds of a+b+c, so the average of their optimal guessing probabilities is at least s=(1+3p)/4. Every pair of distinct Pauli observables is conjugate. Concavity and monotonic decrease of arccos on [0,1] give

$$\arccos\sqrt R\leq2\arccos s.$$

When the right side exceeds pi/2 we retain only R>=0; otherwise taking cosine and squaring yields R >= [2s²-1]_+². Replacing p by its lower bound f gives the formula at the start. This application is our derivation from a published theorem; it is not a claim of a novel general uncertainty relation.

For f=1-e, its high-score expansion is

$$R\geq1-6e+O(e^2).$$

This complements the rational bound and recovers R=1 at perfect score. It does not prove the sharper conjectured slope 3/2.

## Discovery and validation record

The smallest attempted SOS basis {I,U,V} returned a trace normalization near 4/3. That numerical result is not certified optimal and does not rule out other certificates. Enlarging to all 17 reduced words of length at most two returned t approximately 1.040678, with an accuracy warning.

I rounded the candidate matrices to rationals, repaired all polynomial and trace equalities exactly, and added the same positive diagonal shift to both matrices to reach t=1041/1000. Adding the same shift preserves their difference; the trace changes by the known amount. Both resulting matrices have strictly positive exact LDL pivots. No floating-point eigenvalue threshold is used for acceptance.

The independent standard-library verifier accepts the shipped data. Four deliberate mutations reject for their intended reasons. The discovery implementation and the verifier were written separately, but by the same AI: this is not external independent validation. An adversarial reviewer should especially check the controller-score reduction, tensor order in Q, the Choi normalization, and the application of Renes's unsquared fidelity convention.

## Literature sweep and novelty

The sweep covered noncommutative SOS certificates, unitary-moment matrices, complementary-observable recovery, and coherent memory certification.

- [Pironio, Navascués, Acín, *Convergent Relaxations of Polynomial Optimization Problems with Noncommuting Variables*, SIAM J. Optim. 20, 2157–2180 (2010)](https://epubs.siam.org/doi/10.1137/090760155): publisher record inspected. Noncommutative SDP/SOS methodology is established. This certificate is an application; no claim is made that our particular auxiliary hierarchy has been proved convergent.
- [Dykema and Juschenko, *Matrices of unitary moments*](https://arxiv.org/abs/0901.0288): primary abstract inspected. Its tracial-moment setting does not automatically characterize our arbitrary state-weighted, qubit-valued correlations. I did not rely on a theorem from this paper.
- [Renes 2016](https://arxiv.org/html/1605.01420v2): primary definitions, Eq. (3), Theorem 1 and surrounding proof inspected. This is a substantive comparator and the source of the complementary bound above. An explicit recovery map from complementary information is not a new principle.
- Earlier sweep: [An, Cao, Cui, finite-window memory recovery](https://arxiv.org/abs/2608.12803), [Rosset, Buscemi, Liang, faithful memory verification](https://arxiv.org/abs/1710.04710), and [Barizien et al., supermap self-testing](https://arxiv.org/abs/2606.25124) remain necessary comparators. Their existence prevents broad claims that memory localization or internal-operation certification is new.
- The recovery SDP itself is standard: [König, Renner, Schaffner](https://arxiv.org/abs/0807.1338).

This sweep did not locate the particular rational bound or settle whether it is implied by another existing result. Absence from the search is not priority clearance. The potential contribution is the concrete three-branch, addressed-memory guarantee with unrestricted bypass, not the SOS method or recovery principle.

## Next questions

1. External mathematical review of the supplied certificate and its physical interpretation.
2. Find a simpler exact certificate, preferably one that reaches t=1. No result here proves such a certificate exists.
3. Determine whether allowing a reference multiplier or non-polynomial recovery map removes the remaining gap.
4. Compare the score directly with stronger published recovery inequalities; a literature-derived improvement is worthwhile even if it removes a novelty claim.
5. Keep the exact bound and the conjectured sharp line clearly separated in every document.
