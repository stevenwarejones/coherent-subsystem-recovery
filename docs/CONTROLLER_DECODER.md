# Controller-only decoding is minimax-optimal (cannot reach the sharp bound)

This is the operational proof behind the controller-only obstruction summarized in
`FORWARD_RECOVERY.md`. It shows that keeping the fixed two-query stage and replacing the final
Bell-basis decode with an *arbitrary* controller-only channel — even one chosen from the observed
score — cannot beat the fixed-circuit curve `h(f) = max(0, 9f−1)^2/64` on `[7/9, 1]`. Verified by
`proofs/verify_controller_decoder.py`; an independent spectrum re-derivation is in
`tests/verify_forward_curve_independent.py` (scope note at the end).

## 1. Architecture and quantifiers

Prepare two controller qubits in `|++⟩`, apply the same controlled `V` then controlled `U` as the
existing circuit — branches `(I, U, V, UV)` with equal amplitudes — then **discard M** and allow
an arbitrary CPTP channel `Λ` from the 4-dimensional controller to the recovered qubit. `Λ` may use
fresh ancillas, measurements, feed-forward and randomness; it has no further access to `M`, the
reference, the bypass, or an oracle. There is no postselection. For the score-dependent statement,
`Λ = Λ_f` may depend on the scalar score `f` but **not** on information distinguishing the two
opposite-phase instances below. This is a *stronger* decoder allowance than one fixed map.

## 2. The Bell-basis coordinate change is without loss of generality

Apply the fixed Bell-basis change `F` on the two controller qubits first (call them `a, c`), then an
arbitrary `Λ: ac → o`. Because `Λ` is arbitrary and `F` is a fixed unitary, composing `F` into the
front of `Λ` covers exactly the same set of maps — so working in the `a, c` coordinates loses
nothing. `verify_controller_decoder.py` builds `F` and checks `F†F = I` (`BELL_BASIS`).

## 3. The Choi objective

Let `J` be the unnormalized Choi matrix of `Λ`, ordered output `o`, input `a`, input `c`. CP and TP
are `J ⪰ 0` and `Tr_o J = I_{ac}`. For a controller state `σ_{Aac}` the recovered Bell overlap is
```
r(Λ, σ) = ½ Tr(σ^T J),
```
with the **full matrix transpose** under the fixed identification of `A` with the output index `o`
(this follows by expanding `⟨Φ⁺|(id_A ⊗ Λ)(σ)|Φ⁺⟩` with `J_{(o,i),(o',j)} = Λ(|i⟩⟨j|)_{o,o'}`; it is
not a partial transpose). Definitions `E`, `O`, `ρ0`, `ρ̄`, `Q_x` below are exactly those checked in
the verifier.

## 4. Endpoint uniqueness (`f = 1`)

At the ideal instance (`ρ_AM = Φ`, `U = X`, `V = Z`), after the fixed query stage **and the
Bell-basis coordinate change**, tracing `M` gives exactly
```
ρ0 = Φ_{Aa} ⊗ I_c/2
```
(`IDEAL_STATE` in the verifier reconstructs this from the actual query isometry, not from a chosen
formula). With the PSD endpoint slack
```
D0 = ¼ I_{oac} − ρ0^T/2 = ¼ (I − Φ)_{oa} ⊗ I_c ⪰ 0,
```
trace preservation gives `Tr(J)/4 = 1`, so `r = 1` forces `Tr(D0 J) = 0`. For PSD matrices this
puts `range(J) ⊆ ker(D0) = span{|Φ⟩_{oa}} ⊗ H_c`, hence `J = Φ_{oa} ⊗ S_c` with `S_c ⪰ 0`; the
output partial trace then forces `S_c = 2 I_c`, i.e.
```
J = 2 Φ_{oa} ⊗ I_c,   the Choi matrix of Λ = Tr_c.
```
So the final map achieving `r = 1` is unique — the partial trace, no remaining freedom. The
full-rank gauge `I_c/2` is what removes any unused-subspace loophole. This is the standard
noiseless-subsystem / operator-quantum-error-correction rigidity (Kribs–Laflamme–Poulin–Lesosky,
[quant-ph/0504189](https://arxiv.org/abs/quant-ph/0504189), Lemma 2.3 / Theorem 2.5).

## 5. Two-instance obstruction (score-dependent decoding, `7/9 ≤ f ≤ 1`)

Take `ρ_AM = Φ`, `U_± = e^{±it}X`, `V_± = e^{±it}Z`, `x = cos t ∈ [1/2, 1]`, both with score
`f = (5+4x)/9 ∈ [7/9, 1]`. In the Bell-basis controller coordinates the change from the ideal oracle
pair to the phase pair acts (up to a global phase that cancels in the state) by
```
D_± = E ∓ i (sin t /2) O,   E = (1+x)/2 · I_{ac} + (1−x)/2 · Y_a Y_c,   O = X_a X_c + Z_a Z_c,
```
so the average of the two instances' `Aac` states is
```
ρ̄ = (I_A ⊗ E) ρ0 (I_A ⊗ E) + (1−x²)/4 · (I_A ⊗ O) ρ0 (I_A ⊗ O).
```
`ρ̄` is **derived from the physical query isometry** at both phase signs, not posited; the verifier
checks it against the actual isometry at exact rational phases (`MATCHING_PRIMAL` and the
rational-phase reconstructions). With the Hermitian controller matrix
```
Q_x = (1+x)²/16 · I_{ac} + (1−x²)/16 · Y_a Y_c,
```
the `8×8` slack `D_x = I_o ⊗ Q_x − ρ̄^T/2` has the exact spectrum
```
{ 0,  (x+1)/8,  x/4,  (x+1)(2x−1)/8 },   each with multiplicity 2,
```
which is `≥ 0` on `[1/2, 1]` (the only sign-changing factor, `2x−1`, is `≥ 0` there). So `D_x ⪰ 0`,
and for every CPTP `J` (`Tr_o J = I`),
```
½(r_+ + r_-) = ½ Tr(ρ̄^T J) ≤ Tr[(I_o ⊗ Q_x) J] = Tr(Q_x) = (1+x)²/4 = h(f).
```
At least one of the two instances therefore scores `≤ h(f)`. The existing partial-trace decoder
attains `h(f)` on each sign separately (`MATCHING_PRIMAL`), so this is a matching primal/dual pair:
the architecture-wide minimax equals `h(f)` on `[7/9, 1]`. The dual-certificate construction is the
optimal-recovery SDP of Fletcher–Shor–Win ([quant-ph/0606035](https://arxiv.org/abs/quant-ph/0606035),
Section IV) applied to this two-instance problem.

## 6. Why this excludes the sharp implementation for this architecture

On `[7/9, 1)`,
```
(3f−1)/2 − h(f) = 3(1−f)(27f−11)/64 > 0,
```
so this architecture cannot reach the sharp high-score bound even with a score-dependent final map.
The optimal `M`-only recovery is still `R = 1` on both adversarial instances — this is a limitation
of the query-and-discard architecture under unknown phases, not of recovering the subsystem. Below
`x = 1/2` the certificate stops; no claim is made below `f = 7/9`.

## 7. Scope

This concerns the **fixed query stage followed by controller-only decoding**. It is not a lower
bound for all two-query circuits, does not cover a phase-informed map, and does not cover any
recovery that retains access to `M`. Reversing or classically randomizing the query order gives the
same `h`; changing only the four preparation amplitudes (with the query words and Bell decoder
fixed) cannot meet the `f = 1` requirement except as the original circuit up to global phase.
Enlarging the model — an intermediate controller operation between the two calls, additional
queries, coherent order control, or recovery that still touches `M` — is left open in
`OPEN_PROBLEMS.md`.

**Verification scope.** `proofs/verify_controller_decoder.py` reconstructs `ρ0` and the phase-pair
states from the actual query isometry (so `ρ̄` is tied to the physics), and checks the endpoint
kernel, the unique TP Choi, the exact `D_x` characteristic polynomial, and matching primal/dual.
`tests/verify_forward_curve_independent.py` independently recomputes the spectrum of the *supplied*
`ρ̄`/`Q_x` formulas and the gap; it does not re-derive `ρ̄` from the query stage — that derivation
is the primary checker's job. Neither is human review.
