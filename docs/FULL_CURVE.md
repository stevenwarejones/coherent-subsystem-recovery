# The full recovery curve

This consolidates the proof that, under the model of `PROTOCOL.md`,
```
R_min(f) = max{ 1/4, 3f/4, (3f-1)/2 },        0 <= f <= 1,
```
is exact: each segment is a valid lower bound with a matching construction attaining it. The
high-score segment `(3f-1)/2` (the only nontrivial one) is proved by the exact certificate in
`THEOREM.md`; this file supplies the elementary segments and all the matching constructions, so
a reader need not reconstruct them from outside the repository. `proofs/verify_full_curve.py`
checks the arithmetic here; `tests/verify_attainment_independent.py` reconstructs the endpoints
independently.

Notation is that of `PROTOCOL.md`: `L = I + X⊗U + Z⊗V`, `K = L†L/9`, `p = Tr(Kρ)`, `f ≤ p`,
and `R` is the optimal Bell overlap recoverable from M alone.

## Lower bounds

### `R ≥ 1/4` (all `f`)

Discard M and prepare any fixed qubit state. Since `ρ_A = I/2`, the Bell overlap of
`(I/2) ⊗ (fixed)` with `|Φ⁺⟩` is `1/4`. So `R ≥ 1/4` unconditionally.

### `R ≥ 3f/4` — elementary certificate

Pauli orthogonality and unitarity give the exact partial trace
```
Tr_A K = (2/3) I_M.
```
(In `L†L = 3I + X⊗(U+U†) + Z⊗(V+V†) + XZ⊗(U†V−V†U)`, tracing out the qubit A kills every term
whose qubit factor is traceless — `X, Z, XZ` — leaving `Tr_A(3 I⊗I_M)/... = 6 I_M`, so
`Tr_A K = Tr_A(L†L)/9 = (2/3) I_M`.) Therefore
```
Q_0 := (3/2) K = L†L / 6  ⪰ 0,     Tr_A Q_0 = I_M,
```
and the Choi construction of `THEOREM.md` §4 (with `Q_0` in place of the sharp `Q`) yields a
CPTP recovery on M whose Bell acceptance effect is `Q_0/2`, giving
```
R ≥ (1/2) Tr(Q_0 ρ) = (3/4) Tr(K ρ) = (3/4) p ≥ (3/4) f.
```
This is the positive operator `L†L` with the correct normalization — no Gram search, no
numerics.

### `R ≥ (3f-1)/2` — the sharp certificate

Proved in `THEOREM.md` §3–4 and `certificates/sharp_recovery.json`.

The pointwise minimum is the max of the three, and the segments join continuously at `f = 1/3`
(both `1/4`) and `f = 2/3` (both `1/2`).

## Matching constructions (attainability)

Throughout, the common-decoder optimum uses Uhlmann: with `C = T†τT` the compressed
controller-reference matrix, `f_opt = (Tr√C)²/2`, and when `C = (p/2)I` this equals `p`.

### No recoverable information: `f = 1/3`, `R = 1/4`, and all of `[0, 1/3]`

Take M one-dimensional and keep the whole encoded input in the bypass B, so `ρ_AM = I_A/2`.
Choose scalar interventions `U = V = i`. Then `L_0 = I + iX + iZ` has `L_0 L_0† = 3I`, so
`p = 1/3` and `C = Tr_M(LρL†)/9 = I/6`, giving `f_opt = 1/3`. Recovery from M alone is `1/4`
(M carries nothing; `ρ_A = I/2`).

There is also an explicit common decoder: on B apply `W = L_0*/√3 = (I − iX − iZ)/√3` (unitary),
attaining score `1/3`; applying `XW` instead gives score `0`. A probabilistic mixture of these
two CPTP decoders realizes any exact `f ∈ [0, 1/3]` with the same pre-probe state and `R = 1/4`.
The randomness is internal to the common decoder — not access to a branch label — and there is
no postselection.

### The segment `1/3 ≤ f ≤ 2/3`: `R = 3f/4`

Mix the `f = 1/3` endpoint above with the classical-flag endpoint `(f, R) = (2/3, 1/2)` below,
putting the mixing flag **in M** and defining the interventions block-diagonally. Both endpoint
`C` matrices are scalar, so the mixture has `C = (p/2)I` and the allowed optimal decoder attains
`f = p`. With weight `λ` on the classical-flag endpoint,
```
f = (1 + λ)/3,     R = (1 + λ)/4 = 3f/4.
```
Recovery is exactly additive for a block-diagonal state with a readable flag in M: restricting
any CPTP recovery to each input block gives a CPTP recovery there (upper bound by the weighted
block optima), and reading the flag and applying the optimal conditional channel attains it.
Addressed dimension `1 + 4 = 5`.

### The classical-flag endpoint `f = 2/3`, `R = 1/2`

Reference state `ρ_0 = ½[I + (X − Y + Z)/√3]`, scalar interventions
`u = (√3 + i)/2`, `v = (√3 − i)/2`. With `g` ranging over `{I, X, Y, Z}` (a four-valued flag in
M), set `ρ_AM = ¼ ⊕_g g ρ_0 g†`, `U = ⊕_g s_X(g) u`, `V = ⊕_g s_Z(g) v`, where
`g X g† = s_X(g) X`, `g Z g† = s_Z(g) Z`. Direct calculation gives `ρ_A = I/2` (Pauli twirl) and
`C = I/3`, hence `f_opt = 2/3`. The state is separable (classical on the M-flag), so no recovery
from M exceeds Bell overlap `1/2`, and reading the flag and preparing the conjugate of its
associated pure reference state attains `1/2`. Thus `R = 1/2` exactly. The flag is part of the
pre-intervention addressed subsystem, not an extra controller-branch label.

### The segment `2/3 ≤ f ≤ 1`: `R = (3f-1)/2`

The reviewed mixture of the classical-flag endpoint and the Bell endpoint (`U = X, V = Z`,
`f = R = 1`). With probability `λ` on the perfect block, `f = (2 + λ)/3`, `R = (1 + λ)/2 =
(3f−1)/2`, by the same scalar-`C` additivity. Addressed dimension at most six overall.

## Mechanism: why three branches

### Two branches cannot certify `R > 1/2`, even at `f = 1` (general)

For **any** two target qubit unitaries `P_0, P_1`, diagonalize the relative unitary
`P_0† P_1 = Σ_j e^{iθ_j} |e_j⟩⟨e_j|`. Encode `|e_j⟩ ↦ |j⟩_M |j⟩_B`. The reduced state on AM is
classical on M and separable, with equally weighted orthogonal conditional reference states, so
`R = 1/2`. Let the physical branches be `I` and the diagonal `diag(e^{iθ_j})` on M; a single
decoder maps `|j⟩_M|j⟩_B ↦ P_0|e_j⟩_O`, extended to a CPTP channel on all of MB. The two branches
then coherently produce `P_0|ψ⟩` and `P_1|ψ⟩` including their relative phase, so the full target
score is `1` despite no entanglement recoverable from M alone. Same construction for unequal
amplitudes.

Thus, within this fixed-target coherently-controlled unitary-branch family **with an unrestricted
quantum bypass**, two branches cannot certify recovery fidelity above `1/2`. (This is not a
general no-go for all memory-verification protocols; two measurement times in a different
temporal-certification protocol are not two coherent branches of this one.)

### Four Pauli branches: a simple universal forward circuit

For `n` mutually Hilbert–Schmidt-orthogonal target Pauli unitaries, `L_n = Σ_j P_j* ⊗ U_j`,
`K_n = L_n†L_n/n²`. Orthogonality gives `Tr_A K_n = (2/n) I`, so `Q_n = (n/2)K_n` defines a
recovery channel with `R ≥ (n/4) p ≥ (n/4) f`. For all four Paulis this is `R ≥ f`. With ideal
physical Pauli interventions, `L_4 = I + X⊗X − Y⊗Y + Z⊗Z = 4Φ⁺` (the conjugation matters for the
`Y` sign). For the isotropic state `ρ = tΦ⁺ + (1−t)I₄/4`, `f_opt = R = (1+3t)/4`, and this is
sharp on `[1/4, 1]`. Presented as a comparison and intuition, not a novelty claim.

Here `P_j` are the **ideal target** Paulis defining the score; the **physical** interventions
`U_j` may be arbitrary unitaries (Section "elementary certificate" and `THEOREM.md` never assume
they are Pauli). The forward two-query circuit built on this structure is in `FORWARD_RECOVERY.md`.
