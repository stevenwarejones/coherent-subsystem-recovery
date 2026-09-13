# Protocol, access model, and score definitions

This file fixes the exact experiment and the meaning of every symbol. The theorem
statements in `RESULTS.md` are only meaningful relative to this model; changing the
access assumptions changes what is proved.

## Systems

- **Q** — a trusted input qubit.
- **A** — a protected reference qubit, initially maximally entangled with Q. So the
  pre-experiment reference marginal is `rho_A = I/2`.
- **M** — the *addressed* internal subsystem, of arbitrary finite dimension. This is the
  region whose recoverable content the theorem is about.
- **B** — an arbitrary finite-dimensional **quantum bypass**. Its dimension is not bounded
  by that of M. This is the crucial feature of the model: a later decoder may use B.
- **C** — a protected three-level controller.

An arbitrary finite-dimensional encoding isometry sends the input into `M ⊗ B`. Let
`rho_AM` denote the reduced state of reference-plus-addressed-subsystem **before** the
intervention.

## The coherent test

The controller C begins in the uniform superposition `(|0>+|1>+|2>)/sqrt(3)`. Its three
branches coherently apply `I, U, V` **to M alone**, with equal amplitudes. `U` and `V` are
arbitrary finite-dimensional unitaries — not assumed Hermitian, Pauli, commuting, or
anticommuting.

A single subsequent CPTP **decoder** may act on `M ⊗ B`. It may **not** access A, the
controller C, or any additional branch label. There is **no postselection**.

The ideal target on controller-reference-output is
```
|Omega> = (1/sqrt 3) * sum_{a=0,1,2} |a>_C (I_A ⊗ P_a) |Phi+>_{AO},   (P_0,P_1,P_2)=(I,X,Z).
```

## The two scores and the recovery quantity

| Symbol | Meaning | Relation |
|---|---|---|
| `f` | probability of projecting onto the full pure target `|Omega>`, discarding decoder garbage. **This is the experiment-facing score.** | `f <= p` |
| `p` | weight of the actual controller-reference marginal in the target's support projector `P`; `p = Tr(K rho)` with `K = L†L/9`, `L = I + X⊗U + Z⊗V`. | `f <= p` |
| `R` | the largest Bell overlap achievable by a CPTP recovery acting on the **pre-intervention M alone**: `R = max_{D:M->qubit CPTP} <Phi+|(id_A ⊗ D)(rho_AM)|Phi+>`. | squared entanglement fidelity at the maximally mixed input |

`R` is squared entanglement fidelity for the maximally mixed qubit input. It is **not** a
worst-case (diamond-norm) fidelity, not a fraction of stored bits, and not a percentage.
The corresponding average pure-input qubit fidelity is `(2R+1)/3`.

## What the conclusion is, and is not

The theorem is an existence guarantee about recoverable information in M *at the tested
cut*. Coherent implementation of the controlled unitaries — **including their relative
phases** — is a substantive assumption: specifying only the isolated reduced channels
`U(.)U†` and `V(.)V†` is not sufficient (an exact example is in `PRIOR_ART.md`). The trusted
restriction that the probe acts only on M is essential.

It does **not** establish an observer, a semantic record, a unique internal computation, a
realism claim, or continuous retention between tested times. It is **not** fully
device-independent: coherent control, subsystem isolation, the protected systems, and the
target test are all assumptions.
