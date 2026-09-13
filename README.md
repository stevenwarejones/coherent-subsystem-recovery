# Coherent Subsystem Recovery

**Exact sharp recovery bounds for quantum information in an addressed subsystem, from a
three-branch coherent test, with an unrestricted quantum bypass available to the decoder.**

> Status: AI-produced research artifact under adversarial cross-review. Verified in exact
> arithmetic in this repository. **Not reviewed by a human domain expert. No historical
> priority is claimed.** See [`docs/VERIFICATION.md`](docs/VERIFICATION.md),
> [`docs/PRIOR_ART.md`](docs/PRIOR_ART.md), and [`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md).

## The question

A trusted qubit is maximally entangled with a protected reference **A**, then encoded into an
addressed subsystem **M** and an arbitrary finite-dimensional **quantum bypass B**. A protected
three-level controller coherently applies `I, U, V` **to M alone**; a single later decoder may
act on `M ⊗ B` but not on A or the controller. From the coherent-test score `f`, how much
entanglement can be recovered from **M alone** — even though the decoder also holds B?

Full model and score definitions: [`docs/PROTOCOL.md`](docs/PROTOCOL.md).

## The result

For arbitrary finite-dimensional `U, V` and arbitrary finite-dimensional bypass, the minimum
recoverable Bell overlap at score `f` is
```
R_min(f) = max{ 1/4,  3f/4,  (3f - 1)/2 },        0 <= f <= 1,
```
with matching constructions attaining every point (addressed dimension ≤ 6). The nontrivial part
is the high-score segment
```
R >= (3f - 1)/2      for  2/3 <= f <= 1,   sharp,
```
proved by an **exact free-unitary positive-operator certificate** with normalization `t = 1`
(no numerical tolerance): PSD Gram witnesses satisfying, for every pair of unitaries in every
finite dimension, `Q − S = H`, `Tr_A Q = I`, `Q, S ≥ 0`, where `Q` is the Choi matrix of a
CP unital map whose adjoint is the recovery channel on M. The other two segments are elementary. Theorem statements and a
claim ledger: [`docs/RESULTS.md`](docs/RESULTS.md); full proof write-up:
[`docs/THEOREM.md`](docs/THEOREM.md).

**Mechanism.** Two coherent branches are perfectly simulable with only a classical record
(`R = 1/2` even at `f = 1`); three branches certify `R > 1/2` once `f > 2/3`; four Pauli branches
admit a simple forward-only recovery circuit (`R ≥ f`). This is what the resource count buys.

**Implementable decoder (separate addition).** An explicit two-query forward circuit — one
controlled-`U`, one controlled-`V`, no inverses, no knowledge of `U,V` — has the exact worst-case
curve `r_UV ≥ h(f) = max{0,9f−1}²/64` under an extra assumption of reusable controlled forward
access, and this is minimax-optimal for the fixed query stage plus controller-only decoding on
`[7/9,1]`. It is weaker than the existential sharp map (infidelity coefficient `9/4` vs `3/2`);
see [`docs/FORWARD_RECOVERY.md`](docs/FORWARD_RECOVERY.md).

**Forward-query limits (separate addition).** A full-history rational dual certifies a
[robust obstruction](docs/ROBUST_FORWARD_OBSTRUCTION.md) to the sharp high-score guarantee with at
most three active forward calls, including score-selected controllers with fixed or preselected
random order and arbitrary intermediate controller processing. At `f = 0.9` every such controller
has `r ≤ 1692197/2000000 < 0.85`. It does **not** exclude every finite-query algorithm; four active
calls have not been shown sufficient. The [finite-query model](docs/FORWARD_QUERY_MODEL.md) states
the controller and target access explicitly. This does not change the existential sharp curve.

## Verify

```sh
python -m pip install -r requirements.txt      # SymPy + NumPy; the core check needs neither
python run_checks.py
```

The central certificate and the mutation cases use the **standard library only** — no SDP
solver, no network. A second independent route over `python-flint`
(`proofs/verify_alternate.py`) and an independent from-scratch re-derivation
(`tests/verify_independent.py`) are run by CI; see [`docs/VERIFICATION.md`](docs/VERIFICATION.md)
for the five tiers and exactly what each does and does not establish.

## What this is not

Not human-reviewed; not a priority claim; not experimental feasibility; not device-independent;
not a statement about observers, records, or retention over time. The nearest formula in the
literature (Berta–Coles–Wehner 2014) is the *same affine expression for a different score* — an
exact counterexample shows the naive identification fails, and the sharp bound still beats the
strongest translated bound identified in our bounded search (BCW + Renes 2017) on `[2/3, 1)`. Two possible
reductions (supermap identity-comb; Ando/WEP completion) remain **open** and gate any strong
novelty claim: [`docs/OPEN_PROBLEMS.md`](docs/OPEN_PROBLEMS.md).

## Layout

```
certificates/sharp_recovery.json     the shipped rational certificate (t = 1)
proofs/verify_sharp_recovery.py      exact standalone checker (stdlib)
proofs/verify_alternate.py           second route (python-flint, optional)
proofs/verify_full_curve.py          elementary segments + mechanism checks (SymPy)
proofs/verify_forward_recovery.py    forward two-query implementable guarantee (SymPy)
tests/test_certificate_mutations.py  intended-reason corruption cases (stdlib)
tests/verify_independent.py          independent from-scratch exact re-derivation
tests/verify_numeric_instantiation.py   concrete-unitary corroboration (NumPy)
tests/verify_attainment_independent.py   attainability reconstruction (NumPy)
research/                            symbolic attainment + literature-comparison arithmetic
docs/                                PROTOCOL, RESULTS, THEOREM, VERIFICATION, PRIOR_ART, OPEN_PROBLEMS
history/                             the superseded t = 1041/1000 certificate, for provenance
```

Produced by AI systems (Claude and ChatGPT "Astra") working adversarially under the direction of
the repository owner, who is not a physicist. The verifiers, the independent re-derivation, and
this scaffold were assembled and checked by Claude; the theorem and certificate were produced by
Astra. AI assistance and the absence of external expert validation are stated plainly throughout.
