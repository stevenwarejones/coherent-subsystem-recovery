# Verification: what each check does, and what it does not establish

## One command

```sh
python run_checks.py          # SymPy + NumPy (declared in requirements.txt); no solver, no network
```

This runs the six required checks (exact certificate, elementary full-curve segments, the
forward two-query guarantee, the exact forward-circuit curve, the controller-only minimax
obstruction, and the mutation cases) fail-fast, streaming output. It exits nonzero on any
failure and prints its success sentence only after all six pass. The exact certificate and the
mutation cases need only the standard library; the four SymPy/NumPy checks complete the curve.

## The checks, by tier

### Tier 1 — exact, standard library only (the central result)

- **`proofs/verify_sharp_recovery.py`** — reads `certificates/sharp_recovery.json`; checks the
  reduced Gram matrices are positive-definite by exact rational elimination, rebuilds `Q` and
  `S` as free-word polynomials, and checks `Q - S = H` coefficient-by-coefficient and
  `Tr_A Q = I` with words reduced only by adjacent-inverse cancellation. Refuses `python -O`
  (assertions are the proof gates).
- **`tests/test_certificate_mutations.py`** — the unmodified certificate passes; four
  corruptions each reject for their intended reason (`CLAIMED_CONSTANT`, `POSITIVE_DEFINITE`,
  `SOS_IDENTITY`, `PARTIAL_TRACE`); rejected inputs print no success line; `-O` is refused. The
  normalization mutation adds the *same* positive summand to `Q` and `S`, so the difference and
  positivity survive and **only** the `t=1` normalization breaks — it isolates that one check.

### Tier 2 — exact, an independent re-derivation (imports nothing from the shipped scripts)

- **`tests/verify_independent.py`** — parses the certificate; recomputes `g = B R Bᵀ`; checks
  the reduced Grams PD by an independent LDLᵀ; assembles the full 106×106 block Gram matrices and
  checks **they** are PSD directly (so the operator conclusion does not rest on the congruence
  prose); rebuilds `Q - S = H` with `H` from the operator formula, and `Tr_A Q = I`.

### Tier 3 — numerical corroboration (fixed seed; catches convention/orientation bugs)

- **`tests/verify_numeric_instantiation.py`** — instantiates the free-algebra certificate on
  200 random `(U,V)` in dimensions 2–4: `Q−S=H` and `Tr_A Q=I` to ~1e-15, `Q,S` PSD. It then
  **builds the actual recovery channel's Kraus operators from `Q`**, checks trace preservation,
  applies the channel, and confirms the resulting Bell overlap both matches `½Tr(Qρ)` and attains
  `≥ (3p−1)/2` — plus a convention mutation (dropping the Kraus conjugation) that must break TP or
  the overlap. Corroboration of the conventions, not a proof; the identity in Tier 1/2 is the proof.
- **`tests/test_supporting_mutations.py`** — anti-degradation for the SymPy/NumPy supporting
  checks: tampering with a load-bearing value in the literature, forward-independent, or
  attainment scripts must make that script fail with no success line (reproduces the two false
  acceptances noted in review R2, so a demonstration can never pass as verification).
- **`tests/verify_attainment_independent.py`** — rebuilds the attaining constructions from
  scratch: classical endpoint `C = I/3` → `f = 2/3`, `R = 1/2`; perfect endpoint `C = I/2` →
  `f = R = 1`; the flagged mixtures trace `R = (3f−1)/2`.

### Tier 4 — supporting exact checks (SymPy)

- **`proofs/verify_full_curve.py`** — the elementary `3f/4` and `1/4` segments, endpoints,
  interpolation, two-branch and four-branch comparisons.
- **`proofs/verify_forward_recovery.py`** — the forward-only two-query guarantee: the one-square
  free-unitary identity `T†T+20I−4L†L=J†J` (affine tangent), the circuit's isometry/Bell-effect
  conventions on exact instances (incl. complex nonsymmetric unitaries), the phase family, and the
  old selection threshold. See `FORWARD_RECOVERY.md`. Refuses `-O`.
- **`proofs/verify_forward_curve.py`** — the *exact* forward-circuit curve `h(f)=max(0,9f−1)²/64`:
  the overlap identity `L†L=I+B†T+T†B`, the curve and its affine tangent, the switch `f_c` and
  maximum gap, 500 random matrix/circuit cases (dims 1,2,3,4,6), 101 phase instances, reversed
  order, low-score endpoint, and a changed-constant rejection control (SymPy + NumPy).
- **`proofs/verify_controller_decoder.py`** — the controller-only minimax obstruction on
  `[7/9,1]`: endpoint-unique Choi matrix, the exact `D_x` spectrum and PSD interval, matching
  primal/dual, score conversion, rational-phase reconstructions, and two rejection controls
  (SymPy; runs under `-O` via explicit exceptions).
- **`tests/verify_forward_independent.py`** — an independent from-scratch NumPy re-derivation of the
  one-square identity (random `U,V`, dims 1–5), `r_UV ≥ (9p−5)/4`, the phase family, `forward ≤ sharp`.
- **`tests/verify_forward_curve_independent.py`** — an independent re-derivation (imports nothing
  from the two checkers above) of the overlap identity and exact curve, the tangent and switch,
  and the controller `D_x` spectrum, `Tr(Q_x)` and the sharp-minus-`h` gap.
- **`research/check_attainment.py`** — symbolic endpoint/mixture arithmetic for sharpness.
- **`research/check_literature_comparisons.py`** — the two exact counterexamples and the
  comparison table from `PRIOR_ART.md` (BCW failed substitution; independent-observable
  steering counterexample; sharp vs BCW+Renes).

### Tier 5 — optional second route (python-flint)

- **`proofs/verify_alternate.py`** — re-checks the same certificate with a different
  factorization and exact Sylvester determinants. Requires `requirements-optional.txt`.
  Deliberately **not** required for the central result.

## What none of this establishes

- **Not human expert review.** Every check here was written and run by AI (Claude and ChatGPT
  "Astra" across sessions). Passing checks written by the same family of systems is not
  independent scientific review.
- **Not priority or novelty.** See `PRIOR_ART.md`. The searches were bounded; absence of a
  match is not clearance. Two literature reductions remain open (see `OPEN_PROBLEMS.md`).
- **Not experimental feasibility.** The theorem is an idealized existence guarantee under the
  trusted-control assumptions of `PROTOCOL.md`. Robustness to control/phase error, and an
  efficient forward-only implementation of the *sharp* recovery map, are open.
- **Not correctness of the prose steps by machine.** `f <= p`, the Choi normalization, the
  Uhlmann attainability, and the general two-branch construction are mathematical arguments to
  be read (they are elementary, but they are not script output). The certificate proper — the
  free-unitary identity and PSD witnesses — is what the exact checks cover.

Measured runtime: the required `run_checks.py` completes in a few seconds; the numeric
instantiation (Tier 3) takes roughly a minute (200 trials over several dimensions).

## Joint calibration certificate chain

`python proofs/verify_joint_recovery.py` and `python proofs/verify_joint_adversary.py`
use standard-library rational arithmetic, refuse optimized Python, and locate
immutable data under `certificates/` independently of the working directory.
Their success output follows all certificate and supplied mutation gates. Both
are included in `run_checks.py`. The adversary imports arithmetic helpers from
the lower verifier; these are not two independent verification implementations.

`python tests/verify_joint_born.py` additionally needs NumPy and SciPy. It is
floating-point corroboration and is run in a separate CI job; SciPy is not a
core dependency. Exact arithmetic does not validate the physical or statistical
prose automatically. The attached predecessor-calibration comparisons are
provenance only until their missing derivation is supplied.

## Robust forward-query checks

- `proofs/verify_robust_forward_obstruction.py`: standard-library rational
  moment and full-history dual checks, followed by five named mutations. It
  is a required check and refuses optimized Python.
- `tests/verify_robust_forward_bridge.py`: seeded NumPy/SciPy corroboration. It
  prints results without rewriting certificates or checked-in logs.
- `research/check_forward_reductions.py`: exact SymPy identities for the
  fixed-scaling completion examples; their physical and positivity
  interpretations are stated separately in the docs.
- `history/endpoint_forward_queries/verify_certificates.py`: optional historical
  standard-library chain. It requires all four named certificates, rejects
  optimized Python, and defers success output until all mutations pass.

The numerical bridge reuses coefficient builders from the exact verifier; it is
corroboration, not a wholly independent proof implementation. The physical
finite-unitary extension and controller-class arguments are not machine
formalized by any of these checks.
