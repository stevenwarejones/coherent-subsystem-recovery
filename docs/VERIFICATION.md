# Verification: what each check does, and what it does not establish

## One command

```sh
python run_checks.py          # standard library + SymPy; no solver, no network
```

This runs the four required checks (exact certificate, elementary full-curve segments, the
forward two-query guarantee, and the mutation cases) fail-fast, streaming output. It exits
nonzero on any failure and prints its success sentence only after all four pass.

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
  free-unitary identity `T†T+20I−4L†L=J†J`, the circuit's isometry/Bell-effect conventions on
  exact instances (incl. complex nonsymmetric unitaries), the phase-family slope, and the
  selection threshold. See `FORWARD_RECOVERY.md`. Refuses `-O`.
- **`tests/verify_forward_independent.py`** — an independent from-scratch NumPy re-derivation of that identity
  (random `U,V`, dims 1–5), `r_UV ≥ (9p−5)/4`, the phase family, and `forward ≤ sharp`.
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
