# Contributing and claim discipline

This is an AI-produced research artifact under adversarial cross-review, not a
peer-reviewed result. Contributions are welcome; the non-negotiable rule is
**claim discipline**.

## Claim discipline

- Label every statement by evidence type: **exact** (machine-checked exact arithmetic),
  **proved prose** (a read argument), **numerical** (floating-point corroboration),
  **conjecture**, or **literature report**. `docs/RESULTS.md` has the ledger; keep it current.
- Do not describe AI review as peer review. Passing checks written by the same family of
  systems is not independent validation.
- Do not convert a bounded literature search into a priority or novelty claim. `docs/PRIOR_ART.md`
  states what is ceded to prior art; the two open reductions in `docs/OPEN_PROBLEMS.md` gate any
  strong novelty statement in a paper.
- No "first", "unprecedented", "optimal" (beyond the exact sharpness statement that is proved),
  "observer", or "device-independent" language. The theorem is protocol-specific; keep the
  access model attached to every claim.
- The elementary low/middle-score segments are not a second breakthrough; the deep result is the
  high-score certificate. Do not market them as co-equal.

## Verification tiers

`docs/VERIFICATION.md` defines five tiers. The central theorem is Tier 1 (exact, standard
library only). Never make a solver, `python-flint`, or any heavy dependency mandatory for the
core check. New checks must (a) run from a clean checkout via `Path(__file__).resolve()`, (b)
refuse `python -O` if they use assertions as proof gates, and (c) print no success sentence
before all their assertions have passed.

## Certificates

Preserve certificate bytes unless a reviewed change is mathematically necessary. If you add a
mutation test, it must reject for an **intended, named reason** and must include a control that
the unmodified certificate passes in the identical staging layout.

## What not to add

No general Python-source security linter, no opaque wrappers, no sprawling test framework, no
proliferation of dated top-level review files. Keep the root small; put exploration in
`research/` and superseded artifacts in `history/`.
