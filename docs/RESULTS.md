# Results and claim ledger

All symbols are defined in `PROTOCOL.md`. The canonical proof write-up is
`THEOREM.md`; this file states the theorems, their dependencies, and a ledger
separating what is exact-checked from what is prose, numerical, conjectural, or
a literature report.

## Main theorem (sharp high-score bound)

For arbitrary finite-dimensional `U, V` on M and arbitrary finite-dimensional bypass B,
```
R >= (3f - 1) / 2.
```
For every `f in [2/3, 1]` an allowed construction attains equality, so the minimum
recoverability compatible with an exact score `f` on that interval is exactly `(3f-1)/2`.
Because the full curve below is nondecreasing, the same minimum holds under the constraint
*score at least f*. The attaining construction uses addressed dimension at most six.

The proof is an exact free-unitary positive-operator certificate with normalization `t = 1`
(no numerical tolerance, no limiting argument): two positive-semidefinite Gram witnesses over
53 reduced words of length ≤ 3 in `U,U†,V,V†`, satisfying, for **every** pair of unitaries in
**every** finite dimension,
```
Q - S = H,     Tr_A Q = I_M,     Q, S >= 0,
```
where `H = (1/3)[X⊗(U+U†) + Z⊗(V+V†) + XZ⊗(U†V − V†U)]`. `Q` is the Choi matrix of an explicit
CPTP recovery channel on M, giving `R >= (1/2) Tr(Q rho) >= (1/2) Tr(H rho) = (3p-1)/2 >= (3f-1)/2`.

## Full curve (all scores)

```
R_min(f) = max{ 1/4,  3f/4,  (3f-1)/2 },   0 <= f <= 1.
```
The three segments join continuously at `f = 1/3` (both `1/4`) and `f = 2/3` (both `1/2`).
- `R >= 1/4`: discard M, prepare a fixed state; `rho_A = I/2`.
- `R >= 3f/4`: elementary certificate `Q_0 = (3/2)K = L†L/6`, using `Tr_A K = (2/3) I_M`.
- `R >= (3f-1)/2`: the main theorem.
The **only** deep segment is the high-score one; the other two are elementary. Matching
constructions attain every segment (addressed dimension ≤ 6 overall).

## Mechanism (why three branches)

- **Two** fixed-unitary branches can be perfectly simulated (score 1) using only a classical
  addressed record, with `R = 1/2` — no quantum recoverability certified. (General for any two
  target unitaries, via diagonalizing `P_0† P_1`.)
- **Three** branches certify `R > 1/2` once `f > 2/3`.
- **Four** complete Pauli branches admit a simple universal *forward-only* recovery circuit
  with `R >= f`.

An explicit **two-query forward circuit** (one controlled-`U`, one controlled-`V`, no inverses)
gives an *implementable* guarantee `r_UV >= (9f-5)/4` under the extra assumption of reusable
controlled forward access; see `docs/FORWARD_RECOVERY.md`. This is short of the existential sharp
map (`3/2` vs `9/4` slope) and is kept as a separate research addition.

## Claim ledger

| # | Claim | Assumptions | Proof location | Executable support | Independent-review status |
|---|---|---|---|---|---|
| 1 | `R >= (3f-1)/2`, all finite dim | PROTOCOL model | THEOREM.md §3–4 | `proofs/verify_sharp_recovery.py` (exact, stdlib); `proofs/verify_alternate.py` (python-flint) | exact certificate re-derived from scratch in `tests/verify_independent.py`; **no human expert review** |
| 2 | sharpness on `[2/3,1]` | PROTOCOL model | THEOREM.md §5 | `research/check_attainment.py`; `tests/verify_attainment_independent.py` | endpoints + mixture reconstructed independently; the Uhlmann/Bell-overlap steps are prose |
| 3 | `f <= p` (support projector) | PROTOCOL model | THEOREM.md §2 | — | **prose lemma** (elementary; not machine output) |
| 4 | Choi channel gives `R >= (1/2)Tr(Q rho)` | PROTOCOL model | THEOREM.md §4 | `tests/verify_numeric_instantiation.py` (concrete `U,V`) | numeric corroboration only; the identity is the proof |
| 5 | `R >= 3f/4` and `R >= 1/4` | PROTOCOL model | full-curve handoff §2–3 | `proofs/verify_full_curve.py` | elementary; **new segments not yet expert-reviewed** |
| 6 | full curve `R_min(f)` | PROTOCOL model | full-curve handoff | `proofs/verify_full_curve.py` | as above |
| 7 | two-branch obstruction | protocol family | full-curve handoff §5 | `proofs/verify_full_curve.py` (Pauli instance) | general spectral construction is **prose** |
| 8 | four-Pauli `R >= f` | ideal Paulis | full-curve handoff §6 | `proofs/verify_full_curve.py` | comparison/benchmark, **not** a novelty claim |
| 9 | `f <= P_pg` is **false** (BCW is a different score) | — | PRIOR_ART.md | `research/check_literature_comparisons.py` | exact counterexample, verified |
| 10 | sharp bound strictly beats BCW+Renes on `[2/3,1)` | — | PRIOR_ART.md | `research/check_literature_comparisons.py` | exact, verified |
| 11 | forward two-query circuit: `r_UV >= (9f-5)/4`, via `T†T+20I-4L†L=J†J` | PROTOCOL + **reusable controlled forward access** to `U,V` | FORWARD_RECOVERY.md | `proofs/verify_forward_recovery.py`; `tests/verify_forward_independent.py` | exact identity re-derived independently; **extra access assumption**; not the sharp map (slope 9/4 vs 3/2) |
| 12 | selectable implementable guarantee `max{1/4,3f/4,(9f-5)/4}`, switch at `f=5/6` | as row 11 | FORWARD_RECOVERY.md | `proofs/verify_forward_recovery.py` | exact; strictly below the existential sharp curve by `3(1-f)/4` |

Legend of evidence types: **exact** = exact-arithmetic machine check; **prose** = a
mathematical argument checked by reading, not by a script; **numeric** = floating-point
corroboration that would catch a convention error but is not itself a proof; **conjecture**
/ **literature** where stated. No row has been reviewed by a human domain expert.

## Fidelity conventions (do not conflate)

`R` is squared entanglement fidelity at the maximally mixed input; average pure-input qubit
fidelity is `(2R+1)/3`. The appearance of `2/3` in a classical average-fidelity benchmark is
**not** the same object as the coherent-score threshold `f = 2/3`. See
[Nielsen, quant-ph/0205035](https://arxiv.org/abs/quant-ph/0205035).
