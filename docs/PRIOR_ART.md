# Prior-art status

**Nothing here has been reviewed by a human domain expert. No priority is claimed.**

Three literature sweeps were performed (all AI-conducted, bounded, English-language,
primary-text where noted). They did **not** locate a published derivation of the exact
curve `R_min(f) = max{1/4, 3f/4, (3f-1)/2}` for this access model, but absence of a match in
a bounded search is not novelty clearance. The defensible contribution is narrow:

> a sharp, dimension-independent recovery guarantee for the addressed subsystem M alone,
> inferred from the specified three-branch coherent test, even though the final decoder may
> use an arbitrary finite-dimensional quantum bypass.

Everything grander is **ceded to prior art**: recoverability certification, robust
quantum-memory certification, and internal-slot/supermap certification all already exist.

## The closest formula: Berta–Coles–Wehner (2014)

BCW's Eq. (10), specialized to a reference qubit, gives `R_pg = (3 P_pg − 1)/2` — the **same
affine expression** as ours, for *pretty-good guessing and recovery*. So the formula alone
cannot support novelty. But it is a **different score**, and the naive identification fails:

- Exact counterexample (verified in `research/check_literature_comparisons.py`):
  `rho = (3/4)Φ⁺ + (1/4)I₄/4`, `U=X`, `V=Z` gives `p = f_opt = 5/6` while `P_pg = 25/32 < 5/6`.
  So `f ≤ P_pg` is **false**; the coherent score is not the pretty-good guessing probability.
- Stated precisely: this rules out the *direct substitution*, not every conceivable reduction
  from the BCW framework.

A valid, stronger literature-derived comparison does exist (BCW 2014 + Renes 2017, Eq. (12)):
with `s = (1+3f)/4`, `R ≥ 3s² − 3s + 1`. This is legitimate independent corroboration that
recovery → 1. The sharp bound strictly beats it on `[2/3, 1)`:
```
f:      2/3     0.80    0.90    0.95    0.99    1.00
sharp:  0.5000  0.7000  0.8500  0.9250  0.9850  1.0000
BCW+R:  0.4375  0.6175  0.7919  0.8917  0.9777  1.0000
```
At `f = 1−ε`: sharp is `1 − (3/2)ε`; BCW+Renes is `1 − (9/4)ε + (27/16)ε²`. They share the
endpoint and linear scaling but have **different first-order coefficients** (3/2 vs 9/4), so
the literature bound corroborates recoverability but does **not** independently pin the 3/2
slope. (An earlier summary that called these "leading-order agreement" was imprecise; the
correct statement is different first-order coefficients.)

## Coherent control makes relative phases physical

Characterizing only the isolated channels `U(.)U†` and `V(.)V†` is insufficient — relative
phases of the controlled implementations are visible in the test. Exact illustration (verified):
physical branches `(I,X,Z)` vs `(I,−X,Z)` induce identical individual channels, yet give
`C = I/2` vs `I/18`, hence `f_opt = 1` vs `1/9`. This is an assumption to expose, not hide
(Abbott et al., [1810.09826](https://arxiv.org/abs/1810.09826)).

## Comparator map (all AI-inspected; "compared" = specific eqs read in full text)

| Source | Role | Bottom line |
|---|---|---|
| [Berta–Coles–Wehner 2014, 1302.5902](https://arxiv.org/abs/1302.5902) | same affine formula, different score | direct substitution fails (exact c/e); closest comparator |
| [Renes 2016, 1605.01420](https://arxiv.org/abs/1605.01420) / [Renes 2017, 1707.01114](https://arxiv.org/abs/1707.01114) | translated recovery bounds | give the valid BCW+Renes benchmark; sharp bound still stronger |
| [Wang et al. 2022, 2210.11243](https://arxiv.org/abs/2210.11243) | 3-setting steering self-test | formal translation only; arbitrary-dimension 3-observable extension **not** established — keep as caveated, not a certified competitor |
| [Rosset–Buscemi–Liang 2018, 1710.04710](https://arxiv.org/abs/1710.04710) | faithful memory verification | rules out "first memory verification"; different access model (no bypass) |
| [Sekatski et al. 2023, 2304.10408](https://arxiv.org/abs/2304.10408) | tight robust memory certification | rules out "first tight robust memory certification"; DI, different object |
| [Barizien et al. 2026, 2606.25124](https://arxiv.org/abs/2606.25124) | self-testing supermaps / internal slots | closest in motivation; exact identity-comb needs a premise we lack — see OPEN_PROBLEMS |
| [Farenick–Kavruk–Paulsen, 1107.0418](https://arxiv.org/abs/1107.0418) | Ando/WEP completions | reduction not obtained; specific obstruction isolated — see OPEN_PROBLEMS |
| [White et al. 2025, 2107.13934](https://arxiv.org/abs/2107.13934) | unitary-probe process characterization | unitary probes revealing internal quantum structure is **not** new; different accessible data and target |
| [Ziman–Bužek, quant-ph/0612218](https://arxiv.org/abs/quant-ph/0612218); [Harrow, quant-ph/0307091](https://arxiv.org/abs/quant-ph/0307091); [Nielsen–Chuang, quant-ph/9703032](https://arxiv.org/abs/quant-ph/9703032) | programmable processors / coherent teleportation / no-programming | the forward-only decoder uses established ingredients — present it as an application, not a new architecture |
| [Santos et al. 2026, 2601.14191](https://arxiv.org/abs/2601.14191) | two-point-measurement DI memory | temporal causal score ≠ coherent branch score; "two-point" ≠ our "two-branch" |
| [König–Renner–Schaffner, 0807.1338](https://arxiv.org/abs/0807.1338); [Bény–Oreshkov, 0907.5391](https://arxiv.org/abs/0907.5391); [Barnum–Knill, quant-ph/0004088](https://arxiv.org/abs/quant-ph/0004088) | recovery SDP / duality / near-optimal recovery | background; the recovery optimization and its Choi/SDP form are standard |

"Compared" means the indicated definitions/equations were read in full text; it does not mean
every proof in that paper was independently verified. This was not an exhaustive
citation-graph or subscription-database review, and none of it is human expert assessment.
