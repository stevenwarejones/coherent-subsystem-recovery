# Open problems

The exact theorem is proved (subject to human review). These are genuinely open, and the
first two are the **gate on any strong novelty claim in a paper** — not merely "future work."

## 1. The two literature reductions (novelty gate)

Neither is closed. Each has been converted from an open-ended worry into a *specific isolated
obstruction*, which is the honest status:

- **Barizien identity-comb (supermap self-testing).** The exact identity-comb result assumes an
  exact perfect-transmission channel identity plus a paired test; a single scalar coherent
  score at `f < 1` supplies neither. The literal exact result does not reduce to our curve. A
  robust/quantitative version of its proof *might* imply our bound — that is a defined proof
  task, not a missing citation. Our theorem is a recovery statement, not a full supermap
  self-test.
- **Farenick–Kavruk–Paulsen Ando/WEP completion.** The cited transfer theorems move an *already
  feasible* positive completion into an algebra under WEP; they do not establish feasibility
  from arbitrary `U,V`. There is also an exact obstruction to the simplest off-diagonal
  substitution: at `U=X, V=Z`, the block `H_01` has numerical radius `2/3 > 1/2`, so completing
  `Q` with off-diagonal fixed to `H_01` fails even at the perfect point — the certificate
  succeeds precisely by letting `Q`'s off-diagonal differ while keeping `Q − H ≥ 0`. Status:
  **no reduction obtained; the required extra inequality has been isolated** — not "Ando
  excluded" and not "novelty proved."

**Reviewer ask.** Either produce a complete derivation of the same score-to-M-recovery curve
from an existing recovery/supermap/operator-completion theorem, with every access assumption
and normalization matched; or a precise obstruction naming the missing hypothesis. White et
al.'s process-completion approach is a third concrete route worth testing.

## 2. Human expert review

No step here has been reviewed by a human quantum-information specialist. The specific asks:
reconstruct `f ≤ p`, the block identity and Choi orientation/normalization, the Uhlmann
attainability, and the general two-branch spectral construction; and verify the free-unitary
certificate separately from its physical interpretation.

## 3. Efficient forward-only implementation of the sharp map — partially addressed

The certificate specifies the recovery channel algebraically, but it uses adjoint words.

**Resolved for the fixed two-query circuit and the controller-only architecture (see
`docs/FORWARD_RECOVERY.md`).** The explicit two-query forward circuit — one controlled-`U`, one
controlled-`V`, two ancillas, no inverses, no knowledge of `U,V`, no bypass access — has the
*exact* worst-case curve `r_UV ≥ h(f) = max{0,9f−1}^2/64`, with equality attained at every score
(the old affine `(9f−5)/4` is its tangent). Moreover, replacing the final Bell-basis decode with
an *arbitrary controller-only CPTP channel*, even one chosen from the scalar score, gives minimax
exactly `h(f)` on `[7/9,1]` (an `8×8` PSD dual plus endpoint uniqueness; full proof in
`docs/CONTROLLER_DECODER.md`), so **that architecture cannot reach `(3f−1)/2`** — the gap is
`3(1−f)(27f−11)/64 > 0` on `[7/9,1)` (it vanishes at `f = 1`). Reversing or classically randomizing
the query order is also closed; and changing only the four **fixed preparation amplitudes**, with
the specified query words and final Bell decoder retained, cannot meet the `f = 1` requirement
except as the original circuit up to global phase (this is not a statement about arbitrary
preparation plus arbitrary decoding).

**Larger class now addressed, subject to the stated prose review.** The
[robust full-history certificate](ROBUST_FORWARD_OBSTRUCTION.md) excludes
universal sharp recovery with at most three active calls under fixed or
preselected random order, even with arbitrary intermediate controller processing
and score-selected coefficients. The VUVUV slot result is not an exclusion of
all five-call patterns. At least four active calls are necessary in this class;
four have not been shown sufficient.

**Still open.** The [finite-query model](FORWARD_QUERY_MODEL.md) defines the
universal score-independent search and its proposed support-specific complete
moment formulation. Its compactness argument does not cover score-selected
controller families. Neither a universally sharp finite-forward-query recovery
nor an obstruction for every finite count is established. Outcome-dependent
query order, indefinite causal order, extra target operations, and postselection
are outside the stated obstruction. Controlled oracle reuse remains an extra
access assumption beyond the original coherent test.

The [restricted phase family](FORWARD_PHASE_FAMILY.md) has a randomized
three-call recovery meeting the target, so that family cannot establish the
general negative result.

## 4. Robustness to trusted-control error

The theorem assumes exact coherent implementation of `I, U, V` including relative phases. A
defensible corollary: with a calibration budget `η` such that the ideal-model score is at
least `f_obs − η`, `R ≥ (3(f_obs − η) − 1)/2` — immediate substitution, not a new theorem. The
substantive open task is to bound `η` from a physical model of control error / branch leakage /
target-measurement error. An unrestricted leakage model is not covered. Statistical confidence
on the score (a lower confidence limit `f_L`) is a separate, stateable substitution.

## 5. Scope extensions

Infinite-dimensional domains are not claimed. Multi-time / sequential versions, and other
target sets beyond `(I,X,Z)`, are unexplored here.

## Joint calibration for a different two-query recovery

The [balanced recovery](JOINT_CALIBRATION_RECOVERY.md) meets the sharp target
under `q <= (1-f)/12`; a counterexample brackets its universal threshold by
`1/12 <= c_bal <= 522171/4500256`. This does not solve score-only forward
recovery or contradict the existing fixed-query-stage minimax theorem: the
circuit and observed data differ. A useful next question is whether another
accessible controller can handle the adversary without this small-q condition.

## Compilation prerequisite for control-error work

The supplied control-error handoff assumes a weighted LCU/adjoint-access
compiler absent from this checkout. Its full ideal construction and resource
accounting should be integrated and reviewed before promoting the conditional
27 eta_max error budget. Passing its qubit simulation is not a generic compiler
proof. No such hardware-noise guarantee is added in this change.

## Fixed off-diagonal completion scaling

[Two exact examples](COMPLETION_OBSTRUCTION.md) exclude every universal constant
prescription `Q_01 = k H_01`. Coupled, nonlinear, or larger-block completion
reductions remain open; the historical novelty gate remains unresolved.
