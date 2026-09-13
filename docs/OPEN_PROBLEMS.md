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
exactly `h(f)` on `[7/9,1]` (an `8×8` PSD dual plus endpoint uniqueness), so **that architecture
cannot reach `(3f−1)/2`** — the gap is `3(1−f)(27f−11)/64 > 0`. Reversing or classically
randomizing the query order, and merely reweighting the four branch amplitudes, are also closed.

**Still open.** A forward implementation of the *sharp* `(3f−1)/2` guarantee, if one exists, must
change the query stage itself — the next model is an intermediate controller operation *between*
the two calls (the present obstruction fixes the query stage and decodes only from the
controller, so it does not exclude this). The exclusion is also silent on coherent order control
with a different joint decoder, additional queries, and recovery operations that retain access to
`M`; a lower bound over those needs a new argument, and unitary-inversion query lower bounds do
not transfer without a reduction. Everything here carries the **extra access assumption** —
reusable controlled forward access to `U,V`, stronger than one-shot observation.

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
