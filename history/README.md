# History / discovery material (not the current result)

This directory holds superseded artifacts, kept for provenance. **Do not cite these as the
current result**, and do not let their language compete with the sharp theorem.

## `three_branch_1041_certificate/`

The earlier exact certificate proving `R >= (500/1041)(3f-1)` — normalization `t = 1041/1000`,
17 reduced words of length ≤ 2. It was a valid dimension-independent bound, but it is
**superseded** by the sharp `t = 1` certificate in `../certificates/sharp_recovery.json`
(proved bound `R >= (3f-1)/2`, which meets the attaining construction). The `1041/1000`
certificate's only advantage was a mid-fidelity segment that the sharp result now dominates
everywhere; at `f = 1` it reached only `1000/1041 ≈ 0.9606` and never `1`.

Its files still verify on their own (standard library):
```sh
cd history/three_branch_1041_certificate
python3 verify.py
python3 test_mutations.py
```
They are retained so the discovery path is auditable, not as a current claim.
