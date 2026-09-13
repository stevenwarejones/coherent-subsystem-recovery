# A forward-only, two-query implementable recovery guarantee

This is a **separate research addition**. It does not replace or modify the sharp certificate
(`docs/RESULTS.md`, `docs/THEOREM.md`); it addresses the implementability question raised in
`docs/OPEN_PROBLEMS.md` §3 with an explicit circuit, at a slope short of sharp.

## Result

With `f` the coherent-test score and `p` its support-projector probability (`f ≤ p`, same
bridge as the main theorem), an explicit fixed recovery circuit acting on **M alone** attains
```
r_UV  >=  (9p - 5)/4  >=  (9f - 5)/4,
```
equivalently `r_UV ≥ 1 − (9/4)(1 − p)`. Combined with the elementary decoder and the
replacement channel, the *selectable implementable* guarantee is
```
r_selected(f) = max{ 1/4, 3f/4, (9f - 5)/4 },
```
switching from the elementary `3f/4` decoder to this one at `f = 5/6`. These are squared Bell
overlaps, not average input fidelities.

The optimal-recovery theorem is still stronger — its high-score line is `(3f−1)/2`, and the gap
to this implementable circuit is `3(1−f)/4` (e.g. at `f = 0.99`: `0.9775` here vs `0.985`
existential; both reach `1` at `f = 1`). **The sharp theorem's recovery channel has not been
given a forward-query implementation;** this is an explicit lower-slope circuit, not that.

## Circuit

Two controller bits in `|++>`; apply controlled `V` then controlled `U` on the addressed
subsystem, giving branches `I, U, V, UV`; map the four labels by a fixed unitary to the
orthonormal `|β_j> = (I⊗P_j)|Φ+>` with `P_j = (I, X, Z, XZ)`; keep one ancilla as the recovered
qubit and discard the other, the controller, and M. One controlled-`U` call, one controlled-`V`
call, two fresh ancillas, fixed gates, discarding. **No inverse queries, no knowledge of `U,V`,
no access to the bypass B, no postselection.** The circuit description and query count are
independent of `dim(M)`.

The Bell acceptance effect on AM is `T†T/16` with
```
T = I + X⊗U + Z⊗V + XZ⊗UV = (I + X⊗U)(I + Z⊗V).
```
The fourth (XZ) branch is antisymmetric; using `P_j^T` instead of `P_j` flips its sign and is
wrong — the instance checker caught exactly this in the original implementation.

## Proof (one square)

With `B = X⊗U`, `C = Z⊗V` (`B†B = C†C = I`), `L = I + B + C`, `T = (I+B)(I+C)`, and
`J = 3B − I − C − BC`, there is the exact free-unitary identity
```
T†T + 20 I − 4 L†L = J†J.
```
Since `p = Tr(ρ L†L)/9` and `r_UV = Tr(ρ T†T)/16`,
```
r_UV = (9p − 5)/4 + Tr(ρ J†J)/16 ≥ (9p − 5)/4.
```
Only unitarity of `B, C` is used — no commutation, no Hermiticity, no Gram search, no numerics.

## The 9/4 slope is sharp for THIS circuit only

For `ρ = Φ⁺`, `U = e^{it}X`, `V = e^{it}Z` (original target branches fixed):
`p = f_opt = (5 + 4cos t)/9`, `r_UV = (1 + cos t)²/4`, and `(1−r_UV)/(1−f_opt) → 9/4` as `t → 0`.
So no bound `r_UV ≥ 1 − c(1−f)` with `c < 9/4` holds universally for this circuit. This example
has optimal abstract recovery `R = 1` throughout; it bounds *this circuit*, not other circuits,
and is not a proof that the affine bound is this circuit's exact pointwise curve.

## Additional access assumption (do not hide this)

Running the circuit needs **reusable controlled forward access** to `U, V` — the ability to call
them again, coherently, on the same addressed Hilbert space in a fresh recovery run (stationary
implementations, same relative phases, no uncontrolled memory between queries). This is stronger
than merely observing a one-shot intervention: the recovery acts on the pre-probe state (e.g. on
fresh identically prepared trials), and it is **not** a claim that measuring the original test
leaves that state available. The reference and the original bypass remain inaccessible; fresh
ancillas are available.

## Prior art (screened, AI-inspected; not novelty clearance)

Algorithmic recovery is established: Gilyén et al. QSVT Petz recovery
([2006.16924](https://arxiv.org/abs/2006.16924)); simulation of adjoints/Petz maps for unknown
channels ([2602.05828](https://arxiv.org/abs/2602.05828)); query lower bounds for transforming
unknown unitaries ([2405.07625](https://arxiv.org/abs/2405.07625) — those bounds concern inverting
an arbitrary unitary on arbitrary inputs and do not automatically constrain outputting a qubit
correlated with a protected reference; no reduction proved). Coherent teleportation and
programmable processors supply the circuit ingredients. The candidate specific additions are the
score-to-recovery translation and the one-square identity, subject to a targeted comparison with
extraction-channel bounds. **No external or human review; no priority claim.**

## Verification

- `proofs/verify_forward_recovery.py` — SymPy; free-unitary expansion of the identity, isometry
  and Bell-effect conventions on three exact matrix instances (incl. nonsymmetric complex
  unitaries), the phase-family ratio, and the switch threshold. Refuses `python -O`.
- `tests/verify_forward_independent.py` — my from-scratch NumPy re-derivation: the identity over
  random `U,V` in dimensions 1–5, `r ≥ (9p−5)/4`, the phase family, and `forward ≤ sharp`.
