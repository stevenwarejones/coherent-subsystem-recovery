# A forward-only, two-query implementable recovery guarantee

This is a **separate research addition**. It does not replace or modify the sharp certificate
(`docs/RESULTS.md`, `docs/THEOREM.md`); it addresses the implementability question raised in
`docs/OPEN_PROBLEMS.md` §3 with an explicit circuit, at a slope short of sharp.

## Result

With `f` the coherent-test score and `p` its support-projector probability (`f ≤ p`, same
bridge as the main theorem), the explicit fixed two-query recovery circuit acting on **M alone**
has the exact worst-case score-to-recovery curve
```
r_UV  >=  h(f) := max{0, 9f - 1}^2 / 64,     0 <= f <= 1,
```
and for every score `f` there is an allowed physical instance attaining equality — so `h` is
the exact performance of *this fixed circuit* (not the optimum recovery channel, not every
two-query algorithm). The earlier affine bound `(9f−5)/4` is exactly the tangent to `h`: their
difference is `81(1−f)^2/64 ≥ 0`, so `h` improves it by that much for `f ≥ 1/9`. The asymptotic
infidelity coefficient is still `9/4` (`h(f) = 1 − (9/4)(1−f) + O((1−f)^2)`).

Combined with the elementary decoder and the replacement channel, the *selectable implementable*
**guarantee function** is
```
g(f) = max{ 1/4, 3f/4, h(f) },       with  r_selected(f) >= g(f),
```
switching from the elementary `3f/4` decoder to the two-query circuit at
`f_c = (11 + 4√7)/27 = 0.799370564…` (was `5/6` with the affine bound). `g` is a guaranteed
lower bound for the selected implementation, not a proved exact minimax over the selected
family (attaining `h` alone does not establish simultaneous attainment for all three decoders).
These are squared Bell overlaps, not average input fidelities.

The optimal-recovery theorem is still stronger — its high-score line is `(3f−1)/2`. The gap
between the sharp curve and the selected guarantee `g` is piecewise
```
0                          for 0    <= f <= 2/3,
(3f - 2)/4                 for 2/3  <= f <= f_c,
(3f - 1)/2 - (9f - 1)^2/64 for f_c  <= f <= 1,
```
maximal at `(4√7 − 7)/36 = 0.099527…` (down from `1/8`), and zero at `f = 1` (e.g. at
`f = 0.99`: `0.977626…` here vs `0.985` sharp). **The sharp theorem's recovery channel has not
been given a forward-query implementation;** this is an explicit weaker-coefficient circuit,
not that.

### Proof of `h` (overlap identity + Cauchy–Schwarz)

With `B = X⊗U`, `C = Z⊗V` (`B†B = C†C = I`), `L = I + B + C`, `T = (I+B)(I+C)`, using only
unitarity,
```
L†L = I + B†T + T†B.
```
Since `p = Tr(ρ L†L)/9` and `r_UV = Tr(ρ T†T)/16`, and `Tr(ρ B†B) = 1`,
```
9p − 1 = 2 Re Tr(ρ B†T) ≤ 2|Tr(ρ B†T)| ≤ 2 √(Tr(ρ B†B) Tr(ρ T†T)) = 8 √r_UV
```
by Cauchy–Schwarz for `B√ρ, T√ρ` in Hilbert–Schmidt norm (valid for mixed, rank-deficient ρ).
For `p ≥ 1/9`, squaring gives `r_UV ≥ (9p−1)^2/64`; for `p < 1/9`, `r_UV ≥ 0`. Since `h` is
nondecreasing and `f ≤ p`, `r_UV ≥ h(p) ≥ h(f)`. Attainment: `ρ = Φ⁺` on a qubit M with
`U = e^{it}X, V = e^{it}Z` gives `f = p = (5+4cos t)/9` and `r_UV = (1+cos t)^2/4 = (9f−1)^2/64`
on `[1/9,1]`; a CPTP mixture of the two original test decoders (identity with probability `9f`,
and `X`) completes `[0,1/9]` with the pre-probe state and evaluated circuit unchanged. No new
SDP certificate is used. The essential caveat stands: on these examples the optimal `M`-only
recovery is `R = 1`, so `h` is a limitation of *this circuit*, not of recovering the subsystem.

### This is minimax-optimal for the fixed query stage plus controller-only decoding

Replacing the final Bell-basis change and discard with an **arbitrary controller-only CPTP
channel** (an arbitrary map from the 4-dimensional controller to the recovered qubit, allowed
to depend on the scalar score `f` but not on hidden phase information) cannot do better than `h`
on `[7/9,1]`:
```
sup_{Λ_f}  inf_{allowed instances with score f}  r(Λ_f)  =  h(f),     7/9 ≤ f ≤ 1.
```
Two exact facts give this. **Endpoint uniqueness:** at the ideal Bell instance the post-query
controller-reference state is `Φ_{Aa} ⊗ I_c/2`, whose full-rank gauge forces any final channel
achieving `r = 1` to be uniquely the partial trace `Tr_c` (a standard noiseless-subsystem
rigidity; `D_0 = ¼(I−Φ)_{oa}⊗I_c ⪰ 0` and `Tr(D_0 J)=0` pins the Choi matrix `J = 2Φ_{oa}⊗I_c`).
**Two-instance obstruction:** for `U_± = e^{±it}X, V_± = e^{±it}Z` (`x = cos t ∈ [1/2,1]`, both
score `f = (5+4x)/9 ∈ [7/9,1]`), an exact `8×8` PSD dual slack `D_x = I_o⊗Q_x − ρ̄ᵀ/2` with
`Q_x = (1+x)^2/16·I + (1−x^2)/16·Y_aY_c` and spectrum `{0, (x+1)/8, x/4, (x+1)(2x−1)/8}` (each
doubled) gives, for any CPTP Choi `J` (`Tr_o J = I`),
```
½(r_+ + r_-) = ½ Tr(ρ̄ᵀ J) ≤ Tr[(I_o⊗Q_x) J] = Tr(Q_x) = (1+x)^2/4 = h(f),
```
so at least one of the two instances scores `≤ h(f)`, matched by the existing decoder. Since
`(3f−1)/2 − h(f) = 3(1−f)(27f−11)/64 > 0` on `[7/9,1)`, **this architecture cannot reach the
sharp high-score bound even with score-dependent decoding** — a finite analytic obstruction, not
a failed search. It is *not* a lower bound for all two-query circuits: it fixes the query stage
and decodes only from the controller. (Below `x = 1/2` the certificate stops; no claim is made
below `f = 7/9`.)

Two restricted dead ends are also closed exactly: reversing the query order (fourth Pauli `ZX`,
`T_rev = (I+C)(I+B)`) and any classical randomization of the two orders have the same curve `h`
(the phase family saturates both at equal score); and merely reweighting the four branch
amplitudes cannot satisfy the `f = 1` requirement except as the original circuit up to global
phase. Neither excludes coherent order control with a different joint decoder, an intermediate
controller operation between the two calls, adaptive operations, or recovery that retains access
to `M` — those need a new argument (see `OPEN_PROBLEMS.md`).

## Circuit

Prepare the **two controller qubits** in `|++>`; apply controlled `V` on one and controlled `U`
on the other, giving branches `I, U, V, UV`. A fixed basis change acts on **those same two
qubits**, sending the four computational labels to the orthonormal Bell-type states
`|β_j> = (I⊗P_j)|Φ+>` with `P_j = (I, X, Z, XZ)`. Retain **one** of the two qubits as the
recovered output; discard the **other** qubit and M. There is no separate leftover controller
register — the controller qubits *are* the Bell-labelled ancillas — so one must not copy the
labels into a fresh system and then trace the original controller, which would destroy the
coherence. One controlled-`U` call, one controlled-`V` call, two fresh qubits, fixed gates,
discarding. **No inverse queries, no knowledge of `U,V`, no access to the bypass B, no
postselection.** The circuit description and query count are independent of `dim(M)`.

The Bell acceptance effect on AM is `T†T/16` with
```
T = I + X⊗U + Z⊗V + XZ⊗UV = (I + X⊗U)(I + Z⊗V).
```
The fourth (XZ) branch is antisymmetric; using `P_j^T` instead of `P_j` flips its sign and is
wrong — the instance checker caught exactly this in the original implementation.

The `T†T/16` effect is established by a general index calculation, valid in every dimension:
with `β_j(a,c) = (P_j)_{ca}/√2`, contracting the reference and the retained output against the
Bell bra gives environment amplitude `(1/4) Σ_j (P_j)_{ca} W_j`, whose squared norm summed over
the discarded indices is `T†T/16`. The three exact matrix instances in the checker corroborate
this convention (including the `Y`/XZ sign on complex nonsymmetric unitaries); they are not
themselves the arbitrary-dimensional proof.

## The affine tangent (one-square identity)

The affine bound `(9f−5)/4`, the tangent to `h` at `f = 1`, also follows from a one-square
identity, kept here as an independent (weaker) derivation. With `J = 3B − I − C − BC`,
```
T†T + 20 I − 4 L†L = J†J,     so     r_UV = (9p − 5)/4 + Tr(ρ J†J)/16 ≥ (9p − 5)/4.
```
Only unitarity of `B, C` is used. The overlap proof in **Result** above is stronger — it gives
the exact quadratic curve `h`, of which this line is the tangent.

## The exact curve is attained only by these instances; it bounds THIS circuit only

The phase family `ρ = Φ⁺`, `U = e^{it}X`, `V = e^{it}Z` attains `r_UV = (1+cos t)²/4 = (9f−1)²/64`,
so `h` is the exact pointwise curve of this circuit, not merely a lower bound. These instances
have optimal abstract recovery `R = 1` throughout; `h` bounds *this circuit* (and the
controller-only architecture around it), not other circuits or the subsystem's recoverability.

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

- `proofs/verify_forward_recovery.py` — SymPy; the one-square identity (affine tangent), isometry
  and Bell-effect conventions on three exact matrix instances (incl. nonsymmetric complex
  unitaries), the phase-family ratio, and the old switch threshold. Refuses `python -O`.
- `proofs/verify_forward_curve.py` — SymPy + NumPy; the exact overlap identity `L†L=I+B†T+T†B`,
  the curve `h`, the tangent/improvement `81(1−f)²/64`, the switch `f_c` and maximum gap, plus
  500 random matrix/circuit cases (dims 1,2,3,4,6), 101 phase instances, reversed order and the
  low-score endpoint, with a changed-constant rejection control.
- `proofs/verify_controller_decoder.py` — SymPy; the ideal controller state, the endpoint-unique
  Choi matrix, the exact `D_x` spectrum and PSD interval, matching primal/dual and score
  conversion, five rational-phase reconstructions, and two rejection controls. Runs under `-O`
  (explicit exceptions).
- `tests/verify_forward_curve_independent.py` — an independent from-scratch re-derivation
  (imports nothing from the above): the overlap identity and curve over random `U,V` in
  dimensions 1–5, the tangent and switch, and the controller `D_x` spectrum, `Tr(Q_x)` and gap.
- `tests/verify_forward_independent.py` — the earlier independent check of the one-square
  identity, `r ≥ (9p−5)/4`, the phase family, and `forward ≤ sharp`.
- `research/explore_controller_decoder.py` — optional numerical discovery (NumPy + CVXPY,
  `research/requirements.txt`); compares endpoint-pinned and unrestricted final decoders. Not a
  certificate and not part of verification.
