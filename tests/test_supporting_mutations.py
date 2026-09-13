"""Anti-degradation tests for the supporting checks (addresses review R2).

The concern: a supporting script that merely PRINTS a comparison, rather than asserting it,
lets CI treat a demonstration as verification. These tests take temporary copies of the
supporting scripts, tamper with a load-bearing value, and require that the tampered copy
FAILS (nonzero exit) and prints no success line. They leave the repository checkout
unchanged. Two of them reproduce exactly the false acceptances noted in the review.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def must_fail(script_rel, replacements, success_marker):
    src = (ROOT / script_rel).read_text()
    mutated = src
    for old, new in replacements:
        assert old in mutated, (script_rel, "anchor not found:", old)
        mutated = mutated.replace(old, new)
    assert mutated != src, (script_rel, "no change applied")
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / Path(script_rel).name
        p.write_text(mutated)
        r = subprocess.run([sys.executable, str(p)], capture_output=True, text=True, cwd=str(ROOT))
    assert r.returncode != 0, (script_rel, "tampered copy was accepted", r.stdout[-800:])
    assert success_marker not in r.stdout, (script_rel, "printed success despite tampering")
    print(f"PASS mutation: {script_rel} rejects the tampered value (exit {r.returncode}, no success line)")


# R2 reproduction 1: forcing P_pg = 1 must no longer be accepted (was: printed the FALSE line).
must_fail(
    "research/check_literature_comparisons.py",
    [("Ppg = sp.nsimplify((2 * Rpg + 1) / 3)", "Ppg = sp.Integer(1)")],
    "ALL literature-comparison assertions passed",
)

# R2 reproduction 2: zeroing the phase-family predictions must no longer be accepted
# (was: printed disagreements yet still declared PASS).
must_fail(
    "tests/verify_forward_independent.py",
    [("p_pred=(5+4*math.cos(t))/9; r_pred=(1+math.cos(t))**2/4",
      "p_pred=0.0; r_pred=0.0")],
    "PASS independent:",
)

# Extra: corrupting the constructed steering correlation must fail the exact assertion.
must_fail(
    "research/check_literature_comparisons.py",
    [("signs = {\"X\": [1, 1, -1, -1]", "signs = {\"X\": [1, 1, 1, -1]")],
    "ALL literature-comparison assertions passed",
)

# Extra: breaking the attainment C-matrix equality must fail.
must_fail(
    "tests/verify_attainment_independent.py",
    [("assert np.max(np.abs(C_c - I2 / 3)) < 1e-12",
      "assert np.max(np.abs(C_c - I2 / 3)) < 1e-30")],
    "PASS: classical endpoint",
)

print("PASS: all supporting-check tamperings are rejected (no demonstration passes as verification).")
