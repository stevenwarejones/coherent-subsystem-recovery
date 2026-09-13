"""Anti-degradation tests for the supporting checks (review R2, tightened per follow-up).

A supporting script that merely PRINTS a comparison lets CI treat a demonstration as
verification. These tests tamper with a load-bearing value in a temporary copy of a script
and require it to fail FOR THE INTENDED REASON: nonzero exit, no success line, AND the
target assertion's specific marker in stderr. Requiring only "some nonzero exit" would count
an unrelated startup crash or a missing dependency as a successful rejection -- the controls
below show this harness does not. The repository checkout is left unchanged.
"""
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def run_text(script_rel, text):
    """Run `text` as a copy of script_rel (same filename, from ROOT) and return (rc, out, err)."""
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / Path(script_rel).name
        p.write_text(text)
        r = subprocess.run([sys.executable, str(p)], capture_output=True, text=True, cwd=str(ROOT))
        return r.returncode, r.stdout, r.stderr


def original(script_rel):
    return (ROOT / script_rel).read_text()


def mutate(text, replacements):
    for old, new in replacements:
        assert old in text, ("anchor not found", old)
        text = text.replace(old, new)
    return text


def must_fail(script_rel, replacements, success_marker, expected_error):
    """The mutation must fail for its intended reason: nonzero exit, no success line, and the
    specific marker present in stderr."""
    text = mutate(original(script_rel), replacements)
    assert text != original(script_rel), (script_rel, "no change applied")
    rc, out, err = run_text(script_rel, text)
    assert rc != 0, (script_rel, "tampered copy was accepted", out[-800:])
    assert success_marker not in out, (script_rel, "printed success despite tampering")
    assert expected_error in err, (script_rel, "did not fail for the intended reason", expected_error, err[-800:])
    print(f"PASS mutation: {script_rel} -> {expected_error} (exit {rc}, no success line)")


TARGETS = {
    "research/check_literature_comparisons.py": "ALL literature-comparison assertions passed",
    "tests/verify_forward_independent.py": "PASS independent:",
    "tests/verify_attainment_independent.py": "PASS: classical endpoint",
}

# Baseline: every unmodified target passes in this harness's run layout, so a failure below is
# caused by the tampering and not by how the copy is staged.
for rel, marker in TARGETS.items():
    rc, out, err = run_text(rel, original(rel))
    assert rc == 0 and marker in out, ("baseline must pass", rel, err[-800:])
print("PASS baseline: all target scripts pass unmodified in the harness layout")

# Intended-reason mutations (reproducing the two false acceptances the review found, plus two more).
must_fail("research/check_literature_comparisons.py",
          [("Ppg = sp.nsimplify((2 * Rpg + 1) / 3)", "Ppg = sp.Integer(1)")],
          TARGETS["research/check_literature_comparisons.py"], "E-BCW-VALUES")

must_fail("research/check_literature_comparisons.py",
          [('signs = {"X": [1, 1, -1, -1]', 'signs = {"X": [1, 1, 1, -1]')],
          TARGETS["research/check_literature_comparisons.py"], "E-STEERING-CORRELATION")

must_fail("tests/verify_forward_independent.py",
          [("p_pred=(5+4*math.cos(t))/9; r_pred=(1+math.cos(t))**2/4", "p_pred=0.0; r_pred=0.0")],
          TARGETS["tests/verify_forward_independent.py"], "E-FORWARD-PHASE")

# Mathematical mismatch (not an FP-threshold trick): assert C == I/2 when it is really I/3.
must_fail("tests/verify_attainment_independent.py",
          [("assert np.max(np.abs(C_c - I2 / 3)) < 1e-12", "assert np.max(np.abs(C_c - I2 / 2)) < 1e-12")],
          TARGETS["tests/verify_attainment_independent.py"], "E-ATTAINMENT-C")

# Controls: an UNRELATED failure must NOT be mistaken for a successful mutation rejection.
# For each target and its marker, prepend an unrelated crash and confirm the marker is ABSENT
# from stderr -- i.e. must_fail's marker requirement would (correctly) refuse to accept it.
MARKERS = {
    "research/check_literature_comparisons.py": "E-BCW-VALUES",
    "tests/verify_forward_independent.py": "E-FORWARD-PHASE",
    "tests/verify_attainment_independent.py": "E-ATTAINMENT-C",
}
for rel, marker in MARKERS.items():
    for prefix, label in [("raise RuntimeError('UNRELATED_STARTUP')\n", "RuntimeError"),
                          ("assert False, 'UNRELATED_ASSERTION'\n", "AssertionError")]:
        rc, out, err = run_text(rel, prefix + original(rel))
        assert rc != 0, (rel, label, "unrelated failure did not exit nonzero")
        assert marker not in err, (rel, label, "unrelated failure carried the intended marker")
        assert TARGETS[rel] not in out, (rel, label, "unrelated failure printed success")
    print(f"PASS control: unrelated RuntimeError/AssertionError in {rel} lack {marker} (not counted as rejection)")

print("PASS: mutations rejected for their intended reasons; unrelated failures are not counted.")
