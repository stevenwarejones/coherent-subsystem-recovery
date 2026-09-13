#!/usr/bin/env python3
"""Fail-fast, streaming runner for the required checks.

Checks are run in dependency order; output is streamed as each produces it, any
nonzero exit stops the run, and no success sentence is printed unless every
required check has already passed.

Required checks:
    proofs/verify_sharp_recovery.py       exact free-unitary certificate, R >= (3f-1)/2   [stdlib only]
    proofs/verify_full_curve.py           elementary segments completing R_min(f)          [SymPy]
    proofs/verify_forward_recovery.py     forward two-query guarantee r_UV >= (9f-5)/4      [SymPy]
    tests/test_certificate_mutations.py   intended-reason corruption cases                 [stdlib only]

The central exact certificate and the mutation cases need only the standard
library; the two SymPy checks above complete and extend the curve. The
independent re-derivations (SymPy/NumPy), the supporting-check mutations, and the
python-flint alternate route are run by CI and can be run manually
(see docs/VERIFICATION.md); they are invoked separately so a missing optional
dependency never masks a failure of the core checks.
"""
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent

REQUIRED = [
    ("Exact free-unitary certificate (stdlib)", ROOT / "proofs" / "verify_sharp_recovery.py"),
    ("Elementary full-curve segments (SymPy)", ROOT / "proofs" / "verify_full_curve.py"),
    ("Forward-only two-query recovery guarantee (SymPy)", ROOT / "proofs" / "verify_forward_recovery.py"),
    ("Certificate corruption / mutation cases (stdlib)", ROOT / "tests" / "test_certificate_mutations.py"),
]


def run(label, path):
    print(f"\n=== {label}\n=== {path.relative_to(ROOT)}", flush=True)
    proc = subprocess.Popen([sys.executable, "-u", str(path)], cwd=str(ROOT))
    proc.wait()
    if proc.returncode != 0:
        print(f"\nFAILED: {path.relative_to(ROOT)} exited {proc.returncode}", flush=True)
        sys.exit(proc.returncode)


def main():
    for label, path in REQUIRED:
        run(label, path)
    print(
        "\nAll required checks passed: the exact certificate for R >= (3f-1)/2, the "
        "elementary segments completing R_min(f) = max{1/4, 3f/4, (3f-1)/2}, the "
        "forward-only two-query guarantee r_UV >= (9f-5)/4, and the intended-reason "
        "mutation cases.\n"
        "This is exact-arithmetic verification, not human expert review or a priority "
        "claim. See docs/VERIFICATION.md for what these checks do and do not establish."
    )


if __name__ == "__main__":
    main()
