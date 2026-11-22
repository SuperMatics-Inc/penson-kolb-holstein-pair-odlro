# pyright: reportMissingImports=false
#!/usr/bin/env python3
"""
mutants.py -- adversarial mutation audit for the replay supplement.

The replay script carries seven semantic mutants, one for each class of statement the
paper relies on. Each mutant is a deliberately wrong variant of one construction in the
argument; the checks in `replay.py` must refuse it. This script runs each mutant in its
own process and verifies three things:

  1. the process exits with status 1,
  2. its diagnostic starts with "REJECTED:", and
  3. the diagnostic reports the quantity and the sign that the paper says is broken,
     so that a mutant cannot be "rejected" for an unrelated numerical accident.

Exits 0 iff every mutant is rejected as expected, 1 otherwise.

    python mutants.py                 run every mutant through replay.py
    python mutants.py --only NAME     run a single mutant
    python mutants.py --list          list the mutant names and what each breaks
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
REPLAY = os.path.join(HERE, "replay.py")

# name -> (what the mutant corrupts, regex that the rejection must match)
EXPECTED: dict[str, tuple[str, str]] = {
    "omit-odd-phonon-parity": (
        "drop the odd-sublattice phonon parity from the reflection",
        r"min_eig_even=-\d\.\d{6}e-0[12]",
    ),
    "antisymmetric-phonon-reflection": (
        "reflect the oscillator coordinate with the wrong sign",
        r"min_eig_odd=-\d\.\d{6}e-0[12]",
    ),
    "zero-flux-hopping": (
        "replace the pi-flux hopping on the reflected plaquette by real zero-flux hopping",
        r"min_eig_odd=-\d\.\d{6}e-0[12]",
    ),
    "crossing-bond-phonon": (
        "couple one bond oscillator linearly across the reflection plane",
        r"min_eig=-\d\.\d{6}e\+00",
    ),
    "drop-uph-static-term": (
        "delete the static phonon term of the double commutator",
        r"violates matrix double-commutator sum rule by \+\d\.\d{6}e-01.*m_mut\^2=0\.0\d+ > 0",
    ),
    "excessive-phonon-coupling": (
        "omit the Franck-Condon suppression from the dimer variational energy",
        r"m_cert\^2=-\d\.\d+ <= 0; omitting Frank-Condon factor violates Gibbs bound",
    ),
    "above-xy-critical-temperature": (
        "evaluate the monotone floor at a temperature below the ceiling",
        r"m_cert\^2=-\d\.\d+ <= 0; reversed beta-monotonicity gives false positive",
    ),
}


def run_one(name: str) -> bool:
    """Run a single mutant and check exit status and diagnostic."""
    what, pattern = EXPECTED[name]
    env = dict(os.environ)
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    env.setdefault("PYTHONIOENCODING", "utf-8")
    proc = subprocess.run(
        [sys.executable, "-B", REPLAY, "--mutant", name],
        capture_output=True,
        text=True,
        check=False,
        env=env,
        cwd=HERE,
    )
    err = proc.stderr.strip()
    last = (err.splitlines() or [""])[-1]
    ok_status = proc.returncode == 1
    ok_prefix = last.startswith("REJECTED:")
    ok_detail = re.search(pattern, last) is not None
    verdict = "REJECTED" if (ok_status and ok_prefix and ok_detail) else "NOT REJECTED"
    print(f"  [{verdict}] {name}: {what}")
    print(f"      rc={proc.returncode}, diagnostic: {last}")
    if verdict == "REJECTED":
        return True
    # Explain which of the three requirements failed.
    if not ok_status:
        print(f"      expected exit status 1, got {proc.returncode}")
    elif not ok_prefix:
        print("      diagnostic does not start with REJECTED:")
    else:
        print(f"      diagnostic does not report the expected quantity (pattern {pattern!r})")
    return False


def main() -> int:
    parser = argparse.ArgumentParser(description="Adversarial mutation audit for the replay supplement.")
    parser.add_argument("--only", default=None, choices=sorted(EXPECTED), help="Run a single mutant.")
    parser.add_argument("--list", action="store_true", help="List the mutants and what each corrupts.")
    args = parser.parse_args()

    if args.list:
        for name, (what, _) in EXPECTED.items():
            print(f"{name:32s} {what}")
        return 0

    names = [args.only] if args.only else list(EXPECTED)
    print(f"Adversarial mutation audit ({len(names)} mutants)")
    results = [run_one(n) for n in names]
    caught = sum(results)
    print(f"\n{caught}/{len(results)} mutants rejected as expected.")
    if caught != len(results):
        print("MUTATION AUDIT FAILED.", file=sys.stderr)
        return 1
    print("MUTATION AUDIT PASSED.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
