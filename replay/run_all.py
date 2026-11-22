# pyright: reportMissingImports=false
#!/usr/bin/env python3
"""
run_all.py -- run the whole replay supplement and report every exit status.

Runs, in order:
  replay.py            the interval certificates and the finite checks,
  replay.py --mutants  the seven mutants inside the replay script,
  probe_xover.py       the positivity boundary in the Holstein coupling,
  mutants.py           the independent mutation audit,

prints a summary table, and exits 0 only if all four commands exit 0.

    python run_all.py            run everything
    python run_all.py --quiet    show only the summary

The exit status of each command is reported, so a partial run cannot look like a pass.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))


def run(cmd: list[str], quiet: bool) -> tuple[int, float]:
    env = dict(os.environ)
    env.setdefault("PYTHONDONTWRITEBYTECODE", "1")
    env.setdefault("PYTHONIOENCODING", "utf-8")
    t0 = time.time()
    proc = subprocess.run(
        [sys.executable, "-B"] + cmd,
        capture_output=True,
        text=True,
        check=False,
        env=env,
        cwd=HERE,
    )
    dt = time.time() - t0
    if not quiet:
        sys.stdout.write(proc.stdout)
        sys.stdout.write(proc.stderr)
        sys.stdout.flush()
    return proc.returncode, dt


def main() -> int:
    parser = argparse.ArgumentParser(description="Run the whole replay supplement.")
    parser.add_argument("--quiet", action="store_true", help="Print only the summary.")
    args = parser.parse_args()

    steps = [
        ("replay.py (certificates and finite checks)", ["replay.py"]),
        ("replay.py --mutants (in-script mutants)", ["replay.py", "--mutants"]),
        ("probe_xover.py (positivity boundary)", ["probe_xover.py"]),
        ("mutants.py (mutation audit)", ["mutants.py"]),
    ]
    results = []
    for label, cmd in steps:
        rc, dt = run(cmd, args.quiet)
        results.append((label, rc, dt))
        if rc != 0:
            # Keep going so the summary shows every failure at once.
            print(f"\n[FAILED rc={rc}] {label}", file=sys.stderr)

    print("\n" + "=" * 72)
    print("SUMMARY")
    print("=" * 72)
    for label, rc, dt in results:
        status = "ok" if rc == 0 else "FAILED"
        print(f"  {status:6s} rc={rc} {dt:7.2f}s  {label}")
    failed = [label for label, rc, _ in results if rc != 0]
    if failed:
        print(f"\n{len(failed)} of {len(results)} steps failed.", file=sys.stderr)
        return 1
    print(f"\nAll {len(results)} steps passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
