"""Run the formula-audit scripts shipped with this repository."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = [
    "validate_direct_oracle.py",
    "audit_WG1L_eta_form.py",
    "audit_WG2L_eta_form.py",
    "audit_all_ports_eta_form.py",
]


def main() -> int:
    failed = []
    for script in SCRIPTS:
        path = Path(__file__).resolve().parent / script
        print(f"\n=== {script} ===")
        result = subprocess.run([sys.executable, str(path)], cwd=ROOT)
        if result.returncode:
            failed.append((script, result.returncode))

    if failed:
        print("\nFormula audit failed:")
        for script, code in failed:
            print(f"  {script}: exit code {code}")
        return 1

    print("\nAll formula-audit scripts completed successfully.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
