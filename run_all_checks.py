"""Run package tests and the historical convention/enumerative regression suite."""
from pathlib import Path
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
CHECKS: tuple[tuple[str, ...], ...] = (
    ("test_frames.py",),
    ("test_gu_gauge.py",),
    ("bcov_yy.py", "--genus", "3"),
    ("audit_coefficients.py",),
    ("audit_denominators.py", "--genus", "8"),
    ("verify_physical_qseries.py", "--genus", "15", "--degree", "6"),
)


def main() -> None:
    """Fail immediately when any exact regression fails."""
    subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-v"], cwd=ROOT, check=True)
    for call in CHECKS:
        print("\n>>>", sys.executable, *call, flush=True)
        subprocess.run([sys.executable, *call], cwd=ROOT/"bcov", check=True)
    print("\nALL BCOV REGRESSIONS PASSED", flush=True)


if __name__ == "__main__":
    main()
