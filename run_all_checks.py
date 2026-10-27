"""One-command convention, Bernoulli, denominator and enumerative regression suite."""
from pathlib import Path
import subprocess,sys
BASE=Path(__file__).resolve().parent/'bcov'
CHECKS=[
 ['test_frames.py'],
 ['test_gu_gauge.py'],
 ['bcov_yy.py','--genus','3'],
 ['audit_coefficients.py'],
 ['audit_denominators.py','--genus','8'],
 ['verify_physical_qseries.py','--genus','15','--degree','6'],
]
for call in CHECKS:
 print('\n>>>',sys.executable,*call,flush=True)
 subprocess.run([sys.executable,*call],cwd=BASE,check=True)
print('\nALL BCOV v0.2 REGRESSIONS PASSED',flush=True)
