#!/usr/bin/env python3
"""Reproducibility tests, including an isolated no-oracle predictor copy."""
from pathlib import Path
import sys,shutil,tempfile,subprocess,csv,json,hashlib,py_compile
ROOT=Path(__file__).resolve().parent

def run(*args,cwd=ROOT):
    r=subprocess.run([sys.executable,*map(str,args)],cwd=cwd,capture_output=True,text=True,check=True)
    return r.stdout

# Materialized archive contents are available; regenerate predictions exactly.
run('predictor/deep_envelope.py')
run('predictor/gated_recurrence.py')
run('validator/deep_reference.py')
# scorer assumes imports from validator and predictor, handled within its script.
print(run('validator/audit.py').splitlines()[-1])
m=json.loads((ROOT/'results/metrics.json').read_text())
assert m['frozen_gated_exact']==32 and m['frozen_gated_predictions']==32
assert m['frozen_ungated_exact']==28 and m['deep_reference_safe_Q41_Q53']==65
assert m['quotient_defect_activity_accuracy']=='91/91'
assert m['third_defect_activity_accuracy']=='41/41'
with tempfile.TemporaryDirectory() as d:
    d=Path(d)
    shutil.copytree(ROOT/'predictor',d/'predictor',ignore=shutil.ignore_patterns('__pycache__'))
    (d/'data').mkdir()
    for name in ['shallow_calibration_Q1_Q20.csv','calibration_deep_Q20_Q40.csv','forecast_sites_Q41_Q53.csv']:
        shutil.copy(ROOT/'data'/name,d/'data'/name)
    # There is NO validator or reference in this temporary tree.
    assert not (d/'validator').exists()
    print(run('predictor/gated_recurrence.py',cwd=d).strip())
    assert (d/'results/frozen_gated_predictions_Q41_Q53.csv').read_bytes()==(ROOT/'results/frozen_gated_predictions_Q41_Q53.csv').read_bytes()
    print(run('predictor/deep_envelope.py',cwd=d).strip())
    assert (d/'results/predictions_paired_repeated.csv').read_bytes()==(ROOT/'results/predictions_paired_repeated.csv').read_bytes()
py_compile.compile(str(ROOT/'deep_shell_marimo.py'),doraise=True)
print('All exact tests and no-oracle replay PASSED.')
