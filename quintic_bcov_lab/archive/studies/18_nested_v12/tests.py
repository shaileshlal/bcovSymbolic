"""Exact reproduction and genuine predictor/reference isolation."""
import csv
import pathlib
import shutil
import subprocess
import sys
import tempfile
from fractions import Fraction as F
ROOT=pathlib.Path(__file__).resolve().parent

def test_isolation():
 files=['run.py','predict_terminal.py','predict_terminal_values.py','predict_diagonal.py','predict_r6.py']
 with tempfile.TemporaryDirectory() as tmp:
  tmp=pathlib.Path(tmp)
  for p in files:shutil.copy(ROOT/p,tmp/p)
  shutil.copytree(ROOT/'data',tmp/'data')
  (tmp/'results').mkdir()
  z=subprocess.run([sys.executable,str(tmp/'run.py')],cwd=tmp,capture_output=True,text=True)
  assert z.returncode==0,(z.stdout,z.stderr)
  for name in ['terminal_pulse_predictions.csv','terminal_Z_reconstruction.csv','diagonal_pulse_predictions.csv','r6_Z_predictions.csv']:
   assert (tmp/'results'/name).read_bytes()==(ROOT/'results'/name).read_bytes(),name
 print('Predictor-only isolated run: byte-identical predictions, no validator files')

def test_exact():
 z=subprocess.run([sys.executable,str(ROOT/'validation/score.py')],cwd=ROOT,capture_output=True,text=True)
 assert z.returncode==0,(z.stdout,z.stderr)
 print('All exact validation checks passed')

def test_empirical_normalization():
 z=subprocess.run([sys.executable,str(ROOT/'validation/normalize_observations.py')],cwd=ROOT,capture_output=True,text=True)
 assert z.returncode==0,(z.stdout,z.stderr)
 print('Mature-PT integer normalization reproduced')

if __name__=='__main__':
 test_isolation();test_empirical_normalization();test_exact()
