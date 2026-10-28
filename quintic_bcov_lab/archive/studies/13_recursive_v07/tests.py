#!/usr/bin/env python3
"""Tests exact replay with withheld file physically absent, then validates separately."""
import csv,hashlib,subprocess,sys,tempfile,shutil
from pathlib import Path
from validate import check
ROOT=Path(__file__).resolve().parent

def test():
    expected=(ROOT/'results/predictions_Q26_Q30.csv').read_bytes()
    with tempfile.TemporaryDirectory() as tmp:
        d=Path(tmp);(d/'data').mkdir()
        shutil.copy2(ROOT/'recursive_defects.py',d/'recursive_defects.py')
        for name in ['calibration_Q1_Q20.csv','update_Q21_Q25.csv','sites_Q26_Q30.csv']:
            shutil.copy2(ROOT/'data'/name,d/'data'/name)
        assert not (d/'data/reference_only_Q26_Q30.csv').exists()
        subprocess.run([sys.executable,str(d/'recursive_defects.py'),'--calibration',str(d/'data/calibration_Q1_Q20.csv'),
              '--update',str(d/'data/update_Q21_Q25.csv'), '--sites',str(d/'data/sites_Q26_Q30.csv'),
              '--out',str(d/'test_prediction.csv')],check=True)
        actual=(d/'test_prediction.csv').read_bytes()
        assert actual==expected,(hashlib.sha256(actual).hexdigest(),hashlib.sha256(expected).hexdigest())
    metrics=check()
    assert metrics['defect_pulses']=={'zero':11,'one':4,'two':9}
    assert metrics['model_scores']['quadratic_PT']['exact']==15
    assert metrics['model_scores']['cubic_PT']['exact']==11
    assert metrics['model_scores']['linear_times_quadratic_PT']['exact']==10
    assert metrics['model_scores']['mirror_PT']['exact']==17
    assert metrics['defect_interval']['covered']==24
    print('All no-oracle replay, exact pulse reconstruction and benchmark tests passed.')
if __name__=='__main__':test()
