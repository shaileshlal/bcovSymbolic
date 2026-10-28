#!/usr/bin/env python3
"""Exact self-tests, including physically isolated no-reference predictor execution."""
from __future__ import annotations
from pathlib import Path
from fractions import Fraction as F
import tempfile,shutil,subprocess,sys,csv,json
ROOT=Path(__file__).resolve().parent


def data(path):
    with (ROOT/path).open(newline='') as f:return list(csv.DictReader(f))


def main():
    current=data('results/pulse_predictions_Q41_Q53.csv')
    ref={(int(r['degree']),int(r['shell'])):int(r['pt']) for r in data('validator/reference_Q20_Q53.csv') if r['pt']}
    assert len(current)==65
    pred=[r for r in current if r['prediction']]
    assert len(pred)==59
    assert all(int(r['prediction'])==ref[int(r['degree']),int(r['shell'])] for r in pred)
    assert sorted((int(r['degree']),int(r['shell'])) for r in current if not r['prediction'])==[
        (42,10),(43,10),(47,10),(48,10),(52,10),(53,10)]
    assert any(int(z['degree'])==41 and int(z['shell'])==8 and z['prediction']=='534184449206140347572250' for z in pred)
    summary=json.loads((ROOT/'validator/summary.json').read_text())
    assert summary['v10_predicted']==summary['v10_exact']==59
    p=json.loads((ROOT/'results/fitted_pulses.json').read_text())
    prof={(z['residue'],z['front_distance']):(F(z['P']),F(z['Q'])) for z in p['profiles']}
    assert prof[0,2]==tuple(-v for v in prof[2,2])
    assert prof[3,2]==tuple(-v for v in prof[4,2])
    assert prof[1,2]==(F(-1135625),F(28750))
    # No reference or validator files exist in this subprocess sandbox.
    with tempfile.TemporaryDirectory(prefix='quintic_pulse_no_oracle_') as td:
        tmp=Path(td)
        shutil.copytree(ROOT/'data',tmp/'data')
        shutil.copytree(ROOT/'predictor',tmp/'predictor',ignore=shutil.ignore_patterns('__pycache__'))
        assert not (tmp/'validator').exists()
        subprocess.run([sys.executable,str(tmp/'predictor/pulse_profile.py')],cwd=tmp,check=True,capture_output=True,text=True)
        actual=(tmp/'results/pulse_predictions_Q41_Q53.csv').read_bytes()
        expected=(ROOT/'results/pulse_predictions_Q41_Q53.csv').read_bytes()
        assert actual==expected,'no-oracle output differs'
        assert (tmp/'results/fitted_pulses.json').read_bytes()==(ROOT/'results/fitted_pulses.json').read_bytes()
    print('PASS exact 59/59; abstentions 6; antisymmetry; one concrete PT reconstruction; isolated no-oracle replay')

if __name__=='__main__':main()
