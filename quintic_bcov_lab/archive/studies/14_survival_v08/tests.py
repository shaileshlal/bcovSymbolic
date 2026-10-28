#!/usr/bin/env python3
"""End-to-end audit, including a physically isolated, no-reference rollout."""
import os,sys,csv,subprocess,tempfile,json,shutil,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent
PINNED='e180b2d314fe91111c0ea9dabdc5aeb9472c375708ea29903f0d99ac697223d0'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def test():
    assert sha(ROOT/'predictor/recursive_defects.py')==PINNED
    b=json.loads((ROOT/'results/metrics.json').read_text())
    assert (b['teacher']['exact'],b['teacher']['n'])==(106,109)
    assert (b['autonomous']['exact'],b['autonomous']['n'])==(99,109)
    assert b['baselines']['quadratic_PT']['exact']==102
    assert b['baselines']['cubic_PT']['exact']==93
    assert not b['defect_reactivations']
    assert [z['n_active'] for z in b['defect_support_by_block']]==[9,4,2,1,0,0]
    # Physically delete access to the validator and reference for the predictor.
    with tempfile.TemporaryDirectory() as t:
        wd=Path(t)
        shutil.copytree(ROOT/'predictor',wd/'predictor',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        subprocess.run([sys.executable,str(wd/'predictor/rollout.py'),
            '--out',str(wd/'pred.csv')],check=True,cwd=wd,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
        orig=ROOT/'results/autonomous_predictions_Q31_Q53.csv'
        assert (wd/'pred.csv').read_bytes()==orig.read_bytes(),'Predictor changed when reference was absent'
    scored=list(csv.DictReader((ROOT/'results/individual_predictions_scored.csv').open()))
    assert len(scored)==218
    assert sum(int(z['exact']) for z in scored if z['mode']=='teacher')==106
    assert sum(int(z['exact']) for z in scored if z['mode']=='autonomous')==99
    geometric=list(csv.DictReader((ROOT/'results/extremal_geometry_Q35_Q50.csv').open()))
    assert [int(z['leading_PT']) for z in geometric]==[-550,725,925,-1150]
    print('PASS: pinned code, autonomous no-oracle recreation, 109-site scores, defects, geometry')
if __name__=='__main__':test()
