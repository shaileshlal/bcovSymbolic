"""Independent exact boundary projection of high-degree YY coefficients.

Run from the study root: python verify_boundary.py [--all]
Defaults to genus 7 and 9; --all adds 10. Uses no positive-degree GV.
"""
import argparse,json,sys,time
from fractions import Fraction
from pathlib import Path
P=Path(__file__).resolve().parent
sys.path.insert(0,str(P/'quintic_ai_predictor'))
from hybrid.boundary import boundary_system
from engine.bcov_yy import Poly

def main():
 p=argparse.ArgumentParser();p.add_argument('--all',action='store_true');args=p.parse_args()
 data=json.loads((P/'quintic_ai_predictor/reference_amplitudes.json').read_text())['amplitudes']
 d={row['genus']:row for row in data}
 rows=[]
 for g in ([7,9,10] if args.all else [7,9]):
  t=time.monotonic();a=boundary_system(g,Poly.from_data(d[g]['anomaly']),mum_degrees=[0])
  truth=list(map(Fraction,d[g]['ambiguity']))
  assert a.verify(truth)
  free=set(a.free_indices)
  fixed=[j for j in range(len(truth)) if all(vec[j]==0 for vec in a.kernel)]
  assert all(a.particular[j]==truth[j] for j in fixed)
  rows.append({'genus':g,'rank':a.rank,'total':len(truth),'free':sorted(free),
               'high_fixed_from':max(free)+1,'all_fixed_match':True,
               'seconds':round(time.monotonic()-t,3)})
 print(json.dumps(rows,indent=2))
if __name__=='__main__':main()
