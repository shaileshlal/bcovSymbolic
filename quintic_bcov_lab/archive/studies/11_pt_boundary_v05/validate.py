#!/usr/bin/env python3
"""Comparison-only evaluator; NEVER imported by the empirical predictor."""
import csv,json,math
from collections import defaultdict,Counter
from fractions import Fraction
from pathlib import Path
from fit_empirical_transform import fit
from run_experiments import lagrange_quad

def read(path):
 with open(path) as f:return list(csv.DictReader(f))

def med(nums):
 nums=sorted(nums);n=len(nums);return (nums[(n-1)//2]+nums[n//2])/2

def run():
 pred=read('results/pt_predictions.csv')
 ref=read('data/reference_only_Q21_Q30.csv')
 train=read('data/calibration.csv')
 R={(int(r['degree']),int(r['shell'])):r for r in ref}
 T={(int(r['degree']),int(r['shell'])):r for r in train}
 results=[];t_by_d={}
 for p in pred:
  d=int(p['degree']);t=int(p['shell']);truth=int(R[d,t]['pt']);y=int(p['predicted_PT'])
  t_by_d[(d,t)]=y
  results.append(dict(degree=d,shell=t,method=p['model'],predicted_PT=y,true_PT=truth,
                      exact=y==truth,sign=int(y>0)-int(y<0),true_sign=int(truth>0)-int(truth<0),
                      relative_error=abs(y-truth)/max(1,abs(truth))))
 own_quad=[]
 # Construct independent quadratic baseline, phase fitted from calibration (no holdout).
 info=json.loads(Path('results/model_audit.json').read_text());period=info['period'];sign_map={tuple(map(int,k.split(','))):v for k,v in info['signs'].items()}
 for r in ref:
  d=int(r['degree']);t=int(r['shell']);s=d%5
  sample=[(int(q['degree']),int(q['pt'])) for q in train if 6<=int(q['degree'])<=20 and int(q['degree'])%5==s and int(q['shell'])==t]
  assert len(sample)==3
  f=lagrange_quad([q//5 for q,v in sample],[v*sign_map[(s,(q//5)%period)] for q,v in sample])
  baseline=int(f(d//5)*sign_map[(s,(d//5)%period)])
  own_quad.append(baseline==int(r['pt']))
 r1,r2,fitinfo=fit('data/calibration.csv')
 gv_rows=[]
 for d in range(21,31):
  h=int(R[d,0]['gc'])-int(R[d,0]['delta'])
  for t in range(3):
   p=t_by_d[(d,t)]
   if t==0:v=p
   if t==1:v=p-r1(h)*gv_rows[-1]['predicted_GV']
   if t==2:v=p-r1(h-1)*gv_rows[-1]['predicted_GV']-r2(h)*gv_rows[-2]['predicted_GV']
   assert v.denominator==1,(d,t,v)
   v=int(v);truth=int(R[d,t]['gv'])
   gv_rows.append(dict(degree=d,shell=t,predicted_GV=v,true_GV=truth,exact=v==truth))
 c=Counter(r['exact'] for r in results)
 stats={'holdout_sites':len(results),'degree_range':[21,30],
        'independent_PT_quadratic':{'exact':sum(own_quad),'total':len(own_quad)},
        'mirror_PT':{'exact':c[True],'total':len(results),'by_shell':{},'by_model':{}},
        'PT_signs':sum(r['sign']==r['true_sign'] for r in results),
        'PT_median_relative_error':med([r['relative_error'] for r in results]),
        'PT_median_relative_error_misses':med([r['relative_error'] for r in results if not r['exact']]),
        'PT_median_abs_log10_error':med([abs(math.log10(max(1,abs(r['predicted_PT'])))-math.log10(abs(r['true_PT']))) for r in results]),
        'PT_exact_by_degree_block':{},
        'empirical_PT_GV_transform':fitinfo,
        'GV_recovery_using_predicted_PT':{'exact':sum(r['exact'] for r in gv_rows),'total':len(gv_rows),'by_shell':{}},
       }
 for shell in sorted({r['shell'] for r in results}):
  z=[r for r in results if r['shell']==shell]
  stats['mirror_PT']['by_shell'][str(shell)]=dict(exact=sum(r['exact'] for r in z),total=len(z),median_relative_error=med([r['relative_error'] for r in z]))
 for model in ['quadratic','mirror_affine']:
  z=[r for r in results if r['method']==model]
  stats['mirror_PT']['by_model'][model]={'exact':sum(r['exact'] for r in z),'total':len(z)}
 for segment in [(21,25),(26,30)]:
  z=[r for r in results if segment[0]<=r['degree']<=segment[1]]
  stats['PT_exact_by_degree_block'][str(segment)]={'exact':sum(r['exact'] for r in z),'total':len(z)}
 for t in range(3):
  z=[r for r in gv_rows if r['shell']==t]
  stats['GV_recovery_using_predicted_PT']['by_shell'][str(t)]={'exact':sum(r['exact'] for r in z),'total':len(z)}
 Path('results').mkdir(exist_ok=True)
 with open('results/scored_predictions.csv','w',newline='') as f:
  w=csv.DictWriter(f,results[0]);w.writeheader();w.writerows(results)
 with open('results/gv_predictions.csv','w',newline='') as f:
  w=csv.DictWriter(f,gv_rows[0]);w.writeheader();w.writerows(gv_rows)
 Path('results/metrics.json').write_text(json.dumps(stats,indent=2))
 print(json.dumps(stats,indent=2))
 return stats
if __name__=='__main__':run()
