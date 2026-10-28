"""Permutation and matched historic-window controls. Requires run_experiment.py."""
import json, math, itertools, random, sys
from pathlib import Path
import numpy as np
import run_experiment as ex
P=Path(__file__).resolve().parent
fb,checks,records,gv=ex.create_data()
train=list(range(4,10));test=[10,11,12]
alpha=10.
observed=ex.eval_preds(ex.pred_tail('yy',train,test,fb,alpha),fb,test)
rng=random.Random(13579)
perms=[tuple(train)] + rng.sample(list(itertools.permutations(train)),199)
null=[]
for perm in perms:
  mod={g:dict(x) for g,x in fb.items()}
  for g,pg in zip(train,perm): mod[g]['yy1']=fb[pg]['yy1']
  value=ex.eval_preds(ex.pred_tail('yy',train,test,mod,alpha),fb,test)
  null.append({'median':value['median_relative_error'], 'p90':value['p90_relative_error']})
null=np.array([(q['median'],q['p90']) for q in null[1:]])
# Compare to the earlier strong historical baseline, to avoid cherry-picking g4..9.
benchmark={}
for start,end in [(4,9),(5,9),(6,9),(7,10),(8,11)]:
    if end>9: # this is only for retrospective moving-window evaluation
       test_g=[end+1]
    else:test_g=[10,11,12]
    train_g=list(range(start,end+1))
    for mode in ['base2','base3','yy','yy+gv']:
       a={'base2':1,'base3':1,'yy':10,'yy+gv':0.001}[mode]
       pred=ex.pred_tail(mode,train_g,test_g,fb,a)
       # if 4 parameters with 3 observations model may be ill-posed; keep diagnostics.
       benchmark[f'{start}-{end}/{mode}']={'train':train_g,'test':test_g,**{k:v for k,v in ex.eval_preds(pred,fb,test_g).items() if k!='all_relative_errors'}}
result={'permutation_control':{'unshuffled_median':observed['median_relative_error'],
             'unshuffled_p90':observed['p90_relative_error'],
             'null_size':len(null),'null_median_distribution':np.quantile(null[:,0],[0,.05,.25,.50,.75,.95,1]).tolist(),
             'unshuffled_beats_null_median_fraction':float(np.mean(null[:,0]>observed['median_relative_error'])),
             'unshuffled_beats_null_p90_fraction':float(np.mean(null[:,1]>observed['p90_relative_error']))},
        'window_sensitivity':benchmark}
(P/'statistical_controls.json').write_text(json.dumps(result,indent=2)+'\n')
print('Permutation result:',result['permutation_control'])
for k,v in benchmark.items():print(k,'median',v['median_relative_error'],'p90',v['p90_relative_error'])
