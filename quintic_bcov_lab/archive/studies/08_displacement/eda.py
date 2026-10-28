"""Leakage-safe exploratory analysis: label geometry, periodicity, frontier."""
from __future__ import annotations
from pathlib import Path
from collections import Counter
import csv, json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from build_data import residual_dimension, margin10,top_extremal
ROOT=Path(__file__).resolve().parent


def status(g,d):
    m=margin10(g,d)
    if m<0:return 'certified_zero'
    if m==0 and d%5==0 and d>=20:return 'extremal_nonzero'
    return 'unresolved'

windows={}
for end in (12,50,53,54,60,100,200,500,1000):
    stat=Counter(status(g,d) for g in range(4,end+1)
                 for d in range(1,residual_dimension(g)+1))
    windows[str(end)]={k:stat[k] for k in ('certified_zero','extremal_nonzero','unresolved')}

frontier=[]
for g in range(49,101):
    unknown=[d for d in range(1,residual_dimension(g)+1) if status(g,d)=='unresolved']
    top=[d for d in range(1,residual_dimension(g)+1) if status(g,d)=='extremal_nonzero']
    frontier.append({'g':g,'D':residual_dimension(g),'unknown_count':len(unknown),
                     'min_unknown_d':min(unknown) if unknown else None,
                     'first_10_unknown_d':unknown[:10], 'proven_extremal_degrees':top})

# Distinct trends, not a single accidental 'g periodicity': gmod5 is 1 on proven top family,
# rmod5=4 there, while the sign is four-periodic in k=d/5 or twenty-periodic in r.
top=list(csv.DictReader((ROOT/'data'/'extremal_displacements.csv').open()))
mod_tables={}
for k in (2,3,4,5,6,7,8):
    counts={}
    for r in top:
        a=int(r['k'])%k
        counts.setdefault(str(a),set()).add(1 if int(r['n'])>0 else -1)
    mod_tables[str(k)]={key:sorted(v) for key,v in sorted(counts.items())}

trend={
    'counts_by_max_genus':windows,
    'first_unresolved':next(r for r in frontier if r['unknown_count']>0),
    'neighborhood_examples':[r for r in frontier if r['g'] in (50,51,52,53,54,55,56,57,60,75,76,100)],
    'top_family':{
        'n_positive':sum(int(r['n'])>0 for r in top),
        'n_negative':sum(int(r['n'])<0 for r in top),
        'g_mod5_values':sorted(set(int(r['g'])%5 for r in top)),
        'r_mod5_values':sorted(set(int(r['r'])%5 for r in top)),
        'signs_in_order':[1 if int(r['n'])>0 else -1 for r in top],
        'sign_classes_by_modulus':mod_tables,
    },
    'not_trainable_claim':'No reference labels for the interior (positive margin10), so no empirical g-neighborhood law can be validated there.',
}
(ROOT/'results'/'eda_summary.json').write_text(json.dumps(trend,indent=2)+'\n')
with (ROOT/'data'/'frontier_status_summary.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['g','D','unknown_count','min_unknown_d','proven_top_d'])
    w.writeheader()
    for r in frontier:w.writerow({'g':r['g'],'D':r['D'],'unknown_count':r['unknown_count'],
        'min_unknown_d':r['min_unknown_d'],'proven_top_d':','.join(map(str,r['proven_extremal_degrees']))})

# Frontiers at low g with clear data provenance.
fig,ax=plt.subplots(figsize=(9,5.8))
for name,mark in [('certified_zero','.'),('unresolved','o'),('extremal_nonzero','*')]:
    coords=[(g,d) for g in range(45,131)
                   for d in range(15,residual_dimension(g)+1) if status(g,d)==name]
    if not coords:continue
    gg,dd=zip(*coords)
    ax.scatter(gg,dd,marker=mark,s=(9 if mark=='.' else 18 if mark=='o' else 110),
               alpha=(.6 if mark=='.' else .8),label=name)
ax.set(xlabel='genus g',ylabel='curve degree d = r+1',
       title='Residual GV lattice: proven zero vs unresolved vs proven extremal')
ax.legend(loc='upper left',frameon=False)
fig.tight_layout();fig.savefig(ROOT/'figures'/'geometric_boundary.png',dpi=160);plt.close(fig)

fig,ax=plt.subplots(figsize=(9,5))
ks=[int(r['k']) for r in top]
ns=[int(r['n']) for r in top]
ax.plot(ks,ns,marker='.',label='extremal exact integers')
train=[(int(r['k']),int(r['n'])) for r in top if r['split']=='train']
test=[(int(r['k']),int(r['n'])) for r in top if r['split']=='test']
ax.scatter(*zip(*train),s=35,label='train k=4..11')
ax.scatter(*zip(*test),s=40,marker='x',label='held-out k=12..19')
ax.axhline(0,lw=1,linestyle=':')
ax.set(xlabel='k = d/5',ylabel='integer displacement Δn = n_(g,5k)',
       title='Proven extremal quintic GV displacements: periodic sign, quadratic magnitude')
ax.legend(frameon=False);fig.tight_layout()
fig.savefig(ROOT/'figures'/'extremal_pattern.png',dpi=160);plt.close(fig)
print(json.dumps({'counts_1000':windows['1000'],
                  'first_unresolved':trend['first_unresolved'],
                  'top_phase4':mod_tables['4']},indent=2))
