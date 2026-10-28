#!/usr/bin/env python3
"""Static figures generated only after posthoc validation; not training features."""
from pathlib import Path
import json,csv, math
import numpy as np
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parent
metrics=json.loads((ROOT/'results/validated_metrics.json').read_text())
events=list(csv.DictReader(open(ROOT/'results/defect_analysis_only.csv')))

for low,high,label in [(21,25,'earlier'),(26,30,'later')]:
  A=np.full((5,6),np.nan)
  for r in events:
    d=int(r['degree'])
    if low<=d<=high:
      A[int(r['residue']),int(r['shell'])]=int(r['active'])
  fig,ax=plt.subplots(figsize=(7.4,3.4))
  im=ax.imshow(np.ma.masked_invalid(A),aspect='auto',cmap='Greys',vmin=0,vmax=1)
  for i in range(5):
    for j in range(6):
      if not np.isnan(A[i,j]):ax.text(j,i,'×' if A[i,j] else '·',ha='center',va='center',fontsize=14,color='white' if A[i,j] else 'black')
  ax.set(xticks=list(range(6)),yticks=list(range(5)),xlabel='Shell depth t',ylabel='Residue s = degree mod 5',title=f'Nonzero third differences: degrees {low}–{high}')
  ax.invert_yaxis()
  fig.tight_layout();fig.savefig(ROOT/f'figures/defect_support_{label}.png',dpi=170);plt.close(fig)

fig,ax=plt.subplots(figsize=(7.2,3.8))
S=metrics['scores']
shells=np.arange(6)
baseline=S['rolling_quad_second_block']['by_shell']; new=S['new_mirror_bracket_second_block']['by_shell']
ax.bar(shells-0.18,[baseline[str(i)][0]/baseline[str(i)][1] for i in shells],width=.34,label='Rolling quadratic')
ax.bar(shells+0.18,[new[str(i)][0]/new[str(i)][1] for i in shells],width=.34,label='Mirror-curvature defect rule')
ax.set(xlabel='Shell depth t',ylabel='Fraction predicted exactly',ylim=(0,1.08),xticks=shells,title='One-block-ahead exact PT predictions · degrees 26–30')
ax.legend(frameon=False);fig.tight_layout();fig.savefig(ROOT/'figures/exact_accuracy_by_shell.png',dpi=170);plt.close(fig)

D={(int(r['residue']),int(r['shell']),int(r['degree'])>=26):int(r['third_difference']) for r in events}
fig,ax=plt.subplots(figsize=(7.2,4.1))
for s in range(5):
  rows=[]
  for t in range(6):
    a=D.get((s,t,False));b=D.get((s,t,True))
    if a and b:rows.append((t,abs(b/a)))
  if rows:ax.scatter([r[0] for r in rows],[r[1] for r in rows],label=f's={s}',s=36)
ax.set_yscale('log');ax.set(xticks=list(range(1,6)),xlabel='Shell depth t',ylabel='|next defect / previous defect|',title='Observed attenuation for defects persisting to the next block')
ax.legend(frameon=False,ncol=3);fig.tight_layout();fig.savefig(ROOT/'figures/defect_attenuation.png',dpi=170);plt.close(fig)

M=json.loads((ROOT/'results/event_metrics.json').read_text())
labels=['Positional logistic','Full-feature logistic','ExtraTrees','Random forest','Moving front']
keys=['logit_positional','logit_engineered','extratrees_engineered','randomforest_engineered','unit_speed_front']
fig,ax=plt.subplots(figsize=(8,3.8))
ax.bar(labels,[M[k]['accuracy'] for k in keys]);ax.set(ylim=(0,1.08),ylabel='Next-block defect classification accuracy',title='24 residue–shell states, degrees 26–30')
ax.tick_params(axis='x',rotation=25)
fig.tight_layout();fig.savefig(ROOT/'figures/defect_classifiers.png',dpi=170);plt.close(fig)
print('Saved four diagnostic figures + two support maps.')
