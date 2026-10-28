"""Reproducible figures; no seaborn; no manually chosen palettes."""
import csv, json, math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
R=P/'results'; (P/'figures').mkdir(exist_ok=True)
with (R/'metrics.json').open() as f:metrics=json.load(f)
by=metrics['mirror_PT']['by_shell']
ks=sorted(map(int,by))
fig,ax=plt.subplots(figsize=(8.2,4.7))
ax.bar([str(k) for k in ks],[by[str(k)]['exact'] for k in ks])
ax.plot([str(k) for k in ks],[by[str(k)]['total'] for k in ks],marker='o',label='Available holdout coefficients')
ax.set(xlabel='Depth below first nonzero boundary coefficient',ylabel='Exact PT predictions',title='Where the empirical recurrence stops being exact')
ax.legend();fig.tight_layout();fig.savefig(P/'figures/holdout_exact_by_depth.png',dpi=160);plt.close(fig)

def read(p):
 with p.open() as f:return list(csv.DictReader(f))
all_rows=read(P/'data/calibration.csv')+read(P/'data/reference_only_Q21_Q30.csv')
info=json.loads((R/'model_audit.json').read_text())
T={(int(r['degree']),int(r['shell'])):int(r['pt']) for r in all_rows}
phase=lambda q:info['signs'][f'{q%5},{(q//5)%info["period"]}']
fig,ax=plt.subplots(figsize=(8.5,4.9))
for s,t,lab in [(0,4,'s=0, layer 4'),(4,3,'s=4, layer 3'),(3,2,'s=3, layer 2')]:
 qs=sorted(q for q,shell in T if shell==0 and q%5==s and q>=6)
 v=[phase(q)*T[q,t] for q in qs]
 dd=[v[i+2]-2*v[i+1]+v[i] for i in range(3)]
 # Normalize by the eventual (last two equal) curvature.
 ref=dd[-1]
 ax.plot([qs[i+2]//5 for i in range(3)],[d/ref for d in dd],marker='o',label=lab)
ax.axhline(1,linestyle='--',linewidth=1)
ax.set(xlabel='Right endpoint k of three-degree curvature window',ylabel='Curvature / late stabilized curvature',title='Near-boundary PT finite differences stabilize')
ax.legend();fig.tight_layout();fig.savefig(P/'figures/curvature_stabilization.png',dpi=160);plt.close(fig)
print('plots saved')
