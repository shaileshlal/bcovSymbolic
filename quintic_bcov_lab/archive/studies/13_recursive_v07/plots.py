#!/usr/bin/env python3
"""Rebuild figures from already validated integer data."""
import csv,json,math
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
s=json.loads((ROOT/'results/validated_metrics.json').read_text())
rows=list(csv.DictReader((ROOT/'results/defect_shells.csv').open()))
fig,ax=plt.subplots(figsize=(9.6,4.6))
status=np.full((6,5),np.nan)
for r in rows:
    st={'locked':0,'healed':1,'persistent':2}[r['defect_front_state']]
    status[int(r['shell']),int(r['residue'])]=st
im=ax.imshow(np.ma.masked_invalid(status),origin='lower',vmin=0,vmax=2,cmap='viridis',aspect='auto')
ax.set_xticks(range(5),[f'{s} (mod 5)' for s in range(5)])
ax.set_yticks(range(6));ax.set_xlabel('Degree residue');ax.set_ylabel('Shell depth t')
ax.set_title('Defect-pulse states in two forward blocks')
for r in rows:
    s0=int(r['residue']);t=int(r['shell']);stat=r['defect_front_state']
    ax.text(s0,t,('00' if stat=='locked' else '10' if stat=='healed' else '11'),ha='center',va='center',fontsize=9,color='white' if stat=='locked' else 'black')
cb=fig.colorbar(im,ax=ax,shrink=.84,ticks=[0,1,2]);cb.ax.set_yticklabels(['no pulse','pulse then heals','pulse persists'])
fig.tight_layout();fig.savefig(ROOT/'figures/pulse_shell_state.png',dpi=180);plt.close(fig)
mods=['quadratic_PT','cubic_PT','linear_times_quadratic_PT','mirror_PT']
fig,ax=plt.subplots(figsize=(9.4,4.2))
counts=[s['model_scores'][m]['exact'] for m in mods];den=[s['model_scores'][m]['total'] for m in mods]
labels=['Quadratic','Cubic','Linear × leading quadratic','Paired-curvature']
bars=ax.barh(labels,counts,color=plt.rcParams['axes.prop_cycle'].by_key()['color'][0]);ax.invert_yaxis();ax.set_xlim(0,25)
for bar,n,tot in zip(bars,counts,den): ax.text(n+.2,bar.get_y()+bar.get_height()/2,f'{n}/{tot}',va='center')
ax.set_xlabel('Exact integers (forward holdout degrees 26–30)');ax.set_title('Higher polynomial degree does not improve exact prediction')
fig.tight_layout();fig.savefig(ROOT/'figures/model_accuracy.png',dpi=180);plt.close(fig)
