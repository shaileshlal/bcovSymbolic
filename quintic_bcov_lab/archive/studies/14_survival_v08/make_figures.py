#!/usr/bin/env python3
"""Reproducible benchmark graphics; matplotlib only."""
import json,csv
from pathlib import Path
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
REPORT=json.loads((ROOT/'results/metrics.json').read_text())
blocks=['31–35','36–40','41–45','46–50','51–53']
fig,ax=plt.subplots(figsize=(8,4.5))
x=range(len(blocks))
for mode,offset,name in [('teacher',-0.18,'Sequential updates'),('autonomous',0.18,'Autonomous')]:
 vals=[z['exact'] for z in REPORT['by_block'] if z['mode']==mode]
 ax.bar([i+offset for i in x],vals,width=0.36,label=name)
ax.set_xticks(list(x),blocks)
ax.set_ylim(0,27)
ax.set_ylabel('Exact PT predictions per block')
ax.set_xlabel('Curve degree')
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(ROOT/'figures/exact_predictions_by_block.png',dpi=150)
plt.close(fig)

supports=REPORT['defect_support_by_block']
fig,ax=plt.subplots(figsize=(8,4))
xs=[z['degree_start'] for z in supports]
ys=[z['n_active'] for z in supports]
ax.plot(xs,ys,marker='o')
ax.set_xticks(xs,[f"{z['degree_start']}–{z['degree_end']}" for z in supports]);ax.set_ylim(-0.3,10)
ax.set_ylabel('Nonzero third finite differences')
ax.set_xlabel('Curve degree block')
fig.tight_layout()
fig.savefig(ROOT/'figures/defect_front_decay.png',dpi=150)
plt.close(fig)
print('figures saved')
