from pathlib import Path
import csv,json,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from collections import defaultdict
ROOT=Path(__file__).parent

def read(path):
 with open(path,newline='') as f:return list(csv.DictReader(f))
rows=read(ROOT/'results/normalized_scale_defects.csv')
# Single scatter/heatmap for residue 1.
grid={(int(r['shell']),int(r['k'])):int(r['activity']) for r in rows if int(r['residue'])==1}
xs=range(5,11);ys=range(6,11)
fig,ax=plt.subplots(figsize=(7.3,3.9))
for t in ys:
 for k in xs:
  if (t,k) in grid:
   v=grid[t,k]
   ax.scatter(k,t,s=210,marker='s',alpha=0.85 if v else 0.25)
   ax.text(k,t,'1' if v else '0',ha='center',va='center',fontsize=9)
ax.set(xlabel='k = floor(degree/5)',ylabel='Shell depth t',title='Residue 1: normalized-scale defect activity (1=active)')
ax.set_xticks(list(xs));ax.set_yticks(list(ys));ax.grid(alpha=0.16)
fig.tight_layout();fig.savefig(ROOT/'figures/deep_front_residue1.png',dpi=170);plt.close(fig)

scored=read(ROOT/'results/frozen_gated_predictions_Q41_Q53.csv')
base=read(ROOT/'results/predictions_paired_repeated.csv')
ref={(int(r['degree']),int(r['shell'])):int(r['pt']) for r in read(ROOT/'validator/reference_Q20_Q53.csv') if r['pt']}
by=defaultdict(lambda:[0,0,0]);ab=defaultdict(int)
for r in scored:
 t=int(r['shell']);d=int(r['degree']);by[t][0]+=1
 if r['prediction']:
  by[t][1]+=1;by[t][2]+=int(int(r['prediction'])==ref[d,t])
for r in base:
 t=int(r['shell']);d=int(r['degree']);ab[t]+=int(int(r['prediction'])==ref[d,t])
xx=list(range(6,11));fig,ax=plt.subplots(figsize=(7.3,3.9))
ax.bar([x-0.2 for x in xx],[ab[x] for x in xx],width=0.4,label='Ungated, exact')
ax.bar([x+0.2 for x in xx],[by[x][2] for x in xx],width=0.4,label='Gated, exact')
ax.set(xlabel='Shell depth t',ylabel='Exact PT integers',title='Depth holdout Q41–53: integer accuracy')
ax.set_xticks(xx);ax.legend();fig.tight_layout();fig.savefig(ROOT/'figures/exact_deep_by_shell.png',dpi=170);plt.close(fig)
print('Figures generated')
