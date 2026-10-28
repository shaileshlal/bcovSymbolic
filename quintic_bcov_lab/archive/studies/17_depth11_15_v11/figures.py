import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
z=json.loads((ROOT/'results/benchmark.json').read_text())
depths=list(range(11,16))
counts=[z['prediction']['by_depth'][str(x)]['exact'] for x in depths]
eligible=[13,12,10,7,4]
fig,ax=plt.subplots(figsize=(7.3,4.1))
ax.bar(depths,eligible,label='In affine-profile range')
ax.bar(depths,counts,label='Uniquely predicted exact')
ax.set(xlabel='Shell depth t',ylabel='Number of PT sites (of 13)',ylim=(0,14))
ax.set_xticks(depths)
ax.legend()
fig.tight_layout()
fig.savefig(ROOT/'figures/depth11_15_coverage.png',dpi=160)
plt.close(fig)
