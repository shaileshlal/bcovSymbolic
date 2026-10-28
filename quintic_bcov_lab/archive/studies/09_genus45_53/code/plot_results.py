from pathlib import Path
import csv
import matplotlib.pyplot as plt
root=Path(__file__).resolve().parents[1]
rows=list(csv.DictReader((root/'data'/'residue_family_benchmark.csv').open()))
scales={1:8625/2,2:609250,3:609250,4:8625/2}
fig,ax=plt.subplots(figsize=(9.3,5.4))
for s in [1,2,3,4]:
    subset=sorted((r for r in rows if int(r['residue'])==s),key=lambda r:int(r['k']))
    xs=[int(r['k']) for r in subset]
    ys=[abs(int(r['reference_integer']))/scales[s] for r in subset]
    line=ax.plot(xs,ys,label=f'd ≡ {s} (mod 5)',lw=1.7)[0]
    ax.scatter(xs[:3],ys[:3],marker='o',s=35,color=line.get_color())
    ax.scatter(xs[3:],ys[3:],marker='x',s=80,color=line.get_color())
ax.set_xlabel('k in degree d=5k+s')
ax.set_ylabel('|GV charge| / residue-dependent scale')
ax.set_xticks([1,2,3,4,5])
ax.set_title('Quadratic magnitude envelopes; circles trained, crosses held out')
ax.grid(alpha=0.2)
ax.legend()
fig.tight_layout()
fig.savefig(root/'figures'/'periodic_integer_envelope.png',dpi=150)
plt.close(fig)
print('wrote',root/'figures'/'periodic_integer_envelope.png')
