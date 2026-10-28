"""Draw the normalized terminal pulses, without choosing matplotlib colors."""
from pathlib import Path
from collections import defaultdict
import csv
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parent
rows=[]
for name in ('data/terminal_pulse_calibration.csv','results/terminal_pulse_predictions.csv'):
 with (ROOT/name).open() as f:rows+=list(csv.DictReader(f))
by=defaultdict(list)
for z in rows: by[int(z['residue'])].append((int(z['r']),int(z['terminal_pulse'])))
fig,ax=plt.subplots(figsize=(8.5,4.7))
for s,rr in sorted(by.items()):
 rr.sort();scale=max(abs(x) for _,x in rr)
 ax.plot([r for r,v in rr],[v/scale for r,v in rr],marker='o',label=f's={s}')
ax.set(xlabel='Front distance r',ylabel='Pulse / maximum absolute pulse per residue',title='Empirical quadratic terminal pulses')
ax.grid(True,alpha=.2);ax.legend(ncol=5)
fig.tight_layout()
out=ROOT/'figures';out.mkdir(exist_ok=True)
fig.savefig(out/'terminal_pulse_quadratics.png',dpi=160)
plt.close(fig)
print(out/'terminal_pulse_quadratics.png')
