"""Marimo explorer for the quintic PT deep-shell empirical baseline v0.9.

Run: marimo edit deep_shell_marimo.py
The reference data are conditional on the DT/wall-crossing framework and
are used only for visualization/validation. The predictor never imports them.
"""
import marimo

__generated_with = "0.17.0"
app = marimo.App(width="medium")

@app.cell
def _():
    import marimo as mo
    import csv
    import math
    from fractions import Fraction as F
    from pathlib import Path
    from collections import defaultdict
    return mo,csv,math,F,Path,defaultdict

@app.cell
def _(mo):
    mo.md(r"""
    # Quintic PT: recursive defects at shell depths 6–10

    Use \(d=5k+s\). The shallow-data quadratic envelope is
    \[
      W_{s,t}(k)=5k^2+(2s-5)k+c_s-2t,\qquad
      (c_0,\ldots,c_4)=(10,6,4,4,6).
    \]
    The sign-corrected PT scale \(A_{s,t}(k)=\sigma_s(k)\mathrm{PT}(d,t)/W_{s,t}(k)\)
    has an empirical diagonal stabilization front. All rationals below
    are handled exactly; plots use approximate log magnitudes only.

    **Provenance:** exact labels here are obtained from a *conditional*
    DT/wall-crossing validation solver; no modular data enter the empirical predictor.
    """)
    return

@app.cell
def _(mo):
    residue = mo.ui.dropdown(options=[str(x) for x in range(5)],value='1',label='Degree residue s')
    shell = mo.ui.slider(start=6,stop=10,step=1,value=7,label='Shell depth t')
    mo.hstack([residue,shell],gap=3)
    return residue,shell

@app.cell
def _(Path,csv,F,defaultdict):
    root = Path(__file__).resolve().parent
    def load(rel):
        with (root/rel).open(newline='') as f:
            return list(csv.DictReader(f))
    data=load('validator/reference_Q20_Q53.csv')
    front=load('results/normalized_scale_defects.csv')
    third=load('results/third_difference_front.csv')
    preds=load('results/frozen_gated_predictions_Q41_Q53.csv')
    # signed phases in the predictor are learned from the shallow data.
    import sys
    sys.path.insert(0,str(root/'predictor'))
    from recursive_defects import phase_fit
    period,phases=phase_fit(load('data/shallow_calibration_Q1_Q20.csv'))
    cvals=[10,6,4,4,6]
    def envelope(s,t,k):return 5*k*k+(2*s-5)*k+cvals[s]-2*t
    def threshold(s,t):return t-1 if s==0 else t+3-s
    return root,data,front,third,preds,period,phases,envelope,threshold,load

@app.cell
def _(residue,shell,mo,data,F,phases,period,envelope,threshold):
    s=int(residue.value);t=int(shell.value)
    rows=[]
    for r in data:
        d=int(r['degree'])
        if d%5!=s or int(r['shell'])!=t:continue
        k=d//5
        if r['pt']:
            w=envelope(s,t,k)
            a=F(int(r['pt']),phases[s,k%period]*w)
            rows.append({'degree':d,'k':k,'PT':r['pt'],
                         'normalized_A':str(a),'mature':k>=threshold(s,t)})
        else:
            rows.append({'degree':d,'k':k,'PT':'not certified','normalized_A':'','mature':k>=threshold(s,t)})
    mo.md(f'### Residue {s}, depth {t}: predicted last scale change at k = {threshold(s,t)}')
    mo.ui.table(rows,page_size=12)
    return s,t,rows

@app.cell
def _(mo,front,third,s,t):
    v=[r for r in front if int(r['residue'])==s and int(r['shell'])==t]
    u=[r for r in third if int(r['residue'])==s and int(r['shell'])==t]
    mo.md(f'**Scale-change tests:** {sum(int(r["activity"]) == int(r["front_prediction"]) for r in v)}/{len(v)} correct.  '
          f'**Third-difference tests:** {sum(int(r["activity"]) == int(r["front_prediction"]) for r in u)}/{len(u)} correct.')
    return v,u

@app.cell
def _(mo,math,F,rows):
    import matplotlib.pyplot as plt
    pairs=[]
    for r in rows:
        if not r['normalized_A']:continue
        a=F(r['normalized_A'])
        pairs.append((r['k'],a))
    fig,ax=plt.subplots(figsize=(7,3.2))
    if pairs:
        a_last=pairs[-1][1]
        v=[(k,math.log10(abs(a-a_last).numerator)-math.log10(abs(a-a_last).denominator))
           for k,a in pairs if a!=a_last]
        if v:
            ax.scatter([k for k,y in v],[y for k,y in v],label='log10 absolute scale residual')
        for k,a in pairs:
            if a==a_last:
                ax.axvline(k,alpha=0.1)
    ax.set(xlabel='k',ylabel='log10 |A(k)-A(last)|',title='Transient and stabilized normalized scale')
    ax.grid(alpha=0.2)
    fig.tight_layout()
    fig
    return

@app.cell
def _(mo,preds,s,t):
    predrows=[{'degree':r['degree'],'state':r['state'],'PT prediction':r['prediction'] or '(abstain)',
               'anchor degree':r['source_degree']} for r in preds if int(r['degree'])%5==s and int(r['shell'])==t]
    mo.md('### Frozen no-oracle predictor, curve degrees 41–53')
    mo.ui.table(predrows,page_size=10)
    return

if __name__=='__main__':
    app.run()
