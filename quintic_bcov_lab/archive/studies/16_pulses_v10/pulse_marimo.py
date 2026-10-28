"""Interactive marimo notebook for v1.0 empirical PT scale pulses.

Install: pip install marimo matplotlib
Run: marimo edit pulse_marimo.py
The notebook reads the reference ONLY for visualization; the predictor does not.
"""
import marimo

__generated_with = "0.17.0"
app = marimo.App(width="medium")

@app.cell
def _():
    import marimo as mo
    from pathlib import Path
    from fractions import Fraction as F
    import csv
    import json
    import math
    import sys
    import matplotlib.pyplot as plt
    return mo,Path,F,csv,json,math,sys,plt

@app.cell
def _(mo):
    mo.md(r"""
    # Quintic PT pulse hierarchy — empirical v1.0

    Let \(d=5k+s\), \(W_{s,t}(k)=5k^2+(2s-5)k+c_s-2t\).
    The **sign-corrected scale** is
    \[
    A_{s,t}(k)=\frac{\sigma_s(k)\mathrm{PT}(5k+s,t)}{W_{s,t}(k)}.
    \]
    The stabilized value is \(B_{[s],t}\), shared by residue pairs
    \(1\leftrightarrow4\) and \(2\leftrightarrow3\).
    The fitted residual profile is
    \[
    A_{s,t}(k)-B_{[s],t}=(-1)^t(P_{s,r}+Q_{s,r}/W_{s,t}(k)),
    \quad r=K(s,t)-k.
    \]

    **Provenance:** all learner parameters use degree \(\le40\) PT integers.
    The degree 41–53 reference is conditional on the prior DT / wall-crossing
    calculation, and is opened in this notebook only to visualize validation.
    The profile hypothesis was formulated retrospectively.
    """)
    return

@app.cell
def _(mo):
    residue=mo.ui.dropdown(options=[str(i) for i in range(5)],value="1",label="Residue s")
    depth=mo.ui.slider(start=6,stop=10,step=1,value=8,label="Shell depth t")
    mo.hstack([residue,depth],gap=3)
    return residue,depth

@app.cell
def _(Path,csv,json,sys):
    root=Path(__file__).resolve().parent
    sys.path.insert(0,str(root/"predictor"))
    from recursive_defects import read,phase_fit
    from deep_envelope import W
    from pulse_profile import K,group
    cal=read(root/"data/shallow_calibration_Q1_Q20.csv")
    period,phase=phase_fit(cal)
    training=read(root/"data/calibration_deep_Q20_Q40.csv")
    reference=read(root/"validator/reference_Q20_Q53.csv")
    preds=read(root/"results/pulse_predictions_Q41_Q53.csv")
    params=json.loads((root/"results/fitted_pulses.json").read_text())
    summary=json.loads((root/"validator/summary.json").read_text())
    return root,period,phase,W,K,group,training,reference,preds,params,summary

@app.cell
def _(mo,residue,depth,summary):
    s=int(residue.value)
    t=int(depth.value)
    j=summary["by_depth"][str(t)]
    mo.md(f"**Depth {t}:** {j['exact']} exact predictions out of {j['total']} candidate coefficients, with {j['total']-j['predicted']} abstentions.\\n\\n"
          f"Selected residue **{s}**. Separate previously observed values from later conditional reference values below.")
    return s,t

@app.cell
def _(mo,s,t,preds):
    sel=[{'degree':z['degree'],'front_distance':z['front_distance'],
          'predicted_PT':z['prediction'] or 'abstain',
          'source':z['source'] or 'unresolved',
          'support':z['profile_support']} for z in preds if int(z['residue'])==s and int(z['shell'])==t]
    mo.ui.table(sel,page_size=12)
    return sel

@app.cell
def _(mo,s,t,params):
    p=[x for x in params['profiles'] if x['residue']==s]
    b=[x for x in params['B'] if x['shell']==t and x['class']==(0 if s==0 else 1 if s in (1,4) else 2)]
    mo.md('### Exact-rational fitted pulse parameters')
    mo.ui.table([{'r':x['front_distance'],'P':x['P'],'Q':x['Q'],
                  'calibration_support':x['support'],'origin':x['origin']} for x in p],page_size=9)
    if b:mo.md('Stabilized scale B = `'+b[0]['rational']+'` ('+b[0]['source']+').')
    return p,b

@app.cell
def _(mo,plt,F,s,t,K,W,phase,period,reference,training):
    raw=[]
    for z in training+reference:
        if not z.get('pt'):continue
        d=int(z['degree']);tt=int(z['shell']);k=d//5
        if d%5==s and tt==t and d>=25:
            w=W(s,t,k)
            if not w:continue
            A=F(int(z['pt']),phase[s,k%period]*w)
            raw.append((d,k,A,'earlier' if d<=40 else 'held-out reference'))
    seen={}
    for d,k,a,source in raw:seen[d]=(k,a,source)
    ordered=sorted((d,*v) for d,v in seen.items())
    fig,ax=plt.subplots(figsize=(7.0,3.5))
    if ordered:
        last=ordered[-1][2]
        for source in ('earlier','held-out reference'):
            vals=[(k,float(a-last)) for d,k,a,z in ordered if z==source]
            if vals:
                # Log of departure, plotted with robust explicit zero floor.
                x=[k for k,v in vals]
                y=[__import__('math').log10(abs(v)) if v else -2 for k,v in vals]
                ax.scatter(x,y,label=source)
        ax.axvline(K(s,t),linestyle='--',linewidth=1,label='stabilization front')
    ax.set(xlabel='k',ylabel='log10 |A(k)-A(last)| (0 marked -2)',title=f'Residue {s}, depth {t}')
    ax.legend(fontsize=8);ax.grid(alpha=0.2);fig.tight_layout()
    fig
    return

@app.cell
def _(mo,summary):
    mo.md('### Unresolved sites')
    mo.ui.table(summary['unresolved'],page_size=10)
    return

if __name__=='__main__':
    app.run()
