"""Interactive quintic PT shell-pulse notebook. Run: marimo edit recursive_marimo.py.

Uses reference labels ONLY for post-hoc figures and exact identities.
Forecasts were generated in a separate script with no held-out labels.
"""
from __future__ import annotations
import marimo
app=marimo.App(width='full')

@app.cell
def _():
    import marimo as mo
    import matplotlib.pyplot as plt
    import pandas as pd
    import math
    from pathlib import Path
    from recursive_defects import read,groups,phase_fit,evalpoly
    ROOT=Path(__file__).resolve().parent
    return mo,plt,pd,math,ROOT,read,groups,phase_fit,evalpoly

@app.cell
def _(mo):
    mo.md(r"""
    # Recursive defect shells in near-extremal quintic PT integers

    **Purely empirical:** neither modular forms nor DT/PT wall-crossing constraints enter the predictions.
    All training and validation use exact tabulated PT integers; the reference data are independent of the prediction script.

    \[
    T_{s,t}(k)=Q^{(2)}_{s,t}(k)+\sum_{j\ge k_0+3}
       \varepsilon_{s,t}(j)\binom{k-j+2}{2}\mathbf 1_{k\ge j},
    \quad \varepsilon=\Delta^3 T.
    \]
    The finite-support corrections are switched-on **quadratic splines**, not permanent cubic extrapolations.
    """)
    return

@app.cell
def _(mo):
    residue=mo.ui.dropdown(options={f'Residue {s}':s for s in range(5)},value='Residue 2',label='Degree mod 5')
    depth=mo.ui.slider(start=0,stop=5,value=1,step=1,label='Shell depth t')
    mo.hstack([residue,depth]);return residue,depth

@app.cell
def _(ROOT,read,groups,phase_fit):
    cal=read(ROOT/'data/calibration_Q1_Q20.csv')
    upd=read(ROOT/'data/update_Q21_Q25.csv')
    ref=read(ROOT/'data/reference_only_Q26_Q30.csv')
    period,phase=phase_fit(cal)
    G=groups(cal+upd+ref,phase,period)
    return cal,upd,ref,G

@app.cell
def _(residue,depth,G,mo,plt,pd,evalpoly,math):
    s=residue.value;t=depth.value
    seq=[(k,v) for k,v in G.get((s,t),[]) if 5*k+s>=6]
    if len(seq)!=5:
        mo.md('No observed stream at this depth; choose a smaller shell.')
    else:
        k0=seq[0][0];j1=k0+3;j2=k0+4
        x=[v for k,v in seq]
        e1=x[3]-3*x[2]+3*x[1]-x[0]
        e2=x[4]-3*x[3]+3*x[2]-x[1]
        rows=[]
        for k,y in seq:
            base=evalpoly(seq[:3],k)
            b1=e1*math.comb(k-j1+2,2) if k>=j1 else 0
            b2=e2*math.comb(k-j2+2,2) if k>=j2 else 0
            assert base+b1+b2==y
            rows.append(dict(k=k,signed_PT=y,quadratic=base,pulse_1=b1,pulse_2=b2))
        f,ax=plt.subplots(figsize=(7.5,3.5))
        ax.plot([v['k'] for v in rows],[v['signed_PT'] for v in rows],'o-',label='Signed PT')
        ax.plot([v['k'] for v in rows],[v['quadratic'] for v in rows],'s--',label='Initial quadratic')
        ax.set(xlabel='Residue position k',ylabel='Exact integer converted to plotted float',title=f'Residue {s}, shell {t}: pulse corrections')
        ax.legend();f.tight_layout()
        mo.vstack([mo.md(f'Observed first and second pulses: **{e1}**, **{e2}**'),mo.as_html(f),mo.ui.table(pd.DataFrame(rows))])
    return

@app.cell
def _(mo,ROOT,pd):
    forecasts=pd.read_csv(ROOT/'results/forecasts_Q31_Q35.csv',dtype=str)
    score=pd.read_csv(ROOT/'results/predictions_scored_Q26_Q30.csv',dtype=str)
    mo.vstack([mo.md('## Forward-degree holdout (26–30): actual vs competing models'),mo.ui.table(score),
      mo.md('## Prospective one-step forecasts (31–35), **not verified**'),
      mo.ui.table(forecasts[['degree','shell','quadratic_PT','mirror_PT','pulse_state','defect_interval_low','defect_interval_high']])])
    return

if __name__=='__main__':app.run()
