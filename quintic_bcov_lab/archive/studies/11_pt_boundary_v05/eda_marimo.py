"""Marimo notebook: train-only phase and recurrence explorer.

Run `marimo edit eda_marimo.py` with `pip install marimo pandas matplotlib`.
The published validation table is loaded for visualization only, never into
`run_experiments.train_predict()`.
"""
from __future__ import annotations
import marimo
app=marimo.App(width='full')

@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import json
    from pathlib import Path
    import matplotlib.pyplot as plt
    ROOT=Path(__file__).resolve().parent
    return mo,pd,json,ROOT,plt

@app.cell
def _(mo):
    mo.md(r"""
    # Quintic near-Castelnuovo arithmetic — empirical discovery

    Training = PT/GV coefficients at degrees 1–20, held out by **whole degree**
    21–30. Periodic phases, exact finite differences, and mirrored residue
    transfer are learned without modular forms or wall-crossing identities.
    The published comparison data are modularity-conditional, but are never
    loaded by the predictor.
    """)
    return

@app.cell
def _(mo):
    s=mo.ui.dropdown(options={"0 (5k)":0,"1 (5k+1)":1,"2 (5k+2)":2,"3 (5k+3)":3,"4 (5k+4)":4},value="1 (5k+1)",label='Degree residue')
    t=mo.ui.slider(start=0,stop=3,step=1,label='Depth beyond first nonzero shell')
    mo.hstack([s,t],gap=2)
    return s,t

@app.cell
def _(pd,ROOT):
    cal=pd.read_csv(ROOT/'data/calibration.csv')
    cmp=pd.read_csv(ROOT/'data/reference_only_Q21_Q30.csv')
    preds=pd.read_csv(ROOT/'results/scored_predictions.csv')
    return cal,cmp,preds

@app.cell
def _(s,t,cal,cmp,preds,plt,mo,json,ROOT):
    v=int(s.value);j=int(t.value)
    # The sign phase is read from calibration-only model selection.
    audit=json.loads((ROOT/'results/model_audit.json').read_text())
    selected=pd.concat([cal,cmp]) if False else None
    # use marimo notebook's independent plotting imports
    import pandas as pd
    tab=pd.concat([cal,cmp],ignore_index=True)
    tab=tab[(tab['residue']==v)&(tab['shell']==j)].copy()
    if len(tab):
        tab['signed_PT']=[int(r.pt)*int(audit['signs'][f'{v},{(int(r.k))%4}']) for r in tab.itertuples()]
        fig,ax=plt.subplots(figsize=(7,3.8))
        train=tab[tab.degree<=20];test=tab[tab.degree>20]
        ax.plot(train.k,train.signed_PT,marker='o',label='Calibration')
        ax.scatter(test.k,test.signed_PT,marker='x',label='Held-out observed')
        ax.set_title(f'Period-corrected PT: s={v}, shell={j}');ax.set_xlabel('k');ax.legend();fig.tight_layout()
        mo.vstack([mo.as_html(fig),mo.ui.table(tab[['degree','k','signed_PT']])])
    else:
        mo.md('No rows for this layer.')
    return

@app.cell
def _(mo,preds):
    mo.md('## Forward holdout errors by shell')
    mo.ui.table(preds[['degree','shell','method','predicted_PT','true_PT','exact','relative_error']])
    return

@app.cell
def _(pd,ROOT,mo):
    fingerprints=pd.read_csv(ROOT/'results/dt_fingerprint_comparison.csv')
    mo.vstack([mo.md('## Finite differences versus independently listed DT coefficients'),mo.ui.table(fingerprints)])
    return

@app.cell
def _(ROOT,pd,mo):
    future=pd.read_csv(ROOT/'results/forward_Q31_Q40.csv')
    mo.vstack([mo.md('## Recurrence forecasts beyond degree 30 — unverified'),mo.ui.table(future)])
    return

if __name__=='__main__':app.run()
