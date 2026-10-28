"""Marimo explorer for the empirical sparse-defect model.

Run: pip install marimo pandas matplotlib scikit-learn
     marimo edit defect_marimo.py
No modular/PT-DT/GV identity is loaded by the predictor.  Full PT labels in
this GUI are for after-the-fact diagnostics only.
"""
from __future__ import annotations
import marimo
app=marimo.App(width="full")

@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import matplotlib.pyplot as plt
    import json
    from pathlib import Path
    ROOT=Path(__file__).resolve().parent
    return mo,pd,plt,json,ROOT

@app.cell
def _(mo):
    mo.md(r"""
    # Sparse third-difference defects in quintic PT invariants

    The first predictor uses only degrees **1–20**. The second block treats
    degrees **21–25** as a sequential update before forecasting 26–30.
    The labels were published in a modularity-conditional reference table,
    but modularity, wall-crossing and exact DT relations are not used as
    predictive features.

    \[T_{s,t}(k)=\sigma_s(k)\operatorname{PT}(5k+s,m_C+\delta_s+t).\]
    \[C_{s,t}(k)=\Delta^2 T_{s,t}(k),\qquad\varepsilon_{s,t}(k)=\Delta^3 T_{s,t}(k).\]
    """)
    return

@app.cell
def _(mo):
    s=mo.ui.dropdown(options={f"Residue {i}":i for i in range(5)},value="Residue 1",label="Degree mod 5")
    t=mo.ui.slider(start=0,stop=5,value=2,step=1,label="Depth t")
    mo.hstack([s,t])
    return s,t

@app.cell
def _(pd,ROOT):
    observations=pd.read_csv(ROOT/'results/defect_analysis_only.csv',dtype={'third_difference':str})
    predictions=pd.concat([pd.read_csv(ROOT/'results/block_Q21_Q25.csv',dtype=str),pd.read_csv(ROOT/'results/block_Q26_Q30.csv',dtype=str)],ignore_index=True)
    scores=pd.read_csv(ROOT/'results/event_risk_predictions_Q26_Q30.csv')
    future=pd.read_csv(ROOT/'results/future_Q31_Q35.csv',dtype=str)
    return observations,predictions,scores,future

@app.cell
def _(s,t,mo,observations,plt,pd):
    v=s.value;j=t.value
    rows=observations[(observations.residue==v)&(observations.shell==j)].copy()
    if rows.empty:
        mo.md('No published PT data for this residue–shell selection.')
    else:
        fig,ax=plt.subplots(figsize=(7.8,3.5))
        for r in rows.itertuples():
            y=int(r.third_difference)
            ax.scatter(int(r.degree), (1 if y>0 else -1 if y<0 else 0)*__import__('math').log10(1+abs(y)),label=str(r.degree))
        ax.axhline(0,linewidth=.8,color='gray')
        ax.set(xlabel='Curve degree',ylabel='Signed log10(1 + |third difference|)',title=f'Empirical defect sequence · residue {v}, shell {j}')
        fig.tight_layout()
        mo.vstack([mo.as_html(fig),mo.ui.table(rows[['degree','shell','third_difference','active','v5_delta']])])
    return

@app.cell
def _(s,t,mo,predictions,pd):
    v=s.value;j=t.value
    p=predictions[(predictions.residue.astype(int)==v)&(predictions.shell.astype(int)==j)].copy()
    mo.vstack([mo.md('## Causal PT forecasts'),mo.ui.table(p[['degree','shell','predicted_PT','model']])])
    return

@app.cell
def _(mo,scores):
    mo.vstack([mo.md('## Predicted next-block defect activity'),mo.ui.table(scores[['degree','shell','logit_positional_p','logit_engineered_p','unit_speed_front_p','stable_absorbing_p']])])
    return

@app.cell
def _(mo,future):
    mo.vstack([mo.md('## Candidate next-degree predictions (31–35, not independently validated)'),mo.ui.table(future[['degree','shell','predicted_PT','model','verification_status']])])
    return

if __name__=='__main__':app.run()
