"""Empirical quintic PT pulses at shell depths 11--15.

Run `marimo edit depth_marimo.py`. References are read ONLY to display
independent validation, never by the predictive model.
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
    import matplotlib.pyplot as plt
    root=Path(__file__).resolve().parent
    return mo,root,F,csv,json,plt


@app.cell
def _(mo):
    mo.md(r"""
    # Quintic PT: the next five defect shells

    **Shell depth 11--15, empirical v1.1.**

    For \(d=5k+s\), define
    \[
      W_{s,t}(k)=5k^2+(2s-5)k+c_s-2t,
      \quad A_{s,t}(k)=\frac{\sigma_s(k)\,\mathrm{PT}(5k+s,t)}{W_{s,t}(k)}.
    \]

    The working residual is \(A-B=(-1)^t(P_{s,r}+Q_{s,r}/W)\),
    with \(r=K_s(t)-k\). The learned paired-residue model only emits
    an exact integer when its calibration equations uniquely fix the answer.

    **Data provenance:** prediction input uses curve degree at most 40.
    The comparison/reference is constructed separately from a conditional
    DT/PT computation; this notebook uses it for visualization only.
    """)
    return


@app.cell
def _(mo):
    depth=mo.ui.slider(start=11,stop=15,step=1,value=11,label="Shell depth")
    degree=mo.ui.slider(start=41,stop=53,step=1,value=45,label="Curve degree")
    mo.hstack([depth,degree],gap=3)
    return depth,degree


@app.cell
def _(root,csv,json):
    def readcsv(path):
        with path.open() as f:
            return list(csv.DictReader(f))
    preds=readcsv(root/'results/linear_profile_pairs_predictions.csv')
    reference=readcsv(root/'validator/reference_deep_11_15.csv')
    metrics=json.loads((root/'results/benchmark.json').read_text())
    anchors=json.loads((root/'results/anchor_rank.json').read_text())
    closure=json.loads((root/'results/closure_tests.json').read_text())
    return preds,reference,metrics,anchors,closure


@app.cell
def _(mo,depth,degree,preds,reference):
    t=int(depth.value)
    d=int(degree.value)
    p=next(z for z in preds if int(z['degree'])==d and int(z['shell'])==t)
    q=next(z for z in reference if int(z['degree'])==d and int(z['shell'])==t)
    mo.md(f"""### Selected site $(d,t)=({d},{t})$

    - **Front distance:** `{p['r']}`
    - **Blind prediction:** `{p['prediction'] or 'abstain'}`
    - **Conditional reference:** `{q['pt'] if q['pt'] else q['status']}`
    - **Prediction status:** `{'exact' if p['prediction'] and p['prediction']==q['pt'] else 'unidentified or outside model'}`
    """)
    return t,d


@app.cell
def _(mo,t,preds,reference):
    rmap={(int(z['degree']),int(z['shell'])):z for z in reference}
    out=[]
    for p in preds:
        if int(p['shell'])!=t:continue
        d=int(p['degree'])
        q=rmap[d,t]
        out.append(dict(degree=d,r=p['r'],predicted_PT=p['prediction'] or 'abstain',reference_PT=q['pt'],
                        exact=bool(p['prediction'] and p['prediction']==q['pt'])))
    mo.ui.table(out,page_size=13)
    return out


@app.cell
def _(plt,mo,metrics,anchors):
    x=list(range(11,16))
    y=[metrics['prediction']['by_depth'][str(t)]['exact'] for t in x]
    fig,ax=plt.subplots(figsize=(7.1,3.2))
    ax.bar(x,y)
    ax.set(xlabel='Shell depth',ylabel='Exact PT integers',title='Prediction coverage on degrees 41--53',ylim=(0,13.5))
    ax.set_xticks(x)
    mo.as_html(fig)
    return


@app.cell
def _(mo,anchors,closure):
    rows=[{'depth':k,'eligible':v['eligible'],'new_integer_anchors':v['new_independent_anchors']} for k,v in anchors['by_depth'].items()]
    mo.md('### Minimum number of new integer anchors (conditional affine model)')
    mo.ui.table(rows)
    mo.md('The three scale streams have Hankel rank 8 at size 8, with no tested constant-coefficient recurrence of order at most 7. This does **not** rule out finite nonlinear generation.')
    return


if __name__ == "__main__":
    app.run()
