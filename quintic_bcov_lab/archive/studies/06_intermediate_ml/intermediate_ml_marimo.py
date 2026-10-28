# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo>=0.13", "numpy>=1.24", "pandas>=2.0", "scikit-learn>=1.2", "matplotlib>=3.7"]
# ///
"""Exploratory marimo notebook for the *free* quintic YY ambiguity coefficients.

Executed numerical results are also exported by `python run_experiment.py`, which
requires no marimo runtime. All predictions are conditional on the source anomaly
at the target genus (it uses previously fixed lower-genus amplitudes).
"""
import marimo
app = marimo.App(width="full")

@app.cell(hide_code=True)
def _():
    import marimo as mo
    mo.md(r"""
    # Quintic BCOV — restricted intermediate ML

    We only predict the physical coordinates of the homogeneous boundary kernel
    \(\delta a_g=X(1-X)^{m_g}\sum_{r=0}^{D_g-1}c_{g,r}X^r\), equivalently
    \(\delta f_g=\sum_r(-1)^{g-2-r}c_{g,r}x^{r+1}(1-x)^{D_g-1-r}\).

    There are **22** exact rational coefficient labels for \(4\le g\le12\).
    The target-genus physical \(a_{g,j}\), reference-derived zero-pivot
    anchor and midpoint values are **never input features**. The mixed
    nonholomorphic source assumes already-fixed lower genera.

    The main experiment separates **logarithmic magnitude** from **sign**.
    Each row retains the exact rational label, but regression is floating point.
    """)
    return (mo,)

@app.cell
def _():
    import run_experiment as experiment
    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    rows=experiment.load_rows()
    return (experiment,np,plt,pd,rows)

@app.cell
def _(mo):
    split=mo.ui.dropdown(options={"Train 4–8, hold out 9–12":8,"Train 4–9, hold out 10–12":9},
        value="Train 4–8, hold out 9–12",label="Forward-genus holdout")
    split
    return (split,)

@app.cell
def _(mo):
    family=mo.ui.dropdown(options={
       "Genus + residual index":"genus",
       "Mixed source scale + residual index":"mixed",
       "Mixed source scale + genus + residual index":"mixed_plus_genus",
       "Fixed mixed-source slope + genus + residual index":"mixed_residual",
       "Original v0.1 signed-log ExtraTrees":"v01_trees"
       },value="Mixed source scale + genus + residual index",label="Magnitude model")
    family
    return (family,)

@app.cell
def _(experiment,rows,split,family):
    sign_rule="combined_v01" if family.value=="v01_trees" else "majority"
    predicted,params=experiment.run(rows,split.value,family.value,sign_rule)
    metrics=experiment.summarize_predictions(predicted)
    return (predicted,params,metrics,sign_rule)

@app.cell(hide_code=True)
def _(mo,metrics,sign_rule):
    mo.md(rf"""
    ## Out-of-genus evaluation

    - Held-out coefficients: **{metrics['n_coefficients']}** in **{metrics['n_genera']}** genera.
    - Absolute error in \(\log_{{10}}|c|\): **{metrics['MAE_log10_magnitude']:.3f} dex** mean,
      **{metrics['median_abs_log10_magnitude']:.3f} dex** median.
    - Typical multiplicative magnitude factor: **{metrics['median_magnitude_factor']:.2f}**.
    - Coefficient-sign accuracy: **{100*metrics['sign_accuracy_by_coefficient']:.1f}%**
      (only {metrics['n_genera']} independent genus signs; sign policy `{sign_rule}`).
    - Newly appearing residual positions in holdout: **{metrics['n_unseen_position']}**.

    Errors in sign are **not** included in magnitude scores; report both.
    The theoretical boundary kernel is fixed, but this experiment does not
    recover any nontrivial exact rational label.
    """)
    return

@app.cell
def _(mo,pd,predicted):
    display=pd.DataFrame(predicted)[['genus','r','actual_sign','pred_sign','actual_log_mag','pred_log_mag','error_dex','true_exact_c']]
    mo.ui.table(display,label="Predictions against exactly known but held-out residual coefficients")
    return (display,)

@app.cell
def _(mo,np,plt,predicted):
    idx=np.arange(len(predicted))
    real=np.array([a['actual_log_mag'] for a in predicted]);p=np.array([a['pred_log_mag'] for a in predicted])
    fig,ax=plt.subplots(figsize=(8,3.6))
    ax.plot(idx,real,'o-',label='Exact reference log magnitude')
    ax.plot(idx,p,'x--',label='Predicted log magnitude')
    ax.set_xticks(idx,[f"{a['genus']}.{a['r']}" for a in predicted],rotation=60)
    ax.set_xlabel('genus . residual index');ax.set_ylabel(r'$\log_{10}|c_{g,r}|$')
    ax.set_title('Chronological genus holdout — magnitudes only')
    ax.legend();fig.tight_layout()
    mo.mpl.interactive(fig)
    return

@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Source normalization and within-genus shape

    Denote by \(M_{g,r}\) the base-ten logarithm of the largest absolute
    mixed-YY coefficient at \(X^{m_g+1+r}\). A deliberately small model is
    \[
      \widehat{\log_{10}|c_{g,r}|}
      =\alpha+\beta M_{g,r}+\gamma g+\delta r.
    \]

    There is a separate diagnostic: when the **true** \(c_{g,0}\) is supplied,
    approximate \(\log_{10}|c_{g,r}/c_{g,0}|\approx\rho r\).
    This isolates the shape but is *not* blind prediction. Signs are
    empirical and genus-wide in the current archive; no reliable general
    sign law has been established.
    """)
    return

@app.cell
def _(mo,experiment,rows,split,pd,np,plt):
    import json
    from pathlib import Path
    report=json.loads((Path(experiment.__file__).parent/'results'/'benchmark.json').read_text())
    sh=report['shape_tests'][str(split.value)]
    p=sh['predictions']
    fig,ax=plt.subplots(figsize=(7.4,3.1))
    ax.scatter([z['r'] for z in p],[z['true_shape_log10'] for z in p],label='Reference shapes')
    vals=np.arange(1,5);ax.plot(vals,sh['slope_per_index']*vals,label=f"Fitted slope {sh['slope_per_index']:.3f}")
    ax.set_xlabel('r');ax.set_ylabel(r'$\log_{10}|c_{g,r}/c_{g,0}|$')
    ax.set_title('Shape-only diagnostic: true first coefficient supplied')
    ax.legend();fig.tight_layout()
    mo.mpl.interactive(fig)
    return

@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Interpretation

    - This dataset is **nine nontrivial genera**, not 22 independent genus observations.
    - The arithmetic values are retained exactly as rational strings; the ML stage
      estimates decimal logarithms and never certifies a rational reconstruction.
    - High-rank sinusoidal/GV-phase attention would be underdetermined here:
      the first training split contains only residual positions r=0 and r=1.
    - Use `python run_experiment.py` for reproducible CSV predictions, and
      `python tests.py` for the exact data and no-leakage checks.
    """)
    return

if __name__=="__main__":
    app.run()
