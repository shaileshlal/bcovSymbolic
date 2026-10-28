# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo>=0.13", "sympy>=1.12", "numpy>=1.24", "pandas>=2.0", "matplotlib>=3.7"]
# ///
"""Exact quintic DT->PT->GV arithmetic explorer. Open with marimo edit dt_pt_marimo.py."""
import marimo
app=marimo.App(width="full")

@app.cell(hide_code=True)
def _():
    import marimo as mo
    mo.md(r"""
    # Empirical DT/PT extrapolation of the quintic holomorphic ambiguity

    The seven polar rank-zero DT coefficients are the only input to the exact
    modular-form fit. No published degree-21 PT or GV values enter the fit.
    The genus-54 integer is predicted from the rank-zero DT polar coefficient
    **8625** and a single-term wall-crossing factor **37**.

    **Caveat:** This certifies a DT/PT/GV conditional prediction and the BCOV
    homogeneous ambiguity response, not the complete genus-54 B-model amplitude.
    """)
    return (mo,)

@app.cell
def _():
    from pathlib import Path
    import pandas as pd
    import matplotlib.pyplot as plt
    from empirical_benchmark import site,learned,ROWS,FORECAST
    from pt_from_dt_solver import simple_vacuum_pt,pt_from_rankzero_dt
    from modular_core import h_series
    return (Path,pd,plt,site,learned,ROWS,FORECAST,simple_vacuum_pt,pt_from_rankzero_dt,h_series)

@app.cell
def _(mo):
    degree=mo.ui.slider(start=6,stop=40,step=1,value=21,label="Curve degree Q")
    degree
    return (degree,)

@app.cell(hide_code=True)
def _(mo,degree,site,learned):
    record=site(degree.value)
    guess=record['kinematic_multiplier']*learned[record['rankzero_mu']]
    condition=record['wall_inequality_pass'] and not record['nonvacuum_candidate_count']
    mo.md(f"""
    ## Prediction at degree {degree.value}

    - Genus: **{record['genus']}**; charge residue: **{record['rankzero_mu']}**.
    - Signed kinematic multiplier: **{record['kinematic_multiplier']}**.
    - Learned small DT coefficient: **{learned[record['rankzero_mu']]}**.
    - Predicted leading PT/GV: **{guess}**.
    - Single-wall theorem conditions certified by enumeration: **{condition}**.

    A false certification flag means that the *vacuum-only* formula is a
    conjectural extrapolation, not that the integer is necessarily incorrect.
    """)
    return ()

@app.cell(hide_code=True)
def _(mo,pd,ROWS):
    dat=pd.DataFrame([{"Q":a['degree'],"genus":a['genus'],
                       "predicted_PT":a['prediction_PT'],"reference_PT":a['reference_PT_leading'],
                       "split":a['split'],"nonvacuum_possible":a['nonvacuum_summands_possible']}
                      for a in ROWS])
    mo.vstack([mo.md("## Degree-ordered evaluation"),mo.ui.table(dat)])
    return ()

@app.cell
def _(mo,plt,ROWS):
    fig,ax=plt.subplots(figsize=(9,3.3))
    for phase in ["calibration","forward_pre_54","frontier_holdout"]:
        rs=[r for r in ROWS if r['split']==phase]
        ax.scatter([r['degree'] for r in rs], [abs(r['reference_PT_leading']) for r in rs],label=phase)
    ax.set_yscale('log')
    ax.set_xlabel('Q');ax.set_ylabel('absolute leading PT integer');ax.legend()
    ax.grid(True,alpha=.2)
    mo.as_html(fig)
    return ()

@app.cell(hide_code=True)
def _(mo,pd,pt_from_rankzero_dt):
    rows=[]
    for m in range(-53,-48):
        value,vacuum,walls=pt_from_rankzero_dt(21,m)
        rows.append({'m':m,'PT from DT':str(value),'single-term value':str(vacuum),
                     'nonvacuum correction':str(value-vacuum),'candidate count':len(walls)})
    mo.vstack([mo.md('## Degree-21 five-term DT → PT reconstruction'),mo.ui.table(pd.DataFrame(rows))])
    return ()

@app.cell(hide_code=True)
def _(mo,pd,FORECAST):
    mo.vstack([mo.md('## Out-of-table forecasts (not independently validated)'),
               mo.ui.table(pd.DataFrame([{'Q':x['degree'],'genus':x['genus'],
                          'integer':str(x['predicted_PT_leading']),
                          'one-wall exact':x['theorem_single_wall_certified']} for x in FORECAST]))])
    return ()

if __name__=='__main__':
    app.run()
