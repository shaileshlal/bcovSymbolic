"""Marimo exploration of the frozen-v0.7 recursive survival audit (no training here).
Run: marimo edit survival_marimo.py
"""
import marimo
app=marimo.App(width="full")

@app.cell

def _():
    import marimo as mo
    import pandas as pd
    from pathlib import Path
    root=Path(__file__).resolve().parent
    predictions=pd.read_csv(root/'results/individual_predictions_scored.csv',dtype={'predicted_pt':'string','reference_pt':'string','integer_error':'string'})
    defects=pd.read_csv(root/'results/third_differences.csv',dtype={'defect':'string'})
    mo.md('# Quintic PT recursion — degree 31–53 survival\nA comparison of the **unchanged** v0.7 predictions against an independently computed conditional PT reference.')
    return mo,predictions,defects

@app.cell

def _(mo):
    mode=mo.ui.dropdown(options={'Sequential (previous true blocks revealed)':'teacher','Autonomous (only model feedback)':'autonomous'},value='Autonomous (only model feedback)',label='Prediction protocol')
    degree=mo.ui.slider(start=31,stop=53,value=31,label='Curve degree')
    mo.hstack([mode,degree])
    return mode,degree

@app.cell

def _(mo,predictions,mode,degree):
    subset=predictions[(predictions['mode']==mode.value)&(predictions['degree']==degree.value)].copy()
    mo.md(f'### Degree {degree.value}: {int(subset["exact"].sum())} / {len(subset)} exact PT integers')
    mo.ui.table(subset[['shell','predicted_pt','reference_pt','integer_error','exact','mirror_used']])
    return

@app.cell

def _(mo,defects):
    d=defects[(defects['degree']>=26)&(defects['degree']<=53)].copy()
    d=d[d['active']==1]
    mo.md('### Nonzero third-difference defects: the shrinking shell front')
    mo.ui.table(d[['degree','residue','depth','defect']])
    return

if __name__=='__main__':
    app.run()
