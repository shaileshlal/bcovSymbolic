"""Sequential overview of the quintic BCOV project, as a marimo notebook.

Run: marimo edit notebooks/overview_marimo.py
Only reads stage metadata, exact source data and recorded results, never
promotes withheld truth to a predictor feature.
"""
import marimo
app=marimo.App(width="full")

@app.cell

def _():
    import marimo as mo
    import json
    from pathlib import Path
    from fractions import Fraction
    import pandas as pd
    from quintic_bcov_lab import datasets, runner
    return mo,json,Path,Fraction,pd,datasets,runner

@app.cell

def _(mo,runner):
    mo.md("""# Sequential BCOV / PT / GV arithmetic research

    Every stage is preserved separately. Green checkmarks in historical reports
    reflect exact comparisons **within the conditional reference framework**;
    they are not automatic unconditional enumerative proofs.
    """)
    stage = mo.ui.dropdown(options=[s['id'] for s in runner.stages()],value='05_intermediate_dataset',label='Research stage')
    stage
    return (stage,)

@app.cell

def _(stage,runner,mo,Path):
    z=runner.resolve_stage(stage.value)
    root=runner.ROOT/'studies'/z['id']
    mo.md(f"**{z['title']}** — archived version {z['version']}\n\n"
          f"Run: `python -m quintic_bcov_lab stage {z['id']} --action run`  \n"
          f"Test: `python -m quintic_bcov_lab stage {z['id']} --action test`")
    return (root,)

@app.cell

def _(mo,datasets,pd):
    samples=datasets.intermediate_labels()
    df=pd.DataFrame(samples)
    mo.md('### Exact intermediate coefficients (reference-derived labels)')
    mo.ui.table(df)
    return (df,)

@app.cell

def _(mo,root):
    notebooks=sorted(p.name for p in root.rglob('*marimo.py'))
    mo.md('### Original notebook(s) in this stage\n'+('\n'.join('- `'+p+'`' for p in notebooks) if notebooks else 'No original marimo notebook'))
    return ()

if __name__=='__main__':
    app.run()
