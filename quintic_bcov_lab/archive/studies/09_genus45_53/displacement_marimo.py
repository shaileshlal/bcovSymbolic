# /// script
# requires-python = ">=3.10"
# dependencies = ["marimo", "pandas", "matplotlib"]
# ///
"""Genus 45–53 bootstrap and high-genus modular-phase holdouts.

Launch with `marimo edit displacement_marimo.py`.
"""
import marimo
app = marimo.App(width="full")

@app.cell
def _():
    import marimo as mo
    import pandas as pd
    from pathlib import Path
    import matplotlib.pyplot as plt
    base=Path(__file__).resolve().parent
    return mo,pd,base,plt

@app.cell
def _(mo):
    mo.md(r"""
    # GV lattice walker: exact genus 45–53 and modular near-boundary holdouts

    All arithmetic certificates are supplied in the accompanying JSON files.
    The genus 45–53 entries are zero by Castelnuovo except for the known
    geometric exceptional charge $n_{51,20}=175$. The eight high-genus
    holdout labels are **modularity conditional**, not independently proved.
    """)
    return

@app.cell
def _(mo):
    genus = mo.ui.slider(45,53,value=51, label="Genus")
    residue = mo.ui.dropdown(options=[1,2,3,4], value=1, label="Degree mod 5")
    mo.hstack([genus,residue],justify="start")
    return genus,residue

@app.cell
def _(pd,base,genus,mo):
    frame=pd.read_csv(base/'data'/'sequential_g45_g53.csv')
    subset=frame.loc[frame['genus']==genus.value].copy()
    mo.md(f"**Genus {genus.value}: {len(subset)} residual charges**")
    mo.ui.table(subset[['degree','predicted_GV_charge','source']])
    return frame

@app.cell
def _(pd,base,residue,mo):
    rf=pd.read_csv(base/'data'/'residue_family_benchmark.csv')
    one=rf.loc[rf['residue']==residue.value].copy()
    mo.md(f"**Residue {residue.value}: first three training points, next two strict holdouts**")
    mo.ui.table(one[['degree','genus','data_split','reference_integer','model_prediction','naive_signed_quadratic']])
    return one

@app.cell
def _(one,plt,mo):
    fig,ax=plt.subplots(figsize=(7,3.2))
    train=one.loc[one['data_split']=='train']
    hold=one.loc[one['data_split']=='holdout']
    ax.plot(one.k,one.reference_integer,ls=':',label='Published/reference')
    ax.scatter(train.k,train.reference_integer,label='Training',marker='o')
    ax.scatter(hold.k,hold.model_prediction,label='Predicted holdout',marker='x',s=65)
    ax.set_xlabel('k in d=5k+s');ax.set_ylabel('GV lattice integer');ax.legend()
    ax.grid(True,alpha=0.25)
    mo.as_html(fig)
    return

if __name__=='__main__':
    app.run()
