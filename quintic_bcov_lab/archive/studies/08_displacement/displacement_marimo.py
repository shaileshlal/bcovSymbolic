# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo>=0.13", "matplotlib>=3.7", "numpy>=1.24"]
# ///
"""Quintic BCOV residual-GV integer displacement: EDA and exact frontier walker."""
import marimo

__generated_with = "0.19.7"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    mo.md(r"""
    # Quintic BCOV: integer displacements from minimum height

    In the **unrestricted integral GV lattice**, the minimum-height charge is
    zero. The target for a genuine residual site is therefore the exact
    integer \(\Delta n_{g,d}=n_{g,d}\), with \(d=r+1\).

    **Three labels, not two:** proven zeros, proven extremal nonzeros, and
    **unresolved** (never replace these by zero in a training dataset).

    The residual section is
    \[
      \delta f_g(x)=\sum_{r=0}^{D_g-1}(-1)^{g-2-r}c_{g,r}
        x^{r+1}(1-x)^{D_g-r-1},\qquad
      D_g=\left\lfloor\frac{2g-2}{5}\right\rfloor.
    \]
    """)
    return (mo,)


@app.cell
def _():
    import sys, csv, json
    from pathlib import Path
    root=Path(__file__).resolve().parent
    if str(root) not in sys.path:sys.path.insert(0,str(root))
    from build_data import residual_dimension, margin10
    from frontier_solver import classify_charge, first_nonzero_section_taylor
    manifest=json.loads((root/'data'/'dataset_manifest.json').read_text())
    model=json.loads((root/'results'/'model_benchmarks.json').read_text())
    top=list(csv.DictReader((root/'data'/'extremal_displacements.csv').open()))
    return (root, manifest, model, top, residual_dimension, margin10,
            classify_charge, first_nonzero_section_taylor)


@app.cell(hide_code=True)
def _(mo, manifest):
    mo.md(f"""
    ## Label audit through genus 1000

    - **{manifest['total']:,}** residual lattice positions
    - **{manifest['certified_zero']:,}** known zeros by geometric vanishing
    - **{manifest['extremal_nonzero']}** known nonzero extremal values
    - **{manifest['unresolved']:,} unresolved** sites with no provided integer label
    - **{manifest['reference_calibration']}** exact genus-4–12 BCOV residual checks

    The extremal nonzero training set has only eight genera; eight later
    extremal genera are held out. The source theorem, not BCOV data at
    those genera, provides the training and reference integers.
    """)
    return


@app.cell
def _(mo):
    genus=mo.ui.slider(start=49,stop=120,value=54,step=1,label="Explore frontier genus")
    genus
    return (genus,)


@app.cell
def _(mo, genus, residual_dimension, classify_charge):
    g=genus.value
    D=residual_dimension(g)
    inspect=[classify_charge(g,d) for d in range(max(1,D-9),D+1)]
    mo.md(f"### Genus {g}: {D} residual lattice sites")
    mo.ui.table(inspect,label="Last ten residual charges — unresolved labels are not zero")
    return (g,D,inspect)


@app.cell(hide_code=True)
def _(mo, g, D, classify_charge):
    un=[d for d in range(1,D+1) if classify_charge(g,d)['status']=='ABSTAIN_unresolved']
    mo.md(f"""
    **Unresolved charges at genus {g}:** {len(un)}

    First unresolved curve degree: **{un[0] if un else 'none'}**.
    Exact vanishing and extremal data are always applied before any
    learned prediction. An unresolved charge is an explicit abstention.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Periodic extremal family

    For \(d=5k\) at \(g=1+\tfrac52 k(k+1)\), the extremal theorem gives
    a small signed integer. The model-search experiment rediscovers
    \[
       \Delta n_k=25\left(1+\binom{k}{2}\right)
          \left[\cos\left(\frac{\pi k}{2}\right)
               +\sin\left(\frac{\pi k}{2}\right)\right].
    \]
    The bracket equals \((+,+,-,-)\) as \(k\bmod4=(0,1,2,3)\).
    **This is a learned re-expression of a known theorem**, not a new
    evaluation of previously unknown quintic invariants.
    """)
    return


@app.cell
def _(mo, top):
    mo.ui.table([
        {"k":int(t['k']),"g":int(t['g']),"d":int(t['d']),
         "exact_charge":int(t['n']),"split":t['split']}
        for t in top],label="Known extremal integer labels and honest holdouts")
    return


@app.cell
def _(mo, model):
    metrics=model['metrics']
    rows=[{"model":k,"mean absolute integer error":round(v['mae'],7),
           "exact test integers":f"{v['exact_integers']}/8",
           "sign accuracy":round(v['sign_accuracy'],3)} for k,v in metrics.items()]
    mo.ui.table(rows,label="Holdout k=12..19; training k=4..11")
    return


@app.cell
def _(root):
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    img=mpimg.imread(root/'figures'/'geometric_boundary.png')
    fig,ax=plt.subplots(figsize=(9,5))
    ax.imshow(img);ax.axis('off')
    fig
    return


@app.cell
def _(root):
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    img=mpimg.imread(root/'figures'/'extremal_pattern.png')
    fig,ax=plt.subplots(figsize=(9,5))
    ax.imshow(img);ax.axis('off')
    fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Exact lattice-to-section translation

    The previously computed BCOV response is triangular:
    \[
      \Delta n_{g,d}=\sum_rH_{g,d,r}\Delta c_{g,r},\quad
      H_{g,d,d-1}=(-1)^{g-1+d}5^{5d-4g+4}.
    \]
    If the lower charges vanish, the first affected term in the reduced section obeys
    \[
     [x^d]\Delta f_g=5^{4g-4-5d}\Delta n_{g,d}.
    \]
    If \(d=D_g\) is the last free position, there are **no additional terms**.

    Hence \(\Delta f_{51}=175\,5^{100}x^{20}\). At genus 54 the
    result is \(\Delta f_{54}=5^{107} n_{54,21}x^{21}\), with the **integer unknown**.
    """)
    return


@app.cell
def _(mo, classify_charge, first_nonzero_section_taylor):
    places=[(51,20),(52,20),(53,20),(54,21),(55,21),(76,25)]
    records=[]
    for g,d in places:
        p=classify_charge(g,d)
        if p.get('delta_n') not in (None,0):
            p['section_response']=first_nonzero_section_taylor(g,d,p['delta_n'])['coefficient_factor']
        records.append(p)
    mo.ui.table(records,label="Exact physical inputs versus unresolved frontier")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Interpretation and limitations

    The eight exact held-out extremal predictions show that a small Fourier
    phase plus polynomial envelope can encode known arithmetic perfectly.
    They **do not** test the unknown interior sector or independently obtain
    \(n_{54,21}\).

    A large classifier trained on 63,544 zeros and just 16 top-genus labels
    would mostly learn the Castelnuovo support mask. It is not evidence of
    generalization to the missing 135,840 entries. The production policy
    therefore **abstains** at unresolved lattice points until independent
    BPS/modular data are added.
    """)
    return


if __name__ == "__main__":
    app.run()
