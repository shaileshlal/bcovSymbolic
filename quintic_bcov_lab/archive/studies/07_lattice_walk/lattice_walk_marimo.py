"""Interactive read-only GV charge lattice exploration.

Run: marimo edit lattice_walk_marimo.py
Data: results/exact_lattice_audit.json
"""
import marimo
app=marimo.App(width='medium')

@app.cell
def _():
    import marimo as mo
    import json
    from pathlib import Path
    from fractions import Fraction
    path=Path(__file__).parent/'results'/'exact_lattice_audit.json'
    info=json.loads(path.read_text())
    mo.md(r"""# Quintic GV lattice walk

    We work only on the *homogeneous residual* ambiguity
    $\delta a_g=X(1-X)^{m_g}\sum_{r=0}^{D_g-1}c_{g,r}X^r$.
    The integer GV displacement $\Delta n$ maps to a rational coefficient
    displacement $\Delta c=H_g^{-1}\Delta n$.
    """)
    return mo,info,Fraction

@app.cell
def _(mo,info):
    genus=mo.ui.dropdown(options=[str(x['g']) for x in info['records']],value='9',label='Genus')
    degree=mo.ui.slider(start=1,stop=4,step=1,value=2,label='GV charge index (clipped to D)')
    steps=mo.ui.slider(start=-5,stop=5,step=1,value=1,label='Number of GV lattice steps')
    mo.vstack([genus,degree,steps])
    return genus,degree,steps

@app.cell
def _(genus,degree,steps,info,Fraction,mo):
    record=next(x for x in info['records'] if x['g']==int(genus.value))
    D=record['D']
    d=min(D,int(degree.value))
    H=[[Fraction(x) for x in row] for row in record['H']]
    c=[]
    delta_n=[Fraction(steps.value if j==d-1 else 0) for j in range(D)]
    for j in range(D):
        c.append((delta_n[j]-sum((H[j][k]*c[k] for k in range(j)),Fraction(0)))/H[j][j])
    check=[sum((H[j][k]*c[k] for k in range(D)),Fraction(0)) for j in range(D)]
    assert check==delta_n
    monomials=[f'({(-1)**(record["g"]-2-r)*v})x^{r+1}(1-x)^{D-1-r}' for r,v in enumerate(c) if v]
    xm=' + '.join(monomials) if monomials else '0'
    midpoint=sum(((-1)**(record['g']-2-r)*v for r,v in enumerate(c)),Fraction(0))/2**D
    mo.md(fr"""## Genus {record['g']} · D={D}

    **Unit charge change:** $\Delta n_{{g,{d}}}={steps.value}$

    **Exact residual coefficients:** `{[str(x) for x in c]}`

    **Section:** `{xm}`

    **Midpoint shift:** `{midpoint}`

    Independent affine particular verified: **{record['independent_particular_matches_reference_split']}**.

    All section and low-degree GV consistency checks are satisfied by construction.
    Nonzero shifts are not necessarily the physical ambiguity.
    """)
    return

if __name__=='__main__':
    app.run()
