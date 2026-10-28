# /// script
# requires-python = ">=3.11"
# dependencies = ["marimo>=0.13", "matplotlib>=3.7"]
# ///
"""Reactive quintic YY residual ambiguity dataset, exact arithmetic."""
import marimo

__generated_with = "0.19.7"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    mo.md(r"""
    # Quintic BCOV: the *intermediate* holomorphic ambiguity

    **Exact YY coefficients, no neural predictions.** This notebook rebuilds an
    intermediate-only dataset from canonical corrected genus 2–12 amplitudes.
    The homogeneous freedom after conifold gap, orbifold regularity, and the MUM
    constant map is

    \[
      \delta a_g(X)=X(1-X)^{m_g}R_g(X),\quad
      m_g=\Big\lceil\frac{3g-3}{5}\Big\rceil,\quad
      \deg R_g < D_g=\Big\lfloor\frac{2g-2}{5}\Big\rfloor.
    \]

    We use \(x=\psi^{-5}=3125z_{\rm PF}\), \(X=x/(x-1)\), and the
    reduced **section coefficient** \(f_g(x)=(x-1)^{g-1}a_g(X(x))\).
    The common Hodge-period factor is suppressed, not set to one globally.
    """)
    return (mo,)


@app.cell
def _():
    from pathlib import Path
    import json
    from fractions import Fraction
    import sys
    root = Path(__file__).resolve().parent
    if str(root) not in sys.path:
        sys.path.insert(0, str(root))
    from intermediate import build, dimensions, reference_split, kernel_vector, bcov02_manifest
    source_path = root / "data" / "reference_amplitudes.json"
    output_path = root / "data"
    # The actual exact-basis transformation, performed when the notebook runs.
    data, training_rows, audit = build(source_path, output_path)
    genus_data = {d["genus"]: d for d in data}
    return (Path, Fraction, json, root, source_path, output_path, data, training_rows, audit, genus_data, dimensions, reference_split, kernel_vector, bcov02_manifest)


@app.cell(hide_code=True)
def _(mo, data, training_rows, audit):
    mo.md(f"""
    ## Rebuilt dataset

    - **{len(data)} reference genera**: {data[0]['genus']}–{data[-1]['genus']}
    - **{len(training_rows)} intermediate coefficient targets** (each an exact rational)
    - **{sum(x['round_trip_exact'] for x in audit)} / {len(audit)} exact round-trip checks**

    The `anchor` in the dataset is obtained *from the reference polynomial*
    by setting its intermediate pivot entries to zero. This is a canonical
    gauge choice, **not** a new genus-by-genus independent boundary solution.
    Every training label has this explicit provenance.
    """)
    return


@app.cell
def _(mo, genus_data):
    genus = mo.ui.dropdown(
        options={f"g={g} (D={genus_data[g]['D']})": g for g in sorted(genus_data)},
        value="g=9 (D=3)",
        label="Select genus",
    )
    genus
    return (genus,)


@app.cell
def _(genus, genus_data):
    entry = genus_data[genus.value]
    gg = genus.value
    mm = entry['m']
    DD = entry['D']
    KK = entry['K']
    return (entry, gg, mm, DD, KK)


@app.cell(hide_code=True)
def _(mo, gg, mm, DD, KK, entry):
    mo.md(rf"""
    ### Genus {gg}: {DD} undetermined rational numbers

    Original ambiguity: **{KK+1}** YY coefficients, \(j=0,\ldots,{KK}\).
    Canonical residual basis: \(m={mm}\), \(D={DD}\), pivot degrees
    **{', '.join(str(j) for j in entry['pivot_indices']) or '(none)'}**.

    Define the reference-derived anchor \(a_g^\star(X)\) by setting exactly
    those pivot coefficients to zero; then the *physical reference* is

    \[
      a_{{{gg}}}^{{\rm ref}}(X)=a_{{{gg}}}^\star(X)
        +X(1-X)^{{{mm}}}\sum_{{r=0}}^{{{DD-1}}} c_{{{gg},r}}X^r.
    \]

    The un-fixed \(c_{{g,r}}\) are **labels obtained from the known reference**.
    An ML model should predict these, not the boundary-fixed coefficients.
    """)
    return


@app.cell(hide_code=True)
def _(mo, gg, DD):
    mo.md(rf"""
    ## Section reorganization: \(x(1-x)\) basis

    \[
      \delta f_{{{gg}}}(x)=(x-1)^{{{gg-1}}}
        \delta a_{{{gg}}}\!\left(\frac{{x}}{{x-1}}\right)
       =\sum_{{r=0}}^{{{DD-1}}} b_{{{gg},r}}x^{{r+1}}(1-x)^{{{DD}-r-1}},
    \]
    \[
      b_{{{gg},r}}=(-1)^{{{gg}-2-r}}c_{{{gg},r}},\qquad
      \delta f_{{{gg}}}(x)=xQ_{{{gg}}}(x),\quad \deg Q_{{{gg}}}\le {DD-1}.
    \]

    This is a change of basis for the **homogeneous boundary kernel**, not
    a proposed physical symmetry under \(x\leftrightarrow1-x\).
    """)
    return


@app.cell
def _(mo, entry, DD):
    coord_rows = [
        {"r": r,
         "pivot X degree": entry['pivot_indices'][r],
         "physical YY pivot a[g,j]": entry['physical'][entry['pivot_indices'][r]],
         "c[g,r] in X basis": entry['c'][r],
         "b[g,r] in x(1-x) basis": entry['bernstein_x'][r]}
        for r in range(DD)
    ]
    mo.ui.table(coord_rows, label="Exact intermediate coefficients (rational strings)") if coord_rows else mo.md("No residual ambiguity at this genus.")
    return (coord_rows,)


@app.cell
def _(mo, entry, gg, KK, DD):
    raw_rows = [
      {"j":j, "a_g,j (physical)":entry['physical'][j],
       "boundary anchor a*_g,j":entry['anchor'][j],
       "residual delta a_g,j":entry['residual'][j],
       "label status":("INTERMEDIATE" if j in entry['pivot_indices'] else "fixed in canonical anchor")}
      for j in range(KK+1)]
    mo.ui.table(raw_rows, label=f"Full genus-{gg} YY coefficients: original = anchor + residual")
    return (raw_rows,)


@app.cell(hide_code=True)
def _(mo, entry, gg, DD):
    Qpoly = " + ".join(
        f"({v})x^{i}" for i,v in enumerate(entry['xQ']) if v) or "0"
    mo.md(rf"""
    ### Reduced section coefficient as a polynomial in \(x\)

    The exact polynomial \(\delta f_{{{gg}}}(x)=xQ_{{{gg}}}(x)\) is

    ```text
    {Qpoly}
    ```

    The maximum possible degree is **{DD}**. All higher powers cancel
    exactly even though \(a_g(X)\) originally reaches degree \(3g-3\).
    """)
    return


@app.cell
def _(mo, entry, Fraction, gg, DD):
    import matplotlib.pyplot as plt
    if DD:
        # No conversion to floating point until the *normalized* exact values
        xs = [Fraction(k,100) for k in range(101)]
        vals = [sum((Fraction(v) * xx**i for i,v in enumerate(entry['xQ'])),Fraction(0)) for xx in xs]
        norm = max((abs(v) for v in vals), default=Fraction(1)) or Fraction(1)
        fig, axis = plt.subplots(figsize=(7.5, 3.1))
        axis.plot([float(x) for x in xs], [float(v/norm) for v in vals], label="physical residual / max|residual|")
        axis.axhline(0,linewidth=0.5,alpha=0.55)
        axis.set_xlabel("x (normalized conifold coordinate)")
        axis.set_ylabel("normalized residual section")
        axis.set_title(f"Genus {gg}: residual shape (not the full amplitude)")
        axis.legend(loc="best")
        mo.mpl.interactive(fig)
    else:
        mo.md("No residual shape to plot: D=0.")
    return


@app.cell(hide_code=True)
def _(mo, entry, gg, DD, Fraction, kernel_vector, reference_split):
    from fractions import Fraction as _F
    rebuilt = reference_split(gg,[_F(x) for x in entry['physical']])
    checks = {
      'full YY coefficient vector exact': all(str(x)==y for x,y in zip(rebuilt['physical'],entry['physical'])),
      'canonical zero-pivot anchor': all(x==0 for x in (rebuilt['anchor'][j] for j in rebuilt['pivot_indices'])),
      'two independent residual bases agree': [str(x) for x in rebuilt['xQ']]==entry['xQ'],
      'no residual terms above X^(g-1)':all(_F(t)==0 for t in rebuilt['residual'][gg:]),
      'no shift in the MUM constant-map coefficient': rebuilt['residual'][0]==0,
    }
    mo.md("### Exact round-trip audit\n\n"+"\n".join(f"- {'✅' if value else '❌'} {key}" for key,value in checks.items()))
    return (checks,)


@app.cell(hide_code=True)
def _(mo, training_rows, gg):
    mo.md("## ML-ready dataset: residual labels only")
    mo.ui.table([r for r in training_rows if r['genus']==gg],label="Available mixed-YY features and exact intermediate targets")
    return


@app.cell
def _(mo):
    v02_path = mo.ui.text(value="./bcov_dual_corpus_v02",label="Optional local bcov0.2 structural corpus directory",full_width=True)
    v02_path
    return (v02_path,)


@app.cell(hide_code=True)
def _(mo, Path, v02_path, bcov02_manifest):
    observed = bcov02_manifest(Path(v02_path.value))
    mo.md(f"""
    ### v0.2 integration status

    The `bcov0.2` structural-learning archive and the exact-amplitude BCOV
    engine serve different purposes. V0.2 tasks/registry entries can be
    used as **optional structural input features**, but do not supply
    certified physical ambiguity labels by themselves.

    - v0.2 directory found: **{observed['present']}**
    - v0.2 tasks detected: **{observed['tasks']}**
    - structure registries detected: **{observed['registries']}**

    The standalone notebook always builds its verified target labels from
    `data/reference_amplitudes.json`, which is included with this package.
    When you place the v0.2 corpus at the path above, its manifest is checked
    without silently changing numerical labels.
    """)
    return


@app.cell(hide_code=True)
def _(mo, output_path):
    mo.md(f"""
    ## Exports and reproducibility

    Every run regenerates:

    - `data/intermediate_exact.json`: all 11 genera, exact physical, anchor,
      residual and section-basis rational vectors.
    - `data/intermediate_training.csv`: one row **per free intermediate
      coefficient** (22 rows) with optional nonholomorphic features.
    - `data/validation.json`: dimension counts and round-trip certificate flags.

    Source: `data/reference_amplitudes.json`, corrected genus-four convention.
    Rebuild outside marimo with `python intermediate.py`.

    **Important:** since the anchor is defined using reference data, the model
    cannot treat it as an independently available inference-time feature.
    To make future *blind* predictions, replace it with a separately computed
    particular solution of the exact BCOV boundary equations.
    """)
    return


if __name__ == "__main__":
    app.run()
