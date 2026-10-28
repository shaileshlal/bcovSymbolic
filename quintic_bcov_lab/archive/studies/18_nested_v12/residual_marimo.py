"""Quintic PT: nested finite-difference pulses.

Launch with `marimo edit residual_marimo.py` after `pip install marimo matplotlib`.
The predictor is independent of validation/.  Reference tables in this notebook
are displayed solely to audit/compare already-computed empirical predictions.
"""
import marimo

__generated_with = "0.17.0"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import csv
    import json
    import matplotlib.pyplot as plt
    from pathlib import Path
    from fractions import Fraction
    root = Path(__file__).resolve().parent
    def readcsv(name):
        with (root / name).open(newline="") as f:
            return list(csv.DictReader(f))
    return mo, csv, json, plt, root, Fraction, readcsv


@app.cell
def _(mo):
    mo.md(r"""
    # Quintic PT residuals: recursion within finite-difference pulses

    Let \(Z_{s,r}(t)=(-1)^t(\sigma_s(k)\mathrm{PT}(5k+s,t)-W_{s,t}(k)B_{\alpha,t})\).

    Its first *new* defect is \(\varepsilon_{s,r}(t)=\Delta_t^3 Z_{s,r}(t)\).
    The last active pulse lies at \(T_s(r)=2r-1\) for \(s=1\) and \(2r+3\) otherwise.
    The terminal amplitude \(\tau_s(r)=\varepsilon_{s,r}(T_s(r))\) obeys
    \(\Delta_r^3\tau_s=0\) over the observed ranges.

    **Prediction:** 13 PT-derived terminal pulse calibration integers, 18 exact
    later terminal pulse matches. **Limit:** mature scales are observed PT inputs
    at possibly higher degree; this is a conditional shape benchmark, not a blind
    complete PT or BCOV reconstruction.
    """)
    return


@app.cell
def _(mo):
    residue = mo.ui.slider(start=0, stop=4, step=1, value=2, label="Residue s")
    offset = mo.ui.slider(start=0, stop=3, step=1, value=0, label="Distance behind terminal front ℓ")
    mo.hstack([residue, offset], gap=3)
    return residue, offset


@app.cell
def _(readcsv, root, json):
    train = readcsv('data/terminal_pulse_calibration.csv')
    pred = readcsv('results/terminal_pulse_predictions.csv')
    ref = readcsv('validation/terminal_pulse_holdout_truth.csv')
    second_pred = readcsv('results/diagonal_pulse_predictions.csv')
    second_ref = readcsv('validation/diagonal_holdout_truth.csv')
    metrics = json.loads((root / 'results/benchmark.json').read_text())
    return train, pred, ref, second_pred, second_ref, metrics


@app.cell
def _(mo, residue, train, pred, ref):
    s = int(residue.value)
    refs = {(int(z['residue']), int(z['r'])): z for z in ref}
    table = []
    for z in train:
        if int(z['residue']) == s:
            table.append({'r': int(z['r']), 'terminal depth': int(z['t']), 'source': 'calibration',
                          'terminal pulse': z['terminal_pulse'], 'exact': 'input'})
    for z in pred:
        if int(z['residue']) == s:
            actual = refs[s,int(z['r'])]['terminal_pulse']
            table.append({'r': int(z['r']), 'terminal depth': int(z['t']), 'source': 'holdout',
                          'terminal pulse': z['terminal_pulse'], 'exact': z['terminal_pulse'] == actual})
    mo.ui.table(sorted(table,key=lambda q:q['r']),page_size=12)
    return table


@app.cell
def _(plt, mo, table):
    rows = sorted(table,key=lambda q:q['r'])
    xs = [z['r'] for z in rows]
    ys = [int(z['terminal pulse']) for z in rows]
    denom = max(1,max(abs(y) for y in ys))
    fig, ax = plt.subplots(figsize=(7,3.4))
    ax.plot(xs, [y / denom for y in ys], '-o')
    ax.set(xlabel='Front-distance r',ylabel='Terminal pulse / max amplitude',
           title='Quadratic terminal-pulse sequence')
    ax.grid(True,alpha=0.2)
    mo.as_html(fig)
    return


@app.cell
def _(mo, offset, residue, second_pred, second_ref):
    s,ell=int(residue.value),int(offset.value)
    truth={(int(z['residue']),int(z['offset']),int(z['r'])):z['epsilon'] for z in second_ref}
    table=[]
    for z in second_pred:
        if int(z['residue'])==s and int(z['offset'])==ell:
            key=s,ell,int(z['r'])
            actual=truth.get(key)
            table.append({'r':key[2],'predicted epsilon':z['pulse_prediction'],
                          'reference epsilon':actual or 'not certified',
                          'exact':actual is not None and actual==z['pulse_prediction']})
    mo.md(f'### Secondary pulse at residue {s}, offset {ell}')
    mo.ui.table(table)
    return


@app.cell
def _(mo, metrics):
    mo.md(f"""
    ### Validation summary
    - Terminal pulses: **{metrics['terminal_pulses']['exact']}/{metrics['terminal_pulses']['predicted']}** exact.
    - Secondary diagonal pulses: **{metrics['diagonal_pulses']['exact']}/{metrics['diagonal_pulses']['scored']}** exact,
      {metrics['diagonal_pulses']['unscored']} unevaluable.
    - Post-pulse normalized r=6: **{metrics['r6_normalized']['latest_two_affine_exact']}/32** exact.
    - Terminal next-Z reconstruction with three previous normalized values:
      **{metrics['terminal_Z']['exact']}/18** exact.

    All references are conditional on the separate rank-zero DT/PT calculation.
    """)
    return


if __name__ == "__main__":
    app.run()
