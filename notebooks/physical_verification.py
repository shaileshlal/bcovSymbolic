"""Interactive exact physical q-series versus GV multicover verification."""
import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    from importlib.resources import as_file, files
    from bcov.verify_physical_qseries import verify
    return as_file, files, mo, verify


@app.cell
def _(mo):
    mo.md(r"""
    # Physical GW/GV comparison

    Evaluate HKQ atomic amplitudes using periods, the inverse mirror map and
    the physical normalization, then compare positive-degree coefficients
    with the GV multicover formula. Every residual below is rational and exact.
    The constant-map coefficient is reported separately.
    """)
    max_genus = mo.ui.slider(2, 15, value=3, label="Maximum genus", show_value=True)
    degree = mo.ui.slider(1, 6, value=4, label="q degree", show_value=True)
    mo.hstack([max_genus, degree])
    return degree, max_genus


@app.cell
def _(as_file, degree, files, max_genus, verify):
    data = files("bcov").joinpath("data")
    with as_file(data.joinpath("bcov_hkq_atomic_genus_2_15.jsonl.gz")) as atomic:
        with as_file(data.joinpath("GVTab.csv.gz")) as gv:
            results = verify(atomic, gv, degree.value, max_genus.value)
    return (results,)


@app.cell
def _(mo, results):
    mo.ui.table([
        {"genus": row["genus"], "q power": d, "coefficient": c, "residual": residual}
        for row in results
        for d, (c, residual) in enumerate(zip(row["coefficients"], row["GV_residuals"]), 1)
    ])
    return


if __name__ == "__main__":
    app.run()
