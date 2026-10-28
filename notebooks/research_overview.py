"""Explore the integrated 18-stage research archive and exact section geometry."""
import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import json
    import marimo as mo
    from quintic_bcov_lab import backend, datasets, runner
    from quintic_bcov_lab.arithmetic import difference, newton_forecast
    from quintic_bcov_lab.sections import from_reference
    return backend, datasets, difference, from_reference, json, mo, newton_forecast, runner


@app.cell
def _(mo, runner):
    mo.md("""
    # Quintic sequential research

    The 18 original studies retain their predictors, calibration data, validators,
    and saved results. Unified exact analysis uses this repository's BCOV engine.
    Archived benchmarks are conditional research observations, not newly run
    results or proofs of all-genus reconstruction. Abstentions remain abstentions.
    """)
    stage = mo.ui.dropdown(
        options={s["title"]: s["id"] for s in runner.stages()},
        value=runner.stages()[-1]["title"], label="Study",
    )
    stage
    return (stage,)


@app.cell
def _(backend, mo, runner, stage):
    selected = runner.resolve_stage(stage.value)
    status = backend.validate_small()
    mo.md(f"""
    **{selected['id']}** — {selected['title']}

    Active unified backend: `{status['backend']}`. Exact genus-two residual: zero.

    Reproduce this stage in a writable copy:
    ```bash
    quintic-bcov stage {selected['id']} --action test
    quintic-bcov stage {selected['id']} --action run
    ```
    Fresh logs go to `.bcov-research/run_logs/`. The packaged archive stays intact.
    """)
    return


@app.cell
def _(datasets, mo):
    mo.ui.table(datasets.intermediate_labels(), label="22 exact reference-derived residual labels")
    genus = mo.ui.slider(2, 12, value=4, label="Section genus", show_value=True)
    genus
    return (genus,)


@app.cell
def _(datasets, from_reference, genus, mo):
    amplitude = next(a for a in datasets.read_reference()["amplitudes"] if a["genus"] == genus.value)
    section = from_reference(genus.value, amplitude["ambiguity"])
    mo.ui.table([
        {"index": i, "YY residual coefficient": str(c), "reduced Bernstein coefficient": str(b)}
        for i, (c, b) in enumerate(zip(section.coefficients, section.x_bernstein))
    ], label="Reference-derived section; not a blind prediction")
    return


@app.cell
def _(datasets, difference, mo, newton_forecast):
    # Only calibration observations enter this finite-difference illustration.
    calibration = datasets.read_csv("18", "data/terminal_pulse_calibration.csv")
    seed = [int(row["terminal_pulse"]) for row in calibration if row["residue"] == "0"]
    forecast = newton_forecast(seed, order=2, steps=3)
    mo.md(f"""
    ## Terminal-pulse example (residue zero)

    Calibration values: `{seed}`. Second difference: `{difference(seed, 2)}`.

    Quadratic extrapolation to the next three sites: `{[str(x) for x in forecast]}`.
    This illustrates the empirical terminal-pulse hypothesis using calibration
    alone. Mature-scale normalization is still supplied information; this is
    not a complete blind PT or BCOV reconstruction.
    """)
    return


@app.cell
def _(json, mo, runner, stage):
    study = runner.ROOT/"studies"/stage.value
    result_path = study/"results/benchmark.json"
    if result_path.is_file():
        display = mo.json(json.loads(result_path.read_text()))
    else:
        display = mo.md("This stage keeps its benchmark in other files; see its original README.")
    mo.vstack([mo.md("### Supplied benchmark snapshot (not a fresh run)"), display])
    return


if __name__ == "__main__":
    app.run()
