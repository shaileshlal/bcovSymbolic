"""Interactive source-data audit for quintic arithmetic patterns."""
import marimo

__generated_with = "0.25.1"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    from bcov.datasets import load_amplitude
    from bcov.patterns import audit_patterns, pure_x
    return audit_patterns, load_amplitude, mo, pure_x


@app.cell
def _(mo):
    mo.md(r"""
    # Quintic arithmetic laboratory

    Explore exact pure-$X$ coefficients from the supplied genus 2–15 free energies.
    The constant and linear Bernoulli recurrences are source-verified identities;
    they do not determine all holomorphic ambiguities. Higher-depth formulas from
    the separate AI framework are not included until their source is available.
    """)
    genus = mo.ui.slider(2, 15, value=3, label="Genus", show_value=True)
    depth = mo.ui.slider(1, 15, value=4, label="Inclusive jet depth", show_value=True)
    mo.hstack([genus, depth])
    return depth, genus


@app.cell
def _(audit_patterns, depth, genus, load_amplitude, pure_x):
    amplitude = load_amplitude(genus.value)
    report = audit_patterns(amplitude, depth.value)
    hkq = pure_x(amplitude, depth.value, "hkq")
    small_z = pure_x(amplitude, depth.value, "small_z")
    return hkq, report, small_z


@app.cell
def _(hkq, mo, report, small_z):
    mo.ui.table([
        {"power": k, "YY": str(y), "HKQ": str(h), "small-z": str(z), "YY / constant": str(r)}
        for k, (y, h, z, r) in enumerate(zip(report.yy_jet, hkq, small_z, report.ratios))
    ], label="Exact coefficients; small-z is transformed before truncation")
    return


@app.cell
def _(mo, report):
    mo.md(f"""
    **Constant and linear recurrence:** {report.low_x_matches}
    **Full denominator lattice:** {report.denominator_matches}

    | Frame | Least common denominator |
    |---|---:|
    | HKQ | {report.lcd_hkq} |
    | small-z | {report.lcd_small_z} |
    | YY | {report.lcd_yy} |

    YY and HKQ coefficient ratios agree. HKQ and small-z share the full
    denominator lattice, but individual coefficients mix under $X\\mapsto1-X$.
    """)
    return


if __name__ == "__main__":
    app.run()
