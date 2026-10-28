# Experimental results, grouped by *evaluation protocol*

These are the **archived stage benchmarks**, not newly derived theorems.
The exact stage-specific validation code and datasets are bundled.
Several model classes were developed retrospectively during exploratory
analysis; later values held out from parameter fits do not undo that history.

| Stage | Comparison | Archived result | Status |
|---|---|---|---|
| v0.1 | Low-X normalized log cumulants, train g40–53, holdout g54–1000 | ~0.0456% median ratio error | Numerical; exact low-X labels from BCOV + GV vanishings |
| v0.1 | Mixed YY source ExtraTrees, train g2–8, holdout g9–12 | ~80.7% sign accuracy | Genuine held-out genera; not exact ambiguity reconstruction |
| Intermediate dataset | Exact residual YY sections g2–12 | 22 exact rational labels | Reference-derived anchor, exact round trips |
| Intermediate ML | Magnitude source normalization, rolling g9–12 | ~0.12 dex mean magnitude error | Signs not solved |
| GV lattice | Minimum-height residual reconstruction g4–12 | 22/22 | All relevant GV charges vanish; does not predict nonzero displacements |
| Displacement | Extremal phase/quadratic benchmark | 8/8 integer holdouts | Labels generated from a known extremal formula |
| Genus 45–53 | Nonzero primitive charge at (51,20) | 175, exact BCOV sensitivity g52–53 | Uses known input; zero-charge branch comparison |
| DT/PT v0.4 | Degree-residue arithmetic, degrees 21–30 | 10/10 leading integers | DT-charge-aware representation; reference modular conditional |
| PT v0.5 | Near-boundary finite differences, degrees 21–30 | 28/48 PT exact | Feature/model class exploratory |
| Defect v0.6 | Sequential blocks, degrees 21–30 | 31/48 PT exact | Updates after first block; not fixed-window baseline |
| Recursive v0.7 | Degrees 26–30 | 17/24 PT exact | Sparse Newton-pulse model |
| Survival v0.8 | Degrees 31–53 sequential vs autonomous | 106/109 and 99/109 exact | Conditional reference; frozen baseline |
| Deep v0.9 | Depths 6–10 degrees 41–53 | 32/32 exact, 33 abstain | Confidence-gated; conditional reference |
| Pulses v1.0 | Depths 6–10 degrees 41–53 | 59/59 exact, six abstain | Affine profiles, selected retrospectively |
| Depth v1.1 | Depths 11–15 degrees 41–53 | 25/25 exact, 40 abstain | Paired scale model; 8 new anchors for modeled sites |
| Depth v1.1 | Forward degrees 54–70 | 49/49 exact, 36 abstain | Further check, conditional reference |
| Nested v1.2 | Terminal pulse prediction | 18/18 exact from 13 calibration integers | Conditional normalization |
| Nested v1.2 | Secondary-pulse grid | 17/19 evaluable exact | Two explicit misses |
| Nested v1.2 | Stable r=6 residual values | 32/32 exact | Requires mature-scale/normalization information |

**These are not a common leaderboard.** They evaluate different targets, depths,
training windows, abstention policies and supplied exact inputs. In particular,
an empirical predictor need not use modular forms to be evaluated against a
modularity-conditional reference. Treat all of these as numerical baselines.
