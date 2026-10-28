# Data dictionary

`data/intermediate_training.csv` has one row per `(genus, residual_index)` and only genera with `D>0`.

- `genus`: integer `g`
- `residual_index`: `r=0,...,D-1`
- `X_pivot_degree`: `j=m+1+r`; pure-`X` ambiguity reference pivot
- `K=3g-3`, `m=ceil(K/5)`, `D=floor((2g-2)/5)`
- `relative_pivot_degree`: decimal `j/K`, feature only
- `a_pivot_phys`: true YY coefficient `a_{g,j}`, exact rational label, **never** input on a blind test
- `c_X_kernel`: exact rational kernel coordinate in `X^{r+1}(1-X)^m`; **target**
- `b_x_bernstein`: same coordinate in `x^{r+1}(1-x)^{D-r-1}`, **target**
- `label_sign`, `label_signed_log10`: derived target statistics, labels, not model features
- `mixed_terms_at_Xj`, `mixed_v1_terms_at_Xj`, `mixed_v2_terms_at_Xj`, `mixed_v3_terms_at_Xj`: term-count features from the *nonholomorphic anomaly polynomial*, excluding pure X terms
- `mixed_log10_max_at_Xj`: logarithm of the maximum absolute rational coefficient in that X slice, approximate feature

All exact rational values are canonical `"n/d"` or `"n"` strings. Reconstruct them with `fractions.Fraction(value)`.

`data/intermediate_exact.json` contains for each `genus` the exact vectors `physical`, `anchor`, `residual`, `c`, `bernstein_x`, `xQ`, plus the pivot indices and dimensions. All serialized fractions are strings. `anchor` is **reference-derived**; use it for decomposition and validation, **not** as blind input data.
