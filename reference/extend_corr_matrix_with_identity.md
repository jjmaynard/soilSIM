# Extend a Correlation Matrix with Identity Rows/Columns for Missing Names

Adds a 1-on-diagonal, 0-cross-correlation row/column for each name in
`missing_names` not already in `m`'s dimnames - the same "uncorrelated
with everything, safe default" convention
[`build_kssl_fallback_matrix()`](https://jjmaynard.github.io/soilSIM/reference/build_kssl_fallback_matrix.md)
uses for properties without KSSL reference data (see
MULTI_PROPERTY_FUSION_PLAN.md task P2). Used by
[`simulate_cokey_generalized()`](https://jjmaynard.github.io/soilSIM/reference/simulate_cokey_generalized.md)
as a defensive fallback so an unrecognized-but-requested property
degrades to independence rather than erroring on
`m[keep_cols, keep_cols]`. The result stays positive-definite: a
block-diagonal matrix formed from a positive-definite block (`m`) and an
identity block has no negative eigenvalues.

## Usage

``` r
extend_corr_matrix_with_identity(m, missing_names)
```

## Arguments

- m:

  A square, dimnamed correlation matrix.

- missing_names:

  Character vector of names to add (assumed disjoint from
  `rownames(m)`).

## Value

`m` extended to `nrow(m) + length(missing_names)` rows/columns.
