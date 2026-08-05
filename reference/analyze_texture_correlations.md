# Analyze texture (sand/silt/clay) correlations, in both raw and ILR space

Raw Pearson correlations among compositional (simplex-constrained) parts
like sand/silt/clay percentages are spuriously negative due to the
sum-to-100 constraint - `ilr_correlations` (isometric log-ratio space,
via the dependency-free
[`ilr_forward()`](https://jjmaynard.github.io/soilSIM/reference/ilr_forward.md)
in `distributions.R`) avoids that artifact and is the statistically
defensible view for compositional data. Both are reported:
`raw_correlations` for continuity with prior behavior,
`ilr_correlations` as the added analytical value (only computed when all
three fractions are present with at least 5 complete observations).

## Usage

``` r
analyze_texture_correlations(data, texture_properties, methods, config)
```

## Arguments

- data:

  Input data with texture `_r` columns.

- texture_properties:

  Character vector of available texture column names (e.g. a subset/all
  of `c("sandtotal_r","silttotal_r","claytotal_r")`).

- methods:

  Correlation methods (e.g. `c("pearson","spearman")`).

- config:

  Unused; kept for interface compatibility with callers.

## Value

`list(raw_correlations=, ilr_correlations=, n_observations=, note=)`.
