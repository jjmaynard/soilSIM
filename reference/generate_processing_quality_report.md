# Generate SSURGO Processing Quality Report

Blends whichever quality signals are available (validation score, row
retention rate, horizon property completeness, component retention) into
a single weighted `[0, 1]` score and letter grade, renormalizing weights
over the signals actually present.

## Usage

``` r
generate_processing_quality_report(
  original_data,
  processed_data,
  horizon_stats,
  component_stats,
  validation_results
)
```

## Arguments

- original_data:

  Data frame prior to processing.

- processed_data:

  Data frame after processing.

- horizon_stats:

  List possibly containing `property_completeness`.

- component_stats:

  List possibly containing `rows_retained`.

- validation_results:

  List possibly containing `overall_quality$score`.

## Value

List with `overall_quality_score`, `quality_grade`,
`data_retention_rate`, a human-readable `processing_summary`, and the
passed-through `horizon_processing`/`component_processing`/
`validation_summary` inputs.
