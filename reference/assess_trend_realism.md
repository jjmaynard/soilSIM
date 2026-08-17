# Assess Depth-Trend Realism

Reuses the already-real `assess_trend_monotonicity()` and
`assess_realistic_values()` (`gp-modeling.R`) to judge whether a
predicted depth trend is realistic, and counts constraint violations
against `criteria$realistic_ranges` when supplied.

## Usage

``` r
assess_trend_realism(predictions, test_depths, prop, criteria)
```

## Arguments

- predictions:

  GP predictions at `test_depths`.

- test_depths:

  Depths predictions were made at.

- prop:

  Property name (matched against `criteria$realistic_ranges` and
  `assess_realistic_values()`'s built-in ranges).

- criteria:

  List, optionally with `realistic_ranges` (named list of
  `list(min=, max=)`).

## Value

List with `realistic_trend`, `constraint_violations`, `trend_quality`,
`realism_assessment`.
