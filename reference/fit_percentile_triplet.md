# Fit a distribution to a low/representative/high percentile triplet

The single dispatcher
[`extract_property_parameters()`](https://jjmaynard.github.io/soilSIM/reference/extract_property_parameters.md)
(`monte-carlo.R`) calls to turn a SSURGO `_l/_r/_h` triplet into a
family-appropriate fit.

## Usage

``` r
fit_percentile_triplet(
  l,
  r,
  h,
  family,
  lh_probs = c(0.05, 0.95),
  bounds = NULL,
  boundedness = if (is.null(bounds)) "u" else "b"
)
```

## Arguments

- l, r, h:

  Low/representative/high values.

- family:

  One of `"triangular"`, `"uniform"`, `"normal"`, `"lognormal"`,
  `"beta"`, `"metalog"`, `"linear_cdf"`, `"auto"`.

- lh_probs:

  Length-2 vector giving the probabilities `l`/`h` are assumed to
  represent (SSURGO's low/high are professional-judgment bounds, not
  confirmed exact percentiles - this makes that assumption explicit and
  overridable rather than silently hardcoded).

- bounds:

  Optional `c(lower, upper)` physical bounds; required for `"beta"` and
  for `"metalog"` when `boundedness != "u"`.

- boundedness:

  Metalog boundedness; defaults to `"u"` when `bounds` is `NULL`, else
  `"b"`.

## Value

`list(family=, fit=, valid=)`. `valid=FALSE` when any of `l`/`r`/`h` is
non-finite.
