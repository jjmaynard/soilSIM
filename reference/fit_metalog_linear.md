# Fit a metalog distribution via exact linear solve

When the number of interior percentiles (`p` strictly in `(0,1)`) equals
the number of metalog terms, the fit is an EXACTLY-DETERMINED linear
system (`solve(Y, z)`), not an optimization - validated (upstream) to
~1e-12 against `rmetalog::metalog()`. This is the fast path `rmetalog`
itself takes when its solution is already feasible; see
[`check_metalog_feasible()`](https://jjmaynard.github.io/soilSIM/reference/check_metalog_feasible.md)/[`quantile_metalog_with_fallback()`](https://jjmaynard.github.io/soilSIM/reference/quantile_metalog_with_fallback.md)
for the guard `rmetalog`'s own LP feasibility-correction would otherwise
provide.

## Usage

``` r
fit_metalog_linear(
  interior_values,
  interior_probs,
  bounds = NULL,
  boundedness = "u"
)
```

## Arguments

- interior_values, interior_probs:

  Percentile values/probabilities, both strictly interior (excluding
  p=0/p=1 if present).

- bounds:

  `c(lower, upper)`; required unless `boundedness="u"`.

- boundedness:

  One of `"u"` (unbounded), `"sl"` (semi-bounded below), `"su"`
  (semi-bounded above), `"b"` (bounded).

## Value

`list(a=, term=, bounds=, boundedness=)`.
