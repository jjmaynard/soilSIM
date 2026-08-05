# Reduce an `observed_data` entry to an unnamed `c(low, rep, high)` triplet

Accepts a `c(low=,rep=,high=)`/`list(low=,rep=,high=)` triplet, an
unnamed length-3 `c(low,rep,high)` vector, or a raw sample vector
(reduced via its own empirical quantiles at `lh_probs`).

## Usage

``` r
as_lrh_triplet(x, lh_probs = c(0.05, 0.95))
```
