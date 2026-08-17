# Evaluate percentile-reconstruction methods against a known ground-truth distribution

For each replicate: draws a large sample from `true_rng`, computes the
requested percentiles from it (simulating "what we'd observe" as summary
data), reconstructs a distribution from those percentiles using each
method, then scores the reconstruction against the true draws via the
Kolmogorov-Smirnov statistic and relative mean/SD error. Averages scores
across replicates.

## Usage

``` r
validate_percentile_methods_synthetic(
  true_rng,
  percentile_probs = c(0, 0.05, 0.5, 0.95, 1),
  methods = c("linear_cdf", "spline", "kde"),
  n_true = 5000,
  n_sim = 1000,
  n_reps = 20,
  bounds = NULL
)
```

## Arguments

- true_rng:

  A function(n) that draws n samples from the ground-truth distribution.

- percentile_probs:

  Numeric vector of probabilities (0-1) to reveal as "known"
  percentiles.

- methods:

  Character vector of methods to evaluate.

- n_true:

  Number of ground-truth draws per replicate (used both to derive the
  revealed percentiles and as the KS-test reference sample).

- n_sim:

  Number of samples to reconstruct per method per replicate.

- n_reps:

  Number of independent replicates to average over.

- bounds:

  Optional bounds passed through to methods that use them.

## Value

A data frame with one row per method: mean KS statistic, mean absolute
relative mean error, and mean absolute relative SD error across
replicates.
