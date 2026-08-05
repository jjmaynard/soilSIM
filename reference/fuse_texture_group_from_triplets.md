# Fuse two full clay/sand/silt low-rep-high triplets jointly via ILR fusion

Wraps `distributions.R`'s
[`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md)
(once per side) +
[`fuse_bivariate_normal()`](https://jjmaynard.github.io/soilSIM/reference/fuse_bivariate_normal.md) +
[`sample_ilr_posterior()`](https://jjmaynard.github.io/soilSIM/reference/sample_ilr_posterior.md)
for future joint texture-prior updating. Independent per-fraction fusion
(e.g. three separate
[`fuse_beta()`](https://jjmaynard.github.io/soilSIM/reference/fuse_beta.md)
calls) measurably breaks sum-to-100 (documented upstream: up to 10.5
percentage points on realistic synthetic data) - this is why the fusion
happens jointly, in ILR space, instead.

## Usage

``` r
fuse_texture_group_from_triplets(
  prior_triplets,
  lik_triplets,
  z_prior,
  z_lik,
  n_mc = 2000,
  n_samples = 1000,
  total = 100
)
```

## Arguments

- prior_triplets, lik_triplets:

  Named list with elements `clay`, `sand`, `silt`, each a length-3
  numeric vector `c(low, rep, high)`. These keys are positional-role
  placeholders matching
  [`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md)'s
  `low_clay`/`low_sand`/`low_silt` parameter naming (position 1/2/3 of
  the ILR sequential binary partition - see `distributions.R`'s ILR
  section header), not an identity requirement; callers (e.g.
  `monte-carlo.R`'s
  [`fuse_observed_data_into_priors()`](https://jjmaynard.github.io/soilSIM/reference/fuse_observed_data_into_priors.md))
  build them from `composition_groups$texture$members` in configured
  order.

- z_prior, z_lik:

  Standard-normal quantile matching each side's low/high interval (e.g.
  `qnorm(0.95)` for a 5th/95th-percentile low/high).

- n_mc:

  Monte Carlo sample size passed to
  [`estimate_ilr_moments_mc()`](https://jjmaynard.github.io/soilSIM/reference/estimate_ilr_moments_mc.md).

- n_samples:

  Number of posterior composition samples to draw.

- total:

  Composition target sum (default 100).

## Value

`list(posterior_samples = <n_samples x 3 matrix, columns clay/sand/silt, guaranteed sum-to-`total` by construction>, ilr_mu=, ilr_Sigma=)`.
