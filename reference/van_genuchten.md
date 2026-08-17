# Evaluate the van Genuchten Water Retention Curve

Closed-form van Genuchten (1980) volumetric water content at a given
matric potential.

## Usage

``` r
van_genuchten(h, alpha, n, theta_r, theta_s)
```

## Arguments

- h:

  Matric potential (cmH2O, typically negative).

- alpha, n:

  Van Genuchten shape parameters.

- theta_r, theta_s:

  Residual and saturated volumetric water content.

## Value

Numeric vector of volumetric water content, `theta(h)`.

## Examples

``` r
van_genuchten(h = -100, alpha = 0.02, n = 1.3, theta_r = 0.05, theta_s = 0.45)
#> [1] 0.350325
```
