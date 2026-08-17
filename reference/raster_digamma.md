# `fitdistrplus::fitdist()`'s generic `optim()`-based MLE has no vectorized form. The Beta log-likelihood's score/Hessian equations have a known closed form in terms of `digamma()`/`trigamma()` (both vectorized base-R functions), so a fixed-iteration-count Newton-Raphson solver updates ALL cells' `(alpha, beta)` estimates simultaneously per iteration.

[`fitdistrplus::fitdist()`](https://lbbe-software.github.io/fitdistrplus/reference/fitdist.html)'s
generic [`optim()`](https://rdrr.io/r/stats/optim.html)-based MLE has no
vectorized form. The Beta log-likelihood's score/Hessian equations have
a known closed form in terms of
[`digamma()`](https://rdrr.io/r/base/Special.html)/[`trigamma()`](https://rdrr.io/r/base/Special.html)
(both vectorized base-R functions), so a fixed-iteration-count
Newton-Raphson solver updates ALL cells' `(alpha, beta)` estimates
simultaneously per iteration.

## Usage

``` r
raster_digamma(r)
```
