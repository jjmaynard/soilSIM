# KDE sampler: seeds truncated-normal draws around each percentile, then resamples from their kernel density estimate

KDE sampler: seeds truncated-normal draws around each percentile, then
resamples from their kernel density estimate

## Usage

``` r
sim_kde(probs, values, n, sample_size = 10000)
```

## Arguments

- sample_size:

  Number of truncated-normal draws used to build the KDE (higher =
  smoother density estimate, at more compute cost).
