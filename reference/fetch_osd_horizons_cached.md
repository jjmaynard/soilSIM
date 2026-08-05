# Fetch OSD horizon distinctness data for a set of soil series names, with per-series caching

[`soilDB::fetchOSD()`](http://ncss-tech.github.io/soilDB/reference/fetchOSD.md)
is a live network call. Multiple profiles/components in the same
mukey/AOI frequently share the same series name (e.g.
[`simulate_profile_depths_by_mukey()`](https://jjmaynard.github.io/soilSIM/reference/simulate_profile_depths_by_mukey.md)
calls
[`query_osd_distinctness()`](https://jjmaynard.github.io/soilSIM/reference/query_osd_distinctness.md)
once per component), so fetching unconditionally on every call
re-requests identical data redundantly. This wrapper caches each series'
raw `id`/`hzname`/`distinctness` rows on disk - the same
[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)/[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md)
mechanism `R/raster-cache.R` already uses for AOI-keyed raster fetches,
just with a series-name-keyed cache key instead - fetching only the
series not already cached.

## Usage

``` r
fetch_osd_horizons_cached(compnames)
```

## Arguments

- compnames:

  Character vector of soil series names, as would be passed directly to
  [`soilDB::fetchOSD()`](http://ncss-tech.github.io/soilDB/reference/fetchOSD.md).

## Value

A data frame with columns `id`, `hzname`, `distinctness` - one row per
horizon, across all requested series (any already-cached plus any newly
fetched).

## Known quirk preserved (not a bug introduced here)

[`soilDB::fetchOSD()`](http://ncss-tech.github.io/soilDB/reference/fetchOSD.md)
returns its `id` column upper-cased regardless of the requested case
(e.g. requesting `"amador"` returns rows with `id == "AMADOR"`) -
confirmed via a live call before writing this function. Splitting a
combined multi-series fetch into per-series cache entries is done via
case-insensitive matching on `id` to account for this. This casing
difference was already present in
[`soilDB::fetchOSD()`](http://ncss-tech.github.io/soilDB/reference/fetchOSD.md)'s
output before this caching layer existed, and downstream code
([`simulate_and_perturb_soil_profiles()`](https://jjmaynard.github.io/soilSIM/reference/simulate_and_perturb_soil_profiles.md))
was already unaffected by it, since it joins
[`query_osd_distinctness()`](https://jjmaynard.github.io/soilSIM/reference/query_osd_distinctness.md)'s
output by `genhz` only, never by `id`.
