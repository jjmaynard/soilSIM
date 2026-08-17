# Disk Cache for Raster Fusion Fetch Results

[`build_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/build_cache_key.md)/[`cache_get()`](https://jjmaynard.github.io/soilSIM/reference/cache_get.md)/[`cache_set()`](https://jjmaynard.github.io/soilSIM/reference/cache_set.md)/`CACHE_TTL_SECONDS`
are called by
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)/[`run_stage1_fusion_group()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion_group.md)
(`R/raster-fusion.R`) and
[`simulate_ssurgo_mapunit_draws()`](https://jjmaynard.github.io/soilSIM/reference/simulate_ssurgo_mapunit_draws.md)
(`R/ssurgo-simulation.R`) but were genuinely undefined anywhere in the
source bundle they were ported from (`code_ref/reanalysis-platform/`) -
confirmed by a repo-wide grep; `HANDOFF_NOTES.md` lists them as
externals "assumed to already exist" in that bundle's *target* project,
which was never soilSIM.

Rather than inventing a new scheme or adding a dependency
(`memoise`/`R.cache`), this adapts the already-working, already-tested
disk-RDS cache pattern from `R/ssurgo-acquisition.R`'s
[`generate_ssurgo_cache_key()`](https://jjmaynard.github.io/soilSIM/reference/generate_ssurgo_cache_key.md)/[`check_ssurgo_cache()`](https://jjmaynard.github.io/soilSIM/reference/check_ssurgo_cache.md)/
[`cache_ssurgo_data()`](https://jjmaynard.github.io/soilSIM/reference/cache_ssurgo_data.md)
(`digest`-based keying, age/TTL invalidation via file mtime) -
generalized here to the AOI/property-id/depth/kind shape the raster
fusion code expects.

Cache files live under `tools::R_user_dir("soilSIM", "cache")` (the
standard R \>= 4.0 per-package cache location) rather than a
caller-supplied directory, since
[`run_stage1_fusion()`](https://jjmaynard.github.io/soilSIM/reference/run_stage1_fusion.md)'s
own calls (`cache_get(key, CACHE_TTL_SECONDS)`,
`cache_set(key, "ssurgo", data)`) don't thread a directory argument
through at all.

## Usage

``` r
CACHE_TTL_SECONDS
```

## Format

An object of class `numeric` of length 1.
