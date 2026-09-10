# Contributing to soilSIM

## Naming conventions (enforced by `tests/testthat/test-naming-conventions.R`)

soilSIM is a **generic simulation/fusion core** plus **data-source adapters**. Names must make
that split legible and must use one verb per action.

### Verb lexicon

| Verb | Meaning | Returns |
|---|---|---|
| `fetch_` | acquire data from a network service or disk cache | data frame / raster |
| `read_` / `write_` | local file I/O | data / invisible path |
| `fit_` | estimate distribution/model parameters | parameter list or model object |
| `quantile_` | evaluate a quantile function | numeric |
| `simulate_` | draw random realizations (never `sim_`) | array / data frame / result object |
| `fuse_` | Bayesian combination of prior + likelihood (never `bayes_`/`bayesian_`) | fused params / raster / result object |
| `estimate_` | infer a summary quantity from data (correlations, moments) | matrix / numeric |
| `compute_` | deterministic derivation from inputs (never `calculate_`) | numeric / data frame |
| `analyze_` | high-level, multi-step, user-facing wrapper | result object |
| `validate_` | assertion: returns `invisible(TRUE)` or throws (never `check_`) | invisible logical / error |
| `diagnose_` | produce a QA/diagnostic report object (not a pass/fail) | `soilSIM_diagnostics` |
| `default_` | return a package default config/data object (never `get_default_*`) | list |
| `convert_` | unit / representation transform | same type as input |
| `prepare_` | shape raw input into the form an engine expects | data frame |
| `is_` / `has_` / `*_exists` | predicate | logical scalar |

`assess_` is allowed **only** for non-exported single-metric sub-checks.

### Word order

`verb_object_qualifier`, verb first. Not `mukey_draws_lookup` — use `lookup_mukey_draws`.

### Initialisms

Always lower-case: `_gp`, `_aws`, `_vg`, `_ilr`, `_wkt`, `_osd`, `_nrcs`, `_kssl`, `_rfv`.
Domain keys stay literal: `mukey`, `cokey`, `ssurgo`, `solus`.

### Adapter functions

`<verb>_<source>_<object>`, e.g. `fetch_ssurgo_mukey_raster`, `process_ssurgo_horizons`.
Generic-core functions never name a data source.

### Forbidden suffixes

No `_working`, `_compatible`, `_safe`, `_v2`, `_new`, `_old` anywhere. The only allowed capability
suffixes are `_raster` (raster-native counterpart of a scalar function) and `_multiproperty`
(vectorized over a set of properties).

### Tiered exports

- **Tier 1** — pipeline entry points. Full `@examples`, covered by a vignette.
- **Tier 2** — supported building blocks. Exported, `@family`-tagged.
- **Tier 3** — everything else. `@keywords internal`; drop `@export` unless a vignette/test needs it.

## Deprecation

Renamed exports get a shim in `R/deprecated.R` using `lifecycle::deprecate_soft()`. Timeline:
`x.1.0` soft-deprecate → `x.2.0` warn → `x.3.0` remove.

## Result objects

User-facing engines return classed lists (see `R/classes.R`): `soilSIM_simulation`,
`soilSIM_gp_models`, `soilSIM_fusion`, `soilSIM_statistics`, `soilSIM_diagnostics`,
`soilSIM_property_config`. Add a `print`/`summary`/`plot` method rather than expecting users to
walk nested lists.

## Logging

Any new exported function that calls `log_message()` adds `verbose = getOption("ssurgo.verbose", FALSE)`
plus the standard `set_verbose_logging()` / `on.exit()` two-liner as its first body lines.

## Before opening a PR

```r
devtools::document()
devtools::test()
devtools::check(vignettes = TRUE)   # 0 errors / 0 warnings
```
