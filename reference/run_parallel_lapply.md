# Run a Function Over a List/Vector, in Parallel via `future.apply`, with Sequential Fallback

Run a Function Over a List/Vector, in Parallel via `future.apply`, with
Sequential Fallback

## Usage

``` r
run_parallel_lapply(
  X,
  FUN,
  ...,
  n_cores = NULL,
  future_seed = FALSE,
  op_name = "Parallel operation",
  sequential_fallback = NULL,
  catch_errors = TRUE
)
```

## Arguments

- X:

  A list or vector to iterate over.

- FUN:

  A function taking one element of `X` as its first argument, plus
  `...`.

- ...:

  Additional fixed arguments forwarded to `FUN` on every call.

- n_cores:

  Number of parallel workers. `NULL` (default) resolves to
  `max(1, parallel::detectCores() - 1)`, then is clamped to
  `min(n_cores, length(X))` (no point starting more workers than there
  is work). A resolved value `<= 1` skips
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  setup entirely and just runs `lapply(X, FUN, ...)`.

- future_seed:

  Passed through to
  [`future.apply::future_lapply()`](https://future.apply.futureverse.org/reference/future_lapply.html)'s
  `future.seed` argument. Use `TRUE` whenever `FUN` (or anything it
  calls) uses R's RNG, to get `future`'s parallel-safe RNG streams -
  this is more callers than it first appears: besides the obvious cases
  (Monte Carlo draws, profile perturbation),
  [`GPfit::GP_fit()`](https://rdrr.io/pkg/GPfit/man/GP_fit.html)'s
  internal hyperparameter search is a genetic algorithm and genuinely
  uses RNG too (confirmed empirically via `future_lapply()`'s
  "UNRELIABLE VALUE" warning when this was left `FALSE` for a GP-fitting
  caller in this package), even though its fitted result is highly
  stable across runs. Default `FALSE` matches `future_lapply()`'s own
  default; only leave it `FALSE` for `FUN` you've confirmed doesn't
  touch the RNG anywhere in its call chain.

- op_name:

  A short label used in the log message and as
  [`handle_workflow_error()`](https://jjmaynard.github.io/soilSIM/reference/handle_workflow_error.md)'s
  context string if dispatch fails.

- sequential_fallback:

  A zero-argument closure invoked if dispatch fails and
  `catch_errors = TRUE` (the caller pre-binds its own arguments, e.g.
  `function() run_sequential_simulation( simulation_params, ...)`).
  `NULL` (default) falls back to plain `lapply(X, FUN, ...)` instead.
  Ignored when `catch_errors = FALSE`.

- catch_errors:

  If `TRUE` (default), wraps dispatch in
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html) and applies the
  `sequential_fallback` contract above on error. If `FALSE`, dispatch
  errors propagate to the caller's own error handling instead (plan
  setup/restoration via
  [`on.exit()`](https://rdrr.io/r/base/on.exit.html) still happens
  regardless) - needed for callers whose own contract is "no sequential
  fallback, propagate/ convert the error myself" (e.g.
  [`simulate_profile_depths_by_collection_parallel()`](https://jjmaynard.github.io/soilSIM/reference/simulate_profile_depths_by_collection_parallel.md),
  which returns `NULL` on error via its own outer
  [`tryCatch()`](https://rdrr.io/r/base/conditions.html)).

## Value

A list of `FUN(x, ...)` results, one per element of `X`, in `X`'s
order - the same shape
[`parallel::parLapply()`](https://rdrr.io/r/parallel/clusterApply.html)/`mclapply()`
already produced, so callers that previously consumed those results
(e.g. via
[`dplyr::bind_rows()`](https://dplyr.tidyverse.org/reference/bind_rows.html)
or manual concatenation) need no change.
