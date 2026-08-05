# Temporarily raise the package's log verbosity for the duration of a call

Call at the top of a function with its own `verbose` argument,
immediately followed by
`on.exit(options(soil_workflow_log_config = <the returned value>), add = TRUE)`.
When `verbose = TRUE`, raises the shared `soil_workflow_log_config`
option's `log_level` to `raised_level` for the remainder of the calling
function's execution (restored via the caller's own
[`on.exit()`](https://rdrr.io/r/base/on.exit.html)), so every
[`log_message()`](https://jjmaynard.github.io/soilSIM/reference/log_message.md)
call made anywhere during that call - including deep inside shared
helpers that have no `verbose` parameter of their own, and across
package/file boundaries - becomes visible. Only ever *raises* the level,
never lowers it, so nested calls compose correctly regardless of order:
an outer `verbose = TRUE` makes everything below it visible even if an
inner helper's own `verbose` is `FALSE`, and an inner call's own restore
is a harmless no-op when the level was already raised by an outer
caller. When `verbose = FALSE` (the default everywhere), does nothing -
the ambient default (see
[`log_message()`](https://jjmaynard.github.io/soilSIM/reference/log_message.md))
is already quiet.

## Usage

``` r
set_verbose_logging(verbose, raised_level = "INFO")
```

## Arguments

- verbose:

  Logical.

- raised_level:

  Log level to switch to when `verbose = TRUE` (default `"INFO"` -
  matches the package's historical visible-by-default narration;
  `"DEBUG"`-level detail remains available only via an explicit
  `setup_logging(log_level = "DEBUG")` call, not through this per-call
  `verbose` mechanism).

## Value

Invisibly, the log config in effect *before* this call (pass to
[`on.exit()`](https://rdrr.io/r/base/on.exit.html) to restore it).
