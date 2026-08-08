# Shared `future`/`future.apply` Parallel Dispatch Helper

A single internal helper standardizing this package's parallel
processing on `future`/`future.apply`
([`future::multisession`](https://future.futureverse.org/reference/multisession.html),
which behaves identically on Windows and Unix), replacing the
hand-rolled base-`parallel::` boilerplate (OS branching between
[`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html)/`parLapply()`
on Windows and
[`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
elsewhere, explicit
[`library(soilSIM)`](https://jjmaynard.github.io/soilSIM/) on Windows
workers, explicit `soil_workflow_log_config` propagation) that used to
be duplicated across `process_cokeys_parallel()`
(`multivariate-adjustment.R`),
[`run_parallel_simulation()`](https://jjmaynard.github.io/soilSIM/reference/run_parallel_simulation.md)
(`monte-carlo.R`), and
[`maybe_adjust_soil_data_depth_trend()`](https://jjmaynard.github.io/soilSIM/reference/maybe_adjust_soil_data_depth_trend.md)
(`ssurgo-simulation.R`).
