# KSSL Reference Correlation Matrices (Optional Fallback)

Static, pre-computed, genetic-horizon-keyed (O/A/E/B/C/Cr, plus R for
the texture matrix) correlation matrices fit once from KSSL lab data by
an older, unrelated codebase (`code_ref/brdf/property_simulation.R`,
outside this package). These functions let `monte-carlo.R`'s
correlation-structure estimation optionally fall back to this reference
data (via `config$monte_carlo$correlation_fallback = "kssl_global"`)
instead of a plain identity matrix when there isn't enough SSURGO data
to estimate correlations empirically. The underlying matrices are stored
internally as package sysdata (see
`data-raw/build_kssl_reference_correlations.R`)

- not exported/documented as a public dataset, since they're
  implementation detail of this fallback feature, not a general-purpose
  soilSIM dataset.
