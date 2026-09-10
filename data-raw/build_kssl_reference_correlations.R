# One-time, maintainer-run script - NOT part of the package build (data-raw/
# is excluded via .Rbuildignore). Regenerates R/sysdata.rda from the two
# KSSL-derived, genhz-keyed correlation-matrix RDS files bundled alongside
# this script (data-raw/global_cor_matrices.rds,
# data-raw/global_cor_texture_matrices.rds). Those files are static reference
# data fit once from KSSL lab data. See R/kssl-reference-correlations.R for
# how soilSIM consumes them. Re-run this script only if the KSSL fit is ever
# redone (and refresh the two RDS files from the new fit first).

kssl_property_matrices <- readRDS(file.path("data-raw", "global_cor_matrices.rds"))
kssl_texture_matrices <- readRDS(file.path("data-raw", "global_cor_texture_matrices.rds"))

save(
  kssl_property_matrices, kssl_texture_matrices,
  file = file.path("R", "sysdata.rda"),
  compress = "xz"
)
