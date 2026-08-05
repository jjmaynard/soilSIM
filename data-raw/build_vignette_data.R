# One-time, maintainer-run script - NOT part of the package build (data-raw/ is excluded via
# .Rbuildignore). Pre-fetches real data (live SSURGO/OSD/ROSETTA/SOLUS network calls) ONCE for the
# package vignettes and caches the results under inst/extdata/, so the vignettes themselves build
# offline/deterministically (loading the cached real data) while still showing the real live-fetch
# call in their text for transparency. Re-run only if the vignette data needs refreshing.
#
# Run from within a checkout of the soilSIM package directory (working directory = soilSIM/),
# e.g. via `Rscript data-raw/build_vignette_data.R` from inside soilSIM/.

devtools::load_all(".", quiet = TRUE)  # runs .onLoad(), which defensively unsets a stray PROJ_LIB

extdata_dir <- file.path("inst", "extdata")
if (!dir.exists(extdata_dir)) dir.create(extdata_dir, recursive = TRUE)

# Real AOIs already established and live-verified elsewhere in this package (test suite) - reused
# here rather than inventing new ones.
amador_wkt <- "POLYGON((-120.5 38.5, -120.4 38.5, -120.4 38.6, -120.5 38.6, -120.5 38.5))"
salinas_wkt <- "POLYGON((-121.66 36.60, -121.64 36.60, -121.64 36.62, -121.66 36.62, -121.66 36.60))"

save_step <- function(label, path, expr) {
  cat(sprintf("[%s] starting...\n", label))
  result <- tryCatch(expr, error = function(e) {
    cat(sprintf("[%s] FAILED: %s\n", label, conditionMessage(e)))
    NULL
  })
  if (!is.null(result)) {
    saveRDS(result, path)
    cat(sprintf("[%s] saved -> %s\n", label, path))
  }
  result
}

# ---------------------------------------------------------------------------
# 1. SSURGO tabular download for the Amador-area AOI (vignettes 1 and 3)
# ---------------------------------------------------------------------------
ssurgo_amador <- save_step(
  "ssurgo_amador",
  file.path(extdata_dir, "ssurgo_amador.rds"),
  download_and_prepare_ssurgo(
    aoi_wkt = amador_wkt,
    properties = c("sandtotal", "claytotal", "silttotal", "dbovendry", "ph1to1h2o",
                   "cec7", "om", "wthirdbar", "wfifteenbar"),
    max_depth = 150,
    verbose = TRUE
  )
)

# ---------------------------------------------------------------------------
# 2. A real mukey for the Amador-area AOI, and get_aws_data_by_mukey() for it
#    (vignette 2's raw input)
# ---------------------------------------------------------------------------
mukey_info <- tryCatch(
  process_aoi_and_get_mukeys_working(amador_wkt, verbose = TRUE),
  error = function(e) {
    cat("mukey lookup FAILED:", conditionMessage(e), "\n")
    NULL
  }
)

real_mukey <- if (!is.null(mukey_info)) mukey_info$mukey_list[1] else NA_character_
cat("Using real mukey:", real_mukey, "\n")

depth_sim_mukey_data <- if (!is.na(real_mukey)) {
  save_step(
    "depth_sim_mukey_data",
    file.path(extdata_dir, "depth_sim_mukey_data.rds"),
    list(mukey = real_mukey, data = get_aws_data_by_mukey(real_mukey))
  )
} else NULL

# ---------------------------------------------------------------------------
# 3. Full simulate_profile_depths_by_mukey() result for that real mukey
#    (vignette 2's demonstration output - this function always live-fetches
#    internally, both SDA and OSD, so its final output is cached rather than
#    re-run at vignette build time)
# ---------------------------------------------------------------------------
if (!is.na(real_mukey)) {
  save_step(
    "depth_sim_profiles_amador",
    file.path(extdata_dir, "depth_sim_profiles_amador.rds"),
    simulate_profile_depths_by_mukey(real_mukey, n_simulations = 25, seed = 123)
  )
}

# ---------------------------------------------------------------------------
# 4. ROSETTA / AWS data (vignette 3) - built from the real Amador SSURGO download
# ---------------------------------------------------------------------------
if (!is.null(ssurgo_amador)) {
  ssurgo_df <- ssurgo_amador$ssurgo_data

  rosetta_input <- data.frame(
    compname = ssurgo_df$compname,
    cokey = ssurgo_df$cokey,
    hzdept_r = ssurgo_df$hzdept_r,
    hzdepb_r = ssurgo_df$hzdepb_r,
    sand_total = ssurgo_df$sandtotal_r,
    silt_total = ssurgo_df$silttotal_r,
    clay_total = ssurgo_df$claytotal_r,
    bulk_density_third_bar = ssurgo_df$dbovendry_r,
    water_retention_third_bar = ssurgo_df$wthirdbar_r,
    water_retention_15_bar = ssurgo_df$wfifteenbar_r,
    stringsAsFactors = FALSE
  )
  rosetta_input <- rosetta_input[stats::complete.cases(rosetta_input), ]
  rosetta_input <- utils::head(rosetta_input, 15)  # keep the cached example small
  cat("ROSETTA input rows:", nrow(rosetta_input), "\n")

  if (nrow(rosetta_input) > 0) {
    aws_result <- save_step(
      "aws_rosetta_example",
      file.path(extdata_dir, "aws_rosetta_example.rds"),
      list(input = rosetta_input, result = calculate_aws_df(rosetta_input))
    )
  } else {
    cat("No complete ROSETTA input rows available - skipping AWS example.\n")
  }
}

# ---------------------------------------------------------------------------
# 5. Raster fusion: SSURGO x SOLUS for Salinas Valley, single property (vignette 4)
# ---------------------------------------------------------------------------
Sys.unsetenv("PROJ_LIB")
aoi <- terra::vect(salinas_wkt, crs = "epsg:4326")
aoi <- terra::project(aoi, "epsg:5070")

# property_to_sim_column() (R/ssurgo-simulation.R) only recognizes a fixed set of short
# property-id codes (clay/sand/silt/ph/bulk_density/soc/cec/rock_fragments) for the SSURGO
# prior side - "claytotal" (the SSURGO/SOLUS long-form name, correct for solus_variable) isn't
# one of them and errors here.
property_config <- list(id = "clay", solus_variable = "claytotal", dist = "normal")
fusion_clay <- save_step(
  "fusion_clay_salinas",
  file.path(extdata_dir, "fusion_clay_salinas.rds"),
  wrap_nested_rasters(run_stage1_fusion(aoi, property_config, top_depth = 0, bottom_depth = 5))
)

# ---------------------------------------------------------------------------
# 6. Texture-group fusion (clay/sand/silt), jointly via ILR (vignette 4, bonus section)
# ---------------------------------------------------------------------------
composition_groups <- list(texture = list(members = c("clay", "sand", "silt")))
property_configs <- list(
  clay = list(id = "clay", solus_variable = "claytotal", composition_group = "texture"),
  sand = list(id = "sand", solus_variable = "sandtotal", composition_group = "texture"),
  silt = list(id = "silt", solus_variable = "silttotal", composition_group = "texture")
)
fusion_texture <- save_step(
  "fusion_texture_salinas",
  file.path(extdata_dir, "fusion_texture_salinas.rds"),
  wrap_nested_rasters(run_stage1_fusion(
    aoi, property_configs$clay, top_depth = 0, bottom_depth = 5,
    composition_groups = composition_groups, property_configs = property_configs
  ))
)

cat("\nAll vignette data prep steps complete. Check messages above for any FAILED steps.\n")
