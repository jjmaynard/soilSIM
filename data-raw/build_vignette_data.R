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

# Fixed seed so the cached fusion .rds objects regenerate deterministically (opt-in via the `seed` argument). Determinism is still conditional on the
# live SSURGO/SOLUS data itself being unchanged between runs.
VIGNETTE_SEED <- 20260903L

# property_to_sim_column() (R/ssurgo-simulation.R) only recognizes a fixed set of short
# property-id codes (clay/sand/silt/ph/bulk_density/soc/cec/rock_fragments) for the SSURGO
# prior side - "claytotal" (the SSURGO/SOLUS long-form name, correct for solus_variable) isn't
# one of them and errors here.
property_config <- list(id = "clay", solus_variable = "claytotal", dist = "normal")
fusion_clay <- save_step(
  "fusion_clay_salinas",
  file.path(extdata_dir, "fusion_clay_salinas.rds"),
  wrap_nested_rasters(run_stage1_fusion(aoi, property_config, top_depth = 0, bottom_depth = 5,
                                        seed = VIGNETTE_SEED))
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
  # Must be run_stage1_fusion_group() directly, not run_stage1_fusion() - the latter's own
  # dispatch (see its composition_group branch in R/raster-fusion.R) slices its return down
  # to a single requested member (group_result[[property_config$id]]) even when
  # composition_groups/property_configs describe a full group. That silently produced a
  # ONE-member (clay only: posterior/dist/route/...) result here instead of the 3-member
  # clay/sand/silt list the vignette's "Step 2" section actually consumes
  # (fusion_texture$clay/$sand/$silt) - a genuine, previously-undiscovered bug: this script
  # didn't match its own vignette's displayed "real call" comment, which already showed
  # run_stage1_fusion_group() correctly.
  wrap_nested_rasters(run_stage1_fusion_group(
    aoi, "texture", composition_groups, property_configs, top_depth = 0, bottom_depth = 5,
    seed = VIGNETTE_SEED
  ))
)

# ---------------------------------------------------------------------------
# 7. Per-pixel fused ensemble + Saxton-Rawls AWC inputs, Salinas Valley
#    (vignette: raster-fusion-perpixel-awc.Rmd). The bridge's re-marginalization + AWC steps
#    themselves are deterministic and offline, so only the network-dependent parts are cached:
#    (a) the joint SSURGO Monte Carlo ensemble, (b) one SOLUS x SSURGO fused posterior per
#    Saxton-Rawls input property per depth window.
# ---------------------------------------------------------------------------
pp_windows <- list(c(0, 5), c(5, 15), c(15, 30))

# sim-column id -> fetchSOLUS() variable, for the Saxton-Rawls AWC inputs (the ensemble cache is
# trimmed to just these columns since this data set is for the AWC vignette only).
pp_solus <- c(sand_total = "sandtotal", silt_total = "silttotal", clay_total = "claytotal",
              db = "dbovendry", soc = "soc", rfv = "fragvol")

pp_ensemble <- save_step(
  "perpixel_ensemble_salinas",
  file.path(extdata_dir, "perpixel_ensemble_salinas.rds"),
  {
    ens <- extract_mukey_joint_ensemble(aoi, pp_windows, seed = VIGNETTE_SEED)
    if (!is.null(ens)) {
      keep_cols <- intersect(names(pp_solus), ens$properties)
      # Drop tiny map units (< 100 replicates - a component barely overlapping the AOI): with one
      # such mukey `n_kept` = min(...) collapses to a handful of realizations for the whole product.
      keep <- vapply(ens$by_mukey, function(e) nrow(e$windows[[1]]) >= 100L, logical(1))
      ens$by_mukey <- ens$by_mukey[keep]
      # The vignette caps at n_out = 200, so trim each mukey to <= 250 replicates (evenly spaced,
      # row-aligned across windows + replicate_key), and to just the Saxton-Rawls input columns -
      # keeps the cached .rds ~0.3 MB instead of ~3 MB.
      cap <- 250L
      ens$by_mukey <- lapply(ens$by_mukey, function(e) {
        n <- nrow(e$windows[[1]])
        idx <- if (n <= cap) seq_len(n) else unique(round(seq(1, n, length.out = cap)))
        e$windows <- lapply(e$windows, function(m) m[idx, keep_cols, drop = FALSE])
        e$replicate_key <- e$replicate_key[idx, , drop = FALSE]
        e
      })
      ens$properties <- keep_cols
      ens$mukey_raster <- terra::wrap(ens$mukey_raster)  # only the raster field
    }
    ens
  }
)

pp_posteriors <- save_step(
  "perpixel_posteriors_salinas",
  file.path(extdata_dir, "perpixel_posteriors_salinas.rds"),
  {
    # One SSURGO simulation for all pp_solus properties x all pp_windows (was a
    # run_stage1_fusion() loop that re-ran the full simulation once per property x window).
    pp_configs <- stats::setNames(
      lapply(names(pp_solus), function(nm) list(id = nm, solus_variable = pp_solus[[nm]], dist = "auto")),
      names(pp_solus)
    )
    out <- run_stage1_fusion_multi(aoi, pp_configs, pp_windows, seed = VIGNETTE_SEED,
                                   simplify = TRUE)
    # Match the old loop's shape: drop NULL leaves (a failed window), then drop any property left
    # with no windows at all (SOLUS variable unavailable).
    out <- lapply(out, function(pw) pw[!vapply(pw, is.null, logical(1))])
    out <- out[lengths(out) > 0]
    if (length(out) == 0) NULL else wrap_nested_rasters(out)
  }
)

cat("\nAll vignette data prep steps complete. Check messages above for any FAILED steps.\n")
