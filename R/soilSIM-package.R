#' soilSIM: Multi-Source Soil Property Simulation and Bayesian Data Fusion
#'
#' A data-source-agnostic framework for simulating soil properties and fusing
#' estimates from multiple soil data products. Provides percentile-triplet and
#' arbitrary-percentile distribution fitting, correlated Monte Carlo
#' simulation (including compositional texture handling via isometric
#' log-ratio transforms), Gaussian-process depth-trend modeling, descriptive
#' statistics/correlation analysis, and Bayesian updating/fusion in both
#' scalar and raster-native forms. Ships with SSURGO and SOLUS100 data-source
#' adapters wired into the fusion pipeline as prior and likelihood sources
#' respectively; the simulation and fusion core is designed to accommodate
#' additional gridded soil data sources (e.g. HWSD, SoilGrids) as further
#' adapters.
#'
#' @section Overview:
#' The package is split into a generic simulation/fusion **core** - percentile
#' fitting, Monte Carlo simulation, GP depth modeling, and Bayesian fusion -
#' that operates on any low/representative/high triplet or prior/likelihood
#' pair regardless of where the numbers came from, and a set of
#' source-specific **adapters** that acquire, process, and shape a given soil
#' data product's output into that generic form. SSURGO
#' (`ssurgo-acquisition.R`, `ssurgo-processing.R`, `ssurgo-simulation.R`) and
#' SOLUS100 (`solus-simulation.R`) are the two adapters shipped today; new
#' data sources are intended to plug in as additional adapters alongside
#' them, reusing the same core rather than duplicating it.
#'
#' @section Data acquisition & processing:
#' The SSURGO adapter: [download_ssurgo_tabular()]/
#' [download_and_prepare_ssurgo()] acquire and cache tabular SSURGO data for
#' an area of interest; [process_ssurgo_data()] cleans and standardizes it
#' into modeling-ready horizon/component tables;
#' [process_soil_properties_comprehensive()]/[infill_soil_property()] fill
#' missing property values using pedologically informed recovery strategies.
#'
#' @section Statistics & diagnostics:
#' [analyze_soil_statistics()] characterizes processed data (correlations,
#' distribution fitting, outliers, summary statistics);
#' [validate_complete_workflow()]/[generate_validation_report()] provide
#' end-to-end quality assurance across the full simulation workflow.
#'
#' @section Distribution fitting (core):
#' [fit_percentile_triplet()]/[quantile_from_fit()] fit and sample closed-form
#' distributions (Normal, Beta, metalog, triangular, linear-CDF) from a
#' low/representative/high triplet; [simulate_from_percentiles()] generalizes
#' this to an arbitrary number of percentiles via linear-CDF/spline/KDE
#' methods; [ilr_forward()]/[ilr_inverse()] provide the isometric log-ratio
#' transforms used for compositional texture data.
#'
#' @section Monte Carlo simulation (core):
#' [generate_monte_carlo_realizations()] is the master pipeline for correlated
#' Monte Carlo simulation of soil properties from percentile inputs;
#' [simulate_correlated_properties()] is its Cholesky-copula core simulator.
#'
#' @section GP depth modeling & multivariate adjustment (core):
#' [build_stratified_gp_models()]/[simulate_soil_properties()] fit and apply
#' Gaussian-process depth-trend models; [integrate_monte_carlo_with_gp()]
#' combines Monte Carlo output with GP depth trends across a whole dataset
#' while preserving cross-property correlation structure.
#'
#' @section Profile, component & depth simulation:
#' [sim_component_comp()]/[simulate_cokey_generalized()] simulate component
#' composition and per-cokey properties;
#' [simulate_and_perturb_soil_profiles()]/[simulate_profile_depths_by_mukey()]
#' simulate horizon depths and thicknesses for whole soil profiles.
#'
#' @section AWS / Van Genuchten modeling:
#' [van_genuchten()]/[simulate_vg_aws()]/[calculate_aws_df()] estimate
#' available water storage from ROSETTA pedotransfer parameters.
#'
#' @section Bayesian updating (scalar, core):
#' [bayes_update_normal_normal()] and family provide closed-form same-family
#' fusion of a prior and observed data; [bayesian_update()] provides a fully
#' general grid-KDE fusion for arbitrary distributions;
#' [fuse_texture_group_from_triplets()] fuses compositional texture data
#' jointly in ILR space.
#'
#' @section Multi-source raster fusion pipeline:
#' [fuse_property_adaptive()]/[fuse_texture_group()] (in `raster-fusion.R`)
#' are the generic raster-native fusion core, combining prior and likelihood
#' `SpatRaster` inputs cell-by-cell over an area of interest; the SSURGO
#' adapter ([simulate_ssurgo_mapunit_draws()]) supplies the prior side and the
#' SOLUS100 adapter ([fetch_solus_percentiles()]) supplies the likelihood
#' side, with [cache_get()]/[cache_set()] providing disk caching for both.
#' Designed so additional data-source adapters can supply either side without
#' changes to the fusion core.
#'
#' @section Utilities:
#' [is_unsuitable()]/[validate_data_quality()]/[read_soil_data()]/
#' [write_soil_data()] and a shared logging/configuration framework used
#' throughout the package.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats approxfun cor cov density dnorm median na.omit qbeta
#' @importFrom stats qnorm qunif quantile rbeta rnorm runif sd var
#' @importFrom stats IQR approx complete.cases mad pnorm qt shapiro.test
#' @importFrom stats ecdf setNames weighted.mean
#' @importFrom Hmisc rcorr
#' @importFrom Matrix nearPD
#' @importFrom fitdistrplus fitdist gofstat
#' @importFrom dplyr group_by mutate ungroup across all_of summarise_all
#' @importFrom rlang .data :=
#' @importFrom tidyr pivot_longer everything
#' @importFrom jsonlite fromJSON write_json
#' @importFrom readr read_csv read_tsv write_csv
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils head installed.packages read.csv read.delim
## usethis namespace: end
NULL

#' Default value for NULL
#'
#' @param x value to test
#' @param y fallback used when `x` is `NULL`
#' @return `x` if not `NULL`, otherwise `y`
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Package load hook
#'
#' Guards against a documented, real crash mode in the raster fusion pipeline
#' (`R/raster-fusion.R`, `R/raster-cache.R`): a stray system `PROJ_LIB` environment variable
#' pointing at a version-incompatible external GDAL/PROJ installation silently segfaults
#' `terra::resample()` (used to align SSURGO/SOLUS grids before fusion) -
#' `code_ref/reanalysis-platform/HANDOFF_NOTES.md` documents this exact failure and its fix.
#' Idempotent and a no-op when `PROJ_LIB` isn't set - safe to run unconditionally.
#'
#' @param libname,pkgname Standard `.onLoad()` arguments, unused.
#' @noRd
.onLoad <- function(libname, pkgname) {
  if (nzchar(Sys.getenv("PROJ_LIB", unset = ""))) {
    Sys.unsetenv("PROJ_LIB")
  }
}

# Column/variable names referenced via bare NSE inside dplyr verbs
# (mutate()/filter()/group_by()/summarise()/etc.) across the migrated
# mod01/02/06/07/08 sources, which predate this package's later, more
# consistent use of the .data[[...]] pronoun. R CMD check's static analysis
# cannot tell these apart from genuine undefined globals; this is the
# standard way to declare them as intentional NSE symbols rather than
# rewriting every call site to .data[[...]] style (a much larger,
# behavior-risking change out of scope for a migration).
utils::globalVariables(c(
  "adequate_depth_range", "adequate_observations", "adequate_profiles", "chkey",
  "cokey", "compname", "component", "comppct_r", "depth", "depth_bin", "depth_range",
  "fallback_group", "final_group", "fragsize_r", "gp_model_group",
  "has_reskind_restriction", "horizon_suggests_restriction",
  "hz_below_restriction", "hzdepb_r", "hzdept_r", "is_cemented", "lieutex", "mean_val",
  "mean_value", "missing_required_properties", "n_obs", "n_observations",
  "n_profiles", "pct_violations", "prop_available", "prop_value", "property",
  "resdept_l", "resdept_r", "reskind", "restriction_top", "score", "sim_cokey",
  "simulation_number", "soil_group", "taxclname", "taxgrtgroup", "taxorder",
  "taxpartsize", "taxsuborder", "test_group", "texcl", "texture_suggests_restriction",
  "total_obs", "unsuitable_horizon", "value",
  # R/depth-simulation.R (SSURGO horizon-depth simulation)
  "mukey", "id", "hzname", "sim_comppct", "hzdept_l", "hzdept_h", "hzdepb_l",
  "hzdepb_h", "hzthk_l", "rfv_l", "rfv_r", "rfv_h", "genhz", "distinctness",
  "bound_sd", "top", "bottom", "thickness_sd", "Thickness",
  # R/aws-simulation.R (van Genuchten / ROSETTA-based AWS modeling)
  "alpha", "n", "theta_r", "theta_s", "contributing_fraction", "variable",
  # R/property-simulation.R (component-composition / correlated-triangular simulation)
  "Depth", "comppct_l", "comppct_h",
  # R/ssurgo-processing.R::hz_quant_prob_mukey()
  "sand", "silt", "clay", "prob", "txt_class"
))
