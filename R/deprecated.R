#' Deprecated functions in soilSIM
#'
#' These names were renamed in the 0.2.0 naming re-architecture. Each still
#' works and forwards to its replacement, but emits a soft deprecation warning.
#' Timeline: 0.2.0 soft-deprecate, 0.3.0 warn, 0.4.0 remove.
#'
#' @name soilSIM-deprecated
#' @keywords internal
#' @importFrom lifecycle deprecate_soft
NULL

# --- P3: Tier-1 entry-point renames -------------------------------------

#' @rdname soilSIM-deprecated
#' @export
download_and_prepare_ssurgo <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "download_and_prepare_ssurgo()", "fetch_ssurgo_data()")
  fetch_ssurgo_data(...)
}

#' @rdname soilSIM-deprecated
#' @export
infill_soil_data <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "infill_soil_data()", "infill_ssurgo_data()")
  infill_ssurgo_data(...)
}

#' @rdname soilSIM-deprecated
#' @export
generate_monte_carlo_realizations <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "generate_monte_carlo_realizations()", "simulate_monte_carlo()")
  simulate_monte_carlo(...)
}

#' @rdname soilSIM-deprecated
#' @export
build_stratified_gp_models <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "build_stratified_gp_models()", "fit_depth_gp_models()")
  fit_depth_gp_models(...)
}

#' @rdname soilSIM-deprecated
#' @export
adjust_multivariate_depthwise_GP <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "adjust_multivariate_depthwise_GP()", "adjust_simulation_depthwise()")
  adjust_simulation_depthwise(...)
}

#' @rdname soilSIM-deprecated
#' @export
integrate_monte_carlo_with_gp <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "integrate_monte_carlo_with_gp()", "apply_depth_gp_to_simulation()")
  apply_depth_gp_to_simulation(...)
}

#' @rdname soilSIM-deprecated
#' @export
sim_component_comp <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "sim_component_comp()", "simulate_component_composition()")
  simulate_component_composition(...)
}

#' @rdname soilSIM-deprecated
#' @export
calculate_aws_df <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "calculate_aws_df()", "compute_aws()")
  compute_aws(...)
}

#' @rdname soilSIM-deprecated
#' @export
run_stage1_fusion <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "run_stage1_fusion()", "run_fusion()")
  run_fusion(...)
}

#' @rdname soilSIM-deprecated
#' @export
run_stage1_fusion_group <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "run_stage1_fusion_group()", "run_fusion_group()")
  run_fusion_group(...)
}

#' @rdname soilSIM-deprecated
#' @export
run_stage1_fusion_multi <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "run_stage1_fusion_multi()", "run_fusion_multiproperty()")
  run_fusion_multiproperty(...)
}

#' @rdname soilSIM-deprecated
#' @export
validate_complete_workflow <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_complete_workflow()", "diagnose_workflow()")
  diagnose_workflow(...)
}
