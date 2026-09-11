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

# --- P4: systematic renames ---------------------------------------------

#' @rdname soilSIM-deprecated
#' @export
validate_correlation_preservation_diagnostics <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_correlation_preservation_diagnostics()", "diagnose_correlation_preservation()")
  diagnose_correlation_preservation(...)
}

#' @rdname soilSIM-deprecated
#' @export
process_component_data_working_compatible <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "process_component_data_working_compatible()", "process_ssurgo_components()")
  process_ssurgo_components(...)
}

#' @rdname soilSIM-deprecated
#' @export
process_horizon_data_working_compatible <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "process_horizon_data_working_compatible()", "process_ssurgo_horizons()")
  process_ssurgo_horizons(...)
}

#' @rdname soilSIM-deprecated
#' @export
create_ssurgo_property_lookup_working <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "create_ssurgo_property_lookup_working()", "build_ssurgo_property_lookup()")
  build_ssurgo_property_lookup(...)
}

#' @rdname soilSIM-deprecated
#' @export
clean_property_data_ssurgo_compatible <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "clean_property_data_ssurgo_compatible()", "clean_ssurgo_property_data()")
  clean_ssurgo_property_data(...)
}

#' @rdname soilSIM-deprecated
#' @export
validate_within_depth_correlations <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_within_depth_correlations()", "diagnose_within_depth_correlations()")
  diagnose_within_depth_correlations(...)
}

#' @rdname soilSIM-deprecated
#' @export
generate_validation_report_ssurgo <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "generate_validation_report_ssurgo()", "diagnose_ssurgo_download()")
  diagnose_ssurgo_download(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_statistical_analysis_defaults <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_statistical_analysis_defaults()", "default_statistics_config()")
  default_statistics_config(...)
}

#' @rdname soilSIM-deprecated
#' @export
create_infill_compatible_dataset <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "create_infill_compatible_dataset()", "prepare_ssurgo_for_infill()")
  prepare_ssurgo_for_infill(...)
}

#' @rdname soilSIM-deprecated
#' @export
validate_correlation_structures <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_correlation_structures()", "diagnose_correlation_structures()")
  diagnose_correlation_structures(...)
}

#' @rdname soilSIM-deprecated
#' @export
fetch_solus_low_pred_high_multi <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "fetch_solus_low_pred_high_multi()", "fetch_solus_low_pred_high_multiproperty()")
  fetch_solus_low_pred_high_multiproperty(...)
}

#' @rdname soilSIM-deprecated
#' @export
validate_distribution_fidelity <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_distribution_fidelity()", "diagnose_distribution_fidelity()")
  diagnose_distribution_fidelity(...)
}

#' @rdname soilSIM-deprecated
#' @export
calculate_confidence_intervals <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "calculate_confidence_intervals()", "compute_confidence_intervals()")
  compute_confidence_intervals(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_property_contextual_ranges <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_property_contextual_ranges()", "property_contextual_ranges()")
  property_contextual_ranges(...)
}

#' @rdname soilSIM-deprecated
#' @export
validate_soil_science_realism <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_soil_science_realism()", "diagnose_soil_science_realism()")
  diagnose_soil_science_realism(...)
}

#' @rdname soilSIM-deprecated
#' @export
validate_gp_model_performance <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_gp_model_performance()", "diagnose_gp_performance()")
  diagnose_gp_performance(...)
}

#' @rdname soilSIM-deprecated
#' @export
assess_cholesky_decomposition <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "assess_cholesky_decomposition()", "diagnose_cholesky()")
  diagnose_cholesky(...)
}

#' @rdname soilSIM-deprecated
#' @export
fetch_solus_percentiles_multi <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "fetch_solus_percentiles_multi()", "fetch_solus_percentiles_multiproperty()")
  fetch_solus_percentiles_multiproperty(...)
}

#' @rdname soilSIM-deprecated
#' @export
calculate_saxton_rawls_single <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "calculate_saxton_rawls_single()", "compute_saxton_rawls()")
  compute_saxton_rawls(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_appropriate_distributions <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_appropriate_distributions()", "distributions_for_properties()")
  distributions_for_properties(...)
}

#' @rdname soilSIM-deprecated
#' @export
validate_monte_carlo_quality <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_monte_carlo_quality()", "diagnose_simulation()")
  diagnose_simulation(...)
}

#' @rdname soilSIM-deprecated
#' @export
horizon_name_property_infill <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "horizon_name_property_infill()", "infill_property_by_horizon_name()")
  infill_property_by_horizon_name(...)
}

#' @rdname soilSIM-deprecated
#' @export
calculate_summary_statistics <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "calculate_summary_statistics()", "compute_summary_statistics()")
  compute_summary_statistics(...)
}

#' @rdname soilSIM-deprecated
#' @export
assess_property_constraints <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "assess_property_constraints()", "diagnose_property_constraints()")
  diagnose_property_constraints(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_default_property_config <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_default_property_config()", "default_property_config()")
  default_property_config(...)
}

#' @rdname soilSIM-deprecated
#' @export
validate_gp_model_workflow <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_gp_model_workflow()", "diagnose_gp_models()")
  diagnose_gp_models(...)
}

#' @rdname soilSIM-deprecated
#' @export
assess_simulation_coverage <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "assess_simulation_coverage()", "diagnose_simulation_coverage()")
  diagnose_simulation_coverage(...)
}

#' @rdname soilSIM-deprecated
#' @export
assess_depth_trend_realism <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "assess_depth_trend_realism()", "diagnose_depth_trend_realism()")
  diagnose_depth_trend_realism(...)
}

#' @rdname soilSIM-deprecated
#' @export
mukey_texture_draws_lookup <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "mukey_texture_draws_lookup()", "lookup_mukey_texture_draws()")
  lookup_mukey_texture_draws(...)
}

#' @rdname soilSIM-deprecated
#' @export
bayes_update_normal_normal <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "bayes_update_normal_normal()", "fuse_normal_normal()")
  fuse_normal_normal(...)
}

#' @rdname soilSIM-deprecated
#' @export
sim_component_compositions <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "sim_component_compositions()", "simulate_component_compositions_batch()")
  simulate_component_compositions_batch(...)
}

#' @rdname soilSIM-deprecated
#' @export
calculate_complexity_score <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "calculate_complexity_score()", "compute_complexity_score()")
  compute_complexity_score(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_predefined_properties <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_predefined_properties()", "predefined_properties()")
  predefined_properties(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_default_configuration <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_default_configuration()", "default_config()")
  default_config(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_monte_carlo_defaults <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_monte_carlo_defaults()", "default_monte_carlo_config()")
  default_monte_carlo_config(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_available_properties <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_available_properties()", "available_properties()")
  available_properties(...)
}

#' @rdname soilSIM-deprecated
#' @export
validate_gp_predictions <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "validate_gp_predictions()", "diagnose_gp_predictions()")
  diagnose_gp_predictions(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_validation_defaults <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_validation_defaults()", "default_diagnostics_config()")
  default_diagnostics_config(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_rfv_range_category <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_rfv_range_category()", "rfv_range_category()")
  rfv_range_category(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_aws_data_by_mukey <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_aws_data_by_mukey()", "fetch_ssurgo_aws_data()")
  fetch_ssurgo_aws_data(...)
}

#' @rdname soilSIM-deprecated
#' @export
get_default_synonyms <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "get_default_synonyms()", "default_property_synonyms()")
  default_property_synonyms(...)
}

#' @rdname soilSIM-deprecated
#' @export
hz_quant_prob_mukey <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "hz_quant_prob_mukey()", "compute_mukey_horizon_quantiles()")
  compute_mukey_horizon_quantiles(...)
}

#' @rdname soilSIM-deprecated
#' @export
mukey_draws_lookup <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "mukey_draws_lookup()", "lookup_mukey_draws()")
  lookup_mukey_draws(...)
}

#' @rdname soilSIM-deprecated
#' @export
bayesian_update <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "bayesian_update()", "update_prior()")
  update_prior(...)
}

#' @rdname soilSIM-deprecated
#' @export
calculate_mode <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "calculate_mode()", "compute_mode()")
  compute_mode(...)
}

#' @rdname soilSIM-deprecated
#' @export
bayes_fuse <- function(...) {
  lifecycle::deprecate_soft("0.2.0", "bayes_fuse()", "fuse_distribution()")
  fuse_distribution(...)
}

# --- D2: legacy_exceptions renames -------------------------------------

#' @rdname soilSIM-deprecated
#' @export
beta_to_moments <- function(...) {
  lifecycle::deprecate_soft("0.2.1", "beta_to_moments()", "convert_beta_to_moments()")
  convert_beta_to_moments(...)
}

#' @rdname soilSIM-deprecated
#' @export
gamma_to_moments <- function(...) {
  lifecycle::deprecate_soft("0.2.1", "gamma_to_moments()", "convert_gamma_to_moments()")
  convert_gamma_to_moments(...)
}

#' @rdname soilSIM-deprecated
#' @export
moments_to_beta <- function(...) {
  lifecycle::deprecate_soft("0.2.1", "moments_to_beta()", "convert_moments_to_beta()")
  convert_moments_to_beta(...)
}

#' @rdname soilSIM-deprecated
#' @export
moments_to_gamma <- function(...) {
  lifecycle::deprecate_soft("0.2.1", "moments_to_gamma()", "convert_moments_to_gamma()")
  convert_moments_to_gamma(...)
}

#' @rdname soilSIM-deprecated
#' @export
lognormal_to_normal_params <- function(...) {
  lifecycle::deprecate_soft("0.2.1", "lognormal_to_normal_params()", "convert_lognormal_to_normal()")
  convert_lognormal_to_normal(...)
}

#' @rdname soilSIM-deprecated
#' @export
normal_to_lognormal_params <- function(...) {
  lifecycle::deprecate_soft("0.2.1", "normal_to_lognormal_params()", "convert_normal_to_lognormal()")
  convert_normal_to_lognormal(...)
}

#' @rdname soilSIM-deprecated
#' @export
remarginalized_awc <- function(...) {
  lifecycle::deprecate_soft("0.2.1", "remarginalized_awc()", "remarginalize_awc()")
  remarginalize_awc(...)
}
