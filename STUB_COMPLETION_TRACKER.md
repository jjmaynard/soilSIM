# soilSIM Stub Completion Tracker

Tracks implementation of the plan: "Port/build out all soilSIM stub functions + priority utility backfill".
Check items off as they are implemented, tested, and verified (one unit of work at a time, not batched).

## Phase 0 - Gap assessment report
- [x] code/soilSIM_stub_gap_assessment.md written

## Phase 1 - Tier 1 (data-infilling.R)
- [x] cross_component_property_interpolation() implemented
- [x] related_property_estimation() implemented
- [x] Wired into infill_missing_property_data() as Strategy 4/5 (required restructuring infill_soil_property()/infill_missing_property_data() so cross-component strategy runs at whole-dataset scope, not per-cokey - see apply_group_fallback_mean()/process_property_group_fallback())
- [x] Stub-contract test (L53-57) removed
- [x] New direct tests added for both functions
- [x] Regression check against existing infill_soil_property() test done (row 3 now recovers via cross-component interpolation; test updated with explanatory comment)
- [x] devtools::check() clean (verified: 0 errors, 0 warnings, 1 pre-existing NOTE fixed via .Rbuildignore)

## Phase 2 - Tier 2 (gp-modeling.R)
- [x] optimize_gp_hyperparameters() real CV implemented (k_fold_gp_cv() shared helper, compares 3 corr specs)
- [x] apply_local_gp_adjustments_with_correlations() delegates to real function
- [x] apply_nrcs_gp_adjustments_with_correlations() implemented - delegates to already-real apply_nrcs_trend_adjustments() (found during implementation; simpler than originally planned from-scratch build)
- [x] simulate_soil_properties() placeholder wired to generate_monte_carlo_realizations()
- [x] flatten_simulation_array_to_long() adapter written + tested
- [x] Stub-contract tests replaced (L92-105, L124-134, L136-143)
- [x] devtools::check() clean (verified together with Phase 3: 0 errors, 0 warnings, 0 notes)

## Phase 3 - Priority missing utilities
- [x] validate_property_config()
- [x] get_rfv_range_category()
- [x] apply_property_constraints()
- [x] create_validation_config() / add_range_rule() / add_relationship_rule() / apply_validation_rules() (exported)
- [x] summarize_unsuitable_horizons()
- [x] Tests added for all 5
- [x] devtools::check() clean (verified together with Phase 2: 0 errors, 0 warnings, 0 notes)

## Phase 4 - Tier 3 (validation-diagnostics.R, 27 stubs)
- [x] 4a Monte Carlo/distribution cluster (7 fns)
- [x] 4b Correlation-matrix quality cluster (5 fns)
- [x] 4c GP performance cluster (5 fns)
- [x] 4d Cross-property/pedological realism cluster (5 fns)
- [x] 4e Workflow scoring cluster (2 fns) + shared k_fold_gp_cv() helper (reused from Phase 2)
- [x] 4f Report generation cluster (3 fns) + Suggests: rmarkdown, tinytex added
- [x] Combined stub-contract test (L130-145) removed
- [x] Tests written for all 27 (organized into 6 cluster blocks)
- [x] devtools::check() clean (found + fixed 2 real test failures: validate_horizon_characteristics test used mismatched property-name convention vs get_realistic_property_ranges(); render_value_as_lines() mis-detected a classed-list object as a plain list. Also found + fixed an additional untagged stub, assess_simulation_quality(), not in original inventory, while implementing generate_simulation_diagnostics() in Phase 5 - see Phase 5 notes.)

## Phase 5 - Tier 4 (untagged stubs)
- [x] monte-carlo.R: 8 functions implemented (validate_component_data, extract_component_parameters, apply_composition_constraints, assess_component_quality, calculate_simulation_quality_metrics, generate_simulation_diagnostics, summarize_distributions, validate_distribution_setup) + assess_simulation_quality() (untagged, found during this pass, not in original inventory - also implemented for real)
- [x] statistics.R: 7 functions implemented
- [x] Regression check: test-monte-carlo.R quality_assessment/diagnostics shape - confirmed no existing test asserted on old hardcoded values, zero regression risk
- [x] Regression check: test-statistics.R outlier-detection caller shape - confirmed no existing test asserted on old hardcoded values, zero regression risk
- [x] Tests added for all 16 (15 planned + assess_simulation_quality)
- [x] devtools::check() clean (confirmed: Status: OK - 0 errors, 0 warnings, 0 notes)

## Final gate
- [x] Full integrated devtools::check() across all phases: 0 errors, 0 warnings, 0 notes (confirmed via the Phase 5 check run, which covers the fully-merged state of all 5 phases)
