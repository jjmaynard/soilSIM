# soilSIM Stub Completion Tracker

Tracks implementation of the plan: “Port/build out all soilSIM stub
functions + priority utility backfill”. Check items off as they are
implemented, tested, and verified (one unit of work at a time, not
batched).

## Phase 0 - Gap assessment report

code/soilSIM_stub_gap_assessment.md written

## Phase 1 - Tier 1 (data-infilling.R)

cross_component_property_interpolation() implemented

related_property_estimation() implemented

Wired into infill_missing_property_data() as Strategy 4/5 (required
restructuring infill_soil_property()/infill_missing_property_data() so
cross-component strategy runs at whole-dataset scope, not per-cokey -
see apply_group_fallback_mean()/process_property_group_fallback())

Stub-contract test (L53-57) removed

New direct tests added for both functions

Regression check against existing infill_soil_property() test done (row
3 now recovers via cross-component interpolation; test updated with
explanatory comment)

devtools::check() clean (verified: 0 errors, 0 warnings, 1 pre-existing
NOTE fixed via .Rbuildignore)

## Phase 2 - Tier 2 (gp-modeling.R)

optimize_gp_hyperparameters() real CV implemented (k_fold_gp_cv() shared
helper, compares 3 corr specs)

apply_local_gp_adjustments_with_correlations() delegates to real
function

apply_nrcs_gp_adjustments_with_correlations() implemented - delegates to
already-real apply_nrcs_trend_adjustments() (found during
implementation; simpler than originally planned from-scratch build)

simulate_soil_properties() placeholder wired to
generate_monte_carlo_realizations()

flatten_simulation_array_to_long() adapter written + tested

Stub-contract tests replaced (L92-105, L124-134, L136-143)

devtools::check() clean (verified together with Phase 3: 0 errors, 0
warnings, 0 notes)

## Phase 3 - Priority missing utilities

validate_property_config()

get_rfv_range_category()

apply_property_constraints()

create_validation_config() / add_range_rule() / add_relationship_rule()
/ apply_validation_rules() (exported)

summarize_unsuitable_horizons()

Tests added for all 5

devtools::check() clean (verified together with Phase 2: 0 errors, 0
warnings, 0 notes)

## Phase 4 - Tier 3 (validation-diagnostics.R, 27 stubs)

4a Monte Carlo/distribution cluster (7 fns)

4b Correlation-matrix quality cluster (5 fns)

4c GP performance cluster (5 fns)

4d Cross-property/pedological realism cluster (5 fns)

4e Workflow scoring cluster (2 fns) + shared k_fold_gp_cv() helper
(reused from Phase 2)

4f Report generation cluster (3 fns) + Suggests: rmarkdown, tinytex
added

Combined stub-contract test (L130-145) removed

Tests written for all 27 (organized into 6 cluster blocks)

devtools::check() clean (found + fixed 2 real test failures:
validate_horizon_characteristics test used mismatched property-name
convention vs get_realistic_property_ranges(); render_value_as_lines()
mis-detected a classed-list object as a plain list. Also found + fixed
an additional untagged stub, assess_simulation_quality(), not in
original inventory, while implementing generate_simulation_diagnostics()
in Phase 5 - see Phase 5 notes.)

## Phase 5 - Tier 4 (untagged stubs)

monte-carlo.R: 8 functions implemented (validate_component_data,
extract_component_parameters, apply_composition_constraints,
assess_component_quality, calculate_simulation_quality_metrics,
generate_simulation_diagnostics, summarize_distributions,
validate_distribution_setup) + assess_simulation_quality() (untagged,
found during this pass, not in original inventory - also implemented for
real)

statistics.R: 7 functions implemented

Regression check: test-monte-carlo.R quality_assessment/diagnostics
shape - confirmed no existing test asserted on old hardcoded values,
zero regression risk

Regression check: test-statistics.R outlier-detection caller shape -
confirmed no existing test asserted on old hardcoded values, zero
regression risk

Tests added for all 16 (15 planned + assess_simulation_quality)

devtools::check() clean (confirmed: Status: OK - 0 errors, 0 warnings, 0
notes)

## Final gate

Full integrated devtools::check() across all phases: 0 errors, 0
warnings, 0 notes (confirmed via the Phase 5 check run, which covers the
fully-merged state of all 5 phases)
