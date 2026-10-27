# Supplementary Data 1: data dictionary and caption manifest

## Observation-count clarification for release1.0.1

The principal comparison contains 642 site–date observations from four site-seasons at three towers, spanning 489 distinct calendar dates. The separately evaluated Mead hindcast contains 437 site–date observations from three site-seasons at two previously represented towers, spanning 291 distinct calendar dates. Observations on the same calendar date at different towers are separate site–date records, not independent dates or independent replications. The 480 later-year hindcast training observations reuse principal-population observations; the cohorts are not pooled.

Legacy table fields such as `n_days` and cohort-level retained/evaluation-day counts retain their archived names and numeric values. In these summaries they count site–date observations within the stated cohort/role, and may therefore count a calendar date more than once across towers. A count within a single site-season has its corresponding site-specific daily interpretation. Counts of prediction, model, fitting/scoring-reference and fold roles can repeat the same physical observation; 112,240 scored rows are not 112,240 independent observations.

This clarification changes no field name, schema, row, numerical value, model, fold, reference path or aggregation. The 37 publication-table files remain byte-identical to the prior clean public candidate. Historical source manifests and original estimator-bound protocol hashes remain their as-of provenance authorities.


Supplementary Data 1 contains 37 unchanged UTF-8 CSV source tables. Table numbers are scholarly identifiers, not a file count: Tables S1, S2, S5 and S17 each have a/b files; S3–S4, S6–S16 and S18–S30 each have one file (34 supplementary CSVs). Main Tables 1–2 contribute three further source CSVs, giving 37. The supplementary figures are S1–S7; their numbers do not imply seven additional CSVs.

Each entry lists the original caption, physical data-row count (header excluded), column names in source order and units from the source index. Fields and bytes are unchanged, including source missing-value representations. Do not infer a common denominator across labelled estimands or merge principal, hindcast and contextual populations. VUT/CUT denote dependent supplied nighttime-partitioned reference products. B/C denote uniform/conditional retained-core optical representations; first/second train/test indices retain their source meaning. Signed C-minus-B and conditional contrasts below zero favor C; an interaction sign alone does not indicate improvement. RMSE, MAE and bias use g C m^-2 d^-1; squared-error fields use squared units. Relative percentages and common-BB percentage points are distinct. Counts, Boolean/status fields, dates and identifiers retain their header definitions. The machine-readable manifest records exact SHA-256 and size for every member.

## 1. Table 1 — Principal LOYO mean site-season scores

`tables/Table1_primary_LOYO_scores.csv` — 16 data rows; 16 columns.

Columns: `axis`, `model`, `target`, `seed`, `n_site_seasons`, `n_clusters`, `macro_site_season_rmse`, `median_site_season_rmse`, `macro_site_season_mae`, `macro_site_season_bias`, `macro_site_season_nse`, `n_finite_season_nse`, `macro_cluster_rmse`, `fit_reference`, `score_reference`, `arm`.

Units: GPP RMSE/differences:g C m^-2 d^-1; relative changes:percent; counts as labelled.

## 2. Table 1 — Machine-readable source contrasts supporting Table 1 and Figure 4

`tables/Table1_source_contrasts.csv` — 12 data rows; 17 columns.

Columns: `model`, `variant`, `contrast`, `equal_site_season_first_RMSE`, `equal_site_season_second_RMSE`, `equal_site_season_absolute_difference`, `equal_site_season_relative_percent`, `equal_region_pooled_first_RMSE`, `equal_region_pooled_second_RMSE`, `equal_region_pooled_absolute_difference`, `equal_region_pooled_relative_percent`, `relative_effect_magnitude_shrink_percent`, `equal_region_mean_of_site_season_RMSE_first`, `equal_region_mean_of_site_season_RMSE_second`, `equal_region_mean_of_site_season_absolute_difference`, `equal_region_mean_of_site_season_relative_percent`, `definition`.

Units: GPP RMSE/differences:g C m^-2 d^-1; relative changes:percent; counts as labelled.

## 3. Table 2 — Primary displayed conditional comparisons from the saved B/C models

`tables/Table2_BC_crossprediction_primary.csv` — 8 data rows; 25 columns.

Columns: `cohort`, `axis`, `model`, `fit_reference`, `score_reference`, `estimand`, `R_BB`, `R_BC`, `R_CB`, `R_CC`, `training_given_test_B`, `training_given_test_C`, `test_given_training_B`, `test_given_training_C`, `interaction`, `original_pipeline_contrast`, `n_units`, `n_days`, `maximum_identity_residual`, `training_given_test_B_common_BB_pp`, `training_given_test_C_common_BB_pp`, `test_given_training_B_common_BB_pp`, `test_given_training_C_common_BB_pp`, `interaction_common_BB_pp`, `original_pipeline_contrast_common_BB_pp`.

Units: RMSE and signed contrasts: g C m^-2 d^-1; common-BB normalized fields: percentage points; counts as labelled.

## 4. Table S1a — Principal LOYO scores and input contrasts for all five fixed families — scores/metrics

`tables/TableS1a_all_LOYO_models.csv` — 40 data rows; 16 columns.

Columns: `axis`, `model`, `target`, `seed`, `n_site_seasons`, `n_clusters`, `macro_site_season_rmse`, `median_site_season_rmse`, `macro_site_season_mae`, `macro_site_season_bias`, `macro_site_season_nse`, `n_finite_season_nse`, `macro_cluster_rmse`, `fit_reference`, `score_reference`, `arm`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 5. Table S1b — Principal LOYO scores and input contrasts for all five fixed families — input contrasts

`tables/TableS1b_all_LOYO_contrasts.csv` — 40 data rows; 15 columns.

Columns: `axis`, `model`, `variant`, `contrast`, `first_arm`, `second_arm`, `first_RMSE`, `second_RMSE`, `absolute_RMSE_difference`, `relative_percent`, `favourable_site_seasons`, `adverse_site_seasons`, `tie_site_seasons`, `n_site_seasons`, `n_regional_clusters`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 6. Table S2a — Ro1 calendar-proxy held-region absolute scores and input contrasts — scores/metrics

`tables/TableS2a_partial_Ro1_absolute_scores.csv` — 40 data rows; 8 columns.

Columns: `model`, `variant`, `arm`, `LOYO_equal_four_site_season_RMSE`, `partial_Ro1_LOSO_RMSE`, `same_Ro1_LOYO2015_RMSE`, `descriptive_relative_difference_percent`, `comparison_limit`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 7. Table S2b — Ro1 calendar-proxy held-region absolute scores and input contrasts — input contrasts

`tables/TableS2b_partial_Ro1_input_effects.csv` — 30 data rows; 13 columns.

Columns: `axis`, `model`, `fit_reference`, `score_reference`, `contrast`, `metric`, `first_arm`, `second_arm`, `first_RMSE`, `second_RMSE`, `absolute_RMSE_difference`, `relative_percent`, `primary`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 8. Table S3 — Own-reference processing-path by input contrasts

`tables/TableS3_reference_input_interactions.csv` — 20 data rows; 5 columns.

Columns: `model`, `contrast`, `VUT_own_contrast`, `CUT_own_contrast`, `CUT_minus_VUT_own_interaction_RMSE_units`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 9. Table S4 — Principal cohort selection and conservative field-date exclusions

`tables/TableS4_cohort_selection.csv` — 4 data rows; 8 columns.

Columns: `site_id`, `harvest_year`, `old_cases`, `new_selected`, `guard_excluded`, `after_documented_harvest`, `conservative_disputed_sowing_day`, `both_QC`.

Units: Counts/dates/flags as labelled; where present, contrasts retain the column-specific source definitions.

## 10. Table S5a — Broader point-input contextual metrics and regional scores — scores/metrics

`tables/TableS5a_broader_2431_metrics.csv` — 40 data rows; 15 columns.

Columns: `axis`, `model`, `target`, `seed`, `n_site_seasons`, `n_clusters`, `macro_site_season_rmse`, `median_site_season_rmse`, `macro_site_season_mae`, `macro_site_season_bias`, `macro_site_season_nse`, `n_finite_season_nse`, `macro_cluster_rmse`, `fit_reference`, `score_reference`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 11. Table S5b — Broader point-input contextual metrics and regional scores — region-specific scores

`tables/TableS5b_broader_2431_region_scores.csv` — 240 data rows; 15 columns.

Columns: `axis`, `model`, `target`, `seed`, `spatial_cluster`, `n`, `rmse`, `mae`, `bias_pred_minus_obs`, `nse`, `r2_reference_variance`, `pearson_r`, `squared_pearson_r`, `fit_reference`, `score_reference`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 12. Table S6 — Comparison of aggregation estimands

`tables/TableS6_aggregation_estimands.csv` — 40 data rows; 17 columns.

Columns: `model`, `variant`, `contrast`, `equal_site_season_first_RMSE`, `equal_site_season_second_RMSE`, `equal_site_season_absolute_difference`, `equal_site_season_relative_percent`, `equal_region_pooled_first_RMSE`, `equal_region_pooled_second_RMSE`, `equal_region_pooled_absolute_difference`, `equal_region_pooled_relative_percent`, `relative_effect_magnitude_shrink_percent`, `equal_region_mean_of_site_season_RMSE_first`, `equal_region_mean_of_site_season_RMSE_second`, `equal_region_mean_of_site_season_absolute_difference`, `equal_region_mean_of_site_season_relative_percent`, `definition`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 13. Table S7 — All principal LOYO site-season input contrasts

`tables/TableS7_site_season_contrasts.csv` — 160 data rows; 11 columns.

Columns: `model`, `variant`, `site_id`, `harvest_year`, `n`, `contrast`, `first_RMSE`, `second_RMSE`, `absolute_RMSE_difference`, `relative_percent`, `direction`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 14. Table S8 — Signs of principal site-season input contrasts

`tables/TableS8_sign_counts.csv` — 40 data rows; 11 columns.

Columns: `model`, `variant`, `contrast`, `favourable`, `adverse`, `tie`, `n`, `Ne1_2012_delta`, `Ne1_2013_delta`, `Ne3_2013_delta`, `Ro1_2015_delta`.

Units: Counts/dates/flags as labelled; where present, contrasts retain the column-specific source definitions.

## 15. Table S9 — Principal LOYO fold composition

`tables/TableS9_LOYO_fold_composition.csv` — 3 data rows; 12 columns.

Columns: `held_year`, `fold_id`, `n_train`, `n_test`, `training_floor`, `training_surplus_over_floor`, `out_of_test_region_training_cases`, `out_of_test_region_training_percent`, `future_year_training_cases`, `future_year_training_percent`, `training_site_seasons`, `chronological_validation`.

Units: Counts/dates/flags as labelled; where present, contrasts retain the column-specific source definitions.

## 16. Table S10 — Signed prediction-minus-reference bias shifts

`tables/TableS10_bias_shifts.csv` — 10 data rows; 10 columns.

Columns: `model`, `variant`, `A_shared_bias_pred_minus_obs`, `B_uniform_bias_pred_minus_obs`, `C_core_bias_pred_minus_obs`, `B_minus_A_shared_bias_shift`, `C_minus_B_bias_shift`, `C_minus_A_shared_bias_shift`, `bias_shift_recovery_percent`, `interpretation`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 17. Table S11 — Signed site-season contributions to principal RMSE contrasts

`tables/TableS11_Ne1_contribution_decomposition.csv` — 8 data rows; 9 columns.

Columns: `model`, `variant`, `contrast`, `sum_four_deltas`, `macro_delta`, `Ne1_two_season_sum_delta`, `Ne1_signed_contribution_percent`, `Ne3_signed_contribution_percent`, `definition`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 18. Table S12 — Asymmetry of field-date guard by region

`tables/TableS12_regional_guard_asymmetry.csv` — 2 data rows; 7 columns.

Columns: `region`, `old_cases`, `new_selected`, `guard_excluded`, `excluded_percent`, `operation_bound_evidence_available`, `qualification`.

Units: Counts/dates/flags as labelled; where present, contrasts retain the column-specific source definitions.

## 19. Table S13 — Fixed-output omission summaries

`tables/TableS13_fixed_OOF_summary_omissions.csv` — 60 data rows; 20 columns.

Columns: `model`, `reference`, `contrast`, `scenario`, `omitted_site_season`, `n_retained_site_seasons`, `retained_site_seasons`, `n_retained_evaluation_days`, `n_retained_regions`, `retained_regions`, `retained_region_site_season_counts`, `mean_first_RMSE`, `mean_second_RMSE`, `mean_RMSE_difference`, `relative_difference_percent`, `direction`, `same_sign_as_full_mean`, `summary_omission_only_no_retraining`, `omitted_unit_information_removed_from_other_folds_training`, `Ro1_field_calendar_gate_resolved`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 20. Table S14 — Sign stability across fixed-output omission summaries

`tables/TableS14_omission_sign_stability.csv` — 12 data rows; 20 columns.

Columns: `model`, `reference`, `contrast`, `full_four_season_RMSE_difference`, `full_four_season_relative_percent`, `full_mean_direction`, `n_favourable_omission_means`, `n_adverse_omission_means`, `n_tied_omission_means`, `same_sign_survives_all_four_omissions`, `benefit_survives_all_four_omissions`, `sign_reversing_omissions`, `min_omission_relative_percent`, `max_omission_relative_percent`, `max_absolute_change_from_full_RMSE_difference`, `dominant_absolute_unit`, `dominant_unit_signed_net_share_percent`, `maximum_absolute_signed_net_share_percent`, `maximum_absolute_contribution_share_percent`, `scope`.

Units: Counts/dates/flags as labelled; where present, contrasts retain the column-specific source definitions.

## 21. Table S15 — Site-season contributions to fixed-output contrasts

`tables/TableS15_fixed_unit_contributions.csv` — 48 data rows; 9 columns.

Columns: `model`, `reference`, `contrast`, `site_season`, `site_season_RMSE_difference`, `contribution_to_full_mean_RMSE_difference`, `signed_share_of_net_difference_percent`, `share_of_sum_absolute_unit_differences_percent`, `interpretation`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 22. Table S16 — Earlier-year Mead hindcast case counts

`tables/TableS16_T1_case_counts.csv` — 3 data rows; 3 columns.

Columns: `site_id`, `harvest_year`, `n`.

Units: Counts/dates/flags as labelled; where present, contrasts retain the column-specific source definitions.

## 23. Table S17a — Earlier-year hindcast macro and pooled scores — scores/metrics

`tables/TableS17a_T1_all_reference_macro.csv` — 40 data rows; 7 columns.

Columns: `arm`, `model`, `fit_reference`, `score_reference`, `macro_site_season_rmse`, `macro_site_season_mae`, `n_site_seasons`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 24. Table S17b — Earlier-year hindcast macro and pooled scores — pooled scores

`tables/TableS17b_T1_pooled.csv` — 40 data rows; 12 columns.

Columns: `arm`, `model`, `fit_reference`, `score_reference`, `n`, `rmse`, `mae`, `bias_pred_minus_obs`, `nse`, `r2_reference_variance`, `pearson_r`, `squared_pearson_r`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 25. Table S18 — Earlier-year hindcast site-season scores

`tables/TableS18_T1_site_season.csv` — 120 data rows; 14 columns.

Columns: `arm`, `model`, `fit_reference`, `score_reference`, `site_id`, `harvest_year`, `n`, `rmse`, `mae`, `bias_pred_minus_obs`, `nse`, `r2_reference_variance`, `pearson_r`, `squared_pearson_r`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 26. Table S19 — Earlier-year hindcast C-minus-B contrasts

`tables/TableS19_T1_C_minus_B.csv` — 20 data rows; 9 columns.

Columns: `model`, `fit_reference`, `score_reference`, `contrast`, `B_RMSE`, `C_RMSE`, `absolute_difference`, `relative_percent`, `primary`.

Units: GPP error and bias fields:g C m^-2 d^-1; percentages as named; NSE/correlation:dimensionless; counts as labelled.

## 27. Table S20 — Fixed-output squared-bias and centered-error metrics by unit

`tables/TableS20_error_components_by_unit.csv` — 760 data rows; 14 columns.

Columns: `cohort`, `axis`, `scope`, `unit`, `arm`, `model`, `fit_reference`, `score_reference`, `n`, `bias`, `bias_squared`, `centered_mse`, `mse`, `rmse`.

Units: MSE/bias-squared/centered-MSE fields:(g C m^-2 d^-1)^2; RMSE and normalized components:g C m^-2 d^-1; relative components:percentage points; counts as labelled.

## 28. Table S21 — Paired error-component contrasts by unit

`tables/TableS21_paired_component_contrasts.csv` — 380 data rows; 21 columns.

Columns: `cohort`, `axis`, `scope`, `unit`, `model`, `fit_reference`, `score_reference`, `first_arm`, `second_arm`, `contrast`, `n`, `first_rmse`, `second_rmse`, `bias_shift`, `delta_mse`, `delta_bias_squared`, `delta_centered_mse`, `delta_rmse`, `bias_rmse_component`, `centered_rmse_component`, `normalizing_denominator`.

Units: MSE/bias-squared/centered-MSE fields:(g C m^-2 d^-1)^2; RMSE and normalized components:g C m^-2 d^-1; relative components:percentage points; counts as labelled.

## 29. Table S22 — Fixed-output error-component summary estimands

`tables/TableS22_component_summary_estimands.csv` — 180 data rows; 25 columns.

Columns: `cohort`, `axis`, `scope`, `model`, `fit_reference`, `score_reference`, `contrast`, `first_arm`, `second_arm`, `n_units`, `n_days`, `mean_first_rmse`, `mean_second_rmse`, `mean_delta_rmse`, `mean_bias_rmse_component`, `mean_centered_rmse_component`, `mean_delta_mse`, `mean_delta_bias_squared`, `mean_delta_centered_mse`, `delta_rmse_percent`, `bias_rmse_percentage_points`, `centered_rmse_percentage_points`, `favorable_units`, `unfavorable_units`, `unchanged_units`.

Units: MSE/bias-squared/centered-MSE fields:(g C m^-2 d^-1)^2; RMSE and normalized components:g C m^-2 d^-1; relative components:percentage points; counts as labelled.

## 30. Table S23 — Archived monthly conditional-core weight concentration

`tables/TableS23_monthly_weight_concentration.csv` — 84 data rows; 13 columns.

Columns: `cohort`, `site_id`, `year`, `month`, `map_path`, `map_sha256`, `n_cells`, `source_conditional_weight_sum`, `concentration_cell_count`, `max_conditional_weight`, `distance_from_uniform`, `retained_dates_using_map`, `used_by_retained_case_anchor`.

Units: Conditional weights and uniform-distance measure: dimensionless; concentration and positive cell count: cell equivalents/counts.

## 31. Table S24 — Retained-date optical inputs and anchor lineage

`tables/TableS24_retained_daily_input_lineage.csv` — 1,079 data rows; 20 columns.

Columns: `cohort`, `site_id`, `year`, `date`, `T1_role`, `B_NDVI`, `C_NDVI`, `C_minus_B_NDVI`, `absolute_NDVI_difference`, `left_anchor`, `right_anchor`, `left_item_id`, `right_item_id`, `left_map_month`, `right_map_month`, `bracket_days`, `right_anchor_fraction`, `direct_anchor`, `right_anchor_after_case`, `cross_month_bracket`.

Units: NDVI and interpolation fractions: dimensionless; interpolation spans:days; support fields:counts/flags.

## 32. Table S25 — All seven retained site-season input summaries

`tables/TableS25_site_season_input_summary.csv` — 7 data rows; 20 columns.

Columns: `cohort`, `site_id`, `year`, `T1_role`, `n_days`, `mean_C_minus_B_NDVI`, `mean_absolute_NDVI_difference`, `median_absolute_NDVI_difference`, `max_absolute_NDVI_difference`, `n_exact_equal_BC`, `n_difference_above_tolerance`, `active_map_months`, `minimum_active_cell_count`, `maximum_active_cell_count`, `minimum_active_concentration`, `maximum_active_concentration`, `direct_anchor_days`, `right_anchor_after_case_days`, `cross_month_bracket_days`, `maximum_bracket_days`.

Units: NDVI and interpolation fractions: dimensionless; interpolation spans:days; support fields:counts/flags.

## 33. Table S26 — All source scene-anchor statuses and retained-date usage

`tables/TableS26_scene_anchor_usage.csv` — 277 data rows; 10 columns.

Columns: `cohort`, `site_id`, `year`, `anchor`, `item_id`, `map_month`, `shared_anchor_usable`, `retained_dates_using_anchor`, `B_anchor_NDVI`, `C_anchor_NDVI`.

Units: Counts/dates/flags as labelled; where present, contrasts retain the column-specific source definitions.

## 34. Table S27 — All four train/test RMSE cells under explicitly labelled aggregation estimands

`tables/TableS27_BC_cell_scores_estimands.csv` — 960 data rows; 13 columns.

Columns: `cohort`, `axis`, `model`, `fit_reference`, `score_reference`, `train_arm`, `test_arm`, `estimand`, `n_units`, `n_days`, `rmse`, `mae`, `bias_pred_minus_reference`.

Units: GPP RMSE/MAE/bias and signed contrasts: g C m^-2 d^-1. NSE/correlations: dimensionless. *_common_BB_pp fields:100 times the contrast divided by the matching R_BB, reported as common-baseline percentage points; counts/dates as named.

## 35. Table S28 — Conditional representation contrasts and interaction under the labelled estimands

`tables/TableS28_BC_conditional_contrasts_estimands.csv` — 240 data rows; 25 columns.

Columns: `cohort`, `axis`, `model`, `fit_reference`, `score_reference`, `estimand`, `R_BB`, `R_BC`, `R_CB`, `R_CC`, `training_given_test_B`, `training_given_test_C`, `test_given_training_B`, `test_given_training_C`, `interaction`, `original_pipeline_contrast`, `n_units`, `n_days`, `maximum_identity_residual`, `training_given_test_B_common_BB_pp`, `training_given_test_C_common_BB_pp`, `test_given_training_B_common_BB_pp`, `test_given_training_C_common_BB_pp`, `interaction_common_BB_pp`, `original_pipeline_contrast_common_BB_pp`.

Units: GPP RMSE/MAE/bias and signed contrasts: g C m^-2 d^-1. NSE/correlations: dimensionless. *_common_BB_pp fields:100 times the contrast divided by the matching R_BB, reported as common-baseline percentage points; counts/dates as named.

## 36. Table S29 — All train/test cell metrics at native evaluation scopes

`tables/TableS29_BC_cell_scores_native_units.csv` — 1,920 data rows; 20 columns.

Columns: `cohort`, `axis`, `model`, `fit_reference`, `score_reference`, `train_arm`, `test_arm`, `fold_id`, `scope`, `unit`, `n`, `rmse`, `mae`, `bias_pred_minus_reference`, `nse`, `pearson_r`, `squared_pearson_r`, `site_id`, `spatial_cluster`, `harvest_year`.

Units: GPP RMSE/MAE/bias and signed contrasts: g C m^-2 d^-1. NSE/correlations: dimensionless. *_common_BB_pp fields:100 times the contrast divided by the matching R_BB, reported as common-baseline percentage points; counts/dates as named.

## 37. Table S30 — All conditional contrasts and interactions at native evaluation scopes

`tables/TableS30_BC_conditional_contrasts_native_units.csv` — 480 data rows; 26 columns.

Columns: `cohort`, `axis`, `model`, `fit_reference`, `score_reference`, `scope`, `unit`, `R_BB`, `R_BC`, `R_CB`, `R_CC`, `training_given_test_B`, `training_given_test_C`, `test_given_training_B`, `test_given_training_C`, `interaction`, `original_pipeline_contrast`, `n_units`, `n_days`, `maximum_identity_residual`, `training_given_test_B_common_BB_pp`, `training_given_test_C_common_BB_pp`, `test_given_training_B_common_BB_pp`, `test_given_training_C_common_BB_pp`, `interaction_common_BB_pp`, `original_pipeline_contrast_common_BB_pp`.

Units: GPP RMSE/MAE/bias and signed contrasts: g C m^-2 d^-1. NSE/correlations: dimensionless. *_common_BB_pp fields:100 times the contrast divided by the matching R_BB, reported as common-baseline percentage points; counts/dates as named.


The archive preserves numeric source companions, not raw provider data, model objects, full prediction records or repository-deposit approval. Reproduction of model execution additionally requires the separately documented processing/model package and its source access conditions.
