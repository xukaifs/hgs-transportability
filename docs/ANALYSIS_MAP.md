# Analysis map

Version 1.0.0, dated 2026-09-17. No new statistical analyses were run during repository finalization. All paths below are relative to the repository root. Historical paths must be configured locally; scripts are provenance snapshots, not a one-command portable pipeline.

| Manuscript analysis | Main script(s) | Main archived output |
| --- | --- | --- |
| Korean development and temporal evaluation | `code/legacy/KNHANES_HGS_phase2_20260910/select_development.R`; `code/legacy/KNHANES_HGS_phase2_20260910/validate_locked.R` | `example_results/KNHANES_HGS_phase2_20260910/cv_repeat_metrics.csv`; `example_results/KNHANES_HGS_phase2_20260910/summary_validation.csv` |
| Pooled Korean BCT | `code/legacy/KNHANES_pooled_BCT_20260910/run_pooled_BCT.R` | `example_results/KNHANES_pooled_BCT_20260910/pooled_BCT_overall.csv` |
| Korea-to-U.S. frozen transport | `code/legacy/pooled_BCT_US_transport_20260910/run_transport.R`; `code/legacy/BCT_reciprocal_transport_20260911/run_reciprocal.R` | `example_results/BCT_reciprocal_transport_20260911/reciprocal_overall.csv` (Korea to US rows, pooled frozen comparison) |
| U.S. recalibration ladder | `code/legacy/US_BCT_recalibration_ladder_20260910/run_ladder.R` | `example_results/US_BCT_recalibration_ladder_20260910/US_ladder_heldout_overall.csv` |
| Simple-method comparison | `code/legacy/US_simple_methods_ladder_20260910/run_comparison.R` | `example_results/US_simple_methods_ladder_20260910/heldout_normalization_ladder_comparison.csv` |
| Distribution-family sensitivity | `code/legacy/HGS_family_CV_transport_20260911/run_Korea_CV.R`; `code/legacy/HGS_family_CV_transport_20260911/run_US_transport.R` | `example_results/HGS_family_CV_transport_20260911/Korea_OOF_nine_metrics.csv`; `example_results/HGS_family_CV_transport_20260911/US_heldout_nine_metrics.csv` |
| Age/race diagnostics | `code/legacy/HGS_age_race_transport_20260911/run_age_race.R`; `code/legacy/HGS_adjusted_race_age_profile_20260911/run_update.R` | `example_results/HGS_age_race_transport_20260911/age_transport_overall.csv`; `example_results/HGS_adjusted_race_age_profile_20260911/race_age_height_adjusted_means.csv` |
| Multidimensional performance profile | `code/legacy/HGS_adjusted_race_age_profile_20260911/run_update.R` | `example_results/HGS_adjusted_race_age_profile_20260911/performance_profile_with_age.csv` |
| IPW sensitivity | `code/legacy/HGS_family_CV_transport_20260911/run_IPW.R` | `example_results/HGS_family_CV_transport_20260911/IPW_yearly_curves.csv` |
| Reciprocal U.S.-to-Korea transport | `code/legacy/BCT_reciprocal_transport_20260911/run_reciprocal.R` | `example_results/BCT_reciprocal_transport_20260911/reciprocal_overall.csv` (US to Korea rows) |
| HGS-definition sensitivity | `code/sensitivity/hgs_summary_sensitivity.R` | `example_results/sensitivity_hgs_definition.csv` |

Compact frozen parameters are in `models/`. Existing CI tables retain original units and methods; their inventory is in `reproducibility/existing_CI_inventory.csv`. See `REPRODUCIBILITY.md` and `docs/methods_verified.md` for QC, weights, adjustment targets and interpretation limits. Historical study dates, folder names and sessionInfo timestamps are intentionally unchanged. BCT remains the locked primary model; this map does not imply new selection or a composite ranking.
