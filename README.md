# Handgrip-strength normalization and transportability

Reproducibility resources for the study of multidimensional calibration and transportability of handgrip-strength normalization across KNHANES and NHANES.

**Repository:** https://github.com/xukaifs/hgs-transportability  
**Release:** v1.0.0 (2026-09-17)  
**Zenodo DOI:** to be added after archival of the GitHub v1.0.0 release.

**Authors:** Xiaoqing Sun; Taoyang Han; Kai Xu (corresponding author). Full affiliations, contact information and ORCID iDs are in [AUTHORS.md](AUTHORS.md). Citation metadata are in [CITATION.cff](CITATION.cff).

**Analysis status: locked for v1.0.0.** Final bilateral sensitivity is complete. No further statistical expansion is planned for this release. See [REPRODUCIBILITY.md](REPRODUCIBILITY.md) for verified IPW/race methods, seeds, folds, CI provenance and integrity checks.

BCT remains the locked primary model. Gaussian and BCCG are fixed-structure model-family sensitivity analyses. Performance dimensions are reported separately; there is no composite score or post-hoc winner selection.

Alternative HGS summary sensitivity was completed and verified on 2026-09-15. The bilateral mean of per-hand maxima, requiring at least one valid trial per hand, preserved the directions of age and height slopes and overall tail deviations in both frozen transport directions and the U.S. held-out location-scale evaluation. This sensitivity jointly changes the summary definition and eligibility requirement; it does not isolate their separate effects.

## Contents

- `code/legacy/`: unchanged historical R scripts grouped by original analysis directory. These are provenance snapshots, not a single portable pipeline. Some contain historical local paths; configure paths before running. Do not run all scripts indiscriminately.
- `code/sensitivity/hgs_summary_sensitivity.R`: final bilateral-definition sensitivity analysis, with explicit input-root and output-directory arguments.
- `example_results/`: real aggregated study outputs, not synthetic participant data. `sensitivity_hgs_definition.csv` is the six-row final sensitivity table. Original and sensitivity estimates are paired in `example_results/HGS_definition_sensitivity_20260915/sensitivity_vs_primary.csv`.
- `models/`: compact frozen parameter bundles and adaptation constants. No participant-level score files or raw survey data are included.
- `reproducibility/`: manifests, CI inventory, historical locks, fold assignments, integrity checks and R session information.
- `docs/methods_verified.md`: code-verified IPW, race-standardization, model specifications and inferential limitations.
- `docs/ANALYSIS_MAP.md`: mapping from manuscript analyses to archived scripts and aggregate outputs.

## Reproduce the final sensitivity

Requires R, `gamlss`, `gamlss.dist`, `survey`, `splines`, `ggplot2` and the original local analysis input tree. Package versions are recorded in `reproducibility/sessionInfo_20260915.txt`. The historical core currently adds a Windows library path; configure it for another environment.

```text
Rscript code/sensitivity/hgs_summary_sensitivity.R "PATH/TO/Paper 4" "PATH/TO/NEW_OUTPUT"
```

Required private/local inputs beneath that root are the Korean development/validation RDS files, U.S. harmonized QC RDS, historical `KNHANES_HGS_phase2_20260910/core.R`, original pooled source bundles, reciprocal overall CSV and family-comparison U.S. held-out diagnostic CSV. The script checks original maximum HGS against trial-level reconstruction before changing the definition. New source models are fitted with the locked structure; target observations never estimate source-model parameters. Only 2011-2012 U.S. observations estimate U.S. location-scale adaptation.

Raw KNHANES/NHANES data must be obtained separately under the relevant provider conditions. This repository is not a self-contained raw-data redistribution. Exact historical scripts are retained for audit; full portable end-to-end reproduction has not been tested.

## Data and privacy

Participant-level KNHANES and NHANES data are not redistributed. The public repository contains aggregate outputs, compact model/parameter objects, code, fold assignments at the PSU level, and reproducibility metadata. Raw survey files, individual score tables and harmonized participant-level RDS inputs are excluded.

## Citation and archival

Use the metadata in [CITATION.cff](CITATION.cff) when citing the software/reproducibility archive. The exact v1.0.0 release will be archived with Zenodo; once the DOI is minted, the DOI will be added to the repository metadata and manuscript without changing the archived analysis results.

## License

The MIT License applies to repository code and author-generated model/analysis resources. NHANES and KNHANES source data are not redistributed and remain subject to their respective provider terms.
