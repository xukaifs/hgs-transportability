# Reproducibility record

Analysis is closed for v1.0.0. This repository adds no new statistical analysis. BCT remains the locked primary model; model-family comparisons are sensitivities. Version 1.0.0 is dated 2026-09-17.

## Verified final sensitivity

`code/sensitivity/hgs_summary_sensitivity.R` generates `example_results/sensitivity_hgs_definition.csv`. Four source fits use the locked BCT structure and the bilateral mean of per-hand maxima, restricted to both hands having valid measurements. Six scenario/sex rows cover frozen Korea-to-U.S., frozen U.S.-to-Korea and U.S. 2013-2014 held-out evaluation after 2011-2012 location-scale adaptation. All six have zero CDF clipping. The original HGSmax bundles are unchanged. The two definition/eligibility changes are assessed jointly, not separately.

## IPW formula and implementation

The final participation model is a survey-weighted quasibinomial logistic regression:

```r
measured ~ ns(age, 4) + sex + ns(height, 3) + ns(bmi, 3) + year_f
```

`year_f` is categorical survey year. Body weight is used only for completeness screening. Covariates are not inferred from prose: see the archived `run_IPW.R`. Stabilization uses the sex-year weighted participation proportion divided by the predicted probability. A separate sensitivity caps this multiplier at its weighted 99th percentile within sex-year. Source HGS reference models are not refitted for IPW. Missing covariates are excluded rather than imputed. IPW parameter uncertainty is not propagated.

## Race standardization target

Race/ethnicity-specific marginal mean Z scores are standardized to the MEC-weighted empirical age-height distribution of eligible participants of the same sex in the NHANES 2013-2014 held-out sample (2,527 women; 2,371 men under the primary definition). The diagnostic model is `z ~ race + ns(age,4) + ns(height,3)`. Counterfactual race-category model matrices are averaged over that same-sex target. Target-distribution and prior reference/adaptation parameter uncertainty are treated as fixed. This is subgroup calibration, not a causal effect.

## Specification, seeds, folds and uncertainty

- BCT: identity-link mu with age 4-df and height 3-df natural splines; log sigma with the same basis dimensions; constant nu and tau (tau has a log link). Source or fold-specific knots are frozen for scoring.
- Phase 2 seed: 20260910. Pooled family CV seed: 20260911. Two repeats of five PSU-grouped folds, shared across families and sexes; permutation and round-robin assignment with a random offset within year-strata. Exact assignments are retained in `reproducibility/`.
- The final bilateral sensitivity has no stochastic split or resampling.
- Existing CI tables are in `example_results/existing_CI/`; `reproducibility/existing_CI_inventory.csv` identifies columns. These are existing design-based intervals, not newly computed paired method-difference intervals. Some descriptive quantities have no CI. No missing interval has been filled with an invented value.
- Available historical sessionInfo files are preserved with stage prefixes. `sessionInfo_20260915.txt` describes only the final sensitivity environment.

See `docs/methods_verified.md` for exact QC, weights, threshold units and interpretation limits. `code/legacy/` is an audit archive of historical scripts with original local paths. This package has not been validated as a portable one-command pipeline; local harmonized input data and path configuration are required to refit models.

## Integrity and data checks

`MANIFEST.csv` contains byte size, MD5 and SHA256 for every payload file except the manifest itself. Git internals are excluded. `.gitattributes` preserves archived bytes across checkouts. The ZIP has a separate SHA256 sidecar. The ZIP is checked entry-by-entry against the final manifest before delivery.

Run `Rscript reproducibility/verify_archive.R` from the repository root to check model object structure, aggregated-output identifiers and sensitivity directions. Source scripts necessarily contain the variable names `SEQN` and `uid`; those code references are not participant records. Saved PSU fold assignments identify survey clusters, not individual participants. Raw survey files, individual score tables and individual harmonized RDS inputs are not distributed.

## Publication state

Author information is in `AUTHORS.md` and `CITATION.cff`. The public repository is https://github.com/xukaifs/hgs-transportability. Version 1.0.0 is dated 2026-09-17. MIT applies to author-generated repository resources, excluding third-party source data. A Zenodo DOI will be added after archival of the GitHub v1.0.0 release; no DOI is fabricated in this repository.
