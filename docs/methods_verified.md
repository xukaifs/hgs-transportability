# Methods verified directly against analysis code

## HGS definition sensitivity

The alternative outcome is (maximum valid hand 1 + maximum valid hand 2)/2. Both hands must have at least one valid trial; original participant eligibility remains required. Korea retains finite trials between 0 and 100 kg inclusive. U.S. retains nonnegative trials with the original trial effort flag equal to 1, with original examination-status eligibility. A zero trial is not newly excluded. The bilateral mean must be positive. For U.S., the symmetric mean does not require identifying which tested hand was left/right.

Sex-specific Korean 2014-2019 and U.S. 2011-2014 BCT models are refitted for this different outcome and frozen before transport. Primary HGSmax bundles remain unchanged. Four fits only; no CV/model selection. U.S. adaptation uses 2011-2012 weighted mean and population SD of Korean-reference Z; 2013-2014 supplies only evaluation. The pooled factor 1/2 on U.S. MEC weights cancels within a single-cycle weighted mean/SD. Korean weights are annual interview/examination weights divided by six. Results are survey-weighted descriptive points, not new confidence intervals.

## IPW

Source: `code/legacy/HGS_family_CV_transport_20260911/run_IPW.R`.

Actual response model: `svyglm(measured ~ ns(age,4) + sex + ns(height,3) + ns(bmi,3) + factor(year), family=quasibinomial())` using pooled survey PSU, strata and weights. There are no additional socioeconomic, disease, or interaction covariates.

Base: age 20-79, recognized sex, valid design and positive finite height, BMI and body weight. Body weight is an eligibility requirement, **not an additional predictor**. Measured means finite positive HGS with at least one valid trial. Missing predictor/design records are excluded, not imputed.

Stabilized numerator: weighted response proportion in each sex-year. Denominator: pooled-model predicted response probability. Original survey weights are multiplied by this ratio. The capped analysis caps the stabilized ratio at its survey-weighted 99th percentile among measured respondents within sex-year; it does not cap the final survey weights. Uncapped and capped results are both retained. Reference Z is unchanged. This analysis assumes conditional response exchangeability and cannot rule out unobserved selection or protocol changes; response-model uncertainty was not propagated.

## Race/ethnicity standardized marginal means

Source: `code/legacy/HGS_adjusted_race_age_profile_20260911/` R script.

Fit separately for each sex, existing normalization method and calibration: `svyglm(z ~ race + ns(age,4) + ns(height,3))`.

Target is the **same-sex eligible 2013-2014 U.S. held-out analysis sample**, not all ages, not both sexes combined, and not the Korean sample. The complete target comprises 2,527 women and 2,371 men for the primary definition. Each target person's age and height are retained, race is set to each category in turn, and the resulting model matrix is averaged using normalized MEC weights. The marginal mean is the resulting linear contrast of coefficients. Its CI uses the contrast covariance and model residual-df t critical value. The empirical target distribution and previous scoring/calibration parameters are treated as fixed. This is descriptive subgroup calibration, not a causal race effect.

## Frozen structures, thresholds and CV

Primary BCT: identity-link mu = ns(age,4)+ns(height,3); log(sigma) uses the same covariates; nu and log(tau) are intercept-only. Gaussian/NO and BCCG family sensitivities use the same mu/sigma basis dimensions. Knots and boundary knots are learned from each source training set; CV folds never reuse whole-sample basis parameters. Sigma is not literal kg SD for BCT/BCCG.

Phase 2 seed: 20260910. Pooled family sensitivity seed: 20260911. Each uses two repeats of five PSU-grouped folds. Distinct PSUs are permuted within pooled year-strata; a random offset 0-4 and round-robin modulo-five allocation balance fold counts. Assignment is common across sexes and candidate families; a PSU stays wholly in one fold. Exact saved assignments are included. Preserve original row ordering when reconstructing from seeds. Bilateral sensitivity uses no random split or stochastic resampling.

P10 is Z < -1.282; P5 is Z < -1.645. CSV sensitivity p10/p5 are percentages. SD is the survey-weighted population SD, not SE. Age/height slopes are separate marginal weighted linear regressions per 10 years/10 cm, not mutually adjusted effects. Original supplementary CSVs can use proportions; inspect column documentation and originating code before combining.

## Existing CIs and environment

CI-bearing original aggregated tables are copied to `example_results/existing_CI/`; their inventory lists file and columns. Original scale/units are preserved. Survey design accounts for PSU/strata in the original inferential analyses; frozen scoring/adaptation parameter uncertainty is generally excluded. Not every descriptive statistic has a CI. No paired bootstrap or new method-difference CI has been fabricated.

The 2026-09-15 sessionInfo records the environment used for the final sensitivity only. It must not be described as a recovered historical sessionInfo for earlier analyses. MD5 values verify archived file integrity; no historical timestamps or release provenance are invented.
