# Base6 HF Public Analysis — Run Order

This document records the tested execution order for the public analysis scripts used in the Base6 heart-failure prognostic model project.

## 1. Path configuration

Before running the external-validation scripts, set the following environment variables in the R session:

```r
Sys.setenv(
  BASE6_ANALYSIS_DIR = "~/Documents/R语言/NHANES",
  BASE6_EXTERNAL_DATA = "~/Documents/R语言/NHANES/hospital_data.xlsx",
  BASE6_EXTERNAL_OUTPUT_DIR = "~/Documents/R语言/NHANES/external_private_outputs"
)
```

`BASE6_EXTERNAL_DATA` and `BASE6_EXTERNAL_OUTPUT_DIR` point to restricted local resources and must not be committed to GitHub or deposited in Zenodo.

## 2. Standalone model implementation

From the repository root:

```r
source("implementation/Base6_predict_1year.R")
source("implementation/example_usage.R")
```

The example uses fictional input values only.

## 3. External cohort cleaning and Base6 prediction

Run in this order:

```r
source("code/33_Chinese_external_data_cleaning.R")
source("code/34_Chinese_NWR_audit.R")
source("code/35_Chinese_NWR_final_cleaning.R")
source("code/36_Chinese_Base6_external_validation.R")
```

Notes:

- Script 33 is an initial cleaning/audit step and may emit an expected numeric-conversion warning for a raw neutrophil value.
- Scripts 34–35 audit and finalize the NWR derivation.
- Script 36 is an intermediate external-validation script. It predates the final 364-day discrimination/calibration workflow and should not be used as the source of the manuscript-final external performance estimates.

## 4. MAGGIC comparison

```r
source("code/40_MAGGIC_mapping_audit.R")
source("code/41_MAGGIC_score_head_to_head.R")
```

These scripts create the MAGGIC mapping and comparison objects used by later head-to-head analyses.

## 5. GWTG-HF preparation and three-model comparison

Run the verified GWTG-HF cleaning/correction step before the score comparison:

```r
source("code/44_GWTG_verified_corrections.R")
source("code/45_GWTG_score_comparison.R")
```

The verified corrections are based on source-record review and are independent of outcome status or model predictions.

## 6. Paired discrimination comparison

Intermediate 2000-bootstrap analysis:

```r
source("code/46_GWTG_paired_bootstrap.R")
```

Final manuscript robustness analysis with 5000 paired bootstrap resamples:

```r
source("code/46_STEP11N_5000Bootstrap.R")
```

The 5000-bootstrap script is the preferred source for manuscript-final paired AUC and Harrell C-index differences.

## 7. Calibration and post-hoc recalibration

Run:

```r
source("code/46_STEP_11D_B_Base6vsMAGGIC.R")
source("code/46_STEP_11E_calibration.R")
source("code/46_STEP_11F_pohs_hoc_recalbration.R")
source("code/46_STEP_11I_Recalivration_bootstrap.R")
source("code/54_External_Calibration_CI.R")
```

Important:

- `46_STEP_11D_B_Base6vsMAGGIC.R` generates `Chinese_GWTG_STEP11D_CALIBRATION.RData`, which is required by later calibration scripts.
- The post-hoc recalibration results are apparent/updating analyses performed in the external validation cohort and should not be interpreted as independent external validation.
- The public version of `46_STEP_11I_Recalivration_bootstrap.R` must include the naming fix in both helper functions:

```r
apparent <- as.numeric(apparent)
```

This prevents named scalars from producing `NA` values in the final optimism-correction summary.

## 8. Decision curve analysis

Run:

```r
source("code/46_STEP_11G_DCA.R")
source("code/46_SHTP_11H_DCA_bootstrap.R")
source("code/46_STEP11O_DCA_thershold_range_sensitivity.R")
```

The recalibrated Base6 curves are post-hoc sensitivity analyses. The original Base6 curve is the primary external-validation DCA result.

## 9. Horizon sensitivity

```r
source("code/46_STEP11K_360_364_sensitivity.R")
```

This analysis confirms that 360- and 364-day discrimination estimates are similar and that 365-day time-dependent AUC/Brier estimation is not supported because no participant has follow-up beyond 365 days. The primary evaluable 1-year discrimination horizon is therefore 364 days.

## 10. Modernized MAGGIC sensitivity

```r
source("code/46_STEP11L_B6_vs_original_vs_moderinized_MAGGIC.R")
```

This analysis evaluates a modernized RAAS-treatment definition as a sensitivity analysis while retaining the original MAGGIC implementation as the primary benchmark.

## 11. Common-versus-maximum cohort sensitivity

```r
source("code/46_STEP11M_COMMON_MAXIMUM_COHORT.R")
```

This analysis assesses whether restricting all three models to the common n=346 cohort materially changes discrimination estimates.

## 12. Limited subgroup analysis

```r
source("code/46_STEP11P_Subgroup_limited_subgroup.R")
```

Formal subgroup performance is calculated only for subgroup levels with at least 15 deaths.

## 13. Tested status

The following external-validation modules were smoke-tested successfully in the public-code workflow:

- External cohort cleaning and NWR derivation
- Base6 external prediction
- MAGGIC mapping and scoring
- GWTG-HF verified scoring
- Three-model discrimination comparison
- 2000-bootstrap paired discrimination
- 5000-bootstrap paired discrimination
- Calibration
- Post-hoc recalibration
- Recalibration optimism correction
- External calibration confidence intervals
- Decision curve analysis
- Paired-bootstrap DCA
- DCA threshold-range sensitivity
- 360/364-day horizon sensitivity
- Modernized MAGGIC sensitivity
- Common-versus-maximum cohort sensitivity
- Limited subgroup analysis

## 14. Data and privacy

The external validation dataset contains individual-level clinical data and is not distributed in this repository.

Do not commit or archive:

- the source external cohort spreadsheet;
- patient-level `.csv`, `.xlsx`, `.RData`, `.rds`, or similar files;
- intermediate external-validation objects containing patient identifiers;
- private output directories.

The public repository should contain analysis code, model implementation code, metadata, documentation, and non-identifiable example material only.

## 15. Reproducibility scope

The Base6 prediction function and its fictional example can be run directly from the public repository.

Development-cohort analyses can be reproduced using the required NHANES public-use data and the documented analysis workflow.

External-validation scripts are provided for transparency and methodological reproducibility, but they require access to the restricted external cohort and therefore cannot be executed end-to-end from the public repository alone.
