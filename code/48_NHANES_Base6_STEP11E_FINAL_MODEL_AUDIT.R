# ============================================================
# STEP 11E
# FINAL Base6 coefficient provenance / model-freeze audit
#
# Purpose:
#   1. Freeze the exact pre-shrinkage coefficients
#   2. Verify optimism shrinkage
#   3. Generate final shrunken coefficients
#   4. Freeze NWR spline definition
#   5. Freeze LP centering constant
#   6. Freeze 1-year baseline cumulative hazard
#   7. Write the exact prediction formula
#
# IMPORTANT:
#   This step does NOT refit the model.
# ============================================================

rm(list = ls())
gc()

library(dplyr)
library(tibble)


analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
if (!dir.exists(analysis_dir)) {
  dir.create(analysis_dir, recursive = TRUE)
}

cat("============================================\n")
cat("STEP 11E: FINAL Base6 model audit\n")
cat("============================================\n\n")


# ============================================================
# 1. Frozen model specification
# ============================================================

model_name <- "NHANES Base6"

predictors <- c(
  "age80",
  "bmxbmi",
  "lbxrdw",
  "nwr_ns1",
  "nwr_ns2",
  "nwr_ns3",
  "egfr_2021",
  "med_cnt"
)

cat("Model:", model_name, "\n")
cat(
  "Coefficient terms:",
  paste(predictors, collapse = ", "),
  "\n\n"
)


# ============================================================
# 2. Rubin-pooled coefficients BEFORE shrinkage
#
# These are the coefficients obtained after:
#   - 30 multiple imputations
#   - fixed Base6 specification
#   - Rubin pooling
#
# They are NOT yet the final prediction coefficients.
# ============================================================

beta_pre_shrink <- c(
  age80     =  0.038175446,
  bmxbmi    = -0.020972564,
  lbxrdw    =  0.203942007,
  nwr_ns1   =  0.729017806,
  nwr_ns2   =  0.177541263,
  nwr_ns3   =  0.513577990,
  egfr_2021 = -0.009750497,
  med_cnt   =  0.057526831
)


# ============================================================
# 3. Bootstrap-derived global shrinkage factor
#
# Derived from fixed-model internal bootstrap validation.
# ============================================================

shrinkage_factor <- 0.9212384


# ============================================================
# 4. Calculate FINAL shrunken coefficients
# ============================================================

beta_final <- beta_pre_shrink * shrinkage_factor


cat("============================================\n")
cat("A. Shrinkage audit\n")
cat("============================================\n")

cat(
  "Global shrinkage factor =",
  sprintf("%.7f", shrinkage_factor),
  "\n\n"
)

coef_table <- tibble(
  term = names(beta_pre_shrink),

  beta_pre_shrink =
    as.numeric(beta_pre_shrink),

  shrinkage =
    shrinkage_factor,

  beta_final =
    as.numeric(beta_final),

  HR_pre_shrink =
    exp(beta_pre_shrink),

  HR_final =
    exp(beta_final),

  beta_change =
    beta_final -
    beta_pre_shrink,

  percent_shrink =
    100 *
    (
      1 -
        abs(beta_final) /
        abs(beta_pre_shrink)
    )
)

print(
  as.data.frame(coef_table),
  row.names = FALSE,
  digits = 9
)


# ============================================================
# 5. Mathematical consistency check
# ============================================================

reconstructed_beta <- beta_pre_shrink * shrinkage_factor

max_error <- max(
  abs(
    beta_final -
      reconstructed_beta
  )
)

cat("\n============================================\n")
cat("B. Mathematical consistency\n")
cat("============================================\n")

cat(
  "Maximum reconstruction error =",
  format(
    max_error,
    scientific = TRUE
  ),
  "\n"
)

stopifnot(
  max_error < 1e-12
)

cat(
  "PASS: beta_final = beta_pre_shrink × shrinkage_factor\n"
)


# ============================================================
# 6. Freeze FINAL NWR spline
# ============================================================

nwr_internal_knots <- c(
  0.570093457943925,
  0.652173913043478
)

nwr_boundary_knots <- c(
  0.434782608695652,
  0.773607367475292
)

cat("\n============================================\n")
cat("C. Frozen NWR spline\n")
cat("============================================\n")

cat("Natural spline df = 3\n")

cat(
  "Internal knots =",
  paste(
    sprintf(
      "%.15f",
      nwr_internal_knots
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Boundary knots =",
  paste(
    sprintf(
      "%.15f",
      nwr_boundary_knots
    ),
    collapse = ", "
  ),
  "\n"
)


# ============================================================
# 7. Freeze LP centering and baseline hazard
# ============================================================

LP_center <- 4.415976

# Corrected centered development baseline
# at approximately 1 year / 12 months
H0_1y <- 0.05383994

S0_1y <- exp(
  -H0_1y
)

cat("\n============================================\n")
cat("D. Baseline-risk constants\n")
cat("============================================\n")

cat(
  "LP_center =",
  sprintf("%.6f", LP_center),
  "\n"
)

cat(
  "H0_1y =",
  sprintf("%.8f", H0_1y),
  "\n"
)

cat(
  "S0_1y =",
  sprintf("%.8f", S0_1y),
  "\n"
)

cat(
  "Centered reference 1-y risk =",
  sprintf(
    "%.8f",
    1 - S0_1y
  ),
  "\n"
)


# ============================================================
# 8. Print exact FINAL coefficients
# ============================================================

cat("\n============================================\n")
cat("E. FINAL shrunken coefficients\n")
cat("============================================\n")

for (nm in names(beta_final)) {

  cat(
    sprintf(
      "%-12s = %.9f\n",
      nm,
      beta_final[[nm]]
    )
  )
}


# ============================================================
# 9. Exact FINAL linear predictor
# ============================================================

cat("\n============================================\n")
cat("F. FINAL linear predictor\n")
cat("============================================\n\n")

cat(
  "LP =\n",
  sprintf(
    "  %.9f * age80\n",
    beta_final["age80"]
  ),
  sprintf(
    "  %+.9f * bmxbmi\n",
    beta_final["bmxbmi"]
  ),
  sprintf(
    "  %+.9f * lbxrdw\n",
    beta_final["lbxrdw"]
  ),
  sprintf(
    "  %+.9f * nwr_ns1\n",
    beta_final["nwr_ns1"]
  ),
  sprintf(
    "  %+.9f * nwr_ns2\n",
    beta_final["nwr_ns2"]
  ),
  sprintf(
    "  %+.9f * nwr_ns3\n",
    beta_final["nwr_ns3"]
  ),
  sprintf(
    "  %+.9f * egfr_2021\n",
    beta_final["egfr_2021"]
  ),
  sprintf(
    "  %+.9f * med_cnt\n",
    beta_final["med_cnt"]
  ),
  sep = ""
)


# ============================================================
# 10. FINAL prediction function
# ============================================================

predict_base6_1y <- function(
    age80,
    bmxbmi,
    lbxrdw,
    nwr_ns1,
    nwr_ns2,
    nwr_ns3,
    egfr_2021,
    med_cnt
) {

  lp <-
    beta_final["age80"] *
    age80 +

    beta_final["bmxbmi"] *
    bmxbmi +

    beta_final["lbxrdw"] *
    lbxrdw +

    beta_final["nwr_ns1"] *
    nwr_ns1 +

    beta_final["nwr_ns2"] *
    nwr_ns2 +

    beta_final["nwr_ns3"] *
    nwr_ns3 +

    beta_final["egfr_2021"] *
    egfr_2021 +

    beta_final["med_cnt"] *
    med_cnt

  lp_centered <-
    lp -
    LP_center

  risk_1y <-
    1 -
    exp(
      -H0_1y *
        exp(
          lp_centered
        )
    )

  return(
    risk_1y
  )
}


# ============================================================
# 11. Formula statement
# ============================================================

cat("\n============================================\n")
cat("G. FINAL prediction equation\n")
cat("============================================\n\n")

cat(
  "LP_centered = LP - ",
  LP_center,
  "\n\n",
  sep = ""
)

cat(
  "Risk_1y = 1 - exp[-",
  H0_1y,
  " * exp(LP_centered)]\n",
  sep = ""
)

cat(
  "\nEquivalent:\n"
)

cat(
  "Risk_1y = 1 - exp[-",
  H0_1y,
  " * exp(LP - ",
  LP_center,
  ")]\n",
  sep = ""
)


# ============================================================
# 12. Audit of representative LP values
# ============================================================

lp_values <- c(
  low = LP_center - 2,
  center = LP_center,
  high = LP_center + 2
)

risk_values <- 1 -
  exp(
    -H0_1y *
      exp(
        lp_values -
          LP_center
      )
  )

risk_audit <- tibble(
  scenario =
    names(lp_values),

  LP =
    as.numeric(lp_values),

  LP_centered =
    as.numeric(
      lp_values -
        LP_center
    ),

  predicted_1y_risk =
    as.numeric(risk_values)
)

cat("\n============================================\n")
cat("H. Prediction-formula sanity check\n")
cat("============================================\n")

print(
  as.data.frame(risk_audit),
  row.names = FALSE,
  digits = 8
)


# ============================================================
# 13. Model provenance table
# ============================================================

provenance <- tibble(
  component = c(
    "Predictor specification",
    "NWR functional form",
    "Pre-shrink coefficients",
    "Shrinkage factor",
    "Final coefficients",
    "LP centering",
    "1-y baseline cumulative hazard",
    "1-y survival function"
  ),

  source = c(
    "Final Base6 selection",
    "Natural spline df=3; q05/q95 boundaries",
    "30-MI Rubin-pooled fixed Cox model",
    "Fixed-model bootstrap calibration slope",
    "Pre-shrink beta × global shrinkage",
    "Development LP distribution",
    "Corrected centered development baseline",
    "1-exp[-H0*exp(LP-LP_center)]"
  )
)

cat("\n============================================\n")
cat("I. Model provenance\n")
cat("============================================\n")

print(
  as.data.frame(provenance),
  row.names = FALSE
)


# ============================================================
# 14. Save frozen model constants
# ============================================================

save(
  model_name,
  beta_pre_shrink,
  shrinkage_factor,
  beta_final,
  coef_table,
  nwr_internal_knots,
  nwr_boundary_knots,
  LP_center,
  H0_1y,
  S0_1y,
  predict_base6_1y,
  provenance,
  file =
    file.path(
      analysis_dir,
      "NHANES_Base6_STEP11E_FINAL_MODEL_FREEZE.RData"
    )
)

write.csv(
  coef_table,
  file.path(
    analysis_dir,
    "NHANES_Base6_STEP11E_FINAL_COEFFICIENTS.csv"
  ),
  row.names = FALSE
)

write.csv(
  provenance,
  file.path(
    analysis_dir,
    "NHANES_Base6_STEP11E_MODEL_PROVENANCE.csv"
  ),
  row.names = FALSE
)


cat("\n============================================\n")
cat("STEP 11E completed\n")
cat("Final Base6 model constants frozen.\n")
cat("============================================\n")
