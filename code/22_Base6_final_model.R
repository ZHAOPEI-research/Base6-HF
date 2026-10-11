# ============================================================
# STEP 8A
# Final simplified 6-variable HF prognostic model
#
# Predictors:
# Age
# BMI
# RDW
# NWR (natural spline, 3 df)
# eGFR
# Medication count
# ============================================================

rm(list = ls())

library(mice)
library(survival)
library(splines)
library(dplyr)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)


# ============================================================
# 1. Load final data / frozen spline / bootstrap result
# ============================================================

load(
  "NHANES_2011_2018_MICE_FINAL.RData"
)

load(
  "NHANES_2011_2018_STEP6A_GROUP_LASSO_FINAL.RData"
)

load(
  "NHANES_2011_2018_STEP7D_BOOTSTRAP_4567.RData"
)


cat("============================================\n")
cat("STEP 8A: FINAL Base6 model\n")
cat("============================================\n")


# ============================================================
# 2. Bootstrap-derived shrinkage
# ============================================================

shrinkage <-
  final_summary %>%
  dplyr::filter(
    model == "Base6_Med_BMI"
  ) %>%
  dplyr::pull(
    corrected_slope
  )


cat("\n============================================\n")
cat("A. Shrinkage factor\n")
cat("============================================\n")

cat(
  "Base6 optimism-corrected calibration slope =",
  shrinkage,
  "\n"
)


# ============================================================
# 3. Data preparation
# ============================================================

prepare_data <- function(d) {

  d <- d %>%
    mutate(

      # NHANES age top-coded at 80
      age80 =
        pmin(
          ridageyr,
          80
        ),

      # Final NWR
      nwr =
        nwr_model,

      # CKD-EPI 2021
      kappa =
        ifelse(
          female == 1,
          0.7,
          0.9
        ),

      alpha =
        ifelse(
          female == 1,
          -0.241,
          -0.302
        ),

      egfr_2021 =
        142 *

        pmin(
          scr_mgdl / kappa,
          1
        )^alpha *

        pmax(
          scr_mgdl / kappa,
          1
        )^(-1.200) *

        (0.9938^ridageyr) *

        ifelse(
          female == 1,
          1.012,
          1
        )
    )


  # ----------------------------------------------------------
  # Frozen NWR spline
  # ----------------------------------------------------------

  nb <- splines::ns(

    d$nwr,

    knots =
      nwr_knots,

    Boundary.knots =
      nwr_boundary,

    intercept =
      FALSE
  )


  d$nwr_ns1 <- nb[, 1]
  d$nwr_ns2 <- nb[, 2]
  d$nwr_ns3 <- nb[, 3]


  d
}


# ============================================================
# 4. Final Base6 formula
# ============================================================

base6_formula <-

  Surv(
    followup_month,
    death
  ) ~

  age80 +
  bmxbmi +
  lbxrdw +

  nwr_ns1 +
  nwr_ns2 +
  nwr_ns3 +

  egfr_2021 +
  med_cnt


# ============================================================
# 5. Fit in all 30 MI datasets
# ============================================================

M <- imp_final$m

fits <- vector(
  "list",
  M
)

coef_list <- vector(
  "list",
  M
)

vcov_list <- vector(
  "list",
  M
)

completed_list <- vector(
  "list",
  M
)


for (i in 1:M) {

  d <-
    prepare_data(
      complete(
        imp_final,
        i
      )
    )


  completed_list[[i]] <-
    d


  fit <-
    coxph(

      base6_formula,

      data =
        d,

      ties =
        "efron",

      x =
        TRUE,

      y =
        TRUE
    )


  fits[[i]] <-
    fit

  coef_list[[i]] <-
    coef(fit)

  vcov_list[[i]] <-
    vcov(fit)
}


# ============================================================
# 6. Check coefficient consistency
# ============================================================

coef_names <-
  names(
    coef_list[[1]]
  )


same_names <-
  all(
    sapply(
      coef_list,
      function(x) {

        identical(
          names(x),
          coef_names
        )
      }
    )
  )


if (!same_names) {

  stop(
    "Coefficient names differ across imputations."
  )
}


cat("\n============================================\n")
cat("B. Final model coefficient names\n")
cat("============================================\n")

print(
  coef_names
)


# ============================================================
# 7. Rubin pooling
# ============================================================

Q <-
  do.call(
    rbind,
    coef_list
  )


beta_pool <-
  colMeans(
    Q
  )


# Within-imputation variance
U_bar <-
  Reduce(
    "+",
    vcov_list
  ) /
  M


# Between-imputation variance
B_cov <-
  stats::cov(
    Q
  )


# Total Rubin variance
T_cov <-
  U_bar +
  (1 + 1/M) *
  B_cov


SE_pool <-
  sqrt(
    diag(
      T_cov
    )
  )


# ============================================================
# 8. Rubin-pooled coefficients
# ============================================================

pooled_table <- data.frame(

  term =
    coef_names,

  beta =
    as.numeric(
      beta_pool
    ),

  SE =
    as.numeric(
      SE_pool
    ),

  HR =
    exp(
      beta_pool
    ),

  CI_low =
    exp(
      beta_pool -
        1.96 *
        SE_pool
    ),

  CI_high =
    exp(
      beta_pool +
        1.96 *
        SE_pool
    ),

  stringsAsFactors =
    FALSE
)


cat("\n============================================\n")
cat("C. Rubin-pooled Base6 coefficients\n")
cat("============================================\n")

print(
  pooled_table,
  row.names = FALSE
)


# ============================================================
# 9. Coefficient variability across MI
# ============================================================

coef_stability <- data.frame(

  term =
    coef_names,

  pooled_beta =
    beta_pool,

  min_beta =
    apply(
      Q,
      2,
      min
    ),

  max_beta =
    apply(
      Q,
      2,
      max
    ),

  SD_between_MI =
    apply(
      Q,
      2,
      sd
    )
)


cat("\n============================================\n")
cat("D. Coefficient stability across 30 MI\n")
cat("============================================\n")

print(
  coef_stability,
  row.names = FALSE
)


# ============================================================
# 10. Uniform shrinkage
# ============================================================

beta_shrunk <-
  beta_pool *
  shrinkage


shrunk_table <- data.frame(

  term =
    coef_names,

  beta_original =
    as.numeric(
      beta_pool
    ),

  shrinkage =
    shrinkage,

  beta_shrunk =
    as.numeric(
      beta_shrunk
    ),

  HR_shrunk =
    exp(
      beta_shrunk
    ),

  stringsAsFactors =
    FALSE
)


cat("\n============================================\n")
cat("E. Uniformly shrunk coefficients\n")
cat("============================================\n")

print(
  shrunk_table,
  row.names = FALSE
)


# ============================================================
# 11. Freeze model definition
# ============================================================

model_definition <- list(

  model_name =
    "Simplified HF 6-variable model",

  predictors =
    c(
      "Age",
      "BMI",
      "RDW",
      "NWR",
      "eGFR",
      "Medication count"
    ),

  n_predictors =
    6,

  n_coefficients =
    length(
      beta_shrunk
    ),

  age_definition =
    "Age in years, capped at 80 years",

  bmi_definition =
    "BMI, kg/m2",

  rdw_definition =
    "Red cell distribution width",

  nwr_definition =
    "Absolute neutrophil count / total WBC",

  nwr_form =
    "Natural cubic spline with 3 df",

  nwr_knots =
    nwr_knots,

  nwr_boundary =
    nwr_boundary,

  egfr_definition =
    "CKD-EPI 2021 creatinine equation",

  medication_definition =
    "Number of prescription medications",

  shrinkage_factor =
    shrinkage
)


# ============================================================
# 12. Save
# ============================================================

write.csv(

  pooled_table,

  "Base6_Rubin_pooled_coefficients.csv",

  row.names =
    FALSE
)


write.csv(

  coef_stability,

  "Base6_coefficient_stability_MI.csv",

  row.names =
    FALSE
)


write.csv(

  shrunk_table,

  "Base6_shrunk_coefficients.csv",

  row.names =
    FALSE
)


save(

  beta_pool,

  beta_shrunk,

  T_cov,

  pooled_table,

  shrunk_table,

  coef_stability,

  model_definition,

  base6_formula,

  completed_list,

  fits,

  nwr_knots,

  nwr_boundary,

  shrinkage,

  file =
    "NHANES_Base6_FINAL_MODEL.RData"
)


cat("\n============================================\n")
cat("STEP 8A 完成\n")
cat("============================================\n")
