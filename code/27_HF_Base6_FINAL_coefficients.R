# ============================================================
# STEP 8F
# FINAL Base6 pooled coefficients
#
# Final predictors:
# Age + BMI + RDW + NWR spline + eGFR + Medication count
#
# FINAL NWR spline:
# internal knots = 0.5700935, 0.6521739
# boundary knots = 0.4347826, 0.7736074
#
# 使用 STEP 8E 最新bootstrap calibration slope
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
# 1. Load final MI and final bootstrap
# ============================================================

load(
  "NHANES_2011_2018_MICE_FINAL.RData"
)

load(
  "NHANES_Base6_FINAL_BOOTSTRAP.RData"
)


cat("============================================\n")
cat("STEP 8F: FINAL Base6 coefficients\n")
cat("============================================\n")


# ============================================================
# 2. 最新 bootstrap shrinkage factor
# ============================================================

shrinkage_final <-
  final_summary %>%
  dplyr::pull(
    corrected_slope
  )


cat("\n============================================\n")
cat("A. FINAL bootstrap performance\n")
cat("============================================\n")

print(
  final_summary,
  width = Inf
)


cat(
  "\nFINAL shrinkage factor =",
  shrinkage_final,
  "\n"
)


# ============================================================
# 3. 冻结 FINAL spline settings
# ============================================================

cat("\n============================================\n")
cat("B. FINAL NWR spline settings\n")
cat("============================================\n")

cat(
  "Internal knots:",
  paste(
    nwr_knots_final,
    collapse = ", "
  ),
  "\n"
)

cat(
  "Boundary knots:",
  paste(
    nwr_boundary_final,
    collapse = ", "
  ),
  "\n"
)


# ============================================================
# 4. 数据准备函数
# ============================================================

prepare_data <- function(d) {

  d <- d %>%
    mutate(

      # ------------------------------------------------------
      # Age
      # NHANES age >=80 top-coded at 80
      # ------------------------------------------------------

      age80 =
        pmin(
          ridageyr,
          80
        ),


      # ------------------------------------------------------
      # NWR
      # ------------------------------------------------------

      nwr =
        nwr_model,


      # ------------------------------------------------------
      # CKD-EPI 2021 eGFR
      # ------------------------------------------------------

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
  # FINAL frozen NWR spline basis
  # ----------------------------------------------------------

  nb <- splines::ns(

    d$nwr,

    knots =
      nwr_knots_final,

    Boundary.knots =
      nwr_boundary_final,

    intercept =
      FALSE
  )


  d$nwr_ns1 <- nb[, 1]
  d$nwr_ns2 <- nb[, 2]
  d$nwr_ns3 <- nb[, 3]


  return(d)
}


# ============================================================
# 5. FINAL model formula
# ============================================================

base6_formula_final <-

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
# 6. Fit model in all 30 imputations
# ============================================================

M <- imp_final$m


coef_list <- vector(
  "list",
  M
)

vcov_list <- vector(
  "list",
  M
)

fits_final <- vector(
  "list",
  M
)

completed_final <- vector(
  "list",
  M
)


for (i in 1:M) {

  d <- prepare_data(
    complete(
      imp_final,
      i
    )
  )


  completed_final[[i]] <- d


  fit <- coxph(

    base6_formula_final,

    data =
      d,

    ties =
      "efron",

    x =
      TRUE,

    y =
      TRUE
  )


  fits_final[[i]] <- fit

  coef_list[[i]] <-
    coef(
      fit
    )

  vcov_list[[i]] <-
    vcov(
      fit
    )
}


# ============================================================
# 7. 检查系数一致性
# ============================================================

coef_names <- names(
  coef_list[[1]]
)


same_names <- all(
  vapply(

    coef_list,

    function(x) {

      identical(
        names(x),
        coef_names
      )
    },

    logical(1)
  )
)


if (!same_names) {

  stop(
    "Coefficient names differ across imputations."
  )
}


cat("\n============================================\n")
cat("C. Coefficient names\n")
cat("============================================\n")

print(
  coef_names
)


# ============================================================
# 8. Rubin pooling
# ============================================================

Q <- do.call(
  rbind,
  coef_list
)


# pooled beta
beta_pool_final <-
  colMeans(
    Q
  )


# within-imputation covariance
U_bar <-
  Reduce(
    "+",
    vcov_list
  ) /
  M


# between-imputation covariance
B_cov <-
  stats::cov(
    Q
  )


# total variance
T_cov_final <-
  U_bar +
  (1 + 1/M) *
  B_cov


SE_pool <-
  sqrt(
    diag(
      T_cov_final
    )
  )


# ============================================================
# 9. Rubin-pooled coefficient table
# ============================================================

pooled_table_final <- data.frame(

  term =
    coef_names,

  beta =
    as.numeric(
      beta_pool_final
    ),

  SE =
    as.numeric(
      SE_pool
    ),

  HR =
    exp(
      beta_pool_final
    ),

  CI_low =
    exp(
      beta_pool_final -
        1.96 *
        SE_pool
    ),

  CI_high =
    exp(
      beta_pool_final +
        1.96 *
        SE_pool
    ),

  stringsAsFactors =
    FALSE
)


cat("\n============================================\n")
cat("D. FINAL Rubin-pooled coefficients\n")
cat("============================================\n")

print(
  pooled_table_final,
  row.names = FALSE,
  digits = 7
)


# ============================================================
# 10. MI coefficient stability
# ============================================================

coef_stability_final <- data.frame(

  term =
    coef_names,

  pooled_beta =
    as.numeric(
      beta_pool_final
    ),

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
    ),

  stringsAsFactors =
    FALSE
)


cat("\n============================================\n")
cat("E. Coefficient stability across MI\n")
cat("============================================\n")

print(
  coef_stability_final,
  row.names = FALSE,
  digits = 7
)


# ============================================================
# 11. Apply final uniform shrinkage
# ============================================================

beta_shrunk_final <-
  beta_pool_final *
  shrinkage_final


shrunk_table_final <- data.frame(

  term =
    coef_names,

  beta_pooled =
    as.numeric(
      beta_pool_final
    ),

  shrinkage =
    rep(
      shrinkage_final,
      length(
        coef_names
      )
    ),

  beta_final =
    as.numeric(
      beta_shrunk_final
    ),

  exp_beta_final =
    exp(
      beta_shrunk_final
    ),

  stringsAsFactors =
    FALSE
)


cat("\n============================================\n")
cat("F. FINAL shrunk coefficients\n")
cat("============================================\n")

print(
  shrunk_table_final,
  row.names = FALSE,
  digits = 7
)


# ============================================================
# 12. Final model definition
# ============================================================

model_definition_final <- list(

  model_name =
    "Simplified 6-variable HF prognostic model",

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
      beta_shrunk_final
    ),

  age =
    "Age in years, NHANES capped at 80",

  BMI =
    "Body mass index, kg/m^2",

  RDW =
    "Red cell distribution width",

  NWR =
    "Absolute neutrophil count / total WBC",

  NWR_form =
    "Natural cubic spline, 3 df",

  NWR_internal_knots =
    nwr_knots_final,

  NWR_boundary_knots =
    nwr_boundary_final,

  eGFR =
    "CKD-EPI 2021 creatinine equation",

  medication_count =
    "Number of prescription medications",

  shrinkage_factor =
    shrinkage_final,

  corrected_C =
    final_summary$corrected_C
)


# ============================================================
# 13. Save
# ============================================================

write.csv(

  pooled_table_final,

  "Base6_FINAL_Rubin_coefficients.csv",

  row.names =
    FALSE
)


write.csv(

  coef_stability_final,

  "Base6_FINAL_coefficient_stability.csv",

  row.names =
    FALSE
)


write.csv(

  shrunk_table_final,

  "Base6_FINAL_shrunk_coefficients.csv",

  row.names =
    FALSE
)


save(

  beta_pool_final,

  beta_shrunk_final,

  T_cov_final,

  pooled_table_final,

  coef_stability_final,

  shrunk_table_final,

  model_definition_final,

  base6_formula_final,

  completed_final,

  fits_final,

  nwr_knots_final,

  nwr_boundary_final,

  shrinkage_final,

  file =
    "NHANES_Base6_FINAL_COEFFICIENTS_v2.RData"
)


cat("\n============================================\n")
cat("STEP 8F 完成\n")
cat("============================================\n")
