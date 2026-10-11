# ============================================================
# STEP 11D
# NHANES Base6 complex-survey weighted sensitivity analysis
#
# Purpose:
#   Compare the FINAL unweighted Base6 coefficients with
#   coefficients estimated using NHANES complex survey design.
#
# IMPORTANT:
#   - No variable selection
#   - No new spline selection
#   - No new shrinkage
#   - No baseline recalibration
#   - Fixed final NWR spline knots
#   - 30 final multiple imputations
# ============================================================

rm(list = ls())
gc()

library(mice)
library(survival)
library(survey)
library(splines)
library(dplyr)
library(tibble)


analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)

if (!dir.exists(analysis_dir)) {
  dir.create(
    analysis_dir,
    recursive = TRUE
  )
}

setwd(analysis_dir)

options(survey.lonely.psu = "adjust")

cat("============================================\n")
cat("STEP 11D: NHANES survey-weight sensitivity\n")
cat("============================================\n\n")


# ============================================================
# 1. Load FINAL MICE object
# ============================================================

load(
  file.path(
    analysis_dir,
    "NHANES_2011_2018_MICE_FINAL.RData"
  )
)

stopifnot(
  exists("imp_final"),
  inherits(imp_final, "mids"),
  imp_final$m == 30
)

cat("FINAL MICE object loaded.\n")
cat("Number of imputations =", imp_final$m, "\n\n")


# ============================================================
# 2. Frozen FINAL model constants
# ============================================================

# Rubin-pooled coefficients BEFORE shrinkage
beta_unweighted <- c(
  age80     =  0.038175446,
  bmxbmi    = -0.020972564,
  lbxrdw    =  0.203942007,
  nwr_ns1   =  0.729017806,
  nwr_ns2   =  0.177541263,
  nwr_ns3   =  0.513577990,
  egfr_2021 = -0.009750497,
  med_cnt   =  0.057526831
)

# Original optimism-based shrinkage factor
shrinkage_original <- 0.9212384

# FINAL frozen NWR spline
internal_knots <- c(
  0.570093457943925,
  0.652173913043478
)

boundary_knots <- c(
  0.434782608695652,
  0.773607367475292
)


# ============================================================
# 3. Required-variable audit - CORRECTED
# ============================================================

required_vars <- c(
  "seqn",
  "ridageyr",
  "female",
  "death",
  "followup_month",
  "bmxbmi",
  "lbxrdw",
  "nwr_model",
  "scr_mgdl",
  "med_cnt",
  "wtmec8yr",
  "sdmvstra",
  "sdmvpsu"
)

missing_vars <- setdiff(
  required_vars,
  names(imp_final$data)
)

cat("============================================\n")
cat("A. Required-variable check\n")
cat("============================================\n")

if (length(missing_vars) == 0) {

  cat("All required variables are present.\n")

} else {

  cat(
    "Missing variables:",
    paste(
      missing_vars,
      collapse = ", "
    ),
    "\n"
  )

  cat("\nAvailable variables:\n")

  print(
    names(
      imp_final$data
    )
  )

  stop(
    "Required variables missing."
  )
}


# ============================================================
# 4. Survey-design audit
# ============================================================

audit_dat <- imp_final$data


cat("\n============================================\n")
cat("B. Survey-design audit\n")
cat("============================================\n")

cat(
  "N =",
  nrow(audit_dat),
  "\n"
)

cat(
  "Missing WTMEC8YR =",
  sum(
    is.na(
      audit_dat$wtmec8yr
    )
  ),
  "\n"
)

cat(
  "Non-positive WTMEC8YR =",
  sum(
    !is.na(
      audit_dat$wtmec8yr
    ) &
      audit_dat$wtmec8yr <= 0
  ),
  "\n"
)

cat(
  "Missing SDMVSTRA =",
  sum(
    is.na(
      audit_dat$sdmvstra
    )
  ),
  "\n"
)

cat(
  "Missing SDMVPSU =",
  sum(
    is.na(
      audit_dat$sdmvpsu
    )
  ),
  "\n"
)

cat(
  "Unique strata =",
  length(
    unique(
      na.omit(
        audit_dat$sdmvstra
      )
    )
  ),
  "\n"
)

cat(
  "Unique stratum-PSU combinations =",
  nrow(
    unique(
      audit_dat[
        complete.cases(
          audit_dat[
            c(
              "sdmvstra",
              "sdmvpsu"
            )
          ]
        ),
        c(
          "sdmvstra",
          "sdmvpsu"
        )
      ]
    )
  ),
  "\n"
)

cat("\nWTMEC8YR distribution:\n")

print(
  summary(
    audit_dat$wtmec8yr
  )
)


# ============================================================
# 5. CKD-EPI 2021 helper
# ============================================================

calc_egfr_2021 <- function(
    scr,
    age,
    female
) {

  kappa <- ifelse(
    female == 1,
    0.7,
    0.9
  )

  alpha <- ifelse(
    female == 1,
    -0.241,
    -0.302
  )

  142 *
    pmin(
      scr / kappa,
      1
    )^alpha *
    pmax(
      scr / kappa,
      1
    )^(-1.200) *
    0.9938^age *
    ifelse(
      female == 1,
      1.012,
      1
    )
}


# ============================================================
# 6. Construct FINAL analysis dataset
# ============================================================

make_analysis_data <- function(i) {

  d <- mice::complete(
    imp_final,
    action = i
  )


  # ----------------------------------------------------------
  # Age used in final prognostic model
  #
  # NHANES public age is top-coded at 80.
  # Therefore this recreates the same age80 predictor.
  # ----------------------------------------------------------

  d$age80 <- pmin(
    d$ridageyr,
    80
  )


  # ----------------------------------------------------------
  # NWR
  # ----------------------------------------------------------

  d$nwr <-
    d$nwr_model


  # ----------------------------------------------------------
  # eGFR CKD-EPI 2021
  #
  # ridageyr is the age available in NHANES.
  # ----------------------------------------------------------

  d$egfr_2021 <- calc_egfr_2021(

    scr =
      d$scr_mgdl,

    age =
      d$ridageyr,

    female =
      d$female
  )


  # ----------------------------------------------------------
  # Fixed FINAL NWR spline
  # ----------------------------------------------------------

  ns_mat <- splines::ns(

    d$nwr,

    knots =
      internal_knots,

    Boundary.knots =
      boundary_knots,

    intercept =
      FALSE
  )


  d$nwr_ns1 <-
    ns_mat[, 1]

  d$nwr_ns2 <-
    ns_mat[, 2]

  d$nwr_ns3 <-
    ns_mat[, 3]


  # ----------------------------------------------------------
  # Survival outcome
  # ----------------------------------------------------------

  d$time_month <-
    d$followup_month

  d$event <-
    d$death


  d
}


# ============================================================
# 7. First-imputation construction audit
# ============================================================

d1 <- make_analysis_data(1)

cat("\n============================================\n")
cat("C. Derived-variable audit: MI 1\n")
cat("============================================\n")

audit_predictors <- c(
  "age80",
  "bmxbmi",
  "lbxrdw",
  "nwr",
  "egfr_2021",
  "med_cnt",
  "nwr_ns1",
  "nwr_ns2",
  "nwr_ns3",
  "wtmec8yr"
)

for (v in audit_predictors) {

  cat("\n---", v, "---\n")
  print(summary(d1[[v]]))

  cat(
    "Missing:",
    sum(is.na(d1[[v]])),
    "\n"
  )
}


# ============================================================
# 8. Cox formula
# ============================================================

cox_formula <- as.formula(
  paste(
    "Surv(time_month, event) ~",
    paste(
      names(beta_unweighted),
      collapse = " + "
    )
  )
)


# ============================================================
# 9. Survey-weighted Cox across 30 MI
# ============================================================

coef_list <- vector(
  "list",
  imp_final$m
)

vcov_list <- vector(
  "list",
  imp_final$m
)

lp_cor <- rep(
  NA_real_,
  imp_final$m
)

survey_fit_list <- vector(
  "list",
  imp_final$m
)

cat("\n============================================\n")
cat("D. Survey-weighted Cox fitting\n")
cat("============================================\n")

for (i in seq_len(imp_final$m)) {

  cat(
    "MI",
    i,
    "/",
    imp_final$m,
    "\n"
  )

  d <- make_analysis_data(i)

  # Keep complete survey-design + model rows
  vars_needed <- c(
    "time_month",
    "event",
    names(beta_unweighted),
    "wtmec8yr",
    "sdmvstra",
    "sdmvpsu"
  )

  d_fit <- d[
    complete.cases(d[, vars_needed]),
  ]

  if (nrow(d_fit) == 0) {
    stop(
      paste(
        "No complete observations in MI",
        i
      )
    )
  }

  # ----------------------------------------------------------
  # NHANES complex survey design
  # ----------------------------------------------------------

  des <- svydesign(
    ids = ~sdmvpsu,
    strata = ~sdmvstra,
    weights = ~wtmec8yr,
    nest = TRUE,
    data = d_fit
  )

  # ----------------------------------------------------------
  # Survey-weighted Cox
  # ----------------------------------------------------------

  fit <- svycoxph(
    cox_formula,
    design = des
  )

  survey_fit_list[[i]] <- fit

  b <- coef(fit)
  V <- vcov(fit)

  # Ensure identical coefficient ordering
  b <- b[names(beta_unweighted)]

  V <- V[
    names(beta_unweighted),
    names(beta_unweighted),
    drop = FALSE
  ]

  coef_list[[i]] <- b
  vcov_list[[i]] <- V

  # ----------------------------------------------------------
  # LP correlation:
  # unweighted pooled beta vs survey beta
  # ----------------------------------------------------------

  X <- as.matrix(
    d_fit[, names(beta_unweighted)]
  )

  lp_original <- as.numeric(
    X %*% beta_unweighted
  )

  lp_weighted <- as.numeric(
    X %*% b
  )

  lp_cor[i] <- cor(
    lp_original,
    lp_weighted,
    use = "complete.obs",
    method = "pearson"
  )
}


# ============================================================
# 10. Rubin pooling
# ============================================================

terms <- names(beta_unweighted)

Q <- do.call(
  rbind,
  coef_list
)

Qbar <- colMeans(Q)

# Within-imputation covariance
Ubar <- Reduce(
  "+",
  vcov_list
) / length(vcov_list)

# Between-imputation covariance
Bmat <- cov(Q)

m <- imp_final$m

# Rubin total covariance
Tmat <- Ubar +
  (1 + 1 / m) * Bmat

SE_pool <- sqrt(
  diag(Tmat)
)

beta_weighted <- Qbar


# ============================================================
# 11. Survey-weighted pooled coefficients
# ============================================================

weighted_coef_table <- tibble(
  term = terms,
  beta = beta_weighted,
  SE = SE_pool,
  HR = exp(beta_weighted),
  CI_low = exp(
    beta_weighted -
      1.96 * SE_pool
  ),
  CI_high = exp(
    beta_weighted +
      1.96 * SE_pool
  )
)

cat("\n============================================\n")
cat("E. Rubin-pooled survey-weighted coefficients\n")
cat("============================================\n")

print(
  as.data.frame(weighted_coef_table),
  row.names = FALSE
)


# ============================================================
# 12. Compare weighted vs FINAL unweighted coefficients
# ============================================================

comparison <- tibble(
  term = terms,

  beta_unweighted =
    as.numeric(
      beta_unweighted[terms]
    ),

  beta_weighted =
    as.numeric(
      beta_weighted[terms]
    )
) %>%

  mutate(

    absolute_change =
      beta_weighted -
      beta_unweighted,

    absolute_difference =
      abs(absolute_change),

    relative_change_pct =
      ifelse(
        abs(beta_unweighted) > 1e-8,
        100 *
          absolute_change /
          abs(beta_unweighted),
        NA_real_
      ),

    HR_unweighted =
      exp(beta_unweighted),

    HR_weighted =
      exp(beta_weighted),

    same_direction =
      sign(beta_unweighted) ==
      sign(beta_weighted)
  )


cat("\n============================================\n")
cat("F. Weighted vs unweighted coefficients\n")
cat("============================================\n")

print(
  as.data.frame(comparison),
  row.names = FALSE
)


# ============================================================
# 13. Coefficient stability across MI
# ============================================================

coef_stability <- tibble(
  term = terms,
  pooled_beta = beta_weighted,
  min_beta = apply(
    Q,
    2,
    min
  ),
  max_beta = apply(
    Q,
    2,
    max
  ),
  SD_between_MI = apply(
    Q,
    2,
    sd
  )
)

cat("\n============================================\n")
cat("G. Survey coefficient stability across MI\n")
cat("============================================\n")

print(
  as.data.frame(coef_stability),
  row.names = FALSE
)


# ============================================================
# 14. LP correlation
# ============================================================

cat("\n============================================\n")
cat("H. LP correlation\n")
cat("============================================\n")

cat(
  "Correlation between original and survey-weighted LP:\n"
)

print(
  summary(lp_cor)
)

cat(
  "\nMedian correlation =",
  median(
    lp_cor,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Minimum correlation =",
  min(
    lp_cor,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 15. Direction and magnitude audit
# ============================================================

cat("\n============================================\n")
cat("I. Material-shift audit\n")
cat("============================================\n")

cat(
  "Direction changes =",
  sum(
    !comparison$same_direction,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "|relative change| >20% =",
  sum(
    abs(
      comparison$relative_change_pct
    ) > 20,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "|relative change| >50% =",
  sum(
    abs(
      comparison$relative_change_pct
    ) > 50,
    na.rm = TRUE
  ),
  "\n"
)

cat("\nNWR spline coefficients:\n")

print(
  as.data.frame(
    comparison %>%
      filter(
        term %in% c(
          "nwr_ns1",
          "nwr_ns2",
          "nwr_ns3"
        )
      )
  ),
  row.names = FALSE
)


# ============================================================
# 16. OPTIONAL diagnostic:
#     What if original shrinkage were applied?
#
# This is NOT a new final model.
# Only descriptive.
# ============================================================

comparison_shrunk <- comparison %>%
  mutate(
    original_final_beta =
      beta_unweighted *
      shrinkage_original,

    survey_beta_same_shrinkage =
      beta_weighted *
      shrinkage_original
  )

cat("\n============================================\n")
cat("J. Descriptive same-shrinkage comparison\n")
cat("============================================\n")

cat(
  "NOTE: This does NOT constitute re-validation of the ",
  "survey-weighted model.\n"
)

print(
  as.data.frame(
    comparison_shrunk %>%
      select(
        term,
        original_final_beta,
        survey_beta_same_shrinkage
      )
  ),
  row.names = FALSE
)


# ============================================================
# 17. Save
# ============================================================

save(
  weighted_coef_table,
  comparison,
  coef_stability,
  comparison_shrunk,
  beta_weighted,
  SE_pool,
  Tmat,
  Q,
  lp_cor,
  internal_knots,
  boundary_knots,
  file =
    file.path(
      analysis_dir,
      "NHANES_Base6_STEP11D_SURVEY_SENSITIVITY.RData"
    )
)

write.csv(
  comparison,
  file.path(
    analysis_dir,
    "NHANES_Base6_STEP11D_SURVEY_COEFFICIENT_COMPARISON.csv"
  ),
  row.names = FALSE
)

write.csv(
  weighted_coef_table,
  file.path(
    analysis_dir,
    "NHANES_Base6_STEP11D_SURVEY_WEIGHTED_COEFFICIENTS.csv"
  ),
  row.names = FALSE
)


cat("\n============================================\n")
cat("STEP 11D 完成\n")
cat("============================================\n")
