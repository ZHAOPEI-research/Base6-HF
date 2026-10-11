# ============================================================
# STEP 11C
# Paired bootstrap comparison:
# Base6 vs MAGGIC vs GWTG-HF
#
# Outcomes:
#   1. 364-day time-dependent AUC
#   2. Harrell C-index
#
# IMPORTANT:
#   - All three models evaluated in exactly the same patients
#   - Paired bootstrap: same resampled patient indices
#     used for all three models
#   - Data source:
#     Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData
#   - Common cohort object:
#     comparison3
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

library(survival)
library(timeROC)


# Public-repository path configuration
# Restricted external-cohort data and derived individual-level outputs
# must remain outside the public repository.
.base6_external_output_dir <- Sys.getenv("BASE6_EXTERNAL_OUTPUT_DIR")

if (.base6_external_output_dir == "") {
  stop(
    paste0(
      "Set BASE6_EXTERNAL_OUTPUT_DIR to an authorized local directory ",
      "for restricted external-cohort analysis outputs."
    )
  )
}

if (!dir.exists(.base6_external_output_dir)) {
  dir.create(
    .base6_external_output_dir,
    recursive = TRUE
  )
}

setwd(.base6_external_output_dir)

set.seed(20260905)

B <- 2000


# ============================================================
# 1. Clean environment and load latest STEP 11B
# ============================================================

rm(list = setdiff(ls(), c("B")))

load("Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData")

if (!exists("comparison3")) {
  stop(
    "ERROR: comparison3 not found in ",
    "Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData"
  )
}

dat_boot <- comparison3


cat("\n")
cat("============================================\n")
cat("STEP 11C: Paired bootstrap comparison\n")
cat("Base6 vs MAGGIC vs GWTG-HF\n")
cat("============================================\n\n")


# ============================================================
# 2. Cohort audit
# ============================================================

required_cols <- c(
  "ID",
  "status",
  "time",
  "base6_marker",
  "maggic_score",
  "gwtg_score"
)

missing_required <- setdiff(
  required_cols,
  names(dat_boot)
)

if (length(missing_required) > 0) {

  stop(
    paste0(
      "Missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  )
}


cat("============================================\n")
cat("A. Common cohort audit\n")
cat("============================================\n\n")

cat("N =", nrow(dat_boot), "\n")
cat(
  "Deaths =",
  sum(dat_boot$status == 1, na.rm = TRUE),
  "\n"
)
cat(
  "Censored/alive =",
  sum(dat_boot$status == 0, na.rm = TRUE),
  "\n\n"
)

cat("Follow-up time summary (days):\n")
print(summary(dat_boot$time))
cat("\n")

cat("Status distribution:\n")
print(table(dat_boot$status, useNA = "ifany"))
cat("\n")


# ------------------------------------------------------------
# Expected latest STEP 11B cohort
# ------------------------------------------------------------

expected_n <- 346
expected_deaths <- 47

if (nrow(dat_boot) != expected_n) {

  warning(
    paste0(
      "Cohort size differs from latest STEP 11B: ",
      "observed n = ",
      nrow(dat_boot),
      ", expected n = ",
      expected_n
    )
  )
}

if (
  sum(dat_boot$status == 1, na.rm = TRUE) !=
  expected_deaths
) {

  warning(
    paste0(
      "Death count differs from latest STEP 11B: ",
      "observed deaths = ",
      sum(dat_boot$status == 1, na.rm = TRUE),
      ", expected deaths = ",
      expected_deaths
    )
  )
}


# ============================================================
# 3. Missingness and validity audit
# ============================================================

analysis_cols <- c(
  "time",
  "status",
  "base6_marker",
  "maggic_score",
  "gwtg_score"
)

complete_analysis <- complete.cases(
  dat_boot[, analysis_cols]
)

cat("============================================\n")
cat("B. Analysis completeness\n")
cat("============================================\n\n")

cat(
  "Complete rows =",
  sum(complete_analysis),
  "/",
  nrow(dat_boot),
  "\n"
)

cat(
  "Incomplete rows =",
  sum(!complete_analysis),
  "\n\n"
)

if (!all(complete_analysis)) {

  cat("Rows with missing values:\n")

  print(
    dat_boot[
      !complete_analysis,
      required_cols,
      drop = FALSE
    ]
  )

  stop(
    "Missing values found in the common three-model cohort. ",
    "Return to STEP 11B before bootstrap."
  )
}


# ------------------------------------------------------------
# Basic validity checks
# ------------------------------------------------------------

if (!all(dat_boot$status %in% c(0, 1))) {
  stop("ERROR: status contains values other than 0/1.")
}

if (any(dat_boot$time <= 0)) {
  stop("ERROR: follow-up time contains values <= 0.")
}

if (any(dat_boot$time > 365)) {
  warning(
    "Some follow-up times exceed 365 days."
  )
}

if (anyDuplicated(dat_boot$ID) > 0) {

  cat("Duplicated IDs found:\n")

  print(
    dat_boot$ID[
      duplicated(dat_boot$ID)
    ]
  )

  stop(
    "ERROR: duplicate patient IDs found in comparison3."
  )
}


cat("All required variables complete.\n")
cat("No duplicated patient IDs.\n")
cat("Outcome coding valid.\n\n")


# ============================================================
# 4. Marker summaries
# ============================================================

cat("============================================\n")
cat("C. Marker summaries\n")
cat("============================================\n\n")

cat("Base6 marker:\n")
print(summary(dat_boot$base6_marker))
cat("\n")

cat("MAGGIC score:\n")
print(summary(dat_boot$maggic_score))
cat("\n")

cat("GWTG-HF score:\n")
print(summary(dat_boot$gwtg_score))
cat("\n")


# ============================================================
# 5. Metric functions
# ============================================================

calc_cindex <- function(time, status, marker) {

  fit <- survival::concordance(
    survival::Surv(time, status) ~ marker,
    reverse = TRUE
  )

  as.numeric(
    fit$concordance
  )
}


calc_auc_364 <- function(time, status, marker) {

  roc_obj <- timeROC::timeROC(
    T = time,
    delta = status,
    marker = marker,
    cause = 1,
    weighting = "marginal",
    times = 364,
    iid = FALSE
  )

  # timeROC may internally add t = 0 when only one
  # evaluation time is requested.
  # Therefore explicitly extract AUC at t = 364.
  idx <- which.min(
    abs(roc_obj$times - 364)
  )

  auc <- as.numeric(
    roc_obj$AUC[idx]
  )

  if (!is.finite(auc)) {
    stop("AUC at 364 days could not be calculated.")
  }

  auc
}


# ============================================================
# 6. Original point estimates
# ============================================================

auc_base6 <- calc_auc_364(
  dat_boot$time,
  dat_boot$status,
  dat_boot$base6_marker
)

auc_maggic <- calc_auc_364(
  dat_boot$time,
  dat_boot$status,
  dat_boot$maggic_score
)

auc_gwtg <- calc_auc_364(
  dat_boot$time,
  dat_boot$status,
  dat_boot$gwtg_score
)


c_base6 <- calc_cindex(
  dat_boot$time,
  dat_boot$status,
  dat_boot$base6_marker
)

c_maggic <- calc_cindex(
  dat_boot$time,
  dat_boot$status,
  dat_boot$maggic_score
)

c_gwtg <- calc_cindex(
  dat_boot$time,
  dat_boot$status,
  dat_boot$gwtg_score
)


original_metrics <- data.frame(

  model = c(
    "Base6",
    "MAGGIC",
    "GWTG-HF"
  ),

  AUC_364 = c(
    auc_base6,
    auc_maggic,
    auc_gwtg
  ),

  Harrell_C = c(
    c_base6,
    c_maggic,
    c_gwtg
  ),

  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("D. Original point estimates\n")
cat("============================================\n\n")

print(
  original_metrics,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 7. Compare against STEP 11B
# ============================================================

cat("============================================\n")
cat("E. Consistency check with STEP 11B\n")
cat("============================================\n\n")

if (exists("performance_summary")) {

  cat("STEP 11B performance_summary:\n")

  print(
    performance_summary,
    row.names = FALSE
  )

  cat("\n")

} else {

  cat(
    "performance_summary not available; ",
    "using directly recalculated estimates.\n\n"
  )
}


# ------------------------------------------------------------
# Expected values from latest STEP 11B
# ------------------------------------------------------------

expected_auc <- c(
  Base6 = 0.695213,
  MAGGIC = 0.649478,
  GWTG = 0.665513
)

expected_c <- c(
  Base6 = 0.708190,
  MAGGIC = 0.626921,
  GWTG = 0.651800
)

observed_auc <- c(
  Base6 = auc_base6,
  MAGGIC = auc_maggic,
  GWTG = auc_gwtg
)

observed_c <- c(
  Base6 = c_base6,
  MAGGIC = c_maggic,
  GWTG = c_gwtg
)


consistency_check <- data.frame(

  model = c(
    "Base6",
    "MAGGIC",
    "GWTG-HF"
  ),

  AUC_STEP11B = expected_auc,

  AUC_STEP11C = observed_auc,

  AUC_difference =
    observed_auc -
    expected_auc,

  C_STEP11B = expected_c,

  C_STEP11C = observed_c,

  C_difference =
    observed_c -
    expected_c
)


print(
  consistency_check,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 8. Paired bootstrap
# ============================================================

n <- nrow(dat_boot)

boot_results <- matrix(
  NA_real_,
  nrow = B,
  ncol = 6
)

colnames(boot_results) <- c(
  "AUC_Base6",
  "AUC_MAGGIC",
  "AUC_GWTG",
  "C_Base6",
  "C_MAGGIC",
  "C_GWTG"
)


cat("============================================\n")
cat("F. Paired bootstrap running\n")
cat("============================================\n\n")

cat("Bootstrap repetitions =", B, "\n\n")


for (b in seq_len(B)) {

  # ----------------------------------------------------------
  # Same resampled patient indices for all three models
  # ----------------------------------------------------------

  idx <- sample(
    seq_len(n),
    size = n,
    replace = TRUE
  )

  d <- dat_boot[
    idx,
    ,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Skip pathological resamples
  # ----------------------------------------------------------

  if (
    length(unique(d$status)) < 2
  ) {
    next
  }


  # ----------------------------------------------------------
  # Time-dependent AUC
  # ----------------------------------------------------------

  boot_results[
    b,
    "AUC_Base6"
  ] <- tryCatch(

    calc_auc_364(
      d$time,
      d$status,
      d$base6_marker
    ),

    error = function(e) {
      NA_real_
    }
  )


  boot_results[
    b,
    "AUC_MAGGIC"
  ] <- tryCatch(

    calc_auc_364(
      d$time,
      d$status,
      d$maggic_score
    ),

    error = function(e) {
      NA_real_
    }
  )


  boot_results[
    b,
    "AUC_GWTG"
  ] <- tryCatch(

    calc_auc_364(
      d$time,
      d$status,
      d$gwtg_score
    ),

    error = function(e) {
      NA_real_
    }
  )


  # ----------------------------------------------------------
  # Harrell C-index
  # ----------------------------------------------------------

  boot_results[
    b,
    "C_Base6"
  ] <- tryCatch(

    calc_cindex(
      d$time,
      d$status,
      d$base6_marker
    ),

    error = function(e) {
      NA_real_
    }
  )


  boot_results[
    b,
    "C_MAGGIC"
  ] <- tryCatch(

    calc_cindex(
      d$time,
      d$status,
      d$maggic_score
    ),

    error = function(e) {
      NA_real_
    }
  )


  boot_results[
    b,
    "C_GWTG"
  ] <- tryCatch(

    calc_cindex(
      d$time,
      d$status,
      d$gwtg_score
    ),

    error = function(e) {
      NA_real_
    }
  )


  # ----------------------------------------------------------
  # Progress display
  # ----------------------------------------------------------

  if (b %% 200 == 0) {

    cat(
      "Bootstrap",
      b,
      "/",
      B,
      "\n"
    )
  }
}


boot_results <- as.data.frame(
  boot_results
)


cat("\n")
cat("Valid bootstrap runs:\n")

print(
  colSums(
    is.finite(
      as.matrix(
        boot_results
      )
    )
  )
)

cat("\n")


# ============================================================
# 9. Bootstrap CI for each model
# ============================================================

ci_fun <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (
    length(x) == 0
  ) {

    return(
      c(
        lower = NA_real_,
        median = NA_real_,
        upper = NA_real_
      )
    )
  }


  c(

    lower = unname(
      quantile(
        x,
        probs = 0.025,
        na.rm = TRUE
      )
    ),

    median = unname(
      quantile(
        x,
        probs = 0.50,
        na.rm = TRUE
      )
    ),

    upper = unname(
      quantile(
        x,
        probs = 0.975,
        na.rm = TRUE
      )
    )
  )
}


model_ci <- rbind(

  Base6_AUC =
    ci_fun(
      boot_results$AUC_Base6
    ),

  MAGGIC_AUC =
    ci_fun(
      boot_results$AUC_MAGGIC
    ),

  GWTG_AUC =
    ci_fun(
      boot_results$AUC_GWTG
    ),

  Base6_C =
    ci_fun(
      boot_results$C_Base6
    ),

  MAGGIC_C =
    ci_fun(
      boot_results$C_MAGGIC
    ),

  GWTG_C =
    ci_fun(
      boot_results$C_GWTG
    )
)


cat("============================================\n")
cat("G. Bootstrap 95% CI for each model\n")
cat("============================================\n\n")

print(
  round(
    model_ci,
    4
  )
)

cat("\n")


# ============================================================
# 10. Paired bootstrap differences
# ============================================================

boot_diff <- data.frame(

  dAUC_Base6_MAGGIC =
    boot_results$AUC_Base6 -
    boot_results$AUC_MAGGIC,

  dAUC_Base6_GWTG =
    boot_results$AUC_Base6 -
    boot_results$AUC_GWTG,

  dAUC_GWTG_MAGGIC =
    boot_results$AUC_GWTG -
    boot_results$AUC_MAGGIC,


  dC_Base6_MAGGIC =
    boot_results$C_Base6 -
    boot_results$C_MAGGIC,

  dC_Base6_GWTG =
    boot_results$C_Base6 -
    boot_results$C_GWTG,

  dC_GWTG_MAGGIC =
    boot_results$C_GWTG -
    boot_results$C_MAGGIC
)


# ============================================================
# 11. Paired difference CI and bootstrap P-value
# ============================================================

diff_summary_fun <- function(
    x,
    point_estimate
) {

  x <- x[
    is.finite(x)
  ]

  if (
    length(x) == 0
  ) {

    return(
      c(
        point_estimate = point_estimate,
        bootstrap_mean = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p_boot = NA_real_,
        valid_B = 0
      )
    )
  }


  ci <- quantile(
    x,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE
  )


  # ----------------------------------------------------------
  # Two-sided bootstrap tail-probability P value
  # ----------------------------------------------------------

  p_boot <- 2 * min(
    mean(x <= 0),
    mean(x >= 0)
  )

  p_boot <- min(
    p_boot,
    1
  )


  c(

    point_estimate =
      point_estimate,

    bootstrap_mean =
      mean(
        x,
        na.rm = TRUE
      ),

    CI_low =
      unname(
        ci[1]
      ),

    CI_high =
      unname(
        ci[2]
      ),

    p_boot =
      p_boot,

    valid_B =
      length(x)
  )
}


# ------------------------------------------------------------
# Original paired differences
# ------------------------------------------------------------

delta_auc_base6_maggic <-
  auc_base6 -
  auc_maggic

delta_auc_base6_gwtg <-
  auc_base6 -
  auc_gwtg

delta_auc_gwtg_maggic <-
  auc_gwtg -
  auc_maggic


delta_c_base6_maggic <-
  c_base6 -
  c_maggic

delta_c_base6_gwtg <-
  c_base6 -
  c_gwtg

delta_c_gwtg_maggic <-
  c_gwtg -
  c_maggic


diff_summary <- rbind(

  "Base6 - MAGGIC AUC" =
    diff_summary_fun(
      boot_diff$dAUC_Base6_MAGGIC,
      delta_auc_base6_maggic
    ),

  "Base6 - GWTG-HF AUC" =
    diff_summary_fun(
      boot_diff$dAUC_Base6_GWTG,
      delta_auc_base6_gwtg
    ),

  "GWTG-HF - MAGGIC AUC" =
    diff_summary_fun(
      boot_diff$dAUC_GWTG_MAGGIC,
      delta_auc_gwtg_maggic
    ),

  "Base6 - MAGGIC C" =
    diff_summary_fun(
      boot_diff$dC_Base6_MAGGIC,
      delta_c_base6_maggic
    ),

  "Base6 - GWTG-HF C" =
    diff_summary_fun(
      boot_diff$dC_Base6_GWTG,
      delta_c_base6_gwtg
    ),

  "GWTG-HF - MAGGIC C" =
    diff_summary_fun(
      boot_diff$dC_GWTG_MAGGIC,
      delta_c_gwtg_maggic
    )
)


cat("============================================\n")
cat("H. Paired bootstrap differences\n")
cat("============================================\n\n")

print(
  round(
    diff_summary,
    4
  )
)

cat("\n")


# ============================================================
# 12. Publication-style model summary
# ============================================================

final_model_summary <- data.frame(

  model = c(
    "Base6",
    "MAGGIC",
    "GWTG-HF"
  ),

  n = nrow(dat_boot),

  deaths =
    sum(
      dat_boot$status == 1
    ),

  AUC_364 = c(
    auc_base6,
    auc_maggic,
    auc_gwtg
  ),

  AUC_CI_low = c(
    model_ci[
      "Base6_AUC",
      "lower"
    ],
    model_ci[
      "MAGGIC_AUC",
      "lower"
    ],
    model_ci[
      "GWTG_AUC",
      "lower"
    ]
  ),

  AUC_CI_high = c(
    model_ci[
      "Base6_AUC",
      "upper"
    ],
    model_ci[
      "MAGGIC_AUC",
      "upper"
    ],
    model_ci[
      "GWTG_AUC",
      "upper"
    ]
  ),

  Harrell_C = c(
    c_base6,
    c_maggic,
    c_gwtg
  ),

  C_CI_low = c(
    model_ci[
      "Base6_C",
      "lower"
    ],
    model_ci[
      "MAGGIC_C",
      "lower"
    ],
    model_ci[
      "GWTG_C",
      "lower"
    ]
  ),

  C_CI_high = c(
    model_ci[
      "Base6_C",
      "upper"
    ],
    model_ci[
      "MAGGIC_C",
      "upper"
    ],
    model_ci[
      "GWTG_C",
      "upper"
    ]
  ),

  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("I. Final model performance summary\n")
cat("============================================\n\n")

final_model_summary_print <- final_model_summary

num_cols <- sapply(
  final_model_summary_print,
  is.numeric
)

final_model_summary_print[num_cols] <-
  lapply(
    final_model_summary_print[num_cols],
    round,
    digits = 4
  )

print(
  final_model_summary_print,
  row.names = FALSE
)
cat("\n")


# ============================================================
# 13. Publication-style pairwise summary
# ============================================================

final_pairwise_summary <- data.frame(

  comparison =
    rownames(
      diff_summary
    ),

  metric = c(
    "AUC_364",
    "AUC_364",
    "AUC_364",
    "Harrell_C",
    "Harrell_C",
    "Harrell_C"
  ),

  estimate =
    diff_summary[
      ,
      "point_estimate"
    ],

  CI_low =
    diff_summary[
      ,
      "CI_low"
    ],

  CI_high =
    diff_summary[
      ,
      "CI_high"
    ],

  p_boot =
    diff_summary[
      ,
      "p_boot"
    ],

  valid_B =
    diff_summary[
      ,
      "valid_B"
    ],

  row.names = NULL,
  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("J. Final pairwise comparison summary\n")
cat("============================================\n\n")

final_pairwise_summary_print <-
  final_pairwise_summary

num_cols <- sapply(
  final_pairwise_summary_print,
  is.numeric
)

final_pairwise_summary_print[num_cols] <-
  lapply(
    final_pairwise_summary_print[num_cols],
    round,
    digits = 4
  )

print(
  final_pairwise_summary_print,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 14. Save CSV outputs
# ============================================================

write.csv(
  original_metrics,
  "Chinese_GWTG_STEP11C_POINT_ESTIMATES.csv",
  row.names = FALSE
)


write.csv(
  final_model_summary,
  "Chinese_GWTG_STEP11C_MODEL_BOOTSTRAP_CI.csv",
  row.names = FALSE
)


write.csv(
  final_pairwise_summary,
  "Chinese_GWTG_STEP11C_PAIRED_BOOTSTRAP.csv",
  row.names = FALSE
)


write.csv(
  consistency_check,
  "Chinese_GWTG_STEP11C_STEP11B_CONSISTENCY_CHECK.csv",
  row.names = FALSE
)


# ============================================================
# 15. Save complete RData
# ============================================================

save(
  dat_boot,
  original_metrics,
  consistency_check,
  boot_results,
  boot_diff,
  model_ci,
  diff_summary,
  final_model_summary,
  final_pairwise_summary,
  file =
    "Chinese_GWTG_STEP11C_PAIRED_BOOTSTRAP.RData"
)


# ============================================================
# 16. Completion message
# ============================================================

cat("============================================\n")
cat("STEP 11C completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11C_POINT_ESTIMATES.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11C_MODEL_BOOTSTRAP_CI.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11C_PAIRED_BOOTSTRAP.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11C_STEP11B_CONSISTENCY_CHECK.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11C_PAIRED_BOOTSTRAP.RData\n"
)

cat("\n")
cat("Bootstrap repetitions =", B, "\n")
cat("Final common cohort n =", nrow(dat_boot), "\n")
cat(
  "Deaths =",
  sum(dat_boot$status == 1),
  "\n"
)

cat("============================================\n")
