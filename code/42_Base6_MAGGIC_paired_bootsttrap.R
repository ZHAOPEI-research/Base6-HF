# ============================================================
# STEP 10C
# Base6 vs MAGGIC
# Paired bootstrap head-to-head comparison
#
# Chinese external validation cohort
#
# Primary comparison:
#   Base6 vs MAGGIC original
#
# Sensitivity:
#   Base6 vs MAGGIC modernized
#
# Metrics:
#   364-day cumulative/dynamic AUC
#   365-day IPCW Brier
#
# Bootstrap:
#   patient-level paired bootstrap
#   B = 2000
# ============================================================

rm(list = ls())

library(dplyr)
library(survival)

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

load(
  "Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData"
)

set.seed(20260822)

B <- 2000

AUC_TIME <- 364
BRIER_TIME <- 365


cat("============================================\n")
cat("STEP 10C: Base6 vs MAGGIC paired bootstrap\n")
cat("============================================\n")

cat(
  "Common cohort n =",
  nrow(comparison_data),
  "\n"
)

cat(
  "Deaths =",
  sum(comparison_data$status == 1),
  "\n"
)


# ============================================================
# 1. IPCW AUC
# ============================================================

ipcw_auc <- function(
    time,
    status,
    marker,
    t
) {

  censor_event <- 1 - status

  km_censor <- tryCatch(

    survfit(
      Surv(
        time,
        censor_event
      ) ~ 1
    ),

    error = function(e) NULL
  )

  if (is.null(km_censor)) {
    return(NA_real_)
  }


  G_left <- function(x) {

    out <- tryCatch(

      summary(
        km_censor,
        times =
          pmax(
            x - 1e-8,
            0
          ),
        extend = TRUE
      )$surv,

      error = function(e) NA_real_
    )

    as.numeric(out)
  }


  cases <- which(
    status == 1 &
      time <= t
  )

  controls <- which(
    time > t
  )


  if (
    length(cases) == 0 ||
    length(controls) == 0
  ) {
    return(NA_real_)
  }


  case_weights <- numeric(
    length(cases)
  )


  for (j in seq_along(cases)) {

    G_i <- G_left(
      time[
        cases[j]
      ]
    )

    if (
      length(G_i) != 1 ||
      !is.finite(G_i) ||
      G_i <= 0
    ) {
      return(NA_real_)
    }

    case_weights[j] <- 1 / G_i
  }


  numerator <- 0
  denominator <- 0


  for (j in seq_along(cases)) {

    i <- cases[j]

    difference <-
      marker[i] -
      marker[controls]


    pair_score <- ifelse(
      difference > 0,
      1,
      ifelse(
        difference == 0,
        0.5,
        0
      )
    )


    numerator <-
      numerator +
      case_weights[j] *
      sum(
        pair_score
      )


    denominator <-
      denominator +
      case_weights[j] *
      length(
        controls
      )
  }


  if (
    denominator <= 0 ||
    !is.finite(denominator)
  ) {
    return(NA_real_)
  }


  numerator / denominator
}


# ============================================================
# 2. IPCW Brier
# ============================================================

ipcw_brier <- function(
    time,
    status,
    pred_risk,
    t
) {

  censor_event <- 1 - status


  km_censor <- tryCatch(

    survfit(
      Surv(
        time,
        censor_event
      ) ~ 1
    ),

    error = function(e) NULL
  )


  if (is.null(km_censor)) {
    return(NA_real_)
  }


  G_left <- function(x) {

    out <- tryCatch(

      summary(
        km_censor,
        times =
          pmax(
            x - 1e-8,
            0
          ),
        extend = TRUE
      )$surv,

      error = function(e) NA_real_
    )

    as.numeric(out)
  }


  G_t_left <- G_left(t)


  if (
    length(G_t_left) != 1 ||
    !is.finite(G_t_left) ||
    G_t_left <= 0
  ) {
    return(NA_real_)
  }


  contribution <- numeric(
    length(time)
  )


  # ----------------------------------------------------------
  # Event by t
  # ----------------------------------------------------------

  event_rows <- which(
    status == 1 &
      time <= t
  )


  for (i in event_rows) {

    G_i <- G_left(
      time[i]
    )

    if (
      length(G_i) != 1 ||
      !is.finite(G_i) ||
      G_i <= 0
    ) {
      return(NA_real_)
    }


    contribution[i] <-
      (1 - pred_risk[i])^2 /
      G_i
  }


  # ----------------------------------------------------------
  # Known alive through t
  # ----------------------------------------------------------

  alive_rows <- which(
    status == 0 &
      time >= t
  )


  contribution[
    alive_rows
  ] <-
    pred_risk[
      alive_rows
    ]^2 /
    G_t_left


  mean(
    contribution
  )
}


# ============================================================
# 3. Apparent paired comparison
# ============================================================

auc_base6 <- ipcw_auc(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$lp_base6,
  AUC_TIME
)


auc_maggic <- ipcw_auc(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$maggic_score,
  AUC_TIME
)


auc_maggic_modern <- ipcw_auc(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$maggic_score_modern,
  AUC_TIME
)


brier_base6 <- ipcw_brier(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$risk_1y,
  BRIER_TIME
)


brier_maggic <- ipcw_brier(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$risk_maggic_1y,
  BRIER_TIME
)


brier_maggic_modern <- ipcw_brier(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$risk_maggic_modern_1y,
  BRIER_TIME
)


delta_auc_original <-
  auc_base6 -
  auc_maggic


delta_auc_modern <-
  auc_base6 -
  auc_maggic_modern


# IMPORTANT:
#
# For Brier:
# lower = better
#
# Define:
#
# ΔBrier = Base6 - comparator
#
# Therefore:
# negative ΔBrier favors Base6.
# ============================================================

delta_brier_original <-
  brier_base6 -
  brier_maggic


delta_brier_modern <-
  brier_base6 -
  brier_maggic_modern


cat("\n============================================\n")
cat("A. Apparent paired comparison\n")
cat("============================================\n")

cat(
  "Base6 AUC =",
  auc_base6,
  "\n"
)

cat(
  "MAGGIC original AUC =",
  auc_maggic,
  "\n"
)

cat(
  "MAGGIC modernized AUC =",
  auc_maggic_modern,
  "\n"
)

cat(
  "\nΔAUC Base6 - MAGGIC original =",
  delta_auc_original,
  "\n"
)

cat(
  "ΔAUC Base6 - MAGGIC modernized =",
  delta_auc_modern,
  "\n"
)


cat(
  "\nBase6 Brier =",
  brier_base6,
  "\n"
)

cat(
  "MAGGIC original Brier =",
  brier_maggic,
  "\n"
)

cat(
  "MAGGIC modernized Brier =",
  brier_maggic_modern,
  "\n"
)

cat(
  "\nΔBrier Base6 - MAGGIC original =",
  delta_brier_original,
  "\n"
)

cat(
  "ΔBrier Base6 - MAGGIC modernized =",
  delta_brier_modern,
  "\n"
)


# ============================================================
# 4. Bootstrap storage
# ============================================================

boot_results <- vector(
  "list",
  B
)

n <- nrow(
  comparison_data
)


# ============================================================
# 5. Paired bootstrap loop
# ============================================================

for (b in 1:B) {

  if (
    b == 1 ||
    b %% 100 == 0
  ) {

    cat(
      "Bootstrap",
      b,
      "/",
      B,
      "\n"
    )
  }


  idx <- sample(

    seq_len(n),

    size =
      n,

    replace =
      TRUE
  )


  db <- comparison_data[
    idx,
    ,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Validity
  # ----------------------------------------------------------

  if (
    sum(db$status == 1) == 0 ||
    sum(db$time_day > AUC_TIME) == 0
  ) {

    boot_results[[b]] <- data.frame(

      bootstrap =
        b,

      valid =
        FALSE
    )

    next
  }


  # ----------------------------------------------------------
  # Base6
  # ----------------------------------------------------------

  auc_b6 <- ipcw_auc(

    db$time_day,
    db$status,
    db$lp_base6,
    AUC_TIME
  )


  br_b6 <- ipcw_brier(

    db$time_day,
    db$status,
    db$risk_1y,
    BRIER_TIME
  )


  # ----------------------------------------------------------
  # MAGGIC original
  # ----------------------------------------------------------

  auc_mag <- ipcw_auc(

    db$time_day,
    db$status,
    db$maggic_score,
    AUC_TIME
  )


  br_mag <- ipcw_brier(

    db$time_day,
    db$status,
    db$risk_maggic_1y,
    BRIER_TIME
  )


  # ----------------------------------------------------------
  # MAGGIC modernized
  # ----------------------------------------------------------

  auc_mag_mod <- ipcw_auc(

    db$time_day,
    db$status,
    db$maggic_score_modern,
    AUC_TIME
  )


  br_mag_mod <- ipcw_brier(

    db$time_day,
    db$status,
    db$risk_maggic_modern_1y,
    BRIER_TIME
  )


  values <- c(

    auc_b6,
    br_b6,

    auc_mag,
    br_mag,

    auc_mag_mod,
    br_mag_mod
  )


  if (
    any(
      !is.finite(
        values
      )
    )
  ) {

    boot_results[[b]] <- data.frame(

      bootstrap =
        b,

      valid =
        FALSE
    )

    next
  }


  boot_results[[b]] <- data.frame(

    bootstrap =
      b,

    valid =
      TRUE,

    deaths =
      sum(
        db$status == 1
      ),

    AUC_Base6 =
      auc_b6,

    AUC_MAGGIC =
      auc_mag,

    AUC_MAGGIC_modern =
      auc_mag_mod,

    delta_AUC_original =
      auc_b6 -
      auc_mag,

    delta_AUC_modern =
      auc_b6 -
      auc_mag_mod,


    Brier_Base6 =
      br_b6,

    Brier_MAGGIC =
      br_mag,

    Brier_MAGGIC_modern =
      br_mag_mod,

    delta_Brier_original =
      br_b6 -
      br_mag,

    delta_Brier_modern =
      br_b6 -
      br_mag_mod
  )
}


boot_df <- bind_rows(
  boot_results
)


boot_valid <- boot_df %>%
  filter(
    valid
  )


cat("\n============================================\n")
cat("B. Bootstrap completion\n")
cat("============================================\n")

cat(
  "Valid =",
  nrow(
    boot_valid
  ),
  "/",
  B,
  "\n"
)


# ============================================================
# 6. CI helper
# ============================================================

CI <- function(x) {

  x <- x[
    is.finite(x)
  ]

  c(

    median =
      median(x),

    low =
      quantile(
        x,
        0.025,
        names = FALSE
      ),

    high =
      quantile(
        x,
        0.975,
        names = FALSE
      )
  )
}


# ============================================================
# 7. Bootstrap estimates
# ============================================================

ci_auc_base6 <-
  CI(
    boot_valid$AUC_Base6
  )

ci_auc_maggic <-
  CI(
    boot_valid$AUC_MAGGIC
  )

ci_auc_modern <-
  CI(
    boot_valid$AUC_MAGGIC_modern
  )


ci_delta_auc_original <-
  CI(
    boot_valid$delta_AUC_original
  )

ci_delta_auc_modern <-
  CI(
    boot_valid$delta_AUC_modern
  )


ci_brier_base6 <-
  CI(
    boot_valid$Brier_Base6
  )

ci_brier_maggic <-
  CI(
    boot_valid$Brier_MAGGIC
  )

ci_brier_modern <-
  CI(
    boot_valid$Brier_MAGGIC_modern
  )


ci_delta_brier_original <-
  CI(
    boot_valid$delta_Brier_original
  )

ci_delta_brier_modern <-
  CI(
    boot_valid$delta_Brier_modern
  )


# ============================================================
# 8. Approximate paired bootstrap p-values
#
# Two-sided:
# proportion crossing zero.
# ============================================================

bootstrap_p <- function(x) {

  x <- x[
    is.finite(x)
  ]

  p_left <-
    mean(
      x <= 0
    )

  p_right <-
    mean(
      x >= 0
    )

  min(
    1,
    2 *
      min(
        p_left,
        p_right
      )
  )
}


p_delta_auc_original <-
  bootstrap_p(
    boot_valid$delta_AUC_original
  )


p_delta_auc_modern <-
  bootstrap_p(
    boot_valid$delta_AUC_modern
  )


p_delta_brier_original <-
  bootstrap_p(
    boot_valid$delta_Brier_original
  )


p_delta_brier_modern <-
  bootstrap_p(
    boot_valid$delta_Brier_modern
  )


# ============================================================
# 9. Final paired-comparison table
# ============================================================

paired_summary <- data.frame(

  comparison = c(

    "Base6 vs MAGGIC original - AUC",

    "Base6 vs MAGGIC modernized - AUC",

    "Base6 vs MAGGIC original - Brier",

    "Base6 vs MAGGIC modernized - Brier"
  ),

  estimate = c(

    delta_auc_original,

    delta_auc_modern,

    delta_brier_original,

    delta_brier_modern
  ),

  bootstrap_median = c(

    ci_delta_auc_original["median"],

    ci_delta_auc_modern["median"],

    ci_delta_brier_original["median"],

    ci_delta_brier_modern["median"]
  ),

  CI_low = c(

    ci_delta_auc_original["low"],

    ci_delta_auc_modern["low"],

    ci_delta_brier_original["low"],

    ci_delta_brier_modern["low"]
  ),

  CI_high = c(

    ci_delta_auc_original["high"],

    ci_delta_auc_modern["high"],

    ci_delta_brier_original["high"],

    ci_delta_brier_modern["high"]
  ),

  bootstrap_p = c(

    p_delta_auc_original,

    p_delta_auc_modern,

    p_delta_brier_original,

    p_delta_brier_modern
  )
)


cat("\n============================================\n")
cat("C. Paired bootstrap differences\n")
cat("============================================\n")

print(
  paired_summary,
  row.names = FALSE
)


# ============================================================
# 10. Individual model CI table
# ============================================================

model_CI <- data.frame(

  model = c(

    "Base6",

    "MAGGIC_original",

    "MAGGIC_modernized"
  ),

  AUC = c(

    auc_base6,

    auc_maggic,

    auc_maggic_modern
  ),

  AUC_low = c(

    ci_auc_base6["low"],

    ci_auc_maggic["low"],

    ci_auc_modern["low"]
  ),

  AUC_high = c(

    ci_auc_base6["high"],

    ci_auc_maggic["high"],

    ci_auc_modern["high"]
  ),

  Brier = c(

    brier_base6,

    brier_maggic,

    brier_maggic_modern
  ),

  Brier_low = c(

    ci_brier_base6["low"],

    ci_brier_maggic["low"],

    ci_brier_modern["low"]
  ),

  Brier_high = c(

    ci_brier_base6["high"],

    ci_brier_maggic["high"],

    ci_brier_modern["high"]
  )
)


cat("\n============================================\n")
cat("D. Model performance with bootstrap CI\n")
cat("============================================\n")

print(
  model_CI,
  row.names = FALSE
)


# ============================================================
# 11. Direction audit
# ============================================================

cat("\n============================================\n")
cat("E. Interpretation direction\n")
cat("============================================\n")

cat(
  "Positive ΔAUC favors Base6.\n"
)

cat(
  "Negative ΔBrier favors Base6.\n"
)


cat(
  "\nProportion ΔAUC(original) > 0 =",
  mean(
    boot_valid$delta_AUC_original > 0
  ),
  "\n"
)


cat(
  "Proportion ΔBrier(original) < 0 =",
  mean(
    boot_valid$delta_Brier_original < 0
  ),
  "\n"
)


# ============================================================
# 12. Save
# ============================================================

write.csv(

  paired_summary,

  "Chinese_Base6_MAGGIC_paired_bootstrap.csv",

  row.names = FALSE
)


write.csv(

  model_CI,

  "Chinese_Base6_MAGGIC_model_CI.csv",

  row.names = FALSE
)


save(

  boot_df,

  boot_valid,

  paired_summary,

  model_CI,

  file =
    "Chinese_HF_STEP10C_MAGGIC_PAIRED_BOOTSTRAP.RData"
)


cat("\n============================================\n")
cat("STEP 10C 完成\n")
cat("============================================\n")
