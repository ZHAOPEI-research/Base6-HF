# ============================================================
# STEP 11D-B
# 1-year calibration comparison:
# Base6 vs MAGGIC
#
# Outcomes:
#   1. Observed 1-year mortality by Kaplan-Meier
#   2. Mean predicted 1-year risk
#   3. O/P ratio
#   4. Correct IPCW Brier score at 364 days
#   5. Calibration slope
#   6. Quintile calibration
#   7. Quintile calibration MAE
#
# NOTE:
#   GWTG-HF is NOT included in probability calibration because
#   gwtg_score is not itself a calibrated 1-year probability.
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

library(survival)

rm(list = ls())

set.seed(20260905)

# Public-repository path configuration
external_output_dir <- Sys.getenv("BASE6_EXTERNAL_OUTPUT_DIR")

if (external_output_dir == "") {
  stop(
    paste0(
      "BASE6_EXTERNAL_OUTPUT_DIR is not set. ",
      "Set it to an authorized local directory containing restricted ",
      "external-cohort analysis outputs."
    )
  )
}

if (!dir.exists(external_output_dir)) {
  stop(
    paste0(
      "BASE6_EXTERNAL_OUTPUT_DIR does not exist: ",
      external_output_dir
    )
  )
}


# ============================================================
# 1. Load latest STEP 11B common cohort
# ============================================================

load(file.path(external_output_dir, "Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData"))

if (!exists("comparison3")) {
  stop("ERROR: comparison3 not found.")
}

base_data <- comparison3


# ============================================================
# 2. Load MAGGIC predicted 1-year risk
# ============================================================

e_maggic <- new.env()

load(
  file.path(external_output_dir, "Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData"),
  envir = e_maggic
)

if (!exists("comparison_data", envir = e_maggic)) {
  stop(
    "ERROR: comparison_data not found in ",
    "Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData"
  )
}

maggic_data <- get(
  "comparison_data",
  envir = e_maggic
)


# ============================================================
# 3. Basic source audit
# ============================================================

cat("\n")
cat("============================================\n")
cat("STEP 11D-B: Base6 vs MAGGIC calibration\n")
cat("============================================\n\n")

cat("STEP 11B common cohort n =", nrow(base_data), "\n")
cat("MAGGIC comparison_data n =", nrow(maggic_data), "\n\n")


required_base <- c(
  "ID",
  "status",
  "time",
  "base6_marker"
)

required_maggic <- c(
  "ID",
  "risk_maggic_1y"
)


missing_base <- setdiff(
  required_base,
  names(base_data)
)

missing_maggic <- setdiff(
  required_maggic,
  names(maggic_data)
)


if (length(missing_base) > 0) {
  stop(
    paste(
      "Missing from base_data:",
      paste(missing_base, collapse = ", ")
    )
  )
}

if (length(missing_maggic) > 0) {
  stop(
    paste(
      "Missing from maggic_data:",
      paste(missing_maggic, collapse = ", ")
    )
  )
}


# ============================================================
# 4. Merge explicitly by ID
# ============================================================

cal_data <- merge(
  base_data[, required_base],
  maggic_data[, required_maggic],
  by = "ID",
  all = FALSE,
  sort = FALSE
)


# Restore original STEP 11B patient order
cal_data <- cal_data[
  match(
    base_data$ID,
    cal_data$ID
  ),
]

rownames(cal_data) <- NULL


cat("============================================\n")
cat("A. ID alignment audit\n")
cat("============================================\n\n")

cat("Merged n =", nrow(cal_data), "\n")

cat(
  "IDs identical to STEP 11B order =",
  identical(
    cal_data$ID,
    base_data$ID
  ),
  "\n"
)

cat(
  "Unique IDs =",
  length(unique(cal_data$ID)),
  "\n\n"
)


if (nrow(cal_data) != 346) {
  stop(
    paste0(
      "ERROR: merged calibration cohort n = ",
      nrow(cal_data),
      ", expected 346."
    )
  )
}

if (!identical(
  cal_data$ID,
  base_data$ID
)) {
  stop("ERROR: ID alignment failed.")
}


# ============================================================
# 5. Rename Base6 predicted risk clearly
# ============================================================

cal_data$risk_base6_1y <-
  cal_data$base6_marker


# ============================================================
# 6. Data completeness and validity audit
# ============================================================

analysis_cols <- c(
  "time",
  "status",
  "risk_base6_1y",
  "risk_maggic_1y"
)

complete_rows <- complete.cases(
  cal_data[, analysis_cols]
)


cat("============================================\n")
cat("B. Calibration data completeness\n")
cat("============================================\n\n")

cat(
  "Complete rows =",
  sum(complete_rows),
  "/",
  nrow(cal_data),
  "\n"
)

cat(
  "Deaths =",
  sum(cal_data$status == 1),
  "\n\n"
)


if (!all(complete_rows)) {
  stop("ERROR: missing calibration data.")
}


if (!all(cal_data$status %in% c(0, 1))) {
  stop("ERROR: status contains values other than 0/1.")
}

if (any(cal_data$time <= 0)) {
  stop("ERROR: follow-up time contains values <= 0.")
}


for (v in c(
  "risk_base6_1y",
  "risk_maggic_1y"
)) {

  if (
    any(
      cal_data[[v]] < 0 |
      cal_data[[v]] > 1
    )
  ) {

    stop(
      paste(
        "ERROR:",
        v,
        "contains values outside [0,1]."
      )
    )
  }
}


cat("Base6 risk summary:\n")
print(summary(cal_data$risk_base6_1y))
cat("\n")

cat("MAGGIC risk summary:\n")
print(summary(cal_data$risk_maggic_1y))
cat("\n")


# ============================================================
# 7. Kaplan-Meier observed mortality at 364 days
# ============================================================

calc_km_mortality <- function(
    time,
    status,
    t0 = 364
) {

  fit <- survival::survfit(
    survival::Surv(time, status) ~ 1
  )

  s <- summary(
    fit,
    times = t0,
    extend = TRUE
  )$surv

  1 - as.numeric(s)
}


observed_364 <- calc_km_mortality(
  cal_data$time,
  cal_data$status,
  364
)


mean_base6 <- mean(
  cal_data$risk_base6_1y
)

mean_maggic <- mean(
  cal_data$risk_maggic_1y
)


overall_calibration <- data.frame(

  model = c(
    "Base6",
    "MAGGIC"
  ),

  n = nrow(cal_data),

  deaths =
    sum(cal_data$status == 1),

  observed_KM_364 =
    observed_364,

  mean_predicted_364 = c(
    mean_base6,
    mean_maggic
  ),

  O_P = c(
    observed_364 / mean_base6,
    observed_364 / mean_maggic
  ),

  absolute_difference = c(
    abs(observed_364 - mean_base6),
    abs(observed_364 - mean_maggic)
  ),

  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("C. Overall 1-year calibration\n")
cat("============================================\n\n")

print(
  overall_calibration,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 8. Correct IPCW Brier score at 364 days
# ============================================================

calc_brier_364 <- function(
    time,
    status,
    pred,
    t0 = 364
) {

  n <- length(time)


  # ----------------------------------------------------------
  # Estimate censoring survival function G(t)
  #
  # status = 1 means death
  # therefore 1-status means censoring event
  # ----------------------------------------------------------

  censor_fit <- survival::survfit(
    survival::Surv(
      time,
      1 - status
    ) ~ 1
  )


  # ----------------------------------------------------------
  # Helper function for G(t)
  # ----------------------------------------------------------

  get_G <- function(t) {

    g <- summary(
      censor_fit,
      times = t,
      extend = TRUE
    )$surv

    as.numeric(g)
  }


  # ----------------------------------------------------------
  # Initialize all contributions as zero
  #
  # Important:
  # patients censored before/equal t0 contribute zero,
  # but remain in denominator n.
  # ----------------------------------------------------------

  score_terms <- rep(
    0,
    n
  )


  G_t0 <- get_G(t0)


  for (i in seq_len(n)) {

    # --------------------------------------------------------
    # Case 1:
    # Death observed before or at t0
    # --------------------------------------------------------

    if (
      time[i] <= t0 &&
      status[i] == 1
    ) {

      t_left <- max(
        time[i] - 1e-7,
        0
      )

      g <- get_G(
        t_left
      )

      if (
        is.finite(g) &&
        g > 0
      ) {

        score_terms[i] <-
          ((1 - pred[i])^2) / g

      } else {

        score_terms[i] <-
          NA_real_
      }


      # --------------------------------------------------------
      # Case 2:
      # Known alive/event-free beyond t0
      # --------------------------------------------------------

    } else if (
      time[i] > t0
    ) {

      if (
        is.finite(G_t0) &&
        G_t0 > 0
      ) {

        score_terms[i] <-
          (pred[i]^2) / G_t0

      } else {

        score_terms[i] <-
          NA_real_
      }


      # --------------------------------------------------------
      # Case 3:
      # Censored before or at t0
      #
      # contribution remains zero
      # --------------------------------------------------------

    } else {

      score_terms[i] <- 0
    }
  }


  # ----------------------------------------------------------
  # Crucial:
  # divide by FULL original sample size n
  # ----------------------------------------------------------

  sum(
    score_terms,
    na.rm = TRUE
  ) / n
}


# ------------------------------------------------------------
# Calculate model Brier scores
# ------------------------------------------------------------

brier_base6 <- calc_brier_364(
  cal_data$time,
  cal_data$status,
  cal_data$risk_base6_1y,
  t0 = 364
)

brier_maggic <- calc_brier_364(
  cal_data$time,
  cal_data$status,
  cal_data$risk_maggic_1y,
  t0 = 364
)


# ------------------------------------------------------------
# Null model Brier score
# ------------------------------------------------------------

null_risk <- rep(
  observed_364,
  nrow(cal_data)
)

brier_null <- calc_brier_364(
  cal_data$time,
  cal_data$status,
  null_risk,
  t0 = 364
)


brier_summary <- data.frame(

  model = c(
    "Base6",
    "MAGGIC",
    "Null"
  ),

  Brier_364 = c(
    brier_base6,
    brier_maggic,
    brier_null
  ),

  Scaled_Brier = c(
    1 - brier_base6 / brier_null,
    1 - brier_maggic / brier_null,
    0
  ),

  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("D. Corrected Brier score at 364 days\n")
cat("============================================\n\n")

print(
  brier_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 9. Calibration slope
# ============================================================

safe_logit <- function(p) {

  eps <- 1e-6

  p <- pmin(
    pmax(
      p,
      eps
    ),
    1 - eps
  )

  qlogis(p)
}


cal_data$logit_base6 <-
  safe_logit(
    cal_data$risk_base6_1y
  )

cal_data$logit_maggic <-
  safe_logit(
    cal_data$risk_maggic_1y
  )


fit_slope_base6 <- survival::coxph(
  survival::Surv(time, status) ~ logit_base6,
  data = cal_data
)

fit_slope_maggic <- survival::coxph(
  survival::Surv(time, status) ~ logit_maggic,
  data = cal_data
)


extract_slope <- function(fit) {

  s <- summary(fit)

  est <- s$coefficients[
    1,
    "coef"
  ]

  se <- s$coefficients[
    1,
    "se(coef)"
  ]

  c(
    slope = est,
    low = est - 1.96 * se,
    high = est + 1.96 * se
  )
}


slope_base6 <- extract_slope(
  fit_slope_base6
)

slope_maggic <- extract_slope(
  fit_slope_maggic
)


slope_summary <- data.frame(

  model = c(
    "Base6",
    "MAGGIC"
  ),

  calibration_slope = c(
    slope_base6["slope"],
    slope_maggic["slope"]
  ),

  slope_low = c(
    slope_base6["low"],
    slope_maggic["low"]
  ),

  slope_high = c(
    slope_base6["high"],
    slope_maggic["high"]
  ),

  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("E. Calibration slope\n")
cat("============================================\n\n")

print(
  slope_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 10. Quintile calibration
# ============================================================

make_cal_groups <- function(
    data,
    risk_col,
    model_name
) {

  risk <- data[[risk_col]]

  ranks <- rank(
    risk,
    ties.method = "first"
  )

  q <- cut(
    ranks,
    breaks = quantile(
      ranks,
      probs = seq(
        0,
        1,
        0.2
      )
    ),
    include.lowest = TRUE,
    labels = paste0(
      "Q",
      1:5
    )
  )


  out <- lapply(
    levels(q),
    function(g) {

      idx <- which(
        q == g
      )

      d <- data[
        idx,
        ,
        drop = FALSE
      ]


      obs <- calc_km_mortality(
        d$time,
        d$status,
        364
      )


      pred <- mean(
        d[[risk_col]]
      )


      data.frame(

        model =
          model_name,

        group =
          g,

        n =
          nrow(d),

        deaths =
          sum(
            d$status == 1
          ),

        mean_predicted =
          pred,

        observed_KM =
          obs,

        absolute_error =
          abs(
            obs -
              pred
          ),

        stringsAsFactors = FALSE
      )
    }
  )


  do.call(
    rbind,
    out
  )
}


cal_base6_quintile <-
  make_cal_groups(
    cal_data,
    "risk_base6_1y",
    "Base6"
  )

cal_maggic_quintile <-
  make_cal_groups(
    cal_data,
    "risk_maggic_1y",
    "MAGGIC"
  )


calibration_quintiles <- rbind(
  cal_base6_quintile,
  cal_maggic_quintile
)


cat("============================================\n")
cat("F. Quintile calibration\n")
cat("============================================\n\n")

print(
  calibration_quintiles,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 11. Quintile calibration MAE
# ============================================================

cal_mae <- aggregate(
  absolute_error ~ model,
  data = calibration_quintiles,
  FUN = mean
)

names(cal_mae)[2] <-
  "calibration_MAE_quintile"


cat("============================================\n")
cat("G. Quintile calibration MAE\n")
cat("============================================\n\n")

print(
  cal_mae,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 12. Final combined summary
# ============================================================

final_summary <- merge(
  overall_calibration,
  brier_summary[
    brier_summary$model != "Null",
    c(
      "model",
      "Brier_364",
      "Scaled_Brier"
    )
  ],
  by = "model"
)

final_summary <- merge(
  final_summary,
  slope_summary,
  by = "model"
)

final_summary <- merge(
  final_summary,
  cal_mae,
  by = "model"
)


cat("============================================\n")
cat("H. Final calibration summary\n")
cat("============================================\n\n")

print(
  final_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 13. Save CSV outputs
# ============================================================

write.csv(
  final_summary,
  file.path(external_output_dir, "Chinese_GWTG_STEP11D_CALIBRATION_SUMMARY.csv"),
  row.names = FALSE
)

write.csv(
  calibration_quintiles,
  file.path(external_output_dir, "Chinese_GWTG_STEP11D_CALIBRATION_QUINTILES.csv"),
  row.names = FALSE
)

write.csv(
  brier_summary,
  file.path(external_output_dir, "Chinese_GWTG_STEP11D_BRIER.csv"),
  row.names = FALSE
)

write.csv(
  slope_summary,
  file.path(external_output_dir, "Chinese_GWTG_STEP11D_CALIBRATION_SLOPE.csv"),
  row.names = FALSE
)


# ============================================================
# 14. Save RData
# ============================================================

save(
  cal_data,
  overall_calibration,
  brier_summary,
  slope_summary,
  calibration_quintiles,
  cal_mae,
  final_summary,
  file = file.path(
    external_output_dir,
    "Chinese_GWTG_STEP11D_CALIBRATION.RData"
  )
)


# ============================================================
# 15. Completion message
# ============================================================

cat("============================================\n")
cat("STEP 11D-B completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11D_CALIBRATION_SUMMARY.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11D_CALIBRATION_QUINTILES.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11D_BRIER.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11D_CALIBRATION_SLOPE.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11D_CALIBRATION.RData\n"
)

cat("\n")
