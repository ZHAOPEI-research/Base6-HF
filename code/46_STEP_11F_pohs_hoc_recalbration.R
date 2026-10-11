# ============================================================
# STEP 11F
# Post-hoc recalibration of Base6 in external validation cohort
#
# Models compared:
#   1. Base6 original
#   2. Base6 recalibrated: intercept-only
#   3. Base6 recalibrated: intercept + slope
#   4. MAGGIC
#
# Recalibration method:
#
#   z = log[-log(1 - predicted risk)]
#
# Intercept-only:
#   cloglog(p_new) = alpha + z
#
# Intercept + slope:
#   cloglog(p_new) = alpha + beta * z
#
# For survival data, alpha is obtained through recalibration
# of the baseline cumulative hazard using the Cox/Breslow
# framework.
#
# Time horizon:
#   364 days
#
# Input:
#   Chinese_GWTG_STEP11D_CALIBRATION.RData
#
# IMPORTANT:
#   Recalibration is estimated and evaluated in the SAME
#   external cohort. Results are post-hoc/apparent and should
#   not replace the original external validation performance.
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

rm(list = ls())

library(survival)

if (!requireNamespace("timeROC", quietly = TRUE)) {
  stop(
    "Package 'timeROC' is required. ",
    "Install it using install.packages('timeROC')."
  )
}

library(timeROC)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required. ",
    "Install it using install.packages('ggplot2')."
  )
}

library(ggplot2)



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

t0 <- 364


# ============================================================
# 1. Load STEP 11D calibrated comparison dataset
# ============================================================

load(
  "Chinese_GWTG_STEP11D_CALIBRATION.RData"
)


if (!exists("cal_data")) {
  stop(
    "ERROR: cal_data not found in ",
    "Chinese_GWTG_STEP11D_CALIBRATION.RData"
  )
}


cat("\n")
cat("============================================\n")
cat("STEP 11F: Base6 post-hoc recalibration\n")
cat("============================================\n\n")


# ============================================================
# 2. Input audit
# ============================================================

required_cols <- c(
  "ID",
  "time",
  "status",
  "risk_base6_1y",
  "risk_maggic_1y"
)


missing_cols <- setdiff(
  required_cols,
  names(cal_data)
)


if (length(missing_cols) > 0) {

  stop(
    paste0(
      "Missing required columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}


if (
  any(
    !complete.cases(
      cal_data[, required_cols]
    )
  )
) {

  stop(
    "ERROR: missing values found in required variables."
  )
}


if (
  !all(
    cal_data$status %in% c(0, 1)
  )
) {

  stop(
    "ERROR: status must be coded as 0/1."
  )
}


if (
  any(
    cal_data$time <= 0
  )
) {

  stop(
    "ERROR: follow-up time must be > 0."
  )
}


if (
  any(
    cal_data$risk_base6_1y <= 0 |
    cal_data$risk_base6_1y >= 1
  )
) {

  stop(
    "ERROR: Base6 probabilities must lie strictly ",
    "between 0 and 1 for cloglog transformation."
  )
}


cat("============================================\n")
cat("A. Input cohort audit\n")
cat("============================================\n\n")

cat(
  "N =",
  nrow(cal_data),
  "\n"
)

cat(
  "Deaths =",
  sum(cal_data$status == 1),
  "\n"
)

cat(
  "Censored/alive =",
  sum(cal_data$status == 0),
  "\n"
)

cat(
  "Follow-up range =",
  min(cal_data$time),
  "to",
  max(cal_data$time),
  "days\n\n"
)


cat("Original Base6 risk:\n")
print(
  summary(
    cal_data$risk_base6_1y
  )
)

cat("\n")

cat("MAGGIC risk:\n")
print(
  summary(
    cal_data$risk_maggic_1y
  )
)

cat("\n")


# ============================================================
# 3. Complementary log-log transformation
# ============================================================

safe_cloglog_risk <- function(p) {

  eps <- 1e-8

  p <- pmin(
    pmax(
      p,
      eps
    ),
    1 - eps
  )

  log(
    -log(
      1 - p
    )
  )
}


cal_data$z_base6 <- safe_cloglog_risk(
  cal_data$risk_base6_1y
)


cat("============================================\n")
cat("B. Base6 cloglog predictor\n")
cat("============================================\n\n")

print(
  summary(
    cal_data$z_base6
  )
)

cat("\n")


# ============================================================
# 4. Breslow cumulative baseline hazard at t0
#
# eta is the fixed linear predictor.
# ============================================================

breslow_H0_t0 <- function(
    time,
    status,
    eta,
    t0 = 364
) {

  event_times <- sort(
    unique(
      time[
        status == 1 &
          time <= t0
      ]
    )
  )


  if (
    length(event_times) == 0
  ) {

    return(
      0
    )
  }


  H0 <- 0


  for (tt in event_times) {

    d_j <- sum(
      status == 1 &
        time == tt
    )


    risk_set <- (
      time >= tt
    )


    denominator <- sum(
      exp(
        eta[risk_set]
      )
    )


    if (
      !is.finite(denominator) ||
      denominator <= 0
    ) {

      stop(
        paste0(
          "Invalid Breslow denominator at time ",
          tt
        )
      )
    }


    H0 <- H0 +
      d_j / denominator
  }


  as.numeric(
    H0
  )
}


# ============================================================
# 5. INTERCEPT-ONLY recalibration
#
# slope fixed at 1
#
# Original predicted cumulative hazard:
#
#   H_pred = -log(1-p)
#
# Since:
#
#   z = log(H_pred)
#
# recalibration becomes:
#
#   H_new(t0) = H0_ext(t0) * exp(z)
#
# ============================================================

eta_intercept_only <-
  cal_data$z_base6


H0_intercept_only <-
  breslow_H0_t0(
    time =
      cal_data$time,

    status =
      cal_data$status,

    eta =
      eta_intercept_only,

    t0 =
      t0
  )


alpha_intercept_only <-
  log(
    H0_intercept_only
  )


cal_data$risk_base6_recal_intercept <-
  1 - exp(
    -H0_intercept_only *
      exp(
        eta_intercept_only
      )
  )


cat("============================================\n")
cat("C. Intercept-only recalibration\n")
cat("============================================\n\n")

cat(
  "Baseline cumulative hazard at 364 days =",
  H0_intercept_only,
  "\n"
)

cat(
  "Calibration intercept (cloglog scale) =",
  alpha_intercept_only,
  "\n"
)

cat(
  "Slope fixed = 1\n\n"
)

cat(
  "Recalibrated risk summary:\n"
)

print(
  summary(
    cal_data$risk_base6_recal_intercept
  )
)

cat("\n")


# ============================================================
# 6. INTERCEPT + SLOPE recalibration
#
# Cox model:
#
# Surv(time,status) ~ z_base6
#
# beta = recalibration slope
#
# baseline hazard is subsequently re-estimated with
# Breslow estimator.
# ============================================================

fit_full_recal <- survival::coxph(
  survival::Surv(
    time,
    status
  ) ~ z_base6,
  data =
    cal_data,
  ties =
    "breslow",
  x =
    TRUE,
  y =
    TRUE
)


fit_full_summary <- summary(
  fit_full_recal
)


beta_full <-
  as.numeric(
    coef(
      fit_full_recal
    )["z_base6"]
  )


beta_full_se <-
  fit_full_summary$coefficients[
    "z_base6",
    "se(coef)"
  ]


beta_full_low <-
  beta_full -
  1.96 *
  beta_full_se


beta_full_high <-
  beta_full +
  1.96 *
  beta_full_se


eta_full <-
  beta_full *
  cal_data$z_base6


H0_full <-
  breslow_H0_t0(
    time =
      cal_data$time,

    status =
      cal_data$status,

    eta =
      eta_full,

    t0 =
      t0
  )


alpha_full <-
  log(
    H0_full
  )


cal_data$risk_base6_recal_full <-
  1 - exp(
    -H0_full *
      exp(
        eta_full
      )
  )


cat("============================================\n")
cat("D. Intercept + slope recalibration\n")
cat("============================================\n\n")

cat(
  "Calibration slope =",
  beta_full,
  "\n"
)

cat(
  "Slope 95% CI =",
  beta_full_low,
  "to",
  beta_full_high,
  "\n"
)

cat(
  "Baseline cumulative hazard at 364 days =",
  H0_full,
  "\n"
)

cat(
  "Calibration intercept (cloglog scale) =",
  alpha_full,
  "\n\n"
)

cat(
  "Fully recalibrated risk summary:\n"
)

print(
  summary(
    cal_data$risk_base6_recal_full
  )
)

cat("\n")


# ============================================================
# 7. Kaplan-Meier mortality at t0
# ============================================================

calc_km_mortality <- function(
    time,
    status,
    t0 = 364
) {

  fit <- survival::survfit(
    survival::Surv(
      time,
      status
    ) ~ 1
  )


  s <- summary(
    fit,
    times =
      t0,
    extend =
      TRUE
  )$surv


  1 -
    as.numeric(
      s
    )
}


observed_364 <-
  calc_km_mortality(
    cal_data$time,
    cal_data$status,
    t0
  )


cat("============================================\n")
cat("E. Observed 364-day mortality\n")
cat("============================================\n\n")

cat(
  "Observed KM mortality =",
  observed_364,
  "\n\n"
)


# ============================================================
# 8. Correct IPCW Brier score
# ============================================================

calc_brier_364 <- function(
    time,
    status,
    pred,
    t0 = 364
) {

  n <- length(
    time
  )


  censor_fit <- survival::survfit(
    survival::Surv(
      time,
      1 - status
    ) ~ 1
  )


  get_G <- function(t) {

    g <- summary(
      censor_fit,
      times =
        t,
      extend =
        TRUE
    )$surv

    as.numeric(
      g
    )
  }


  score_terms <- rep(
    0,
    n
  )


  G_t0 <- get_G(
    t0
  )


  for (i in seq_len(n)) {

    # --------------------------------------------------------
    # Observed death before/equal t0
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
      # Known event-free at t0
      # --------------------------------------------------------

    } else if (
      time[i] > t0
    ) {

      if (
        is.finite(G_t0) &&
        G_t0 > 0
      ) {

        score_terms[i] <-
          (pred[i]^2) /
          G_t0

      } else {

        score_terms[i] <-
          NA_real_
      }


      # --------------------------------------------------------
      # Censored before/equal t0
      # contribution = zero
      # --------------------------------------------------------

    } else {

      score_terms[i] <- 0
    }
  }


  sum(
    score_terms,
    na.rm = TRUE
  ) / n
}


# ============================================================
# 9. Overall calibration metrics
# ============================================================

risk_cols <- c(
  "risk_base6_1y",
  "risk_base6_recal_intercept",
  "risk_base6_recal_full",
  "risk_maggic_1y"
)

model_names <- c(
  "Base6 original",
  "Base6 recalibrated: intercept only",
  "Base6 recalibrated: intercept + slope",
  "MAGGIC"
)

overall_list <- vector(
  "list",
  length(risk_cols)
)

for (i in seq_along(risk_cols)) {

  pred <- cal_data[[risk_cols[i]]]

  mean_pred <- mean(pred)

  brier <- calc_brier_364(
    time = cal_data$time,
    status = cal_data$status,
    pred = pred,
    t0 = t0
  )

  overall_list[[i]] <- data.frame(
    model = model_names[i],
    n = nrow(cal_data),
    deaths = sum(cal_data$status == 1),
    observed_KM_364 = observed_364,
    mean_predicted_364 = mean_pred,
    O_P = observed_364 / mean_pred,
    absolute_difference = abs(observed_364 - mean_pred),
    Brier_364 = brier,
    stringsAsFactors = FALSE
  )
}

overall_summary <- do.call(
  rbind,
  overall_list
)

rownames(overall_summary) <- NULL


# ============================================================
# 10. Null Brier score and scaled Brier
# ============================================================

null_pred <- rep(
  observed_364,
  nrow(
    cal_data
  )
)


brier_null <-
  calc_brier_364(
    time =
      cal_data$time,

    status =
      cal_data$status,

    pred =
      null_pred,

    t0 =
      t0
  )


overall_summary$Scaled_Brier <-
  1 -
  overall_summary$Brier_364 /
  brier_null


cat("============================================\n")
cat("F. Overall recalibration performance\n")
cat("============================================\n\n")

print(
  overall_summary,
  row.names =
    FALSE
)

cat("\n")

cat(
  "Null Brier =",
  brier_null,
  "\n\n"
)


# ============================================================
# 11. Quintile calibration function
# ============================================================

make_cal_groups <- function(
    data,
    risk_col,
    model_name,
    t0 = 364
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
      probs = seq(0, 1, 0.2)
    ),
    include.lowest = TRUE,
    labels = paste0("Q", 1:5)
  )

  out <- lapply(
    levels(q),
    function(g) {

      idx <- which(q == g)

      d <- data[
        idx,
        ,
        drop = FALSE
      ]

      obs <- calc_km_mortality(
        d$time,
        d$status,
        t0
      )

      pred <- mean(
        d[[risk_col]]
      )

      data.frame(
        model = model_name,
        group = g,
        n = nrow(d),
        deaths = sum(d$status == 1),
        mean_predicted = pred,
        observed_KM = obs,
        absolute_error = abs(obs - pred),
        stringsAsFactors = FALSE
      )
    }
  )

  do.call(
    rbind,
    out
  )
}

# ============================================================
# 12. Quintile calibration for all models
# ============================================================

q_original <-
  make_cal_groups(
    cal_data,
    "risk_base6_1y",
    "Base6 original",
    t0
  )


q_intercept <-
  make_cal_groups(
    cal_data,
    "risk_base6_recal_intercept",
    "Base6 recalibrated: intercept only",
    t0
  )


q_full <-
  make_cal_groups(
    cal_data,
    "risk_base6_recal_full",
    "Base6 recalibrated: intercept + slope",
    t0
  )


q_maggic <-
  make_cal_groups(
    cal_data,
    "risk_maggic_1y",
    "MAGGIC",
    t0
  )


calibration_quintiles <- rbind(
  q_original,
  q_intercept,
  q_full,
  q_maggic
)


# ============================================================
# 13. Quintile calibration MAE
# ============================================================

mae_summary <- aggregate(
  absolute_error ~ model,
  data =
    calibration_quintiles,
  FUN =
    mean
)


names(
  mae_summary
)[2] <-
  "calibration_MAE_quintile"


cat("============================================\n")
cat("G. Quintile calibration MAE\n")
cat("============================================\n\n")

print(
  mae_summary,
  row.names =
    FALSE
)

cat("\n")


# ============================================================
# 14. Merge overall + MAE
# ============================================================

final_recal_summary <- merge(
  overall_summary,
  mae_summary,
  by =
    "model",
  sort =
    FALSE
)


# restore logical display order

desired_order <- c(
  "Base6 original",
  "Base6 recalibrated: intercept only",
  "Base6 recalibrated: intercept + slope",
  "MAGGIC"
)


final_recal_summary$model <- factor(
  final_recal_summary$model,
  levels =
    desired_order
)


final_recal_summary <-
  final_recal_summary[
    order(
      final_recal_summary$model
    ),
  ]


final_recal_summary$model <-
  as.character(
    final_recal_summary$model
  )


rownames(
  final_recal_summary
) <- NULL


cat("============================================\n")
cat("H. Final recalibration summary\n")
cat("============================================\n\n")

print(
  final_recal_summary,
  row.names =
    FALSE
)

cat("\n")


# ============================================================
# 15. Discrimination audit
#
# Recalibration should not materially alter ranking if the
# recalibration slope is positive.
# ============================================================

calc_auc_364 <- function(
    time,
    status,
    marker,
    t0 = 364
) {

  roc <- timeROC::timeROC(
    T =
      time,

    delta =
      status,

    marker =
      marker,

    cause =
      1,

    weighting =
      "marginal",

    times =
      t0,

    iid =
      FALSE
  )


  idx <- which.min(
    abs(
      roc$times -
        t0
    )
  )


  as.numeric(
    roc$AUC[
      idx
    ]
  )
}


calc_cindex <- function(
    time,
    status,
    marker
) {

  fit <- survival::concordance(
    survival::Surv(
      time,
      status
    ) ~ marker,
    reverse =
      TRUE
  )


  as.numeric(
    fit$concordance
  )
}


discrimination_summary <- data.frame(

  model = c(
    "Base6 original",
    "Base6 recalibrated: intercept only",
    "Base6 recalibrated: intercept + slope"
  ),

  AUC_364 = c(

    calc_auc_364(
      cal_data$time,
      cal_data$status,
      cal_data$risk_base6_1y,
      t0
    ),

    calc_auc_364(
      cal_data$time,
      cal_data$status,
      cal_data$risk_base6_recal_intercept,
      t0
    ),

    calc_auc_364(
      cal_data$time,
      cal_data$status,
      cal_data$risk_base6_recal_full,
      t0
    )
  ),

  Harrell_C = c(

    calc_cindex(
      cal_data$time,
      cal_data$status,
      cal_data$risk_base6_1y
    ),

    calc_cindex(
      cal_data$time,
      cal_data$status,
      cal_data$risk_base6_recal_intercept
    ),

    calc_cindex(
      cal_data$time,
      cal_data$status,
      cal_data$risk_base6_recal_full
    )
  ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("I. Discrimination invariance audit\n")
cat("============================================\n\n")

print(
  discrimination_summary,
  row.names =
    FALSE
)

cat("\n")


# ============================================================
# 16. Recalibration parameters
# ============================================================

recalibration_parameters <- data.frame(

  model = c(
    "Intercept-only",
    "Intercept + slope"
  ),

  intercept_cloglog = c(
    alpha_intercept_only,
    alpha_full
  ),

  slope = c(
    1,
    beta_full
  ),

  slope_CI_low = c(
    NA_real_,
    beta_full_low
  ),

  slope_CI_high = c(
    NA_real_,
    beta_full_high
  ),

  baseline_cumulative_hazard_364 = c(
    H0_intercept_only,
    H0_full
  ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("J. Recalibration parameters\n")
cat("============================================\n\n")

print(
  recalibration_parameters,
  row.names =
    FALSE
)

cat("\n")


# ============================================================
# 17. Prepare calibration plot data
#
# For clarity in the main figure:
#   Original Base6
#   Intercept-only recalibrated Base6
#   Full recalibrated Base6
#
# MAGGIC is retained in tables and earlier STEP 11E figure.
# ============================================================

plot_data <- rbind(
  q_original,
  q_intercept,
  q_full
)


plot_data$model <- factor(
  plot_data$model,
  levels = c(
    "Base6 original",
    "Base6 recalibrated: intercept only",
    "Base6 recalibrated: intercept + slope"
  )
)


# ============================================================
# 18. Common plotting range
# ============================================================

plot_max <- max(
  plot_data$mean_predicted,
  plot_data$observed_KM,
  na.rm =
    TRUE
)


plot_max <- ceiling(
  plot_max *
    10
) / 10


plot_max <- max(
  0.4,
  plot_max
)


# ============================================================
# 19. Recalibration figure
# ============================================================

p_recal <- ggplot(
  plot_data,
  aes(
    x =
      mean_predicted,

    y =
      observed_KM,

    group =
      model,

    shape =
      model,

    linetype =
      model
  )
) +

  geom_abline(
    intercept =
      0,

    slope =
      1,

    linetype =
      "dashed",

    linewidth =
      0.7
  ) +

  geom_line(
    linewidth =
      0.8
  ) +

  geom_point(
    size =
      3
  ) +

  scale_shape_manual(
    values = c(
      "Base6 original" = 16,
      "Base6 recalibrated: intercept only" = 17,
      "Base6 recalibrated: intercept + slope" = 15
    )
  ) +

  scale_linetype_manual(
    values = c(
      "Base6 original" = "solid",
      "Base6 recalibrated: intercept only" = "longdash",
      "Base6 recalibrated: intercept + slope" = "dotdash"
    )
  ) +

  scale_x_continuous(
    name =
      "Mean predicted 1-year mortality risk",

    limits =
      c(
        0,
        plot_max
      ),

    breaks =
      seq(
        0,
        plot_max,
        0.1
      ),

    expand =
      c(
        0,
        0
      )
  ) +

  scale_y_continuous(
    name =
      "Observed 1-year mortality (Kaplan–Meier)",

    limits =
      c(
        0,
        plot_max
      ),

    breaks =
      seq(
        0,
        plot_max,
        0.1
      ),

    expand =
      c(
        0,
        0
      )
  ) +

  labs(
    title =
      "Post-hoc recalibration of Base6",

    subtitle =
      "External validation cohort; 364-day mortality",

    shape =
      "Prediction",

    linetype =
      "Prediction",

    caption =
      paste0(
        "Recalibration parameters were estimated in the same ",
        "external-validation cohort; recalibrated performance ",
        "therefore represents apparent post-hoc performance."
      )
  ) +

  theme_classic(
    base_size =
      12
  ) +

  theme(

    plot.title =
      element_text(
        face =
          "bold",
        size =
          14
      ),

    legend.position =
      "top",

    legend.title =
      element_text(
        face =
          "bold"
      ),

    plot.caption =
      element_text(
        hjust =
          0,
        size =
          9
      ),

    aspect.ratio =
      1
  )


print(
  p_recal
)


# ============================================================
# 20. Save CSV outputs
# ============================================================

write.csv(
  recalibration_parameters,
  "Chinese_GWTG_STEP11F_RECALIBRATION_PARAMETERS.csv",
  row.names =
    FALSE
)


write.csv(
  final_recal_summary,
  "Chinese_GWTG_STEP11F_RECALIBRATION_SUMMARY.csv",
  row.names =
    FALSE
)


write.csv(
  calibration_quintiles,
  "Chinese_GWTG_STEP11F_RECALIBRATION_QUINTILES.csv",
  row.names =
    FALSE
)


write.csv(
  discrimination_summary,
  "Chinese_GWTG_STEP11F_DISCRIMINATION_INVARIANCE.csv",
  row.names =
    FALSE
)


# ============================================================
# 21. Save figures
# ============================================================

ggsave(
  filename =
    "Chinese_GWTG_STEP11F_RECALIBRATION.png",

  plot =
    p_recal,

  width =
    7.5,

  height =
    7.0,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


ggsave(
  filename =
    "Chinese_GWTG_STEP11F_RECALIBRATION.pdf",

  plot =
    p_recal,

  width =
    7.5,

  height =
    7.0,

  units =
    "in",

  device =
    cairo_pdf
)


# ============================================================
# 22. Save RData
# ============================================================

save(
  cal_data,
  fit_full_recal,
  recalibration_parameters,
  overall_summary,
  final_recal_summary,
  calibration_quintiles,
  mae_summary,
  discrimination_summary,
  p_recal,
  file =
    "Chinese_GWTG_STEP11F_RECALIBRATION.RData"
)


# ============================================================
# 23. Final audit
# ============================================================

cat("============================================\n")
cat("K. Final key comparison\n")
cat("============================================\n\n")


print(
  final_recal_summary[
    ,
    c(
      "model",
      "observed_KM_364",
      "mean_predicted_364",
      "O_P",
      "absolute_difference",
      "Brier_364",
      "Scaled_Brier",
      "calibration_MAE_quintile"
    )
  ],
  row.names =
    FALSE
)

cat("\n")


# ============================================================
# 24. Completion
# ============================================================

cat("============================================\n")
cat("STEP 11F completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11F_RECALIBRATION_PARAMETERS.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11F_RECALIBRATION_SUMMARY.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11F_RECALIBRATION_QUINTILES.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11F_DISCRIMINATION_INVARIANCE.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11F_RECALIBRATION.png\n"
)

cat(
  "6. Chinese_GWTG_STEP11F_RECALIBRATION.pdf\n"
)

cat(
  "7. Chinese_GWTG_STEP11F_RECALIBRATION.RData\n"
)

cat("\n")
cat(
  "IMPORTANT: recalibrated results are post-hoc/apparent ",
  "because recalibration was fitted in this validation cohort.\n"
)

cat("============================================\n")
