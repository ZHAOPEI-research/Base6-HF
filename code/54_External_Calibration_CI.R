rm(list = ls())
gc()

library(survival)
library(dplyr)

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

# ============================================================
# 1. Load external calibration data
# ============================================================

e <- new.env()

load(
  "Chinese_GWTG_STEP11F_RECALIBRATION.RData",
  envir = e
)

cal_data <- e$cal_data

cat("============================================\n")
cat("STEP 54: External calibration 95% CI\n")
cat("============================================\n")

cat("n =", nrow(cal_data), "\n")
cat("Columns:\n")
print(names(cal_data))


# ============================================================
# 2. Helper to detect columns
# ============================================================

pick_col <- function(dat, candidates, label) {

  hit <- candidates[
    candidates %in% names(dat)
  ]

  if (length(hit) == 0) {

    stop(
      paste0(
        "Could not identify column for ",
        label,
        ". Available columns:\n",
        paste(
          names(dat),
          collapse = ", "
        )
      )
    )
  }

  hit[1]
}


time_col <- pick_col(
  cal_data,
  c(
    "time",
    "followup_days",
    "followup_day",
    "followup"
  ),
  "follow-up time"
)

status_col <- pick_col(
  cal_data,
  c(
    "status",
    "death"
  ),
  "death status"
)

base6_col <- pick_col(
  cal_data,
  c(
    "risk_base6_1y",
    "risk_base6",
    "base6_risk",
    "Base6_risk",
    "base6_pred",
    "base6_marker",
    "pred_base6",
    "pred_base6_original"
  ),
  "Base6 predicted risk"
)

maggic_col <- pick_col(
  cal_data,
  c(
    "risk_maggic_1y",
    "risk_maggic",
    "maggic_risk",
    "MAGGIC_risk",
    "maggic_pred",
    "pred_maggic",
    "MAGGIC_pred"
  ),
  "MAGGIC predicted risk"
)

cat("\nDetected columns:\n")
cat("time   =", time_col, "\n")
cat("status =", status_col, "\n")
cat("Base6  =", base6_col, "\n")
cat("MAGGIC =", maggic_col, "\n")


# ============================================================
# 3. Standardize analysis dataset
# ============================================================

dat <- data.frame(
  time =
    cal_data[[time_col]],

  status =
    cal_data[[status_col]],

  risk_base6 =
    cal_data[[base6_col]],

  risk_maggic =
    cal_data[[maggic_col]]
)

stopifnot(
  nrow(dat) == 346
)

cat("\nDeaths =", sum(dat$status), "\n")

if (sum(dat$status) != 47) {
  stop("Death count is not 47.")
}


# ============================================================
# 4. KM risk at 364 days
# ============================================================

km_risk <- function(time, status, t = 364) {

  fit <- survfit(
    Surv(time, status) ~ 1
  )

  s <- summary(
    fit,
    times = t,
    extend = TRUE
  )$surv

  1 - as.numeric(s)
}


# ============================================================
# 5. IPCW Brier score
# ============================================================

ipcw_brier <- function(
    time,
    status,
    pred_risk,
    t = 364
) {

  censor_status <- 1 - status

  km_cens <- survfit(
    Surv(time, censor_status) ~ 1
  )

  G_at <- function(s) {

    out <- summary(
      km_cens,
      times = s,
      extend = TRUE
    )$surv

    as.numeric(out)
  }

  G_t <- G_at(t)

  if (
    length(G_t) == 0 ||
    !is.finite(G_t) ||
    G_t <= 0
  ) {
    return(NA_real_)
  }

  n <- length(time)

  w <- numeric(n)
  y <- numeric(n)

  for (i in seq_len(n)) {

    if (
      time[i] <= t &&
      status[i] == 1
    ) {

      s_i <- max(
        time[i] - 1e-8,
        0
      )

      G_i <- G_at(s_i)

      if (
        length(G_i) > 0 &&
        is.finite(G_i) &&
        G_i > 0
      ) {
        w[i] <- 1 / G_i
        y[i] <- 1
      }

    } else if (
      time[i] > t
    ) {

      w[i] <- 1 / G_t
      y[i] <- 0

    } else {

      w[i] <- 0
      y[i] <- 0
    }
  }

  mean(
    w *
      (y - pred_risk)^2,
    na.rm = TRUE
  )
}


# ============================================================
# 6. Calibration metrics
# ============================================================

get_metrics <- function(
    d,
    risk_col,
    t = 364
) {

  obs <- km_risk(
    d$time,
    d$status,
    t
  )

  pred <- mean(
    d[[risk_col]],
    na.rm = TRUE
  )

  op <- obs / pred

  abs_diff <- abs(
    obs - pred
  )

  brier <- ipcw_brier(
    time = d$time,
    status = d$status,
    pred_risk = d[[risk_col]],
    t = t
  )

  # Null prediction = observed KM risk
  null_pred <- rep(
    obs,
    nrow(d)
  )

  null_brier <- ipcw_brier(
    time = d$time,
    status = d$status,
    pred_risk = null_pred,
    t = t
  )

  scaled_brier <-
    1 -
    brier / null_brier

  c(
    observed_risk = obs,
    mean_predicted = pred,
    OP_ratio = op,
    absolute_difference = abs_diff,
    Brier = brier,
    Scaled_Brier = scaled_brier
  )
}


# ============================================================
# 7. Apparent point estimates
# ============================================================

point_base6 <- get_metrics(
  dat,
  "risk_base6"
)

point_maggic <- get_metrics(
  dat,
  "risk_maggic"
)

cat("\n============================================\n")
cat("A. Apparent calibration metrics\n")
cat("============================================\n")

print(
  rbind(
    Base6 = point_base6,
    MAGGIC = point_maggic
  )
)


# ============================================================
# 8. Audit against frozen manuscript values
# ============================================================

expected <- data.frame(
  model = c(
    "Base6",
    "MAGGIC"
  ),

  observed_risk = c(
    0.1475792,
    0.1475792
  ),

  mean_predicted = c(
    0.08554724,
    0.19616474
  ),

  OP_ratio = c(
    1.7251192,
    0.7523227
  ),

  absolute_difference = c(
    0.062031947,
    0.048585552
  ),

  Brier = c(
    0.1236033,
    0.1239257
  ),

  Scaled_Brier = c(
    0.01741787,
    0.01485505
  )
)

calculated <- data.frame(
  model = c(
    "Base6",
    "MAGGIC"
  ),

  rbind(
    point_base6,
    point_maggic
  ),
  row.names = NULL
)

audit <- calculated %>%
  left_join(
    expected,
    by = "model",
    suffix = c(
      "_new",
      "_expected"
    )
  ) %>%
  mutate(
    diff_observed =
      observed_risk_new -
      observed_risk_expected,

    diff_mean_pred =
      mean_predicted_new -
      mean_predicted_expected,

    diff_OP =
      OP_ratio_new -
      OP_ratio_expected,

    diff_absdiff =
      absolute_difference_new -
      absolute_difference_expected,

    diff_Brier =
      Brier_new -
      Brier_expected,

    diff_scaled =
      Scaled_Brier_new -
      Scaled_Brier_expected
  )

cat("\n============================================\n")
cat("B. Point-estimate audit\n")
cat("============================================\n")

print(
  audit %>%
    select(
      model,
      starts_with("diff_")
    ),
  row.names = FALSE
)


# ============================================================
# 9. Bootstrap
# ============================================================

B <- 2000

set.seed(20261004)

boot_base6 <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(point_base6)
)

boot_maggic <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(point_maggic)
)

colnames(boot_base6) <-
  names(point_base6)

colnames(boot_maggic) <-
  names(point_maggic)

n <- nrow(dat)

for (b in seq_len(B)) {

  ix <- sample.int(
    n,
    size = n,
    replace = TRUE
  )

  db <- dat[
    ix,
    ,
    drop = FALSE
  ]

  boot_base6[b, ] <-
    tryCatch(
      get_metrics(
        db,
        "risk_base6"
      ),
      error =
        function(e)
          rep(
            NA_real_,
            length(point_base6)
          )
    )

  boot_maggic[b, ] <-
    tryCatch(
      get_metrics(
        db,
        "risk_maggic"
      ),
      error =
        function(e)
          rep(
            NA_real_,
            length(point_maggic)
          )
    )

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


# ============================================================
# 10. Percentile CI
# ============================================================

get_ci <- function(
    boot_mat,
    point,
    model_name
) {

  data.frame(
    Model =
      model_name,

    Metric =
      names(point),

    Estimate =
      as.numeric(point),

    CI_low =
      apply(
        boot_mat,
        2,
        quantile,
        probs = 0.025,
        na.rm = TRUE
      ),

    CI_high =
      apply(
        boot_mat,
        2,
        quantile,
        probs = 0.975,
        na.rm = TRUE
      ),

    Valid_B =
      apply(
        boot_mat,
        2,
        function(x)
          sum(
            is.finite(x)
          )
      ),

    stringsAsFactors = FALSE
  )
}


ci_base6 <- get_ci(
  boot_base6,
  point_base6,
  "Base6 original"
)

ci_maggic <- get_ci(
  boot_maggic,
  point_maggic,
  "MAGGIC"
)

calibration_CI <- bind_rows(
  ci_base6,
  ci_maggic
)


cat("\n============================================\n")
cat("C. External calibration 95% CI\n")
cat("============================================\n")

print(
  calibration_CI,
  row.names = FALSE
)


# ============================================================
# 11. Manuscript-ready summary
# ============================================================

manuscript_table <- calibration_CI %>%
  mutate(
    Result =
      sprintf(
        "%.3f (%.3f–%.3f)",
        Estimate,
        CI_low,
        CI_high
      )
  ) %>%
  select(
    Model,
    Metric,
    Result,
    Valid_B
  )

cat("\n============================================\n")
cat("D. Manuscript-ready values\n")
cat("============================================\n")

print(
  manuscript_table,
  row.names = FALSE
)


# ============================================================
# 12. Save
# ============================================================

write.csv(
  calibration_CI,
  "External_Calibration_95CI_2000boot.csv",
  row.names = FALSE
)

write.csv(
  audit,
  "External_Calibration_95CI_AUDIT.csv",
  row.names = FALSE
)

save(
  dat,
  point_base6,
  point_maggic,
  boot_base6,
  boot_maggic,
  calibration_CI,
  audit,
  file =
    "External_Calibration_95CI_2000BOOT.RData"
)

cat("\nSTEP 54 completed\n")
