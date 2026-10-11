# ============================================================
# STEP 11I
# Bootstrap optimism correction for Base6 recalibration
#
# Recalibration strategies:
#   1. Intercept-only
#   2. Intercept + slope
#
# Bootstrap workflow:
#   - resample patients
#   - refit recalibration in bootstrap sample
#   - evaluate in bootstrap sample (training/apparent)
#   - evaluate same fitted recalibration in original cohort
#   - estimate optimism
#
# Metrics:
#   1. Brier score at 364 days
#   2. Quintile calibration MAE
#   3. O/P ratio
#   4. |log(O/P)| calibration error
#
# Input:
#   Chinese_GWTG_STEP11F_RECALIBRATION.RData
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

rm(list = ls())

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

set.seed(20260905)

t0 <- 364
B <- 1000


# ============================================================
# 1. Load STEP 11F
# ============================================================

load(
  "Chinese_GWTG_STEP11F_RECALIBRATION.RData"
)


if (!exists("cal_data")) {

  stop(
    "ERROR: cal_data not found in STEP 11F RData."
  )
}


cat("\n")
cat("============================================\n")
cat("STEP 11I: Bootstrap optimism correction\n")
cat("============================================\n\n")


# ============================================================
# 2. Input audit
# ============================================================

required_cols <- c(
  "ID",
  "time",
  "status",
  "risk_base6_1y"
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
    "ERROR: missing required data."
  )
}


cat("============================================\n")
cat("A. Cohort audit\n")
cat("============================================\n\n")

cat(
  "N =",
  nrow(cal_data),
  "\n"
)

cat(
  "Deaths =",
  sum(
    cal_data$status == 1
  ),
  "\n"
)

cat(
  "Bootstrap repetitions =",
  B,
  "\n\n"
)


# ============================================================
# 3. cloglog transformation
# ============================================================

safe_cloglog <- function(p) {

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


cal_data$z_base6 <-
  safe_cloglog(
    cal_data$risk_base6_1y
  )


# ============================================================
# 4. Breslow cumulative baseline hazard
# ============================================================

breslow_H0 <- function(
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
      NA_real_
    )
  }


  H0 <- 0


  for (tt in event_times) {

    d_j <- sum(
      status == 1 &
        time == tt
    )


    risk_set <-
      time >= tt


    denom <- sum(
      exp(
        eta[
          risk_set
        ]
      )
    )


    if (
      !is.finite(denom) ||
      denom <= 0
    ) {

      return(
        NA_real_
      )
    }


    H0 <-
      H0 +
      d_j / denom
  }


  as.numeric(
    H0
  )
}


# ============================================================
# 5. Fit intercept-only recalibration
# ============================================================

fit_intercept_only <- function(
    data,
    t0 = 364
) {

  eta <- data$z_base6


  H0 <- breslow_H0(
    time = data$time,
    status = data$status,
    eta = eta,
    t0 = t0
  )


  if (
    !is.finite(H0) ||
    H0 <= 0
  ) {

    return(
      NULL
    )
  }


  list(
    H0 = H0,
    beta = 1
  )
}


# ============================================================
# 6. Fit full recalibration
# ============================================================

fit_full_recalibration <- function(
    data,
    t0 = 364
) {

  fit <- tryCatch(

    survival::coxph(
      survival::Surv(
        time,
        status
      ) ~ z_base6,
      data = data,
      ties = "breslow"
    ),

    error = function(e) NULL
  )


  if (is.null(fit)) {

    return(
      NULL
    )
  }


  beta <- as.numeric(
    coef(fit)[
      "z_base6"
    ]
  )


  if (
    !is.finite(beta)
  ) {

    return(
      NULL
    )
  }


  eta <-
    beta *
    data$z_base6


  H0 <- breslow_H0(
    time = data$time,
    status = data$status,
    eta = eta,
    t0 = t0
  )


  if (
    !is.finite(H0) ||
    H0 <= 0
  ) {

    return(
      NULL
    )
  }


  list(
    H0 = H0,
    beta = beta
  )
}


# ============================================================
# 7. Apply recalibration parameters
# ============================================================

predict_recalibrated <- function(
    data,
    fit_object
) {

  eta <-
    fit_object$beta *
    data$z_base6


  pred <-
    1 -
    exp(
      -fit_object$H0 *
        exp(
          eta
        )
    )


  pmin(
    pmax(
      pred,
      0
    ),
    1
  )
}


# ============================================================
# 8. KM observed mortality
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
    times = t0,
    extend = TRUE
  )$surv


  1 -
    as.numeric(
      s
    )
}


# ============================================================
# 9. Correct IPCW Brier score
# ============================================================

calc_brier <- function(
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
      times = t,
      extend = TRUE
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


  if (
    !is.finite(G_t0) ||
    G_t0 <= 0
  ) {

    return(
      NA_real_
    )
  }


  for (i in seq_len(n)) {

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
          ((1 - pred[i])^2) /
          g

      } else {

        return(
          NA_real_
        )
      }


    } else if (
      time[i] > t0
    ) {

      score_terms[i] <-
        pred[i]^2 /
        G_t0


    } else {

      score_terms[i] <- 0
    }
  }


  sum(
    score_terms
  ) / n
}


# ============================================================
# 10. Quintile calibration MAE
# ============================================================

calc_quintile_mae <- function(
    data,
    pred,
    t0 = 364
) {

  ranks <- rank(
    pred,
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
    labels = FALSE
  )


  errors <- rep(
    NA_real_,
    5
  )


  for (g in 1:5) {

    idx <- which(
      q == g
    )


    if (
      length(idx) == 0
    ) {

      next
    }


    obs <- calc_km_mortality(
      data$time[idx],
      data$status[idx],
      t0
    )


    mean_pred <- mean(
      pred[idx]
    )


    errors[g] <-
      abs(
        obs -
          mean_pred
      )
  }


  mean(
    errors,
    na.rm = TRUE
  )
}


# ============================================================
# 11. Evaluation function
# ============================================================

evaluate_predictions <- function(
    data,
    pred,
    t0 = 364
) {

  obs <- calc_km_mortality(
    data$time,
    data$status,
    t0
  )


  mean_pred <- mean(
    pred
  )


  O_P <-
    obs /
    mean_pred


  brier <- calc_brier(
    time = data$time,
    status = data$status,
    pred = pred,
    t0 = t0
  )


  mae <- calc_quintile_mae(
    data = data,
    pred = pred,
    t0 = t0
  )


  c(

    observed =
      obs,

    mean_predicted =
      mean_pred,

    O_P =
      O_P,

    abs_log_OP =
      abs(
        log(
          O_P
        )
      ),

    Brier =
      brier,

    MAE =
      mae
  )
}


# ============================================================
# 12. Apparent performance in original cohort
# ============================================================

fit_original_intercept <-
  fit_intercept_only(
    cal_data,
    t0
  )


fit_original_full <-
  fit_full_recalibration(
    cal_data,
    t0
  )


pred_original_intercept <-
  predict_recalibrated(
    cal_data,
    fit_original_intercept
  )


pred_original_full <-
  predict_recalibrated(
    cal_data,
    fit_original_full
  )


apparent_intercept <-
  evaluate_predictions(
    cal_data,
    pred_original_intercept,
    t0
  )


apparent_full <-
  evaluate_predictions(
    cal_data,
    pred_original_full,
    t0
  )


cat("============================================\n")
cat("B. Apparent recalibrated performance\n")
cat("============================================\n\n")

print(
  rbind(
    Intercept_only =
      apparent_intercept,

    Intercept_plus_slope =
      apparent_full
  )
)

cat("\n")


# ============================================================
# 13. Bootstrap optimism loop
# ============================================================

metric_names <- c(
  "O_P",
  "abs_log_OP",
  "Brier",
  "MAE"
)


boot_intercept <- matrix(
  NA_real_,
  nrow = B,
  ncol = 8
)


boot_full <- matrix(
  NA_real_,
  nrow = B,
  ncol = 9
)


colnames(
  boot_intercept
) <- c(
  "train_OP",
  "test_OP",
  "train_abslogOP",
  "test_abslogOP",
  "train_Brier",
  "test_Brier",
  "train_MAE",
  "test_MAE"
)


colnames(
  boot_full
) <- c(
  "beta",
  "train_OP",
  "test_OP",
  "train_abslogOP",
  "test_abslogOP",
  "train_Brier",
  "test_Brier",
  "train_MAE",
  "test_MAE"
)


n <- nrow(
  cal_data
)


cat("============================================\n")
cat("C. Bootstrap optimism correction running\n")
cat("============================================\n\n")


for (b in seq_len(B)) {

  idx <- sample(
    seq_len(n),
    size = n,
    replace = TRUE
  )


  boot_data <- cal_data[
    idx,
    ,
    drop = FALSE
  ]


  if (
    length(
      unique(
        boot_data$status
      )
    ) < 2
  ) {

    next
  }


  # ----------------------------------------------------------
  # Intercept-only recalibration
  # ----------------------------------------------------------

  fit_i <- fit_intercept_only(
    boot_data,
    t0
  )


  if (!is.null(fit_i)) {

    pred_train_i <-
      predict_recalibrated(
        boot_data,
        fit_i
      )

    pred_test_i <-
      predict_recalibrated(
        cal_data,
        fit_i
      )


    perf_train_i <-
      evaluate_predictions(
        boot_data,
        pred_train_i,
        t0
      )

    perf_test_i <-
      evaluate_predictions(
        cal_data,
        pred_test_i,
        t0
      )


    boot_intercept[
      b,
    ] <- c(

      perf_train_i[
        "O_P"
      ],

      perf_test_i[
        "O_P"
      ],

      perf_train_i[
        "abs_log_OP"
      ],

      perf_test_i[
        "abs_log_OP"
      ],

      perf_train_i[
        "Brier"
      ],

      perf_test_i[
        "Brier"
      ],

      perf_train_i[
        "MAE"
      ],

      perf_test_i[
        "MAE"
      ]
    )
  }


  # ----------------------------------------------------------
  # Full recalibration
  # ----------------------------------------------------------

  fit_f <- fit_full_recalibration(
    boot_data,
    t0
  )


  if (!is.null(fit_f)) {

    pred_train_f <-
      predict_recalibrated(
        boot_data,
        fit_f
      )

    pred_test_f <-
      predict_recalibrated(
        cal_data,
        fit_f
      )


    perf_train_f <-
      evaluate_predictions(
        boot_data,
        pred_train_f,
        t0
      )

    perf_test_f <-
      evaluate_predictions(
        cal_data,
        pred_test_f,
        t0
      )


    boot_full[
      b,
    ] <- c(

      fit_f$beta,

      perf_train_f[
        "O_P"
      ],

      perf_test_f[
        "O_P"
      ],

      perf_train_f[
        "abs_log_OP"
      ],

      perf_test_f[
        "abs_log_OP"
      ],

      perf_train_f[
        "Brier"
      ],

      perf_test_f[
        "Brier"
      ],

      perf_train_f[
        "MAE"
      ],

      perf_test_f[
        "MAE"
      ]
    )
  }


  if (
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
}


boot_intercept <-
  as.data.frame(
    boot_intercept
  )


boot_full <-
  as.data.frame(
    boot_full
  )


# ============================================================
# 14. Valid bootstrap audit
# ============================================================

valid_intercept <- complete.cases(
  boot_intercept
)

valid_full <- complete.cases(
  boot_full
)


cat("\n")
cat("============================================\n")
cat("D. Valid bootstrap runs\n")
cat("============================================\n\n")

cat(
  "Intercept-only valid B =",
  sum(
    valid_intercept
  ),
  "\n"
)

cat(
  "Full recalibration valid B =",
  sum(
    valid_full
  ),
  "\n\n"
)


# ============================================================
# 15. Optimism correction helper
# ============================================================

correct_error_metric <- function(
    apparent,
    train,
    test
) {

  apparent <- as.numeric(apparent)

  valid <-
    is.finite(train) &
    is.finite(test)


  # Error metric:
  # lower is better.
  #
  # optimism = test error - training error
  #
  # corrected error =
  # apparent error + optimism

  optimism <- mean(
    test[valid] -
      train[valid]
  )


  corrected <-
    apparent +
    optimism


  c(
    apparent =
      apparent,

    mean_train =
      mean(
        train[valid]
      ),

    mean_test =
      mean(
        test[valid]
      ),

    optimism =
      optimism,

    corrected =
      corrected,

    valid_B =
      sum(
        valid
      )
  )
}


correct_raw_OP <- function(
    apparent,
    train,
    test
) {

  apparent <- as.numeric(apparent)

  valid <-
    is.finite(train) &
    is.finite(test)


  # Raw O/P is corrected using the mean
  # train-test bootstrap bias.
  #
  # corrected = apparent - mean(train-test)

  bias <- mean(
    train[valid] -
      test[valid]
  )


  corrected <-
    apparent -
    bias


  c(
    apparent =
      apparent,

    mean_train =
      mean(
        train[valid]
      ),

    mean_test =
      mean(
        test[valid]
      ),

    bootstrap_bias =
      bias,

    corrected =
      corrected,

    valid_B =
      sum(
        valid
      )
  )
}


# ============================================================
# 16. Intercept-only optimism correction
# ============================================================

summary_intercept_OP <-
  correct_raw_OP(

    apparent =
      apparent_intercept[
        "O_P"
      ],

    train =
      boot_intercept$train_OP,

    test =
      boot_intercept$test_OP
  )


summary_intercept_abslogOP <-
  correct_error_metric(

    apparent =
      apparent_intercept[
        "abs_log_OP"
      ],

    train =
      boot_intercept$train_abslogOP,

    test =
      boot_intercept$test_abslogOP
  )


summary_intercept_Brier <-
  correct_error_metric(

    apparent =
      apparent_intercept[
        "Brier"
      ],

    train =
      boot_intercept$train_Brier,

    test =
      boot_intercept$test_Brier
  )


summary_intercept_MAE <-
  correct_error_metric(

    apparent =
      apparent_intercept[
        "MAE"
      ],

    train =
      boot_intercept$train_MAE,

    test =
      boot_intercept$test_MAE
  )


# ============================================================
# 17. Full recalibration optimism correction
# ============================================================

summary_full_OP <-
  correct_raw_OP(

    apparent =
      apparent_full[
        "O_P"
      ],

    train =
      boot_full$train_OP,

    test =
      boot_full$test_OP
  )


summary_full_abslogOP <-
  correct_error_metric(

    apparent =
      apparent_full[
        "abs_log_OP"
      ],

    train =
      boot_full$train_abslogOP,

    test =
      boot_full$test_abslogOP
  )


summary_full_Brier <-
  correct_error_metric(

    apparent =
      apparent_full[
        "Brier"
      ],

    train =
      boot_full$train_Brier,

    test =
      boot_full$test_Brier
  )


summary_full_MAE <-
  correct_error_metric(

    apparent =
      apparent_full[
        "MAE"
      ],

    train =
      boot_full$train_MAE,

    test =
      boot_full$test_MAE
  )


# ============================================================
# 18. Final optimism-corrected summary
# ============================================================

optimism_summary <- data.frame(

  model = rep(
    c(
      "Base6 recalibrated: intercept only",
      "Base6 recalibrated: intercept + slope"
    ),
    each = 4
  ),

  metric = rep(
    c(
      "O_P",
      "abs_log_OP",
      "Brier_364",
      "Quintile_MAE"
    ),
    times = 2
  ),

  apparent = c(

    summary_intercept_OP[
      "apparent"
    ],

    summary_intercept_abslogOP[
      "apparent"
    ],

    summary_intercept_Brier[
      "apparent"
    ],

    summary_intercept_MAE[
      "apparent"
    ],

    summary_full_OP[
      "apparent"
    ],

    summary_full_abslogOP[
      "apparent"
    ],

    summary_full_Brier[
      "apparent"
    ],

    summary_full_MAE[
      "apparent"
    ]
  ),

  corrected = c(

    summary_intercept_OP[
      "corrected"
    ],

    summary_intercept_abslogOP[
      "corrected"
    ],

    summary_intercept_Brier[
      "corrected"
    ],

    summary_intercept_MAE[
      "corrected"
    ],

    summary_full_OP[
      "corrected"
    ],

    summary_full_abslogOP[
      "corrected"
    ],

    summary_full_Brier[
      "corrected"
    ],

    summary_full_MAE[
      "corrected"
    ]
  ),

  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("E. Optimism-corrected recalibration performance\n")
cat("============================================\n\n")

print(
  optimism_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 19. Bootstrap distribution of full recalibration slope
# ============================================================

beta_valid <- boot_full$beta[
  is.finite(
    boot_full$beta
  )
]


beta_boot_summary <- data.frame(

  original_beta =
    fit_original_full$beta,

  bootstrap_mean =
    mean(
      beta_valid
    ),

  CI_low =
    unname(
      quantile(
        beta_valid,
        0.025
      )
    ),

  CI_high =
    unname(
      quantile(
        beta_valid,
        0.975
      )
    ),

  valid_B =
    length(
      beta_valid
    )
)


cat("============================================\n")
cat("F. Bootstrap recalibration slope\n")
cat("============================================\n\n")

print(
  beta_boot_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 20. Save outputs
# ============================================================

write.csv(
  optimism_summary,
  "Chinese_GWTG_STEP11I_OPTIMISM_CORRECTED_RECALIBRATION.csv",
  row.names = FALSE
)


write.csv(
  beta_boot_summary,
  "Chinese_GWTG_STEP11I_RECALIBRATION_SLOPE_BOOTSTRAP.csv",
  row.names = FALSE
)


save(
  boot_intercept,
  boot_full,
  apparent_intercept,
  apparent_full,
  optimism_summary,
  beta_boot_summary,
  file =
    "Chinese_GWTG_STEP11I_RECALIBRATION_OPTIMISM.RData"
)


cat("============================================\n")
cat("STEP 11I completed\n")
cat("============================================\n\n")

cat("Saved:\n")
cat("1. Chinese_GWTG_STEP11I_OPTIMISM_CORRECTED_RECALIBRATION.csv\n")
cat("2. Chinese_GWTG_STEP11I_RECALIBRATION_SLOPE_BOOTSTRAP.csv\n")
cat("3. Chinese_GWTG_STEP11I_RECALIBRATION_OPTIMISM.RData\n")
cat("============================================\n")
