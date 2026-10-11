# ============================================================
# STEP 11K
# Prediction-horizon sensitivity analysis
#
# Horizons:
#   360 days
#   364 days (primary analysis)
#   365 days
#
# Models:
#   Base6
#   MAGGIC
#   GWTG-HF
#
# Outcomes:
#   1. Time-dependent AUC
#   2. Paired differences in AUC
#   3. Kaplan-Meier observed mortality
#   4. Mean predicted 1-year mortality
#   5. O/P ratio
#   6. IPCW Brier score
#
# IMPORTANT:
#   Maximum follow-up is 365 days.
#   Therefore discrimination/Brier at exactly 365 days
#   may not be estimable because no subject has T > 365.
#
# Inputs:
#   Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData
#   Chinese_GWTG_STEP11D_CALIBRATION.RData
#
# Bootstrap:
#   B = 2000
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

rm(list = ls())

library(survival)

if (!requireNamespace("timeROC", quietly = TRUE)) {
  stop(
    "Package 'timeROC' is required. ",
    "Install with install.packages('timeROC')."
  )
}

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

set.seed(20260906)

horizons <- c(
  360,
  364,
  365
)

B <- 2000


cat("\n")
cat("============================================\n")
cat("STEP 11K: Horizon sensitivity analysis\n")
cat("============================================\n\n")


# ============================================================
# 1. Load latest three-model common cohort
# ============================================================

e11b <- new.env()

load(
  "Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData",
  envir = e11b
)


if (!exists("comparison3", envir = e11b)) {

  stop(
    "ERROR: comparison3 not found in STEP 11B RData."
  )
}


comparison3 <- get(
  "comparison3",
  envir = e11b
)


# ============================================================
# 2. Load probability calibration data
# ============================================================

e11d <- new.env()

load(
  "Chinese_GWTG_STEP11D_CALIBRATION.RData",
  envir = e11d
)


if (!exists("cal_data", envir = e11d)) {

  stop(
    "ERROR: cal_data not found in STEP 11D RData."
  )
}


cal_data <- get(
  "cal_data",
  envir = e11d
)


# ============================================================
# 3. Source-data audit
# ============================================================

required_comparison <- c(
  "ID",
  "status",
  "time",
  "base6_marker",
  "maggic_score",
  "gwtg_score"
)


missing_comparison <- setdiff(
  required_comparison,
  names(comparison3)
)


if (length(missing_comparison) > 0) {

  stop(
    paste0(
      "Missing from comparison3: ",
      paste(
        missing_comparison,
        collapse = ", "
      )
    )
  )
}


required_cal <- c(
  "ID",
  "risk_base6_1y",
  "risk_maggic_1y"
)


missing_cal <- setdiff(
  required_cal,
  names(cal_data)
)


if (length(missing_cal) > 0) {

  stop(
    paste0(
      "Missing from cal_data: ",
      paste(
        missing_cal,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 4. Merge probability risks into common cohort
# ============================================================

risk_data <- cal_data[
  ,
  required_cal,
  drop = FALSE
]


analysis_data <- merge(
  comparison3,
  risk_data,
  by = "ID",
  all.x = TRUE,
  sort = FALSE
)


# Restore STEP 11B order
analysis_data <- analysis_data[
  match(
    comparison3$ID,
    analysis_data$ID
  ),
  ,
  drop = FALSE
]


rownames(analysis_data) <- NULL


# ============================================================
# 5. Alignment audit
# ============================================================

cat("============================================\n")
cat("A. Cohort alignment audit\n")
cat("============================================\n\n")


cat(
  "N =",
  nrow(analysis_data),
  "\n"
)

cat(
  "Deaths =",
  sum(
    analysis_data$status == 1
  ),
  "\n"
)

cat(
  "IDs identical to STEP 11B =",
  identical(
    analysis_data$ID,
    comparison3$ID
  ),
  "\n"
)

cat(
  "Follow-up min =",
  min(
    analysis_data$time
  ),
  "\n"
)

cat(
  "Follow-up max =",
  max(
    analysis_data$time
  ),
  "\n\n"
)


if (nrow(analysis_data) != 346) {

  stop(
    "ERROR: analysis cohort is not n=346."
  )
}


if (
  !identical(
    analysis_data$ID,
    comparison3$ID
  )
) {

  stop(
    "ERROR: ID alignment failed."
  )
}


analysis_cols <- c(
  "status",
  "time",
  "base6_marker",
  "maggic_score",
  "gwtg_score",
  "risk_base6_1y",
  "risk_maggic_1y"
)


if (
  any(
    !complete.cases(
      analysis_data[
        ,
        analysis_cols,
        drop = FALSE
      ]
    )
  )
) {

  stop(
    "ERROR: missing analysis data after merge."
  )
}


# ============================================================
# 6. Check Base6 risk consistency
# ============================================================

max_base6_diff <- max(
  abs(
    analysis_data$base6_marker -
      analysis_data$risk_base6_1y
  )
)


cat(
  "Maximum difference between STEP11B Base6 marker ",
  "and STEP11D Base6 risk = ",
  max_base6_diff,
  "\n\n",
  sep = ""
)


if (
  max_base6_diff > 1e-10
) {

  warning(
    "Base6 marker differs between STEP11B and STEP11D."
  )
}


# ============================================================
# 7. Horizon estimability audit
# ============================================================

estimability <- data.frame(

  horizon = horizons,

  n_followup_gt_horizon =
    sapply(
      horizons,
      function(tt) {

        sum(
          analysis_data$time > tt
        )
      }
    ),

  n_followup_ge_horizon =
    sapply(
      horizons,
      function(tt) {

        sum(
          analysis_data$time >= tt
        )
      }
    ),

  n_events_by_horizon =
    sapply(
      horizons,
      function(tt) {

        sum(
          analysis_data$status == 1 &
            analysis_data$time <= tt
        )
      }
    ),

  stringsAsFactors = FALSE
)


estimability$AUC_estimable <-
  estimability$n_followup_gt_horizon > 0


cat("============================================\n")
cat("B. Horizon estimability audit\n")
cat("============================================\n\n")

print(
  estimability,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 8. Kaplan-Meier mortality
# ============================================================

calc_km_mortality <- function(
    time,
    status,
    t0
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
# 9. Robust time-dependent AUC
# ============================================================

calc_auc <- function(
    time,
    status,
    marker,
    t0
) {

  # cumulative/dynamic AUC requires subjects known
  # event-free beyond the evaluation horizon

  if (
    sum(
      time > t0
    ) == 0
  ) {

    return(
      NA_real_
    )
  }


  out <- tryCatch(

    timeROC::timeROC(
      T = time,
      delta = status,
      marker = marker,
      cause = 1,
      weighting = "marginal",
      times = t0,
      iid = FALSE
    ),

    error = function(e) NULL
  )


  if (is.null(out)) {

    return(
      NA_real_
    )
  }


  idx <- which.min(
    abs(
      out$times -
        t0
    )
  )


  auc <- as.numeric(
    out$AUC[
      idx
    ]
  )


  if (!is.finite(auc)) {

    return(
      NA_real_
    )
  }


  auc
}


# ============================================================
# 10. Correct IPCW Brier score
# ============================================================

calc_brier <- function(
    time,
    status,
    pred,
    t0
) {

  n <- length(
    time
  )


  # At terminal maximum follow-up there are no known
  # event-free observations beyond t0.
  # Do not force an unstable terminal Brier estimate.

  if (
    sum(
      time > t0
    ) == 0
  ) {

    return(
      NA_real_
    )
  }


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


  score_terms <- rep(
    0,
    n
  )


  for (i in seq_len(n)) {

    # --------------------------------------------------------
    # Observed death by horizon
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
        !is.finite(g) ||
        g <= 0
      ) {

        return(
          NA_real_
        )
      }


      score_terms[i] <-
        ((1 - pred[i])^2) /
        g


      # --------------------------------------------------------
      # Known event-free beyond horizon
      # --------------------------------------------------------

    } else if (
      time[i] > t0
    ) {

      score_terms[i] <-
        pred[i]^2 /
        G_t0


      # --------------------------------------------------------
      # Censored at/before horizon:
      # zero contribution, remains in denominator
      # --------------------------------------------------------

    } else {

      score_terms[i] <- 0
    }
  }


  sum(
    score_terms
  ) / n
}


# ============================================================
# 11. Point-estimate discrimination by horizon
# ============================================================

discrimination_list <- list()

counter <- 0


for (tt in horizons) {

  counter <- counter + 1


  discrimination_list[[counter]] <- data.frame(

    horizon =
      tt,

    model =
      "Base6",

    AUC =
      calc_auc(
        analysis_data$time,
        analysis_data$status,
        analysis_data$base6_marker,
        tt
      ),

    stringsAsFactors =
      FALSE
  )


  counter <- counter + 1


  discrimination_list[[counter]] <- data.frame(

    horizon =
      tt,

    model =
      "MAGGIC",

    AUC =
      calc_auc(
        analysis_data$time,
        analysis_data$status,
        analysis_data$maggic_score,
        tt
      ),

    stringsAsFactors =
      FALSE
  )


  counter <- counter + 1


  discrimination_list[[counter]] <- data.frame(

    horizon =
      tt,

    model =
      "GWTG-HF",

    AUC =
      calc_auc(
        analysis_data$time,
        analysis_data$status,
        analysis_data$gwtg_score,
        tt
      ),

    stringsAsFactors =
      FALSE
  )
}


discrimination_horizon <-
  do.call(
    rbind,
    discrimination_list
  )


rownames(
  discrimination_horizon
) <- NULL


cat("============================================\n")
cat("C. AUC by prediction horizon\n")
cat("============================================\n\n")

print(
  discrimination_horizon,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 12. Calibration / Brier by horizon
#
# NOTE:
# risk_base6_1y and risk_maggic_1y remain the fixed
# nominal 1-year predicted probabilities.
#
# This evaluates whether conclusions are sensitive to
# defining the observed 1-year endpoint at 360/364/365 days.
# ============================================================

calibration_list <- list()

counter <- 0


for (tt in horizons) {

  observed <- calc_km_mortality(
    analysis_data$time,
    analysis_data$status,
    tt
  )


  # ----------------------------------------------------------
  # Base6
  # ----------------------------------------------------------

  counter <- counter + 1


  pred_base6 <- mean(
    analysis_data$risk_base6_1y
  )


  calibration_list[[counter]] <- data.frame(

    horizon =
      tt,

    model =
      "Base6",

    observed_KM =
      observed,

    mean_predicted =
      pred_base6,

    O_P =
      observed /
      pred_base6,

    absolute_difference =
      abs(
        observed -
          pred_base6
      ),

    Brier =
      calc_brier(
        analysis_data$time,
        analysis_data$status,
        analysis_data$risk_base6_1y,
        tt
      ),

    stringsAsFactors =
      FALSE
  )


  # ----------------------------------------------------------
  # MAGGIC
  # ----------------------------------------------------------

  counter <- counter + 1


  pred_maggic <- mean(
    analysis_data$risk_maggic_1y
  )


  calibration_list[[counter]] <- data.frame(

    horizon =
      tt,

    model =
      "MAGGIC",

    observed_KM =
      observed,

    mean_predicted =
      pred_maggic,

    O_P =
      observed /
      pred_maggic,

    absolute_difference =
      abs(
        observed -
          pred_maggic
      ),

    Brier =
      calc_brier(
        analysis_data$time,
        analysis_data$status,
        analysis_data$risk_maggic_1y,
        tt
      ),

    stringsAsFactors =
      FALSE
  )
}


calibration_horizon <-
  do.call(
    rbind,
    calibration_list
  )


rownames(
  calibration_horizon
) <- NULL


cat("============================================\n")
cat("D. Calibration/Brier by prediction horizon\n")
cat("============================================\n\n")

print(
  calibration_horizon,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 13. Point-estimate AUC differences
# ============================================================

get_auc_value <- function(
    horizon_value,
    model_name
) {

  x <- discrimination_horizon$AUC[
    discrimination_horizon$horizon ==
      horizon_value &
      discrimination_horizon$model ==
      model_name
  ]


  if (length(x) != 1) {

    return(
      NA_real_
    )
  }


  as.numeric(
    x
  )
}


pairwise_point_list <- list()

counter <- 0


for (tt in horizons) {

  auc_base6 <- get_auc_value(
    tt,
    "Base6"
  )

  auc_maggic <- get_auc_value(
    tt,
    "MAGGIC"
  )

  auc_gwtg <- get_auc_value(
    tt,
    "GWTG-HF"
  )


  counter <- counter + 1


  pairwise_point_list[[counter]] <-
    data.frame(

      horizon =
        tt,

      comparison =
        "Base6 - MAGGIC",

      delta_AUC =
        auc_base6 -
        auc_maggic,

      stringsAsFactors =
        FALSE
    )


  counter <- counter + 1


  pairwise_point_list[[counter]] <-
    data.frame(

      horizon =
        tt,

      comparison =
        "Base6 - GWTG-HF",

      delta_AUC =
        auc_base6 -
        auc_gwtg,

      stringsAsFactors =
        FALSE
    )


  counter <- counter + 1


  pairwise_point_list[[counter]] <-
    data.frame(

      horizon =
        tt,

      comparison =
        "GWTG-HF - MAGGIC",

      delta_AUC =
        auc_gwtg -
        auc_maggic,

      stringsAsFactors =
        FALSE
    )
}


pairwise_auc_point <-
  do.call(
    rbind,
    pairwise_point_list
  )


rownames(
  pairwise_auc_point
) <- NULL


cat("============================================\n")
cat("E. Point-estimate AUC differences\n")
cat("============================================\n\n")

print(
  pairwise_auc_point,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 14. Paired bootstrap AUC sensitivity
#
# 365 days will generally remain NA because no patients
# have follow-up >365 days.
# ============================================================

cat("============================================\n")
cat("F. Paired bootstrap AUC sensitivity running\n")
cat("============================================\n\n")

cat(
  "Bootstrap repetitions =",
  B,
  "\n\n"
)


n <- nrow(
  analysis_data
)


boot_rows <- list()

boot_counter <- 0


for (b in seq_len(B)) {

  idx <- sample(
    seq_len(n),
    size = n,
    replace = TRUE
  )


  d <- analysis_data[
    idx,
    ,
    drop = FALSE
  ]


  if (
    length(
      unique(
        d$status
      )
    ) < 2
  ) {

    next
  }


  for (tt in horizons) {

    auc_base6 <- calc_auc(
      d$time,
      d$status,
      d$base6_marker,
      tt
    )


    auc_maggic <- calc_auc(
      d$time,
      d$status,
      d$maggic_score,
      tt
    )


    auc_gwtg <- calc_auc(
      d$time,
      d$status,
      d$gwtg_score,
      tt
    )


    boot_counter <- boot_counter + 1


    boot_rows[[boot_counter]] <-
      data.frame(

        bootstrap_id =
          b,

        horizon =
          tt,

        AUC_Base6 =
          auc_base6,

        AUC_MAGGIC =
          auc_maggic,

        AUC_GWTG =
          auc_gwtg,

        dAUC_Base6_MAGGIC =
          auc_base6 -
          auc_maggic,

        dAUC_Base6_GWTG =
          auc_base6 -
          auc_gwtg,

        dAUC_GWTG_MAGGIC =
          auc_gwtg -
          auc_maggic,

        stringsAsFactors =
          FALSE
      )
  }


  if (
    b %% 200 == 0
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


boot_horizon <-
  do.call(
    rbind,
    boot_rows
  )


rownames(
  boot_horizon
) <- NULL


# ============================================================
# 15. Bootstrap summary helper
# ============================================================

summarize_boot <- function(
    x
) {

  x <- x[
    is.finite(
      x
    )
  ]


  if (
    length(x) == 0
  ) {

    return(
      c(
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


  p_boot <- 2 * min(
    mean(
      x <= 0
    ),
    mean(
      x >= 0
    )
  )


  p_boot <- min(
    p_boot,
    1
  )


  c(
    bootstrap_mean =
      mean(
        x
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
      length(
        x
      )
  )
}


# ============================================================
# 16. Pairwise bootstrap AUC differences
# ============================================================

comparison_specs <- list(

  "Base6 - MAGGIC" =
    "dAUC_Base6_MAGGIC",

  "Base6 - GWTG-HF" =
    "dAUC_Base6_GWTG",

  "GWTG-HF - MAGGIC" =
    "dAUC_GWTG_MAGGIC"
)


pairwise_boot_list <- list()

counter <- 0


for (tt in horizons) {

  for (
    comp_name in names(
      comparison_specs
    )
  ) {

    col_name <-
      comparison_specs[[comp_name]]


    x <- boot_horizon[
      boot_horizon$horizon == tt,
      col_name
    ]


    stats <- summarize_boot(
      x
    )


    point_estimate <-
      pairwise_auc_point$delta_AUC[
        pairwise_auc_point$horizon == tt &
          pairwise_auc_point$comparison == comp_name
      ]


    counter <- counter + 1


    pairwise_boot_list[[counter]] <-
      data.frame(

        horizon =
          tt,

        comparison =
          comp_name,

        point_estimate =
          ifelse(
            length(
              point_estimate
            ) == 1,
            point_estimate,
            NA_real_
          ),

        bootstrap_mean =
          stats[
            "bootstrap_mean"
          ],

        CI_low =
          stats[
            "CI_low"
          ],

        CI_high =
          stats[
            "CI_high"
          ],

        p_boot =
          stats[
            "p_boot"
          ],

        valid_B =
          stats[
            "valid_B"
          ],

        stringsAsFactors =
          FALSE
      )
  }
}


pairwise_auc_bootstrap <-
  do.call(
    rbind,
    pairwise_boot_list
  )


rownames(
  pairwise_auc_bootstrap
) <- NULL


cat("\n")
cat("============================================\n")
cat("G. Paired bootstrap AUC sensitivity\n")
cat("============================================\n\n")

print(
  pairwise_auc_bootstrap,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 17. Horizon consistency summary
# ============================================================

horizon_summary <- merge(
  discrimination_horizon,
  calibration_horizon,
  by = c(
    "horizon",
    "model"
  ),
  all = TRUE,
  sort = FALSE
)


horizon_summary <- horizon_summary[
  order(
    horizon_summary$horizon,
    match(
      horizon_summary$model,
      c(
        "Base6",
        "MAGGIC",
        "GWTG-HF"
      )
    )
  ),
  ,
  drop = FALSE
]


rownames(
  horizon_summary
) <- NULL


cat("============================================\n")
cat("H. Final horizon sensitivity summary\n")
cat("============================================\n\n")

print(
  horizon_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 18. Explicit interpretation flags
# ============================================================

terminal_horizon_note <- data.frame(

  horizon =
    horizons,

  AUC_estimable =
    estimability$AUC_estimable,

  interpretation =
    ifelse(
      estimability$AUC_estimable,
      "Full discrimination analysis estimable",
      paste0(
        "Terminal horizon: no follow-up beyond ",
        horizons,
        " days; AUC/Brier not reliably estimable"
      )
    ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("I. Estimability interpretation\n")
cat("============================================\n\n")

print(
  terminal_horizon_note,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 19. Save CSV outputs
# ============================================================

write.csv(
  estimability,
  "Chinese_GWTG_STEP11K_HORIZON_ESTIMABILITY.csv",
  row.names = FALSE
)


write.csv(
  discrimination_horizon,
  "Chinese_GWTG_STEP11K_AUC_BY_HORIZON.csv",
  row.names = FALSE
)


write.csv(
  calibration_horizon,
  "Chinese_GWTG_STEP11K_CALIBRATION_BY_HORIZON.csv",
  row.names = FALSE
)


write.csv(
  pairwise_auc_bootstrap,
  "Chinese_GWTG_STEP11K_PAIRED_AUC_BOOTSTRAP.csv",
  row.names = FALSE
)


write.csv(
  horizon_summary,
  "Chinese_GWTG_STEP11K_HORIZON_SUMMARY.csv",
  row.names = FALSE
)


write.csv(
  terminal_horizon_note,
  "Chinese_GWTG_STEP11K_HORIZON_INTERPRETATION.csv",
  row.names = FALSE
)


# ============================================================
# 20. Save RData
# ============================================================

save(
  analysis_data,
  estimability,
  discrimination_horizon,
  calibration_horizon,
  pairwise_auc_point,
  boot_horizon,
  pairwise_auc_bootstrap,
  horizon_summary,
  terminal_horizon_note,
  file =
    "Chinese_GWTG_STEP11K_HORIZON_SENSITIVITY.RData"
)


# ============================================================
# 21. Completion
# ============================================================

cat("============================================\n")
cat("STEP 11K completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11K_HORIZON_ESTIMABILITY.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11K_AUC_BY_HORIZON.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11K_CALIBRATION_BY_HORIZON.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11K_PAIRED_AUC_BOOTSTRAP.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11K_HORIZON_SUMMARY.csv\n"
)

cat(
  "6. Chinese_GWTG_STEP11K_HORIZON_INTERPRETATION.csv\n"
)

cat(
  "7. Chinese_GWTG_STEP11K_HORIZON_SENSITIVITY.RData\n"
)

cat("\n")

cat(
  "IMPORTANT: 365-day discrimination/Brier may be NA ",
  "because the maximum observed follow-up is exactly ",
  "365 days. This is expected and supports the use of ",
  "364 days as the primary evaluable 1-year horizon.\n"
)

cat("============================================\n")
