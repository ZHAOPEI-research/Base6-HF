# ============================================================
# STEP 11L
# Sensitivity analysis:
# Base6 vs original MAGGIC vs modernized MAGGIC
#
# Purpose:
#   Assess whether updating MAGGIC medication handling for the
#   contemporary treatment era materially changes performance.
#
# Models:
#   1. Base6
#   2. Original MAGGIC
#   3. Modernized MAGGIC
#
# Outcomes:
#   - 364-day time-dependent AUC
#   - Harrell C-index
#   - paired bootstrap differences
#   - overall calibration
#   - IPCW Brier score
#   - quintile calibration MAE
#   - selected-threshold DCA
#
# Bootstrap:
#   B = 2000
#
# Inputs:
#   Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData
#   Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData
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

t0 <- 364
B <- 2000

selected_thresholds <- c(
  0.05,
  0.10,
  0.15,
  0.20,
  0.25,
  0.30
)


cat("\n")
cat("============================================\n")
cat("STEP 11L: Modernized MAGGIC sensitivity\n")
cat("============================================\n\n")


# ============================================================
# 1. Load current common cohort
# ============================================================

e11b <- new.env()

load(
  "Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData",
  envir = e11b
)

if (!exists("comparison3", envir = e11b)) {
  stop("ERROR: comparison3 not found in STEP 11B.")
}

common3 <- get(
  "comparison3",
  envir = e11b
)


# ============================================================
# 2. Load MAGGIC dataset
# ============================================================

em <- new.env()

load(
  "Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData",
  envir = em
)

if (!exists("comparison_data", envir = em)) {
  stop(
    "ERROR: comparison_data not found in ",
    "Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData"
  )
}

maggic_source <- get(
  "comparison_data",
  envir = em
)


# ============================================================
# 3. Required columns
# ============================================================

required_common <- c(
  "ID",
  "status",
  "time",
  "base6_marker"
)

required_maggic <- c(
  "ID",
  "maggic_score",
  "risk_maggic_1y",
  "maggic_score_modern",
  "risk_maggic_modern_1y"
)


missing_common <- setdiff(
  required_common,
  names(common3)
)

missing_maggic <- setdiff(
  required_maggic,
  names(maggic_source)
)


if (length(missing_common) > 0) {
  stop(
    paste(
      "Missing from common3:",
      paste(missing_common, collapse = ", ")
    )
  )
}

if (length(missing_maggic) > 0) {
  stop(
    paste(
      "Missing from MAGGIC source:",
      paste(missing_maggic, collapse = ", ")
    )
  )
}


# ============================================================
# 4. Merge explicitly by ID
# ============================================================

analysis_data <- merge(
  common3[, required_common],
  maggic_source[, required_maggic],
  by = "ID",
  all = FALSE,
  sort = FALSE
)


analysis_data <- analysis_data[
  match(
    common3$ID,
    analysis_data$ID
  ),
  ,
  drop = FALSE
]

rownames(analysis_data) <- NULL


cat("============================================\n")
cat("A. Cohort alignment audit\n")
cat("============================================\n\n")

cat("N =", nrow(analysis_data), "\n")

cat(
  "Deaths =",
  sum(analysis_data$status == 1),
  "\n"
)

cat(
  "IDs identical to STEP11B =",
  identical(
    analysis_data$ID,
    common3$ID
  ),
  "\n"
)

cat(
  "Unique IDs =",
  length(unique(analysis_data$ID)),
  "\n\n"
)


if (nrow(analysis_data) != 346) {
  stop(
    paste0(
      "ERROR: expected n=346, obtained n=",
      nrow(analysis_data)
    )
  )
}

if (!identical(
  analysis_data$ID,
  common3$ID
)) {
  stop("ERROR: ID alignment failed.")
}


# ============================================================
# 5. Complete-data audit
# ============================================================

required_analysis <- c(
  "time",
  "status",
  "base6_marker",
  "maggic_score",
  "risk_maggic_1y",
  "maggic_score_modern",
  "risk_maggic_modern_1y"
)


if (
  any(
    !complete.cases(
      analysis_data[, required_analysis]
    )
  )
) {
  stop(
    "ERROR: missing data in common sensitivity cohort."
  )
}


for (v in c(
  "base6_marker",
  "risk_maggic_1y",
  "risk_maggic_modern_1y"
)) {

  if (
    any(
      analysis_data[[v]] < 0 |
      analysis_data[[v]] > 1
    )
  ) {
    stop(
      paste0(
        "ERROR: ",
        v,
        " contains values outside [0,1]."
      )
    )
  }
}


cat("Original MAGGIC risk:\n")
print(summary(analysis_data$risk_maggic_1y))

cat("\nModernized MAGGIC risk:\n")
print(summary(analysis_data$risk_maggic_modern_1y))

cat("\nOriginal MAGGIC score:\n")
print(summary(analysis_data$maggic_score))

cat("\nModernized MAGGIC score:\n")
print(summary(analysis_data$maggic_score_modern))

cat("\n")


# ============================================================
# 6. Audit how many patients actually changed
# ============================================================

analysis_data$score_change <-
  analysis_data$maggic_score_modern -
  analysis_data$maggic_score

analysis_data$risk_change <-
  analysis_data$risk_maggic_modern_1y -
  analysis_data$risk_maggic_1y


change_audit <- data.frame(

  n = nrow(analysis_data),

  score_changed = sum(
    analysis_data$score_change != 0
  ),

  risk_changed = sum(
    abs(
      analysis_data$risk_change
    ) > 1e-12
  ),

  mean_score_change = mean(
    analysis_data$score_change
  ),

  mean_risk_change = mean(
    analysis_data$risk_change
  ),

  max_abs_score_change = max(
    abs(
      analysis_data$score_change
    )
  ),

  max_abs_risk_change = max(
    abs(
      analysis_data$risk_change
    )
  )
)


cat("============================================\n")
cat("B. Modernization impact audit\n")
cat("============================================\n\n")

print(
  change_audit,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 7. AUC function
# ============================================================

calc_auc <- function(
    time,
    status,
    marker,
    t0 = 364
) {

  roc <- timeROC::timeROC(
    T = time,
    delta = status,
    marker = marker,
    cause = 1,
    weighting = "marginal",
    times = t0,
    iid = FALSE
  )


  idx <- which.min(
    abs(
      roc$times - t0
    )
  )


  auc <- as.numeric(
    roc$AUC[idx]
  )


  if (!is.finite(auc)) {
    return(NA_real_)
  }


  auc
}


# ============================================================
# 8. Harrell C function
# ============================================================

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
    reverse = TRUE
  )

  as.numeric(
    fit$concordance
  )
}


# ============================================================
# 9. Point-estimate discrimination
# ============================================================

auc_base6 <- calc_auc(
  analysis_data$time,
  analysis_data$status,
  analysis_data$base6_marker,
  t0
)

auc_original <- calc_auc(
  analysis_data$time,
  analysis_data$status,
  analysis_data$maggic_score,
  t0
)

auc_modern <- calc_auc(
  analysis_data$time,
  analysis_data$status,
  analysis_data$maggic_score_modern,
  t0
)


c_base6 <- calc_cindex(
  analysis_data$time,
  analysis_data$status,
  analysis_data$base6_marker
)

c_original <- calc_cindex(
  analysis_data$time,
  analysis_data$status,
  analysis_data$maggic_score
)

c_modern <- calc_cindex(
  analysis_data$time,
  analysis_data$status,
  analysis_data$maggic_score_modern
)


discrimination_point <- data.frame(

  model = c(
    "Base6",
    "MAGGIC original",
    "MAGGIC modernized"
  ),

  AUC_364 = c(
    auc_base6,
    auc_original,
    auc_modern
  ),

  Harrell_C = c(
    c_base6,
    c_original,
    c_modern
  ),

  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("C. Discrimination point estimates\n")
cat("============================================\n\n")

print(
  discrimination_point,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 10. Paired bootstrap discrimination
# ============================================================

boot_discrimination <- matrix(
  NA_real_,
  nrow = B,
  ncol = 6
)


colnames(boot_discrimination) <- c(
  "AUC_Base6",
  "AUC_MAGGIC_original",
  "AUC_MAGGIC_modern",
  "C_Base6",
  "C_MAGGIC_original",
  "C_MAGGIC_modern"
)


n <- nrow(analysis_data)


cat("============================================\n")
cat("D. Paired bootstrap running\n")
cat("============================================\n\n")

cat("B =", B, "\n\n")


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


  boot_discrimination[
    b,
    "AUC_Base6"
  ] <- tryCatch(
    calc_auc(
      d$time,
      d$status,
      d$base6_marker,
      t0
    ),
    error = function(e) NA_real_
  )


  boot_discrimination[
    b,
    "AUC_MAGGIC_original"
  ] <- tryCatch(
    calc_auc(
      d$time,
      d$status,
      d$maggic_score,
      t0
    ),
    error = function(e) NA_real_
  )


  boot_discrimination[
    b,
    "AUC_MAGGIC_modern"
  ] <- tryCatch(
    calc_auc(
      d$time,
      d$status,
      d$maggic_score_modern,
      t0
    ),
    error = function(e) NA_real_
  )


  boot_discrimination[
    b,
    "C_Base6"
  ] <- tryCatch(
    calc_cindex(
      d$time,
      d$status,
      d$base6_marker
    ),
    error = function(e) NA_real_
  )


  boot_discrimination[
    b,
    "C_MAGGIC_original"
  ] <- tryCatch(
    calc_cindex(
      d$time,
      d$status,
      d$maggic_score
    ),
    error = function(e) NA_real_
  )


  boot_discrimination[
    b,
    "C_MAGGIC_modern"
  ] <- tryCatch(
    calc_cindex(
      d$time,
      d$status,
      d$maggic_score_modern
    ),
    error = function(e) NA_real_
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


boot_discrimination <- as.data.frame(
  boot_discrimination
)


cat("\nValid bootstrap runs:\n")

print(
  colSums(
    is.finite(
      as.matrix(
        boot_discrimination
      )
    )
  )
)

cat("\n")


# ============================================================
# 11. Bootstrap CI helper
# ============================================================

ci_fun <- function(x) {

  x <- x[
    is.finite(x)
  ]


  if (length(x) == 0) {
    return(
      c(
        low = NA_real_,
        high = NA_real_
      )
    )
  }


  c(

    low = unname(
      quantile(
        x,
        0.025
      )
    ),

    high = unname(
      quantile(
        x,
        0.975
      )
    )
  )
}


# ============================================================
# 12. Model-specific CIs
# ============================================================

model_ci <- rbind(

  Base6_AUC =
    ci_fun(
      boot_discrimination$AUC_Base6
    ),

  MAGGIC_original_AUC =
    ci_fun(
      boot_discrimination$AUC_MAGGIC_original
    ),

  MAGGIC_modern_AUC =
    ci_fun(
      boot_discrimination$AUC_MAGGIC_modern
    ),

  Base6_C =
    ci_fun(
      boot_discrimination$C_Base6
    ),

  MAGGIC_original_C =
    ci_fun(
      boot_discrimination$C_MAGGIC_original
    ),

  MAGGIC_modern_C =
    ci_fun(
      boot_discrimination$C_MAGGIC_modern
    )
)


# ============================================================
# 13. Pairwise differences
# ============================================================

boot_diff <- data.frame(

  dAUC_Base6_original =
    boot_discrimination$AUC_Base6 -
    boot_discrimination$AUC_MAGGIC_original,

  dAUC_Base6_modern =
    boot_discrimination$AUC_Base6 -
    boot_discrimination$AUC_MAGGIC_modern,

  dAUC_modern_original =
    boot_discrimination$AUC_MAGGIC_modern -
    boot_discrimination$AUC_MAGGIC_original,


  dC_Base6_original =
    boot_discrimination$C_Base6 -
    boot_discrimination$C_MAGGIC_original,

  dC_Base6_modern =
    boot_discrimination$C_Base6 -
    boot_discrimination$C_MAGGIC_modern,

  dC_modern_original =
    boot_discrimination$C_MAGGIC_modern -
    boot_discrimination$C_MAGGIC_original
)


summarize_difference <- function(
    x,
    point
) {

  x <- x[
    is.finite(x)
  ]


  if (length(x) == 0) {

    return(
      c(
        point_estimate = point,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p_boot = NA_real_,
        valid_B = 0
      )
    )
  }


  ci <- quantile(
    x,
    c(
      0.025,
      0.975
    )
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

    point_estimate =
      point,

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


pairwise_summary <- rbind(

  "Base6 - MAGGIC original AUC" =
    summarize_difference(
      boot_diff$dAUC_Base6_original,
      auc_base6 - auc_original
    ),

  "Base6 - MAGGIC modernized AUC" =
    summarize_difference(
      boot_diff$dAUC_Base6_modern,
      auc_base6 - auc_modern
    ),

  "MAGGIC modernized - original AUC" =
    summarize_difference(
      boot_diff$dAUC_modern_original,
      auc_modern - auc_original
    ),

  "Base6 - MAGGIC original C" =
    summarize_difference(
      boot_diff$dC_Base6_original,
      c_base6 - c_original
    ),

  "Base6 - MAGGIC modernized C" =
    summarize_difference(
      boot_diff$dC_Base6_modern,
      c_base6 - c_modern
    ),

  "MAGGIC modernized - original C" =
    summarize_difference(
      boot_diff$dC_modern_original,
      c_modern - c_original
    )
)


cat("============================================\n")
cat("E. Model bootstrap CI\n")
cat("============================================\n\n")

print(
  round(
    model_ci,
    4
  )
)

cat("\n")


cat("============================================\n")
cat("F. Paired discrimination differences\n")
cat("============================================\n\n")

print(
  round(
    pairwise_summary,
    4
  )
)

cat("\n")


# ============================================================
# 14. KM observed mortality
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


  surv <- summary(
    fit,
    times = t0,
    extend = TRUE
  )$surv


  1 - as.numeric(
    surv
  )
}


observed_364 <- calc_km_mortality(
  analysis_data$time,
  analysis_data$status,
  t0
)


# ============================================================
# 15. IPCW Brier
# ============================================================

calc_brier <- function(
    time,
    status,
    pred,
    t0 = 364
) {

  n <- length(time)


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

    as.numeric(g)
  }


  G_t0 <- get_G(t0)


  if (
    !is.finite(G_t0) ||
    G_t0 <= 0
  ) {
    return(NA_real_)
  }


  score_terms <- rep(
    0,
    n
  )


  for (i in seq_len(n)) {

    if (
      time[i] <= t0 &&
      status[i] == 1
    ) {

      g <- get_G(
        max(
          time[i] - 1e-7,
          0
        )
      )


      if (
        !is.finite(g) ||
        g <= 0
      ) {
        return(NA_real_)
      }


      score_terms[i] <-
        ((1 - pred[i])^2) /
        g


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
# 16. Quintile MAE
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


  group <- cut(
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
      group == g
    )


    obs <- calc_km_mortality(
      data$time[idx],
      data$status[idx],
      t0
    )


    pred_mean <- mean(
      pred[idx]
    )


    errors[g] <- abs(
      obs -
        pred_mean
    )
  }


  mean(
    errors,
    na.rm = TRUE
  )
}


# ============================================================
# 17. Calibration comparison
# ============================================================

prediction_map <- list(

  "Base6" =
    analysis_data$base6_marker,

  "MAGGIC original" =
    analysis_data$risk_maggic_1y,

  "MAGGIC modernized" =
    analysis_data$risk_maggic_modern_1y
)


calibration_list <- list()

counter <- 0


for (
  model_name in names(
    prediction_map
  )
) {

  pred <- prediction_map[[model_name]]


  counter <- counter + 1


  calibration_list[[counter]] <- data.frame(

    model =
      model_name,

    n =
      nrow(
        analysis_data
      ),

    deaths =
      sum(
        analysis_data$status == 1
      ),

    observed_KM_364 =
      observed_364,

    mean_predicted =
      mean(
        pred
      ),

    O_P =
      observed_364 /
      mean(
        pred
      ),

    absolute_difference =
      abs(
        observed_364 -
          mean(
            pred
          )
      ),

    Brier_364 =
      calc_brier(
        analysis_data$time,
        analysis_data$status,
        pred,
        t0
      ),

    quintile_MAE =
      calc_quintile_mae(
        analysis_data,
        pred,
        t0
      ),

    stringsAsFactors =
      FALSE
  )
}


calibration_summary <- do.call(
  rbind,
  calibration_list
)

rownames(
  calibration_summary
) <- NULL


cat("============================================\n")
cat("G. Calibration sensitivity\n")
cat("============================================\n\n")

print(
  calibration_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 18. DCA helper functions
# ============================================================

build_censoring_function <- function(
    time,
    status
) {

  censor_fit <- survival::survfit(
    survival::Surv(
      time,
      1 - status
    ) ~ 1
  )


  function(t) {

    out <- summary(
      censor_fit,
      times = t,
      extend = TRUE
    )$surv

    as.numeric(out)
  }
}


calc_net_benefit <- function(
    data,
    pred,
    threshold,
    t0 = 364
) {

  n <- nrow(data)

  get_G <- build_censoring_function(
    data$time,
    data$status
  )


  G_t0 <- get_G(t0)


  if (
    !is.finite(G_t0) ||
    G_t0 <= 0
  ) {
    return(NA_real_)
  }


  w_event <- rep(
    0,
    n
  )

  w_nonevent <- rep(
    0,
    n
  )


  for (i in seq_len(n)) {

    if (
      data$time[i] <= t0 &&
      data$status[i] == 1
    ) {

      g <- get_G(
        max(
          data$time[i] - 1e-7,
          0
        )
      )


      if (
        is.finite(g) &&
        g > 0
      ) {

        w_event[i] <- 1 / g
      }


    } else if (
      data$time[i] > t0
    ) {

      w_nonevent[i] <-
        1 /
        G_t0
    }
  }


  treat <- pred >= threshold


  TP <- sum(
    w_event[
      treat
    ]
  ) / n


  FP <- sum(
    w_nonevent[
      treat
    ]
  ) / n


  TP -
    FP *
    threshold /
    (1 - threshold)
}


# ============================================================
# 19. Selected-threshold DCA
# ============================================================

dca_rows <- list()

counter <- 0


for (
  model_name in names(
    prediction_map
  )
) {

  pred <- prediction_map[[model_name]]


  for (pt in selected_thresholds) {

    counter <- counter + 1


    dca_rows[[counter]] <- data.frame(

      model =
        model_name,

      threshold =
        pt,

      net_benefit =
        calc_net_benefit(
          analysis_data,
          pred,
          pt,
          t0
        ),

      stringsAsFactors =
        FALSE
    )
  }
}


dca_selected <- do.call(
  rbind,
  dca_rows
)


rownames(
  dca_selected
) <- NULL


cat("============================================\n")
cat("H. Selected-threshold DCA\n")
cat("============================================\n\n")

print(
  dca_selected,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 20. Modernized minus original DCA
# ============================================================

dca_original <- dca_selected[
  dca_selected$model ==
    "MAGGIC original",
  c(
    "threshold",
    "net_benefit"
  )
]

names(
  dca_original
)[2] <- "NB_original"


dca_modern <- dca_selected[
  dca_selected$model ==
    "MAGGIC modernized",
  c(
    "threshold",
    "net_benefit"
  )
]

names(
  dca_modern
)[2] <- "NB_modern"


dca_change <- merge(
  dca_original,
  dca_modern,
  by = "threshold"
)


dca_change$delta_NB_modern_minus_original <-
  dca_change$NB_modern -
  dca_change$NB_original


cat("============================================\n")
cat("I. Modernized minus original MAGGIC DCA\n")
cat("============================================\n\n")

print(
  dca_change,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 21. Final publication-style sensitivity table
# ============================================================

final_summary <- merge(
  discrimination_point,
  calibration_summary,
  by = "model",
  all = TRUE,
  sort = FALSE
)


final_summary <- final_summary[
  match(
    c(
      "Base6",
      "MAGGIC original",
      "MAGGIC modernized"
    ),
    final_summary$model
  ),
  ,
  drop = FALSE
]


rownames(
  final_summary
) <- NULL


cat("============================================\n")
cat("J. Final modernized MAGGIC sensitivity summary\n")
cat("============================================\n\n")

print(
  final_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 22. Save outputs
# ============================================================

write.csv(
  change_audit,
  "Chinese_GWTG_STEP11L_MODERN_MAGGIC_CHANGE_AUDIT.csv",
  row.names = FALSE
)


write.csv(
  discrimination_point,
  "Chinese_GWTG_STEP11L_DISCRIMINATION.csv",
  row.names = FALSE
)


write.csv(
  data.frame(
    comparison =
      rownames(
        pairwise_summary
      ),
    pairwise_summary,
    row.names = NULL
  ),
  "Chinese_GWTG_STEP11L_PAIRED_BOOTSTRAP.csv",
  row.names = FALSE
)


write.csv(
  calibration_summary,
  "Chinese_GWTG_STEP11L_CALIBRATION.csv",
  row.names = FALSE
)


write.csv(
  dca_selected,
  "Chinese_GWTG_STEP11L_DCA_SELECTED_THRESHOLDS.csv",
  row.names = FALSE
)


write.csv(
  dca_change,
  "Chinese_GWTG_STEP11L_DCA_MODERN_MINUS_ORIGINAL.csv",
  row.names = FALSE
)


write.csv(
  final_summary,
  "Chinese_GWTG_STEP11L_FINAL_SUMMARY.csv",
  row.names = FALSE
)


save(
  analysis_data,
  change_audit,
  discrimination_point,
  boot_discrimination,
  boot_diff,
  model_ci,
  pairwise_summary,
  calibration_summary,
  dca_selected,
  dca_change,
  final_summary,
  file =
    "Chinese_GWTG_STEP11L_MODERN_MAGGIC_SENSITIVITY.RData"
)


# ============================================================
# 23. Completion
# ============================================================

cat("============================================\n")
cat("STEP 11L completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11L_MODERN_MAGGIC_CHANGE_AUDIT.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11L_DISCRIMINATION.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11L_PAIRED_BOOTSTRAP.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11L_CALIBRATION.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11L_DCA_SELECTED_THRESHOLDS.csv\n"
)

cat(
  "6. Chinese_GWTG_STEP11L_DCA_MODERN_MINUS_ORIGINAL.csv\n"
)

cat(
  "7. Chinese_GWTG_STEP11L_FINAL_SUMMARY.csv\n"
)

cat(
  "8. Chinese_GWTG_STEP11L_MODERN_MAGGIC_SENSITIVITY.RData\n"
)

cat("============================================\n")
