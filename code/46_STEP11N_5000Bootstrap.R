# ============================================================
# STEP 11N
# 5000-bootstrap robustness analysis
#
# Purpose:
#   Reassess key pairwise discrimination differences using
#   a larger number of paired bootstrap resamples.
#
# Primary comparisons:
#   1. Base6 vs MAGGIC
#   2. Base6 vs GWTG-HF
#
# Secondary:
#   3. GWTG-HF vs MAGGIC
#
# Metrics:
#   - 364-day time-dependent AUC
#   - Harrell C-index
#
# Design:
#   Paired bootstrap on the SAME common n=346 cohort.
#   The same resampled patient indices are used for all models
#   within each bootstrap replicate.
#
# Bootstrap:
#   B = 5000
#
# Input:
#   Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData
#
# Output:
#   Chinese_GWTG_STEP11N_5000_BOOTSTRAP_*.csv
#   Chinese_GWTG_STEP11N_5000_BOOTSTRAP.RData
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

rm(list = ls())

library(survival)

if (!requireNamespace("timeROC", quietly = TRUE)) {
  stop(
    "Package 'timeROC' is required."
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
B <- 5000


cat("\n")
cat("============================================\n")
cat("STEP 11N: 5000-bootstrap robustness\n")
cat("============================================\n\n")


# ============================================================
# 1. Load common three-model cohort
# ============================================================

e11b <- new.env()

load(
  "Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData",
  envir = e11b
)


if (
  !exists(
    "comparison3",
    envir = e11b
  )
) {

  stop(
    "ERROR: comparison3 not found."
  )
}


dat <- get(
  "comparison3",
  envir = e11b
)


required_cols <- c(
  "ID",
  "status",
  "time",
  "base6_marker",
  "maggic_score",
  "gwtg_score"
)


missing_cols <- setdiff(
  required_cols,
  names(dat)
)


if (
  length(
    missing_cols
  ) > 0
) {

  stop(
    paste0(
      "ERROR: comparison3 missing: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}


# ============================================================
# 2. Cohort audit
# ============================================================

cat("============================================\n")
cat("A. Cohort audit\n")
cat("============================================\n\n")

cat(
  "N =",
  nrow(dat),
  "\n"
)

cat(
  "Deaths =",
  sum(
    dat$status == 1
  ),
  "\n"
)

cat(
  "Censored/alive =",
  sum(
    dat$status == 0
  ),
  "\n"
)

cat(
  "Unique IDs =",
  length(
    unique(
      dat$ID
    )
  ),
  "\n"
)

cat(
  "Missing analysis rows =",
  sum(
    !complete.cases(
      dat[
        ,
        required_cols,
        drop = FALSE
      ]
    )
  ),
  "\n\n"
)


if (
  nrow(dat) != 346
) {

  stop(
    "ERROR: expected n=346."
  )
}


if (
  any(
    duplicated(
      dat$ID
    )
  )
) {

  stop(
    "ERROR: duplicated IDs."
  )
}


if (
  any(
    !complete.cases(
      dat[
        ,
        required_cols,
        drop = FALSE
      ]
    )
  )
) {

  stop(
    "ERROR: missing values in analysis dataset."
  )
}


# ============================================================
# 3. Metric functions
# ============================================================

calc_auc <- function(
    time,
    status,
    marker,
    t0 = 364
) {

  roc <- tryCatch(

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


  if (
    is.null(
      roc
    )
  ) {

    return(
      NA_real_
    )
  }


  idx <- which.min(
    abs(
      roc$times -
        t0
    )
  )


  auc <- as.numeric(
    roc$AUC[idx]
  )


  if (
    !is.finite(
      auc
    )
  ) {

    return(
      NA_real_
    )
  }


  auc
}


calc_cindex <- function(
    time,
    status,
    marker
) {

  fit <- tryCatch(

    survival::concordance(
      survival::Surv(
        time,
        status
      ) ~ marker,
      reverse = TRUE
    ),

    error = function(e) NULL
  )


  if (
    is.null(
      fit
    )
  ) {

    return(
      NA_real_
    )
  }


  cindex <- as.numeric(
    fit$concordance
  )


  if (
    !is.finite(
      cindex
    )
  ) {

    return(
      NA_real_
    )
  }


  cindex
}


# ============================================================
# 4. Point estimates
# ============================================================

auc_base6 <- calc_auc(
  dat$time,
  dat$status,
  dat$base6_marker,
  t0
)


auc_maggic <- calc_auc(
  dat$time,
  dat$status,
  dat$maggic_score,
  t0
)


auc_gwtg <- calc_auc(
  dat$time,
  dat$status,
  dat$gwtg_score,
  t0
)


c_base6 <- calc_cindex(
  dat$time,
  dat$status,
  dat$base6_marker
)


c_maggic <- calc_cindex(
  dat$time,
  dat$status,
  dat$maggic_score
)


c_gwtg <- calc_cindex(
  dat$time,
  dat$status,
  dat$gwtg_score
)


point_model_summary <- data.frame(

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
cat("B. Point-estimate discrimination\n")
cat("============================================\n\n")

print(
  point_model_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 5. Point pairwise differences
# ============================================================

point_pairwise <- data.frame(

  comparison = c(
    "Base6 - MAGGIC",
    "Base6 - GWTG-HF",
    "GWTG-HF - MAGGIC"
  ),

  delta_AUC = c(
    auc_base6 -
      auc_maggic,

    auc_base6 -
      auc_gwtg,

    auc_gwtg -
      auc_maggic
  ),

  delta_C = c(
    c_base6 -
      c_maggic,

    c_base6 -
      c_gwtg,

    c_gwtg -
      c_maggic
  ),

  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("C. Point-estimate pairwise differences\n")
cat("============================================\n\n")

print(
  point_pairwise,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 6. Allocate bootstrap storage
# ============================================================

boot_results <- matrix(
  NA_real_,
  nrow = B,
  ncol = 12
)


colnames(
  boot_results
) <- c(

  "AUC_Base6",
  "AUC_MAGGIC",
  "AUC_GWTG",

  "C_Base6",
  "C_MAGGIC",
  "C_GWTG",

  "dAUC_Base6_MAGGIC",
  "dAUC_Base6_GWTG",
  "dAUC_GWTG_MAGGIC",

  "dC_Base6_MAGGIC",
  "dC_Base6_GWTG",
  "dC_GWTG_MAGGIC"
)


# ============================================================
# 7. Paired bootstrap
# ============================================================

n <- nrow(
  dat
)


cat("============================================\n")
cat("D. 5000 paired bootstrap running\n")
cat("============================================\n\n")

cat(
  "Bootstrap repetitions =",
  B,
  "\n\n"
)


for (
  b in seq_len(
    B
  )
) {

  idx <- sample(
    seq_len(
      n
    ),
    size = n,
    replace = TRUE
  )


  d <- dat[
    idx,
    ,
    drop = FALSE
  ]


  # ----------------------------------------------------------
  # Safety check:
  # bootstrap sample must contain both outcome states
  # ----------------------------------------------------------

  if (
    length(
      unique(
        d$status
      )
    ) < 2
  ) {

    next
  }


  # ----------------------------------------------------------
  # AUC
  # ----------------------------------------------------------

  auc_b <- calc_auc(
    d$time,
    d$status,
    d$base6_marker,
    t0
  )


  auc_m <- calc_auc(
    d$time,
    d$status,
    d$maggic_score,
    t0
  )


  auc_g <- calc_auc(
    d$time,
    d$status,
    d$gwtg_score,
    t0
  )


  # ----------------------------------------------------------
  # C-index
  # ----------------------------------------------------------

  c_b <- calc_cindex(
    d$time,
    d$status,
    d$base6_marker
  )


  c_m <- calc_cindex(
    d$time,
    d$status,
    d$maggic_score
  )


  c_g <- calc_cindex(
    d$time,
    d$status,
    d$gwtg_score
  )


  # ----------------------------------------------------------
  # Store model values
  # ----------------------------------------------------------

  boot_results[
    b,
    "AUC_Base6"
  ] <- auc_b


  boot_results[
    b,
    "AUC_MAGGIC"
  ] <- auc_m


  boot_results[
    b,
    "AUC_GWTG"
  ] <- auc_g


  boot_results[
    b,
    "C_Base6"
  ] <- c_b


  boot_results[
    b,
    "C_MAGGIC"
  ] <- c_m


  boot_results[
    b,
    "C_GWTG"
  ] <- c_g


  # ----------------------------------------------------------
  # Store paired differences
  # ----------------------------------------------------------

  boot_results[
    b,
    "dAUC_Base6_MAGGIC"
  ] <- auc_b -
    auc_m


  boot_results[
    b,
    "dAUC_Base6_GWTG"
  ] <- auc_b -
    auc_g


  boot_results[
    b,
    "dAUC_GWTG_MAGGIC"
  ] <- auc_g -
    auc_m


  boot_results[
    b,
    "dC_Base6_MAGGIC"
  ] <- c_b -
    c_m


  boot_results[
    b,
    "dC_Base6_GWTG"
  ] <- c_b -
    c_g


  boot_results[
    b,
    "dC_GWTG_MAGGIC"
  ] <- c_g -
    c_m


  if (
    b %% 250 == 0
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


boot_results <- as.data.frame(
  boot_results
)


# ============================================================
# 8. Valid bootstrap audit
# ============================================================

valid_counts <- data.frame(

  variable =
    names(
      boot_results
    ),

  valid_B =
    vapply(
      boot_results,
      function(x) {

        sum(
          is.finite(
            x
          )
        )
      },
      numeric(1)
    ),

  stringsAsFactors =
    FALSE
)


cat("\n")
cat("============================================\n")
cat("E. Valid bootstrap runs\n")
cat("============================================\n\n")

print(
  valid_counts,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 9. Bootstrap summary helpers
# ============================================================

summarize_metric <- function(
    x,
    point_estimate
) {

  x <- x[
    is.finite(
      x
    )
  ]


  if (
    length(
      x
    ) == 0
  ) {

    return(
      c(
        point_estimate = point_estimate,
        bootstrap_mean = NA_real_,
        bootstrap_median = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
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


  c(

    point_estimate =
      point_estimate,

    bootstrap_mean =
      mean(
        x
      ),

    bootstrap_median =
      median(
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

    valid_B =
      length(
        x
      )
  )
}


summarize_difference <- function(
    x,
    point_estimate
) {

  x <- x[
    is.finite(
      x
    )
  ]


  if (
    length(
      x
    ) == 0
  ) {

    return(
      c(
        point_estimate = point_estimate,
        bootstrap_mean = NA_real_,
        bootstrap_median = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p_boot = NA_real_,
        probability_positive = NA_real_,
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

    point_estimate =
      point_estimate,

    bootstrap_mean =
      mean(
        x
      ),

    bootstrap_median =
      median(
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

    probability_positive =
      mean(
        x > 0
      ),

    valid_B =
      length(
        x
      )
  )
}


# ============================================================
# 10. Model-specific 5000-bootstrap CIs
# ============================================================

model_summary <- rbind(

  cbind(
    model = "Base6",
    metric = "AUC_364",
    as.data.frame(
      t(
        summarize_metric(
          boot_results$AUC_Base6,
          auc_base6
        )
      )
    )
  ),

  cbind(
    model = "MAGGIC",
    metric = "AUC_364",
    as.data.frame(
      t(
        summarize_metric(
          boot_results$AUC_MAGGIC,
          auc_maggic
        )
      )
    )
  ),

  cbind(
    model = "GWTG-HF",
    metric = "AUC_364",
    as.data.frame(
      t(
        summarize_metric(
          boot_results$AUC_GWTG,
          auc_gwtg
        )
      )
    )
  ),

  cbind(
    model = "Base6",
    metric = "Harrell_C",
    as.data.frame(
      t(
        summarize_metric(
          boot_results$C_Base6,
          c_base6
        )
      )
    )
  ),

  cbind(
    model = "MAGGIC",
    metric = "Harrell_C",
    as.data.frame(
      t(
        summarize_metric(
          boot_results$C_MAGGIC,
          c_maggic
        )
      )
    )
  ),

  cbind(
    model = "GWTG-HF",
    metric = "Harrell_C",
    as.data.frame(
      t(
        summarize_metric(
          boot_results$C_GWTG,
          c_gwtg
        )
      )
    )
  )
)


rownames(
  model_summary
) <- NULL


cat("============================================\n")
cat("F. Model-specific 5000-bootstrap CI\n")
cat("============================================\n\n")

print(
  model_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 11. Pairwise 5000-bootstrap differences
# ============================================================

pairwise_summary <- rbind(

  cbind(

    comparison =
      "Base6 - MAGGIC",

    metric =
      "AUC_364",

    as.data.frame(
      t(
        summarize_difference(
          boot_results$dAUC_Base6_MAGGIC,
          auc_base6 -
            auc_maggic
        )
      )
    )
  ),


  cbind(

    comparison =
      "Base6 - GWTG-HF",

    metric =
      "AUC_364",

    as.data.frame(
      t(
        summarize_difference(
          boot_results$dAUC_Base6_GWTG,
          auc_base6 -
            auc_gwtg
        )
      )
    )
  ),


  cbind(

    comparison =
      "GWTG-HF - MAGGIC",

    metric =
      "AUC_364",

    as.data.frame(
      t(
        summarize_difference(
          boot_results$dAUC_GWTG_MAGGIC,
          auc_gwtg -
            auc_maggic
        )
      )
    )
  ),


  cbind(

    comparison =
      "Base6 - MAGGIC",

    metric =
      "Harrell_C",

    as.data.frame(
      t(
        summarize_difference(
          boot_results$dC_Base6_MAGGIC,
          c_base6 -
            c_maggic
        )
      )
    )
  ),


  cbind(

    comparison =
      "Base6 - GWTG-HF",

    metric =
      "Harrell_C",

    as.data.frame(
      t(
        summarize_difference(
          boot_results$dC_Base6_GWTG,
          c_base6 -
            c_gwtg
        )
      )
    )
  ),


  cbind(

    comparison =
      "GWTG-HF - MAGGIC",

    metric =
      "Harrell_C",

    as.data.frame(
      t(
        summarize_difference(
          boot_results$dC_GWTG_MAGGIC,
          c_gwtg -
            c_maggic
        )
      )
    )
  )
)


rownames(
  pairwise_summary
) <- NULL


cat("============================================\n")
cat("G. Pairwise 5000-bootstrap differences\n")
cat("============================================\n\n")

print(
  pairwise_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 12. Focused robustness table
# ============================================================

focused_summary <- pairwise_summary[
  pairwise_summary$comparison %in%
    c(
      "Base6 - MAGGIC",
      "Base6 - GWTG-HF"
    ),
  ,
  drop = FALSE
]


focused_summary$CI_excludes_zero <-
  with(
    focused_summary,
    CI_low > 0 |
      CI_high < 0
  )


focused_summary$direction <-
  ifelse(
    focused_summary$point_estimate > 0,
    "Base6 higher",
    ifelse(
      focused_summary$point_estimate < 0,
      "Base6 lower",
      "No difference"
    )
  )


cat("============================================\n")
cat("H. Focused robustness summary\n")
cat("============================================\n\n")

print(
  focused_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 13. Compare with STEP 11C 2000-bootstrap results
# ============================================================

comparison_2000_5000 <- NULL


if (
  file.exists(
    "Chinese_GWTG_STEP11C_PAIRED_BOOTSTRAP.RData"
  )
) {

  e11c <- new.env()


  load(
    "Chinese_GWTG_STEP11C_PAIRED_BOOTSTRAP.RData",
    envir = e11c
  )


  if (
    exists(
      "final_pairwise_summary",
      envir = e11c
    )
  ) {

    old <- get(
      "final_pairwise_summary",
      envir = e11c
    )


    old_key <- paste(
      old$comparison,
      old$metric,
      sep = " | "
    )


    new_key <- paste(
      pairwise_summary$comparison,
      pairwise_summary$metric,
      sep = " | "
    )


    comparison_2000_5000 <- data.frame(

      comparison =
        pairwise_summary$comparison,

      metric =
        pairwise_summary$metric,

      estimate_5000 =
        pairwise_summary$point_estimate,

      CI_low_5000 =
        pairwise_summary$CI_low,

      CI_high_5000 =
        pairwise_summary$CI_high,

      p_5000 =
        pairwise_summary$p_boot,

      estimate_2000 =
        NA_real_,

      CI_low_2000 =
        NA_real_,

      CI_high_2000 =
        NA_real_,

      p_2000 =
        NA_real_,

      stringsAsFactors =
        FALSE
    )


    for (
      i in seq_len(
        nrow(
          comparison_2000_5000
        )
      )
    ) {

      key_i <- paste(
        comparison_2000_5000$comparison[i],
        comparison_2000_5000$metric[i],
        sep = " | "
      )


      j <- match(
        key_i,
        old_key
      )


      if (
        !is.na(
          j
        )
      ) {

        comparison_2000_5000$estimate_2000[i] <-
          old$estimate[j]


        comparison_2000_5000$CI_low_2000[i] <-
          old$CI_low[j]


        comparison_2000_5000$CI_high_2000[i] <-
          old$CI_high[j]


        comparison_2000_5000$p_2000[i] <-
          old$p_boot[j]
      }
    }


    cat("============================================\n")
    cat("I. 2000 vs 5000 bootstrap comparison\n")
    cat("============================================\n\n")

    print(
      comparison_2000_5000,
      row.names = FALSE
    )

    cat("\n")
  }
}


# ============================================================
# 14. Save outputs
# ============================================================

write.csv(
  point_model_summary,
  "Chinese_GWTG_STEP11N_POINT_MODEL_PERFORMANCE.csv",
  row.names = FALSE
)


write.csv(
  point_pairwise,
  "Chinese_GWTG_STEP11N_POINT_PAIRWISE.csv",
  row.names = FALSE
)


write.csv(
  valid_counts,
  "Chinese_GWTG_STEP11N_VALID_BOOTSTRAP_COUNTS.csv",
  row.names = FALSE
)


write.csv(
  model_summary,
  "Chinese_GWTG_STEP11N_MODEL_5000_BOOTSTRAP.csv",
  row.names = FALSE
)


write.csv(
  pairwise_summary,
  "Chinese_GWTG_STEP11N_PAIRWISE_5000_BOOTSTRAP.csv",
  row.names = FALSE
)


write.csv(
  focused_summary,
  "Chinese_GWTG_STEP11N_FOCUSED_ROBUSTNESS.csv",
  row.names = FALSE
)


if (
  !is.null(
    comparison_2000_5000
  )
) {

  write.csv(
    comparison_2000_5000,
    "Chinese_GWTG_STEP11N_2000_VS_5000.csv",
    row.names = FALSE
  )
}


# ============================================================
# 15. Save complete RData
# ============================================================

save(
  dat,
  point_model_summary,
  point_pairwise,
  boot_results,
  valid_counts,
  model_summary,
  pairwise_summary,
  focused_summary,
  comparison_2000_5000,
  file =
    "Chinese_GWTG_STEP11N_5000_BOOTSTRAP.RData"
)


# ============================================================
# 16. Completion
# ============================================================

cat("============================================\n")
cat("STEP 11N completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11N_POINT_MODEL_PERFORMANCE.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11N_POINT_PAIRWISE.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11N_VALID_BOOTSTRAP_COUNTS.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11N_MODEL_5000_BOOTSTRAP.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11N_PAIRWISE_5000_BOOTSTRAP.csv\n"
)

cat(
  "6. Chinese_GWTG_STEP11N_FOCUSED_ROBUSTNESS.csv\n"
)

if (
  !is.null(
    comparison_2000_5000
  )
) {

  cat(
    "7. Chinese_GWTG_STEP11N_2000_VS_5000.csv\n"
  )
}

cat(
  "8. Chinese_GWTG_STEP11N_5000_BOOTSTRAP.RData\n"
)

cat("============================================\n")
