# ============================================================
# STEP 11M
# Sensitivity analysis:
# Common cohort vs maximum available cohort
#
# Purpose:
#   Assess whether restricting all models to the common
#   three-model cohort (n = 346) materially changes
#   performance estimates.
#
# Models:
#   Base6
#   MAGGIC
#   GWTG-HF
#
# Metrics:
#   - 364-day time-dependent AUC
#   - Harrell C-index
#
# Key design:
#   GWTG-HF maximum cohort = 348
#   GWTG-HF common cohort  = 346
#
#   The common cohort is nested within the maximum cohort.
#   Therefore GWTG cohort sensitivity uses a nested paired
#   bootstrap from the 348-patient maximum cohort.
#
# Inputs:
#   Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData
#   Chinese_HF_EXTERNAL_STEP9B_VALIDATION.RData
#   Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData
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
B <- 2000


cat("\n")
cat("============================================\n")
cat("STEP 11M: Common vs maximum cohort sensitivity\n")
cat("============================================\n\n")


# ============================================================
# 1. Helper functions
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


  if (is.null(roc)) {
    return(NA_real_)
  }


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


  if (is.null(fit)) {
    return(NA_real_)
  }


  as.numeric(
    fit$concordance
  )
}


calc_perf <- function(
    data
) {

  c(

    AUC =
      calc_auc(
        data$time,
        data$status,
        data$marker,
        t0
      ),

    C =
      calc_cindex(
        data$time,
        data$status,
        data$marker
      )
  )
}


summarize_boot_diff <- function(
    x,
    point_estimate
) {

  x <- x[
    is.finite(x)
  ]


  if (length(x) == 0) {

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
      mean(x),

    CI_low =
      unname(ci[1]),

    CI_high =
      unname(ci[2]),

    p_boot =
      p_boot,

    valid_B =
      length(x)
  )
}


# ============================================================
# 2. Load STEP 11B
# ============================================================

e11b <- new.env()

load(
  "Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData",
  envir = e11b
)


if (!exists("comparison3", envir = e11b)) {
  stop("ERROR: comparison3 not found.")
}


if (!exists("gwtg_score_data", envir = e11b)) {
  stop("ERROR: gwtg_score_data not found.")
}


comparison3 <- get(
  "comparison3",
  envir = e11b
)


gwtg_score_data <- get(
  "gwtg_score_data",
  envir = e11b
)


# ============================================================
# 3. Validate GWTG maximum cohort
# ============================================================

required_gwtg <- c(
  "ID",
  "status",
  "time",
  "gwtg_score"
)


missing_gwtg <- setdiff(
  required_gwtg,
  names(gwtg_score_data)
)


if (length(missing_gwtg) > 0) {

  stop(
    paste0(
      "ERROR: gwtg_score_data missing: ",
      paste(
        missing_gwtg,
        collapse = ", "
      )
    )
  )
}


cat("============================================\n")
cat("A. GWTG scored-data audit\n")
cat("============================================\n\n")

cat(
  "GWTG maximum N =",
  nrow(gwtg_score_data),
  "\n"
)

cat(
  "Unique IDs =",
  length(
    unique(
      gwtg_score_data$ID
    )
  ),
  "\n"
)

cat(
  "Deaths =",
  sum(
    gwtg_score_data$status == 1
  ),
  "\n"
)

cat(
  "Missing GWTG score =",
  sum(
    is.na(
      gwtg_score_data$gwtg_score
    )
  ),
  "\n"
)

cat(
  "Score range =",
  min(
    gwtg_score_data$gwtg_score
  ),
  "to",
  max(
    gwtg_score_data$gwtg_score
  ),
  "\n\n"
)


if (
  nrow(gwtg_score_data) != 348
) {
  stop(
    "ERROR: expected GWTG maximum cohort n=348."
  )
}


if (
  any(
    duplicated(
      gwtg_score_data$ID
    )
  )
) {
  stop(
    "ERROR: duplicated IDs in gwtg_score_data."
  )
}


# ============================================================
# 4. Validate common cohort against GWTG maximum cohort
# ============================================================

gwtg_check <- merge(

  comparison3[
    ,
    c(
      "ID",
      "gwtg_score"
    ),
    drop = FALSE
  ],

  gwtg_score_data[
    ,
    c(
      "ID",
      "gwtg_score"
    ),
    drop = FALSE
  ],

  by = "ID",

  suffixes = c(
    "_common",
    "_maximum"
  )
)


gwtg_check$diff <-
  gwtg_check$gwtg_score_common -
  gwtg_check$gwtg_score_maximum


cat("============================================\n")
cat("B. GWTG score consistency check\n")
cat("============================================\n\n")

cat(
  "Matched common IDs =",
  nrow(gwtg_check),
  "\n"
)

cat(
  "Maximum absolute score difference =",
  max(
    abs(
      gwtg_check$diff
    )
  ),
  "\n"
)

cat(
  "Exact matches =",
  sum(
    gwtg_check$diff == 0
  ),
  "/",
  nrow(gwtg_check),
  "\n\n"
)


if (
  nrow(gwtg_check) != 346
) {
  stop(
    "ERROR: expected 346 matched GWTG common IDs."
  )
}


if (
  max(
    abs(
      gwtg_check$diff
    )
  ) != 0
) {
  stop(
    "ERROR: GWTG scores differ between maximum and common cohorts."
  )
}


# ============================================================
# 5. Construct GWTG datasets
# ============================================================

gwtg_max <- data.frame(

  ID =
    gwtg_score_data$ID,

  status =
    gwtg_score_data$status,

  time =
    gwtg_score_data$time,

  marker =
    gwtg_score_data$gwtg_score,

  stringsAsFactors =
    FALSE
)


gwtg_common <- data.frame(

  ID =
    comparison3$ID,

  status =
    comparison3$status,

  time =
    comparison3$time,

  marker =
    comparison3$gwtg_score,

  stringsAsFactors =
    FALSE
)


# ============================================================
# 6. Load Base6 maximum available cohort
# ============================================================

e_base6 <- new.env()

load(
  "Chinese_HF_EXTERNAL_STEP9B_VALIDATION.RData",
  envir = e_base6
)


if (!exists("external_data", envir = e_base6)) {
  stop(
    "ERROR: external_data not found in STEP9B."
  )
}


base6_source <- get(
  "external_data",
  envir = e_base6
)


required_base6 <- c(
  "ID",
  "status",
  "time_day",
  "risk_1y"
)


if (
  !all(
    required_base6 %in%
    names(base6_source)
  )
) {
  stop(
    "ERROR: Base6 source missing required columns."
  )
}


base6_max <- data.frame(

  ID =
    base6_source$ID,

  status =
    base6_source$status,

  time =
    base6_source$time_day,

  marker =
    base6_source$risk_1y,

  stringsAsFactors =
    FALSE
)


base6_max <- base6_max[
  complete.cases(base6_max),
  ,
  drop = FALSE
]


base6_common <- data.frame(

  ID =
    comparison3$ID,

  status =
    comparison3$status,

  time =
    comparison3$time,

  marker =
    comparison3$base6_marker,

  stringsAsFactors =
    FALSE
)


# ============================================================
# 7. Load MAGGIC maximum available cohort
# ============================================================

e_maggic <- new.env()

load(
  "Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData",
  envir = e_maggic
)


# Prefer maggic_data because it is expected to retain
# the maximum available MAGGIC cohort.
if (
  exists(
    "maggic_data",
    envir = e_maggic
  )
) {

  maggic_source <- get(
    "maggic_data",
    envir = e_maggic
  )

  maggic_source_name <- "maggic_data"

} else if (
  exists(
    "comparison_data",
    envir = e_maggic
  )
) {

  maggic_source <- get(
    "comparison_data",
    envir = e_maggic
  )

  maggic_source_name <- "comparison_data"

} else {

  stop(
    "ERROR: neither maggic_data nor comparison_data found."
  )
}


required_maggic <- c(
  "ID",
  "status",
  "time_day",
  "maggic_score"
)


if (
  !all(
    required_maggic %in%
    names(maggic_source)
  )
) {
  stop(
    paste0(
      "ERROR: ",
      maggic_source_name,
      " missing required columns."
    )
  )
}


maggic_max <- data.frame(

  ID =
    maggic_source$ID,

  status =
    maggic_source$status,

  time =
    maggic_source$time_day,

  marker =
    maggic_source$maggic_score,

  stringsAsFactors =
    FALSE
)


maggic_max <- maggic_max[
  complete.cases(maggic_max),
  ,
  drop = FALSE
]


maggic_common <- data.frame(

  ID =
    comparison3$ID,

  status =
    comparison3$status,

  time =
    comparison3$time,

  marker =
    comparison3$maggic_score,

  stringsAsFactors =
    FALSE
)


# ============================================================
# 8. Cohort audit
# ============================================================

cohort_audit <- data.frame(

  model = c(
    "Base6",
    "MAGGIC",
    "GWTG-HF"
  ),

  maximum_source = c(
    "external_data",
    maggic_source_name,
    "gwtg_score_data"
  ),

  n_maximum = c(
    nrow(base6_max),
    nrow(maggic_max),
    nrow(gwtg_max)
  ),

  deaths_maximum = c(
    sum(base6_max$status == 1),
    sum(maggic_max$status == 1),
    sum(gwtg_max$status == 1)
  ),

  n_common = c(
    nrow(base6_common),
    nrow(maggic_common),
    nrow(gwtg_common)
  ),

  deaths_common = c(
    sum(base6_common$status == 1),
    sum(maggic_common$status == 1),
    sum(gwtg_common$status == 1)
  ),

  excluded_n = c(
    nrow(base6_max) -
      nrow(base6_common),

    nrow(maggic_max) -
      nrow(maggic_common),

    nrow(gwtg_max) -
      nrow(gwtg_common)
  ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("C. Cohort audit\n")
cat("============================================\n\n")

print(
  cohort_audit,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 9. Point-estimate performance
# ============================================================

make_perf_row <- function(
    data,
    model_name,
    cohort_name
) {

  perf <- calc_perf(
    data
  )


  data.frame(

    model =
      model_name,

    cohort =
      cohort_name,

    n =
      nrow(data),

    deaths =
      sum(
        data$status == 1
      ),

    AUC_364 =
      unname(
        perf["AUC"]
      ),

    Harrell_C =
      unname(
        perf["C"]
      ),

    stringsAsFactors =
      FALSE
  )
}


performance_point <- rbind(

  make_perf_row(
    base6_max,
    "Base6",
    "Maximum available"
  ),

  make_perf_row(
    base6_common,
    "Base6",
    "Common n=346"
  ),

  make_perf_row(
    maggic_max,
    "MAGGIC",
    "Maximum available"
  ),

  make_perf_row(
    maggic_common,
    "MAGGIC",
    "Common n=346"
  ),

  make_perf_row(
    gwtg_max,
    "GWTG-HF",
    "Maximum available"
  ),

  make_perf_row(
    gwtg_common,
    "GWTG-HF",
    "Common n=346"
  )
)


rownames(
  performance_point
) <- NULL


cat("============================================\n")
cat("D. Point-estimate performance\n")
cat("============================================\n\n")

print(
  performance_point,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 10. Point differences
# maximum - common
# ============================================================

get_metric <- function(
    model_name,
    cohort_name,
    metric_name
) {

  x <- performance_point[
    performance_point$model == model_name &
      performance_point$cohort == cohort_name,
    metric_name
  ]


  if (length(x) != 1) {
    return(NA_real_)
  }


  as.numeric(x)
}


point_difference <- data.frame(

  model = c(
    "Base6",
    "MAGGIC",
    "GWTG-HF"
  ),

  delta_AUC_max_minus_common = c(

    get_metric(
      "Base6",
      "Maximum available",
      "AUC_364"
    ) -
      get_metric(
        "Base6",
        "Common n=346",
        "AUC_364"
      ),

    get_metric(
      "MAGGIC",
      "Maximum available",
      "AUC_364"
    ) -
      get_metric(
        "MAGGIC",
        "Common n=346",
        "AUC_364"
      ),

    get_metric(
      "GWTG-HF",
      "Maximum available",
      "AUC_364"
    ) -
      get_metric(
        "GWTG-HF",
        "Common n=346",
        "AUC_364"
      )
  ),

  delta_C_max_minus_common = c(

    get_metric(
      "Base6",
      "Maximum available",
      "Harrell_C"
    ) -
      get_metric(
        "Base6",
        "Common n=346",
        "Harrell_C"
      ),

    get_metric(
      "MAGGIC",
      "Maximum available",
      "Harrell_C"
    ) -
      get_metric(
        "MAGGIC",
        "Common n=346",
        "Harrell_C"
      ),

    get_metric(
      "GWTG-HF",
      "Maximum available",
      "Harrell_C"
    ) -
      get_metric(
        "GWTG-HF",
        "Common n=346",
        "Harrell_C"
      )
  ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("E. Point change: maximum minus common\n")
cat("============================================\n\n")

print(
  point_difference,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 11. Nested bootstrap
#
# Resample from maximum cohort.
# In each bootstrap sample:
#   1. calculate performance in full resampled cohort
#   2. retain rows whose original IDs belong to common cohort
#   3. calculate performance in nested common subset
#
# This preserves dependence between cohorts.
# ============================================================

nested_bootstrap <- function(
    max_data,
    common_ids,
    B = 2000
) {

  n_max <- nrow(
    max_data
  )


  out <- matrix(
    NA_real_,
    nrow = B,
    ncol = 6
  )


  colnames(out) <- c(
    "AUC_max",
    "AUC_common",
    "dAUC",
    "C_max",
    "C_common",
    "dC"
  )


  for (b in seq_len(B)) {

    idx <- sample(
      seq_len(n_max),
      size = n_max,
      replace = TRUE
    )


    d_max <- max_data[
      idx,
      ,
      drop = FALSE
    ]


    d_common <- d_max[
      d_max$ID %in%
        common_ids,
      ,
      drop = FALSE
    ]


    if (
      nrow(d_common) < 2
    ) {
      next
    }


    if (
      length(
        unique(
          d_max$status
        )
      ) < 2
    ) {
      next
    }


    if (
      length(
        unique(
          d_common$status
        )
      ) < 2
    ) {
      next
    }


    auc_max <- calc_auc(
      d_max$time,
      d_max$status,
      d_max$marker,
      t0
    )


    auc_common <- calc_auc(
      d_common$time,
      d_common$status,
      d_common$marker,
      t0
    )


    c_max <- calc_cindex(
      d_max$time,
      d_max$status,
      d_max$marker
    )


    c_common <- calc_cindex(
      d_common$time,
      d_common$status,
      d_common$marker
    )


    out[
      b,
      "AUC_max"
    ] <- auc_max


    out[
      b,
      "AUC_common"
    ] <- auc_common


    out[
      b,
      "dAUC"
    ] <- auc_max -
      auc_common


    out[
      b,
      "C_max"
    ] <- c_max


    out[
      b,
      "C_common"
    ] <- c_common


    out[
      b,
      "dC"
    ] <- c_max -
      c_common


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


  as.data.frame(
    out
  )
}


# ============================================================
# 12. Run bootstrap only if cohorts actually differ
# ============================================================

cat("============================================\n")
cat("F. Nested bootstrap cohort sensitivity\n")
cat("============================================\n\n")


boot_results <- list()


run_or_skip <- function(
    model_name,
    max_data,
    common_data
) {

  same_ids <- identical(
    sort(
      max_data$ID
    ),
    sort(
      common_data$ID
    )
  )


  if (same_ids) {

    cat(
      model_name,
      ": maximum cohort identical to common cohort; ",
      "bootstrap difference = 0 by construction.\n",
      sep = ""
    )

    return(NULL)
  }


  cat(
    "Running ",
    model_name,
    " nested bootstrap...\n",
    sep = ""
  )


  nested_bootstrap(
    max_data,
    common_data$ID,
    B
  )
}


boot_results$Base6 <- run_or_skip(
  "Base6",
  base6_max,
  base6_common
)


boot_results$MAGGIC <- run_or_skip(
  "MAGGIC",
  maggic_max,
  maggic_common
)


boot_results$GWTG_HF <- run_or_skip(
  "GWTG-HF",
  gwtg_max,
  gwtg_common
)


cat("\n")


# ============================================================
# 13. Bootstrap summary
# ============================================================

make_summary <- function(
    model_name,
    boot_df,
    point_auc,
    point_c
) {

  if (is.null(boot_df)) {

    return(
      data.frame(

        model =
          model_name,

        delta_AUC =
          point_auc,

        bootstrap_mean_AUC_diff =
          0,

        AUC_CI_low =
          0,

        AUC_CI_high =
          0,

        AUC_p_boot =
          NA_real_,

        delta_C =
          point_c,

        bootstrap_mean_C_diff =
          0,

        C_CI_low =
          0,

        C_CI_high =
          0,

        C_p_boot =
          NA_real_,

        valid_B =
          NA_integer_,

        stringsAsFactors =
          FALSE
      )
    )
  }


  auc_stats <- summarize_boot_diff(
    boot_df$dAUC,
    point_auc
  )


  c_stats <- summarize_boot_diff(
    boot_df$dC,
    point_c
  )


  data.frame(

    model =
      model_name,

    delta_AUC =
      unname(
        auc_stats["point_estimate"]
      ),

    bootstrap_mean_AUC_diff =
      unname(
        auc_stats["bootstrap_mean"]
      ),

    AUC_CI_low =
      unname(
        auc_stats["CI_low"]
      ),

    AUC_CI_high =
      unname(
        auc_stats["CI_high"]
      ),

    AUC_p_boot =
      unname(
        auc_stats["p_boot"]
      ),

    delta_C =
      unname(
        c_stats["point_estimate"]
      ),

    bootstrap_mean_C_diff =
      unname(
        c_stats["bootstrap_mean"]
      ),

    C_CI_low =
      unname(
        c_stats["CI_low"]
      ),

    C_CI_high =
      unname(
        c_stats["CI_high"]
      ),

    C_p_boot =
      unname(
        c_stats["p_boot"]
      ),

    valid_B =
      min(
        unname(
          auc_stats["valid_B"]
        ),
        unname(
          c_stats["valid_B"]
        )
      ),

    stringsAsFactors =
      FALSE
  )
}


boot_summary <- rbind(

  make_summary(

    "Base6",

    boot_results$Base6,

    point_difference$delta_AUC_max_minus_common[
      point_difference$model ==
        "Base6"
    ],

    point_difference$delta_C_max_minus_common[
      point_difference$model ==
        "Base6"
    ]
  ),


  make_summary(

    "MAGGIC",

    boot_results$MAGGIC,

    point_difference$delta_AUC_max_minus_common[
      point_difference$model ==
        "MAGGIC"
    ],

    point_difference$delta_C_max_minus_common[
      point_difference$model ==
        "MAGGIC"
    ]
  ),


  make_summary(

    "GWTG-HF",

    boot_results$GWTG_HF,

    point_difference$delta_AUC_max_minus_common[
      point_difference$model ==
        "GWTG-HF"
    ],

    point_difference$delta_C_max_minus_common[
      point_difference$model ==
        "GWTG-HF"
    ]
  )
)


rownames(
  boot_summary
) <- NULL


cat("============================================\n")
cat("G. Final cohort-restriction sensitivity\n")
cat("============================================\n\n")

print(
  boot_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 14. Excluded-patient audit
# ============================================================

excluded_base6 <- base6_max[
  !base6_max$ID %in%
    base6_common$ID,
  ,
  drop = FALSE
]


excluded_maggic <- maggic_max[
  !maggic_max$ID %in%
    maggic_common$ID,
  ,
  drop = FALSE
]


excluded_gwtg <- gwtg_max[
  !gwtg_max$ID %in%
    gwtg_common$ID,
  ,
  drop = FALSE
]


excluded_summary <- data.frame(

  model = c(
    "Base6",
    "MAGGIC",
    "GWTG-HF"
  ),

  excluded_n = c(
    nrow(excluded_base6),
    nrow(excluded_maggic),
    nrow(excluded_gwtg)
  ),

  excluded_deaths = c(
    sum(
      excluded_base6$status == 1
    ),
    sum(
      excluded_maggic$status == 1
    ),
    sum(
      excluded_gwtg$status == 1
    )
  ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("H. Excluded-patient audit\n")
cat("============================================\n\n")

print(
  excluded_summary,
  row.names = FALSE
)

cat("\n")


if (
  nrow(
    excluded_gwtg
  ) > 0
) {

  cat(
    "GWTG-HF patients excluded by common-cohort restriction:\n"
  )

  print(
    excluded_gwtg
  )

  cat("\n")
}


# ============================================================
# 15. Save outputs
# ============================================================

write.csv(
  gwtg_check,
  "Chinese_GWTG_STEP11M_GWTG_SCORE_CONSISTENCY.csv",
  row.names = FALSE
)


write.csv(
  cohort_audit,
  "Chinese_GWTG_STEP11M_COHORT_AUDIT.csv",
  row.names = FALSE
)


write.csv(
  performance_point,
  "Chinese_GWTG_STEP11M_PERFORMANCE_BY_COHORT.csv",
  row.names = FALSE
)


write.csv(
  point_difference,
  "Chinese_GWTG_STEP11M_POINT_DIFFERENCE.csv",
  row.names = FALSE
)


write.csv(
  boot_summary,
  "Chinese_GWTG_STEP11M_NESTED_BOOTSTRAP.csv",
  row.names = FALSE
)


write.csv(
  excluded_summary,
  "Chinese_GWTG_STEP11M_EXCLUDED_SUMMARY.csv",
  row.names = FALSE
)


if (
  nrow(
    excluded_gwtg
  ) > 0
) {

  write.csv(
    excluded_gwtg,
    "Chinese_GWTG_STEP11M_GWTG_EXCLUDED_PATIENTS.csv",
    row.names = FALSE
  )
}


# ============================================================
# 16. Save RData
# ============================================================

save(
  comparison3,
  gwtg_score_data,
  gwtg_check,
  base6_max,
  base6_common,
  maggic_max,
  maggic_common,
  gwtg_max,
  gwtg_common,
  cohort_audit,
  performance_point,
  point_difference,
  boot_results,
  boot_summary,
  excluded_summary,
  excluded_gwtg,
  file =
    "Chinese_GWTG_STEP11M_COHORT_SENSITIVITY.RData"
)


# ============================================================
# 17. Completion
# ============================================================

cat("============================================\n")
cat("STEP 11M completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11M_GWTG_SCORE_CONSISTENCY.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11M_COHORT_AUDIT.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11M_PERFORMANCE_BY_COHORT.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11M_POINT_DIFFERENCE.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11M_NESTED_BOOTSTRAP.csv\n"
)

cat(
  "6. Chinese_GWTG_STEP11M_EXCLUDED_SUMMARY.csv\n"
)

cat(
  "7. Chinese_GWTG_STEP11M_GWTG_EXCLUDED_PATIENTS.csv ",
  "(if applicable)\n"
)

cat(
  "8. Chinese_GWTG_STEP11M_COHORT_SENSITIVITY.RData\n"
)

cat("============================================\n")
