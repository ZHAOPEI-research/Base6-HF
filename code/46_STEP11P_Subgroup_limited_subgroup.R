# ============================================================
# STEP 11P
# Subgroup event-count audit + limited subgroup analysis
#
# Purpose:
#   Explore whether discrimination patterns are broadly
#   consistent across clinically relevant subgroups, while
#   avoiding over-interpretation in small-event strata.
#
# Candidate subgroups:
#   - Sex
#   - Age <70 vs >=70
#   - LVEF <50 vs >=50
#   - CKD
#   - DM
#   - AF
#
# Rule:
#   Formal subgroup performance is only calculated when
#   subgroup event count >= 15.
#
# Metrics:
#   - 364-day time-dependent AUC
#   - Harrell C-index
#
# Input:
#   Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData
#
# For subgroup variables, use gwtg_score_data from STEP11B.
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

rm(list = ls())

library(survival)

if (!requireNamespace("timeROC", quietly = TRUE)) {
  stop("Package 'timeROC' is required.")
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

t0 <- 364
min_events <- 15


cat("\n")
cat("============================================\n")
cat("STEP 11P: Subgroup event-count audit\n")
cat("============================================\n\n")


# ============================================================
# 1. Load STEP11B
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
# 2. Validate required columns
# ============================================================

required_common <- c(
  "ID",
  "status",
  "time",
  "base6_marker",
  "maggic_score",
  "gwtg_score"
)


missing_common <- setdiff(
  required_common,
  names(comparison3)
)


if (length(missing_common) > 0) {
  stop(
    paste0(
      "comparison3 missing: ",
      paste(missing_common, collapse = ", ")
    )
  )
}


required_subgroup <- c(
  "ID",
  "sex",
  "age",
  "LVEF",
  "CKD",
  "DM",
  "AF"
)


missing_subgroup <- setdiff(
  required_subgroup,
  names(gwtg_score_data)
)


if (length(missing_subgroup) > 0) {
  stop(
    paste0(
      "gwtg_score_data missing: ",
      paste(missing_subgroup, collapse = ", ")
    )
  )
}


# ============================================================
# 3. Merge subgroup variables into common n=346 cohort
# ============================================================

subgroup_vars <- gwtg_score_data[
  ,
  required_subgroup,
  drop = FALSE
]


dat <- merge(
  comparison3,
  subgroup_vars,
  by = "ID",
  all.x = TRUE,
  sort = FALSE
)


# restore original common-cohort order
dat <- dat[
  match(
    comparison3$ID,
    dat$ID
  ),
  ,
  drop = FALSE
]


rownames(dat) <- NULL


if (nrow(dat) != nrow(comparison3)) {
  stop("ERROR: subgroup merge changed cohort size.")
}


cat("============================================\n")
cat("A. Common cohort audit\n")
cat("============================================\n\n")

cat("N =", nrow(dat), "\n")
cat("Deaths =", sum(dat$status == 1), "\n")
cat("Alive/censored =", sum(dat$status == 0), "\n\n")


# ============================================================
# 4. Standardize subgroup variables
# ============================================================

# ------------------------------------------------------------
# Sex
# Based on source coding:
#   sex = 1 male
#   sex = 2 female
# ------------------------------------------------------------

dat$sex_group <- ifelse(
  dat$sex == 1,
  "Male",
  ifelse(
    dat$sex == 2,
    "Female",
    NA
  )
)


# ------------------------------------------------------------
# Age
# ------------------------------------------------------------

dat$age_group <- ifelse(
  is.na(dat$age),
  NA,
  ifelse(
    dat$age < 70,
    "<70",
    ">=70"
  )
)


# ------------------------------------------------------------
# LVEF
# ------------------------------------------------------------

dat$lvef_group <- ifelse(
  is.na(dat$LVEF),
  NA,
  ifelse(
    dat$LVEF < 50,
    "<50%",
    ">=50%"
  )
)


# ------------------------------------------------------------
# CKD
# ------------------------------------------------------------

dat$ckd_group <- ifelse(
  is.na(dat$CKD),
  NA,
  ifelse(
    dat$CKD == 1,
    "CKD",
    "No CKD"
  )
)


# ------------------------------------------------------------
# DM
# ------------------------------------------------------------

dat$dm_group <- ifelse(
  is.na(dat$DM),
  NA,
  ifelse(
    dat$DM == 1,
    "DM",
    "No DM"
  )
)


# ------------------------------------------------------------
# AF
# AF may be stored as character in source
# ------------------------------------------------------------

af_num <- suppressWarnings(
  as.numeric(
    as.character(
      dat$AF
    )
  )
)


dat$af_group <- ifelse(
  is.na(af_num),
  NA,
  ifelse(
    af_num == 1,
    "AF",
    "No AF"
  )
)


# ============================================================
# 5. Metric functions
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


  out <- as.numeric(
    fit$concordance
  )


  if (!is.finite(out)) {
    return(NA_real_)
  }


  out
}


# ============================================================
# 6. Define subgroup structure
# ============================================================

subgroup_definitions <- list(

  Sex = "sex_group",

  Age = "age_group",

  LVEF = "lvef_group",

  CKD = "ckd_group",

  DM = "dm_group",

  AF = "af_group"
)


# ============================================================
# 7. Event-count audit
# ============================================================

audit_rows <- list()

counter <- 0


for (sg_name in names(subgroup_definitions)) {

  sg_col <- subgroup_definitions[[sg_name]]

  levels_present <- unique(
    dat[[sg_col]]
  )

  levels_present <- levels_present[
    !is.na(levels_present)
  ]


  for (lv in levels_present) {

    d <- dat[
      dat[[sg_col]] == lv &
        !is.na(dat[[sg_col]]),
      ,
      drop = FALSE
    ]


    counter <- counter + 1


    audit_rows[[counter]] <- data.frame(

      subgroup =
        sg_name,

      level =
        lv,

      n =
        nrow(d),

      deaths =
        sum(
          d$status == 1,
          na.rm = TRUE
        ),

      alive_censored =
        sum(
          d$status == 0,
          na.rm = TRUE
        ),

      event_rate =
        mean(
          d$status == 1,
          na.rm = TRUE
        ),

      eligible_for_formal_analysis =
        sum(
          d$status == 1,
          na.rm = TRUE
        ) >= min_events,

      stringsAsFactors =
        FALSE
    )
  }
}


event_audit <- do.call(
  rbind,
  audit_rows
)


rownames(event_audit) <- NULL


cat("============================================\n")
cat("B. Subgroup event-count audit\n")
cat("============================================\n\n")

print(
  event_audit,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 8. Formal subgroup discrimination
#
# Only if events >= min_events
# ============================================================

perf_rows <- list()

counter <- 0


for (i in seq_len(nrow(event_audit))) {

  sg_name <- event_audit$subgroup[i]
  lv <- event_audit$level[i]


  if (
    !event_audit$eligible_for_formal_analysis[i]
  ) {
    next
  }


  sg_col <- subgroup_definitions[[sg_name]]


  d <- dat[
    dat[[sg_col]] == lv &
      !is.na(dat[[sg_col]]),
    ,
    drop = FALSE
  ]


  auc_base6 <- calc_auc(
    d$time,
    d$status,
    d$base6_marker,
    t0
  )


  auc_maggic <- calc_auc(
    d$time,
    d$status,
    d$maggic_score,
    t0
  )


  auc_gwtg <- calc_auc(
    d$time,
    d$status,
    d$gwtg_score,
    t0
  )


  c_base6 <- calc_cindex(
    d$time,
    d$status,
    d$base6_marker
  )


  c_maggic <- calc_cindex(
    d$time,
    d$status,
    d$maggic_score
  )


  c_gwtg <- calc_cindex(
    d$time,
    d$status,
    d$gwtg_score
  )


  counter <- counter + 1


  perf_rows[[counter]] <- data.frame(

    subgroup =
      sg_name,

    level =
      lv,

    n =
      nrow(d),

    deaths =
      sum(d$status == 1),

    AUC_Base6 =
      auc_base6,

    AUC_MAGGIC =
      auc_maggic,

    AUC_GWTG =
      auc_gwtg,

    delta_AUC_Base6_MAGGIC =
      auc_base6 -
      auc_maggic,

    delta_AUC_Base6_GWTG =
      auc_base6 -
      auc_gwtg,

    C_Base6 =
      c_base6,

    C_MAGGIC =
      c_maggic,

    C_GWTG =
      c_gwtg,

    delta_C_Base6_MAGGIC =
      c_base6 -
      c_maggic,

    delta_C_Base6_GWTG =
      c_base6 -
      c_gwtg,

    stringsAsFactors =
      FALSE
  )
}


if (length(perf_rows) > 0) {

  subgroup_performance <- do.call(
    rbind,
    perf_rows
  )

  rownames(
    subgroup_performance
  ) <- NULL

} else {

  subgroup_performance <- data.frame()
}


cat("============================================\n")
cat("C. Eligible subgroup discrimination\n")
cat("============================================\n\n")

print(
  subgroup_performance,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 9. Directional consistency summary
# ============================================================

if (nrow(subgroup_performance) > 0) {

  directional_summary <- data.frame(

    subgroup =
      subgroup_performance$subgroup,

    level =
      subgroup_performance$level,

    n =
      subgroup_performance$n,

    deaths =
      subgroup_performance$deaths,

    Base6_AUC_higher_than_MAGGIC =
      subgroup_performance$delta_AUC_Base6_MAGGIC > 0,

    Base6_AUC_higher_than_GWTG =
      subgroup_performance$delta_AUC_Base6_GWTG > 0,

    Base6_C_higher_than_MAGGIC =
      subgroup_performance$delta_C_Base6_MAGGIC > 0,

    Base6_C_higher_than_GWTG =
      subgroup_performance$delta_C_Base6_GWTG > 0,

    stringsAsFactors =
      FALSE
  )


  overall_direction <- data.frame(

    comparison = c(
      "Base6 AUC > MAGGIC",
      "Base6 AUC > GWTG-HF",
      "Base6 C > MAGGIC",
      "Base6 C > GWTG-HF"
    ),

    eligible_subgroups = rep(
      nrow(subgroup_performance),
      4
    ),

    positive_subgroups = c(
      sum(
        directional_summary$Base6_AUC_higher_than_MAGGIC,
        na.rm = TRUE
      ),

      sum(
        directional_summary$Base6_AUC_higher_than_GWTG,
        na.rm = TRUE
      ),

      sum(
        directional_summary$Base6_C_higher_than_MAGGIC,
        na.rm = TRUE
      ),

      sum(
        directional_summary$Base6_C_higher_than_GWTG,
        na.rm = TRUE
      )
    ),

    stringsAsFactors =
      FALSE
  )


  overall_direction$proportion_positive <-
    overall_direction$positive_subgroups /
    overall_direction$eligible_subgroups

} else {

  directional_summary <- data.frame()

  overall_direction <- data.frame()
}


cat("============================================\n")
cat("D. Directional consistency\n")
cat("============================================\n\n")

print(
  overall_direction,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 10. Missingness audit for subgroup variables
# ============================================================

missingness_audit <- data.frame(

  subgroup = names(subgroup_definitions),

  variable = unname(
    unlist(
      subgroup_definitions
    )
  ),

  missing_n = vapply(
    subgroup_definitions,
    function(x) {
      sum(
        is.na(
          dat[[x]]
        )
      )
    },
    numeric(1)
  ),

  stringsAsFactors =
    FALSE
)


missingness_audit$missing_pct <-
  100 *
  missingness_audit$missing_n /
  nrow(dat)


cat("============================================\n")
cat("E. Subgroup variable missingness\n")
cat("============================================\n\n")

print(
  missingness_audit,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 11. Optional descriptive flag
#
# No significance tests are performed here.
# ============================================================

if (nrow(subgroup_performance) > 0) {

  subgroup_performance$interpretation_flag <-
    ifelse(
      subgroup_performance$delta_C_Base6_MAGGIC > 0 &
        subgroup_performance$delta_C_Base6_GWTG > 0,
      "Base6 C higher than both",
      ifelse(
        subgroup_performance$delta_C_Base6_MAGGIC > 0 |
          subgroup_performance$delta_C_Base6_GWTG > 0,
        "Base6 C higher than one comparator",
        "No directional C advantage"
      )
    )
}


# ============================================================
# 12. Save outputs
# ============================================================

write.csv(
  event_audit,
  "Chinese_GWTG_STEP11P_SUBGROUP_EVENT_AUDIT.csv",
  row.names = FALSE
)


write.csv(
  subgroup_performance,
  "Chinese_GWTG_STEP11P_SUBGROUP_PERFORMANCE.csv",
  row.names = FALSE
)


write.csv(
  directional_summary,
  "Chinese_GWTG_STEP11P_DIRECTIONAL_DETAIL.csv",
  row.names = FALSE
)


write.csv(
  overall_direction,
  "Chinese_GWTG_STEP11P_DIRECTIONAL_SUMMARY.csv",
  row.names = FALSE
)


write.csv(
  missingness_audit,
  "Chinese_GWTG_STEP11P_SUBGROUP_MISSINGNESS.csv",
  row.names = FALSE
)


# ============================================================
# 13. Save RData
# ============================================================

save(
  dat,
  event_audit,
  subgroup_performance,
  directional_summary,
  overall_direction,
  missingness_audit,
  min_events,
  file =
    "Chinese_GWTG_STEP11P_SUBGROUP_ANALYSIS.RData"
)


# ============================================================
# 14. Completion
# ============================================================

cat("============================================\n")
cat("STEP 11P completed\n")
cat("============================================\n\n")

cat(
  "Rule used: subgroup performance calculated only when ",
  "deaths >= ",
  min_events,
  ".\n\n",
  sep = ""
)

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11P_SUBGROUP_EVENT_AUDIT.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11P_SUBGROUP_PERFORMANCE.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11P_DIRECTIONAL_DETAIL.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11P_DIRECTIONAL_SUMMARY.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11P_SUBGROUP_MISSINGNESS.csv\n"
)

cat(
  "6. Chinese_GWTG_STEP11P_SUBGROUP_ANALYSIS.RData\n"
)

cat("============================================\n")
