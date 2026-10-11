# ============================================================
# STEP 9B-1
# Fix external 1-year time-dependent AUC
#
# Absolute-risk horizon remains frozen at 1 year.
# Discrimination evaluated at day 364 because follow-up
# is administratively censored at day 365.
# ============================================================

rm(list = ls())

library(dplyr)
library(survival)
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

load(
  "Chinese_HF_EXTERNAL_STEP9B_VALIDATION.RData"
)

cat("============================================\n")
cat("STEP 9B-1: External AUC audit\n")
cat("============================================\n")


# ============================================================
# A. Follow-up boundary
# ============================================================

cat("\n============================================\n")
cat("A. Follow-up endpoint audit\n")
cat("============================================\n")

cat(
  "Maximum follow-up =",
  max(external_data$time_day),
  "days\n"
)

cat(
  "Alive/censored at day 365 =",
  sum(
    external_data$status == 0 &
      external_data$time_day == 365
  ),
  "\n"
)

cat(
  "Deaths exactly at day 365 =",
  sum(
    external_data$status == 1 &
      external_data$time_day == 365
  ),
  "\n"
)


# ============================================================
# B. Controls available for timeROC
# ============================================================

check_times <- c(
  360,
  364,
  365
)

control_check <- data.frame(

  time =
    check_times,

  cases =
    sapply(
      check_times,
      function(t) {
        sum(
          external_data$status == 1 &
            external_data$time_day <= t
        )
      }
    ),

  controls =
    sapply(
      check_times,
      function(t) {
        sum(
          external_data$time_day > t
        )
      }
    ),

  censored_before =
    sapply(
      check_times,
      function(t) {
        sum(
          external_data$status == 0 &
            external_data$time_day <= t
        )
      }
    )
)


cat("\n============================================\n")
cat("B. Cases and controls by time\n")
cat("============================================\n")

print(
  control_check,
  row.names = FALSE
)


# ============================================================
# C. Time-dependent AUC at 360 and 364 days
# ============================================================

roc_fix <- timeROC::timeROC(

  T =
    external_data$time_day,

  delta =
    external_data$status,

  marker =
    external_data$lp_base6,

  cause =
    1,

  weighting =
    "marginal",

  times =
    c(
      360,
      364
    ),

  iid =
    TRUE
)


auc_fix <- data.frame(

  time_day =
    c(
      360,
      364
    ),

  AUC =
    as.numeric(
      roc_fix$AUC
    )
)


cat("\n============================================\n")
cat("C. Near-1-year time-dependent AUC\n")
cat("============================================\n")

print(
  auc_fix,
  row.names = FALSE
)


# ============================================================
# D. 95% CI
# ============================================================

auc_ci <- tryCatch(

  confint(
    roc_fix
  ),

  error =
    function(e) NULL
)


cat("\n============================================\n")
cat("D. AUC confidence intervals\n")
cat("============================================\n")

if (
  is.null(
    auc_ci
  )
) {

  cat(
    "confint() unavailable.\n"
  )

} else {

  print(
    auc_ci
  )
}


# ============================================================
# E. Harrell C as supplementary discrimination
# ============================================================

c_obj <- survival::concordance(

  Surv(
    time_day,
    status
  ) ~
    lp_base6,

  data =
    external_data,

  reverse =
    TRUE
)


harrell_C <-
  as.numeric(
    c_obj$concordance
  )


harrell_SE <-
  sqrt(
    c_obj$var
  )


cat("\n============================================\n")
cat("E. External Harrell C-index\n")
cat("============================================\n")

cat(
  "C-index =",
  harrell_C,
  "\n"
)

cat(
  "Approximate 95% CI =",
  harrell_C -
    1.96 *
    harrell_SE,
  "to",
  harrell_C +
    1.96 *
    harrell_SE,
  "\n"
)


# ============================================================
# F. Save
# ============================================================

write.csv(

  auc_fix,

  "Chinese_Base6_AUC_360_364days.csv",

  row.names = FALSE
)


save(

  roc_fix,

  auc_fix,

  auc_ci,

  harrell_C,

  harrell_SE,

  control_check,

  file =
    "Chinese_HF_EXTERNAL_AUC_FIX.RData"
)


cat("\n============================================\n")
cat("STEP 9B-1 完成\n")
cat("============================================\n")
