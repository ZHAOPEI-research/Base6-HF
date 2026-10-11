rm(list = ls())
gc()

library(survival)
library(timeROC)
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
# 1. Load original STEP11P results
# ============================================================

e <- new.env()

load(
  "Chinese_GWTG_STEP11P_SUBGROUP_ANALYSIS.RData",
  envir = e
)

dat <- e$dat
subgroup_performance <- e$subgroup_performance
min_events <- e$min_events

cat("============================================\n")
cat("STEP 52: Subgroup discrimination 95% CI\n")
cat("============================================\n")

cat("External cohort n =", nrow(dat), "\n")
cat("Deaths =", sum(dat$status), "\n")
cat("Minimum subgroup deaths =", min_events, "\n")


# ============================================================
# 2. Subgroup-variable map
# ============================================================

subgroup_var_map <- c(
  Sex  = "sex_group",
  Age  = "age_group",
  LVEF = "lvef_group",
  CKD  = "ckd_group",
  DM   = "dm_group",
  AF   = "af_group"
)


# ============================================================
# 3. Harrell C
# ============================================================

get_cindex <- function(
    time,
    status,
    marker
) {

  if (
    length(unique(status)) < 2 ||
    sum(status == 1) == 0
  ) {
    return(NA_real_)
  }

  z <- tryCatch(
    survival::concordance(
      Surv(time, status) ~ marker,
      reverse = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(z)) {
    return(NA_real_)
  }

  out <- as.numeric(z$concordance)

  if (!is.finite(out)) {
    return(NA_real_)
  }

  out
}


# ============================================================
# 4. 364-day AUC
# ============================================================

get_auc364 <- function(
    time,
    status,
    marker
) {

  # Need deaths by 364 days and survivors beyond 364 days
  n_case <- sum(
    status == 1 &
      time <= 364
  )

  n_control <- sum(
    time > 364
  )

  if (
    n_case == 0 ||
    n_control == 0
  ) {
    return(NA_real_)
  }

  z <- tryCatch(
    timeROC(
      T = time,
      delta = status,
      marker = marker,
      cause = 1,
      weighting = "marginal",
      times = 364,
      iid = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(z)) {
    return(NA_real_)
  }

  # timeROC may include t = 0 automatically
  auc <- as.numeric(
    tail(
      z$AUC,
      1
    )
  )

  if (!is.finite(auc)) {
    return(NA_real_)
  }

  auc
}


# ============================================================
# 5. Model markers
# ============================================================

model_markers <- c(
  Base6  = "base6_marker",
  MAGGIC = "maggic_score",
  GWTG   = "gwtg_score"
)


# ============================================================
# 6. Recalculate original subgroup point estimates
#    before bootstrap
# ============================================================

point_results <- list()

counter <- 1

for (r in seq_len(nrow(subgroup_performance))) {

  sg <- subgroup_performance$subgroup[r]
  lv <- subgroup_performance$level[r]

  sg_var <- subgroup_var_map[[sg]]

  dsg <- dat[
    dat[[sg_var]] == lv,
    ,
    drop = FALSE
  ]

  row_out <- data.frame(
    subgroup = sg,
    level = lv,
    n = nrow(dsg),
    deaths = sum(dsg$status),
    stringsAsFactors = FALSE
  )

  for (model in names(model_markers)) {

    marker_var <- model_markers[[model]]

    row_out[[paste0("AUC_", model)]] <-
      get_auc364(
        time = dsg$time,
        status = dsg$status,
        marker = dsg[[marker_var]]
      )

    row_out[[paste0("C_", model)]] <-
      get_cindex(
        time = dsg$time,
        status = dsg$status,
        marker = dsg[[marker_var]]
      )
  }

  point_results[[counter]] <- row_out

  counter <- counter + 1
}

point_results <- bind_rows(point_results)


cat("\n============================================\n")
cat("A. Recalculated point estimates\n")
cat("============================================\n")

print(
  point_results,
  row.names = FALSE
)


# ============================================================
# 7. Audit against original STEP11P values
# ============================================================

audit <- point_results %>%
  left_join(
    subgroup_performance %>%
      select(
        subgroup,
        level,
        AUC_Base6,
        AUC_MAGGIC,
        AUC_GWTG,
        C_Base6,
        C_MAGGIC,
        C_GWTG
      ),
    by = c(
      "subgroup",
      "level"
    ),
    suffix = c(
      "_new",
      "_original"
    )
  ) %>%
  mutate(
    diff_AUC_Base6 =
      AUC_Base6_new - AUC_Base6_original,

    diff_AUC_MAGGIC =
      AUC_MAGGIC_new - AUC_MAGGIC_original,

    diff_AUC_GWTG =
      AUC_GWTG_new - AUC_GWTG_original,

    diff_C_Base6 =
      C_Base6_new - C_Base6_original,

    diff_C_MAGGIC =
      C_MAGGIC_new - C_MAGGIC_original,

    diff_C_GWTG =
      C_GWTG_new - C_GWTG_original
  )


cat("\n============================================\n")
cat("B. Point-estimate audit\n")
cat("============================================\n")

print(
  audit %>%
    select(
      subgroup,
      level,
      starts_with("diff_")
    ),
  row.names = FALSE
)


# ============================================================
# 8. Bootstrap settings
# ============================================================

B <- 2000

set.seed(20261004)

boot_results <- list()

counter <- 1


# ============================================================
# 9. Bootstrap within each subgroup
# ============================================================

for (r in seq_len(nrow(subgroup_performance))) {

  sg <- subgroup_performance$subgroup[r]
  lv <- subgroup_performance$level[r]

  sg_var <- subgroup_var_map[[sg]]

  dsg <- dat[
    dat[[sg_var]] == lv,
    ,
    drop = FALSE
  ]

  n_sg <- nrow(dsg)
  deaths_sg <- sum(dsg$status)

  cat(
    "\n--------------------------------------------\n"
  )

  cat(
    "Subgroup:",
    sg,
    "/",
    lv,
    "| n =",
    n_sg,
    "| deaths =",
    deaths_sg,
    "\n"
  )

  # safety rule
  if (deaths_sg < min_events) {

    cat(
      "Skipped: fewer than",
      min_events,
      "deaths\n"
    )

    next
  }

  for (b in seq_len(B)) {

    # Same patient resample used for all 3 models
    ix <- sample.int(
      n_sg,
      size = n_sg,
      replace = TRUE
    )

    db <- dsg[
      ix,
      ,
      drop = FALSE
    ]

    metrics <- list()

    for (model in names(model_markers)) {

      marker_var <- model_markers[[model]]

      metrics[[paste0("AUC_", model)]] <-
        get_auc364(
          time = db$time,
          status = db$status,
          marker = db[[marker_var]]
        )

      metrics[[paste0("C_", model)]] <-
        get_cindex(
          time = db$time,
          status = db$status,
          marker = db[[marker_var]]
        )
    }

    boot_results[[counter]] <-
      data.frame(
        subgroup = sg,
        level = lv,
        bootstrap = b,

        AUC_Base6 =
          metrics$AUC_Base6,

        AUC_MAGGIC =
          metrics$AUC_MAGGIC,

        AUC_GWTG =
          metrics$AUC_GWTG,

        C_Base6 =
          metrics$C_Base6,

        C_MAGGIC =
          metrics$C_MAGGIC,

        C_GWTG =
          metrics$C_GWTG,

        stringsAsFactors = FALSE
      )

    counter <- counter + 1

    if (b %% 500 == 0) {

      cat(
        "Bootstrap",
        b,
        "/",
        B,
        "\n"
      )
    }
  }
}


boot_results <- bind_rows(
  boot_results
)


# ============================================================
# 10. Add paired differences
# ============================================================

boot_results <- boot_results %>%
  mutate(

    delta_AUC_Base6_MAGGIC =
      AUC_Base6 - AUC_MAGGIC,

    delta_AUC_Base6_GWTG =
      AUC_Base6 - AUC_GWTG,

    delta_C_Base6_MAGGIC =
      C_Base6 - C_MAGGIC,

    delta_C_Base6_GWTG =
      C_Base6 - C_GWTG
  )


# ============================================================
# 11. Helper for percentile CI
# ============================================================

q025 <- function(x) {

  quantile(
    x,
    probs = 0.025,
    na.rm = TRUE,
    names = FALSE
  )
}


q975 <- function(x) {

  quantile(
    x,
    probs = 0.975,
    na.rm = TRUE,
    names = FALSE
  )
}


# ============================================================
# 12. Metric-specific CIs
# ============================================================

ci_summary <- boot_results %>%
  group_by(
    subgroup,
    level
  ) %>%
  summarise(

    AUC_Base6_low =
      q025(AUC_Base6),

    AUC_Base6_high =
      q975(AUC_Base6),

    AUC_MAGGIC_low =
      q025(AUC_MAGGIC),

    AUC_MAGGIC_high =
      q975(AUC_MAGGIC),

    AUC_GWTG_low =
      q025(AUC_GWTG),

    AUC_GWTG_high =
      q975(AUC_GWTG),


    C_Base6_low =
      q025(C_Base6),

    C_Base6_high =
      q975(C_Base6),

    C_MAGGIC_low =
      q025(C_MAGGIC),

    C_MAGGIC_high =
      q975(C_MAGGIC),

    C_GWTG_low =
      q025(C_GWTG),

    C_GWTG_high =
      q975(C_GWTG),


    valid_AUC_Base6 =
      sum(is.finite(AUC_Base6)),

    valid_AUC_MAGGIC =
      sum(is.finite(AUC_MAGGIC)),

    valid_AUC_GWTG =
      sum(is.finite(AUC_GWTG)),

    valid_C_Base6 =
      sum(is.finite(C_Base6)),

    valid_C_MAGGIC =
      sum(is.finite(C_MAGGIC)),

    valid_C_GWTG =
      sum(is.finite(C_GWTG)),

    .groups = "drop"
  )


# ============================================================
# 13. Paired difference CIs
#
# Exploratory only; no subgroup P values are generated.
# ============================================================

difference_CI <- boot_results %>%
  group_by(
    subgroup,
    level
  ) %>%
  summarise(

    delta_AUC_Base6_MAGGIC_low =
      q025(delta_AUC_Base6_MAGGIC),

    delta_AUC_Base6_MAGGIC_high =
      q975(delta_AUC_Base6_MAGGIC),

    delta_AUC_Base6_GWTG_low =
      q025(delta_AUC_Base6_GWTG),

    delta_AUC_Base6_GWTG_high =
      q975(delta_AUC_Base6_GWTG),

    delta_C_Base6_MAGGIC_low =
      q025(delta_C_Base6_MAGGIC),

    delta_C_Base6_MAGGIC_high =
      q975(delta_C_Base6_MAGGIC),

    delta_C_Base6_GWTG_low =
      q025(delta_C_Base6_GWTG),

    delta_C_Base6_GWTG_high =
      q975(delta_C_Base6_GWTG),

    .groups = "drop"
  )


# ============================================================
# 14. Merge point estimates + CIs
# ============================================================

final_subgroup_CI <- subgroup_performance %>%
  left_join(
    ci_summary,
    by = c(
      "subgroup",
      "level"
    )
  ) %>%
  left_join(
    difference_CI,
    by = c(
      "subgroup",
      "level"
    )
  )


# ============================================================
# 15. Manuscript-ready C-index table
# ============================================================

C_table <- final_subgroup_CI %>%
  transmute(

    Subgroup =
      paste(
        subgroup,
        level,
        sep = ": "
      ),

    `n / deaths` =
      paste0(
        n,
        " / ",
        deaths
      ),

    Base6 =
      sprintf(
        "%.3f (%.3f–%.3f)",
        C_Base6,
        C_Base6_low,
        C_Base6_high
      ),

    MAGGIC =
      sprintf(
        "%.3f (%.3f–%.3f)",
        C_MAGGIC,
        C_MAGGIC_low,
        C_MAGGIC_high
      ),

    `GWTG-HF` =
      sprintf(
        "%.3f (%.3f–%.3f)",
        C_GWTG,
        C_GWTG_low,
        C_GWTG_high
      )
  )


# ============================================================
# 16. Manuscript-ready AUC table
# ============================================================

AUC_table <- final_subgroup_CI %>%
  transmute(

    Subgroup =
      paste(
        subgroup,
        level,
        sep = ": "
      ),

    `n / deaths` =
      paste0(
        n,
        " / ",
        deaths
      ),

    Base6 =
      sprintf(
        "%.3f (%.3f–%.3f)",
        AUC_Base6,
        AUC_Base6_low,
        AUC_Base6_high
      ),

    MAGGIC =
      sprintf(
        "%.3f (%.3f–%.3f)",
        AUC_MAGGIC,
        AUC_MAGGIC_low,
        AUC_MAGGIC_high
      ),

    `GWTG-HF` =
      sprintf(
        "%.3f (%.3f–%.3f)",
        AUC_GWTG,
        AUC_GWTG_low,
        AUC_GWTG_high
      )
  )


# ============================================================
# 17. Print
# ============================================================

cat("\n============================================\n")
cat("C. Subgroup Harrell C with 95% CI\n")
cat("============================================\n")

print(
  C_table,
  row.names = FALSE
)


cat("\n============================================\n")
cat("D. Subgroup 364-day AUC with 95% CI\n")
cat("============================================\n")

print(
  AUC_table,
  row.names = FALSE
)


cat("\n============================================\n")
cat("E. Valid bootstrap replicates\n")
cat("============================================\n")

print(
  ci_summary %>%
    select(
      subgroup,
      level,
      starts_with("valid_")
    ),
  row.names = FALSE
)


# ============================================================
# 18. Save
# ============================================================

write.csv(
  final_subgroup_CI,
  "Chinese_GWTG_STEP11P_SUBGROUP_CI_FULL.csv",
  row.names = FALSE
)

write.csv(
  C_table,
  "Chinese_GWTG_STEP11P_SUBGROUP_CI_CINDEX.csv",
  row.names = FALSE
)

write.csv(
  AUC_table,
  "Chinese_GWTG_STEP11P_SUBGROUP_CI_AUC364.csv",
  row.names = FALSE
)

write.csv(
  difference_CI,
  "Chinese_GWTG_STEP11P_SUBGROUP_DIFFERENCE_CI.csv",
  row.names = FALSE
)

write.csv(
  audit,
  "Chinese_GWTG_STEP11P_SUBGROUP_CI_AUDIT.csv",
  row.names = FALSE
)

save(
  boot_results,
  final_subgroup_CI,
  C_table,
  AUC_table,
  difference_CI,
  audit,
  file =
    "Chinese_GWTG_STEP11P_SUBGROUP_CI_2000BOOT.RData"
)

cat("\nSTEP 52 completed\n")
