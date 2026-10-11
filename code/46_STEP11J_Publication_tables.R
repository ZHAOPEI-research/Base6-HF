# ============================================================
# STEP 11J
# Publication-ready model comparison tables
#
# This step DOES NOT recalculate model performance.
# It consolidates validated results from STEP 11C-11I.
#
# Inputs:
#   Chinese_GWTG_STEP11C_PAIRED_BOOTSTRAP.RData
#   Chinese_GWTG_STEP11D_CALIBRATION.RData
#   Chinese_GWTG_STEP11F_RECALIBRATION.RData
#   Chinese_GWTG_STEP11H_DCA_PAIRED_BOOTSTRAP.RData
#   Chinese_GWTG_STEP11I_RECALIBRATION_OPTIMISM_FIXED.RData
#
# Outputs:
#   Publication-ready CSV tables
#   Master RData
# ============================================================


# ============================================================
# 0. Clean environment
# ============================================================

rm(list = ls())


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

cat("\n")
cat("============================================\n")
cat("STEP 11J: Publication-ready result tables\n")
cat("============================================\n\n")


# ============================================================
# 1. Helper functions
# ============================================================

fmt_num <- function(x, digits = 3) {

  ifelse(
    is.na(x),
    "",
    sprintf(
      paste0("%.", digits, "f"),
      x
    )
  )
}


fmt_ci <- function(
    estimate,
    low,
    high,
    digits = 3
) {

  ifelse(
    is.na(estimate) |
      is.na(low) |
      is.na(high),
    "",
    paste0(
      sprintf(
        paste0("%.", digits, "f"),
        estimate
      ),
      " (",
      sprintf(
        paste0("%.", digits, "f"),
        low
      ),
      "–",
      sprintf(
        paste0("%.", digits, "f"),
        high
      ),
      ")"
    )
  )
}


fmt_p <- function(p) {

  ifelse(
    is.na(p),
    "",
    ifelse(
      p < 0.001,
      "<0.001",
      sprintf(
        "%.3f",
        p
      )
    )
  )
}


check_file <- function(f) {

  if (!file.exists(f)) {

    stop(
      paste0(
        "Required file not found: ",
        f
      )
    )
  }
}


# ============================================================
# 2. Define input files
# ============================================================

file_11c <-
  "Chinese_GWTG_STEP11C_PAIRED_BOOTSTRAP.RData"

file_11d <-
  "Chinese_GWTG_STEP11D_CALIBRATION.RData"

file_11f <-
  "Chinese_GWTG_STEP11F_RECALIBRATION.RData"

file_11h <-
  "Chinese_GWTG_STEP11H_DCA_PAIRED_BOOTSTRAP.RData"

file_11i <-
  "Chinese_GWTG_STEP11I_RECALIBRATION_OPTIMISM_FIXED.RData"


input_files <- c(
  file_11c,
  file_11d,
  file_11f,
  file_11h,
  file_11i
)


for (f in input_files) {
  check_file(f)
}


cat("All required input files found.\n\n")


# ============================================================
# 3. Load each RData into a separate environment
# ============================================================

e11c <- new.env()
e11d <- new.env()
e11f <- new.env()
e11h <- new.env()
e11i <- new.env()


load(
  file_11c,
  envir = e11c
)

load(
  file_11d,
  envir = e11d
)

load(
  file_11f,
  envir = e11f
)

load(
  file_11h,
  envir = e11h
)

load(
  file_11i,
  envir = e11i
)


# ============================================================
# 4. Required-object audit
# ============================================================

required_objects <- list(

  STEP11C = c(
    "final_model_summary",
    "final_pairwise_summary"
  ),

  STEP11D = c(
    "final_summary"
  ),

  STEP11F = c(
    "final_recal_summary",
    "recalibration_parameters"
  ),

  STEP11H = c(
    "dca_pairwise",
    "dca_pairwise_selected"
  ),

  STEP11I = c(
    "optimism_summary",
    "beta_boot_summary_fixed"
  )
)


env_map <- list(
  STEP11C = e11c,
  STEP11D = e11d,
  STEP11F = e11f,
  STEP11H = e11h,
  STEP11I = e11i
)


for (step_name in names(required_objects)) {

  env <- env_map[[step_name]]

  needed <- required_objects[[step_name]]

  missing <- needed[
    !vapply(
      needed,
      exists,
      logical(1),
      envir = env
    )
  ]


  if (length(missing) > 0) {

    stop(
      paste0(
        step_name,
        " missing objects: ",
        paste(
          missing,
          collapse = ", "
        )
      )
    )
  }
}


cat("All required result objects found.\n\n")


# ============================================================
# 5. Retrieve objects
# ============================================================

perf <-
  get(
    "final_model_summary",
    envir = e11c
  )


pairwise <-
  get(
    "final_pairwise_summary",
    envir = e11c
  )


calibration <-
  get(
    "final_summary",
    envir = e11d
  )


recal <-
  get(
    "final_recal_summary",
    envir = e11f
  )


recal_params <-
  get(
    "recalibration_parameters",
    envir = e11f
  )


dca_pairwise <-
  get(
    "dca_pairwise",
    envir = e11h
  )


dca_selected <-
  get(
    "dca_pairwise_selected",
    envir = e11h
  )


optimism <-
  get(
    "optimism_summary",
    envir = e11i
  )


slope_boot <-
  get(
    "beta_boot_summary_fixed",
    envir = e11i
  )


# ============================================================
# 6. TABLE 1
# Three-model discrimination
# ============================================================

table1_discrimination_raw <- data.frame(

  Model =
    perf$model,

  N =
    perf$n,

  Deaths =
    perf$deaths,

  AUC_364 =
    perf$AUC_364,

  AUC_CI_low =
    perf$AUC_CI_low,

  AUC_CI_high =
    perf$AUC_CI_high,

  Harrell_C =
    perf$Harrell_C,

  C_CI_low =
    perf$C_CI_low,

  C_CI_high =
    perf$C_CI_high,

  stringsAsFactors = FALSE
)


table1_discrimination <- data.frame(

  Model =
    perf$model,

  N =
    perf$n,

  Deaths =
    perf$deaths,

  `364-day AUC (95% CI)` =
    fmt_ci(
      perf$AUC_364,
      perf$AUC_CI_low,
      perf$AUC_CI_high
    ),

  `Harrell C-index (95% CI)` =
    fmt_ci(
      perf$Harrell_C,
      perf$C_CI_low,
      perf$C_CI_high
    ),

  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("TABLE 1. Three-model discrimination\n")
cat("============================================\n\n")

print(
  table1_discrimination,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 7. TABLE 2
# Pairwise discrimination comparison
# ============================================================

table2_pairwise_raw <-
  pairwise


table2_pairwise <- data.frame(

  Comparison =
    pairwise$comparison,

  Metric =
    pairwise$metric,

  Difference =
    fmt_num(
      pairwise$estimate,
      3
    ),

  `95% CI` =
    ifelse(
      is.na(pairwise$CI_low),
      "",
      paste0(
        fmt_num(
          pairwise$CI_low,
          3
        ),
        "–",
        fmt_num(
          pairwise$CI_high,
          3
        )
      )
    ),

  `Bootstrap P` =
    fmt_p(
      pairwise$p_boot
    ),

  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("TABLE 2. Pairwise discrimination\n")
cat("============================================\n\n")

print(
  table2_pairwise,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 8. TABLE 3
# Original probability calibration
#
# Base6 vs MAGGIC only
# ============================================================

calibration <- calibration[
  calibration$model %in%
    c(
      "Base6",
      "MAGGIC"
    ),
  ,
  drop = FALSE
]


calibration <- calibration[
  match(
    c(
      "Base6",
      "MAGGIC"
    ),
    calibration$model
  ),
  ,
  drop = FALSE
]


table3_calibration_raw <- calibration


table3_calibration <- data.frame(

  Model =
    calibration$model,

  N =
    calibration$n,

  Deaths =
    calibration$deaths,

  `Observed KM risk` =
    fmt_num(
      calibration$observed_KM_364,
      3
    ),

  `Mean predicted risk` =
    fmt_num(
      calibration$mean_predicted_364,
      3
    ),

  `O/P ratio` =
    fmt_num(
      calibration$O_P,
      3
    ),

  `Brier score` =
    fmt_num(
      calibration$Brier_364,
      4
    ),

  `Scaled Brier` =
    fmt_num(
      calibration$Scaled_Brier,
      3
    ),

  `Quintile MAE` =
    fmt_num(
      calibration$calibration_MAE_quintile,
      3
    ),

  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("TABLE 3. Original model calibration\n")
cat("============================================\n\n")

print(
  table3_calibration,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 9. Extract optimism-corrected recalibration metrics
# ============================================================

get_corrected_metric <- function(
    model_name,
    metric_name
) {

  x <- optimism$corrected[
    optimism$model == model_name &
      optimism$metric == metric_name
  ]


  if (length(x) != 1) {
    return(NA_real_)
  }


  as.numeric(x)
}


intercept_corrected_OP <-
  get_corrected_metric(
    "Base6 recalibrated: intercept only",
    "O_P"
  )


intercept_corrected_Brier <-
  get_corrected_metric(
    "Base6 recalibrated: intercept only",
    "Brier_364"
  )


full_corrected_OP <-
  get_corrected_metric(
    "Base6 recalibrated: intercept + slope",
    "O_P"
  )


full_corrected_Brier <-
  get_corrected_metric(
    "Base6 recalibrated: intercept + slope",
    "Brier_364"
  )


# ============================================================
# 10. TABLE 4
# Base6 post-hoc recalibration
# ============================================================

desired_recal_order <- c(
  "Base6 original",
  "Base6 recalibrated: intercept only",
  "Base6 recalibrated: intercept + slope"
)


recal_base6 <- recal[
  recal$model %in%
    desired_recal_order,
  ,
  drop = FALSE
]


recal_base6 <- recal_base6[
  match(
    desired_recal_order,
    recal_base6$model
  ),
  ,
  drop = FALSE
]


recal_base6$optimism_corrected_OP <-
  c(
    NA_real_,
    intercept_corrected_OP,
    full_corrected_OP
  )


recal_base6$optimism_corrected_Brier <-
  c(
    NA_real_,
    intercept_corrected_Brier,
    full_corrected_Brier
  )


table4_recalibration_raw <- recal_base6


table4_recalibration <- data.frame(

  Model =
    recal_base6$model,

  `Mean predicted risk` =
    fmt_num(
      recal_base6$mean_predicted_364,
      3
    ),

  `Observed KM risk` =
    fmt_num(
      recal_base6$observed_KM_364,
      3
    ),

  `Apparent O/P` =
    fmt_num(
      recal_base6$O_P,
      3
    ),

  `Optimism-corrected O/P` =
    fmt_num(
      recal_base6$optimism_corrected_OP,
      3
    ),

  `Apparent Brier` =
    fmt_num(
      recal_base6$Brier_364,
      4
    ),

  `Optimism-corrected Brier` =
    fmt_num(
      recal_base6$optimism_corrected_Brier,
      4
    ),

  `Quintile MAE (apparent)` =
    fmt_num(
      recal_base6$calibration_MAE_quintile,
      3
    ),

  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("TABLE 4. Base6 recalibration\n")
cat("============================================\n\n")

print(
  table4_recalibration,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 11. TABLE 4B
# Recalibration parameters
# ============================================================

table4b_recalibration_parameters <- data.frame(

  Strategy =
    recal_params$model,

  `cloglog intercept` =
    fmt_num(
      recal_params$intercept_cloglog,
      3
    ),

  Slope =
    fmt_num(
      recal_params$slope,
      3
    ),

  `Slope 95% CI` =
    ifelse(
      is.na(
        recal_params$slope_CI_low
      ),
      "",
      paste0(
        fmt_num(
          recal_params$slope_CI_low,
          3
        ),
        "–",
        fmt_num(
          recal_params$slope_CI_high,
          3
        )
      )
    ),

  `Bootstrap slope 95% CI` =
    c(
      "",
      paste0(
        fmt_num(
          slope_boot$CI_low,
          3
        ),
        "–",
        fmt_num(
          slope_boot$CI_high,
          3
        )
      )
    ),

  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ============================================================
# 12. TABLE 5
# Main DCA:
# Base6 original vs MAGGIC
# ============================================================

selected_thresholds <- c(
  0.05,
  0.10,
  0.15,
  0.20,
  0.25,
  0.30
)


is_selected_threshold <- function(x) {

  vapply(
    x,
    function(xx) {

      any(
        abs(
          xx -
            selected_thresholds
        ) < 1e-8
      )
    },
    logical(1)
  )
}


dca_main <- dca_pairwise[
  dca_pairwise$comparison ==
    "Base6 original - MAGGIC" &
    is_selected_threshold(
      dca_pairwise$threshold
    ),
  ,
  drop = FALSE
]


dca_main <- dca_main[
  order(
    dca_main$threshold
  ),
  ,
  drop = FALSE
]


table5_dca_main_raw <- dca_main


table5_dca_main <- data.frame(

  `Threshold probability` =
    paste0(
      round(
        dca_main$threshold *
          100
      ),
      "%"
    ),

  `Δ Net benefit (Base6 - MAGGIC)` =
    fmt_num(
      dca_main$delta_NB,
      4
    ),

  `95% CI` =
    paste0(
      fmt_num(
        dca_main$CI_low,
        4
      ),
      "–",
      fmt_num(
        dca_main$CI_high,
        4
      )
    ),

  `Bootstrap P` =
    fmt_p(
      dca_main$p_boot
    ),

  `P(Base6 > MAGGIC)` =
    fmt_num(
      dca_main$probability_A_better,
      3
    ),

  check.names = FALSE,
  stringsAsFactors = FALSE
)


cat("============================================\n")
cat("TABLE 5. DCA main pairwise comparison\n")
cat("============================================\n\n")

print(
  table5_dca_main,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 13. TABLE 6
# Post-hoc DCA sensitivity:
# recalibrated Base6 vs MAGGIC
# ============================================================

dca_recal_sensitivity <- dca_pairwise[
  dca_pairwise$comparison %in%
    c(
      "Base6 recalibrated intercept - MAGGIC",
      "Base6 recalibrated full - MAGGIC"
    ) &
    is_selected_threshold(
      dca_pairwise$threshold
    ),
  ,
  drop = FALSE
]


dca_recal_sensitivity <-
  dca_recal_sensitivity[
    order(
      dca_recal_sensitivity$comparison,
      dca_recal_sensitivity$threshold
    ),
    ,
    drop = FALSE
  ]


table6_dca_sensitivity_raw <-
  dca_recal_sensitivity


table6_dca_sensitivity <- data.frame(

  Comparison =
    dca_recal_sensitivity$comparison,

  Threshold =
    paste0(
      round(
        dca_recal_sensitivity$threshold *
          100
      ),
      "%"
    ),

  `Δ Net benefit` =
    fmt_num(
      dca_recal_sensitivity$delta_NB,
      4
    ),

  `95% CI` =
    paste0(
      fmt_num(
        dca_recal_sensitivity$CI_low,
        4
      ),
      "–",
      fmt_num(
        dca_recal_sensitivity$CI_high,
        4
      )
    ),

  `Bootstrap P` =
    fmt_p(
      dca_recal_sensitivity$p_boot
    ),

  check.names = FALSE,
  stringsAsFactors = FALSE
)


# ============================================================
# 14. MASTER summary
#
# This intentionally contains only metrics directly comparable
# across the three original models.
#
# Calibration fields for GWTG-HF are left unavailable because
# gwtg_score is not a directly calibrated 1-year probability.
# ============================================================

master_summary <- merge(
  table1_discrimination_raw,
  table3_calibration_raw[
    ,
    c(
      "model",
      "mean_predicted_364",
      "observed_KM_364",
      "O_P",
      "Brier_364",
      "Scaled_Brier",
      "calibration_MAE_quintile"
    )
  ],
  by.x =
    "Model",
  by.y =
    "model",
  all.x =
    TRUE,
  sort =
    FALSE
)


master_summary <- master_summary[
  match(
    c(
      "Base6",
      "MAGGIC",
      "GWTG-HF"
    ),
    master_summary$Model
  ),
  ,
  drop = FALSE
]


cat("============================================\n")
cat("MASTER SUMMARY\n")
cat("============================================\n\n")

print(
  master_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 15. Save RAW tables
# ============================================================

write.csv(
  table1_discrimination_raw,
  "Chinese_GWTG_STEP11J_TABLE1_DISCRIMINATION_RAW.csv",
  row.names = FALSE
)

write.csv(
  table2_pairwise_raw,
  "Chinese_GWTG_STEP11J_TABLE2_PAIRWISE_RAW.csv",
  row.names = FALSE
)

write.csv(
  table3_calibration_raw,
  "Chinese_GWTG_STEP11J_TABLE3_CALIBRATION_RAW.csv",
  row.names = FALSE
)

write.csv(
  table4_recalibration_raw,
  "Chinese_GWTG_STEP11J_TABLE4_RECALIBRATION_RAW.csv",
  row.names = FALSE
)

write.csv(
  table5_dca_main_raw,
  "Chinese_GWTG_STEP11J_TABLE5_DCA_MAIN_RAW.csv",
  row.names = FALSE
)

write.csv(
  master_summary,
  "Chinese_GWTG_STEP11J_MASTER_SUMMARY_RAW.csv",
  row.names = FALSE
)


# ============================================================
# 16. Save publication-ready tables
# ============================================================

write.csv(
  table1_discrimination,
  "Chinese_GWTG_STEP11J_TABLE1_DISCRIMINATION.csv",
  row.names = FALSE
)

write.csv(
  table2_pairwise,
  "Chinese_GWTG_STEP11J_TABLE2_PAIRWISE.csv",
  row.names = FALSE
)

write.csv(
  table3_calibration,
  "Chinese_GWTG_STEP11J_TABLE3_CALIBRATION.csv",
  row.names = FALSE
)

write.csv(
  table4_recalibration,
  "Chinese_GWTG_STEP11J_TABLE4_RECALIBRATION.csv",
  row.names = FALSE
)

write.csv(
  table4b_recalibration_parameters,
  "Chinese_GWTG_STEP11J_TABLE4B_RECALIBRATION_PARAMETERS.csv",
  row.names = FALSE
)

write.csv(
  table5_dca_main,
  "Chinese_GWTG_STEP11J_TABLE5_DCA_MAIN.csv",
  row.names = FALSE
)

write.csv(
  table6_dca_sensitivity,
  "Chinese_GWTG_STEP11J_TABLE6_DCA_RECAL_SENSITIVITY.csv",
  row.names = FALSE
)


# ============================================================
# 17. Save complete RData
# ============================================================

save(
  table1_discrimination_raw,
  table1_discrimination,
  table2_pairwise_raw,
  table2_pairwise,
  table3_calibration_raw,
  table3_calibration,
  table4_recalibration_raw,
  table4_recalibration,
  table4b_recalibration_parameters,
  table5_dca_main_raw,
  table5_dca_main,
  table6_dca_sensitivity_raw,
  table6_dca_sensitivity,
  master_summary,
  file =
    "Chinese_GWTG_STEP11J_PUBLICATION_TABLES.RData"
)


# ============================================================
# 18. Completion
# ============================================================

cat("============================================\n")
cat("STEP 11J completed\n")
cat("============================================\n\n")

cat("Publication-ready outputs:\n")

cat(
  "1. TABLE1_DISCRIMINATION\n"
)

cat(
  "2. TABLE2_PAIRWISE\n"
)

cat(
  "3. TABLE3_CALIBRATION\n"
)

cat(
  "4. TABLE4_RECALIBRATION\n"
)

cat(
  "5. TABLE4B_RECALIBRATION_PARAMETERS\n"
)

cat(
  "6. TABLE5_DCA_MAIN\n"
)

cat(
  "7. TABLE6_DCA_RECAL_SENSITIVITY\n"
)

cat(
  "8. MASTER_SUMMARY_RAW\n"
)

cat(
  "9. Chinese_GWTG_STEP11J_PUBLICATION_TABLES.RData\n"
)

cat("\n")
cat(
  "NOTE: GWTG-HF calibration fields are intentionally ",
  "left unavailable because the GWTG-HF score is not a ",
  "directly calibrated 1-year probability.\n"
)

cat("============================================\n")
