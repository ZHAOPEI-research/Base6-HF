# ============================================================
# STEP 11Q
# Final robustness / sensitivity synthesis
#
# Purpose:
#   Summarize STEP 11K–11P into manuscript-ready tables.
#
# Components:
#   11K Horizon sensitivity
#   11L Modernized MAGGIC sensitivity
#   11M Common vs maximum cohort
#   11N 5000-bootstrap robustness
#   11O DCA threshold-range sensitivity
#   11P Subgroup directional consistency
#
# No new statistical inference is performed.
#
# Outputs:
#   - Overall robustness summary
#   - Key numerical sensitivity table
#   - Manuscript interpretation table
#   - Supplement-ready CSVs
#   - Final RData
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
cat("STEP 11Q: Final robustness synthesis\n")
cat("============================================\n\n")


# ============================================================
# 1. Helper functions
# ============================================================

safe_load <- function(file) {

  if (!file.exists(file)) {
    stop(
      paste0(
        "Required file not found: ",
        file
      )
    )
  }


  e <- new.env()

  load(
    file,
    envir = e
  )

  e
}


get_required <- function(
    env,
    object_name,
    file_label
) {

  if (
    !exists(
      object_name,
      envir = env
    )
  ) {

    stop(
      paste0(
        "Object '",
        object_name,
        "' not found in ",
        file_label
      )
    )
  }


  get(
    object_name,
    envir = env
  )
}


fmt3 <- function(x) {

  ifelse(
    is.na(x),
    NA_character_,
    sprintf(
      "%.3f",
      x
    )
  )
}


fmt4 <- function(x) {

  ifelse(
    is.na(x),
    NA_character_,
    sprintf(
      "%.4f",
      x
    )
  )
}


fmt_ci <- function(
    low,
    high,
    digits = 3
) {

  if (
    is.na(low) ||
    is.na(high)
  ) {
    return(NA_character_)
  }


  fmt <- paste0(
    "%.",
    digits,
    "f"
  )


  paste0(
    sprintf(fmt, low),
    " to ",
    sprintf(fmt, high)
  )
}


fmt_p <- function(p) {

  if (is.na(p)) {
    return(NA_character_)
  }


  if (p < 0.001) {
    return("<0.001")
  }


  sprintf(
    "%.3f",
    p
  )
}


# ============================================================
# 2. Load STEP 11K
# ============================================================

e11k <- safe_load(
  "Chinese_GWTG_STEP11K_HORIZON_SENSITIVITY.RData"
)


# Try common expected names
obj11k <- ls(
  envir = e11k
)


cat("STEP 11K objects:\n")
print(obj11k)
cat("\n")


# ============================================================
# 3. Load STEP 11L
# ============================================================

e11l <- safe_load(
  "Chinese_GWTG_STEP11L_MODERN_MAGGIC_SENSITIVITY.RData"
)


obj11l <- ls(
  envir = e11l
)


cat("STEP 11L objects:\n")
print(obj11l)
cat("\n")


# ============================================================
# 4. Load STEP 11M
# ============================================================

e11m <- safe_load(
  "Chinese_GWTG_STEP11M_COHORT_SENSITIVITY.RData"
)


cohort_audit <- get_required(
  e11m,
  "cohort_audit",
  "STEP 11M"
)


performance_point_11m <- get_required(
  e11m,
  "performance_point",
  "STEP 11M"
)


point_difference_11m <- get_required(
  e11m,
  "point_difference",
  "STEP 11M"
)


boot_summary_11m <- get_required(
  e11m,
  "boot_summary",
  "STEP 11M"
)


# ============================================================
# 5. Load STEP 11N
# ============================================================

e11n <- safe_load(
  "Chinese_GWTG_STEP11N_5000_BOOTSTRAP.RData"
)


point_model_summary_11n <- get_required(
  e11n,
  "point_model_summary",
  "STEP 11N"
)


pairwise_summary_11n <- get_required(
  e11n,
  "pairwise_summary",
  "STEP 11N"
)


focused_summary_11n <- get_required(
  e11n,
  "focused_summary",
  "STEP 11N"
)


# ============================================================
# 6. Load STEP 11O
# ============================================================

e11o <- safe_load(
  "Chinese_GWTG_STEP11O_DCA_RANGE_SENSITIVITY.RData"
)


point_interval_summary_11o <- get_required(
  e11o,
  "point_interval_summary",
  "STEP 11O"
)


focused_interval_summary_11o <- get_required(
  e11o,
  "focused_interval_summary",
  "STEP 11O"
)


dominance_summary_11o <- get_required(
  e11o,
  "dominance_summary",
  "STEP 11O"
)


# ============================================================
# 7. Load STEP 11P
# ============================================================

e11p <- safe_load(
  "Chinese_GWTG_STEP11P_SUBGROUP_ANALYSIS.RData"
)


event_audit_11p <- get_required(
  e11p,
  "event_audit",
  "STEP 11P"
)


subgroup_performance_11p <- get_required(
  e11p,
  "subgroup_performance",
  "STEP 11P"
)


overall_direction_11p <- get_required(
  e11p,
  "overall_direction",
  "STEP 11P"
)


# ============================================================
# 8. Build STEP 11N key robustness table
# ============================================================

key_11n <- pairwise_summary_11n[
  pairwise_summary_11n$comparison %in%
    c(
      "Base6 - MAGGIC",
      "Base6 - GWTG-HF"
    ),
  ,
  drop = FALSE
]


key_11n_table <- data.frame(

  sensitivity =
    "5000-bootstrap discrimination robustness",

  comparison =
    key_11n$comparison,

  metric =
    key_11n$metric,

  estimate =
    key_11n$point_estimate,

  CI =
    mapply(
      fmt_ci,
      key_11n$CI_low,
      key_11n$CI_high
    ),

  p_value =
    vapply(
      key_11n$p_boot,
      fmt_p,
      character(1)
    ),

  probability_positive =
    key_11n$probability_positive,

  interpretation =
    ifelse(
      key_11n$CI_low > 0 &
        key_11n$comparison == "Base6 - MAGGIC" &
        key_11n$metric == "Harrell_C",
      "Base6 showed a statistically supported higher Harrell C-index than MAGGIC.",
      ifelse(
        key_11n$point_estimate > 0,
        "Point estimate favored Base6, but the confidence interval included zero.",
        "No directional advantage for Base6."
      )
    ),

  stringsAsFactors =
    FALSE
)


# ============================================================
# 9. Build STEP 11M cohort-restriction table
# ============================================================

key_11m_table <- data.frame(

  sensitivity =
    "Common vs maximum available cohort",

  comparison =
    point_difference_11m$model,

  metric =
    c(
      "AUC/C",
      "AUC/C",
      "AUC/C"
    ),

  estimate =
    paste0(
      "ΔAUC=",
      fmt4(
        point_difference_11m$delta_AUC_max_minus_common
      ),
      "; ΔC=",
      fmt4(
        point_difference_11m$delta_C_max_minus_common
      )
    ),

  CI =
    NA_character_,

  p_value =
    NA_character_,

  probability_positive =
    NA_real_,

  interpretation =
    ifelse(
      abs(
        point_difference_11m$delta_AUC_max_minus_common
      ) < 0.001 &
        abs(
          point_difference_11m$delta_C_max_minus_common
        ) < 0.001,
      "Performance was virtually unchanged after common-cohort restriction.",
      "Some performance change was observed after cohort restriction."
    ),

  stringsAsFactors =
    FALSE
)


# ============================================================
# 10. Build STEP 11O interval DCA table
# ============================================================

key_11o <- focused_interval_summary_11o


key_11o_table <- data.frame(

  sensitivity =
    "DCA threshold-range sensitivity",

  comparison =
    paste0(
      key_11o$range,
      ": ",
      key_11o$comparison
    ),

  metric =
    "Standardized integrated ΔNB",

  estimate =
    key_11o$point_estimate,

  CI =
    mapply(
      fmt_ci,
      key_11o$CI_low,
      key_11o$CI_high,
      MoreArgs = list(
        digits = 4
      )
    ),

  p_value =
    vapply(
      key_11o$p_boot,
      fmt_p,
      character(1)
    ),

  probability_positive =
    key_11o$probability_positive,

  interpretation =
    ifelse(
      grepl(
        "Base6 original",
        key_11o$comparison,
        fixed = TRUE
      ),
      "Original Base6 showed no consistent integrated net-benefit advantage over MAGGIC.",
      "Recalibrated Base6 showed positive integrated net-benefit differences, but confidence intervals included zero."
    ),

  stringsAsFactors =
    FALSE
)


# ============================================================
# 11. Build STEP 11P subgroup consistency table
# ============================================================

key_11p_table <- data.frame(

  sensitivity =
    "Exploratory subgroup consistency",

  comparison =
    overall_direction_11p$comparison,

  metric =
    "Directional consistency",

  estimate =
    paste0(
      overall_direction_11p$positive_subgroups,
      "/",
      overall_direction_11p$eligible_subgroups,
      " (",
      sprintf(
        "%.0f%%",
        100 *
          overall_direction_11p$proportion_positive
      ),
      ")"
    ),

  CI =
    NA_character_,

  p_value =
    NA_character_,

  probability_positive =
    NA_real_,

  interpretation =
    ifelse(
      overall_direction_11p$comparison ==
        "Base6 C > MAGGIC",
      "Base6 had a higher C-index than MAGGIC in all evaluable subgroups.",
      ifelse(
        overall_direction_11p$proportion_positive >= 0.8,
        "Directional consistency was high across evaluable subgroups.",
        "Directional consistency was moderate across evaluable subgroups."
      )
    ),

  stringsAsFactors =
    FALSE
)


# ============================================================
# 12. Manually encode STEP 11K key horizon results
#
# Based on completed STEP 11K:
#   360 and 364 days virtually identical.
#   365-day timeROC not estimable because no subject had
#   follow-up strictly beyond 365 days.
# ============================================================

key_11k_table <- data.frame(

  sensitivity =
    "Prediction-horizon sensitivity",

  comparison = c(
    "360 vs 364 days",
    "365-day estimability"
  ),

  metric = c(
    "Discrimination/calibration",
    "time-dependent AUC"
  ),

  estimate = c(
    "Virtually unchanged",
    "Not estimable"
  ),

  CI = c(
    NA_character_,
    NA_character_
  ),

  p_value = c(
    NA_character_,
    NA_character_
  ),

  probability_positive = c(
    NA_real_,
    NA_real_
  ),

  interpretation = c(
    "Results at 360 and 364 days were essentially identical, supporting robustness to the exact 1-year horizon definition.",
    "Exact 365-day time-dependent AUC was not estimable because no patient had follow-up strictly beyond 365 days."
  ),

  stringsAsFactors =
    FALSE
)


# ============================================================
# 13. Manually encode STEP 11L key modernized-MAGGIC results
#
# Based on completed STEP 11L.
# ============================================================

key_11l_table <- data.frame(

  sensitivity =
    "Modernized MAGGIC treatment coding",

  comparison = c(
    "Base6 - modernized MAGGIC AUC",
    "Base6 - modernized MAGGIC C",
    "Modernized vs original MAGGIC"
  ),

  metric = c(
    "ΔAUC",
    "ΔC",
    "Overall performance"
  ),

  estimate = c(
    0.0422,
    0.0721,
    NA_real_
  ),

  CI = c(
    "-0.0389 to 0.1300",
    "approximately 0.0000 to 0.1486",
    NA_character_
  ),

  p_value = c(
    "0.315",
    "0.051",
    NA_character_
  ),

  probability_positive = c(
    NA_real_,
    NA_real_,
    NA_real_
  ),

  interpretation = c(
    "Base6 retained a higher AUC point estimate than modernized MAGGIC, without statistical significance.",
    "The Base6 C-index remained higher than modernized MAGGIC, with borderline evidence after modernization.",
    "Modernized MAGGIC modestly improved discrimination and calibration, but did not materially alter the overall comparative pattern."
  ),

  stringsAsFactors =
    FALSE
)


# ============================================================
# 14. Combine overall robustness summary
# ============================================================

robustness_summary <- rbind(
  key_11k_table,
  key_11l_table,
  key_11m_table,
  key_11n_table,
  key_11o_table,
  key_11p_table
)


rownames(
  robustness_summary
) <- NULL


cat("============================================\n")
cat("A. Overall robustness summary\n")
cat("============================================\n\n")

print(
  robustness_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 15. Build compact manuscript-ready synthesis
# ============================================================

manuscript_synthesis <- data.frame(

  analysis = c(
    "Horizon sensitivity",
    "Modernized MAGGIC",
    "Common-cohort restriction",
    "5000-bootstrap discrimination",
    "Threshold-range DCA",
    "Exploratory subgroup analysis"
  ),

  key_finding = c(

    paste0(
      "Performance estimates at 360 and 364 days were essentially identical; ",
      "exact 365-day AUC was not estimable because no participant had follow-up beyond 365 days."
    ),

    paste0(
      "Modernized treatment coding modestly improved MAGGIC performance, ",
      "but Base6 retained higher discrimination point estimates."
    ),

    paste0(
      "Restricting analyses to the common n=346 cohort changed AUC and C-index ",
      "by less than 0.001 for MAGGIC and GWTG-HF."
    ),

    paste0(
      "With 5000 paired bootstrap replicates, Base6 had a higher Harrell C-index ",
      "than MAGGIC (ΔC=0.0813, 95% CI 0.0069 to 0.1563; P=0.0276); ",
      "other pairwise differences were not statistically significant."
    ),

    paste0(
      "Original Base6 showed no consistent integrated net-benefit advantage over MAGGIC. ",
      "After recalibration, integrated ΔNB was positive across all examined ranges, ",
      "although confidence intervals included zero."
    ),

    paste0(
      "Among strata with at least 15 deaths, Base6 C-index exceeded MAGGIC in 10/10 ",
      "evaluable subgroups and GWTG-HF in 9/10; subgroup analyses were descriptive."
    )
  ),

  manuscript_role = c(
    "Supplement",
    "Supplement",
    "Supplement",
    "Main robustness result / Supplement",
    "Supplement",
    "Exploratory Supplement"
  ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("B. Manuscript-ready synthesis\n")
cat("============================================\n\n")

print(
  manuscript_synthesis,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 16. Build one-line robustness conclusions
# ============================================================

robustness_conclusions <- data.frame(

  domain = c(
    "Discrimination",
    "Calibration",
    "Clinical utility",
    "Cohort restriction",
    "Subgroups",
    "Overall"
  ),

  conclusion = c(

    paste0(
      "Base6 consistently showed the highest discrimination point estimates; ",
      "its Harrell C-index was significantly higher than MAGGIC in the 5000-bootstrap analysis."
    ),

    paste0(
      "The main Base6 limitation was underprediction of absolute risk; ",
      "simple intercept-only recalibration substantially improved calibration."
    ),

    paste0(
      "Recalibration improved net-benefit patterns across clinically relevant threshold ranges, ",
      "but interval-level confidence intervals versus MAGGIC included zero."
    ),

    paste0(
      "Using the common three-model cohort had negligible impact on performance estimates."
    ),

    paste0(
      "Exploratory subgroup findings were directionally consistent, particularly for C-index versus MAGGIC, ",
      "but limited event counts preclude strong subgroup claims."
    ),

    paste0(
      "Sensitivity analyses supported the robustness of the main comparative conclusions without demonstrating ",
      "universal superiority of Base6 across every metric or clinical scenario."
    )
  ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("C. Final robustness conclusions\n")
cat("============================================\n\n")

print(
  robustness_conclusions,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 17. Supplement-ready compact numerical table
# ============================================================

supplement_numeric <- data.frame(

  analysis = c(

    "11N Base6 vs MAGGIC",
    "11N Base6 vs GWTG-HF",

    "11O Recalibrated intercept vs MAGGIC, 5%-20%",
    "11O Recalibrated intercept vs MAGGIC, 5%-25%",
    "11O Recalibrated intercept vs MAGGIC, 10%-30%",

    "11M MAGGIC maximum vs common",
    "11M GWTG-HF maximum vs common",

    "11P Base6 C > MAGGIC",
    "11P Base6 C > GWTG-HF"
  ),

  endpoint = c(

    "Harrell C difference",
    "Harrell C difference",

    "Standardized integrated ΔNB",
    "Standardized integrated ΔNB",
    "Standardized integrated ΔNB",

    "ΔAUC / ΔC",
    "ΔAUC / ΔC",

    "Directional subgroup consistency",
    "Directional subgroup consistency"
  ),

  result = c(

    "0.0813 (95% CI 0.0069 to 0.1563), P=0.0276",

    "0.0564 (95% CI -0.0222 to 0.1384), P=0.1436",

    "0.00836 (95% CI -0.00613 to 0.02303)",

    "0.01367 (95% CI -0.00209 to 0.02865)",

    "0.01673 (95% CI -0.00435 to 0.03576)",

    paste0(
      fmt4(
        point_difference_11m$delta_AUC_max_minus_common[
          point_difference_11m$model == "MAGGIC"
        ]
      ),
      " / ",
      fmt4(
        point_difference_11m$delta_C_max_minus_common[
          point_difference_11m$model == "MAGGIC"
        ]
      )
    ),

    paste0(
      fmt4(
        point_difference_11m$delta_AUC_max_minus_common[
          point_difference_11m$model == "GWTG-HF"
        ]
      ),
      " / ",
      fmt4(
        point_difference_11m$delta_C_max_minus_common[
          point_difference_11m$model == "GWTG-HF"
        ]
      )
    ),

    "10/10",

    "9/10"
  ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("D. Supplement-ready numerical summary\n")
cat("============================================\n\n")

print(
  supplement_numeric,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 18. Save CSV outputs
# ============================================================

write.csv(
  robustness_summary,
  "Chinese_GWTG_STEP11Q_ROBUSTNESS_SUMMARY.csv",
  row.names = FALSE
)


write.csv(
  manuscript_synthesis,
  "Chinese_GWTG_STEP11Q_MANUSCRIPT_SYNTHESIS.csv",
  row.names = FALSE
)


write.csv(
  robustness_conclusions,
  "Chinese_GWTG_STEP11Q_ROBUSTNESS_CONCLUSIONS.csv",
  row.names = FALSE
)


write.csv(
  supplement_numeric,
  "Chinese_GWTG_STEP11Q_SUPPLEMENT_NUMERIC.csv",
  row.names = FALSE
)


# ============================================================
# 19. Save complete RData
# ============================================================

save(
  robustness_summary,
  manuscript_synthesis,
  robustness_conclusions,
  supplement_numeric,
  cohort_audit,
  performance_point_11m,
  point_difference_11m,
  boot_summary_11m,
  point_model_summary_11n,
  pairwise_summary_11n,
  focused_summary_11n,
  point_interval_summary_11o,
  focused_interval_summary_11o,
  dominance_summary_11o,
  event_audit_11p,
  subgroup_performance_11p,
  overall_direction_11p,
  file =
    "Chinese_GWTG_STEP11Q_FINAL_ROBUSTNESS_SYNTHESIS.RData"
)


# ============================================================
# 20. Completion
# ============================================================

cat("============================================\n")
cat("STEP 11Q completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11Q_ROBUSTNESS_SUMMARY.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11Q_MANUSCRIPT_SYNTHESIS.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11Q_ROBUSTNESS_CONCLUSIONS.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11Q_SUPPLEMENT_NUMERIC.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11Q_FINAL_ROBUSTNESS_SYNTHESIS.RData\n"
)

cat("============================================\n")
