# ============================================================
# STEP 11H
# Paired bootstrap differences in Decision Curve Analysis
#
# Comparisons:
#   1. Base6 original vs MAGGIC
#   2. Base6 recalibrated: intercept only vs MAGGIC
#   3. Base6 recalibrated: intercept + slope vs MAGGIC
#
# Uses the SAME bootstrap resamples generated in STEP 11G.
#
# Input:
#   Chinese_GWTG_STEP11G_DCA.RData
#
# Outputs:
#   - threshold-specific delta net benefit
#   - bootstrap 95% CI
#   - two-sided bootstrap P value
# ============================================================


# ============================================================
# 0. Settings
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

load(
  "Chinese_GWTG_STEP11G_DCA.RData"
)


cat("\n")
cat("============================================\n")
cat("STEP 11H: Paired bootstrap DCA differences\n")
cat("============================================\n\n")


# ============================================================
# 1. Input audit
# ============================================================

if (!exists("dca_boot")) {
  stop(
    "ERROR: dca_boot not found in STEP 11G RData."
  )
}

if (!exists("dca_final")) {
  stop(
    "ERROR: dca_final not found in STEP 11G RData."
  )
}


required_boot_cols <- c(
  "model",
  "threshold",
  "net_benefit",
  "bootstrap_id"
)

missing_boot_cols <- setdiff(
  required_boot_cols,
  names(dca_boot)
)

if (length(missing_boot_cols) > 0) {

  stop(
    paste0(
      "Missing columns in dca_boot: ",
      paste(
        missing_boot_cols,
        collapse = ", "
      )
    )
  )
}


cat("============================================\n")
cat("A. Bootstrap audit\n")
cat("============================================\n\n")

cat(
  "Bootstrap rows =",
  nrow(dca_boot),
  "\n"
)

cat(
  "Valid bootstrap IDs =",
  length(
    unique(
      dca_boot$bootstrap_id
    )
  ),
  "\n"
)

cat(
  "Threshold range =",
  min(dca_boot$threshold),
  "to",
  max(dca_boot$threshold),
  "\n\n"
)

cat("Models:\n")

print(
  unique(
    dca_boot$model
  )
)

cat("\n")


# ============================================================
# 2. Define comparisons
# ============================================================

comparison_map <- list(

  "Base6 original - MAGGIC" = c(
    "Base6 original",
    "MAGGIC"
  ),

  "Base6 recalibrated intercept - MAGGIC" = c(
    "Base6 recalibrated: intercept only",
    "MAGGIC"
  ),

  "Base6 recalibrated full - MAGGIC" = c(
    "Base6 recalibrated: intercept + slope",
    "MAGGIC"
  )
)


thresholds <- sort(
  unique(
    dca_boot$threshold
  )
)


# ============================================================
# 3. Function to calculate paired bootstrap difference
# ============================================================

calc_paired_nb_difference <- function(
    boot_data,
    model_a,
    model_b,
    threshold_value
) {

  a <- boot_data[
    boot_data$model == model_a &
      abs(
        boot_data$threshold -
          threshold_value
      ) < 1e-10,
    c(
      "bootstrap_id",
      "net_benefit"
    )
  ]

  b <- boot_data[
    boot_data$model == model_b &
      abs(
        boot_data$threshold -
          threshold_value
      ) < 1e-10,
    c(
      "bootstrap_id",
      "net_benefit"
    )
  ]


  names(a)[2] <- "NB_A"
  names(b)[2] <- "NB_B"


  paired <- merge(
    a,
    b,
    by = "bootstrap_id"
  )


  paired$delta_NB <-
    paired$NB_A -
    paired$NB_B


  paired <- paired[
    is.finite(
      paired$delta_NB
    ),
    ,
    drop = FALSE
  ]


  if (nrow(paired) == 0) {

    return(
      c(
        bootstrap_mean = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_,
        p_boot = NA_real_,
        probability_A_better = NA_real_,
        valid_B = 0
      )
    )
  }


  ci <- quantile(
    paired$delta_NB,
    probs = c(
      0.025,
      0.975
    ),
    na.rm = TRUE
  )


  p_boot <- 2 * min(
    mean(
      paired$delta_NB <= 0
    ),
    mean(
      paired$delta_NB >= 0
    )
  )

  p_boot <- min(
    p_boot,
    1
  )


  c(

    bootstrap_mean =
      mean(
        paired$delta_NB
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

    probability_A_better =
      mean(
        paired$delta_NB > 0
      ),

    valid_B =
      nrow(
        paired
      )
  )
}


# ============================================================
# 4. Original point-estimate difference function
# ============================================================

get_point_delta <- function(
    final_data,
    model_a,
    model_b,
    threshold_value
) {

  a <- final_data$net_benefit[
    final_data$model == model_a &
      abs(
        final_data$threshold -
          threshold_value
      ) < 1e-10
  ]

  b <- final_data$net_benefit[
    final_data$model == model_b &
      abs(
        final_data$threshold -
          threshold_value
      ) < 1e-10
  ]


  if (
    length(a) != 1 ||
    length(b) != 1
  ) {

    return(
      NA_real_
    )
  }


  as.numeric(
    a - b
  )
}


# ============================================================
# 5. Run all comparisons
# ============================================================

results_list <- list()

counter <- 0


for (comparison_name in names(comparison_map)) {

  models <- comparison_map[[comparison_name]]

  model_a <- models[1]
  model_b <- models[2]


  for (pt in thresholds) {

    counter <- counter + 1


    stats <- calc_paired_nb_difference(
      boot_data = dca_boot,
      model_a = model_a,
      model_b = model_b,
      threshold_value = pt
    )


    point_delta <- get_point_delta(
      final_data = dca_final,
      model_a = model_a,
      model_b = model_b,
      threshold_value = pt
    )


    results_list[[counter]] <- data.frame(

      comparison = comparison_name,

      model_A = model_a,

      model_B = model_b,

      threshold = pt,

      delta_NB = point_delta,

      bootstrap_mean = stats["bootstrap_mean"],

      CI_low = stats["CI_low"],

      CI_high = stats["CI_high"],

      p_boot = stats["p_boot"],

      probability_A_better =
        stats["probability_A_better"],

      valid_B = stats["valid_B"],

      stringsAsFactors = FALSE
    )
  }
}


dca_pairwise <- do.call(
  rbind,
  results_list
)

rownames(dca_pairwise) <- NULL


# ============================================================
# 6. Selected clinically interpretable thresholds
# ============================================================

selected_thresholds <- c(
  0.05,
  0.10,
  0.15,
  0.20,
  0.25,
  0.30
)


selected_flag <- sapply(
  dca_pairwise$threshold,
  function(x) {

    any(
      abs(
        x -
          selected_thresholds
      ) < 1e-8
    )
  }
)


dca_pairwise_selected <-
  dca_pairwise[
    selected_flag,
    ,
    drop = FALSE
  ]


cat("============================================\n")
cat("B. Selected-threshold paired differences\n")
cat("============================================\n\n")

print(
  dca_pairwise_selected[
    ,
    c(
      "comparison",
      "threshold",
      "delta_NB",
      "CI_low",
      "CI_high",
      "p_boot",
      "probability_A_better",
      "valid_B"
    )
  ],
  row.names = FALSE
)

cat("\n")


# ============================================================
# 7. Original Base6 vs MAGGIC only
# ============================================================

main_comparison <- dca_pairwise[
  dca_pairwise$comparison ==
    "Base6 original - MAGGIC",
  ,
  drop = FALSE
]


cat("============================================\n")
cat("C. Main analysis: Base6 original vs MAGGIC\n")
cat("============================================\n\n")

print(
  main_comparison[
    sapply(
      main_comparison$threshold,
      function(x) {

        any(
          abs(
            x -
              selected_thresholds
          ) < 1e-8
        )
      }
    ),
    c(
      "threshold",
      "delta_NB",
      "CI_low",
      "CI_high",
      "p_boot",
      "probability_A_better"
    )
  ],
  row.names = FALSE
)

cat("\n")


# ============================================================
# 8. Significant / supported threshold intervals
# ============================================================

dca_pairwise$CI_excludes_zero <-
  with(
    dca_pairwise,
    CI_low > 0 |
      CI_high < 0
  )


dca_pairwise$direction <- ifelse(
  dca_pairwise$CI_low > 0,
  "Model A higher",
  ifelse(
    dca_pairwise$CI_high < 0,
    "Model B higher",
    "No clear difference"
  )
)


cat("============================================\n")
cat("D. Thresholds with CI excluding zero\n")
cat("============================================\n\n")

supported <- dca_pairwise[
  dca_pairwise$CI_excludes_zero,
  c(
    "comparison",
    "threshold",
    "delta_NB",
    "CI_low",
    "CI_high",
    "p_boot",
    "direction"
  )
]


if (nrow(supported) == 0) {

  cat(
    "No threshold had a bootstrap 95% CI ",
    "excluding zero.\n\n"
  )

} else {

  print(
    supported,
    row.names = FALSE
  )

  cat("\n")
}


# ============================================================
# 9. Save
# ============================================================

write.csv(
  dca_pairwise,
  "Chinese_GWTG_STEP11H_DCA_PAIRED_DIFFERENCES.csv",
  row.names = FALSE
)


write.csv(
  dca_pairwise_selected,
  "Chinese_GWTG_STEP11H_DCA_SELECTED_THRESHOLDS.csv",
  row.names = FALSE
)


save(
  dca_pairwise,
  dca_pairwise_selected,
  main_comparison,
  supported,
  file =
    "Chinese_GWTG_STEP11H_DCA_PAIRED_BOOTSTRAP.RData"
)


cat("============================================\n")
cat("STEP 11H completed\n")
cat("============================================\n\n")

cat("Saved:\n")
cat("1. Chinese_GWTG_STEP11H_DCA_PAIRED_DIFFERENCES.csv\n")
cat("2. Chinese_GWTG_STEP11H_DCA_SELECTED_THRESHOLDS.csv\n")
cat("3. Chinese_GWTG_STEP11H_DCA_PAIRED_BOOTSTRAP.RData\n")
cat("============================================\n")
