# ============================================================
# STEP 11O
# DCA threshold-range sensitivity analysis
#
# Purpose:
#   Evaluate clinical utility over clinically relevant
#   threshold RANGES rather than focusing on isolated
#   threshold-specific significance.
#
# Threshold ranges:
#   1. 5% - 20%
#   2. 5% - 25%
#   3. 10% - 30%
#
# Models:
#   - Base6 original
#   - MAGGIC
#   - Base6 recalibrated intercept
#   - Base6 recalibrated full
#
# Main outcomes:
#   - Mean net benefit over threshold range
#   - Standardized area under NB curve
#   - Mean / standardized integrated delta-NB vs MAGGIC
#   - Paired bootstrap CI and P
#   - Proportion of thresholds with:
#       NB > MAGGIC
#       NB > Treat All
#       NB > Treat None
#
# Input:
#   Chinese_GWTG_STEP11G_DCA.RData
#
# Uses existing STEP 11G bootstrap results.
# NO DCA bootstrap is rerun.
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
cat("STEP 11O: DCA threshold-range sensitivity\n")
cat("============================================\n\n")


# ============================================================
# 1. Load STEP 11G
# ============================================================

input_file <-
  "Chinese_GWTG_STEP11G_DCA.RData"


if (!file.exists(input_file)) {

  stop(
    paste0(
      "Required file not found: ",
      input_file
    )
  )
}


e11g <- new.env()

load(
  input_file,
  envir = e11g
)


required_objects <- c(
  "dca_final",
  "dca_boot"
)


missing_objects <- required_objects[
  !vapply(
    required_objects,
    exists,
    logical(1),
    envir = e11g
  )
]


if (length(missing_objects) > 0) {

  stop(
    paste0(
      "STEP 11G missing: ",
      paste(
        missing_objects,
        collapse = ", "
      )
    )
  )
}


dca_final <- get(
  "dca_final",
  envir = e11g
)


dca_boot <- get(
  "dca_boot",
  envir = e11g
)


# ============================================================
# 2. Column-name detector
# ============================================================

find_column <- function(
    data,
    candidates,
    label
) {

  found <- candidates[
    candidates %in%
      names(data)
  ]


  if (length(found) == 0) {

    stop(
      paste0(
        "Cannot identify ",
        label,
        " column.\nAvailable columns: ",
        paste(
          names(data),
          collapse = ", "
        )
      )
    )
  }


  found[1]
}


model_col_final <- find_column(
  dca_final,
  c(
    "model",
    "Model",
    "model_name"
  ),
  "model"
)


threshold_col_final <- find_column(
  dca_final,
  c(
    "threshold",
    "Threshold",
    "pt"
  ),
  "threshold"
)


nb_col_final <- find_column(
  dca_final,
  c(
    "net_benefit",
    "NB",
    "nb",
    "NetBenefit"
  ),
  "net benefit"
)


model_col_boot <- find_column(
  dca_boot,
  c(
    "model",
    "Model",
    "model_name"
  ),
  "bootstrap model"
)


threshold_col_boot <- find_column(
  dca_boot,
  c(
    "threshold",
    "Threshold",
    "pt"
  ),
  "bootstrap threshold"
)


nb_col_boot <- find_column(
  dca_boot,
  c(
    "net_benefit",
    "NB",
    "nb",
    "NetBenefit"
  ),
  "bootstrap net benefit"
)


bootstrap_col <- find_column(
  dca_boot,
  c(
    "bootstrap_id",
    "boot_id",
    "bootstrap",
    "replicate",
    "b"
  ),
  "bootstrap ID"
)


# ============================================================
# 3. Standardize column names
# ============================================================

dca_point <- data.frame(

  model =
    as.character(
      dca_final[[model_col_final]]
    ),

  threshold =
    as.numeric(
      dca_final[[threshold_col_final]]
    ),

  net_benefit =
    as.numeric(
      dca_final[[nb_col_final]]
    ),

  stringsAsFactors =
    FALSE
)


dca_boot_std <- data.frame(

  bootstrap_id =
    dca_boot[[bootstrap_col]],

  model =
    as.character(
      dca_boot[[model_col_boot]]
    ),

  threshold =
    as.numeric(
      dca_boot[[threshold_col_boot]]
    ),

  net_benefit =
    as.numeric(
      dca_boot[[nb_col_boot]]
    ),

  stringsAsFactors =
    FALSE
)


# ============================================================
# 4. Model-name audit
# ============================================================

cat("============================================\n")
cat("A. DCA object audit\n")
cat("============================================\n\n")

cat("Point-estimate models:\n")

print(
  unique(
    dca_point$model
  )
)

cat("\nBootstrap models:\n")

print(
  unique(
    dca_boot_std$model
  )
)

cat("\n")

cat(
  "Point thresholds:",
  min(dca_point$threshold),
  "to",
  max(dca_point$threshold),
  "\n"
)

cat(
  "Bootstrap thresholds:",
  min(dca_boot_std$threshold),
  "to",
  max(dca_boot_std$threshold),
  "\n"
)

cat(
  "Bootstrap replicates =",
  length(
    unique(
      dca_boot_std$bootstrap_id
    )
  ),
  "\n\n"
)


# ============================================================
# 5. Identify required model names
# ============================================================

available_models <-
  unique(
    dca_point$model
  )


find_model_name <- function(
    exact_candidates,
    regex_candidates,
    label
) {

  # Exact match first
  for (x in exact_candidates) {

    if (x %in% available_models) {
      return(x)
    }
  }


  # Regex fallback
  for (pat in regex_candidates) {

    hit <- available_models[
      grepl(
        pat,
        available_models,
        ignore.case = TRUE
      )
    ]


    if (length(hit) == 1) {
      return(hit)
    }
  }


  stop(
    paste0(
      "Could not uniquely identify model: ",
      label,
      "\nAvailable models: ",
      paste(
        available_models,
        collapse = "; "
      )
    )
  )
}


model_base6 <- find_model_name(

  exact_candidates = c(
    "Base6 original",
    "Base6"
  ),

  regex_candidates = c(
    "^Base6 original$"
  ),

  label = "Base6 original"
)


model_maggic <- find_model_name(

  exact_candidates = c(
    "MAGGIC"
  ),

  regex_candidates = c(
    "^MAGGIC$"
  ),

  label = "MAGGIC"
)


model_recal_intercept <- find_model_name(

  exact_candidates = c(
    "Base6 recalibrated intercept",
    "Base6 recalibrated: intercept only"
  ),

  regex_candidates = c(
    "recalibrated.*intercept"
  ),

  label = "Base6 recalibrated intercept"
)


model_recal_full <- find_model_name(

  exact_candidates = c(
    "Base6 recalibrated full",
    "Base6 recalibrated: intercept + slope"
  ),

  regex_candidates = c(
    "recalibrated.*full",
    "intercept.*slope"
  ),

  label = "Base6 recalibrated full"
)


analysis_models <- c(
  model_base6,
  model_maggic,
  model_recal_intercept,
  model_recal_full
)


cat("============================================\n")
cat("B. Models used in STEP 11O\n")
cat("============================================\n\n")

print(
  data.frame(
    role = c(
      "Base6 original",
      "MAGGIC",
      "Base6 recalibrated intercept",
      "Base6 recalibrated full"
    ),
    object_model_name =
      analysis_models,
    stringsAsFactors = FALSE
  ),
  row.names = FALSE
)

cat("\n")


# ============================================================
# 6. Restrict to required models
# ============================================================

dca_point <- dca_point[
  dca_point$model %in%
    analysis_models,
  ,
  drop = FALSE
]


dca_boot_std <- dca_boot_std[
  dca_boot_std$model %in%
    analysis_models,
  ,
  drop = FALSE
]


# ============================================================
# 7. Define threshold ranges
# ============================================================

range_definitions <- data.frame(

  range_name = c(
    "5%-20%",
    "5%-25%",
    "10%-30%"
  ),

  lower = c(
    0.05,
    0.05,
    0.10
  ),

  upper = c(
    0.20,
    0.25,
    0.30
  ),

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("C. Prespecified sensitivity ranges\n")
cat("============================================\n\n")

print(
  range_definitions,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 8. Trapezoidal integration helper
# ============================================================

trapz_manual <- function(
    x,
    y
) {

  valid <-
    is.finite(x) &
    is.finite(y)


  x <- x[valid]
  y <- y[valid]


  if (length(x) < 2) {
    return(NA_real_)
  }


  ord <- order(x)

  x <- x[ord]
  y <- y[ord]


  sum(
    diff(x) *
      (
        head(y, -1) +
          tail(y, -1)
      ) /
      2
  )
}


# ============================================================
# 9. Interval-summary function
# ============================================================

summarize_interval <- function(
    data,
    lower,
    upper
) {

  d <- data[
    data$threshold >=
      lower - 1e-10 &
      data$threshold <=
      upper + 1e-10,
    ,
    drop = FALSE
  ]


  d <- d[
    order(
      d$threshold
    ),
    ,
    drop = FALSE
  ]


  if (nrow(d) == 0) {

    return(
      c(
        n_thresholds = 0,
        mean_NB = NA_real_,
        area_NB = NA_real_,
        standardized_area_NB = NA_real_
      )
    )
  }


  area <- trapz_manual(
    d$threshold,
    d$net_benefit
  )


  standardized_area <-
    area /
    (upper - lower)


  c(

    n_thresholds =
      nrow(d),

    mean_NB =
      mean(
        d$net_benefit
      ),

    area_NB =
      area,

    standardized_area_NB =
      standardized_area
  )
}


# ============================================================
# 10. Point-estimate interval summary
# ============================================================

point_interval_rows <- list()

counter <- 0


for (
  r in seq_len(
    nrow(
      range_definitions
    )
  )
) {

  lower <-
    range_definitions$lower[r]

  upper <-
    range_definitions$upper[r]

  range_name <-
    range_definitions$range_name[r]


  for (
    model_name in analysis_models
  ) {

    d <- dca_point[
      dca_point$model ==
        model_name,
      ,
      drop = FALSE
    ]


    stats <- summarize_interval(
      d,
      lower,
      upper
    )


    counter <- counter + 1


    point_interval_rows[[counter]] <-
      data.frame(

        range =
          range_name,

        lower =
          lower,

        upper =
          upper,

        model =
          model_name,

        n_thresholds =
          unname(
            stats["n_thresholds"]
          ),

        mean_NB =
          unname(
            stats["mean_NB"]
          ),

        area_NB =
          unname(
            stats["area_NB"]
          ),

        standardized_area_NB =
          unname(
            stats[
              "standardized_area_NB"
            ]
          ),

        stringsAsFactors =
          FALSE
      )
  }
}


point_interval_summary <- do.call(
  rbind,
  point_interval_rows
)


rownames(
  point_interval_summary
) <- NULL


cat("============================================\n")
cat("D. Point-estimate interval net benefit\n")
cat("============================================\n\n")

print(
  point_interval_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 11. Threshold-level reference strategies
#
# Treat-none:
#   NB = 0
#
# Treat-all:
#   NB = observed risk -
#        (1-observed risk)*pt/(1-pt)
#
# Use observed 364-day KM mortality from prior analysis.
# ============================================================

observed_risk <- 0.1475792


reference_grid <- data.frame(

  threshold =
    sort(
      unique(
        dca_point$threshold
      )
    ),

  stringsAsFactors =
    FALSE
)


reference_grid$NB_treat_none <-
  0


reference_grid$NB_treat_all <-
  observed_risk -
  (
    1 -
      observed_risk
  ) *
  reference_grid$threshold /
  (
    1 -
      reference_grid$threshold
  )


# ============================================================
# 12. Threshold dominance proportions
# ============================================================

dominance_rows <- list()

counter <- 0


for (
  r in seq_len(
    nrow(
      range_definitions
    )
  )
) {

  lower <-
    range_definitions$lower[r]

  upper <-
    range_definitions$upper[r]

  range_name <-
    range_definitions$range_name[r]


  maggic_d <- dca_point[
    dca_point$model ==
      model_maggic &
      dca_point$threshold >=
      lower - 1e-10 &
      dca_point$threshold <=
      upper + 1e-10,
    c(
      "threshold",
      "net_benefit"
    ),
    drop = FALSE
  ]


  names(
    maggic_d
  )[2] <- "NB_MAGGIC"


  for (
    model_name in analysis_models
  ) {

    d <- dca_point[
      dca_point$model ==
        model_name &
        dca_point$threshold >=
        lower - 1e-10 &
        dca_point$threshold <=
        upper + 1e-10,
      c(
        "threshold",
        "net_benefit"
      ),
      drop = FALSE
    ]


    names(
      d
    )[2] <- "NB_model"


    merged <- merge(
      d,
      maggic_d,
      by = "threshold",
      all = FALSE
    )


    merged <- merge(
      merged,
      reference_grid,
      by = "threshold",
      all.x = TRUE
    )


    counter <- counter + 1


    dominance_rows[[counter]] <-
      data.frame(

        range =
          range_name,

        model =
          model_name,

        n_thresholds =
          nrow(
            merged
          ),

        proportion_above_MAGGIC =
          if (
            model_name ==
            model_maggic
          ) {

            NA_real_

          } else {

            mean(
              merged$NB_model >
                merged$NB_MAGGIC
            )
          },

        proportion_above_treat_all =
          mean(
            merged$NB_model >
              merged$NB_treat_all
          ),

        proportion_above_treat_none =
          mean(
            merged$NB_model >
              merged$NB_treat_none
          ),

        stringsAsFactors =
          FALSE
      )
  }
}


dominance_summary <- do.call(
  rbind,
  dominance_rows
)


rownames(
  dominance_summary
) <- NULL


cat("============================================\n")
cat("E. Threshold dominance proportions\n")
cat("============================================\n\n")

print(
  dominance_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 13. Bootstrap interval summaries
# ============================================================

bootstrap_ids <- sort(
  unique(
    dca_boot_std$bootstrap_id
  )
)


boot_interval_rows <- list()

counter <- 0


cat("============================================\n")
cat("F. Bootstrap interval aggregation\n")
cat("============================================\n\n")


for (
  r in seq_len(
    nrow(
      range_definitions
    )
  )
) {

  lower <-
    range_definitions$lower[r]

  upper <-
    range_definitions$upper[r]

  range_name <-
    range_definitions$range_name[r]


  for (
    b in bootstrap_ids
  ) {

    db <- dca_boot_std[
      dca_boot_std$bootstrap_id ==
        b,
      ,
      drop = FALSE
    ]


    for (
      model_name in analysis_models
    ) {

      d <- db[
        db$model ==
          model_name,
        ,
        drop = FALSE
      ]


      stats <- summarize_interval(
        d,
        lower,
        upper
      )


      counter <- counter + 1


      boot_interval_rows[[counter]] <-
        data.frame(

          bootstrap_id =
            b,

          range =
            range_name,

          lower =
            lower,

          upper =
            upper,

          model =
            model_name,

          mean_NB =
            unname(
              stats["mean_NB"]
            ),

          area_NB =
            unname(
              stats["area_NB"]
            ),

          standardized_area_NB =
            unname(
              stats[
                "standardized_area_NB"
              ]
            ),

          stringsAsFactors =
            FALSE
        )
    }
  }


  cat(
    "Completed range ",
    range_name,
    "\n",
    sep = ""
  )
}


boot_interval_summary <- do.call(
  rbind,
  boot_interval_rows
)


rownames(
  boot_interval_summary
) <- NULL


cat("\n")


# ============================================================
# 14. Bootstrap distribution helper
# ============================================================

summarize_boot_distribution <- function(
    x,
    point_estimate
) {

  x <- x[
    is.finite(
      x
    )
  ]


  if (length(x) == 0) {

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
# 15. Model-vs-MAGGIC interval differences
# ============================================================

comparison_models <- c(
  model_base6,
  model_recal_intercept,
  model_recal_full
)


interval_difference_rows <- list()

counter <- 0


for (
  range_name in range_definitions$range_name
) {

  # ----------------------------------------------------------
  # MAGGIC bootstrap rows
  # ----------------------------------------------------------

  mag_boot <- boot_interval_summary[
    boot_interval_summary$range ==
      range_name &
      boot_interval_summary$model ==
      model_maggic,
    ,
    drop = FALSE
  ]


  mag_point <- point_interval_summary[
    point_interval_summary$range ==
      range_name &
      point_interval_summary$model ==
      model_maggic,
    ,
    drop = FALSE
  ]


  for (
    model_name in comparison_models
  ) {

    model_boot <- boot_interval_summary[
      boot_interval_summary$range ==
        range_name &
        boot_interval_summary$model ==
        model_name,
      ,
      drop = FALSE
    ]


    model_point <- point_interval_summary[
      point_interval_summary$range ==
        range_name &
        point_interval_summary$model ==
        model_name,
      ,
      drop = FALSE
    ]


    paired <- merge(

      model_boot[
        ,
        c(
          "bootstrap_id",
          "mean_NB",
          "standardized_area_NB"
        ),
        drop = FALSE
      ],

      mag_boot[
        ,
        c(
          "bootstrap_id",
          "mean_NB",
          "standardized_area_NB"
        ),
        drop = FALSE
      ],

      by =
        "bootstrap_id",

      suffixes =
        c(
          "_model",
          "_maggic"
        )
    )


    # --------------------------------------------------------
    # Mean NB difference
    # --------------------------------------------------------

    delta_mean_boot <-
      paired$mean_NB_model -
      paired$mean_NB_maggic


    delta_mean_point <-
      model_point$mean_NB -
      mag_point$mean_NB


    stats_mean <- summarize_boot_distribution(
      delta_mean_boot,
      delta_mean_point
    )


    counter <- counter + 1


    interval_difference_rows[[counter]] <-
      data.frame(

        range =
          range_name,

        comparison =
          paste0(
            model_name,
            " - MAGGIC"
          ),

        summary_metric =
          "Mean net benefit difference",

        point_estimate =
          unname(
            stats_mean[
              "point_estimate"
            ]
          ),

        bootstrap_mean =
          unname(
            stats_mean[
              "bootstrap_mean"
            ]
          ),

        bootstrap_median =
          unname(
            stats_mean[
              "bootstrap_median"
            ]
          ),

        CI_low =
          unname(
            stats_mean[
              "CI_low"
            ]
          ),

        CI_high =
          unname(
            stats_mean[
              "CI_high"
            ]
          ),

        p_boot =
          unname(
            stats_mean[
              "p_boot"
            ]
          ),

        probability_positive =
          unname(
            stats_mean[
              "probability_positive"
            ]
          ),

        valid_B =
          unname(
            stats_mean[
              "valid_B"
            ]
          ),

        stringsAsFactors =
          FALSE
      )


    # --------------------------------------------------------
    # Standardized area difference
    # --------------------------------------------------------

    delta_area_boot <-
      paired$standardized_area_NB_model -
      paired$standardized_area_NB_maggic


    delta_area_point <-
      model_point$standardized_area_NB -
      mag_point$standardized_area_NB


    stats_area <- summarize_boot_distribution(
      delta_area_boot,
      delta_area_point
    )


    counter <- counter + 1


    interval_difference_rows[[counter]] <-
      data.frame(

        range =
          range_name,

        comparison =
          paste0(
            model_name,
            " - MAGGIC"
          ),

        summary_metric =
          "Standardized integrated NB difference",

        point_estimate =
          unname(
            stats_area[
              "point_estimate"
            ]
          ),

        bootstrap_mean =
          unname(
            stats_area[
              "bootstrap_mean"
            ]
          ),

        bootstrap_median =
          unname(
            stats_area[
              "bootstrap_median"
            ]
          ),

        CI_low =
          unname(
            stats_area[
              "CI_low"
            ]
          ),

        CI_high =
          unname(
            stats_area[
              "CI_high"
            ]
          ),

        p_boot =
          unname(
            stats_area[
              "p_boot"
            ]
          ),

        probability_positive =
          unname(
            stats_area[
              "probability_positive"
            ]
          ),

        valid_B =
          unname(
            stats_area[
              "valid_B"
            ]
          ),

        stringsAsFactors =
          FALSE
      )
  }
}


interval_difference_summary <- do.call(
  rbind,
  interval_difference_rows
)


rownames(
  interval_difference_summary
) <- NULL


cat("============================================\n")
cat("G. Interval-level paired differences vs MAGGIC\n")
cat("============================================\n\n")

print(
  interval_difference_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 16. Focused summary:
# standardized integrated NB only
# ============================================================

focused_interval_summary <-
  interval_difference_summary[
    interval_difference_summary$summary_metric ==
      "Standardized integrated NB difference",
    ,
    drop = FALSE
  ]


focused_interval_summary$CI_excludes_zero <-
  with(
    focused_interval_summary,
    CI_low > 0 |
      CI_high < 0
  )


focused_interval_summary$direction <-
  ifelse(
    focused_interval_summary$point_estimate > 0,
    "Model higher",
    ifelse(
      focused_interval_summary$point_estimate < 0,
      "MAGGIC higher",
      "No difference"
    )
  )


cat("============================================\n")
cat("H. Focused interval DCA robustness\n")
cat("============================================\n\n")

print(
  focused_interval_summary,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 17. Threshold-wise ΔNB summary WITHOUT threshold-wise P values
# ============================================================

threshold_delta_rows <- list()

counter <- 0


for (
  model_name in comparison_models
) {

  model_d <- dca_point[
    dca_point$model ==
      model_name,
    c(
      "threshold",
      "net_benefit"
    ),
    drop = FALSE
  ]


  names(
    model_d
  )[2] <- "NB_model"


  mag_d <- dca_point[
    dca_point$model ==
      model_maggic,
    c(
      "threshold",
      "net_benefit"
    ),
    drop = FALSE
  ]


  names(
    mag_d
  )[2] <- "NB_MAGGIC"


  merged <- merge(
    model_d,
    mag_d,
    by = "threshold"
  )


  merged$delta_NB <-
    merged$NB_model -
    merged$NB_MAGGIC


  merged$model <-
    model_name


  counter <- counter + 1

  threshold_delta_rows[[counter]] <-
    merged
}


threshold_delta <- do.call(
  rbind,
  threshold_delta_rows
)


rownames(
  threshold_delta
) <- NULL


# ============================================================
# 18. Sensitivity conclusion flags
# ============================================================

sensitivity_flags <- data.frame(

  range =
    focused_interval_summary$range,

  comparison =
    focused_interval_summary$comparison,

  delta_integrated_NB =
    focused_interval_summary$point_estimate,

  CI_low =
    focused_interval_summary$CI_low,

  CI_high =
    focused_interval_summary$CI_high,

  p_boot =
    focused_interval_summary$p_boot,

  probability_model_better =
    focused_interval_summary$probability_positive,

  statistically_supported =
    focused_interval_summary$CI_excludes_zero,

  stringsAsFactors =
    FALSE
)


cat("============================================\n")
cat("I. Final DCA sensitivity flags\n")
cat("============================================\n\n")

print(
  sensitivity_flags,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 19. Save outputs
# ============================================================

write.csv(
  point_interval_summary,
  "Chinese_GWTG_STEP11O_INTERVAL_NET_BENEFIT.csv",
  row.names = FALSE
)


write.csv(
  dominance_summary,
  "Chinese_GWTG_STEP11O_THRESHOLD_DOMINANCE.csv",
  row.names = FALSE
)


write.csv(
  interval_difference_summary,
  "Chinese_GWTG_STEP11O_INTERVAL_PAIRED_DIFFERENCES.csv",
  row.names = FALSE
)


write.csv(
  focused_interval_summary,
  "Chinese_GWTG_STEP11O_FOCUSED_INTERVAL_DCA.csv",
  row.names = FALSE
)


write.csv(
  threshold_delta,
  "Chinese_GWTG_STEP11O_THRESHOLD_DELTA_NB.csv",
  row.names = FALSE
)


write.csv(
  sensitivity_flags,
  "Chinese_GWTG_STEP11O_SENSITIVITY_FLAGS.csv",
  row.names = FALSE
)


# ============================================================
# 20. Save RData
# ============================================================

save(
  dca_point,
  dca_boot_std,
  range_definitions,
  point_interval_summary,
  reference_grid,
  dominance_summary,
  boot_interval_summary,
  interval_difference_summary,
  focused_interval_summary,
  threshold_delta,
  sensitivity_flags,
  file =
    "Chinese_GWTG_STEP11O_DCA_RANGE_SENSITIVITY.RData"
)


# ============================================================
# 21. Completion
# ============================================================

cat("============================================\n")
cat("STEP 11O completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11O_INTERVAL_NET_BENEFIT.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11O_THRESHOLD_DOMINANCE.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11O_INTERVAL_PAIRED_DIFFERENCES.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11O_FOCUSED_INTERVAL_DCA.csv\n"
)

cat(
  "5. Chinese_GWTG_STEP11O_THRESHOLD_DELTA_NB.csv\n"
)

cat(
  "6. Chinese_GWTG_STEP11O_SENSITIVITY_FLAGS.csv\n"
)

cat(
  "7. Chinese_GWTG_STEP11O_DCA_RANGE_SENSITIVITY.RData\n"
)

cat("\n")

cat(
  "NOTE: STEP 11O evaluates whole threshold ranges. ",
  "Threshold-specific P values are intentionally not used ",
  "for interpretation.\n"
)

cat("============================================\n")
