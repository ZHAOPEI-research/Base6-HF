# ============================================================
# STEP 11G
# Decision Curve Analysis (DCA)
#
# Endpoint:
#   364-day all-cause mortality
#
# Models:
#   1. Base6 original
#   2. MAGGIC
#   3. Base6 recalibrated: intercept only
#   4. Base6 recalibrated: intercept + slope
#
# Method:
#   IPCW time-dependent net benefit
#
# Main analysis:
#   Base6 original vs MAGGIC
#
# Sensitivity analysis:
#   + recalibrated Base6 variants
#
# Threshold range:
#   0.02 to 0.30
#
# Input:
#   Chinese_GWTG_STEP11F_RECALIBRATION.RData
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

rm(list = ls())

library(survival)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required. ",
    "Install using install.packages('ggplot2')."
  )
}

library(ggplot2)


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

set.seed(20260905)

t0 <- 364

thresholds <- seq(
  0.02,
  0.30,
  by = 0.01
)

B <- 1000


# ============================================================
# 1. Load STEP 11F
# ============================================================

load(
  "Chinese_GWTG_STEP11F_RECALIBRATION.RData"
)

if (!exists("cal_data")) {
  stop(
    "ERROR: cal_data not found in ",
    "Chinese_GWTG_STEP11F_RECALIBRATION.RData"
  )
}


cat("\n")
cat("============================================\n")
cat("STEP 11G: Decision Curve Analysis\n")
cat("============================================\n\n")


# ============================================================
# 2. Input audit
# ============================================================

required_cols <- c(
  "ID",
  "time",
  "status",
  "risk_base6_1y",
  "risk_maggic_1y",
  "risk_base6_recal_intercept",
  "risk_base6_recal_full"
)

missing_cols <- setdiff(
  required_cols,
  names(cal_data)
)

if (length(missing_cols) > 0) {

  stop(
    paste0(
      "Missing required columns: ",
      paste(
        missing_cols,
        collapse = ", "
      )
    )
  )
}


if (
  any(
    !complete.cases(
      cal_data[, required_cols]
    )
  )
) {
  stop("ERROR: missing values in required variables.")
}


if (!all(cal_data$status %in% c(0, 1))) {
  stop("ERROR: status must be coded 0/1.")
}


cat("============================================\n")
cat("A. Cohort audit\n")
cat("============================================\n\n")

cat("N =", nrow(cal_data), "\n")

cat(
  "Deaths =",
  sum(cal_data$status == 1),
  "\n"
)

cat(
  "Follow-up range =",
  min(cal_data$time),
  "to",
  max(cal_data$time),
  "days\n\n"
)

cat(
  "Threshold range =",
  min(thresholds),
  "to",
  max(thresholds),
  "\n\n"
)


# ============================================================
# 3. Censoring survival helper
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


  get_G <- function(t) {

    if (length(t) == 0) {
      return(numeric(0))
    }


    sapply(
      t,
      function(tt) {

        out <- summary(
          censor_fit,
          times = tt,
          extend = TRUE
        )$surv

        as.numeric(out)
      }
    )
  }


  get_G
}


# ============================================================
# 4. IPCW net benefit function
# ============================================================

calc_dca_ipcw <- function(
    data,
    pred,
    thresholds,
    t0 = 364
) {

  n <- nrow(data)

  time <- data$time
  status <- data$status

  get_G <- build_censoring_function(
    time = time,
    status = status
  )


  G_t0 <- get_G(t0)

  if (
    !is.finite(G_t0) ||
    G_t0 <= 0
  ) {
    stop("ERROR: invalid censoring survival at t0.")
  }


  # ----------------------------------------------------------
  # IPCW event weights
  #
  # Death <= t0:
  #   1 / G(T_i-)
  #
  # Event-free > t0:
  #   1 / G(t0)
  #
  # Censored <= t0:
  #   0
  # ----------------------------------------------------------

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
      time[i] <= t0 &&
      status[i] == 1
    ) {

      t_left <- max(
        time[i] - 1e-7,
        0
      )

      g <- get_G(t_left)

      if (
        is.finite(g) &&
        g > 0
      ) {

        w_event[i] <- 1 / g
      }


    } else if (
      time[i] > t0
    ) {

      w_nonevent[i] <- 1 / G_t0
    }
  }


  # ----------------------------------------------------------
  # Weighted event/non-event prevalence
  # ----------------------------------------------------------

  weighted_events <-
    sum(w_event) / n

  weighted_nonevents <-
    sum(w_nonevent) / n


  out <- lapply(
    thresholds,
    function(pt) {

      treat <- pred >= pt


      TP_weighted <-
        sum(
          w_event[treat]
        ) / n


      FP_weighted <-
        sum(
          w_nonevent[treat]
        ) / n


      odds_threshold <-
        pt / (1 - pt)


      nb_model <-
        TP_weighted -
        FP_weighted *
        odds_threshold


      nb_all <-
        weighted_events -
        weighted_nonevents *
        odds_threshold


      nb_none <- 0


      data.frame(

        threshold =
          pt,

        net_benefit =
          nb_model,

        net_benefit_all =
          nb_all,

        net_benefit_none =
          nb_none,

        weighted_event_rate =
          weighted_events,

        weighted_nonevent_rate =
          weighted_nonevents,

        stringsAsFactors =
          FALSE
      )
    }
  )


  do.call(
    rbind,
    out
  )
}


# ============================================================
# 5. Point-estimate DCA for each model
# ============================================================

risk_map <- list(

  "Base6 original" =
    "risk_base6_1y",

  "MAGGIC" =
    "risk_maggic_1y",

  "Base6 recalibrated: intercept only" =
    "risk_base6_recal_intercept",

  "Base6 recalibrated: intercept + slope" =
    "risk_base6_recal_full"
)


dca_list <- list()


for (model_name in names(risk_map)) {

  risk_col <- risk_map[[model_name]]

  tmp <- calc_dca_ipcw(
    data = cal_data,
    pred = cal_data[[risk_col]],
    thresholds = thresholds,
    t0 = t0
  )

  tmp$model <- model_name

  dca_list[[model_name]] <- tmp
}


dca_point <- do.call(
  rbind,
  dca_list
)

rownames(dca_point) <- NULL


cat("============================================\n")
cat("B. Point-estimate DCA complete\n")
cat("============================================\n\n")

cat(
  "Rows =",
  nrow(dca_point),
  "\n\n"
)


# ============================================================
# 6. Bootstrap DCA
# ============================================================

cat("============================================\n")
cat("C. Bootstrap DCA running\n")
cat("============================================\n\n")

cat(
  "Bootstrap repetitions =",
  B,
  "\n\n"
)


boot_results <- vector(
  "list",
  B
)


n <- nrow(cal_data)


for (b in seq_len(B)) {

  idx <- sample(
    seq_len(n),
    size = n,
    replace = TRUE
  )


  d <- cal_data[
    idx,
    ,
    drop = FALSE
  ]


  if (
    length(unique(d$status)) < 2
  ) {
    next
  }


  boot_model_list <- list()


  for (model_name in names(risk_map)) {

    risk_col <- risk_map[[model_name]]


    tmp <- tryCatch(

      calc_dca_ipcw(
        data = d,
        pred = d[[risk_col]],
        thresholds = thresholds,
        t0 = t0
      ),

      error = function(e) NULL
    )


    if (!is.null(tmp)) {

      tmp$model <- model_name
      tmp$bootstrap_id <- b

      boot_model_list[[model_name]] <- tmp
    }
  }


  if (length(boot_model_list) > 0) {

    boot_results[[b]] <- do.call(
      rbind,
      boot_model_list
    )
  }


  if (b %% 100 == 0) {

    cat(
      "Bootstrap",
      b,
      "/",
      B,
      "\n"
    )
  }
}


boot_results <- boot_results[
  !sapply(
    boot_results,
    is.null
  )
]


dca_boot <- do.call(
  rbind,
  boot_results
)

rownames(dca_boot) <- NULL


cat("\n")

cat(
  "Valid bootstrap runs =",
  length(
    unique(
      dca_boot$bootstrap_id
    )
  ),
  "\n\n"
)


# ============================================================
# 7. Bootstrap CI for each model / threshold
# ============================================================

ci_list <- list()

counter <- 0


for (model_name in names(risk_map)) {

  for (pt in thresholds) {

    x <- dca_boot$net_benefit[
      dca_boot$model == model_name &
        abs(
          dca_boot$threshold - pt
        ) < 1e-12
    ]


    x <- x[
      is.finite(x)
    ]


    if (length(x) == 0) {
      next
    }


    counter <- counter + 1


    ci_list[[counter]] <- data.frame(

      model =
        model_name,

      threshold =
        pt,

      NB_boot_mean =
        mean(x),

      NB_CI_low =
        unname(
          quantile(
            x,
            0.025
          )
        ),

      NB_CI_high =
        unname(
          quantile(
            x,
            0.975
          )
        ),

      valid_B =
        length(x),

      stringsAsFactors =
        FALSE
    )
  }
}


dca_ci <- do.call(
  rbind,
  ci_list
)

rownames(dca_ci) <- NULL


# ============================================================
# 8. Merge point estimates with bootstrap CI
# ============================================================

dca_final <- merge(
  dca_point,
  dca_ci,
  by = c(
    "model",
    "threshold"
  ),
  all.x = TRUE,
  sort = FALSE
)


# restore order

model_order <- names(
  risk_map
)


dca_final$model <- factor(
  dca_final$model,
  levels = model_order
)


dca_final <- dca_final[
  order(
    dca_final$model,
    dca_final$threshold
  ),
]


dca_final$model <- as.character(
  dca_final$model
)

rownames(dca_final) <- NULL


# ============================================================
# 9. Treat-all and treat-none curves
# ============================================================

reference_curves <- unique(
  dca_point[
    ,
    c(
      "threshold",
      "net_benefit_all",
      "net_benefit_none"
    )
  ]
)


reference_curves <- reference_curves[
  order(
    reference_curves$threshold
  ),
]


# ============================================================
# 10. Main-analysis dataset
#
# Original Base6 vs MAGGIC
# ============================================================

dca_main <- dca_final[
  dca_final$model %in% c(
    "Base6 original",
    "MAGGIC"
  ),
]


dca_main$model <- factor(
  dca_main$model,
  levels = c(
    "Base6 original",
    "MAGGIC"
  )
)


# ============================================================
# 11. Improved main DCA plot
# Base6 original vs MAGGIC
# ============================================================

dca_main <- dca_final[
  dca_final$model %in% c(
    "Base6 original",
    "MAGGIC"
  ),
]

dca_main$model <- factor(
  dca_main$model,
  levels = c(
    "Base6 original",
    "MAGGIC"
  ),
  labels = c(
    "Base6",
    "MAGGIC"
  )
)


p_main <- ggplot() +

  # ----------------------------------------------------------
# Treat none
# ----------------------------------------------------------

geom_hline(
  yintercept = 0,
  linewidth = 0.8,
  colour = "grey35"
) +

  # ----------------------------------------------------------
# Treat all
# ----------------------------------------------------------

geom_line(
  data = reference_curves,
  aes(
    x = threshold,
    y = net_benefit_all
  ),
  linewidth = 0.9,
  linetype = "dotted",
  colour = "grey45"
) +

  # ----------------------------------------------------------
# Bootstrap CI ribbons
# ----------------------------------------------------------

geom_ribbon(
  data = dca_main,
  aes(
    x = threshold,
    ymin = NB_CI_low,
    ymax = NB_CI_high,
    fill = model,
    group = model
  ),
  alpha = 0.12,
  colour = NA
) +

  # ----------------------------------------------------------
# Model curves
# ----------------------------------------------------------

geom_line(
  data = dca_main,
  aes(
    x = threshold,
    y = net_benefit,
    colour = model,
    linetype = model
  ),
  linewidth = 1.25
) +

  scale_colour_manual(
    values = c(
      "Base6" = "#0072B2",
      "MAGGIC" = "#D55E00"
    )
  ) +

  scale_fill_manual(
    values = c(
      "Base6" = "#0072B2",
      "MAGGIC" = "#D55E00"
    )
  ) +

  scale_linetype_manual(
    values = c(
      "Base6" = "solid",
      "MAGGIC" = "longdash"
    )
  ) +

  scale_x_continuous(
    name = "Threshold probability",
    limits = c(
      min(thresholds),
      max(thresholds)
    ),
    breaks = seq(
      0.05,
      0.30,
      0.05
    ),
    labels = function(x) {
      paste0(
        round(x * 100),
        "%"
      )
    }
  ) +

  ylab(
    "Net benefit"
  ) +

  labs(
    title =
      "Decision curve analysis for 1-year mortality",

    subtitle =
      "Original Base6 versus MAGGIC at 364 days",

    colour =
      "Model",

    fill =
      "Model",

    linetype =
      "Model",

    caption =
      paste0(
        "Grey horizontal line: treat none; ",
        "grey dotted curve: treat all. ",
        "Shaded bands indicate bootstrap 95% CIs."
      )
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    plot.title =
      element_text(
        face = "bold",
        size = 14
      ),

    legend.position =
      "top",

    legend.title =
      element_text(
        face = "bold"
      ),

    legend.key.width =
      grid::unit(
        1.5,
        "cm"
      ),

    plot.caption =
      element_text(
        hjust = 0,
        size = 9
      )
  )


print(
  p_main
)


# ============================================================
# 12. Improved sensitivity DCA plot
# ============================================================

dca_sensitivity <- dca_final


# ------------------------------------------------------------
# Shorter display labels
# ------------------------------------------------------------

display_labels <- c(

  "Base6 original" =
    "Base6",

  "MAGGIC" =
    "MAGGIC",

  "Base6 recalibrated: intercept only" =
    "Base6 recalibrated (intercept)",

  "Base6 recalibrated: intercept + slope" =
    "Base6 recalibrated (full)"
)


dca_sensitivity$model_display <-
  display_labels[
    dca_sensitivity$model
  ]


dca_sensitivity$model_display <- factor(
  dca_sensitivity$model_display,
  levels = c(
    "Base6",
    "MAGGIC",
    "Base6 recalibrated (intercept)",
    "Base6 recalibrated (full)"
  )
)


p_sensitivity <- ggplot() +

  # ----------------------------------------------------------
# Treat none
# ----------------------------------------------------------

geom_hline(
  yintercept = 0,
  linewidth = 0.8,
  colour = "grey35"
) +

  # ----------------------------------------------------------
# Treat all
# ----------------------------------------------------------

geom_line(
  data = reference_curves,
  aes(
    x = threshold,
    y = net_benefit_all
  ),
  linewidth = 0.9,
  linetype = "dotted",
  colour = "grey45"
) +

  # ----------------------------------------------------------
# Model curves
# ----------------------------------------------------------

geom_line(
  data = dca_sensitivity,
  aes(
    x = threshold,
    y = net_benefit,
    colour = model_display,
    linetype = model_display
  ),
  linewidth = 1.15
) +

  scale_colour_manual(
    values = c(
      "Base6" =
        "#0072B2",

      "MAGGIC" =
        "#D55E00",

      "Base6 recalibrated (intercept)" =
        "#009E73",

      "Base6 recalibrated (full)" =
        "#CC79A7"
    )
  ) +

  scale_linetype_manual(
    values = c(
      "Base6" =
        "solid",

      "MAGGIC" =
        "longdash",

      "Base6 recalibrated (intercept)" =
        "dotdash",

      "Base6 recalibrated (full)" =
        "twodash"
    )
  ) +

  scale_x_continuous(
    name =
      "Threshold probability",

    limits =
      c(
        min(thresholds),
        max(thresholds)
      ),

    breaks =
      seq(
        0.05,
        0.30,
        0.05
      ),

    labels =
      function(x) {

        paste0(
          round(
            x * 100
          ),
          "%"
        )
      }
  ) +

  ylab(
    "Net benefit"
  ) +

  labs(
    title =
      "Decision curve sensitivity analysis",

    subtitle =
      "Original and post-hoc recalibrated Base6 versus MAGGIC",

    colour =
      "Prediction strategy",

    linetype =
      "Prediction strategy",

    caption =
      paste0(
        "Grey horizontal line: treat none; ",
        "grey dotted curve: treat all. ",
        "Recalibrated Base6 curves are post-hoc sensitivity analyses."
      )
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    plot.title =
      element_text(
        face =
          "bold",
        size =
          14
      ),

    legend.position =
      "top",

    legend.title =
      element_text(
        face =
          "bold"
      ),

    legend.text =
      element_text(
        size =
          9
      ),

    legend.key.width =
      grid::unit(
        1.4,
        "cm"
      ),

    plot.caption =
      element_text(
        hjust =
          0,
        size =
          9
      )
  )


print(
  p_sensitivity
)


# ============================================================
# 13. Selected clinically interpretable thresholds
# ============================================================

selected_thresholds <- c(
  0.05,
  0.10,
  0.15,
  0.20,
  0.25,
  0.30
)


selected_dca <- dca_final[
  sapply(
    dca_final$threshold,
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
    "model",
    "threshold",
    "net_benefit",
    "NB_CI_low",
    "NB_CI_high",
    "net_benefit_all",
    "net_benefit_none"
  )
]


selected_dca$model <- factor(
  selected_dca$model,
  levels = model_order
)


selected_dca <- selected_dca[
  order(
    selected_dca$threshold,
    selected_dca$model
  ),
]


selected_dca$model <- as.character(
  selected_dca$model
)


cat("============================================\n")
cat("D. Selected-threshold net benefit\n")
cat("============================================\n\n")

print(
  selected_dca,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 14. Net benefit advantage vs MAGGIC
# ============================================================

base6_main <- dca_final[
  dca_final$model ==
    "Base6 original",
  c(
    "threshold",
    "net_benefit"
  )
]


names(
  base6_main
)[2] <-
  "NB_Base6"


maggic_main <- dca_final[
  dca_final$model ==
    "MAGGIC",
  c(
    "threshold",
    "net_benefit"
  )
]


names(
  maggic_main
)[2] <-
  "NB_MAGGIC"


nb_difference <- merge(
  base6_main,
  maggic_main,
  by = "threshold"
)


nb_difference$delta_NB_Base6_minus_MAGGIC <-
  nb_difference$NB_Base6 -
  nb_difference$NB_MAGGIC


cat("============================================\n")
cat("E. Base6 minus MAGGIC net benefit\n")
cat("============================================\n\n")

print(
  nb_difference[
    nb_difference$threshold %in%
      selected_thresholds,
  ],
  row.names = FALSE
)

cat("\n")


# ============================================================
# 15. Save CSV outputs
# ============================================================

write.csv(
  dca_final,
  "Chinese_GWTG_STEP11G_DCA_FULL.csv",
  row.names = FALSE
)


write.csv(
  selected_dca,
  "Chinese_GWTG_STEP11G_DCA_SELECTED_THRESHOLDS.csv",
  row.names = FALSE
)


write.csv(
  nb_difference,
  "Chinese_GWTG_STEP11G_DCA_BASE6_MINUS_MAGGIC.csv",
  row.names = FALSE
)


# ============================================================
# 16. Save figures
# ============================================================

ggsave(
  filename =
    "Chinese_GWTG_STEP11G_DCA_MAIN.png",

  plot =
    p_main,

  width =
    7.5,

  height =
    6.5,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


ggsave(
  filename =
    "Chinese_GWTG_STEP11G_DCA_MAIN.pdf",

  plot =
    p_main,

  width =
    7.5,

  height =
    6.5,

  units =
    "in",

  device =
    cairo_pdf
)


ggsave(
  filename =
    "Chinese_GWTG_STEP11G_DCA_SENSITIVITY.png",

  plot =
    p_sensitivity,

  width =
    8.0,

  height =
    6.5,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


ggsave(
  filename =
    "Chinese_GWTG_STEP11G_DCA_SENSITIVITY.pdf",

  plot =
    p_sensitivity,

  width =
    8.0,

  height =
    6.5,

  units =
    "in",

  device =
    cairo_pdf
)


# ============================================================
# 17. Save RData
# ============================================================

save(
  dca_point,
  dca_boot,
  dca_ci,
  dca_final,
  reference_curves,
  selected_dca,
  nb_difference,
  p_main,
  p_sensitivity,
  file =
    "Chinese_GWTG_STEP11G_DCA.RData"
)


# ============================================================
# 18. Final audit
# ============================================================

cat("============================================\n")
cat("F. Final DCA audit\n")
cat("============================================\n\n")

cat(
  "Threshold range =",
  min(thresholds),
  "to",
  max(thresholds),
  "\n"
)

cat(
  "Number of thresholds =",
  length(thresholds),
  "\n"
)

cat(
  "Bootstrap repetitions requested =",
  B,
  "\n"
)

cat(
  "Valid bootstrap runs =",
  length(
    unique(
      dca_boot$bootstrap_id
    )
  ),
  "\n\n"
)


# ============================================================
# 19. Completion message
# ============================================================

cat("============================================\n")
cat("STEP 11G completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11G_DCA_FULL.csv\n"
)

cat(
  "2. Chinese_GWTG_STEP11G_DCA_SELECTED_THRESHOLDS.csv\n"
)

cat(
  "3. Chinese_GWTG_STEP11G_DCA_BASE6_MINUS_MAGGIC.csv\n"
)

cat(
  "4. Chinese_GWTG_STEP11G_DCA_MAIN.png\n"
)

cat(
  "5. Chinese_GWTG_STEP11G_DCA_MAIN.pdf\n"
)

cat(
  "6. Chinese_GWTG_STEP11G_DCA_SENSITIVITY.png\n"
)

cat(
  "7. Chinese_GWTG_STEP11G_DCA_SENSITIVITY.pdf\n"
)

cat(
  "8. Chinese_GWTG_STEP11G_DCA.RData\n"
)

cat("\n")

cat(
  "NOTE: recalibrated Base6 curves are post-hoc ",
  "sensitivity analyses and should not be presented as ",
  "independent external validation.\n"
)

cat("============================================\n")
