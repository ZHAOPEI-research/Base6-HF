# ============================================================
# STEP 11E
# Calibration plot:
# Base6 vs MAGGIC
#
# Time horizon:
#   364 days
#
# Plot components:
#   1. Ideal 45-degree calibration line
#   2. Quintile-based mean predicted mortality
#   3. Kaplan-Meier observed mortality
#   4. 95% CI for observed KM mortality
#   5. Base6 vs MAGGIC comparison
#
# Input:
#   Chinese_GWTG_STEP11D_CALIBRATION.RData
#
# Outputs:
#   PNG
#   PDF
#   CSV
#   RData
# ============================================================


# ============================================================
# 0. Packages and settings
# ============================================================

rm(list = ls())

library(survival)

if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop(
    "Package 'ggplot2' is required. ",
    "Please install it first using install.packages('ggplot2')."
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


# ============================================================
# 1. Load STEP 11D results
# ============================================================

load(
  "Chinese_GWTG_STEP11D_CALIBRATION.RData"
)


if (!exists("cal_data")) {
  stop(
    "ERROR: cal_data not found in ",
    "Chinese_GWTG_STEP11D_CALIBRATION.RData"
  )
}

if (!exists("final_summary")) {
  stop(
    "ERROR: final_summary not found in ",
    "Chinese_GWTG_STEP11D_CALIBRATION.RData"
  )
}


cat("\n")
cat("============================================\n")
cat("STEP 11E: Calibration plot\n")
cat("Base6 vs MAGGIC\n")
cat("============================================\n\n")


# ============================================================
# 2. Input data audit
# ============================================================

required_cols <- c(
  "ID",
  "time",
  "status",
  "risk_base6_1y",
  "risk_maggic_1y"
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


cat("============================================\n")
cat("A. Input audit\n")
cat("============================================\n\n")

cat(
  "N =",
  nrow(cal_data),
  "\n"
)

cat(
  "Deaths =",
  sum(
    cal_data$status == 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Follow-up range =",
  min(cal_data$time),
  "to",
  max(cal_data$time),
  "days\n\n"
)


if (
  any(
    !complete.cases(
      cal_data[, required_cols]
    )
  )
) {

  stop(
    "ERROR: missing values detected in calibration input."
  )
}


if (
  any(
    cal_data$risk_base6_1y < 0 |
    cal_data$risk_base6_1y > 1
  )
) {

  stop(
    "ERROR: Base6 predicted risk outside [0,1]."
  )
}


if (
  any(
    cal_data$risk_maggic_1y < 0 |
    cal_data$risk_maggic_1y > 1
  )
) {

  stop(
    "ERROR: MAGGIC predicted risk outside [0,1]."
  )
}


cat("Input audit passed.\n\n")


# ============================================================
# 3. Function:
# Kaplan-Meier mortality + 95% CI at t0
# ============================================================

get_km_mortality_ci <- function(
    time,
    status,
    t0 = 364
) {

  fit <- survival::survfit(
    survival::Surv(
      time,
      status
    ) ~ 1,
    conf.type = "log-log"
  )


  s <- summary(
    fit,
    times = t0,
    extend = TRUE
  )


  surv_est <- as.numeric(
    s$surv
  )

  surv_low <- as.numeric(
    s$lower
  )

  surv_high <- as.numeric(
    s$upper
  )


  # ----------------------------------------------------------
  # Convert survival probabilities to mortality probabilities
  #
  # Mortality = 1 - Survival
  #
  # Because the transformation reverses bounds:
  #
  # mortality lower = 1 - survival upper
  # mortality upper = 1 - survival lower
  # ----------------------------------------------------------

  mortality <- 1 - surv_est

  mortality_low <- 1 - surv_high

  mortality_high <- 1 - surv_low


  c(
    observed_KM =
      mortality,

    CI_low =
      mortality_low,

    CI_high =
      mortality_high
  )
}


# ============================================================
# 4. Function:
# construct quintile calibration groups
# ============================================================

make_calibration_plot_data <- function(
    data,
    risk_col,
    model_name,
    t0 = 364
) {

  risk <- data[[risk_col]]


  # ----------------------------------------------------------
  # Rank-based quintiles
  #
  # This guarantees near-equal group sizes even when prediction
  # values contain ties.
  # ----------------------------------------------------------

  r <- rank(
    risk,
    ties.method = "first"
  )


  breaks <- quantile(
    r,
    probs = seq(
      0,
      1,
      by = 0.2
    ),
    names = FALSE
  )


  cal_group <- cut(
    r,
    breaks = breaks,
    include.lowest = TRUE,
    labels = paste0(
      "Q",
      1:5
    )
  )


  out <- lapply(
    levels(cal_group),
    function(g) {

      idx <- which(
        cal_group == g
      )


      d <- data[
        idx,
        ,
        drop = FALSE
      ]


      km <- get_km_mortality_ci(
        time =
          d$time,

        status =
          d$status,

        t0 =
          t0
      )


      data.frame(

        model =
          model_name,

        group =
          g,

        n =
          nrow(d),

        deaths =
          sum(
            d$status == 1,
            na.rm = TRUE
          ),

        mean_predicted =
          mean(
            d[[risk_col]],
            na.rm = TRUE
          ),

        median_predicted =
          median(
            d[[risk_col]],
            na.rm = TRUE
          ),

        observed_KM =
          km["observed_KM"],

        CI_low =
          km["CI_low"],

        CI_high =
          km["CI_high"],

        absolute_error =
          abs(
            km["observed_KM"] -
              mean(
                d[[risk_col]],
                na.rm = TRUE
              )
          ),

        stringsAsFactors = FALSE
      )
    }
  )


  result <- do.call(
    rbind,
    out
  )

  rownames(result) <- NULL

  result
}


# ============================================================
# 5. Build Base6 calibration data
# ============================================================

plot_base6 <- make_calibration_plot_data(
  data =
    cal_data,

  risk_col =
    "risk_base6_1y",

  model_name =
    "Base6",

  t0 =
    t0
)


# ============================================================
# 6. Build MAGGIC calibration data
# ============================================================

plot_maggic <- make_calibration_plot_data(
  data =
    cal_data,

  risk_col =
    "risk_maggic_1y",

  model_name =
    "MAGGIC",

  t0 =
    t0
)


# ============================================================
# 7. Combine plot data
# ============================================================

calibration_plot_data <- rbind(
  plot_base6,
  plot_maggic
)


calibration_plot_data$model <- factor(
  calibration_plot_data$model,
  levels = c(
    "Base6",
    "MAGGIC"
  )
)


calibration_plot_data$group <- factor(
  calibration_plot_data$group,
  levels = paste0(
    "Q",
    1:5
  )
)


cat("============================================\n")
cat("B. Calibration plot data\n")
cat("============================================\n\n")

print(
  calibration_plot_data,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 8. Calibration MAE check
# ============================================================

mae_check <- aggregate(
  absolute_error ~ model,
  data =
    calibration_plot_data,
  FUN =
    mean
)

names(mae_check)[2] <-
  "calibration_MAE"


cat("============================================\n")
cat("C. Calibration MAE check\n")
cat("============================================\n\n")

print(
  mae_check,
  row.names = FALSE
)

cat("\n")


# ============================================================
# 9. Extract model summary annotations
# ============================================================

base6_row <- final_summary[
  final_summary$model == "Base6",
  ,
  drop = FALSE
]

maggic_row <- final_summary[
  final_summary$model == "MAGGIC",
  ,
  drop = FALSE
]


if (
  nrow(base6_row) != 1 ||
  nrow(maggic_row) != 1
) {

  stop(
    "ERROR: Base6/MAGGIC rows not uniquely found ",
    "in final_summary."
  )
}


annotation_text <- paste0(

  "Base6: slope ",
  sprintf(
    "%.2f",
    base6_row$calibration_slope
  ),
  "; Brier ",
  sprintf(
    "%.3f",
    base6_row$Brier_364
  ),

  "     ",

  "MAGGIC: slope ",
  sprintf(
    "%.2f",
    maggic_row$calibration_slope
  ),
  "; Brier ",
  sprintf(
    "%.3f",
    maggic_row$Brier_364
  )
)


cat("============================================\n")
cat("D. Figure annotation\n")
cat("============================================\n\n")

cat(
  annotation_text,
  "\n\n"
)


# ============================================================
# 10. Determine common plotting range
# ============================================================

plot_max <- max(
  calibration_plot_data$mean_predicted,
  calibration_plot_data$CI_high,
  na.rm = TRUE
)


# Give some visual headroom
plot_max <- min(
  1,
  ceiling(
    plot_max * 10
  ) / 10
)


# Avoid too narrow plotting range
plot_max <- max(
  plot_max,
  0.4
)


cat(
  "Common axis maximum =",
  plot_max,
  "\n\n"
)


# ============================================================
# 11. Main combined calibration plot
# ============================================================

p_combined <- ggplot(
  calibration_plot_data,
  aes(
    x = mean_predicted,
    y = observed_KM,
    group = model,
    shape = model,
    linetype = model
  )
) +

  # ----------------------------------------------------------
# Ideal calibration
# ----------------------------------------------------------

geom_abline(
  intercept = 0,
  slope = 1,
  linewidth = 0.7,
  linetype = "dashed"
) +


  # ----------------------------------------------------------
# 95% CI of observed KM mortality
# ----------------------------------------------------------

geom_errorbar(
  aes(
    ymin = CI_low,
    ymax = CI_high
  ),
  width = 0.012,
  linewidth = 0.6
) +


  # ----------------------------------------------------------
# Connect calibration groups
# ----------------------------------------------------------

geom_line(
  linewidth = 0.8
) +


  # ----------------------------------------------------------
# Observed calibration points
# ----------------------------------------------------------

geom_point(
  size = 3
) +


  # ----------------------------------------------------------
# Axes
# ----------------------------------------------------------

scale_x_continuous(
  name =
    "Mean predicted 1-year mortality risk",

  limits =
    c(
      0,
      plot_max
    ),

  breaks =
    seq(
      0,
      plot_max,
      by = 0.1
    ),

  expand =
    c(
      0,
      0
    )
) +

  scale_y_continuous(
    name =
      "Observed 1-year mortality (Kaplan–Meier)",

    limits =
      c(
        0,
        plot_max
      ),

    breaks =
      seq(
        0,
        plot_max,
        by = 0.1
      ),

    expand =
      c(
        0,
        0
      )
  ) +


  # ----------------------------------------------------------
# Model style
# ----------------------------------------------------------

scale_shape_manual(
  values = c(
    "Base6" = 16,
    "MAGGIC" = 17
  )
) +

  scale_linetype_manual(
    values = c(
      "Base6" = "solid",
      "MAGGIC" = "longdash"
    )
  ) +


  # ----------------------------------------------------------
# Titles
# ----------------------------------------------------------

labs(
  title =
    "Calibration of 1-year mortality prediction",

  subtitle =
    "Base6 versus MAGGIC in the common external-validation cohort",

  caption =
    paste0(
      "Points represent quintiles of predicted risk; ",
      "vertical bars indicate 95% CIs for Kaplan–Meier ",
      "observed mortality at 364 days.\n",
      annotation_text
    ),

  shape =
    "Model",

  linetype =
    "Model"
) +


  # ----------------------------------------------------------
# Publication-style theme
# ----------------------------------------------------------

theme_classic(
  base_size = 12
) +

  theme(

    plot.title = element_text(
      size = 14,
      face = "bold"
    ),

    plot.subtitle = element_text(
      size = 11
    ),

    plot.caption = element_text(
      size = 9,
      hjust = 0
    ),

    axis.title = element_text(
      size = 11
    ),

    axis.text = element_text(
      size = 10
    ),

    legend.position =
      "top",

    legend.title =
      element_text(
        face = "bold"
      ),

    aspect.ratio =
      1
  )


# ============================================================
# 12. Display combined plot
# ============================================================

print(
  p_combined
)


# ============================================================
# 13. Faceted version
#
# Useful as supplementary figure / visual audit
# ============================================================

p_facet <- ggplot(
  calibration_plot_data,
  aes(
    x =
      mean_predicted,

    y =
      observed_KM
  )
) +

  geom_abline(
    intercept = 0,
    slope = 1,
    linewidth = 0.7,
    linetype = "dashed"
  ) +

  geom_errorbar(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    width = 0.012,
    linewidth = 0.6
  ) +

  geom_line(
    aes(
      group = 1
    ),
    linewidth = 0.8
  ) +

  geom_point(
    size = 3
  ) +

  facet_wrap(
    ~ model,
    nrow = 1
  ) +

  scale_x_continuous(
    name =
      "Mean predicted 1-year mortality risk",

    limits =
      c(
        0,
        plot_max
      ),

    breaks =
      seq(
        0,
        plot_max,
        by = 0.1
      ),

    expand =
      c(
        0,
        0
      )
  ) +

  scale_y_continuous(
    name =
      "Observed 1-year mortality (Kaplan–Meier)",

    limits =
      c(
        0,
        plot_max
      ),

    breaks =
      seq(
        0,
        plot_max,
        by = 0.1
      ),

    expand =
      c(
        0,
        0
      )
  ) +

  labs(
    title =
      "Calibration of 1-year mortality prediction",

    subtitle =
      "Quintile-based external calibration at 364 days",

    caption =
      paste0(
        "Vertical bars indicate 95% CIs for ",
        "Kaplan–Meier observed mortality."
      )
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(

    plot.title =
      element_text(
        size = 14,
        face = "bold"
      ),

    strip.text =
      element_text(
        size = 12,
        face = "bold"
      ),

    aspect.ratio =
      1,

    plot.caption =
      element_text(
        size = 9,
        hjust = 0
      )
  )


# ============================================================
# 14. Display faceted plot
# ============================================================

print(
  p_facet
)


# ============================================================
# 15. Save combined figure - PNG
# ============================================================

ggsave(
  filename =
    "Chinese_GWTG_STEP11E_CALIBRATION_BASE6_MAGGIC.png",

  plot =
    p_combined,

  width =
    7.2,

  height =
    7.0,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


# ============================================================
# 16. Save combined figure - PDF
# ============================================================

ggsave(
  filename =
    "Chinese_GWTG_STEP11E_CALIBRATION_BASE6_MAGGIC.pdf",

  plot =
    p_combined,

  width =
    7.2,

  height =
    7.0,

  units =
    "in",

  device =
    cairo_pdf
)


# ============================================================
# 17. Save faceted figure
# ============================================================

ggsave(
  filename =
    "Chinese_GWTG_STEP11E_CALIBRATION_FACET.png",

  plot =
    p_facet,

  width =
    9.0,

  height =
    5.0,

  units =
    "in",

  dpi =
    600,

  bg =
    "white"
)


ggsave(
  filename =
    "Chinese_GWTG_STEP11E_CALIBRATION_FACET.pdf",

  plot =
    p_facet,

  width =
    9.0,

  height =
    5.0,

  units =
    "in",

  device =
    cairo_pdf
)


# ============================================================
# 18. Save calibration plotting data
# ============================================================

write.csv(
  calibration_plot_data,
  "Chinese_GWTG_STEP11E_CALIBRATION_PLOT_DATA.csv",
  row.names = FALSE
)


write.csv(
  mae_check,
  "Chinese_GWTG_STEP11E_CALIBRATION_MAE_CHECK.csv",
  row.names = FALSE
)


# ============================================================
# 19. Save RData
# ============================================================

save(
  calibration_plot_data,
  mae_check,
  p_combined,
  p_facet,
  final_summary,
  file =
    "Chinese_GWTG_STEP11E_CALIBRATION_PLOT.RData"
)


# ============================================================
# 20. Final output audit
# ============================================================

cat("============================================\n")
cat("E. Final plot-data audit\n")
cat("============================================\n\n")


cat("Base6:\n")

print(
  calibration_plot_data[
    calibration_plot_data$model == "Base6",
    c(
      "group",
      "n",
      "deaths",
      "mean_predicted",
      "observed_KM",
      "CI_low",
      "CI_high",
      "absolute_error"
    )
  ],
  row.names = FALSE
)

cat("\n")


cat("MAGGIC:\n")

print(
  calibration_plot_data[
    calibration_plot_data$model == "MAGGIC",
    c(
      "group",
      "n",
      "deaths",
      "mean_predicted",
      "observed_KM",
      "CI_low",
      "CI_high",
      "absolute_error"
    )
  ],
  row.names = FALSE
)

cat("\n")


# ============================================================
# 21. Completion message
# ============================================================

cat("============================================\n")
cat("STEP 11E completed\n")
cat("============================================\n\n")

cat("Saved files:\n")

cat(
  "1. Chinese_GWTG_STEP11E_CALIBRATION_BASE6_MAGGIC.png\n"
)

cat(
  "2. Chinese_GWTG_STEP11E_CALIBRATION_BASE6_MAGGIC.pdf\n"
)

cat(
  "3. Chinese_GWTG_STEP11E_CALIBRATION_FACET.png\n"
)

cat(
  "4. Chinese_GWTG_STEP11E_CALIBRATION_FACET.pdf\n"
)

cat(
  "5. Chinese_GWTG_STEP11E_CALIBRATION_PLOT_DATA.csv\n"
)

cat(
  "6. Chinese_GWTG_STEP11E_CALIBRATION_MAE_CHECK.csv\n"
)

cat(
  "7. Chinese_GWTG_STEP11E_CALIBRATION_PLOT.RData\n"
)

cat("\n")
cat("STEP 11E successfully completed.\n")
cat("============================================\n")
