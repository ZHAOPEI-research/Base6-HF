# ============================================================
# STEP 8I
# FINAL Base6 supplementary calibration metrics
#
# 1. Null-model IPCW Brier score
# 2. Scaled Brier score
# 3. Calibration slope
# 4. Calibration-in-the-large
# ============================================================

rm(list = ls())

library(mice)
library(survival)
library(splines)
library(dplyr)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load(
  "NHANES_2011_2018_MICE_FINAL.RData"
)

load(
  "NHANES_Base6_FINAL_COEFFICIENTS_v2.RData"
)

load(
  "NHANES_Base6_FINAL_ABSOLUTE_RISK_FIXED.RData"
)

load(
  "NHANES_Base6_FINAL_PERFORMANCE.RData"
)

cat("============================================\n")
cat("STEP 8I: Calibration and scaled Brier\n")
cat("============================================\n")


# ============================================================
# 1. Target times
# ============================================================

target_times <- c(
  12,
  36,
  60
)

target_labels <- c(
  "1-year",
  "3-year",
  "5-year"
)


# ============================================================
# 2. Data preparation
# ============================================================

prepare_data <- function(d) {

  d <- d %>%
    mutate(

      age80 =
        pmin(
          ridageyr,
          80
        ),

      nwr =
        nwr_model,

      kappa =
        ifelse(
          female == 1,
          0.7,
          0.9
        ),

      alpha =
        ifelse(
          female == 1,
          -0.241,
          -0.302
        ),

      egfr_2021 =
        142 *
        pmin(scr_mgdl / kappa, 1)^alpha *
        pmax(scr_mgdl / kappa, 1)^(-1.200) *
        (0.9938^ridageyr) *
        ifelse(
          female == 1,
          1.012,
          1
        )
    )


  nb <- splines::ns(

    d$nwr,

    knots =
      nwr_knots_final,

    Boundary.knots =
      nwr_boundary_final,

    intercept =
      FALSE
  )


  d$nwr_ns1 <- nb[, 1]
  d$nwr_ns2 <- nb[, 2]
  d$nwr_ns3 <- nb[, 3]


  d
}


# ============================================================
# 3. LP
# ============================================================

calculate_lp <- function(d) {

  as.numeric(

    beta_shrunk_final["age80"] *
      d$age80 +

      beta_shrunk_final["bmxbmi"] *
      d$bmxbmi +

      beta_shrunk_final["lbxrdw"] *
      d$lbxrdw +

      beta_shrunk_final["nwr_ns1"] *
      d$nwr_ns1 +

      beta_shrunk_final["nwr_ns2"] *
      d$nwr_ns2 +

      beta_shrunk_final["nwr_ns3"] *
      d$nwr_ns3 +

      beta_shrunk_final["egfr_2021"] *
      d$egfr_2021 +

      beta_shrunk_final["med_cnt"] *
      d$med_cnt
  )
}


# ============================================================
# 4. Risk
# ============================================================

calculate_risk <- function(
    lp,
    t
) {

  H0 <-
    baseline_final$H0_mean[
      baseline_final$time_month == t
    ]


  1 -
    exp(
      -H0 *
        exp(
          lp -
            LP_center
        )
    )
}


# ============================================================
# 5. IPCW Brier
# ============================================================

ipcw_brier <- function(
    time,
    status,
    pred_risk,
    t
) {

  censor_status <-
    1 -
    status


  km_cens <-
    survfit(
      Surv(
        time,
        censor_status
      ) ~ 1
    )


  G_at <- function(s) {

    as.numeric(
      summary(
        km_cens,
        times = s,
        extend = TRUE
      )$surv
    )
  }


  G_t <-
    G_at(t)


  contribution <-
    numeric(
      length(time)
    )


  for (i in seq_along(time)) {

    if (
      status[i] == 1 &&
      time[i] <= t
    ) {

      G_i <-
        G_at(
          max(
            time[i] -
              1e-8,
            0
          )
        )


      if (
        is.finite(G_i) &&
        G_i > 0
      ) {

        contribution[i] <-
          (1 - pred_risk[i])^2 /
          G_i
      }


    } else if (
      time[i] > t
    ) {

      if (
        is.finite(G_t) &&
        G_t > 0
      ) {

        contribution[i] <-
          pred_risk[i]^2 /
          G_t
      }
    }
  }


  mean(
    contribution
  )
}


# ============================================================
# 6. Results
# ============================================================

metric_list <- list()

counter <- 1


for (m in 1:imp_final$m) {

  cat(
    "MI ",
    m,
    "/",
    imp_final$m,
    "\n",
    sep = ""
  )


  d <- prepare_data(
    complete(
      imp_final,
      m
    )
  )


  d$lp <-
    calculate_lp(d)


  d$lp_centered <-
    d$lp -
    LP_center


  for (j in seq_along(target_times)) {

    t <-
      target_times[j]


    risk <-
      calculate_risk(
        d$lp,
        t
      )


    # ========================================================
    # Model Brier
    # ========================================================

    brier_model <-
      ipcw_brier(

        d$followup_month,

        d$death,

        risk,

        t
      )


    # ========================================================
    # Null prediction:
    # KM marginal event risk at time t
    # ========================================================

    sf_null <-
      survfit(

        Surv(
          followup_month,
          death
        ) ~ 1,

        data = d
      )


    S_t <-
      summary(

        sf_null,

        times = t,

        extend = TRUE

      )$surv


    null_risk <-
      rep(
        1 -
          as.numeric(S_t),
        nrow(d)
      )


    brier_null <-
      ipcw_brier(

        d$followup_month,

        d$death,

        null_risk,

        t
      )


    scaled_brier <-
      1 -
      brier_model /
      brier_null


    # ========================================================
    # Calibration slope
    #
    # Ideal = 1
    # ========================================================

    slope_fit <-
      coxph(

        Surv(
          followup_month,
          death
        ) ~
          lp_centered,

        data = d,

        ties = "efron"
      )


    calibration_slope <-
      as.numeric(
        coef(
          slope_fit
        )[1]
      )


    # ========================================================
    # Calibration-in-the-large
    #
    # Fixed LP as offset;
    # estimate multiplicative recalibration of baseline hazard.
    #
    # log baseline multiplier ideal = 0
    # ========================================================

    # Observed KM survival
    observed_risk <-
      1 -
      as.numeric(S_t)


    mean_predicted <-
      mean(
        risk
      )


    metric_list[[counter]] <-
      data.frame(

        imputation =
          m,

        time_month =
          t,

        time_label =
          target_labels[j],

        Brier_model =
          brier_model,

        Brier_null =
          brier_null,

        Scaled_Brier =
          scaled_brier,

        calibration_slope =
          calibration_slope,

        mean_predicted =
          mean_predicted,

        observed_KM =
          observed_risk,

        calibration_difference =
          mean_predicted -
          observed_risk
      )


    counter <-
      counter + 1
  }
}


metrics_raw <-
  dplyr::bind_rows(
    metric_list
  )


# ============================================================
# 7. Summary
# ============================================================

metrics_summary <-
  metrics_raw %>%

  dplyr::group_by(
    time_month,
    time_label
  ) %>%

  dplyr::summarise(

    Brier_model =
      median(
        Brier_model
      ),

    Brier_null =
      median(
        Brier_null
      ),

    Scaled_Brier =
      median(
        Scaled_Brier
      ),

    calibration_slope =
      median(
        calibration_slope
      ),

    mean_predicted =
      median(
        mean_predicted
      ),

    observed_KM =
      median(
        observed_KM
      ),

    calibration_difference =
      median(
        calibration_difference
      ),

    .groups =
      "drop"
  )


cat("\n============================================\n")
cat("A. Scaled Brier and overall calibration\n")
cat("============================================\n")

print(
  metrics_summary,
  width = Inf
)


# ============================================================
# 8. Calibration ratio
# ============================================================

metrics_summary <- metrics_summary %>%
  mutate(

    observed_predicted_ratio =
      observed_KM /
      mean_predicted
  )


cat("\n============================================\n")
cat("B. Observed / predicted risk\n")
cat("============================================\n")

print(
  metrics_summary %>%
    dplyr::select(
      time_label,
      mean_predicted,
      observed_KM,
      calibration_difference,
      observed_predicted_ratio
    ),
  width = Inf
)


# ============================================================
# 9. Save
# ============================================================

write.csv(

  metrics_raw,

  "Base6_FINAL_calibration_metrics_raw.csv",

  row.names = FALSE
)


write.csv(

  metrics_summary,

  "Base6_FINAL_calibration_metrics_summary.csv",

  row.names = FALSE
)


save(

  metrics_raw,

  metrics_summary,

  file =
    "NHANES_Base6_FINAL_CALIBRATION_METRICS.RData"
)


cat("\n============================================\n")
cat("STEP 8I 完成\n")
cat("============================================\n")
