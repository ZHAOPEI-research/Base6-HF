# ============================================================
# STEP 8H
# FINAL Base6 fixed-time performance
#
# 1. IPCW time-dependent AUC: 1 / 3 / 5 years
# 2. IPCW Brier score
# 3. 3-year calibration by risk deciles
#
# 使用已经冻结的：
# beta
# spline knots
# LP center
# baseline hazard
# ============================================================

rm(list = ls())

library(mice)
library(survival)
library(splines)
library(dplyr)
library(timeROC)

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

cat("============================================\n")
cat("STEP 8H: FINAL performance assessment\n")
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

        pmin(
          scr_mgdl / kappa,
          1
        )^alpha *

        pmax(
          scr_mgdl / kappa,
          1
        )^(-1.200) *

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
# 3. Linear predictor
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
# 4. Absolute risk from final model
# ============================================================

calculate_risk <- function(
    lp,
    t
) {

  H0 <-
    baseline_final$H0_mean[
      baseline_final$time_month ==
        t
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
# 5. IPCW Brier score
#
# Y(t)=1 if death occurs by t
#
# Event before t:
# weight = 1/G(T-)
#
# Event-free beyond t:
# weight = 1/G(t)
#
# Censored before t contributes 0
# ============================================================

ipcw_brier <- function(
    time,
    status,
    pred_risk,
    t
) {

  # ----------------------------------------------------------
  # KM for censoring distribution:
  # censoring is treated as the "event"
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # G(s) helper
  # ----------------------------------------------------------

  G_at <- function(s) {

    z <- summary(
      km_cens,
      times = s,
      extend = TRUE
    )$surv

    as.numeric(z)
  }


  # G(t)
  G_t <-
    G_at(t)


  contribution <-
    numeric(
      length(time)
    )


  for (i in seq_along(time)) {

    # --------------------------------------------------------
    # Event by t
    # --------------------------------------------------------

    if (
      status[i] == 1 &&
      time[i] <= t
    ) {

      # G(T_i -)
      ti_left <-
        max(
          time[i] -
            1e-8,
          0
        )

      G_i <-
        G_at(
          ti_left
        )


      if (
        is.finite(G_i) &&
        G_i > 0
      ) {

        contribution[i] <-
          (1 - pred_risk[i])^2 /
          G_i
      }


      # --------------------------------------------------------
      # Known event-free at t
      # --------------------------------------------------------

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


      # --------------------------------------------------------
      # Censored before t:
      # contribution remains zero
      # --------------------------------------------------------

    }
  }


  mean(
    contribution
  )
}


# ============================================================
# 6. Evaluation in all 30 imputations
# ============================================================

performance_list <- list()

perf_counter <- 1

calibration_list <- list()


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


  d$lp_final <-
    calculate_lp(d)


  # ==========================================================
  # Time-dependent AUC
  # ==========================================================

  roc_obj <- timeROC::timeROC(

    T =
      d$followup_month,

    delta =
      d$death,

    marker =
      d$lp_final,

    cause =
      1,

    weighting =
      "marginal",

    times =
      target_times,

    iid =
      FALSE
  )


  # ==========================================================
  # Fixed times
  # ==========================================================

  for (j in seq_along(target_times)) {

    t <-
      target_times[j]


    pred_risk <-
      calculate_risk(
        d$lp_final,
        t
      )


    brier <-
      ipcw_brier(

        time =
          d$followup_month,

        status =
          d$death,

        pred_risk =
          pred_risk,

        t =
          t
      )


    performance_list[[perf_counter]] <-
      data.frame(

        imputation =
          m,

        time_month =
          t,

        time_label =
          target_labels[j],

        AUC =
          roc_obj$AUC[j],

        Brier =
          brier,

        mean_predicted_risk =
          mean(
            pred_risk
          ),

        median_predicted_risk =
          median(
            pred_risk
          )
      )


    perf_counter <-
      perf_counter +
      1
  }


  # ==========================================================
  # 3-year calibration
  # ==========================================================

  risk3 <-
    calculate_risk(
      d$lp_final,
      36
    )


  # deciles based on current MI predictions
  cuts <-
    unique(
      quantile(
        risk3,
        probs =
          seq(
            0,
            1,
            0.1
          ),
        na.rm =
          TRUE
      )
    )


  # Safety if duplicate cuts
  if (
    length(cuts) >= 4
  ) {

    d$cal_group <-
      cut(

        risk3,

        breaks =
          cuts,

        include.lowest =
          TRUE,

        labels =
          FALSE
      )


    d$risk3 <-
      risk3


    groups <-
      sort(
        unique(
          d$cal_group
        )
      )


    for (g in groups) {

      dg <-
        d[
          d$cal_group ==
            g,
        ]


      sf <-
        survfit(

          Surv(
            followup_month,
            death
          ) ~ 1,

          data =
            dg
        )


      s36 <-
        summary(

          sf,

          times =
            36,

          extend =
            TRUE

        )$surv


      calibration_list[[
        length(
          calibration_list
        ) +
          1
      ]] <-
        data.frame(

          imputation =
            m,

          group =
            g,

          n =
            nrow(
              dg
            ),

          mean_predicted =
            mean(
              dg$risk3
            ),

          observed_KM =
            1 -
            as.numeric(
              s36
            )
        )
    }
  }
}


# ============================================================
# 7. Performance results
# ============================================================

performance_raw <-
  dplyr::bind_rows(
    performance_list
  )


performance_summary <-
  performance_raw %>%

  dplyr::group_by(
    time_month,
    time_label
  ) %>%

  dplyr::summarise(

    AUC_median =
      median(
        AUC
      ),

    AUC_min =
      min(
        AUC
      ),

    AUC_max =
      max(
        AUC
      ),

    Brier_median =
      median(
        Brier
      ),

    Brier_min =
      min(
        Brier
      ),

    Brier_max =
      max(
        Brier
      ),

    mean_predicted_risk =
      median(
        mean_predicted_risk
      ),

    .groups =
      "drop"
  )


cat("\n============================================\n")
cat("A. Fixed-time discrimination and Brier score\n")
cat("============================================\n")

print(
  performance_summary,
  width = Inf
)


# ============================================================
# 8. Pool 3-year calibration groups
# ============================================================

calibration_raw <-
  dplyr::bind_rows(
    calibration_list
  )


calibration_summary <-
  calibration_raw %>%

  dplyr::group_by(
    group
  ) %>%

  dplyr::summarise(

    n_median =
      median(
        n
      ),

    predicted =
      median(
        mean_predicted
      ),

    observed =
      median(
        observed_KM
      ),

    observed_min =
      min(
        observed_KM
      ),

    observed_max =
      max(
        observed_KM
      ),

    .groups =
      "drop"
  )


cat("\n============================================\n")
cat("B. 3-year calibration by decile\n")
cat("============================================\n")

print(
  calibration_summary,
  n = Inf
)


# ============================================================
# 9. Calibration errors
# ============================================================

calibration_summary <-
  calibration_summary %>%
  dplyr::mutate(

    absolute_error =
      abs(
        predicted -
          observed
      )
  )


calibration_metrics <- data.frame(

  MAE_decile =
    mean(
      calibration_summary$
        absolute_error
    ),

  max_error_decile =
    max(
      calibration_summary$
        absolute_error
    )
)


cat("\n============================================\n")
cat("C. 3-year calibration error\n")
cat("============================================\n")

print(
  calibration_metrics,
  row.names = FALSE
)


# ============================================================
# 10. Plot 3-year calibration
# ============================================================

plot(

  calibration_summary$predicted,

  calibration_summary$observed,

  xlim =
    c(
      0,
      max(
        calibration_summary$predicted,
        calibration_summary$observed
      )
    ),

  ylim =
    c(
      0,
      max(
        calibration_summary$predicted,
        calibration_summary$observed
      )
    ),

  xlab =
    "Predicted 3-year mortality",

  ylab =
    "Observed 3-year mortality (Kaplan-Meier)",

  main =
    "Base6 3-year calibration",

  pch =
    19
)


abline(
  0,
  1,
  lty = 2
)


# ============================================================
# 11. Save
# ============================================================

write.csv(

  performance_raw,

  "Base6_FINAL_fixed_time_performance_raw.csv",

  row.names =
    FALSE
)


write.csv(

  performance_summary,

  "Base6_FINAL_fixed_time_performance_summary.csv",

  row.names =
    FALSE
)


write.csv(

  calibration_raw,

  "Base6_FINAL_3year_calibration_raw.csv",

  row.names =
    FALSE
)


write.csv(

  calibration_summary,

  "Base6_FINAL_3year_calibration_summary.csv",

  row.names =
    FALSE
)


save(

  performance_raw,

  performance_summary,

  calibration_raw,

  calibration_summary,

  calibration_metrics,

  file =
    "NHANES_Base6_FINAL_PERFORMANCE.RData"
)


cat("\n============================================\n")
cat("STEP 8H 完成\n")
cat("============================================\n")
