# ============================================================
# STEP 8J
# Bootstrap optimism correction for:
# 1-year / 3-year / 5-year
#
# FINAL frozen Base6 model
#
# Metrics:
# - time-dependent AUC
# - IPCW Brier score
#
# 200 bootstrap × 30 MI
#
# 模型结构、beta、spline、LP center、baseline hazard
# 全部保持冻结
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


# ============================================================
# 1. Load final frozen model
# ============================================================

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
cat("STEP 8J: Bootstrap fixed-time validation\n")
cat("============================================\n")


# ============================================================
# 2. Settings
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

B <- 200

set.seed(
  20260822
)


# ============================================================
# 3. Data preparation
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


  return(d)
}


# ============================================================
# 4. Final LP
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
# 5. Absolute risk
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
# 6. IPCW Brier
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

    out <- summary(
      km_cens,
      times = s,
      extend = TRUE
    )$surv

    as.numeric(out)
  }


  G_t <-
    G_at(t)


  contribution <-
    numeric(
      length(time)
    )


  for (i in seq_along(time)) {

    # --------------------------------------------------------
    # death by t
    # --------------------------------------------------------

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
          (1 -
             pred_risk[i])^2 /
          G_i
      }


      # --------------------------------------------------------
      # event-free at t
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
    }
  }


  mean(
    contribution
  )
}


# ============================================================
# 7. time-dependent AUC helper
# ============================================================

get_auc <- function(
    time,
    status,
    marker,
    times
) {

  roc <- timeROC::timeROC(

    T =
      time,

    delta =
      status,

    marker =
      marker,

    cause =
      1,

    weighting =
      "marginal",

    times =
      times,

    iid =
      FALSE
  )


  as.numeric(
    roc$AUC
  )
}


# ============================================================
# 8. Storage
# ============================================================

bootstrap_results <- list()

counter <- 1


# ============================================================
# 9. MI loop
# ============================================================

for (m in 1:imp_final$m) {

  cat(
    "\nMI ",
    m,
    "/",
    imp_final$m,
    "\n",
    sep = ""
  )


  d_original <- prepare_data(
    complete(
      imp_final,
      m
    )
  )


  n <-
    nrow(
      d_original
    )


  d_original$lp <-
    calculate_lp(
      d_original
    )


  # ==========================================================
  # Apparent performance
  # ==========================================================

  auc_app <- get_auc(

    time =
      d_original$followup_month,

    status =
      d_original$death,

    marker =
      d_original$lp,

    times =
      target_times
  )


  brier_app <- sapply(

    target_times,

    function(t) {

      risk <-
        calculate_risk(
          d_original$lp,
          t
        )


      ipcw_brier(

        time =
          d_original$followup_month,

        status =
          d_original$death,

        pred_risk =
          risk,

        t =
          t
      )
    }
  )


  # ==========================================================
  # Bootstrap loop
  # ==========================================================

  for (b in 1:B) {

    ix <- sample.int(

      n,

      size =
        n,

      replace =
        TRUE
    )


    d_boot <-
      d_original[
        ix,
        ,
        drop = FALSE
      ]


    # --------------------------------------------------------
    # safeguard
    # --------------------------------------------------------

    if (
      length(
        unique(
          d_boot$death
        )
      ) < 2
    ) {

      next
    }


    # ========================================================
    # Refit SAME model structure in bootstrap sample
    #
    # Important:
    # spline settings frozen
    # ========================================================

    fit_boot <- tryCatch(

      suppressWarnings(

        coxph(

          Surv(
            followup_month,
            death
          ) ~

            age80 +
            bmxbmi +
            lbxrdw +

            nwr_ns1 +
            nwr_ns2 +
            nwr_ns3 +

            egfr_2021 +
            med_cnt,

          data =
            d_boot,

          ties =
            "efron",

          x =
            TRUE,

          y =
            TRUE
        )
      ),

      error =
        function(e)
          NULL
    )


    if (
      is.null(
        fit_boot
      )
    ) {

      next
    }


    # ========================================================
    # Bootstrap model LP
    #
    # train
    # ========================================================

    lp_train <- predict(

      fit_boot,

      newdata =
        d_boot,

      type =
        "lp"
    )


    # ========================================================
    # Bootstrap model LP
    #
    # test on original sample
    # ========================================================

    lp_test <- predict(

      fit_boot,

      newdata =
        d_original,

      type =
        "lp"
    )


    # ========================================================
    # AUC train/test
    # ========================================================

    auc_train <- tryCatch(

      get_auc(

        time =
          d_boot$followup_month,

        status =
          d_boot$death,

        marker =
          lp_train,

        times =
          target_times
      ),

      error =
        function(e)
          rep(
            NA_real_,
            length(
              target_times
            )
          )
    )


    auc_test <- tryCatch(

      get_auc(

        time =
          d_original$followup_month,

        status =
          d_original$death,

        marker =
          lp_test,

        times =
          target_times
      ),

      error =
        function(e)
          rep(
            NA_real_,
            length(
              target_times
            )
          )
    )


    # ========================================================
    # Brier
    #
    # Here absolute risks must be calibrated to each
    # bootstrap Cox model.
    #
    # Therefore estimate bootstrap baseline hazard.
    # ========================================================

    bh_boot <- basehaz(

      fit_boot,

      centered =
        TRUE
    )


    # --------------------------------------------------------
    # Cox predict(type="lp") is centered consistently
    # with the fit object.
    #
    # baseline hazard from centered=TRUE corresponds
    # to the model's centered LP.
    # --------------------------------------------------------

    get_H0_boot <- function(t) {

      ix2 <-
        which(
          bh_boot$time <= t
        )


      if (
        length(ix2) == 0
      ) {

        return(0)
      }


      tail(
        bh_boot$hazard[
          ix2
        ],
        1
      )
    }


    H0_boot <- sapply(

      target_times,

      get_H0_boot
    )


    # ========================================================
    # Bootstrap train/test risk
    # ========================================================

    for (j in seq_along(target_times)) {

      t <-
        target_times[j]


      # ------------------------------------------------------
      # Risk based on bootstrap model
      # ------------------------------------------------------

      risk_train <-
        1 -
        exp(
          -H0_boot[j] *
            exp(
              lp_train
            )
        )


      risk_test <-
        1 -
        exp(
          -H0_boot[j] *
            exp(
              lp_test
            )
        )


      brier_train <-
        ipcw_brier(

          time =
            d_boot$followup_month,

          status =
            d_boot$death,

          pred_risk =
            risk_train,

          t =
            t
        )


      brier_test <-
        ipcw_brier(

          time =
            d_original$followup_month,

          status =
            d_original$death,

          pred_risk =
            risk_test,

          t =
            t
        )


      bootstrap_results[[counter]] <-
        data.frame(

          imputation =
            m,

          bootstrap =
            b,

          time_month =
            t,

          time_label =
            target_labels[j],

          apparent_AUC =
            auc_app[j],

          AUC_train =
            auc_train[j],

          AUC_test =
            auc_test[j],

          optimism_AUC =
            auc_train[j] -
            auc_test[j],

          apparent_Brier =
            brier_app[j],

          Brier_train =
            brier_train,

          Brier_test =
            brier_test,

          optimism_Brier =
            brier_train -
            brier_test
        )


      counter <-
        counter +
        1
    }
  }
}


# ============================================================
# 10. Combine results
# ============================================================

bootstrap_results <-
  dplyr::bind_rows(
    bootstrap_results
  )


cat("\n============================================\n")
cat("A. Bootstrap completion\n")
cat("============================================\n")

print(
  bootstrap_results %>%
    dplyr::count(
      time_label
    )
)


# ============================================================
# 11. MI-level optimism correction
# ============================================================

mi_summary <-
  bootstrap_results %>%

  dplyr::group_by(
    imputation,
    time_month,
    time_label
  ) %>%

  dplyr::summarise(

    apparent_AUC =
      first(
        apparent_AUC
      ),

    optimism_AUC =
      mean(
        optimism_AUC,
        na.rm =
          TRUE
      ),

    corrected_AUC =
      apparent_AUC -
      optimism_AUC,


    apparent_Brier =
      first(
        apparent_Brier
      ),

    optimism_Brier =
      mean(
        optimism_Brier,
        na.rm =
          TRUE
      ),

    corrected_Brier =
      apparent_Brier -
      optimism_Brier,

    B_valid =
      sum(
        is.finite(
          optimism_AUC
        ) &
          is.finite(
            optimism_Brier
          )
      ),

    .groups =
      "drop"
  )


# ============================================================
# 12. Final pooled summary across MI
# ============================================================

final_summary <-
  mi_summary %>%

  dplyr::group_by(
    time_month,
    time_label
  ) %>%

  dplyr::summarise(

    apparent_AUC =
      median(
        apparent_AUC
      ),

    optimism_AUC =
      median(
        optimism_AUC
      ),

    corrected_AUC =
      median(
        corrected_AUC
      ),

    min_corrected_AUC =
      min(
        corrected_AUC
      ),

    max_corrected_AUC =
      max(
        corrected_AUC
      ),


    apparent_Brier =
      median(
        apparent_Brier
      ),

    optimism_Brier =
      median(
        optimism_Brier
      ),

    corrected_Brier =
      median(
        corrected_Brier
      ),

    min_corrected_Brier =
      min(
        corrected_Brier
      ),

    max_corrected_Brier =
      max(
        corrected_Brier
      ),

    median_B_valid =
      median(
        B_valid
      ),

    .groups =
      "drop"
  )


cat("\n============================================\n")
cat("B. Optimism-corrected fixed-time performance\n")
cat("============================================\n")

print(
  final_summary,
  width = Inf
)


# ============================================================
# 13. Numerical audit
# ============================================================

cat("\n============================================\n")
cat("C. Numerical audit\n")
cat("============================================\n")

cat(
  "Non-finite optimism AUC:",
  sum(
    !is.finite(
      bootstrap_results$optimism_AUC
    )
  ),
  "\n"
)

cat(
  "Non-finite optimism Brier:",
  sum(
    !is.finite(
      bootstrap_results$optimism_Brier
    )
  ),
  "\n"
)


cat(
  "|AUC optimism| > 0.15:",
  sum(
    abs(
      bootstrap_results$optimism_AUC
    ) > 0.15,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "|Brier optimism| > 0.15:",
  sum(
    abs(
      bootstrap_results$optimism_Brier
    ) > 0.15,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 14. Save
# ============================================================

write.csv(

  bootstrap_results,

  "Base6_FINAL_bootstrap_fixed_time_raw.csv",

  row.names =
    FALSE
)


write.csv(

  mi_summary,

  "Base6_FINAL_bootstrap_fixed_time_MI_summary.csv",

  row.names =
    FALSE
)


write.csv(

  final_summary,

  "Base6_FINAL_bootstrap_fixed_time_summary.csv",

  row.names =
    FALSE
)


save(

  bootstrap_results,

  mi_summary,

  final_summary,

  file =
    "NHANES_Base6_FINAL_BOOTSTRAP_FIXED_TIME.RData"
)


cat("\n============================================\n")
cat("STEP 8J 完成\n")
cat("============================================\n")
