# ============================================================
# STEP 8G
# FINAL Base6 baseline survival
#
# 使用：
# 1. Rubin pooled coefficients
# 2. bootstrap uniform shrinkage
# 3. frozen NWR spline
#
# 在每个MI中：
# 用固定 shrunk LP 作为 offset
# 重新估计 Cox baseline cumulative hazard
#
# 最后跨30个MI汇总 baseline hazard
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


# ============================================================
# 1. Load final model
# ============================================================

load(
  "NHANES_2011_2018_MICE_FINAL.RData"
)

load(
  "NHANES_Base6_FINAL_COEFFICIENTS_v2.RData"
)


cat("============================================\n")
cat("STEP 8G: FINAL baseline survival\n")
cat("============================================\n")


# ============================================================
# 2. Target prediction times
# ============================================================

target_months <- c(
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


  # ----------------------------------------------------------
  # FINAL frozen NWR spline
  # ----------------------------------------------------------

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
# 4. Fixed shrunk LP
# ============================================================

calculate_lp <- function(d) {

  lp <-

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


  as.numeric(lp)
}


# ============================================================
# 5. Helper: cumulative baseline hazard at time t
# ============================================================

get_H0_at_time <- function(
    bh,
    t
) {

  ix <- which(
    bh$time <= t
  )


  if (
    length(ix) == 0
  ) {

    return(0)
  }


  tail(
    bh$hazard[ix],
    1
  )
}


# ============================================================
# 6. Risk-set support
# ============================================================

d_support <- prepare_data(
  complete(
    imp_final,
    1
  )
)


risk_support <- data.frame(

  time_month =
    target_months,

  time_label =
    target_labels,

  n_at_risk =
    sapply(
      target_months,
      function(t) {

        sum(
          d_support$followup_month >= t
        )
      }
    ),

  deaths_by_time =
    sapply(
      target_months,
      function(t) {

        sum(
          d_support$death == 1 &
            d_support$followup_month <= t
        )
      }
    )
)


cat("\n============================================\n")
cat("A. Prediction-time support\n")
cat("============================================\n")

print(
  risk_support,
  row.names = FALSE
)


# ============================================================
# 7. Baseline hazard in each MI
#
# Important:
# β is FIXED at shrunk pooled β.
#
# Cox:
# Surv ~ offset(LP)
#
# This re-estimates H0 after coefficient shrinkage.
# ============================================================

baseline_MI <- list()

baseline_curves <- list()


for (i in 1:imp_final$m) {

  cat(
    "MI ",
    i,
    "/",
    imp_final$m,
    "\n",
    sep = ""
  )


  d <- prepare_data(
    complete(
      imp_final,
      i
    )
  )


  d$lp_final <-
    calculate_lp(d)


  # ----------------------------------------------------------
  # Fixed coefficient model through offset
  # ----------------------------------------------------------

  fit_offset <- coxph(

    Surv(
      followup_month,
      death
    ) ~
      offset(
        lp_final
      ),

    data =
      d,

    ties =
      "efron"
  )


  bh <- basehaz(

    fit_offset,

    centered =
      FALSE
  )


  bh$imputation <-
    i


  baseline_curves[[i]] <-
    bh


  # ----------------------------------------------------------
  # Target-time baseline hazards
  # ----------------------------------------------------------

  H0 <- sapply(

    target_months,

    function(t) {

      get_H0_at_time(
        bh,
        t
      )
    }
  )


  baseline_MI[[i]] <-
    data.frame(

      imputation =
        i,

      time_month =
        target_months,

      time_label =
        target_labels,

      H0 =
        H0,

      S0 =
        exp(
          -H0
        )
    )
}


baseline_MI <-
  bind_rows(
    baseline_MI
  )


baseline_curves <-
  bind_rows(
    baseline_curves
  )


# ============================================================
# 8. Pool baseline hazard across 30 MI
#
# 使用平均 cumulative baseline hazard
# 同时报告MI间范围
# ============================================================

baseline_final <- baseline_MI %>%

  group_by(
    time_month,
    time_label
  ) %>%

  summarise(

    H0_mean =
      mean(H0),

    H0_median =
      median(H0),

    H0_min =
      min(H0),

    H0_max =
      max(H0),

    S0_from_mean_H0 =
      exp(
        -H0_mean
      ),

    S0_mean =
      mean(S0),

    S0_min =
      min(S0),

    S0_max =
      max(S0),

    .groups =
      "drop"
  )


cat("\n============================================\n")
cat("B. FINAL pooled baseline survival\n")
cat("============================================\n")

print(
  baseline_final,
  width = Inf
)


# ============================================================
# 9. MI variability
# ============================================================

cat("\n============================================\n")
cat("C. Baseline survival by MI\n")
cat("============================================\n")

print(
  as.data.frame(baseline_MI),
  row.names = FALSE
)


# ============================================================
# 10. Example risk calculation
#
# Use pooled mean H0
#
# risk(t) =
# 1 - exp[-H0(t) * exp(LP)]
# ============================================================

predict_absolute_risk <- function(
    age,
    bmi,
    rdw,
    nwr,
    egfr,
    med_count,
    time_month
) {

  # ----------------------------------------------------------
  # Match NHANES development definition
  # ----------------------------------------------------------

  age_use <-
    pmin(
      age,
      80
    )


  nb <- splines::ns(

    nwr,

    knots =
      nwr_knots_final,

    Boundary.knots =
      nwr_boundary_final,

    intercept =
      FALSE
  )


  lp <-

    beta_shrunk_final["age80"] *
    age_use +

    beta_shrunk_final["bmxbmi"] *
    bmi +

    beta_shrunk_final["lbxrdw"] *
    rdw +

    beta_shrunk_final["nwr_ns1"] *
    nb[, 1] +

    beta_shrunk_final["nwr_ns2"] *
    nb[, 2] +

    beta_shrunk_final["nwr_ns3"] *
    nb[, 3] +

    beta_shrunk_final["egfr_2021"] *
    egfr +

    beta_shrunk_final["med_cnt"] *
    med_count


  row <- baseline_final[
    baseline_final$time_month ==
      time_month,
  ]


  if (
    nrow(row) != 1
  ) {

    stop(
      "time_month must be one of 12, 36, 60."
    )
  }


  H0 <-
    row$H0_mean


  survival_probability <-
    exp(
      -H0 *
        exp(lp)
    )


  mortality_risk <-
    1 -
    survival_probability


  data.frame(

    time_month =
      time_month,

    linear_predictor =
      as.numeric(lp),

    survival_probability =
      as.numeric(
        survival_probability
      ),

    mortality_risk =
      as.numeric(
        mortality_risk
      )
  )
}


# ============================================================
# 11. Example patient
#
# 这里只检查函数能否运行
# 不是临床阈值
# ============================================================

example_patient <- predict_absolute_risk(

  age =
    70,

  bmi =
    30,

  rdw =
    14,

  nwr =
    0.60,

  egfr =
    60,

  med_count =
    7,

  time_month =
    36
)


cat("\n============================================\n")
cat("D. Example prediction\n")
cat("============================================\n")

print(
  example_patient
)


# ============================================================
# 12. Save
# ============================================================

write.csv(

  risk_support,

  "Base6_prediction_time_support.csv",

  row.names =
    FALSE
)


write.csv(

  baseline_MI,

  "Base6_baseline_survival_by_MI.csv",

  row.names =
    FALSE
)


write.csv(

  baseline_final,

  "Base6_FINAL_baseline_survival.csv",

  row.names =
    FALSE
)


save(

  baseline_MI,

  baseline_final,

  baseline_curves,

  risk_support,

  predict_absolute_risk,

  beta_shrunk_final,

  nwr_knots_final,

  nwr_boundary_final,

  file =
    "NHANES_Base6_FINAL_ABSOLUTE_RISK.RData"
)


cat("\n============================================\n")
cat("STEP 8G 完成\n")
cat("============================================\n")
