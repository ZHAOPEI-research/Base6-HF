# ============================================================
# STEP 8G FIXED
# FINAL absolute risk prediction
#
# 修复：
# 手工计算 Breslow baseline cumulative hazard
# 明确使用 centered linear predictor
#
# Final:
# Risk(t) =
# 1 - exp[-H0_centered(t) * exp(LP - LP_center)]
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

cat("============================================\n")
cat("STEP 8G FIXED: Absolute risk\n")
cat("============================================\n")


# ============================================================
# 1. Prediction times
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


  return(d)
}


# ============================================================
# 3. Final LP
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
# 4. First calculate LP in all MI datasets
# ============================================================

completed_list <- vector(
  "list",
  imp_final$m
)

lp_means <- numeric(
  imp_final$m
)


for (i in 1:imp_final$m) {

  d <- prepare_data(
    complete(
      imp_final,
      i
    )
  )

  d$lp_final <-
    calculate_lp(d)


  completed_list[[i]] <-
    d


  lp_means[i] <-
    mean(
      d$lp_final
    )
}


# ============================================================
# 5. Freeze ONE common LP center
#
# 未来外部验证/计算器必须使用同一个值
# ============================================================

LP_center <- mean(
  lp_means
)


cat("\n============================================\n")
cat("A. Linear predictor centering\n")
cat("============================================\n")

cat(
  "MI-specific LP means:\n"
)

print(
  summary(
    lp_means
  )
)

cat(
  "\nFINAL frozen LP center =",
  LP_center,
  "\n"
)


# ============================================================
# 6. Manual Breslow cumulative baseline hazard
# ============================================================

breslow_H0 <- function(
    time,
    status,
    lp,
    center
) {

  lp_c <-
    lp -
    center


  event_times <-
    sort(
      unique(
        time[
          status == 1
        ]
      )
    )


  H_increment <-
    numeric(
      length(
        event_times
      )
    )


  for (j in seq_along(event_times)) {

    tj <-
      event_times[j]


    # Number of events exactly at tj
    dj <-
      sum(
        status == 1 &
          time == tj
      )


    # Risk set
    risk_set <-
      time >= tj


    denom <-
      sum(
        exp(
          lp_c[
            risk_set
          ]
        )
      )


    H_increment[j] <-
      dj /
      denom
  }


  data.frame(

    time =
      event_times,

    hazard_increment =
      H_increment,

    H0 =
      cumsum(
        H_increment
      )
  )
}


# ============================================================
# 7. Helper: H0 at target time
# ============================================================

get_H0 <- function(
    bh,
    t
) {

  ix <-
    which(
      bh$time <= t
    )


  if (
    length(ix) == 0
  ) {

    return(0)
  }


  bh$H0[
    max(ix)
  ]
}


# ============================================================
# 8. Calculate H0 for each MI
# ============================================================

baseline_MI <- list()

baseline_curve_MI <- list()


for (i in 1:imp_final$m) {

  cat(
    "MI ",
    i,
    "/",
    imp_final$m,
    "\n",
    sep = ""
  )


  d <-
    completed_list[[i]]


  bh <-
    breslow_H0(

      time =
        d$followup_month,

      status =
        d$death,

      lp =
        d$lp_final,

      center =
        LP_center
    )


  bh$imputation <-
    i


  baseline_curve_MI[[i]] <-
    bh


  H0_values <-
    sapply(

      target_months,

      function(t) {

        get_H0(
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

      H0_centered =
        H0_values,

      S0_centered =
        exp(
          -H0_values
        )
    )
}


baseline_MI <-
  dplyr::bind_rows(
    baseline_MI
  )


# ============================================================
# 9. Pool baseline hazard
# ============================================================

baseline_final <-
  baseline_MI %>%

  dplyr::group_by(
    time_month,
    time_label
  ) %>%

  dplyr::summarise(

    H0_mean =
      mean(
        H0_centered
      ),

    H0_median =
      median(
        H0_centered
      ),

    H0_min =
      min(
        H0_centered
      ),

    H0_max =
      max(
        H0_centered
      ),

    S0_mean_H =
      exp(
        -H0_mean
      ),

    S0_mean =
      mean(
        S0_centered
      ),

    S0_min =
      min(
        S0_centered
      ),

    S0_max =
      max(
        S0_centered
      ),

    .groups =
      "drop"
  )


cat("\n============================================\n")
cat("B. Corrected centered baseline survival\n")
cat("============================================\n")

print(
  baseline_final,
  width = Inf
)


# ============================================================
# 10. Final absolute-risk function
# ============================================================

predict_absolute_risk <- function(
    age,
    bmi,
    rdw,
    nwr,
    egfr,
    med_count
) {

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


  lp_c <-
    as.numeric(
      lp -
        LP_center
    )


  risk <- sapply(

    1:nrow(
      baseline_final
    ),

    function(j) {

      H0 <-
        baseline_final$H0_mean[j]


      1 -
        exp(
          -H0 *
            exp(
              lp_c
            )
        )
    }
  )


  data.frame(

    time =
      baseline_final$time_label,

    LP =
      as.numeric(lp),

    LP_centered =
      lp_c,

    mortality_risk =
      as.numeric(risk)
  )
}


# ============================================================
# 11. Example patient
# ============================================================

example_patient <-
  predict_absolute_risk(

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
      7
  )


cat("\n============================================\n")
cat("C. Example patient - CORRECTED\n")
cat("============================================\n")

print(
  example_patient,
  row.names = FALSE
)


# ============================================================
# 12. Sanity check:
# risk distribution in development cohort
# ============================================================

d1 <-
  completed_list[[1]]


d1$lp_centered <-
  d1$lp_final -
  LP_center


risk_distribution <- list()


for (j in 1:nrow(baseline_final)) {

  H0 <-
    baseline_final$H0_mean[j]


  risk <-
    1 -
    exp(
      -H0 *
        exp(
          d1$lp_centered
        )
    )


  risk_distribution[[j]] <-
    data.frame(

      time =
        baseline_final$time_label[j],

      min =
        min(risk),

      p05 =
        quantile(
          risk,
          0.05
        ),

      p25 =
        quantile(
          risk,
          0.25
        ),

      median =
        median(risk),

      p75 =
        quantile(
          risk,
          0.75
        ),

      p95 =
        quantile(
          risk,
          0.95
        ),

      max =
        max(risk)
    )
}


risk_distribution <-
  bind_rows(
    risk_distribution
  )


cat("\n============================================\n")
cat("D. Development-cohort predicted-risk distribution\n")
cat("============================================\n")

print(
  risk_distribution,
  row.names = FALSE
)


# ============================================================
# 13. Observed crude event proportions
# only sanity check, not calibration analysis
# ============================================================

cat("\n============================================\n")
cat("E. Observed event counts\n")
cat("============================================\n")


d_support <-
  completed_list[[1]]


observed_support <- data.frame(

  time =
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


print(
  observed_support,
  row.names = FALSE
)


# ============================================================
# 14. Save
# ============================================================

write.csv(

  baseline_MI,

  "Base6_FINAL_centered_baseline_MI.csv",

  row.names = FALSE
)


write.csv(

  baseline_final,

  "Base6_FINAL_centered_baseline.csv",

  row.names = FALSE
)


write.csv(

  risk_distribution,

  "Base6_FINAL_predicted_risk_distribution.csv",

  row.names = FALSE
)


save(

  LP_center,

  baseline_MI,

  baseline_final,

  baseline_curve_MI,

  predict_absolute_risk,

  risk_distribution,

  beta_shrunk_final,

  nwr_knots_final,

  nwr_boundary_final,

  file =
    "NHANES_Base6_FINAL_ABSOLUTE_RISK_FIXED.RData"
)


cat("\n============================================\n")
cat("STEP 8G FIXED 完成\n")
cat("============================================\n")
