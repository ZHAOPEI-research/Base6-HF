# ============================================================
# STEP 8E
# FINAL Base6 bootstrap validation
#
# Final model:
# Age + BMI + RDW + NWR spline + eGFR + Medication count
#
# NWR:
# 3-df natural spline
# internal knots = 1/3 and 2/3 quantiles
# boundary knots = 5th and 95th percentiles
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

cat("============================================\n")
cat("STEP 8E: FINAL Base6 bootstrap\n")
cat("============================================\n")


# ============================================================
# 1. Freeze FINAL spline settings
# ============================================================

d0 <- complete(
  imp_final,
  1
)

nwr0 <- d0$nwr_model


nwr_knots_final <- as.numeric(
  quantile(
    nwr0,
    probs = c(
      1/3,
      2/3
    ),
    na.rm = TRUE
  )
)


nwr_boundary_final <- as.numeric(
  quantile(
    nwr0,
    probs = c(
      0.05,
      0.95
    ),
    na.rm = TRUE
  )
)


cat("\n============================================\n")
cat("A. FINAL NWR spline settings\n")
cat("============================================\n")

cat(
  "Internal knots:",
  paste(
    round(
      nwr_knots_final,
      7
    ),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Boundary knots:",
  paste(
    round(
      nwr_boundary_final,
      7
    ),
    collapse = ", "
  ),
  "\n"
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
# 3. FINAL Base6 formula
# ============================================================

base6_formula <-

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
  med_cnt


# ============================================================
# 4. C-index
# ============================================================

get_cindex <- function(
    time,
    status,
    lp
) {

  z <- survival::concordance(

    Surv(
      time,
      status
    ) ~ lp,

    reverse = TRUE
  )

  as.numeric(
    z$concordance
  )
}


# ============================================================
# 5. Calibration slope
# ============================================================

get_slope <- function(
    data,
    lp
) {

  tmp <- data.frame(

    time =
      data$followup_month,

    status =
      data$death,

    lp =
      as.numeric(lp)
  )


  fit <- coxph(

    Surv(
      time,
      status
    ) ~ lp,

    data =
      tmp,

    ties =
      "efron"
  )


  as.numeric(
    coef(fit)[1]
  )
}


# ============================================================
# 6. Bootstrap settings
# ============================================================

B <- 200

set.seed(
  20260822
)


bootstrap_results <- list()

counter <- 1


# ============================================================
# 7. 30 MI
# ============================================================

for (m in 1:imp_final$m) {

  cat(
    "MI ",
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


  n <- nrow(
    d_original
  )


  # ----------------------------------------------------------
  # Apparent model
  # ----------------------------------------------------------

  fit_original <- coxph(

    base6_formula,

    data =
      d_original,

    ties =
      "efron",

    x =
      TRUE,

    y =
      TRUE
  )


  lp_original <- predict(

    fit_original,

    newdata =
      d_original,

    type =
      "lp"
  )


  apparent_C <- get_cindex(

    d_original$followup_month,

    d_original$death,

    lp_original
  )


  apparent_slope <- get_slope(

    d_original,

    lp_original
  )


  # ==========================================================
  # Bootstrap
  # ==========================================================

  for (b in 1:B) {

    ix <- sample.int(

      n,

      size =
        n,

      replace =
        TRUE
    )


    d_boot <- d_original[
      ix,
      ,
      drop = FALSE
    ]


    if (
      length(
        unique(
          d_boot$death
        )
      ) < 2
    ) {
      next
    }


    fit_boot <- tryCatch(

      suppressWarnings(

        coxph(

          base6_formula,

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


    # --------------------------------------------------------
    # Train performance
    # --------------------------------------------------------

    lp_train <- predict(

      fit_boot,

      newdata =
        d_boot,

      type =
        "lp"
    )


    C_train <- get_cindex(

      d_boot$followup_month,

      d_boot$death,

      lp_train
    )


    slope_train <- get_slope(

      d_boot,

      lp_train
    )


    # --------------------------------------------------------
    # Test on original data
    # --------------------------------------------------------

    lp_test <- predict(

      fit_boot,

      newdata =
        d_original,

      type =
        "lp"
    )


    C_test <- get_cindex(

      d_original$followup_month,

      d_original$death,

      lp_test
    )


    slope_test <- get_slope(

      d_original,

      lp_test
    )


    bootstrap_results[[counter]] <-
      data.frame(

        imputation =
          m,

        bootstrap =
          b,

        apparent_C =
          apparent_C,

        optimism_C =
          C_train -
          C_test,

        apparent_slope =
          apparent_slope,

        optimism_slope =
          slope_train -
          slope_test
      )


    counter <-
      counter +
      1
  }
}


bootstrap_results <- bind_rows(
  bootstrap_results
)


# ============================================================
# 8. Completion
# ============================================================

cat("\n============================================\n")
cat("B. Bootstrap completion\n")
cat("============================================\n")

cat(
  "Valid bootstrap runs:",
  nrow(
    bootstrap_results
  ),
  "/",
  imp_final$m *
    B,
  "\n"
)


# ============================================================
# 9. MI-level summary
# ============================================================

mi_summary <- bootstrap_results %>%

  group_by(
    imputation
  ) %>%

  summarise(

    apparent_C =
      first(
        apparent_C
      ),

    optimism_C =
      mean(
        optimism_C,
        na.rm = TRUE
      ),

    corrected_C =
      apparent_C -
      optimism_C,

    apparent_slope =
      first(
        apparent_slope
      ),

    optimism_slope =
      mean(
        optimism_slope,
        na.rm = TRUE
      ),

    corrected_slope =
      apparent_slope -
      optimism_slope,

    B_valid =
      n(),

    .groups =
      "drop"
  )


# ============================================================
# 10. Final summary
# ============================================================

final_summary <- mi_summary %>%

  summarise(

    apparent_C =
      median(
        apparent_C
      ),

    optimism_C =
      median(
        optimism_C
      ),

    corrected_C =
      median(
        corrected_C
      ),

    min_corrected_C =
      min(
        corrected_C
      ),

    max_corrected_C =
      max(
        corrected_C
      ),

    apparent_slope =
      median(
        apparent_slope
      ),

    optimism_slope =
      median(
        optimism_slope
      ),

    corrected_slope =
      median(
        corrected_slope
      ),

    min_corrected_slope =
      min(
        corrected_slope
      ),

    max_corrected_slope =
      max(
        corrected_slope
      )
  )


cat("\n============================================\n")
cat("C. FINAL Base6 internal validation\n")
cat("============================================\n")

print(
  final_summary
)


# ============================================================
# 11. Numerical sanity audit
# ============================================================

cat("\n============================================\n")
cat("D. Bootstrap sanity audit\n")
cat("============================================\n")

cat(
  "Non-finite optimism C:",
  sum(
    !is.finite(
      bootstrap_results$optimism_C
    )
  ),
  "\n"
)

cat(
  "Non-finite optimism slope:",
  sum(
    !is.finite(
      bootstrap_results$optimism_slope
    )
  ),
  "\n"
)

cat(
  "|optimism C| > 0.15:",
  sum(
    abs(
      bootstrap_results$optimism_C
    ) > 0.15,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "|optimism slope| > 1:",
  sum(
    abs(
      bootstrap_results$optimism_slope
    ) > 1,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 12. Save
# ============================================================

write.csv(

  mi_summary,

  "Base6_FINAL_bootstrap_MI_summary.csv",

  row.names =
    FALSE
)


write.csv(

  final_summary,

  "Base6_FINAL_bootstrap_summary.csv",

  row.names =
    FALSE
)


save(

  bootstrap_results,

  mi_summary,

  final_summary,

  base6_formula,

  nwr_knots_final,

  nwr_boundary_final,

  file =
    "NHANES_Base6_FINAL_BOOTSTRAP.RData"
)


cat("\n============================================\n")
cat("STEP 8E 完成\n")
cat("============================================\n")
