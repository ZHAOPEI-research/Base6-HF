# ============================================================
# STEP 7D
# Bootstrap comparison of nested simplified models
#
# Base4:
# Age + RDW + NWR + eGFR
#
# Base5:
# Base4 + Medication count
#
# Base6:
# Base5 + BMI
#
# Candidate7:
# Base6 + Albumin
#
# 200 bootstrap × 30 MI
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
  "NHANES_2011_2018_STEP6A_GROUP_LASSO_FINAL.RData"
)

cat("============================================\n")
cat("STEP 7D: Bootstrap 4/5/6/7-variable models\n")
cat("============================================\n")

B <- 200

set.seed(20260822)


# ============================================================
# 1. Data preparation
# ============================================================

prepare_data <- function(d) {

  d <- d %>%
    mutate(

      age80 =
        pmin(ridageyr, 80),

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


  # Frozen NWR spline
  nb <- splines::ns(

    d$nwr,

    knots =
      nwr_knots,

    Boundary.knots =
      nwr_boundary,

    intercept =
      FALSE
  )


  d$nwr_ns1 <- nb[, 1]
  d$nwr_ns2 <- nb[, 2]
  d$nwr_ns3 <- nb[, 3]

  d
}


# ============================================================
# 2. Models
# ============================================================

models <- list(

  Base4 =
    Surv(
      followup_month,
      death
    ) ~
    age80 +
    lbxrdw +
    nwr_ns1 +
    nwr_ns2 +
    nwr_ns3 +
    egfr_2021,


  Base5_Med =
    Surv(
      followup_month,
      death
    ) ~
    age80 +
    lbxrdw +
    nwr_ns1 +
    nwr_ns2 +
    nwr_ns3 +
    egfr_2021 +
    med_cnt,


  Base6_Med_BMI =
    Surv(
      followup_month,
      death
    ) ~
    age80 +
    lbxrdw +
    nwr_ns1 +
    nwr_ns2 +
    nwr_ns3 +
    egfr_2021 +
    med_cnt +
    bmxbmi,


  Candidate7 =
    Surv(
      followup_month,
      death
    ) ~
    age80 +
    lbxrdw +
    nwr_ns1 +
    nwr_ns2 +
    nwr_ns3 +
    egfr_2021 +
    med_cnt +
    bmxbmi +
    lbxsal
)


# ============================================================
# 3. C-index
# ============================================================

get_cindex <- function(
    time,
    status,
    lp
) {

  z <-
    survival::concordance(

      Surv(
        time,
        status
      ) ~ lp,

      reverse =
        TRUE
    )

  as.numeric(
    z$concordance
  )
}


# ============================================================
# 4. Calibration slope
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


  fit <-
    coxph(

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
# 5. Results
# ============================================================

bootstrap_results <- list()

counter <- 1


# ============================================================
# 6. MI loop
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


  d_original <-
    prepare_data(
      complete(
        imp_final,
        m
      )
    )


  n <- nrow(
    d_original
  )


  for (model_name in names(models)) {

    cat(
      "  ",
      model_name,
      "\n",
      sep = ""
    )


    form <-
      models[[model_name]]


    # ========================================================
    # Original / apparent model
    # ========================================================

    fit_original <-
      coxph(

        form,

        data =
          d_original,

        ties =
          "efron",

        x =
          TRUE,

        y =
          TRUE
      )


    lp_original <-
      predict(
        fit_original,
        newdata =
          d_original,
        type =
          "lp"
      )


    apparent_C <-
      get_cindex(

        d_original$followup_month,

        d_original$death,

        lp_original
      )


    apparent_slope <-
      get_slope(
        d_original,
        lp_original
      )


    # ========================================================
    # Bootstrap
    # ========================================================

    for (b in 1:B) {

      ix <-
        sample.int(

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


      if (
        length(
          unique(
            d_boot$death
          )
        ) < 2
      ) {

        next
      }


      fit_boot <-
        tryCatch(

          coxph(

            form,

            data =
              d_boot,

            ties =
              "efron",

            x =
              TRUE,

            y =
              TRUE
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


      # ------------------------------------------------------
      # Bootstrap training performance
      # ------------------------------------------------------

      lp_train <-
        predict(

          fit_boot,

          newdata =
            d_boot,

          type =
            "lp"
        )


      C_train <-
        get_cindex(

          d_boot$followup_month,

          d_boot$death,

          lp_train
        )


      slope_train <-
        get_slope(
          d_boot,
          lp_train
        )


      # ------------------------------------------------------
      # Bootstrap model tested on original data
      # ------------------------------------------------------

      lp_test <-
        predict(

          fit_boot,

          newdata =
            d_original,

          type =
            "lp"
        )


      C_test <-
        get_cindex(

          d_original$followup_month,

          d_original$death,

          lp_test
        )


      slope_test <-
        get_slope(
          d_original,
          lp_test
        )


      # ------------------------------------------------------
      # Save
      # ------------------------------------------------------

      bootstrap_results[[counter]] <-
        data.frame(

          imputation =
            m,

          model =
            model_name,

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
        counter + 1
    }
  }
}


# ============================================================
# 7. Combine
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
      model
    )
)


# ============================================================
# 8. MI-level summary
# ============================================================

mi_summary <-
  bootstrap_results %>%

  dplyr::group_by(
    imputation,
    model
  ) %>%

  dplyr::summarise(

    apparent_C =
      first(
        apparent_C
      ),

    optimism_C =
      mean(
        optimism_C,
        na.rm =
          TRUE
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
        na.rm =
          TRUE
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
# 9. Overall summary
# ============================================================

final_summary <-
  mi_summary %>%

  dplyr::group_by(
    model
  ) %>%

  dplyr::summarise(

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
      ),

    .groups =
      "drop"
  )


cat("\n============================================\n")
cat("B. Optimism-corrected 4/5/6/7-model comparison\n")
cat("============================================\n")

print(
  final_summary,
  n = Inf
)


# ============================================================
# 10. Incremental corrected C
# ============================================================

wide <-
  final_summary %>%

  dplyr::select(
    model,
    corrected_C,
    corrected_slope
  ) %>%

  tidyr::pivot_wider(

    names_from =
      model,

    values_from =
      c(
        corrected_C,
        corrected_slope
      )
  )


incremental <- data.frame(

  comparison =
    c(
      "Base5 vs Base4",
      "Base6 vs Base5",
      "Candidate7 vs Base6"
    ),

  delta_corrected_C =
    c(

      wide$corrected_C_Base5_Med -
        wide$corrected_C_Base4,

      wide$corrected_C_Base6_Med_BMI -
        wide$corrected_C_Base5_Med,

      wide$corrected_C_Candidate7 -
        wide$corrected_C_Base6_Med_BMI
    )
)


cat("\n============================================\n")
cat("C. Incremental optimism-corrected C-index\n")
cat("============================================\n")

print(
  incremental,
  row.names = FALSE
)


# ============================================================
# 11. Save
# ============================================================

write.csv(

  bootstrap_results,

  "HF_bootstrap_4567_raw.csv",

  row.names =
    FALSE
)


write.csv(

  final_summary,

  "HF_bootstrap_4567_summary.csv",

  row.names =
    FALSE
)


write.csv(

  incremental,

  "HF_bootstrap_4567_incremental.csv",

  row.names =
    FALSE
)


save(

  bootstrap_results,

  mi_summary,

  final_summary,

  incremental,

  models,

  file =
    "NHANES_2011_2018_STEP7D_BOOTSTRAP_4567.RData"
)


cat("\n============================================\n")
cat("STEP 7D 完成\n")
cat("============================================\n")
