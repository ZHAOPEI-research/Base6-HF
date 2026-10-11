# ============================================================
# STEP 7A
# Bootstrap internal validation of fixed candidate models
#
# Core6 / Core9 / Extended10
#
# 200 bootstrap samples × 30 MI datasets
#
# Outputs:
# 1. Apparent C-index
# 2. Bootstrap optimism
# 3. Optimism-corrected C-index
# 4. Calibration slope optimism
#
# 注意：
# 这是固定模型结构的第一阶段内部验证。
# 最终主模型还会进行 bootstrap -> MI 的完整验证。
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
cat("STEP 7A: Bootstrap internal validation\n")
cat("============================================\n")


# ============================================================
# 1. 参数
# ============================================================

B <- 200

set.seed(
  20260822
)


# ============================================================
# 2. 数据准备函数
# ============================================================

prepare_data <- function(d) {

  d <- d %>%
    mutate(

      age80 =
        pmin(
          ridageyr,
          80
        ),

      pulse_rate =
        bpxpls,

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


  # 固定 NWR spline basis
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
# 3. 三个模型
# ============================================================

formulas <- list(

  Core6 =

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

    lbxsbu +
    lbxsal,


  Core9 =

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

    lbxsbu +
    lbxsal +

    pulse_rate +
    lbxsnasi +
    egfr_2021,


  Extended10 =

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

    lbxsbu +
    lbxsal +

    pulse_rate +
    lbxsnasi +
    egfr_2021 +

    med_cnt
)


# ============================================================
# 4. Harrell C 函数
# ============================================================

get_cindex <- function(time, status, lp) {

  out <- survival::concordance(

    Surv(
      time,
      status
    ) ~
      lp,

    reverse =
      TRUE
  )

  as.numeric(
    out$concordance
  )
}


# ============================================================
# 5. Calibration slope函数
#
# 在测试数据：
# Surv ~ LP
#
# 理想 slope = 1
# ============================================================

get_calibration_slope <- function(
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
    ) ~
      lp,

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
# 6. 储存结果
# ============================================================

all_results <- list()

result_counter <- 1


# ============================================================
# 7. 循环30个MI
# ============================================================

for (
  m in 1:imp_final$m
) {

  cat(
    "\nMI ",
    m,
    "/",
    imp_final$m,
    "\n",
    sep = ""
  )


  # ----------------------------------------------------------
  # 当前完整数据集
  # ----------------------------------------------------------

  d_original <-
    prepare_data(
      complete(
        imp_final,
        m
      )
    )


  n <-
    nrow(
      d_original
    )


  # ==========================================================
  # 每个候选模型
  # ==========================================================

  for (
    model_name in names(
      formulas
    )
  ) {

    cat(
      "  ",
      model_name,
      "\n",
      sep = ""
    )


    form <-
      formulas[[
        model_name
      ]]


    # ========================================================
    # Apparent model
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


    # 原始模型在原始数据上的 calibration slope 应≈1
    apparent_slope <-
      get_calibration_slope(

        d_original,

        lp_original
      )


    # ========================================================
    # Bootstrap循环
    # ========================================================

    for (
      b in 1:B
    ) {

      # ------------------------------------------------------
      # Bootstrap patient sample
      # ------------------------------------------------------

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


      # ------------------------------------------------------
      # 确保bootstrap样本有事件和非事件
      # 极低概率防护
      # ------------------------------------------------------

      if (
        length(
          unique(
            d_boot$death
          )
        ) <
        2
      ) {

        next
      }


      # ------------------------------------------------------
      # 在bootstrap sample重新fit
      # ------------------------------------------------------

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

          error = function(e)
            NULL
        )


      if (
        is.null(
          fit_boot
        )
      ) {

        next
      }


      # ======================================================
      # Bootstrap model在bootstrap数据：
      # training/apparent-bootstrap performance
      # ======================================================

      lp_boot_train <-
        predict(

          fit_boot,

          newdata =
            d_boot,

          type =
            "lp"
        )


      C_boot_train <-
        get_cindex(

          d_boot$followup_month,

          d_boot$death,

          lp_boot_train
        )


      slope_boot_train <-
        get_calibration_slope(

          d_boot,

          lp_boot_train
        )


      # ======================================================
      # 同一bootstrap model应用回原始MI数据：
      # test performance
      # ======================================================

      lp_boot_test <-
        predict(

          fit_boot,

          newdata =
            d_original,

          type =
            "lp"
        )


      C_boot_test <-
        get_cindex(

          d_original$followup_month,

          d_original$death,

          lp_boot_test
        )


      slope_boot_test <-
        get_calibration_slope(

          d_original,

          lp_boot_test
        )


      # ======================================================
      # optimism
      # ======================================================

      optimism_C <-
        C_boot_train -
        C_boot_test


      optimism_slope <-
        slope_boot_train -
        slope_boot_test


      all_results[[
        result_counter
      ]] <-
        data.frame(

          imputation =
            m,

          model =
            model_name,

          bootstrap =
            b,

          apparent_C =
            apparent_C,

          apparent_slope =
            apparent_slope,

          C_boot_train =
            C_boot_train,

          C_boot_test =
            C_boot_test,

          optimism_C =
            optimism_C,

          slope_boot_train =
            slope_boot_train,

          slope_boot_test =
            slope_boot_test,

          optimism_slope =
            optimism_slope
        )


      result_counter <-
        result_counter +
        1
    }
  }
}


# ============================================================
# 8. 合并结果
# ============================================================

bootstrap_results <-
  dplyr::bind_rows(
    all_results
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
# 9. 每个MI × 模型 optimism
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

    mean_optimism_C =
      mean(
        optimism_C,
        na.rm =
          TRUE
      ),

    corrected_C =
      apparent_C -
      mean_optimism_C,


    apparent_slope =
      first(
        apparent_slope
      ),

    mean_optimism_slope =
      mean(
        optimism_slope,
        na.rm =
          TRUE
      ),

    corrected_slope =
      apparent_slope -
      mean_optimism_slope,

    B_valid =
      dplyr::n(),

    .groups =
      "drop"
  )


# ============================================================
# 10. 跨30个MI汇总
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
        mean_optimism_C
      ),

    corrected_C =
      median(
        corrected_C
      ),

    corrected_C_min =
      min(
        corrected_C
      ),

    corrected_C_max =
      max(
        corrected_C
      ),


    apparent_slope =
      median(
        apparent_slope
      ),

    optimism_slope =
      median(
        mean_optimism_slope
      ),

    corrected_slope =
      median(
        corrected_slope
      ),

    corrected_slope_min =
      min(
        corrected_slope
      ),

    corrected_slope_max =
      max(
        corrected_slope
      ),

    .groups =
      "drop"
  )


cat("\n============================================\n")
cat("B. Optimism-corrected performance\n")
cat("============================================\n")

print(
  final_summary,
  n = Inf
)


# ============================================================
# 11. 简化模型相对Extended10的corrected C损失
# ============================================================

comparison <- final_summary %>%

  dplyr::select(
    model,
    corrected_C,
    corrected_slope
  )


extended_C <-
  comparison$corrected_C[
    comparison$model ==
      "Extended10"
  ]


comparison <-
  comparison %>%
  dplyr::mutate(

    delta_C_vs_Extended10 =
      corrected_C -
      extended_C
  )


cat("\n============================================\n")
cat("C. Corrected performance comparison\n")
cat("============================================\n")

print(
  comparison,
  n = Inf
)


# ============================================================
# 12. 保存
# ============================================================

write.csv(

  bootstrap_results,

  "HF_bootstrap_fixed_models_raw.csv",

  row.names =
    FALSE
)


write.csv(

  mi_summary,

  "HF_bootstrap_fixed_models_MI_summary.csv",

  row.names =
    FALSE
)


write.csv(

  final_summary,

  "HF_bootstrap_fixed_models_final_summary.csv",

  row.names =
    FALSE
)


save(

  bootstrap_results,

  mi_summary,

  final_summary,

  comparison,

  formulas,

  nwr_knots,

  nwr_boundary,

  file =
    "NHANES_2011_2018_STEP7A_BOOTSTRAP_FIXED.RData"
)


cat("\n============================================\n")
cat("STEP 7A 完成\n")
cat("============================================\n")
