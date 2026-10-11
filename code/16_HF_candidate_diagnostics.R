# ============================================================
# STEP 6C
# Candidate model diagnostics
#
# 1. 重新构造 Core6 / Core9 / Extended10
# 2. Medication-count effect
# 3. Proportional hazards assumption
#
# 独立运行，不依赖 model_fits 中间对象
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
# 1. 读取最终MI数据
# ============================================================

load(
  "NHANES_2011_2018_MICE_FINAL.RData"
)

# 读取Step 6A冻结的NWR spline knots
load(
  "NHANES_2011_2018_STEP6A_GROUP_LASSO_FINAL.RData"
)

cat("============================================\n")
cat("STEP 6C: Candidate model diagnostics\n")
cat("============================================\n")


# ============================================================
# 2. 构造正式分析数据
# ============================================================

prepare_model_data <- function(d) {

  d <- d %>%
    mutate(

      # NHANES age >=80 top-coded at 80
      age80 =
        pmin(
          ridageyr,
          80
        ),

      # adult pulse
      pulse_rate =
        bpxpls,

      # final NWR
      nwr =
        nwr_model,

      # CKD-EPI 2021
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
  # 固定NWR spline basis
  # ----------------------------------------------------------

  nb <- ns(
    d$nwr,
    knots =
      nwr_knots,
    Boundary.knots =
      nwr_boundary,
    intercept =
      FALSE
  )


  d$nwr_ns1 <-
    nb[, 1]

  d$nwr_ns2 <-
    nb[, 2]

  d$nwr_ns3 <-
    nb[, 3]


  return(d)
}


# ============================================================
# 3. 30个MI数据
# ============================================================

completed_list_final <- lapply(

  1:imp_final$m,

  function(i) {

    prepare_model_data(
      complete(
        imp_final,
        i
      )
    )
  }
)


cat(
  "\nCompleted datasets:",
  length(
    completed_list_final
  ),
  "\n"
)


# ============================================================
# 4. 三个候选模型
# ============================================================

formulas <- list(

  # ----------------------------------------------------------
  # Core6
  # Age + BMI + RDW + NWR + BUN + Albumin
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Core9
  # Core6 + Pulse + Sodium + eGFR
  # ----------------------------------------------------------

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


  # ----------------------------------------------------------
  # Extended10
  # Core9 + medication count
  # ----------------------------------------------------------

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
# 5. 重新拟合全部模型
# ============================================================

model_fits <- list()


for (
  i in seq_along(
    completed_list_final
  )
) {

  d <-
    completed_list_final[[i]]


  for (
    model_name in names(
      formulas
    )
  ) {

    fit <-
      coxph(

        formulas[[model_name]],

        data =
          d,

        ties =
          "efron",

        x =
          TRUE,

        y =
          TRUE,

        model =
          TRUE
      )


    model_fits[[
      paste0(
        model_name,
        "_MI",
        i
      )
    ]] <-
      fit
  }
}


cat("\n============================================\n")
cat("A. Model reconstruction check\n")
cat("============================================\n")


for (
  model_name in names(
    formulas
  )
) {

  f <-
    model_fits[[
      paste0(
        model_name,
        "_MI1"
      )
    ]]

  cat(
    model_name,
    ": ",
    length(
      coef(f)
    ),
    " coefficients\n",
    sep = ""
  )
}


cat(
  "\nTotal fitted models:",
  length(
    model_fits
  ),
  "\n"
)


# ============================================================
# 6. Medication-count effect
# ============================================================

med_effect <- lapply(

  1:imp_final$m,

  function(i) {

    fit <-
      model_fits[[
        paste0(
          "Extended10_MI",
          i
        )
      ]]


    b <-
      coef(fit)[
        "med_cnt"
      ]


    se <-
      sqrt(
        vcov(fit)[
          "med_cnt",
          "med_cnt"
        ]
      )


    data.frame(

      imputation =
        i,

      beta =
        b,

      SE =
        se,

      HR =
        exp(b),

      HR_low =
        exp(
          b -
            1.96 *
            se
        ),

      HR_high =
        exp(
          b +
            1.96 *
            se
        )
    )
  }
) %>%
  dplyr::bind_rows()


cat("\n============================================\n")
cat("B. Medication-count effect across MI\n")
cat("============================================\n")


med_summary <-
  med_effect %>%
  dplyr::summarise(

    median_beta =
      median(
        beta
      ),

    min_beta =
      min(
        beta
      ),

    max_beta =
      max(
        beta
      ),

    median_HR =
      median(
        HR
      ),

    min_HR =
      min(
        HR
      ),

    max_HR =
      max(
        HR
      )
  )


print(
  med_summary
)


# ============================================================
# 7. PH assumption
# ============================================================

ph_results <-
  list()

counter <-
  1


for (
  model_name in names(
    formulas
  )
) {

  cat(
    "\nPH testing:",
    model_name,
    "\n"
  )


  for (
    i in 1:imp_final$m
  ) {

    fit <-
      model_fits[[
        paste0(
          model_name,
          "_MI",
          i
        )
      ]]


    zph <-
      cox.zph(
        fit,
        transform =
          "km",
        terms =
          TRUE,
        singledf =
          TRUE
      )


    tab <-
      as.data.frame(
        zph$table
      )


    tab$term <-
      rownames(
        tab
      )


    tab$model <-
      model_name


    tab$imputation <-
      i


    rownames(
      tab
    ) <-
      NULL


    ph_results[[counter]] <-
      tab


    counter <-
      counter +
      1
  }
}


ph_results <-
  dplyr::bind_rows(
    ph_results
  )


# ============================================================
# 8. 自动识别P值列
# ============================================================

p_col <-
  grep(

    "^p$|pvalue|p.value",

    names(
      ph_results
    ),

    ignore.case =
      TRUE,

    value =
      TRUE
  )


if (
  length(
    p_col
  ) ==
  0
) {

  stop(
    "无法识别 cox.zph 的 p-value 列。"
  )
}


p_col <-
  p_col[1]


cat(
  "\nDetected p-value column:",
  p_col,
  "\n"
)


# ============================================================
# 9. Predictor-specific PH tests
# ============================================================

ph_summary <-
  ph_results %>%

  dplyr::filter(
    term !=
      "GLOBAL"
  ) %>%

  dplyr::group_by(
    model,
    term
  ) %>%

  dplyr::summarise(

    median_p =
      median(
        .data[[p_col]],
        na.rm =
          TRUE
      ),

    min_p =
      min(
        .data[[p_col]],
        na.rm =
          TRUE
      ),

    max_p =
      max(
        .data[[p_col]],
        na.rm =
          TRUE
      ),

    pct_p_lt_005 =
      mean(
        .data[[p_col]] <
          0.05,
        na.rm =
          TRUE
      ) *
      100,

    .groups =
      "drop"
  ) %>%

  dplyr::arrange(
    model,
    median_p
  )


cat("\n============================================\n")
cat("C. PH assumption by predictor\n")
cat("============================================\n")


print(
  ph_summary,
  n = Inf
)


# ============================================================
# 10. GLOBAL PH test
# ============================================================

ph_global <-
  ph_results %>%

  dplyr::filter(
    term ==
      "GLOBAL"
  ) %>%

  dplyr::group_by(
    model
  ) %>%

  dplyr::summarise(

    median_p =
      median(
        .data[[p_col]],
        na.rm =
          TRUE
      ),

    min_p =
      min(
        .data[[p_col]],
        na.rm =
          TRUE
      ),

    max_p =
      max(
        .data[[p_col]],
        na.rm =
          TRUE
      ),

    pct_p_lt_005 =
      mean(
        .data[[p_col]] <
          0.05,
        na.rm =
          TRUE
      ) *
      100,

    .groups =
      "drop"
  )


cat("\n============================================\n")
cat("D. GLOBAL PH test\n")
cat("============================================\n")


print(
  ph_global
)


# ============================================================
# 11. 保存
# ============================================================

write.csv(
  med_effect,
  "HF_medication_count_effect_MI.csv",
  row.names =
    FALSE
)

write.csv(
  ph_summary,
  "HF_candidate_models_PH_summary.csv",
  row.names =
    FALSE
)

write.csv(
  ph_global,
  "HF_candidate_models_PH_global.csv",
  row.names =
    FALSE
)


save(

  completed_list_final,

  formulas,

  model_fits,

  med_effect,

  med_summary,

  ph_results,

  ph_summary,

  ph_global,

  nwr_knots,

  nwr_boundary,

  file =
    "NHANES_2011_2018_STEP6C_DIAGNOSTICS.RData"
)


cat("\n============================================\n")
cat("STEP 6C 完成\n")
cat("============================================\n")
