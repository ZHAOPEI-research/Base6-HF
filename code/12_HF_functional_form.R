# ============================================================
# STEP 5C
# Simplified HF prognostic model
#
# 寻找简化的连续变量函数形式
#
# NWR:
# linear vs log vs quadratic vs spline
#
# SBP:
# linear vs quadratic vs spline
#
# 目的：
# 在预测性能与模型简洁性之间折中
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
  "NHANES_2011_2018_STEP4D_MICE50_FINAL.RData"
)

cat("============================================\n")
cat("STEP 5C: Functional form simplification\n")
cat("============================================\n")


# ============================================================
# 1. 构造30个completed datasets
# ============================================================

completed_list <- vector(
  "list",
  imp50$m
)


for (i in 1:imp50$m) {

  d <- complete(
    imp50,
    i
  )

  d <- d %>%
    mutate(

      pulse_rate = bpxpls,

      nwr =
        lbdneno /
        lbxwbcsi,

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
        ),

      age80 =
        pmin(
          ridageyr,
          80
        ),

      # ----------------------------------------
      # NWR simplified transformations
      # ----------------------------------------

      log_nwr =
        log(
          nwr
        ),

      nwr2 =
        nwr^2,

      # ----------------------------------------
      # SBP quadratic
      #
      # 先中心化，提高数值稳定性
      # ----------------------------------------

      sbp_c =
        sbp_mean - 130,

      sbp_c2 =
        sbp_c^2
    )

  completed_list[[i]] <- d
}


# ============================================================
# 2. 除目标变量以外的固定候选预测因子
# ============================================================

other_predictors <- c(

  "age80",
  "female",

  "pulse_rate",
  "bmxbmi",

  "lbxhgb",
  "lbxrdw",

  "egfr_2021",
  "lbxsbu",

  "lbxsal",
  "lbxsnasi",

  "med_cnt",

  "aceiarb_use",
  "bb_use",
  "mra_use"
)


other_rhs <-
  paste(
    other_predictors,
    collapse = " + "
  )


# ============================================================
# 3. NWR模型公式
# ============================================================

nwr_forms <- list(

  linear =
    as.formula(
      paste(
        "Surv(followup_month, death) ~",
        "nwr + sbp_mean +",
        other_rhs
      )
    ),

  log =
    as.formula(
      paste(
        "Surv(followup_month, death) ~",
        "log_nwr + sbp_mean +",
        other_rhs
      )
    ),

  quadratic =
    as.formula(
      paste(
        "Surv(followup_month, death) ~",
        "nwr + nwr2 + sbp_mean +",
        other_rhs
      )
    ),

  spline =
    as.formula(
      paste(
        "Surv(followup_month, death) ~",
        "ns(nwr, df=3) + sbp_mean +",
        other_rhs
      )
    )
)


# ============================================================
# 4. SBP模型公式
#
# NWR这里先统一使用最灵活的spline，
# 避免SBP函数形式比较被NWR错误指定干扰
# ============================================================

sbp_forms <- list(

  linear =
    as.formula(
      paste(
        "Surv(followup_month, death) ~",
        "sbp_mean + ns(nwr, df=3) +",
        other_rhs
      )
    ),

  quadratic =
    as.formula(
      paste(
        "Surv(followup_month, death) ~",
        "sbp_c + sbp_c2 + ns(nwr, df=3) +",
        other_rhs
      )
    ),

  spline =
    as.formula(
      paste(
        "Surv(followup_month, death) ~",
        "ns(sbp_mean, df=3) + ns(nwr, df=3) +",
        other_rhs
      )
    )
)


# ============================================================
# 5. 通用比较函数
# ============================================================

compare_forms <- function(
    forms,
    data_list,
    variable_name
) {

  out <- list()

  k <- 1


  for (
    i in seq_along(
      data_list
    )
  ) {

    d <-
      data_list[[i]]


    for (
      form_name in names(
        forms
      )
    ) {

      fit <-
        coxph(
          forms[[form_name]],
          data = d,
          ties = "efron",
          x = TRUE
        )


      out[[k]] <-
        data.frame(

          variable =
            variable_name,

          imputation =
            i,

          form =
            form_name,

          logLik =
            as.numeric(
              logLik(
                fit
              )
            ),

          df =
            attr(
              logLik(
                fit
              ),
              "df"
            ),

          AIC =
            AIC(
              fit
            )
        )


      k <-
        k + 1
    }
  }


  bind_rows(
    out
  )
}


# ============================================================
# 6. NWR比较
# ============================================================

cat("\nTesting NWR functional forms...\n")

nwr_results <-
  compare_forms(
    nwr_forms,
    completed_list,
    "nwr"
  )


# ============================================================
# 7. SBP比较
# ============================================================

cat("Testing SBP functional forms...\n")

sbp_results <-
  compare_forms(
    sbp_forms,
    completed_list,
    "sbp_mean"
  )


form_results <-
  bind_rows(
    nwr_results,
    sbp_results
  )


# ============================================================
# 8. 每个插补数据集计算相对于最佳AIC的ΔAIC
# ============================================================

form_results <-
  form_results %>%
  group_by(
    variable,
    imputation
  ) %>%
  mutate(

    best_AIC =
      min(
        AIC
      ),

    delta_from_best =
      AIC -
      best_AIC
  ) %>%
  ungroup()


# ============================================================
# 9. 汇总
# ============================================================

form_summary <-
  form_results %>%
  group_by(
    variable,
    form
  ) %>%
  summarise(

    median_AIC =
      median(
        AIC
      ),

    median_delta =
      median(
        delta_from_best
      ),

    p25_delta =
      quantile(
        delta_from_best,
        0.25
      ),

    p75_delta =
      quantile(
        delta_from_best,
        0.75
      ),

    pct_best =
      mean(
        delta_from_best <
          0.000001
      ) *
      100,

    pct_within_2 =
      mean(
        delta_from_best <= 2
      ) *
      100,

    pct_within_4 =
      mean(
        delta_from_best <= 4
      ) *
      100,

    .groups = "drop"
  ) %>%
  arrange(
    variable,
    median_delta
  )


cat("\n============================================\n")
cat("A. Functional form comparison\n")
cat("============================================\n")

print(
  form_summary,
  n = Inf
)


# ============================================================
# 10. NWR详细比较
# ============================================================

cat("\n============================================\n")
cat("B. NWR forms\n")
cat("============================================\n")

print(
  form_summary %>%
    filter(
      variable == "nwr"
    ),
  n = Inf
)


# ============================================================
# 11. SBP详细比较
# ============================================================

cat("\n============================================\n")
cat("C. SBP forms\n")
cat("============================================\n")

print(
  form_summary %>%
    filter(
      variable == "sbp_mean"
    ),
  n = Inf
)


# ============================================================
# 12. 简化选择规则
#
# 这只是辅助，不自动替代研究判断
# ============================================================

simple_form_decision <-
  data.frame(

    variable =
      c(
        "nwr",
        "sbp_mean"
      ),

    preferred_if_similar =
      c(
        "log",
        "linear"
      ),

    principle =
      c(

        paste0(
          "Prefer log(NWR) if median ΔAIC ",
          "is within ~2 of spline"
        ),

        paste0(
          "Prefer linear SBP if median ΔAIC ",
          "is within ~2 of best nonlinear model"
        )
      )
  )


cat("\n============================================\n")
cat("D. Simplification principle\n")
cat("============================================\n")

print(
  simple_form_decision,
  row.names = FALSE
)


# ============================================================
# 13. 保存
# ============================================================

write.csv(
  form_results,
  "HF_functional_forms_all_MI.csv",
  row.names = FALSE
)

write.csv(
  form_summary,
  "HF_functional_forms_summary.csv",
  row.names = FALSE
)


save(
  completed_list,
  form_results,
  form_summary,
  file =
    "NHANES_2011_2018_STEP5C_FUNCTIONAL_FORM.RData"
)


cat("\n============================================\n")
cat("STEP 5C 完成\n")
cat("============================================\n")
