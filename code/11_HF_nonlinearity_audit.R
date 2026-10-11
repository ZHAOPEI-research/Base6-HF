# ============================================================
# STEP 5B
# Simplified HF prognostic model
# Continuous predictor nonlinearity audit
#
# Primary development analysis:
# ordinary Cox model
#
# NHANES survey-weighted model以后作为敏感性分析
#
# 每个连续变量：
#   linear model
#       vs
#   natural cubic spline df=3
#
# 在30个MI数据集分别比较
# ============================================================

rm(list = ls())

library(survival)
library(splines)
library(dplyr)
library(mice)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load(
  "NHANES_2011_2018_STEP4D_MICE50_FINAL.RData"
)

cat("============================================\n")
cat("STEP 5B: Nonlinearity audit\n")
cat("============================================\n")


# ============================================================
# 1. 从30个插补数据重新构造正式分析变量
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

      # ----------------------------------------
      # Pulse rate
      # ----------------------------------------

      pulse_rate =
        bpxpls,


      # ----------------------------------------
      # NWR
      # ----------------------------------------

      nwr =
        lbdneno /
        lbxwbcsi,


      # ----------------------------------------
      # CKD-EPI 2021 eGFR
      # ----------------------------------------

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


      # ----------------------------------------
      # NHANES age为80岁封顶
      #
      # 原始数据已经如此；
      # 明确保存一份 predictor
      # ----------------------------------------

      age80 =
        pmin(
          ridageyr,
          80
        )
    )

  completed_list[[i]] <- d
}


# ============================================================
# 2. 候选预测因子
# ============================================================

continuous_vars <- c(

  "age80",
  "sbp_mean",
  "pulse_rate",
  "bmxbmi",

  "lbxhgb",
  "lbxrdw",
  "nwr",

  "egfr_2021",
  "lbxsbu",

  "lbxsal",
  "lbxsnasi",

  "med_cnt"
)


binary_vars <- c(

  "female",
  "aceiarb_use",
  "bb_use",
  "mra_use"
)


all_predictors <- c(
  continuous_vars,
  binary_vars
)


# ============================================================
# 3. Cox公式生成函数
# ============================================================

make_formula <- function(
    target,
    nonlinear = FALSE
) {

  others <-
    setdiff(
      all_predictors,
      target
    )


  if (
    nonlinear
  ) {

    target_term <-
      paste0(
        "ns(",
        target,
        ", df=3)"
      )

  } else {

    target_term <-
      target
  }


  rhs <-
    paste(
      c(
        target_term,
        others
      ),
      collapse = " + "
    )


  as.formula(
    paste0(
      "Surv(followup_month, death) ~ ",
      rhs
    )
  )
}


# ============================================================
# 4. 每个插补数据集：
# linear vs spline
# ============================================================

nonlinear_results <- list()

counter <- 1


for (
  v in continuous_vars
) {

  cat(
    "\nTesting:",
    v,
    "\n"
  )


  form_linear <-
    make_formula(
      v,
      nonlinear = FALSE
    )


  form_spline <-
    make_formula(
      v,
      nonlinear = TRUE
    )


  for (
    i in seq_along(
      completed_list
    )
  ) {

    d <-
      completed_list[[i]]


    # ----------------------------------------
    # Linear model
    # ----------------------------------------

    fit_linear <-
      coxph(
        form_linear,
        data = d,
        ties = "efron",
        x = TRUE
      )


    # ----------------------------------------
    # Spline model
    # ----------------------------------------

    fit_spline <-
      coxph(
        form_spline,
        data = d,
        ties = "efron",
        x = TRUE
      )


    # ----------------------------------------
    # Likelihood ratio comparison
    # ----------------------------------------

    lr <-
      2 *
      (
        as.numeric(
          logLik(
            fit_spline
          )
        ) -
          as.numeric(
            logLik(
              fit_linear
            )
          )
      )


    df_diff <-
      attr(
        logLik(
          fit_spline
        ),
        "df"
      ) -
      attr(
        logLik(
          fit_linear
        ),
        "df"
      )


    p_nonlin <-
      pchisq(
        lr,
        df = df_diff,
        lower.tail = FALSE
      )


    nonlinear_results[[counter]] <-
      data.frame(

        variable =
          v,

        imputation =
          i,

        logLik_linear =
          as.numeric(
            logLik(
              fit_linear
            )
          ),

        logLik_spline =
          as.numeric(
            logLik(
              fit_spline
            )
          ),

        LR =
          lr,

        df =
          df_diff,

        p_nonlin =
          p_nonlin,

        AIC_linear =
          AIC(
            fit_linear
          ),

        AIC_spline =
          AIC(
            fit_spline
          ),

        delta_AIC =
          AIC(
            fit_spline
          ) -
          AIC(
            fit_linear
          )
      )


    counter <-
      counter + 1
  }
}


nonlinear_results <-
  bind_rows(
    nonlinear_results
  )


# ============================================================
# 5. 跨30次插补汇总
# ============================================================

nonlinear_summary <-
  nonlinear_results %>%
  group_by(
    variable
  ) %>%
  summarise(

    median_LR =
      median(
        LR
      ),

    median_p =
      median(
        p_nonlin
      ),

    min_p =
      min(
        p_nonlin
      ),

    max_p =
      max(
        p_nonlin
      ),

    pct_p_lt_005 =
      mean(
        p_nonlin < 0.05
      ) *
      100,

    median_delta_AIC =
      median(
        delta_AIC
      ),

    min_delta_AIC =
      min(
        delta_AIC
      ),

    max_delta_AIC =
      max(
        delta_AIC
      ),

    .groups = "drop"
  ) %>%
  arrange(
    median_p
  )


cat("\n============================================\n")
cat("A. Nonlinearity summary across 30 MI datasets\n")
cat("============================================\n")

print(
  nonlinear_summary,
  n = Inf
)


# ============================================================
# 6. 简单诊断分类
#
# 注意：
# 只是审计标签
# 不是变量选择规则
# ============================================================

nonlinear_summary <-
  nonlinear_summary %>%
  mutate(

    nonlinear_signal =
      case_when(

        pct_p_lt_005 >= 80 &
          median_delta_AIC <= -2 ~

          "Strong",

        pct_p_lt_005 >= 50 ~

          "Possible",

        TRUE ~

          "Weak"
      )
  )


cat("\n============================================\n")
cat("B. Nonlinearity signal\n")
cat("============================================\n")

print(
  nonlinear_summary %>%
    select(
      variable,
      median_p,
      pct_p_lt_005,
      median_delta_AIC,
      nonlinear_signal
    ),
  n = Inf
)


# ============================================================
# 7. 第一个插补数据集：
# 为每个变量生成risk-shape预测数据
#
# 用于后续画图
# ============================================================

d1 <-
  completed_list[[1]]


# ----------------------------------------
# 其他连续变量固定中位数
# 二分类变量固定最常见值
# ----------------------------------------

reference_values <- list()

for (
  v in continuous_vars
) {

  reference_values[[v]] <-
    median(
      d1[[v]],
      na.rm = TRUE
    )
}


for (
  v in binary_vars
) {

  reference_values[[v]] <-
    as.numeric(
      names(
        which.max(
          table(
            d1[[v]]
          )
        )
      )
    )
}


curve_list <- list()


for (
  v in continuous_vars
) {

  cat(
    "Creating curve:",
    v,
    "\n"
  )


  # ----------------------------------------
  # 用第1-99百分位避免极端尾部支配图形
  # ----------------------------------------

  q1 <-
    quantile(
      d1[[v]],
      0.01,
      na.rm = TRUE
    )


  q99 <-
    quantile(
      d1[[v]],
      0.99,
      na.rm = TRUE
    )


  grid <-
    seq(
      q1,
      q99,
      length.out = 100
    )


  newdat <-
    as.data.frame(
      lapply(
        all_predictors,
        function(x)
          rep(
            reference_values[[x]],
            100
          )
      )
    )


  names(
    newdat
  ) <-
    all_predictors


  newdat[[v]] <-
    grid


  fit <-
    coxph(

      make_formula(
        v,
        nonlinear = TRUE
      ),

      data =
        d1,

      ties =
        "efron",

      x =
        TRUE
    )


  lp <-
    predict(
      fit,
      newdata =
        newdat,
      type =
        "lp",
      se.fit =
        TRUE
    )


  # ----------------------------------------
  # 以中位数作为HR=1参考
  # ----------------------------------------

  ref_dat <-
    newdat[
      1,
      ,
      drop = FALSE
    ]


  ref_dat[[v]] <-
    median(
      d1[[v]],
      na.rm = TRUE
    )


  ref_lp <-
    as.numeric(
      predict(
        fit,
        newdata =
          ref_dat,
        type =
          "lp"
      )
    )


  curve_list[[v]] <-
    data.frame(

      variable =
        v,

      value =
        grid,

      logHR =
        lp$fit -
        ref_lp,

      HR =
        exp(
          lp$fit -
            ref_lp
        ),

      lower =
        exp(
          lp$fit -
            ref_lp -
            1.96 *
            lp$se.fit
        ),

      upper =
        exp(
          lp$fit -
            ref_lp +
            1.96 *
            lp$se.fit
        )
    )
}


curve_data <-
  bind_rows(
    curve_list
  )


# ============================================================
# 8. 年龄top-code特别报告
# ============================================================

cat("\n============================================\n")
cat("C. Age top-coding audit\n")
cat("============================================\n")


age80_n <-
  sum(
    d1$age80 == 80
  )


cat(
  "Age==80 (includes all true ages >=80):",
  age80_n,
  "/",
  nrow(
    d1
  ),
  "(",
  round(
    age80_n /
      nrow(
        d1
      ) *
      100,
    1
  ),
  "%)\n"
)


# ============================================================
# 9. 保存
# ============================================================

write.csv(
  nonlinear_results,
  "HF_nonlinearity_all_MI.csv",
  row.names = FALSE
)

write.csv(
  nonlinear_summary,
  "HF_nonlinearity_summary.csv",
  row.names = FALSE
)

write.csv(
  curve_data,
  "HF_nonlinearity_curves.csv",
  row.names = FALSE
)


save(
  completed_list,
  nonlinear_results,
  nonlinear_summary,
  curve_data,
  file =
    "NHANES_2011_2018_STEP5B_NONLINEARITY.RData"
)


cat("\n============================================\n")
cat("STEP 5B 完成\n")
cat("============================================\n")
