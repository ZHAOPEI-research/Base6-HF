# ============================================================
# STEP 6A FINAL
# Simplified HF prognostic model
#
# Group-LASSO Cox variable selection stability
# across 30 multiply imputed datasets
#
# grpreg version: 3.6.0
#
# SBP quadratic = one group
# NWR spline    = one group
# ============================================================

rm(list = ls())

library(mice)
library(survival)
library(splines)
library(dplyr)
library(grpreg)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load(
  "NHANES_2011_2018_MICE_FINAL.RData"
)

cat("============================================\n")
cat("STEP 6A FINAL: Group-LASSO Cox\n")
cat("============================================\n")

cat(
  "grpreg version:",
  as.character(packageVersion("grpreg")),
  "\n"
)


# ============================================================
# 1. 固定 NWR spline knots
#
# 必须冻结：
# 后续30个MI、bootstrap、外部验证均使用同一组knots
# ============================================================

d0 <- complete(
  imp_final,
  1
)

d0$nwr <- d0$nwr_model


# ------------------------------------------------------------
# ns(df=3) 相当于3个basis
#
# 为保证未来外部验证可复现，
# 显式冻结内部knots和Boundary knots
# ------------------------------------------------------------

nwr_knots <- as.numeric(
  quantile(
    d0$nwr,
    probs = c(
      1/3,
      2/3
    ),
    na.rm = TRUE
  )
)

nwr_boundary <- range(
  d0$nwr,
  na.rm = TRUE
)


cat("\n============================================\n")
cat("A. Frozen NWR spline settings\n")
cat("============================================\n")

cat(
  "Internal knots:",
  paste(
    round(nwr_knots, 6),
    collapse = ", "
  ),
  "\n"
)

cat(
  "Boundary knots:",
  paste(
    round(nwr_boundary, 6),
    collapse = ", "
  ),
  "\n"
)


# ============================================================
# 2. 构建单个MI分析数据
# ============================================================

prepare_data <- function(d) {

  d <- d %>%
    mutate(

      # Outcome
      time =
        followup_month,

      status =
        death,


      # NHANES age top-coded at 80
      age80 =
        pmin(
          ridageyr,
          80
        ),


      # Adult pulse
      pulse_rate =
        bpxpls,


      # NWR
      nwr =
        nwr_model,


      # ------------------------------------------------------
      # SBP quadratic
      # 中心化不改变模型拟合，
      # 但改善一次项与二次项的数值表现
      # ------------------------------------------------------

      sbp_c =
        sbp_mean - 130,

      sbp_c2 =
        sbp_c^2,


      # ------------------------------------------------------
      # CKD-EPI 2021 eGFR
      # ------------------------------------------------------

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


  d$nwr_ns1 <- nb[, 1]
  d$nwr_ns2 <- nb[, 2]
  d$nwr_ns3 <- nb[, 3]


  return(d)
}


# ============================================================
# 3. 生成30个MI数据集
# ============================================================

completed_list <- lapply(

  1:imp_final$m,

  function(i) {

    prepare_data(
      complete(
        imp_final,
        i
      )
    )
  }
)


cat(
  "\nCompleted datasets:",
  length(completed_list),
  "\n"
)


# ============================================================
# 4. Predictor design matrix
# ============================================================

x_vars <- c(

  # group 1
  "age80",

  # group 2
  "female",

  # group 3: SBP
  "sbp_c",
  "sbp_c2",

  # group 4
  "pulse_rate",

  # group 5
  "bmxbmi",

  # group 6
  "lbxhgb",

  # group 7
  "lbxrdw",

  # group 8: NWR
  "nwr_ns1",
  "nwr_ns2",
  "nwr_ns3",

  # group 9
  "egfr_2021",

  # group 10
  "lbxsbu",

  # group 11
  "lbxsal",

  # group 12
  "lbxsnasi",

  # group 13
  "med_cnt",

  # group 14
  "aceiarb_use",

  # group 15
  "bb_use",

  # group 16
  "mra_use"
)


# ============================================================
# 5. Group definition
# ============================================================

group <- c(

  1,

  2,

  3, 3,

  4,

  5,

  6,

  7,

  8, 8, 8,

  9,

  10,

  11,

  12,

  13,

  14,

  15,

  16
)


group_labels <- c(

  "Age",
  "Sex",
  "SBP",
  "Pulse",
  "BMI",
  "Hemoglobin",
  "RDW",
  "NWR",
  "eGFR",
  "BUN",
  "Albumin",
  "Sodium",
  "Medication_count",
  "ACEI_ARB_ARNI",
  "Beta_blocker",
  "MRA"
)


group_dictionary <- data.frame(

  group =
    1:16,

  predictor =
    group_labels,

  stringsAsFactors =
    FALSE
)


cat("\n============================================\n")
cat("B. Group definition\n")
cat("============================================\n")

print(
  group_dictionary,
  row.names = FALSE
)


# ============================================================
# 6. Design matrix检查
# ============================================================

X0 <- as.matrix(
  completed_list[[1]][
    x_vars
  ]
)


cat("\n============================================\n")
cat("C. Design matrix check\n")
cat("============================================\n")

cat(
  "Rows:",
  nrow(X0),
  "\n"
)

cat(
  "Columns:",
  ncol(X0),
  "\n"
)

cat(
  "Group vector length:",
  length(group),
  "\n"
)


if (
  ncol(X0) !=
  length(group)
) {

  stop(
    "ERROR: X columns != group length."
  )
}


if (
  any(
    !is.finite(X0)
  )
) {

  stop(
    "ERROR: Design matrix contains NA/Inf."
  )
}


# ============================================================
# 7. 固定10-fold CV
#
# 30个MI使用完全相同的fold assignment
#
# 并按照死亡状态分层
# ============================================================

set.seed(
  20260822
)


status0 <-
  completed_list[[1]]$status


fold_id <-
  rep(
    NA_integer_,
    length(status0)
  )


for (
  s in c(
    0,
    1
  )
) {

  ix <-
    which(
      status0 == s
    )


  fold_id[ix] <-
    sample(
      rep(
        1:10,
        length.out =
          length(ix)
      )
    )
}


cat("\n============================================\n")
cat("D. CV fold check\n")
cat("============================================\n")

print(
  table(
    fold_id,
    status0
  )
)


# ============================================================
# 8. 存储结果
# ============================================================

selection_min <- matrix(

  0,

  nrow =
    imp_final$m,

  ncol =
    length(
      group_labels
    ),

  dimnames =
    list(
      paste0(
        "MI",
        1:imp_final$m
      ),
      group_labels
    )
)


selection_1se <-
  selection_min


lambda_table <- data.frame(

  imputation =
    1:imp_final$m,

  lambda_min =
    NA_real_,

  lambda_1se =
    NA_real_,

  groups_min =
    NA_integer_,

  groups_1se =
    NA_integer_
)


cv_models <- vector(
  "list",
  imp_final$m
)


# ============================================================
# 9. 辅助函数：
# 把predict(type="groups")转换为group index
# ============================================================

extract_group_indices <- function(x) {

  if (
    length(x) == 0
  ) {

    return(
      integer(0)
    )
  }


  # predict(type="groups")可能返回factor
  # levels即原始group labels

  z <-
    suppressWarnings(
      as.integer(
        as.character(x)
      )
    )


  # 如果不是数字标签，则尝试factor内部code
  if (
    any(
      is.na(z)
    )
  ) {

    z <-
      as.integer(x)
  }


  z <-
    unique(
      z[
        !is.na(z)
      ]
    )


  z <-
    z[
      z >= 1 &
        z <= length(
          group_labels
        )
    ]


  sort(z)
}


# ============================================================
# 10. 30个MI分别跑Group-LASSO
# ============================================================

for (
  i in 1:imp_final$m
) {

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


  X <-
    as.matrix(
      d[
        x_vars
      ]
    )


  y <-
    Surv(
      d$time,
      d$status
    )


  # ----------------------------------------------------------
  # Group-LASSO Cox
  # ----------------------------------------------------------

  cvfit <- cv.grpsurv(

    X =
      X,

    y =
      y,

    group =
      group,

    penalty =
      "grLasso",

    nfolds =
      10,

    fold =
      fold_id,

    seed =
      20260822,

    se =
      "quick",

    trace =
      FALSE
  )


  cv_models[[i]] <-
    cvfit


  # ==========================================================
  # lambda.min
  # ==========================================================

  idx_min <-
    which.min(
      cvfit$cve
    )


  lambda_min <-
    cvfit$lambda[
      idx_min
    ]


  selected_min_raw <-
    predict(

      cvfit,

      lambda =
        lambda_min,

      type =
        "groups"
    )


  selected_groups_min <-
    extract_group_indices(
      selected_min_raw
    )


  if (
    length(
      selected_groups_min
    ) > 0
  ) {

    selection_min[
      i,
      selected_groups_min
    ] <- 1
  }


  # ==========================================================
  # lambda.1se
  #
  # largest lambda satisfying:
  #
  # CV error <= min CV error + SE(min)
  # ==========================================================

  threshold <-
    cvfit$cve[
      idx_min
    ] +
    cvfit$cvse[
      idx_min
    ]


  eligible <-
    which(
      cvfit$cve <=
        threshold
    )


  if (
    length(
      eligible
    ) == 0
  ) {

    stop(
      paste0(
        "MI ",
        i,
        ": no eligible lambda for 1SE."
      )
    )
  }


  lambda_1se <-
    max(
      cvfit$lambda[
        eligible
      ]
    )


  selected_1se_raw <-
    predict(

      cvfit,

      lambda =
        lambda_1se,

      type =
        "groups"
    )


  selected_groups_1se <-
    extract_group_indices(
      selected_1se_raw
    )


  if (
    length(
      selected_groups_1se
    ) > 0
  ) {

    selection_1se[
      i,
      selected_groups_1se
    ] <- 1
  }


  # ==========================================================
  # 保存
  # ==========================================================

  lambda_table$lambda_min[i] <-
    lambda_min


  lambda_table$lambda_1se[i] <-
    lambda_1se


  lambda_table$groups_min[i] <-
    length(
      selected_groups_min
    )


  lambda_table$groups_1se[i] <-
    length(
      selected_groups_1se
    )
}


# ============================================================
# 11. Selection frequency
# ============================================================

selection_summary <- data.frame(

  predictor =
    group_labels,

  selected_min_n =
    colSums(
      selection_min
    ),

  selected_min_pct =
    round(
      colMeans(
        selection_min
      ) *
        100,
      1
    ),

  selected_1se_n =
    colSums(
      selection_1se
    ),

  selected_1se_pct =
    round(
      colMeans(
        selection_1se
      ) *
        100,
      1
    )
)


selection_summary <-
  selection_summary %>%
  arrange(
    desc(
      selected_1se_pct
    ),
    desc(
      selected_min_pct
    )
  )


cat("\n============================================\n")
cat("E. Group-LASSO selection stability\n")
cat("============================================\n")

print(
  selection_summary,
  row.names = FALSE
)


# ============================================================
# 12. 模型大小
# ============================================================

cat("\n============================================\n")
cat("F. Selected predictor count\n")
cat("============================================\n")


cat("\nlambda.min:\n")

print(
  summary(
    lambda_table$
      groups_min
  )
)


cat("\nlambda.1se:\n")

print(
  summary(
    lambda_table$
      groups_1se
  )
)


cat(
  "\nMedian selected at lambda.min:",
  median(
    lambda_table$
      groups_min
  ),
  "\n"
)


cat(
  "Median selected at lambda.1se:",
  median(
    lambda_table$
      groups_1se
  ),
  "\n"
)


# ============================================================
# 13. 稳定性阈值
# ============================================================

stable_50 <-
  selection_summary %>%
  filter(
    selected_1se_pct >= 50
  )


stable_70 <-
  selection_summary %>%
  filter(
    selected_1se_pct >= 70
  )


stable_90 <-
  selection_summary %>%
  filter(
    selected_1se_pct >= 90
  )


cat("\n============================================\n")
cat("G. Selected >=50% at lambda.1se\n")
cat("============================================\n")

print(
  stable_50,
  row.names = FALSE
)


cat("\n============================================\n")
cat("H. Selected >=70% at lambda.1se\n")
cat("============================================\n")

print(
  stable_70,
  row.names = FALSE
)


cat("\n============================================\n")
cat("I. Selected >=90% at lambda.1se\n")
cat("============================================\n")

print(
  stable_90,
  row.names = FALSE
)


# ============================================================
# 14. 选中组合的频率
# ============================================================

selection_pattern_1se <- apply(

  selection_1se,

  1,

  function(x) {

    selected <-
      group_labels[
        x == 1
      ]


    if (
      length(
        selected
      ) == 0
    ) {

      return(
        "None"
      )
    }


    paste(
      selected,
      collapse = " + "
    )
  }
)


pattern_table <- sort(
  table(
    selection_pattern_1se
  ),
  decreasing = TRUE
)


cat("\n============================================\n")
cat("J. Most common lambda.1se models\n")
cat("============================================\n")

print(
  pattern_table
)


# ============================================================
# 15. 保存
# ============================================================

write.csv(
  selection_summary,
  "HF_groupLASSO_selection_stability_FINAL.csv",
  row.names = FALSE
)


write.csv(
  lambda_table,
  "HF_groupLASSO_lambda_by_MI_FINAL.csv",
  row.names = FALSE
)


write.csv(
  as.data.frame(
    pattern_table
  ),
  "HF_groupLASSO_model_patterns_FINAL.csv",
  row.names = FALSE
)


save(
  completed_list,
  nwr_knots,
  nwr_boundary,
  group,
  group_dictionary,
  fold_id,
  selection_min,
  selection_1se,
  selection_summary,
  lambda_table,
  pattern_table,
  cv_models,
  file =
    "NHANES_2011_2018_STEP6A_GROUP_LASSO_FINAL.RData"
)


cat("\n============================================\n")
cat("STEP 6A FINAL 完成\n")
cat("============================================\n")
