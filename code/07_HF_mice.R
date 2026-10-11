# ============================================================
# STEP 4B revised
# NHANES HF cohort 多重插补
#
# 原则：
# 1. 原始测量变量插补
# 2. NWR / eGFR 在插补后重新计算
# 3. outcome / follow-up / survey design 不插补
# ============================================================

rm(list = ls())

library(dplyr)
library(mice)
library(survival)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load("NHANES_2011_2018_STEP3B_WEIGHT_MISSING.RData")


cat("============================================\n")
cat("STEP 4B revised: Multiple Imputation\n")
cat("============================================\n")


# ============================================================
# 1. 工作数据
# ============================================================

mi_dat <- hf_base2


# ============================================================
# 2. Survival variables
# ============================================================

mi_dat <- mi_dat %>%
  mutate(

    death =
      as.integer(
        mortstat == 1
      ),

    followup_month =
      permth_exm
  )


# ============================================================
# 3. Nelson-Aalen cumulative hazard
# ============================================================

mi_dat$nelson_aalen <-
  mice::nelsonaalen(
    mi_dat,
    timevar = "followup_month",
    statusvar = "death"
  )


cat("\n死亡状态：\n")

print(
  table(
    mi_dat$death
  )
)


cat("\n随访时间：\n")

print(
  summary(
    mi_dat$followup_month
  )
)


cat("\nNelson-Aalen：\n")

print(
  summary(
    mi_dat$nelson_aalen
  )
)


# ============================================================
# 4. 类型
# ============================================================

mi_dat$cycle <-
  factor(
    mi_dat$cycle
  )

mi_dat$female <-
  as.integer(
    mi_dat$female
  )

mi_dat$death <-
  as.integer(
    mi_dat$death
  )

mi_dat$aceiarb_use <-
  as.integer(
    mi_dat$aceiarb_use
  )

mi_dat$bb_use <-
  as.integer(
    mi_dat$bb_use
  )

mi_dat$mra_use <-
  as.integer(
    mi_dat$mra_use
  )


# ============================================================
# 5. MICE输入变量
#
# 注意：
# 不直接放 nwr
# 不直接放 egfr_2021
#
# 改放：
# lbdneno
# lbxwbcsi
# scr_mgdl
# ============================================================

mice_vars <- c(

  "seqn",

  "cycle",
  "ridageyr",
  "female",

  "sbp_mean",
  "bpxpls",
  "bmxbmi",

  "lbxhgb",
  "lbxrdw",

  # NWR组成变量
  "lbdneno",
  "lbxwbcsi",

  # eGFR组成变量
  "scr_mgdl",

  "lbxsal",
  "lbxsbu",
  "lbxsnasi",

  "med_cnt",
  "aceiarb_use",
  "bb_use",
  "mra_use",

  "death",
  "followup_month",
  "nelson_aalen",

  "wtmec8yr",
  "sdmvstra",
  "sdmvpsu"
)


mice_vars <-
  intersect(
    mice_vars,
    names(mi_dat)
  )


mi_input <-
  mi_dat %>%
  select(
    all_of(
      mice_vars
    )
  )


# ============================================================
# 6. Missingness
# ============================================================

cat("\n============================================\n")
cat("A. MICE input missingness\n")
cat("============================================\n")

missing_before <-
  data.frame(

    variable =
      names(
        mi_input
      ),

    missing =
      sapply(
        mi_input,
        function(x)
          sum(
            is.na(x)
          )
      ),

    pct =
      round(
        sapply(
          mi_input,
          function(x)
            mean(
              is.na(x)
            ) *
            100
        ),
        2
      )
  )

print(
  missing_before[
    missing_before$missing > 0,
  ],
  row.names = FALSE
)


# ============================================================
# 7. MICE初始化
# ============================================================

ini <-
  mice(
    mi_input,
    maxit = 0,
    printFlag = FALSE
  )

method <-
  ini$method

pred <-
  ini$predictorMatrix


# ============================================================
# 8. 指定插补变量
# ============================================================

method[] <- ""


impute_continuous <- c(

  "sbp_mean",
  "bpxpls",
  "bmxbmi",

  "lbxhgb",
  "lbxrdw",

  "lbdneno",
  "lbxwbcsi",

  "scr_mgdl",

  "lbxsal",
  "lbxsbu",
  "lbxsnasi"
)


impute_continuous <-
  intersect(
    impute_continuous,
    names(
      method
    )
  )


method[
  impute_continuous
] <- "pmm"


cat("\n============================================\n")
cat("B. Imputation methods\n")
cat("============================================\n")

print(
  method[
    method != ""
  ]
)


# ============================================================
# 9. Predictor matrix
# ============================================================

pred[,] <- 1

diag(pred) <- 0


# ============================================================
# 10. SEQN完全排除
# ============================================================

if (
  "seqn" %in%
  colnames(pred)
) {

  pred[, "seqn"] <- 0
  pred["seqn", ] <- 0
}


# ============================================================
# 11. 不插补变量
# ============================================================

do_not_impute <- c(

  "cycle",
  "ridageyr",
  "female",

  "med_cnt",
  "aceiarb_use",
  "bb_use",
  "mra_use",

  "death",
  "followup_month",
  "nelson_aalen",

  "wtmec8yr",
  "sdmvstra",
  "sdmvpsu"
)


do_not_impute <-
  intersect(
    do_not_impute,
    names(
      method
    )
  )


method[
  do_not_impute
] <- ""


# ============================================================
# 12. Survival auxiliary information
#
# death + Nelson-Aalen作为predictors
# followup_month不再重复参与
# ============================================================

if (
  "followup_month" %in%
  colnames(pred)
) {

  pred[
    impute_continuous,
    "followup_month"
  ] <- 0
}


if (
  "death" %in%
  colnames(pred)
) {

  pred[
    impute_continuous,
    "death"
  ] <- 1
}


if (
  "nelson_aalen" %in%
  colnames(pred)
) {

  pred[
    impute_continuous,
    "nelson_aalen"
  ] <- 1
}


# ============================================================
# 13. Survey design辅助变量
# ============================================================

survey_aux <- c(
  "cycle",
  "wtmec8yr",
  "sdmvstra",
  "sdmvpsu"
)


survey_aux <-
  intersect(
    survey_aux,
    colnames(
      pred
    )
  )


for (
  v in survey_aux
) {

  pred[
    impute_continuous,
    v
  ] <- 1
}


# ============================================================
# 14. 正式MICE
# ============================================================

set.seed(
  20260819
)


imp <-
  mice(

    mi_input,

    m = 30,

    maxit = 20,

    method =
      method,

    predictorMatrix =
      pred,

    seed =
      20260819,

    printFlag =
      TRUE
  )


# ============================================================
# 15. 查看logged events
# ============================================================

cat("\n============================================\n")
cat("C. Logged events\n")
cat("============================================\n")

if (
  is.null(
    imp$loggedEvents
  )
) {

  cat(
    "No logged events.\n"
  )

} else {

  print(
    imp$loggedEvents
  )
}


# ============================================================
# 16. 第1个完整数据集
# ============================================================

completed1 <-
  complete(
    imp,
    1
  )


# ============================================================
# 17. 重新计算NWR
# ============================================================

completed1 <-
  completed1 %>%
  mutate(

    nwr =
      case_when(

        !is.na(lbdneno) &
          !is.na(lbxwbcsi) &
          lbxwbcsi > 0 ~

          lbdneno /
          lbxwbcsi,

        TRUE ~
          NA_real_
      )
  )


# ============================================================
# 18. 重新计算CKD-EPI 2021 eGFR
# ============================================================

completed1 <-
  completed1 %>%
  mutate(

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
        scr_mgdl /
          kappa,
        1
      )^alpha *

      pmax(
        scr_mgdl /
          kappa,
        1
      )^(-1.200) *

      (0.9938^ridageyr) *

      ifelse(
        female == 1,
        1.012,
        1
      )
  )


# ============================================================
# 19. 检查插补后是否仍有缺失
# ============================================================

cat("\n============================================\n")
cat("D. Completed dataset missingness\n")
cat("============================================\n")


check_vars <- c(
  impute_continuous,
  "nwr",
  "egfr_2021"
)


print(
  sapply(
    completed1[
      check_vars
    ],
    function(x)
      sum(
        is.na(x)
      )
  )
)


# ============================================================
# 20. 数学一致性检查
# ============================================================

cat("\n============================================\n")
cat("E. Derived-variable consistency\n")
cat("============================================\n")


cat(
  "NWR range:\n"
)

print(
  summary(
    completed1$nwr
  )
)


cat(
  "\nNWR >1:",
  sum(
    completed1$nwr > 1,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\neGFR range:\n"
)

print(
  summary(
    completed1$egfr_2021
  )
)


# ============================================================
# 21. 插补值分布
# ============================================================

cat("\n============================================\n")
cat("F. Imputed values summary\n")
cat("============================================\n")


for (
  v in impute_continuous
) {

  cat(
    "\n--- ",
    v,
    " ---\n",
    sep = ""
  )

  print(
    summary(
      as.numeric(
        unlist(
          imp$imp[[v]]
        )
      )
    )
  )
}


# ============================================================
# 22. 保存
# ============================================================

save(
  imp,
  mi_input,
  method,
  pred,
  completed1,
  file =
    "NHANES_2011_2018_STEP4B_MICE_REVISED.RData"
)


cat("\n============================================\n")
cat("STEP 4B revised 完成\n")
cat("============================================\n")
