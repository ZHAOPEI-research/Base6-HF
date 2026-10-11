# ============================================================
# STEP 3B
# HF队列随访修正 + 缺失模式 + Survey weight 审计
# ============================================================

rm(list = ls())

library(dplyr)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load("NHANES_2011_2018_STEP3A_HF_BASE.RData")


cat("============================================\n")
cat("STEP 3B\n")
cat("随访、缺失与 Survey weight 审计\n")
cat("============================================\n")


# ============================================================
# 1. 先检查 permth_exm == 0
# ============================================================

cat("\n============================================\n")
cat("A. PERMTH_EXM == 0 检查\n")
cat("============================================\n")

permth_zero <- dat %>%
  filter(
    ridageyr >= 20,
    hf == 1,
    !is.na(permth_exm),
    permth_exm == 0
  ) %>%
  select(
    seqn,
    cycle,
    ridageyr,
    sex,
    mortstat,
    permth_int,
    permth_exm,
    wtmec2yr
  )

print(permth_zero)

cat(
  "\nHF患者中 permth_exm == 0 人数:",
  nrow(permth_zero),
  "\n"
)


# ============================================================
# 2. 重建 HF 基础队列
#
# 改为 permth_exm >= 0
# ============================================================

hf_base2 <- dat %>%
  filter(
    ridageyr >= 20,
    hf == 1,
    !is.na(permth_exm),
    permth_exm >= 0,
    !is.na(wtmec2yr),
    wtmec2yr > 0
  )


cat("\n============================================\n")
cat("B. 修正后的 HF base cohort\n")
cat("============================================\n")

cat(
  "HF base n =",
  nrow(hf_base2),
  "\n"
)

print(
  hf_base2 %>%
    count(cycle)
)

cat("\n死亡状态：\n")

print(
  table(
    hf_base2$mortstat,
    useNA = "ifany"
  )
)


# ============================================================
# 3. 创建 2011-2018 MEC 8-year weight
#
# NHANES官方规则：
# 4个连续2-year cycles
# WTMEC8YR = WTMEC2YR / 4
# ============================================================

hf_base2 <- hf_base2 %>%
  mutate(
    wtmec8yr =
      wtmec2yr / 4
  )


cat("\n============================================\n")
cat("C. WTMEC8YR 检查\n")
cat("============================================\n")

print(
  hf_base2 %>%
    group_by(cycle) %>%
    summarise(
      n = n(),
      min_wt = min(wtmec8yr),
      median_wt = median(wtmec8yr),
      max_wt = max(wtmec8yr),
      sum_wt = sum(wtmec8yr),
      .groups = "drop"
    )
)


# ============================================================
# 4. Main-model candidate variables
#
# 暂时不包含：
# walk_tol
# hsCRP
# fasting glucose / TG / LDL
#
# 因为这些存在结构性或subsample问题
# ============================================================

main_candidate_vars <- c(

  # outcome / survey
  "seqn",
  "cycle",
  "mortstat",
  "permth_exm",
  "wtmec8yr",
  "sdmvstra",
  "sdmvpsu",

  # demographic
  "ridageyr",
  "female",

  # examination
  "sbp_mean",
  "dbp_mean",
  "bpxpls",
  "bmxbmi",

  # CBC
  "lbxhgb",
  "lbxrdw",
  "lbxwbcsi",
  "lbdneno",
  "nwr",

  # standard biochemistry
  "scr_mgdl",
  "egfr_2021",
  "lbxsal",
  "lbxsbu",
  "lbxsnasi",

  # medication
  "med_cnt",
  "aceiarb_use",
  "bb_use",
  "mra_use"
)

main_candidate_vars <-
  intersect(
    main_candidate_vars,
    names(hf_base2)
  )


# ============================================================
# 5. 缺失率
# ============================================================

missing_main <- data.frame(

  variable =
    main_candidate_vars,

  n_missing =
    sapply(
      hf_base2[
        main_candidate_vars
      ],
      function(x)
        sum(is.na(x))
    ),

  pct_missing =
    round(
      sapply(
        hf_base2[
          main_candidate_vars
        ],
        function(x)
          mean(is.na(x)) * 100
      ),
      2
    )
)

missing_main <-
  missing_main %>%
  arrange(
    desc(
      pct_missing
    )
  )


cat("\n============================================\n")
cat("D. Main candidate 缺失率\n")
cat("============================================\n")

print(
  missing_main,
  row.names = FALSE
)


# ============================================================
# 6. 每个人缺几个关键变量
# ============================================================

model_predictors <- c(
  "ridageyr",
  "female",
  "sbp_mean",
  "bpxpls",
  "bmxbmi",
  "lbxhgb",
  "lbxrdw",
  "nwr",
  "egfr_2021",
  "lbxsal",
  "lbxsbu",
  "lbxsnasi"
)

model_predictors <-
  intersect(
    model_predictors,
    names(hf_base2)
  )


hf_base2 <-
  hf_base2 %>%
  mutate(

    n_missing_model =
      rowSums(
        is.na(
          as.data.frame(
            select(
              .,
              all_of(
                model_predictors
              )
            )
          )
        )
      )
  )


cat("\n============================================\n")
cat("E. 每人缺失变量数量\n")
cat("============================================\n")

print(
  table(
    hf_base2$n_missing_model,
    useNA = "ifany"
  )
)


# ============================================================
# 7. Complete-case数量
#
# 这里只统计，不作为最终方案
# ============================================================

cc_flag <-
  complete.cases(
    hf_base2[
      model_predictors
    ]
  )


cat("\n============================================\n")
cat("F. Complete-case 情况\n")
cat("============================================\n")

cat(
  "完整病例:",
  sum(cc_flag),
  "\n"
)

cat(
  "非完整病例:",
  sum(!cc_flag),
  "\n"
)

cat(
  "完整病例比例:",
  round(
    mean(cc_flag) * 100,
    1
  ),
  "%\n"
)


# ============================================================
# 8. Complete-case vs incomplete基本比较
#
# 查看缺失是否明显集中在某些患者
# ============================================================

hf_base2 <- hf_base2 %>%
  mutate(
    complete_main =
      as.integer(
        cc_flag
      )
  )


missing_compare <-
  hf_base2 %>%
  group_by(
    complete_main
  ) %>%
  summarise(

    n = n(),

    age_mean =
      mean(
        ridageyr,
        na.rm = TRUE
      ),

    female_pct =
      mean(
        female == 1,
        na.rm = TRUE
      ) * 100,

    death_pct =
      mean(
        mortstat == 1,
        na.rm = TRUE
      ) * 100,

    median_followup =
      median(
        permth_exm,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


cat("\n============================================\n")
cat("G. Complete vs incomplete\n")
cat("============================================\n")

print(
  missing_compare
)


# ============================================================
# 9. Walk_tol单独审计
# ============================================================

cat("\n============================================\n")
cat("H. Walk_tol 子样本\n")
cat("============================================\n")

walk_sub <- hf_base2 %>%
  filter(
    !is.na(
      walk_tol
    )
  )

cat(
  "HF总人数:",
  nrow(
    hf_base2
  ),
  "\n"
)

cat(
  "Walk_tol可用:",
  nrow(
    walk_sub
  ),
  "\n"
)

cat(
  "可用比例:",
  round(
    nrow(walk_sub) /
      nrow(hf_base2) *
      100,
    1
  ),
  "%\n"
)

print(
  table(
    walk_sub$walk_tol
  )
)


# ============================================================
# 10. Fasting变量提醒
# ============================================================

cat("\n============================================\n")
cat("I. Fasting变量可用人数\n")
cat("============================================\n")

fasting_vars <- c(
  "lbxglu",
  "lbxtr",
  "lbdldl"
)

fasting_vars <-
  intersect(
    fasting_vars,
    names(hf_base2)
  )


for (
  v in fasting_vars
) {

  cat(
    v,
    ": ",
    sum(
      !is.na(
        hf_base2[[v]]
      )
    ),
    "/",
    nrow(
      hf_base2
    ),
    "\n",
    sep = ""
  )
}


# ============================================================
# 11. 保存
# ============================================================

save(
  hf_base2,
  missing_main,
  missing_compare,
  walk_sub,
  permth_zero,
  file =
    "NHANES_2011_2018_STEP3B_WEIGHT_MISSING.RData"
)


cat("\n============================================\n")
cat("STEP 3B 完成\n")
cat("============================================\n")
