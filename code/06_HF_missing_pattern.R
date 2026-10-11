# ============================================================
# STEP 4A
# HF cohort 缺失模式与分析变量分层
# ============================================================

rm(list = ls())

library(dplyr)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load("NHANES_2011_2018_STEP3B_WEIGHT_MISSING.RData")


cat("============================================\n")
cat("STEP 4A: Missingness pattern audit\n")
cat("============================================\n")


# ============================================================
# 1. 定义暂定主模型变量
#
# 注意：
# 这里只用于缺失审计
# 并不是说最后一定全部进入预测模型
# ============================================================

main_vars <- c(

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
  "lbxsnasi",

  "med_cnt",
  "aceiarb_use",
  "bb_use",
  "mra_use"
)

main_vars <- intersect(
  main_vars,
  names(hf_base2)
)


# ============================================================
# 2. 为每个变量建立 missing flag
# ============================================================

miss_dat <- hf_base2

for (v in main_vars) {

  new_name <- paste0(
    "miss_",
    v
  )

  miss_dat[[new_name]] <-
    as.integer(
      is.na(
        miss_dat[[v]]
      )
    )
}


# ============================================================
# 3. 每人缺失数量
# ============================================================

miss_names <- paste0(
  "miss_",
  main_vars
)

miss_dat <- miss_dat %>%
  mutate(

    n_missing_main =
      rowSums(
        across(
          all_of(
            miss_names
          )
        )
      )
  )


cat("\n============================================\n")
cat("A. 每人主变量缺失数\n")
cat("============================================\n")

print(
  table(
    miss_dat$n_missing_main
  )
)


# ============================================================
# 4. 最常见缺失组合
# ============================================================

missing_pattern <- miss_dat %>%
  transmute(

    pattern = apply(
      as.data.frame(
        select(
          .,
          all_of(
            miss_names
          )
        )
      ),
      1,
      paste,
      collapse = ""
    )

  ) %>%
  count(
    pattern,
    sort = TRUE
  ) %>%
  mutate(

    pct =
      round(
        n /
          sum(n) *
          100,
        1
      )
  )


cat("\n============================================\n")
cat("B. 最常见缺失模式\n")
cat("============================================\n")

print(
  head(
    missing_pattern,
    20
  )
)


cat("\n变量顺序：\n")

print(
  data.frame(
    position =
      seq_along(
        main_vars
      ),
    variable =
      main_vars
  )
)


# ============================================================
# 5. 查看实验室变量是否成簇缺失
# ============================================================

lab_vars <- c(
  "lbxhgb",
  "lbxrdw",
  "nwr",
  "egfr_2021",
  "lbxsal",
  "lbxsbu",
  "lbxsnasi"
)

lab_vars <- intersect(
  lab_vars,
  names(hf_base2)
)


lab_missing <- hf_base2 %>%
  mutate(

    lab_missing_n =
      rowSums(
        is.na(
          as.data.frame(
            select(
              .,
              all_of(
                lab_vars
              )
            )
          )
        )
      )
  )


cat("\n============================================\n")
cat("C. 实验室变量缺失聚集情况\n")
cat("============================================\n")

print(
  table(
    lab_missing$lab_missing_n
  )
)


# ============================================================
# 6. 各周期缺失
# ============================================================

cat("\n============================================\n")
cat("D. 各周期主变量缺失率\n")
cat("============================================\n")

cycle_missing <- miss_dat %>%
  group_by(
    cycle
  ) %>%
  summarise(

    n = n(),

    across(
      all_of(
        miss_names
      ),
      ~ round(
        mean(.x) *
          100,
        1
      )
    ),

    .groups = "drop"
  )

print(
  cycle_missing
)


# ============================================================
# 7. Complete vs incomplete
# 更详细审计
# ============================================================

miss_dat <- miss_dat %>%
  mutate(

    incomplete =
      as.integer(
        n_missing_main > 0
      )
  )


cat("\n============================================\n")
cat("E. Complete / incomplete 基线比较\n")
cat("============================================\n")

missing_compare2 <- miss_dat %>%
  group_by(
    incomplete
  ) %>%
  summarise(

    n = n(),

    age_mean =
      mean(
        ridageyr,
        na.rm = TRUE
      ),

    age_sd =
      sd(
        ridageyr,
        na.rm = TRUE
      ),

    female_pct =
      mean(
        female == 1,
        na.rm = TRUE
      ) *
      100,

    death_pct =
      mean(
        mortstat == 1,
        na.rm = TRUE
      ) *
      100,

    followup_median =
      median(
        permth_exm,
        na.rm = TRUE
      ),

    wtmec_median =
      median(
        wtmec8yr,
        na.rm = TRUE
      ),

    .groups = "drop"
  )

print(
  missing_compare2
)


# ============================================================
# 8. 缺失与死亡的关系
# ============================================================

cat("\n============================================\n")
cat("F. 各变量缺失 vs mortality\n")
cat("============================================\n")

death_missing <- do.call(
  rbind,
  lapply(
    main_vars,
    function(v) {

      miss <-
        is.na(
          hf_base2[[v]]
        )

      data.frame(

        variable =
          v,

        death_pct_complete =
          round(
            mean(
              hf_base2$mortstat[
                !miss
              ] ==
                1,
              na.rm = TRUE
            ) *
              100,
            1
          ),

        death_pct_missing =
          ifelse(
            sum(miss) > 0,

            round(
              mean(
                hf_base2$mortstat[
                  miss
                ] ==
                  1,
                na.rm = TRUE
              ) *
                100,
              1
            ),

            NA
          ),

        missing_n =
          sum(miss)
      )
    }
  )
)

print(
  death_missing
)


# ============================================================
# 9. 缺失与年龄的关系
# ============================================================

cat("\n============================================\n")
cat("G. 各变量缺失 vs age\n")
cat("============================================\n")

age_missing <- do.call(
  rbind,
  lapply(
    main_vars,
    function(v) {

      miss <-
        is.na(
          hf_base2[[v]]
        )

      data.frame(

        variable =
          v,

        age_complete =
          round(
            mean(
              hf_base2$ridageyr[
                !miss
              ],
              na.rm = TRUE
            ),
            1
          ),

        age_missing =
          ifelse(
            sum(miss) > 0,

            round(
              mean(
                hf_base2$ridageyr[
                  miss
                ],
                na.rm = TRUE
              ),
              1
            ),

            NA
          ),

        missing_n =
          sum(miss)
      )
    }
  )
)

print(
  age_missing
)


# ============================================================
# 10. 缺失与性别关系
# ============================================================

cat("\n============================================\n")
cat("H. 各变量缺失 vs sex\n")
cat("============================================\n")

sex_missing <- do.call(
  rbind,
  lapply(
    main_vars,
    function(v) {

      miss <-
        is.na(
          hf_base2[[v]]
        )

      data.frame(

        variable =
          v,

        female_pct_complete =
          round(
            mean(
              hf_base2$female[
                !miss
              ] ==
                1,
              na.rm = TRUE
            ) *
              100,
            1
          ),

        female_pct_missing =
          ifelse(
            sum(miss) > 0,

            round(
              mean(
                hf_base2$female[
                  miss
                ] ==
                  1,
                na.rm = TRUE
              ) *
                100,
              1
            ),

            NA
          ),

        missing_n =
          sum(miss)
      )
    }
  )
)

print(
  sex_missing
)


# ============================================================
# 11. 分析变量层级
# ============================================================

analysis_tiers <- list(

  # ----------------------------------------------------------
  # Tier 1:
  # 可用于全739人主分析框架
  # ----------------------------------------------------------

  main_full_sample = c(

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
    "lbxsnasi",

    "med_cnt",
    "aceiarb_use",
    "bb_use",
    "mra_use"
  ),


  # ----------------------------------------------------------
  # Tier 2:
  # 功能状态扩展模型
  # 因缺失约54%
  # ----------------------------------------------------------

  functional_extension = c(

    "walk_tol"
  ),


  # ----------------------------------------------------------
  # Tier 3:
  # 仅2015-2018 / 特殊周期变量
  # ----------------------------------------------------------

  limited_cycles = c(

    "lbxhscrp"
  ),


  # ----------------------------------------------------------
  # Tier 4:
  # fasting subsample
  # 需要对应 subsample weight
  # ----------------------------------------------------------

  fasting_subsample = c(

    "lbxglu",
    "lbxtr",
    "lbdldl"
  )
)


cat("\n============================================\n")
cat("I. 分析变量分层\n")
cat("============================================\n")

print(
  analysis_tiers
)


# ============================================================
# 12. 保存
# ============================================================

save(
  miss_dat,
  missing_pattern,
  cycle_missing,
  missing_compare2,
  death_missing,
  age_missing,
  sex_missing,
  analysis_tiers,
  file =
    "NHANES_2011_2018_STEP4A_MISSING_PATTERN.RData"
)


cat("\n============================================\n")
cat("STEP 4A 完成\n")
cat("============================================\n")
