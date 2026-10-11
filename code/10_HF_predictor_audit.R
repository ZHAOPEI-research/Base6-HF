# ============================================================
# STEP 5A
# Simplified HF prognostic model
# Candidate predictor audit
#
# 目的：
# 1. 构建每个插补数据集中的正式候选预测变量
# 2. 明确 BPXPLS = pulse rate
# 3. 检查分布
# 4. 检查相关性
# 5. 检查共线性 / 冗余
#
# 暂不进行变量选择
# ============================================================

rm(list = ls())

library(mice)
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
cat("STEP 5A: Predictor audit\n")
cat("============================================\n")


# ============================================================
# 1. 建立30个completed datasets
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


  # ----------------------------------------------------------
  # Pulse rate
  #
  # BPXPLS = 60-second pulse in adult NHANES
  # 不命名为 BPXCHR / heart rate
  # ----------------------------------------------------------

  d <- d %>%
    mutate(

      pulse_rate =
        bpxpls
    )


  # ----------------------------------------------------------
  # NWR
  # ----------------------------------------------------------

  d <- d %>%
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


  # ----------------------------------------------------------
  # CKD-EPI 2021 eGFR
  # ----------------------------------------------------------

  d <- d %>%
    mutate(

      egfr_kappa =
        ifelse(
          female == 1,
          0.7,
          0.9
        ),

      egfr_alpha =
        ifelse(
          female == 1,
          -0.241,
          -0.302
        ),

      egfr_2021 =

        142 *

        pmin(
          scr_mgdl / egfr_kappa,
          1
        )^egfr_alpha *

        pmax(
          scr_mgdl / egfr_kappa,
          1
        )^(-1.200) *

        (0.9938^ridageyr) *

        ifelse(
          female == 1,
          1.012,
          1
        )
    )


  completed_list[[i]] <- d
}


cat(
  "Completed datasets:",
  length(completed_list),
  "\n"
)


# ============================================================
# 2. 正式候选预测变量池
# ============================================================

candidate_vars <- c(

  # Demographic
  "ridageyr",
  "female",

  # Vital signs
  "sbp_mean",
  "pulse_rate",

  # Anthropometric
  "bmxbmi",

  # Hematology
  "lbxhgb",
  "lbxrdw",
  "nwr",

  # Renal
  "egfr_2021",
  "lbxsbu",

  # Biochemistry
  "lbxsal",
  "lbxsnasi",

  # Treatment / disease burden
  "med_cnt",
  "aceiarb_use",
  "bb_use",
  "mra_use"
)


cat("\n============================================\n")
cat("A. Candidate predictors\n")
cat("============================================\n")

print(
  candidate_vars
)


# ============================================================
# 3. 变量名称映射
# ============================================================

variable_dictionary <- data.frame(

  variable =
    candidate_vars,

  label = c(

    "Age",
    "Female sex",

    "Mean systolic blood pressure",
    "Pulse rate",

    "Body mass index",

    "Hemoglobin",
    "Red cell distribution width",
    "Neutrophil-to-WBC ratio",

    "eGFR (CKD-EPI 2021)",
    "Blood urea nitrogen",

    "Serum albumin",
    "Serum sodium",

    "Number of prescription medications",
    "ACEI/ARB/ARNI use",
    "Beta-blocker use",
    "MRA use"
  ),

  stringsAsFactors = FALSE
)


print(
  variable_dictionary,
  row.names = FALSE
)


# ============================================================
# 4. 第1个插补数据集分布
#
# 先用于临床/数据范围审计
# ============================================================

d1 <- completed_list[[1]]


continuous_vars <- c(

  "ridageyr",
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


cat("\n============================================\n")
cat("B. Continuous predictor distributions\n")
cat("============================================\n")


distribution_table <- do.call(

  rbind,

  lapply(
    continuous_vars,
    function(v) {

      x <- d1[[v]]

      data.frame(

        variable = v,

        n =
          sum(
            !is.na(x)
          ),

        mean =
          mean(
            x,
            na.rm = TRUE
          ),

        sd =
          sd(
            x,
            na.rm = TRUE
          ),

        min =
          min(
            x,
            na.rm = TRUE
          ),

        p01 =
          as.numeric(
            quantile(
              x,
              0.01,
              na.rm = TRUE
            )
          ),

        p25 =
          as.numeric(
            quantile(
              x,
              0.25,
              na.rm = TRUE
            )
          ),

        median =
          median(
            x,
            na.rm = TRUE
          ),

        p75 =
          as.numeric(
            quantile(
              x,
              0.75,
              na.rm = TRUE
            )
          ),

        p99 =
          as.numeric(
            quantile(
              x,
              0.99,
              na.rm = TRUE
            )
          ),

        max =
          max(
            x,
            na.rm = TRUE
          )
      )
    }
  )
)


print(
  distribution_table,
  row.names = FALSE
)


# ============================================================
# 5. Binary predictors
# ============================================================

binary_vars <- c(
  "female",
  "aceiarb_use",
  "bb_use",
  "mra_use"
)


cat("\n============================================\n")
cat("C. Binary predictor distributions\n")
cat("============================================\n")


binary_table <- do.call(

  rbind,

  lapply(
    binary_vars,
    function(v) {

      data.frame(

        variable =
          v,

        n =
          length(
            d1[[v]]
          ),

        n_yes =
          sum(
            d1[[v]] == 1,
            na.rm = TRUE
          ),

        pct_yes =
          round(
            mean(
              d1[[v]] == 1,
              na.rm = TRUE
            ) *
              100,
            1
          )
      )
    }
  )
)


print(
  binary_table,
  row.names = FALSE
)


# ============================================================
# 6. 30个插补数据集分别计算Spearman相关
#
# 不把30份数据直接堆起来计算相关性
# ============================================================

cor_vars <- c(

  "ridageyr",
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


cor_list <- vector(
  "list",
  imp50$m
)


for (i in 1:imp50$m) {

  d <- completed_list[[i]]

  cor_list[[i]] <-
    cor(
      d[
        cor_vars
      ],
      method = "spearman",
      use = "pairwise.complete.obs"
    )
}


# ============================================================
# 7. 每一对变量总结30次相关系数
# ============================================================

cor_pairs <- expand.grid(

  var1 =
    cor_vars,

  var2 =
    cor_vars,

  stringsAsFactors = FALSE
)


# 只保留上三角
cor_pairs <- cor_pairs[
  match(
    cor_pairs$var1,
    cor_vars
  ) <
    match(
      cor_pairs$var2,
      cor_vars
    ),
]


cor_summary <- do.call(

  rbind,

  lapply(
    seq_len(
      nrow(
        cor_pairs
      )
    ),
    function(j) {

      v1 <-
        cor_pairs$var1[j]

      v2 <-
        cor_pairs$var2[j]


      r_values <-
        sapply(
          cor_list,
          function(cm) {

            cm[
              v1,
              v2
            ]
          }
        )


      data.frame(

        var1 =
          v1,

        var2 =
          v2,

        median_r =
          median(
            r_values,
            na.rm = TRUE
          ),

        min_r =
          min(
            r_values,
            na.rm = TRUE
          ),

        max_r =
          max(
            r_values,
            na.rm = TRUE
          ),

        median_abs_r =
          median(
            abs(
              r_values
            ),
            na.rm = TRUE
          )
      )
    }
  )
)


cor_summary <-
  cor_summary %>%
  arrange(
    desc(
      median_abs_r
    )
  )


cat("\n============================================\n")
cat("D. Strongest Spearman correlations\n")
cat("============================================\n")

print(
  head(
    cor_summary,
    25
  ),
  row.names = FALSE
)


# ============================================================
# 8. 标记 |r| >= 0.70
# ============================================================

cat("\n============================================\n")
cat("E. Potentially redundant pairs |r| >= 0.70\n")
cat("============================================\n")


high_cor <- cor_summary %>%
  filter(
    median_abs_r >= 0.70
  )


if (
  nrow(
    high_cor
  ) == 0
) {

  cat(
    "No predictor pairs with median |r| >= 0.70\n"
  )

} else {

  print(
    high_cor,
    row.names = FALSE
  )
}


# ============================================================
# 9. VIF
#
# 这里只用于冗余审计
# 不代表最终模型
#
# 对30个插补数据分别计算
# ============================================================

# 手工VIF函数
# 避免依赖car包

calculate_vif <- function(data, vars) {

  result <- data.frame(

    variable =
      vars,

    vif =
      NA_real_
  )


  for (
    j in seq_along(
      vars
    )
  ) {

    outcome <-
      vars[j]

    predictors <-
      vars[
        -j
      ]


    formula_vif <-
      as.formula(
        paste(
          outcome,
          "~",
          paste(
            predictors,
            collapse = "+"
          )
        )
      )


    fit <-
      lm(
        formula_vif,
        data = data
      )


    r2 <-
      summary(
        fit
      )$r.squared


    result$vif[j] <-
      1 /
      (1 - r2)
  }


  return(
    result
  )
}


vif_vars <- candidate_vars


vif_all <- vector(
  "list",
  imp50$m
)


for (
  i in 1:imp50$m
) {

  vif_i <-
    calculate_vif(
      completed_list[[i]],
      vif_vars
    )

  vif_i$imputation <-
    i

  vif_all[[i]] <-
    vif_i
}


vif_all <-
  bind_rows(
    vif_all
  )


vif_summary <- vif_all %>%
  group_by(
    variable
  ) %>%
  summarise(

    median_vif =
      median(
        vif,
        na.rm = TRUE
      ),

    min_vif =
      min(
        vif,
        na.rm = TRUE
      ),

    max_vif =
      max(
        vif,
        na.rm = TRUE
      ),

    .groups = "drop"

  ) %>%
  arrange(
    desc(
      median_vif
    )
  )


cat("\n============================================\n")
cat("F. VIF summary across 30 imputations\n")
cat("============================================\n")


print(
  vif_summary,
  row.names = FALSE
)


# ============================================================
# 10. 标记明显共线性
# ============================================================

cat("\n============================================\n")
cat("G. Variables with median VIF >=5\n")
cat("============================================\n")


high_vif <-
  vif_summary %>%
  filter(
    median_vif >= 5
  )


if (
  nrow(
    high_vif
  ) == 0
) {

  cat(
    "None.\n"
  )

} else {

  print(
    high_vif,
    row.names = FALSE
  )
}


# ============================================================
# 11. 事件数 / 候选参数数量
# ============================================================

cat("\n============================================\n")
cat("H. Sample size / events / candidate predictors\n")
cat("============================================\n")


cat(
  "HF patients:",
  nrow(d1),
  "\n"
)

cat(
  "Deaths:",
  sum(
    d1$death == 1
  ),
  "\n"
)

cat(
  "Candidate predictors:",
  length(
    candidate_vars
  ),
  "\n"
)

cat(
  "Events per candidate predictor:",
  round(
    sum(
      d1$death == 1
    ) /
      length(
        candidate_vars
      ),
    1
  ),
  "\n"
)


# ============================================================
# 12. 保存
# ============================================================

write.csv(
  distribution_table,
  "HF_candidate_distribution.csv",
  row.names = FALSE
)

write.csv(
  binary_table,
  "HF_candidate_binary_distribution.csv",
  row.names = FALSE
)

write.csv(
  cor_summary,
  "HF_candidate_correlations_MI.csv",
  row.names = FALSE
)

write.csv(
  vif_summary,
  "HF_candidate_VIF_MI.csv",
  row.names = FALSE
)


save(
  completed_list,
  candidate_vars,
  variable_dictionary,
  distribution_table,
  binary_table,
  cor_summary,
  vif_summary,
  file =
    "NHANES_2011_2018_STEP5A_PREDICTOR_AUDIT.RData"
)


cat("\n============================================\n")
cat("STEP 5A 完成\n")
cat("============================================\n")
