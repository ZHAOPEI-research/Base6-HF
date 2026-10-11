# ============================================================
# STEP 9A-1
# Chinese external cohort: NWR data-quality audit
#
# 目的：
# 1. 找出 neutrophil 无法转 numeric 的原始值
# 2. 找出 NWR > 1
# 3. 找出 ANC/WBC 与 neu% 明显不一致的患者
# 4. 判断是否可以用 neu%/100 作为最终 NWR
#
# 暂时不修改任何患者数据
# ============================================================

rm(list = ls())

library(dplyr)

# Public-repository path configuration
# Restricted external-cohort data and derived individual-level outputs
# must remain outside the public repository.
.base6_external_output_dir <- Sys.getenv("BASE6_EXTERNAL_OUTPUT_DIR")

if (.base6_external_output_dir == "") {
  stop(
    paste0(
      "Set BASE6_EXTERNAL_OUTPUT_DIR to an authorized local directory ",
      "for restricted external-cohort analysis outputs."
    )
  )
}

if (!dir.exists(.base6_external_output_dir)) {
  dir.create(
    .base6_external_output_dir,
    recursive = TRUE
  )
}

setwd(.base6_external_output_dir)

load(
  "Chinese_HF_EXTERNAL_STEP9A_CLEAN.RData"
)

cat("============================================\n")
cat("STEP 9A-1: Chinese NWR audit\n")
cat("============================================\n")


# ============================================================
# 1. 原始变量类型
# ============================================================

cat("\n============================================\n")
cat("A. Raw variable classes\n")
cat("============================================\n")

cat(
  "WBC class:",
  class(chinese_data_raw$WBC),
  "\n"
)

cat(
  "neutrophil class:",
  class(chinese_data_raw$neutrophil),
  "\n"
)

cat(
  "neu% class:",
  class(chinese_data_raw$`neu%`),
  "\n"
)


# ============================================================
# 2. 找出 neutrophil 不能直接转换为数字的值
# ============================================================

neu_raw_chr <-
  trimws(
    as.character(
      chinese_data_raw$neutrophil
    )
  )


neu_numeric <-
  suppressWarnings(
    as.numeric(
      neu_raw_chr
    )
  )


bad_conversion <-
  chinese_data_raw %>%
  mutate(

    row_id =
      dplyr::row_number(),

    neutrophil_raw =
      neu_raw_chr,

    neutrophil_numeric =
      neu_numeric

  ) %>%
  filter(

    !is.na(neutrophil_raw) &
      neutrophil_raw != "" &
      is.na(neutrophil_numeric)

  ) %>%
  select(

    row_id,
    ID,

    WBC,
    `neu%`,
    neutrophil,

    neutrophil_raw
  )


cat("\n============================================\n")
cat("B. Neutrophil values failing numeric conversion\n")
cat("============================================\n")

cat(
  "Number of problematic rows:",
  nrow(bad_conversion),
  "\n"
)

print(
  as.data.frame(
    bad_conversion
  ),
  row.names = FALSE
)


# ============================================================
# 3. 所有 NWR > 1 的病例
# ============================================================

nwr_over1 <-
  chinese_validation %>%
  mutate(
    row_id =
      dplyr::row_number()
  ) %>%
  filter(
    !is.na(nwr) &
      nwr > 1
  ) %>%
  select(

    row_id,
    ID,
    status,
    time_day,

    WBC,
    neutrophil,
    `neu%`,

    WBC_clean,
    neutrophil_clean,

    nwr,
    nwr_from_neu_pct,
    nwr_diff
  )


cat("\n============================================\n")
cat("C. NWR > 1 cases\n")
cat("============================================\n")

print(
  as.data.frame(
    nwr_over1
  ),
  row.names = FALSE
)


# ============================================================
# 4. ANC/WBC 与 neu%/100 差异 >0.05
# ============================================================

nwr_discordant <-
  chinese_validation %>%
  mutate(
    row_id =
      dplyr::row_number()
  ) %>%
  filter(

    !is.na(nwr_diff) &
      abs(nwr_diff) > 0.05

  ) %>%
  select(

    row_id,
    ID,

    WBC,
    neutrophil,
    `neu%`,

    nwr,
    nwr_from_neu_pct,
    nwr_diff
  ) %>%

  arrange(
    desc(
      abs(nwr_diff)
    )
  )


cat("\n============================================\n")
cat("D. ANC/WBC vs neu% discordance >0.05\n")
cat("============================================\n")

cat(
  "Number:",
  nrow(nwr_discordant),
  "\n"
)

print(
  as.data.frame(
    nwr_discordant
  ),
  row.names = FALSE
)


# ============================================================
# 5. 检查 neu% 本身异常值
#
# 百分比原则上应在 0-100
# 极低值也单独展示
# ============================================================

neu_pct_suspicious <-
  chinese_validation %>%
  mutate(
    row_id =
      dplyr::row_number()
  ) %>%
  filter(

    !is.na(neu_pct_clean) &

      (
        neu_pct_clean < 20 |
          neu_pct_clean > 95
      )

  ) %>%
  select(

    row_id,
    ID,

    WBC,
    neutrophil,
    `neu%`,

    nwr,
    nwr_from_neu_pct,
    nwr_diff
  )


cat("\n============================================\n")
cat("E. Suspicious neutrophil percentages\n")
cat("============================================\n")

print(
  as.data.frame(
    neu_pct_suspicious
  ),
  row.names = FALSE
)


# ============================================================
# 6. 只在两种方法都合理的病例中评价一致性
#
# NWR理论上：
#
# absolute neutrophil / WBC
# ~= neutrophil percentage / 100
# ============================================================

valid_compare <-
  chinese_validation %>%
  filter(

    !is.na(nwr),

    !is.na(nwr_from_neu_pct),

    nwr >= 0,
    nwr <= 1,

    nwr_from_neu_pct >= 0,
    nwr_from_neu_pct <= 1
  )


cat("\n============================================\n")
cat("F. Agreement in plausible observations\n")
cat("============================================\n")


cat(
  "n =",
  nrow(valid_compare),
  "\n"
)


cat(
  "Pearson correlation =",
  cor(
    valid_compare$nwr,
    valid_compare$nwr_from_neu_pct,
    method = "pearson"
  ),
  "\n"
)


cat(
  "Spearman correlation =",
  cor(
    valid_compare$nwr,
    valid_compare$nwr_from_neu_pct,
    method = "spearman"
  ),
  "\n"
)


cat("\nAbsolute difference:\n")

print(
  summary(
    abs(
      valid_compare$nwr -
        valid_compare$nwr_from_neu_pct
    )
  )
)


cat(
  "\n|difference| <= 0.01:",
  sum(
    abs(
      valid_compare$nwr -
        valid_compare$nwr_from_neu_pct
    ) <= 0.01
  ),
  "/",
  nrow(valid_compare),
  "\n"
)


cat(
  "|difference| <= 0.02:",
  sum(
    abs(
      valid_compare$nwr -
        valid_compare$nwr_from_neu_pct
    ) <= 0.02
  ),
  "/",
  nrow(valid_compare),
  "\n"
)


# ============================================================
# 7. 对照开发队列 NWR 范围
#
# Frozen spline:
# boundary:
# 0.4347826, 0.7736074
#
# 原NHANES开发样本大约：
# 0.088 - 0.968
# ============================================================

cat("\n============================================\n")
cat("G. External NWR distribution comparison\n")
cat("============================================\n")


cat("\nANC/WBC NWR:\n")

print(
  quantile(
    chinese_validation$nwr,
    probs = c(
      0,
      0.01,
      0.025,
      0.05,
      0.25,
      0.5,
      0.75,
      0.95,
      0.975,
      0.99,
      1
    ),
    na.rm = TRUE
  )
)


cat("\nneu% / 100:\n")

print(
  quantile(
    chinese_validation$nwr_from_neu_pct,
    probs = c(
      0,
      0.01,
      0.025,
      0.05,
      0.25,
      0.5,
      0.75,
      0.95,
      0.975,
      0.99,
      1
    ),
    na.rm = TRUE
  )
)


# ============================================================
# 8. 列出最极端的10个 ANC/WBC NWR
# ============================================================

extreme_nwr <-
  chinese_validation %>%
  mutate(
    row_id =
      dplyr::row_number()
  ) %>%
  arrange(
    desc(nwr)
  ) %>%
  select(

    row_id,
    ID,

    WBC,
    neutrophil,
    `neu%`,

    nwr,
    nwr_from_neu_pct,
    nwr_diff

  ) %>%
  head(10)


cat("\n============================================\n")
cat("H. Ten highest ANC/WBC NWR values\n")
cat("============================================\n")

print(
  as.data.frame(
    extreme_nwr
  ),
  row.names = FALSE
)


# ============================================================
# 9. 保存审计结果
# ============================================================

write.csv(
  bad_conversion,
  "Chinese_NWR_bad_conversion.csv",
  row.names = FALSE
)

write.csv(
  nwr_over1,
  "Chinese_NWR_over1.csv",
  row.names = FALSE
)

write.csv(
  nwr_discordant,
  "Chinese_NWR_discordant.csv",
  row.names = FALSE
)

write.csv(
  neu_pct_suspicious,
  "Chinese_NWR_suspicious_neu_pct.csv",
  row.names = FALSE
)


save(

  bad_conversion,
  nwr_over1,
  nwr_discordant,
  neu_pct_suspicious,
  valid_compare,
  extreme_nwr,

  file =
    "Chinese_HF_EXTERNAL_NWR_AUDIT.RData"
)


cat("\n============================================\n")
cat("STEP 9A-1 完成\n")
cat("============================================\n")
