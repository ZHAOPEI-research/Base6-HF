# ============================================================
# STEP 3A
# 构造基础 HF analysis cohort
# ============================================================

rm(list = ls())

library(dplyr)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load("NHANES_2011_2018_STEP2_CLEAN.RData")

cat("============================================\n")
cat("STEP 3A: HF基础队列构造\n")
cat("============================================\n")


# ============================================================
# 1. CKD-EPI 2021 eGFR
# ============================================================

dat <- dat %>%
  mutate(

    kappa = case_when(
      female == 1 ~ 0.7,
      female == 0 ~ 0.9,
      TRUE ~ NA_real_
    ),

    alpha = case_when(
      female == 1 ~ -0.241,
      female == 0 ~ -0.302,
      TRUE ~ NA_real_
    ),

    egfr_2021 = case_when(

      !is.na(scr_mgdl) &
        scr_mgdl > 0 &
        !is.na(ridageyr) &
        !is.na(female) ~

        142 *
        pmin(scr_mgdl / kappa, 1)^alpha *
        pmax(scr_mgdl / kappa, 1)^(-1.200) *
        (0.9938^ridageyr) *
        ifelse(female == 1, 1.012, 1),

      TRUE ~ NA_real_
    )
  )


# ============================================================
# 2. 基础纳排流程
# ============================================================

flow <- data.frame(
  step = character(),
  n = integer()
)

add_flow <- function(label, x) {
  flow <<- rbind(
    flow,
    data.frame(
      step = label,
      n = nrow(x)
    )
  )
}


# 全体死亡链接 eligible 数据
cohort0 <- dat

add_flow(
  "Eligible mortality-linked participants",
  cohort0
)


# 年龄 >=20
cohort1 <- cohort0 %>%
  filter(
    ridageyr >= 20
  )

add_flow(
  "Age >=20 years",
  cohort1
)


# 自报 CHF
cohort2 <- cohort1 %>%
  filter(
    hf == 1
  )

add_flow(
  "Self-reported physician-diagnosed HF",
  cohort2
)


# MEC examination follow-up time available
cohort3 <- cohort2 %>%
  filter(
    !is.na(permth_exm),
    permth_exm > 0
  )

add_flow(
  "Valid MEC-based follow-up",
  cohort3
)


# Positive MEC weight
hf_base <- cohort3 %>%
  filter(
    !is.na(wtmec2yr),
    wtmec2yr > 0
  )

add_flow(
  "Positive MEC examination weight",
  hf_base
)


cat("\n============================================\n")
cat("A. Cohort flow\n")
cat("============================================\n")

print(flow)


cat("\n============================================\n")
cat("B. HF base cohort by cycle\n")
cat("============================================\n")

print(
  hf_base %>%
    count(cycle)
)


cat("\n============================================\n")
cat("C. Death status\n")
cat("============================================\n")

print(
  table(
    hf_base$mortstat,
    useNA = "ifany"
  )
)


# ============================================================
# 3. eGFR 检查
# ============================================================

cat("\n============================================\n")
cat("D. eGFR 2021\n")
cat("============================================\n")

print(
  summary(
    hf_base$egfr_2021
  )
)

cat(
  "\neGFR missing:",
  sum(
    is.na(hf_base$egfr_2021)
  ),
  "\n"
)

cat(
  "eGFR <15:",
  sum(
    hf_base$egfr_2021 < 15,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "eGFR 15-29:",
  sum(
    hf_base$egfr_2021 >= 15 &
      hf_base$egfr_2021 < 30,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "eGFR 30-59:",
  sum(
    hf_base$egfr_2021 >= 30 &
      hf_base$egfr_2021 < 60,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "eGFR >=60:",
  sum(
    hf_base$egfr_2021 >= 60,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 4. HF cohort 内关键变量缺失
# ============================================================

key_vars <- c(
  "ridageyr",
  "female",
  "sbp_mean",
  "dbp_mean",
  "bpxpls",
  "bmxbmi",
  "lbxhgb",
  "lbxrdw",
  "lbxwbcsi",
  "lbdneno",
  "nwr",
  "scr_mgdl",
  "egfr_2021",
  "lbxsal",
  "lbxsbu",
  "lbxsnasi",
  "walk_tol",
  "med_cnt",
  "aceiarb_use",
  "bb_use",
  "mra_use",
  "mortstat",
  "permth_exm",
  "wtmec2yr"
)

key_vars <- intersect(
  key_vars,
  names(hf_base)
)

missing_hf <- data.frame(
  variable = key_vars,
  n_missing = sapply(
    hf_base[key_vars],
    function(x) sum(is.na(x))
  ),
  pct_missing = round(
    sapply(
      hf_base[key_vars],
      function(x) mean(is.na(x)) * 100
    ),
    1
  )
) %>%
  arrange(
    desc(pct_missing)
  )

cat("\n============================================\n")
cat("E. HF cohort key-variable missingness\n")
cat("============================================\n")

print(
  missing_hf,
  row.names = FALSE
)


# ============================================================
# 5. 极端值个案
# ============================================================

extreme_cases <- hf_base %>%
  filter(
    flag_wbc_extreme == 1 |
      flag_lymph_extreme == 1 |
      flag_creatinine_extreme == 1 |
      flag_sodium_extreme == 1
  ) %>%
  select(
    seqn,
    cycle,
    ridageyr,
    sex,
    mortstat,
    lbxwbcsi,
    lbdneno,
    lbdlymno,
    nwr,
    lbxscr,
    egfr_2021,
    lbxsnasi,
    lbxhgb,
    lbxrdw
  )

cat("\n============================================\n")
cat("F. HF cohort extreme-value cases\n")
cat("============================================\n")

print(extreme_cases)


# ============================================================
# 6. Walk_tol 在正式 HF 基础队列中的分布
# ============================================================

cat("\n============================================\n")
cat("G. Walk_tol in HF cohort\n")
cat("============================================\n")

print(
  table(
    hf_base$walk_tol,
    useNA = "ifany"
  )
)


# ============================================================
# 7. 保存
# ============================================================

save(
  dat,
  hf_base,
  flow,
  missing_hf,
  extreme_cases,
  file =
    "NHANES_2011_2018_STEP3A_HF_BASE.RData"
)

cat("\n============================================\n")
cat("STEP 3A 完成\n")
cat("============================================\n")
