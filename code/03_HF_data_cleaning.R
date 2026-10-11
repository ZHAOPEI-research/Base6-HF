# ============================================================
# STEP 2
# NHANES 2011-2018
# 类型统一 + 问卷无效值处理 + 基础衍生变量
# ============================================================

rm(list = ls())

library(dplyr)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load("NHANES_2011_2018_RAW_HFstudy.RData")

dat <- final_all

cat("============================================\n")
cat("STEP 2: 数据编码清洗\n")
cat("============================================\n")

cat("原始样本量:", nrow(dat), "\n")


# ============================================================
# 1. 工具函数
# ============================================================

# ------------------------------------------------------------
# 连续变量：character -> numeric
# ------------------------------------------------------------

to_num <- function(x) {

  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}


# ------------------------------------------------------------
# Yes / No 型变量
#
# Yes = 1
# No  = 0
# Refused / Don't know / Missing = NA
# ------------------------------------------------------------

yn_clean <- function(x) {

  x <- trimws(
    as.character(x)
  )

  case_when(

    x == "Yes" ~ 1L,

    x == "No" ~ 0L,

    x %in% c(
      "Refused",
      "Don't know",
      "Don't Know",
      "Unknown",
      ""
    ) ~ NA_integer_,

    is.na(x) ~ NA_integer_,

    TRUE ~ NA_integer_
  )
}


# ============================================================
# 2. 连续变量统一 numeric
# ============================================================

numeric_vars <- c(

  # survey
  "sddsrvyr",
  "wtint2yr",
  "wtmec2yr",
  "sdmvstra",
  "sdmvpsu",

  # demographic
  "ridageyr",
  "indfmpir",

  # BP
  "bpxsy1",
  "bpxsy2",
  "bpxsy3",
  "bpxdi1",
  "bpxdi2",
  "bpxdi3",
  "bpxpls",

  # BMI
  "bmxbmi",

  # CBC
  "lbxhgb",
  "lbxrdw",
  "lbxwbcsi",
  "lbdneno",
  "lbdlymno",

  # glucose
  "lbxglu",
  "lbdglusi",
  "lbxgh",

  # biochemistry
  "lbxscr",
  "lbdscrsi",
  "lbxsal",
  "lbdsalsi",
  "lbxsbu",
  "lbdsbusi",
  "lbxsnasi",

  # kidney
  "urdact",

  # inflammation
  "lbxhscrp",

  # lipid
  "lbdhdd",
  "lbdhddsi",
  "lbxtr",
  "lbdtrsi",
  "lbdldl",
  "lbdldnsi",
  "lbxtc",
  "lbdtcsi",

  # activity
  "pad680",

  # medication
  "med_cnt",

  # mortality
  "eligstat",
  "mortstat",
  "ucod_leading",
  "permth_int",
  "permth_exm"
)

numeric_vars <- intersect(
  numeric_vars,
  names(dat)
)

dat <- dat %>%
  mutate(
    across(
      all_of(numeric_vars),
      to_num
    )
  )


# ============================================================
# 3. Sex
# ============================================================

dat <- dat %>%
  mutate(

    sex = case_when(

      riagendr == "Male" ~ "Male",

      riagendr == "Female" ~ "Female",

      TRUE ~ NA_character_
    ),

    female = case_when(

      riagendr == "Female" ~ 1L,

      riagendr == "Male" ~ 0L,

      TRUE ~ NA_integer_
    )
  )


# ============================================================
# 4. HF / cardiovascular disease
# ============================================================

dat <- dat %>%
  mutate(

    hf = yn_clean(mcq160b),

    chd = yn_clean(mcq160c),

    angina = yn_clean(mcq160d),

    mi = yn_clean(mcq160e),

    stroke = yn_clean(mcq160f)
  )


# ============================================================
# 5. Hypertension
# ============================================================

dat <- dat %>%
  mutate(

    hypertension_dx =
      yn_clean(
        bpq020
      ),

    hypertension_med =
      yn_clean(
        bpq050a
      )
  )


# ============================================================
# 6. Diabetes
# ============================================================

# DIQ010:
# Yes
# No
# Borderline
# Don't know

dat <- dat %>%
  mutate(

    diabetes = case_when(

      diq010 == "Yes" ~ 1L,

      diq010 == "No" ~ 0L,

      # Borderline 不等于确诊 diabetes
      diq010 == "Borderline" ~ 0L,

      diq010 %in% c(
        "Don't know",
        "Don't Know",
        "Refused"
      ) ~ NA_integer_,

      TRUE ~ NA_integer_
    ),

    diabetes_borderline =
      case_when(

        diq010 == "Borderline" ~ 1L,

        diq010 %in% c(
          "Yes",
          "No"
        ) ~ 0L,

        TRUE ~ NA_integer_
      ),

    insulin_use =
      yn_clean(
        diq050
      ),

    diabetes_med =
      yn_clean(
        diq070
      )
  )


# ============================================================
# 7. Kidney disease
# ============================================================

dat <- dat %>%
  mutate(

    kidney_disease =
      yn_clean(
        kiq022
      ),

    dialysis =
      yn_clean(
        kiq025
      )
  )


# ============================================================
# 8. Smoking
# ============================================================

dat <- dat %>%
  mutate(

    ever_smoker =
      yn_clean(
        smq020
      ),

    smoking_status =
      case_when(

        smq020 == "No" ~
          "Never",

        smq020 == "Yes" &
          smq040 == "Not at all" ~
          "Former",

        smq020 == "Yes" &
          smq040 %in%
          c(
            "Every day",
            "Some days"
          ) ~
          "Current",

        TRUE ~
          NA_character_
      )
  )


# ============================================================
# 9. PAQ
# ============================================================

dat <- dat %>%
  mutate(

    vigorous_work =
      yn_clean(
        paq605
      ),

    moderate_work =
      yn_clean(
        paq620
      ),

    vigorous_recreation =
      yn_clean(
        paq650
      ),

    moderate_recreation =
      yn_clean(
        paq665
      )
  )



# ============================================================
# 10. PAD680
# Minutes sedentary activity
#
# NHANES有效范围：0–1380 min/day
# 7777 = Refused
# 9999 = Don't know
# ============================================================

dat <- dat %>%
  mutate(

    # 保留原始值，方便审计
    pad680_raw = pad680,

    # 特殊编码标记
    pad680_invalid_reason = case_when(

      pad680 == 7777 ~ "Refused",

      pad680 == 9999 ~ "Don't know",

      !is.na(pad680) &
        (pad680 < 0 | pad680 > 1380) ~
        "Out of valid range",

      TRUE ~ NA_character_
    ),

    # 正式用于分析的变量
    sedentary_min = case_when(

      is.na(pad680) ~
        NA_real_,

      pad680 %in% c(
        7777,
        9999
      ) ~
        NA_real_,

      pad680 >= 0 &
        pad680 <= 1380 ~
        pad680,

      TRUE ~
        NA_real_
    )
  )




cat("\n============================================\n")
cat("PAD680 清洗检查\n")
cat("============================================\n")

cat("\n原始 PAD680 特殊值：\n")

print(
  table(
    dat$pad680[
      dat$pad680 > 1380
    ],
    useNA = "ifany"
  )
)

cat("\n无效原因：\n")

print(
  table(
    dat$pad680_invalid_reason,
    useNA = "ifany"
  )
)

cat("\n清洗后 sedentary_min：\n")

print(
  summary(
    dat$sedentary_min
  )
)

# ============================================================
# 11. PFQ061B
# walking limitation
#
# 不直接命名 NYHA
# ============================================================

dat <- dat %>%
  mutate(

    walk_tol = case_when(

      pfq061b ==
        "No difficulty" ~ 1,

      pfq061b ==
        "Some difficulty" ~ 2,

      pfq061b ==
        "Much difficulty" ~ 3,

      pfq061b ==
        "Unable to do" ~ 4,

      # 无法判断是身体限制还是主观不做
      pfq061b ==
        "Do not do this activity" ~ NA_real_,

      pfq061b %in%
        c(
          "Don't know",
          "Don't Know",
          "Refused"
        ) ~ NA_real_,

      TRUE ~ NA_real_
    ),

    walk_limitation =
      case_when(

        walk_tol == 1 ~ 0L,

        walk_tol %in%
          c(
            2,
            3,
            4
          ) ~ 1L,

        TRUE ~ NA_integer_
      )
  )


# ============================================================
# 12. NWR
#
# neutrophil absolute count / WBC
# ============================================================

dat <- dat %>%
  mutate(

    nwr = case_when(

      !is.na(lbdneno) &
        !is.na(lbxwbcsi) &
        lbxwbcsi > 0 ~
        lbdneno / lbxwbcsi,

      TRUE ~
        NA_real_
    )
  )


# ============================================================
# 13. 平均 SBP / DBP
#
# 使用实际存在的血压读数平均
# 不要求必须三次完整
# ============================================================

dat <- dat %>%
  rowwise() %>%
  mutate(

    sbp_mean = {

      x <- c(
        bpxsy1,
        bpxsy2,
        bpxsy3
      )

      if (
        all(
          is.na(x)
        )
      ) {

        NA_real_

      } else {

        mean(
          x,
          na.rm = TRUE
        )
      }
    },

    dbp_mean = {

      x <- c(
        bpxdi1,
        bpxdi2,
        bpxdi3
      )

      if (
        all(
          is.na(x)
        )
      ) {

        NA_real_

      } else {

        mean(
          x,
          na.rm = TRUE
        )
      }
    }
  ) %>%
  ungroup()


# ============================================================
# 14. Creatinine consistency
#
# 不删除极值
# 仅确认两种单位关系
# ============================================================

dat <- dat %>%
  mutate(

    scr_mgdl =
      lbxscr,

    scr_umol =
      lbdscrsi,

    scr_conversion_ratio =
      case_when(

        !is.na(scr_mgdl) &
          scr_mgdl > 0 &
          !is.na(scr_umol) ~
          scr_umol /
          scr_mgdl,

        TRUE ~
          NA_real_
      )
  )


# ============================================================
# 15. 暂时标记可疑极端值
#
# 这里只 FLAG，不删除
# ============================================================

dat <- dat %>%
  mutate(

    flag_wbc_extreme =
      case_when(

        is.na(lbxwbcsi) ~ NA_integer_,

        lbxwbcsi > 100 ~ 1L,

        TRUE ~ 0L
      ),

    flag_lymph_extreme =
      case_when(

        is.na(lbdlymno) ~ NA_integer_,

        lbdlymno > 100 ~ 1L,

        TRUE ~ 0L
      ),

    flag_creatinine_extreme =
      case_when(

        is.na(lbxscr) ~ NA_integer_,

        lbxscr > 10 ~ 1L,

        TRUE ~ 0L
      ),

    flag_sodium_extreme =
      case_when(

        is.na(lbxsnasi) ~ NA_integer_,

        lbxsnasi < 110 |
          lbxsnasi > 160 ~ 1L,

        TRUE ~ 0L
      )
  )


# ============================================================
# 16. 检查 HF
# ============================================================

cat("\n============================================\n")
cat("A. HF 清洗后分布\n")
cat("============================================\n")

print(
  table(
    dat$hf,
    useNA = "ifany"
  )
)

cat("\n按周期：\n")

print(
  dat %>%
    count(
      cycle,
      hf
    )
)


# ============================================================
# 17. ≥20岁 HF 人数
# ============================================================

cat("\n============================================\n")
cat("B. ≥20岁 HF 候选队列\n")
cat("============================================\n")

hf_age_check <- dat %>%
  filter(
    ridageyr >= 20
  ) %>%
  count(
    cycle,
    hf
  )

print(hf_age_check)


# ============================================================
# 18. Walk_tol
# ============================================================

cat("\n============================================\n")
cat("C. Walk_tol 分布\n")
cat("============================================\n")

print(
  table(
    dat$walk_tol,
    useNA = "ifany"
  )
)

cat("\nHF患者中的 Walk_tol：\n")

print(
  dat %>%
    filter(
      hf == 1,
      ridageyr >= 20
    ) %>%
    count(
      walk_tol
    )
)


# ============================================================
# 19. NWR
# ============================================================

cat("\n============================================\n")
cat("D. NWR 检查\n")
cat("============================================\n")

print(
  summary(
    dat$nwr
  )
)

cat(
  "NWR > 1：",
  sum(
    dat$nwr > 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "NWR < 0：",
  sum(
    dat$nwr < 0,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 20. 极端值 FLAG
# ============================================================

cat("\n============================================\n")
cat("E. 可疑极端值数量\n")
cat("============================================\n")

cat(
  "WBC >100:",
  sum(
    dat$flag_wbc_extreme == 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Lymphocyte >100:",
  sum(
    dat$flag_lymph_extreme == 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Creatinine >10:",
  sum(
    dat$flag_creatinine_extreme == 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Na <110 or >160:",
  sum(
    dat$flag_sodium_extreme == 1,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 21. MEC weight
# ============================================================

cat("\n============================================\n")
cat("F. MEC Weight\n")
cat("============================================\n")

print(
  dat %>%
    group_by(
      cycle
    ) %>%
    summarise(

      n = n(),

      wtmec_zero =
        sum(
          wtmec2yr <= 0,
          na.rm = TRUE
        ),

      wtmec_positive =
        sum(
          wtmec2yr > 0,
          na.rm = TRUE
        ),

      .groups = "drop"
    )
)


# ============================================================
# 22. 初步 HF candidate
#
# 这里只建立候选对象
# 尚未做最终纳排
# ============================================================

hf_candidate <- dat %>%
  filter(

    ridageyr >= 20,

    hf == 1
  )


cat("\n============================================\n")
cat("G. HF candidate\n")
cat("============================================\n")

cat(
  "HF candidate n =",
  nrow(hf_candidate),
  "\n"
)

print(
  hf_candidate %>%
    count(
      cycle
    )
)


# ============================================================
# 23. 保存
# ============================================================

save(
  dat,
  hf_candidate,
  file =
    "NHANES_2011_2018_STEP2_CLEAN.RData"
)

cat("\n============================================\n")
cat("STEP 2 完成\n")
cat("============================================\n")

cat(
  "保存：NHANES_2011_2018_STEP2_CLEAN.RData\n"
)
