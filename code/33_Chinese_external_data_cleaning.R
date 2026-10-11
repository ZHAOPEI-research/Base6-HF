# ============================================================
# STEP 9A
# Chinese external validation data cleaning
#
# 中国住院HF队列
# 外部随访最长12个月
#
# 目标：
# 1. Excel读取
# 2. 结局/随访时间审计
# 3. NWR
# 4. CKD-EPI 2021 eGFR
# 5. medication count
# 6. Base6变量范围和缺失审计
#
# 暂不计算模型预测
# ============================================================

rm(list = ls())

library(readxl)
library(dplyr)




cat("============================================\n")
cat("STEP 9A: Chinese external data cleaning\n")
cat("============================================\n")


# ============================================================
# 1. Read Excel
# ============================================================

# External hospital cohort data are not distributed with this repository.
# Authorized users should provide the local path to the deidentified dataset.

external_data_path <- Sys.getenv("BASE6_EXTERNAL_DATA")

if (external_data_path == "") {
  stop(
    paste0(
      "External-cohort data are not included in this public repository. ",
      "Set the environment variable BASE6_EXTERNAL_DATA to the local path ",
      "of the authorized deidentified hospital dataset."
    )
  )
}

chinese_data <- read_excel(
  external_data_path
)


cat(
  "\nOriginal sample size =",
  nrow(chinese_data),
  "\n"
)


cat(
  "Number of variables =",
  ncol(chinese_data),
  "\n"
)


# 保留原始数据
chinese_data_raw <- chinese_data

# 分析副本
chinese_validation <- chinese_data


# ============================================================
# 2. Required-variable check
# ============================================================

required_vars <- c(

  "ID",
  "status",
  "time",

  "sex",
  "age",

  "BMI",
  "WBC",
  "neutrophil",
  "RDW",

  "CR",

  "ACEI/ARB/ANRi",
  "β-blocker",
  "MRA",
  "dapagliflozin",
  "diuretic",
  "digoxin",
  "anti_platelet",
  "statin",
  "anticoagulants",
  "nitrates",
  "CCB",
  "PPI",
  "renal_drug",
  "NASIDs",
  "febuxostat",
  "agents_electrolyte",
  "other_drug_cnt"
)


missing_vars <- setdiff(
  required_vars,
  names(chinese_validation)
)


cat("\n============================================\n")
cat("A. Required-variable check\n")
cat("============================================\n")


if (length(missing_vars) == 0) {

  cat(
    "All required variables are present.\n"
  )

} else {

  cat(
    "Missing variables:\n"
  )

  print(
    missing_vars
  )

  stop(
    "Required variables are missing."
  )
}


# ============================================================
# 3. Outcome audit
#
# time单位 = days
# 中国队列最长随访12个月
# ============================================================

cat("\n============================================\n")
cat("B. Outcome and follow-up\n")
cat("============================================\n")


cat("\nstatus:\n")

print(
  table(
    chinese_validation$status,
    useNA = "ifany"
  )
)


cat("\nFollow-up time (days):\n")

print(
  summary(
    chinese_validation$time
  )
)


cat(
  "\nMissing follow-up:",
  sum(
    is.na(
      chinese_validation$time
    )
  ),
  "\n"
)


cat(
  "Missing status:",
  sum(
    is.na(
      chinese_validation$status
    )
  ),
  "\n"
)


# ============================================================
# 4. Derived follow-up variables
#
# 原始time_day永久保留
#
# 后续1-year validation：
# target = 365 days
#
# time_month仅用于描述，
# 不作为主要外部validation时间尺度
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    time_day =
      as.numeric(
        time
      ),

    time_month =
      time_day /
      30.4375
  )


cat("\nDerived follow-up months:\n")

print(
  summary(
    chinese_validation$time_month
  )
)


cat(
  "\nPatients with follow-up >=365 days:",
  sum(
    chinese_validation$time_day >= 365,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 5. Sex coding audit
# ============================================================

cat("\n============================================\n")
cat("C. Sex coding\n")
cat("============================================\n")


print(
  table(
    chinese_validation$sex,
    useNA = "ifany"
  )
)


# ============================================================
# 6. Construct female
#
# 暂支持：
# 1 male / 2 female
# 男/女
# M/F
# Male/Female
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    female =
      case_when(

        sex %in% c(
          2,
          "2",
          "Female",
          "female",
          "F",
          "女"
        ) ~ 1,

        sex %in% c(
          1,
          "1",
          "Male",
          "male",
          "M",
          "男"
        ) ~ 0,

        TRUE ~ NA_real_
      )
  )


cat("\nDerived female:\n")

print(
  table(
    chinese_validation$female,
    useNA = "ifany"
  )
)


# ============================================================
# 7. Age transformation
#
# 开发数据NHANES年龄>=80全部编码80
#
# 外部应用必须保持一致：
# age_model = min(age, 80)
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    age_raw =
      as.numeric(
        age
      ),

    age_model =
      pmin(
        age_raw,
        80
      )
  )


cat("\n============================================\n")
cat("D. Age audit\n")
cat("============================================\n")


cat("\nOriginal age:\n")

print(
  summary(
    chinese_validation$age_raw
  )
)


cat(
  "\nAge >80:",
  sum(
    chinese_validation$age_raw > 80,
    na.rm = TRUE
  ),
  "\n"
)


cat("\nModel age after top-coding:\n")

print(
  summary(
    chinese_validation$age_model
  )
)


# ============================================================
# 8. Creatinine conversion
#
# Chinese CR = μmol/L
#
# 1 mg/dL = 88.4 μmol/L
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    scr_umol =
      as.numeric(
        CR
      ),

    scr_mgdl =
      scr_umol /
      88.4
  )


cat("\n============================================\n")
cat("E. Creatinine conversion\n")
cat("============================================\n")


cat("\nCR (umol/L):\n")

print(
  summary(
    chinese_validation$scr_umol
  )
)


cat("\nScr (mg/dL):\n")

print(
  summary(
    chinese_validation$scr_mgdl
  )
)


# ============================================================
# 9. CKD-EPI 2021 eGFR
#
# 为严格复现NHANES开发过程：
# 使用 age_model
#
# female:
# kappa = 0.7
# alpha = -0.241
#
# male:
# kappa = 0.9
# alpha = -0.302
# ============================================================

chinese_validation <- chinese_validation %>%
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
        scr_mgdl / kappa,
        1
      )^alpha *

      pmax(
        scr_mgdl / kappa,
        1
      )^(-1.200) *

      (0.9938^age_model) *

      ifelse(
        female == 1,
        1.012,
        1
      )
  )


cat("\n============================================\n")
cat("F. CKD-EPI 2021 eGFR\n")
cat("============================================\n")


print(
  summary(
    chinese_validation$egfr_2021
  )
)


cat(
  "\neGFR missing:",
  sum(
    is.na(
      chinese_validation$egfr_2021
    )
  ),
  "\n"
)


cat(
  "eGFR <15:",
  sum(
    chinese_validation$egfr_2021 < 15,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "eGFR 15-29:",
  sum(
    chinese_validation$egfr_2021 >= 15 &
      chinese_validation$egfr_2021 < 30,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "eGFR 30-59:",
  sum(
    chinese_validation$egfr_2021 >= 30 &
      chinese_validation$egfr_2021 < 60,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "eGFR >=60:",
  sum(
    chinese_validation$egfr_2021 >= 60,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 10. Construct NWR
#
# NWR =
# absolute neutrophil count / total WBC
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    WBC_clean =
      as.numeric(
        WBC
      ),

    neutrophil_clean =
      as.numeric(
        neutrophil
      ),

    nwr =
      ifelse(

        !is.na(WBC_clean) &
          !is.na(neutrophil_clean) &
          WBC_clean > 0,

        neutrophil_clean /
          WBC_clean,

        NA_real_
      )
  )


cat("\n============================================\n")
cat("G. NWR\n")
cat("============================================\n")


print(
  summary(
    chinese_validation$nwr
  )
)


cat(
  "\nNWR <=0:",
  sum(
    chinese_validation$nwr <= 0,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "NWR >1:",
  sum(
    chinese_validation$nwr > 1,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 11. Compare NWR with neu%
#
# 如果 neu% 为 60、70 这种百分数：
#
# neu% / 100 ≈ NWR
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    neu_pct_clean =
      as.numeric(
        `neu%`
      ),

    nwr_from_neu_pct =
      neu_pct_clean /
      100,

    nwr_diff =
      nwr -
      nwr_from_neu_pct
  )


cat("\n============================================\n")
cat("H. NWR vs neu%\n")
cat("============================================\n")


cat("\nneu%:\n")

print(
  summary(
    chinese_validation$neu_pct_clean
  )
)


cat("\nNWR - neu%/100:\n")

print(
  summary(
    chinese_validation$nwr_diff
  )
)


cat(
  "\n|difference| >0.05:",
  sum(
    abs(
      chinese_validation$nwr_diff
    ) > 0.05,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 12. Medication variables
#
# Important:
#
# ACEI/ARB/ANRi includes
# ACEI/ARB and ARNi
#
# Therefore:
#
# ACEI/ARB and ARNi ARE NOT counted again.
#
# Each listed drug variable represents:
# 0 = no
# 1 = one drug from this category
#
# other_drug_cnt:
# actual number of additional drugs
# ============================================================

med_binary_vars <- c(

  "ACEI/ARB/ANRi",

  "β-blocker",

  "MRA",

  "dapagliflozin",

  "diuretic",

  "digoxin",

  "anti_platelet",

  "statin",

  "anticoagulants",

  "nitrates",

  "CCB",

  "PPI",

  "renal_drug",

  "NASIDs",

  "febuxostat",

  "agents_electrolyte"
)


# ============================================================
# 13. Medication coding audit
# ============================================================

cat("\n============================================\n")
cat("I. Medication coding audit\n")
cat("============================================\n")


for (v in med_binary_vars) {

  cat(
    "\n--- ",
    v,
    " ---\n",
    sep = ""
  )

  print(
    table(
      chinese_validation[[v]],
      useNA = "ifany"
    )
  )
}


cat("\nother_drug_cnt:\n")

print(
  summary(
    chinese_validation$other_drug_cnt
  )
)


# ============================================================
# 14. Invalid medication values
# ============================================================

invalid_med <- do.call(

  rbind,

  lapply(

    med_binary_vars,

    function(v) {

      x <-
        chinese_validation[[v]]


      bad <-
        unique(
          x[
            !is.na(x) &
              !x %in% c(
                0,
                1
              )
          ]
        )


      data.frame(

        variable =
          v,

        invalid_values =
          if (
            length(bad) == 0
          ) {

            ""

          } else {

            paste(
              bad,
              collapse = ", "
            )
          },

        stringsAsFactors =
          FALSE
      )
    }
  )
)


cat("\n============================================\n")
cat("J. Non-binary medication values\n")
cat("============================================\n")


print(
  invalid_med,
  row.names = FALSE
)


# ============================================================
# 15. Medication missing count
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    med_missing_n =
      rowSums(

        is.na(

          across(

            all_of(
              med_binary_vars
            )
          )
        )
      ) +

      as.integer(
        is.na(
          other_drug_cnt
        )
      )
  )


cat("\n============================================\n")
cat("K. Medication missingness\n")
cat("============================================\n")


print(
  table(
    chinese_validation$med_missing_n,
    useNA = "ifany"
  )
)


# ============================================================
# 16. Construct med_cnt
#
# 非缺失情况下：
#
# med_cnt =
# sum(binary drug variables)
# +
# other_drug_cnt
#
# 如果任何药物信息缺失，
# 暂时把med_cnt设为NA
#
# 不直接把missing视为0。
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    med_binary_sum =
      rowSums(
        across(
          all_of(
            med_binary_vars
          )
        ),
        na.rm = FALSE
      ),

    med_cnt =
      ifelse(

        med_missing_n == 0,

        med_binary_sum +
          other_drug_cnt,

        NA_real_
      )
  )


cat("\n============================================\n")
cat("L. Derived medication count\n")
cat("============================================\n")


print(
  summary(
    chinese_validation$med_cnt
  )
)


cat(
  "\nmed_cnt missing:",
  sum(
    is.na(
      chinese_validation$med_cnt
    )
  ),
  "\n"
)


# ============================================================
# 17. Main Base6 variables
# ============================================================

base6_vars <- c(

  "age_model",
  "BMI",
  "RDW",
  "nwr",
  "egfr_2021",
  "med_cnt"
)


# ============================================================
# 18. Predictor distributions
# ============================================================

cat("\n============================================\n")
cat("M. Base6 predictor distributions\n")
cat("============================================\n")


base6_distribution <- do.call(

  rbind,

  lapply(

    base6_vars,

    function(v) {

      x <-
        as.numeric(
          chinese_validation[[v]]
        )


      if (
        sum(
          !is.na(x)
        ) == 0
      ) {

        return(

          data.frame(

            variable =
              v,

            n =
              0,

            missing =
              length(x),

            pct_missing =
              100,

            min =
              NA,

            p01 =
              NA,

            p25 =
              NA,

            median =
              NA,

            p75 =
              NA,

            p99 =
              NA,

            max =
              NA
          )
        )
      }


      data.frame(

        variable =
          v,

        n =
          sum(
            !is.na(x)
          ),

        missing =
          sum(
            is.na(x)
          ),

        pct_missing =
          round(
            mean(
              is.na(x)
            ) *
              100,
            2
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
  base6_distribution,
  row.names = FALSE
)


# ============================================================
# 19. Complete-case audit
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    base6_complete =
      complete.cases(
        across(
          all_of(
            base6_vars
          )
        )
      ),

    outcome_complete =
      !is.na(status) &
      !is.na(time_day)
  )


cat("\n============================================\n")
cat("N. Complete-case audit\n")
cat("============================================\n")


cat(
  "Total n =",
  nrow(
    chinese_validation
  ),
  "\n"
)


cat(
  "Outcome complete =",
  sum(
    chinese_validation$outcome_complete
  ),
  "\n"
)


cat(
  "Base6 complete =",
  sum(
    chinese_validation$base6_complete
  ),
  "\n"
)


cat(
  "Base6 complete proportion =",
  round(
    mean(
      chinese_validation$base6_complete
    ) *
      100,
    1
  ),
  "%\n"
)


cat(
  "Complete predictor + outcome =",
  sum(
    chinese_validation$base6_complete &
      chinese_validation$outcome_complete
  ),
  "\n"
)


# ============================================================
# 20. 1-year support
# ============================================================

cat("\n============================================\n")
cat("O. 1-year validation support\n")
cat("============================================\n")


cat(
  "Follow-up >=365 days:",
  sum(
    chinese_validation$time_day >= 365,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Deaths within 365 days:",
  sum(
    chinese_validation$status == 1 &
      chinese_validation$time_day <= 365,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Censored before 365 days:",
  sum(
    chinese_validation$status == 0 &
      chinese_validation$time_day < 365,
    na.rm = TRUE
  ),
  "\n"
)

# ============================================================
# 21. Save
# ============================================================

# Outputs derived from the restricted external cohort must be
# written outside the public repository.
#
# Authorized users should specify a local output directory using:
# BASE6_EXTERNAL_OUTPUT_DIR

external_output_dir <- Sys.getenv(
  "BASE6_EXTERNAL_OUTPUT_DIR"
)

if (external_output_dir == "") {
  stop(
    paste0(
      "External-cohort outputs are not written to the public repository. ",
      "Set the environment variable BASE6_EXTERNAL_OUTPUT_DIR to an ",
      "authorized local directory for restricted analysis outputs."
    )
  )
}

if (!dir.exists(external_output_dir)) {
  dir.create(
    external_output_dir,
    recursive = TRUE
  )
}


# Aggregate audit outputs

write.csv(
  base6_distribution,
  file.path(
    external_output_dir,
    "Chinese_Base6_predictor_audit.csv"
  ),
  row.names = FALSE
)


write.csv(
  invalid_med,
  file.path(
    external_output_dir,
    "Chinese_medication_coding_audit.csv"
  ),
  row.names = FALSE
)


# Restricted individual-level analysis objects
# These files must not be uploaded to the public repository.

save(
  chinese_data_raw,
  chinese_validation,
  base6_distribution,
  invalid_med,
  med_binary_vars,
  file = file.path(
    external_output_dir,
    "Chinese_HF_EXTERNAL_STEP9A_CLEAN.RData"
  )
)


cat("\n============================================\n")
cat("STEP 9A 完成\n")
cat("============================================\n")
