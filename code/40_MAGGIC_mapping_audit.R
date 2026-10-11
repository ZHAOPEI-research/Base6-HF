# ============================================================
# STEP 10A
# MAGGIC variable mapping and coding audit
# Chinese external HF cohort
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
  "Chinese_HF_EXTERNAL_STEP9A_FINAL.RData"
)

cat("============================================\n")
cat("STEP 10A: MAGGIC mapping audit\n")
cat("============================================\n")


# ============================================================
# 1. Required MAGGIC variables
# ============================================================

maggic_vars <- c(
  "age",
  "sex",
  "LVEF",
  "NYHA",
  "BMI",
  "CR",
  "SBP",
  "DM",
  "COPD",
  "smoke",
  "Duration_HF（day）",
  "β-blocker",
  "ACEI/ARB"
)

missing_vars <- setdiff(
  maggic_vars,
  names(chinese_validation)
)


cat("\n============================================\n")
cat("A. Required-variable check\n")
cat("============================================\n")

if (length(missing_vars) == 0) {

  cat("All 13 MAGGIC variables are present.\n")

} else {

  print(missing_vars)

  stop(
    "MAGGIC variables missing."
  )
}


# ============================================================
# 2. Continuous variables
# ============================================================

cat("\n============================================\n")
cat("B. Continuous-variable audit\n")
cat("============================================\n")


continuous_vars <- c(
  "age",
  "LVEF",
  "BMI",
  "CR",
  "SBP",
  "Duration_HF（day）"
)


for (v in continuous_vars) {

  cat(
    "\n--- ",
    v,
    " ---\n",
    sep = ""
  )

  print(
    summary(
      chinese_validation[[v]]
    )
  )

  cat(
    "Missing:",
    sum(
      is.na(
        chinese_validation[[v]]
      )
    ),
    "\n"
  )
}


# ============================================================
# 3. Sex
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
# 4. NYHA
# ============================================================

cat("\n============================================\n")
cat("D. NYHA coding\n")
cat("============================================\n")

print(
  table(
    chinese_validation$NYHA,
    useNA = "ifany"
  )
)


# ============================================================
# 5. Diabetes
# ============================================================

cat("\n============================================\n")
cat("E. Diabetes coding\n")
cat("============================================\n")

print(
  table(
    chinese_validation$DM,
    useNA = "ifany"
  )
)


# ============================================================
# 6. COPD
# ============================================================

cat("\n============================================\n")
cat("F. COPD coding\n")
cat("============================================\n")

print(
  table(
    chinese_validation$COPD,
    useNA = "ifany"
  )
)


# ============================================================
# 7. Smoking
# ============================================================

cat("\n============================================\n")
cat("G. Smoking coding\n")
cat("============================================\n")

print(
  table(
    chinese_validation$smoke,
    useNA = "ifany"
  )
)


# ============================================================
# 8. Beta blocker
# ============================================================

cat("\n============================================\n")
cat("H. Beta-blocker coding\n")
cat("============================================\n")

print(
  table(
    chinese_validation$`β-blocker`,
    useNA = "ifany"
  )
)


# ============================================================
# 9. ACEI/ARB
# ============================================================

cat("\n============================================\n")
cat("I. ACEI/ARB coding\n")
cat("============================================\n")

print(
  table(
    chinese_validation$`ACEI/ARB`,
    useNA = "ifany"
  )
)


# ============================================================
# 10. ARNI context
#
# Important for interpretation:
# MAGGIC uses ACEI/ARB, not ARNI.
# ============================================================

cat("\n============================================\n")
cat("J. RAAS medication context\n")
cat("============================================\n")

cat("\nACEI/ARB:\n")

print(
  table(
    chinese_validation$`ACEI/ARB`,
    useNA = "ifany"
  )
)

cat("\nARNi:\n")

print(
  table(
    chinese_validation$ARNi,
    useNA = "ifany"
  )
)

cat("\nACEI/ARB/ANRi:\n")

print(
  table(
    chinese_validation$`ACEI/ARB/ANRi`,
    useNA = "ifany"
  )
)


# ============================================================
# 11. HF duration
#
# MAGGIC distinguishes HF duration
# <18 months vs >=18 months.
#
# 18 months approximated as 548 days
# for audit.
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    hf_duration_days =
      as.numeric(
        `Duration_HF（day）`
      ),

    hf_duration_ge18m =
      case_when(

        is.na(hf_duration_days) ~
          NA_real_,

        hf_duration_days >= 548 ~
          1,

        hf_duration_days < 548 ~
          0
      )
  )


cat("\n============================================\n")
cat("K. HF duration >=18 months\n")
cat("============================================\n")

print(
  table(
    chinese_validation$hf_duration_ge18m,
    useNA = "ifany"
  )
)


cat("\nHF duration days:\n")

print(
  summary(
    chinese_validation$hf_duration_days
  )
)


# ============================================================
# 12. Creatinine unit check
#
# MAGGIC literature generally uses μmol/L
# in scoring implementations.
# Our Chinese CR is already μmol/L.
# ============================================================

cat("\n============================================\n")
cat("L. Creatinine audit\n")
cat("============================================\n")

print(
  summary(
    chinese_validation$CR
  )
)

cat(
  "Missing CR:",
  sum(
    is.na(
      chinese_validation$CR
    )
  ),
  "\n"
)


# ============================================================
# 13. Complete MAGGIC cases
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    maggic_complete =
      complete.cases(
        across(
          all_of(
            maggic_vars
          )
        )
      )
  )


cat("\n============================================\n")
cat("M. MAGGIC complete-case availability\n")
cat("============================================\n")

cat(
  "Total n =",
  nrow(
    chinese_validation
  ),
  "\n"
)

cat(
  "Complete MAGGIC =",
  sum(
    chinese_validation$maggic_complete
  ),
  "\n"
)

cat(
  "Complete proportion =",
  round(
    mean(
      chinese_validation$maggic_complete
    ) *
      100,
    2
  ),
  "%\n"
)


# ============================================================
# 14. Missingness by MAGGIC variable
# ============================================================

missing_summary <- data.frame(

  variable =
    maggic_vars,

  missing =
    sapply(

      maggic_vars,

      function(v) {

        sum(
          is.na(
            chinese_validation[[v]]
          )
        )
      }
    )
)


cat("\n============================================\n")
cat("N. MAGGIC variable missingness\n")
cat("============================================\n")

print(
  missing_summary,
  row.names = FALSE
)


# ============================================================
# 15. Save audit
# ============================================================

write.csv(

  missing_summary,

  "Chinese_MAGGIC_missingness_audit.csv",

  row.names = FALSE
)


save(

  chinese_validation,

  maggic_vars,

  missing_summary,

  file =
    "Chinese_HF_STEP10A_MAGGIC_AUDIT.RData"
)


cat("\n============================================\n")
cat("STEP 10A 完成\n")
cat("============================================\n")
