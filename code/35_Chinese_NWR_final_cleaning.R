# ============================================================
# STEP 9A-2
# FINAL Chinese external NWR cleaning
#
# Principle:
# 1. Primary source = absolute neutrophil / WBC
# 2. neu% / 100 only used as fallback
# 3. No outcome-dependent/manual patient correction
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
cat("STEP 9A-2: FINAL Chinese NWR cleaning\n")
cat("============================================\n")


# ============================================================
# 1. Parse raw neutrophil values
#
# Handles examples such as:
# Decimal-comma numeric strings are normalized before conversion.
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    neutrophil_raw =
      trimws(
        as.character(
          neutrophil
        )
      ),

    neutrophil_final_numeric =
      suppressWarnings(
        as.numeric(
          gsub(
            ",",
            ".",
            neutrophil_raw,
            fixed = TRUE
          )
        )
      ),

    WBC_final_numeric =
      as.numeric(
        WBC
      ),

    neu_pct_final_numeric =
      as.numeric(
        `neu%`
      )
  )


# ============================================================
# 2. Construct two candidate NWR measures
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    # Development-model definition
    nwr_anc =
      ifelse(
        !is.na(neutrophil_final_numeric) &
          !is.na(WBC_final_numeric) &
          WBC_final_numeric > 0,
        neutrophil_final_numeric /
          WBC_final_numeric,
        NA_real_
      ),

    # QC / fallback definition
    nwr_pct =
      ifelse(
        !is.na(neu_pct_final_numeric) &
          neu_pct_final_numeric >= 0 &
          neu_pct_final_numeric <= 100,
        neu_pct_final_numeric /
          100,
        NA_real_
      )
  )


# ============================================================
# 3. Validity flags
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    anc_valid =
      !is.na(nwr_anc) &
      nwr_anc > 0 &
      nwr_anc <= 1,

    pct_valid =
      !is.na(nwr_pct) &
      nwr_pct > 0 &
      nwr_pct <= 1
  )


# ============================================================
# 4. FINAL NWR
#
# Priority:
#
# valid ANC/WBC
#    ↓
# valid neu% / 100
#    ↓
# NA
#
# This deliberately does NOT replace a plausible ANC/WBC
# merely because neu% disagrees with it.
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    nwr_final =
      case_when(

        anc_valid ~
          nwr_anc,

        !anc_valid &
          pct_valid ~
          nwr_pct,

        TRUE ~
          NA_real_
      ),

    nwr_source =
      case_when(

        anc_valid ~
          "ANC/WBC",

        !anc_valid &
          pct_valid ~
          "neu% fallback",

        TRUE ~
          "missing"
      )
  )


# ============================================================
# 5. Difference for QC
# ============================================================

chinese_validation <- chinese_validation %>%
  mutate(

    nwr_method_diff =
      nwr_anc -
      nwr_pct
  )


# ============================================================
# 6. Source summary
# ============================================================

cat("\n============================================\n")
cat("A. FINAL NWR source\n")
cat("============================================\n")

print(
  table(
    chinese_validation$nwr_source,
    useNA = "ifany"
  )
)


# ============================================================
# 7. Show fallback cases
# ============================================================

fallback_cases <- chinese_validation %>%
  mutate(
    row_id =
      row_number()
  ) %>%
  filter(
    nwr_source ==
      "neu% fallback"
  ) %>%
  select(

    row_id,
    ID,

    WBC,
    neutrophil,
    `neu%`,

    neutrophil_final_numeric,

    nwr_anc,
    nwr_pct,
    nwr_final,
    nwr_source
  )


cat("\n============================================\n")
cat("B. Cases using neu% fallback\n")
cat("============================================\n")

print(
  as.data.frame(
    fallback_cases
  ),
  row.names = FALSE
)


# ============================================================
# 8. Cases where both sources valid but disagree
#
# We still use ANC/WBC because it matches
# the development predictor definition.
# ============================================================

discordant_valid <- chinese_validation %>%
  mutate(
    row_id =
      row_number()
  ) %>%
  filter(

    anc_valid,
    pct_valid,

    abs(
      nwr_method_diff
    ) > 0.05
  ) %>%
  select(

    row_id,
    ID,

    WBC,
    neutrophil,
    `neu%`,

    nwr_anc,
    nwr_pct,
    nwr_method_diff,
    nwr_final
  ) %>%
  arrange(
    desc(
      abs(
        nwr_method_diff
      )
    )
  )


cat("\n============================================\n")
cat("C. Valid but discordant ANC/WBC vs neu%\n")
cat("============================================\n")

cat(
  "Number:",
  nrow(
    discordant_valid
  ),
  "\n"
)

print(
  as.data.frame(
    discordant_valid
  ),
  row.names = FALSE
)


# ============================================================
# 9. FINAL distribution
# ============================================================

cat("\n============================================\n")
cat("D. FINAL NWR distribution\n")
cat("============================================\n")

print(
  summary(
    chinese_validation$nwr_final
  )
)


cat("\nQuantiles:\n")

print(
  quantile(
    chinese_validation$nwr_final,
    probs = c(
      0,
      0.01,
      0.025,
      0.05,
      0.25,
      0.50,
      0.75,
      0.95,
      0.975,
      0.99,
      1
    ),
    na.rm = TRUE
  )
)


cat(
  "\nNWR <=0:",
  sum(
    chinese_validation$nwr_final <= 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "NWR >1:",
  sum(
    chinese_validation$nwr_final > 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "NWR missing:",
  sum(
    is.na(
      chinese_validation$nwr_final
    )
  ),
  "\n"
)


# ============================================================
# 10. Compare to frozen NHANES spline boundaries
# ============================================================

nwr_boundary_lower <-
  0.434782608695652

nwr_boundary_upper <-
  0.773607367475292


cat("\n============================================\n")
cat("E. External NWR vs development spline boundaries\n")
cat("============================================\n")


cat(
  "Below lower boundary:",
  sum(
    chinese_validation$nwr_final <
      nwr_boundary_lower,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Within boundaries:",
  sum(
    chinese_validation$nwr_final >=
      nwr_boundary_lower &
      chinese_validation$nwr_final <=
      nwr_boundary_upper,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Above upper boundary:",
  sum(
    chinese_validation$nwr_final >
      nwr_boundary_upper,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 11. Replace analysis NWR
#
# From now on:
# nwr = final cleaned NWR
# ============================================================

chinese_validation$nwr <-
  chinese_validation$nwr_final


# ============================================================
# 12. Recalculate Base6 complete cases
# ============================================================

base6_vars_final <- c(

  "age_model",
  "BMI",
  "RDW",
  "nwr",
  "egfr_2021",
  "med_cnt"
)


chinese_validation <- chinese_validation %>%
  mutate(

    base6_complete_final =
      complete.cases(
        across(
          all_of(
            base6_vars_final
          )
        )
      )
  )


cat("\n============================================\n")
cat("F. FINAL Base6 complete cases\n")
cat("============================================\n")


cat(
  "Total n =",
  nrow(
    chinese_validation
  ),
  "\n"
)


cat(
  "Complete Base6 =",
  sum(
    chinese_validation$base6_complete_final
  ),
  "\n"
)


cat(
  "Complete proportion =",
  round(
    mean(
      chinese_validation$base6_complete_final
    ) *
      100,
    2
  ),
  "%\n"
)


# ============================================================
# 13. Missing-variable detail
# ============================================================

missing_final <- data.frame(

  variable =
    base6_vars_final,

  missing =
    sapply(
      base6_vars_final,
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
cat("G. Remaining predictor missingness\n")
cat("============================================\n")


print(
  missing_final,
  row.names = FALSE
)


# ============================================================
# 14. Save
# ============================================================

write.csv(

  fallback_cases,

  "Chinese_NWR_FINAL_fallback_cases.csv",

  row.names = FALSE
)


write.csv(

  discordant_valid,

  "Chinese_NWR_FINAL_discordant_valid.csv",

  row.names = FALSE
)


save(

  chinese_validation,

  fallback_cases,

  discordant_valid,

  missing_final,

  base6_vars_final,

  file =
    "Chinese_HF_EXTERNAL_STEP9A_FINAL.RData"
)


cat("\n============================================\n")
cat("STEP 9A-2 完成\n")
cat("============================================\n")
