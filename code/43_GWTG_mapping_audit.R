# ============================================================
# STEP 11A
# GWTG-HF variable mapping / unit / range audit
#
# Chinese hospitalized HF external cohort
#
# Variables:
#   age
#   SBP
#   BUN
#   Plus = pulse / heart rate
#   Na
#   COPD
#   race = all non-Black
#
# IMPORTANT:
# Do NOT calculate score yet.
# First determine BUN units and data validity.
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
cat("STEP 11A: GWTG-HF mapping audit\n")
cat("============================================\n")


# ============================================================
# 1. Required variables
# ============================================================

gwtg_vars <- c(
  "age",
  "SBP",
  "BUN",
  "Plus",
  "Na",
  "COPD"
)


missing_vars <- setdiff(
  gwtg_vars,
  names(chinese_validation)
)


cat("\n============================================\n")
cat("A. Required-variable check\n")
cat("============================================\n")

if (length(missing_vars) == 0) {

  cat(
    "All required GWTG-HF variables are present.\n"
  )

} else {

  print(
    missing_vars
  )

  stop(
    "Required GWTG-HF variables missing."
  )
}


# ============================================================
# 2. Variable classes
# ============================================================

cat("\n============================================\n")
cat("B. Variable classes\n")
cat("============================================\n")

for (v in gwtg_vars) {

  cat(
    v,
    ":",
    class(
      chinese_validation[[v]]
    ),
    "\n"
  )
}


# ============================================================
# 3. Convert to audit numeric variables
# ============================================================

gwtg_audit <- chinese_validation %>%
  mutate(

    gwtg_age =
      suppressWarnings(
        as.numeric(age)
      ),

    gwtg_sbp =
      suppressWarnings(
        as.numeric(SBP)
      ),

    gwtg_bun_raw =
      suppressWarnings(
        as.numeric(BUN)
      ),

    gwtg_hr =
      suppressWarnings(
        as.numeric(Plus)
      ),

    gwtg_na =
      suppressWarnings(
        as.numeric(Na)
      ),

    gwtg_copd =
      suppressWarnings(
        as.numeric(COPD)
      )
  )


# ============================================================
# 4. Continuous-variable audit
# ============================================================

continuous_gwtg <- c(
  "gwtg_age",
  "gwtg_sbp",
  "gwtg_bun_raw",
  "gwtg_hr",
  "gwtg_na"
)


cat("\n============================================\n")
cat("C. Continuous-variable distributions\n")
cat("============================================\n")


for (v in continuous_gwtg) {

  cat(
    "\n--- ",
    v,
    " ---\n",
    sep = ""
  )

  print(
    summary(
      gwtg_audit[[v]]
    )
  )

  cat(
    "Missing:",
    sum(
      is.na(
        gwtg_audit[[v]]
      )
    ),
    "\n"
  )
}


# ============================================================
# 5. Detailed quantiles
# ============================================================

distribution_table <- do.call(

  rbind,

  lapply(

    continuous_gwtg,

    function(v) {

      x <- gwtg_audit[[v]]

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

        p05 =
          as.numeric(
            quantile(
              x,
              0.05,
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

        p95 =
          as.numeric(
            quantile(
              x,
              0.95,
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


cat("\n============================================\n")
cat("D. Detailed distributions\n")
cat("============================================\n")

print(
  distribution_table,
  row.names = FALSE
)


# ============================================================
# 6. COPD
# ============================================================

cat("\n============================================\n")
cat("E. COPD coding\n")
cat("============================================\n")

print(
  table(
    gwtg_audit$gwtg_copd,
    useNA = "ifany"
  )
)


# ============================================================
# 7. BUN unit diagnostics
#
# Original GWTG-HF uses BUN in mg/dL.
#
# Typical possibilities:
#
# Chinese hospital field BUN may be:
#
#   mmol/L:
#       often approximately 2.5 - 10 mmol/L
#
#   mg/dL:
#       often approximately 7 - 25 mg/dL
#
# Conversion:
#
#   BUN mg/dL = BUN mmol/L * 2.801
#
# DO NOT convert automatically yet.
# ============================================================

cat("\n============================================\n")
cat("F. BUN unit diagnostic\n")
cat("============================================\n")

cat("\nRaw BUN:\n")

print(
  summary(
    gwtg_audit$gwtg_bun_raw
  )
)


cat(
  "\nBUN raw < 3:",
  sum(
    gwtg_audit$gwtg_bun_raw < 3,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "BUN raw 3-15:",
  sum(
    gwtg_audit$gwtg_bun_raw >= 3 &
      gwtg_audit$gwtg_bun_raw <= 15,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "BUN raw >15:",
  sum(
    gwtg_audit$gwtg_bun_raw > 15,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "BUN raw >50:",
  sum(
    gwtg_audit$gwtg_bun_raw > 50,
    na.rm = TRUE
  ),
  "\n"
)


# Candidate conversion if raw is mmol/L

gwtg_audit$gwtg_bun_if_mgdl <-
  gwtg_audit$gwtg_bun_raw *
  2.801


cat("\nIf raw BUN is mmol/L, converted mg/dL would be:\n")

print(
  summary(
    gwtg_audit$gwtg_bun_if_mgdl
  )
)


# ============================================================
# 8. Range against original GWTG-HF score boundaries
# ============================================================

cat("\n============================================\n")
cat("G. Original-score range audit\n")
cat("============================================\n")


cat(
  "Age <19:",
  sum(
    gwtg_audit$gwtg_age < 19,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Age >110:",
  sum(
    gwtg_audit$gwtg_age > 110,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nSBP <50:",
  sum(
    gwtg_audit$gwtg_sbp < 50,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "SBP >=200:",
  sum(
    gwtg_audit$gwtg_sbp >= 200,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nHR <79:",
  sum(
    gwtg_audit$gwtg_hr < 79,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "HR >=105:",
  sum(
    gwtg_audit$gwtg_hr >= 105,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nNa <130:",
  sum(
    gwtg_audit$gwtg_na < 130,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Na >=139:",
  sum(
    gwtg_audit$gwtg_na >= 139,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 9. Potentially implausible values
# ============================================================

implausible <- gwtg_audit %>%
  mutate(
    row_id =
      row_number()
  ) %>%
  filter(

    gwtg_age < 18 |
      gwtg_age > 110 |

      gwtg_sbp < 50 |
      gwtg_sbp > 260 |

      gwtg_hr < 25 |
      gwtg_hr > 220 |

      gwtg_na < 110 |
      gwtg_na > 165 |

      gwtg_bun_raw < 0
  ) %>%
  select(

    row_id,
    ID,

    gwtg_age,
    gwtg_sbp,
    gwtg_bun_raw,
    gwtg_hr,
    gwtg_na,
    gwtg_copd
  )


cat("\n============================================\n")
cat("H. Potentially implausible cases\n")
cat("============================================\n")

cat(
  "Number:",
  nrow(
    implausible
  ),
  "\n"
)

print(
  as.data.frame(
    implausible
  ),
  row.names = FALSE
)


# ============================================================
# 10. Complete-case availability
# ============================================================

gwtg_complete_vars <- c(
  "gwtg_age",
  "gwtg_sbp",
  "gwtg_bun_raw",
  "gwtg_hr",
  "gwtg_na",
  "gwtg_copd"
)


gwtg_audit <- gwtg_audit %>%
  mutate(

    gwtg_complete =
      complete.cases(
        across(
          all_of(
            gwtg_complete_vars
          )
        )
      )
  )


cat("\n============================================\n")
cat("I. GWTG-HF complete cases\n")
cat("============================================\n")

cat(
  "Total n =",
  nrow(
    gwtg_audit
  ),
  "\n"
)

cat(
  "Complete n =",
  sum(
    gwtg_audit$gwtg_complete
  ),
  "\n"
)

cat(
  "Complete proportion =",
  round(
    mean(
      gwtg_audit$gwtg_complete
    ) *
      100,
    2
  ),
  "%\n"
)


# ============================================================
# 11. Missingness
# ============================================================

missing_summary <- data.frame(

  variable =
    gwtg_complete_vars,

  missing =
    sapply(

      gwtg_complete_vars,

      function(v) {

        sum(
          is.na(
            gwtg_audit[[v]]
          )
        )
      }
    )
)


cat("\n============================================\n")
cat("J. Missingness\n")
cat("============================================\n")

print(
  missing_summary,
  row.names = FALSE
)


# ============================================================
# 12. Save
# ============================================================

write.csv(

  distribution_table,

  "Chinese_GWTG_variable_distribution.csv",

  row.names = FALSE
)


write.csv(

  implausible,

  "Chinese_GWTG_implausible_cases.csv",

  row.names = FALSE
)


save(

  gwtg_audit,

  distribution_table,

  missing_summary,

  implausible,

  file =
    "Chinese_HF_STEP11A_GWTG_AUDIT.RData"
)


cat("\n============================================\n")
cat("STEP 11A 完成\n")
cat("============================================\n")
