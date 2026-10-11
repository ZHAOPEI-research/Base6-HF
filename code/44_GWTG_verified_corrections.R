# ============================================================
# STEP 11A-2
# Verified correction of Na, BUN and pulse
#
# IMPORTANT:
#   Corrections are based on source-record verification.
#
#   Do NOT automatically convert physiologically extreme
#   but verified measurements to NA.
#
# ============================================================

rm(list = ls())
gc()

library(readxl)
library(dplyr)
library(tibble)




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

cat("============================================\n")
cat("STEP 11A-2: Verified data corrections\n")
cat("============================================\n\n")


# ============================================================
# 1. Read Chinese external dataset
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

dat <- read_excel(
  external_data_path
)


dat <- as.data.frame(dat)

cat(
  "N =",
  nrow(dat),
  "\n"
)

stopifnot(
  all(
    c(
      "ID",
      "Na",
      "BUN",
      "Plus"
    ) %in%
      names(dat)
  )
)


# ============================================================
# 2. Convert fields to numeric safely
# ============================================================

dat$Na <- suppressWarnings(
  as.numeric(
    dat$Na
  )
)

dat$BUN <- suppressWarnings(
  as.numeric(
    dat$BUN
  )
)

dat$Plus <- suppressWarnings(
  as.numeric(
    dat$Plus
  )
)


# ============================================================
# 3. Source-record verified values
#
# IMPORTANT:
# Patient-level correction/verification values are restricted and are
# intentionally not embedded in the public repository.
#
# Authorized users should provide a local CSV with the columns:
#   ID, variable, verified_value, note
#
# Expected variable names are: Na, BUN, Plus
# ============================================================

external_corrections_path <- Sys.getenv("BASE6_EXTERNAL_CORRECTIONS")

if (external_corrections_path == "") {
  stop(
    paste0(
      "Restricted source-verification values are not distributed with ",
      "this public repository. Set BASE6_EXTERNAL_CORRECTIONS to the ",
      "local path of the authorized correction/verification CSV."
    )
  )
}

if (!file.exists(external_corrections_path)) {
  stop(
    paste0(
      "BASE6_EXTERNAL_CORRECTIONS file does not exist: ",
      external_corrections_path
    )
  )
}

verified_values <- read.csv(
  external_corrections_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_correction_cols <- c(
  "ID",
  "variable",
  "verified_value",
  "note"
)

missing_correction_cols <- setdiff(
  required_correction_cols,
  names(verified_values)
)

if (length(missing_correction_cols) > 0) {
  stop(
    paste0(
      "Missing required correction columns: ",
      paste(missing_correction_cols, collapse = ", ")
    )
  )
}

verified_values <- verified_values[
  ,
  required_correction_cols,
  drop = FALSE
]

if (!all(verified_values$variable %in% c("Na", "BUN", "Plus"))) {
  stop("Correction file contains unsupported variable names.")
}

verified_values$verified_value <- suppressWarnings(
  as.numeric(verified_values$verified_value)
)

if (anyNA(verified_values$verified_value)) {
  stop("verified_value contains non-numeric or missing values.")
}

cat("\n============================================\n")
cat("A. Restricted verification dictionary loaded\n")
cat("============================================\n")
cat("Number of verification records =", nrow(verified_values), "\n")
cat("Variables =", paste(sort(unique(verified_values$variable)), collapse = ", "), "\n")
cat("Patient-level verification values are not printed.\n")

# ============================================================
# 4. Verify IDs exist
# ============================================================

missing_ids <- setdiff(
  unique(verified_values$ID),
  dat$ID
)

cat("\n============================================\n")
cat("B. Study-ID check\n")
cat("============================================\n")

if (length(missing_ids) == 0) {

  cat("PASS: all study IDs found.\n")

} else {

  cat("Missing study IDs:\n")
  print(missing_ids)

  stop("Some verified study IDs were not found.")
}


# ============================================================
# 5. Audit CURRENT Excel values
# ============================================================

verification_audit <- verified_values %>%

  rowwise() %>%

  mutate(

    current_value = {

      idx <- which(dat$ID == ID)

      as.numeric(
        dat[
          idx,
          variable
        ]
      )
    },

    difference =
      current_value -
      verified_value,

    matches_source =
      abs(
        current_value -
          verified_value
      ) < 1e-8

  ) %>%

  ungroup()


cat("\n============================================\n")
cat("C. Source-record verification\n")
cat("============================================\n")

print(
  as.data.frame(verification_audit),
  row.names = FALSE
)


# ============================================================
# 6. Stop if corrected Excel is inconsistent
# ============================================================

if (
  any(
    !verification_audit$matches_source
  )
) {

  cat("\nERROR: Some Excel values do not match source records:\n")

  print(
    as.data.frame(
      verification_audit %>%
        filter(
          !matches_source
        )
    ),
    row.names = FALSE
  )

  stop(
    "The authorized external dataset does not match verified source records."
  )

} else {

  cat(
    "\nPASS: all manually corrected Excel values match source records.\n"
  )
}


# ============================================================
# 7. No further correction required
#
# The Excel file itself is now the corrected master dataset.
# ============================================================

cat("\nNo additional value replacement performed.\n")
cat("Verified extreme measurements are retained as observed.\n")

# ============================================================
# 8. IMPORTANT:
#    Do NOT remove verified extreme values
# ============================================================

cat("\n============================================\n")
cat("E. Extreme-value audit AFTER verification\n")
cat("============================================\n")


cat("\nNa summary:\n")

print(
  summary(
    dat$Na
  )
)


cat("\nNa < 110 mmol/L:\n")

print(
  dat[
    !is.na(dat$Na) &
      dat$Na < 110,
    c(
      "ID",
      "Na"
    )
  ]
)


cat("\nBUN summary (mmol/L):\n")

print(
  summary(
    dat$BUN
  )
)


cat("\nPulse summary:\n")

print(
  summary(
    dat$Plus
  )
)


cat("\nPulse < 30 /min:\n")

print(
  dat[
    !is.na(dat$Plus) &
      dat$Plus < 30,
    c(
      "ID",
      "Plus"
    )
  ]
)


# ============================================================
# 9. Create GWTG variables
# ============================================================

gwtg_clean <- dat


# ------------------------------------------------------------
# BUN:
# Chinese dataset uses mmol/L
# Convert to mg/dL for original GWTG-HF score
#
# 1 mmol/L urea nitrogen ≈ 2.801 mg/dL
# ------------------------------------------------------------

gwtg_clean$BUN_mgdl <-
  gwtg_clean$BUN *
  2.801


# ------------------------------------------------------------
# Sodium
#
# Retain source-verified extreme sodium values.
# Original GWTG scoring naturally assigns:
# Na <=130 -> 4 points.
# ------------------------------------------------------------

gwtg_clean$Na_gwtg <-
  gwtg_clean$Na


# ------------------------------------------------------------
# Heart rate
# ------------------------------------------------------------

gwtg_clean$HR_gwtg <-
  gwtg_clean$Plus


# ============================================================
# 10. Final GWTG input audit
# ============================================================

cat("\n============================================\n")
cat("F. Final GWTG input distributions\n")
cat("============================================\n")


cat("\nBUN mmol/L:\n")

print(
  summary(
    gwtg_clean$BUN
  )
)


cat("\nBUN mg/dL:\n")

print(
  summary(
    gwtg_clean$BUN_mgdl
  )
)


cat("\nNa mmol/L:\n")

print(
  summary(
    gwtg_clean$Na_gwtg
  )
)


cat("\nHR /min:\n")

print(
  summary(
    gwtg_clean$HR_gwtg
  )
)


# ============================================================
# 11. Check remaining missing data
# ============================================================

gwtg_input_vars <- c(
  "age",
  "SBP",
  "BUN_mgdl",
  "HR_gwtg",
  "Na_gwtg",
  "COPD"
)

cat("\n============================================\n")
cat("G. GWTG missing-data audit\n")
cat("============================================\n")

for (
  v in gwtg_input_vars
) {

  cat(
    sprintf(
      "%-12s missing = %d\n",
      v,
      sum(
        is.na(
          gwtg_clean[[v]]
        )
      )
    )
  )
}


gwtg_clean$gwtg_complete <-
  complete.cases(
    gwtg_clean[
      gwtg_input_vars
    ]
  )


cat(
  "\nComplete GWTG cases =",
  sum(
    gwtg_clean$gwtg_complete
  ),
  "/",
  nrow(
    gwtg_clean
  ),
  "\n"
)

cat(
  "Completeness =",
  round(
    100 *
      mean(
        gwtg_clean$gwtg_complete
      ),
    2
  ),
  "%\n"
)


# ============================================================
# 12. Save verified GWTG-clean dataset
# ============================================================

save(
  gwtg_clean,
  verification_audit,
  file =
    "Chinese_GWTG_STEP11A2_VERIFIED_CLEAN.RData"
)

write.csv(
  verification_audit,
  "Chinese_GWTG_STEP11A2_VERIFICATION_AUDIT.csv",
  row.names = FALSE
)

cat("\n============================================\n")
cat("STEP 11A-2 completed\n")
cat("Verified extreme values retained.\n")
cat("All source-record corrections confirmed.\n")
cat("GWTG complete cases =",
    sum(gwtg_clean$gwtg_complete),
    "/",
    nrow(gwtg_clean),
    "\n")
cat("============================================\n")
