# ============================================================
# STEP 11D-A
# Audit and recover 1-year predicted risks for calibration
#
# Goal:
#   Base6 vs MAGGIC 1-year calibration
#
# NOTE:
#   GWTG-HF is retained for discrimination comparison only,
#   because gwtg_score is not itself a 1-year predicted
#   probability.
# ============================================================


# ============================================================
# 0. Clean environment
# ============================================================

rm(list = ls())



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

# ============================================================
# 1. Load latest STEP 11B common cohort
# ============================================================

load("Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData")

if (!exists("comparison3")) {
  stop("ERROR: comparison3 not found.")
}

cal_base <- comparison3


cat("\n")
cat("============================================\n")
cat("STEP 11D-A: Calibration data audit\n")
cat("============================================\n\n")

cat("Common cohort n =", nrow(cal_base), "\n")
cat(
  "Deaths =",
  sum(cal_base$status == 1, na.rm = TRUE),
  "\n\n"
)

cat("Columns in common cohort:\n")
print(names(cal_base))
cat("\n")


# ============================================================
# 2. Confirm Base6 1-year risk
# ============================================================

if (!"base6_marker" %in% names(cal_base)) {
  stop("ERROR: base6_marker not found.")
}

cat("============================================\n")
cat("A. Base6 predicted 1-year risk\n")
cat("============================================\n\n")

print(
  summary(
    cal_base$base6_marker
  )
)

cat(
  "\nMissing Base6 risk =",
  sum(is.na(cal_base$base6_marker)),
  "\n"
)

cat(
  "Outside [0,1] =",
  sum(
    cal_base$base6_marker < 0 |
      cal_base$base6_marker > 1,
    na.rm = TRUE
  ),
  "\n\n"
)


# ============================================================
# 3. Search available RData files for MAGGIC 1-year risk
# ============================================================

rdata_files <- list.files(
  pattern = "\\.(RData|rdata)$",
  ignore.case = TRUE,
  full.names = TRUE
)

cat("============================================\n")
cat("B. RData files found\n")
cat("============================================\n\n")

print(
  basename(
    rdata_files
  )
)

cat("\n")


# ============================================================
# 4. Inspect RData files in temporary environments
# ============================================================

candidate_list <- list()
candidate_counter <- 0

for (f in rdata_files) {

  e <- new.env(parent = emptyenv())

  load(
    f,
    envir = e
  )

  obj_names <- ls(envir = e)

  for (obj_name in obj_names) {

    obj <- get(
      obj_name,
      envir = e
    )

    if (
      is.data.frame(obj) ||
      inherits(obj, "tbl_df")
    ) {

      nm <- names(obj)

      if (
        all(
          c(
            "ID",
            "risk_maggic_1y"
          ) %in% nm
        )
      ) {

        candidate_counter <- candidate_counter + 1

        candidate_list[[candidate_counter]] <- data.frame(
          file = basename(f),
          object = obj_name,
          n = nrow(obj),

          risk_missing =
            sum(
              is.na(
                obj$risk_maggic_1y
              )
            ),

          risk_min =
            suppressWarnings(
              min(
                obj$risk_maggic_1y,
                na.rm = TRUE
              )
            ),

          risk_median =
            suppressWarnings(
              median(
                obj$risk_maggic_1y,
                na.rm = TRUE
              )
            ),

          risk_max =
            suppressWarnings(
              max(
                obj$risk_maggic_1y,
                na.rm = TRUE
              )
            ),

          stringsAsFactors = FALSE
        )
      }
    }
  }
}


# ============================================================
# 5. Print candidate MAGGIC datasets
# ============================================================

cat("============================================\n")
cat("C. Candidate MAGGIC risk datasets\n")
cat("============================================\n\n")

if (length(candidate_list) == 0) {

  cat(
    "No dataframe containing both ID and ",
    "risk_maggic_1y was found.\n\n"
  )

} else {

  candidate_summary <- do.call(
    rbind,
    candidate_list
  )

  rownames(candidate_summary) <- NULL

  print(
    candidate_summary,
    row.names = FALSE
  )

  cat("\n")
}


# ============================================================
# 6. Search for all MAGGIC-related columns
# ============================================================

cat("============================================\n")
cat("D. Objects containing MAGGIC-related columns\n")
cat("============================================\n\n")

for (f in rdata_files) {

  e <- new.env(parent = emptyenv())

  load(
    f,
    envir = e
  )

  for (obj_name in ls(envir = e)) {

    obj <- get(
      obj_name,
      envir = e
    )

    if (
      is.data.frame(obj) ||
      inherits(obj, "tbl_df")
    ) {

      maggic_cols <- grep(
        "maggic|risk.*1y|1y.*risk",
        names(obj),
        value = TRUE,
        ignore.case = TRUE
      )

      if (length(maggic_cols) > 0) {

        cat(
          "File:",
          basename(f),
          "\n"
        )

        cat(
          "Object:",
          obj_name,
          "\n"
        )

        cat(
          "N:",
          nrow(obj),
          "\n"
        )

        cat(
          "Columns:",
          paste(
            maggic_cols,
            collapse = ", "
          ),
          "\n\n"
        )
      }
    }
  }
}


cat("============================================\n")
cat("STEP 11D-A completed\n")
cat("============================================\n")
