# ============================================================
# STEP 11B
# GWTG-HF score + three-model comparison
#
# Uses:
#   - VERIFIED Chinese GWTG-clean dataset
#   - Existing Base6 external predictions
#   - Existing MAGGIC score data
#
# Important:
#   GWTG-HF was originally developed for IN-HOSPITAL mortality.
#   Here it is used as a secondary comparator for ~1-year
#   mortality discrimination only.
# ============================================================

rm(list = ls())
gc()

library(dplyr)
library(tibble)
library(survival)
library(timeROC)

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
cat("STEP 11B: GWTG-HF scoring and comparison\n")
cat("============================================\n\n")


# ============================================================
# 1. Load verified GWTG-clean data
# ============================================================

load(
  "Chinese_GWTG_STEP11A2_VERIFIED_CLEAN.RData"
)

stopifnot(
  exists("gwtg_clean"),
  nrow(gwtg_clean) == 348
)

cat("Verified GWTG data loaded.\n")
cat("N =", nrow(gwtg_clean), "\n\n")


# ============================================================
# 2. GWTG-HF component scoring functions
# ============================================================

score_age <- function(x) {

  dplyr::case_when(
    x < 20  ~ 0,
    x < 30  ~ 3,
    x < 40  ~ 6,
    x < 50  ~ 8,
    x < 60  ~ 11,
    x < 70  ~ 14,
    x < 80  ~ 17,
    x < 90  ~ 19,
    x < 100 ~ 22,
    x < 110 ~ 25,
    TRUE    ~ 28
  )
}


score_hr <- function(x) {

  dplyr::case_when(
    x < 80  ~ 0,
    x < 85  ~ 1,
    x < 90  ~ 3,
    x < 95  ~ 4,
    x < 100 ~ 5,
    x < 105 ~ 6,
    TRUE    ~ 8
  )
}


score_sbp <- function(x) {

  dplyr::case_when(
    x < 60  ~ 28,
    x < 70  ~ 26,
    x < 80  ~ 24,
    x < 90  ~ 23,
    x < 100 ~ 21,
    x < 110 ~ 19,
    x < 120 ~ 17,
    x < 130 ~ 15,
    x < 140 ~ 13,
    x < 150 ~ 11,
    x < 160 ~ 9,
    x < 170 ~ 8,
    x < 180 ~ 6,
    x < 190 ~ 4,
    x < 200 ~ 2,
    TRUE    ~ 0
  )
}


score_na <- function(x) {

  dplyr::case_when(
    x <= 130 ~ 4,
    x <= 133 ~ 3,
    x <= 136 ~ 2,
    x <= 138 ~ 1,
    TRUE     ~ 0
  )
}


score_bun <- function(x) {

  dplyr::case_when(
    x < 10  ~ 0,
    x < 20  ~ 2,
    x < 30  ~ 4,
    x < 40  ~ 6,
    x < 50  ~ 8,
    x < 60  ~ 9,
    x < 70  ~ 11,
    x < 80  ~ 13,
    x < 90  ~ 15,
    x < 100 ~ 17,
    x < 110 ~ 19,
    x < 120 ~ 21,
    x < 130 ~ 23,
    x < 140 ~ 25,
    x < 150 ~ 27,
    TRUE    ~ 28
  )
}


# ============================================================
# 3. Calculate GWTG-HF score
#
# Chinese cohort:
#   all patients are nonblack -> +3 uniformly
#
# COPD:
#   assumed 0/1 in column COPD
# ============================================================

gwtg_score_data <- gwtg_clean %>%

  mutate(

    gwtg_age =
      score_age(age),

    gwtg_hr =
      score_hr(HR_gwtg),

    gwtg_sbp =
      score_sbp(SBP),

    gwtg_bun =
      score_bun(BUN_mgdl),

    gwtg_na =
      score_na(Na_gwtg),

    gwtg_nonblack =
      3,

    gwtg_copd =
      ifelse(
        COPD == 1,
        2,
        0
      ),

    gwtg_score =
      gwtg_age +
      gwtg_hr +
      gwtg_sbp +
      gwtg_bun +
      gwtg_na +
      gwtg_nonblack +
      gwtg_copd
  )


# ============================================================
# 4. Score completeness audit
# ============================================================

cat("============================================\n")
cat("A. GWTG score completeness\n")
cat("============================================\n")

cat(
  "Missing GWTG score =",
  sum(is.na(gwtg_score_data$gwtg_score)),
  "\n"
)

cat(
  "Complete GWTG score =",
  sum(!is.na(gwtg_score_data$gwtg_score)),
  "/",
  nrow(gwtg_score_data),
  "\n"
)

cat("\nScore distribution:\n")

print(
  summary(
    gwtg_score_data$gwtg_score
  )
)


# ============================================================
# 5. Component audit
# ============================================================

cat("\n============================================\n")
cat("B. GWTG component summaries\n")
cat("============================================\n")

component_summary <- gwtg_score_data %>%

  summarise(

    age_min = min(gwtg_age, na.rm = TRUE),
    age_max = max(gwtg_age, na.rm = TRUE),

    hr_min = min(gwtg_hr, na.rm = TRUE),
    hr_max = max(gwtg_hr, na.rm = TRUE),

    sbp_min = min(gwtg_sbp, na.rm = TRUE),
    sbp_max = max(gwtg_sbp, na.rm = TRUE),

    bun_min = min(gwtg_bun, na.rm = TRUE),
    bun_max = max(gwtg_bun, na.rm = TRUE),

    na_min = min(gwtg_na, na.rm = TRUE),
    na_max = max(gwtg_na, na.rm = TRUE),

    copd_min = min(gwtg_copd, na.rm = TRUE),
    copd_max = max(gwtg_copd, na.rm = TRUE),

    score_min = min(gwtg_score, na.rm = TRUE),
    score_max = max(gwtg_score, na.rm = TRUE)
  )

print(component_summary)


# ============================================================
# 6. Extreme-sodium scoring audit
#
# Patient identifiers and individual-level source values are not printed
# in the public workflow. Audit the minimum observed sodium value only.
# ============================================================

cat("\n============================================\n")
cat("C. Extreme sodium score audit\n")
cat("============================================\n")

extreme_na_audit <- gwtg_score_data %>%
  arrange(Na_gwtg) %>%
  slice(1) %>%
  transmute(
    Na_gwtg,
    sodium_points = gwtg_na,
    gwtg_score
  )

print(extreme_na_audit)

# ============================================================
# 7. Load Base6 external-validation results
# ============================================================

base6_file <- "Chinese_HF_EXTERNAL_STEP9B_VALIDATION.RData"

stopifnot(
  file.exists(base6_file)
)

base6_env <- new.env()

load(
  base6_file,
  envir = base6_env
)

cat("\n============================================\n")
cat("D0. Base6 file contents\n")
cat("============================================\n")

cat(
  "Loaded:",
  base6_file,
  "\n"
)

print(
  ls(base6_env)
)


# ============================================================
# 8. Load MAGGIC comparison results
# ============================================================

maggic_file <- "Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData"

stopifnot(
  file.exists(maggic_file)
)

maggic_env <- new.env()

load(
  maggic_file,
  envir = maggic_env
)

cat("\nMAGGIC file contents:\n")

print(
  ls(maggic_env)
)


# ============================================================
# 9. Helper: find data frame in an environment
# ============================================================

get_dataframes <- function(env) {

  nms <- ls(env)

  nms[
    sapply(
      nms,
      function(x) {
        is.data.frame(
          get(
            x,
            envir = env
          )
        )
      }
    )
  ]
}


base6_df_names <-
  get_dataframes(base6_env)

maggic_df_names <-
  get_dataframes(maggic_env)


cat("\nBase6 data frames:\n")
print(base6_df_names)

cat("\nMAGGIC data frames:\n")
print(maggic_df_names)


# ============================================================
# 10. Display columns to make matching transparent
# ============================================================

cat("\n============================================\n")
cat("D1. Base6 dataframe columns\n")
cat("============================================\n")

for (nm in base6_df_names) {

  obj <- get(
    nm,
    envir = base6_env
  )

  cat("\n---", nm, "---\n")

  print(
    names(obj)
  )
}


cat("\n============================================\n")
cat("D2. MAGGIC dataframe columns\n")
cat("============================================\n")

for (nm in maggic_df_names) {

  obj <- get(
    nm,
    envir = maggic_env
  )

  cat("\n---", nm, "---\n")

  print(
    names(obj)
  )
}


# ============================================================
# 11. Automatically locate Base6 dataset + prediction marker
# ============================================================

base6_marker_candidates <- c(
  "lp_final",
  "LP_final",
  "lp",
  "LP",
  "base6_lp",
  "pred_lp",
  "lp_shrunk",
  "LP_shrunk",
  "risk_1y",
  "pred_1y",
  "pred_risk_1y",
  "pred_risk"
)

base6_obj <- NULL
base6_obj_name <- NULL
base6_col <- NULL


for (nm in base6_df_names) {

  obj <- get(
    nm,
    envir = base6_env
  )

  # Must have patient ID
  if (!"ID" %in% names(obj)) {
    next
  }

  hit <- base6_marker_candidates[
    base6_marker_candidates %in%
      names(obj)
  ]

  if (length(hit) > 0) {

    base6_obj <-
      obj

    base6_obj_name <-
      nm

    base6_col <-
      hit[1]

    break
  }
}


if (is.null(base6_obj)) {

  stop(
    "Base6 dataframe found, but prediction/LP column ",
    "could not be identified. See D1 column output."
  )
}


cat("\n============================================\n")
cat("D3. Base6 object identified\n")
cat("============================================\n")

cat(
  "Object =",
  base6_obj_name,
  "\n"
)

cat(
  "Marker =",
  base6_col,
  "\n"
)

cat(
  "N =",
  nrow(base6_obj),
  "\n"
)


# ============================================================
# 12. Automatically locate MAGGIC dataset
# ============================================================

maggic_marker_candidates <- c(
  "maggic_score",
  "MAGGIC_score",
  "score_maggic",
  "maggic_original",
  "maggic_score_original"
)

maggic_obj <- NULL
maggic_obj_name <- NULL
maggic_col <- NULL


for (nm in maggic_df_names) {

  obj <- get(
    nm,
    envir = maggic_env
  )

  if (!"ID" %in% names(obj)) {
    next
  }

  hit <- maggic_marker_candidates[
    maggic_marker_candidates %in%
      names(obj)
  ]

  if (length(hit) > 0) {

    maggic_obj <-
      obj

    maggic_obj_name <-
      nm

    maggic_col <-
      hit[1]

    break
  }
}


if (is.null(maggic_obj)) {

  stop(
    "MAGGIC dataframe found, but score column ",
    "could not be identified. See D2 column output."
  )
}


cat("\n============================================\n")
cat("D4. MAGGIC object identified\n")
cat("============================================\n")

cat(
  "Object =",
  maggic_obj_name,
  "\n"
)

cat(
  "Marker =",
  maggic_col,
  "\n"
)

cat(
  "N =",
  nrow(maggic_obj),
  "\n"
)


# ============================================================
# 13. Outcome variables
#
# Use outcome/time from current corrected master dataset,
# rather than relying on old saved objects.
# ============================================================

outcome_data <- gwtg_score_data %>%
  transmute(
    ID = ID,
    status = status,
    time = time
  )


# ============================================================
# 14. Build common three-model cohort
# ============================================================

base6_compact <- base6_obj %>%
  transmute(
    ID = ID,
    base6_marker =
      .data[[base6_col]]
  )


maggic_compact <- maggic_obj %>%
  transmute(
    ID = ID,
    maggic_score =
      .data[[maggic_col]]
  )


gwtg_compact <- gwtg_score_data %>%
  transmute(
    ID = ID,
    gwtg_score = gwtg_score
  )


comparison3 <- outcome_data %>%

  inner_join(
    base6_compact,
    by = "ID"
  ) %>%

  inner_join(
    maggic_compact,
    by = "ID"
  ) %>%

  inner_join(
    gwtg_compact,
    by = "ID"
  ) %>%

  filter(
    complete.cases(
      status,
      time,
      base6_marker,
      maggic_score,
      gwtg_score
    )
  )


cat("\n============================================\n")
cat("D. Common three-model cohort\n")
cat("============================================\n")

cat(
  "Common n =",
  nrow(comparison3),
  "\n"
)

cat(
  "Deaths =",
  sum(
    comparison3$status == 1
  ),
  "\n"
)

cat(
  "Censored/alive =",
  sum(
    comparison3$status == 0
  ),
  "\n"
)

cat("\nMarker summaries:\n")

print(
  summary(
    comparison3[
      c(
        "base6_marker",
        "maggic_score",
        "gwtg_score"
      )
    ]
  )
)


# ============================================================
# 15. Helper: 364-day time-dependent AUC
# ============================================================

get_auc364 <- function(marker) {

  roc <- timeROC(
    T = comparison3$time,
    delta = comparison3$status,
    marker = marker,
    cause = 1,
    weighting = "marginal",
    times = 364,
    iid = TRUE
  )

  idx <- which(
    roc$times == 364
  )

  if (length(idx) == 0) {

    return(
      tail(
        roc$AUC,
        1
      )
    )
  }

  roc$AUC[idx]
}


auc_base6 <- get_auc364(
  comparison3$base6_marker
)

auc_maggic <- get_auc364(
  comparison3$maggic_score
)

auc_gwtg <- get_auc364(
  comparison3$gwtg_score
)


# ============================================================
# 16. Helper: Harrell C
# ============================================================

get_cindex <- function(marker) {

  fit <- survival::coxph(
    Surv(time, status) ~ marker,
    data = comparison3
  )

  as.numeric(
    summary(fit)$concordance[1]
  )
}


c_base6 <- get_cindex(
  comparison3$base6_marker
)

c_maggic <- get_cindex(
  comparison3$maggic_score
)

c_gwtg <- get_cindex(
  comparison3$gwtg_score
)


# ============================================================
# 17. Performance summary
# ============================================================

performance_summary <- tibble(
  model = c(
    "Base6",
    "MAGGIC",
    "GWTG-HF"
  ),

  AUC_364 = c(
    auc_base6,
    auc_maggic,
    auc_gwtg
  ),

  Harrell_C = c(
    c_base6,
    c_maggic,
    c_gwtg
  )
)


cat("\n============================================\n")
cat("E. Three-model discrimination\n")
cat("============================================\n")

print(
  as.data.frame(
    performance_summary
  ),
  row.names = FALSE,
  digits = 6
)


# ============================================================
# 18. Pairwise point-estimate differences
# ============================================================

pairwise_summary <- tibble(

  comparison = c(
    "Base6 - MAGGIC",
    "Base6 - GWTG-HF",
    "GWTG-HF - MAGGIC"
  ),

  delta_AUC = c(
    auc_base6 - auc_maggic,
    auc_base6 - auc_gwtg,
    auc_gwtg - auc_maggic
  ),

  delta_C = c(
    c_base6 - c_maggic,
    c_base6 - c_gwtg,
    c_gwtg - c_maggic
  )
)


cat("\n============================================\n")
cat("F. Pairwise point-estimate differences\n")
cat("============================================\n")

print(
  as.data.frame(
    pairwise_summary
  ),
  row.names = FALSE,
  digits = 6
)


# ============================================================
# 19. Save
# ============================================================

save(
  gwtg_score_data,
  comparison3,
  performance_summary,
  pairwise_summary,
  auc_base6,
  auc_maggic,
  auc_gwtg,
  c_base6,
  c_maggic,
  c_gwtg,
  base6_obj_name,
  base6_col,
  maggic_obj_name,
  maggic_col,
  file =
    "Chinese_GWTG_STEP11B_SCORE_COMPARISON.RData"
)

write.csv(
  performance_summary,
  "Chinese_GWTG_STEP11B_PERFORMANCE.csv",
  row.names = FALSE
)

write.csv(
  pairwise_summary,
  "Chinese_GWTG_STEP11B_PAIRWISE_POINT_ESTIMATES.csv",
  row.names = FALSE
)

cat("\n============================================\n")
cat("STEP 11B completed\n")
cat("============================================\n")
