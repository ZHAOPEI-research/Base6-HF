# ============================================================
# STEP 10B
# MAGGIC original integer score
# + 1-year mortality prediction
# + head-to-head comparison with FINAL Base6
#
# Chinese hospitalized HF cohort
# ============================================================

rm(list = ls())

library(dplyr)
library(survival)

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
# 1. Load data
# ============================================================

load(
  "Chinese_HF_STEP10A_MAGGIC_AUDIT.RData"
)

load(
  "Chinese_HF_EXTERNAL_STEP9B_VALIDATION.RData"
)


cat("============================================\n")
cat("STEP 10B: MAGGIC score + head-to-head\n")
cat("============================================\n")


# ============================================================
# 2. Prepare MAGGIC dataset
# ============================================================

maggic_data <- chinese_validation %>%
  mutate(

    age_maggic =
      as.numeric(age),

    ef_maggic =
      as.numeric(LVEF),

    sbp_maggic =
      as.numeric(SBP),

    bmi_maggic =
      as.numeric(BMI),

    cr_maggic =
      as.numeric(CR),

    nyha_maggic =
      as.numeric(NYHA),

    duration_days =
      suppressWarnings(
        as.numeric(
          gsub(
            ",",
            ".",
            as.character(
              `Duration_HF（day）`
            ),
            fixed = TRUE
          )
        )
      ),

    # --------------------------------------------------------
    # sex:
    # 1 = male
    # 2 = female
    # --------------------------------------------------------

    male_maggic =
      ifelse(
        sex == 1,
        1,
        0
      ),

    smoker_maggic =
      as.numeric(smoke),

    dm_maggic =
      as.numeric(DM),

    copd_maggic =
      as.numeric(COPD),

    bb_use_maggic =
      as.numeric(
        `β-blocker`
      ),

    aceiarb_use_maggic =
      as.numeric(
        `ACEI/ARB`
      ),

    raas_modern_use =
      as.numeric(
        `ACEI/ARB/ANRi`
      )
  )


# ============================================================
# 3. Main EF points
#
# Original MAGGIC:
#
# EF <20    = 7
# 20-24     = 6
# 25-29     = 5
# 30-34     = 3
# 35-39     = 2
# >=40      = 0
# ============================================================

maggic_data <- maggic_data %>%
  mutate(

    pts_ef =
      case_when(

        ef_maggic < 20 ~ 7,

        ef_maggic >= 20 &
          ef_maggic < 25 ~ 6,

        ef_maggic >= 25 &
          ef_maggic < 30 ~ 5,

        ef_maggic >= 30 &
          ef_maggic < 35 ~ 3,

        ef_maggic >= 35 &
          ef_maggic < 40 ~ 2,

        ef_maggic >= 40 ~ 0,

        TRUE ~ NA_real_
      )
  )


# ============================================================
# 4. Age category
# ============================================================

maggic_data <- maggic_data %>%
  mutate(

    age_group =
      case_when(

        age_maggic < 55 ~ 1,

        age_maggic >= 55 &
          age_maggic < 60 ~ 2,

        age_maggic >= 60 &
          age_maggic < 65 ~ 3,

        age_maggic >= 65 &
          age_maggic < 70 ~ 4,

        age_maggic >= 70 &
          age_maggic < 75 ~ 5,

        age_maggic >= 75 &
          age_maggic < 80 ~ 6,

        age_maggic >= 80 ~ 7,

        TRUE ~ NA_real_
      ),

    ef_age_group =
      case_when(

        ef_maggic < 30 ~ 1,

        ef_maggic >= 30 &
          ef_maggic < 40 ~ 2,

        ef_maggic >= 40 ~ 3,

        TRUE ~ NA_real_
      )
  )


# ============================================================
# 5. Age × EF points
#
#                 Age
# EF<30:     0 1 2 4 6 8 10
# EF30-39:   0 2 4 6 8 10 13
# EF>=40:    0 3 5 7 9 12 15
# ============================================================

age_points_matrix <- matrix(

  c(

    # EF <30
    0, 1, 2, 4, 6, 8, 10,

    # EF 30-39
    0, 2, 4, 6, 8, 10, 13,

    # EF >=40
    0, 3, 5, 7, 9, 12, 15
  ),

  nrow = 3,

  byrow = TRUE
)


maggic_data$pts_age <- mapply(

  function(
    ef_group,
    age_group
  ) {

    if (
      is.na(ef_group) ||
      is.na(age_group)
    ) {
      return(NA_real_)
    }

    age_points_matrix[
      ef_group,
      age_group
    ]
  },

  maggic_data$ef_age_group,
  maggic_data$age_group
)


# ============================================================
# 6. SBP category
# ============================================================

maggic_data <- maggic_data %>%
  mutate(

    sbp_group =
      case_when(

        sbp_maggic < 110 ~ 1,

        sbp_maggic >= 110 &
          sbp_maggic < 120 ~ 2,

        sbp_maggic >= 120 &
          sbp_maggic < 130 ~ 3,

        sbp_maggic >= 130 &
          sbp_maggic < 140 ~ 4,

        sbp_maggic >= 140 &
          sbp_maggic < 150 ~ 5,

        sbp_maggic >= 150 ~ 6,

        TRUE ~ NA_real_
      )
  )


# ============================================================
# 7. SBP × EF points
#
#                 SBP
#             <110 110 120 130 140 >=150
#
# EF<30        5   4   3   2   1   0
# EF30-39      3   2   1   1   0   0
# EF>=40       2   1   1   0   0   0
# ============================================================

sbp_points_matrix <- matrix(

  c(

    5, 4, 3, 2, 1, 0,

    3, 2, 1, 1, 0, 0,

    2, 1, 1, 0, 0, 0
  ),

  nrow = 3,

  byrow = TRUE
)


maggic_data$pts_sbp <- mapply(

  function(
    ef_group,
    sbp_group
  ) {

    if (
      is.na(ef_group) ||
      is.na(sbp_group)
    ) {
      return(NA_real_)
    }

    sbp_points_matrix[
      ef_group,
      sbp_group
    ]
  },

  maggic_data$ef_age_group,
  maggic_data$sbp_group
)


# ============================================================
# 8. BMI points
#
# Original score:
#
# <15      6
# 15-19    5
# 20-24    3
# 25-29    2
# >=30     0
# ============================================================

maggic_data <- maggic_data %>%
  mutate(

    pts_bmi =
      case_when(

        bmi_maggic < 15 ~ 6,

        bmi_maggic >= 15 &
          bmi_maggic < 20 ~ 5,

        bmi_maggic >= 20 &
          bmi_maggic < 25 ~ 3,

        bmi_maggic >= 25 &
          bmi_maggic < 30 ~ 2,

        bmi_maggic >= 30 ~ 0,

        TRUE ~ NA_real_
      )
  )


# ============================================================
# 9. Creatinine points
#
# µmol/L
# ============================================================

maggic_data <- maggic_data %>%
  mutate(

    pts_cr =
      case_when(

        cr_maggic < 90 ~ 0,

        cr_maggic >= 90 &
          cr_maggic < 110 ~ 1,

        cr_maggic >= 110 &
          cr_maggic < 130 ~ 2,

        cr_maggic >= 130 &
          cr_maggic < 150 ~ 3,

        cr_maggic >= 150 &
          cr_maggic < 170 ~ 4,

        cr_maggic >= 170 &
          cr_maggic < 210 ~ 5,

        cr_maggic >= 210 &
          cr_maggic < 250 ~ 6,

        cr_maggic >= 250 ~ 8,

        TRUE ~ NA_real_
      )
  )


# ============================================================
# 10. NYHA points
#
# I   0
# II  2
# III 6
# IV  8
# ============================================================

maggic_data <- maggic_data %>%
  mutate(

    pts_nyha =
      case_when(

        nyha_maggic == 1 ~ 0,

        nyha_maggic == 2 ~ 2,

        nyha_maggic == 3 ~ 6,

        nyha_maggic == 4 ~ 8,

        TRUE ~ NA_real_
      )
  )


# ============================================================
# 11. Binary clinical points
#
# IMPORTANT:
#
# HF first diagnosed within past 18 months = +2
#
# NOT >=18 months.
# ============================================================

maggic_data <- maggic_data %>%
  mutate(

    pts_male =
      male_maggic * 1,

    pts_smoker =
      smoker_maggic * 1,

    pts_dm =
      dm_maggic * 3,

    pts_copd =
      copd_maggic * 2,

    hf_recent18 =
      ifelse(
        duration_days < 548,
        1,
        0
      ),

    pts_hf_recent =
      hf_recent18 * 2,

    pts_no_bb =
      ifelse(
        bb_use_maggic == 0,
        3,
        0
      ),

    # Strict original MAGGIC
    pts_no_aceiarb =
      ifelse(
        aceiarb_use_maggic == 0,
        1,
        0
      ),

    # Modernized sensitivity analysis only
    pts_no_raas_modern =
      ifelse(
        raas_modern_use == 0,
        1,
        0
      )
  )


# ============================================================
# 12. Total ORIGINAL MAGGIC integer score
# ============================================================

score_components <- c(

  "pts_ef",
  "pts_age",
  "pts_sbp",
  "pts_bmi",
  "pts_cr",
  "pts_nyha",

  "pts_male",
  "pts_smoker",
  "pts_dm",
  "pts_copd",
  "pts_hf_recent",
  "pts_no_bb",
  "pts_no_aceiarb"
)


maggic_data$maggic_score <- rowSums(

  maggic_data[
    score_components
  ],

  na.rm = FALSE
)


# ============================================================
# 13. Modernized RAAS sensitivity score
#
# Only difference:
#
# ACEI/ARB/ARNI counts as RAAS therapy.
#
# NOT the primary analysis.
# ============================================================

maggic_data$maggic_score_modern <-

  maggic_data$maggic_score -

  maggic_data$pts_no_aceiarb +

  maggic_data$pts_no_raas_modern


# ============================================================
# 14. Original MAGGIC Table 4:
# integer score -> 1-year mortality probability
# ============================================================

score_lookup <- data.frame(

  score = 0:50,

  risk_1y = c(

    0.015,
    0.016,
    0.018,
    0.020,
    0.022,
    0.024,
    0.027,
    0.029,
    0.032,
    0.036,

    0.039,
    0.043,
    0.048,
    0.052,
    0.058,
    0.063,
    0.070,
    0.077,
    0.084,
    0.093,

    0.102,
    0.111,
    0.122,
    0.134,
    0.147,
    0.160,
    0.175,
    0.191,
    0.209,
    0.227,

    0.248,
    0.269,
    0.292,
    0.316,
    0.342,
    0.369,
    0.398,
    0.427,
    0.458,
    0.490,

    0.523,
    0.557,
    0.591,
    0.625,
    0.659,
    0.692,
    0.725,
    0.757,
    0.787,
    0.816,

    0.842
  )
)


# ============================================================
# 15. Score range audit
# ============================================================

cat("\n============================================\n")
cat("A. MAGGIC score range\n")
cat("============================================\n")

print(
  summary(
    maggic_data$maggic_score
  )
)

cat(
  "\nScores <0:",
  sum(
    maggic_data$maggic_score < 0,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Scores >50:",
  sum(
    maggic_data$maggic_score > 50,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Missing scores:",
  sum(
    is.na(
      maggic_data$maggic_score
    )
  ),
  "\n"
)


# ============================================================
# 16. Map score to original 1-year mortality
# ============================================================

maggic_data$risk_maggic_1y <-
  score_lookup$risk_1y[
    match(
      maggic_data$maggic_score,
      score_lookup$score
    )
  ]


maggic_data$risk_maggic_modern_1y <-
  score_lookup$risk_1y[
    match(
      maggic_data$maggic_score_modern,
      score_lookup$score
    )
  ]


cat("\n============================================\n")
cat("B. MAGGIC score and predicted risk\n")
cat("============================================\n")

cat("\nMAGGIC score:\n")

print(
  summary(
    maggic_data$maggic_score
  )
)

cat("\nMAGGIC original 1-year predicted risk:\n")

print(
  summary(
    maggic_data$risk_maggic_1y
  )
)

cat(
  "\nMean predicted risk =",
  mean(
    maggic_data$risk_maggic_1y,
    na.rm = TRUE
  ),
  "\n"
)


# ============================================================
# 17. Score component summary
# ============================================================

component_summary <- data.frame(

  component =
    score_components,

  mean_points =
    sapply(
      score_components,
      function(v) {
        mean(
          maggic_data[[v]],
          na.rm = TRUE
        )
      }
    ),

  max_points =
    sapply(
      score_components,
      function(v) {
        max(
          maggic_data[[v]],
          na.rm = TRUE
        )
      }
    )
)


cat("\n============================================\n")
cat("C. MAGGIC score components\n")
cat("============================================\n")

print(
  component_summary,
  row.names = FALSE
)


# ============================================================
# 18. HF duration rule audit
# ============================================================

cat("\n============================================\n")
cat("D. HF diagnosed within past 18 months\n")
cat("============================================\n")

print(
  table(
    maggic_data$hf_recent18,
    useNA = "ifany"
  )
)

cat(
  "\nRecent HF (<548 d) receives +2 points.\n"
)


# ============================================================
# 19. Original vs modernized RAAS score
# ============================================================

cat("\n============================================\n")
cat("E. Original vs modernized RAAS definition\n")
cat("============================================\n")

cat(
  "Patients with changed MAGGIC score:",
  sum(
    maggic_data$maggic_score !=
      maggic_data$maggic_score_modern
  ),
  "\n"
)

cat("\nOriginal score:\n")

print(
  summary(
    maggic_data$maggic_score
  )
)

cat("\nModernized score:\n")

print(
  summary(
    maggic_data$maggic_score_modern
  )
)


# ============================================================
# 20. Merge with Base6 external predictions
#
# Head-to-head MUST use same patients.
# Base6 complete n previously = 346.
# ============================================================

base6_compare <- external_data %>%
  select(
    ID,
    time_day,
    status,
    lp_base6,
    risk_1y
  )


comparison_data <- maggic_data %>%

  select(

    ID,

    maggic_score,
    risk_maggic_1y,

    maggic_score_modern,
    risk_maggic_modern_1y

  ) %>%

  inner_join(
    base6_compare,
    by = "ID"
  )


cat("\n============================================\n")
cat("F. Common head-to-head cohort\n")
cat("============================================\n")

cat(
  "Common n =",
  nrow(
    comparison_data
  ),
  "\n"
)

cat(
  "Deaths =",
  sum(
    comparison_data$status == 1
  ),
  "\n"
)


# ============================================================
# 21. IPCW AUC helper
#
# 364 days because all administrative follow-up <=365.
# ============================================================

ipcw_auc <- function(
    time,
    status,
    marker,
    t
) {

  censor_event <- 1 - status

  km_censor <- survfit(
    Surv(
      time,
      censor_event
    ) ~ 1
  )


  G_left <- function(x) {

    as.numeric(
      summary(
        km_censor,
        times =
          pmax(
            x - 1e-8,
            0
          ),
        extend = TRUE
      )$surv
    )
  }


  cases <- which(
    status == 1 &
      time <= t
  )

  controls <- which(
    time > t
  )


  if (
    length(cases) == 0 ||
    length(controls) == 0
  ) {
    return(NA_real_)
  }


  case_weights <- sapply(

    cases,

    function(i) {

      1 /
        G_left(
          time[i]
        )
    }
  )


  numerator <- 0
  denominator <- 0


  for (j in seq_along(cases)) {

    i <- cases[j]

    diff <-
      marker[i] -
      marker[controls]


    cc <- ifelse(
      diff > 0,
      1,
      ifelse(
        diff == 0,
        0.5,
        0
      )
    )


    numerator <-
      numerator +
      case_weights[j] *
      sum(cc)


    denominator <-
      denominator +
      case_weights[j] *
      length(
        controls
      )
  }


  numerator /
    denominator
}


# ============================================================
# 22. Correct IPCW Brier helper at 365 days
# ============================================================

ipcw_brier <- function(
    time,
    status,
    pred_risk,
    t
) {

  censor_event <- 1 - status


  km_censor <- survfit(
    Surv(
      time,
      censor_event
    ) ~ 1
  )


  G_left <- function(x) {

    as.numeric(
      summary(
        km_censor,
        times =
          pmax(
            x - 1e-8,
            0
          ),
        extend = TRUE
      )$surv
    )
  }


  G_t_left <-
    G_left(t)


  contribution <- numeric(
    length(time)
  )


  event_rows <- which(
    status == 1 &
      time <= t
  )


  for (i in event_rows) {

    G_i <-
      G_left(
        time[i]
      )


    contribution[i] <-
      (1 - pred_risk[i])^2 /
      G_i
  }


  alive_rows <- which(
    status == 0 &
      time >= t
  )


  contribution[
    alive_rows
  ] <-
    pred_risk[
      alive_rows
    ]^2 /
    G_t_left


  mean(
    contribution
  )
}


# ============================================================
# 23. Head-to-head AUC
# ============================================================

auc_base6 <- ipcw_auc(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$lp_base6,
  364
)


auc_maggic <- ipcw_auc(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$maggic_score,
  364
)


auc_maggic_modern <- ipcw_auc(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$maggic_score_modern,
  364
)


cat("\n============================================\n")
cat("G. Head-to-head 364-day AUC\n")
cat("============================================\n")

cat(
  "Base6 =",
  auc_base6,
  "\n"
)

cat(
  "MAGGIC original =",
  auc_maggic,
  "\n"
)

cat(
  "MAGGIC modernized sensitivity =",
  auc_maggic_modern,
  "\n"
)


# ============================================================
# 24. Head-to-head Brier
# ============================================================

brier_base6 <- ipcw_brier(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$risk_1y,
  365
)


brier_maggic <- ipcw_brier(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$risk_maggic_1y,
  365
)


brier_maggic_modern <- ipcw_brier(

  comparison_data$time_day,
  comparison_data$status,
  comparison_data$risk_maggic_modern_1y,
  365
)


cat("\n============================================\n")
cat("H. Head-to-head 365-day Brier\n")
cat("============================================\n")

cat(
  "Base6 =",
  brier_base6,
  "\n"
)

cat(
  "MAGGIC original =",
  brier_maggic,
  "\n"
)

cat(
  "MAGGIC modernized sensitivity =",
  brier_maggic_modern,
  "\n"
)


# ============================================================
# 25. Observed risk in common cohort
# ============================================================

sf_common <- survfit(

  Surv(
    time_day,
    status
  ) ~ 1,

  data =
    comparison_data
)


observed_common <- 1 -
  as.numeric(

    summary(
      sf_common,
      times = 365,
      extend = TRUE
    )$surv
  )


cat("\n============================================\n")
cat("I. Overall calibration in common cohort\n")
cat("============================================\n")

cat(
  "Observed KM =",
  observed_common,
  "\n"
)

cat(
  "Base6 mean predicted =",
  mean(
    comparison_data$risk_1y
  ),
  "\n"
)

cat(
  "MAGGIC mean predicted =",
  mean(
    comparison_data$risk_maggic_1y
  ),
  "\n"
)

cat(
  "MAGGIC modernized mean predicted =",
  mean(
    comparison_data$risk_maggic_modern_1y
  ),
  "\n"
)


cat(
  "\nBase6 O/P =",
  observed_common /
    mean(
      comparison_data$risk_1y
    ),
  "\n"
)

cat(
  "MAGGIC original O/P =",
  observed_common /
    mean(
      comparison_data$risk_maggic_1y
    ),
  "\n"
)

cat(
  "MAGGIC modernized O/P =",
  observed_common /
    mean(
      comparison_data$risk_maggic_modern_1y
    ),
  "\n"
)


# ============================================================
# 26. Final preliminary comparison
# ============================================================

comparison_summary <- data.frame(

  model = c(
    "Base6",
    "MAGGIC_original",
    "MAGGIC_modernized"
  ),

  n =
    nrow(
      comparison_data
    ),

  AUC_364 = c(
    auc_base6,
    auc_maggic,
    auc_maggic_modern
  ),

  Brier_365 = c(
    brier_base6,
    brier_maggic,
    brier_maggic_modern
  ),

  mean_predicted = c(

    mean(
      comparison_data$risk_1y
    ),

    mean(
      comparison_data$risk_maggic_1y
    ),

    mean(
      comparison_data$risk_maggic_modern_1y
    )
  ),

  observed_KM =
    observed_common,

  O_P = c(

    observed_common /
      mean(
        comparison_data$risk_1y
      ),

    observed_common /
      mean(
        comparison_data$risk_maggic_1y
      ),

    observed_common /
      mean(
        comparison_data$risk_maggic_modern_1y
      )
  )
)


cat("\n============================================\n")
cat("J. PRELIMINARY HEAD-TO-HEAD SUMMARY\n")
cat("============================================\n")

print(
  comparison_summary,
  row.names = FALSE
)


# ============================================================
# 27. Save
# ============================================================

write.csv(

  comparison_summary,

  "Chinese_Base6_MAGGIC_preliminary_comparison.csv",

  row.names = FALSE
)


write.csv(

  component_summary,

  "Chinese_MAGGIC_score_components.csv",

  row.names = FALSE
)


save(

  maggic_data,

  comparison_data,

  comparison_summary,

  component_summary,

  score_lookup,

  file =
    "Chinese_HF_STEP10B_MAGGIC_COMPARISON.RData"
)


cat("\n============================================\n")
cat("STEP 10B 完成\n")
cat("============================================\n")
