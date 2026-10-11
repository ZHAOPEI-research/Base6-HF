# ============================================================
# STEP 6B
# Candidate simplified model comparison
#
# Core-6
# Core-9
# Extended-10
#
# 这一步只比较 apparent fit / discrimination
# 不是最终内部验证
# ============================================================

rm(list = ls())

library(mice)
library(survival)
library(splines)
library(dplyr)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load(
  "NHANES_2011_2018_STEP6A_GROUP_LASSO_FINAL.RData"
)

load(
  "NHANES_2011_2018_MICE_FINAL.RData"
)

cat("============================================\n")
cat("STEP 6B: Candidate model comparison\n")
cat("============================================\n")


# ============================================================
# 1. 使用Step 6A已经冻结的NWR knots
# ============================================================

cat("\nNWR internal knots:\n")
print(nwr_knots)

cat("\nNWR boundary knots:\n")
print(nwr_boundary)


# ============================================================
# 2. 每个MI构造最终变量
# ============================================================

prepare_model_data <- function(d) {

  d <- d %>%
    mutate(

      age80 =
        pmin(
          ridageyr,
          80
        ),

      pulse_rate =
        bpxpls,

      nwr =
        nwr_model,

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

        (0.9938^ridageyr) *

        ifelse(
          female == 1,
          1.012,
          1
        )
    )


  # ----------------------------------------------------------
  # 固定NWR spline basis
  # ----------------------------------------------------------

  nb <- ns(

    d$nwr,

    knots =
      nwr_knots,

    Boundary.knots =
      nwr_boundary,

    intercept =
      FALSE
  )


  d$nwr_ns1 <- nb[, 1]
  d$nwr_ns2 <- nb[, 2]
  d$nwr_ns3 <- nb[, 3]


  d
}


completed_list_final <- lapply(

  1:imp_final$m,

  function(i) {

    prepare_model_data(
      complete(
        imp_final,
        i
      )
    )
  }
)


# ============================================================
# 3. 三个候选模型
# ============================================================

formulas <- list(

  # ----------------------------------------------------------
  # Core-6
  # 6 predictors
  # 实际8个coefficients
  # ----------------------------------------------------------

  Core6 =
    Surv(
      followup_month,
      death
    ) ~
    age80 +
    bmxbmi +
    lbxrdw +
    nwr_ns1 +
    nwr_ns2 +
    nwr_ns3 +
    lbxsbu +
    lbxsal,


  # ----------------------------------------------------------
  # Core-9
  # Core6 + Pulse + Na + eGFR
  # ----------------------------------------------------------

  Core9 =
    Surv(
      followup_month,
      death
    ) ~
    age80 +
    bmxbmi +
    lbxrdw +
    nwr_ns1 +
    nwr_ns2 +
    nwr_ns3 +
    lbxsbu +
    lbxsal +
    pulse_rate +
    lbxsnasi +
    egfr_2021,


  # ----------------------------------------------------------
  # Extended-10
  # Core9 + medication count
  # ----------------------------------------------------------

  Extended10 =
    Surv(
      followup_month,
      death
    ) ~
    age80 +
    bmxbmi +
    lbxrdw +
    nwr_ns1 +
    nwr_ns2 +
    nwr_ns3 +
    lbxsbu +
    lbxsal +
    pulse_rate +
    lbxsnasi +
    egfr_2021 +
    med_cnt
)


# ============================================================
# 4. 检查模型parameter数量
# ============================================================

cat("\n============================================\n")
cat("A. Candidate model dimensions\n")
cat("============================================\n")

for (model_name in names(formulas)) {

  fit_test <-
    coxph(
      formulas[[model_name]],
      data =
        completed_list_final[[1]],
      ties =
        "efron"
    )

  cat(
    model_name,
    ": ",
    length(
      coef(fit_test)
    ),
    " coefficients\n",
    sep = ""
  )
}


# ============================================================
# 5. 在30个MI分别拟合
# ============================================================

model_results <- list()

model_fits <- list()

counter <- 1


for (i in seq_along(completed_list_final)) {

  d <- completed_list_final[[i]]


  for (model_name in names(formulas)) {

    fit <- coxph(

      formulas[[model_name]],

      data =
        d,

      ties =
        "efron",

      x =
        TRUE,

      y =
        TRUE
    )


    # --------------------------------------------------------
    # Harrell's C-index
    #
    # 对 coxph object:
    # survival::concordance() 会自动处理风险方向
    # 不需要 reverse = TRUE
    # --------------------------------------------------------

    conc <- concordance(
      fit
    )

    c_index <- as.numeric(
      conc$concordance
    )


    model_results[[counter]] <-
      data.frame(

        imputation =
          i,

        model =
          model_name,

        n =
          fit$n,

        events =
          fit$nevent,

        parameters =
          length(
            coef(
              fit
            )
          ),

        logLik =
          as.numeric(
            logLik(
              fit
            )
          ),

        AIC =
          AIC(
            fit
          ),

        C_index =
          c_index
      )


    model_fits[[
      paste0(
        model_name,
        "_MI",
        i
      )
    ]] <- fit


    counter <-
      counter + 1
  }
}


model_results <-
  bind_rows(
    model_results
  )


# ============================================================
# 6. 每个MI内相对最佳模型计算delta AIC
# ============================================================

model_results <-
  model_results %>%
  group_by(
    imputation
  ) %>%
  mutate(

    min_AIC =
      min(
        AIC
      ),

    delta_AIC =
      AIC -
      min_AIC,

    max_C =
      max(
        C_index
      ),

    delta_C =
      max_C -
      C_index
  ) %>%
  ungroup()


# ============================================================
# 7. 汇总
# ============================================================

model_summary <-
  model_results %>%
  group_by(
    model
  ) %>%
  summarise(

    parameters =
      first(
        parameters
      ),

    median_AIC =
      median(
        AIC
      ),

    median_delta_AIC =
      median(
        delta_AIC
      ),

    p25_delta_AIC =
      quantile(
        delta_AIC,
        0.25
      ),

    p75_delta_AIC =
      quantile(
        delta_AIC,
        0.75
      ),

    median_C =
      median(
        C_index
      ),

    min_C =
      min(
        C_index
      ),

    max_C =
      max(
        C_index
      ),

    median_delta_C =
      median(
        delta_C
      ),

    .groups =
      "drop"
  ) %>%
  arrange(
    median_delta_AIC
  )


cat("\n============================================\n")
cat("B. Candidate-model apparent comparison\n")
cat("============================================\n")

print(
  model_summary,
  n = Inf
)


# ============================================================
# 8. Core6 vs Core9
# ============================================================

core6_vs_9 <- model_results %>%
  select(
    imputation,
    model,
    AIC,
    C_index
  ) %>%
  tidyr::pivot_wider(

    names_from =
      model,

    values_from =
      c(
        AIC,
        C_index
      )
  ) %>%
  mutate(

    AIC_gain_Core9 =
      AIC_Core6 -
      AIC_Core9,

    C_gain_Core9 =
      C_index_Core9 -
      C_index_Core6,

    AIC_gain_Extended10 =
      AIC_Core6 -
      AIC_Extended10,

    C_gain_Extended10 =
      C_index_Extended10 -
      C_index_Core6
  )


cat("\n============================================\n")
cat("C. Added value over Core-6\n")
cat("============================================\n")

cat("\nCore9 vs Core6:\n")

print(
  summary(
    core6_vs_9$
      AIC_gain_Core9
  )
)

cat("\nC-index gain:\n")

print(
  summary(
    core6_vs_9$
      C_gain_Core9
  )
)


cat("\nExtended10 vs Core6:\n")

print(
  summary(
    core6_vs_9$
      AIC_gain_Extended10
  )
)

cat("\nC-index gain:\n")

print(
  summary(
    core6_vs_9$
      C_gain_Extended10
  )
)


# ============================================================
# 9. Medication count added value
#
# Extended10 vs Core9
# ============================================================

core6_vs_9 <- core6_vs_9 %>%
  mutate(

    AIC_gain_medication =
      AIC_Core9 -
      AIC_Extended10,

    C_gain_medication =
      C_index_Extended10 -
      C_index_Core9
  )


cat("\n============================================\n")
cat("D. Medication-count added value\n")
cat("============================================\n")

cat("\nAIC improvement from Medication_count:\n")

print(
  summary(
    core6_vs_9$
      AIC_gain_medication
  )
)


cat("\nC-index improvement from Medication_count:\n")

print(
  summary(
    core6_vs_9$
      C_gain_medication
  )
)


# ============================================================
# 10. 保存
# ============================================================

write.csv(
  model_results,
  "HF_candidate_models_all_MI.csv",
  row.names = FALSE
)

write.csv(
  model_summary,
  "HF_candidate_models_summary.csv",
  row.names = FALSE
)

write.csv(
  core6_vs_9,
  "HF_candidate_models_added_value.csv",
  row.names = FALSE
)


save(
  completed_list_final,
  formulas,
  model_results,
  model_summary,
  core6_vs_9,
  model_fits,
  nwr_knots,
  nwr_boundary,
  file =
    "NHANES_2011_2018_STEP6B_CANDIDATE_MODELS.RData"
)


cat("\n============================================\n")
cat("STEP 6B 完成\n")
cat("============================================\n")
