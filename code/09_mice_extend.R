# ============================================================
# STEP 4D
# 继续 MICE 迭代至50次
# ============================================================

rm(list = ls())

library(mice)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load("NHANES_2011_2018_STEP4B_MICE_REVISED.RData")

cat("============================================\n")
cat("STEP 4D: Extend MICE iterations\n")
cat("============================================\n")

cat(
  "当前迭代次数:",
  imp$iteration,
  "\n"
)

# ------------------------------------------------------------
# 继续增加30次
# 原来20次 → 总计50次
# ------------------------------------------------------------

imp50 <- mice.mids(
  imp,
  maxit = 30,
  printFlag = TRUE
)

cat(
  "\n新的迭代次数:",
  imp50$iteration,
  "\n"
)


# ============================================================
# 1. Logged events
# ============================================================

cat("\n============================================\n")
cat("A. Logged events\n")
cat("============================================\n")

if (is.null(imp50$loggedEvents)) {

  cat("No logged events.\n")

} else {

  print(
    imp50$loggedEvents
  )
}


# ============================================================
# 2. Convergence
# ============================================================

cat("\n============================================\n")
cat("B. Convergence after 50 iterations\n")
cat("============================================\n")

conv50 <- convergence(
  imp50
)

# 只显示真正被插补的变量
imp_vars <- names(
  imp50$method[
    imp50$method != ""
  ]
)

conv50_imp <- conv50[
  conv50$vrb %in% imp_vars,
]


# ------------------------------------------------------------
# 每个变量只看最后一次迭代
# ------------------------------------------------------------

conv50_final <- do.call(
  rbind,
  lapply(
    split(
      conv50_imp,
      conv50_imp$vrb
    ),
    function(x) {

      x[
        which.max(x$.it),
        ,
        drop = FALSE
      ]
    }
  )
)

rownames(conv50_final) <- NULL

print(
  conv50_final
)


# ============================================================
# 3. 标记 PSRF > 1.10
# ============================================================

cat("\n============================================\n")
cat("C. Variables with final PSRF >1.10\n")
cat("============================================\n")

bad_psrf <- conv50_final[
  !is.na(conv50_final$psrf) &
    conv50_final$psrf > 1.10,
]

if (nrow(bad_psrf) == 0) {

  cat(
    "None. All final PSRF <= 1.10\n"
  )

} else {

  print(
    bad_psrf
  )
}


# ============================================================
# 4. 最后10次迭代范围
# ============================================================

cat("\n============================================\n")
cat("D. Last 10 iterations PSRF range\n")
cat("============================================\n")

last10_summary <- do.call(
  rbind,
  lapply(
    split(
      conv50_imp,
      conv50_imp$vrb
    ),
    function(x) {

      x <- x[
        order(
          x$.it
        ),
      ]

      z <- tail(
        x,
        10
      )

      data.frame(

        variable =
          z$vrb[1],

        min_psrf =
          min(
            z$psrf,
            na.rm = TRUE
          ),

        median_psrf =
          median(
            z$psrf,
            na.rm = TRUE
          ),

        max_psrf =
          max(
            z$psrf,
            na.rm = TRUE
          ),

        final_psrf =
          tail(
            z$psrf,
            1
          ),

        final_ac =
          tail(
            z$ac,
            1
          )
      )
    }
  )
)

rownames(last10_summary) <- NULL

print(
  last10_summary
)


# ============================================================
# 5. Trace plots
# ============================================================

cat("\n============================================\n")
cat("E. Trace plots\n")
cat("============================================\n")

for (v in imp_vars) {

  print(
    plot(
      imp50,
      v
    )
  )
}


# ============================================================
# 6. 保存
# ============================================================

write.csv(
  conv50_final,
  "MICE50_final_convergence.csv",
  row.names = FALSE
)

write.csv(
  last10_summary,
  "MICE50_last10_convergence.csv",
  row.names = FALSE
)

save(
  imp50,
  conv50,
  conv50_final,
  last10_summary,
  file =
    "NHANES_2011_2018_STEP4D_MICE50.RData"
)

cat("\n============================================\n")
cat("STEP 4D 完成\n")
cat("============================================\n")





# ============================================================
# STEP 4E
# Final MICE v3
#
# 改动：
# 1. NWR 直接作为预测变量进行 PMM 插补
# 2. 不再插补 lbdneno / lbxwbcsi
# 3. Scr仍PMM插补，eGFR随后重算
# 4. 30 imputations × 50 iterations
# ============================================================

rm(list = ls())

library(dplyr)
library(mice)
library(survival)

analysis_dir <- Sys.getenv(
  "BASE6_ANALYSIS_DIR",
  unset = getwd()
)
setwd(analysis_dir)

load("NHANES_2011_2018_STEP3B_WEIGHT_MISSING.RData")

cat("============================================\n")
cat("STEP 4E: Final MICE v3\n")
cat("============================================\n")


# ============================================================
# 1. 工作数据
# ============================================================

mi_dat <- hf_base2 %>%
  mutate(

    death =
      as.integer(
        mortstat == 1
      ),

    followup_month =
      permth_exm,

    # NWR使用原始观测值
    nwr_model =
      nwr
  )


# ============================================================
# 2. Nelson-Aalen
# ============================================================

mi_dat$nelson_aalen <-
  mice::nelsonaalen(
    mi_dat,
    timevar = "followup_month",
    statusvar = "death"
  )


# ============================================================
# 3. 类型
# ============================================================

mi_dat$cycle <- factor(
  mi_dat$cycle
)

mi_dat$female <- as.integer(
  mi_dat$female
)

mi_dat$death <- as.integer(
  mi_dat$death
)

mi_dat$aceiarb_use <- as.integer(
  mi_dat$aceiarb_use
)

mi_dat$bb_use <- as.integer(
  mi_dat$bb_use
)

mi_dat$mra_use <- as.integer(
  mi_dat$mra_use
)


# ============================================================
# 4. MICE变量
# ============================================================

mice_vars <- c(

  "seqn",

  "cycle",
  "ridageyr",
  "female",

  "sbp_mean",
  "bpxpls",
  "bmxbmi",

  "lbxhgb",
  "lbxrdw",

  # 直接插补NWR
  "nwr_model",

  # eGFR基础变量
  "scr_mgdl",

  "lbxsal",
  "lbxsbu",
  "lbxsnasi",

  "med_cnt",
  "aceiarb_use",
  "bb_use",
  "mra_use",

  "death",
  "followup_month",
  "nelson_aalen",

  "wtmec8yr",
  "sdmvstra",
  "sdmvpsu"
)


mi_input <- mi_dat %>%
  select(
    all_of(
      mice_vars
    )
  )


# ============================================================
# 5. 缺失检查
# ============================================================

cat("\n============================================\n")
cat("A. Missingness\n")
cat("============================================\n")

missing_table <- data.frame(

  variable =
    names(
      mi_input
    ),

  n_missing =
    sapply(
      mi_input,
      function(x)
        sum(
          is.na(x)
        )
    ),

  pct_missing =
    round(
      sapply(
        mi_input,
        function(x)
          mean(
            is.na(x)
          ) *
          100
      ),
      2
    )
)

print(
  missing_table[
    missing_table$n_missing > 0,
  ],
  row.names = FALSE
)


# ============================================================
# 6. 初始化
# ============================================================

ini <- mice(
  mi_input,
  maxit = 0,
  printFlag = FALSE
)

method <- ini$method
pred <- ini$predictorMatrix

method[] <- ""


# ============================================================
# 7. 被插补变量
# ============================================================

impute_vars <- c(

  "sbp_mean",
  "bpxpls",
  "bmxbmi",

  "lbxhgb",
  "lbxrdw",

  "nwr_model",

  "scr_mgdl",

  "lbxsal",
  "lbxsbu",
  "lbxsnasi"
)


method[
  impute_vars
] <- "pmm"


cat("\n============================================\n")
cat("B. Methods\n")
cat("============================================\n")

print(
  method[
    method != ""
  ]
)


# ============================================================
# 8. Predictor matrix
# ============================================================

pred[,] <- 1
diag(pred) <- 0


# SEQN不用作预测器
pred[, "seqn"] <- 0
pred["seqn", ] <- 0


# ============================================================
# 9. Survival信息
# ============================================================

# 不重复使用原始follow-up time
pred[
  impute_vars,
  "followup_month"
] <- 0


# outcome + Nelson-Aalen作为辅助预测
pred[
  impute_vars,
  "death"
] <- 1

pred[
  impute_vars,
  "nelson_aalen"
] <- 1


# ============================================================
# 10. Survey design辅助预测
# ============================================================

survey_aux <- c(
  "cycle",
  "wtmec8yr",
  "sdmvstra",
  "sdmvpsu"
)

for (v in survey_aux) {

  pred[
    impute_vars,
    v
  ] <- 1
}


# ============================================================
# 11. 正式MICE
# ============================================================

set.seed(
  20260820
)

imp_final <- mice(

  mi_input,

  m = 30,

  maxit = 50,

  method =
    method,

  predictorMatrix =
    pred,

  seed =
    20260820,

  printFlag =
    TRUE
)


# ============================================================
# 12. Logged events
# ============================================================

cat("\n============================================\n")
cat("C. Logged events\n")
cat("============================================\n")

if (
  is.null(
    imp_final$loggedEvents
  )
) {

  cat(
    "No logged events.\n"
  )

} else {

  print(
    imp_final$loggedEvents
  )
}


# ============================================================
# 13. Convergence
# ============================================================

cat("\n============================================\n")
cat("D. Final convergence\n")
cat("============================================\n")

conv_final <- convergence(
  imp_final
)


imp_vars_actual <- names(
  imp_final$method[
    imp_final$method != ""
  ]
)


conv_imp <- conv_final[
  as.character(
    conv_final$vrb
  ) %in%
    imp_vars_actual,
]


conv_last <- do.call(

  rbind,

  lapply(
    split(
      conv_imp,
      as.character(
        conv_imp$vrb
      ),
      drop = TRUE
    ),

    function(x) {

      x <- x[
        order(
          x$.it
        ),
      ]

      tail(
        x,
        1
      )
    }
  )
)


rownames(
  conv_last
) <- NULL

print(
  conv_last
)


# ============================================================
# 14. 检查30个数据集中的NWR
# ============================================================

cat("\n============================================\n")
cat("E. NWR validity across 30 imputations\n")
cat("============================================\n")


nwr_validity <- do.call(

  rbind,

  lapply(
    1:imp_final$m,

    function(i) {

      d <- complete(
        imp_final,
        i
      )

      x <-
        d$nwr_model


      data.frame(

        imputation =
          i,

        min =
          min(
            x
          ),

        max =
          max(
            x
          ),

        n_le_0 =
          sum(
            x <= 0
          ),

        n_ge_1 =
          sum(
            x >= 1
          )
      )
    }
  )
)


print(
  nwr_validity
)


cat(
  "\nTotal NWR <=0:",
  sum(
    nwr_validity$n_le_0
  ),
  "\n"
)

cat(
  "Total NWR >=1:",
  sum(
    nwr_validity$n_ge_1
  ),
  "\n"
)


# ============================================================
# 15. 第一个completed dataset构造正式衍生变量
# ============================================================

d1 <- complete(
  imp_final,
  1
)

d1 <- d1 %>%
  mutate(

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


cat("\n============================================\n")
cat("F. Derived predictors check\n")
cat("============================================\n")

cat("\nNWR:\n")
print(
  summary(
    d1$nwr
  )
)

cat("\neGFR:\n")
print(
  summary(
    d1$egfr_2021
  )
)


# ============================================================
# 16. 保存
# ============================================================

write.csv(
  conv_last,
  "MICE_FINAL_convergence.csv",
  row.names = FALSE
)

write.csv(
  nwr_validity,
  "MICE_FINAL_NWR_validity.csv",
  row.names = FALSE
)

save(
  imp_final,
  method,
  pred,
  conv_final,
  conv_last,
  nwr_validity,
  file =
    "NHANES_2011_2018_MICE_FINAL.RData"
)


cat("\n============================================\n")
cat("STEP 4E 完成\n")
cat("============================================\n")
